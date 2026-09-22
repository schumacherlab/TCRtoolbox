import collections
import json
import subprocess
from pathlib import Path
from datetime import datetime
from typing import TextIO

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from tcr_toolbox.utils.plot_utils import startfig


def run_mixmhcpred(
    patient_df: pd.DataFrame,
    patient_hla_df: pd.DataFrame,
    patient: str,
    run_dir: Path,
    kmer_list: list,
    gene_id_col: str,
    aa_seq_col: str,
    hla_df_patient_col: str,
    mixmhcpred_exe: str,
    allele_1_col: str = "allele_1_mix_enc",
    allele_2_col: str = "allele_2_mix_enc",
    variant: bool = False,
    minigene_aa_var_start_col: str = None,
    minigene_aa_var_end_col: str = None,
    var_type_for_kmer_overlap_col: str = None,
):
    """Generate the minimal set of variant-overlapping minimal peptides for one patient/donor and predict their binding confidence to the patient's HLA class I alleles with MixMHCpred3.0.

    `kmer_list`-length peptides are tiled across each antigen's amino acid
    sequence; when `variant` is False, all such kmers are kept. When `variant`
    is True, minimal peptides are defined as containing at least one mutated
    residue, for SNVs, in-frame insertions, and frameshifting indels; for
    in-frame deletions, minimal peptides are instead required to contain both
    amino acids forming the neojunction created by the deletion.

    Because the same peptide sequence can come from more than one minigene
    (which may have different mRNA expression), identical peptides are
    predicted only once; their mapping back to every source minigene is saved
    to JSON files in `run_dir/mixmhcpred/`, so `parse_mixmhcpred_out_txt` can
    later attribute expression correctly to each duplicate. The deduplicated
    peptides are written to a fasta, and their binding confidence (%Rank) to
    all of the patient's HLA class I alleles (from `patient_hla_df`) is
    predicted using MixMHCpred3.0 (`mixmhcpred_exe`); the prediction output
    and logs are saved under `run_dir`.

    Parameters
    ----------
    patient_df : pd.DataFrame
        One row per antigen for this patient, containing at least
        `gene_id_col` and `aa_seq_col` (and, if `variant` is True,
        `minigene_aa_var_start_col`, `minigene_aa_var_end_col`, and
        `var_type_for_kmer_overlap_col`).
    patient_hla_df : pd.DataFrame
        This is an example of a correctly encoded `patient_hla_df` for two
        patients:

        | patient     | gene  | allele_1_mix_enc | allele_2_mix_enc |
        |-------------|-------|------------------|------------------|
        | patient01   | HLA-A | A0201            | A2402            |
        | patient01   | HLA-B | B0702            | B4901            |
        | patient01   | HLA-C | C0727            | C0763            |
        | patient02   | HLA-A | A0201            | A2601            |
        | patient02   | HLA-B | B1553            | B4037            |
        | patient02   | HLA-C | C0304            | C0304            |
    patient : str
        This patient's identifier; must match a value in
        `patient_hla_df[hla_df_patient_col]`, and is used to name output files.
    run_dir : Path
        Directory under which a `mixmhcpred/` subdirectory (with the query
        fasta, logs, and MixMHCpred output) is created.
    kmer_list : list
        Peptide lengths to tile across each antigen sequence (e.g. `[8, 9, 10,
        11]` for the minimal-peptide 8-11-mer definition).
    gene_id_col : str
        Column in `patient_df` holding each antigen's name/gene ID, used as
        the kmer name prefix.
    aa_seq_col : str
        Column in `patient_df` holding each antigen's full amino acid
        sequence.
    hla_df_patient_col : str
        Column in `patient_hla_df` identifying the patient each row's HLA
        alleles belong to. 'patient' in the patient_hla_df argument table
        example.
    mixmhcpred_exe : str
        Path to the MixMHCpred3.0 executable.
    allele_1_col : str, default = "allele_1_mix_enc"
        Column in `patient_hla_df` holding the patient's first HLA class I
        allele.
    allele_2_col : str, default = "allele_2_mix_enc"
        Column in `patient_hla_df` holding the patient's second HLA class I
        allele.
    variant : bool, default = False
        If True, restrict kmers per antigen to those overlapping its annotated
        variant (see above); requires `minigene_aa_var_start_col`,
        `minigene_aa_var_end_col`, and `var_type_for_kmer_overlap_col` to be
        set.
    minigene_aa_var_start_col : str, default = None
        Column in `patient_df` holding each antigen's variant start position
        (see `var_start_end_for_kmer_overlap`); required when `variant` is
        True.
    minigene_aa_var_end_col : str, default = None
        Column in `patient_df` holding each antigen's variant end position;
        required when `variant` is True.
    var_type_for_kmer_overlap_col : str, default = None
        Column in `patient_df` holding each antigen's variant type
        (`"inframe_del"` or `"missense_or_inframe_ins_or_frameshift"`);
        required when `variant` is True.

    Returns
    -------
    Path
        Path to the MixMHCpred output `.txt` file, named
        `<min(kmer_list)>-<max(kmer_list)>_kmers_<patient>_mixmhcpred.txt`.

    Notes
    -----
    Each kmer written to the intermediate FASTA (and therefore each row of
    the eventual MixMHCpred output) is named `<name>_<k>-<i>`, where `name`
    is that row's `gene_id_col` value (e.g. the minigene/antigen name),
    `k` is the kmer length (one of `kmer_list`), and `i` is the kmer's
    0-based start position within that antigen's amino acid sequence
    (`aa_seq_col`). `parse_mixmhcpred_out_txt` parses this name back apart
    to recover `gene_id_col` and re-attribute duplicate kmers to every
    source minigene that produced them.

    Raises
    ------
    ValueError
        If `gene_id_col` or `aa_seq_col` contains any NaN values for this
        patient.
    Exception
        If `variant` is True and an antigen's `var_type_for_kmer_overlap_col`
        is neither `"inframe_del"` nor
        `"missense_or_inframe_ins_or_frameshift"`.
    """
    patient_df = patient_df.copy()
    if patient_df[gene_id_col].isna().any():
        raise ValueError(f"NaN values found in {gene_id_col} for patient {patient}")

    if patient_df[aa_seq_col].isna().any():
        raise ValueError(f"NaN values found in {aa_seq_col} for patient {patient}")

    mixmhcpred_dir = run_dir / "mixmhcpred"
    mixmhcpred_dir.mkdir(exist_ok=True)
    fasta_file = mixmhcpred_dir / f"{min(kmer_list)}-{max(kmer_list)}_kmers_{patient}.fasta"
    patient_df.loc[:, f"{aa_seq_col}_len"] = patient_df.loc[:, aa_seq_col].str.len()

    kmer_duplicate_name_dict: dict = collections.defaultdict(list)
    kmer_to_first_duplicate_name: dict = dict()

    with open(fasta_file, "w") as fasta_out:
        log_dir = fasta_file.parent / "logs"
        log_dir.mkdir(exist_ok=True)
        out_txt_file = fasta_file.with_suffix(".txt").with_stem(fasta_file.stem + "_mixmhcpred")
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = log_dir / f"{fasta_file.stem}_{timestamp}.log"

        kmer_cache: set = set()

        if variant:
            for name, aa_seq, aa_seq_len, var_start, var_end, var_type in zip(
                patient_df[gene_id_col],
                patient_df[aa_seq_col],
                patient_df[f"{aa_seq_col}_len"],
                patient_df[minigene_aa_var_start_col],
                patient_df[minigene_aa_var_end_col],
                patient_df[var_type_for_kmer_overlap_col],
            ):
                for k in kmer_list:
                    if aa_seq_len < k:
                        print(f"Warning: sequence {name} is shorter ({aa_seq_len}) than kmer length {k}. Skipping kmer generation.")
                        continue
                    for i in np.arange(0, aa_seq_len - k + 1):
                        kmer = aa_seq[i : i + k]
                        kmer_start = i
                        kmer_end = i + k - 1

                        if var_type == "inframe_del":
                            if not (kmer_start <= var_start and var_end <= kmer_end):
                                continue
                        elif var_type == "missense_or_inframe_ins_or_frameshift":
                            if not (kmer_start < var_end and var_start <= kmer_end):
                                continue
                        else:
                            raise Exception(f"{var_type} is not supported! Use: inframe_del or missense_or_inframe_ins_or_frameshift")

                        if kmer in kmer_cache:
                            kmer_duplicate_name_dict[kmer_to_first_duplicate_name[kmer]].append(f"{name}_{k}-{i}")
                            continue

                        print(f">{name}_{k}-{i}", file=fasta_out)
                        print(str(kmer), file=fasta_out)
                        kmer_to_first_duplicate_name[kmer] = f"{name}_{k}-{i}"
                        kmer_cache.add(kmer)

        else:
            for name, aa_seq, aa_seq_len in zip(patient_df[gene_id_col], patient_df[aa_seq_col], patient_df[f"{aa_seq_col}_len"]):
                for k in kmer_list:
                    if aa_seq_len < k:
                        print(f"Warning: sequence {name} is shorter ({aa_seq_len}) than kmer length {k}. Skipping kmer generation.")
                        continue
                    for i in np.arange(0, aa_seq_len - k + 1):
                        kmer = aa_seq[i : i + k]
                        if kmer in kmer_cache:
                            kmer_duplicate_name_dict[kmer_to_first_duplicate_name[kmer]].append(f"{name}_{k}-{i}")
                            continue

                        print(f">{name}_{k}-{i}", file=fasta_out)
                        print(str(kmer), file=fasta_out)
                        kmer_to_first_duplicate_name[kmer] = f"{name}_{k}-{i}"
                        kmer_cache.add(kmer)

    kmer_json_file = mixmhcpred_dir / f"{patient}_duplicate_kmer_name_dict.json"
    with open(kmer_json_file, "w") as json_out:
        json.dump(dict(kmer_duplicate_name_dict), json_out, indent=2)

    kmer_json_file = mixmhcpred_dir / f"{patient}_kmer_to_first_duplicate_name.json"
    with open(kmer_json_file, "w") as json_out:
        json.dump(dict(kmer_to_first_duplicate_name), json_out, indent=2)

    hla_list = (
        patient_hla_df.loc[patient_hla_df[hla_df_patient_col] == patient, allele_1_col].dropna().tolist()
        + patient_hla_df.loc[patient_hla_df[hla_df_patient_col] == patient, allele_2_col].dropna().tolist()
    )
    hla_str = ",".join(hla_list)
    print(hla_str)

    with open(log_file, "w") as f:
        result = subprocess.run([mixmhcpred_exe, "-a", hla_str, "-i", fasta_file, "-o", out_txt_file], stdout=f, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        raise RuntimeError(
            f"MixMHCpred3.0 failed for patient {patient} (exit code {result.returncode}). See log file for details: {log_file}"
        )

    return out_txt_file


def run_mixmhcpred_patient_df_dict(
    patient_df_dict: dict,
    patient_hla_df: pd.DataFrame,
    run_dir: Path,
    kmer_list: list,
    gene_id_col: str,
    aa_seq_col: str,
    hla_df_patient_col: str,
    mixmhcpred_exe: str,
    allele_1_col: str = "allele_1_mix_enc",
    allele_2_col: str = "allele_2_mix_enc",
    variant: bool = False,
    minigene_aa_var_start_col: str = None,
    minigene_aa_var_end_col: str = None,
    var_type_for_kmer_overlap_col: str = None,
) -> list:
    """Run `run_mixmhcpred` for every patient in `patient_df_dict`.

    Parameters
    ----------
    patient_df_dict : dict
        Mapping of patient identifier to that patient's antigen dataframe (see
        `run_mixmhcpred`'s `patient_df`); each is processed independently.
    patient_hla_df : pd.DataFrame
        HLA typing table shared across all patients (see `run_mixmhcpred`).
    run_dir : Path
        Directory under which each patient's `mixmhcpred/` output is created.
    kmer_list : list
        Peptide lengths to tile across each antigen sequence, passed through
        to `run_mixmhcpred` for every patient.
    gene_id_col : str
        Column holding each antigen's name/gene ID, passed through for every
        patient.
    aa_seq_col : str
        Column holding each antigen's amino acid sequence, passed through for
        every patient.
    hla_df_patient_col : str
        Column in `patient_hla_df` identifying which patient each row's
        alleles belong to.
    mixmhcpred_exe : str
        Path to the MixMHCpred3.0 executable.
    allele_1_col : str, default = "allele_1_mix_enc"
        Column in `patient_hla_df` holding each patient's first HLA class I
        allele.
    allele_2_col : str, default = "allele_2_mix_enc"
        Column in `patient_hla_df` holding each patient's second HLA class I
        allele.
    variant : bool, default = False
        Passed through to `run_mixmhcpred`; if True, restrict kmers per
        antigen to those overlapping its annotated variant.
    minigene_aa_var_start_col : str, default = None
        Passed through to `run_mixmhcpred`; required when `variant` is True.
    minigene_aa_var_end_col : str, default = None
        Passed through to `run_mixmhcpred`; required when `variant` is True.
    var_type_for_kmer_overlap_col : str, default = None
        Passed through to `run_mixmhcpred`; required when `variant` is True.

    Returns
    -------
    list[Path]
        The MixMHCpred output `.txt` path returned by `run_mixmhcpred` for
        each patient, in `patient_df_dict` iteration order.

    Raises
    ------
    ValueError
        Propagated from `run_mixmhcpred` if any patient's `gene_id_col` or
        `aa_seq_col` contains NaN values.
    Exception
        Propagated from `run_mixmhcpred` if `variant` is True and an antigen's
        `var_type_for_kmer_overlap_col` is not a supported variant type.
    """
    mixmhcpred_out_txt_list = []
    patient_df_dict = patient_df_dict.copy()
    patient_hla_df = patient_hla_df.copy()

    for patient, patient_df in patient_df_dict.items():
        print(f"Processing: {patient}")
        out_txt_file = run_mixmhcpred(
            patient_df=patient_df,
            patient_hla_df=patient_hla_df,
            patient=patient,
            run_dir=run_dir,
            kmer_list=kmer_list,
            gene_id_col=gene_id_col,
            aa_seq_col=aa_seq_col,
            hla_df_patient_col=hla_df_patient_col,
            allele_1_col=allele_1_col,
            allele_2_col=allele_2_col,
            mixmhcpred_exe=mixmhcpred_exe,
            variant=variant,
            minigene_aa_var_start_col=minigene_aa_var_start_col,
            minigene_aa_var_end_col=minigene_aa_var_end_col,
            var_type_for_kmer_overlap_col=var_type_for_kmer_overlap_col,
        )

        mixmhcpred_out_txt_list.append(out_txt_file)

    return mixmhcpred_out_txt_list


def parse_mixmhcpred_out_txt(
    mixmhcpred_out_txt: Path,
    patient_df: pd.DataFrame,
    patient: str,
    gene_id_col: str,
    mrna_expression_col: str,
    keep_top_number_peptides: int,
    log_file: TextIO,
    include_cols_patient_df_list: list = None,
    rank_by_binning: bool = False,
    mrna_sorting_weight: float = 0.2,
    percent_rank_threshold: float = 2.0,
    hit_col: str = None,
    confident_binder_percent_rank_threshold: float = None,
) -> pd.DataFrame:
    """Rank a patient's MixMHCpred predictions and select the top minimal peptides for screening.

    MixMHCpred predictions are first mapped back to every source minigene,
    re-expanding any peptide `run_mixmhcpred` deduplicated. Peptides with
    %Rank_bestAllele at or above `percent_rank_threshold` (2% by default) are
    then filtered out to reduce the minimal-peptide sorting set size. Protein
    mRNA expression (`mrna_expression_col`) and %Rank_bestAllele are converted
    to percentile ranks, where 1.0 represents the highest expression and
    lowest (best) %Rank_bestAllele. A combined sort score then weights the
    mRNA percentile rank by `mrna_sorting_weight` and the %Rank_bestAllele
    percentile rank by the remainder, `1 - mrna_sorting_weight` (20%/80% with
    the default `mrna_sorting_weight` of 0.2), and minimal peptides are
    sorted by this score. Duplicate minimal peptide amino acid
    sequences are then removed, retaining the first-ranked sequence. The top
    `keep_top_number_peptides` minimal peptides, as based on weighted
    predicted MHC class I binding and RNA expression, are selected for
    screening.

    Parameters
    ----------
    mixmhcpred_out_txt : Path
        Path to the MixMHCpred output `.txt` file for `patient` (as returned
        by `run_mixmhcpred`); its directory must also contain the two extra
        JSON files `run_mixmhcpred` saved alongside it,
        `<patient>_duplicate_kmer_name_dict.json` and
        `<patient>_kmer_to_first_duplicate_name.json`.
    patient_df : pd.DataFrame
        This patient's antigen dataframe (as passed to `run_mixmhcpred`), used
        to look up each peptide's `gene_id_col`, `mrna_expression_col`,
        `hit_col`, and `include_cols_patient_df_list` values.
    patient : str
        This patient's identifier; used in output/log file names and plot
        titles.
    gene_id_col : str
        Column (in both `patient_df` and, derived from kmer names, the
        MixMHCpred output) identifying each peptide's source antigen.
    mrna_expression_col : str
        Column in `patient_df` holding each antigen's mRNA expression value
        (e.g. TPM); looked up for each peptide and converted to a percentile
        rank.
    keep_top_number_peptides : int
        Number of top-sort-scoring peptides to keep for `selected_peptides_df`.
    log_file : TextIO
        Open file handle that diagnostic messages are printed to.
    include_cols_patient_df_list : list, default = None
        Extra `patient_df` meta columns to add to the output alongside
        `gene_id_col`/`mrna_expression_col`/`hit_col`.
    rank_by_binning : bool, default = False
        If True, percentile-rank `%Rank_bestAllele` and `mrna_expression_col`
        by quantile-binning into `min(n_unique_%Rank, n_unique_mRNA)` bins
        first, instead of ranking every value directly.
    mrna_sorting_weight : float, default = 0.2
        Weight given to the mRNA-expression percentile rank in the combined
        sort score; `%Rank_bestAllele` gets weight `1 - mrna_sorting_weight`.
    percent_rank_threshold : float, default = 2.0
        `%Rank_bestAllele` cutoff (exclusive); peptides at or above this are
        dropped before ranking/selection.
    hit_col : str, default = None
        Optional boolean column in `patient_df` flagging known "hit" recognized
        peptides, used only for the extra diagnostics described above; if None,
        those diagnostics are skipped.
    confident_binder_percent_rank_threshold : float, default = None
        `%Rank_bestAllele` cutoff for calling a `hit_col` peptide a confident
        binder. Only needed when `hit_col` is given; if left at the default
        None, one unrelated log message that always runs will just report a
        count of 0 instead of a real number, which is harmless.

    Returns
    -------
    pd.DataFrame
        The top `keep_top_number_peptides` rows by `sort_score`, reindexed
        from 0, with `%Rank_bestAllele_scaled`, `{mrna_expression_col}_scaled`,
        and `sort_score` columns added.
    """
    patient_df = patient_df.copy()

    include_patient_df_col_list = []
    if include_cols_patient_df_list:
        for col in include_cols_patient_df_list:
            include_patient_df_col_list.append(col)

    mixmhcpred_out_txt = Path(mixmhcpred_out_txt)

    json_file = mixmhcpred_out_txt.parent / f"{patient}_duplicate_kmer_name_dict.json"
    with open(json_file, "r") as json_in:
        kmer_duplicate_name_dict = json.load(json_in)

    json_file = mixmhcpred_out_txt.parent / f"{patient}_kmer_to_first_duplicate_name.json"
    with open(json_file, "r") as json_in:
        kmer_to_first_duplicate_name = json.load(json_in)

    mix_df = pd.read_csv(mixmhcpred_out_txt, delimiter="\t", header=11)
    mix_df.loc[:, "kmer_name"] = mix_df.loc[:, "Peptide"].map(kmer_to_first_duplicate_name)

    print(f"Any kmer_name is missing after mapping kmer aa sequence to kmer_name: {mix_df.loc[:, 'kmer_name'].isna().any()}", file=log_file, flush=True)
    print(f"Number of predicted peptides: {mix_df.shape[0]}", file=log_file)

    mix_df["dup_list"] = mix_df["kmer_name"].map(kmer_duplicate_name_dict)
    base = mix_df.loc[mix_df["dup_list"].isna()].drop(columns="dup_list")
    base_with_dup = mix_df.loc[mix_df["dup_list"].notna()].drop(columns="dup_list")

    expanded = mix_df.loc[mix_df["dup_list"].notna()].explode("dup_list")

    expanded["kmer_name"] = expanded["dup_list"]
    expanded.drop(columns="dup_list", inplace=True)

    mix_df = pd.concat([base, base_with_dup, expanded], ignore_index=True)
    print(f"Number of predicted kmers after adding duplicate kmer names: {mix_df.shape[0]}", file=log_file)

    mix_df.loc[:, gene_id_col] = mix_df.loc[:, "kmer_name"].str.split(r"_(\d+)-(\d+)", expand=True)[0]

    include_patient_df_col_list.append(gene_id_col)
    include_patient_df_col_list.append(mrna_expression_col)

    if hit_col:
        include_patient_df_col_list.append(hit_col)

    mix_df = mix_df.loc[mix_df["%Rank_bestAllele"] < percent_rank_threshold, :]
    merge_df = pd.merge(mix_df, patient_df.loc[:, include_patient_df_col_list], on=gene_id_col, how="left")

    print(f"Number of predicted kmers below {percent_rank_threshold} percent_rank_threshold: {merge_df.shape[0]}", file=log_file)
    print(
        f"Number of predicted kmers below {confident_binder_percent_rank_threshold} percent_rank_threshold: {merge_df.loc[merge_df['%Rank_bestAllele'] < confident_binder_percent_rank_threshold, :].shape[0]}",
        file=log_file,
    )

    print(f"Number of unique %Rank_bestAllele values: {len(merge_df['%Rank_bestAllele'].unique())}", file=log_file)
    print(f"Number of unique {mrna_expression_col} values: {len(merge_df[f'{mrna_expression_col}'].unique())}", file=log_file)
    if rank_by_binning:
        number_of_bins = min([len(merge_df[f"{mrna_expression_col}"].unique()), len(merge_df["%Rank_bestAllele"].unique())])
        print(f"Number of bins: {number_of_bins}", file=log_file)
        merge_df[f"{mrna_expression_col}_binned"] = pd.qcut(merge_df[mrna_expression_col], q=number_of_bins, duplicates="drop", labels=False)
        merge_df["%Rank_bestAllele_binned"] = pd.qcut(-merge_df["%Rank_bestAllele"], q=number_of_bins, duplicates="drop", labels=False)
        merge_df["%Rank_bestAllele_scaled"] = merge_df["%Rank_bestAllele_binned"] / merge_df["%Rank_bestAllele_binned"].max()
        merge_df[f"{mrna_expression_col}_scaled"] = merge_df[f"{mrna_expression_col}_binned"] / merge_df[f"{mrna_expression_col}_binned"].max()

    if not rank_by_binning:
        merge_df["%Rank_bestAllele_scaled"] = merge_df["%Rank_bestAllele"].rank(pct=True, method="max", ascending=False)
        merge_df[f"{mrna_expression_col}_scaled"] = merge_df[mrna_expression_col].rank(pct=True, method="max")

    print(f"Number of unique ranks %Rank_bestAllele_scaled: {len(merge_df['%Rank_bestAllele_scaled'].unique())}", file=log_file)
    print(f"Number of unique ranks {mrna_expression_col}_scaled: {len(merge_df[f'{mrna_expression_col}_scaled'].unique())}", file=log_file)

    merge_df.loc[:, "sort_score"] = (1 - mrna_sorting_weight) * merge_df["%Rank_bestAllele_scaled"] + mrna_sorting_weight * merge_df[f"{mrna_expression_col}_scaled"]

    merge_df.sort_values("sort_score", inplace=True, ascending=False)
    merge_df.drop_duplicates(
        "Peptide", keep="first", inplace=True
    )  # This will cause some underestimation of the number of recovered minigene hits because you cannot match on peptide sequence
    merge_df_before_selection = merge_df.copy()
    merge_df_before_selection.to_csv(mixmhcpred_out_txt.parent / f"{patient}_parsed_and_sorted_peptides_df.csv", index=False)

    ax, fig, gs = startfig(5, 5)
    ax.hist(merge_df.loc[:, "%Rank_bestAllele_scaled"], np.arange(0, 1.05, 0.05))
    ax.set_xlabel("%Rank_bestAllele_scaled", fontsize=7)
    ax.set_ylabel("Frequency", fontsize=7)
    ax.set_title(patient, fontsize=7)
    ax.tick_params("both", labelsize=7)
    ax.set_xlim(0, 1.0)
    fig.tight_layout()
    fig.savefig(mixmhcpred_out_txt.parent / "outs" / f"{patient}_percent_rank_mixmhcpred_filtered_percent_rank_bestAllele_scaled_hist.pdf")
    plt.close()

    ax, fig, gs = startfig(5, 5)
    ax.hist(merge_df.loc[:, f"{mrna_expression_col}_scaled"], np.arange(0, 1.05, 0.05))
    ax.set_xlabel(f"{mrna_expression_col}_scaled", fontsize=7)
    ax.set_ylabel("Frequency", fontsize=7)
    ax.set_title(patient, fontsize=7)
    ax.tick_params("both", labelsize=7)
    ax.set_xlim(0, 1.0)
    fig.tight_layout()
    fig.savefig(mixmhcpred_out_txt.parent / "outs" / f"{patient}_percent_rank_mixmhcpred_filtered_{mrna_expression_col}_scaled_hist.pdf")
    plt.close()

    print(f"Number of kmers after dropping Peptide aa sequence duplicates: {merge_df.shape[0]}", file=log_file)

    selected_peptides_df = merge_df.iloc[:keep_top_number_peptides].copy()
    selected_peptides_df.reset_index(inplace=True, drop=True)

    print("", file=log_file)
    print(f"Spearmanr sort_score with mrna_expression_col_scaled: {round(spearmanr(merge_df['sort_score'], merge_df[f'{mrna_expression_col}_scaled'])[0], 2)}", file=log_file)
    print(f"Spearmanr sort_score with %Rank_bestAllele_scaled: {round(spearmanr(merge_df['sort_score'], merge_df['%Rank_bestAllele_scaled'])[0], 2)}", file=log_file)
    print(
        f"Spearmanr sort_score with mrna_expression_col_scaled for "
        f"{keep_top_number_peptides} selected peptides: "
        f"{round(spearmanr(selected_peptides_df['sort_score'], selected_peptides_df[f'{mrna_expression_col}_scaled'])[0], 2)}",
        file=log_file,
    )
    print(
        f"Spearmanr sort_score with %Rank_bestAllele_scaled for "
        f"{keep_top_number_peptides} selected peptides: "
        f"{round(spearmanr(selected_peptides_df['sort_score'], selected_peptides_df['%Rank_bestAllele_scaled'])[0], 2)}",
        file=log_file,
    )

    frac_unique_muts = len(set(patient_df[gene_id_col].unique()).intersection(set(selected_peptides_df[gene_id_col].unique()))) / len(patient_df[gene_id_col].unique())
    print(f"Maximal %rankBestAllele (higher means less confident): {selected_peptides_df.loc[:, '%Rank_bestAllele'].max()}", file=log_file)

    if merge_df["kmer_name"].duplicated().any():
        raise Exception("Duplicate kmer_names should not exist!", merge_df["kmer_name"].value_counts())
    outside_percent_rank_series = merge_df.loc[~merge_df["kmer_name"].isin(selected_peptides_df["kmer_name"]), "%Rank_bestAllele"]
    if outside_percent_rank_series.empty:
        print(f"Minimal %rankBestAllele outside selection: NaN, all peptides below {percent_rank_threshold}% rank were selected", file=log_file)
    else:
        min_rank_idx = outside_percent_rank_series.idxmin()
        print(
            f"Minimal %rankBestAllele outside selection: {merge_df.loc[min_rank_idx, '%Rank_bestAllele']}, "
            f"{mrna_expression_col}: {merge_df.loc[min_rank_idx, mrna_expression_col]}, sort_score: {merge_df.loc[min_rank_idx, 'sort_score']},\n lowest sort_score selection: {selected_peptides_df['sort_score'].min()}",
            file=log_file,
        )

    print(f"Percentage of unique {gene_id_col} names in {selected_peptides_df.shape[0]} selected peptides set: {round(frac_unique_muts * 100, 1)}%", file=log_file)
    print(f"Top-30 {gene_id_col} counts: {selected_peptides_df[gene_id_col].value_counts()[:30]}", file=log_file)

    if hit_col:
        print(patient, file=log_file)
        hit_df = merge_df.loc[merge_df[hit_col], :].copy()

        for nth_best_binder_gene_id_col in range(1, 5):
            n_smallest_per_gene = hit_df.groupby(gene_id_col)["%Rank_bestAllele"].apply(
                lambda x: x.nsmallest(nth_best_binder_gene_id_col).iloc[-1] if len(x) >= nth_best_binder_gene_id_col else np.nan
            )

            print(f"top-{nth_best_binder_gene_id_col} %rankBestAllele peptide per minigene:", file=log_file)
            for gene_id, value in n_smallest_per_gene.items():
                print(f"{gene_id}: {value}", file=log_file)

        hit_df = hit_df.loc[hit_df["%Rank_bestAllele"] < confident_binder_percent_rank_threshold, :].copy()
        hit_df.loc[:, "intersect_col"] = hit_df["Peptide"] + "_" + hit_df[gene_id_col]

        print(f"Number of minigene hit derived peptides below {confident_binder_percent_rank_threshold} %rankBestAllele: {hit_df.shape[0]}", file=log_file)
        print(f"Number of unique minigene hit names below {confident_binder_percent_rank_threshold} %rankBestAllele: {len(hit_df[gene_id_col].unique())}", file=log_file)

        selected_peptides_df.loc[:, "intersect_col"] = selected_peptides_df["Peptide"] + "_" + selected_peptides_df[gene_id_col]
        print(
            f"Number of minigene hit derived peptides below {confident_binder_percent_rank_threshold} %rankBestAllele that are selected: {len(set(hit_df['intersect_col'].unique()).intersection(set(selected_peptides_df['intersect_col'].unique())))}",
            file=log_file,
        )
        if hit_df.empty:
            raise ValueError(
                f"No confidently predicted binders (below {confident_binder_percent_rank_threshold} %rankBestAllele) found for patient {patient}."
            )
        frac_peptides_below_threshold_in_all_minigene_hits = len(set(hit_df["intersect_col"].unique()).intersection(set(selected_peptides_df["intersect_col"].unique()))) / len(
            set(hit_df["intersect_col"].unique())
        )
        print(f"Fraction of peptides below {confident_binder_percent_rank_threshold} in all minigene hits recovered: {frac_peptides_below_threshold_in_all_minigene_hits}", file=log_file)
        selected_peptides_df.drop("intersect_col", axis=1, inplace=True)
        print("\n")

    # Compare on kmer_name (unique per row, asserted above) rather than full-row tuples:
    # a tuple containing NaN (e.g. a missing mRNA expression value) never equals itself,
    # so a selected row with any NaN column would otherwise be miscolored as unselected.
    selected_kmer_names = set(selected_peptides_df["kmer_name"])
    color_list = ["green" if kmer_name in selected_kmer_names else "lightgrey" for kmer_name in merge_df_before_selection["kmer_name"]]
    ax, fig, gs = startfig(6, 6)
    ax.scatter(
        merge_df_before_selection.loc[:, "%Rank_bestAllele"], merge_df_before_selection.loc[:, f"{mrna_expression_col}"], color=color_list, s=5, alpha=0.45, edgecolors="none"
    )
    ax.set_xlabel("%Rank_bestAllele", fontsize=7)
    ax.set_ylabel(f"{mrna_expression_col}", fontsize=7)
    ax.set_title(patient, fontsize=7)
    ax.tick_params("both", labelsize=7)
    ax.set_xlim(-5, 80)
    fig.tight_layout()
    fig.savefig(mixmhcpred_out_txt.parent / "outs" / f"{patient}_percent_rank_mixmhcpred_vs_{mrna_expression_col}_scatter.pdf")
    plt.close()

    ax, fig, gs = startfig(6, 6)
    ax.scatter(
        merge_df_before_selection.loc[:, "%Rank_bestAllele"], merge_df_before_selection.loc[:, f"{mrna_expression_col}"], color=color_list, s=5, alpha=0.45, edgecolors="none"
    )
    ax.set_xlabel("%Rank_bestAllele", fontsize=7)
    ax.set_ylabel(f"{mrna_expression_col}", fontsize=7)
    ax.set_title(patient, fontsize=7)
    ax.tick_params("both", labelsize=7)
    ax.set_xlim(-0.1, 2)
    fig.tight_layout()
    fig.savefig(mixmhcpred_out_txt.parent / "outs" / f"{patient}_percent_rank_mixmhcpred_vs_{mrna_expression_col}_scatter_zoom.pdf")
    plt.close()

    return selected_peptides_df


def parse_mixmhcpred_out_txt_list(
    mixmhcpred_out_txt_list: list,
    patient_df_dict: dict,
    gene_id_col: str,
    mrna_expression_col: str,
    keep_top_number_peptides: int | dict,
    include_cols_patient_df_list: list = None,
    rank_by_binning: bool = False,
    mrna_sorting_weight: float = 0.2,
    percent_rank_threshold: float = 2.0,
    hit_col: str = None,
    confident_binder_percent_rank_threshold: float = None,
):
    """Run `parse_mixmhcpred_out_txt` for every patient's MixMHCpred output and save each patient's selection to its own CSV.

    For each file in `mixmhcpred_out_txt_list`, the patient name is read
    straight from that file's name (as written by `run_mixmhcpred`). The
    user sets how many top-ranked minimal peptides to select for screening
    via `keep_top_number_peptides`, either as a single number applied to
    every patient, or as a dict giving each patient their own number.
    Each patient's selected peptides get a `"patient"` column and are saved
    to `<patient>_selected_peptides_df.csv` next to that patient's input
    file. Diagnostic logging for every patient is collected into one shared
    `outs/parse_mixmhcpred_out_txt.log` file.

    Parameters
    ----------
    mixmhcpred_out_txt_list : list
        MixMHCpred output `.txt` paths (e.g. from
        `run_mixmhcpred_patient_df_dict`), one per patient; patient names
        themselves must not contain underscores, since the patient name is
        read back out of the file name.
    patient_df_dict : dict
        Mapping of patient identifier to that patient's antigen dataframe (see
        `parse_mixmhcpred_out_txt`'s `patient_df`).
    gene_id_col : str
        Column identifying each peptide's source antigen, passed through for
        every patient.
    mrna_expression_col : str
        Column holding each antigen's mRNA expression value, passed through
        for every patient.
    keep_top_number_peptides : int | dict
        Number of top-sort-scoring peptides to keep per patient: either a
        single int applied to all patients, or a dict mapping patient name to
        a per-patient count.
    include_cols_patient_df_list : list, default = None
        Passed through to `parse_mixmhcpred_out_txt` for every patient.
    rank_by_binning : bool, default = False
        Passed through to `parse_mixmhcpred_out_txt` for every patient.
    mrna_sorting_weight : float, default = 0.2
        Passed through to `parse_mixmhcpred_out_txt` for every patient.
    percent_rank_threshold : float, default = 2.0
        Passed through to `parse_mixmhcpred_out_txt` for every patient.
    hit_col : str, default = None
        Passed through to `parse_mixmhcpred_out_txt` for every patient.
    confident_binder_percent_rank_threshold : float, default = None
        Passed through to `parse_mixmhcpred_out_txt` for every patient.

    Returns
    -------
    list[Path]
        Path to the written `<patient>_selected_peptides_df.csv` for each
        entry in `mixmhcpred_out_txt_list`, in the same order (not the
        dataframes themselves).
    """
    patient_df_dict = patient_df_dict.copy()
    selected_peptide_df_file_list = []
    outs_dir = mixmhcpred_out_txt_list[0].parent / "outs"
    outs_dir.mkdir(exist_ok=True)
    log_file_path = outs_dir / "parse_mixmhcpred_out_txt.log"
    with open(log_file_path, "w") as log_file:
        for mixmhcpred_out_txt in mixmhcpred_out_txt_list:
            print(f"Parsing: {mixmhcpred_out_txt}", file=log_file)
            patient = mixmhcpred_out_txt.stem.split("_")[2]

            if isinstance(keep_top_number_peptides, dict):
                keep_top_number_peptides_patient = keep_top_number_peptides[patient]
            elif isinstance(keep_top_number_peptides, int):
                keep_top_number_peptides_patient = keep_top_number_peptides
            else:
                raise TypeError(f"keep_top_number_peptides must be int or dict, got {type(keep_top_number_peptides).__name__}")

            patient_df = patient_df_dict[patient]
            selected_peptides_df = parse_mixmhcpred_out_txt(
                mixmhcpred_out_txt=mixmhcpred_out_txt,
                patient_df=patient_df,
                patient=patient,
                gene_id_col=gene_id_col,
                mrna_expression_col=mrna_expression_col,
                keep_top_number_peptides=keep_top_number_peptides_patient,
                log_file=log_file,
                include_cols_patient_df_list=include_cols_patient_df_list,
                rank_by_binning=rank_by_binning,
                mrna_sorting_weight=mrna_sorting_weight,
                percent_rank_threshold=percent_rank_threshold,
                hit_col=hit_col,
                confident_binder_percent_rank_threshold=confident_binder_percent_rank_threshold,
            )
            selected_peptides_df.loc[:, "patient"] = patient
            selected_peptides_df.to_csv(mixmhcpred_out_txt.parent / f"{patient}_selected_peptides_df.csv")
            selected_peptide_df_file_list.append(mixmhcpred_out_txt.parent / f"{patient}_selected_peptides_df.csv")
            print("\n\n", file=log_file)

    return selected_peptide_df_file_list
