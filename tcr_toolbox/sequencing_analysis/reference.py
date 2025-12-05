from __future__ import annotations

import collections
import json
import os
from re import search
from typing import Union

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from dotenv import load_dotenv
import pysam

from tcr_toolbox.tcr_assembly.constants import tcr_golden_gate_constant_nt_seq_dict
from tcr_toolbox.utils.constants import p20_cloning_sites, p20_pMX_trim_seqs

load_dotenv()
tcr_toolbox_data_path = os.getenv("tcr_toolbox_data_path")


def get_indices(group):
    return group.index.tolist()


def generate_assembly_nt_refs(
    tcr_refs_df: pd.DataFrame,
    epitope_barcode_refs: Union[pd.DatFrame, str, os.PathLike[str], None] = None,
    tcr_name_col_name: str = "name",
    epitope_name_col_name: str = "name",
    fasta_alpha_out_fname: Union[str, os.PathLike[str], None] = None,
    fasta_beta_out_fname: Union[str, os.PathLike[str], None] = None,
    fasta_beta_epitope_plate_seq_out_fname: Union[str, os.PathLike[str], None] = None,
    fasta_epitope_out_fname: Union[str, os.PathLike[str], None] = None,
    trimmed_beta_model_tcr_dict: dict = None,
    model_epitope_dict: dict = None,
    alpha_order_col_name: str = "cdr3j_alpha_nt_order_primers",
    beta_order_col_name: str = "cdr3j_beta_nt_order_primers",
    trim_assembly_primers_from_cdr3j: bool = True,
    v_alpha_col: str = "TRAV_IMGT_allele_collapsed",
    v_beta_col: str = "TRBV_IMGT_allele_collapsed",
    add_v_gene_to_duplicate_filter: bool = True,
    epitope_order_col_name: str = "epitope_order",
    read_length: int = 150,
    epitope_barcode_length: int = 18,
    gtf: bool = True,
    verbose: int = 0,
    trim_constant_seq: bool = True,
):
    """
    Write reference FASTA (.fa) files for TCR alpha, TCR beta, epitope, or combined PAIR-scan TCR alpha & beta & epitope plate alignment.

    The function optionally trims constant sequences, assembly primers, and adjusts for sequencing read length
    and epitope barcode length. Duplicate TCR sequences can be filtered. If duplicate TCR sequences occur, a dictionary mapping kept
    TCR  names to filtered TCR names duplicates is written as `first_idx_map_duplicate_idx_dict.json`.

    Parameters
    ----------
    tcr_refs_df : pd.DataFrame
        DataFrame containing all the TCRs that need to be written to the reference .fa file.
    epitope_barcode_refs : Union[None, pd.DataFrame, str, os.PathLike[str]], default  = None
        DataFrame or fasta file containing antigen/epitope-encoding-minigene + barcode nucleotide order sequences that need to be written to the \
        reference .fa file
    tcr_name_col_name : str, default='name'
        Column in `tcr_refs_df` containing TCR names.
    epitope_name_col_name : str, default='name'
        Column in `epitope_barcode_refs` containing epitope names.
    fasta_alpha_out_fname : Union[str, os.PathLike[str], None], default=None
        Output file path for TCR alpha reference FASTA.
    fasta_beta_out_fname : Union[str, os.PathLike[str], None], default=None
        Output file path for TCR beta reference FASTA.
    fasta_beta_epitope_plate_seq_out_fname : Union[str, os.PathLike[str], None], default=None
        Output file path for plate-based TCR beta and epitope reference FASTA.
    fasta_epitope_out_fname : Union[str, os.PathLike[str], None], default=None
        Output file path for epitope reference FASTA.
    trimmed_beta_model_tcr_dict : dict, default=None
        Dictionary of beta chain sequences from commonly used model TCRs to include in reference files.
    model_epitope_dict : dict, default=None
        Dictionary of commonly used model epitopes to include in reference files.
    alpha_order_col_name : str, default='cdr3j_alpha_nt_order_primers'
        Column in `tcr_refs_df` containing CDR3-J alpha nucleotide sequences.
    beta_order_col_name : str, default='cdr3j_beta_nt_order_primers'
        Column in `tcr_refs_df` containing CDR3-J beta nucleotide sequences.
    trim_assembly_primers_from_cdr3j : bool, default=True
        Whether to trim assembly ortho primers from CDR3-J sequences.
    v_alpha_col : str, default='TRAV_IMGT_allele_collapsed'
        Column in `tcr_refs_df` containing TRAV gene annotations.
    v_beta_col : str, default='TRBV_IMGT_allele_collapsed'
        Column in `tcr_refs_df` containing TRBV gene annotations.
    add_v_gene_to_duplicate_filter : bool, default=True
        Whether to include V genes when filtering duplicate TCR sequences.
    epitope_order_col_name : str, default='epitope_order'
        Column containing epitope CDS and barcode sequences.
    read_length : int, default=150
        Sequencing read length for trimming references.
    epitope_barcode_length : int, default=18
        Epitope barcode length for trimming and alignment purposes.
    gtf : bool, default=True
        Whether to generate corresponding GTF annotation files.
    verbose : int, default=0
        If 1, prints amino acid translations of first 5 reference sequences for verification.
    trim_constant_seq : bool, default=True
        Whether to trim constant sequences from TCR reference sequences.

    Notes
    -----
    Duplicate TCR sequences are handled by keeping one representative reference sequence and storing other
    duplicates names in a JSON dictionary:
    ```
    first_idx_map_duplicate_idx_dict = {
        'Kept_TCR_Name': ['Duplicate_TCR_1_name', 'Duplicate_TCR_2_name', ...],
        ...
    }
    ```
    This dictionary allows consistent alignment and reference usage while preserving information about filtered duplicates.

    The function can also optionally add model TCR beta chains and model epitopes for validation or benchmarking purposes.


    Examples
    --------
    Generate standardized TCR assembly QC TCR alpha and beta reference .fa files:

    >>> generate_assembly_nt_refs(
    ...     tcr_refs_df=tcr_refs_df,
    ...     tcr_name_col_name="name",
    ...     fasta_alpha_out_fname='references/run_name_150bp_alpha.fa',
    ...     fasta_beta_out_fname='references/run_name_150bp_beta.fa',
    ...     alpha_order_col_name="cdr3j_alpha_nt_order_primers",
    ...     beta_order_col_name="cdr3j_beta_nt_order_primers",
    ...     trim_assembly_primers_from_cdr3j=True,
    ...     v_alpha_col="TRAV",
    ...     v_beta_col="TRBV",
    ...     add_v_gene_to_duplicate_filter=True,
    ...     read_length=150,
    ...     gtf=True,
    ...     verbose=1,
    ...     trim_constant_seq=True,
    ... )

    Output:
    tcr_ref_df shape: 2688
    Removing TRAV + CDR3 alpha duplicates...
    Shape before TRAV + CDR3 alpha duplicate removal: 2688
    Shape after duplicate removal: 2135
    Removing TRBV + CDR3 beta duplicates...
    Shape before TRBV + CDR3 beta duplicate removal: 2688
    Shape after duplicate removal: 1537

    Translation of first 5 reconstructed alpha chains:
    ['METLLGLLILWLQLQWVSSKQEVTQIPAALSVPEGENLVLNCSFTDSAIYNLQWFRQDPGKGLTSLLLIQSSQREQTSGRLNASLDKSSGRSTLYIAASQPGDSATYLCAVRPSDKLIFGTGTRLQVFPNI',
    'MKTFAGFSFLFLWLQLDCMSRGEDVEQSLFLSVREGDSSVINCTYTDSSSTYLYWYKQEPGAGLQLLTYIFSNMDMKQDQRLTVLLNKKDKHLSLRIADTQTGDSAIYFCAAAEGGFKTIFGAGTRLFVKANI', ...]

    Translation of first 5 reconstructed beta chains:
    ['MGTSLLCWMALCLLGADHADTGVSQNPRHKITKRGQNVTFRCDPISEHNRLYWYRQTLGQGPEFLTYFQNEAQLEKSRLLSDRFSAERPKGSFSTLEIQRTEQGDSAMYLCASELANEQFFGPGTRLTVLE',
    'MSISLLCCAAFPLLWAGPVNAGVTQTPKFRILKIGQSMTLQCTQDMNHNYMYWYRQDPGMGLKLIYYSVGAGITDKGEVPNGYNVSRSTTEDFPLRLELAAPSQTSVYFCASSTPGQGSYEQYFGPGTRLTVTE', ...]

    Finished writing references!


    Generate PAIR-scan plate-based single-TCR and epitope reference .fa files:

    >>> generate_assembly_nt_refs(
    ...     tcr_refs_df = tcr_refs_df,
    ...     epitope_barcode_refs = epitope_barcode_refs,
    ...     tcr_name_col_name = 'name',
    ...     epitope_name_col_name = 'name',
    ...     fasta_alpha_out_fname = 'references/run_name_150bp_alpha.fa',
    ...     fasta_beta_out_fname = 'references/run_name_150bp_beta.fa',
    ...     fasta_beta_epitope_plate_seq_out_fname = 'references/run_name_150bp_pairscan_beta_epi_plate.fa',
    ...     fasta_epitope_out_fname = 'references/run_name_150bp_epi.fa',
    ...     trimmed_beta_model_tcr_dict = None,
    ...     model_epitope_dict = None,
    ...     alpha_order_col_name = 'cdr3j_alpha_nt_order_primers',
    ...     beta_order_col_name = 'cdr3j_beta_nt_order_primers',
    ...     trim_assembly_primers_from_cdr3j = True,
    ...     v_alpha_col = 'TRAV_IMGT_allele_collapsed',
    ...     v_beta_col = 'TRBV_IMGT_allele_collapsed',
    ...     add_v_gene_to_duplicate_filter = True,
    ...     epitope_order_col_name = 'sequence',
    ...     epitope_barcode_length = 18,
    ...     read_length = 150,
    ...     gtf = True,
    ...     verbose = 1
    ... )
    tcr_ref_df shape: 983
    Removing TRAV + CDR3 alpha duplicates...
    Shape before TRAV + CDR3 alpha duplicate removal: 983
    Shape after duplicate removal: 909
    Removing TRBV + CDR3 beta duplicates...
    Shape before TRBV + CDR3 beta duplicate removal: 983
    Shape after duplicate removal: 850
    epitope_barcode_refs shape: 2778

    Translation of first 5 epitope minigene sequences:
    ['EPPEVGSDCTTIHYNDMCNSSCMGGMNRRPI*', 'EPPEVGSDCTTIHYNDMCNSSCMGGMNRRPI*', ...]

    Translation of first 5 reconstructed alpha chains:
    ['MACPGFLWALVISTCLEFSMAQTVTQSQPEMSVQEAETVTLSCTYDTSESDYYLFWYKQPPSRQMILVIRQEAYKQQNATENRFSVNFQKAAKSFSLKISDSQLGDAAMYFCAYRRGRSGGSEKLVFGKGTKLTVNPYI',
    'MKTFAGFSFLFLWLQLDCMSRGEDVEQSLFLSVREGDSSVINCTYTDSSSTYLYWYKQEPGAGLQLLTYIFSNMDMKQDQRLTVLLNKKDKHLSLRIADTQTGDSAIYFCAEKGSGGGADGLTFGKGTHLIIQPYI', ...]

    Translation of first 5 reconstructed beta chains:
    ['MGTSLLCWMALCLLGADHADTGVSQNPRHKITKRGQNVTFRCDPISEHNRLYWYRQTLGQGPEFLTYFQNEAQLEKSRLLSDRFSAERPKGSFSTLEIQRTEQGDSAMYLCASSPQGSTGELFFGEGSRLTVLE',
    'MGCRLLCCAVLCLLGAVPIDTEVTQTPKHLVMGMTNKKSLKCEQHMGHRAMYWYKQKAKKPPELMFVYSYEKLSINESVPSRFSPECPNSSLLNLHLHALQPEDSALYLCASSQTPGQTGSPLHFGNGTRLTVTE', ...]

    Finished writing references!
    """
    tcr_refs_df = tcr_refs_df.copy()
    if not tcr_refs_df.index.is_unique:
        print("WARNING: tcr_refs_df index was not resetted! Resetting index...")
        tcr_refs_df.reset_index(drop=True, inplace=True)

    if tcr_refs_df.loc[:, alpha_order_col_name].isna().any() or tcr_refs_df.loc[:, beta_order_col_name].isna().any():
        print("TCRs without CDR3-J order sequences found. Please verify that these are wells that on purpose did not get TCRs assigned.\nRemoving these TCRs!")
        print("Shape before CDR3-J order sequence NaN removal:", tcr_refs_df.shape[0])
        if tcr_refs_df.loc[:, alpha_order_col_name].isna().any():
            tcr_refs_df = tcr_refs_df.loc[~tcr_refs_df.loc[:, alpha_order_col_name].isna(), :]
        if tcr_refs_df.loc[:, beta_order_col_name].isna().any():
            tcr_refs_df = tcr_refs_df.loc[~tcr_refs_df.loc[:, beta_order_col_name].isna(), :]
        print("Shape after CDR3-J order sequence NaN removal:", tcr_refs_df.shape[0])
    tcr_refs_df = tcr_refs_df.copy()
    if tcr_refs_df.loc[:, alpha_order_col_name].isna().any() or tcr_refs_df.loc[:, beta_order_col_name].isna().any():
        print("TCRs without CDR3-J order sequences found. Please verify that these are wells that on purpose did not get TCRs assigned.\nRemoving these TCRs!")
        print("Shape before CDR3-J order sequence NaN removal:", tcr_refs_df.shape[0])
        if tcr_refs_df.loc[:, alpha_order_col_name].isna().any():
            tcr_refs_df = tcr_refs_df.loc[~tcr_refs_df.loc[:, alpha_order_col_name].isna(), :]
        if tcr_refs_df.loc[:, beta_order_col_name].isna().any():
            tcr_refs_df = tcr_refs_df.loc[~tcr_refs_df.loc[:, beta_order_col_name].isna(), :]
        print("Shape after CDR3-J order sequence NaN removal:", tcr_refs_df.shape[0])

    alpha_tcr_refs_df = tcr_refs_df.copy()
    beta_tcr_refs_df = tcr_refs_df.copy()

    if model_epitope_dict is None:
        model_epitope_dict = {}
    with open(
        os.path.join(
            tcr_toolbox_data_path,
            "tcr_toolbox_datasets",
            "tcr_assembly",
            "references",
            "v_gene_stock_sequences",
            "cys104_trimmed_v_gene_stock_nt_dict.json",
        )
    ) as infile:
        cys104_trimmed_v_gene_stock_nt_dict = json.load(infile)

    print("tcr_ref_df shape:", tcr_refs_df.shape[0])
    if fasta_alpha_out_fname:
        if trim_assembly_primers_from_cdr3j:
            # plate Fw/Rev + Fw/Rev orthoprimers are together 40 bp:
            alpha_tcr_refs_df[alpha_order_col_name] = alpha_tcr_refs_df[alpha_order_col_name].str[40:-40]

        if add_v_gene_to_duplicate_filter:
            alpha_tcr_refs_df["TRAV_CDR3_alpha_filter"] = alpha_tcr_refs_df[v_alpha_col] + "_" + alpha_tcr_refs_df[alpha_order_col_name]

        if add_v_gene_to_duplicate_filter and alpha_tcr_refs_df["TRAV_CDR3_alpha_filter"].duplicated().any():
            duplicate_indices = alpha_tcr_refs_df.groupby("TRAV_CDR3_alpha_filter").apply(get_indices)[
                alpha_tcr_refs_df.groupby("TRAV_CDR3_alpha_filter").apply(get_indices).apply(len) > 1
            ]
            alpha_first_idx_map_duplicate_idx_dict = {}
            for idx_list in duplicate_indices.values:
                alpha_first_idx_map_duplicate_idx_dict[alpha_tcr_refs_df.loc[idx_list[0], tcr_name_col_name]] = [
                    alpha_tcr_refs_df.loc[idx, tcr_name_col_name] for idx in idx_list[1:]
                ]
            file_base = os.path.join(os.path.dirname(fasta_alpha_out_fname), os.path.splitext(os.path.basename(fasta_alpha_out_fname))[0])
            with open(file_base + "_first_idx_map_duplicate_idx_dict.json", "w") as json_file:
                json.dump(alpha_first_idx_map_duplicate_idx_dict, json_file, indent=4)
            print("Removing TRAV + CDR3 alpha duplicates...")
            print("Shape before TRAV + CDR3 alpha duplicate removal:", alpha_tcr_refs_df.shape[0])
            alpha_tcr_refs_df.drop_duplicates("TRAV_CDR3_alpha_filter", inplace=True)
            alpha_tcr_refs_df.reset_index(drop=True, inplace=True)
            print("Shape after duplicate removal:", alpha_tcr_refs_df.shape[0])
        elif alpha_tcr_refs_df[alpha_order_col_name].duplicated().any():
            duplicate_indices = alpha_tcr_refs_df.groupby(alpha_order_col_name).apply(get_indices)[
                alpha_tcr_refs_df.groupby(alpha_order_col_name).apply(get_indices).apply(len) > 1
            ]
            alpha_first_idx_map_duplicate_idx_dict = {}
            for idx_list in duplicate_indices.values:
                alpha_first_idx_map_duplicate_idx_dict[alpha_tcr_refs_df.loc[idx_list[0], tcr_name_col_name]] = [
                    alpha_tcr_refs_df.loc[idx, tcr_name_col_name] for idx in idx_list[1:]
                ]
            file_base = os.path.join(os.path.dirname(fasta_alpha_out_fname), os.path.splitext(os.path.basename(fasta_alpha_out_fname))[0])

            with open(file_base + "_first_idx_map_duplicate_idx_dict.json", "w") as json_file:
                json.dump(alpha_first_idx_map_duplicate_idx_dict, json_file, indent=4)
            print("Removing CDR3 alpha duplicates...")
            print("Shape before CDR3 alpha duplicate removal:", alpha_tcr_refs_df.shape[0])
            alpha_tcr_refs_df.drop_duplicates(alpha_order_col_name, inplace=True)
            alpha_tcr_refs_df.reset_index(drop=True, inplace=True)
            print("Shape after duplicate removal:", alpha_tcr_refs_df.shape[0])

    if fasta_beta_out_fname or fasta_beta_epitope_plate_seq_out_fname:
        if trim_assembly_primers_from_cdr3j:
            # plate Fw/Rev + Fw/Rev orthoprimers are together 40 bp:
            beta_tcr_refs_df[beta_order_col_name] = beta_tcr_refs_df[beta_order_col_name].str[40:-40]

        if add_v_gene_to_duplicate_filter:
            beta_tcr_refs_df["TRBV_CDR3_beta_filter"] = beta_tcr_refs_df[v_beta_col] + "_" + beta_tcr_refs_df[beta_order_col_name]

        if add_v_gene_to_duplicate_filter and beta_tcr_refs_df["TRBV_CDR3_beta_filter"].duplicated().any():
            duplicate_indices = beta_tcr_refs_df.groupby("TRBV_CDR3_beta_filter").apply(get_indices)[
                beta_tcr_refs_df.groupby("TRBV_CDR3_beta_filter").apply(get_indices).apply(len) > 1
            ]
            beta_first_idx_map_duplicate_idx_dict = {}
            for idx_list in duplicate_indices.values:
                beta_first_idx_map_duplicate_idx_dict[beta_tcr_refs_df.loc[idx_list[0], tcr_name_col_name]] = [beta_tcr_refs_df.loc[idx, tcr_name_col_name] for idx in idx_list[1:]]
            if fasta_beta_out_fname:
                file_base = os.path.join(os.path.dirname(fasta_beta_out_fname), os.path.splitext(os.path.basename(fasta_beta_out_fname))[0])
            elif fasta_beta_epitope_plate_seq_out_fname:
                file_base = os.path.join(os.path.dirname(fasta_beta_epitope_plate_seq_out_fname), os.path.splitext(os.path.basename(fasta_beta_epitope_plate_seq_out_fname))[0])
            with open(file_base + "_first_idx_map_duplicate_idx_dict.json", "w") as json_file:
                json.dump(beta_first_idx_map_duplicate_idx_dict, json_file, indent=4)
            print("Removing TRBV + CDR3 beta duplicates...")
            print("Shape before TRBV + CDR3 beta duplicate removal:", beta_tcr_refs_df.shape[0])
            beta_tcr_refs_df.drop_duplicates("TRBV_CDR3_beta_filter", inplace=True)
            beta_tcr_refs_df.reset_index(drop=True, inplace=True)
            print("Shape after duplicate removal:", beta_tcr_refs_df.shape[0])
        elif beta_tcr_refs_df[beta_order_col_name].duplicated().any():
            duplicate_indices = beta_tcr_refs_df.groupby(beta_order_col_name).apply(get_indices)[beta_tcr_refs_df.groupby(beta_order_col_name).apply(get_indices).apply(len) > 1]
            beta_first_idx_map_duplicate_idx_dict = {}
            for idx_list in duplicate_indices.values:
                beta_first_idx_map_duplicate_idx_dict[beta_tcr_refs_df.loc[idx_list[0], tcr_name_col_name]] = [beta_tcr_refs_df.loc[idx, tcr_name_col_name] for idx in idx_list[1:]]
            if fasta_beta_out_fname:
                file_base = os.path.join(os.path.dirname(fasta_beta_out_fname), os.path.splitext(os.path.basename(fasta_beta_out_fname))[0])
            elif fasta_beta_epitope_plate_seq_out_fname:
                file_base = os.path.join(os.path.dirname(fasta_beta_epitope_plate_seq_out_fname), os.path.splitext(os.path.basename(fasta_beta_epitope_plate_seq_out_fname))[0])
            with open(file_base + "_first_idx_map_duplicate_idx_dict.json", "w") as json_file:
                json.dump(beta_first_idx_map_duplicate_idx_dict, json_file, indent=4)
            print("Removing CDR3 beta duplicates...")
            print("Shape before CDR3 beta duplicate removal:", beta_tcr_refs_df.shape[0])
            beta_tcr_refs_df.drop_duplicates(beta_order_col_name, inplace=True)
            beta_tcr_refs_df.reset_index(drop=True, inplace=True)
            print("Shape after duplicate removal:", beta_tcr_refs_df.shape[0])

    if isinstance(epitope_barcode_refs, pd.DataFrame):
        epitope_barcode_refs_df = epitope_barcode_refs.copy()
    elif isinstance(epitope_barcode_refs, (str, os.PathLike)):
        epitope_barcode_refs_seq_df_construct_list = []
        epitope_barcode_refs_name_df_construct_list = []
        with pysam.FastxFile(epitope_barcode_refs, "r") as fasta_in:
            for entry in fasta_in:
                epitope_barcode_refs_name_df_construct_list.append(entry.name)
                epitope_barcode_refs_seq_df_construct_list.append(entry.sequence)

        epitope_barcode_refs_df = pd.DataFrame(
            {epitope_name_col_name: epitope_barcode_refs_name_df_construct_list, epitope_order_col_name: epitope_barcode_refs_seq_df_construct_list}
        )

    if epitope_barcode_refs is not None:
        print("epitope_barcode_refs_df shape:", epitope_barcode_refs_df.shape[0])
        if epitope_barcode_refs_df[epitope_order_col_name].duplicated().any():
            print("Removing epitope duplicates...")
            print("Shape before epitope duplicate removal:", epitope_barcode_refs_df.shape[0])
            epitope_barcode_refs_df.drop_duplicates(epitope_order_col_name, inplace=True)
            epitope_barcode_refs_df.reset_index(drop=True, inplace=True)
            print("Shape after duplicate removal:", epitope_barcode_refs_df.shape[0])

    print("\n\n")
    if epitope_barcode_refs_df is not None:
        epitope_nt_ref_list = []
        epitope_nt_dict = collections.defaultdict(str)

        for epitope_idx in epitope_barcode_refs_df.index:
            epitope_nt_dict[epitope_idx] = epitope_barcode_refs_df.loc[epitope_idx, epitope_order_col_name][
                search(p20_cloning_sites["end_5_BsmBI_p20"], epitope_barcode_refs_df.loc[epitope_idx, epitope_order_col_name]).span()[1] : search(
                    p20_cloning_sites["end_3_BsmBI_p20"], epitope_barcode_refs_df.loc[epitope_idx, epitope_order_col_name]
                ).span()[0]
            ]

            if not len(epitope_nt_dict[epitope_idx][:-epitope_barcode_length]) % 3 == 0:
                raise Exception("Peptide coding nucleotide sequence is not a multiple of 3!")

            if not (epitope_nt_dict[epitope_idx][-(3 + epitope_barcode_length) : -epitope_barcode_length] == "TAA") | (
                epitope_nt_dict[epitope_idx][-(3 + epitope_barcode_length) : -epitope_barcode_length] == "TAG"
            ):
                raise Exception(
                    "Input barcode length is incorrect!", epitope_nt_dict[epitope_idx], epitope_nt_dict[epitope_idx][-(3 + epitope_barcode_length) : -epitope_barcode_length]
                )

            # The minimal length to align to with bwa is 36 nt,
            # as epitopes might be very similar, ideally we only align to the non-similar barcode nt sequences:
            if epitope_barcode_length < 36:
                epitope_nt_ref_list.append(epitope_nt_dict[epitope_idx][-36:])
            else:
                epitope_nt_ref_list.append(epitope_nt_dict[epitope_idx][-epitope_barcode_length:])

        epitope_names = epitope_barcode_refs[epitope_name_col_name].values.tolist()

        if verbose == 1:
            print("Translation of first 5 epitope sequences:")
            print([str(Seq(epitope[:-epitope_barcode_length]).translate()) for epitope in epitope_nt_dict.values()[:5]])
            print("\n\n")

        if model_epitope_dict:
            epitope_nt_ref_list += list(model_epitope_dict.values())
            epitope_names += list(model_epitope_dict.keys())

        if len(epitope_names) != len(epitope_nt_ref_list):
            raise Exception("Length epitope names list is not equal to length epitope sequence list!")

        with open(fasta_epitope_out_fname, mode="w") as fasta_out:
            for name, epi_seq in zip(epitope_names, epitope_nt_ref_list):
                fasta_out.write(">" + str(name) + "\n")
                fasta_out.write(epi_seq + "\n")

        if gtf:
            with open(os.path.splitext(fasta_epitope_out_fname)[0] + ".gtf", mode="w") as gtf_out:
                for name, tcr_seq in zip(epitope_names, epitope_nt_ref_list):
                    gtf_out.write(
                        str(name)
                        + "\thavana\ttranscript\t"
                        + "1\t"
                        + str(len(tcr_seq))
                        + '\t.\t+\t.\tgene_id "'
                        + str(name)
                        + '"; '
                        + 'gene_version "1"; gene_name "'
                        + str(name)
                        + '";\n'
                    )
                    gtf_out.write(
                        str(name)
                        + "\thavana\texon\t"
                        + "1\t"
                        + str(len(tcr_seq))
                        + '\t.\t+\t.\tgene_id "'
                        + str(name)
                        + '"; '
                        + 'gene_version "1"; exon_number "1"; gene_name "'
                        + str(name)
                        + '";\n'
                    )

    if fasta_alpha_out_fname:
        alpha_nt_ref_list = []

        # Check whether correct start and end amino acids are in alpha order after trimming the golden gate sites by 9 and -10:
        if not all([True if Seq(cdr3_seq).translate()[0] == "C" else False for cdr3_seq in alpha_tcr_refs_df[alpha_order_col_name].str[9:-10]]) & all(
            [True if Seq(cdr3_seq).translate()[-1] == "I" else False for cdr3_seq in alpha_tcr_refs_df[alpha_order_col_name].str[9:-10]]
        ):
            raise Exception("Alpha order column is in wrong format!")

        for tcr_idx in alpha_tcr_refs_df.index:
            alpha_tcr_refs_df.loc[tcr_idx, v_alpha_col] = alpha_tcr_refs_df.loc[tcr_idx, v_alpha_col].replace("/", "_")
            if "*" in alpha_tcr_refs_df.loc[tcr_idx, v_alpha_col]:
                # We only have *01 alleles in our V gene stocks. The + '*01' should be changed once we also have other alleles:
                if alpha_tcr_refs_df.loc[tcr_idx, v_alpha_col].split("*")[1] != "01":
                    print("Warning: non-*01 allele detected in:", v_alpha_col)
                alpha_nt_ref_list.append(
                    cys104_trimmed_v_gene_stock_nt_dict[alpha_tcr_refs_df.loc[tcr_idx, v_alpha_col].split("*")[0] + "*01"]
                    + alpha_tcr_refs_df.loc[tcr_idx, alpha_order_col_name][9:-10]
                )
            elif "*" not in alpha_tcr_refs_df.loc[tcr_idx, v_alpha_col]:
                alpha_nt_ref_list.append(
                    cys104_trimmed_v_gene_stock_nt_dict[alpha_tcr_refs_df.loc[tcr_idx, v_alpha_col] + "*01"] + alpha_tcr_refs_df.loc[tcr_idx, alpha_order_col_name][9:-10]
                )
            else:
                raise Exception("TRAV encoding is broken!")

        # Check whether reconstructed nt sequences are a multiple of 3:
        if not all([False if len(tcr_seq) % 3 != 0 else True for tcr_seq in alpha_nt_ref_list]):
            tmp_incorrect_len_alpha_dict = collections.defaultdict(str)
            for i, tcr_seq in enumerate(alpha_nt_ref_list):
                if len(tcr_seq) % 3 != 0:
                    tmp_incorrect_len_alpha_dict[i] = tcr_seq
            raise Exception("Reconstructed alpha chain is not a multiple of 3!", tmp_incorrect_len_alpha_dict)

        alpha_names = alpha_tcr_refs_df[tcr_name_col_name].values.tolist()

        if len(alpha_names) != len(alpha_nt_ref_list):
            raise Exception("Length alpha names list is not equal to length alpha sequence list!")

        if verbose == 1:
            print("Translation of first 5 reconstructed alpha chains:")
            print([str(Seq(tcr).translate()) for tcr in alpha_nt_ref_list[:5]])
            print("\n\n")

        # If you have an alpha primer, implement alpha trimming here similar as done for fasta_beta_fname:
        if trim_constant_seq:
            alpha_nt_ref_list = [tcr_seq[-(read_length - len(p20_pMX_trim_seqs["pcr_pMX_rev_muTCRa"])) :] for tcr_seq in alpha_nt_ref_list].copy()

        with open(fasta_alpha_out_fname, mode="w") as fasta_out:
            for name, tcr_seq in zip(alpha_names, alpha_nt_ref_list):
                fasta_out.write(">" + str(name) + "\n")
                fasta_out.write(tcr_seq + "\n")

        if gtf:
            with open(os.path.splitext(fasta_alpha_out_fname)[0] + ".gtf", mode="w") as gtf_out:
                for name, tcr_seq in zip(alpha_names, alpha_nt_ref_list):
                    gtf_out.write(
                        str(name)
                        + "\thavana\ttranscript\t"
                        + "1\t"
                        + str(len(tcr_seq))
                        + '\t.\t+\t.\tgene_id "'
                        + str(name)
                        + '"; '
                        + 'gene_version "1"; gene_name "'
                        + str(name)
                        + '";\n'
                    )
                    gtf_out.write(
                        str(name)
                        + "\thavana\texon\t"
                        + "1\t"
                        + str(len(tcr_seq))
                        + '\t.\t+\t.\tgene_id "'
                        + str(name)
                        + '"; '
                        + 'gene_version "1"; exon_number "1"; gene_name "'
                        + str(name)
                        + '";\n'
                    )

    if fasta_beta_out_fname or fasta_beta_epitope_plate_seq_out_fname:
        beta_nt_ref_list = []

        # Check whether correct start and end amino acids are in beta order after trimming the golden gate sites by 8 and -9:
        if not all([True if Seq(cdr3_seq).translate()[0] == "C" else False for cdr3_seq in beta_tcr_refs_df[beta_order_col_name].str[8:-9]]) & all(
            [True if Seq(cdr3_seq).translate()[-1] == "E" else False for cdr3_seq in beta_tcr_refs_df[beta_order_col_name].str[8:-9]]
        ):
            raise Exception("Beta order column is in wrong format!")

        for tcr_idx in beta_tcr_refs_df.index:
            beta_tcr_refs_df.loc[tcr_idx, v_beta_col] = beta_tcr_refs_df.loc[tcr_idx, v_beta_col].replace("/", "_")
            if "*" in beta_tcr_refs_df.loc[tcr_idx, v_beta_col]:
                if beta_tcr_refs_df.loc[tcr_idx, v_beta_col].split("*")[1] != "01":
                    print("Warning: non-*01 allele detected in:", v_beta_col)
                # We only have *01 alleles in our V gene stocks. The + '*01' should be changed once we also have other alleles:
                beta_nt_ref_list.append(
                    cys104_trimmed_v_gene_stock_nt_dict[beta_tcr_refs_df.loc[tcr_idx, v_beta_col].split("*")[0] + "*01"] + beta_tcr_refs_df.loc[tcr_idx, beta_order_col_name][8:-9]
                )
            elif "*" not in beta_tcr_refs_df.loc[tcr_idx, v_beta_col]:
                beta_nt_ref_list.append(
                    cys104_trimmed_v_gene_stock_nt_dict[beta_tcr_refs_df.loc[tcr_idx, v_beta_col] + "*01"] + beta_tcr_refs_df.loc[tcr_idx, beta_order_col_name][8:-9]
                )
            else:
                raise Exception("TRBV encoding is broken!")

        # Check whether reconstructed nt sequences are a multiple of 3:
        if not all([False if len(tcr_seq) % 3 != 0 else True for tcr_seq in beta_nt_ref_list]):
            tmp_incorrect_len_beta_dict = collections.defaultdict(str)
            for i, tcr_seq in enumerate(beta_nt_ref_list):
                if len(tcr_seq) % 3 != 0:
                    tmp_incorrect_len_beta_dict[i] = tcr_seq
            raise Exception("Reconstructed beta chain is not a multiple of 3!", tmp_incorrect_len_beta_dict)

        beta_names = beta_tcr_refs_df[tcr_name_col_name].values.tolist()

        if verbose == 1:
            print("Translation of first 5 reconstructed beta chains:")
            print([str(Seq(tcr).translate()) for tcr in beta_nt_ref_list[:5]])
            print("\n\n")

        if len(beta_names) != len(beta_nt_ref_list):
            raise Exception("Length beta names list is not equal to length beta sequence list!")

        beta_nt_ref_list_untrimmed = beta_nt_ref_list.copy()
        if trim_constant_seq:
            beta_nt_ref_list = [tcr_seq[-(read_length - len(p20_pMX_trim_seqs["pcr_pMX_rev_muTCRb"])) :] for tcr_seq in beta_nt_ref_list].copy()

        if fasta_beta_out_fname:
            with open(fasta_beta_out_fname, mode="w") as fasta_out:
                for name, tcr_seq in zip(beta_names, beta_nt_ref_list):
                    fasta_out.write(">" + str(name) + "\n")
                    fasta_out.write(tcr_seq + "\n")

            if gtf:
                with open(os.path.splitext(fasta_beta_out_fname)[0] + ".gtf", mode="w") as gtf_out:
                    for name, tcr_seq in zip(beta_names, beta_nt_ref_list):
                        gtf_out.write(
                            str(name)
                            + "\thavana\ttranscript\t"
                            + "1\t"
                            + str(len(tcr_seq))
                            + '\t.\t+\t.\tgene_id "'
                            + str(name)
                            + '"; '
                            + 'gene_version "1"; gene_name "'
                            + str(name)
                            + '";\n'
                        )
                        gtf_out.write(
                            str(name)
                            + "\thavana\texon\t"
                            + "1\t"
                            + str(len(tcr_seq))
                            + '\t.\t+\t.\tgene_id "'
                            + str(name)
                            + '"; '
                            + 'gene_version "1"; exon_number "1"; gene_name "'
                            + str(name)
                            + '";\n'
                        )

        if fasta_beta_epitope_plate_seq_out_fname:
            # Plate-based is single end 150 bp sequencing.
            # Barcode + UMI + RT priming site + small piece of muTCRb constant chain is 45 bp
            # but Refs contain 3 bp from muTCRb.
            # So, the TCR ref str need to be trimmed to the remaining length that can be covered by the remaining length
            # by doing tcr_seq[-(read_length - 42):]
            beta_nt_ref_list = [tcr_seq[-(read_length - 42) :] for tcr_seq in beta_nt_ref_list_untrimmed].copy()
            beta_names = ["tcr_" + tcr_name for tcr_name in beta_names].copy()

            # I now assumed that only for plate-seq experiments we want to be able to quickly add control model TCRs
            # That is why the length of all model TCRs need to be equal to (read_length-45), assuming read_length = 150 bp
            if trimmed_beta_model_tcr_dict:
                if not all([len(tcr_seq) == (read_length - 45) for tcr_seq in trimmed_beta_model_tcr_dict.values()]):
                    raise Exception("Model TCRs are not all trimmed to a length of (read_length-45)")

                else:
                    beta_nt_ref_list += list(trimmed_beta_model_tcr_dict.values())
                    beta_names += ["tcr_" + tcr_name for tcr_name in list(trimmed_beta_model_tcr_dict.keys())]

            # You only want to add model epitopes for aligning plate-based data that needs trimmed beta TCRs:
            if epitope_barcode_refs is not None:
                beta_nt_ref_list += epitope_nt_ref_list
                beta_names += ["epi_" + epitope_name for epitope_name in epitope_names]

            if model_epitope_dict and epitope_barcode_refs is not None:
                beta_nt_ref_list += list(model_epitope_dict.values())
                beta_names += ["epi_" + epitope_name for epitope_name in list(model_epitope_dict.keys())]

            if len(beta_names) != len(beta_nt_ref_list):
                raise Exception("Length beta names list is not equal to length beta sequence list!")

            with open(fasta_beta_epitope_plate_seq_out_fname, mode="w") as fasta_out:
                for name, tcr_seq in zip(beta_names, beta_nt_ref_list):
                    fasta_out.write(">" + str(name) + "\n")
                    fasta_out.write(tcr_seq + "\n")

            # you want to always write the gtf file for plate-based sequence alignment:
            with open(os.path.splitext(fasta_beta_epitope_plate_seq_out_fname)[0] + ".gtf", mode="w") as gtf_out:
                for name, tcr_seq in zip(beta_names, beta_nt_ref_list):
                    gtf_out.write(
                        str(name)
                        + "\thavana\ttranscript\t"
                        + "1\t"
                        + str(len(tcr_seq))
                        + '\t.\t+\t.\tgene_id "'
                        + str(name)
                        + '"; '
                        + 'gene_version "1"; gene_name "'
                        + str(name)
                        + '";\n'
                    )
                    gtf_out.write(
                        str(name)
                        + "\thavana\texon\t"
                        + "1\t"
                        + str(len(tcr_seq))
                        + '\t.\t+\t.\tgene_id "'
                        + str(name)
                        + '"; '
                        + 'gene_version "1"; exon_number "1"; gene_name "'
                        + str(name)
                        + '";\n'
                    )

    return print("Finished writing references!")


def generate_assembly_nanopore_nt_refs(
    tcr_refs_df: pd.DataFrame,
    tcr_name_col_name: str = "name",
    fasta_out_fname: Union[str, os.PathLike[str], None] = None,
    alpha_order_col_name: str = "cdr3j_alpha_nt_order_primers",
    beta_order_col_name: str = "cdr3j_beta_nt_order_primers",
    trim_assembly_primers_from_cdr3j: bool = True,
    include_mu_constant_beta_and_p2a: bool = False,
    v_alpha_col: str = "TRAV_IMGT_allele_collapsed",
    v_beta_col: str = "TRBV_IMGT_allele_collapsed",
    add_number_of_negatives_by_mutating_full_refs: int = None,
    add_number_of_negatives_by_mutating_cdr3j_in_refs: int = None,
    add_number_of_negatives_by_mutating_v_gene_in_refs: int = None,
    min_nt_diff_negative_ref_seqs: int = None,
    max_nt_diff_negative_ref_seqs: int = None,
    verbose: int = 0,
):
    """
    Generate full-length TCR reference sequences for nanopore sequencing QC and write them to a FASTA file.

    This function reconstructs full-length TCR alpha and beta nucleotide sequences from
    provided CDR3-J order sequences and V gene annotations, optionally trims assembly primers,
    and can append mu constant beta + P2A sequences. Duplicate full-length TCRs are removed,
    and optionally, negative TCR reference sequences can be generated by introducing mutations
    in full TCRs, CDR3-J regions, or V genes.

    Parameters
    ----------
    tcr_refs_df : pd.DataFrame
        DataFrame containing TCR information for reference generation.
    tcr_name_col_name : str, default 'name'
        Column in `tcr_refs_df` containing TCR names.
    fasta_out_fname : Union[str, os.PathLike[str], None], default None
        Path to output FASTA file for full-length TCR nanopore reference sequences.
    alpha_order_col_name : str, default 'cdr3j_alpha_nt_order_primers'
        Column in `tcr_refs_df` containing CDR3-J alpha nucleotide sequences.
    beta_order_col_name : str, default 'cdr3j_beta_nt_order_primers'
        Column in `tcr_refs_df` containing CDR3-J beta nucleotide sequences.
    trim_assembly_primers_from_cdr3j : bool, default True
        Whether to trim 40 bp assembly orthoprimers from CDR3-J sequences.
    include_mu_constant_beta_and_p2a : bool, default False
        If True, appends mu constant beta and P2A sequence between beta and alpha chains.
    v_alpha_col : str, default 'TRAV_IMGT_allele_collapsed'
        Column containing TRAV gene annotations.
    v_beta_col : str, default 'TRBV_IMGT_allele_collapsed'
        Column containing TRBV gene annotations.
    add_number_of_negatives_by_mutating_full_refs : int, default None
        Number of negative sequences to generate by mutating full TCR sequences.
    add_number_of_negatives_by_mutating_cdr3j_in_refs : int, default None
        Number of negative sequences to generate by mutating CDR3-J regions in TCR sequences.
    add_number_of_negatives_by_mutating_v_gene_in_refs : int, default None
        Number of negative sequences to generate by mutating V gene regions in TCR sequences.
    min_nt_diff_negative_ref_seqs : int, default None
        Minimum number of nucleotide differences for generated negative sequences.
    max_nt_diff_negative_ref_seqs : int, default None
        Maximum number of nucleotide differences for generated negative sequences.
    verbose : int, default 0
        If 1, prints nucleotide reference sequence amino acid translations.

    Raises
    ------
    Exception
        - If alpha or beta order columns are in incorrect format (wrong start/end amino acids).
        - If reconstructed sequences are not multiples of 3.
        - If alpha and beta reference lists or TCR names become out of sync.
        - If more negative sequences are requested than available references.

    Notes
    -----
    - TCRs without CDR3-J order sequences are removed.
    - Duplicate full-length TCRs are identified based on paired alpha + beta sequences
      and removed. A JSON file mapping retained TCR names to duplicates is saved:
    ```
    first_idx_map_duplicate_idx_dict = {
        'Kept_TCR_Name': ['Duplicate_TCR_1_name', 'Duplicate_TCR_2_name', ...],
        ...
    }
    ```
    - Negative reference sequences (if requested) introduce mutations in selected regions
      and are appended to the reference FASTA.

    Outputs
    -------
    - Writes a FASTA file (`fasta_out_fname`) with full-length TCR sequences.
    - Writes a BED file (`.bed`) with the start and end positions of each sequence.
    - Optionally writes a `_first_idx_map_duplicate_idx_dict.json` mapping file for duplicates.

    Example
    -------
    >>> generate_assembly_nanopore_nt_refs(
    ...     tcr_refs_df=tcr_refs_df,
    ...     tcr_name_col_name="name",
    ...     fasta_out_fname="/references/run_name_nanopore.fa",
    ...     alpha_order_col_name="cdr3j_alpha_nt_order_primers",
    ...     beta_order_col_name="cdr3j_beta_nt_order_primers",
    ...     trim_assembly_primers_from_cdr3j=True,
    ...     include_mu_constant_beta_and_p2a=True,
    ...     v_alpha_col="TRAV_IMGT_allele_collapsed",
    ...     v_beta_col="TRBV_IMGT_allele_collapsed",
    ...     add_number_of_negatives_by_mutating_full_refs=None,
    ...     add_number_of_negatives_by_mutating_cdr3j_in_refs=None,
    ...     add_number_of_negatives_by_mutating_v_gene_in_refs=None,
    ...     min_nt_diff_negative_ref_seqs=None,
    ...     max_nt_diff_negative_ref_seqs=None,
    ...     verbose=0,
    ... )

    Output:
    TCRs without CDR3-J order sequences found. Please verify that these are wells that on purpose did not get TCRs assigned.
    Removing these TCRs!
    Shape before CDR3-J order sequence NaN removal: 7680
    Shape after CDR3-J order sequence NaN removal: 7600
    Removing paired alpha and beta full-length TCR nt duplicates...
    Shape before TCR duplicate removal: 7600
    Shape after TCR duplicate removal: 3693
    """

    rng = np.random.default_rng()

    tcr_refs_df = tcr_refs_df.copy()
    if tcr_refs_df.loc[:, alpha_order_col_name].isna().any() or tcr_refs_df.loc[:, beta_order_col_name].isna().any():
        print("TCRs without CDR3-J order sequences found. Please verify that these are wells that on purpose did not get TCRs assigned.\nRemoving these TCRs!")
        print("Shape before CDR3-J order sequence NaN removal:", tcr_refs_df.shape[0])
        if tcr_refs_df.loc[:, alpha_order_col_name].isna().any():
            tcr_refs_df = tcr_refs_df.loc[~tcr_refs_df.loc[:, alpha_order_col_name].isna(), :]
        if tcr_refs_df.loc[:, beta_order_col_name].isna().any():
            tcr_refs_df = tcr_refs_df.loc[~tcr_refs_df.loc[:, beta_order_col_name].isna(), :]
        print("Shape after CDR3-J order sequence NaN removal:", tcr_refs_df.shape[0])

    with open(
        os.path.join(
            tcr_toolbox_data_path,
            "tcr_toolbox_datasets",
            "tcr_assembly",
            "references",
            "v_gene_stock_sequences",
            "cys104_trimmed_v_gene_stock_nt_dict.json",
        )
    ) as infile:
        cys104_trimmed_v_gene_stock_nt_dict = json.load(infile)

    if trim_assembly_primers_from_cdr3j:
        # plate Fw/Rev + Fw/Rev orthoprimers are together 40 bp:
        tcr_refs_df.loc[:, alpha_order_col_name] = tcr_refs_df.loc[:, alpha_order_col_name].str[40:-40]
        tcr_refs_df.loc[:, beta_order_col_name] = tcr_refs_df.loc[:, beta_order_col_name].str[40:-40]

    tcr_refs_df.loc[:, "TCR_filter"] = (
        tcr_refs_df.loc[:, v_alpha_col] + "_" + tcr_refs_df.loc[:, alpha_order_col_name] + "__" + tcr_refs_df.loc[:, v_beta_col] + "_" + tcr_refs_df.loc[:, beta_order_col_name]
    )

    if tcr_refs_df.loc[:, "TCR_filter"].duplicated().any():
        duplicate_indices = tcr_refs_df.groupby("TCR_filter").apply(get_indices)[tcr_refs_df.groupby("TCR_filter").apply(get_indices).apply(len) > 1]
        first_idx_map_duplicate_idx_dict = {}
        for idx_list in duplicate_indices.values:
            first_idx_map_duplicate_idx_dict[tcr_refs_df.loc[idx_list[0], tcr_name_col_name]] = [tcr_refs_df.loc[idx, tcr_name_col_name] for idx in idx_list[1:]]
        file_base = os.path.join(os.path.dirname(fasta_out_fname), os.path.splitext(os.path.basename(fasta_out_fname))[0])

        with open(file_base + "_first_idx_map_duplicate_idx_dict.json", "w") as json_file:
            json.dump(first_idx_map_duplicate_idx_dict, json_file, indent=4)
        print("Removing paired alpha and beta full-length TCR nt duplicates...")
        print("Shape before TCR duplicate removal:", tcr_refs_df.shape[0])
        tcr_refs_df.drop_duplicates("TCR_filter", inplace=True, keep="first")
        tcr_refs_df.reset_index(drop=True, inplace=True)
        print("Shape after TCR duplicate removal:", tcr_refs_df.shape[0])

    alpha_nt_ref_list = []
    alpha_cys104_pos_list = []
    # Check whether correct start and end amino acids are in alpha order after trimming the golden gate sites by 9 and -10:
    if not all([True if Seq(cdr3_seq).translate()[0] == "C" else False for cdr3_seq in tcr_refs_df[alpha_order_col_name].str[9:-10]]) & all(
        [True if Seq(cdr3_seq).translate()[-1] == "I" else False for cdr3_seq in tcr_refs_df[alpha_order_col_name].str[9:-10]]
    ):
        raise Exception("Alpha order column is in wrong format!")

    for tcr_idx in tcr_refs_df.index:
        tcr_refs_df.loc[tcr_idx, v_alpha_col] = tcr_refs_df.loc[tcr_idx, v_alpha_col].replace("/", "_")
        if "*" in tcr_refs_df.loc[tcr_idx, v_alpha_col]:
            # We only have *01 alleles in our V gene stocks. The + '*01' should be changed once we also have other alleles:
            if tcr_refs_df.loc[tcr_idx, v_alpha_col].split("*")[1] != "01":
                print("Warning: non-*01 allele detected in:", v_alpha_col)
            alpha_nt_ref_list.append(
                cys104_trimmed_v_gene_stock_nt_dict[tcr_refs_df.loc[tcr_idx, v_alpha_col].split("*")[0] + "*01"] + tcr_refs_df.loc[tcr_idx, alpha_order_col_name][9:-10]
            )
            alpha_cys104_pos_list.append(len(cys104_trimmed_v_gene_stock_nt_dict[tcr_refs_df.loc[tcr_idx, v_alpha_col].split("*")[0] + "*01"]))
        elif "*" not in tcr_refs_df.loc[tcr_idx, v_alpha_col]:
            alpha_nt_ref_list.append(cys104_trimmed_v_gene_stock_nt_dict[tcr_refs_df.loc[tcr_idx, v_alpha_col] + "*01"] + tcr_refs_df.loc[tcr_idx, alpha_order_col_name][9:-10])
            alpha_cys104_pos_list.append(len(cys104_trimmed_v_gene_stock_nt_dict[tcr_refs_df.loc[tcr_idx, v_alpha_col] + "*01"]))
        else:
            raise Exception("TRAV encoding is broken!")

        # Check whether reconstructed nt sequences are a multiple of 3:
    if not all([False if len(tcr_seq) % 3 != 0 else True for tcr_seq in alpha_nt_ref_list]):
        tmp_incorrect_len_alpha_dict = collections.defaultdict(str)
        for i, tcr_seq in enumerate(alpha_nt_ref_list):
            if len(tcr_seq) % 3 != 0:
                tmp_incorrect_len_alpha_dict[i] = tcr_seq
        raise Exception("Reconstructed alpha chain is not a multiple of 3!", tmp_incorrect_len_alpha_dict)

    if verbose == 1:
        print("Translation of first 5 reconstructed alpha chains:")
        print([str(Seq(tcr).translate()) for tcr in alpha_nt_ref_list[:5]])
        print("\n\n")

    beta_nt_ref_list = []
    beta_cys104_pos_list = []
    # Check whether correct start and end amino acids are in beta order after trimming the golden gate sites by 8 and -9:
    if not all([True if Seq(cdr3_seq).translate()[0] == "C" else False for cdr3_seq in tcr_refs_df[beta_order_col_name].str[8:-9]]) & all(
        [True if Seq(cdr3_seq).translate()[-1] == "E" else False for cdr3_seq in tcr_refs_df[beta_order_col_name].str[8:-9]]
    ):
        raise Exception("Beta order column is in wrong format!")

    for tcr_idx in tcr_refs_df.index:
        tcr_refs_df.loc[tcr_idx, v_beta_col] = tcr_refs_df.loc[tcr_idx, v_beta_col].replace("/", "_")
        if "*" in tcr_refs_df.loc[tcr_idx, v_beta_col]:
            if tcr_refs_df.loc[tcr_idx, v_beta_col].split("*")[1] != "01":
                print("Warning: non-*01 allele detected in:", v_beta_col)
            # We only have *01 alleles in our V gene stocks. The + '*01' should be changed once we also have other alleles:
            beta_nt_ref_list.append(
                cys104_trimmed_v_gene_stock_nt_dict[tcr_refs_df.loc[tcr_idx, v_beta_col].split("*")[0] + "*01"] + tcr_refs_df.loc[tcr_idx, beta_order_col_name][8:-9]
            )
            beta_cys104_pos_list.append(len(cys104_trimmed_v_gene_stock_nt_dict[tcr_refs_df.loc[tcr_idx, v_beta_col].split("*")[0] + "*01"]))
        elif "*" not in tcr_refs_df.loc[tcr_idx, v_beta_col]:
            beta_nt_ref_list.append(cys104_trimmed_v_gene_stock_nt_dict[tcr_refs_df.loc[tcr_idx, v_beta_col] + "*01"] + tcr_refs_df.loc[tcr_idx, beta_order_col_name][8:-9])
            beta_cys104_pos_list.append(len(cys104_trimmed_v_gene_stock_nt_dict[tcr_refs_df.loc[tcr_idx, v_beta_col] + "*01"]))
        else:
            raise Exception("TRBV encoding is broken!")

    # Check whether reconstructed nt sequences are a multiple of 3:
    if not all([False if len(tcr_seq) % 3 != 0 else True for tcr_seq in beta_nt_ref_list]):
        tmp_incorrect_len_beta_dict = collections.defaultdict(str)
        for i, tcr_seq in enumerate(beta_nt_ref_list):
            if len(tcr_seq) % 3 != 0:
                tmp_incorrect_len_beta_dict[i] = tcr_seq
        raise Exception("Reconstructed beta chain is not a multiple of 3!", tmp_incorrect_len_beta_dict)

    if verbose == 1:
        print("Translation of first 5 reconstructed beta chains:")
        print([str(Seq(tcr).translate()) for tcr in beta_nt_ref_list[:5]])
        print("\n\n")

    if not len(alpha_nt_ref_list) == len(beta_nt_ref_list):
        raise Exception("Alpha and beta nt reference lists became out of sync!")

    tcr_names = tcr_refs_df.loc[:, tcr_name_col_name]
    if not len(tcr_names) == len(alpha_nt_ref_list):
        raise Exception("tcr_refs_df became out of sync with alpha and beta nucleotide refs lists!")

    tcr_nt_refs_list = []
    for alpha_ref, beta_ref in zip(alpha_nt_ref_list, beta_nt_ref_list):
        if include_mu_constant_beta_and_p2a:
            tcr_nt_refs_list.append(beta_ref + tcr_golden_gate_constant_nt_seq_dict["muTRBC_P2A_nt"] + alpha_ref)
        else:
            tcr_nt_refs_list.append(beta_ref + alpha_ref)

    if verbose == 1:
        print("Translation of first 5 final TCR references:")
        print([str(Seq(tcr).translate()) for tcr in tcr_nt_refs_list[:5]])
        print("\n\n")

    if add_number_of_negatives_by_mutating_full_refs or add_number_of_negatives_by_mutating_cdr3j_in_refs or add_number_of_negatives_by_mutating_v_gene_in_refs:
        tcr_names = tcr_names.to_numpy()
        number_non_negative_tcrs = len(tcr_nt_refs_list)

    if add_number_of_negatives_by_mutating_full_refs:
        if add_number_of_negatives_by_mutating_full_refs > number_non_negative_tcrs:
            raise Exception("More negative TCR ref sequences cannot be generated than TCR ref sequences by mutating full TCR refs!")

        tcr_ref_sample_idxs = rng.choice(number_non_negative_tcrs, size=add_number_of_negatives_by_mutating_full_refs, replace=False)
        for tcr_seq, name in zip(np.array(tcr_nt_refs_list)[tcr_ref_sample_idxs], tcr_names[tcr_ref_sample_idxs]):
            tcr_seq_arr = np.array(list(tcr_seq))
            pos_to_change_arr = rng.choice(len(tcr_seq_arr), size=rng.choice(np.arange(min_nt_diff_negative_ref_seqs, max_nt_diff_negative_ref_seqs + 1)), replace=False)
            for pos in pos_to_change_arr:
                tcr_seq_arr[pos] = rng.choice(list("ATGC"))
            tcr_neg_seq = "".join(tcr_seq_arr)

            if tcr_neg_seq not in tcr_nt_refs_list and (name + "_full_n") not in tcr_names:
                tcr_nt_refs_list.append(tcr_neg_seq)
                tcr_names = np.append(tcr_names, (name + "_full_n"))

    if add_number_of_negatives_by_mutating_cdr3j_in_refs:
        if add_number_of_negatives_by_mutating_cdr3j_in_refs > number_non_negative_tcrs:
            raise Exception("More negative TCR ref sequences cannot be generated than TCR ref sequences by mutating CDR3-J region in TCR refs!")
        muTRBC_len = len(tcr_golden_gate_constant_nt_seq_dict["muTRBC_P2A_nt"])
        tcr_ref_sample_idxs = rng.choice(number_non_negative_tcrs, size=add_number_of_negatives_by_mutating_cdr3j_in_refs, replace=False)

        for tcr_seq, name, beta_cys104_pos, alpha_cys104_pos, alpha_ref, beta_ref in zip(
            np.array(tcr_nt_refs_list)[tcr_ref_sample_idxs],
            tcr_names[tcr_ref_sample_idxs],
            np.array(beta_cys104_pos_list)[tcr_ref_sample_idxs],
            np.array(alpha_cys104_pos_list)[tcr_ref_sample_idxs],
            np.array(alpha_nt_ref_list)[tcr_ref_sample_idxs],
            np.array(beta_nt_ref_list)[tcr_ref_sample_idxs],
        ):
            tcr_seq_arr = np.array(list(tcr_seq))
            alpha_bool = rng.choice([True, False])
            if alpha_bool:
                alpha_start = len(beta_ref) + muTRBC_len
                low = alpha_start + alpha_cys104_pos
                high = alpha_start + len(alpha_ref)
                pos_to_change_arr = rng.choice(np.arange(low, high), size=rng.choice(np.arange(min_nt_diff_negative_ref_seqs, max_nt_diff_negative_ref_seqs + 1)), replace=False)
                for pos in pos_to_change_arr:
                    tcr_seq_arr[pos] = rng.choice(list("ATGC"))
                tcr_neg_seq = "".join(tcr_seq_arr)

                if tcr_neg_seq not in tcr_nt_refs_list and (name + "_cdr3j_n") not in tcr_names:
                    tcr_nt_refs_list.append(tcr_neg_seq)
                    tcr_names = np.append(tcr_names, (name + "_cdr3j_n"))
            else:
                low = beta_cys104_pos
                high = len(beta_ref)
                pos_to_change_arr = rng.choice(np.arange(low, high), size=rng.choice(np.arange(min_nt_diff_negative_ref_seqs, max_nt_diff_negative_ref_seqs + 1)), replace=False)
                for pos in pos_to_change_arr:
                    tcr_seq_arr[pos] = rng.choice(list("ATGC"))
                tcr_neg_seq = "".join(tcr_seq_arr)

                if tcr_neg_seq not in tcr_nt_refs_list and (name + "_cdr3j_n") not in tcr_names:
                    tcr_nt_refs_list.append(tcr_neg_seq)
                    tcr_names = np.append(tcr_names, (name + "_cdr3j_n"))

    if add_number_of_negatives_by_mutating_v_gene_in_refs:
        if add_number_of_negatives_by_mutating_cdr3j_in_refs > number_non_negative_tcrs:
            raise Exception("More negative TCR ref sequences cannot be generated than TCR ref sequences by mutating V gene region in TCR refs!")

        muTRBC_len = len(tcr_golden_gate_constant_nt_seq_dict["muTRBC_P2A_nt"])
        tcr_ref_sample_idxs = rng.choice(number_non_negative_tcrs, size=add_number_of_negatives_by_mutating_v_gene_in_refs, replace=False)

        for tcr_seq, name, beta_cys104_pos, alpha_cys104_pos, alpha_ref, beta_ref in zip(
            np.array(tcr_nt_refs_list)[tcr_ref_sample_idxs],
            tcr_names[tcr_ref_sample_idxs],
            np.array(beta_cys104_pos_list)[tcr_ref_sample_idxs],
            np.array(alpha_cys104_pos_list)[tcr_ref_sample_idxs],
            np.array(alpha_nt_ref_list)[tcr_ref_sample_idxs],
            np.array(beta_nt_ref_list)[tcr_ref_sample_idxs],
        ):
            tcr_seq_arr = np.array(list(tcr_seq))
            alpha_bool = rng.choice([True, False])
            if alpha_bool:
                alpha_start = len(beta_ref) + muTRBC_len
                low = alpha_start
                high = alpha_start + alpha_cys104_pos
                pos_to_change_arr = rng.choice(np.arange(low, high), size=rng.choice(np.arange(min_nt_diff_negative_ref_seqs, max_nt_diff_negative_ref_seqs + 1)), replace=False)
                for pos in pos_to_change_arr:
                    tcr_seq_arr[pos] = rng.choice(list("ATGC"))
                tcr_neg_seq = "".join(tcr_seq_arr)

                if tcr_neg_seq not in tcr_nt_refs_list and (name + "_v_n") not in tcr_names:
                    tcr_nt_refs_list.append(tcr_neg_seq)
                    tcr_names = np.append(tcr_names, (name + "_v_n"))
            else:
                low = 0
                high = beta_cys104_pos
                pos_to_change_arr = rng.choice(np.arange(low, high), size=rng.choice(np.arange(min_nt_diff_negative_ref_seqs, max_nt_diff_negative_ref_seqs + 1)), replace=False)
                for pos in pos_to_change_arr:
                    tcr_seq_arr[pos] = rng.choice(list("ATGC"))
                tcr_neg_seq = "".join(tcr_seq_arr)

                if tcr_neg_seq not in tcr_nt_refs_list and (name + "_v_n") not in tcr_names:
                    tcr_nt_refs_list.append(tcr_neg_seq)
                    tcr_names = np.append(tcr_names, (name + "_v_n"))

    bed_out_fname = os.path.join(os.path.dirname(fasta_out_fname), os.path.splitext(os.path.basename(fasta_out_fname))[0] + ".bed")
    with open(fasta_out_fname, mode="w") as fasta_out, open(bed_out_fname, mode="w") as bed_out:
        for name, tcr_seq in zip(tcr_names, tcr_nt_refs_list):
            fasta_out.write(">" + str(name) + "\n")
            fasta_out.write(tcr_seq + "\n")
            bed_out.write(str(name) + "\t" + str(0) + "\t" + str(len(tcr_seq)) + "\t" + str(name) + "\n")
