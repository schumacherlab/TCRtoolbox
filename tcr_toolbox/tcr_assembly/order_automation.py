import collections
import datetime
import glob
import json
import math
import os
import re
from re import search
from typing import Literal, Union

import numpy as np
import pandas as pd
from Bio import Restriction, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from dotenv import load_dotenv
from pydna.amplify import Anneal
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer

from tcr_toolbox.tcr_assembly.constants import (
    alt_v_genes_beta,
    insert_alpha,
    insert_beta,
    j_gene_res119_rest_orphan_dict_alpha,
    j_gene_res119_rest_orphan_dict_beta,
    sites_to_check,
    snapgene_feature_dict,
)
from tcr_toolbox.utils.codon_optimizer import codon_optimize, make_search_strings
from tcr_toolbox.utils.constants import aa_to_all_codons_list
from tcr_toolbox.utils.logger import log_function_to_file
from tcr_toolbox.utils.utils import (
    generate_random_dna_sequence,
    levenshtein_distance,
    levenshtein_ratio,
    write_echo_dispense_csv,
    write_hamilton_pipetting_excel_sheet,
    write_idot_dispense_csv,
    add_well_coordinates,
)
import warnings
from Bio import BiopythonDeprecationWarning

warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)


load_dotenv()
tcr_toolbox_data_path = os.getenv("tcr_toolbox_data_path")
v_gene_barcode_tracing_path = os.getenv("v_gene_barcode_tracing_path")


def return_assembly_run_number(run_dirs_path: str | os.PathLike, run_name: str):
    """
    Generate a new, uniquely numbered run directory name and path for an assembly run.

    This function scans an existing directory containing run folders (named in the
    format `r<number>_<run_name>`), determines the next available run number,
    and constructs a new directory name and full path. It ensures that the provided
    `run_name` does not already exist among existing runs.

    Parameters
    ----------
    run_dirs_path : str or os.PathLike
        Path to the parent directory containing previous run directories.
    run_name : str
        Descriptive name for the new run. Must not duplicate an existing run name.

    Returns
    -------
    numbered_run_name : str
        The new run directory name in the format `r<number>_<run_name>`.
    run_path : str
        The full path to the new run directory.

    Raises
    ------
    ValueError
        If a run with the given `run_name` already exists in `run_dirs_path`.

    Notes
    -----
    - Existing runs must follow the naming convention `r<number>_<run_name>`.
    - The function prints the new run number to stdout.
    """

    run_dir_list = [file for file in os.listdir(run_dirs_path) if os.path.isdir(os.path.join(run_dirs_path, file))]

    run_name_list = [re.split(r"r\d+_", run_dir, maxsplit=1)[-1] for run_dir in run_dir_list]

    if run_name in run_name_list:
        raise ValueError(f"Run name '{run_name}' already exists! Cannot assign a new run number.")

    run_number_list = [int(re.match(r"r(\d+)_", run_dir).group(1)) for run_dir in run_dir_list if re.match(r"r\d+_", run_dir)]
    next_run_number = max(run_number_list, default=0) + 1

    numbered_run_name = f"r{next_run_number}_{run_name}"
    run_path = os.path.join(run_dirs_path, numbered_run_name)

    print("Your run_number is:", f"r{next_run_number}")
    return numbered_run_name, run_path


def init_oligo_order_run_dir(run_name: str):
    """
    Initialize a new oligo order (TCR assembly) run directory and its required subdirectories.

    This function creates a new, uniquely numbered directory for a TCR assembly run
    inside the `tcr_toolbox_tcr_assembly_runs` folder located under the
    `tcr_toolbox_data_path` environment variable. It also creates the required
    subdirectory structure used throughout the oligo ordering and TCR
    assembly pipeline.

    Each assembly run corresponds to a set of TCRs for which V genes are
    to be pre-mixed together in a single V gene dispensing run.

    Parameters
    ----------
    run_name : str
        Descriptive name for the new oligo order run. Must not duplicate an existing run name.

    Returns
    -------
    numbered_run_name : str
        The new run directory name in the format `r<number>_<run_name>`.
    run_path : str
        The absolute path to the created run directory.

    Raises
    ------
    EnvironmentError
        If the environment variable `tcr_toolbox_data_path` is not set.
    FileNotFoundError
        If the expected `tcr_toolbox_tcr_assembly_runs` directory does not exist.
    FileExistsError
        If a directory for this run already exists.

    Notes
    -----
    The following subdirectories are created under the new run directory:
        - `cdr3j_oligo_order_sheets`
        - `plate_sheets`
        - `v_genes_premix_dispense/`
            - `v_gene_plate_maps/`
            - `v_gene_stock/`
        - `pydna_logs`
        - `pydna_ligation_dicts`
    """
    if not tcr_toolbox_data_path:
        raise EnvironmentError(
            "Environment variable 'tcr_toolbox_data_path' is not set. Did you create or edit your .env file with the correct path?"
        )

    run_dirs_path = os.path.join(tcr_toolbox_data_path, "tcr_toolbox_tcr_assembly_runs")

    if not os.path.isdir(run_dirs_path):
        raise FileNotFoundError(
            f"{run_dirs_path} does not exist!\n"
            "Possible causes:\n"
            "  1) You did not download the accompanying tcr_toolbox data package.\n"
            "  2) The 'tcr_toolbox_data_path' in your .env file is incorrect.\n"
            "  3) You renamed the 'tcr_toolbox_tcr_assembly_run directory.\n"
            "     Directories in the tcr_toolbox data package should not be renamed!"
            f"Current 'tcr_toolbox_data_path': {tcr_toolbox_data_path}"
        )

    numbered_run_name, run_path = return_assembly_run_number(run_dirs_path=run_dirs_path, run_name=run_name)
    if os.path.isdir(run_path):
        raise FileExistsError(f"{run_path} already exist!\n\nDid you already prepare the assembly for this TCR library?")
    os.mkdir(run_path)
    os.mkdir(os.path.join(run_path, "cdr3j_oligo_order_sheets"))
    os.mkdir(os.path.join(run_path, "plate_sheets"))
    os.mkdir(os.path.join(run_path, "v_genes_premix_dispense"))  # will be removed after running pydna simulation
    os.mkdir(os.path.join(run_path, "v_genes_premix_dispense", "v_gene_plate_maps"))  # will be removed after running pydna simulation
    os.mkdir(os.path.join(run_path, "v_genes_premix_dispense", "v_gene_stock"))  # will be removed after running pydna simulation
    os.mkdir(os.path.join(run_path, "pydna_logs"))
    os.mkdir(os.path.join(run_path, "pydna_ligation_dicts"))

    return numbered_run_name, run_path


def validate_oligo_order_run_dir(run_path: str | os.PathLike):
    """Check that a run_path is a valid standardized oligo_order run directory."""
    required_subdirs = ["cdr3j_oligo_order_sheets", "plate_sheets", "pydna_ligation_dicts", "pydna_logs"]

    if not os.path.isdir(run_path):
        raise FileNotFoundError(
            f"{run_path} does not exist!\n\n"
            "To run the full assembly, you first need to execute `run_mode = oligo_order`."
            "This step will perform the following:\n\n"
            "1. Prepare the CDR3-J oligo pool order.\n"
            "2. Identify which V gene stock tubes need to be replaced if they have too little remaining volume "
            "while the oligo pool is being manufactured.\n\n"
            "Please complete these steps before running `run_mode = run_assembly`."
            "It could also be that run_path does not exist because you aciddentally double quoted in the input json."
        )

    missing = [sub for sub in required_subdirs if not os.path.isdir(os.path.join(run_path, sub))]
    if missing:
        raise FileNotFoundError(f"The following required subdirectories are missing from {run_path}:\n  " + "\n  ".join(missing))

    return True


def expand_oligo_order_run_dir_to_run_assembly(run_path: str | os.PathLike):
    """
    Expand an existing oligo order run directory into a full TCR assembly run directory.

    This function adds the additional subdirectories required for performing a full
    TCR assembly and sequencing quality analysis. It assumes that the provided
    directory path corresponds to a valid oligo order run created by
    `init_oligo_order_run_dir`.

    The newly added directories include structures needed for V gene pre-mix
    dispensing and sequencing quality analysis.

    Parameters
    ----------
    run_path : str or os.PathLike
        Path to the oligo order run directory to be expanded.

    Returns
    -------
    numbered_run_name : str
        The name of the expanded run directory (e.g., `r5_my_tcr_run`).
    run_path : str
        The absolute path to the expanded run directory.

    Raises
    ------
    FileExistsError
        If any of the new directories to be created already exist.
    FileNotFoundError
        If the specified `run_path` does not exist.

    Notes
    -----
    The following subdirectories are created under `run_path`:
        - `v_genes_premix_dispense/`
            - `v_gene_plate_maps/`
            - `v_gene_stock/`
        - `sequencing_quality_analysis/`
            - `references/`
    """
    os.mkdir(os.path.join(run_path, "v_genes_premix_dispense"))
    os.mkdir(os.path.join(run_path, "v_genes_premix_dispense", "v_gene_plate_maps"))
    os.mkdir(os.path.join(run_path, "v_genes_premix_dispense", "v_gene_stock"))
    os.mkdir(os.path.join(run_path, "sequencing_quality_analysis"))
    os.mkdir(os.path.join(run_path, "sequencing_quality_analysis", "references"))

    numbered_run_name = os.path.basename(run_path)

    return numbered_run_name, run_path


def translate_v_gene_input_to_imgt(
    run_tcr_df: pd.DataFrame,
    v_alpha_col: str,
    v_beta_col: str,
    vdj_translation_dict_json: Union[str, os.PathLike[str]] = os.path.join(
        tcr_toolbox_data_path,
        "tcr_toolbox_datasets",
        "tcr_reconstruction",
        "VDJ_gene_sequences",
        "after_benchmark",
        "functional_with_L-PART1+V-EXON_after_benchmark",
        "20230803_vdj_translation_dict_aa.json",
    ),
) -> pd.DataFrame:
    """
    Translate V gene identifiers in a TCR dataset to standardized IMGT nomenclature.

    This function converts the TRAV (V-alpha) and TRBV (V-beta) gene identifiers
    in the input DataFrame to their corresponding IMGT-standardized names using
    a provided translation dictionary in JSON format.

    It also validates that all V genes can be translated, raising an error if
    unrecognized entries are found. If translations alter gene names (beyond
    allele differences), concise logs are printed summarizing the changes.

    Parameters
    ----------
    run_tcr_df : pandas.DataFrame
        Input DataFrame containing TCR records with V-alpha and V-beta columns.
    v_alpha_col : str
        Name of the column in `run_tcr_df` containing TRAV (V-alpha) gene identifiers.
    v_beta_col : str
        Name of the column in `run_tcr_df` containing TRBV (V-beta) gene identifiers.
    vdj_translation_dict_json : str or os.PathLike, optional
        Path to the JSON file containing the VDJ translation dictionary.
        Defaults to:
        `tcr_toolbox_datasets/tcr_reconstruction/VDJ_gene_sequences/after_benchmark/functional_with_L-PART1+V-EXON_after_benchmark/20230803_vdj_translation_dict_aa.json`

    Returns
    -------
    pandas.DataFrame
        A copy of the input DataFrame with the TRAV and TRBV columns
        translated to IMGT-standardized gene names.

    Raises
    ------
    ValueError
        If any TRAV or TRBV gene cannot be translated using the provided dictionary.
    FileNotFoundError
        If the specified JSON translation dictionary file does not exist.
    """

    run_tcr_df = run_tcr_df.copy()
    run_tcr_df[v_alpha_col] = run_tcr_df[v_alpha_col].str.replace("_", "/").str.strip()
    run_tcr_df[v_beta_col] = run_tcr_df[v_beta_col].str.replace("_", "/").str.strip()

    with open(vdj_translation_dict_json, "r") as f:
        vdj_translation_dict_aa = json.load(f)

    trav_dict = {k: v[0] for k, v in vdj_translation_dict_aa["TRAV"].items()}
    trbv_dict = {k: v[0] for k, v in vdj_translation_dict_aa["TRBV"].items()}

    original_v_alpha_series = run_tcr_df.loc[:, v_alpha_col].copy()
    original_v_beta_series = run_tcr_df.loc[:, v_beta_col].copy()

    run_tcr_df[v_alpha_col] = run_tcr_df[v_alpha_col].map(trav_dict)
    run_tcr_df[v_beta_col] = run_tcr_df[v_beta_col].map(trbv_dict)

    if run_tcr_df[v_alpha_col].isna().any():
        missing = original_v_alpha_series.loc[run_tcr_df[v_alpha_col].isna()].unique()
        raise ValueError(f"Untranslatable TRAV genes found: {missing}")
    if run_tcr_df[v_beta_col].isna().any():
        missing = original_v_beta_series.loc[run_tcr_df[v_beta_col].isna()].unique()
        raise ValueError(f"Untranslatable TRBV genes found: {missing}")

    fixed_v_alpha_dict = collections.defaultdict(set)
    fixed_v_beta_dict = collections.defaultdict(set)

    for orig, trans in zip(original_v_alpha_series, run_tcr_df[v_alpha_col]):
        if orig != trans:
            # Skip allele-only differences
            if orig.split("*")[0] == trans.split("*")[0]:
                continue
            fixed_v_alpha_dict[trans].add(orig)

    for orig, trans in zip(original_v_beta_series, run_tcr_df[v_beta_col]):
        if orig != trans:
            # Skip allele-only differences
            if orig.split("*")[0] == trans.split("*")[0]:
                continue
            fixed_v_beta_dict[trans].add(orig)

    # --- Print concise logs ---
    if fixed_v_alpha_dict:
        print("TRAV genes translated to IMGT:")
        for k, v in fixed_v_alpha_dict.items():
            print(f"  {k}: {sorted(v)}")
    if fixed_v_beta_dict:
        print("TRBV genes translated to IMGT:")
        for k, v in fixed_v_beta_dict.items():
            print(f"  {k}: {sorted(v)}")

    return run_tcr_df


def add_collapsed_v_alleles_to_01_col(tcr_df: pd.DataFrame, v_alpha_col: str, v_beta_col: str):
    tcr_df["TRAV_IMGT_allele_collapsed"] = tcr_df[v_alpha_col].str.split("*").str[0] + "*01"
    tcr_df["TRBV_IMGT_allele_collapsed"] = tcr_df[v_beta_col].str.split("*").str[0] + "*01"
    return tcr_df


def list_v_alpha_and_v_beta_dna_files_dir(
    v_gene_dir: Union[str, os.PathLike[str]] = os.path.join(
        tcr_toolbox_data_path,
        "tcr_toolbox_datasets",
        "tcr_assembly",
        "plasmids",
        "TRV_plasmid_stocks_ordered_at_Twist",
    ),
):
    v_alphas_list = [os.path.basename(v_gene_file).split(".dna")[0] for v_gene_file in glob.glob(v_gene_dir + "/TRAV*.dna")]
    v_betas_list = [os.path.basename(v_gene_file).split(".dna")[0] for v_gene_file in glob.glob(v_gene_dir + "/TRBV*.dna")]

    if not set([v_alpha.split("*")[1] for v_alpha in v_alphas_list]) == {"01"}:
        raise NotImplementedError("Non-*01 allele detected in V alpha gene stock plasmid folder. Should we implement non-*01 alleles?")
    if not set([v_beta.split("*")[1] for v_beta in v_betas_list]) == {"01"}:
        raise NotImplementedError("Non-*01 allele detected in V beta gene stock plasmid folder. Should we implement non-*01 alleles?")

    return v_alphas_list, v_betas_list


def filter_tcr_df_on_v_genes_in_stock(
    run_path: Union[str, os.PathLike[str]], numbered_run_name: str, tcr_df: pd.DataFrame, v_alpha_col: str, v_beta_col: str, tcr_df_name: str = ""
):
    """Remove TCRs from TCR DataFrame of which we do not have their V beta and/or alpha genes in stock. Write .xlsx
    files that contain which V beta and alpha genes and alleles we do not have in our V gene vector stock.


    Parameters
    ----------
    run_path : Union[str, os.PathLike[str]]
        String path to the standardised assembly run directory.
    numbered_run_name : str
        Name of the assembly run.
    tcr_df : pd.DataFrame
        DataFrame containing TCRs that need to be filtered. For example, a re-formated (e.g., V genes names re-formatted
        to IMGT format) 10X Genomics VDJ output file or MiXCR output file. tcr_df needs to contain at least a CDR3 alpha
        and beta amino acid sequence, TRAV, TRBV, TRAJ, and TRBJ annotation.
    v_alpha_col : str
        Column that stores TCR V alpha gene + allele. V alpha gene + allele needs to be in IMGT format. Currently, only
        *01 alleles are implemented because we do not have *02 alleles in our V alpha gene vector stock. If you have
        V alpha gene alleles other than *01, a solution is currently to collapse these non-*01 alleles to allele *01.
    v_beta_col : str
        Column that stores TCR V beta gene + allele. V beta gene + allele needs to be in IMGT format. Currently, only
        *01 alleles are implemented because we do not have *02 alleles in our V beta gene vector stock. If you have
        V beta gene alleles other than *01, a solution is currently to collapse these non-*01 alleles to allele *01.
    tcr_df_name : str = ''
        Unique name for writing TCR DataFrame .xlsx files that contain which V beta and alpha genes and alleles we do
        not have in our vector stock. Required if there are separate .csv files that store the TCRs that need to be
        processed together in an assembly run, as for each TCR DataFrame we want to separately write in the
        assembly run directory a uniquely named .xlsx file.

    Returns
    -------
    tcr_df : pd.DataFrame
        Filtered TCR DataFrame that only contains TCRs of which we have their V genes in our V gene vector stock.

    Notes
    -----
    Function only works when v_alpha_col and v_beta_col are in IMGT format.
    """

    v_alphas_list, v_betas_list = list_v_alpha_and_v_beta_dna_files_dir(
        v_gene_dir=os.path.join(
            tcr_toolbox_data_path,
            "tcr_toolbox_datasets",
            "tcr_assembly",
            "plasmids",
            "TRV_plasmid_stocks_ordered_at_Twist",
        )
    )

    tcr_df[v_alpha_col] = tcr_df[v_alpha_col].str.replace("/", "_").str.strip()
    tcr_df[v_beta_col] = tcr_df[v_beta_col].str.replace("/", "_").str.strip()

    v_alphas_alleles_not_in_stock = pd.DataFrame(tcr_df[~tcr_df[v_alpha_col].isin(v_alphas_list)][v_alpha_col].value_counts())
    v_alphas_alleles_not_in_stock.rename({v_alpha_col: "count"}, axis=1, inplace=True)

    v_betas_alleles_not_in_stock = pd.DataFrame(tcr_df[~tcr_df[v_beta_col].isin(v_betas_list)][v_beta_col].value_counts())
    v_betas_alleles_not_in_stock.rename({v_beta_col: "count"}, axis=1, inplace=True)

    tcr_df = add_collapsed_v_alleles_to_01_col(tcr_df=tcr_df, v_alpha_col=v_alpha_col, v_beta_col=v_beta_col)

    v_alphas_not_in_stock = pd.DataFrame(tcr_df[~tcr_df["TRAV_IMGT_allele_collapsed"].isin(v_alphas_list)]["TRAV_IMGT_allele_collapsed"].value_counts())
    v_alphas_not_in_stock.rename({"TRAV_IMGT_allele_collapsed": "count"}, axis=1, inplace=True)

    v_betas_not_in_stock = pd.DataFrame(tcr_df[~tcr_df["TRBV_IMGT_allele_collapsed"].isin(v_betas_list)]["TRBV_IMGT_allele_collapsed"].value_counts())
    v_betas_not_in_stock.rename({"TRBV_IMGT_allele_collapsed": "count"}, axis=1, inplace=True)

    base_path = os.path.join(run_path, "v_genes_premix_dispense", "v_gene_stock")
    if tcr_df_name:
        suffix = f"{numbered_run_name}_{tcr_df_name}_"
    else:
        suffix = f"{numbered_run_name}_"

    v_alphas_alleles_not_in_stock.to_excel(os.path.join(base_path, f"{suffix}v_alphas_alleles_not_in_stock.xlsx"))
    v_betas_alleles_not_in_stock.to_excel(os.path.join(base_path, f"{suffix}v_betas_alleles_not_in_stock.xlsx"))
    v_alphas_not_in_stock.to_excel(os.path.join(base_path, f"{suffix}v_alphas_not_in_stock.xlsx"))
    v_betas_not_in_stock.to_excel(os.path.join(base_path, f"{suffix}v_betas_not_in_stock.xlsx"))

    print("Number of TCRs before filtering:", tcr_df.shape[0])
    tcr_df = tcr_df[tcr_df["TRAV_IMGT_allele_collapsed"].isin(v_alphas_list)].copy()
    tcr_df = tcr_df[tcr_df["TRBV_IMGT_allele_collapsed"].isin(v_betas_list)].copy()
    print("Number of TCRs after filtering:", tcr_df.shape[0])

    tcr_df.reset_index(inplace=True, drop=True)

    return tcr_df


def count_run_v_gene_usage(run_path: Union[str, os.PathLike[str]], numbered_run_name: str, run_tcr_df: pd.DataFrame):
    """
    Only works when v_alpha_col and v_beta_col are in IMGT format.
    """

    base_path = os.path.join(run_path, "v_genes_premix_dispense", "v_gene_stock")
    v_alphas_usage_df = pd.DataFrame(run_tcr_df["TRAV_IMGT_allele_collapsed"].value_counts()).copy()
    v_alphas_usage_df.rename({"TRAV_IMGT_allele_collapsed": "count"}, axis=1, inplace=True)
    v_alphas_usage_df.to_excel(os.path.join(base_path, f"{numbered_run_name}_v_alphas_usage_count.xlsx"))

    v_betas_usage_df = pd.DataFrame(run_tcr_df["TRBV_IMGT_allele_collapsed"].value_counts()).copy()
    v_betas_usage_df.rename({"TRBV_IMGT_allele_collapsed": "count"}, axis=1, inplace=True)
    v_betas_usage_df.to_excel(os.path.join(base_path, f"{numbered_run_name}_v_betas_usage_count.xlsx"))

    return v_alphas_usage_df, v_betas_usage_df


def concatenate_alpha(cdr3_alpha: str, j_gene_alpha: str, idx: int, filter_cannot_be_codon_optimized: bool = True, steepness: float = 4.0, idx_seed: int = None):
    """Concatenate alpha subunits into a single string."""
    if idx_seed is None:
        idx_seed = idx

    # remove Cysteine 104 from CDR3
    if cdr3_alpha[0] != "C":
        raise ValueError(f"First amino acid in CDR3 is not Cysteine for row {idx, cdr3_alpha}.")

    # remove C
    cdr3_alpha_no_c = cdr3_alpha[1:]
    # add J119 + rest + orphan residue
    cdr3_alpha_j119 = cdr3_alpha_no_c + j_gene_res119_rest_orphan_dict_alpha[j_gene_alpha]

    # codon optimize
    cdr3_alpha_j119_nt_opt = codon_optimize(sequence_aa=cdr3_alpha_j119, enzymes=["BbsI"], sites_to_check=sites_to_check, if_fail="return_none", idx=idx_seed, steepness=steepness)

    if cdr3_alpha_j119_nt_opt is None:
        if filter_cannot_be_codon_optimized:
            print(f"{cdr3_alpha} CDR3-J alpha aa sequence of TCR: {idx} cannot be codon optimized!")
            # print(f"Removing TCR: {idx}\n")
        else:
            print(f"{cdr3_alpha} CDR3-J alpha aa sequence of TCR: {idx} cannot be codon optimized!")
            print(f"As filter_cannot_be_codon_optimized is not set to True, TCR: {idx} will not be removed!\n")
        return filter_cannot_be_codon_optimized, None

    # 2nt spacer, Overhang (Cys104+preceding base)
    prepend = insert_alpha["2nt_spacer_1"] + insert_alpha["cysteine_overhang"]
    # muTRAC partial,	2nt spacer
    append = insert_alpha["mu_trac_patial_overhang"] + insert_alpha["2nt_spacer_2"]
    full_sequence = prepend + cdr3_alpha_j119_nt_opt + append

    # check for restriction sites
    search_string, restriction_enzyme_sites, restriction_edge_sites = make_search_strings(enzymes=["BbsI"])
    if search(search_string, full_sequence):
        if filter_cannot_be_codon_optimized:
            print(f"{cdr3_alpha} CDR3-J alpha aa sequence of TCR: {idx} contains restriction site(s)!")
            print(f"Removing TCR: {idx}\n")
        else:
            print(f"{cdr3_alpha} CDR3-J alpha aa sequence of TCR: {idx} contains restriction site(s)!")
            print(f"As filter_cannot_be_codon_optimized is not set to True, TCR: {idx} will not be removed!\n")
        return filter_cannot_be_codon_optimized, None

    full_sequence_with_restriction = insert_alpha["bbsi_1"] + full_sequence + insert_alpha["bbsi_2"]

    return False, full_sequence_with_restriction


def G_is_succeeding_nt_cys104_beta(cdr3_beta_aa: str):
    if cdr3_beta_aa[0] != "C":
        raise ValueError(f"First amino acid in CDR3 is not Cysteine for row {cdr3_beta_aa}.")

    return any([codon[0] == "G" for codon in aa_to_all_codons_list[cdr3_beta_aa[1]]])


def A_is_succeeding_nt_cys104_beta(cdr3_beta_aa: str):
    if cdr3_beta_aa[0] != "C":
        raise ValueError(f"First amino acid in CDR3 is not Cysteine for row {cdr3_beta_aa}.")

    return any([codon[0] == "A" for codon in aa_to_all_codons_list[cdr3_beta_aa[1]]])


def cys104_beta_can_be_cloned(cdr3_beta_aa: str):
    return G_is_succeeding_nt_cys104_beta(cdr3_beta_aa) or A_is_succeeding_nt_cys104_beta(cdr3_beta_aa)


def return_optimal_codon_starting_with_G_nt(aa: str):
    codon_list = [codon for codon in aa_to_all_codons_list[aa] if codon[0] == "G"]
    if not codon_list:
        raise Exception("Amino acid:", aa, "is not encoded by codons that start with an G nucleotide!")

    # We can select the first G nt starting codon as the most optimal codon,
    # as codons are ordered on their usage frequency:
    codon = codon_list[0]
    return codon


def return_optimal_codon_starting_with_A_nt(aa: str):
    codon_list = [codon for codon in aa_to_all_codons_list[aa] if codon[0] == "A"]
    if not codon_list:
        raise Exception("Amino acid:", aa, "is not encoded by codons that start with an A nucleotide!")

    # We can select the first A nt starting codon as the most optimal codon,
    # as codons are ordered on their usage frequency:
    codon = codon_list[0]
    return codon


def concatenate_beta(
    cdr3_beta: str,
    v_gene_beta: str,
    j_gene_beta: str,
    idx: int,
    filter_cannot_be_codon_optimized: bool = True,
    filter_succeeding_nt_cys_104_beta: bool = True,
    steepness: float = 4.0,
    idx_seed: int = None,
):
    """Concatenate beta subunits into a single string."""
    if idx_seed is None:
        idx_seed = idx

    append_is_alt_v_beta_gene = "_alt" if v_gene_beta in alt_v_genes_beta else ""
    alt_v_beta_gene_bool = v_gene_beta in alt_v_genes_beta

    # remove Cysteine 104 from CDR3
    if cdr3_beta[0] != "C":
        raise ValueError(f"First amino acid in CDR3 is not Cysteine for row {idx, cdr3_beta}.")

    # add J119 + rest + orphan residue
    cdr3_beta_j119 = cdr3_beta + j_gene_res119_rest_orphan_dict_beta[j_gene_beta]
    # for 4-base overhang cloning, we do not want to codon optimize the first two codons and the last codon
    # remove first two and last amino acid from cdr3_beta_j119
    cdr3_beta_j119 = cdr3_beta_j119[2:-1]

    cdr3_beta_j119_nt_opt = codon_optimize(sequence_aa=cdr3_beta_j119, enzymes=["BbsI"], sites_to_check=sites_to_check, if_fail="return_none", idx=idx_seed, steepness=steepness)

    if cdr3_beta_j119_nt_opt is None:
        if filter_cannot_be_codon_optimized:
            print(f"{cdr3_beta} CDR3-J beta aa sequence of TCR: {idx} cannot be codon optimized!")
            # print(f"Removing TCR: {idx}\n")
        else:
            print(f"{cdr3_beta} CDR3-J beta aa sequence of TCR: {idx} cannot be codon optimized!")
            print(f"As filter_cannot_be_codon_optimized is not set to True, TCR: {idx} will not be removed!\n")
        return filter_cannot_be_codon_optimized, None

    succeeding_g_nt = G_is_succeeding_nt_cys104_beta(cdr3_beta)
    succeeding_a_nt = A_is_succeeding_nt_cys104_beta(cdr3_beta)
    if succeeding_g_nt and succeeding_a_nt:
        raise Exception(
            f"Either G_is_succeeding_nt_cys104_beta {succeeding_g_nt} or A_is_succeeding_nt_cys104_beta {succeeding_a_nt} function is broken, as both return true for this CDR3 beta!"
        )

    elif succeeding_g_nt:
        if not append_is_alt_v_beta_gene:
            second_codon_starting_with_G = return_optimal_codon_starting_with_G_nt(aa=cdr3_beta[1])
            prepend = insert_beta["bbsi_1"] + insert_beta["2nt_spacer_1"] + insert_beta["cysteine_overhang"] + second_codon_starting_with_G[1:]
        else:
            if filter_succeeding_nt_cys_104_beta:
                print(
                    cdr3_beta,
                    f"CDR3 beta Cys104 succeeding nucleotide is G {succeeding_g_nt} or A {succeeding_a_nt}\nbut TRBV20-1 or TRBV29-1 is {alt_v_beta_gene_bool}",
                    v_gene_beta,
                    "\nRemoving TCR:",
                    idx,
                    "\n",
                )
                return filter_succeeding_nt_cys_104_beta, None
            else:
                print(
                    cdr3_beta,
                    f"CDR3 beta Cys104 succeeding nucleotide is G {succeeding_g_nt} or A {succeeding_a_nt}\nbut TRBV20-1 or TRBV29-1 is {alt_v_beta_gene_bool}",
                    v_gene_beta,
                    f"\nAs filtering_succeeding_nt_cys104_beta is set to {filter_succeeding_nt_cys_104_beta},\nusing S amino acid that has A as starting nucleotide for 2nd amino acid CDR3-J beta:",
                    idx,
                    "\n",
                )
                second_codon_starting_with_A = return_optimal_codon_starting_with_A_nt(aa="S")  # First S codon is AGC
                prepend = insert_beta["bbsi_1"] + insert_beta["2nt_spacer_1"] + insert_beta["cysteine_overhang" + append_is_alt_v_beta_gene] + second_codon_starting_with_A[1:]

    elif succeeding_a_nt:
        if append_is_alt_v_beta_gene:
            second_codon_starting_with_A = return_optimal_codon_starting_with_A_nt(aa=cdr3_beta[1])
            prepend = insert_beta["bbsi_1"] + insert_beta["2nt_spacer_1"] + insert_beta["cysteine_overhang" + append_is_alt_v_beta_gene] + second_codon_starting_with_A[1:]
        else:
            if filter_succeeding_nt_cys_104_beta:
                print(
                    cdr3_beta,
                    f"CDR3 beta Cys104 succeeding nucleotide is G {succeeding_g_nt} or A {succeeding_a_nt}\n but TRBV20-1 or TRBV29-1 is {alt_v_beta_gene_bool}",
                    v_gene_beta,
                    "\nRemoving TCR:",
                    idx,
                    "\n",
                )
                return filter_succeeding_nt_cys_104_beta, None
            else:
                print(
                    cdr3_beta,
                    f"CDR3 beta Cys104 succeeding nucleotide is G {succeeding_g_nt} or A {succeeding_a_nt}\nbut TRBV20-1 or TRBV29-1 is {alt_v_beta_gene_bool}",
                    v_gene_beta,
                    f"\nAs filtering_succeeding_nt_cys104_beta is set to {filter_succeeding_nt_cys_104_beta},\nusing A amino acid that has G as starting nucleotide as 2nd amino acid for CDR3-J beta:",
                    idx,
                    "\n",
                )
                second_codon_starting_with_G = return_optimal_codon_starting_with_G_nt(aa="A")  # First A codon is GCC
                prepend = insert_beta["bbsi_1"] + insert_beta["2nt_spacer_1"] + insert_beta["cysteine_overhang"] + second_codon_starting_with_G[1:]

    else:
        if filter_succeeding_nt_cys_104_beta:
            print(
                cdr3_beta,
                f"CDR3 beta Cys104 succeeding nucleotide is G {succeeding_g_nt} or A {succeeding_a_nt}!\nand either an G or A succeeding nucleotide is required.\nRemoving TCR:",
                idx,
                "\n",
            )
            return filter_succeeding_nt_cys_104_beta, None
        else:
            if append_is_alt_v_beta_gene:
                print(
                    cdr3_beta,
                    f"CDR3 beta Cys104 succeeding nucleotide is G {succeeding_g_nt} or A {succeeding_a_nt}\nand TRBV20-1 or TRBV29-1 is {alt_v_beta_gene_bool}",
                    v_gene_beta,
                    f"\nAs filtering_succeeding_nt_cys104_beta is set to {filter_succeeding_nt_cys_104_beta},\nusing S amino acid that has A as starting nucleotide for 2nd amino acid CDR3-J beta:",
                    idx,
                    "\n",
                )
                second_codon_starting_with_A = return_optimal_codon_starting_with_A_nt(aa="S")  # First S codon is AGC
                prepend = insert_beta["bbsi_1"] + insert_beta["2nt_spacer_1"] + insert_beta["cysteine_overhang"] + second_codon_starting_with_A[1:]
            else:
                print(
                    cdr3_beta,
                    f"CDR3 beta Cys104 succeeding nucleotide is G {succeeding_g_nt} or A {succeeding_a_nt}\nand TRBV20-1 or TRBV29-1 is {alt_v_beta_gene_bool}",
                    v_gene_beta,
                    f"\nAs filtering_succeeding_nt_cys104_beta is set to {filter_succeeding_nt_cys_104_beta},\nusing A amino acid that has G as starting nucleotide as 2nd amino acid for CDR3-J beta:",
                    idx,
                    "\n",
                )
                second_codon_starting_with_G = return_optimal_codon_starting_with_G_nt(aa="A")  # First A codon is GCC
                prepend = insert_beta["bbsi_1"] + insert_beta["2nt_spacer_1"] + insert_beta["cysteine_overhang"] + second_codon_starting_with_G[1:]

    # last codon_j is GAG, which makes GAGG with the 4th nt, then 2nt spacer, and then BbsI
    append = insert_beta["last_codon_j"] + insert_beta["4th_nt_overhang"] + insert_beta["2nt_spacer_2"] + insert_beta["bbsi_2"]

    full_sequence = prepend + cdr3_beta_j119_nt_opt + append

    # check for restriction sites in full_sequence without the intended BbsI sites
    search_string, restriction_enzyme_sites, restriction_edge_sites = make_search_strings(enzymes=["BbsI"])
    if search(search_string, full_sequence[len(insert_beta["bbsi_1"]) : -len(insert_beta["bbsi_2"])]):
        if filter_cannot_be_codon_optimized:
            print(f"{cdr3_beta} CDR3-J beta aa sequence of TCR: {idx} contains restriction site!")
            print(f"Removing TCR: {idx}\n")
        else:
            print(f"{cdr3_beta} CDR3-J beta aa sequence of TCR: {idx} contains restriction site!")
            print(f"As filter_cannot_be_codon_optimized is not set to True, TCR: {idx} will not be removed!\n")
        return filter_cannot_be_codon_optimized, None

    return False, full_sequence


def add_cdr3j_alpha_beta_order_columns(
    tcr_df_v_gene_filtered: pd.DataFrame,
    filter_cannot_be_codon_optimized: bool = True,
    filter_succeeding_nt_cys_104_beta: bool = True,
    allow_duplicates: bool = False,
    steepness: float = 4,
):
    """Add column containing CDR3-J alpha and beta order nucleotide sequences to tcr_df_v_gene_filtered.

    CDR3-J alpha and beta order nucleotide sequences are generated by: (i) codon optimizing their amino acid sequences
    and (ii) adding the Cys104 + muCαlpha/muCbeta 4-base overhangs and BbsI recognition site nucleotide sequences.


    Parameters
    ----------
    tcr_df_v_gene_filtered : pd.DataFrame
        Filtered TCR DataFrame that only contains TCRs of which we have their V genes in our V gene vector stock. Has
        to be generated by filter_tcr_df_on_v_genes_in_stock.
    filter_cannot_be_codon_optimized : bool, default = True
        Filter TCRs of which the CDR3-J alpha and/or beta cannot be codon optimized.
    filter_succeeding_nt_cys_104_beta : bool, default = False
        Filter TCRs of which the combination of the Cys104 beta succeeding nucleotide and the V beta gene results in an
        incorrect Cys104 succeeding amino acid. If False, force incorrect Cys104 amino acid. For example, if the beta
        Cys104 succeeding nucleotide is A and the V beta gene is not TRBV20-1 or TRBV29-1 that have a Cys104 with a
        succeeding A nucleotide, filter_succeeding_nt_cys_104_beta = False will force an incorrect S amino acid that
        starts with an A nucleotide as the 2nd amino acid of the CDR3-J beta.
    allow_duplicates : bool, default = False
        Allow duplicates in the CDR3-J alpha and beta order nucleotide sequences. If False, the function will
        optimize the CDR3-J alpha and beta order nucleotide sequences until a unique sequence is found. This
        can take a long time if there are many duplicates. If True, the function will not check for duplicates and
        will return the first CDR3-J alpha and beta order nucleotide sequences that are generated. This is much faster
        but may result in duplicates in the final TCR DataFrame.
    steepness : float, default = 4.0
        Steepness of the codon optimization. A higher steepness will result in a more optimal codon usage. A lower
        steepness will result in a less optimal codon usage. The steepness is used in the codon_optimize function.
        The steepness is used to determine the steepness of the sigmoid function that is used to calculate the
        codon usage frequency. A minimal steepness of 1 is required, 4 is recommended.
        Steepness 4 has been tested for 100x the same tcr and converges in about 1 minute to generate 100 unique nt alpha and beta sequences.

    Returns
    -------
    tcr_df_v_gene_filtered : pd.DataFrame
        Input tcr_df_v_gene_filtered DataFrame with added columns containing CDR3-J alpha and beta order nucleotide
        sequences.
    """
    tcr_df_v_gene_filtered = tcr_df_v_gene_filtered.copy()
    tcr_df_v_gene_filtered["cdr3j_alpha_nt_order"] = None
    tcr_df_v_gene_filtered["cdr3j_beta_nt_order"] = None
    for idx, cdr3_alpha, cdr3_beta, TRBV, TRAJ, TRBJ in zip(
        tcr_df_v_gene_filtered.index,
        tcr_df_v_gene_filtered["cdr3_alpha_aa"],
        tcr_df_v_gene_filtered["cdr3_beta_aa"],
        tcr_df_v_gene_filtered["TRBV_IMGT_allele_collapsed"],
        tcr_df_v_gene_filtered["TRAJ_IMGT"],
        tcr_df_v_gene_filtered["TRBJ_IMGT"],
    ):
        # remove allele from V gene if present
        TRBV = TRBV.split("*")[0]
        TRAJ = TRAJ.split("*")[0]
        TRBJ = TRBJ.split("*")[0]

        if allow_duplicates:
            alpha_filter_bool, cdr3_alpha_nt = concatenate_alpha(
                cdr3_alpha, TRAJ, idx, filter_cannot_be_codon_optimized=filter_cannot_be_codon_optimized, steepness=steepness, idx_seed=idx
            )

            beta_filter_bool, cdr3_beta_nt = concatenate_beta(
                cdr3_beta,
                TRBV,
                TRBJ,
                idx,
                filter_cannot_be_codon_optimized=filter_cannot_be_codon_optimized,
                filter_succeeding_nt_cys_104_beta=filter_succeeding_nt_cys_104_beta,
                steepness=steepness,
                idx_seed=idx,
            )
        else:
            idx_seed = idx
            duplicate_a = True
            duplicate_b = True
            while duplicate_a or duplicate_b:
                if duplicate_a:
                    alpha_filter_bool, cdr3_alpha_nt = concatenate_alpha(
                        cdr3_alpha, TRAJ, idx, filter_cannot_be_codon_optimized=filter_cannot_be_codon_optimized, steepness=steepness, idx_seed=idx_seed
                    )
                    if (cdr3_alpha_nt not in tcr_df_v_gene_filtered["cdr3j_alpha_nt_order"].values) and (cdr3_alpha_nt is not None):
                        duplicate_a = False
                        if idx != idx_seed:
                            print(f"Unique CDR3-J alpha nt seq made for TCR: {idx}!\n")
                    else:
                        print(f"Duplicate CDR3-J alpha nt seq made for TCR: {idx}, retrying!")
                if duplicate_b:
                    beta_filter_bool, cdr3_beta_nt = concatenate_beta(
                        cdr3_beta,
                        TRBV,
                        TRBJ,
                        idx,
                        filter_cannot_be_codon_optimized=filter_cannot_be_codon_optimized,
                        filter_succeeding_nt_cys_104_beta=filter_succeeding_nt_cys_104_beta,
                        steepness=steepness,
                        idx_seed=idx_seed,
                    )
                    if (cdr3_beta_nt not in tcr_df_v_gene_filtered["cdr3j_beta_nt_order"].values) and (cdr3_beta_nt is not None):
                        duplicate_b = False
                        if idx != idx_seed:
                            print(f"Unique CDR3-J beta nt seq made for TCR: {idx}!\n")
                    else:
                        print(f"Duplicate CDR3-J beta nt seq made for TCR: {idx}, retrying!")

                idx_seed += 1
                if idx_seed > 10000:
                    if duplicate_a and not duplicate_b:
                        print(f"Duplicate CDR3-J alpha cannot be optimized away for TCR: {idx}! Keeping duplicate CDR3-J alpha!")
                        duplicate_a = False
                    elif duplicate_b and not duplicate_a:
                        raise Exception("Duplicate CDR3-J beta cannot be optimized away for TCR:", idx, "!")

        if (beta_filter_bool and cdr3_beta_nt is not None) or (alpha_filter_bool and cdr3_alpha_nt is not None):
            raise Exception("filter_bool is True but CDR3-J beta is not None!\nThis should not occur!")

        filter_bool = alpha_filter_bool or beta_filter_bool

        if filter_bool:
            tcr_df_v_gene_filtered.drop(idx, axis=0, inplace=True)
        else:
            tcr_df_v_gene_filtered.loc[idx, "cdr3j_alpha_nt_order"] = cdr3_alpha_nt
            tcr_df_v_gene_filtered.loc[idx, "cdr3j_beta_nt_order"] = cdr3_beta_nt

    tcr_df_v_gene_filtered.reset_index(inplace=True, drop=True)
    return tcr_df_v_gene_filtered


def assign_tcrs_to_wells_from_layout(df_plate: pd.DataFrame, tcr_assignment_plate_df_stack: pd.DataFrame) -> pd.DataFrame:
    """
    Assign TCRs to wells following a predefined plate layout template. New layouts can be defined in
    tcr_toolbox_data_path + `/tcr_toolbox_datasets/tcr_assembly/tcr_assignment_plate_layout`.
    """

    if not df_plate.shape[0] == tcr_assignment_plate_df_stack.shape[0]:
        raise Exception("Plate dataframe does not have the same shape as the tcr_assignment_plate_stack dataframe!")

    df_plate = df_plate.reset_index(drop=True)

    df_merged = pd.merge(tcr_assignment_plate_df_stack, df_plate, how="left", left_on="tcr_idx", right_index=True)
    df_merged.loc[:, "well"] = df_merged.loc[:, "well_letter"] + df_merged.loc[:, "well_number"].astype(str)
    df_merged.drop(["well_letter", "well_number"], axis=1, inplace=True)
    df_merged.insert(0, "well", df_merged.pop("well"))
    df_merged = add_well_coordinates(df_merged)
    df_merged.drop(["row", "column"], axis=1, inplace=True)

    return df_merged


def add_standardized_tcr_assembly_well_plate_layout(run_tcr_df: pd.DataFrame, grouping_col: str, well_plate_size: int = 384):
    """
    Assign TCRs to 384-well plates by group. TCRs from different groups never share the same plate.
    Raises an error if the last plate of a group would be less than half-filled.

    Parameters
    ----------
    run_tcr_df : pd.DataFrame
        DataFrame containing TCRs and their group assignments.
    grouping_col : str
        Column name for TCR group assignment (required).
    well_plate_size : int, default=384
        Maximum wells per plate.

    Returns
    -------
    plates_df_dict : defaultdict(pd.DataFrame)
        Dictionary of plates where keys are plate numbers and values are plate DataFrames containing TCRs assigned to wells of that plate.
    """
    if not grouping_col:
        raise ValueError("run_tcr_df needs to have a TCR grouping column!")

    if "well" in run_tcr_df.columns:
        raise Exception("Well column already exists in input DataFrame! Were TCRs already assigned to wells?")

    if well_plate_size == 384:
        tcr_assignment_plate_df_stack_384 = pd.read_excel(
            os.path.join(
                tcr_toolbox_data_path,
                "tcr_toolbox_datasets",
                "tcr_assembly",
                "tcr_assignment_plate_layout",
                "tcr_assignment_plate_layout_stack_384.xlsx",
            )
        )
    else:
        raise NotImplementedError(f"{well_plate_size}-well plates have not been implementend yet!")

    plates_df_dict = collections.defaultdict(pd.DataFrame)
    plate_counter = 1
    plate_assignments = []

    for group in run_tcr_df.loc[:, grouping_col].unique():
        group_df = run_tcr_df.loc[run_tcr_df.loc[:, grouping_col] == group, :].copy()
        n_tcrs = len(group_df)
        full_plates, remainder = divmod(n_tcrs, well_plate_size)

        if remainder and remainder < well_plate_size // 2:
            raise ValueError(
                f"Group '{group}' would have a last plate with {remainder} TCRs (< half of {well_plate_size})!\nWe only assemble at least half-filled last plates per group.\nPlease remove plate by removing its TCRs or add TCRS until plate is at least half-filled."
            )

        group_df = group_df.reset_index(drop=True)
        num_plates = full_plates + (1 if remainder else 0)
        plates_for_group = np.repeat(np.arange(plate_counter, plate_counter + num_plates), well_plate_size)[:n_tcrs]
        group_df.loc[:, "plate_number"] = plates_for_group
        plate_counter += num_plates
        plate_assignments.append(group_df)

    all_df = pd.concat(plate_assignments, ignore_index=True)

    for plate_num, df_plate in all_df.groupby("plate_number", sort=True):
        df_plate = df_plate.reset_index(drop=True).reindex(range(well_plate_size))
        plates_df_dict[plate_num] = assign_tcrs_to_wells_from_layout(df_plate, tcr_assignment_plate_df_stack_384)

    plate_to_group = {plate_num: df.loc[:, grouping_col].iloc[0] for plate_num, df in plates_df_dict.items() if not df.empty}

    for plate_num in sorted(plates_df_dict.keys()):
        print(f"Writing {plate_to_group[plate_num]} to plate: {plate_num}")

    return plates_df_dict


def add_plate_ortho_primer_combinations_to_cdr3j_seqs(plates_df_dict: collections.defaultdict[pd.DataFrame], remove_longer_than_200_nt: bool = True):
    """Add 384-well plate sub-pool-specific Fw + Rev ortho primers and well-specific ortho primers to CDR3-J + 4-base
    overhangs + BbsI recognition sites order nucleotide sequences.

    Plates are ordered together in oligo pools of maximal 9 plates. plate sub-pool-specific Fw + Rev ortho primers are
    used to first amplify plate-specific oligo's from the full ordered oligo pool. Afterwards, well-specific
    ortho primers are used to amplify well-specific CDR3-J alpha and beta oligo's from the amplified
    plate-specific-sub-pool.

    Parameters
    ----------
    plates_df_dict : collections.defaultdict[pd.DataFrame]
        Dictionary where keys are plate numbers and values are DataFrames that store plate data. Has to be generated
        by add_standardized_tcr_assembly_well_plate_layout.
    remove_longer_than_200_nt : bool, default = True
        If True, order oligos that are longer than 200 nucleotides, after ortho primers are added, will be removed.
        This is needed because Twist determines its oligo order cost category based on maximal oligo length
        in the ordered pool.

    Returns
    -------
    plates_df_dict : collections.defaultdict[pd.DataFrame]
        Dictionary where keys are plate numbers and values are DataFrames that store plate data.
        'cdr3j_alpha_nt_order_primers' and 'cdr3j_beta_nt_order_primers' columns have been added to each DataFrame that
        contain CDR3-J order nucleotide sequences + required Fw + Rev ortho primer sequences.

    Notes
    -----
    Standardized plate sub-pool-specific Fw + Rev ortho primer combinations can be found in:
    /tcr_toolbox_datasets/tcr_assembly/final_ortho_primer_combination_source_target_dispense_sheets/plate_pool_primers_df.xlsx

    Standardized well-specific Fw + Rev ortho primer combinations can be found in:
    /tcr_toolbox_datasets/tcr_assembly/final_ortho_primer_combination_source_target_dispense_sheets/final_ortho_primer_combination_assignment_df_384_short_correct_rev.xlsx
    """

    if plates_df_dict[1].shape[0] == 384:
        final_ortho_combination_assignment_df = pd.read_excel(
            os.path.join(
                tcr_toolbox_data_path,
                "tcr_toolbox_datasets",
                "tcr_assembly",
                "final_ortho_primer_combination_source_target_dispense_sheets",
                "final_ortho_primer_combination_assignment_df_384_short_correct_rev.xlsx",
            ),
            index_col=0,
        )
    else:
        raise NotImplementedError(f"{plates_df_dict[1].shape[0]}-well plate orthoprimer combination assignment has not been implemented yet!")

    plate_primers_df = pd.read_excel(
        os.path.join(
            tcr_toolbox_data_path,
            "tcr_toolbox_datasets",
            "tcr_assembly",
            "final_ortho_primer_combination_source_target_dispense_sheets",
            "plate_pool_primers_df.xlsx",
        ),
        index_col=0,
    )

    run_sub_pool = 1
    for i, plate in enumerate(plates_df_dict.keys()):
        print("run_sub_pool:", run_sub_pool)
        print("plate:", plate)
        plates_df_dict[plate] = pd.merge(plates_df_dict[plate], final_ortho_combination_assignment_df, how="left", on="well")

        if (i != 0) and (i % plate_primers_df.shape[0] == 0):
            run_sub_pool += 1

        plate_primer_int = i % plate_primers_df.shape[0]
        plates_df_dict[plate]["run_sub_pool"] = run_sub_pool

        for idx in plates_df_dict[plate].index:
            if not (isinstance(plates_df_dict[plate].loc[idx]["cdr3j_alpha_nt_order"], str) or isinstance(plates_df_dict[plate].loc[idx]["cdr3j_beta_nt_order"], str)):
                if math.isnan(plates_df_dict[plate].loc[idx]["cdr3j_alpha_nt_order"]) or math.isnan(plates_df_dict[plate].loc[idx]["cdr3j_beta_nt_order"]):
                    plates_df_dict[plate].loc[idx, "cdr3j_alpha_nt_order_primers"] = np.nan
                    plates_df_dict[plate].loc[idx, "cdr3j_beta_nt_order_primers"] = np.nan
                else:
                    raise Exception("CDR3-J alpha or beta type can only be NaN if not a str!")

            else:
                if remove_longer_than_200_nt and (
                    (len(plates_df_dict[plate].loc[idx, "cdr3j_alpha_nt_order"]) > (200 - 80)) or (len(plates_df_dict[plate].loc[idx, "cdr3j_beta_nt_order"]) > (200 - 80))
                ):
                    print(
                        "TCR idx:",
                        idx,
                        "+ primers is longer than 200 nt. Removed from plate!",
                        "If you do not want to remove this TCR idx, set remove_longer_than_200_nt to False.",
                    )
                    plates_df_dict[plate].loc[idx, "cdr3j_alpha_nt_order_primers"] = np.nan
                    plates_df_dict[plate].loc[idx, "cdr3j_beta_nt_order_primers"] = np.nan

                else:
                    plate_primer_fw = plate_primers_df.loc[plate_primer_int, "fw"]
                    ortho_primer_fw = plates_df_dict[plate].loc[idx, "ortho_primer_combination"].split("-")[0]
                    cdr3j_alpha = plates_df_dict[plate].loc[idx, "cdr3j_alpha_nt_order"]
                    ortho_primer_rev = str(Dseq(plates_df_dict[plate].loc[idx, "ortho_primer_combination"].split("-")[1]).reverse_complement())
                    plate_primer_rev = str(Dseq(plate_primers_df.loc[plate_primer_int, "rv_correct"]).reverse_complement())
                    cdr3j_beta = plates_df_dict[plate].loc[idx, "cdr3j_beta_nt_order"]

                    plates_df_dict[plate].loc[idx, "cdr3j_alpha_nt_order_primers"] = plate_primer_fw + ortho_primer_fw + cdr3j_alpha + ortho_primer_rev + plate_primer_rev

                    plates_df_dict[plate].loc[idx, "cdr3j_beta_nt_order_primers"] = plate_primer_fw + ortho_primer_fw + cdr3j_beta + ortho_primer_rev + plate_primer_rev

    return plates_df_dict


def make_tcr_assembly_run_cdr3j_seq_order_sheets(plates_df_dict: pd.DataFrame, run_path: Union[str, os.PathLike[str]], numbered_run_name: str, name_annotation_cols: list = None):
    """Write Twist CDR3-J oligo pool order sheet(s). One order .xlsx sheet is written per 9 plates.

    Parameters
    ----------
    plates_df_dict : collections.defaultdict[pd.DataFrame]
        Dictionary where keys are plate numbers and values are DataFrames that store plate data. Has to be generated by
        add_plate_ortho_primer_combinations_to_cdr3j_seqs.
    run_path : Union[str, os.PathLike[str]],
        String path to the standardized assembly run directory.
    numbered_run_name :  str
        Name of the assembly run.
    name_annotation_cols : list = []
        List of annotations columns that need to be used to annotate oligo names. Annotations are added to oligo name
        strings in the order of name_annotation_cols list with a "_" in between.

    Returns
    -------
    plates_df_dict : collections.defaultdict[pd.DataFrame]
        Dictionary where keys are plate numbers and values are DataFrames that store plate data. Empty wells that are
        dropped from each plate DataFrame because these are not included in the oligo pool order sheet(s).

    """
    run_sub_pool_dict = collections.defaultdict(int)
    for plate in plates_df_dict.keys():
        print("Preparing plate:", plate)
        run_sub_pool_dict[plate] = plates_df_dict[plate]["run_sub_pool"].unique()[0]

        plates_df_dict[plate]["plate"] = plate
        plates_df_dict[plate]["name"] = (
            plates_df_dict[plate]["run_sub_pool"].astype(str)
            + "_"
            + plates_df_dict[plate]["plate"].astype(str)
            + "_"
            + plates_df_dict[plate]["well"].astype(str)
            + "_"
            + plates_df_dict[plate]["tcr_idx"].astype(str)
        )
        for annotation_col in name_annotation_cols:
            plates_df_dict[plate]["name"] = plates_df_dict[plate]["name"] + "_" + plates_df_dict[plate][annotation_col].astype(str)

        print(f"Writing plate sheet {plate}")
        plates_df_dict[plate].to_excel(os.path.join(run_path, "plate_sheets", numbered_run_name + "_plate_" + str(plate) + ".xlsx"))

    run_sub_pool_order_df_dict = collections.defaultdict(pd.DataFrame)
    for plate, run_sub_pool in run_sub_pool_dict.items():
        if run_sub_pool not in run_sub_pool_order_df_dict.keys():
            print("Processing run_sub_pool:", run_sub_pool, "...")
            run_sub_pool_order_df_dict[run_sub_pool] = pd.DataFrame({"name": [], "sequence": []})

        if not (plates_df_dict[plate].index.equals(pd.Index(np.arange(384)))):
            raise Exception("Plate df index was not reset!")

        if plates_df_dict[plate]["cdr3j_alpha_nt_order_primers"].isna().any() or plates_df_dict[plate]["cdr3j_beta_nt_order_primers"].isna().any():
            print("Removing NaN CDR3-J seqs from plate:", plate, "!", "If you did not expect these NaNs, plate requires manual inspection!")
            plates_df_dict[plate].dropna(subset=["cdr3j_alpha_nt_order_primers", "cdr3j_beta_nt_order_primers"], inplace=True)
            plates_df_dict[plate].reset_index(drop=True, inplace=True)

        tmp_alpha_plate_df = plates_df_dict[plate].loc[:, ["name", "cdr3j_alpha_nt_order_primers"]].copy()
        tmp_alpha_plate_df["name"] = tmp_alpha_plate_df["name"] + "_alpha"
        tmp_alpha_plate_df.rename({"cdr3j_alpha_nt_order_primers": "sequence"}, axis=1, inplace=True)

        tmp_beta_plate_df = plates_df_dict[plate].loc[:, ["name", "cdr3j_beta_nt_order_primers"]].copy()
        tmp_beta_plate_df["name"] = tmp_beta_plate_df["name"] + "_beta"
        tmp_beta_plate_df.rename({"cdr3j_beta_nt_order_primers": "sequence"}, axis=1, inplace=True)

        run_sub_pool_order_df_dict[run_sub_pool] = pd.concat([run_sub_pool_order_df_dict[run_sub_pool], tmp_alpha_plate_df, tmp_beta_plate_df])

    for idx, sequence in enumerate(run_sub_pool_order_df_dict[run_sub_pool]["sequence"]):
        if len(sequence) > 200:
            print(
                "Sequence idx",
                idx,
                "is longer than 200 nt!",
                "Sequence length:",
                len(sequence),
                "As having oligo sequences longer than 200 nt will substantially increase oligo pool price, you might want to order this TCR separately and remove TCR sequence from oligo pool by setting remove_longer_than_200_nt = True in add_plate_ortho_primer_combinations_to_cdr3j_seqs!",
            )

    for run_sub_pool in run_sub_pool_order_df_dict.keys():
        run_sub_pool_order_df_dict[run_sub_pool].reset_index(inplace=True, drop=True)
        plates_list = run_sub_pool_order_df_dict[run_sub_pool]["name"].str.split("_").str[1].unique().tolist()
        plates_list = [int(plate) for plate in plates_list]

        print("Writing run sub pool", run_sub_pool, ".xlsx sheet!")
        run_sub_pool_order_df_dict[run_sub_pool].to_excel(
            os.path.join(
                run_path,
                "cdr3j_oligo_order_sheets",
                numbered_run_name + "_run_sub_pool_" + str(run_sub_pool) + "_plates_" + str(min(plates_list)) + "-" + str(max(plates_list)) + ".xlsx",
            ),
            index=False,
        )

    return plates_df_dict


def make_tcr_96_well_plate_eblock_idt_order_sheets(alpha_order: pd.Series, beta_order: pd.Series, identifiers: pd.Series, wells: pd.Series, save_path: str, verbose: bool = False):
    """Make order sheets for IDT or Twist."""
    from tcr_toolbox.tcr_assembly.constants import sites_to_check
    from tcr_toolbox.utils.codon_optimizer import make_search_strings

    (search_string, restriction_enzyme_sites, restriction_edge_sites) = make_search_strings(enzymes=["BbsI"], sites_to_check=sites_to_check)

    if len(alpha_order) != len(beta_order) or len(alpha_order) != len(identifiers) or len(alpha_order) != len(wells):
        raise ValueError("Length of alpha_order, beta_order, identifiers and wells is not the same.")

    # loop over sequences and add random sequences to increase above minimal order length of 300
    sequences = []
    for idx, (alpha, beta, identifier, well) in enumerate(zip(alpha_order, beta_order, identifiers, wells)):
        if not pd.isna(alpha) and not pd.isna(beta):
            # loop until we have correct random inserts
            correct_sequence = False
            i = 0
            while not correct_sequence:
                random_before_alpha = generate_random_dna_sequence(length=6)
                random_after_beta = generate_random_dna_sequence(length=6)
                alpha_sequence = random_before_alpha + alpha
                beta_sequence = beta + random_after_beta
                middle_len = 301 - (len(alpha_sequence) + len(beta_sequence))

                # if the middle length is negative, we remove the random sequences from the ends
                if middle_len < 1:
                    alpha_sequence = alpha
                    beta_sequence = beta
                    middle_len = 301 - (len(alpha_sequence) + len(beta_sequence))
                    no_restriction_edge_alpha = True
                    no_restriction_edge_beta = True
                else:
                    no_restriction_edge_alpha = search(restriction_edge_sites, random_before_alpha) is None and search(search_string, random_before_alpha) is None
                    no_restriction_edge_beta = search(restriction_edge_sites, random_after_beta) is None and search(search_string, random_after_beta) is None

                random_between_alpha_beta = generate_random_dna_sequence(length=middle_len)
                no_restriction_middle = search(restriction_edge_sites, random_between_alpha_beta) is None and search(search_string, random_between_alpha_beta) is None
                sequence = alpha_sequence + random_between_alpha_beta + beta_sequence

                if verbose and i % 100 == 1:
                    print(search(restriction_edge_sites, random_before_alpha), search(search_string, alpha))
                    print(random_before_alpha, alpha, random_between_alpha_beta, beta, random_after_beta)
                    print(no_restriction_edge_alpha, no_restriction_edge_beta, no_restriction_middle)
                if no_restriction_edge_alpha and no_restriction_edge_beta and no_restriction_middle:
                    sequences.append(sequence)
                    correct_sequence = True
                    if verbose == 2:
                        print(f"{idx} correct sequence, iterations {i}")
                i += 1
        else:
            sequences.append(None)

    name = identifiers
    well_position = wells

    df = pd.DataFrame({"Well Position": well_position, "Name": name, "Sequence": sequences})
    # save at save_path in Excel format
    if save_path:
        df.to_excel(save_path, index=False)
    else:
        return df


def substract_v_usage_from_v_rack_csv(v_rack_most_recent_df_dict: dict, v_usage_df: pd.DataFrame, run_mode: str = "oligo_order"):
    """
    Parameters
    ----------
    v_rack_most_recent_df_dict : dict
        Dictionary that contains the most recent v-rack with 2d barcode tubes
    v_usage_df : pd.DataFrame
        DataFrame that contains the log of v-gene usage

    Returns
    -------
    v_rack_most_recent_df_dict: dict
        Dictionary that contains Va and Vb adapted v-rack with 2d barcode tubes


    """
    for row in v_usage_df.iterrows():
        v_ab = row[1]["source_ID"]
        va_or_vb = "Va" if v_ab == "FK00701986" else "Vb" if v_ab == "FK00701992" else None
        if va_or_vb is None:
            raise Exception(f"v_alpha_or_v_beta is not TRAV or TRBV! {va_or_vb}{v_ab}")
        vol = row[1]["Vol"]
        v_gene = row[1]["matcher_ID"]
        print("v_alpha_or_v_beta:", va_or_vb)
        print("vol:", vol)
        print("v_gene:", v_gene)

        # check if v_gene is in v_rack_most_recent_df
        if v_gene not in v_rack_most_recent_df_dict[va_or_vb]["TRAV_or_TRBV"].unique():
            print(f"WARNING: V gene {v_gene} is not in v_rack_most_recent_df!, attemting replace _ with /")
            v_gene = v_gene.replace("_", "/")
            if v_gene not in v_rack_most_recent_df_dict[va_or_vb]["TRAV_or_TRBV"].unique():
                raise Exception(f"V gene {v_gene} is not in v_rack_most_recent_df!")

        ul_before = v_rack_most_recent_df_dict[va_or_vb].loc[v_rack_most_recent_df_dict[va_or_vb]["TRAV_or_TRBV"] == v_gene, "uL_in_tube"]
        v_rack_most_recent_df_dict[va_or_vb].loc[v_rack_most_recent_df_dict[va_or_vb]["TRAV_or_TRBV"] == v_gene, "uL_in_tube"] -= vol
        ul_after = v_rack_most_recent_df_dict[va_or_vb].loc[v_rack_most_recent_df_dict[va_or_vb]["TRAV_or_TRBV"] == v_gene, "uL_in_tube"]

        if run_mode == "run_assembly":
            print("ul before:", ul_before.values)
            print("ul after:", ul_after.values)
            print()

    return v_rack_most_recent_df_dict


def substract_v_usage_from_v_rack_trc(v_rack_most_recent_df_dict: dict, v_usage_tcr_path: str, run_mode: str = "oligo_order"):
    with open(v_usage_tcr_path, "r", encoding="latin-1") as f:
        v_usage_lines = f.readlines()

    # get va_or_vb = 'Va' if v_ab == 'FK00701986' else 'Vb' if v_ab == 'FK00701992' else None
    for line in v_usage_lines:
        if "<Scan V-gene rack barcodes>,   Closed with: <OK (1)>," in line:
            if "FK00701986" not in line and "FK00701992" not in line:
                raise Exception("scanned V-gene rack barcode is not FK00701986 or FK00701992!")

    # Different format .trc inputs after "> channel" str split:
    # pipette_step = " 1: TRAV_rack, D1, 5 uL "
    # pipette_step = " 1: TRAV_rack, D3 21µL "
    pattern = r"([A-H]\d+)[, ]+\s*(\d+)\s*(?:uL|µL|ul|μL)\b"
    for line in v_usage_lines:
        # check aspirate to get the volumes pipetted out of the v-gene 2d barcode racks
        if "Channel Aspirate (Single Step) - complete" in line:  # we do not care about "Head Aspirate" (water dispense step)
            if "channel" not in line:  # the pipette steps contain 'channel'
                raise Exception("aspirate line does not contain channel!")
            pipetted_steps = line.split("> channel")[1:]
            for pipette_step in pipetted_steps:
                match = re.search(pattern, pipette_step)
                if not match:
                    raise ValueError(f"Could not parse well and volume from: '{pipette_step}'")

                well = match.group(1)
                vol = int(match.group(2))
                va_or_vb = "Va" if "TRAV" in pipette_step else "Vb" if "TRBV" in pipette_step else None
                if va_or_vb is None:
                    raise Exception("va_or_vb is not Va or Vb!")

                ul_before = v_rack_most_recent_df_dict[va_or_vb].loc[v_rack_most_recent_df_dict[va_or_vb]["V_rack_well"] == well, "uL_in_tube"]
                v_rack_most_recent_df_dict[va_or_vb].loc[v_rack_most_recent_df_dict[va_or_vb]["V_rack_well"] == well, "uL_in_tube"] -= vol
                ul_after = v_rack_most_recent_df_dict[va_or_vb].loc[v_rack_most_recent_df_dict[va_or_vb]["V_rack_well"] == well, "uL_in_tube"]
                if run_mode == "run_assembly":
                    print("ul before:", ul_before.values)
                    print("ul after:", ul_after.values)
                    print()

        # this step pipettes the v-gene into the idot source plate, we don't do anything with this at the moment
        elif "Channel Dispense (Single Step) - complete" in line:
            # "1: idot_V_source_plate, A2, 5 uL"
            pass

    return v_rack_most_recent_df_dict


def replace_tube_in_v_rack(v_rack_most_recent_df_dict: dict, tube_df: pd.DataFrame, va_vb: str, run_mode: str = "oligo_order"):
    """
    Parameters
    ----------
    v_rack_most_recent_df_dict : dict
        Dictionary that contains the most recent v-rack with 2d barcode tubes
    tube_df : pd.DataFrame
        DataFrame that contains the tube that should be replaced
    va_vb : str
        String that indicates if the tube is a Va or Vb tube

    Returns
    -------
    v_rack_most_recent_df_dict: dict
        Dictionary that contains Va and Vb adapted v-rack with 2d barcode tubes
    """
    if not len(tube_df) == 1:
        raise Exception("New tube scan log has more than one tube!")
    if str(tube_df["TRAV_or_TRBV"][0]) not in v_rack_most_recent_df_dict[va_vb]["TRAV_or_TRBV"].unique():
        raise Exception(f"New tube scan log has TRAV_or_TRBV ({str(tube_df['TRAV_or_TRBV'][0])}) that are not in new_full_V_rack_hamilton_log!")
    if not all(tube_df.columns == v_rack_most_recent_df_dict[va_vb].columns):
        raise Exception("New tube scan log does not have the same columns as new_full_V_rack_hamilton_log!")

    barcodes_before = v_rack_most_recent_df_dict[va_vb]["V_tube_barcode"].unique()
    v_rack_loc = v_rack_most_recent_df_dict[va_vb]["TRAV_or_TRBV"] == str(tube_df["TRAV_or_TRBV"][0])
    if run_mode == "run_assembly":
        for col in [  # 'V_rack_name', # 'V_rack_barcode', #'V_rack_well', # 'TRAV_or_TRBV'
            "V_tube_barcode",
            "uL_in_tube",
        ]:
            print("col:", col)
            print("vrack value:", v_rack_most_recent_df_dict[va_vb].loc[v_rack_loc, col].values)
            print("tube value:", tube_df[col][0])
            v_rack_most_recent_df_dict[va_vb].loc[v_rack_loc, col] = tube_df[col][0]
            print("rack value after replacement:", v_rack_most_recent_df_dict[va_vb].loc[v_rack_loc, col].values)
            print()
        barcodes_after = v_rack_most_recent_df_dict[va_vb]["V_tube_barcode"].unique()

        print(f"Difference (removed tube) of tube barcodes before and after: {set(barcodes_before).difference(set(barcodes_after))}")
        print(f"Difference (added tube) of tube barcodes before and after: {set(barcodes_after).difference(set(barcodes_before))}")
        print(f"tmp tube barcode: {tube_df['V_tube_barcode'][0]}")
        print()

    return v_rack_most_recent_df_dict


def track_v_barcode_tubes_and_v_usage(hamilton_log_mode: Literal[".csv", ".trc"] = ".trc", run_mode: str = "oligo_order"):
    """
    Parameters
    hamilton_log_mode : Literal['.csv', '.trc'], default = '.trc'

    Returns
    -------
    v_rack_most_recent_df_dict: pd.DataFrame
        DataFrame that contains the most recent v-rack with 2d barcode tubes.

    Writes
    -------
    v-rack with 2d barcode tubes .csv file.

    """

    # 1. read most recent files of barcode tracing
    idot_V_source_usage_hamilton_log_path = os.path.join(v_gene_barcode_tracing_path, "idot_V_source_usage_hamilton_log")
    new_full_V_rack_hamilton_log_path = os.path.join(v_gene_barcode_tracing_path, "new_full_V_rack_hamilton_log")
    new_V_barcode_tube_scan_log_path = os.path.join(v_gene_barcode_tracing_path, "new_V_barcode_tube_scan_log")
    V_rack_position_scan_log_path = os.path.join(v_gene_barcode_tracing_path, "V_rack_position_scan_log")

    # RACKS: get most recent v_rack
    v_rack_lists = {}
    v_rack_most_recent_df_dict = {}
    rack_file_list = os.listdir(new_full_V_rack_hamilton_log_path)
    # files all have (YYYYMMDD_HHMMSS_*.csv), so first remove _ and then sort by number (old[0] to new[-1])
    for va_vb in ["Va", "Vb"]:
        rack_file_list_va_vb = [file_name for file_name in rack_file_list if va_vb in file_name]
        rack_file_list_va_vb = sorted(rack_file_list_va_vb, key=lambda x: x.split("_")[0] + x.split("_")[1])
        v_rack_lists[va_vb] = rack_file_list_va_vb
        v_rack_most_recent_df_dict[va_vb] = pd.read_csv(os.path.join(new_full_V_rack_hamilton_log_path, v_rack_lists[va_vb][-1]))

    most_recent_va_rack = v_rack_lists["Va"][-1]
    most_recent_vb_rack = v_rack_lists["Vb"][-1]
    most_recent_va_datetime = most_recent_va_rack.split("_")[0]
    most_recent_vb_datetime = most_recent_vb_rack.split("_")[0]
    if run_mode == "run_assembly":
        print("Most recent v-racks:", most_recent_va_rack, most_recent_vb_rack, "\n")

    # # check if the most recent alpha and beta date is the same
    if most_recent_va_datetime != most_recent_vb_datetime:
        raise Exception("Most recent Va and Vb rack scan dates are not the same!")
    else:
        most_recent_rack_datetime = most_recent_va_datetime

    # replace old tubes with newly scanned tubes in v_rack and update pipetted volumes
    tube_scan_list = os.listdir(new_V_barcode_tube_scan_log_path)
    tube_scan_list = [file_name for file_name in tube_scan_list if file_name.endswith(".csv")]
    v_usage_log_list = os.listdir(idot_V_source_usage_hamilton_log_path)
    v_usage_log_date_file_dict = {}
    v_usage_log_list = [file_name for file_name in v_usage_log_list if file_name.endswith(hamilton_log_mode)]
    if hamilton_log_mode == ".trc":
        # get all timestamps from the files and add them to the name
        for file_name in v_usage_log_list:
            # read the first line: "2024-07-02 15:47:58> SYSTEM : Analyze method - ..."
            with open(os.path.join(idot_V_source_usage_hamilton_log_path, file_name), "r", encoding="latin-1") as f:
                first_line = f.readline()
            # Match both: "2025-01-28 09:54:52>" and "2025-08-05 13:58:13.807"
            match = re.match(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}(?:\.\d{3})?)", first_line)
            if not match:
                raise Exception(f"First line of file does not contain a valid date! {file_name}")
            date_time = match.group(1)
            date_time_no_ms = date_time.split(".")[0]
            date, time = date_time_no_ms.split(" ")
            date = date.replace("-", "")
            time = time.replace(":", "")
            new_file_name = f"{date}_{time}_V_usage_{file_name}"
            v_usage_log_date_file_dict[new_file_name] = file_name
        v_usage_log_list = list(v_usage_log_date_file_dict.keys())

    combined_tube_scan_v_usage_list = tube_scan_list + v_usage_log_list
    tube_replaced_or_pipetted = False
    for file_name in combined_tube_scan_v_usage_list:
        parts = file_name.split("_")
        if len(parts) < 2:
            raise Exception(f"Filename does not contain enough parts to extract a date! {file_name}")
        date = parts[0] + parts[1]
        if not date.isdigit():
            raise Exception(f"Current file log does not contain a valid date! {file_name}")

    combined_tube_scan_v_usage_list_sorted = sorted(combined_tube_scan_v_usage_list, key=lambda x: x.split("_")[0] + x.split("_")[1])
    for path in combined_tube_scan_v_usage_list_sorted:
        if "Va-tube" in path:
            va_vb = "Va"
        elif "Vb-tube" in path:
            va_vb = "Vb"
        elif "V_usage" in path:
            va_vb = None
        else:
            raise Exception(f"current file log does not contain Va or Vb! {path}")

        if "tube" in path:
            tmp_tube_time = path.split("_")[0] + path.split("_")[1]
            if tmp_tube_time > most_recent_rack_datetime:
                path_new_tube_scan = os.path.join(new_V_barcode_tube_scan_log_path, path)
                tmp_tube = pd.read_csv(path_new_tube_scan)
                if run_mode == "run_assembly":
                    print(f"######## current tube:{str(tmp_tube['TRAV_or_TRBV'][0])} ({path}) ########")
                v_rack_most_recent_df_dict = replace_tube_in_v_rack(v_rack_most_recent_df_dict=v_rack_most_recent_df_dict, tube_df=tmp_tube, va_vb=va_vb, run_mode=run_mode)
                tube_replaced_or_pipetted = True

        elif "V_usage" in path:
            tmp_v_usage_time = path.split("_")[0] + path.split("_")[1]
            if tmp_v_usage_time > most_recent_rack_datetime:
                if run_mode == "run_assembly":
                    print(f"######## substracting usage from {path} ########")
                path_v_gene_usage = os.path.join(idot_V_source_usage_hamilton_log_path, path)

                if hamilton_log_mode == ".csv":
                    tmp_v_usage_df = pd.read_csv(path_v_gene_usage)
                    v_rack_most_recent_df_dict = substract_v_usage_from_v_rack_csv(
                        v_rack_most_recent_df_dict=v_rack_most_recent_df_dict, v_usage_df=tmp_v_usage_df, run_mode=run_mode
                    )
                elif hamilton_log_mode == ".trc":
                    v_usage_tcr_path = os.path.join(idot_V_source_usage_hamilton_log_path, v_usage_log_date_file_dict[path])
                    v_rack_most_recent_df_dict = substract_v_usage_from_v_rack_trc(
                        v_rack_most_recent_df_dict=v_rack_most_recent_df_dict, v_usage_tcr_path=v_usage_tcr_path, run_mode=run_mode
                    )
                tube_replaced_or_pipetted = True
        else:
            raise Exception(f"current file log does not contain tube or pipetting! {path}")

    # compare barcodes and well location with the most recent scan of the v-rack:
    rack_scan_list_sorted = sorted(os.listdir(V_rack_position_scan_log_path), key=lambda x: x.split("_")[0] + x.split("_")[1])
    rack_scan_list_vb = [file_name for file_name in rack_scan_list_sorted if "Vb" in file_name]
    rack_scan_list_va = [file_name for file_name in rack_scan_list_sorted if "Va" in file_name]
    most_recent_vb_scan_df = pd.read_csv(os.path.join(V_rack_position_scan_log_path, rack_scan_list_vb[-1]))
    most_recent_va_scan_df = pd.read_csv(os.path.join(V_rack_position_scan_log_path, rack_scan_list_va[-1]))

    for va_vb in ["Va", "Vb"]:
        v_rack = v_rack_most_recent_df_dict[va_vb]
        v_scan = most_recent_va_scan_df if va_vb == "Va" else most_recent_vb_scan_df
        v_scan_cutoff = "E06" if va_vb == "Va" else "D07"
        index_cutoff = v_scan[v_scan["V_rack_well"] == v_scan_cutoff].index.tolist()[0]
        v_scan = v_scan.loc[:index_cutoff]

        if not all(v_scan["V_tube_barcode"].values == v_rack["V_tube_barcode"].values):
            # # check if barcodes match
            raise Exception(
                f"v_tube_barcodesscan barcodes {v_scan['V_tube_barcode'].values}rack barcodes {v_rack['V_tube_barcode'].values}Va tube barcodes do not match between most recent scan and v_rack_most_recent_df!"
            )

        for well in v_scan["V_rack_well"].unique():
            well_scan = v_scan[v_scan["V_rack_well"] == well]
            well_rack = v_rack[v_rack["V_rack_well"] == well]

            if not all(well_scan["V_tube_barcode"].values == well_rack["V_tube_barcode"].values):
                raise Exception(
                    f"well: {well}well scan: {well_scan['V_tube_barcode'].values}well rack: {well_rack['V_tube_barcode'].values}V tube barcodes do not match between most recent scan and v_rack_most_recent_df for well {well}!"
                )

    if tube_replaced_or_pipetted:
        va_rack_path_barcode = v_rack_lists["Va"][-1].split("_")[-1]
        vb_rack_path_barcode = v_rack_lists["Vb"][-1].split("_")[-1]

        now = datetime.datetime.now()
        now = now.strftime("%Y%m%d_%H%M%S")
        va_new_rack_path = os.path.join(new_full_V_rack_hamilton_log_path, f"{now}_Va_rack_checkpoint_{va_rack_path_barcode}")
        vb_new_rack_path = os.path.join(new_full_V_rack_hamilton_log_path, f"{now}_Vb_rack_checkpoint_{vb_rack_path_barcode}")

        va_rack = v_rack_most_recent_df_dict["Va"]
        vb_rack = v_rack_most_recent_df_dict["Vb"]
        if run_mode == "run_assembly":
            va_rack.to_csv(va_new_rack_path, index=False)
            vb_rack.to_csv(vb_new_rack_path, index=False)
            print("New v-rack saved:", va_new_rack_path, vb_new_rack_path)

    else:
        if run_mode == "run_assembly":
            print("No new tubes or pipetting done, not writing new v-rack.")

    return v_rack_most_recent_df_dict


def make_idot_v_gene_premixing_dispense_csv(
    plates_df_dict: dict,
    source_name: str,
    run_path: Union[str, os.PathLike[str]],
    transfer_volume_nl: float,
    numbered_run_name: str = "",
    max_idot_vol_nl: float = 70_000,
    v_gene_hamilton_ul: float = 21.11,
    water_hamilton_ul: float = 58.88,
    min_v_gene_stock_tube_vol_ul: float = 50.00,
    run_mode: str = "oligo_order",
    skip_barcoded_v_gene_stock_tube_vol_update: bool = False,
):
    """
    Generate the iDOT V gene premixing dispense instruction files for a TCR assembly run.

    This function integrates all plate-level DataFrames from a TCR assembly run to
    calculate V gene usage, assign wells on an iDOT source plate, and generate
    dispense instructions for both iDOT dispenser and Hamilton pipetting robotics systems.

    It ensures that:
    - Each V gene has sufficient source volume available.
    - The total transferred volumes match Hamilton and iDOT limits.
    - The final dispense instructions reflect all TRAV and TRBV genes used
      across the run.

    The function outputs:
    - An Excel file describing the V gene iDOT source plate layout.
    - A `.csv` file containing iDOT dispense instructions.
    - A Hamilton pipetting worklist (`.xlsx` and `.csv`) for transferring
      V genes from barcode-tracked racks to the iDOT source plate.

    Parameters
    ----------
    plates_df_dict : dict
        Dictionary mapping plate identifiers to DataFrames containing TCR assembly
        plate data. Generated by
        `make_tcr_assembly_run_cdr3j_seq_order_sheets()`.
    source_name : str
        Name of the iDOT 384-well source plate containing V genes.
    run_path : str or os.PathLike
        Path to the standardized assembly run directory.
    transfer_volume_nl : float
        Volume of V gene to dispense (in nanoliters) per transfer.
    numbered_run_name : str, optional
        Name of the assembly run in the format `rX_<run_name>`. Default is empty.
    max_idot_vol_nl : float, default=70_000
        Maximum total volume (in nanoliters) that can be held in a 384-well iDOT
        source plate well.
    v_gene_hamilton_ul : float, default=21.11
        µL volume of V gene transferred by the Hamilton robot from the V-rack
        to the iDOT 96-well source plate.
    water_hamilton_ul : float, default=58.88
        µL volume of nuclease-free water added by the Hamilton to dilute the V genes
        in the iDOT source plate.
    min_v_gene_stock_tube_vol_ul : float, default=50.0
        Minimum allowable remaining V gene stock tube volume (in microliters)
        after all transfers are performed.
    run_mode : {'oligo_order', 'run_assembly', 'simulation'}, default='oligo_order'
        Defines whether the function is running in planning mode (`oligo_order`),
        live execution mode (`run_assembly`), or dry-run (`simulation`).
    skip_barcoded_v_gene_stock_tube_vol_update : bool, default=False
        If True, skips updating barcoded V gene stock tube volumes and bypasses checks
        for sufficient remaining volume.

    Returns
    -------
    run_tcr_df : pandas.DataFrame
        Combined DataFrame containing all TCRs from the assembly run after merging
        all plates in `plates_df_dict` and removing empty wells.

    Raises
    ------
    Exception
        If NaN values are found in critical sequence columns,
        if volume constraints are violated, or
        if V gene usage counts mismatch between intermediate and final outputs.
    FileNotFoundError
        If required reference layout or tracking files are missing.
    NotImplementedError
        If `run_mode` is not one of ('oligo_order', 'run_assembly', 'simulation').

    Notes
    -----
    - The function performs extensive internal validation to ensure data and
      volume consistency between iDOT and Hamilton systems.
    - In `oligo_order` mode, warnings about low-volume V tubes are printed.
      In `run_assembly` mode, these conditions raise exceptions instead.
    - The minimum remaining volume in any V tube after transfer is defined internally
      (currently 50 µL).
    - The Hamilton and iDOT setup assumes a dead volume of 10 µL per source well.

    Outputs
    -------
    The following files are generated under `run_path/v_genes_premix_dispense/`:
        - `<source_name>.xlsx`: iDOT source plate layout.
        - `<numbered_run_name>_idot_dispense.csv`: iDOT dispense instructions.
        - `<timestamp>_<source_name>_hamilton_v-rack_to_idot-source_worklist.xlsx`:
          Hamilton pipetting worklist for V gene transfer.

    """
    minimum_vol = min_v_gene_stock_tube_vol_ul

    run_tcr_df = pd.concat([plates_df_dict[plate] for plate in plates_df_dict.keys()])

    if run_tcr_df["cdr3j_alpha_nt_order_primers"].isna().any() or run_tcr_df["cdr3j_beta_nt_order_primers"].isna().any():
        raise Exception("NaN CDR3-J + ortho primer order sequences should not occur, as these were removed by make_tcr_assembly_run_cdr3j_seq_order_sheets from plates_df_dict!")

    if (v_gene_hamilton_ul + water_hamilton_ul) != ((max_idot_vol_nl / 1000) + 10):
        raise Exception(
            "max_idot_vol_nl does not match what is added by the Hamilton to the idot source wells!",
            "v_gene_hamilton_ul + water_hamilton_ul needs to be equal to (max_idot_vol_nl/1000) + 10 uL dead volume!",
        )

    v_alphas_usage_df, v_betas_usage_df = count_run_v_gene_usage(run_path=run_path, numbered_run_name=numbered_run_name, run_tcr_df=run_tcr_df)

    v_alphas_usage_df["TRAV_or_TRBV"] = v_alphas_usage_df.index
    v_alphas_usage_df = v_alphas_usage_df.reset_index(drop=True)
    v_alphas_usage_df["number_of_idot_wells_needed"] = v_alphas_usage_df["count"].multiply(transfer_volume_nl).div(max_idot_vol_nl).apply(np.ceil)

    v_betas_usage_df["TRAV_or_TRBV"] = v_betas_usage_df.index
    v_betas_usage_df = v_betas_usage_df.reset_index(drop=True)
    v_betas_usage_df["number_of_idot_wells_needed"] = v_betas_usage_df["count"].multiply(transfer_volume_nl).div(max_idot_vol_nl).apply(np.ceil)

    v_both_usage_df = pd.concat([v_alphas_usage_df, v_betas_usage_df])
    v_both_usage_check_df = v_both_usage_df.copy()
    v_both_usage_df.reset_index(drop=True, inplace=True)
    v_both_usage_df = (
        v_both_usage_df.loc[v_both_usage_df.index.repeat(v_both_usage_df["number_of_idot_wells_needed"])]
        .drop(["count", "number_of_idot_wells_needed"], axis=1)
        .reset_index(drop=True)
    )

    if v_both_usage_df.shape[0] > 96:
        raise Exception("Not yet implemented that both V gene run usage does not fit on one idot source plate!", v_both_usage_df.shape[0])

    v_gene_assignment_idot_source_plate_df = pd.read_excel(
        os.path.join(
            tcr_toolbox_data_path,
            "tcr_toolbox_datasets",
            "tcr_assembly",
            "v_gene_assignment_idot_source_plate_layout",
            "v_gene_assignment_idot_source_plate_layout_stack.xlsx",
        )
    )
    v_gene_assignment_idot_source_plate_df["well"] = v_gene_assignment_idot_source_plate_df["well_letter"] + v_gene_assignment_idot_source_plate_df["well_number"].astype(str)
    v_gene_assignment_idot_source_plate_df.drop(["well_letter", "well_number"], axis=1, inplace=True)

    v_gene_idot_source_plate = (
        pd.merge(v_gene_assignment_idot_source_plate_df, v_both_usage_df, how="left", left_index=True, right_index=True).drop("v_gene", axis=1).reset_index(drop=True)
    )
    v_gene_idot_source_plate = v_gene_idot_source_plate[~v_gene_idot_source_plate.loc[:, "TRAV_or_TRBV"].isna()].copy()
    v_gene_idot_source_plate["plate"] = source_name

    v_gene_idot_source_plate.to_excel(os.path.join(run_path, "v_genes_premix_dispense", source_name + ".xlsx"))

    run_tcr_df_v_gene_premix = run_tcr_df.loc[:, ["well", "tcr_idx", "TRAV_IMGT_allele_collapsed", "TRBV_IMGT_allele_collapsed", "plate"]].copy()
    run_tcr_df_v_gene_premix = pd.melt(
        run_tcr_df_v_gene_premix,
        id_vars=["well", "tcr_idx", "plate"],
        value_vars=["TRAV_IMGT_allele_collapsed", "TRBV_IMGT_allele_collapsed"],
        var_name="TRAV_or_TRBV_iterator",
        value_name="TRAV_or_TRBV",
    )

    run_tcr_df_v_gene_premix["TRAV_or_TRBV_iterator"] = run_tcr_df_v_gene_premix["TRAV_or_TRBV_iterator"].map(
        {"TRAV_IMGT_allele_collapsed": "TRAV", "TRBV_IMGT_allele_collapsed": "TRBV"}
    )
    if run_tcr_df_v_gene_premix.shape[0] != run_tcr_df.shape[0] * 2:
        raise Exception("V gene melt operation went wrong")

    dispense_df = write_idot_dispense_csv(
        target_plates_list=run_tcr_df_v_gene_premix["plate"].unique().tolist(),
        source_df_plates_col_name="plate",
        target_df_plates_col_name="plate",
        source_df=v_gene_idot_source_plate,
        target_df=run_tcr_df_v_gene_premix,
        source_match_col_name="TRAV_or_TRBV",
        target_match_col_name="TRAV_or_TRBV",
        source_df_well_col_name="well",
        target_df_well_col_name="well",
        csv_write_path=os.path.join(run_path, "v_genes_premix_dispense", numbered_run_name + "_idot_dispense.csv"),
        target_source_subset_iter_list=run_tcr_df_v_gene_premix["TRAV_or_TRBV_iterator"].unique().tolist(),
        target_source_subset_iter_col_name="TRAV_or_TRBV_iterator",
        # volume_col_name: str = '',
        transfer_volume_nl=transfer_volume_nl,
        max_source_plate_volume_nl=max_idot_vol_nl,
    )

    if skip_barcoded_v_gene_stock_tube_vol_update:
        run_mode = "simulation"

    now = datetime.datetime.now()
    current_date = now.strftime("%Y%m%d_%H%M%S")
    new_full_V_rack_hamilton_log_path = os.path.join(v_gene_barcode_tracing_path, "new_full_V_rack_hamilton_log", "log")

    if run_mode == "run_assembly":
        va_new_rack_path_log = os.path.join(new_full_V_rack_hamilton_log_path, f"{current_date}_v_rack_update_log.txt")
        # track_v_barcode_tubes_and_v_usage is called inside the wrap_and_log_to_file for logging to file
        v_rack_most_recent_df_dict = log_function_to_file(
            func=track_v_barcode_tubes_and_v_usage,
            print_log_file_path=va_new_rack_path_log,
            # function params:
            hamilton_log_mode=".trc",
            run_mode=run_mode,
        )
    elif run_mode in ("oligo_order", "simulation"):
        # track_v_barcode_tubes_and_v_usage is called inside the wrap_and_log_to_file for logging to file
        v_rack_most_recent_df_dict = track_v_barcode_tubes_and_v_usage(hamilton_log_mode=".trc", run_mode=run_mode)
    else:
        raise NotImplementedError("Other run_modes than oligo_order and run_assembly are not implemented yet!")

    v_rack_most_recent_df_va_vb = pd.concat([v_rack_most_recent_df_dict["Va"].copy().reset_index(drop=True), v_rack_most_recent_df_dict["Vb"].copy().reset_index(drop=True)])
    v_rack_most_recent_df_va_vb["TRAV_or_TRBV"] = v_rack_most_recent_df_va_vb["TRAV_or_TRBV"].str.replace("/", "_")
    # idot_V_source_transfer_worklist = v_gene_barcode_tracing_path + 'idot_V_source_transfer_worklist'
    v_rack_to_idot_instructions_path = os.path.join(run_path, "v_genes_premix_dispense", f"{current_date}_{source_name}_hamilton_v-rack_to_idot-source_worklist.xlsx")

    v_rack_to_idot_worklist_df = write_hamilton_pipetting_excel_sheet(
        target_plates_list=[source_name],  # --> idot plate name
        source_df_plates_col_name="V_rack_barcode",  # --> alpha rack of beta rack barcode
        target_df_plates_col_name="plate",  # --> barcode voor hamilton dilutions (run_id)
        source_df=v_rack_most_recent_df_va_vb,  # --> concat alpha en beta rack met identifier voor alpha en beta rack (col_name)
        target_df=v_gene_idot_source_plate,  # --> 96 wells voor idot v-gene plate
        source_match_col_name="TRAV_or_TRBV",  # --> TRAV_or_TRBV (v-gene col)
        target_match_col_name="TRAV_or_TRBV",
        source_df_well_col_name="V_rack_well",
        target_df_well_col_name="well",
        max_source_plate_volume_ul=max(v_rack_most_recent_df_va_vb["uL_in_tube"]),  # --> max vol rack tubes (of all tubes in rack)
        excel_write_path=v_rack_to_idot_instructions_path,  # --> path to write excel sheet
        target_annotation_col_names_list=None,  # --> columns to also take (annotations)
        # volume_col_name: str = '',                    # --> target df can code how much vol (optional, this or
        transfer_volume_ul=v_gene_hamilton_ul,  # --> transfer vol of v-gene (optional, this or volume_col_name)
        # volume_solvent_col_name: str = '',            # --> optional
        volume_solvent_ul=water_hamilton_ul,  # --> optional
        target_plate_format=96,
        also_write_csv=True,
    )
    if skip_barcoded_v_gene_stock_tube_vol_update:
        # check how much volume is left in the v-rack tubes after the transfer is hypothetically done (only for warning)
        v_gene_vols_substract_in_run = v_rack_to_idot_worklist_df[["matcher_ID", "Vol"]].groupby("matcher_ID").sum()

        v_gene_not_enough_volume_msg_list = []

        # Iterate over v_genes and their required transfer volumes
        for v_gene, vol in zip(v_gene_vols_substract_in_run.index.to_list(), v_gene_vols_substract_in_run["Vol"].to_list()):
            # Current volume in the rack before transfer
            vol_before = v_rack_most_recent_df_va_vb.loc[v_rack_most_recent_df_va_vb["TRAV_or_TRBV"] == v_gene, "uL_in_tube"].iloc[0]

            # Compute the remaining volume after this transfer
            vol_after = vol_before - vol

            # Check if this would drop below the minimum allowed
            if vol_after < minimum_vol:
                total_needed_vol = v_gene_idot_source_plate.loc[v_gene_idot_source_plate["TRAV_or_TRBV"] == v_gene, :].shape[0] * v_gene_hamilton_ul
                tmp_msg = (
                    f"{v_gene} has less than {minimum_vol} µL left in the tube after transfer!\n"
                    f"Total needed {v_gene} volume is: {total_needed_vol} µL.\n"
                    f"Current {v_gene} batch volume is: {vol_before} µL.\n"
                    f"During CDR3-J oligo pool manufacturing, the current batch of {v_gene} needs "
                    "to be replaced\nwith a new batch before this assembly "
                    f"run can be performed!\n"
                )
                v_gene_not_enough_volume_msg_list.append(tmp_msg)

        # Handle exceptions as before
        if v_gene_not_enough_volume_msg_list:
            if run_mode in ("oligo_order", "simulation"):
                formatted_msgs = [f"WARNING: {msg}" for msg in v_gene_not_enough_volume_msg_list]
                print("\n".join(formatted_msgs))
            elif run_mode == "run_assembly":
                formatted_msgs = [f"ERROR: {msg}" for msg in v_gene_not_enough_volume_msg_list]
                reminder_msg = "Reminder: Run run_mode = oligo_order first to identify tubes that require replacement before performing the wet-lab run_assembly mode!"
                formatted_msgs.append(f"ERROR: {reminder_msg}")
                raise Exception("\n".join(formatted_msgs))

    v_both_usage_check_df.set_index("TRAV_or_TRBV", inplace=True)
    dispense_df_value_counts = dispense_df.loc[:, "Liquid Name"].value_counts()[dispense_df.loc[:, "Liquid Name"].value_counts().index.str.startswith(("TRAV", "TRBV"))]

    if not set(v_both_usage_check_df.index) == set(dispense_df_value_counts.index):
        raise Exception("Unique set of V genes do not overlap between run both usage sheet and final idot dispense csv!")

    if not all(v_both_usage_check_df.loc[dispense_df_value_counts.index, "count"] == dispense_df_value_counts):
        raise Exception("Usage counts do not overlap between run both usage sheet and final idot dispense csv!")

    return run_tcr_df, dispense_df, dispense_df_value_counts, v_both_usage_check_df


def make_echo_v_gene_premixing_dispense_csv(
    plates_df_dict: dict,
    source_name: str,
    run_path: Union[str, os.PathLike[str]],
    transfer_volume_nl: float,
    numbered_run_name: str = "",
    max_echo_vol_nl: int = 50,
    run_mode: str = "oligo_order",
    skip_barcoded_v_gene_tube_vol_update: bool = False,
):
    """Write echo V gene premixing instruction csv file.

    This function likely will not be used anymore because we are not acquiring an Echo dispenser anymore.
    However, the function might be useful for collaborators that have an Echo dispenser.

    Parameters
    ----------
    plates_df_dict : dict
        Dictionary where keys are plate numbers and values are DataFrames that store plate data. Has to be generated by
        make_tcr_assembly_run_cdr3j_seq_order_sheets.
    source_name : str
        Name of the V gene 384-well echo source plate from which is dispensed.
    run_path: Union[str, os.PathLike[str]],
        String path to the standardized assembly run directory that was initialized with
        init_tcr_assembly_run_dir.
    transfer_volume_nl : float = 0.00
        V gene volume to dispense in nL.
    numbered_run_name : str
        Name of the assembly run.
    max_echo_vol_nl : int = 50
        Maximum 384-well echo plate source well volume that can be dispensed from.

    Returns
    -------
    run_tcr_df : pd.DataFrame
        DataFrame that contains all TCRs of the assembly run that was generated by concatenating all plate DataFrames in
        the plates_df_dict and by removing empty wells. This DataFrame is only used for manual inspect and not
        used by any functions after this function.

    Examples
    --------
    >>> source_name = 'cd4_source_1' # you should be able to write this label on the echo source plate
    >>> run_tcr_df = make_echo_v_gene_premixing_dispense_csv(plates_df_dict = plates_df_dict,
    ...                                                      source_name = source_name,
    ...                                                      run_path = run_path,
    ...                                                      transfer_volume_nl = 1.00 * 1000, # Echo needs nanoliter instructions
    ...                                                      numbered_run_name = numbered_run_name,
    ...                                                      max_echo_vol_nl = 50 * 1000 # Echo needs nanoliter instructions
    ...                                                     )
    """
    print(
        "WARNING: make_idot_v_gene_premixing_dispense_csv() only generates the wet-lab tested Echo dispense .CSV file — it does NOT produce the Hamilton pipetting .XLSX worklist for preparing the Echo 384-well source plate."
    )
    run_tcr_df = pd.concat([plates_df_dict[plate] for plate in plates_df_dict.keys()])

    # It took approximately 20 min to pre-mix a full 384-well plate on the Echo
    print("Estimation of hours needed to pre-mix run V genes on Echo:", round((run_tcr_df.shape[0] / 384) * 20 / 60, 2))
    if run_tcr_df.shape[0] > len(plates_df_dict.keys()) * 384:
        raise Exception("This should not occur!")

    if run_tcr_df["cdr3j_alpha_nt_order_primers"].isna().any() or run_tcr_df["cdr3j_beta_nt_order_primers"].isna().any():
        raise Exception("NaN CDR3-J + ortho primer order sequences should not occur, as these were removed by make_tcr_assembly_run_cdr3j_seq_order_sheets from plates_df_dict!")

    v_alphas_usage_df, v_betas_usage_df = count_run_v_gene_usage(run_path=run_path, numbered_run_name=numbered_run_name, run_tcr_df=run_tcr_df)

    v_alphas_usage_df.reset_index(inplace=True)
    v_alphas_usage_df.rename({"index": "TRAV_or_TRBV"}, axis=1, inplace=True)
    v_alphas_usage_df["number_of_echo_wells_needed"] = v_alphas_usage_df["count"].div((max_echo_vol_nl / 1000)).apply(np.ceil)

    v_betas_usage_df.reset_index(inplace=True)
    v_betas_usage_df.rename({"index": "TRAV_or_TRBV"}, axis=1, inplace=True)
    v_betas_usage_df["number_of_echo_wells_needed"] = v_betas_usage_df["count"].div((max_echo_vol_nl / 1000)).apply(np.ceil)

    v_both_usage_df = pd.concat([v_alphas_usage_df, v_betas_usage_df])
    v_both_usage_check_df = v_both_usage_df.copy()
    v_both_usage_df.reset_index(drop=True, inplace=True)
    v_both_usage_df = (
        v_both_usage_df.loc[v_both_usage_df.index.repeat(v_both_usage_df["number_of_echo_wells_needed"])]
        .drop(["count", "number_of_echo_wells_needed"], axis=1)
        .reset_index(drop=True)
    )

    if v_both_usage_df.shape[0] > 384:
        raise Exception("Not yet implemented that both V gene run usage does not fit on one echo source plate!", v_both_usage_df.shape[0])

    v_gene_assignment_echo_source_plate_df = pd.read_excel(
        os.path.join(
            tcr_toolbox_data_path,
            "tcr_toolbox_datasets",
            "tcr_assembly",
            "v_gene_assignment_echo_source_plate_layout",
            "v_gene_assignment_echo_source_plate_layout_stack.xlsx",
        )
    )
    v_gene_assignment_echo_source_plate_df["well"] = v_gene_assignment_echo_source_plate_df["well_letter"] + v_gene_assignment_echo_source_plate_df["well_number"].astype(str)
    v_gene_assignment_echo_source_plate_df.drop(["well_letter", "well_number"], axis=1, inplace=True)

    v_gene_echo_source_plate = (
        pd.merge(v_gene_assignment_echo_source_plate_df, v_both_usage_df, how="left", left_index=True, right_index=True).drop("v_gene", axis=1).reset_index(drop=True)
    )

    v_gene_echo_source_plate.to_excel(os.path.join(run_path, "v_genes_premix_dispense", source_name + ".xlsx"))

    v_gene_echo_source_plate["plate"] = source_name

    run_tcr_df_v_gene_premix = run_tcr_df.loc[:, ["well", "tcr_idx", "TRAV_IMGT_allele_collapsed", "TRBV_IMGT_allele_collapsed", "plate"]].copy()
    run_tcr_df_v_gene_premix = pd.melt(
        run_tcr_df_v_gene_premix,
        id_vars=["well", "tcr_idx", "plate"],
        value_vars=["TRAV_IMGT_allele_collapsed", "TRBV_IMGT_allele_collapsed"],
        var_name="TRAV_or_TRBV_iterator",
        value_name="TRAV_or_TRBV",
    )

    run_tcr_df_v_gene_premix["TRAV_or_TRBV_iterator"] = run_tcr_df_v_gene_premix["TRAV_or_TRBV_iterator"].map(
        {"TRAV_IMGT_allele_collapsed": "TRAV", "TRBV_IMGT_allele_collapsed": "TRBV"}
    )
    if run_tcr_df_v_gene_premix.shape[0] != run_tcr_df.shape[0] * 2:
        raise Exception("V gene melt operation went wrong")

    dispense_df = write_echo_dispense_csv(
        target_plates_list=run_tcr_df_v_gene_premix["plate"].unique().tolist(),
        source_df_plates_col_name="plate",
        target_df_plates_col_name="plate",
        source_df=v_gene_echo_source_plate,
        target_df=run_tcr_df_v_gene_premix,
        source_match_col_name="TRAV_or_TRBV",
        target_match_col_name="TRAV_or_TRBV",
        source_df_well_col_name="well",
        target_df_well_col_name="well",
        csv_write_path=os.path.join(run_path, "v_genes_premix_dispense", numbered_run_name + "_echo_dispense.csv"),
        target_source_subset_iter_list=run_tcr_df_v_gene_premix["TRAV_or_TRBV_iterator"].unique().tolist(),
        target_source_subset_iter_col_name="TRAV_or_TRBV_iterator",
        # volume_col_name: str = '',
        transfer_volume_nl=transfer_volume_nl,
        max_source_plate_volume_nl=max_echo_vol_nl,
    )

    v_both_usage_check_df.set_index("TRAV_or_TRBV", inplace=True)

    if not set(v_both_usage_check_df.index) == set(dispense_df["Sample ID"].value_counts().index):
        raise Exception("Unique set of V genes do not overlap between run both usage sheet and final echo dispense csv!")

    if not all(v_both_usage_check_df.loc[dispense_df["Sample ID"].value_counts().index, :]["count"] == dispense_df["Sample ID"].value_counts()):
        raise Exception("Usage counts do not overlap between run both usage sheet and final echo dispense csv!")

    return run_tcr_df


def pydna_plate_sub_pool_amp_pcr(run_sub_pool_file: Union[str, os.PathLike[str]], run_sub_pool: int):
    """
    Simulate plate sub-pool amplification PCR for CDR3-J alpha and beta oligos in a 384-well plate.

    This function uses plate-specific forward (Fw) and reverse (Rev) orthogonal primers
    from an oligo pool order `.xlsx` sheet and simulates PCR amplification using `pydna`.
    One PCR product per alpha and beta oligo per well is expected; otherwise, an Exception is raised,
    indicating a problem with prerequisite functions.

    Prerequisite functions
    ----------------------
    0. init_tcr_assembly_run_dir
    1. filter_tcr_df_on_v_genes_in_stock
    2. add_cdr3j_alpha_beta_order_columns
    3. concatenate_separate_tcr_df_csvs_into_run_tcr_df
       (Only run if separate CSVs store TCRs for a single assembly run.)
    4. add_standardized_tcr_assembly_well_plate_layout
    5. add_plate_ortho_primer_combinations_to_cdr3j_seqs
    6. make_tcr_assembly_run_cdr3j_seq_order_sheets
    7. make_idot_v_gene_premixing_dispense_csv

    Run after this function
    -----------------------
    pydna_cdr3j_ortho_primer_amp_pcr

    Parameters
    ----------
    run_sub_pool_file : Union[str, os.PathLike[str]]
        Path to the oligo pool order `.xlsx` sheet.
    run_sub_pool : int
        Sub-pool number of the oligo order.
        - If < 10 plates in the run, use 1.
        - If ≥ 10 plates, the first 9 plates correspond to `run_sub_pool = 1`,
          the next 9 plates to `run_sub_pool = 2`.

    Returns
    -------
    plate_alpha_pcr_products_dict : collections.defaultdict[lambda: collections.defaultdict[list]]
        Nested dictionary representing CDR3-J alpha PCR products organized by plates and wells.
        The first level of keys represents plates, while the second level represents wells within each plate.
        Each well key maps to a list of :class:`pydna.amplicon.Amplicon` objects representing the PCR product(s).
    plate_beta_pcr_products_dict : collections.defaultdict[lambda: collections.defaultdict[list]]
        Nested dictionary representing CDR3-J beta PCR products organized by plates and wells.
        The first level of keys represents plates, while the second level represents wells within each plate.
        Each well key maps to a list of :class:`pydna.amplicon.Amplicon` objects representing the PCR product(s).

    Notes
    -----
    Pseudo-example of the plate_alpha_pcr_products_dict or plate_beta_pcr_products_dict dictionaries:
    {
        '1': {
            'A1': [Amplicon(180), Amplicon(160)], # if more than one product is amplified, function will raise an Exception.
            'A2': [Amplicon(186)],
            ...
        },
        '2': {
            'A1': [], # if less than one product is amplified, function will raise an Exception.
            'A2': [Amplicon(174)],
            ...
        },
        ...
    }
    """
    plate_alpha_pcr_products_dict = collections.defaultdict(lambda: collections.defaultdict(list))
    plate_beta_pcr_products_dict = collections.defaultdict(lambda: collections.defaultdict(list))

    plate_primers_df = pd.read_excel(
        os.path.join(
            tcr_toolbox_data_path,
            "tcr_toolbox_datasets",
            "tcr_assembly",
            "final_ortho_primer_combination_source_target_dispense_sheets",
            "plate_pool_primers_df.xlsx",
        ),
        index_col=0,
    )

    run_sub_pool_df = pd.read_excel(run_sub_pool_file)
    run_sub_pool_df["plate"] = run_sub_pool_df["name"].str.split("_").str[1].astype(int)
    run_sub_pool_df["well"] = run_sub_pool_df["name"].str.split("_").str[2]
    run_sub_pool_df["alpha_or_beta"] = run_sub_pool_df["name"].str.split("_").str[-1]

    seen_plate_fw_primers_list = []
    seen_plate_rev_primers_list = []

    for plate in run_sub_pool_df["plate"].unique():
        print("Running plate PCR:", plate)
        tmp_plate_df = run_sub_pool_df[run_sub_pool_df["plate"] == plate].copy()

        # All ortho primers have length 20 nt:
        if (tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "alpha", "sequence"].str[0:20].unique().shape[0] > 1) or (
            tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "alpha", "sequence"].str[-20:].unique().shape[0] > 1
        ):
            raise Exception("CDR3-J alpha nt sequences contain more than one unique plate sub-pool fw or rev primer!")

        if (tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "beta", "sequence"].str[0:20].unique().shape[0] > 1) or (
            tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "beta", "sequence"].str[-20:].unique().shape[0] > 1
        ):
            raise Exception("CDR3-J beta nt sequences contain more than one unique plate sub-pool fw or rev primer!")

        if (
            tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "alpha", "sequence"].str[0:20].tolist()
            != tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "beta", "sequence"].str[0:20].tolist()
        ):
            raise Exception("Fw CDR3-J alpha and beta primers do not match!")

        if (
            tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "alpha", "sequence"].str[-20:].tolist()
            != tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "beta", "sequence"].str[-20:].tolist()
        ):
            raise Exception("Rev CDR3-J alpha and beta primers do not match!")

        seen_plate_fw_primers_list.append(tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "alpha", "sequence"].str[0:20].values[0])
        seen_plate_fw_primers_list.append(tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "alpha", "sequence"].str[-20:].values[0])

        if run_sub_pool == 1:
            plate_idx = plate - 1
        else:
            plate_idx = plate - (1 + (plate_primers_df.shape[0] * (run_sub_pool - 1)))

        primers = {}
        primers["plate_primer_fw"] = Primer(Seq(plate_primers_df.loc[plate_idx, "fw"]))
        primers["plate_primer_rv"] = Primer(Seq(plate_primers_df.loc[plate_idx, "rv_correct"]))

        if tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "alpha", "well"].tolist() != tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "beta", "well"].tolist():
            raise Exception("TCR well annotation of CDR3-J alpha and beta oligos does not match!")

        for tcr_well, seq in zip(tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "alpha", "well"], tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "alpha", "sequence"]):
            plate_alpha_pcr_products_dict[plate][tcr_well] = Anneal((primers["plate_primer_fw"], primers["plate_primer_rv"]), Dseqrecord(Seq(seq))).products
            if len(plate_alpha_pcr_products_dict[plate][tcr_well]) != 1:
                raise Exception("tcr well:", tcr_well, "has more or less than 1 CDR3-J alpha plate sub-pool PCR amplicon!", plate_alpha_pcr_products_dict[plate][tcr_well])

        for tcr_well, seq in zip(tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "beta", "well"], tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "beta", "sequence"]):
            plate_beta_pcr_products_dict[plate][tcr_well] = Anneal((primers["plate_primer_fw"], primers["plate_primer_rv"]), Dseqrecord(Seq(seq))).products
            if len(plate_beta_pcr_products_dict[plate][tcr_well]) != 1:
                raise Exception("tcr well:", tcr_well, "has more or less than 1 CDR3-J alpha plate sub-pool PCR amplicon!", plate_beta_pcr_products_dict[plate][tcr_well])

    if len(seen_plate_fw_primers_list) != len(set(seen_plate_fw_primers_list)):
        raise Exception("There are duplicate plate sub-pool fw primers!", seen_plate_fw_primers_list)

    if len(seen_plate_rev_primers_list) != len(set(seen_plate_rev_primers_list)):
        raise Exception("There are duplicate plate sub-pool fw primers!", seen_plate_fw_primers_list)

    return plate_alpha_pcr_products_dict, plate_beta_pcr_products_dict


def pydna_cdr3j_ortho_primer_amp_pcr(
    run_sub_pool_file: Union[str, os.PathLike[str]],
    plate_alpha_pcr_products_dict: collections.defaultdict[lambda: collections.defaultdict[list]],
    plate_beta_pcr_products_dict: collections.defaultdict[lambda: collections.defaultdict[list]],
):
    """
    Simulate well-specific Fw + Rev orthogonal primer PCR on all amplified plate sub-pool
    CDR3-J alpha and beta oligos.

    If multiple plate sub-pool CDR3-J alpha and/or beta oligos are amplified by well-specific
    Fw + Rev orthogonal primers, all resulting amplicons are recorded in the
    `plate_alpha_pcr_products_dict` and/or `plate_beta_pcr_products_dict`. This allows
    the `pydna_golden_gate_tcr_assembly` simulation to test whether incorrect Golden Gate
    products could be assembled.

    Parameters
    ----------
    run_sub_pool_file : Union[str, os.PathLike[str]]
        Path to the oligo pool order `.xlsx` sheet.
    plate_alpha_pcr_products_dict : collections.defaultdict[lambda: collections.defaultdict[list]]
        Nested dictionary representing CDR3-J alpha PCR products organized by plates and wells.
        The first level of keys represents plates, while the second level represents wells within each plate.
        Each well key maps to a list of :class:`pydna.amplicon.Amplicon` objects representing the PCR product(s).
        Must be generated by `pydna_plate_sub_pool_amp_pcr`.
    plate_beta_pcr_products_dict : collections.defaultdict[lambda: collections.defaultdict[list]]
        Nested dictionary representing CDR3-J beta PCR products organized by plates and wells.
        The first level of keys represents plates, while the second level represents wells within each plate.
        Each well key maps to a list of :class:`pydna.amplicon.Amplicon` objects representing the PCR product(s).
        Must be generated by `pydna_plate_sub_pool_amp_pcr`.

    Returns
    -------
    cdr3_alpha_product_dict : collections.defaultdict[lambda: collections.defaultdict[lambda: collections.defaultdict[list]]]
        Nested dictionary representing CDR3-J alpha PCR products organized by plates and wells.
        - First-level keys: plates
        - Second-level keys: wells of the well-specific Fw + Rev primer combination used for PCR
        - Third-level keys: wells whose CDR3-J alpha oligos are amplified by the second-level well-specific primers
        Each third-level key maps to a list of :class:`pydna.amplicon.Amplicon` objects.
        If a single third-level key equals the second-level key and contains one product, amplification was successful.
    cdr3_beta_product_dict : collections.defaultdict[lambda: collections.defaultdict[lambda: collections.defaultdict[list]]]
        Nested dictionary representing CDR3-J beta PCR products organized by plates and wells.
        - First-level keys: plates
        - Second-level keys: wells of the well-specific Fw + Rev primer combination used for PCR
        - Third-level keys: wells whose CDR3-J beta oligos are amplified by the second-level well-specific primers
        Each third-level key maps to a list of :class:`pydna.amplicon.Amplicon` objects.
        If a single third-level key equals the second-level key and contains one product, amplification was successful.

    Examples
    --------
    >>> cdr3_alpha_product_dict, cdr3_beta_product_dict = pydna_cdr3j_ortho_primer_amp_pcr(
    ...     run_sub_pool_file=run_sub_pool_file,
    ...     plate_alpha_pcr_products_dict=plate_alpha_pcr_products_dict,
    ...     plate_beta_pcr_products_dict=plate_beta_pcr_products_dict
    ... )

    Notes
    -----
    Pseudo-example of `cdr3_alpha_product_dict` or `cdr3_beta_product_dict`:
    {
        '1': {
              'A1': {
                     'A1': [Amplicon(140)] # A1 Fw + Rev well-specific ortho primers correctly amplify a single product from the plate 1 sub-pool CDR3-J alpha/beta oligos of well A1.
              },
              'A2': {
                     'A2': [Amplicon(146)],
                     'P2': [Amplicon(137)] # A2 Fw + Rev well-specific ortho primers incorrectly amplify a product from the plate 1 sub-pool CDR3-J alpha/beta oligos of well P2.
              },
              ...
        },
        '2': {
              'A1': {
                     'A1': [Amplicon(137)]
              },
              'A2': {
                     'A2': [Amplicon(134), Amplicon(104)], # A2 Fw + Rev well-specific ortho primers incorrectly amplify two products from the plate 2 sub-pool CDR3-J alpha/beta oligos of well A2.
              },
              ...
        },
        ...
    }
    """

    if not list(plate_alpha_pcr_products_dict.keys()) == list(plate_beta_pcr_products_dict.keys()):
        raise Exception("Plate keys in alpha and beta plate product dict are out of sync!!")

    final_ortho_combination_assignment_df = pd.read_excel(
        os.path.join(
            tcr_toolbox_data_path,
            "tcr_toolbox_datasets",
            "tcr_assembly",
            "final_ortho_primer_combination_source_target_dispense_sheets",
            "final_ortho_primer_combination_assignment_df_384_short_correct_rev.xlsx",
        ),
        index_col=0,
    )
    run_sub_pool_df = pd.read_excel(run_sub_pool_file)
    run_sub_pool_df["plate"] = run_sub_pool_df["name"].str.split("_").str[1].astype(int)
    run_sub_pool_df["well"] = run_sub_pool_df["name"].str.split("_").str[2]
    run_sub_pool_df["alpha_or_beta"] = run_sub_pool_df["name"].str.split("_").str[-1]

    cdr3_alpha_product_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))
    cdr3_beta_product_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))

    for plate in run_sub_pool_df["plate"].unique():
        print("Running CDR3-J ortho primer PCR plate:", plate)
        tmp_plate_df = run_sub_pool_df[run_sub_pool_df["plate"] == plate].copy()

        if tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "alpha", "well"].tolist() != tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "beta", "well"].tolist():
            raise Exception("TCR well annotation of CDR3-J alpha and beta oligos does not match!")

        for ortho_well in tmp_plate_df.loc[tmp_plate_df["alpha_or_beta"] == "alpha", "well"]:
            primers = {}
            primers["tcr_fw"] = Primer(
                Seq(
                    final_ortho_combination_assignment_df.loc[final_ortho_combination_assignment_df["well"] == ortho_well, "ortho_primer_combination"]
                    .str.split("-")
                    .str[0]
                    .values[0]
                )
            )
            primers["tcr_rv"] = Primer(
                Seq(
                    final_ortho_combination_assignment_df.loc[final_ortho_combination_assignment_df["well"] == ortho_well, "ortho_primer_combination"]
                    .str.split("-")
                    .str[1]
                    .values[0]
                )
            )

            for tcr_well in plate_alpha_pcr_products_dict[plate].keys():
                # alpha:
                if len(plate_alpha_pcr_products_dict[plate][tcr_well]) > 1:
                    # This code block was useful to allow simulation of PCR on multiple plate sub-pool amplicons
                    # However, it will now never be called anymore because we do not allow anymore multiple
                    # plate sub-pool CDR3-J alpha amplicons in the pydna_plate_sub_pool_amp_pcr simulation function.
                    # We should, however, keep this code block to keep the function independent of choices in the
                    # pydna_plate_sub_pool_amp_pcr function.
                    tmp_alpha_list = [Anneal((primers["tcr_fw"], primers["tcr_rv"]), sub_amplicon).products for sub_amplicon in plate_alpha_pcr_products_dict[plate][tcr_well]]
                    tmp_alpha_list = [sub_amplicon for nested_amplicon_list in tmp_alpha_list for sub_amplicon in nested_amplicon_list]
                    if len(tmp_alpha_list) > 0:
                        cdr3_alpha_product_dict[plate][ortho_well][tcr_well] = tmp_alpha_list
                else:
                    tmp_product = Anneal((primers["tcr_fw"], primers["tcr_rv"]), plate_alpha_pcr_products_dict[plate][tcr_well][0]).products
                    if tmp_product:
                        cdr3_alpha_product_dict[plate][ortho_well][tcr_well] = tmp_product

                # beta:
                if len(plate_beta_pcr_products_dict[plate][tcr_well]) > 1:
                    # This code block was useful to allow simulation of PCR on multiple plate sub-pool amplicons
                    # However, it will now never be called anymore because we do not allow anymore multiple
                    # plate sub-pool CDR3-J alpha amplicons in the pydna_plate_sub_pool_amp_pcr simulation function.
                    # We should, however, keep this code block to keep the function independent of choices in the
                    # pydna_plate_sub_pool_amp_pcr function.
                    tmp_beta_list = [Anneal((primers["tcr_fw"], primers["tcr_rv"]), sub_amplicon).products for sub_amplicon in plate_beta_pcr_products_dict[plate][tcr_well]]
                    tmp_beta_list = [sub_amplicon for nested_amplicon_list in tmp_beta_list for sub_amplicon in nested_amplicon_list]
                    if len(tmp_beta_list) > 0:
                        cdr3_beta_product_dict[plate][ortho_well][tcr_well] = tmp_beta_list
                else:
                    tmp_product = Anneal((primers["tcr_fw"], primers["tcr_rv"]), plate_beta_pcr_products_dict[plate][tcr_well][0]).products
                    if tmp_product:
                        cdr3_beta_product_dict[plate][ortho_well][tcr_well] = tmp_product

    return cdr3_alpha_product_dict, cdr3_beta_product_dict


def pydna_golden_gate_tcr_assembly(
    v_genes_premix_dispense_run_file: Union[str, os.PathLike[str]],
    v_genes_source_run_files_list: list[Union[str, os.PathLike[str]]],
    cdr3_alpha_product_dict: collections.defaultdict[lambda: collections.defaultdict[lambda: collections.defaultdict[list]]],
    cdr3_beta_product_dict: collections.defaultdict[lambda: collections.defaultdict[lambda: collections.defaultdict[list]]],
    run_sub_pool_file: Union[str, os.PathLike[str]],
    run_path: Union[str, os.PathLike[str]],
    numbered_run_name: str,
    run_sub_pool: int,
    run_mode: str = "oligo_order",
):
    """
    Simulate TCR Golden Gate assembly for a TCR assembly run and generate a Golden Gate
    TCR assembly log `.txt` file.

    The log file will be empty if all TCRs are assembled correctly. Otherwise, any
    assembly errors are logged with the corresponding plate and well identifiers. The
    function also returns a `ligation_products_dict` that can be used to compare the
    amino acid translation of assembled TCRs with the expected TCR sequences using
    `check_translation_assembled_tcr_ligation_products`.

    Parameters
    ----------
    v_genes_premix_dispense_run_file : Union[str, os.PathLike[str]]
        Path to the V gene premixing dispense instruction .csv file.
    v_genes_source_run_files_list : list[Union[str, os.PathLike[str]]]
        List of paths to V gene premixing source plate layout files.
    cdr3_alpha_product_dict : collections.defaultdict[lambda: collections.defaultdict[lambda: collections.defaultdict[list]]]
        Nested dictionary representing CDR3-J alpha PCR products organized by plates and wells.
        The first-level keys represent plates, the second-level keys represent the wells of the
        well-specific Fw + Rev ortho primer combinations used for PCR amplification, and the third-level
        keys represent wells whose CDR3-J alpha oligos have one or more PCR products amplified by the
        second-level primers. The third-level keys map to lists of :class:`pydna.amplicon.Amplicon` objects.
        Must be generated by `pydna_cdr3j_ortho_primer_amp_pcr`.
    cdr3_beta_product_dict : collections.defaultdict[lambda: collections.defaultdict[lambda: collections.defaultdict[list]]]
        Nested dictionary representing CDR3-J beta PCR products organized by plates and wells.
        The first-level keys represent plates, the second-level keys represent the wells of the
        well-specific Fw + Rev ortho primer combinations used for PCR amplification, and the third-level
        keys represent wells whose CDR3-J beta oligos have one or more PCR products amplified by the
        second-level primers. The third-level keys map to lists of :class:`pydna.amplicon.Amplicon` objects.
        Must be generated by `pydna_cdr3j_ortho_primer_amp_pcr`.
    run_sub_pool_file : Union[str, os.PathLike[str]]
        Path to the oligo pool order `.xlsx` sheet.
    run_path : Union[str, os.PathLike[str]]
        Path to the standardized assembly run directory initialized with `init_tcr_assembly_run_dir`.
    numbered_run_name : str
        Name of the assembly run.
    run_sub_pool : int
        Sub-pool number of the oligo order sheet.
    run_mode : str, default="oligo_order"
        Defines whether the function runs in planning (`oligo_order`), wet-lab execution (`run_assembly`),
        or dry-run (`simulation`) mode.

    Returns
    -------
    ligation_products_dict : collections.defaultdict[lambda: collections.defaultdict[lambda: collections.defaultdict[list]]]
        Nested dictionary representing circular Golden Gate TCR ligation products.
        - First-level keys: plates
        - Second-level keys: wells of the well-specific Fw + Rev primer combination used for PCR
        - Third-level keys: wells whose assigned CDR3-J alpha and beta oligos produce PCR products that can be ligated
          with the second-level well-specific V gene pre-mix.
        Values are :class:`pydna.dseq.Dseq` objects representing circular Golden Gate ligation products. A third-level
        key is only present if a single circular product can be ligated from the associated PCR products. The log file
        records cases where multiple PCR products exist for a third-level well.

    Examples
    --------
    >>> ligation_products_dict = pydna_golden_gate_tcr_assembly(
    ...     v_genes_premix_dispense_run_file=run_path + '/v_genes_premix_dispense/' + numbered_run_name + '_echo_dispense.csv',
    ...     v_genes_source_run_files_list=[run_path + '/v_genes_premix_dispense/' + source_name + '.xlsx'],
    ...     cdr3_alpha_product_dict=cdr3_alpha_product_dict,
    ...     cdr3_beta_product_dict=cdr3_beta_product_dict,
    ...     run_sub_pool_file=run_sub_pool_file,
    ...     run_path=run_path,
    ...     numbered_run_name=numbered_run_name,
    ...     run_sub_pool=run_sub_pool
    ... )

    Notes
    -----
    Pseudo-example of `ligation_products_dict`:
    {
        '1': {
              'A1': {
                     'A1': Dseq(o6985) # If a circular ligation product can be assembled for any of the well-specific CDR3-J amplification product(s) + well-specific pre-mixed V genes, that well will have a key that contains a :class:`pydna.dseq.Dseq` object storing the circular ligation product.
                           CCAA..GGCG
                           GGTT..CCGC
              },
              'A2': {
                     'A2': Dseq(o6985) # If a circular ligation product can be assembled for any of the well-specific CDR3-J amplification product(s) + well-specific pre-mixed V genes, that well will have a key that contains a :class:`pydna.dseq.Dseq` object storing the circular ligation product.
                           CCAA..GGCG
                           GGTT..CCGC
              },
              ...
        },
        '2': {
              'A1': {  # If no circular ligation product can be assembled for any of the well-specific CDR3-J amplification product(s) + well-specific pre-mixed V genes, Fw + Rev ortho primer well values will be empty.
              },
              'A2': {
                     'A2': Dseq(o6985) # If a circular ligation product can be assembled for any of the well-specific CDR3-J amplification product(s) + well-specific pre-mixed V genes, that well will have a key that contains a :class:`pydna.dseq.Dseq` object storing the circular ligation product.
                           CCAA..GGCG
                           GGTT..CCGC
              },
              ...
        },
        ...
    }
    """
    if not list(cdr3_alpha_product_dict.keys()) == list(cdr3_beta_product_dict.keys()):
        raise Exception("Plate keys in alpha and beta product dict are out of sync!!")

    run_sub_pool_df = pd.read_excel(run_sub_pool_file)
    run_sub_pool_df["plate"] = run_sub_pool_df["name"].str.split("_").str[1].astype(int)

    # Only the first line without header = 3 and on_bad_lines = 'skip' is needed when an echo v gene premix .csv would be read:
    column_names = ["Source Well", "Target Well", "Volume [uL]", "Liquid Name", "Additional Volume Per Source Well", "a", "b", "c"]
    v_genes_dispense_run_df = pd.read_csv(v_genes_premix_dispense_run_file, names=column_names, header=None, on_bad_lines="skip")
    v_genes_dispense_run_df.rename({"Volume [uL]": "Transfer Volume"}, axis=1, inplace=True)

    header_idx_list = []
    for well, idx in zip(v_genes_dispense_run_df.loc[:, "Source Well"], v_genes_dispense_run_df.loc[:, "Source Well"].index):
        if well == "V_gene_pre-mixing":
            header_idx_list.append(idx)

    index_ranges_to_remove = [list(np.arange(header_idx, header_idx + 3)) for header_idx in header_idx_list]
    v_genes_dispense_run_df["Source Plate Name"] = pd.Series([np.nan] * len(v_genes_dispense_run_df), dtype=object)
    v_genes_dispense_run_df["Target Plate Name"] = pd.Series([np.nan] * len(v_genes_dispense_run_df), dtype=object)

    for current_plate_indices, next_plate_indices in zip(index_ranges_to_remove, index_ranges_to_remove[1:] + [None]):
        source_plate = v_genes_dispense_run_df.iloc[current_plate_indices[1], 1]
        target_plate = int(v_genes_dispense_run_df.iloc[current_plate_indices[1], 5].split("_")[1])

        if next_plate_indices:
            v_genes_dispense_run_df.loc[current_plate_indices[2] + 1 : next_plate_indices[0] - 1, "Source Plate Name"] = source_plate
            v_genes_dispense_run_df.loc[current_plate_indices[2] + 1 : next_plate_indices[0] - 1, "Target Plate Name"] = target_plate
        else:
            v_genes_dispense_run_df.loc[current_plate_indices[2] + 1 : v_genes_dispense_run_df.shape[0], "Source Plate Name"] = source_plate
            v_genes_dispense_run_df.loc[current_plate_indices[2] + 1 : v_genes_dispense_run_df.shape[0], "Target Plate Name"] = target_plate

    indices_to_remove = [idx for sublist in index_ranges_to_remove for idx in sublist]
    v_genes_dispense_run_df.drop(indices_to_remove, axis=0, inplace=True)
    v_genes_dispense_run_df.reset_index(inplace=True, drop=True)

    v_genes_dispense_run_df.drop([column for column in v_genes_dispense_run_df.columns if column in ["a", "b", "c"]], axis=1, inplace=True)
    v_genes_dispense_run_df.drop(["Additional Volume Per Source Well"], axis=1, inplace=True)
    v_genes_dispense_run_df.reset_index(inplace=True, drop=True)

    v_genes_dispense_source_df_dict = collections.defaultdict(pd.DataFrame)
    for source_df_file in v_genes_source_run_files_list:
        v_genes_dispense_source_df_dict[os.path.basename(source_df_file).split(".")[0]] = pd.read_excel(source_df_file, index_col=0)

    pMX_S1_Kana_2_record = list(
        SeqIO.parse(
            os.path.join(
                tcr_toolbox_data_path,
                "tcr_toolbox_datasets",
                "tcr_assembly",
                "plasmids",
                "pMX_S1_Kana_2_cloning.dna",
            ),
            "snapgene",
        ).records
    )
    pMX_S1_Kana_2_record = Dseq(str(pMX_S1_Kana_2_record[0].seq), circular=True)

    pMX_muTCRB_record = list(
        SeqIO.parse(
            os.path.join(
                tcr_toolbox_data_path,
                "tcr_toolbox_datasets",
                "tcr_assembly",
                "plasmids",
                "muTRBC.dna",
            ),
            "snapgene",
        ).records
    )
    pMX_muTCRB_record = Dseq(str(pMX_muTCRB_record[0].seq), circular=True)

    ligation_products_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))

    with open(os.path.join(run_path, "pydna_logs", numbered_run_name + "_run_sub_pool_" + str(run_sub_pool) + "_golden_gate_log_" + run_mode + ".txt"), "w") as fout:
        for plate in run_sub_pool_df["plate"].unique():
            fout.write("\n\n-------------------------------------------------------------------------------------------------------\n")
            fout.write("Plate: " + str(plate) + "\n")

            tmp_v_genes_plate_df = v_genes_dispense_run_df[v_genes_dispense_run_df["Target Plate Name"] == plate].copy()

            if not list(cdr3_alpha_product_dict[plate].keys()) == list(cdr3_beta_product_dict[plate].keys()):
                raise Exception("Ortho primer combination keys in alpha and beta plate product dict are out of sync!!")

            for ortho_well in cdr3_alpha_product_dict[plate].keys():
                tmp_error_str = ""

                if (len(cdr3_alpha_product_dict[plate][ortho_well].keys()) > 2) | (len(cdr3_beta_product_dict[plate][ortho_well].keys()) > 2):
                    tmp_error_str += "CDR3 ortho PCR amplified more than one product!\n"

                # TRAV:
                v_gene_well_list = tmp_v_genes_plate_df.loc[
                    (tmp_v_genes_plate_df["Target Well"] == ortho_well) & (tmp_v_genes_plate_df["Liquid Name"].str.startswith("TRAV")), "Source Well"
                ].values.tolist()

                source_name_list = tmp_v_genes_plate_df.loc[
                    (tmp_v_genes_plate_df["Target Well"] == ortho_well) & (tmp_v_genes_plate_df["Liquid Name"].str.startswith("TRAV")), "Source Plate Name"
                ].values.tolist()
                if len(source_name_list) > 1:
                    raise Exception("Multiple source plates match TRAV ortho well. This is not allowed!", source_name_list)

                v_genes_dispense_source_df = v_genes_dispense_source_df_dict[source_name_list[0]]

                tmp_v_gene_test_list = []
                for v_gene_well in v_gene_well_list:
                    tmp_v_gene_test_list.append(v_genes_dispense_source_df.loc[v_genes_dispense_source_df["well"] == v_gene_well, "TRAV_or_TRBV"].values[0])
                if len(set(tmp_v_gene_test_list)) != 1:
                    raise Exception("Not all V gene source wells of ortho primer combination:", ortho_well, "contain the same V gene!")

                trav_name_str_tmp = v_genes_dispense_source_df.loc[v_genes_dispense_source_df["well"] == v_gene_well_list[0], "TRAV_or_TRBV"].values[0] + ".dna"

                # TRBV:
                v_gene_well_list = tmp_v_genes_plate_df.loc[
                    (tmp_v_genes_plate_df["Target Well"] == ortho_well) & (tmp_v_genes_plate_df["Liquid Name"].str.startswith("TRBV")), "Source Well"
                ].values.tolist()

                source_name_list = tmp_v_genes_plate_df.loc[
                    (tmp_v_genes_plate_df["Target Well"] == ortho_well) & (tmp_v_genes_plate_df["Liquid Name"].str.startswith("TRBV")), "Source Plate Name"
                ].values.tolist()
                if len(source_name_list) > 1:
                    raise Exception("Multiple source plates match TRBV ortho well. This is not allowed!", source_name_list)

                v_genes_dispense_source_df = v_genes_dispense_source_df_dict[source_name_list[0]]

                tmp_v_gene_test_list = []
                for v_gene_well in v_gene_well_list:
                    tmp_v_gene_test_list.append(v_genes_dispense_source_df.loc[v_genes_dispense_source_df["well"] == v_gene_well, "TRAV_or_TRBV"].values[0])
                if len(set(tmp_v_gene_test_list)) != 1:
                    raise Exception("Not all V gene source wells of ortho primer combination:", ortho_well, "contain the same V gene!")

                trbv_name_str_tmp = v_genes_dispense_source_df.loc[v_genes_dispense_source_df["well"] == v_gene_well_list[0], "TRAV_or_TRBV"].values[0] + ".dna"

                trav_name_str_tmp = trav_name_str_tmp.replace("/", "_")
                trbv_name_str_tmp = trbv_name_str_tmp.replace("/", "_")

                trav = list(
                    SeqIO.parse(
                        os.path.join(
                            tcr_toolbox_data_path,
                            "tcr_toolbox_datasets",
                            "tcr_assembly",
                            "plasmids",
                            "TRV_plasmid_stocks_ordered_at_Twist",
                            trav_name_str_tmp,
                        ),
                        "snapgene",
                    ).records
                )
                trav = Dseq(str(trav[0].seq), circular=True)

                trbv = list(
                    SeqIO.parse(
                        os.path.join(
                            tcr_toolbox_data_path,
                            "tcr_toolbox_datasets",
                            "tcr_assembly",
                            "plasmids",
                            "TRV_plasmid_stocks_ordered_at_Twist",
                            trbv_name_str_tmp,
                        ),
                        "snapgene",
                    ).records
                )
                trbv = Dseq(str(trbv[0].seq), circular=True)

                for tcr_well in cdr3_alpha_product_dict[plate][ortho_well].keys():
                    if not list(cdr3_alpha_product_dict[plate][ortho_well].keys()) == list(cdr3_beta_product_dict[plate][ortho_well].keys()):
                        raise Exception("TCR well keys in alpha and beta ortho primer combination plate dict are out of sync!!")

                    number_alpha_bbsi_cut_products = len(cdr3_alpha_product_dict[plate][ortho_well][tcr_well][0].seq.cut(Restriction.BbsI))

                    number_beta_bbsi_cut_products = len(cdr3_beta_product_dict[plate][ortho_well][tcr_well][0].seq.cut(Restriction.BbsI))
                    if (number_alpha_bbsi_cut_products == 3) & (number_beta_bbsi_cut_products == 3):
                        try:
                            ligation_product = (
                                trbv.cut(Restriction.BbsI)[1]
                                + cdr3_beta_product_dict[plate][ortho_well][tcr_well][0].seq.cut(Restriction.BbsI)[1]
                                + pMX_muTCRB_record.cut(Restriction.BbsI)[1]
                                + trav.cut(Restriction.BbsI)[1]
                                + cdr3_alpha_product_dict[plate][ortho_well][tcr_well][0].seq.cut(Restriction.BbsI)[1]
                                + pMX_S1_Kana_2_record.cut(Restriction.BbsI)[0]
                            ).looped()

                            if ligation_product.circular:
                                ligation_products_dict[plate][ortho_well][tcr_well] = ligation_product
                            else:
                                tmp_error_str += "BbsI restriction products compatible sticky ends did not ligate into a circular vector!\n"
                        except Exception as err:
                            tmp_error_str += "TCR well: " + str(tcr_well) + "\n"
                            tmp_error_str += str(err) + " " + str(trav_name_str_tmp) + " " + str(trbv_name_str_tmp) + "\n"

                    else:
                        tmp_error_str += (
                            "TCR well: "
                            + str(tcr_well)
                            + " has more or less than 3 BbsI restriction products!\n"
                            + str(cdr3_alpha_product_dict[plate][ortho_well][tcr_well][0].seq.cut(Restriction.BbsI))
                            + "\n"
                            + str(cdr3_beta_product_dict[plate][ortho_well][tcr_well][0].seq.cut(Restriction.BbsI))
                            + "\n"
                        )
                        tmp_error_str += "Running multiple BbsI cut product ligation simulation!\n"

                        if number_alpha_bbsi_cut_products != number_beta_bbsi_cut_products:
                            raise Exception(
                                f"{plate}, {ortho_well} have different numbers of alpha and beta BbsI cut products!\n alpha {cdr3_alpha_product_dict[plate][ortho_well][tcr_well][0].seq} has {number_alpha_bbsi_cut_products}\nand beta {cdr3_beta_product_dict[plate][ortho_well][tcr_well][0].seq} has {number_beta_bbsi_cut_products} products.\n We do not support ligation simulation with non-matching numbers of alpha and beta BbsI cut products!"
                            )

                        for i, (cdr3_alpha_restriction_product, cdr3_beta_restriction_product) in enumerate(
                            zip(
                                cdr3_alpha_product_dict[plate][ortho_well][tcr_well][0].seq.cut(Restriction.BbsI),
                                cdr3_beta_product_dict[plate][ortho_well][tcr_well][0].seq.cut(Restriction.BbsI),
                            )
                        ):
                            try:
                                ligation_product = (
                                    trbv.cut(Restriction.BbsI)[1]
                                    + cdr3_beta_restriction_product
                                    + pMX_muTCRB_record.cut(Restriction.BbsI)[1]
                                    + trav.cut(Restriction.BbsI)[1]
                                    + cdr3_alpha_restriction_product
                                    + pMX_S1_Kana_2_record.cut(Restriction.BbsI)[0]
                                ).looped()

                                if ligation_product.circular:
                                    ligation_products_dict[plate][ortho_well][tcr_well] = ligation_product
                                else:
                                    tmp_error_str += "BbsI restriction products with compatible sticky ends did not ligate into a circular vector!\n"
                            except Exception as err:
                                tmp_error_str += "TCR well: " + str(tcr_well) + "\n"
                                tmp_error_str += str(err) + " " + str(trav_name_str_tmp) + " " + str(trbv_name_str_tmp) + "\n"

                if len(ligation_products_dict[plate][ortho_well].keys()) == 0:
                    tmp_error_str += "No Golden Gate product assembled at all!\n"

                if tmp_error_str:
                    fout.write("Ortho primer combination: " + str(ortho_well) + "\n")
                    fout.write(tmp_error_str)
                    fout.write("\n")

    tmp_dict = {int(plate): {well: {well: str(inner[well])} for well, inner in wells.items()} for plate, wells in ligation_products_dict.items()}

    with open(os.path.join(run_path, "pydna_ligation_dicts", numbered_run_name + "_ligation_dict_" + str(run_sub_pool) + ".json"), "w", encoding="utf-8") as fout:
        json.dump(tmp_dict, fout)

    return ligation_products_dict


def pydna_golden_gate_eblock_tcr_assembly(
    eblock_df_dict: collections.defaultdict[pd.DataFrame],
    group_df_dict: collections.defaultdict[pd.DataFrame],
    run_path: Union[str, os.PathLike[str]],
    numbered_run_name: str,
    run_mode: str = "oligo_order",
):
    if not len(set(group_df_dict.keys()).difference(set([name.split("_order")[0] for name in eblock_df_dict.keys()]))) == 0:
        raise Exception("group_df_dict keys do not equal eblock_df_dict keys!")

    pMX_S1_Kana_2_record = list(
        SeqIO.parse(
            os.path.join(
                tcr_toolbox_data_path,
                "tcr_toolbox_datasets",
                "tcr_assembly",
                "plasmids",
                "pMX_S1_Kana_2_cloning.dna",
            ),
            "snapgene",
        ).records
    )
    pMX_S1_Kana_2_record = Dseq(str(pMX_S1_Kana_2_record[0].seq), circular=True)

    pMX_muTCRB_record = list(
        SeqIO.parse(
            os.path.join(
                tcr_toolbox_data_path,
                "tcr_toolbox_datasets",
                "tcr_assembly",
                "plasmids",
                "muTRBC.dna",
            ),
            "snapgene",
        ).records
    )
    pMX_muTCRB_record = Dseq(str(pMX_muTCRB_record[0].seq), circular=True)

    ligation_products_dict = collections.defaultdict(lambda: collections.defaultdict(list))

    with open(os.path.join(run_path, "pydna_logs", numbered_run_name + "golden_gate_log_" + run_mode + ".txt"), "w") as fout:
        for plate in eblock_df_dict.keys():
            fout.write("\n\n-------------------------------------------------------------------------------------------------------\n")
            fout.write("Plate: " + str(plate) + "\n")

            tmp_eblock_plate_df = eblock_df_dict[plate].copy()
            tmp_eblock_plate_df.dropna(subset=["Name"], inplace=True)

            tmp_v_genes_plate_df = group_df_dict[plate.split("_order")[0]].copy()
            tmp_v_genes_plate_df.dropna(subset=["TCR name"], inplace=True)

            if not all(tmp_eblock_plate_df["Name"] == tmp_v_genes_plate_df["TCR name"]):
                raise Exception("TCR names in eblock order sheet are out of sync with TCR names in V gene csv file!")

            if not all(tmp_eblock_plate_df["Well Position"] == tmp_v_genes_plate_df["well"]):
                raise Exception("Wells in eblock order sheet are out of sync with wells in V gene csv file!")

            for well, sequence in zip(tmp_eblock_plate_df["Well Position"], tmp_eblock_plate_df["Sequence"]):
                tmp_error_str = ""

                trav_name_str_tmp = tmp_v_genes_plate_df.loc[tmp_v_genes_plate_df["well"] == well, "TRAV"].values[0]
                if "*" in trav_name_str_tmp:
                    if "01" in trav_name_str_tmp.split("*")[1]:
                        raise Exception("We do not currently support non-*01 alleles yet!")
                    trav_name_str_tmp = trav_name_str_tmp.split("*")[0]
                trav_name_str_tmp = trav_name_str_tmp + "*01.dna"

                trbv_name_str_tmp = tmp_v_genes_plate_df.loc[tmp_v_genes_plate_df["well"] == well, "TRBV"].values[0]
                if "*" in trbv_name_str_tmp:
                    if "01" in trbv_name_str_tmp.split("*")[1]:
                        raise Exception("We do not currently support non-*01 alleles yet!")
                    trbv_name_str_tmp = trbv_name_str_tmp.split("*")[0]
                trbv_name_str_tmp = trbv_name_str_tmp + "*01.dna"

                trav_name_str_tmp = trav_name_str_tmp.replace("/", "_")
                trbv_name_str_tmp = trbv_name_str_tmp.replace("/", "_")

                trav = list(
                    SeqIO.parse(
                        tcr_toolbox_data_path
                        + "/tcr_toolbox_datasets/tcr_assembly/plasmids/TRV_plasmid_stocks_ordered_at_Twist/"
                        + trav_name_str_tmp,
                        "snapgene",
                    ).records
                )
                trav = Dseq(str(trav[0].seq), circular=True)

                trbv = list(
                    SeqIO.parse(
                        tcr_toolbox_data_path
                        + "/tcr_toolbox_datasets/tcr_assembly/plasmids/TRV_plasmid_stocks_ordered_at_Twist/"
                        + trbv_name_str_tmp,
                        "snapgene",
                    ).records
                )
                trbv = Dseq(str(trbv[0].seq), circular=True)

                if len(Dseq(sequence).cut(Restriction.BbsI)) == 5:
                    try:
                        ligation_product = (
                            trbv.cut(Restriction.BbsI)[1]
                            + Dseq(sequence).cut(Restriction.BbsI)[3]
                            + pMX_muTCRB_record.cut(Restriction.BbsI)[1]
                            + trav.cut(Restriction.BbsI)[1]
                            + Dseq(sequence).cut(Restriction.BbsI)[1]
                            + pMX_S1_Kana_2_record.cut(Restriction.BbsI)[0]
                        ).looped()

                        if ligation_product.circular:
                            ligation_products_dict[plate][well] = ligation_product
                        else:
                            tmp_error_str += "BbsI restriction products compatible sticky ends did not ligate into a circular vector!\n"
                    except Exception as err:
                        tmp_error_str += str(err) + " " + str(trav_name_str_tmp) + " " + str(trbv_name_str_tmp) + "\n"
                else:
                    tmp_error_str += "Has more or less than 5 BbsI restriction products!\n" + str(Dseq(sequence).cut(Restriction.BbsI)) + "\n"

                if not ligation_products_dict[plate][well]:
                    tmp_error_str += "No Golden Gate product assembled at all!\n"

                if tmp_error_str:
                    fout.write("Well: " + str(well) + "\n")
                    fout.write(tmp_error_str)
                    fout.write("\n")

    return ligation_products_dict


def check_translation_assembled_tcr_ligation_products(
    plates_df_dict: pd.DataFrame,
    ligation_products_dict: collections.defaultdict[lambda: collections.defaultdict[lambda: collections.defaultdict[list]]],
    run_path: Union[str, os.PathLike[str]],
    numbered_run_name: str,
    run_sub_pool: int,
    run_mode: str = "oligo_order",
):
    """
    Compare the amino acid translation of simulated TCR Golden Gate ligation products 
    in wells with the reconstructed TCR amino acid sequence that should be assembled 
    in those wells. 

    The function writes a ligation reconstruction translation comparison log `.txt` file. 
    The log file will be empty if all simulated TCR ligation product translations match 
    the reconstructed TCR sequences. Any discrepancies are logged with plate and well 
    identifiers.

    Prerequisite function(s)
    ------------------------
    0. init_tcr_assembly_run_dir
    1. filter_tcr_df_on_v_genes_in_stock
    2. add_cdr3j_alpha_beta_order_columns
    3. concatenate_separate_tcr_df_csvs_into_run_tcr_df
       Run only if there are separate `.csv` files storing TCRs for the assembly run.
    4. add_standardized_tcr_assembly_well_plate_layout
    5. add_plate_ortho_primer_combinations_to_cdr3j_seqs
    6. make_tcr_assembly_run_cdr3j_seq_order_sheets
    7. make_idot_v_gene_premixing_dispense_csv
    8. pydna_plate_sub_pool_amp_pcr
    9. pydna_cdr3j_ortho_primer_amp_pcr
    10. pydna_golden_gate_tcr_assembly
    11. reconstruct_full_tcr

    Parameters
    ----------
    plates_df_dict : dict
        Dictionary where keys are plate numbers and values are DataFrames storing plate data. 
        Must be generated by `make_tcr_assembly_run_cdr3j_seq_order_sheets` with a 
        'reconstructed_tcrab_full_aa' column added via `reconstruct_full_tcr`. 
    ligation_products_dict : collections.defaultdict[lambda: collections.defaultdict[lambda: collections.defaultdict[list]]]
        Nested dictionary of circular Golden Gate ligation products (output from 
        `pydna_golden_gate_tcr_assembly`).
    run_path : Union[str, os.PathLike[str]]
        Path to the standardized assembly run directory initialized with `init_tcr_assembly_run_dir`.
    numbered_run_name : str
        Name of the assembly run. 
    run_sub_pool : int
        Sub-pool number of the oligo order sheet.
    run_mode : str, default="oligo_order"
        Defines whether the function runs in planning (`oligo_order`), wet-lab execution 
        (`run_assembly`), or dry-run (`simulation`) mode.

    Returns
    -------
    correct_translation_ligation_products_dict : collections.defaultdict[lambda: collections.defaultdict[lambda: collections.defaultdict[list]]]
        Nested dictionary indicating whether ligation products match the reconstructed TCR sequences.
        - First-level keys: plates
        - Second-level keys: wells of the well-specific Fw + Rev ortho primer combination used for PCR
        - Third-level keys: wells whose assigned CDR3-J alpha and beta oligos produced PCR products that 
          were used to generate a circular ligation product
        Values are `True` if the ligation product matches the reconstructed TCR sequence, `False` otherwise.
        A third-level key only exists if a ligation product could be formed from the PCR products and pre-mixed 
        V genes.

    Examples
    --------
    >>> correct_translation_ligation_products_dict = \
    ...     check_translation_assembled_tcr_ligation_products(
    ...         plates_df_dict=plates_df_dict,
    ...         ligation_products_dict=ligation_products_dict,
    ...         run_path=run_path,
    ...         numbered_run_name=numbered_run_name,
    ...         run_sub_pool=run_sub_pool
    ...     )

    Notes
    -----
    Pseudo-example of `correct_translation_ligation_products_dict`:
    {
        '1': {
              'A1': {
                     'A1': True  # True if ligation product in ligation_products_dict input dictionary is not equal to reconstructed TCR of well 'A1'.
              },
              'A2': {
                     'A2': True  # True if ligation product in ligation_products_dict input dictionary is not equal to reconstructed TCR of well 'A2'.
              },
              ...
        },
        '2': {
              'A1': {  # No third level well key and value if no ligation product could be ligated and stored in ligation_products_dict input dictonary.
              },
              'A2': {
                     'A2': False   # False if ligation product in ligation_products_dict input dictionary is not equal to reconstructed TCR of well 'A2'.
              },
              ...
        },
        ...
    }
    """

    five_end = "CCAAGCCGCC"
    three_end = "GTCGACGATAAAATAAAAGATTTTA"

    correct_translation_ligation_products_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))

    with open(
        os.path.join(run_path, "pydna_logs", numbered_run_name + "_run_sub_pool_" + str(run_sub_pool) + "_ligation_translation_vs_reconstruction_log_" + run_mode + ".txt"), "w"
    ) as fout:
        for plate in ligation_products_dict.keys():
            fout.write("\n\n-------------------------------------------------------------------------------------------------------\n")
            fout.write("Plate: " + str(plate) + "\n")

            for ortho_well in ligation_products_dict[plate].keys():
                tmp_error_str = ""

                if len(ligation_products_dict[plate][ortho_well].keys()) == 0:
                    tmp_error_str += "Ortho primer combination well: " + str(ortho_well) + " has no assembled Golden Gate product at all!\n"
                    fout.write(tmp_error_str)
                    fout.write("\n")
                    continue

                if not isinstance(plates_df_dict[plate].loc[plates_df_dict[plate]["well"] == ortho_well, "reconstructed_tcrab_full_aa"].values[0], str):
                    tmp_error_str += "Ortho primer combination well: " + str(ortho_well) + " reference TCRab was not reconstructed!\n"
                    fout.write(tmp_error_str)
                    fout.write("\n")
                    continue

                for tcr_well in ligation_products_dict[plate][ortho_well].keys():
                    if len(ligation_products_dict[plate][ortho_well][tcr_well]) == 0:
                        tmp_error_str += f"{plate}, {ortho_well} does not have a ligation product!\n"
                        continue

                    start = ligation_products_dict[plate][ortho_well][tcr_well].find(five_end) + len(five_end)
                    end = ligation_products_dict[plate][ortho_well][tcr_well].find(three_end)

                    if ligation_products_dict[plate][ortho_well][tcr_well][start:end].translate(to_stop=True)[0] != "M":
                        ligation_product_tmp = ligation_products_dict[plate][ortho_well][tcr_well][start + 1 : end].translate(to_stop=True)
                    else:
                        ligation_product_tmp = ligation_products_dict[plate][ortho_well][tcr_well][start:end].translate(to_stop=True)

                    correct_translation_ligation_products_dict[plate][ortho_well][tcr_well] = (
                        str(ligation_product_tmp) == plates_df_dict[plate].loc[plates_df_dict[plate]["well"] == ortho_well, "reconstructed_tcrab_full_aa"].values[0]
                    )

                    if not correct_translation_ligation_products_dict[plate][ortho_well][tcr_well]:
                        tmp_error_str += "Ligation product translation does not match the full reconstructed TCRab sequence of ortho primer combination: " + str(ortho_well) + "!\n"
                        tmp_error_str += "TCR well: " + str(tcr_well) + ", tcr_well == ortho_well: " + str(tcr_well == ortho_well) + "\n"
                        tmp_error_str += (
                            "Levenshtein Ratio between Ligation translation and Reconstruction aa: "
                            + str(
                                round(
                                    levenshtein_ratio(
                                        str(ligation_product_tmp),
                                        str(plates_df_dict[plate].loc[plates_df_dict[plate]["well"] == ortho_well, "reconstructed_tcrab_full_aa"].values[0]),
                                    ),
                                    4,
                                )
                            )
                            + "\n"
                        )
                        tmp_error_str += (
                            "Levenshtein dist. between Ligation translation and Reconstruction aa: "
                            + str(
                                round(
                                    levenshtein_distance(
                                        str(ligation_product_tmp),
                                        str(plates_df_dict[plate].loc[plates_df_dict[plate]["well"] == ortho_well, "reconstructed_tcrab_full_aa"].values[0]),
                                    ),
                                    2,
                                )
                            )
                            + "\n"
                        )
                        tmp_error_str += "Ligation product translation aa: " + str(ligation_product_tmp) + "\n"
                        tmp_error_str += (
                            "Reconstruction aa: " + str(plates_df_dict[plate].loc[plates_df_dict[plate]["well"] == ortho_well, "reconstructed_tcrab_full_aa"].values[0]) + "\n"
                        )

                if tmp_error_str:
                    fout.write(tmp_error_str)
                    fout.write("\n")

    return correct_translation_ligation_products_dict


def check_translation_assembled_eblock_tcr_ligation_products(
    group_df_dict: pd.DataFrame,
    ligation_products_dict: collections.defaultdict[lambda: collections.defaultdict[list]],
    run_path: Union[str, os.PathLike[str]],
    numbered_run_name: str,
    run_mode: str = "oligo_order",
):
    five_end = "CCAAGCCGCC"
    three_end = "GTCGACGATAAAATAAAAGATTTTA"

    correct_translation_ligation_products_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))

    with open(os.path.join(run_path, "pydna_logs" + numbered_run_name + "_ligation_translation_vs_reconstruction_log_" + run_mode + ".txt"), "w") as fout:
        for plate in ligation_products_dict.keys():
            group_plate = plate.split("_order")[0]
            fout.write("\n\n-------------------------------------------------------------------------------------------------------\n")
            fout.write("Plate: " + str(group_plate) + "\n")

            for well in ligation_products_dict[plate].keys():
                tmp_error_str = ""

                if not ligation_products_dict[plate][well]:
                    tmp_error_str += "Well: " + str(well) + " has no assembled Golden Gate product at all!\n"
                    fout.write(tmp_error_str)
                    fout.write("\n")
                    continue

                if not isinstance(group_df_dict[group_plate].loc[group_df_dict[group_plate]["well"] == well, "reconstructed_tcrab_full_aa"].values[0], str):
                    tmp_error_str += "Well: " + str(well) + " reference TCRab was not reconstructed!\n"
                    fout.write(tmp_error_str)
                    fout.write("\n")
                    continue

                start = ligation_products_dict[plate][well].find(five_end) + len(five_end)
                end = ligation_products_dict[plate][well].find(three_end)

                if ligation_products_dict[plate][well][start:end].translate(to_stop=True)[0] != "M":
                    ligation_product_tmp = ligation_products_dict[plate][well][start + 1 : end].translate(to_stop=True)
                else:
                    ligation_product_tmp = ligation_products_dict[plate][well][start:end].translate(to_stop=True)

                correct_translation_ligation_products_dict[group_plate][well] = (
                    str(ligation_product_tmp) == group_df_dict[group_plate].loc[group_df_dict[group_plate]["well"] == well, "reconstructed_tcrab_full_aa"].values[0]
                )

                if not correct_translation_ligation_products_dict[group_plate][well]:
                    tmp_error_str += "Ligation product translation does not match the full reconstructed TCRab sequence of well: " + str(well) + "!\n"
                    tmp_error_str += (
                        "Levenshtein Ratio between Ligation translation and Reconstruction aa: "
                        + str(
                            round(
                                levenshtein_ratio(
                                    str(ligation_product_tmp),
                                    str(group_df_dict[group_plate].loc[group_df_dict[group_plate]["well"] == well, "reconstructed_tcrab_full_aa"].values[0]),
                                ),
                                4,
                            )
                        )
                        + "\n"
                    )
                    tmp_error_str += (
                        "Levenshtein dist. between Ligation translation and Reconstruction aa: "
                        + str(
                            round(
                                levenshtein_distance(
                                    str(ligation_product_tmp),
                                    str(group_df_dict[group_plate].loc[group_df_dict[group_plate]["well"] == well, "reconstructed_tcrab_full_aa"].values[0]),
                                ),
                                2,
                            )
                        )
                        + "\n"
                    )
                    tmp_error_str += "Ligation product translation aa: " + str(ligation_product_tmp) + "\n"
                    tmp_error_str += (
                        "Reconstruction aa: " + str(group_df_dict[group_plate].loc[group_df_dict[group_plate]["well"] == well, "reconstructed_tcrab_full_aa"].values[0]) + "\n"
                    )

                if tmp_error_str:
                    fout.write(tmp_error_str)
                    fout.write("\n")

    return correct_translation_ligation_products_dict


def generate_snapgene_files(tcr_list: list, output_path: str, run_path: Union[str, os.PathLike[str]], run_name: str):
    """
    Write Golden Gate nucleotide products of a selection of TCRs to SnapGene-readable
    GenBank files for Sanger sequencing alignment.

    Parameters
    ----------
    tcr_list : list
        List of TCRs to validate. Entries can have the following formats:
        - Short format: plate_well identifier, e.g., `20_C14`
        - Long format: full identifier with annotation, e.g., `3_20_C14_121_WTDNNCYLA`
          In the long format, only the plate and well information (second and third elements)
          will be extracted.
    output_path : str
        Path to the output folder where SnapGene `.gb` files will be saved.
    run_path : Union[str, os.PathLike[str]]
        Path to the standardized assembly run directory.

    Notes
    -----
    - The function reads all `*_ligation_dict_*.json` files from
      `run_path/pydna_ligation_dicts/` and merges them into a single dictionary.
    - Only TCRs present in `tcr_list` are processed and written to GenBank files.
    - Each output GenBank file includes features defined in `snapgene_feature_dict`
      with circular DNA topology.

    Raises
    ------
    ValueError
        If any string in `tcr_list` does not match the expected short or long format.
    """

    run_name = os.path.basename(os.path.normpath(run_path))

    # use format xx_plate_well_xx_xx or plate_well, other formats are not accepted
    for i, tcr in enumerate(tcr_list):
        # Check for short format: number_letter+number (e.g., 20_C14)
        if re.match(r"^\d+_[A-Za-z]+\d+$", tcr):
            pass
        # Check for long format: number_number_letter+number_number_anything (e.g., 3_20_C14_121_WTDNNCYLA)
        elif re.match(r"^\d+_\d+_[A-Za-z]+\d+_\d+_.*$", tcr):
            tcr_list[i] = "_".join(tcr.split("_")[1:3])
        else:
            raise ValueError(f"String '{tcr}' does not match expected format")

    # create empty defaultdict
    tcr_validation_dict = collections.defaultdict(list)
    tmp_dict = collections.defaultdict(list)

    # read all ligation_product.json dictionaries for each subpool and make one dictionary out of it
    for json_path in glob.glob(os.path.join(run_path, "pydna_ligation_dicts", run_name + "_ligation_dict_*.json")):
        with open(json_path, "r", encoding="utf-8") as read_file:
            tmp_dict_json = json.load(read_file)

        # convert plate numbers back to int
        tmp_dict_json = {int(plate): {well: {well: Dseqrecord(inner[well], circular=True)} for well, inner in wells.items()} for plate, wells in tmp_dict_json.items()}
        # concat
        tmp_dict = tmp_dict | tmp_dict_json

    # create empty defaultdict
    tcr_validation_dict = collections.defaultdict(list)

    # append well, plate and sequence to empty defaultdict if it is in list.
    for plate, well in [tcr.split("_") for tcr in tcr_list]:
        tcr_validation_dict[plate + "_" + well].append(tmp_dict[int(plate)][well][well])

    # write genbank file
    for id, tcr_list in tcr_validation_dict.items():
        tcr = tcr_list[0]
        for name, seq in snapgene_feature_dict.items():
            start = tcr.seq.find(seq)
            if start != -1:
                end = start + len(seq)
                feature = SeqFeature(FeatureLocation(start, end), type="misc_feature", qualifiers={"label": name})
                tcr.features.append(feature)

        tcr.annotations["molecule_type"] = "DNA"
        tcr.annotations["topology"] = "circular"

        SeqIO.write(tcr, os.path.join(output_path, run_name + "_" + id + ".gb"), "genbank")
