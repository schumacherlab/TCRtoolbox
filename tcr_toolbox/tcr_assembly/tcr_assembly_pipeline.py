import collections
import glob
import os
import re
import shutil
import ast

import pandas as pd
from typing import Union
from dotenv import load_dotenv

from tcr_toolbox.sequencing_analysis.reference import generate_assembly_nanopore_nt_refs, generate_assembly_nt_refs
from tcr_toolbox.tcr_assembly.order_automation import (
    add_cdr3j_alpha_beta_order_columns,
    add_plate_ortho_primer_combinations_to_cdr3j_seqs,
    add_standardized_tcr_assembly_well_plate_layout,
    check_translation_assembled_tcr_ligation_products,
    expand_oligo_order_run_dir_to_run_assembly,
    filter_tcr_df_on_v_genes_in_stock,
    init_oligo_order_run_dir,
    make_idot_v_gene_premixing_dispense_csv,
    make_tcr_assembly_run_cdr3j_seq_order_sheets,
    pydna_cdr3j_ortho_primer_amp_pcr,
    pydna_golden_gate_tcr_assembly,
    pydna_plate_sub_pool_amp_pcr,
    translate_v_gene_input_to_imgt,
    validate_oligo_order_run_dir,
)
from tcr_toolbox.tcr_reconstruction.constants import pydna_tcr_golden_gate_constant_aa_seq_dict
from tcr_toolbox.tcr_reconstruction.reconstruction_tcr_assembly import reconstruct_tcrs_assembly
from tcr_toolbox.utils.plot_utils import plot_v_genes_per_well

load_dotenv()
tcr_toolbox_data_path = os.getenv("tcr_toolbox_data_path")


def print_run_mode_info(run_mode: str):
    """
    Prints a clear informational block about the pipeline mode.

    Parameters
    ----------
    mode : str
        The run mode, currently "simulation", "oligo_order" or "run_assembly" are implemented.
    """
    border = "=" * 80
    print(border)
    print(f"RUN MODE SELECTED: {run_mode.upper()}")
    print(border)

    if run_mode == "oligo_order":
        print(
            "In 'oligo_order' mode, the pipeline performs the following steps:\n"
            "  1. Initializes the standardized assembly run directory.\n"
            "  2. Assigns TCRs to plates.\n"
            "  3. Prepares the CDR3-J oligo pools.\n"
            "  4. Identifies which V gene stock tubes with too little volume\n"
            "     need to be replaced during oligo pool manufacturing.\n"
            "  5. Simulates in silico whether the assembly yields correct\n"
            "     final ligation products using pydna.\n\n"
            "This mode should be run first before run_mode = `run_assembly`."
        )
    elif run_mode == "run_assembly":
        print(
            "In 'run_assembly' mode, the pipeline performs the following steps:\n"
            "  1. Writes V gene premixing instruction files for robotic assembly.\n"
            "     to perform the assembly in the wet-lab.\n"
            "      Files are written to standardized assembly run directory initialized in the oligo_order mode.\n\n"
            "This mode should be run only after:\n"
            "  - The CDR3-J oligo pool has been manufactured.\n"
            "  - V gene stock tubes requiring replacement have been replaced."
        )
    elif run_mode == "simulation":
        print(
            "In 'simulation' mode, the full TCR assembly pipeline is run, after which\n"
            "the standardized assembly run directory is deleted.\n\n"
            "This mode is useful to check whether:\n"
            "  1. CDR3-J amino acid sequences can be codon optimized and ligated to V genes\n"
            "     in the Golden Gate assembly.\n"
            "  2. Unique V genes are in stock and, for those in stock, whether any tubes\n"
            "     need to be replaced to have enough volume.\n"
            "  3. The last plate per TCR group can at least be half filled.\n\n"
            "After checking all of the above, the TCR library input can be re-designed and re-evaluated.\n\n"
            "NOTE: When running in 'simulation' mode, the pipeline log file is written in the current directory."
        )
    else:
        raise NotImplementedError(f"There is no {run_mode}. Supported run_modes = `simulation`, `oligo_order`, `run_assembly`.")

    print(border)


def init_run_dir(run_mode: str, run_name: str | None = None, run_path: str | os.PathLike | None = None, **kwargs):
    if not tcr_toolbox_data_path:
        raise EnvironmentError(
            "Environment variable 'tcr_toolbox_data_path' is not set. Did you create or edit your .env file with the correct path?"
        )

    implemented_run_modes = ["simulation", "oligo_order", "run_assembly"]

    if run_mode not in implemented_run_modes:
        raise NotImplementedError(f"There is no {run_mode}. Supported run_modes = `simulation`, `oligo_order`, `run_assembly`.")

    run_path_local = run_path
    numbered_run_name = None
    if run_mode in ("oligo_order", "simulation"):
        if not run_path:
            numbered_run_name, run_path_local = init_oligo_order_run_dir(run_name=run_name)
            if run_mode == "oligo_order":
                log_file = os.path.join(run_path_local, f"{numbered_run_name}_oligo_order.log")

            if run_mode == "simulation":
                log_file = os.path.join(os.getcwd(), f"{numbered_run_name}_simulation.log")  # In simulation mode we cannot write the log file to the standardized run directory

            print("\n")
        else:
            raise Exception(
                f"A run_path was provided ({run_path}), but run_mode='{run_mode}' always initializes a new standardized run directory.\nOnly run_mode='run_assembly' is allowed to use an existing run_path, which must be a directory previously created by run_mode='oligo_order'."
            )

    elif run_mode == "run_assembly":
        if not run_path:
            raise Exception(
                "A path of a standardized run directory initialized by run_mode = 'oligo_order' needs to be provided in\nrun_mode = 'run_assembly' via the argument run_path."
            )
        if validate_oligo_order_run_dir(run_path):  # raises if invalid
            numbered_run_name, run_path_local = expand_oligo_order_run_dir_to_run_assembly(run_path)
            log_file = os.path.join(run_path_local, numbered_run_name + "_run_assembly.log")
            print("\n")
    else:
        raise NotImplementedError(f"There is no {run_mode}. Supported run_modes = `simulation`, `oligo_order`, `run_assembly`.")

    return numbered_run_name, run_path_local, log_file


def run_tcr_assembly(run_mode: str, numbered_run_name: str, run_path: str | None = None, **kwargs):
    print_run_mode_info(run_mode=run_mode)

    try:
        run_path_local = run_tcr_assembly_pipeline(run_mode=run_mode, numbered_run_name=numbered_run_name, run_path=run_path, **kwargs)

        return run_path_local

    except Exception:
        if run_mode == "simulation":
            shutil.rmtree(run_path)
        elif run_mode == "oligo_order":
            shutil.rmtree(run_path)
        elif run_mode == "run_assembly":
            shutil.rmtree(os.path.join(run_path, "v_genes_premix_dispense"))
            shutil.rmtree(os.path.join(run_path, "sequencing_quality_analysis"))
            for pydna_log in glob.glob(os.path.join(run_path, "pydna_logs")):
                if pydna_log.endswith("_run_assembly.txt"):
                    os.remove(pydna_log)
        raise


def run_tcr_assembly_pipeline(
    run_mode: str,
    run_tcr_csv: Union[str, os.PathLike[str]],
    numbered_run_name: Union[str, None],
    run_path: Union[str, os.PathLike[str], None] = None,
    grouping_col: str = "group",
    oligo_name_annotation_col_list: list | None = None,
    filter_cannot_be_codon_optimized: bool = True,
    filter_succeeding_nt_cys_104_beta: bool = False,
    allow_cdr3j_nt_duplicates: bool = False,
    steepness: float = 4,
    well_plate_size: int = 384,
    remove_longer_than_200_nt_oligos: bool = False,
    v_gene_premix_dispenser: str = "idot",
    v_gene_transfer_volume_nl: int = 150,
    max_v_gene_source_well_volume_nl: int = 70_000,
    v_gene_hamilton_ul: int = 21.12,
    water_hamilton_ul: int = 58.88,
    skip_barcoded_v_gene_stock_tube_vol_update: bool = False,
    min_v_gene_stock_tube_vol_ul: float = 50.00,
    generate_illumina_refs: bool = True,
    epitope_barcode_refs: Union[pd.DataFrame, str, os.PathLike[str], None] = None,
    epitope_name_col_name: str = "name",
    epitope_order_col_name: str = "sequence",
    trimmed_beta_model_tcr_dict: dict = None,
    tcr_read_length: int = 150,
    model_epitope_dict: dict = None,
    epitope_barcode_length: int = 18,
    generate_nanopore_ref: bool = True,
    add_number_of_negatives_by_mutating_full_refs: int = None,
    add_number_of_negatives_by_mutating_cdr3j_in_refs: int = None,
    add_number_of_negatives_by_mutating_v_gene_in_refs: int = None,
    min_nt_diff_negative_ref_seqs: int = None,
    max_nt_diff_negative_ref_seqs: int = None,
) -> Union[str, os.PathLike[str]]:
    """
    Run the complete TCR assembly pipeline, including oligo ordering, in silico simulation,
    V gene premix dispense file generation, and sequencing reference generation.

    Parameters
    ----------
    run_mode : str
    Specifies the mode of the pipeline. Options are:
    - 'simulation':
        The full TCR assembly pipeline is run in silico, after which the
        standardized assembly run directory is deleted. This mode checks:
          1. Whether CDR3-J amino acid sequences can be codon optimized and
             ligated to V genes via Golden Gate assembly.
          2. That all required V genes are in stock and identifies any tubes
             that need replacement to ensure sufficient volume.
          3. That the last plate per TCR group can be at least half filled.
        Notes:
          - Use this mode to check whether your run_tcr_csv input passes all assembly
            checks or whether your input needs to be re-designed.
          - The simulation log will be written in the directory in which the
            tcr-toolbox run-tcr-assembly CLI is run.
    - 'oligo_order':
        The pipeline performs the following steps:
          1. Initializes the standardized assembly run directory.
          2. Assigns TCRs to wells in plates based on the grouping column.
          3. Prepares CDR3-J oligo pools for synthesis.
          4. Identifies V gene stock tubes with insufficient volume
             that require replacement for oligo pool manufacturing.
          5. Simulates in silico whether assembly yields correct
             ligation products using `pydna`.
        Notes:
          - This mode should be run first before run_mode = `run_assembly`.
          - No V gene premix files for robotic assembly are generated in this mode.
    - 'run_assembly':
        The pipeline performs the following steps:
          1. Writes V gene premixing instruction files for robotic assembly
             to perform the assembly in the wet-lab.
             These files are written to the standardized assembly run directory
             that was initialized in 'oligo_order' mode.
        Notes:
          - This mode should be run only after:
              * The CDR3-J oligo pool has been manufactured.
              * V gene stock tubes requiring replacement have been replaced.

    run_tcr_csv : str or os.PathLike
        Path to the input CSV containing the following required columns:
        - "TRAV_IMGT": TRAV gene in IMGT format (10X Genomics format works).
        - "TRAJ_IMGT": TRAJ gene in IMGT format (10X Genomics format works).
        - "cdr3_alpha_aa": CDR3 alpha amino acid sequence.
        - "TRBV_IMGT": TRBV gene in IMGT format (10X Genomics format works).
        - "TRBJ_IMGT": TRBJ gene in IMGT format (10X Genomics format works).
        - "cdr3_beta_aa": CDR3 beta amino acid sequence.
        - "custom_name": Custom TCR name. For example, a TCR with custom name 'ywe_151_293' assigned to
           plate 2, well K14 will have the standardized assembly name:
           1_2_K14_253_ywe_151_293 → [oligo subpool]_[plate number]_[well]_[well number]_[custom name].
        - "group": Used to assign TCRs to plates in groups. Groups must fill at least half a plate
            if they cannot fully fill the last plate.

    numbered_run_name : str or None
        Numbered run name for this run, used to name output files and directories.

    run_path : str or os.PathLike, optional
        Base standardized run directory for outputs. Required for all modes.
        Deleted automatically in `simulation` mode.

    grouping_col : str, default='group'
        Column name used to group TCRs into plates.

    oligo_name_annotation_col_list : list of str or None
        Columns from the input CSV used to annotate oligo names. Defaults to ['custom_name'].

    filter_cannot_be_codon_optimized : bool, default=True
        Whether to remove TCR sequences that cannot be codon-optimized.

    filter_succeeding_nt_cys_104_beta : bool, default=False
        Whether to filter CDR3β sequences based on the expected second amino acid after
        the conserved Cys104.

        Details:
        - In Golden Gate assembly, BbsI generates four-base sticky overhangs that join TRBV
        and CDR3β-Jβ fragments.
        - The first three nucleotides of the overhang come from the conserved Cys104 codon (TGT).
        - The fourth nucleotide of the overhang is contributed by the final nucleotide of the
        TRBV codon immediately before the Cys104 codon.
        - In human TRBV genes, this nucleotide is highly conserved:
            * G for most TRBV genes
            * A for TRBV20-1 and TRBV29-1
        (based on sequence analysis of ~89,000 TCRs).
        - Consequently, the second amino acid in the CDR3β sequence is typically:
            * Alanine (A) when the nucleotide is G
            * Serine (S) when the nucleotide is A
        Sequences not matching these patterns may represent sequencing artifacts.

        Behavior:
        - If True: sequences where the second amino acid does not provide the expected nucleotide
        for the Cys104 overhang are removed.
        - If False: sequences are retained, and the second amino acid is replaced by the V gene
        compatible amino acid to produce the correct 4-base overhang.

    allow_cdr3j_nt_duplicates : bool, default=False
        Whether to allow duplicate CDR3-J nucleotide sequences. Set to False if you want to distinguish
        all TCRs in your library by only sequencing the CDR3-J of a single chain (alpha or beta).
        CDR3-J nucleotide sequence duplicates are made unique by picking unique codons with the codon
        optimizer (i.e., codon optimizer induced diversity).

    steepness : float, default = 4
        Controls the bias toward more frequently used codons during codon optimization.
        - Higher values (e.g., 4) strongly favor more optimal codons in Homo sapiens genes.
        - Lower values (≥1) make codon choice more uniform among synonymous codons.
        - Must be ≥ 1 in `make_codon_sampling_probability_distribution_with_prior`.
        - Can be gradually reduced during optimization to relax codon bias when constraints
        (GC content, CAI, restriction sites) are not met.

    well_plate_size : int, default=384
        Plate size.

    remove_longer_than_200_nt_oligos : bool, default=False
        Remove oligos longer than 200 nucleotides.
        Twist oligo pools have pricing tiers based on oligo length and pool size.

    v_gene_premix_dispenser : str, default='idot'
        Dispenser type for V gene premix. Options: 'idot', 'echo'.
        `echo` has been successfully used in wet-lab but is not actively
        maintained.

    v_gene_transfer_volume_nl : int, default=150
        Volume (nL) of V gene stock to transfer to each assembly well.

    max_v_gene_source_well_volume_nl : int, default=70000
        Maximum allowed volume per V gene source well (nL).

    v_gene_hamilton_ul : float, default=21.12
        Hamilton transfer volume for V gene stock (µL).

    water_hamilton_ul : float, default=58.88
        Hamilton water volume for diluting transferred V gene stock (µL).

    skip_barcoded_v_gene_stock_tube_vol_update : bool, default=False
        Skip updating V gene stock tube volumes for barcoded tubes.
        Set to `True` if you do not want to use the barcoded V gene stock volume
        tracking system of this package.

    min_v_gene_stock_tube_vol_ul : float, default=50.0
        Minimum volume to maintain in V gene stock tube (µL).

    generate_illumina_refs : bool, default=True
        Generate Illumina reference FASTA files. Reference files are only generated in
        `run_assembly` mode.

    epitope_barcode_refs : Union[None, pd.DataFrame, str, os.PathLike[str], None], default  = None
        CSV or FASTA file containing antigen/epitope-encoding-minigene + barcode nucleotide
        order sequences that need to be written to the reference .fa file Include if you want to write a
        combined TCR beta + epitope nucleotide sequence reference for e.g., PAIR-scan.

    epitope_name_col_name : str, default='name'
        Column name for epitope identifiers in epitope CSV.

    epitope_order_col_name : str, default='sequence'
        Column name for epitope nucleotide sequences in epitope CSV.

    trimmed_beta_model_tcr_dict : dict or None
        Optionally include trimmed beta chain nucleotides of common model TCRs in reference files
        that are stored in this dictionary.

    tcr_read_length : int, default=150
        Illumina sequencing read length.

    model_epitope_dict : dict or None
        Optional include epitope minigenes of common model epitopes in reference files that are
        stored in this dictionary.

    epitope_barcode_length : int, default=18
        Length of epitope barcode sequences.

    generate_nanopore_ref : bool, default=True
        Generate Nanopore reference FASTA files. Reference files are only generated in `run_assembly`
        mode.

    add_number_of_negatives_by_mutating_full_refs : int or None
        Add negative control sequences by mutating full TCR references.

    add_number_of_negatives_by_mutating_cdr3j_in_refs : int or None
        Add negative control sequences by mutating CDR3-J regions.

    add_number_of_negatives_by_mutating_v_gene_in_refs : int or None
        Add negative control sequences by mutating V genes.

    min_nt_diff_negative_ref_seqs : int or None
        Minimum nucleotide difference for negative reference sequences.

    max_nt_diff_negative_ref_seqs : int or None
        Maximum nucleotide difference for negative reference sequences.

    Returns
    -------
    str or os.PathLike
        Path to the run directory for 'oligo_order' and 'run_assembly' modes.
        In 'simulation' mode, returns a status message indicating completion.

    Notes
    -----
    - Integrates multiple subfunctions including:
    `translate_v_gene_input_to_imgt`, `filter_tcr_df_on_v_genes_in_stock`,
    `add_cdr3j_alpha_beta_order_columns`, `add_standardized_tcr_assembly_well_plate_layout`,
    `add_plate_ortho_primer_combinations_to_cdr3j_seqs`,
    `make_tcr_assembly_run_cdr3j_seq_order_sheets`,
    `make_idot_v_gene_premixing_dispense_csv`,
    `reconstruct_tcrs_assembly`, `pydna_plate_sub_pool_amp_pcr`,
    `pydna_cdr3j_ortho_primer_amp_pcr`, `pydna_golden_gate_tcr_assembly`,
    `check_translation_assembled_tcr_ligation_products`,
    `generate_assembly_nt_refs`, and `generate_assembly_nanopore_nt_refs`.
    - Validates input CSV for required columns and NaN values.
    - Automatically handles TCR plate assignment, oligo naming, V gene premix generation,
    and `pydna` simulation logging.
    - Generates sequencing reference files compatible with Illumina and Nanopore workflows.
    """
    required_col_list = ["TRAV_IMGT", "TRAJ_IMGT", "cdr3_alpha_aa", "TRBV_IMGT", "TRBJ_IMGT", "cdr3_beta_aa", "custom_name", grouping_col]

    if oligo_name_annotation_col_list is None:
        oligo_name_annotation_col_list = ["custom_name"]

    run_tcr_df = pd.read_csv(run_tcr_csv)

    missing_cols = [col for col in required_col_list if col not in run_tcr_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    run_tcr_df = run_tcr_df.loc[:, required_col_list]

    nan_columns = [col for col in run_tcr_df.columns if run_tcr_df[col].isna().any()]

    if nan_columns:
        error_msg = "The following required columns contain NaN values:\n" + "\n".join(f"- {col}: {run_tcr_df[col].isna().sum()} NaN(s)" for col in nan_columns)
        raise ValueError(error_msg)

    run_tcr_df = translate_v_gene_input_to_imgt(
        run_tcr_df=run_tcr_df,
        v_alpha_col="TRAV_IMGT",
        v_beta_col="TRBV_IMGT",
        vdj_translation_dict_json=os.path.join(
            tcr_toolbox_data_path,
            "tcr_toolbox_datasets",
            "tcr_reconstruction",
            "VDJ_gene_sequences",
            "after_benchmark",
            "functional_with_L-PART1+V-EXON_after_benchmark",
            "20230803_vdj_translation_dict_aa.json",
        ),
    )

    if run_mode in ("oligo_order", "simulation"):
        print("\nFiltering run_tcr_df for TRAV and TRBV genes available in stock tubes:")
        run_tcr_df = filter_tcr_df_on_v_genes_in_stock(
            run_path=run_path,
            numbered_run_name=numbered_run_name,
            tcr_df=run_tcr_df,
            v_alpha_col="TRAV_IMGT",
            v_beta_col="TRBV_IMGT",
            tcr_df_name="",  # since we only have one individual tcr_df we don't need a name here
        )

        print("\n\nGenerating CDR3-J alpha and beta order nucleotide sequences:")
        run_tcr_df = add_cdr3j_alpha_beta_order_columns(
            tcr_df_v_gene_filtered=run_tcr_df,
            filter_cannot_be_codon_optimized=filter_cannot_be_codon_optimized,
            filter_succeeding_nt_cys_104_beta=filter_succeeding_nt_cys_104_beta,
            allow_duplicates=allow_cdr3j_nt_duplicates,
            steepness=steepness,
        )

        print("\nAssigning TCRs to wells in plates...")
        plates_df_dict = add_standardized_tcr_assembly_well_plate_layout(run_tcr_df=run_tcr_df, grouping_col=grouping_col, well_plate_size=well_plate_size)

        print("\n\nAdding well-specific orthoprimer combinations to CDR3-J order nucleotide sequences:")
        plates_df_dict = add_plate_ortho_primer_combinations_to_cdr3j_seqs(plates_df_dict=plates_df_dict, remove_longer_than_200_nt=remove_longer_than_200_nt_oligos)

        print("\n\nMaking CDR3-J nucleotide oligo pool order sheets:")
        plates_df_dict = make_tcr_assembly_run_cdr3j_seq_order_sheets(
            plates_df_dict=plates_df_dict, run_path=run_path, numbered_run_name=numbered_run_name, name_annotation_cols=oligo_name_annotation_col_list
        )

    if run_mode == "run_assembly":
        plate_sheet_dir = os.path.join(run_path, "plate_sheets")
        plates_df_dict = collections.defaultdict(pd.DataFrame)
        for plate_xlsx in glob.glob(os.path.join(plate_sheet_dir, "*.xlsx")):
            filename = os.path.basename(plate_xlsx)
            stem = os.path.splitext(filename)[0]
            plate_num = int(stem.split("_")[-1])
            plates_df_dict[plate_num] = pd.read_excel(plate_xlsx, index_col=0)
            plates_df_dict[plate_num]["well_coord"] = plates_df_dict[plate_num]["well_coord"].apply(lambda x: tuple(ast.literal_eval(x)) if isinstance(x, str) else x)
            plates_df_dict[plate_num].dropna(subset=["cdr3j_alpha_nt_order_primers", "cdr3j_beta_nt_order_primers"], inplace=True)
            plates_df_dict[plate_num].reset_index(drop=True, inplace=True)

        tcr_refs_df = pd.concat([tcr_df for tcr_df in plates_df_dict.values()])
        tcr_refs_df.reset_index(drop=True, inplace=True)
        if tcr_refs_df.loc[:, "name"].duplicated().any():
            raise Exception("There are duplicate TCR names. Did you manually edit the .xlsx plate sheets?")
        if tcr_refs_df.loc[:, "name"].isna().any():
            raise Exception("There are missing TCR names. Did you manually edit the .xlsx plate sheets?")

        print("\n")
        for plate in plates_df_dict.keys():
            print("Plotting V gene plate map:", plate)
            plot_v_genes_per_well(
                plates_df_dict[plate],
                save_path=os.path.join(run_path, "v_genes_premix_dispense", "v_gene_plate_maps", numbered_run_name + "_plate_" + str(plate) + ".pdf"),
                plate_title="plate_" + str(plate),
                well_384_or_96=well_plate_size,
                TRAV_col_name="TRAV_IMGT_allele_collapsed",
                TRBV_col_name="TRBV_IMGT_allele_collapsed",
            )

        print("\n\nWriting v gene premix dispense instruction files:")
    elif run_mode in ("oligo_order", "simulation"):
        print("\n\nChecking whether V gene stock tubes have enough remaing volume for this assembly:")
    source_name = re.search(r"(r\d+)_", numbered_run_name).group(1) + "_source_1"
    if v_gene_premix_dispenser == "idot":
        run_tcr_df = make_idot_v_gene_premixing_dispense_csv(
            plates_df_dict=plates_df_dict,
            source_name=source_name,
            run_path=run_path,
            transfer_volume_nl=v_gene_transfer_volume_nl,
            numbered_run_name=numbered_run_name,
            max_idot_vol_nl=max_v_gene_source_well_volume_nl,
            v_gene_hamilton_ul=v_gene_hamilton_ul,
            water_hamilton_ul=water_hamilton_ul,
            min_v_gene_stock_tube_vol_ul=min_v_gene_stock_tube_vol_ul,
            run_mode=run_mode,
            skip_barcoded_v_gene_stock_tube_vol_update=skip_barcoded_v_gene_stock_tube_vol_update,
        )
    elif v_gene_premix_dispenser == "echo":
        raise NotImplementedError(
            "`make_echo_v_gene_premixing_dispense_csv` function needs to be updated to\ncorrectly update the volumes of V gene stock tubes in our V gene stock database."
        )
    else:
        raise NotImplementedError(f"v_gene_premix_dispenser {v_gene_premix_dispenser} has not been implemented yet.")

    print("\n\nTCR assembly preparation run finished!\nStarting in sillico pydna simulation ligation product check of assembly run.")
    print("\nReconstructing TCR amino acid sequence:")
    plates_df_dict = reconstruct_tcrs_assembly(
        include_leader=True,
        include_constant=True,
        mouse_or_human="mouse",
        constant_beta=pydna_tcr_golden_gate_constant_aa_seq_dict["muTRBC_aa"],
        constant_alpha=pydna_tcr_golden_gate_constant_aa_seq_dict["muTRAC_aa"],
        exclude_c_fw=False,
        plates_df_dict=plates_df_dict,
        verbose=False,
    )

    base_dir = os.path.join(run_path, "cdr3j_oligo_order_sheets")
    pattern = os.path.join(base_dir, f"{numbered_run_name}_run_sub_pool_*.xlsx")
    files = glob.glob(pattern)

    sub_pools = sorted({int(re.search(r"_run_sub_pool_(\d+)_", f).group(1)) for f in files})

    if not sub_pools:
        raise FileNotFoundError(f"No subpool files found matching pattern: {pattern}")

    for run_sub_pool in sub_pools:
        print("\nPerforming pydna simulation for run_sub_pool:", run_sub_pool)
        matching_files = glob.glob(os.path.join(base_dir, f"{numbered_run_name}_run_sub_pool_{run_sub_pool}_plates_*.xlsx"))

        if len(matching_files) != 1:
            raise Exception(f"Expected exactly one .xlsx match for subpool {run_sub_pool}, but found {len(matching_files)}.")

        run_sub_pool_file = matching_files[0]
        print(f"Processing subpool {run_sub_pool}: {os.path.basename(run_sub_pool_file)}")

        plate_alpha_pcr_products_dict, plate_beta_pcr_products_dict = pydna_plate_sub_pool_amp_pcr(run_sub_pool_file=run_sub_pool_file, run_sub_pool=run_sub_pool)

        cdr3_alpha_product_dict, cdr3_beta_product_dict = pydna_cdr3j_ortho_primer_amp_pcr(
            run_sub_pool_file=run_sub_pool_file, plate_alpha_pcr_products_dict=plate_alpha_pcr_products_dict, plate_beta_pcr_products_dict=plate_beta_pcr_products_dict
        )

        print("Simulating Golden Gate:")
        ligation_products_dict = pydna_golden_gate_tcr_assembly(
            v_genes_premix_dispense_run_file=os.path.join(run_path, "v_genes_premix_dispense", numbered_run_name + "_idot_dispense.csv"),
            v_genes_source_run_files_list=[os.path.join(run_path, "v_genes_premix_dispense", source_name + ".xlsx")],
            cdr3_alpha_product_dict=cdr3_alpha_product_dict,
            cdr3_beta_product_dict=cdr3_beta_product_dict,
            run_sub_pool_file=run_sub_pool_file,
            run_path=run_path,
            numbered_run_name=numbered_run_name,
            run_sub_pool=run_sub_pool,
            run_mode=run_mode,
        )

        print("Comparing amino acid translation of ligation product to reconstructed TCR amino acid sequence:")
        correct_translation_ligation_products_dict = check_translation_assembled_tcr_ligation_products(
            plates_df_dict=plates_df_dict,
            ligation_products_dict=ligation_products_dict,
            run_path=run_path,
            numbered_run_name=numbered_run_name,
            run_sub_pool=run_sub_pool,
            run_mode=run_mode,
        )

    print("pydna simulation finished.\nIf you simulated in run_mode = 'oligo_order' or 'run_assembly', you check pydna logs in the pydna_logs directory.\n")

    if run_mode == "oligo_order":
        shutil.rmtree(os.path.join(run_path, "v_genes_premix_dispense"))
        print("Finished oligo_order mode!")
        return run_path

    if run_mode == "simulation":
        if any(fname.endswith(".log") for fname in os.listdir(run_path)):
            print(
                f"WARNING: A pipeline log file was found in '{run_path}'.\nIn simulation run_mode, the pipeline log file should be written to a directory other than the standardized run directory ('{run_path}'), as this directory will be deleted after the simulation."
            )
        shutil.rmtree(run_path)
        return "Finished simulation mode! Simulation results can be found in the pipeline .log file."

    if generate_illumina_refs:
        print("\nGenerating Illumina reference .fa files:")
        if epitope_barcode_refs:
            ext = os.path.splitext(epitope_barcode_refs)[1].lower()
            if ext == ".csv":
                epitope_barcode_refs = pd.read_csv(epitope_barcode_refs)
            elif ext == ".fa" or ext == ".fasta":
                # Keep as-is (FASTA file path)
                pass
            else:
                raise FileNotFoundError(f"{epitope_barcode_refs} must be a CSV or FASTA file.")
        else:
            epitope_barcode_refs = None

        generate_assembly_nt_refs(
            tcr_refs_df=tcr_refs_df,
            epitope_barcode_refs=epitope_barcode_refs,  # None
            tcr_name_col_name="name",
            epitope_name_col_name=epitope_name_col_name,  # "name"
            fasta_alpha_out_fname=os.path.join(run_path, "sequencing_quality_analysis", "references", f"{numbered_run_name}_{tcr_read_length}bp_alpha.fa"),  # None
            fasta_beta_out_fname=os.path.join(run_path, "sequencing_quality_analysis", "references", f"{numbered_run_name}_{tcr_read_length}bp_beta.fa"),
            fasta_beta_epitope_plate_seq_out_fname=os.path.join(
                run_path, "sequencing_quality_analysis", "references", f"{numbered_run_name}_{tcr_read_length}bp_pairscan_beta_epi_plate.fa"
            ),
            fasta_epitope_out_fname=os.path.join(run_path, "sequencing_quality_analysis", "references", f"{numbered_run_name}_{tcr_read_length}bp_epi.fa"),
            trimmed_beta_model_tcr_dict=trimmed_beta_model_tcr_dict,  # None
            model_epitope_dict=model_epitope_dict,  # None
            alpha_order_col_name="cdr3j_alpha_nt_order_primers",
            beta_order_col_name="cdr3j_beta_nt_order_primers",
            trim_assembly_primers_from_cdr3j=True,
            v_alpha_col="TRAV_IMGT_allele_collapsed",
            v_beta_col="TRBV_IMGT_allele_collapsed",
            add_v_gene_to_duplicate_filter=True,
            epitope_order_col_name=epitope_order_col_name,  # "sequence"
            epitope_barcode_length=epitope_barcode_length,  # 18
            read_length=tcr_read_length,  # 150
            gtf=True,
            verbose=1,
        )

    if generate_nanopore_ref:
        print("\nGenerating nanopore reference .fa file:")
        generate_assembly_nanopore_nt_refs(
            tcr_refs_df=tcr_refs_df,
            tcr_name_col_name="name",
            fasta_out_fname=os.path.join(run_path, "sequencing_quality_analysis", "references", f"{numbered_run_name}_full_tcr_nanopore.fa"),
            alpha_order_col_name="cdr3j_alpha_nt_order_primers",
            beta_order_col_name="cdr3j_beta_nt_order_primers",
            trim_assembly_primers_from_cdr3j=True,
            include_mu_constant_beta_and_p2a=True,
            v_alpha_col="TRAV_IMGT_allele_collapsed",
            v_beta_col="TRBV_IMGT_allele_collapsed",
            add_number_of_negatives_by_mutating_full_refs=add_number_of_negatives_by_mutating_full_refs,  # None
            add_number_of_negatives_by_mutating_cdr3j_in_refs=add_number_of_negatives_by_mutating_cdr3j_in_refs,  # None
            add_number_of_negatives_by_mutating_v_gene_in_refs=add_number_of_negatives_by_mutating_v_gene_in_refs,  # None
            min_nt_diff_negative_ref_seqs=min_nt_diff_negative_ref_seqs,  # None
            max_nt_diff_negative_ref_seqs=max_nt_diff_negative_ref_seqs,  # None
            verbose=1,
        )
        tcr_refs_df.to_csv(os.path.join(run_path, "sequencing_quality_analysis", "references", "tcr_refs_df.csv"), index=False)

    print("Finished run_assembly mode!")
    return run_path
