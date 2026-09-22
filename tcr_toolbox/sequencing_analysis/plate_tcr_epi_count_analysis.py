import collections
import itertools
import json
import math
import os
import re
from pathlib import Path
from typing import Tuple, Union

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import scipy as sp
from adjustText import adjust_text
from dotenv import load_dotenv
from scipy.stats import binom, spearmanr

from tcr_toolbox.sequencing_analysis.utils import collapse_duplicate_tcr_names_reference_file, read_gdna_counts_csv, remove_dna_barcode_from_string
from tcr_toolbox.utils.plot_utils import startfig, view_color_dict
from tcr_toolbox.utils.stat_utils import fdrcorrection, normalized_entropy, bootstrap_spearman_ci
from tcr_toolbox.utils.utils import ensure_list

load_dotenv()
tcr_toolbox_data_path = os.getenv("tcr_toolbox_data_path")
if tcr_toolbox_data_path is None:
    raise EnvironmentError("The 'tcr_toolbox_data_path' environment variable is not set (checked .env and the process environment).")


def init_umi_count_analysis_dir(project_dir: str | os.PathLike):
    path = Path(project_dir)
    (path / "preprocessing").mkdir(parents=True, exist_ok=False)
    (path / "outs").mkdir(parents=True, exist_ok=False)
    (path / "outs" / "plate_maps").mkdir(parents=True, exist_ok=False)


def read_umi_count_tsv(
    project_dir: str | os.PathLike,
    plate_name_split_location: int = 2,
    barcode_xlsx_file: Union[str, os.PathLike[str]] = tcr_toolbox_data_path + "/tcr_toolbox_datasets/pair_scan_luna_plate_seq/barcodes/final_tcr_barcodes_03_02_2022-15-46.xlsx",
) -> dict:

    project_dir = Path(project_dir)

    barcodes_df = pd.read_excel(barcode_xlsx_file)
    barcodes_df["Barcode"] = barcodes_df["Name"].str.split("-").str[1]

    counts_df_dict = {}
    for count_tsv in (project_dir / "counts").glob("*_counts.tsv"):
        plate_name = count_tsv.stem.split("_")[plate_name_split_location]
        print("Reading plate:", plate_name)
        counts_df_dict[plate_name] = pd.read_csv(count_tsv, sep="\t")
        counts_df_dict[plate_name]["well"] = counts_df_dict[plate_name]["cell"].map(dict(zip(barcodes_df["Barcode"], barcodes_df["Well position"])))
    return counts_df_dict


def filter_UMI_counts_df_dict(
    project_dir: str | os.PathLike,
    counts_df_dict: dict,
    epi_umi_threshold: int | dict,
    tcr_umi_threshold: int | dict,
    min_start_fraction_of_top_umi_count: float | dict,
    min_end_fraction_of_top_umi_count: float | dict,
    num_bins: int = 50,
):
    """
    Well-specific UMI filtering. Filtering removes low-abundance transcripts and non-top transcripts based on
    their UMI count ratio relative to the top transcript per well. Filtering is performed separately for TCR and
    Ag transcripts. The minimum ratio required to keep a non-top transcript is not fixed: it follows a straight
    line in log-log space between (top transcript count = UMI threshold, ratio = min_start_fraction_of_top_umi_count)
    and (top transcript count = 99th percentile, ratio = min_end_fraction_of_top_umi_count), so wells with a higher
    top transcript count are held to a lower minimum ratio, since at a high top count even a smaller fraction still
    corresponds to a large absolute UMI count that can be called with confidence. For each plate, a diagnostic QC
    PDF of the filtering thresholds is written to project_dir/preprocessing.

    Parameters
    ----------
    project_dir : str | os.PathLike
        String path to plate-based TCR and epitope sequencing project directory. Diagnostic QC filtering plots are
        written to the "preprocessing" subdirectory of this path.
    counts_df_dict : dict
        Dictionary where keys represent plate names and values represent DataFrames that store 384-well plate UMI
        counts, as returned by :func:`read_umi_count_tsv`.
    epi_umi_threshold : int | dict
        Minimum UMI count required to keep an epitope transcript in a well. If a dict, keys must match the plate
        names in counts_df_dict and values are the per-plate threshold.
    tcr_umi_threshold : int | dict
        Minimum UMI count required to keep a TCR transcript in a well. If a dict, keys must match the plate names
        in counts_df_dict and values are the per-plate threshold.
    min_start_fraction_of_top_umi_count : float | dict
        Ratio at the low end (top transcript count = UMI threshold) of the diagonal threshold line. If a dict,
        keys must match the plate names in counts_df_dict.
    min_end_fraction_of_top_umi_count : float | dict
        Ratio at the high end (top transcript count = 99th percentile) of the diagonal threshold line. If a dict,
        keys must match the plate names in counts_df_dict.
    num_bins : int, default = 50
        Number of log-spaced bins used for the top UMI transcript count histograms in the diagnostic QC plot.

    Returns
    -------
    counts_df_dict : dict
        Dictionary where keys represent plate names and values represent DataFrames that store the filtered
        384-well plate UMI counts (only transcripts that passed the top/non-top filtering criteria above).
    """
    counts_df_dict = counts_df_dict.copy()
    out_dir = Path(project_dir) / "preprocessing"
    out_dir.mkdir(parents=True, exist_ok=True)

    def make_diagonal_threshold_function(x1, y1, x2, y2):
        """
        Create a power-law threshold function in log-log space.
        """

        logx1 = np.log10(x1)
        logx2 = np.log10(x2)

        logy1 = np.log10(y1)
        logy2 = np.log10(y2)

        if logx2 == logx1:
            raise ValueError(f"Cannot compute a diagonal threshold: x1 ({x1}) and x2 ({x2}) resolve to the same value in log space.")
        slope = (logy2 - logy1) / (logx2 - logx1)

        def threshold(x):

            x = np.asarray(x)

            return 10 ** (logy1 + slope * (np.log10(x) - logx1))

        return threshold

    for plate_name in counts_df_dict.keys():
        print("Filtering plate:", plate_name)

        if isinstance(epi_umi_threshold, dict):
            epi_umi_threshold_int = epi_umi_threshold[plate_name]
        else:
            epi_umi_threshold_int = epi_umi_threshold

        if isinstance(tcr_umi_threshold, dict):
            tcr_umi_threshold_int = tcr_umi_threshold[plate_name]
        else:
            tcr_umi_threshold_int = tcr_umi_threshold

        if isinstance(min_start_fraction_of_top_umi_count, dict):
            min_start_fraction_of_top_umi_count_float = min_start_fraction_of_top_umi_count[plate_name]
        else:
            min_start_fraction_of_top_umi_count_float = min_start_fraction_of_top_umi_count

        if isinstance(min_end_fraction_of_top_umi_count, dict):
            min_end_fraction_of_top_umi_count_float = min_end_fraction_of_top_umi_count[plate_name]
        else:
            min_end_fraction_of_top_umi_count_float = min_end_fraction_of_top_umi_count

        df = counts_df_dict[plate_name].copy()
        df["epi_or_tcr"] = df["gene"].str.split("_").str[0]
        df["top_count"] = df.groupby(["well", "epi_or_tcr"])["count"].transform("max")

        df = df.loc[df["top_count"] > 3, :].copy()
        df["ratio_to_top"] = df["count"] / df["top_count"]
        secondary_df = df.loc[df["ratio_to_top"] < 1, :].copy()

        depth_df = df.groupby(["well", "epi_or_tcr"])["top_count"].max().reset_index()
        epi_depth = depth_df.loc[depth_df["epi_or_tcr"] == "epi", "top_count"]
        tcr_depth = depth_df.loc[depth_df["epi_or_tcr"] == "tcr", "top_count"]

        xmin = 1
        epi_xmax = max(epi_depth.max(), 1) * 1.25
        tcr_xmax = max(tcr_depth.max(), 1) * 1.25

        epi_x99 = np.percentile(epi_depth, 99)
        tcr_x99 = np.percentile(tcr_depth, 99)

        epi_bins = np.logspace(np.log10(xmin), np.log10(epi_xmax), num_bins)
        tcr_bins = np.logspace(np.log10(xmin), np.log10(tcr_xmax), num_bins)

        epi_threshold_fn = make_diagonal_threshold_function(
            x1=epi_umi_threshold_int, y1=min_start_fraction_of_top_umi_count_float, x2=epi_x99, y2=min_end_fraction_of_top_umi_count_float
        )
        tcr_threshold_fn = make_diagonal_threshold_function(
            x1=tcr_umi_threshold_int, y1=min_start_fraction_of_top_umi_count_float, x2=tcr_x99, y2=min_end_fraction_of_top_umi_count_float
        )

        fig = plt.figure(figsize=(6.5, 4.5), constrained_layout=True)
        gs = fig.add_gridspec(2, 2, height_ratios=[1, 4])
        ax_hist_epi = fig.add_subplot(gs[0, 0])
        ax_hist_tcr = fig.add_subplot(gs[0, 1])
        ax_epi = fig.add_subplot(gs[1, 0])
        ax_tcr = fig.add_subplot(gs[1, 1], sharey=ax_epi)

        ax_hist_epi.hist(epi_depth, bins=epi_bins, color="lightgrey")
        ax_hist_tcr.hist(tcr_depth, bins=tcr_bins, color="lightgrey")

        for ax, thr, title, xmax in [(ax_hist_epi, epi_umi_threshold_int, "epi", epi_xmax), (ax_hist_tcr, tcr_umi_threshold_int, "tcr", tcr_xmax)]:
            ax.set_xscale("log")
            ax.set_xlim(xmin, xmax)
            ax.axvline(thr, color="red", linestyle="--")
            ax.set_title(title, fontsize=7)
            ax.set_ylabel("Number of wells", fontsize=7)
            ax.tick_params("both", labelsize=7)

        plot_configs = [(ax_epi, "epi", epi_umi_threshold_int, epi_threshold_fn, epi_xmax), (ax_tcr, "tcr", tcr_umi_threshold_int, tcr_threshold_fn, tcr_xmax)]

        scatters = []
        for ax, group, thr, threshold_fn, xmax in plot_configs:
            sub = secondary_df.loc[secondary_df["epi_or_tcr"] == group, :].copy()

            sc = ax.scatter(sub["top_count"], sub["ratio_to_top"], c=sub["count"], cmap="viridis", s=5, alpha=0.3, norm=mpl.colors.LogNorm(), edgecolors="none")
            scatters.append(sc)

            ax.axvline(thr, color="red", linestyle="--")

            x_line = np.logspace(np.log10(thr), np.log10(xmax), 500)
            y_line = threshold_fn(x_line)
            ax.plot(x_line, y_line, color="red", linestyle="--", linewidth=1)

            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(xmin, xmax)
            ax.set_xlabel("Top UMI transcript count per well", fontsize=7)
            ax.set_ylabel("Fraction non-top of top UMI count", fontsize=7)
            ax.tick_params("both", labelsize=7)

        cbar_epi = fig.colorbar(scatters[0], ax=ax_epi, location="bottom", pad=0.1, aspect=40)
        cbar_tcr = fig.colorbar(scatters[1], ax=ax_tcr, location="bottom", pad=0.1, aspect=40)
        cbar_epi.solids.set_alpha(1)
        cbar_tcr.solids.set_alpha(1)

        for cbar in [cbar_epi, cbar_tcr]:
            cbar.set_label("non-top UMI count", fontsize=7)
            cbar.ax.tick_params(labelsize=7)

        fig.savefig(out_dir / (f"{plate_name}_top_UMI_count_hist_and_fraction_non_top_of_top_scatter.pdf"))
        plt.close(fig)

        epi_df = df.loc[df["epi_or_tcr"] == "epi", :].copy()
        tcr_df = df.loc[df["epi_or_tcr"] == "tcr", :].copy()

        epi_primary_df = epi_df[epi_df["ratio_to_top"] == 1.0]
        epi_secondary_df = epi_df[epi_df["ratio_to_top"] < 1.0]

        tcr_primary_df = tcr_df[tcr_df["ratio_to_top"] == 1.0]
        tcr_secondary_df = tcr_df[tcr_df["ratio_to_top"] < 1.0]

        epi_primary_kept_df = epi_primary_df[epi_primary_df["count"] >= epi_umi_threshold_int].copy()

        epi_secondary_kept_df = epi_secondary_df[
            (epi_secondary_df["count"] >= epi_umi_threshold_int) & (epi_secondary_df["ratio_to_top"] >= epi_threshold_fn(epi_secondary_df["top_count"]))
        ].copy()

        tcr_primary_kept_df = tcr_primary_df[tcr_primary_df["count"] >= tcr_umi_threshold_int].copy()

        tcr_secondary_kept_df = tcr_secondary_df[
            (tcr_secondary_df["count"] >= tcr_umi_threshold_int) & (tcr_secondary_df["ratio_to_top"] >= tcr_threshold_fn(tcr_secondary_df["top_count"]))
        ].copy()

        epi_filter_df = pd.concat([epi_primary_kept_df, epi_secondary_kept_df]).reset_index(drop=True)
        tcr_filter_df = pd.concat([tcr_primary_kept_df, tcr_secondary_kept_df]).reset_index(drop=True)

        filtered_df = pd.concat([epi_filter_df, tcr_filter_df]).reset_index(drop=True)

        counts_df_dict[plate_name] = filtered_df

        print("Number of unique epis:", filtered_df.loc[filtered_df["epi_or_tcr"] == "epi", "gene"].nunique())
        print("Number of unique tcrs:", filtered_df.loc[filtered_df["epi_or_tcr"] == "tcr", "gene"].nunique())

    return counts_df_dict


def validate_validating_pair_dict(validating_pair_dict: dict) -> None:
    if len(validating_pair_dict) == 0:
        raise ValueError("validating_pair_dict cannot be an empty dictionary")
    if not all(isinstance(v, bool) for v in validating_pair_dict.values()):
        raise ValueError("validating_pair_dict values must all be boolean (True or False)")
    if not all(pair.split("-")[0].startswith("epi_") for pair in validating_pair_dict):
        raise ValueError('First element of the "-" split in validating_pair_dict keys does not start with "epi_". Did you encode pairs as epitope-TCR?')
    return None


def plot_detected_tcr_epi_plate_map(
    project_dir: str | os.PathLike,
    counts_df_dict: dict,
    validating_pair_dict: dict = None,
    epitope_barcode_length: int = 18,
    detected_fontsize: int = 4,
    barcode_xlsx_file: str | os.PathLike = tcr_toolbox_data_path + "/tcr_toolbox_datasets/pair_scan_luna_plate_seq/barcodes/final_tcr_barcodes_03_02_2022-15-46.xlsx",
):
    """
    Write a .pdf 384-well plate map plot to outs/plate_maps showing detected TCR and epitope transcripts per well.

    For every well, all detected transcript reference names (with the DNA barcode stripped from epitope names) are
    printed as text inside that well's cell of the plate map. If validating_pair_dict is provided, any well in
    which both members of a known-reactive (True) TCR-epitope pair are detected is additionally highlighted in
    light green.


    Parameters
    ----------
    project_dir : str | os.PathLike
        String path to plate-based TCR and epitope sequencing project directory. Plate map PDFs are written to the
        "outs/plate_maps" subdirectory of this path.
    counts_df_dict : dict
        UMI count filtered dictionary where keys represent plate names and values represent DataFrames that store
        filtered 384-well plate UMI counts, as returned by :func:`filter_UMI_counts_df_dict`.
    validating_pair_dict : dict, optional
        Dictionary where keys are epitope-TCR pair strings ("<epitope_reference_name>-<tcr_reference_name>") and
        values are booleans indicating whether that pair was confirmed reactive (True) in an independent 1-to-1
        validation assay. Must pass :func:`validate_validating_pair_dict`. If provided, wells in which a pair with
        value True is co-detected are highlighted in the plate map.
    epitope_barcode_length : int, default = 18
        Epitope barcode length. Used to remove the DNA barcode from epitope reference name strings before plotting.
    detected_fontsize : int, default = 4
        Plotting fontsize of detected transcripts per well. Adjust when transcript reference names are long and do
        not fit in a well.
    barcode_xlsx_file : str | os.PathLike, default = tcr_toolbox_data_path + '/tcr_toolbox_datasets/pair_scan_luna_plate_seq/barcodes/final_tcr_barcodes_03_02_2022-15-46.xlsx'
        DNA Barcode to 384-well plate well name map .xlsx file.

    Returns
    -------
    None
        Writes one "{plate_name}_plate_map.pdf" file per plate to project_dir/outs/plate_maps.
    """
    barcodes_df = pd.read_excel(barcode_xlsx_file)
    barcodes_df["Barcode"] = barcodes_df["Name"].str.split("-").str[1]

    letter_to_number_map = dict(zip([letter for letter in "ABCDEFGHIJKLMNOP"], [number for number in np.arange(len("ABCDEFGHIJKLMNOP"))]))
    for plate_name in counts_df_dict.keys():
        print("Plotting plate:", plate_name)

        count_df = counts_df_dict[plate_name].copy()
        array_plate_df = pd.DataFrame([], index=letter_to_number_map.keys(), columns=np.arange(24))
        array_plate = np.zeros((16, 24))

        for well in barcodes_df["Well position"]:
            well_gene_detect_list = []
            i = int(letter_to_number_map[re.split(r"\d", well)[0]])
            j = int(re.split("[A-P]", well)[1]) - 1

            if count_df[count_df["well"] == well].shape[0] != 0:
                well_gene_detect_list = count_df[count_df["well"] == well]["gene"].unique().tolist()
                well_gene_detect_list = [
                    remove_dna_barcode_from_string(string=gene, barcode_length=epitope_barcode_length) if gene.startswith("epi_") else gene for gene in well_gene_detect_list
                ]

            if well_gene_detect_list:
                array_plate_df.iloc[i, j] = well_gene_detect_list

            if validating_pair_dict is not None:
                validate_validating_pair_dict(validating_pair_dict)

                for pair, is_validating in validating_pair_dict.items():
                    if is_validating and (pair.split("-")[0] in well_gene_detect_list) and (pair.split("-")[1] in well_gene_detect_list):
                        array_plate[i, j] = 1
                        break

        from matplotlib.colors import ListedColormap

        cmap = ListedColormap(["white", "lightgreen"])
        ax, fig, gs = startfig(w=32, h=55)
        ax.matshow(array_plate, cmap=cmap)
        ax.set_xticks(np.arange(array_plate.shape[1]), labels=[str(int(column) + 1) for column in np.arange(array_plate.shape[1])], fontsize=8)
        ax.set_yticks(np.arange(array_plate.shape[0]), labels=letter_to_number_map.keys(), fontsize=8)
        ax.xaxis.set_ticks_position("top")
        ax.yaxis.set_ticks_position("left")

        for i in range(array_plate_df.shape[0]):
            for j in range(array_plate_df.shape[1]):
                if isinstance(array_plate_df.iloc[i, j], list):
                    plot_str_tmp = "".join([str(gene) + "\n" for gene in array_plate_df.iloc[i, j]])[:-1]
                    text = ax.text(j, i, plot_str_tmp, ha="center", va="center", color="black", fontsize=detected_fontsize)
                else:
                    continue

        ax.set_title(plate_name, pad=25)
        fig.tight_layout()
        fig.savefig(Path(project_dir) / "outs" / "plate_maps" / f"{plate_name}_plate_map.pdf")
        plt.close()


def write_tcr_epi_pair_well_count_files(
    project_dir: str | os.PathLike,
    plates_list: list,
    plates_list_name: str,
    counts_df_dict: dict,
    collapse_epitope_barcode: bool,
    collapse_custom_ref_names_dict: dict = None,
    collapse_tcr_technical_duplicates: bool = False,
    reference_file: bool = None,
    epitope_barcode_length: int = 18,
    plates_rows_for_pair_counting_dict: dict = None,
    plates_rows_for_normalization_dict: dict = None,
    validating_pair_dict: dict = None,
):
    """
    Count detected TCR-epitope pairs per well across a set of plates and write TCR-epitope pair count files to
    outs/<plates_list_name>.

    If plates_rows_for_pair_counting_dict and plates_rows_for_normalization_dict are provided, a subset of rows
    on each plate can instead be reserved for computing per-transcript detection counts used for normalization
    (e.g., wells where T cells and B cells were sorted separately rather than as co-culture doublets), while the
    remaining rows are used for pair counting.

    Files written to outs/<plates_list_name>:
        - pair_count_matrix_df_<plates_list_name>.xlsx/.csv: Epitope-by-TCR matrix of well counts.
        - pair_counts_sorted_<plates_list_name>.xlsx/.csv: detected pairs sorted by well count.
        - tcr_co_occurrence_counts_<plates_list_name>.xlsx/.csv (and _total_well_normalized variant): TCR-TCR
          well co-occurrence counts (and well-count-normalized fractions).
        - epitopes_co_occurrence_counts_<plates_list_name>.xlsx/.csv (and _total_well_normalized variant):
          Epitope-Epitope well co-occurrence counts (and well-count-normalized fractions).
        - transcript_tuple_well_counter_<plates_list_name>.json: well counts per unique combination of detected
          transcripts (used by :func:`add_tcr_epi_pairs_only_supported_by_co_occurrence_bool_col`).
        - summary_statistics_<plates_list_name>.txt: run- and plate-level summary statistics (see below).
        - norm_tcr_detect_count_<plates_list_name>.csv / norm_epi_detect_count_<plates_list_name>.csv: per-transcript
          detection counts in the normalization rows, if plates_rows_for_normalization_dict is provided.

    Files written per plate to outs/plate_maps:
        - {plate_name}_tcrs_per_well.xlsx / {plate_name}_epitopes_per_well.xlsx: unique detected TCRs/epitopes
          per well.
        - {plate_name}_wells_per_pair.xlsx: wells in which each detected pair was found.
        - {plate_name}_pairs_per_well.xlsx: detected pairs per well.
        - {plate_name}_reactive_pair_in_well_bool.xlsx: per-well boolean of whether a known-reactive pair from
          validating_pair_dict was detected, if validating_pair_dict is provided.

    Parameters
    ----------
    project_dir : str | os.PathLike
        String path to plate-based TCR and epitope sequencing project directory.
    plates_list : list
        List of plate names to be counted. Each name must be a key in counts_df_dict.
    plates_list_name : str
        Name of the list of plates to be counted, used to name the outs/<plates_list_name> output directory and
        output files. For example, 'run-4-1' for run-4-1 plates.
    counts_df_dict : dict
        UMI count filtered dictionary where keys represent plate names and values represent DataFrames that store
        filtered 384-well plate UMI counts, as returned by :func:`filter_UMI_counts_df_dict`.
    collapse_epitope_barcode : bool
        If True, remove the DNA barcode from epitope reference name strings before counting. Used to collapse
        distinct DNA barcode variants of the same epitope construct into a single epitope reference name.
    collapse_custom_ref_names_dict : dict, optional
        Dictionary mapping specific transcript reference names to a replacement reference name, applied to both
        TCR and epitope transcripts after barcode/duplicate collapsing. Use this to merge transcript reference
        names that should be treated as identical but are not already merged by collapse_epitope_barcode or
        collapse_tcr_technical_duplicates.
    collapse_tcr_technical_duplicates : bool, default = False
        If True, collapse TCR reference names that are technical duplicates of one another (e.g., re-synthesized
        or re-cloned versions of the same TCR) using the mapping generated by
        :func:`tcr_toolbox.sequencing_analysis.utils.collapse_duplicate_tcr_names_reference_file` from
        reference_file.
    reference_file : bool, optional
        Path to the reference file passed to
        :func:`tcr_toolbox.sequencing_analysis.utils.collapse_duplicate_tcr_names_reference_file` to build the TCR
        technical duplicate collapsing map. Required when collapse_tcr_technical_duplicates is True.
    epitope_barcode_length : int, default = 18
        Epitope barcode length. Used to remove the DNA barcode from epitope reference name strings when
        collapse_epitope_barcode is True.
    plates_rows_for_pair_counting_dict : dict, optional
        Dictionary with keys representing plates and associated values representing the letters of the rows that
        should be used for TCR-epitope pair counting. Use this argument if, for example, into only a subset of
        rows T cell-B cell conjugates were sorted from a co-culture. Must be provided together with
        plates_rows_for_normalization_dict, and the two row sets for a given plate must not overlap.
    plates_rows_for_normalization_dict : dict, optional
        Dictionary with keys representing plates and associated values representing the letters of the rows that
        should be used for per-transcript detection count normalization. Use this argument if, for example, into a
        subset of rows Jurkat T cells and B cells were separately sorted (i.e., no co-culture selection was
        performed). Must be provided together with plates_rows_for_pair_counting_dict.
    validating_pair_dict : dict, optional
        Dictionary where keys are epitope-TCR pair strings ("<epitope_reference_name>-<tcr_reference_name>") and
        values are booleans indicating whether that pair was validated to be reactive (True) or non-reactive (False)
        in an independent 1-to-1 validation assay. Must pass :func:`validate_validating_pair_dict`. If provided,
        a per-well boolean of whether a True pair was detected is written to {plate_name}_reactive_pair_in_well_bool.xlsx
        and summarized in the run summary statistics.

    Returns
    -------
    None
        Writes the pair count, co-occurrence, and summary statistics files described above.

    Examples
    --------
    Use a subset of rows for TCR-epitope pair counting and a subset of rows for normalization (e.g., rows in which
    Jurkat T cell and B cell singlets were separately sorted):

    >>> plates_list = ['3-1-1', '3-1-2', '3-1-3', '3-1-4']
    >>> plates_rows_for_pair_counting_dict = {
    ...     '3-1-1': ['A', 'B', 'C', 'D'],
    ...     '3-1-2': ['A', 'B', 'C', 'D', 'E', 'F'],
    ...     '3-1-3': ['A', 'B', 'C', 'D', 'E'],
    ...     '3-1-4': ['A', 'B', 'C', 'D', 'E', 'F']
    ... }
    >>> plates_rows_for_normalization_dict = {
    ...      '3-1-1': [],
    ...      '3-1-2': ['G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P'],
    ...      '3-1-3': ['G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P'],
    ...      '3-1-4': ['G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
    ...  }
    >>> write_tcr_epi_pair_well_count_files(project_dir=project_dir,
    ...                                     plates_list=plates_list,
    ...                                     plates_list_name='run-3-1',
    ...                                     counts_df_dict=counts_df_dict,
    ...                                     collapse_epitope_barcode=True,
    ...                                     epitope_barcode_length=18,
    ...                                     plates_rows_for_pair_counting_dict=plates_rows_for_pair_counting_dict,
    ...                                     plates_rows_for_normalization_dict=plates_rows_for_normalization_dict,
    ...                                     )

    Use all rows for TCR-epitope pair counting:

    >>> plates_list = ['4-2-1', '4-2-2', '4-2-3']
    >>> write_tcr_epi_pair_well_count_files(project_dir=project_dir,
    ...                                     plates_list=plates_list,
    ...                                     plates_list_name='run-4-2',
    ...                                     counts_df_dict=counts_df_dict,
    ...                                     collapse_epitope_barcode=True,
    ...                                     epitope_barcode_length=18,
    ...                                     )
    """
    detected_pair_dict = collections.defaultdict(int)
    well_detect_pair_dict = collections.defaultdict(lambda: collections.defaultdict(list))
    pair_detect_well_dict = collections.defaultdict(lambda: collections.defaultdict(list))
    tcr_detect_well_dict = collections.defaultdict(lambda: collections.defaultdict(list))
    epi_detect_well_dict = collections.defaultdict(lambda: collections.defaultdict(list))
    transcript_tuple_well_counter = collections.defaultdict(int)
    number_of_tcr_per_well_dict = collections.defaultdict(list)
    number_of_epi_per_well_dict = collections.defaultdict(list)
    tcr_co_occur_dict = collections.defaultdict(list)
    tcr_co_occur_plate_dict = collections.defaultdict(lambda: collections.defaultdict(list))
    epi_co_occur_dict = collections.defaultdict(list)
    epi_co_occur_plate_dict = collections.defaultdict(lambda: collections.defaultdict(list))
    plate_well_counter_dict = collections.defaultdict(int)
    well_counter = 0
    well_with_detected_pair_counter = 0
    well_with_detected_tcr_counter = 0
    well_with_detected_epi_counter = 0

    if plates_rows_for_pair_counting_dict and not plates_rows_for_normalization_dict:
        raise Exception("When providing plates_rows_for_pair_counting_dict, plates_rows_for_normalization_dict must also be provided!")

    if plates_rows_for_normalization_dict and not plates_rows_for_pair_counting_dict:
        raise Exception("When providing plates_rows_for_normalization_dict, plates_rows_for_pair_counting_dict must also be provided!")

    if plates_rows_for_normalization_dict:
        norm_tcr_detect_count_dict = collections.defaultdict(int)
        norm_epi_detect_count_dict = collections.defaultdict(int)

    if validating_pair_dict is not None:
        validate_validating_pair_dict(validating_pair_dict)

        detected_prior_reactive_pair_well_dict = collections.defaultdict(lambda: collections.defaultdict(bool))

    if collapse_tcr_technical_duplicates:
        duplicate_tcr_collapse_dict = collapse_duplicate_tcr_names_reference_file(reference_file=reference_file)

    plates_list_outs = Path(project_dir) / "outs" / plates_list_name
    plates_list_outs.mkdir(exist_ok=False)
    plate_map_outs = Path(project_dir) / "outs" / "plate_maps"
    summary_stat_write_str = "## " + plates_list_name + "\n\n"

    # This reorder_target_by_well_list code block needs to become a small util function...
    # but that can only happen once order_automation_pydna_mm_2 is merged into main:
    letters = "ABCDEFGHIJKLMNOP"
    col_numbers = np.arange(1, 25)
    reorder_target_by_well_list = []
    for letter in letters:
        for col_number in col_numbers:
            reorder_target_by_well_list.append(letter + str(col_number))
    sorter_index = dict(zip(reorder_target_by_well_list, range(len(reorder_target_by_well_list))))

    for plate_name in plates_list:
        print("Processing plate:", plate_name)
        count_df = counts_df_dict[plate_name].copy()

        if plates_rows_for_pair_counting_dict and plates_rows_for_normalization_dict:
            if (plate_name not in plates_rows_for_pair_counting_dict.keys()) or (plate_name not in plates_rows_for_normalization_dict.keys()):
                raise Exception(
                    plate_name, "does not have a key in both plates_rows_for_pair_counting_dict and plates_rows_for_normalization_dict while this is required when these are used!"
                )

            if not plates_rows_for_pair_counting_dict[plate_name] and not plates_rows_for_normalization_dict[plate_name]:
                raise Exception(
                    "Either plates_rows_for_pair_counting_dict:",
                    plates_rows_for_pair_counting_dict[plate_name],
                    "or plates_rows_for_normalization_dict:",
                    plates_rows_for_normalization_dict[plate_name],
                    "together do not store at least one row for plate:",
                    plate_name,
                    "If a plate is listed in the provided plates_list and a plates_rows_for_pair_counting_dict and plates_rows_for_normalization_dict is provided as argument, these dicts together need to store at least one row from the plate that should be used for file writing!",
                )

            if set(plates_rows_for_pair_counting_dict[plate_name]).intersection(set(plates_rows_for_normalization_dict[plate_name])):
                raise Exception(
                    "There are overlapping rows in plates_rows_for_pair_counting_dict and plates_rows_for_normalization_dict for plate:",
                    plate_name,
                    "while this is not allowed. Use a row from the same plate either for TCR-epitope pairs counting or for normalization!",
                )

            if plates_rows_for_normalization_dict[plate_name]:
                norm_count_df = count_df[count_df["well"].str.startswith(tuple(plates_rows_for_normalization_dict[plate_name]))].copy()
                if norm_count_df.empty:
                    raise Exception(
                        "norm_count_df is empty while this should not be possible when plates_rows_for_normalization_dict contains rows for plate:",
                        plate_name,
                        "Did you use capital letters to store the row names in the list?",
                    )

                for well in norm_count_df["well"].unique():
                    epi_transcripts = norm_count_df[norm_count_df["well"] == well]["gene"][norm_count_df[norm_count_df["well"] == well]["gene"].str.startswith("epi_")]
                    tcr_transcripts = norm_count_df[norm_count_df["well"] == well]["gene"][norm_count_df[norm_count_df["well"] == well]["gene"].str.startswith("tcr_")]

                    if collapse_epitope_barcode:
                        epi_transcripts = [remove_dna_barcode_from_string(string=transcript, barcode_length=epitope_barcode_length) for transcript in epi_transcripts]
                    if collapse_tcr_technical_duplicates:
                        tcr_transcripts = [duplicate_tcr_collapse_dict[transcript] for transcript in tcr_transcripts]

                    if collapse_custom_ref_names_dict:
                        epi_transcripts = [
                            collapse_custom_ref_names_dict[transcript] if transcript in collapse_custom_ref_names_dict.keys() else transcript for transcript in epi_transcripts
                        ]
                        tcr_transcripts = [
                            collapse_custom_ref_names_dict[transcript] if transcript in collapse_custom_ref_names_dict.keys() else transcript for transcript in tcr_transcripts
                        ]

                    for tcr_transcript in tcr_transcripts:
                        norm_tcr_detect_count_dict[tcr_transcript] += 1

                    for epi_transcript in epi_transcripts:
                        norm_epi_detect_count_dict[epi_transcript] += 1

                norm_tcr_detect_count_df = pd.DataFrame(pd.Series(norm_tcr_detect_count_dict), columns=["count"])
                norm_tcr_detect_count_df = norm_tcr_detect_count_df.sort_values("count", ascending=False)
                norm_tcr_detect_count_df.to_csv(plates_list_outs / f"norm_tcr_detect_count_{plates_list_name}.csv", index_label="reference_name")

                norm_epi_detect_count_df = pd.DataFrame(pd.Series(norm_epi_detect_count_dict), columns=["count"])
                norm_epi_detect_count_df = norm_epi_detect_count_df.sort_values("count", ascending=False)
                norm_epi_detect_count_df.to_csv(plates_list_outs / f"norm_epi_detect_count_{plates_list_name}.csv", index_label="reference_name")

            if plates_rows_for_pair_counting_dict[plate_name]:
                # Subset to rows that should be used for TCR-epitope pairs counting:
                count_df = count_df[count_df["well"].str.startswith(tuple(plates_rows_for_pair_counting_dict[plate_name]))].copy()
                if count_df.empty:
                    raise Exception(
                        "count_df is empty but that should not be possible when plates_rows_for_pair_counting_dict contains rows for plate:",
                        plate_name,
                        "Did you use capital letters to store the row names in the list?",
                    )
            else:
                # If users provide a plates_rows_for_pair_counting_dict and this dict does not store any rows for this plate_name
                # (e.g., when the plate is used for normalization only), we do not want to write any TCR-epitope pairs
                # count files for this plate using the code below and continue to the next plate by calling continue:
                continue

        for well in count_df["well"].unique():
            well_counter += 1
            plate_well_counter_dict[plate_name] += 1

            epi_transcripts = count_df[count_df["well"] == well]["gene"][count_df[count_df["well"] == well]["gene"].str.startswith("epi_")]
            tcr_transcripts = count_df[count_df["well"] == well]["gene"][count_df[count_df["well"] == well]["gene"].str.startswith("tcr_")]

            if collapse_epitope_barcode:
                epi_transcripts = [remove_dna_barcode_from_string(string=transcript, barcode_length=epitope_barcode_length) for transcript in epi_transcripts]

            if collapse_tcr_technical_duplicates:
                tcr_transcripts = [duplicate_tcr_collapse_dict[transcript] for transcript in tcr_transcripts]

            if collapse_custom_ref_names_dict:
                epi_transcripts = [
                    collapse_custom_ref_names_dict[transcript] if transcript in collapse_custom_ref_names_dict.keys() else transcript for transcript in epi_transcripts
                ]
                tcr_transcripts = [
                    collapse_custom_ref_names_dict[transcript] if transcript in collapse_custom_ref_names_dict.keys() else transcript for transcript in tcr_transcripts
                ]

            for tcr_transcript in tcr_transcripts:
                tcr_detect_well_dict[plate_name][well].append(tcr_transcript)
            for epi_transcript in epi_transcripts:
                epi_detect_well_dict[plate_name][well].append(epi_transcript)

            if tcr_detect_well_dict[plate_name][well] and epi_detect_well_dict[plate_name][well]:
                well_with_detected_pair_counter += 1
            elif tcr_detect_well_dict[plate_name][well] and not epi_detect_well_dict[plate_name][well]:
                well_with_detected_tcr_counter += 1
            elif not tcr_detect_well_dict[plate_name][well] and epi_detect_well_dict[plate_name][well]:
                well_with_detected_epi_counter += 1

            transcript_tuple_well_counter[tuple(list(np.unique(tcr_transcripts)) + list(np.unique(epi_transcripts)))] += 1

            if len(np.unique(tcr_transcripts)) > 1:
                tcr_co_occur_dict["".join([tcr_transcript + "-" for tcr_transcript in np.unique(tcr_transcripts)])[:-1]].append(1)
                tcr_co_occur_plate_dict[plate_name]["".join([tcr_transcript + "-" for tcr_transcript in np.unique(tcr_transcripts)])[:-1]].append(1)
            elif len(np.unique(tcr_transcripts)) == 1:
                tcr_co_occur_dict[np.unique(tcr_transcripts)[0]].append(0)
                tcr_co_occur_plate_dict[plate_name][np.unique(tcr_transcripts)[0]].append(0)

            if len(np.unique(epi_transcripts)) > 1:
                epi_co_occur_dict["".join([epi_transcript + "-" for epi_transcript in np.unique(epi_transcripts)])[:-1]].append(1)
                epi_co_occur_plate_dict[plate_name]["".join([epi_transcript + "-" for epi_transcript in np.unique(epi_transcripts)])[:-1]].append(1)
            elif len(np.unique(epi_transcripts)) == 1:
                epi_co_occur_dict[np.unique(epi_transcripts)[0]].append(0)
                epi_co_occur_plate_dict[plate_name][np.unique(epi_transcripts)[0]].append(0)

            tcr_detect_well_dict[plate_name][well] = np.unique(tcr_detect_well_dict[plate_name][well])
            epi_detect_well_dict[plate_name][well] = np.unique(epi_detect_well_dict[plate_name][well])
            number_of_tcr_per_well_dict[plate_name].append(len(tcr_detect_well_dict[plate_name][well]))
            number_of_epi_per_well_dict[plate_name].append(len(epi_detect_well_dict[plate_name][well]))

            if validating_pair_dict is not None:
                detected_prior_reactive_pair_well_dict[plate_name][well] = False

            for epi_transcript, tcr_transcript in itertools.product(epi_transcripts, tcr_transcripts):
                detected_pair_dict[epi_transcript + "-" + tcr_transcript] += 1
                well_detect_pair_dict[plate_name][epi_transcript + "-" + tcr_transcript].append(well)
                pair_detect_well_dict[plate_name][well].append(epi_transcript + "-" + tcr_transcript)

                if validating_pair_dict is not None and validating_pair_dict.get(epi_transcript + "-" + tcr_transcript) is True:
                    detected_prior_reactive_pair_well_dict[plate_name][well] = True

        tcr_detect_well_df = pd.DataFrame([(well, tcr) for well, tcr in tcr_detect_well_dict[plate_name].items()], columns=["well", "tcr"])
        tcr_detect_well_df["well_rank"] = tcr_detect_well_df["well"].map(sorter_index)
        tcr_detect_well_df.sort_values(["well_rank"], inplace=True)
        tcr_detect_well_df.drop("well_rank", axis=1, inplace=True)
        tcr_detect_well_df.reset_index(inplace=True, drop=True)
        tcr_detect_well_df.to_excel(plate_map_outs / f"{plate_name}_tcrs_per_well.xlsx")

        epi_detect_well_df = pd.DataFrame([(well, epi) for well, epi in epi_detect_well_dict[plate_name].items()], columns=["well", "epi"])
        epi_detect_well_df["well_rank"] = epi_detect_well_df["well"].map(sorter_index)
        epi_detect_well_df.sort_values(["well_rank"], inplace=True)
        epi_detect_well_df.drop("well_rank", axis=1, inplace=True)
        epi_detect_well_df.reset_index(inplace=True, drop=True)
        epi_detect_well_df.to_excel(plate_map_outs / f"{plate_name}_epitopes_per_well.xlsx")

        pd.DataFrame([(pair, wells) for pair, wells in well_detect_pair_dict[plate_name].items()], columns=["pair", "wells"]).to_excel(
            plate_map_outs / f"{plate_name}_wells_per_pair.xlsx"
        )

        pair_detect_well_df = pd.DataFrame([(well, pairs) for well, pairs in pair_detect_well_dict[plate_name].items()], columns=["well", "pairs"])
        pair_detect_well_df["well_rank"] = pair_detect_well_df["well"].map(sorter_index)
        pair_detect_well_df.sort_values(["well_rank"], inplace=True)
        pair_detect_well_df.drop("well_rank", axis=1, inplace=True)
        pair_detect_well_df.reset_index(inplace=True, drop=True)
        pair_detect_well_df.to_excel(plate_map_outs / f"{plate_name}_pairs_per_well.xlsx")

    pd.Series(detected_pair_dict, index=detected_pair_dict.keys()).sort_values(ascending=False).to_frame().rename({0: "well_count"}, axis=1).to_excel(
        plates_list_outs / f"pair_counts_sorted_{plates_list_name}.xlsx"
    )
    pd.Series(detected_pair_dict, index=detected_pair_dict.keys()).sort_values(ascending=False).to_frame().rename({0: "well_count"}, axis=1).to_csv(
        plates_list_outs / f"pair_counts_sorted_{plates_list_name}.csv"
    )

    print("Writing", plates_list_name, "summary statistics...")
    tcr_averages_list = []
    tcr_medians_list = []
    epi_averages_list = []
    epi_medians_list = []
    for plate_name in plates_list:
        summary_stat_write_str += "# " + plate_name + "\n"
        summary_stat_write_str += "Total number of co-culture sort wells counted for the following plate statistics: " + str(plate_well_counter_dict[plate_name]) + "\n"
        summary_stat_write_str += "Average # unique TCRs/well: " + str(round(np.average([number for number in number_of_tcr_per_well_dict[plate_name]]), 2)) + "\n"
        tcr_averages_list.append(np.average([number for number in number_of_tcr_per_well_dict[plate_name]]))

        summary_stat_write_str += "Average # unique epitopes/well: " + str(round(np.average([number for number in number_of_epi_per_well_dict[plate_name]]), 2)) + "\n"
        epi_averages_list.append(np.average([number for number in number_of_epi_per_well_dict[plate_name]]))

        summary_stat_write_str += "Median # unique TCRs/well: " + str(np.median([number for number in number_of_tcr_per_well_dict[plate_name]])) + "\n"
        tcr_medians_list.append(np.median([number for number in number_of_tcr_per_well_dict[plate_name]]))

        summary_stat_write_str += "Median # unique epitopes/well: " + str(np.median([number for number in number_of_epi_per_well_dict[plate_name]])) + "\n"
        epi_medians_list.append(np.median([number for number in number_of_epi_per_well_dict[plate_name]]))

        tcr_co_occur_plate_df = pd.DataFrame(pd.Series({tcr: np.sum(tcr_list) for tcr, tcr_list in tcr_co_occur_plate_dict[plate_name].items()}), columns=["count"])
        summary_stat_write_str += "TCR co-occurrence count average: " + str(round(tcr_co_occur_plate_df["count"].mean(), 6)) + "\n"
        epi_co_occur_plate_df = pd.DataFrame(pd.Series({epi: np.sum(epi_list) for epi, epi_list in epi_co_occur_plate_dict[plate_name].items()}), columns=["count"])
        summary_stat_write_str += "Epitopes co-occurrence count average: " + str(round(epi_co_occur_plate_df["count"].mean(), 6)) + "\n"

        if validating_pair_dict is not None:
            if plates_rows_for_pair_counting_dict and not plates_rows_for_pair_counting_dict[plate_name]:
                summary_stat_write_str += "\n"
                continue

            else:
                detected_prior_reactive_pair_well_df = pd.DataFrame(
                    [(well, reactive_bool) for well, reactive_bool in detected_prior_reactive_pair_well_dict[plate_name].items()], columns=["well", "reactive_bool"]
                )
                if not detected_prior_reactive_pair_well_df.shape[0] == plate_well_counter_dict[plate_name]:
                    raise Exception("Number of counted wells in detected_prior_reactive_pair_well_df is not equal to number of counted wells in plate_well_counter_dict!")

                summary_stat_write_str += (
                    "Fraction of sorted wells with reactive TCR-epitope pairs: "
                    + str(round(detected_prior_reactive_pair_well_df["reactive_bool"].sum() / detected_prior_reactive_pair_well_df.shape[0], 2))
                    + "\n"
                )
                detected_prior_reactive_pair_well_df["well_rank"] = detected_prior_reactive_pair_well_df["well"].map(sorter_index)
                detected_prior_reactive_pair_well_df.sort_values(["well_rank"], inplace=True)
                detected_prior_reactive_pair_well_df.drop("well_rank", axis=1, inplace=True)
                detected_prior_reactive_pair_well_df.reset_index(inplace=True, drop=True)
                detected_prior_reactive_pair_well_df.to_excel(plate_map_outs / f"{plate_name}_reactive_pair_in_well_bool.xlsx")
        summary_stat_write_str += "\n"

    summary_stat_write_str += "\n"
    print("Total # of co-culture sort wells counted:", well_counter)
    summary_stat_write_str += "Total # of co-culture sort wells counted for following aggregate run statistics: " + str(well_counter) + "\n"
    summary_stat_write_str += "Total # of co-culture sort wells with at least one detected TCR-Epitope pair: " + str(well_with_detected_pair_counter) + "\n"
    summary_stat_write_str += (
        "Fraction total co-culture sort wells with at least one detected TCR-Epitope pair of total counted sort wells: "
        + str(round(well_with_detected_pair_counter / well_counter, 2))
        + "\n"
    )
    summary_stat_write_str += "Total # of co-culture sort wells with at least one detected TCR but not an Epitope: " + str(well_with_detected_tcr_counter) + "\n"
    summary_stat_write_str += "Total # of co-culture sort wells with at least one detected Epitope but not a TCR: " + str(well_with_detected_epi_counter) + "\n"
    summary_stat_write_str += "Mean average # of TCRs/well/plate: " + str(round(np.average([number for number in tcr_averages_list if not math.isnan(number)]), 2)) + "\n"
    summary_stat_write_str += "Mean average # of epitopes/well/plate: " + str(round(np.average([number for number in epi_averages_list if not math.isnan(number)]), 2)) + "\n"
    summary_stat_write_str += "Mean median # of TCRs/well/plate: " + str(round(np.average([number for number in tcr_medians_list if not math.isnan(number)]), 2)) + "\n"
    summary_stat_write_str += "Mean median # of epitopes/well/plate: " + str(round(np.average([number for number in epi_medians_list if not math.isnan(number)]), 2)) + "\n"

    tcr_co_occur_df = pd.DataFrame(pd.Series({tcr: np.sum(tcr_list) for tcr, tcr_list in tcr_co_occur_dict.items()}), columns=["count"])
    tcr_co_occur_df = tcr_co_occur_df.sort_values("count", ascending=False)
    tcr_co_occur_df.to_excel(plates_list_outs / f"tcr_co_occurrence_counts_{plates_list_name}.xlsx", index_label="co_occurring_tcr")
    tcr_co_occur_df.to_csv(plates_list_outs / f"tcr_co_occurrence_counts_{plates_list_name}.csv", index_label="co_occurring_tcr")
    summary_stat_write_str += (
        "Average TCR co-occurrence count per well: " + str(round(np.average([count for tcr_list in tcr_co_occur_dict.values() for count in tcr_list]), 6)) + "\n"
    )

    summary_stat_write_str += "TCR co-occurrence count average: " + str(round(tcr_co_occur_df["count"].mean(), 6)) + "\n"
    total_number_of_detected_tcr_wells = np.sum(len(number_tcr_list) for number_tcr_list in number_of_tcr_per_well_dict.values())
    if not total_number_of_detected_tcr_wells == well_counter:
        raise Exception("Total number of counted wells in number_of_tcr_per_well_dict is not equal to number of counted wells in well_counter!")
    tcr_co_occur_df["count"] = tcr_co_occur_df["count"].div(well_counter)
    tcr_co_occur_df.to_excel(plates_list_outs / f"tcr_co_occurrence_counts_total_well_normalized_{plates_list_name}.xlsx", index_label="co_occurring_tcr")
    tcr_co_occur_df.to_csv(plates_list_outs / f"tcr_co_occurrence_counts_total_well_normalized_{plates_list_name}.csv", index_label="co_occurring_tcr")

    summary_stat_write_str += "Average TCR co-occurrence fraction of all wells: " + str(round(tcr_co_occur_df["count"].mean(), 6)) + "\n"

    epi_co_occur_df = pd.DataFrame(pd.Series({epi: np.sum(epi_list) for epi, epi_list in epi_co_occur_dict.items()}), columns=["count"])
    epi_co_occur_df = epi_co_occur_df.sort_values("count", ascending=False)
    epi_co_occur_df.to_excel(plates_list_outs / f"epitopes_co_occurrence_counts_{plates_list_name}.xlsx", index_label="co_occurring_epitopes")
    epi_co_occur_df.to_csv(plates_list_outs / f"epitopes_co_occurrence_counts_{plates_list_name}.csv", index_label="co_occurring_epitopes")
    summary_stat_write_str += (
        "Average epitopes co-occurrence count per well: " + str(round(np.average([count for epi_list in epi_co_occur_dict.values() for count in epi_list]), 6)) + "\n"
    )

    summary_stat_write_str += "Epitope co-occurrence count average: " + str(round(epi_co_occur_df["count"].mean(), 6)) + "\n"
    total_number_of_detected_epi_wells = np.sum(len(number_epi_list) for number_epi_list in number_of_epi_per_well_dict.values())
    if not total_number_of_detected_epi_wells == well_counter:
        raise Exception("Total number of counted wells in number_of_epi_per_well_dict is not equal to number of counted wells in well_counter!")
    epi_co_occur_df["count"] = epi_co_occur_df["count"].div(well_counter)
    epi_co_occur_df.to_excel(plates_list_outs / f"epitopes_co_occurrence_counts_total_well_normalized_{plates_list_name}.xlsx", index_label="co_occurring_epitopes")
    epi_co_occur_df.to_csv(plates_list_outs / f"epitopes_co_occurrence_counts_total_well_normalized_{plates_list_name}.csv", index_label="co_occurring_epitopes")

    summary_stat_write_str += "Average epitope co-occurrence fraction of all wells: " + str(round(epi_co_occur_df["count"].mean(), 6)) + "\n"

    # tuple keys are not supported in json format:
    transcript_tuple_well_counter = {str(tuple_key): count for tuple_key, count in transcript_tuple_well_counter.items()}
    with open(plates_list_outs / f"transcript_tuple_well_counter_{plates_list_name}.json", "w") as json_file:
        json.dump(transcript_tuple_well_counter, json_file)

    pair_count_tuples_list = [(pair.split("-")[0], pair.split("-")[1], well_count) for pair, well_count in detected_pair_dict.items()]
    pair_count_matrix_df = pd.DataFrame(pair_count_tuples_list, columns=["Epitope", "TCR", "Count"])
    pair_count_matrix_df = pair_count_matrix_df.pivot(index="Epitope", columns="TCR", values="Count").fillna(0).copy()
    pair_count_matrix_df.to_excel(plates_list_outs / f"pair_count_matrix_df_{plates_list_name}.xlsx")
    pair_count_matrix_df.to_csv(plates_list_outs / f"pair_count_matrix_df_{plates_list_name}.csv")

    with open(plates_list_outs / f"summary_statistics_{plates_list_name}.txt", "w") as summary_out_txt:
        summary_out_txt.write(summary_stat_write_str)


def read_tcr_epi_pair_counts_csv(pair_count_matrix_csv: str | os.PathLike) -> pd.DataFrame:
    """Read in pair_count_matrix csv file that was generated by
    :func:`tcr_toolbox.sequencing_analysis.write_tcr_epi_pair_well_count_files`

    Parameters
    ----------
    pair_count_matrix_csv : str | os.PathLike
        String path to pair_count_matrix csv file. pair_count_matrix_df needs to be
        written to a .csv file by :func:`tcr_toolbox.sequencing_analysis.plate_tcr_epi_count_analysis.write_tcr_epi_pair_well_count_files`.

    Returns
    -------
    pair_count_matrix_df
        DataFrame with TCR-epitope pair well counts.
    """
    pair_count_matrix_df = pd.read_csv(pair_count_matrix_csv, index_col=0)

    if not pair_count_matrix_df.index[0].startswith("epi_"):
        raise Exception(
            "pair_count_matrix_df does not store Epitope reference names in the index!", "Did you generate the pair_count_matrix_csv file with write_tcr_epi_pair_well_count_files?"
        )

    if not pair_count_matrix_df.columns[0].startswith("tcr_"):
        raise Exception(
            "pair_count_matrix_df does not store TCR reference names in the index!", "Did you generate the pair_count_matrix_csv file with write_tcr_epi_pair_well_count_files?"
        )

    return pair_count_matrix_df


def calculate_tcr_epi_pair_exp_prob_matrix(
    pair_count_matrix_df: pd.DataFrame,
    normalized_epitope_bulk_counts_df: pd.DataFrame = None,
    normalized_tcr_bulk_counts_df: pd.DataFrame = None,
    bulk_epitope_count_col: str = "",
    bulk_tcr_count_col: str = "",
    normalize_using_plate: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate expected TCR-epitope pair frequencies in case of random pairing.

    Expected pair frequency is calculated as each Epitope's normalized baseline bulk library frequency multiplied
    by each TCR's normalized baseline bulk library frequency

    Parameters
    ----------
    pair_count_matrix_df : pd.DataFrame
        Epitope-by-TCR matrix that contains detected TCR-Epitope pair well counts. pair_count_matrix_df needs to be
        written to a .csv file by :func:`write_tcr_epi_pair_well_count_files`
        and read by :func:`read_tcr_epi_pair_counts_csv`.
    normalized_epitope_bulk_counts_df : pd.DataFrame, optional
        DataFrame that contains epitope-barcode bulk gDNA-sequencing read counts, total-read-count normalized by
        :func:`tcr_toolbox.utils.stat_utils.normalize_counts_df_by_total_counts`. Required (and only used) when
        normalize_using_plate is False.
    normalized_tcr_bulk_counts_df : pd.DataFrame, optional
        DataFrame that contains TCR beta bulk gDNA-sequencing read counts, total-read-count normalized by
        :func:`tcr_toolbox.utils.stat_utils.normalize_counts_df_by_total_counts`. Required (and only used) when
        normalize_using_plate is False.
    bulk_epitope_count_col : str
        Column of normalized_epitope_bulk_counts_df that contains the total-read-count-normalized epitope-barcode read
        counts. Required (and only used) when normalize_using_plate is False.
    bulk_tcr_count_col : str
        Column of normalized_tcr_bulk_counts_df that contains the total-read-count-normalized TCR beta read counts.
        Required (and only used) when normalize_using_plate is False.
    normalize_using_plate : bool, default = False
        If True, do not use bulk TCR and bulk Epitope read counts and normalize using total TCR and Epitope well
        counts in sorted plates (pair_count_matrix_df row/column sums) instead. Use this argument when bulk data is
        not available.

    Returns
    -------
    exp_pair_prob_matrix_df : pd.DataFrame
        Epitope-by-TCR matrix of expected pair frequencies (Epitope bulk frequency x TCR bulk frequency).
    pair_count_matrix_df : pd.DataFrame
        The input pair_count_matrix_df. When normalize_using_plate is False, this is subset to only the Epitopes
        and TCRs shared with normalized_epitope_bulk_counts_df/normalized_tcr_bulk_counts_df, so its index/columns
        match exp_pair_prob_matrix_df exactly; when normalize_using_plate is True, it is returned unchanged.
    """

    pair_count_matrix_df = pair_count_matrix_df.copy()

    if not pair_count_matrix_df.index[0].startswith("epi_"):
        raise Exception(
            "pair_count_matrix_df does not store Epitope reference names in the index!", "Did you generate the pair_count_matrix_csv file with write_tcr_epi_pair_well_count_files?"
        )

    if not pair_count_matrix_df.columns[0].startswith("tcr_"):
        raise Exception(
            "pair_count_matrix_df does not store TCR reference names in the index!", "Did you generate the pair_count_matrix_csv file with write_tcr_epi_pair_well_count_files?"
        )

    total_called_pairs = pair_count_matrix_df.sum().sum()
    print("Total called TCR-epitope pair counts:", total_called_pairs)

    # Normalize using plate sort data:
    if normalize_using_plate:
        print("Calculate expected TCR-epitope pair probability matrix using pair_count_matrix_df plate sequencing data...")
        norm_row_sums_series = pair_count_matrix_df.sum(1) / total_called_pairs
        norm_col_sums_series = pair_count_matrix_df.sum(0) / total_called_pairs

    # Normalize using bulk sorts:
    else:
        print("Calculating expected TCR-epitope pair probability matrix using TCR and epitope gDNA bulk sequencing data...")
        normalized_epitope_bulk_counts_df = normalized_epitope_bulk_counts_df.copy()
        normalized_tcr_bulk_counts_df = normalized_tcr_bulk_counts_df.copy()
        if not normalized_epitope_bulk_counts_df.index[0].startswith("epi_"):
            print('Adding "epi_" prefix to normalized_epitope_bulk_counts_df reference names!')
            normalized_epitope_bulk_counts_df.index = ["epi_" + ref_name for ref_name in normalized_epitope_bulk_counts_df.index]

        if not normalized_tcr_bulk_counts_df.index[0].startswith("tcr_"):
            print('Adding "tcr_" prefix to normalized_epitope_bulk_counts_df reference names!')
            normalized_tcr_bulk_counts_df.index = ["tcr_" + ref_name for ref_name in normalized_tcr_bulk_counts_df.index]

        if len(set(normalized_epitope_bulk_counts_df.index).difference(set(pair_count_matrix_df.index))) > 0:
            print(
                "Epitope reference name(s) in normalized_epitope_bulk_counts_df that are missing in pair_count_matrix_df:",
                set(normalized_epitope_bulk_counts_df.index).difference(set(pair_count_matrix_df.index)),
            )
            print(
                "Epitope reference name(s) in pair_count_matrix_df that are missing in normalized_epitope_bulk_counts_df:",
                set(pair_count_matrix_df.index).difference(set(normalized_epitope_bulk_counts_df.index)),
            )

        if len(set(normalized_tcr_bulk_counts_df.index).difference(set(pair_count_matrix_df.columns))) > 0:
            print(
                "TCR reference name(s) in normalized_tcr_bulk_counts_df that are missing in pair_count_matrix_df:",
                set(normalized_tcr_bulk_counts_df.index).difference(set(pair_count_matrix_df.columns)),
            )
            print(
                "TCR reference name(s) in pair_count_matrix_df that are missing in normalized_tcr_bulk_counts_df:",
                set(pair_count_matrix_df.columns).difference(set(normalized_tcr_bulk_counts_df.index)),
            )

        intersect_epitopes_list = list(set(pair_count_matrix_df.index).intersection(normalized_epitope_bulk_counts_df.index))

        intersect_tcr_list = list(set(pair_count_matrix_df.columns).intersection(normalized_tcr_bulk_counts_df.index))

        normalized_epitope_bulk_counts_df = normalized_epitope_bulk_counts_df.loc[intersect_epitopes_list, :]
        normalized_tcr_bulk_counts_df = normalized_tcr_bulk_counts_df.loc[intersect_tcr_list, :]

        pair_count_matrix_df = pair_count_matrix_df.loc[intersect_epitopes_list, intersect_tcr_list]

        norm_row_sums_series = normalized_epitope_bulk_counts_df[bulk_epitope_count_col]
        norm_col_sums_series = normalized_tcr_bulk_counts_df[bulk_tcr_count_col]

    exp_pair_prob_matrix_arr = np.outer(norm_row_sums_series, norm_col_sums_series)
    exp_pair_prob_matrix_df = pd.DataFrame(exp_pair_prob_matrix_arr, index=norm_row_sums_series.index, columns=norm_col_sums_series.index)
    exp_pair_prob_matrix_df.index.name = "Epitope"
    exp_pair_prob_matrix_df.columns.name = "TCR"

    if exp_pair_prob_matrix_df.max().max() > 1:
        raise Exception(
            "Maximal value of exp_pair_prob_matrix_df is larger than 1!",
            "Did you normalize both normalized_epitope_bulk_counts_df and normalized_tcr_bulk_counts_df using tcr_toolbox.utils.stat_utils.normalize_counts_df_by_total_counts?",
        )

    return exp_pair_prob_matrix_df, pair_count_matrix_df


def obs_tcr_epi_count_binom_sf_test(
    project_dir: str | os.PathLike,
    pair_count_matrix_df: pd.DataFrame,
    exp_pair_prob_matrix_df: pd.DataFrame,
    plates_list_name: str,
    validating_pair_dict: dict = None,
    fdr_correction: bool = True,
    alpha: float = 0.05,
    normalize_total_called_pairs: bool = True,
    normalize_total_wells_with_detected_pair: bool = False,
    sort_fold_change: bool = False,
    sort_pvalue: bool = True,
) -> pd.DataFrame:
    """
    Estimate the probability of observing k sorted pair counts in n sorted wells with an expected pair frequency p as:
    Enrichment P value: P(X >= k observed pairs) for X ~ Binomial(n, p)
    Computed as: P(X≥k ) = binom.sf(k−1, n,p), where p is the pair's expected probability (from exp_pair_prob_matrix_df).
    Result is written to a returned pair_fold_changes_df that also reports the fold enrichment as
    observed pair frequency / expected pair frequency.

    Parameters
    ----------
    project_dir : str | os.PathLike
        String path to plate-based TCR and epitope sequencing project directory. Used to read the
        summary_statistics_<plates_list_name>.txt file written by
        :func:`write_tcr_epi_pair_well_count_files` when normalize_total_wells_with_detected_pair is True.
    pair_count_matrix_df : pd.DataFrame
        Epitope-by-TCR matrix that contains detected TCR-Epitope pair well counts. pair_count_matrix_df needs to be
        written to a .csv file by :func:`write_tcr_epi_pair_well_count_files`
        and read by :func:`read_tcr_epi_pair_counts_csv`.
    exp_pair_prob_matrix_df : pd.DataFrame
        Epitope-by-TCR matrix that contains expected TCR-epitope pair frequencies. Needs to be generated by
        :func:`calculate_tcr_epi_pair_exp_prob_matrix`.
    plates_list_name : str
        Name of the plates_list used when writing the pair count files with
        :func:`write_tcr_epi_pair_well_count_files` (e.g., 'run-4-1'). Only used to locate the
        summary_statistics_<plates_list_name>.txt file when normalize_total_wells_with_detected_pair is True.
    validating_pair_dict : dict, optional
        Dictionary where keys are epitope-TCR pair strings ("<epitope_reference_name>-<tcr_reference_name>") and
        values are booleans indicating whether that pair was confirmed reactive (True) or tested and not confirmed
        (False) in an independent 1-to-1 validation assay. Must pass :func:`validate_validating_pair_dict`. If
        provided, a "validating_pair" column is added to pair_fold_changes_df (True/False/None, where None means
        the pair was not tested in 1-to-1 co-culture).
    fdr_correction : bool, default = True
        If True, use Benjamini/Hochberg multiple hypothesis correction to adjust P values for controlling the
        false discovery rate (FDR) in independent tests.
    alpha : float, default = 0.05
        Family-wise error rate used for Benjamini/Hochberg multiple hypothesis correction.
    normalize_total_called_pairs : bool, default = True
        If True, use the total number of called TCR-epitope pairs (sum of pair_count_matrix_df) as the number of
        binomial trials N. Mutually exclusive with normalize_total_wells_with_detected_pair.
    normalize_total_wells_with_detected_pair : bool, default = False
        If True, use the total number of co-culture sort wells with at least one detected TCR-epitope pair (read
        from the summary_statistics_<plates_list_name>.txt file written by
        :func:`write_tcr_epi_pair_well_count_files`) as the number of binomial trials N. Mutually exclusive with
        normalize_total_called_pairs.
    sort_fold_change : bool, default = False
        If True, sort rows of pair_fold_changes_df by Observed pair fraction/Expected pair fraction fold changes.
    sort_pvalue : bool, default = True
        If True, sort rows of pair_fold_changes_df by their calculated P values. If fdr_correction = True, then rows
        will be sorted by their adjusted P values.

    Returns
    -------
    pair_fold_changes_df
        DataFrame that can be used to identify significantly enriched TCR-epitope pairs:
            - Epitope: Epitope from TCR-epitope pair.
            - TCR: TCR from TCR-epitope pair.
            - Fold_change: Observed TCR-epitope pair frequency/Expected TCR-epitope pair frequency fold change.
            - Pair: TCR-epitope pair string.
            - Pearson_residual: Pearson residual of the observed vs. expected well counts.
            - validating_pair: True/False/None status from validating_pair_dict, or None if pair was not measured and stored
              in validating_pair_dict.
            - BH_Reject: True if a hypothesis is rejected based on alpha after Benjamini/Hochberg multiple
              correction, False if not. Only present if fdr_correction = True.
            - P_value_corrected: P values adjusted for multiple hypothesis testing to limit false discovery rate.
              Only present if fdr_correction = True.
            - P_value: Probability (P value) that the probability of k observed TCR-epitope well counts in n sorted
              wells (i.e., observed TCR-epitope pair frequency) is greater than the expected TCR-epitope pair
              frequency p.
            - Exp_pair_prob: expected TCR-epitope pair fraction calculated using
              :func:`calculate_tcr_epi_pair_exp_prob_matrix`.
        Only pairs with a non-zero observed well count (Fold_change > 0) are retained.
    """
    pair_count_matrix_df = pair_count_matrix_df.copy()
    if not pair_count_matrix_df.index.name == "Epitope":
        pair_count_matrix_df.index.name = "Epitope"
    if not pair_count_matrix_df.columns.name == "TCR":
        pair_count_matrix_df.columns.name = "TCR"

    exp_pair_prob_matrix_df = exp_pair_prob_matrix_df.copy()

    if not all(pair_count_matrix_df.index == exp_pair_prob_matrix_df.index):
        raise Exception("pair_count_matrix_df epitope index is not equal to exp_pair_prob_matrix_df epitope index!")

    if not all(pair_count_matrix_df.columns == exp_pair_prob_matrix_df.columns):
        raise Exception("pair_count_matrix_df TCR column is not equal to exp_pair_prob_matrix_df TCR column!")

    if normalize_total_called_pairs and normalize_total_wells_with_detected_pair:
        raise ValueError("normalize_total_called_pairs and normalize_total_wells_with_detected_pair are mutually exclusive - set only one to True.")
    if not normalize_total_called_pairs and not normalize_total_wells_with_detected_pair:
        raise ValueError("Either normalize_total_called_pairs or normalize_total_wells_with_detected_pair must be True.")

    if normalize_total_called_pairs:
        total_called_pairs = pair_count_matrix_df.sum().sum()
        print("Total called # TCR-epitope pairs:", total_called_pairs)
        normalized_pair_count_matrix_df = pair_count_matrix_df / total_called_pairs
        N = total_called_pairs
    elif normalize_total_wells_with_detected_pair:
        with open(Path(project_dir) / "outs" / plates_list_name / f"summary_statistics_{plates_list_name}.txt", "r") as text_in:
            for line in text_in:
                if line.startswith("Total # of co-culture sort wells with at least one detected TCR-Epitope pair:"):
                    total_wells_with_detected_pair = int(line.split(":")[1])
        normalized_pair_count_matrix_df = pair_count_matrix_df / total_wells_with_detected_pair
        N = total_wells_with_detected_pair

    E = N * exp_pair_prob_matrix_df.to_numpy()
    pearson_residual_arr = (pair_count_matrix_df.to_numpy() - E) / np.sqrt(E)

    pearson_residual_df = pd.DataFrame(pearson_residual_arr, index=pair_count_matrix_df.index, columns=pair_count_matrix_df.columns)
    pearson_residual_df.index.name = "Epitope"
    pearson_residual_df.columns.name = "TCR"

    pair_fold_changes_df = normalized_pair_count_matrix_df / exp_pair_prob_matrix_df
    pair_fold_changes_df.reset_index(inplace=True)
    pair_fold_changes_df = pd.melt(pair_fold_changes_df, id_vars=["Epitope"], value_name="Fold_change")
    pair_fold_changes_df["Pair"] = pair_fold_changes_df["Epitope"] + "-" + pair_fold_changes_df["TCR"]

    pearson_residual_df_reset = pearson_residual_df.reset_index().melt(id_vars=["Epitope"], var_name="TCR", value_name="Pearson_residual")
    pearson_residual_df_reset["Pair"] = pearson_residual_df_reset["Epitope"] + "-" + pearson_residual_df_reset["TCR"]
    pearson_residual_df_reset.drop(columns=["Epitope", "TCR"], inplace=True)
    pair_fold_changes_df = pair_fold_changes_df.merge(pearson_residual_df_reset, on="Pair")

    if validating_pair_dict is not None:
        validate_validating_pair_dict(validating_pair_dict)

        validating_pairs = []
        for pair in pair_fold_changes_df["Pair"]:
            if pair not in validating_pair_dict:
                validating_pairs.append(None)  # None encodes pair not tested in 1-to-1 co-culture
            elif validating_pair_dict[pair]:
                validating_pairs.append(True)
            else:
                validating_pairs.append(False)
        pair_fold_changes_df["validating_pair"] = validating_pairs
    else:
        pair_fold_changes_df["validating_pair"] = None

    # Binom test obs count greater than exp is equal to binom.sf(k-1, n, p):
    # Note that we have to use the non-total_called_pairs normalized original pair_count_matrix_df:
    if normalize_total_called_pairs:
        pvals_arr = binom.sf(pair_count_matrix_df - np.ones_like(pair_count_matrix_df), total_called_pairs, exp_pair_prob_matrix_df)

    if normalize_total_wells_with_detected_pair:
        pvals_arr = binom.sf(pair_count_matrix_df - np.ones_like(pair_count_matrix_df), total_wells_with_detected_pair, exp_pair_prob_matrix_df)
    pvals_df = pd.DataFrame(pvals_arr, index=pair_count_matrix_df.index, columns=pair_count_matrix_df.columns)

    if fdr_correction:
        # NaN/zero counts in pair_count_matrix_df need to be filled to P values = 1, otherwise all pvals_corrected
        # will be NaN: https://github.com/statsmodels/statsmodels/issues/2899
        # If you do not want NaN/zero counts to be included in multiple hypothesis counts, a mask should be used as
        # in above statmodels github link.
        pvals = pvals_df.fillna(1).to_numpy().flatten()
        reject, pvals_corrected = fdrcorrection(pvals, alpha=alpha, method="indep", is_sorted=False)
        reject_df = pd.DataFrame(reject.reshape(pvals_df.shape), index=pvals_df.index, columns=pvals_df.columns)
        pvals_corrected_df = pd.DataFrame(pvals_corrected.reshape(pvals_df.shape), index=pvals_df.index, columns=pvals_df.columns)

        reject_df.reset_index(inplace=True)
        reject_df = pd.melt(reject_df, id_vars=["Epitope"], value_name="BH_Reject")
        reject_df["Pair"] = reject_df["Epitope"] + "-" + reject_df["TCR"]

        pvals_corrected_df.reset_index(inplace=True)
        pvals_corrected_df = pd.melt(pvals_corrected_df, id_vars=["Epitope"], value_name="P_value_corrected")
        pvals_corrected_df["Pair"] = pvals_corrected_df["Epitope"] + "-" + pvals_corrected_df["TCR"]

        pair_fold_changes_df = pair_fold_changes_df.merge(reject_df, on="Pair", suffixes=("", "_left"))
        pair_fold_changes_df.drop(["Epitope_left", "TCR_left"], axis=1, inplace=True)

        pair_fold_changes_df = pair_fold_changes_df.merge(pvals_corrected_df, on="Pair", suffixes=("", "_left"))
        pair_fold_changes_df.drop(["Epitope_left", "TCR_left"], axis=1, inplace=True)

    pvals_df.reset_index(inplace=True)
    pvals_df = pd.melt(pvals_df, id_vars=["Epitope"], value_name="P_value")
    pvals_df["Pair"] = pvals_df["Epitope"] + "-" + pvals_df["TCR"]

    pair_fold_changes_df = pair_fold_changes_df.merge(pvals_df, on="Pair", suffixes=("", "_left"))
    pair_fold_changes_df.drop(["Epitope_left", "TCR_left"], axis=1, inplace=True)

    exp_pair_prob_matrix_df.reset_index(inplace=True)
    exp_pair_prob_matrix_df = pd.melt(exp_pair_prob_matrix_df, id_vars=["Epitope"], value_name="Exp_pair_prob")
    exp_pair_prob_matrix_df["Pair"] = exp_pair_prob_matrix_df["Epitope"] + "-" + exp_pair_prob_matrix_df["TCR"]

    pair_fold_changes_df = pair_fold_changes_df.merge(exp_pair_prob_matrix_df, on="Pair", suffixes=("", "_left"))
    pair_fold_changes_df.drop(["Epitope_left", "TCR_left"], axis=1, inplace=True)

    print(
        pair_fold_changes_df[pair_fold_changes_df["Fold_change"] > 0].shape[0],
        "TCR-epitope pairs have non-zero well counts of",
        pair_fold_changes_df.shape[0],
        "tested TCR-epitope pairs",
    )
    pair_fold_changes_df = pair_fold_changes_df[pair_fold_changes_df["Fold_change"] > 0]

    if sort_fold_change:
        pair_fold_changes_df = pair_fold_changes_df.sort_values("Fold_change", ascending=False)
    elif sort_pvalue:
        if fdr_correction:
            pair_fold_changes_df = pair_fold_changes_df.sort_values(["P_value_corrected", "P_value"], ascending=True)
        else:
            pair_fold_changes_df = pair_fold_changes_df.sort_values("P_value", ascending=True)

    pair_fold_changes_df.reset_index(inplace=True, drop=True)

    return pair_fold_changes_df


def check_for_reactive_partner(partner_list, reactive_partner_list) -> bool:
    return any(partner in partner_list for partner in reactive_partner_list)


def add_reactive_epi_and_tcr_co_occurrence_bool_col(
    pair_fold_changes_df: pd.DataFrame,
    epi_co_occurrence_csv: str | os.PathLike,
    tcr_co_occurrence_csv: str | os.PathLike,
    validating_pair_dict: dict,
    epi_co_occur_threshold: Union[int, float],
    tcr_co_occur_threshold: Union[int, float],
) -> pd.DataFrame:
    """
    Flag detected TCRs and epitopes that frequently co-occur in the same wells as known-reactive positive control
    TCRs and epitopes.

    A detected epitope gets flagged (epi_co_occurs_with_reactive_epi = True) when it co-occurs with a known-reactive
    epitope (from validating_pair_dict) in more wells than epi_co_occur_threshold, and is itself not one of the
    known-reactive epitopes. Analogously, a detected TCR gets flagged (tcr_co_occurs_with_reactive_tcr = True) when
    it co-occurs with a known-reactive TCR in more wells than tcr_co_occur_threshold. This is used to identify
    TCR-epitope pairs whose apparent significance may instead be explained by frequent well co-occurrence with a
    true reactive pair, rather than by the pair itself being reactive.

    Parameters
    ----------
    pair_fold_changes_df : pd.DataFrame
        pair_fold_changes_df DataFrame generated by :func:`obs_tcr_epi_count_binom_sf_test`.
    epi_co_occurrence_csv : str | os.PathLike
        String path to the epitopes_co_occurrence_counts_<plates_list_name>.csv file written by
        :func:`write_tcr_epi_pair_well_count_files`, listing well co-occurrence counts between pairs of epitopes.
    tcr_co_occurrence_csv : str | os.PathLike
        String path to the tcr_co_occurrence_counts_<plates_list_name>.csv file written by
        :func:`write_tcr_epi_pair_well_count_files`, listing well co-occurrence counts between pairs of TCRs.
    validating_pair_dict : dict
        Dictionary where keys are epitope-TCR pair strings ("<epitope_reference_name>-<tcr_reference_name>") and
        values are booleans indicating whether that pair was confirmed reactive (True) or non-reactive (False) in an
        independent 1-to-1 validation assay. Must pass :func:`validate_validating_pair_dict`. The epitopes and TCRs of
        reactive pairs with value True are used as the reference set of "reactive" transcripts for co-occurrence flagging.
    epi_co_occur_threshold : int | float
        Minimum well co-occurrence count with a known-reactive epitope (strictly greater than) required to flag a
        detected epitope as epi_co_occurs_with_reactive_epi.
    tcr_co_occur_threshold : int | float
        Minimum well co-occurrence count with a known-reactive TCR (strictly greater than) required to flag a
        detected TCR as tcr_co_occurs_with_reactive_tcr.

    Returns
    -------
    pair_fold_changes_df
        Input DataFrame with two additional boolean columns added:
            - epi_co_occurs_with_reactive_epi: True if the pair's Epitope frequently co-occurs with a
              known-reactive epitope.
            - tcr_co_occurs_with_reactive_tcr: True if the pair's TCR frequently co-occurs with a known-reactive
              TCR.
    """
    pair_fold_changes_df = pair_fold_changes_df.copy()
    validate_validating_pair_dict(validating_pair_dict)

    reactive_epi_list = np.unique([pair.split("-")[0] for pair in validating_pair_dict if validating_pair_dict[pair] is True])
    reactive_tcr_list = np.unique([pair.split("-")[1] for pair in validating_pair_dict if validating_pair_dict[pair] is True])

    epi_co_occur_df = read_gdna_counts_csv(count_file=epi_co_occurrence_csv, count_file_ref_name_col="co_occurring_epitopes")
    epi_co_occur_df["co_occurring_epitopes_lists"] = epi_co_occur_df.index.str.split("-")
    epi_co_occur_df["co_occurs_with_reactive"] = epi_co_occur_df["co_occurring_epitopes_lists"].apply(lambda x: check_for_reactive_partner(x, reactive_epi_list))
    epi_co_occur_df = epi_co_occur_df[(epi_co_occur_df["co_occurs_with_reactive"]) & (epi_co_occur_df["count"] > epi_co_occur_threshold)]
    unique_co_occur_epi_list = list(set(epi_co_occur_df.explode("co_occurring_epitopes_lists")["co_occurring_epitopes_lists"].tolist()))
    unique_co_occur_epi_list = [epi for epi in unique_co_occur_epi_list if epi not in reactive_epi_list]

    tcr_co_occur_df = read_gdna_counts_csv(count_file=tcr_co_occurrence_csv, count_file_ref_name_col="co_occurring_tcr")
    tcr_co_occur_df["co_occurring_tcr_lists"] = tcr_co_occur_df.index.str.split("-")
    tcr_co_occur_df["co_occurs_with_reactive"] = tcr_co_occur_df["co_occurring_tcr_lists"].apply(lambda x: check_for_reactive_partner(x, reactive_tcr_list))
    tcr_co_occur_df = tcr_co_occur_df[(tcr_co_occur_df["co_occurs_with_reactive"]) & (tcr_co_occur_df["count"] > tcr_co_occur_threshold)]
    unique_co_occur_tcr_list = list(set(tcr_co_occur_df.explode("co_occurring_tcr_lists")["co_occurring_tcr_lists"].tolist()))
    unique_co_occur_tcr_list = [tcr for tcr in unique_co_occur_tcr_list if tcr not in reactive_tcr_list]

    pair_fold_changes_df["epi_co_occurs_with_reactive_epi"] = [epi in unique_co_occur_epi_list for epi in pair_fold_changes_df["Epitope"]]

    pair_fold_changes_df["tcr_co_occurs_with_reactive_tcr"] = [tcr in unique_co_occur_tcr_list for tcr in pair_fold_changes_df["TCR"]]

    return pair_fold_changes_df


def subset_tuple_counter(counter: Union[collections.defaultdict, dict], transcript: str) -> dict:
    return {tuple_key: count for tuple_key, count in counter.items() if transcript in tuple_key}


def count_transcript_with_prefix_in_tuple_counter(counter: Union[collections.Counter, dict], prefix: str) -> collections.Counter:
    prefix_counter = collections.Counter()
    for tuple_key, count in counter.items():
        for transcript in tuple_key:
            if transcript.startswith(prefix):
                prefix_counter.update([transcript] * count)
    return prefix_counter


def is_any_transcript_ratio_above_threshold(counter: Union[collections.Counter, dict], transcript: str, transcript_ratio: float) -> bool:
    # If this ever becomes slow, a more efficient solution would probably be to call a counter.most_common and
    # take the most abundant other_transcript which should be at index 2 of the counter.
    for other_transcript in counter:
        if other_transcript != transcript and counter[other_transcript] / counter[transcript] >= transcript_ratio:
            return True
    return False


def add_tcr_epi_pairs_only_supported_by_co_occurrence_bool_col(
    project_dir: str | os.PathLike, plates_list_name: str, pair_fold_changes_df: pd.DataFrame, epi_transcript_ratio: float = 0.40, tcr_transcript_ratio: float = 0.40
) -> pd.DataFrame:
    """
    Epi-Epi and TCR-TCR co-occurrence filtering removes TCR-Ag pairs with insufficient independent supporting wells.
    For each TCR, a TCR-Ag pair is removed if a second Ag co-occurs in these wells above a set threshold
    (epi_transcript_ratio) of its supporting wells, ensuring that only TCR-Ag pairs with sufficient independent well
    evidence are retained. TCR-TCR co-occurrence is filtered in the same manner using the tcr_transcript_ratio.

    In addition, for every significant (BH_Reject) TCR or epitope, this function flags
    TCR_significant_with_more_than_one_epi / epi_significant_with_more_than_one_TCR = True when that TCR/epitope
    remains significantly paired with more than one epitope/TCR after excluding co-occurrence-only-supported
    pairs, which indicates a TCR (or epitope) that is cross-detected with multiple partners.

    Parameters
    ----------
    project_dir : str | os.PathLike
        String path to plate-based TCR and epitope sequencing project directory. Used to read
        outs/<plates_list_name>/transcript_tuple_well_counter_<plates_list_name>.json.
    plates_list_name : str
        Name of the plates_list used when writing the pair count files with
        :func:`write_tcr_epi_pair_well_count_files` (e.g., 'run-4-1').
    pair_fold_changes_df : pd.DataFrame
        DataFrame that can be used to identify significantly enriched TCR-epitope pairs, as returned by
        :func:`obs_tcr_epi_count_binom_sf_test` (optionally passed through
        :func:`add_reactive_epi_and_tcr_co_occurrence_bool_col`). Must contain "TCR", "Epitope", and "BH_Reject"
        columns.
    epi_transcript_ratio : float, default = 0.40
        Minimum fraction of a TCR-epitope pair's supporting wells in which a competing epitope must also be detected
        to flag that pair as epi_only_supported_by_co_occurrence. Ranges from 0.0
        (most stringent: any co-occurrence at all with a competing epitope is enough to flag the pair) to
        1.0 (least stringent: a competing epitope must be present in every one of the pair's supporting wells to
        flag it).
    tcr_transcript_ratio : float, default = 0.40
        Minimum fraction of a TCR-epitope pair's supporting wells in which a competing TCR must also be detected
        to flag that pair as tcr_only_supported_by_co_occurrence. Ranges from 0.0
        (most stringent: any co-occurrence at all with a competing TCR is enough to flag the pair) to
        1.0 (least stringent: a competing TCR must be present in every one of the pair's supporting wells to
        flag it).

    Returns
    -------
    pair_fold_changes_df
        Input DataFrame with four additional boolean columns added:
            - epi_only_supported_by_co_occurrence: True if the pair's Epitope is likely only detected due to
              co-occurrence with another epitope.
            - tcr_only_supported_by_co_occurrence: True if the pair's TCR is likely only detected due to
              co-occurrence with another TCR.
            - TCR_significant_with_more_than_one_epi: True if the pair's TCR remains significantly paired with more
              than one epitope after excluding co-occurrence-only-supported pairs.
            - epi_significant_with_more_than_one_TCR: True if the pair's Epitope remains significantly paired with
              more than one TCR after excluding co-occurrence-only-supported pairs.

    """

    with open(Path(project_dir) / "outs" / plates_list_name / f"transcript_tuple_well_counter_{plates_list_name}.json", "r") as json_file:
        transcript_tuple_well_counter = json.load(json_file)

    # tuple keys are not supported in json --> we wrote the tuple keys as str and convert back to tuple:
    transcript_tuple_well_counter = {tuple(eval(tuple_key)): count for tuple_key, count in transcript_tuple_well_counter.items()}

    pair_fold_changes_df["tcr_only_supported_by_co_occurrence"] = False
    pair_fold_changes_df["epi_only_supported_by_co_occurrence"] = False

    for tcr in pair_fold_changes_df.loc[:, "TCR"].unique():
        tcr_epi_counter = count_transcript_with_prefix_in_tuple_counter(subset_tuple_counter(transcript_tuple_well_counter, tcr), prefix="epi_")
        tcr_epi_counter_list = tcr_epi_counter.most_common()
        if len(tcr_epi_counter_list) == 1:
            continue

        for epi_count_tuple in tcr_epi_counter_list:
            tmp_subset_epi_counter = subset_tuple_counter(subset_tuple_counter(transcript_tuple_well_counter, tcr), epi_count_tuple[0])
            tmp_subset_epi_epi_co_occur_counter = count_transcript_with_prefix_in_tuple_counter(tmp_subset_epi_counter, prefix="epi_")
            if is_any_transcript_ratio_above_threshold(tmp_subset_epi_epi_co_occur_counter, epi_count_tuple[0], transcript_ratio=epi_transcript_ratio):
                pair_fold_changes_df.loc[(pair_fold_changes_df["TCR"] == tcr) & (pair_fold_changes_df["Epitope"] == epi_count_tuple[0]), "epi_only_supported_by_co_occurrence"] = (
                    True
                )

    for epi in pair_fold_changes_df.loc[:, "Epitope"].unique():
        epi_tcr_counter = count_transcript_with_prefix_in_tuple_counter(subset_tuple_counter(transcript_tuple_well_counter, epi), prefix="tcr_")
        epi_tcr_counter_list = epi_tcr_counter.most_common()
        if len(epi_tcr_counter_list) == 1:
            continue

        for tcr_count_tuple in epi_tcr_counter_list:
            tmp_subset_tcr_counter = subset_tuple_counter(subset_tuple_counter(transcript_tuple_well_counter, epi), tcr_count_tuple[0])
            tmp_subset_tcr_tcr_co_occur_counter = count_transcript_with_prefix_in_tuple_counter(tmp_subset_tcr_counter, prefix="tcr_")
            if is_any_transcript_ratio_above_threshold(tmp_subset_tcr_tcr_co_occur_counter, tcr_count_tuple[0], transcript_ratio=tcr_transcript_ratio):
                pair_fold_changes_df.loc[(pair_fold_changes_df["Epitope"] == epi) & (pair_fold_changes_df["TCR"] == tcr_count_tuple[0]), "tcr_only_supported_by_co_occurrence"] = (
                    True
                )

    pair_fold_changes_df["TCR_significant_with_more_than_one_epi"] = False
    pair_fold_changes_df["epi_significant_with_more_than_one_TCR"] = False

    for tcr in pair_fold_changes_df.loc[pair_fold_changes_df["BH_Reject"], "TCR"].unique():
        if (
            len(
                pair_fold_changes_df.loc[
                    (pair_fold_changes_df["TCR"] == tcr) & ~pair_fold_changes_df["epi_only_supported_by_co_occurrence"] & pair_fold_changes_df["BH_Reject"], "Epitope"
                ].unique()
            )
            > 1
        ):
            pair_fold_changes_df.loc[pair_fold_changes_df["TCR"] == tcr, "TCR_significant_with_more_than_one_epi"] = True

    for epi in pair_fold_changes_df.loc[pair_fold_changes_df["BH_Reject"], "Epitope"].unique():
        if (
            len(
                pair_fold_changes_df.loc[
                    (pair_fold_changes_df["Epitope"] == epi) & ~pair_fold_changes_df["tcr_only_supported_by_co_occurrence"] & pair_fold_changes_df["BH_Reject"], "TCR"
                ].unique()
            )
            > 1
        ):
            pair_fold_changes_df.loc[pair_fold_changes_df["Epitope"] == epi, "epi_significant_with_more_than_one_TCR"] = True

    return pair_fold_changes_df


def _get_ma_color_alpha_and_annotation(
    pair: str,
    reject: bool,
    validated: bool | None,
    epi_co_occur_bool: bool,
    tcr_co_occur_bool: bool,
    custom_pair_color_dict: dict,
    annotate: bool,
    ax: mpl.axes.Axes,
    exp_pair_prob: float,
    fold_change: float,
    annotation_fontsize: float,
    coloring_mode: str = "significance",
) -> tuple[str, float, mpl.text.Text | None]:
    """
    coloring_mode options:
        - "significance": #8dc63f if reject else grey
        - "pair_dict": color from pair_color_dict if present else grey
        - "co_occurrence": colored by co-occurrence support
        - "co_occurrence_with_pair_dict": co-occurrence coloring with pair_dict override
        - "validated": colored by validating_pair status
        - "reactive_co_occurrence": colored by reactive co-occurrence with pair_dict override
        - "significance_reactive_co_occurrence": colored by significance and reactive co-occurrence
    """
    annotate_dot = False

    if coloring_mode == "significance":
        if reject:
            color, alpha = "#97C654", 1.0
            annotate_dot = True
        else:
            color, alpha = "grey", 0.2

    elif coloring_mode == "pair_dict":
        if pair in custom_pair_color_dict:
            color, alpha = custom_pair_color_dict[pair], 1.0
            annotate_dot = True
        else:
            color, alpha = "grey", 0.2

    elif coloring_mode == "co_occurrence":
        if epi_co_occur_bool and tcr_co_occur_bool and reject:
            color, alpha = "darkred", 1.0
            annotate_dot = True
        elif epi_co_occur_bool and reject:
            color, alpha = "red", 1.0
            annotate_dot = True
        elif tcr_co_occur_bool and reject:
            color, alpha = "black", 1.0
            annotate_dot = True
        elif reject:
            color, alpha = "#97C654", 1.0
            annotate_dot = True
        else:
            color, alpha = "grey", 0.1

    elif coloring_mode == "co_occurrence_with_pair_dict":
        if pair in custom_pair_color_dict:
            color, alpha = custom_pair_color_dict[pair], 1.0
            annotate_dot = True
        elif epi_co_occur_bool and tcr_co_occur_bool and reject:
            color, alpha = "darkred", 1.0
            annotate_dot = True
        elif epi_co_occur_bool and reject:
            color, alpha = "red", 1.0
            annotate_dot = True
        elif tcr_co_occur_bool and reject:
            color, alpha = "black", 1.0
            annotate_dot = True
        elif reject:
            color, alpha = "#97C654", 1.0
            annotate_dot = True
        else:
            color, alpha = "grey", 0.1

    elif coloring_mode == "validated":
        if pair in custom_pair_color_dict:
            color, alpha = custom_pair_color_dict[pair], 1.0
            annotate_dot = True
        elif pd.notna(validated) and validated:
            color, alpha = "#97C654", 1.0
            annotate_dot = True
        elif pd.notna(validated) and not validated:
            color, alpha = "#f7941d", 1.0
            annotate_dot = True
        elif pd.isna(validated) and reject:
            color, alpha = "#5C6C46", 1.0
            annotate_dot = True
        else:
            color, alpha = "grey", 0.1

    elif coloring_mode == "reactive_co_occurrence":
        if pair in custom_pair_color_dict:
            color, alpha = custom_pair_color_dict[pair], 1.0
            annotate_dot = True
        elif epi_co_occur_bool and tcr_co_occur_bool:
            color, alpha = "darkred", 1.0
            annotate_dot = reject
        elif epi_co_occur_bool:
            color, alpha = "red", 1.0
            annotate_dot = reject
        elif tcr_co_occur_bool:
            color, alpha = "black", 1.0
            annotate_dot = reject
        else:
            color, alpha = "grey", 0.1

    elif coloring_mode == "significance_reactive_co_occurrence":
        if epi_co_occur_bool and tcr_co_occur_bool:
            color, alpha = "darkred", 1.0
            annotate_dot = reject
        elif epi_co_occur_bool:
            color, alpha = "red", 1.0
            annotate_dot = reject
        elif tcr_co_occur_bool:
            color, alpha = "black", 1.0
            annotate_dot = reject
        elif reject:
            color, alpha = "#8dc63f", 1.0
            annotate_dot = True
        else:
            color, alpha = "grey", 0.1

    else:
        raise ValueError(f"Unknown coloring_mode: {coloring_mode}")

    text = ax.text(exp_pair_prob, fold_change, pair, fontsize=annotation_fontsize, zorder=8) if annotate and annotate_dot else None

    return color, alpha, text


def _plot_ma_scatter(
    plot_df: pd.DataFrame,
    ax: mpl.axes.Axes,
    fig: mpl.figure.Figure,
    xticks_list: list,
    xlim_min: float,
    xlim_max: float,
    ylim_max: float,
    title: str,
    annotate: bool,
    annotation_fontsize: float,
    pair_color_dict: dict,
    save_path: Path,
    coloring_mode: str = "significance",
    **kwargs,
) -> None:

    if (plot_df["Fold_change"] > ylim_max).any():
        print(f"Warning: ylim_max of {ylim_max} is smaller than the Fold_change of {(plot_df['Fold_change'] > ylim_max).sum()} pair(s)")
    if (plot_df["Fold_change"] < -1).any():
        print(f"Warning: ylim_min of -1 is larger than the Fold_change of {(plot_df['Fold_change'] < -1).sum()} pair(s)")
    if (plot_df["Exp_pair_prob"] > xlim_max).any():
        print(f"Warning: xlim_max of {xlim_max} is smaller than the Exp_pair_prob of {(plot_df['Exp_pair_prob'] > xlim_max).sum()} pair(s)")
    if (plot_df["Exp_pair_prob"] < xlim_min).any():
        print(f"Warning: xlim_min of {xlim_min} is larger than the Exp_pair_prob of {(plot_df['Exp_pair_prob'] < xlim_min).sum()} pair(s)")

    dot_color_list = []
    dot_alpha_list = []
    text_list = [] if annotate else None

    for pair, reject, validated, epi_co_occur_bool, tcr_co_occur_bool, exp_pair_prob, fold_change in zip(
        plot_df["Pair"],
        plot_df["BH_Reject"],
        plot_df["validating_pair"],
        plot_df["epi_co_occurs_with_reactive_epi"]
        if coloring_mode in ("reactive_co_occurrence", "significance_reactive_co_occurrence")
        else plot_df["epi_only_supported_by_co_occurrence"],
        plot_df["tcr_co_occurs_with_reactive_tcr"]
        if coloring_mode in ("reactive_co_occurrence", "significance_reactive_co_occurrence")
        else plot_df["tcr_only_supported_by_co_occurrence"],
        plot_df["Exp_pair_prob"],
        plot_df["Fold_change"],
    ):
        color, alpha, text = _get_ma_color_alpha_and_annotation(
            pair, reject, validated, epi_co_occur_bool, tcr_co_occur_bool, pair_color_dict, annotate, ax, exp_pair_prob, fold_change, annotation_fontsize, coloring_mode
        )
        dot_color_list.append(color)
        dot_alpha_list.append(alpha)
        if text is not None:
            text_list.append(text)

    ax.grid(zorder=-4)
    ax.scatter(plot_df["Exp_pair_prob"], plot_df["Fold_change"], c=dot_color_list, alpha=dot_alpha_list, edgecolors="none", zorder=4, s=8)
    if annotate:
        adjust_text(text_list, ax=ax, arrowprops=dict(arrowstyle="-", color="k", zorder=4, lw=0.25), **kwargs)

    ax.set_xticks(xticks_list)
    ax.set_xticklabels([round(x, 5) for x in ax.get_xticks()], fontsize=7, rotation=90)
    ax.set_ylabel("Observed/Expected pair fraction", fontsize=7)
    ax.set_xlabel("Expected pair fraction (when random)", fontsize=7)
    ax.set_ylim(-1, ylim_max)
    ax.set_xlim(xlim_min, xlim_max)
    ax.set_title(title, fontsize=7)
    ax.tick_params(top=False, right=False)
    ax.tick_params("both", labelsize=7)
    plt.tight_layout()
    fig.savefig(save_path)
    plt.close()


def _prepare_ma_plot_df(
    pair_fold_changes_df: pd.DataFrame,
    co_occurrence_filter: bool = False,
    cross_detected_filter: bool = False,
    sort_by_pair_color_dict: dict = None,
    sort_by_reactive_co_occurrence: bool = False,
) -> pd.DataFrame:
    plot_df = pair_fold_changes_df.copy()

    if co_occurrence_filter:
        plot_df = plot_df.loc[~((plot_df["epi_only_supported_by_co_occurrence"]) | (plot_df["tcr_only_supported_by_co_occurrence"]))]

    plot_df.sort_values("P_value_corrected", inplace=True)

    if cross_detected_filter:
        plot_df.drop_duplicates("TCR", keep="first", inplace=True)

    if sort_by_pair_color_dict:
        tmp_plot_df = plot_df[plot_df["Pair"].isin(list(sort_by_pair_color_dict.keys()))]
        plot_df = plot_df.drop(tmp_plot_df.index)
        plot_df = pd.concat([plot_df, tmp_plot_df])
    elif sort_by_reactive_co_occurrence:
        tmp_plot_df = plot_df.loc[plot_df["BH_Reject"] & ((plot_df["epi_co_occurs_with_reactive_epi"]) | (plot_df["tcr_co_occurs_with_reactive_tcr"])), :]
        plot_df = plot_df.drop(tmp_plot_df.index)
        plot_df = pd.concat([plot_df, tmp_plot_df])
    else:
        tmp_plot_df = plot_df.loc[plot_df["BH_Reject"], :].copy()
        plot_df = plot_df.drop(tmp_plot_df.index)
        plot_df = pd.concat([plot_df, tmp_plot_df])

    return plot_df


def plot_tcr_epi_pair_ma_plot(
    pair_fold_changes_df: pd.DataFrame,
    title: str,
    xlim_max: float,
    xticks_list: list,
    xlim_min: float,
    ylim_max: float,
    save_dir: str | os.PathLike,
    custom_pair_color_dict: dict = None,
    plot_filter_only_supported_co_occurrence: bool = True,
    plot_reactive_pair_co_occurrence: bool = False,
    annotate: bool = True,
    annotation_fontsize: float = 2,
    w: float = 6,
    h: float = 6,
    **kwargs,
) -> None:
    """
    Generate MA-like plots comparing each detected TCR-epitope pair's expected pair frequency (in case of random
    pairing) with its fold enrichment (observed/expected pair frequency), colored and filtered under several
    schemes, and write each as a .pdf to save_dir.

    Output files: MA_plot_<filters>_colored[_annotated].pdf (up to 10 files). Filter/annotation tags in the
    filenames mean the following:
        - only_supported_by_co_occurrence: pair dots are removed if the co-occurrence filter found insufficient
          independent well evidence for that pair (a competing epitope/TCR co-occurred in >= threshold fraction
          of the pair's supporting wells).
        - cross_detected: pair dots are removed if their TCR is also part of another, more significant pair with a
          different Ag; only the most significant pair per TCR is kept.
        - significance: pair dots are colored if P_value_corrected (also called "PAIR-Scan confidence in a pair
          being reactive") is below a user-defined FDR significance threshold (alpha, set in
          :func:`obs_tcr_epi_count_binom_sf_test`).
        - validated: pair dots are colored by individual pair validation co-culture result (via validating_pair_dict).
        - _annotated: pair dots are labeled with text describing the pair.

    Parameters
    ----------
    pair_fold_changes_df : pd.DataFrame
        DataFrame with per-pair Observed/Expected fold changes, P values, and co-occurrence/validation annotations,
        as returned by :func:`add_tcr_epi_pairs_only_supported_by_co_occurrence_bool_col`.
    title : str
        Title used for every generated plot.
    xlim_max : float
        Maximum x-axis (Expected pair fraction) limit.
    xticks_list : list
        Tick positions for the x-axis (Expected pair fraction).
    xlim_min : float
        Minimum x-axis (Expected pair fraction) limit.
    ylim_max : float
        Maximum y-axis (Observed/Expected pair fraction fold change) limit. The y-axis minimum is fixed at -1.
    save_dir : str | os.PathLike
        Directory the generated .pdf plots (and, if custom_pair_color_dict is provided, a color legend .pdf) are
        written to.
    custom_pair_color_dict : dict, optional
        Dictionary mapping specific TCR-epitope pair strings to a matplotlib color. When provided, additional
        "pair_dict"-colored (and, if plot_filter_only_supported_co_occurrence, "co_occurrence_with_pair_dict"-
        colored) plots are generated, and pairs in this dict are always drawn on top and annotated.
    plot_filter_only_supported_co_occurrence : bool, default = True
        If True, also generate the co-occurrence-colored and co-occurrence-/cross-detected-filtered plot variants
        (plots 3-8 in the implementation). Requires pair_fold_changes_df to have the
        *_only_supported_by_co_occurrence columns.
    plot_reactive_pair_co_occurrence : bool, default = False
        If True, also generate the reactive-co-occurrence-colored plot variants (plots 9-10). Requires
        custom_pair_color_dict to be provided and pair_fold_changes_df to have the
        epi_co_occurs_with_reactive_epi / tcr_co_occurs_with_reactive_tcr columns from
        :func:`add_reactive_epi_and_tcr_co_occurrence_bool_col`.
    annotate : bool, default = True
        If True, annotate colored (non-grey) dots with their pair string, using adjustText to avoid overlaps.
    annotation_fontsize : float, default = 2
        Fontsize used for pair annotations when annotate is True.
    w : float, default = 6
        Figure width, passed to :func:`tcr_toolbox.utils.plot_utils.startfig`.
    h : float, default = 6
        Figure height, passed to :func:`tcr_toolbox.utils.plot_utils.startfig`.
    **kwargs
        Additional keyword arguments forwarded to adjustText.adjust_text when annotate is True.

    Returns
    -------
    None
        Writes one or more "MA_plot_*.pdf" files (and, if applicable, a pair-color legend .pdf) to save_dir.
    """

    save_dir = Path(save_dir)
    if custom_pair_color_dict is None:
        custom_pair_color_dict = {}
    annotated_str = "_annotated" if annotate else ""

    # Plot 1: significance colored
    plot_df = _prepare_ma_plot_df(pair_fold_changes_df)
    ax, fig, gs = startfig(w, h)
    _plot_ma_scatter(
        plot_df,
        ax,
        fig,
        xticks_list,
        xlim_min,
        xlim_max,
        ylim_max,
        title,
        annotate,
        annotation_fontsize,
        custom_pair_color_dict,
        save_dir / f"MA_plot_significance_colored{annotated_str}.pdf",
        coloring_mode="significance",
        **kwargs,
    )

    # Plot 2: pair_dict colored
    if custom_pair_color_dict:
        view_color_dict(custom_pair_color_dict, label_length_scaler=3.25, save_path=save_dir / "MA_plot_pair_dict_color_legend.pdf")
        plot_df = _prepare_ma_plot_df(pair_fold_changes_df, sort_by_pair_color_dict=custom_pair_color_dict)
        ax, fig, gs = startfig(w, h)
        _plot_ma_scatter(
            plot_df,
            ax,
            fig,
            xticks_list,
            xlim_min,
            xlim_max,
            ylim_max,
            title,
            annotate,
            annotation_fontsize,
            custom_pair_color_dict,
            save_dir / f"MA_plot_pair_dict_colored{annotated_str}.pdf",
            coloring_mode="pair_dict",
            **kwargs,
        )

    if plot_filter_only_supported_co_occurrence:
        # Plot 3: co-occurrence colored
        plot_df = _prepare_ma_plot_df(pair_fold_changes_df)
        ax, fig, gs = startfig(w, h)
        _plot_ma_scatter(
            plot_df,
            ax,
            fig,
            xticks_list,
            xlim_min,
            xlim_max,
            ylim_max,
            title,
            annotate,
            annotation_fontsize,
            custom_pair_color_dict,
            save_dir / f"MA_plot_significance_and_only_supported_by_co_occurrence_colored{annotated_str}.pdf",
            coloring_mode="co_occurrence",
            **kwargs,
        )

        # Plot 4: co-occurrence with pair_dict colored
        if custom_pair_color_dict:
            plot_df = _prepare_ma_plot_df(pair_fold_changes_df)
            ax, fig, gs = startfig(w, h)
            _plot_ma_scatter(
                plot_df,
                ax,
                fig,
                xticks_list,
                xlim_min,
                xlim_max,
                ylim_max,
                title,
                annotate,
                annotation_fontsize,
                custom_pair_color_dict,
                save_dir / f"MA_plot_significance_pair_and_only_supported_by_co_occurrence_colored{annotated_str}.pdf",
                coloring_mode="co_occurrence_with_pair_dict",
                **kwargs,
            )

        # Plot 5: co-occurrence filtered, significance colored
        plot_df = _prepare_ma_plot_df(pair_fold_changes_df, co_occurrence_filter=True)
        ax, fig, gs = startfig(w, h)
        _plot_ma_scatter(
            plot_df,
            ax,
            fig,
            xticks_list,
            xlim_min,
            xlim_max,
            ylim_max,
            title,
            annotate,
            annotation_fontsize,
            custom_pair_color_dict,
            save_dir / f"MA_plot_only_supported_by_co_occurrence_filtered_and_significance_colored{annotated_str}.pdf",
            coloring_mode="significance",
            **kwargs,
        )

        # Plot 6: co-occurrence and cross-detected filtered, significance colored
        plot_df = _prepare_ma_plot_df(pair_fold_changes_df, co_occurrence_filter=True, cross_detected_filter=True)
        ax, fig, gs = startfig(w, h)
        _plot_ma_scatter(
            plot_df,
            ax,
            fig,
            xticks_list,
            xlim_min,
            xlim_max,
            ylim_max,
            title,
            annotate,
            annotation_fontsize,
            custom_pair_color_dict,
            save_dir / f"MA_plot_cross_detected_and_only_supported_by_co_occurrence_filtered_and_significance_colored{annotated_str}.pdf",
            coloring_mode="significance",
            **kwargs,
        )

        # Plot 7: co-occurrence and cross-detected filtered, pair_dict colored
        if custom_pair_color_dict:
            ax, fig, gs = startfig(w, h)
            _plot_ma_scatter(
                plot_df,
                ax,
                fig,
                xticks_list,
                xlim_min,
                xlim_max,
                ylim_max,
                title,
                annotate,
                annotation_fontsize,
                custom_pair_color_dict,
                save_dir / f"MA_plot_cross_detected_and_only_supported_by_co_occurrence_filtered_and_pair_dict_colored{annotated_str}.pdf",
                coloring_mode="pair_dict",
                **kwargs,
            )

            # Plot 8: pair_dict colored (repeated with co-occurrence filter context)
            plot_df = _prepare_ma_plot_df(pair_fold_changes_df, co_occurrence_filter=True, sort_by_pair_color_dict=custom_pair_color_dict)
            ax, fig, gs = startfig(w, h)
            _plot_ma_scatter(
                plot_df,
                ax,
                fig,
                xticks_list,
                xlim_min,
                xlim_max,
                ylim_max,
                title,
                annotate,
                annotation_fontsize,
                custom_pair_color_dict,
                save_dir / f"MA_plot_only_supported_by_co_occurrence_filtered_and_pair_dict_colored{annotated_str}.pdf",
                coloring_mode="pair_dict",
                **kwargs,
            )

    # Plot 9: reactive co-occurrence colored
    if plot_reactive_pair_co_occurrence:
        if not custom_pair_color_dict:
            raise Exception("Plot reactive pair co-occurrence requires custom_pair_color_dict!")

        plot_df = _prepare_ma_plot_df(pair_fold_changes_df, sort_by_reactive_co_occurrence=True)
        ax, fig, gs = startfig(w, h)
        _plot_ma_scatter(
            plot_df,
            ax,
            fig,
            xticks_list,
            xlim_min,
            xlim_max,
            ylim_max,
            title,
            annotate,
            annotation_fontsize,
            custom_pair_color_dict,
            save_dir / f"MA_plot_pair_dict_and_reactive_co_occurrence_colored{annotated_str}.pdf",
            coloring_mode="reactive_co_occurrence",
            **kwargs,
        )

        # Plot 10: significance and reactive co-occurrence colored
        plot_df = _prepare_ma_plot_df(pair_fold_changes_df)
        ax, fig, gs = startfig(w, h)
        _plot_ma_scatter(
            plot_df,
            ax,
            fig,
            xticks_list,
            xlim_min,
            xlim_max,
            ylim_max,
            title,
            annotate,
            annotation_fontsize,
            custom_pair_color_dict,
            save_dir / f"MA_plot_significance_and_reactive_co_occurrence_colored{annotated_str}.pdf",
            coloring_mode="significance_reactive_co_occurrence",
            **kwargs,
        )

    # Plot 11: co-occurrence and cross-detected filtered, validated colored
    plot_df = _prepare_ma_plot_df(pair_fold_changes_df, co_occurrence_filter=True, cross_detected_filter=True)
    ax, fig, gs = startfig(w, h)
    _plot_ma_scatter(
        plot_df,
        ax,
        fig,
        xticks_list,
        xlim_min,
        xlim_max,
        ylim_max,
        title,
        annotate,
        annotation_fontsize,
        custom_pair_color_dict,
        save_dir / f"MA_plot_cross_detected_and_only_supported_by_co_occurrence_filtered_and_significance_and_validated_colored{annotated_str}.pdf",
        coloring_mode="validated",
        **kwargs,
    )


def _get_pvalue_rank_color_alpha_and_annotation(
    reject: bool,
    validated: bool | None,
    pair: str,
    custom_pair_color_dict: dict,
    annotate: bool,
    ax: mpl.axes.Axes,
    p_rank: int,
    p_corrected: float,
    annotation_fontsize: float,
    use_validating_colors: bool = True,
) -> tuple[str, float, mpl.text.Text | None]:
    if use_validating_colors:
        if pair in custom_pair_color_dict:
            color, alpha = custom_pair_color_dict[pair], 1.0
        elif pd.notna(validated) and validated:
            color, alpha = "#97C654", 1.0
        elif pd.notna(validated) and not validated:
            color, alpha = "#f7941d", 1.0
        elif pd.isna(validated) and reject:
            color, alpha = "#5C6C46", 1.0
        else:
            color, alpha = "grey", 0.20
    else:
        if reject:
            color, alpha = "#97C654", 1.0
        else:
            color, alpha = "grey", 0.20

    text = ax.text(p_rank, p_corrected, pair, fontsize=annotation_fontsize, zorder=8) if annotate and alpha == 1.0 else None

    return color, alpha, text


def _plot_pvalue_rank_scatter(
    plot_df: pd.DataFrame,
    ax: mpl.axes.Axes,
    fig: mpl.figure.Figure,
    xlim_min: float,
    xlim_max: float,
    ylim_min: float,
    ylim_max: float,
    title: str,
    annotate: bool,
    annotation_fontsize: float,
    custom_pair_color_dict: dict,
    save_path: Path,
    adj_pval_alpha: float | None = None,
    use_validating_colors=True,
) -> None:
    color_list = []
    alpha_list = []
    text_list = [] if annotate else None

    for reject, validated, pair, p_rank, p_corrected in zip(
        plot_df["BH_Reject"], plot_df["validating_pair"], plot_df["Pair"], np.arange(0, len(plot_df["P_value_corrected"])), -np.log10(plot_df["P_value_corrected"])
    ):
        color, a, text = _get_pvalue_rank_color_alpha_and_annotation(
            reject, validated, pair, custom_pair_color_dict, annotate, ax, p_rank, p_corrected, annotation_fontsize, use_validating_colors
        )
        color_list.append(color)
        alpha_list.append(a)
        if text is not None:
            text_list.append(text)

    ax.scatter(np.arange(0, len(plot_df["P_value_corrected"])), -np.log10(plot_df["P_value_corrected"]), edgecolor="none", color=color_list, alpha=alpha_list, zorder=4, s=10)
    if annotate:
        adjust_text(text_list, ax=ax, arrowprops=dict(arrowstyle="-", color="k", zorder=4, lw=0.25))

    if adj_pval_alpha:
        ax.axhline(-np.log10(adj_pval_alpha), color="black", linewidth=0.5, zorder=6)

    ax.set_xlim(xlim_min, xlim_max)
    ax.set_ylim(ylim_min, ylim_max)
    ax.set_ylabel("-Log10(adjusted P value)", fontsize=7)
    ax.set_xlabel("adjusted P value rank", fontsize=7)
    ax.set_title(title, fontsize=7)
    ax.tick_params(top=False, right=False)
    ax.tick_params("both", labelsize=7)
    fig.tight_layout()
    fig.savefig(save_path)
    plt.close()


def _prepare_pvalue_rank_plot_df(
    pair_fold_changes_df: pd.DataFrame, filter_significant_with_more_than_one_epi: bool, ylim_max: float, xlim_max: float, co_occurrence_filter: bool = False
) -> pd.DataFrame:
    plot_df = pair_fold_changes_df.copy()
    if co_occurrence_filter:
        plot_df = plot_df.loc[~plot_df["tcr_only_supported_by_co_occurrence"] & ~plot_df["epi_only_supported_by_co_occurrence"], :]
    plot_df.sort_values("P_value_corrected", inplace=True)
    if filter_significant_with_more_than_one_epi:
        plot_df.drop_duplicates("TCR", keep="first", inplace=True)
    if (plot_df["P_value_corrected"] == 0).any():
        pseudo_count = plot_df.loc[plot_df["P_value_corrected"] > 0, "P_value_corrected"].min() * 0.1
        plot_df.loc[:, "P_value_corrected"] += pseudo_count
    if -np.log10(plot_df["P_value_corrected"]).min() > ylim_max:
        print(f"Warning: ylim_min is larger than minimal corrected log10(P value): {-np.log10(plot_df['P_value_corrected']).min()}")
    if np.max(np.arange(0, len(plot_df["P_value_corrected"]))) > xlim_max:
        print("Warning: xlim_max is smaller than maximal corrected log10(P value) rank:", np.max(np.arange(0, len(plot_df["P_value_corrected"]))))
    return plot_df


def plot_pvalue_corrected_rank_vs_pvalue_corrected(
    pair_fold_changes_df: pd.DataFrame,
    save_dir: str | os.PathLike,
    title: str,
    xlim_max: float,
    xlim_min: float,
    ylim_max: float,
    ylim_min: float,
    adj_pval_alpha: float | None = None,
    w: float = 6,
    h: float = 6,
    filter_significant_with_more_than_one_epi: bool = True,
    annotate: bool = True,
    annotation_fontsize: float = 2,
    custom_pair_color_dict: dict = None,
):
    """
    Generate rank plots comparing each detected TCR-epitope pair's -log10(adjusted P value) with its rank (lowest
    adjusted P value has the highest rank), colored and filtered under the same underlying schemes as the MA plots
    (see :func:`plot_tcr_epi_pair_ma_plot`), and write each as a .pdf to save_dir.

    Output files: p_value_corrected_vs_p_value_rank[_cross_detected_filtered][_co_occurrence_filtered]
    [_only_significance_colored][_annotated].pdf (up to 3 files). Filter/annotation tags in the filenames mean the
    following:
        - cross_detected_filtered: pair dots are removed if their TCR is also part of another, more significant
          pair with a different Ag; only the most significant pair per TCR is kept. Present by default (controlled
          by filter_significant_with_more_than_one_epi).
        - co_occurrence_filtered: pair dots are removed if the co-occurrence filter found insufficient independent
          well evidence for that pair (a competing epitope/TCR co-occurred in >= threshold fraction of the pair's
          supporting wells).
        - validated (default, no tag in filename): pair dots are colored by individual pair validation co-culture
          result (annotated in pair_fold_changes_df via validating_pair_dict).
        - only_significance_colored: pair dots are colored only by whether P_value_corrected (also called
          "PAIR-Scan confidence in a pair being reactive") is below a user-defined FDR significance
          threshold (alpha, set in :func:`obs_tcr_epi_count_binom_sf_test`), instead of by validation status.
        - _annotated: pair dots are labeled with text describing the pair.

    If any adjusted P values are exactly 0, a small pseudo-count (10% of the smallest non-zero adjusted P value) is
    added to all adjusted P values before taking the log so that they remain plottable.

    Parameters
    ----------
    pair_fold_changes_df : pd.DataFrame
        DataFrame with per-pair Observed/Expected fold changes, P values, and co-occurrence/validation annotations,
        as returned by :func:`add_tcr_epi_pairs_only_supported_by_co_occurrence_bool_col`.
    save_dir : str | os.PathLike
        Directory the generated .pdf plots are written to.
    title : str
        Title used for every generated plot.
    xlim_max : float
        Maximum x-axis (adjusted P value rank) limit.
    xlim_min : float
        Minimum x-axis (adjusted P value rank) limit.
    ylim_max : float
        Maximum y-axis (-log10 adjusted P value) limit.
    ylim_min : float
        Minimum y-axis (-log10 adjusted P value) limit.
    adj_pval_alpha : float, optional
        If provided, draw a horizontal reference line at -log10(adj_pval_alpha) on every plot (e.g., the
        significance threshold used elsewhere in the pipeline).
    w : float, default = 6
        Figure width, passed to :func:`tcr_toolbox.utils.plot_utils.startfig`.
    h : float, default = 6
        Figure height, passed to :func:`tcr_toolbox.utils.plot_utils.startfig`.
    filter_significant_with_more_than_one_epi : bool, default = True
        If True, keep only the most significant pair per TCR (duplicate TCRs dropped) before plotting, and append
        "_cross_detected_filtered" to output filenames.
    annotate : bool, default = True
        If True, annotate colored (non-grey/full-alpha) dots with their pair string, using adjustText to avoid
        overlaps.
    annotation_fontsize : float, default = 2
        Fontsize used for pair annotations when annotate is True.
    custom_pair_color_dict : dict, optional
        Dictionary mapping specific TCR-epitope pair strings to a matplotlib color, used instead of the default
        validation-based coloring for those pairs.

    Returns
    -------
    None
        Writes up to three "p_value_corrected_vs_p_value_rank*.pdf" files to save_dir.
    """
    save_dir = Path(save_dir)
    if custom_pair_color_dict is None:
        custom_pair_color_dict = {}

    cross_detected_filter_str = "_cross_detected_filtered" if filter_significant_with_more_than_one_epi else ""
    annotated_str = "_annotated" if annotate else ""

    # Plot 1: no co-occurrence filter
    plot_df = _prepare_pvalue_rank_plot_df(pair_fold_changes_df, filter_significant_with_more_than_one_epi, ylim_max, xlim_max, co_occurrence_filter=False)
    ax, fig, gs = startfig(w, h)
    _plot_pvalue_rank_scatter(
        plot_df=plot_df,
        ax=ax,
        fig=fig,
        adj_pval_alpha=adj_pval_alpha,
        xlim_min=xlim_min,
        xlim_max=xlim_max,
        ylim_min=ylim_min,
        ylim_max=ylim_max,
        title=title,
        annotate=annotate,
        annotation_fontsize=annotation_fontsize,
        custom_pair_color_dict=custom_pair_color_dict,
        save_path=save_dir / f"p_value_corrected_vs_p_value_rank{cross_detected_filter_str}{annotated_str}.pdf",
        use_validating_colors=True,
    )

    # Plot 2: with co-occurrence filter
    plot_df = _prepare_pvalue_rank_plot_df(pair_fold_changes_df, filter_significant_with_more_than_one_epi, ylim_max, xlim_max, co_occurrence_filter=True)

    ax, fig, gs = startfig(w, h)
    _plot_pvalue_rank_scatter(
        plot_df=plot_df,
        ax=ax,
        fig=fig,
        adj_pval_alpha=adj_pval_alpha,
        xlim_min=xlim_min,
        xlim_max=xlim_max,
        ylim_min=ylim_min,
        ylim_max=ylim_max,
        title=title,
        annotate=annotate,
        annotation_fontsize=annotation_fontsize,
        custom_pair_color_dict=custom_pair_color_dict,
        save_path=save_dir / f"p_value_corrected_vs_p_value_rank{cross_detected_filter_str}_co_occurrence_filtered{annotated_str}.pdf",
        use_validating_colors=True,
    )

    # Plot 3: co-occurrence filter, significance only coloring
    plot_df = _prepare_pvalue_rank_plot_df(pair_fold_changes_df, filter_significant_with_more_than_one_epi, ylim_max, xlim_max, co_occurrence_filter=True)
    ax, fig, gs = startfig(w, h)
    _plot_pvalue_rank_scatter(
        plot_df=plot_df,
        ax=ax,
        fig=fig,
        adj_pval_alpha=adj_pval_alpha,
        xlim_min=xlim_min,
        xlim_max=xlim_max,
        ylim_min=ylim_min,
        ylim_max=ylim_max,
        title=title,
        annotate=annotate,
        annotation_fontsize=annotation_fontsize,
        custom_pair_color_dict=custom_pair_color_dict,
        save_path=save_dir / f"p_value_corrected_vs_p_value_rank{cross_detected_filter_str}_co_occurrence_filtered_only_significance_colored{annotated_str}.pdf",
        use_validating_colors=False,
    )


def add_pair_entropy_over_plate_col(pair_fold_changes_df: pd.DataFrame, project_dir: Union[str, os.PathLike[str]]):

    pair_plate_counts = collections.defaultdict(lambda: collections.defaultdict(int))

    plate_name_list = []
    plate_files = sorted((Path(project_dir) / "outs" / "plate_maps").glob("*pairs_per_well.xlsx"))

    for pairs_well_xlsx in plate_files:
        plate_name = pairs_well_xlsx.stem.split("_pairs_per_well")[0]
        plate_name_list.append(plate_name)
        plate_df = pd.read_excel(pairs_well_xlsx, index_col=0)
        plate_df["pairs"] = plate_df["pairs"].apply(ensure_list)
        for pair_list in plate_df["pairs"]:
            if not pair_list:
                continue
            for pair in pair_list:
                pair_plate_counts[pair][plate_name] += 1

    if not len(set(plate_name_list)) == len(plate_name_list):
        raise Exception("Duplicate plate names detected in plate_map directory!")

    if len(plate_name_list) < 2:
        raise ValueError("Need at least two plates to compute entropy")

    pair_entropy = {}

    for pair, plate_dict in pair_plate_counts.items():
        counts = [plate_dict.get(p, 0) for p in plate_name_list]
        pair_entropy[pair] = normalized_entropy(counts)

    pair_fold_changes_df["pair_plate_entropy"] = pair_fold_changes_df["Pair"].map(pair_entropy).fillna(0)

    if pair_fold_changes_df["pair_plate_entropy"].isna().any():
        print("Warning: some pairs missing from plate maps")

    return pair_fold_changes_df


def filter_reactive_pairs(pair_fold_changes_df: pd.DataFrame):
    reactive_pair_df = pair_fold_changes_df.copy()
    reactive_pair_df.sort_values("P_value_corrected", ascending=True, inplace=True)

    if reactive_pair_df["Pair"].duplicated().any():
        raise Exception("There are duplicate Pairs in the pair_fold_changes_df Pair column!")

    reactive_pair_df = reactive_pair_df.loc[reactive_pair_df["BH_Reject"], :].copy()
    reactive_pair_not_co_occurrence_filtered_df = reactive_pair_df.copy()
    reactive_pair_df = reactive_pair_df.loc[~(reactive_pair_df["tcr_only_supported_by_co_occurrence"] | reactive_pair_df["epi_only_supported_by_co_occurrence"]), :].copy()

    cross_detected_idx = reactive_pair_df.loc[
        reactive_pair_df["TCR_significant_with_more_than_one_epi"]
        & reactive_pair_df.loc[reactive_pair_df["TCR_significant_with_more_than_one_epi"], "TCR"].duplicated(keep="first"),
        :,
    ].index
    cross_detected_idx_co_occurrence = reactive_pair_not_co_occurrence_filtered_df.loc[
        reactive_pair_not_co_occurrence_filtered_df["TCR_significant_with_more_than_one_epi"]
        & reactive_pair_not_co_occurrence_filtered_df.loc[reactive_pair_not_co_occurrence_filtered_df["TCR_significant_with_more_than_one_epi"], "TCR"].duplicated(keep="first"),
        :,
    ].index
    cross_detected_df = reactive_pair_df.loc[reactive_pair_df["TCR_significant_with_more_than_one_epi"], :]
    reactive_pair_df.drop(index=cross_detected_idx, inplace=True)
    reactive_pair_not_co_occurrence_filtered_df.drop(index=cross_detected_idx_co_occurrence, inplace=True)
    if reactive_pair_df["TCR"].duplicated().any():
        raise Exception("Duplicate/cross-detected TCRs are still present in reactive pairs filtered pair_fold_changes_df!")

    filtered_pair_df = pair_fold_changes_df.copy()
    filtered_pair_df.drop(index=reactive_pair_df.index, inplace=True)

    return reactive_pair_df, reactive_pair_not_co_occurrence_filtered_df, cross_detected_df, filtered_pair_df


def eval_performance_run(pair_fold_changes_df: pd.DataFrame, validating_pair_dict: dict):
    """
    Calculate sensitivity and precision of a PAIR-Scan screen against a set of validated reactive pairs in
    validating_pair_dict.

    Parameters
    ----------
    pair_fold_changes_df : pd.DataFrame
        DataFrame with one row per unique detected TCR-epitope pair (no duplicate "Pair" values), as returned by
        :func:`add_tcr_epi_pairs_only_supported_by_co_occurrence_bool_col`.
    validating_pair_dict : dict
        Dictionary where keys are epitope-TCR pair strings ("<epitope_reference_name>-<tcr_reference_name>") and
        values are booleans indicating whether that pair was confirmed reactive (True) or non-reactive (False) in
        an independent 1-to-1 validation assay. Must pass :func:`validate_validating_pair_dict`.

    Returns
    -------
    sensitivity : float
        Fraction of validated reactive pairs recovered in the pair hit set.
    sensitivity_without_co_occurrence_filter : float
        Fraction of validated reactive pairs recovered in the pair hit set before the co-occurrence
        filter is applied.
    precision : float
        Fraction of the pair hit set that is a validated reactive pair.
    """

    validate_validating_pair_dict(validating_pair_dict)

    eval_df = pair_fold_changes_df.copy()
    if eval_df["Pair"].duplicated().any():
        raise Exception("There are duplicate Pairs in the pair_fold_changes_df Pair column!")

    validating_pairs_set = {pair for pair, is_validating in validating_pair_dict.items() if is_validating is True}

    reactive_pair_df, reactive_pair_not_co_occurrence_filtered_df, _, _ = filter_reactive_pairs(pair_fold_changes_df=eval_df)

    detected_pairs_set = set(reactive_pair_df.loc[:, "Pair"].unique())
    significant_pairs_without_co_occurrence_filter_set = set(reactive_pair_not_co_occurrence_filtered_df.loc[:, "Pair"].unique())

    sensitivity = len(validating_pairs_set.intersection(detected_pairs_set)) / len(validating_pairs_set)
    sensitivity_without_co_occurrence_filter = len(validating_pairs_set.intersection(significant_pairs_without_co_occurrence_filter_set)) / len(validating_pairs_set)
    precision = len(detected_pairs_set.intersection(validating_pairs_set)) / len(detected_pairs_set)

    return sensitivity, sensitivity_without_co_occurrence_filter, precision


def write_reactive_pair_results(pair_fold_changes_df: pd.DataFrame, project_dir: Union[str, os.PathLike[str]], plates_list_name: str, validating_pair_dict: dict = None):
    """
    Write the screen pair hit set, related intermediate DataFrames, and a human-readable results summary to
    outs/<plates_list_name>.

    Runs :func:`filter_reactive_pairs` on pair_fold_changes_df and writes its four outputs to .xlsx files. If
    validating_pair_dict is provided, also runs :func:`eval_performance_run` and writes a detailed breakdown of
    which validated reactive pairs were recovered, missed, or explained by co-occurrence
    filtering/cross-detection to reactive_pair_result_<plates_list_name>.txt, alongside the run's sensitivity,
    sensitivity_without_co_occurrence_filter, and precision.

    Output files written to outs/<plates_list_name>:
        - reactive_pair_df_<plates_list_name>.xlsx: the pair hit set.
        - reactive_pair_not_co_occurrence_filtered_df_<plates_list_name>.xlsx: significant pair hit set before the
          co-occurrence filter.
        - cross_detected_df_<plates_list_name>.xlsx: pairs dropped for having a cross-detected TCR.
        - filtered_pair_df_<plates_list_name>.xlsx: all remaining (non-hit) pairs.
        - reactive_pair_result_<plates_list_name>.txt: human-readable summary (sensitivity/precision when
          validating_pair_dict is provided, identified hits, and cross-detected pairs).

    Parameters
    ----------
    pair_fold_changes_df : pd.DataFrame
        DataFrame with one row per unique detected TCR-epitope pair (no duplicate "Pair" values), as returned by
        :func:`add_tcr_epi_pairs_only_supported_by_co_occurrence_bool_col`.
    project_dir : str | os.PathLike
        String path to plate-based TCR and epitope sequencing project directory. Output files are written to
        project_dir/outs/<plates_list_name>.
    plates_list_name : str
        Name used for the output directory and output files (e.g., 'run-4-1', or a combined-screen name such as
        "screenA_vs_screenB_stouffer" when called from :func:`plot_and_combine_screen_pvals_stouffer`).
    validating_pair_dict : dict, optional
        Dictionary where keys are epitope-TCR pair strings ("<epitope_reference_name>-<tcr_reference_name>") and
        values are booleans indicating whether that pair was confirmed reactive (True) or non-reactive (False) in
        an independent 1-to-1 validation assay. Must pass :func:`validate_validating_pair_dict`. If provided,
        sensitivity/precision and detailed missing/newly-identified pair breakdowns are included in the results
        text file.

    Returns
    -------
    None
        Writes the files described above to project_dir/outs/<plates_list_name>.
    """
    plates_list_outs = Path(project_dir) / "outs" / plates_list_name

    pair_fold_changes_df = pair_fold_changes_df.copy()

    if pair_fold_changes_df["Pair"].duplicated().any():
        raise Exception("There are duplicate Pairs in the pair_fold_changes_df Pair column!")

    total_detected_pairs_set = set(pair_fold_changes_df["Pair"].unique())
    total_detected_epi_set = set(pair_fold_changes_df["Epitope"].unique())
    total_detected_tcr_set = set(pair_fold_changes_df["TCR"].unique())

    if validating_pair_dict is not None:
        validate_validating_pair_dict(validating_pair_dict)
        validating_pairs_set = {pair for pair, is_validating in validating_pair_dict.items() if is_validating is True}

    reactive_pair_df, reactive_pair_not_co_occurrence_filtered_df, cross_detected_df, filtered_pair_df = filter_reactive_pairs(pair_fold_changes_df=pair_fold_changes_df)
    reactive_pair_df.to_excel(plates_list_outs / f"reactive_pair_df_{plates_list_name}.xlsx")
    reactive_pair_not_co_occurrence_filtered_df.to_excel(plates_list_outs / f"reactive_pair_not_co_occurrence_filtered_df_{plates_list_name}.xlsx")
    cross_detected_df.to_excel(plates_list_outs / f"cross_detected_df_{plates_list_name}.xlsx")
    filtered_pair_df.to_excel(plates_list_outs / f"filtered_pair_df_{plates_list_name}.xlsx")

    if validating_pair_dict is not None:
        sensitivity, sensitivity_without_co_occurrence_filter, precision = eval_performance_run(
            pair_fold_changes_df=pair_fold_changes_df, validating_pair_dict=validating_pair_dict
        )

    reactive_pairs_set = set(reactive_pair_df.loc[:, "Pair"].unique())
    reactive_pairs_without_co_occurrence_filter_set = set(reactive_pair_not_co_occurrence_filtered_df.loc[:, "Pair"].unique())
    diff_reactive_pairs_with_and_without_co_occurrence_filter_set = reactive_pairs_without_co_occurrence_filter_set.difference(reactive_pairs_set)
    not_significant_pairs_set = set(filtered_pair_df.loc[~filtered_pair_df["BH_Reject"], "Pair"].unique())
    filtered_pair_with_co_occurrence_df = filtered_pair_df.loc[
        (filtered_pair_df["tcr_only_supported_by_co_occurrence"] | filtered_pair_df["epi_only_supported_by_co_occurrence"]), :
    ].copy()
    filtered_pair_not_significant_with_co_occurrence_df = filtered_pair_with_co_occurrence_df.loc[~filtered_pair_with_co_occurrence_df["BH_Reject"], :]
    not_significant_pairs_with_co_occurrence_set = set(filtered_pair_not_significant_with_co_occurrence_df.loc[:, "Pair"].unique())

    if validating_pair_dict is not None:
        missing_pairs_set = validating_pairs_set.difference(reactive_pairs_set)
        missing_pairs_not_significant_set = missing_pairs_set.intersection(not_significant_pairs_set)
        missing_pairs_significant_with_co_occurrence_set = missing_pairs_set.intersection(diff_reactive_pairs_with_and_without_co_occurrence_filter_set)
        missing_pairs_not_significant_with_co_occurrence_set = missing_pairs_set.intersection(not_significant_pairs_with_co_occurrence_set)
        missing_pairs_never_detected_set = validating_pairs_set.difference(total_detected_pairs_set)
        if missing_pairs_never_detected_set:
            missing_pairs_without_detected_epi_list = []
            missing_pairs_without_detected_tcr_list = []
            for pair in missing_pairs_never_detected_set:
                epi_tmp = pair.split("-")[0]
                tcr_tmp = pair.split("-")[1]
                if epi_tmp not in total_detected_epi_set:
                    missing_pairs_without_detected_epi_list.append(pair)
                if tcr_tmp not in total_detected_tcr_set:
                    missing_pairs_without_detected_tcr_list.append(pair)

        new_reactive_pairs_set = reactive_pairs_set.difference(validating_pairs_set)

    with open(plates_list_outs / f"reactive_pair_result_{plates_list_name}.txt", "w") as reactive_pair_result_out_txt:
        if validating_pair_dict is not None:
            reactive_pair_result_out_txt.write("Sensitivity: " + str(round(sensitivity, 2)) + "\n")
            reactive_pair_result_out_txt.write("Sensitivity without co-occurrence filtering: " + str(round(sensitivity_without_co_occurrence_filter, 2)) + "\n")
            reactive_pair_result_out_txt.write("Precision: " + str(round(precision, 2)) + "\n\n")

        reactive_pair_result_out_txt.write("Identified reactive pairs (sorted in from most-to-least significant):\n")
        for i, pair in enumerate(reactive_pair_df["Pair"]):
            reactive_pair_result_out_txt.write(str(i + 1) + ". " + pair + "\n")
        reactive_pair_result_out_txt.write("\n")

        if validating_pair_dict is not None:
            reactive_pair_result_out_txt.write("Missing previous reactive pairs:\n")
            for i, pair in enumerate(missing_pairs_set):
                reactive_pair_result_out_txt.write(str(i + 1) + ". " + pair + "\n")
            reactive_pair_result_out_txt.write("\n")

            reactive_pair_result_out_txt.write("Missing previous reactive pairs that are not significant:\n")
            for i, pair in enumerate(missing_pairs_not_significant_set):
                reactive_pair_result_out_txt.write(str(i + 1) + ". " + pair + "\n")
            reactive_pair_result_out_txt.write("\n")

            reactive_pair_result_out_txt.write("Missing previous reactive pairs that are significant but were co-occurrence filtered:\n")
            for i, pair in enumerate(missing_pairs_significant_with_co_occurrence_set):
                reactive_pair_result_out_txt.write(str(i + 1) + ". " + pair + "\n")
            reactive_pair_result_out_txt.write("\n")

            reactive_pair_result_out_txt.write("Missing previous reactive pairs that are not significant and have co-occurrence:\n")
            for i, pair in enumerate(missing_pairs_not_significant_with_co_occurrence_set):
                reactive_pair_result_out_txt.write(str(i + 1) + ". " + pair + "\n")
            reactive_pair_result_out_txt.write("\n")

            reactive_pair_result_out_txt.write("Missing previous reactive pairs that were never detected:\n")
            for i, pair in enumerate(missing_pairs_never_detected_set):
                reactive_pair_result_out_txt.write(str(i + 1) + ". " + pair + "\n")
            reactive_pair_result_out_txt.write("\n")

            if missing_pairs_never_detected_set:
                reactive_pair_result_out_txt.write("Missing previous reactive pairs that were never detected with a never detected Epitope:\n")
                for i, pair in enumerate(missing_pairs_without_detected_epi_list):
                    reactive_pair_result_out_txt.write(str(i + 1) + ". " + pair + "\n")
                reactive_pair_result_out_txt.write("\n")

                reactive_pair_result_out_txt.write("Missing previous reactive pairs that were never detected with a never detected TCR:\n")
                for i, pair in enumerate(missing_pairs_without_detected_tcr_list):
                    reactive_pair_result_out_txt.write(str(i + 1) + ". " + pair + "\n")
                reactive_pair_result_out_txt.write("\n")

            reactive_pair_result_out_txt.write("Newly identified reactive pairs:\n")
            for i, pair in enumerate(new_reactive_pairs_set):
                reactive_pair_result_out_txt.write(str(i + 1) + ". " + pair + "\n")
            reactive_pair_result_out_txt.write("\n")

        reactive_pair_result_out_txt.write(
            "Identified significant non-co-occurring pairs with cross-detected TCR (sorted from most-to-least significant; only the most significant pair of each unique TCR was kept as identified reactive pair):\n"
        )
        for i, (pair, p_val_corrected) in enumerate(zip(cross_detected_df["Pair"], cross_detected_df["P_value_corrected"])):
            reactive_pair_result_out_txt.write(str(i + 1) + ". " + pair + ", " + str(p_val_corrected) + "\n")
        reactive_pair_result_out_txt.write("\n")


def _make_screen_comparison_scatter(
    merged_pair_fold_changes_df: pd.DataFrame,
    screen_1: str,
    screen_2: str,
    adj_pval_alpha: float,
    zoom_adj_pval_alpha: float,
    zoom_log_fold: float,
    w: float,
    h: float,
    output_path: Path,
    *,
    annotate: bool = False,
    zoom: bool = False,
    validation_coloring: bool = True,
    only_plot_validating: bool = False,
):
    output_path = Path(output_path)
    merged_pair_fold_changes_df = merged_pair_fold_changes_df.copy()

    if only_plot_validating:
        mask = pd.to_numeric(merged_pair_fold_changes_df["validating_pair"], errors="coerce") == 1
        merged_pair_fold_changes_df = merged_pair_fold_changes_df.loc[mask, :]
    x = merged_pair_fold_changes_df["P_value_corrected_" + screen_1]
    y = merged_pair_fold_changes_df["P_value_corrected_" + screen_2]

    sres = spearmanr(x, y)
    ci_low, ci_high = bootstrap_spearman_ci(x, y)

    ax, fig, gs = startfig(w, h)
    log_pval_1 = -np.log10(merged_pair_fold_changes_df["P_value_corrected_" + screen_1])
    log_pval_2 = -np.log10(merged_pair_fold_changes_df["P_value_corrected_" + screen_2])
    log_zoom = -np.log10(zoom_adj_pval_alpha)

    color_list = []
    text_list = []

    for idx, (pval_1, pval_2, pair, reject, validated, rej1, rej2) in enumerate(
        zip(
            log_pval_1,
            log_pval_2,
            merged_pair_fold_changes_df["Pair"],
            merged_pair_fold_changes_df["P_value_corrected_combined"] < adj_pval_alpha,
            merged_pair_fold_changes_df["validating_pair"],
            merged_pair_fold_changes_df["P_value_corrected_" + screen_1] < adj_pval_alpha,
            merged_pair_fold_changes_df["P_value_corrected_" + screen_2] < adj_pval_alpha,
        )
    ):
        reject_any = rej1 or rej2

        if validation_coloring:
            if validated is True:
                color = "#97C654"
            elif validated is False:
                color = "#f7941d"
            elif pd.isna(validated) and reject:
                color = "#5C6C46"
            elif reject_any:
                color = "lightblue"
            else:
                color = "lightgrey"
        else:
            if reject:
                color = "#97C654"
            elif reject_any:
                color = "lightblue"
            else:
                color = "lightgrey"

        color_list.append(color)
        if annotate and color != "lightgrey":
            text_list.append(ax.text(pval_1, pval_2, pair, fontsize=2, zorder=8))

    ax.scatter(log_pval_1, log_pval_2, color=color_list, s=3)
    ax.set_xlabel(f"-Log10(adjusted P value {screen_1})", fontsize=7)
    ax.set_ylabel(f"-Log10(adjusted P value {screen_2})", fontsize=7)
    ax.set_title((f"R={sres.statistic:.2f} [{ci_low:.2f}-{ci_high:.2f}], P={sres.pvalue:.2e}"), fontsize=7)
    ax.tick_params(top=False, right=False)
    ax.tick_params("both", labelsize=7)

    if zoom:
        ax.set_xlim(-0.75, log_zoom)
        ax.set_ylim(-0.75, log_zoom)
    else:
        max_lim = max(log_pval_1.max(), log_pval_2.max()) * 1.05
        ax.set_xlim(-5, max_lim)
        ax.set_ylim(-5, max_lim)

    if annotate and text_list:
        adjust_text(text_list, ax=ax, arrowprops=dict(arrowstyle="-", color="k", zorder=4, lw=0.25))

    fig.tight_layout()
    fig.savefig(output_path.with_name(output_path.stem + "_adj_pval" + output_path.suffix))
    plt.close()

    merged_pair_fold_changes_df = merged_pair_fold_changes_df.copy()
    merged_pair_fold_changes_df = merged_pair_fold_changes_df.loc[merged_pair_fold_changes_df["P_value_corrected_combined"] < adj_pval_alpha, :]

    ax, fig, gs = startfig(w, h)
    fold_change_1 = merged_pair_fold_changes_df["Fold_change_" + screen_1]
    fold_change_2 = merged_pair_fold_changes_df["Fold_change_" + screen_2]

    sres = spearmanr(fold_change_1, fold_change_2)
    ci_low, ci_high = bootstrap_spearman_ci(fold_change_1, fold_change_2)

    color_list = []
    text_list = []

    for idx, (fold_1, fold_2, pair, validated) in enumerate(zip(fold_change_1, fold_change_2, merged_pair_fold_changes_df["Pair"], merged_pair_fold_changes_df["validating_pair"])):
        if validation_coloring:
            if validated is True:
                color = "#97C654"
            elif validated is False:
                color = "#f7941d"
            elif pd.isna(validated):
                color = "#5C6C46"
            else:
                color = "lightgrey"
        else:
            color = "#97C654"

        color_list.append(color)
        if annotate and color != "lightgrey":
            text_list.append(ax.text(fold_1, fold_2, pair, fontsize=2, zorder=8))

    ax.scatter(fold_change_1, fold_change_2, color=color_list, s=3)
    ax.set_xlabel(f"Fold change {screen_1}", fontsize=7)
    ax.set_ylabel(f"Fold change {screen_2}", fontsize=7)
    ax.set_title((f"R={sres.statistic:.2f} [{ci_low:.2f}-{ci_high:.2f}], P={sres.pvalue:.2e}"), fontsize=7)
    ax.tick_params(top=False, right=False)
    ax.tick_params("both", labelsize=7)

    if zoom:
        ax.set_xlim(-0.75, zoom_log_fold)
        ax.set_ylim(-0.75, zoom_log_fold)
    else:
        max_lim = max(fold_change_1.max(), fold_change_2.max()) * 1.05
        ax.set_xlim(-(max_lim * 0.03), max_lim)
        ax.set_ylim(-(max_lim * 0.03), max_lim)

    if annotate and text_list:
        adjust_text(text_list, ax=ax, arrowprops=dict(arrowstyle="-", color="k", zorder=4, lw=0.25))

    fig.tight_layout()
    fig.savefig(output_path.with_name(output_path.stem + "_fold_change" + output_path.suffix))
    plt.close()


def plot_and_combine_screen_pvals_stouffer(
    screen_output_run_dir_dict: dict,
    adj_pval_alpha: float,
    project_dir: str | os.PathLike,
    zoom_adj_pval_alpha: float,
    zoom_log_fold: float,
    validating_pair_dict: dict = None,
    subset_epi_list: list = None,
    only_plot_validating: bool = False,
    w: float = 5.5,
    h: float = 5.5,
):
    """
    Combine BH-adjusted P values for detected TCR-epitope pairs across multiple independent screens using
    Stouffer's weighted Z-score method, plot pairwise screen comparisons, and write the combined pair hit set.

    For each screen in screen_output_run_dir_dict, this function reads its pair_fold_changes_df_<run>.xlsx and
    pair_count_matrix_df_<run>.csv, drops pairs only supported by co-occurrence, keeps only the most significant
    pair per TCR, and (optionally) subsets to a list of epitopes. The screens' P_value_corrected and Fold_change
    columns are then outer-merged on "Pair", and each screen's adjusted P values are combined per pair with
    scipy.stats.combine_pvalues(method="stouffer"), weighted by that screen's total called pair count, into a
    single P_value_corrected_combined.

    Output files: for every pair of screens, six scatter plots (annotated/non-annotated, zoomed/non-zoomed,
    validation-colored/significance-only-colored) comparing -log10(adjusted P value) and, for jointly significant
    pairs, fold change, are written to outs/<screen_1>_vs_..._stouffer[_<epitopes>], along with their Spearman
    correlation and bootstrapped confidence interval. The combined, deduplicated (one pair per TCR) pair hit set is
    then written via :func:`write_reactive_pair_results`.

    Parameters
    ----------
    screen_output_run_dir_dict : dict
        Dictionary mapping a screen/run name to the string path of its outs/<plates_list_name> output directory,
        containing pair_fold_changes_df_<dir_name>.xlsx and pair_count_matrix_df_<dir_name>.csv.
    adj_pval_alpha : float
        Significance threshold applied to P_value_corrected_combined (and to each screen's own adjusted P value) to
        determine BH_Reject/coloring in the comparison plots.
    project_dir : str | os.PathLike
        String path to the plate-based TCR and epitope sequencing project directory. The combined-screen output
        directory outs/<screen_1>_vs_..._stouffer[_<epitopes>] is created under this path.
    zoom_adj_pval_alpha : float
        Adjusted P value used to set the axis limits of the "zoomed" -log10(P value) comparison plots (limits are
        set to [-0.75, -log10(zoom_adj_pval_alpha)]).
    zoom_log_fold : float
        Fold-change axis limit used for the "zoomed" fold-change comparison plots.
    validating_pair_dict : dict, optional
        Dictionary where keys are epitope-TCR pair strings ("<epitope_reference_name>-<tcr_reference_name>") and
        values are booleans indicating whether that pair was confirmed reactive (True) or non-reactive (False) in
        an independent 1-to-1 validation assay. Must pass :func:`validate_validating_pair_dict`. Passed through to
        :func:`write_reactive_pair_results` for the combined result.
    subset_epi_list : list, optional
        List of epitope reference names to subset every screen's pairs to before combining. If provided, the
        combined output directory and result files are suffixed with the joined epitope names.
    only_plot_validating : bool, default = False
        If True, restrict the pairwise screen comparison scatter plots to only pairs with a truthy
        "validating_pair" value.
    w : float, default = 5.5
        Figure width used for the pairwise screen comparison scatter plots.
    h : float, default = 5.5
        Figure height used for the pairwise screen comparison scatter plots.

    Returns
    -------
    None
        Writes pairwise screen comparison .pdf plots and merged_pair_fold_changes_df.xlsx to
        outs/<screen_1>_vs_..._stouffer[_<epitopes>], and writes the combined pair hit set via
        :func:`write_reactive_pair_results`.
    """
    screen_name_list = list(screen_output_run_dir_dict.keys())

    pair_fold_changes_df_dict = {}
    pair_count_matrix_dict = {}
    for screen in screen_name_list:
        pair_fold_changes_df_dict[screen] = pd.read_excel(
            Path(screen_output_run_dir_dict[screen]) / f"pair_fold_changes_df_{Path(screen_output_run_dir_dict[screen]).name}.xlsx", index_col=0
        )
        pair_fold_changes_df_dict[screen]["validating_pair"] = pair_fold_changes_df_dict[screen]["validating_pair"].map({1.0: True, 0.0: False})

        pair_fold_changes_df_dict[screen] = (
            pair_fold_changes_df_dict[screen]
            .loc[~(pair_fold_changes_df_dict[screen]["tcr_only_supported_by_co_occurrence"] | pair_fold_changes_df_dict[screen]["epi_only_supported_by_co_occurrence"]), :]
            .copy()
        )
        pair_fold_changes_df_dict[screen].drop_duplicates("TCR", keep="first", inplace=True)
        if subset_epi_list is not None:
            pair_fold_changes_df_dict[screen] = pair_fold_changes_df_dict[screen].loc[pair_fold_changes_df_dict[screen].loc[:, "Epitope"].isin(subset_epi_list), :].copy()
        pair_count_matrix_dict[screen] = pd.read_csv(
            Path(screen_output_run_dir_dict[screen]) / f"pair_count_matrix_df_{Path(screen_output_run_dir_dict[screen]).name}.csv", index_col=0
        )

    merged_pair_fold_changes_df = pair_fold_changes_df_dict[screen_name_list[0]]

    keep_columns_list = ["Pair", "P_value_corrected", "validating_pair", "Fold_change"]

    merged_pair_fold_changes_df = merged_pair_fold_changes_df.loc[:, keep_columns_list]

    merged_pair_fold_changes_df.rename(dict(zip(keep_columns_list[1:], [column + "_" + screen_name_list[0] for column in keep_columns_list[1:]])), axis=1, inplace=True)

    for screen in screen_name_list[1:]:
        pair_fold_changes_df_dict[screen] = pair_fold_changes_df_dict[screen].loc[:, keep_columns_list]
        pair_fold_changes_df_dict[screen].rename(dict(zip(keep_columns_list[1:], [column + "_" + screen for column in keep_columns_list[1:]])), axis=1, inplace=True)
        merged_pair_fold_changes_df = merged_pair_fold_changes_df.merge(pair_fold_changes_df_dict[screen], on="Pair", how="outer", suffixes=("", ""))

    merged_pair_fold_changes_df.reset_index(drop=True, inplace=True)

    def get_validating_status(row):
        non_nan = row.dropna()
        if non_nan.empty:
            return np.nan
        if non_nan.nunique() > 1:
            raise ValueError(f"Conflicting validating_pair values across screens for pair: {row.name}")
        return non_nan.iloc[0]

    merged_pair_fold_changes_df["validating_pair"] = merged_pair_fold_changes_df.loc[:, ["validating_pair_" + screen for screen in screen_name_list]].apply(
        get_validating_status, axis=1
    )
    merged_pair_fold_changes_df.drop(["validating_pair" + "_" + screen for screen in screen_name_list], axis=1, inplace=True)

    weights_dict = {}
    for screen in screen_name_list:
        merged_pair_fold_changes_df["P_value_corrected_" + screen] = merged_pair_fold_changes_df["P_value_corrected_" + screen].fillna(1.0)
        if (merged_pair_fold_changes_df["P_value_corrected_" + screen] == 0).any():
            pseudo_count = merged_pair_fold_changes_df.loc[merged_pair_fold_changes_df.loc[:, "P_value_corrected_" + screen] > 0, "P_value_corrected_" + screen].min() * 0.1
            merged_pair_fold_changes_df["P_value_corrected_" + screen] += pseudo_count

        merged_pair_fold_changes_df["Fold_change_" + screen] = merged_pair_fold_changes_df["Fold_change_" + screen].fillna(0)
        weights_dict[screen] = pair_count_matrix_dict[screen].sum().sum()

    combined_pval_list = []
    for idx in merged_pair_fold_changes_df.index:
        pvals_corrected = list(merged_pair_fold_changes_df.loc[idx, ["P_value_corrected_" + screen for screen in screen_name_list]].values)
        combined_pval_list.append(sp.stats.combine_pvalues(pvals_corrected, method="stouffer", weights=[weights_dict[screen] for screen in screen_name_list])[1])

    merged_pair_fold_changes_df["P_value_corrected_combined"] = combined_pval_list

    output_dir = Path(project_dir) / "outs"
    comparison_name = "_vs_".join(screen_name_list) + "_stouffer"
    if subset_epi_list is not None:
        comparison_name += "_" + "-".join(subset_epi_list)

    screen_comparison_path = output_dir / comparison_name
    screen_comparison_path.mkdir(exist_ok=False)

    for screen_1, screen_2 in itertools.combinations(screen_name_list, 2):
        base = str(screen_comparison_path / f"{screen_1}_vs_{screen_2}")
        plot_kwargs = dict(
            merged_pair_fold_changes_df=merged_pair_fold_changes_df,
            screen_1=screen_1,
            screen_2=screen_2,
            adj_pval_alpha=adj_pval_alpha,
            zoom_adj_pval_alpha=zoom_adj_pval_alpha,
            zoom_log_fold=zoom_log_fold,
            w=w,
            h=h,
        )
        _make_screen_comparison_scatter(
            **plot_kwargs, annotate=True, zoom=False, validation_coloring=True, only_plot_validating=only_plot_validating, output_path=base + "_annotated.pdf"
        )
        _make_screen_comparison_scatter(
            **plot_kwargs, annotate=True, zoom=True, validation_coloring=True, only_plot_validating=only_plot_validating, output_path=base + "_zoom_annotated.pdf"
        )
        _make_screen_comparison_scatter(**plot_kwargs, annotate=False, zoom=False, validation_coloring=True, only_plot_validating=only_plot_validating, output_path=base + ".pdf")
        _make_screen_comparison_scatter(
            **plot_kwargs, annotate=False, zoom=True, validation_coloring=True, only_plot_validating=only_plot_validating, output_path=base + "_zoom.pdf"
        )
        _make_screen_comparison_scatter(
            **plot_kwargs, annotate=False, zoom=False, validation_coloring=False, only_plot_validating=only_plot_validating, output_path=base + "_only_significance_colored.pdf"
        )
        _make_screen_comparison_scatter(
            **plot_kwargs, annotate=False, zoom=True, validation_coloring=False, only_plot_validating=only_plot_validating, output_path=base + "_zoom_only_significance_colored.pdf"
        )

    merged_pair_fold_changes_df.drop(["P_value_corrected_" + screen for screen in screen_name_list], axis=1, inplace=True)
    merged_pair_fold_changes_df.rename({"P_value_corrected_combined": "P_value_corrected"}, axis=1, inplace=True)
    merged_pair_fold_changes_df.sort_values("P_value_corrected", ascending=True, inplace=True)
    merged_pair_fold_changes_df["Epitope"] = merged_pair_fold_changes_df["Pair"].str.split("-").str[0]
    merged_pair_fold_changes_df["TCR"] = merged_pair_fold_changes_df["Pair"].str.split("-").str[1]
    merged_pair_fold_changes_df.drop_duplicates("TCR", keep="first", inplace=True)
    merged_pair_fold_changes_df.reset_index(drop=True, inplace=True)
    merged_pair_fold_changes_df["tcr_only_supported_by_co_occurrence"] = False
    merged_pair_fold_changes_df["epi_only_supported_by_co_occurrence"] = False
    merged_pair_fold_changes_df["TCR_significant_with_more_than_one_epi"] = False
    merged_pair_fold_changes_df["epi_significant_with_more_than_one_TCR"] = False
    merged_pair_fold_changes_df["BH_Reject"] = merged_pair_fold_changes_df["P_value_corrected"] < adj_pval_alpha
    merged_pair_fold_changes_df.to_excel(screen_comparison_path / "merged_pair_fold_changes_df.xlsx")

    if subset_epi_list is None:
        write_reactive_pair_results(
            pair_fold_changes_df=merged_pair_fold_changes_df,
            project_dir=project_dir,
            plates_list_name="_vs_".join(screen_name_list) + "_stouffer",
            validating_pair_dict=validating_pair_dict,
        )
    else:
        write_reactive_pair_results(
            pair_fold_changes_df=merged_pair_fold_changes_df,
            project_dir=project_dir,
            plates_list_name="_vs_".join(screen_name_list) + "_stouffer_" + "-".join(subset_epi_list),
            validating_pair_dict=validating_pair_dict,
        )
