import collections
import itertools
import os
import io
import random
from typing import Literal, Union

import numpy as np
import pandas as pd
from dotenv import load_dotenv

from tcr_toolbox.utils.constants import aa_to_codon

load_dotenv()
tcr_toolbox_data_path = os.getenv("tcr_toolbox_data_path")


def reverse_translate(seq):
    """Reverse translate a protein sequence into a DNA sequence using the most common nucleotide 3-mers."""
    seq = seq.upper()
    return "".join([aa_to_codon[seq[i]] for i in range(0, len(seq))])


def generate_random_dna_sequence(length):
    return "".join(random.choice("ACTG") for _ in range(length))


def interactive_well_from_df(df, well_384_or_96, hover_columns, out_path, color_column):
    """Create an interactive plotly figure from a dataframe with well names as index and columns as interactive columns.

    Args:
        df (pandas df): pandas dataframe with well names as index and columns as interactive columns
        hover_columns (list): list of column names to be interactive inside the wells
        out_path (str): path to save the html file (e.g. "my_plot.html")
        color_column (str): color the wells based on this column
    Returns:
        fig: plotly figure"""
    import plotly.express as px

    df = add_wells_plate_layout(plate_df=df, well_384_or_96=well_384_or_96)

    # Create a scatter plot with custom data and layout
    fig = px.scatter(df, x="column", y="row", color=color_column, hover_data=hover_columns, width=1200, height=600)

    # Customize the axes labels and ranges
    fig.update_xaxes(title="Column", tickmode="array", tickvals=list(range(1, 25)), range=[0, 25])
    fig.update_yaxes(title="Row", tickmode="array", tickvals=list(range(0, 17)), range=[0, 17])
    fig.update_traces(marker=dict(size=12))
    fig.update_yaxes(autorange="reversed")  # reverse x axis
    fig.write_html(out_path)


def levenshtein_distance(seq1, seq2):
    """
    Calculates the levenshtein distance between two sequences.
    """
    # initialize matrix
    matrix = [[0 for x in range(len(seq2) + 1)] for x in range(len(seq1) + 1)]
    # initialize first row
    for i in range(len(seq1) + 1):
        matrix[i][0] = i
    # initialize first column
    for j in range(len(seq2) + 1):
        matrix[0][j] = j
    # fill in matrix
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            # if the characters are the same, there is no cost
            if seq1[i - 1] == seq2[j - 1]:
                cost = 0
            # if the characters are not the same, there is a cost of 1
            else:
                cost = 1
            # calculate the minimum cost
            matrix[i][j] = min(
                matrix[i - 1][j] + 1,  # deletion
                matrix[i][j - 1] + 1,  # insertion
                matrix[i - 1][j - 1] + cost,
            )  # substitution
    # return the bottom right element of the matrix
    return matrix[-1][-1]


def levenshtein_ratio(seq1, seq2):
    max_len = max(len(seq1), len(seq2))
    return (max_len - levenshtein_distance(seq1, seq2)) / max_len


def add_well_coordinates(plate_df):
    if "well" not in plate_df.columns:
        raise "Please add well column to add well coordinates"

    plate_df["well_coord"] = plate_df["well"].apply(lambda x: (ord(x[0]) - 64, int(x[1:])))
    plate_df["row"] = plate_df["well_coord"].apply(lambda x: chr(x[0] + 64))
    plate_df["column"] = plate_df["well_coord"].apply(lambda x: x[1])

    return plate_df


def add_wells_plate_layout(plate_df, well_384_or_96=96):
    # make 96 well layout (A1 to A12 and B1 to B12 etc, untill the end of the plate)
    if well_384_or_96 == 96:
        wells = [x + str(y) for x in "ABCDEFGH" for y in range(1, 13)]
    elif well_384_or_96 == 384:
        wells = [x + str(y) for x in "ABCDEFGHIJKLMNOP" for y in range(1, 25)]
    else:
        raise ValueError("well_384_or_96 must be either 96 or 384")

    # add well column
    plate_df["well"] = wells[: plate_df.shape[0]]

    # if no entries: make empty columns
    if plate_df.shape[0] < 96:
        for i in range(plate_df.shape[0], 96):
            plate_df.loc[i] = [np.nan for x in range(plate_df.shape[1])]
            plate_df.loc[i, "well"] = wells[i]

    # Create a column for the well coordinates
    plate_df = add_well_coordinates(plate_df)

    return plate_df


def write_hamilton_pipetting_excel_sheet(
    target_plates_list: list,
    source_df_plates_col_name: str,
    target_df_plates_col_name: str,
    source_df: pd.DataFrame,
    target_df: pd.DataFrame,
    source_match_col_name: str,
    target_match_col_name: str,
    source_df_well_col_name: str,
    target_df_well_col_name: str,
    max_source_plate_volume_ul: float,
    excel_write_path: Union[str, os.PathLike[str]],
    target_source_subset_iter_list: list = None,
    target_source_subset_iter_col_name: str = "",
    target_annotation_col_names_list: list = None,
    volume_col_name: str = "",
    transfer_volume_ul: float = 0.00,
    volume_solvent_col_name: str = "",
    volume_solvent_ul: float = 0.00,
    target_plate_format: int = Literal[96, 384],
    also_write_csv: bool = False,
):
    """
    Write a Hamilton pipetting instruction .xlsx sheet based on source and target plate DataFrames.

    This function generates a detailed pipetting instruction table for Hamilton STAR robots. It
    matches samples between source and target plates, ensures volume constraints, handles
    multiple source subsets, and optionally adds annotation columns or solvent volumes. The
    resulting pipetting instructions are written to an Excel .xlsx file and optionally a CSV.

    Parameters
    ----------
    target_plates_list : list
        List of target plate names that should be pipetted into.
    source_df_plates_col_name : str
        Column name in the source DataFrame containing source plate names.
    target_df_plates_col_name : str
        Column name in the target DataFrame containing target plate names.
    source_df : pd.DataFrame
        Source DataFrame storing sample contents of source plate wells.
    target_df : pd.DataFrame
        Target DataFrame specifying samples to be pipetted into each target well.
    source_match_col_name : str
        Column in source DataFrame containing sample names for matching wells between source and target.
    target_match_col_name : str
        Column in target DataFrame containing sample names for matching wells between source and target.
    source_df_well_col_name : str
        Column in source DataFrame containing well identifiers.
    target_df_well_col_name : str
        Column in target DataFrame containing well identifiers.
    max_source_plate_volume_ul : float
        Maximum allowed volume in a source well (µL).
    excel_write_path : Union[str, os.PathLike[str]]
        Path to write the Hamilton .xlsx instruction sheet.
    target_source_subset_iter_list : list, optional
        Required when multiple different samples from different source subsets are pipetted into the same target wells.
        Indicates the order in which source subsets are used. Example: ['Fw', 'Rev'].
    target_source_subset_iter_col_name : str, optional
        Column in the target DataFrame specifying the source subset for each sample. Names must match those in
        `target_source_subset_iter_list`.
    target_annotation_col_names_list : list, optional
        List of column names from the target DataFrame to include as annotations in the pipetting sheet.
    volume_col_name : str, optional
        Column in the target DataFrame specifying variable transfer volumes for each well.
    transfer_volume_ul : float, optional
        Constant transfer volume (µL) for all wells if `volume_col_name` is not provided.
    volume_solvent_col_name : str, optional
        Column in the target DataFrame specifying solvent volume for each well.
    volume_solvent_ul : float, optional
        Constant solvent volume (µL) for all wells if `volume_solvent_col_name` is not provided.
    target_plate_format : int, optional
        Format of target plate: 96 or 384 wells.
    also_write_csv : bool, default False
        If True, also write a CSV file alongside the Excel file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the Hamilton pipetting instructions that were written to Excel.

    Writes
    ------
    Hamilton pipetting instruction .xlsx file (and optionally CSV).

    Notes
    -----
    - Usage example can be found in `make_idot_v_gene_premixing_dispense_csv`.
    - Unique samples are not allowed to be stored in multiple source plates. For example, the same
      primer (e.g., "Fw PacBio257") cannot be stored in multiple wells across different source plates.
    - Ensures correct pipetting order for 384-well plates by separating even and odd rows to allow
      simultaneous pipetting using 8-channel Hamilton STAR pipettes.
    - Handles source well exhaustion by checking `max_source_plate_volume_ul` and switching to
      additional wells if needed.

    Raises
    ------
    Exception
        - If a sample is found in multiple source plates.
        - If no matching source plate or well is found for a target sample.
        - If there are insufficient source wells to transfer the requested volume.
        - If lengths of target subsets iterator and target plate list are out of sync.
    """

    # Needed because otherwise the dfs in the notebook also get re-ordered below
    source_df = source_df.copy()
    target_df = target_df.copy()
    letters = "ABCDEFGHIJKLMNOP"
    col_numbers = np.arange(1, 25)

    # remove leading 0 if present in well coordinates
    source_df[source_df_well_col_name] = remove_leading_zeros_well_coordinates(source_df[source_df_well_col_name])
    target_df[target_df_well_col_name] = remove_leading_zeros_well_coordinates(target_df[target_df_well_col_name])

    if target_plate_format == 384:
        # Eight Hamilton STAR 96 channels can only simultaneously pipette into even or uneven rows of the same column
        # Re-order target_df to allow for simultaneous even and uneven row pipetting in the same columns:
        even_letters = letters[::2]
        uneven_letters = letters[1::2]
        reorder_target_by_well_list = [letter + str(col_number) for col_number, letter in itertools.product(col_numbers, even_letters + uneven_letters)]
    elif target_plate_format == 96:
        reorder_target_by_well_list = []
        for col_number in col_numbers:
            for letter in letters:
                reorder_target_by_well_list.append(letter + str(col_number))
    else:
        raise ValueError("target_plate_format must be either 96 or 384")

    sorter_index = dict(zip(reorder_target_by_well_list, range(len(reorder_target_by_well_list))))
    target_df["well_rank"] = target_df[target_df_well_col_name].map(sorter_index)

    if not target_source_subset_iter_col_name:
        target_df.sort_values(["well_rank"], inplace=True)
    else:
        target_df.sort_values([target_source_subset_iter_col_name, "well_rank"], inplace=True)

    # Initialize the columns of the Excel pipetting sheet in a dispense_dict:
    dispense_dict = {"source_ID": [], "source_well_ID": [], "target_ID": [], "target_well_ID": [], "Vol": [], "Vol_solvent": [], "matcher_ID": []}

    # For each pipetting instruction, add annotation columns from target_annotation_col_names_list to dispense_dict:
    if target_annotation_col_names_list:
        for annotation_col in target_annotation_col_names_list:
            dispense_dict[annotation_col] = []

    # Initialize the iterator lists for generating the pipetting instruction sheet:
    target_df_subsets_iterator_list = []
    target_plate_names_zip_list = []

    if not target_source_subset_iter_list:
        for target_plate in target_plates_list:
            target_df_subsets_iterator_list.append(target_df[target_df[target_df_plates_col_name] == target_plate].copy())
        target_plate_names_zip_list = target_plates_list

    else:
        for target_plate in target_plates_list:
            for source_subset in target_source_subset_iter_list:
                # If we do not pipette from the current source subset at all into the current target plate,
                # continue to next source subset:
                if target_df[(target_df[target_df_plates_col_name] == target_plate) & (target_df[target_source_subset_iter_col_name] == source_subset)].shape[0] == 0:
                    continue

                target_df_subsets_iterator_list.append(
                    target_df[(target_df[target_df_plates_col_name] == target_plate) & (target_df[target_source_subset_iter_col_name] == source_subset)].copy()
                )
                target_plate_names_zip_list.append(target_plate)

    if len(target_df_subsets_iterator_list) != len(target_plate_names_zip_list):
        raise Exception("Lengths of target_df_subset_iterator_list and target_plate_names_zip_list are out of sync!")

    source_well_counter = collections.defaultdict(int)
    for tmp_target_df, target_plate in zip(target_df_subsets_iterator_list, target_plate_names_zip_list):
        for well in tmp_target_df[target_df_well_col_name]:
            # sample_matcher stores the name of the Sample that needs to be pipetted into the current target well,
            # and this Sample name is used to check whether there is a matching source well that contains this Sample:
            sample_matcher = tmp_target_df[tmp_target_df[target_df_well_col_name] == well][target_match_col_name].values[0]

            source_plate_list = source_df[source_df[source_match_col_name] == sample_matcher][source_df_plates_col_name].unique().tolist()

            if len(source_plate_list) > 1:
                raise Exception("Current Sample", sample_matcher, "is stored in source wells in multiple source plates! This is not allowed!")
            elif len(source_plate_list) == 0:
                raise Exception(
                    "No matching source plates found for this Sample!",
                    sample_matcher,
                    well,
                    target_plate,
                    "Are you sure Sample names have the same format in source and target df?",
                )
            source_plate = source_plate_list[0]

            matching_source_well_list = source_df[source_df[source_match_col_name] == sample_matcher][source_df_well_col_name].to_list()

            if not matching_source_well_list:
                raise Exception("No matching source wells found for this Sample! Are you sure Sample names have the same format in source and target df?")

            for idx, source_well in enumerate(matching_source_well_list):
                if source_well_counter[source_well] * transfer_volume_ul >= max_source_plate_volume_ul:
                    # If there are no source wells left for a specific Sample, more source wells should be prepared:
                    if len(matching_source_well_list[idx:]) == 1:
                        raise Exception(
                            "Not enough enough source wells for transferring required amount of µL!",
                            "Target well:",
                            well,
                            "Sample:",
                            sample_matcher,
                            "Last emptied Source well:",
                            source_well,
                        )
                    # If there are still other source wells left that also contain the Sample,
                    # pipette from next available source well:
                    continue
                else:
                    source_well_counter[source_well] += 1
                    dispense_dict["source_well_ID"].append(source_well)
                    dispense_dict["matcher_ID"].append(sample_matcher)
                    break  # exit the for loop if a well was used for pipetting

            # Only if we could pipette from a source well, we want to append plates and target well names to...
            # dispense dict value lists:
            dispense_dict["target_ID"].append(target_plate)
            dispense_dict["source_ID"].append(source_plate)
            dispense_dict["target_well_ID"].append(well)

            if not volume_col_name:
                dispense_dict["Vol"].append(transfer_volume_ul)
            else:
                dispense_dict["Vol"].append(tmp_target_df[tmp_target_df[target_df_well_col_name] == well][volume_col_name].values[0])

            if not volume_solvent_col_name:
                dispense_dict["Vol_solvent"].append(volume_solvent_ul)
            else:
                dispense_dict["Vol_solvent"].append(tmp_target_df[tmp_target_df[target_df_well_col_name] == well][volume_solvent_col_name].values[0])

            # Annotate pipetting instructions with user-defined target plate annotation columns:
            if target_annotation_col_names_list:
                for target_annotation_col in target_annotation_col_names_list:
                    dispense_dict[target_annotation_col].append(tmp_target_df[tmp_target_df[target_df_well_col_name] == well][target_annotation_col].values[0])

    # We do not have to check for NaNs in rows because all value lists in the dispense dict need to be of equal length
    # otherwise pd.DataFrame.from_dict will return a ValueError:
    dispense_df = pd.DataFrame.from_dict(dispense_dict)
    dispense_df.to_excel(excel_write_path)
    if also_write_csv:
        dispense_df.to_csv(excel_write_path.split(".xlsx")[0] + ".csv")

    return dispense_df


def write_echo_dispense_csv(
    target_plates_list: list,
    source_df_plates_col_name: str,
    target_df_plates_col_name: str,
    source_df: pd.DataFrame,
    target_df: pd.DataFrame,
    source_match_col_name: str,
    target_match_col_name: str,
    source_df_well_col_name: str,
    target_df_well_col_name: str,
    csv_write_path: Union[str, os.PathLike[str]],
    target_source_subset_iter_list: list = None,
    target_source_subset_iter_col_name: str = "",
    volume_col_name: str = "",
    transfer_volume_nl: float = 0.00,  # echo needs nanoliter instructions
    max_source_plate_volume_nl: int = 50 * 1000,  # echo needs nanoliter instructions
):
    """Write an echo dispense .csv sheet:

    Parameters
    ----------
    target_plates_list : list
        List of target plates names that should be dispensed into.
    source_df_plates_col_name : str
        Column in the source DataFrame containing source plate names.
    target_df_plates_col_name : str
        Column in the target DataFrame containing target plate names.
    source_df : pd.DataFrame
        Source DataFrame that stores Sample well contents of source plate(s).
    target_df : pd.DataFrame
        Target DataFrame that stores what Samples should be dispensed into each well of target plate(s).
    source_match_col_name : str
        Column in source DataFrame that containing Sample names for matching wells between source and target DataFrame.
    target_match_col_name : str
        Column in the target DataFrame that containing Sample names for matching wells between source and target DataFrame.
    source_df_well_col_name : str
        Column in source DataFrame containing well names.
    target_df_well_col_name : str
        Column in target DataFrame containing well names.
    max_source_plate_volume_nl : float
        Maximum nanoliter volume allowed in a source plate well (nL).
    csv_write_path : Union[str, os.PathLike[str]]
        csv write File Path.
    target_source_subset_iter_list : list, optional
        Required when multiple different Samples from different source subsets are pipetted into the same target well(s).
        Indicates the order in which those source subsets are used to pipette target plate(s). For example,
        when there are two separate "Fw" and "Rev" primers source subsets, target_source_subset_iter_list = ['Fw', 'Rev']
        indicates that first the Fw and next the Rev primers are pipetted.
        Source subsets may or may not be equal to source plates.
        Names of source subsets in this list should exactly match the names in the target_source_subset_iter_col_name in
        the target DataFrame.
    target_source_subset_iter_col_name : str, optional
        Required when multiple different Samples from different source subsets are pipetted into the same target well(s).
        Column in the target DataFrame that specifies, for each Sample, in which of the source subsets that Sample is stored.
        Source subsets may or may not be equal to source plates.
        Names of source subsets in this column should exactly match the names in the target_source_subset_iter_list.
    volume_col_name : str, optional
        Column in the target DataFrame that stores how much µL transfer volume should be transferred to each target well.
    transfer_volume_nl : float, optional
        If constant nL transfer volume is transferred to each target well, transfer volume in nanoliters (nL) that should be transferred.

    Returns
    -------
    pd.DataFrame
        DataFrame that contains the dispense instruction .csv sheet that was written to the .csv file.

    Writes
    -------
    Echo dispense instruction .csv file.

    Notes
    -------
    Usage example can be found in tcr_toolbox.tcr_assembly.order_automation.make_echo_v_gene_premixing_dispense_csv.
    Unique Samples are not allowed to be stored in separate plates. For example, the same Fw primer PacBio257 is not
    allowed to be stored in multiple source wells in two different source plates.
    """

    # Needed because otherwise the dfs in the notebook also get re-ordered below
    source_df = source_df.copy()
    target_df = target_df.copy()

    letters = "ABCDEFGHIJKLMNOP"
    col_numbers = np.arange(1, 25)

    reorder_target_by_well_list = []
    for letter in letters:
        for col_number in col_numbers:
            reorder_target_by_well_list.append(letter + str(col_number))

    sorter_index = dict(zip(reorder_target_by_well_list, range(len(reorder_target_by_well_list))))
    target_df["well_rank"] = target_df[target_df_well_col_name].map(sorter_index)

    if not target_source_subset_iter_col_name:
        target_df.sort_values(["well_rank"], inplace=True)
    else:
        target_df.sort_values([target_source_subset_iter_col_name, "well_rank"], inplace=True)

    # Initialize the columns of the dispense csv in a dispense_dict:
    dispense_dict = {
        "Source Plate Name": [],  # source_ID
        "Source Well": [],  # source_well_ID
        "Destination Plate Name": [],  # target_ID
        "Destination Well": [],  # target_well_ID
        "Transfer Volume": [],  # Vol
        "Sample ID": [],
    }  # :TRAV_TRBV

    # Initialize the iterator lists for generating the dispensing instruction sheet:
    target_df_subsets_iterator_list = []
    target_plate_names_zip_list = []

    if not target_source_subset_iter_list:
        for target_plate in target_plates_list:
            target_df_subsets_iterator_list.append(target_df[target_df[target_df_plates_col_name] == target_plate].copy())
        target_plate_names_zip_list = target_plates_list

    else:
        for target_plate in target_plates_list:
            for source_subset in target_source_subset_iter_list:
                # If we do not pipette from the current source subset at all into the current target plate,
                # continue to next source subset:
                if target_df[(target_df[target_df_plates_col_name] == target_plate) & (target_df[target_source_subset_iter_col_name] == source_subset)].shape[0] == 0:
                    continue

                target_df_subsets_iterator_list.append(
                    target_df[(target_df[target_df_plates_col_name] == target_plate) & (target_df[target_source_subset_iter_col_name] == source_subset)].copy()
                )
                target_plate_names_zip_list.append(target_plate)

    if len(target_df_subsets_iterator_list) != len(target_plate_names_zip_list):
        raise Exception("Lengths of target_df_subset_iterator_list and target_plate_names_zip_list are out of sync!")

    source_well_counter = collections.defaultdict(int)
    for tmp_target_df, target_plate in zip(target_df_subsets_iterator_list, target_plate_names_zip_list):
        for well in tmp_target_df[target_df_well_col_name]:
            # sample_matcher stores the name of the Sample that needs to be dispensed into the current target well,
            # and this Sample name is used to check whether there is a matching source well that contains this Sample:
            sample_matcher = tmp_target_df[tmp_target_df[target_df_well_col_name] == well][target_match_col_name].values[0]

            source_plate_list = source_df[source_df[source_match_col_name] == sample_matcher][source_df_plates_col_name].unique().tolist()

            if len(source_plate_list) > 1:
                raise Exception("Current Sample", sample_matcher, "is stored in source wells divided over multiple source plates! This is not allowed!")
            source_plate = source_plate_list[0]

            matching_source_well_list = source_df[source_df[source_match_col_name] == sample_matcher][source_df_well_col_name].to_list()

            if not matching_source_well_list:
                raise Exception("No matching source wells found for this Sample! Are you sure Sample names have the same format in source and target df?")

            for idx, source_well in enumerate(matching_source_well_list):
                if source_well_counter[source_well] * transfer_volume_nl >= max_source_plate_volume_nl:
                    # If there are no source wells left for a specific Sample, more source wells should be prepared:
                    if len(matching_source_well_list[idx:]) == 1:
                        raise Exception(
                            "Not enough enough source wells for transferring required amount of nL!",
                            "Target well:",
                            well,
                            "Sample:",
                            sample_matcher,
                            "Last emptied Source well:",
                            source_well,
                        )
                    # If there are still other source wells left that also contain the Sample,
                    # dispense from next available source well:
                    continue
                else:
                    source_well_counter[source_well] += 1
                    dispense_dict["Source Well"].append(source_well)
                    dispense_dict["Sample ID"].append(sample_matcher)
                    break  # exit the for loop if a well was used for dispensing

            # Only if we could dispense from a source well, we want to append plates and target well names to...
            # dispense dict value lists:
            dispense_dict["Destination Plate Name"].append(target_plate)
            dispense_dict["Source Plate Name"].append(source_plate)
            dispense_dict["Destination Well"].append(well)

            if not volume_col_name:
                dispense_dict["Transfer Volume"].append(transfer_volume_nl)
            else:
                dispense_dict["Transfer Volume"].append(tmp_target_df[tmp_target_df[target_df_well_col_name] == well][volume_col_name].values[0])

    # We do not have to check for NaNs in rows because all value lists in the dispense dict need to be of equal length
    # otherwise pd.DataFrame.from_dict will return a ValueError:
    dispense_df = pd.DataFrame.from_dict(dispense_dict)
    dispense_df.to_csv(csv_write_path)

    return dispense_df


def write_idot_dispense_csv(
    target_plates_list: list,
    source_df_plates_col_name: str,
    target_df_plates_col_name: str,
    source_df: pd.DataFrame,
    target_df: pd.DataFrame,
    source_match_col_name: str,
    target_match_col_name: str,
    source_df_well_col_name: str,
    target_df_well_col_name: str,
    csv_write_path: Union[str, os.PathLike[str]],
    target_source_subset_iter_list: list = None,
    target_source_subset_iter_col_name: str = "",
    volume_col_name: str = "",
    transfer_volume_nl: float = 0.00,
    max_source_plate_volume_nl: int = 68_000,
):
    """
    Write an I.DOT dispense instruction .csv sheet based on source and target plate DataFrames.

    This function generates a dispense instruction table for I.DOT liquid handlers. It matches
    samples between source and target plates, ensures source well volume constraints, handles
    multiple source subsets, and writes the resulting instructions in the required I.DOT CSV format.

    Parameters
    ----------
    target_plates_list : list
        List of target plate names that should be dispensed into.
    source_df_plates_col_name : str
        Column name in the source DataFrame containing source plate names.
    target_df_plates_col_name : str
        Column name in the target DataFrame containing target plate names.
    source_df : pd.DataFrame
        Source DataFrame storing sample well contents of source plates.
    target_df : pd.DataFrame
        Target DataFrame specifying samples to be dispensed into each target well.
    source_match_col_name : str
        Column in source DataFrame containing sample names for matching wells between source and target.
    target_match_col_name : str
        Column in target DataFrame containing sample names for matching wells between source and target.
    source_df_well_col_name : str
        Column in source DataFrame containing well identifiers.
    target_df_well_col_name : str
        Column in target DataFrame containing well identifiers.
    csv_write_path : Union[str, os.PathLike[str]]
        Path to write the I.DOT dispense CSV instruction file.
    target_source_subset_iter_list : list, optional
        Required when multiple different samples from different source subsets are dispensed into the same target wells.
        Indicates the order in which source subsets are used. Example: ['Fw', 'Rev'].
    target_source_subset_iter_col_name : str, optional
        Column in the target DataFrame specifying the source subset for each sample. Names must match those in
        `target_source_subset_iter_list`.
    volume_col_name : str, optional
        Column in the target DataFrame specifying variable transfer volumes for each well.
    transfer_volume_nl : float, optional
        Constant transfer volume (nL) for all wells if `volume_col_name` is not provided.
    max_source_plate_volume_nl : int, default 68000
        Maximum allowed volume in a source well (nL).

    Returns
    -------
    pd.DataFrame
        DataFrame containing the I.DOT dispense instructions that were written to the CSV file.

    Writes
    ------
    I.DOT dispense instruction CSV file in the format required for I.DOT liquid handlers.

    Notes
    -----
    - Unique samples are not allowed to be stored in multiple source plates. For example, the same
      primer cannot be stored in multiple wells across different source plates.
    - Ensures correct dispensing order by sorting wells and handling multiple source subsets.
    - Checks that no source well exceeds `max_source_plate_volume_nl` and switches to other wells if needed.
    - Only supports one source plate in the final CSV.
    - Adds template headers for I.DOT format and adjusts additional required fields.
    - Raises exceptions if sample matching fails, if source wells are insufficient, or if final CSV indices are inconsistent.
    """

    # Needed because otherwise the dfs in the notebook also get re-ordered below
    source_df = source_df.copy()
    target_df = target_df.copy()

    letters = "ABCDEFGHIJKLMNOP"
    col_numbers = np.arange(1, 25)

    reorder_target_by_well_list = []
    for letter in letters:
        for col_number in col_numbers:
            reorder_target_by_well_list.append(letter + str(col_number))

    sorter_index = dict(zip(reorder_target_by_well_list, range(len(reorder_target_by_well_list))))
    target_df["well_rank"] = target_df[target_df_well_col_name].map(sorter_index)

    if not target_source_subset_iter_col_name:
        target_df.sort_values(["well_rank"], inplace=True)
    else:
        target_df.sort_values([target_source_subset_iter_col_name, "well_rank"], inplace=True)

    # Initialize the columns of dispense .csv in a dispense_dict:
    dispense_dict = {
        "Source Plate Name": [],  # source_ID
        "Source Well": [],  # source_well_ID
        "Target Plate Name": [],  # target_ID
        "Target Well": [],  # target_well_ID
        "Transfer Volume": [],  # Vol
        "Liquid Name": [],
    }  # :TRAV_TRBV

    # Initialize the iterator lists for generating the dispensing instruction sheet:
    target_df_subsets_iterator_list = []
    target_plate_names_zip_list = []

    if not target_source_subset_iter_list:
        for target_plate in target_plates_list:
            target_df_subsets_iterator_list.append(target_df[target_df[target_df_plates_col_name] == target_plate].copy())
        target_plate_names_zip_list = target_plates_list

    else:
        for target_plate in target_plates_list:
            for source_subset in target_source_subset_iter_list:
                # If we do not pipette from the current source subset at all into the current target plate,
                # continue to next source subset:
                if target_df[(target_df[target_df_plates_col_name] == target_plate) & (target_df[target_source_subset_iter_col_name] == source_subset)].shape[0] == 0:
                    continue

                target_df_subsets_iterator_list.append(
                    target_df[(target_df[target_df_plates_col_name] == target_plate) & (target_df[target_source_subset_iter_col_name] == source_subset)].copy()
                )
                target_plate_names_zip_list.append(target_plate)

    if len(target_df_subsets_iterator_list) != len(target_plate_names_zip_list):
        raise Exception("Lengths of target_df_subset_iterator_list and target_plate_names_zip_list are out of sync!")

    source_well_counter = collections.defaultdict(int)
    for tmp_target_df, target_plate in zip(target_df_subsets_iterator_list, target_plate_names_zip_list):
        for well in tmp_target_df[target_df_well_col_name]:
            # sample_matcher stores the name of the Sample that needs to be dispensed into the current target well,
            # and this Sample name is used to check whether there is a matching source well that contains this Sample:
            sample_matcher = tmp_target_df[tmp_target_df[target_df_well_col_name] == well][target_match_col_name].values[0]

            source_plate_list = source_df[source_df[source_match_col_name] == sample_matcher][source_df_plates_col_name].unique().tolist()

            if len(source_plate_list) > 1:
                raise Exception("Current Sample", sample_matcher, "is stored in source wells divided over multiple source plates! This is not allowed!")
            source_plate = source_plate_list[0]

            matching_source_well_list = source_df[source_df[source_match_col_name] == sample_matcher][source_df_well_col_name].to_list()

            if not matching_source_well_list:
                raise Exception("No matching source wells found for this Sample! Are you sure Sample names have the same format in source and target df?")

            for idx, source_well in enumerate(matching_source_well_list):
                if source_well_counter[source_well] * transfer_volume_nl >= max_source_plate_volume_nl:
                    # If there are no source wells left for a specific Sample, more source wells should be prepared:
                    if len(matching_source_well_list[idx:]) == 1:
                        raise Exception(
                            "Not enough enough source wells for transferring required amount of nL!",
                            "Target plate:",
                            target_plate,
                            "Target well:",
                            well,
                            "Sample:",
                            sample_matcher,
                            "Last emptied Source well:",
                            source_well,
                            "Last emptied Source well count:",
                            source_well_counter[source_well],
                            "Volume last emptied Source well:",
                            source_well_counter[source_well] * transfer_volume_nl,
                        )
                    # If there are still other source wells left that also contain the Sample,
                    # dispense from next available source well:
                    continue
                else:
                    source_well_counter[source_well] += 1
                    dispense_dict["Source Well"].append(source_well)
                    dispense_dict["Liquid Name"].append(sample_matcher)
                    break  # exit the for loop if a well was used for dispensing

            # Only if we could dispense from a source well, we want to append plates and target well names to...
            # dispense dict value lists:
            dispense_dict["Target Plate Name"].append(target_plate)
            dispense_dict["Source Plate Name"].append(source_plate)
            dispense_dict["Target Well"].append(well)

            if not volume_col_name:
                dispense_dict["Transfer Volume"].append(transfer_volume_nl)
            else:
                dispense_dict["Transfer Volume"].append(tmp_target_df[tmp_target_df[target_df_well_col_name] == well][volume_col_name].values[0])

    # We do not have to check for NaNs in rows because all value lists in the dispense dict need to be of equal length
    # otherwise pd.DataFrame.from_dict will return a ValueError:
    dispense_df = pd.DataFrame.from_dict(dispense_dict)

    # Everything above the following lines is the same as in the write_echo_dispense_csv function
    # The following lines are only needed to convert to the I.DOT .csv file format
    dispense_df["Transfer Volume"] = "{:.2E}".format(transfer_volume_nl / 1000)
    if len(dispense_df["Source Plate Name"].unique()) > 1:
        raise Exception("There is more than one source plate in dispense_df, while we currently support only one source plate!")

    header_df = pd.read_csv(
        os.path.join(
            tcr_toolbox_data_path,
            "tcr_toolbox_datasets",
            "tcr_assembly",
            "v_gene_assignment_idot_source_plate_layout",
            "idot_header_template.csv",
        )
    )
    final_dispense_df = pd.DataFrame()
    tcr_idx = 0
    for target_plate in dispense_df["Target Plate Name"].unique():
        tmp_plate_dispense_df = dispense_df.loc[dispense_df["Target Plate Name"] == target_plate, ["Source Well", "Target Well", "Transfer Volume", "Liquid Name"]].copy()
        tmp_plate_dispense_df.rename({"Transfer Volume": "Volume [uL]"}, axis=1, inplace=True)
        tmp_plate_dispense_df.reset_index(drop=True, inplace=True)
        tmp_plate_dispense_df["Additional Volume Per Source Well"] = "{:.2E}".format(0)
        tmp_plate_dispense_df = pd.concat([header_df, tmp_plate_dispense_df])
        tmp_plate_dispense_df.reset_index(drop=True, inplace=True)
        tmp_plate_dispense_df.iloc[1, 1] = dispense_df["Source Plate Name"].unique()[0]
        tmp_plate_dispense_df.iloc[1, 5] = "plate_" + str(target_plate)

        if final_dispense_df.empty:
            tmp_plate_dispense_df.drop(3, axis=0, inplace=True)
            tmp_plate_dispense_df.reset_index(drop=True, inplace=True)
            final_dispense_df = tmp_plate_dispense_df.copy()
            tcr_idx = tmp_plate_dispense_df.shape[0]
        else:
            tmp_plate_dispense_df.drop(3, axis=0, inplace=True)
            tmp_plate_dispense_df.reset_index(drop=True, inplace=True)
            tmp_plate_dispense_df.index = [idx + tcr_idx for idx in tmp_plate_dispense_df.index]
            final_dispense_df = pd.concat([final_dispense_df, tmp_plate_dispense_df]).copy()
            tcr_idx += tmp_plate_dispense_df.shape[0]

    if final_dispense_df.index.duplicated().any():
        raise Exception("There are duplicate indices in final_dispense_df index! This means tcr_idx counting is broken!")

    if not all(final_dispense_df.index == np.arange(0, final_dispense_df.shape[0])):
        raise Exception("final_dispense_df index is not in order! This means tcr_idx counting is broken!")

    # if final_dispense_df.shape[0] != target_df.shape[0] + (target_df['plate'].unique().shape[0] * 3):
    #    raise Exception('final_dispense_df has incorrect number of rows!')

    final_dispense_df.to_csv(csv_write_path, header=False, index=False)

    return final_dispense_df


def remove_leading_zeros_well_coordinates(well_column: Union[list, pd.Series]):
    """
    Rack scanners sometimes have hardcoded that the well annotation is A01, A02, A03, ..., A12, B01, B02, B03.
    Some machines cannot read wells like this and require A1, A2, A3, ..., A12, B1, B2, B3. etc.
    This function solves that issue.

    Parameters
    ----------
    well_column : Union[list, pd.Series]
        well_column_adapted that contains the well column that needs to be cleaned

    Returns
    -------
    well_column_adapted: list
        A list that contains the well column that has been cleaned
    """

    char_list, remainder = zip(*[(well[0], int(well[1:])) for well in well_column])
    remainder = [str(int(rem)) for rem in remainder]
    if "0" in remainder:
        raise Exception("A remainder is 0, this is not allowed! as it suggests a well is named A00, A000, etc. Which is ambiguous!")

    well_column_adapted = [char + rem for char, rem in zip(char_list, remainder)]

    return well_column_adapted


def add_tcr_name_to_tcr_df(tcr_df: pd.DataFrame, grouping_col: str):
    tcr_df["tcr_name"] = tcr_df.groupby(grouping_col).cumcount() + 1  # Starts numbering from 1
    tcr_df["tcr_name"] = tcr_df[grouping_col] + "_" + tcr_df["tcr_name"].astype(str)

    if tcr_df["tcr_name"].duplicated().any():
        raise Exception("There are duplicate TCR names in the tcr_df!")

    return tcr_df


def stream_stderr(pipe: io.IOBase, logfile: str | os.PathLike, is_bytes: bool = True):
    if is_bytes:
        for line in iter(pipe.readline, b""):
            logfile.write(line.decode())
            logfile.flush()
    else:
        for line in iter(pipe.readline, ""):
            logfile.write(line)
            logfile.flush()

    pipe.close()
