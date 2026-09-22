import os
import pandas as pd
import re
import pysam
from collections import defaultdict


def remove_dna_barcode_from_string(string: str, barcode_length: int) -> str:
    string = re.sub(rf"[ACGT]{{{barcode_length}}}", "", string)  # Remove ACGT barcode
    string = re.sub(r"_\d+_\d+$", "", string)  # Remove trailing "_num_num"
    string = string.replace("-", "_")  # Replace hyphens with underscores so that only TCR and Epi are separated by "-"
    string = re.sub(r"^_+", "", string)  # Remove leading underscores
    string = re.sub(r"_+$", "", string)  # Remove trailing underscores
    string = re.sub(r"(?<=[ABC]\d{4})_\d+$", "", string)  # Remove technical duplicate suffix behind HLA code (e.g., A0201) if it exists
    string = string.replace("__", "_")  # Collapse 2x underscores to 1x underscore
    return string


def collapse_duplicate_tcr_names_reference_file(reference_file: str | os.PathLike) -> dict:
    # Read all TCR names from FASTA headers
    tcr_name_list = []
    with pysam.FastxFile(reference_file, "r") as fasta_in:
        for entry in fasta_in:
            tcr_name_list.append(entry.name)

    # Regex to extract prefix, well, mid number, and patient
    pattern = re.compile(r"^tcr_(\d+_\d+)_([A-Z]\d+)_([0-9]+)_(.+)$")

    # Group entries by (prefix, mid_num, patient)
    grouped = defaultdict(list)

    def well_to_coord(well: str):
        row = ord(well[0].upper()) - 64  # 1-based row index
        col = int(well[1:])  # 1-based column
        return (row, col)

    # Assign entries to groups
    for tcr in tcr_name_list:
        match = pattern.match(tcr)
        if not match:
            continue
        prefix, well, mid_num, patient = match.groups()
        coord = well_to_coord(well)
        grouped[(prefix, mid_num, patient)].append((tcr, coord))

    # Build dictionary: key = entry, value = canonical lowest well
    canonical_dict = {}
    for key, group in grouped.items():
        group_sorted = sorted(group, key=lambda x: (x[1][0], x[1][1]))  # plate order
        canonical_tcr = group_sorted[0][0]
        for tcr_entry, _ in group:
            canonical_dict[tcr_entry] = canonical_tcr

    return canonical_dict


def overlap_between_ref_and_count_names(
    count_file: str | os.PathLike, count_file_ref_name_col: str, reference_file: str | os.PathLike, library_ref_id_list: list = None, print_diff_names: bool = True
):
    """
    Print # of overlapping reference names between reference .fa and counts csv/tsv file.

    # diff, # union, and, if print_diff_names = True, diff reference names are also printed. This function should be used for
    checking whether the generated library cell line covers all unique transcripts of the library.

    Parameters
    ----------
    count_file : str | os.PathLike
        String path to count file.
    count_file_ref_name_col : str
        Column that contains reference names in count file.
    reference_file : str | os.PathLike
        String path to reference .fa file.
    library_ref_id_list: list
        Only include reference names that match a reference ID.
    print_diff_names : bool, default = True
        If True, print non-overlapping difference between count file and reference file reference names.

    Examples
    --------
    >>> counts_dir = project_dir / "counts"
    >>> for count_file in counts_dir.glob(*.csv'):
    >>>     print(count_file)
    >>>     overlap_between_ref_and_count_names(count_file = count_file,
    ...                                         count_file_ref_name_col = 'reference_name',
    ...                                         reference_file = '/path/to/project/references/beta.fa',
    ...                                         print_diff_names = True
    ...                                         )
    >>>     print('\n')
    /path/to/project/data/run1/TCR/counts/sample1_S13_R1_001_counts.csv
    # intersecting: 92
    # diff: 1
    # union: 93
    Diff: {'102_198'}

    /path/to/project/data/run1/TCR/counts/sample2_S14_R1_001_counts.csv
    # intersecting: 92
    # diff: 1
    # union: 93
    Diff: {'102_198'}

    >>> overlap_between_ref_and_count_names(count_file = counts_dir / 'sample3_S8_R1_001_counts.tsv',
    ...                                     count_file_ref_name_col = 'gene',
    ...                                     reference_file = '/path/to/project/references/beta_plate.fa',
    ...                                     print_diff_names = True
    ...                                     )
    # intersecting: 291
    # diff: 7
    # union: 298
    Diff: {'tcr_DMF5', 'epi_GLC', 'tcr_102_198', 'tcr_DMF4', 'tcr_33_57', 'epi_p20_stuffer', 'tcr_1G4'}
    """
    with pysam.FastxFile(reference_file, "r") as fasta_in:
        if not library_ref_id_list:
            ref_names_set = {entry.name for entry in fasta_in}
        else:
            ref_names_set = set()
            ref_ids = set(library_ref_id_list)
            for entry in fasta_in:
                if any(ref_id in entry.name for ref_id in ref_ids):
                    ref_names_set.add(entry.name)

    # Only python engine allows regex separators:
    count_df = pd.read_csv(count_file, sep=r"\t|,", engine="python", usecols=[count_file_ref_name_col])
    if len(count_df[count_file_ref_name_col].unique()) != count_df.shape[0]:
        raise Exception("There are duplicate reference name in the count file!")
    count_names_set = set(count_df[count_file_ref_name_col])

    intersecting = ref_names_set.intersection(count_names_set)
    print("# intersecting:", len(intersecting))
    print("# diff:", len(ref_names_set.difference(count_names_set)))
    print("# union:", len(ref_names_set.union(count_names_set)))

    sensitivity = len(intersecting) / len(ref_names_set) * 100
    print(f"Sensitivity: {sensitivity:.2f}% ({len(intersecting)}/{len(ref_names_set)} reference names found in counts)")

    if print_diff_names:
        print("refs not in count names:", ref_names_set.difference(count_names_set))
        print("count names not in ref subset:", count_names_set.difference(ref_names_set))


def read_gdna_counts_csv(
    count_file: str | os.PathLike,
    count_file_ref_name_col: str,
    collapse_epitope_barcode: bool = False,
    collapse_custom_ref_names_dict: dict = None,
    collapse_tcr_technical_duplicates: bool = False,
    reference_file=None,
    epitope_barcode_length: int = 18,
) -> pd.DataFrame:
    """
    Read gDNA read count .csv file that was generated by :func:`tcr_toolbox.sequencing_analysis.count_aligned_gdna`

    Parameters
    ----------
    count_file : str | os.PathLike
        String path to count file.
    count_file_ref_name_col : str
        Column that contains reference names in count file.
    collapse_epitope_barcode : bool, default = False
        If True, remove DNA barcode from epitope reference name string. Used to collapse DNA barcodes counts
        to their corresponding epitope reference name.
    epitope_barcode_length : int, required if collapse_epitope_barcode = True
        Epitope barcode length. Used to remove the DNA barcode from epitope reference name to collapse epitopes per barcode.

    Returns
    -------
    count_df
        DataFrame with read counts per reference name.
    """
    # Only python engine allows regex separators:
    count_df = pd.read_csv(count_file, sep=r"\t|,", engine="python")
    count_df["reference_name_working"] = count_df[count_file_ref_name_col]

    if collapse_epitope_barcode:
        count_df["reference_name_working"] = count_df["reference_name_working"].apply(lambda x: remove_dna_barcode_from_string(x, epitope_barcode_length))

    if collapse_tcr_technical_duplicates:
        duplicate_tcr_collapse_dict = collapse_duplicate_tcr_names_reference_file(reference_file=reference_file)

        if duplicate_tcr_collapse_dict:
            first_key = next(iter(duplicate_tcr_collapse_dict))
            if first_key.startswith("tcr_"):
                duplicate_tcr_collapse_dict = {k.removeprefix("tcr_"): v.removeprefix("tcr_") for k, v in duplicate_tcr_collapse_dict.items()}

            count_df["reference_name_working"] = count_df["reference_name_working"].map(duplicate_tcr_collapse_dict)

            if count_df["reference_name_working"].isna().any():
                missing = count_df.loc[count_df["reference_name_working"].isna(), count_file_ref_name_col].unique()
                raise Exception(f"Missing TCR names in collapse map: {missing}")

    if collapse_custom_ref_names_dict is not None:
        count_df["reference_name_working"] = count_df["reference_name_working"].replace(collapse_custom_ref_names_dict)

    count_df = count_df.groupby("reference_name_working", dropna=False).sum()
    count_df.index.name = "reference_name"
    count_df.drop(columns=[count_file_ref_name_col], inplace=True, errors="ignore")

    return count_df
