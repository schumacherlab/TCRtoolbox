import pandas
from tcr_toolbox.tcr_reconstruction.reconstruction_utils import reconstruct_vdj, reconstruct_full_tcr
from typing import Union
from tcr_toolbox.utils.logger import init_logger
from dotenv import load_dotenv
import os

def reconstruct_tcrs_simple(
    include_leader: bool,
    include_constant: bool,
    mouse_or_human: Union[str, None],
    constant_beta: Union[str, None],
    constant_alpha: Union[str, None],
    exclude_c_fw: bool,
    dataframe: pandas.DataFrame = None,
    dataset_file_path: str = None,
    output_file_path: str = None,
    verbose=False,
):
    """Reconstruct TCRs from a dataset file.

    Parameters
    ----------
    include_leader : bool
        Whether to include the leader sequence in the reconstruction.
    include_constant : bool
        Whether to include the constant sequence in the reconstruction.
    mouse_or_human : str
        Whether to use murine or human constant sequences.
    constant_beta : str
        Constant sequence for beta chain.
    constant_alpha : str
        Constant sequence for alpha chain.
    exclude_c_fw : bool
        Whether to exclude the C region framework from the reconstruction.
    dataframe : pandas.DataFrame
        pandas dataframe containing the dataset. Either provide this or dataset_file_path.
    dataset_file_path : str
        Path to the dataset file. Either provide this or a dataframe.
    output_file_path : str
        Path to the output file.

    Returns
    -------
    list
        List of reconstructed TCRs.
    """
    load_dotenv()
    tcr_toolbox_data_path = os.getenv("tcr_toolbox_data_path")

    # logger = init_logger("TCR_reconstruction.log", level_msg="INFO")
    tcr_toolbox_dir = tcr_toolbox_data_path + "/tcr_toolbox_datasets"
    translation_dict_aa = (
        tcr_toolbox_dir + "/tcr_reconstruction/VDJ_gene_sequences/after_benchmark/functional_with_L-PART1+V-EXON_after_benchmark/20230803_vdj_translation_dict_aa.json"
    )
    translation_dict_nt = (
        tcr_toolbox_dir + "/tcr_reconstruction/VDJ_gene_sequences/after_benchmark/functional_with_L-PART1+V-EXON_after_benchmark/20230803_vdj_translation_dict_nt.json"
    )

    # either df or dataset_file_path should be provided
    if dataframe is not None:
        dataset = dataframe
    elif dataset_file_path:
        # try different delimiters untill the data loads properly
        dataset = pandas.read_csv(dataset_file_path, delimiter=";", low_memory=False)
        # remove the index, as it might cause errors when looping over the dataset
        dataset = dataset.reset_index(drop=True)

        # if only one column, try other delimiters
        if dataset.shape[1] == 1:
            dataset = pandas.read_csv(dataset_file_path, delimiter="\t", low_memory=False)
        if dataset.shape[1] == 1:
            dataset = pandas.read_csv(dataset_file_path, low_memory=False)
        if dataset.shape[1] == 1:
            raise ValueError("Could not read the dataset file properly. Please check the file format.")
    else:
        raise ValueError("Either dataframe or dataset_file_path should be provided.")

    if not all([x in dataset.columns for x in ["TRAV", "TRAJ", "TRBV", "TRBJ"]]):
        raise ValueError("Columns TRAV, TRAJ, TRBV, TRBJ are not present in the dataframe.")

    dataset = dataset.reset_index(drop=True)
    if "TRBD" not in dataset.columns:
        dataset["TRBD"] = None

    for vdj in ["TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ"]:
        (dataset[vdj + "_imgt_aa"], dataset[vdj + "_seq_aa"], dataset[vdj + "_imgt_nt"], dataset[vdj + "_seq_nt"]) = reconstruct_vdj(
            dataset, vdj, translation_dict_nt, translation_dict_aa, verbose=verbose
        )
    total_len = len(dataset)
    for nt_or_aa in ["aa", "nt"]:
        for vdj in ["TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ"]:
            count = sum(dataset[vdj + "_seq_" + nt_or_aa].notna())
            count_original = sum(dataset[vdj].notna())
            if verbose:
                print("{0} imputed: {1} / {4} total annotations ({2}) ({3})".format(vdj, count, total_len, nt_or_aa.upper(), count_original))

    # beta
    dataset["full_seq_reconstruct_beta_aa"] = reconstruct_full_tcr(
        df=dataset,
        v_column_nt="TRBV_seq_nt",
        v_column_aa="TRBV_seq_aa",
        j_column_nt="TRBJ_seq_nt",
        j_column_aa="TRBJ_seq_aa",
        cdr3_column_aa="cdr3_beta_aa",
        include_leader=include_leader,
        exclude_C_F_W=exclude_c_fw,
        include_constant=include_constant,
        constant_sequence=constant_beta,
        human_or_mouse_constant=mouse_or_human,
    )
    # alpha
    dataset["full_seq_reconstruct_alpha_aa"] = reconstruct_full_tcr(
        df=dataset,
        v_column_nt="TRAV_seq_nt",
        v_column_aa="TRAV_seq_aa",
        j_column_nt="TRAJ_seq_nt",
        j_column_aa="TRAJ_seq_aa",
        cdr3_column_aa="cdr3_alpha_aa",
        include_leader=include_leader,
        exclude_C_F_W=exclude_c_fw,
        include_constant=include_constant,
        constant_sequence=constant_alpha,
        human_or_mouse_constant=mouse_or_human,
    )
    if verbose:
        print(
            "Could reconstruct full BETA TCR for {0} entries of total {1} CDR3b entries".format(
                sum(dataset["full_seq_reconstruct_beta_aa"].notna()), sum(dataset["cdr3_beta_aa"].notna())
            )
        )

        print(
            "Could  reconstruct full ALPHA TCR for {0} entries of total {1} CDR3a entries".format(
                sum(dataset["full_seq_reconstruct_alpha_aa"].notna()), sum(dataset["cdr3_alpha_aa"].notna())
            )
        )

    if output_file_path:
        dataset.to_csv(output_file_path)
    else:
        return dataset
