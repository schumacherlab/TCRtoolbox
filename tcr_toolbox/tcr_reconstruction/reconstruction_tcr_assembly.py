import pandas
import numpy as np
from tcr_toolbox.tcr_reconstruction.constants import pydna_tcr_golden_gate_constant_aa_seq_dict
from tcr_toolbox.tcr_reconstruction.reconstruction_utils import reconstruct_vdj, reconstruct_full_tcr
import os
from dotenv import load_dotenv

load_dotenv()
tcr_toolbox_data_path = os.getenv("tcr_toolbox_data_path")


def reconstruct_tcrs_assembly(
    include_leader: bool = True,
    include_constant: bool = True,
    mouse_or_human: str = "mouse",
    constant_beta: str = pydna_tcr_golden_gate_constant_aa_seq_dict["muTRBC_aa"],
    constant_alpha: str = pydna_tcr_golden_gate_constant_aa_seq_dict["muTRAC_aa"],
    exclude_c_fw: bool = False,
    plates_df_dict: dict = None,  # was dataset
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
    translation_dict_aa = tcr_toolbox_data_path + "/tcr_toolbox_datasets/tcr_reconstruction/VDJ_gene_sequences/after_benchmark/functional_with_L-PART1+V-EXON_after_benchmark/20230803_vdj_translation_dict_aa.json"

    translation_dict_nt = tcr_toolbox_data_path + "/tcr_toolbox_datasets/tcr_reconstruction/VDJ_gene_sequences/after_benchmark/functional_with_L-PART1+V-EXON_after_benchmark/20230803_vdj_translation_dict_nt.json"

    for plate in plates_df_dict.keys():
        tmp_plate_df = plates_df_dict[plate].copy()

        # Fix to be able to use Bjørn's function:
        tmp_plate_df["TRAV"] = tmp_plate_df["TRAV_IMGT_allele_collapsed"].str.split("*").str[0].str.replace("_", "/")
        tmp_plate_df["TRBV"] = tmp_plate_df["TRBV_IMGT_allele_collapsed"].str.split("*").str[0].str.replace("_", "/")
        tmp_plate_df["TRAJ"] = tmp_plate_df["TRAJ_IMGT"].str.split("*").str[0].str.replace("_", "/")
        tmp_plate_df["TRBJ"] = tmp_plate_df["TRBJ_IMGT"].str.split("*").str[0].str.replace("_", "/")

        if not all([column in tmp_plate_df.columns for column in ["TRAV", "TRAJ", "TRBV", "TRBJ"]]):
            raise ValueError("Columns TRAV, TRAJ, TRBV, TRBJ are not present in the dataframe.")

        tmp_plate_df = tmp_plate_df.reset_index(drop=True)
        # if not TRBD present, add it
        if "TRBD" not in tmp_plate_df.columns:
            tmp_plate_df["TRBD"] = None

        for vdj in ["TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ"]:
            (tmp_plate_df[vdj + "_imgt_aa"], tmp_plate_df[vdj + "_seq_aa"], tmp_plate_df[vdj + "_imgt_nt"], tmp_plate_df[vdj + "_seq_nt"]) = reconstruct_vdj(
                tmp_plate_df, vdj, translation_dict_nt, translation_dict_aa, verbose=verbose
            )
        total_len = len(tmp_plate_df)
        for nt_or_aa in ["aa", "nt"]:
            for vdj in ["TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ"]:
                count = sum(tmp_plate_df[vdj + "_seq_" + nt_or_aa].notna())
                count_original = sum(tmp_plate_df[vdj].notna())
                if verbose:
                    print("{0} imputed: {1} / {4} total annotations ({2}) ({3})".format(vdj, count, total_len, nt_or_aa.upper(), count_original))

        # beta
        tmp_plate_df["full_seq_reconstruct_beta_aa"] = reconstruct_full_tcr(
            df=tmp_plate_df,
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
        tmp_plate_df["full_seq_reconstruct_alpha_aa"] = reconstruct_full_tcr(
            df=tmp_plate_df,
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
                    sum(tmp_plate_df["full_seq_reconstruct_beta_aa"].notna()), sum(tmp_plate_df["cdr3_beta_aa"].notna())
                )
            )

            print(
                "Could  reconstruct full ALPHA TCR for {0} entries of total {1} CDR3a entries".format(
                    sum(tmp_plate_df["full_seq_reconstruct_alpha_aa"].notna()), sum(tmp_plate_df["cdr3_alpha_aa"].notna())
                )
            )

        if tmp_plate_df.index.duplicated().any():
            raise Exception("plate df index contains duplicates!")

        reconstructed_tcrab_full_list = []
        for idx in tmp_plate_df.index:
            if pandas.isna(tmp_plate_df.loc[idx, "full_seq_reconstruct_alpha_aa"]) or pandas.isna(tmp_plate_df.loc[idx, "full_seq_reconstruct_beta_aa"]):
                reconstructed_tcrab_full_list.append(np.nan)
            else:
                reconstructed_tcrab_full_list.append(
                    tmp_plate_df.loc[idx, "full_seq_reconstruct_beta_aa"]
                    + pydna_tcr_golden_gate_constant_aa_seq_dict["P2A_aa"]
                    + tmp_plate_df.loc[idx, "full_seq_reconstruct_alpha_aa"]
                    + pydna_tcr_golden_gate_constant_aa_seq_dict["T2A_aa"]
                    + pydna_tcr_golden_gate_constant_aa_seq_dict["puroR_aa"]
                )

        plates_df_dict[plate]["reconstructed_tcrab_full_aa"] = reconstructed_tcrab_full_list

    return plates_df_dict
