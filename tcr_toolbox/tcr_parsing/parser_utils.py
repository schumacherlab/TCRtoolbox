import logging
import pandas as pd

def calculate_clonotype_frequencies(
    tcr_meta_df: pd.DataFrame,
    donor_col: str,
    detected_chains_col: str = "detected_chains",
    clonotype_nt_col: str = "clonotype_nt",
    clonotype_aa_col: str = "clonotype_aa",
    unique_clonotype_nt_col: str = "unique_clonotype_nt",
    unique_clonotype_aa_col: str = "unique_clonotype_aa",
    clonotype_count_nt_col: str = "clonotype_count_nt",
    clonotype_count_aa_col: str = "clonotype_count_aa",
    clonotype_freq_nt_col: str = "clonotype_frequency_nt",
    clonotype_freq_aa_col: str = "clonotype_frequency_aa",
    unique_clonotype_count_nt_col: str = "unique_clonotype_count_nt",
    unique_clonotype_count_aa_col: str = "unique_clonotype_count_aa",
    unique_clonotype_freq_nt_col: str = "unique_clonotype_frequency_nt",
    unique_clonotype_freq_aa_col: str = "unique_clonotype_frequency_aa",
):
    tcr_meta_df = tcr_meta_df.copy()
    tcr_meta_df[donor_col] = tcr_meta_df[donor_col].astype(str)
    # dual_chain_names_list = ["dual_beta", "dual_alpha", "dual_alpha_and_dual_beta"]

    # dual_mask = tcr_meta_df[detected_chains_col].isin(dual_chain_names_list)
    # non_dual_tcr_meta_df = tcr_meta_df.loc[~dual_mask].copy()

    # dual_tcr_meta_df = (tcr_meta_df.loc[dual_mask].drop_duplicates(
    #             subset=[donor_col, "barcode"], keep="first"
    #         ).copy())

    # freq_calculate_df = pd.concat([non_dual_tcr_meta_df, dual_tcr_meta_df], ignore_index=True)
    freq_calculate_df = tcr_meta_df.copy()
    total_cells_per_donor = freq_calculate_df.groupby(donor_col)["barcode"].nunique().rename("total_cells").reset_index()

    counts_nt_df = freq_calculate_df.groupby([donor_col, clonotype_nt_col])["barcode"].nunique().reset_index(name=clonotype_count_nt_col)
    counts_nt_df = counts_nt_df.merge(total_cells_per_donor, on=donor_col, how="left")
    counts_nt_df[clonotype_freq_nt_col] = counts_nt_df[clonotype_count_nt_col] / counts_nt_df["total_cells"]

    tcr_meta_df = tcr_meta_df.merge(counts_nt_df[[donor_col, clonotype_nt_col, clonotype_count_nt_col, clonotype_freq_nt_col]], on=[donor_col, clonotype_nt_col], how="left")

    counts_aa_df = freq_calculate_df.groupby([donor_col, clonotype_aa_col])["barcode"].nunique().reset_index(name=clonotype_count_aa_col)
    counts_aa_df = counts_aa_df.merge(total_cells_per_donor, on=donor_col, how="left")
    counts_aa_df[clonotype_freq_aa_col] = counts_aa_df[clonotype_count_aa_col] / counts_aa_df["total_cells"]

    tcr_meta_df = tcr_meta_df.merge(counts_aa_df[[donor_col, clonotype_aa_col, clonotype_count_aa_col, clonotype_freq_aa_col]], on=[donor_col, clonotype_aa_col], how="left")

    counts_nt_df = freq_calculate_df.groupby([donor_col, unique_clonotype_nt_col])["barcode"].nunique().reset_index(name=unique_clonotype_count_nt_col)
    counts_nt_df = counts_nt_df.merge(total_cells_per_donor, on=donor_col, how="left")
    counts_nt_df[unique_clonotype_freq_nt_col] = counts_nt_df[unique_clonotype_count_nt_col] / counts_nt_df["total_cells"]

    tcr_meta_df = tcr_meta_df.merge(
        counts_nt_df[[donor_col, unique_clonotype_nt_col, unique_clonotype_count_nt_col, unique_clonotype_freq_nt_col]], on=[donor_col, unique_clonotype_nt_col], how="left"
    )

    counts_aa_df = freq_calculate_df.groupby([donor_col, unique_clonotype_aa_col])["barcode"].nunique().reset_index(name=unique_clonotype_count_aa_col)
    counts_aa_df = counts_aa_df.merge(total_cells_per_donor, on=donor_col, how="left")
    counts_aa_df[unique_clonotype_freq_aa_col] = counts_aa_df[unique_clonotype_count_aa_col] / counts_aa_df["total_cells"]

    tcr_meta_df = tcr_meta_df.merge(
        counts_aa_df[[donor_col, unique_clonotype_aa_col, unique_clonotype_count_aa_col, unique_clonotype_freq_aa_col]], on=[donor_col, unique_clonotype_aa_col], how="left"
    )

    tcr_meta_df.sort_values(unique_clonotype_freq_aa_col, inplace=True, ascending=False)

    return tcr_meta_df


def drop_clonotype_duplicates(
    tcr_meta_df: pd.DataFrame,
    also_drop_clonotype_beta_aa_duplicates: bool = False,
    donor_col: str = "patient",
    clonotype_freq_col: str = "unique_clonotype_frequency_aa",
    clonotype_col: str = "unique_clonotype_aa",
    clonotype_beta_col: str = "unique_clonotype_beta_aa",
    clonotype_shared_with_col: str = "shared_with",
):
    tcr_meta_df = tcr_meta_df.copy()
    tcr_meta_df[donor_col] = tcr_meta_df[donor_col].astype(str)
    # Drop duplicates within each patient group
    print("Number of TCRs before dropping duplicate clonotypes:", tcr_meta_df.shape[0])
    tcr_meta_df = tcr_meta_df.sort_values(clonotype_freq_col, ascending=False).drop_duplicates(subset=[donor_col, clonotype_col], keep="first").reset_index(drop=True)

    if also_drop_clonotype_beta_aa_duplicates:
        tcr_meta_df = tcr_meta_df.sort_values(clonotype_freq_col, ascending=False).drop_duplicates(subset=[donor_col, clonotype_beta_col], keep="first").reset_index(drop=True)

    print("Number of TCRs after dropping duplicate clonotypes:", tcr_meta_df.shape[0])
    # Identify shared clonotypes across patients
    clonotype_patient_map = tcr_meta_df.groupby(clonotype_col)[donor_col].apply(lambda donors: "-".join(sorted(donors.unique())) if len(donors) > 1 else None)

    tcr_meta_df[clonotype_shared_with_col] = tcr_meta_df[clonotype_col].map(clonotype_patient_map)

    return tcr_meta_df
