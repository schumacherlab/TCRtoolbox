import pandas as pd


def calculate_clonotype_frequencies(
    tcr_meta_df: pd.DataFrame,
    donor_col: str,
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
    """Count cell barcodes per clonotype and compute clonotype frequencies within each donor.

    For each of `clonotype_nt`, `clonotype_aa`, `unique_clonotype_nt`, and `unique_clonotype_aa`,
    counts the number of unique cell barcodes sharing that clonotype string within each `donor_col`
    group, and divides by the total number of unique barcodes in that donor to get the frequencies
    of these clonotype strings.

    Parameters
    ----------
    tcr_meta_df : DataFrame
        A dataframe that contains cell barcodes (rows) annotated with clonotype strings, CD8_and_CD4 calls,
        and a donor group.
    donor_col : str
        Column identifying the donor/patient (or e.g. donor + coreceptor) to compute frequencies
        within.
    clonotype_nt_col, clonotype_aa_col : str, default="clonotype_nt", "clonotype_aa"
        Columns holding the (possibly dual-chain) clonotype strings.
    unique_clonotype_nt_col, unique_clonotype_aa_col : str, default="unique_clonotype_nt", "unique_clonotype_aa"
        Columns holding the single_pair clonotype strings.
    clonotype_count_nt_col, clonotype_count_aa_col, clonotype_freq_nt_col, clonotype_freq_aa_col : str
        Output column names for the `clonotype_nt` and `clonotype_aa` counts and frequencies.
    unique_clonotype_count_nt_col, unique_clonotype_count_aa_col, unique_clonotype_freq_nt_col, unique_clonotype_freq_aa_col : str
        Output column names for the `unique_clonotype_nt`and `unique_clonotype_aa` counts and
        frequencies.

    Returns
    -------
    DataFrame
        `tcr_meta_df` with the 8 count/frequency columns above added, sorted by
        `unique_clonotype_freq_aa_col` descending.
    """
    tcr_meta_df = tcr_meta_df.copy()
    tcr_meta_df[donor_col] = tcr_meta_df[donor_col].astype(str)

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
    clonotype_freq_col: str = "clonotype_frequency_aa",
    clonotype_col: str = "unique_clonotype_aa",
    clonotype_beta_col: str = "unique_clonotype_beta_aa",
    clonotype_shared_with_col: str = "shared_with",
):
    """Keep only the first, most frequent occurrence of each duplicate clonotype per donor.

    Within each `donor_col` group, keeps only the first, most frequent occurrence of each
    duplicate `clonotype_col` clonotype string, dropping the rest (duplicate `clonotype_col`
    strings can have different `clonotype_freq_col` values if these originate from different
    dual chain clonotypes). If `also_drop_clonotype_beta_aa_duplicates` is set, a second such
    pass then keeps only the first, most frequent occurrence of each duplicate `clonotype_beta_col`
    beta chain within each donor. Finally, adds `clonotype_shared_with_col`, listing for
    clonotype strings found in more than one donor the sorted, "-" joined set of all donors
    carrying that clonotype (including the row's own donor), or `None` if the clonotype is only
    found in one donor. This column is useful to identify public clonotypes in your cohort.

    Parameters
    ----------
    tcr_meta_df : DataFrame
        A dataframe that contains cell barcodes (rows) annotated with clonotype strings, CD8_and_CD4 calls,
        and a donor group that are further annotated with clonotype_count and frequency columns by
        `calculate_clonotype_frequencies`.
    also_drop_clonotype_beta_aa_duplicates : bool, default=False
        If True, also keep only the first, most frequent occurrence of each duplicate
        `clonotype_beta_col` beta chain within each donor, after the primary `clonotype_col` step.
    donor_col : str, default="patient"
        Column identifying the donor/patient group within which duplicates are resolved and
        cross-donor sharing is detected; cast to str.
    clonotype_freq_col : str, default="clonotype_frequency_aa"
        Column used to rank duplicates; the highest-value row is kept.
    clonotype_col : str, default="unique_clonotype_aa"
        Column whose duplicate occurrences are resolved (one row kept per distinct value, per donor).
    clonotype_beta_col : str, default="unique_clonotype_beta_aa"
        Column used for the optional second, beta-chain-only pass.
    clonotype_shared_with_col : str, default="shared_with"
        Output column name for the cross-donor sharing list.

    Returns
    -------
    DataFrame
        `tcr_meta_df` with duplicate clonotypes dropped and `clonotype_shared_with_col` added.
    """
    tcr_meta_df = tcr_meta_df.copy()
    tcr_meta_df[donor_col] = tcr_meta_df[donor_col].astype(str)
    # Drop duplicates within each donor group
    print("Number of TCRs before dropping duplicate clonotypes:", tcr_meta_df.shape[0])
    tcr_meta_df = tcr_meta_df.sort_values(clonotype_freq_col, ascending=False).drop_duplicates(subset=[donor_col, clonotype_col], keep="first").reset_index(drop=True)

    if also_drop_clonotype_beta_aa_duplicates:
        tcr_meta_df = tcr_meta_df.sort_values(clonotype_freq_col, ascending=False).drop_duplicates(subset=[donor_col, clonotype_beta_col], keep="first").reset_index(drop=True)

    print("Number of TCRs after dropping duplicate clonotypes:", tcr_meta_df.shape[0])
    # Identify shared clonotypes across donors
    clonotype_patient_map = tcr_meta_df.groupby(clonotype_col)[donor_col].apply(lambda donors: "-".join(sorted(donors.unique())) if len(donors) > 1 else None)

    tcr_meta_df[clonotype_shared_with_col] = tcr_meta_df[clonotype_col].map(clonotype_patient_map)

    tcr_meta_df = tcr_meta_df.sort_values([donor_col, clonotype_freq_col], ascending=[True, False]).reset_index(drop=True)

    return tcr_meta_df
