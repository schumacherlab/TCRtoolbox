import collections
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tcr_toolbox.utils.plot_utils import startfig


def compute_umi_ratio(barcode_df: pd.DataFrame, umi_col: str):
    sorted_barcode_df = barcode_df.sort_values(umi_col, ascending=False)
    return sorted_barcode_df[umi_col].iloc[1] / sorted_barcode_df[umi_col].iloc[0]


def construct_clonotype_string(barcode_df: pd.DataFrame, chain_prefix: str, cdr3_nt_col: str, cdr3_aa_col: str):
    nt = [f"{row[chain_prefix + 'V']}_{row[chain_prefix + 'J']}_{row[cdr3_nt_col]}" for _, row in barcode_df.iterrows()]
    aa = [f"{row[chain_prefix + 'V']}_{row[chain_prefix + 'J']}_{row[cdr3_aa_col]}" for _, row in barcode_df.iterrows()]
    return "_".join(nt), "_".join(aa)


def process_dual_alpha_and_dual_beta_barcode(barcode: str, dual_alpha_and_dual_beta_contig_df: pd.DataFrame, umi_ratio_threshold: float = None):
    barcode_df = dual_alpha_and_dual_beta_contig_df.loc[dual_alpha_and_dual_beta_contig_df.loc[:, "barcode"] == barcode, :].copy()

    barcode_df.loc[:, "tmp_clonotype_alpha_nt"] = barcode_df.loc[:, "TRAV"] + "_" + barcode_df.loc[:, "TRAJ"] + "_" + barcode_df.loc[:, "cdr3_alpha_nt"]

    barcode_df.loc[:, "tmp_clonotype_beta_nt"] = barcode_df.loc[:, "TRBV"] + "_" + barcode_df.loc[:, "TRBJ"] + "_" + barcode_df.loc[:, "cdr3_beta_nt"]

    alpha_clonotype_aa_to_indices = barcode_df.groupby("tmp_clonotype_alpha_nt").groups
    alpha_barcode_df = barcode_df.drop_duplicates("tmp_clonotype_alpha_nt").copy()
    beta_clonotype_aa_to_indices = barcode_df.groupby("tmp_clonotype_beta_nt").groups
    beta_barcode_df = barcode_df.drop_duplicates("tmp_clonotype_beta_nt").copy()

    alpha_drop = None
    beta_drop = None

    if (alpha_barcode_df.shape[0] != 2) or (beta_barcode_df.shape[0] != 2):
        raise Exception("Both alpha_barcode_df and beta_barcode_df should have 2 rows in dual_alpha_and_dual_beta barcodes!", alpha_barcode_df.shape[0], beta_barcode_df.shape[0])

    alpha_umi_ratio = compute_umi_ratio(alpha_barcode_df, "umis_a")
    if umi_ratio_threshold and alpha_umi_ratio < umi_ratio_threshold:
        sorted_alpha_df = alpha_barcode_df.sort_values("umis_a", ascending=False)
        to_drop_key = sorted_alpha_df.iloc[1]["tmp_clonotype_alpha_nt"]
        alpha_drop = alpha_clonotype_aa_to_indices[to_drop_key]
        alpha_barcode_df = sorted_alpha_df.iloc[:1].copy()

    beta_umi_ratio = compute_umi_ratio(beta_barcode_df, "umis_b")
    if umi_ratio_threshold and beta_umi_ratio < umi_ratio_threshold:
        sorted_beta_df = beta_barcode_df.sort_values("umis_b", ascending=False)
        to_drop_key = sorted_beta_df.iloc[1]["tmp_clonotype_beta_nt"]
        beta_drop = beta_clonotype_aa_to_indices[to_drop_key]
        beta_barcode_df = sorted_beta_df.iloc[:1].copy()

    alpha_barcode_df = alpha_barcode_df.sort_values("umis_a", ascending=False)
    beta_barcode_df = beta_barcode_df.sort_values("umis_b", ascending=False)

    clonotype_alpha_nt, clonotype_alpha_aa = construct_clonotype_string(alpha_barcode_df, "TRA", "cdr3_alpha_nt", "cdr3_alpha_aa")
    clonotype_beta_nt, clonotype_beta_aa = construct_clonotype_string(beta_barcode_df, "TRB", "cdr3_beta_nt", "cdr3_beta_aa")

    return {
        "barcode": barcode,
        "alpha_umi_ratio": alpha_umi_ratio,
        "beta_umi_ratio": beta_umi_ratio,
        "alpha_drop": alpha_drop,
        "beta_drop": beta_drop,
        "clonotype_nt": clonotype_alpha_nt + "__" + clonotype_beta_nt,
        "clonotype_aa": clonotype_alpha_aa + "__" + clonotype_beta_aa,
    }


def process_dual_chain_barcode(barcode: str, dual_chain_contig_df: pd.DataFrame, umi_col: str, umi_ratio_threshold: float = None):
    barcode_df = dual_chain_contig_df.loc[dual_chain_contig_df.loc[:, "barcode"] == barcode, :].copy()
    umi_ratio = compute_umi_ratio(barcode_df, umi_col)
    drop_index = None
    if umi_ratio_threshold and umi_ratio < umi_ratio_threshold:
        sorted_barcode_df = barcode_df.sort_values(umi_col, ascending=False)
        drop_index = sorted_barcode_df.index[1]
        barcode_df = sorted_barcode_df.iloc[:1].copy()

    barcode_df = barcode_df.sort_values(umi_col, ascending=False)

    if umi_col.endswith("_a"):
        clonotype_alpha_nt, clonotype_alpha_aa = construct_clonotype_string(barcode_df, "TRA", "cdr3_alpha_nt", "cdr3_alpha_aa")
        if len(barcode_df.loc[:, "cdr3_beta_nt"].unique()) > 1:
            raise Exception("Something went wrong with extracting dual_alpha chain barcodes")
        clonotype_beta_nt = barcode_df.loc[:, "TRBV"].unique()[0] + "_" + barcode_df.loc[:, "TRBJ"].unique()[0] + "_" + barcode_df.loc[:, "cdr3_beta_nt"].unique()[0]
        clonotype_beta_aa = barcode_df.loc[:, "TRBV"].unique()[0] + "_" + barcode_df.loc[:, "TRBJ"].unique()[0] + "_" + barcode_df.loc[:, "cdr3_beta_aa"].unique()[0]
    elif umi_col.endswith("_b"):
        clonotype_beta_nt, clonotype_beta_aa = construct_clonotype_string(barcode_df, "TRB", "cdr3_beta_nt", "cdr3_beta_aa")
        if len(barcode_df.loc[:, "cdr3_alpha_nt"].unique()) > 1:
            raise Exception("Something went wrong with extracting dual_beta chain barcodes")
        clonotype_alpha_nt = barcode_df.loc[:, "TRAV"].unique()[0] + "_" + barcode_df.loc[:, "TRAJ"].unique()[0] + "_" + barcode_df.loc[:, "cdr3_alpha_nt"].unique()[0]
        clonotype_alpha_aa = barcode_df.loc[:, "TRAV"].unique()[0] + "_" + barcode_df.loc[:, "TRAJ"].unique()[0] + "_" + barcode_df.loc[:, "cdr3_alpha_aa"].unique()[0]
    else:
        raise Exception("umi_col name needs to end with chain suffix:", umi_col)

    return {
        "barcode": barcode,
        "umi_ratio": umi_ratio,
        "drop_index": drop_index,
        "clonotype_nt": clonotype_alpha_nt + "__" + clonotype_beta_nt,
        "clonotype_aa": clonotype_alpha_aa + "__" + clonotype_beta_aa,
    }


def process_dual_alpha_and_dual_beta_barcode_list(
    contig_df: pd.DataFrame, barcodes_list: list, umi_ratio_threshold: float, alpha_umi_ratio_list: list, beta_umi_ratio_list: list, counter: collections.Counter, threads: int = 8
):
    dual_alpha_and_dual_beta_contig_df = contig_df.loc[contig_df.loc[:, "barcode"].isin(barcodes_list), :].copy()
    results = []

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(process_dual_alpha_and_dual_beta_barcode, barcode, dual_alpha_and_dual_beta_contig_df, umi_ratio_threshold) for barcode in barcodes_list]
        for f in as_completed(futures):
            results.append(f.result())

    for res in results:
        alpha_umi_ratio_list.append(res["alpha_umi_ratio"])
        beta_umi_ratio_list.append(res["beta_umi_ratio"])

        drop_indices = set()
        if res["alpha_drop"] is not None:
            drop_indices.update(res["alpha_drop"])
            counter.update(["umi_count_ratio_dual_alpha_removed_from_dual_alpha_and_dual_beta"])
        if res["beta_drop"] is not None:
            drop_indices.update(res["beta_drop"])
            counter.update(["umi_count_ratio_dual_beta_removed_from_dual_alpha_and_dual_beta"])

        if drop_indices:
            contig_df.drop(index=list(drop_indices), inplace=True)

        if res["alpha_drop"] is not None and res["beta_drop"] is not None:
            continue

        elif res["alpha_drop"] is None and res["beta_drop"] is None:
            contig_df.loc[contig_df["barcode"] == res["barcode"], ["clonotype_nt", "clonotype_aa", "detected_chains"]] = (
                res["clonotype_nt"],
                res["clonotype_aa"],
                "dual_alpha_and_dual_beta",
            )
            continue
        elif res["alpha_drop"] is None and res["beta_drop"] is not None:
            contig_df.loc[contig_df["barcode"] == res["barcode"], ["clonotype_nt", "clonotype_aa", "detected_chains"]] = (res["clonotype_nt"], res["clonotype_aa"], "dual_alpha")
            continue
        elif res["alpha_drop"] is not None and res["beta_drop"] is None:
            contig_df.loc[contig_df["barcode"] == res["barcode"], ["clonotype_nt", "clonotype_aa", "detected_chains"]] = (res["clonotype_nt"], res["clonotype_aa"], "dual_beta")
            continue
        else:
            raise Exception("Failed to test None status of res['alpha_drop'] or res['beta_drop'].")

    return contig_df, counter


def process_dual_chain_barcode_list(
    contig_df: pd.DataFrame,
    barcodes_list: list,
    umi_col: str,
    umi_ratio_list: list,
    counter: collections.Counter,
    dual_label: str,
    drop_label: str,
    umi_ratio_threshold: float = None,
    threads: int = 8,
):
    dual_chain_contig_df = contig_df.loc[contig_df.loc[:, "barcode"].isin(barcodes_list)].copy()
    results = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        future_to_barcode = {
            executor.submit(
                process_dual_chain_barcode,
                barcode,
                dual_chain_contig_df,
                umi_col,
                umi_ratio_threshold
            ): barcode
            for barcode in barcodes_list
        }

        for f in as_completed(future_to_barcode):
            barcode = future_to_barcode[f]

            try:
                results.append(f.result())

            except Exception as e:
                raise RuntimeError(f"Failed barcode: {barcode}") from e

    dropped_idx_list = []
    for res in results:
        umi_ratio_list.append(res["umi_ratio"])
        if res["drop_index"] is not None:
            dropped_idx_list.append(dropped_idx_list)
            contig_df.drop(res["drop_index"], inplace=True)
            counter.update(list(drop_label))
            continue

        contig_df.loc[contig_df["barcode"] == res["barcode"], ["clonotype_nt", "clonotype_aa", "dual_chain_umi_count_ratio", "detected_chains"]] = (
            res["clonotype_nt"],
            res["clonotype_aa"],
            res["umi_ratio"],
            dual_label,
        )

    return contig_df, counter


def parse_10x_from_contig(
    contig_csv: Union[pd.DataFrame, str, os.PathLike[str]],
    contig_csv_name: str,
    log_file: Union[str, os.PathLike[str]],
    meta_data_column_list: list = None,
    filter_on_quality: bool = True,
    keep_highest_read_only: bool = False,
    remove_dual_alpha_and_dual_beta_barcodes: bool = False,
    remove_dual_beta_barcodes: bool = False,
    remove_dual_alpha_barcodes: bool = False,
    remove_single_chain_barcodes: bool = True,
    remove_only_alpha_or_beta_chain_barcodes: bool = True,
    dual_chain_umi_count_ratio_threshold: float = None,
    threads: int = 8,
):
    """Parse a 10x ``filtered_contig_annotations.csv`` to call and count TCR alpha + beta clonotypes.

    Briefly: contigs are filtered on Cell Ranger's ``high_confidence``, ``is_cell``,
    ``full_length``, and ``productive`` flags. Barcodes with only a single chain, or with
    multiple chains that are all alpha or all beta but not both, are removed. Dual TCR alpha and
    dual TCR beta chains are filtered per barcode, separately for the alpha and beta chains, based
    on the ratio of the non-top chain's UMI count to the top chain's UMI count
    (``dual_chain_umi_count_ratio_threshold``). Dual chains that pass this filter are ordered
    deterministically by UMI count when building the clonotype string, to prevent incorrectly
    counting dual chain clonotypes twice. Amino acid and nucleotide clonotypes are defined by
    TRAV/TRAJ/TRBV/TRBJ gene usage and the CDR3alpha/CDR3beta amino acid or nucleotide sequences.
    To support screening cells expressing dual alpha and/or dual beta chains, each such barcode is
    expanded into the 2 (one dual chain) or 4 (two dual chains) possible alpha+beta combinations
    (unique_clonotype) as separate rows, each carrying the same parent dual-chain
    ``clonotype_nt``/``clonotype_aa`` string (clonotype counts and frequencies are computed
    later, by :func:`~tcr_toolbox.tcr_parsing.parser_utils.calculate_clonotype_frequencies`).
    See the example notebook for the full column reference (``clonotype_nt``, ``clonotype_aa``,
    ``detected_chains``,``dual_chain_umi_count_ratio``, ``unique_clonotype_nt``, ``unique_clonotype_aa``,
    and their alpha/beta splits).

    Parameters
    ----------
    contig_csv : DataFrame or path
        A ``filtered_contig_annotations.csv`` path, or an already-loaded dataframe.
    contig_csv_name : str
        Sample/run identifier, stored in the ``data_origin`` column.
    log_file : path
        Path to append parsing log lines to.
    meta_data_column_list : list, optional
        Extra columns from ``contig_csv`` to carry through per chain, suffixed ``_a``/``_b``.
    filter_on_quality : bool, default=True
        Restrict to rows where ``high_confidence``, ``is_cell``, ``full_length``, and
        ``productive`` are all true.
    keep_highest_read_only : bool, default=False
        Keep only the highest-``reads`` contig per barcode per chain.
    remove_dual_alpha_and_dual_beta_barcodes : bool, default=False
        Drop dual_alpha_and_dual_beta barcodes instead of resolving/keeping them.
    remove_dual_beta_barcodes : bool, default=False
        Drop dual_beta barcodes instead of resolving/keeping them.
    remove_dual_alpha_barcodes : bool, default=False
        Drop dual_alpha barcodes instead of resolving/keeping them.
    remove_single_chain_barcodes : bool, default=True
        Drop barcodes with only 1 detected chain.
    remove_only_alpha_or_beta_chain_barcodes : bool, default=True
        Drop barcodes where all detected chains are the same type.
    dual_chain_umi_count_ratio_threshold : float, optional
        Ratio of the non-top chain's UMI count to the top chain's UMI count below which the
        non-top chain is filtered. For example, a dual_alpha cell barcode then becomes a
        single_pair barcode (see Notes).
    threads : int, default=8
        Worker processes used when resolving dual-chain barcodes.

    Returns
    -------
    DataFrame or tuple
        ``contig_df``, or - if any dual-chain barcodes were kept - a ``(contig_df,
        alpha_dual_chain_umi_count_ratio_list, beta_dual_chain_umi_count_ratio_list)`` tuple.

    Notes
    -----
    Barcode classification, counted in ``log_file``:

    - ``single_chain``: barcode has only 1 detected chain (whether alpha or beta).
    - ``only_alpha_or_beta_chains``: barcode has only a single chain type detected.
    - ``single_pair`` (not counted under its own key; reflected in the final
      ``detected_chains`` distribution): barcode has 1 detected alpha and beta pair.
    - ``3 rows``: barcode has 3 detected chains.
    - ``dual_alpha``: of those 3 chains, 2 are alpha and 1 is beta.
    - ``dual_beta``: of those 3 chains, 2 are beta and 1 is alpha.
    - ``nt_duplicate_dual_alpha_or_nt_duplicate_dual_beta``: of those 3 chains, the seemingly
      duplicated alpha or beta rows are identical nt clonotypes rather than a true dual pair
      (so they collapse to 1 unique alpha/beta clonotype, not 2). The barcode cannot be
      resolved into a dual_alpha or dual_beta pair, so it is dropped from the output.
    - ``4 rows``: barcode has 4 detected chains.
    - ``dual_alpha_and_dual_beta``: of those 4 chains, 2 are alpha and 2 are beta.
    - ``4_rows_due_to_3_alpha_or_beta``: of those 4 chains, not exactly 2 alpha + 2 beta (e.g.
      3 alpha + 1 beta). The barcode cannot be resolved into a dual_alpha_and_dual_beta pair,
      so it is dropped from the output.

    If ``dual_chain_umi_count_ratio_threshold`` is set, the non-top chain's UMI count is
    removed when its ratio to the top chain's UMI count falls below that threshold
    (``umi_count_ratio_dual_alpha_removed`` / ``umi_count_ratio_dual_beta_removed``),
    converting dual_alpha/dual_beta to single_pair, and dual_alpha_and_dual_beta to
    single_pair/dual_alpha/dual_beta depending on which chain(s) were removed
    (``..._removed_from_dual_alpha_and_dual_beta``). Dual chains that pass the filter are
    ordered by UMI count, descending, before their V/J/CDR3 annotations are joined into the
    clonotype string, so the same underlying dual-chain pair always yields the same
    ``clonotype_nt``/``clonotype_aa`` regardless of contig row order.
    """
    if not remove_dual_alpha_and_dual_beta_barcodes and (remove_dual_beta_barcodes or remove_dual_alpha_barcodes):
        raise Exception("It is not allowed to set remove_dual_alpha_and_dual_beta_barcodes to False and remove_dual_beta_barcodes or remove_dual_alpha_barcodes to True!")

    if isinstance(contig_csv, pd.DataFrame):
        contig_df = contig_csv.copy()
    else:
        contig_df = pd.read_csv(contig_csv)

    if filter_on_quality:
        contig_df = contig_df.loc[contig_df.loc[:, "high_confidence"], :].copy()
        contig_df = contig_df.loc[contig_df.loc[:, "is_cell"], :].copy()
        contig_df = contig_df.loc[contig_df.loc[:, "full_length"], :].copy()
        contig_df = contig_df.loc[contig_df.loc[:, "productive"], :].copy()

    ambiguous_counter = collections.Counter()
    remove_barcodes_list = []

    # check number of rows per barcode
    vc = contig_df.loc[:, "barcode"].value_counts()
    # take only barcodes with max 4 rows
    contig_df = contig_df.loc[contig_df.loc[:, "barcode"].isin(vc[vc <= 4].index), :].copy()

    dual_alpha_barcode_list = []
    dual_beta_barcode_list = []
    dual_alpha_and_dual_beta_barcode_list = []

    for barcode, barcode_df in contig_df.groupby("barcode"):
        barcode_df = barcode_df.loc[barcode_df.loc[:, "chain"].isin({"TRA", "TRB"}), :]
        if barcode_df.empty:
            ambiguous_counter.update(["no_TRA_or_TRB_chains"])
            remove_barcodes_list.append(barcode)
            continue

        chains = barcode_df.loc[:, "chain"]
        unique_chains = chains.unique()
        n_rows = barcode_df.shape[0]

        barcode_df.loc[:, "tmp_clonotype_nt"] = barcode_df.loc[:, "v_gene"] + "_" + barcode_df.loc[:, "j_gene"] + "_" + barcode_df.loc[:, "cdr3_nt"]
        alpha_chain_count = barcode_df.loc[barcode_df.loc[:, "chain"] == "TRA", "tmp_clonotype_nt"].value_counts()
        beta_chain_count = barcode_df.loc[barcode_df.loc[:, "chain"] == "TRB", "tmp_clonotype_nt"].value_counts()

        if n_rows == 1:
            ambiguous_counter.update(["single_chain"])
            if remove_single_chain_barcodes:
                remove_barcodes_list.append(barcode)
            continue

        if len(unique_chains) == 1:
            ambiguous_counter.update(["only_alpha_or_beta_chains"])
            if remove_only_alpha_or_beta_chain_barcodes:
                remove_barcodes_list.append(barcode)
            continue

        if n_rows == 4:
            ambiguous_counter.update(["4 rows"])
            if alpha_chain_count.shape[0] == 2 and beta_chain_count.shape[0] == 2:
                ambiguous_counter.update(["dual_alpha_and_dual_beta"])
                if remove_dual_alpha_and_dual_beta_barcodes:
                    remove_barcodes_list.append(barcode)
                else:
                    dual_alpha_and_dual_beta_barcode_list.append(barcode)
            else:
                ambiguous_counter.update(["4_rows_due_to_3_alpha_or_beta"])
                remove_barcodes_list.append(barcode)
            continue

        if n_rows == 3:
            ambiguous_counter.update(["3 rows"])
            if beta_chain_count.shape[0] == 2:
                ambiguous_counter.update(["dual_beta"])
                if remove_dual_beta_barcodes:
                    remove_barcodes_list.append(barcode)
                else:
                    dual_beta_barcode_list.append(barcode)
            elif alpha_chain_count.shape[0] == 2:
                ambiguous_counter.update(["dual_alpha"])
                if remove_dual_alpha_barcodes:
                    remove_barcodes_list.append(barcode)
                else:
                    dual_alpha_barcode_list.append(barcode)
            else:
                ambiguous_counter.update(["nt_duplicate_dual_alpha_or_nt_duplicate_dual_beta"])
                remove_barcodes_list.append(barcode)

    # for each barcode if only an alpha or beta present or 2 beta, remove entry
    if (
        remove_dual_alpha_and_dual_beta_barcodes
        or remove_dual_alpha_barcodes
        or remove_dual_beta_barcodes
        or remove_single_chain_barcodes
        or remove_only_alpha_or_beta_chain_barcodes
    ):
        removing_parts = []

        if remove_dual_alpha_and_dual_beta_barcodes:
            removing_parts.append("dual_alpha_and_dual_beta")
        if remove_dual_alpha_barcodes:
            removing_parts.append("dual_alpha")
        if remove_dual_beta_barcodes:
            removing_parts.append("dual_beta")
        if remove_single_chain_barcodes:
            removing_parts.append("single_chain")
        if remove_only_alpha_or_beta_chain_barcodes:
            removing_parts.append("only_alpha_or_beta_chain")

        removing = ", ".join(removing_parts)

        with open(log_file, "a") as log_out:
            print(ambiguous_counter, file=log_out)
            print(contig_df.shape[0], f"before removing {removing} rows", file=log_out)
            contig_df = contig_df.loc[~contig_df.loc[:, "barcode"].isin(remove_barcodes_list), :].copy()
            print(contig_df.shape[0], f"after removing {removing} rows", file=log_out)

    selected_barcode_set = set(contig_df.loc[:, "barcode"].values)

    # split for tra and trb
    alpha_contig_df = contig_df.loc[contig_df.loc[:, "chain"] == "TRA", :].copy()
    # keep only the contig for each barcode with the highest read count
    if keep_highest_read_only:
        alpha_contig_df = alpha_contig_df.sort_values(by="reads", ascending=False)
        alpha_contig_df = alpha_contig_df.drop_duplicates(subset="barcode", keep="first")

    if meta_data_column_list is not None:
        alpha_contig_df = alpha_contig_df.loc[:, ["barcode", "v_gene", "d_gene", "j_gene", "cdr3", "cdr3_nt", "umis"] + meta_data_column_list]
    else:
        alpha_contig_df = alpha_contig_df.loc[:, ["barcode", "v_gene", "d_gene", "j_gene", "cdr3", "cdr3_nt", "umis"]]

    alpha_contig_df = alpha_contig_df.rename(columns={"v_gene": "TRAV", "j_gene": "TRAJ", "cdr3": "cdr3_alpha_aa", "cdr3_nt": "cdr3_alpha_nt"})
    if meta_data_column_list is not None:
        meta_rename_dict_tmp = dict(zip(meta_data_column_list, [meta_col + "_a" for meta_col in meta_data_column_list]))
        alpha_contig_df = alpha_contig_df.rename(columns=meta_rename_dict_tmp)

    beta_contig_df = contig_df.loc[contig_df.loc[:, "chain"] == "TRB", :].copy()
    if keep_highest_read_only:
        beta_contig_df = beta_contig_df.sort_values(by="reads", ascending=False)
        beta_contig_df = beta_contig_df.drop_duplicates(subset="barcode", keep="first")

    if meta_data_column_list is not None:
        beta_contig_df = beta_contig_df.loc[:, ["barcode", "v_gene", "d_gene", "j_gene", "cdr3", "cdr3_nt", "umis"] + meta_data_column_list]
    else:
        beta_contig_df = beta_contig_df.loc[:, ["barcode", "v_gene", "d_gene", "j_gene", "cdr3", "cdr3_nt", "umis"]]

    beta_contig_df = beta_contig_df.rename(columns={"v_gene": "TRBV", "d_gene": "TRBD", "j_gene": "TRBJ", "cdr3": "cdr3_beta_aa", "cdr3_nt": "cdr3_beta_nt"})
    if meta_data_column_list is not None:
        meta_rename_dict_tmp = dict(zip(meta_data_column_list, [meta_col + "_b" for meta_col in meta_data_column_list]))
        beta_contig_df = beta_contig_df.rename(columns=meta_rename_dict_tmp)

    # combine beta and alpha
    contig_df = beta_contig_df.merge(alpha_contig_df, on="barcode", how="outer", suffixes=("_b", "_a"))

    with open(log_file, "a") as log_out:
        print(contig_df.shape[0], "TCRs before dual chain filtering", file=log_out)

    contig_df.loc[:, "data_origin"] = contig_csv_name
    contig_df.loc[:, "clonotype_nt"] = None
    contig_df.loc[:, "clonotype_aa"] = None
    contig_df.loc[:, "detected_chains"] = "single_pair"
    contig_df.loc[:, "dual_chain_umi_count_ratio"] = None

    if not remove_dual_alpha_barcodes or not remove_dual_beta_barcodes or not dual_alpha_and_dual_beta_barcode_list:
        dual_chain_removed_counter = collections.Counter()
        alpha_dual_chain_umi_count_ratio_list = []
        beta_dual_chain_umi_count_ratio_list = []

    if dual_alpha_and_dual_beta_barcode_list and not (remove_dual_alpha_and_dual_beta_barcodes and remove_dual_alpha_barcodes and remove_dual_beta_barcodes):
        contig_df, dual_chain_removed_counter = process_dual_alpha_and_dual_beta_barcode_list(
            contig_df,
            dual_alpha_and_dual_beta_barcode_list,
            dual_chain_umi_count_ratio_threshold,
            alpha_dual_chain_umi_count_ratio_list,
            beta_dual_chain_umi_count_ratio_list,
            dual_chain_removed_counter,
            threads=threads,
        )

    if not remove_dual_alpha_barcodes and dual_alpha_barcode_list:
        contig_df, dual_chain_removed_counter = process_dual_chain_barcode_list(
            contig_df,
            dual_alpha_barcode_list,
            umi_col="umis_a",
            umi_ratio_list=alpha_dual_chain_umi_count_ratio_list,
            counter=dual_chain_removed_counter,
            dual_label="dual_alpha",
            drop_label=["umi_count_ratio_dual_alpha_removed"],
            umi_ratio_threshold=dual_chain_umi_count_ratio_threshold,
            threads=threads,
        )

    if not remove_dual_beta_barcodes and dual_beta_barcode_list:
        contig_df, dual_chain_removed_counter = process_dual_chain_barcode_list(
            contig_df,
            dual_beta_barcode_list,
            umi_col="umis_b",
            umi_ratio_list=beta_dual_chain_umi_count_ratio_list,
            counter=dual_chain_removed_counter,
            dual_label="dual_beta",
            drop_label=["umi_count_ratio_dual_beta_removed"],
            umi_ratio_threshold=dual_chain_umi_count_ratio_threshold,
            threads=threads,
        )

    if not (remove_dual_alpha_barcodes or remove_dual_beta_barcodes or remove_dual_alpha_and_dual_beta_barcodes):
        with open(log_file, "a") as log_out:
            print(dual_chain_removed_counter, file=log_out)

    vc = contig_df.loc[:, "barcode"].value_counts()
    single_pair_contig_df = contig_df.loc[contig_df["barcode"].isin(vc[vc == 1].index), :].copy()
    if not (single_pair_contig_df.loc[:, "detected_chains"].unique() == "single_pair").all():
        raise Exception("Dual chain and single_pair annotations are incorrect or dual chain barcode filtering using dual_chain_umi_count_ratio_threshold went wrong!")

    single_pair_contig_df.loc[:, "clonotype_nt"] = (
        single_pair_contig_df.loc[:, "TRAV"]
        + "_"
        + single_pair_contig_df.loc[:, "TRAJ"]
        + "_"
        + single_pair_contig_df.loc[:, "cdr3_alpha_nt"]
        + "__"
        + single_pair_contig_df.loc[:, "TRBV"]
        + "_"
        + single_pair_contig_df.loc[:, "TRBJ"]
        + "_"
        + single_pair_contig_df.loc[:, "cdr3_beta_nt"]
    )

    single_pair_contig_df.loc[:, "clonotype_aa"] = (
        single_pair_contig_df.loc[:, "TRAV"]
        + "_"
        + single_pair_contig_df.loc[:, "TRAJ"]
        + "_"
        + single_pair_contig_df.loc[:, "cdr3_alpha_aa"]
        + "__"
        + single_pair_contig_df.loc[:, "TRBV"]
        + "_"
        + single_pair_contig_df.loc[:, "TRBJ"]
        + "_"
        + single_pair_contig_df.loc[:, "cdr3_beta_aa"]
    )
    contig_df.update(single_pair_contig_df)
    contig_df.loc[:, "unique_clonotype_nt"] = (
        contig_df.loc[:, "TRAV"]
        + "_"
        + contig_df.loc[:, "TRAJ"]
        + "_"
        + contig_df.loc[:, "cdr3_alpha_nt"]
        + "__"
        + contig_df.loc[:, "TRBV"]
        + "_"
        + contig_df.loc[:, "TRBJ"]
        + "_"
        + contig_df.loc[:, "cdr3_beta_nt"]
    )

    contig_df.loc[:, "unique_clonotype_aa"] = (
        contig_df.loc[:, "TRAV"]
        + "_"
        + contig_df.loc[:, "TRAJ"]
        + "_"
        + contig_df.loc[:, "cdr3_alpha_aa"]
        + "__"
        + contig_df.loc[:, "TRBV"]
        + "_"
        + contig_df.loc[:, "TRBJ"]
        + "_"
        + contig_df.loc[:, "cdr3_beta_aa"]
    )

    contig_df.loc[:, "clonotype_alpha_nt"] = contig_df.loc[:, "clonotype_nt"].str.split("__").str[0]
    contig_df.loc[:, "clonotype_beta_nt"] = contig_df.loc[:, "clonotype_nt"].str.split("__").str[1]
    contig_df.loc[:, "clonotype_alpha_aa"] = contig_df.loc[:, "clonotype_aa"].str.split("__").str[0]
    contig_df.loc[:, "clonotype_beta_aa"] = contig_df.loc[:, "clonotype_aa"].str.split("__").str[1]

    contig_df.loc[:, "unique_clonotype_alpha_nt"] = contig_df.loc[:, "unique_clonotype_nt"].str.split("__").str[0]
    contig_df.loc[:, "unique_clonotype_beta_nt"] = contig_df.loc[:, "unique_clonotype_nt"].str.split("__").str[1]
    contig_df.loc[:, "unique_clonotype_alpha_aa"] = contig_df.loc[:, "unique_clonotype_aa"].str.split("__").str[0]
    contig_df.loc[:, "unique_clonotype_beta_aa"] = contig_df.loc[:, "unique_clonotype_aa"].str.split("__").str[1]

    if not len(selected_barcode_set.intersection(set(contig_df.loc[:, "barcode"].values))) == len(selected_barcode_set):
        raise Exception("contig_df barcodes were lost while generating clonotype columns!")
    contig_df.reset_index(inplace=True, drop=True)
    with open(log_file, "a") as log_out:
        print(f"{contig_df.shape[0]} final selected TCR barcodes:", file=log_out)
        print(f"{contig_df.loc[:, 'detected_chains'].value_counts(normalize=True, dropna=False).to_string(header=False)}\n", file=log_out)

    if not remove_dual_alpha_barcodes or not remove_dual_beta_barcodes:
        return (contig_df, alpha_dual_chain_umi_count_ratio_list, beta_dual_chain_umi_count_ratio_list)
    else:
        return contig_df


def parse_10x_contig_list(
    gem_10x_run_name_list: list,
    contig_file_list: list,
    outs_dir: Union[str, os.PathLike[str]] = None,
    meta_data_column_list: list = None,
    filter_on_quality: bool = True,
    keep_highest_read_only: bool = False,
    remove_dual_alpha_and_dual_beta_barcodes: bool = False,
    remove_dual_beta_barcodes: bool = False,
    remove_dual_alpha_barcodes: bool = False,
    remove_single_chain_barcodes: bool = True,
    remove_only_alpha_or_beta_chain_barcodes: bool = True,
    dual_chain_umi_count_ratio_threshold: float = None,
    threads: int = 8,
) -> pd.DataFrame:
    """Run :func:`parse_10x_from_contig` over multiple 10x runs/samples and combine the results.

    Barcodes are prefixed with ``"{gem_10x_run_name}_"`` to stay unique across samples, then
    all per-sample dataframes are concatenated into ``combined_contig_df``. If any dual-chain
    barcodes were kept, also plots the per-barcode dual-chain UMI count ratio QC histogram to
    ``outs_dir / "umi_count_ratio_least_abundant_most_abundant_dual_chain_hist.pdf"``.

    Parameters
    ----------
    gem_10x_run_name_list : list
        Sample/run identifiers, in the same order as ``contig_file_list``.
    contig_file_list : list
        Paths to each run's ``filtered_contig_annotations.csv``.
    outs_dir : path, default=None
        Output directory; the parsing log and UMI ratio QC histogram are written here.
    meta_data_column_list : list, optional
        Passed through to :func:`parse_10x_from_contig`.
    filter_on_quality : bool, default=True
        Passed through to :func:`parse_10x_from_contig`.
    keep_highest_read_only : bool, default=False
        Passed through to :func:`parse_10x_from_contig`.
    remove_dual_alpha_and_dual_beta_barcodes : bool, default=False
        Passed through to :func:`parse_10x_from_contig`.
    remove_dual_beta_barcodes : bool, default=False
        Passed through to :func:`parse_10x_from_contig`.
    remove_dual_alpha_barcodes : bool, default=False
        Passed through to :func:`parse_10x_from_contig`.
    remove_single_chain_barcodes : bool, default=True
        Passed through to :func:`parse_10x_from_contig`.
    remove_only_alpha_or_beta_chain_barcodes : bool, default=True
        Passed through to :func:`parse_10x_from_contig`.
    dual_chain_umi_count_ratio_threshold : float, optional
        Passed through to :func:`parse_10x_from_contig`; also used as the histogram's
        threshold line.
    threads : int, default=8
        Passed through to :func:`parse_10x_from_contig`.

    Returns
    -------
    DataFrame
        The concatenated ``combined_contig_df`` across all runs, with run-prefixed barcodes.
    """
    log_file = outs_dir / "parse_10x_from_contig.log"

    parsed_contig_df_dict = {}
    alpha_dual_chain_umi_count_ratio_list_dict = {}
    beta_dual_chain_umi_count_ratio_list_dict = {}
    for contig_file, gem_10x_run_name in zip(contig_file_list, gem_10x_run_name_list):
        with open(log_file, "a") as log_out:
            print(gem_10x_run_name, file=log_out)
        contig_df = pd.read_csv(contig_file)
        if not remove_dual_alpha_barcodes or not remove_dual_beta_barcodes:
            (parsed_contig_df_dict[gem_10x_run_name], alpha_dual_chain_umi_count_ratio_list_dict[gem_10x_run_name], beta_dual_chain_umi_count_ratio_list_dict[gem_10x_run_name]) = (
                parse_10x_from_contig(
                    contig_csv=contig_df,
                    contig_csv_name=gem_10x_run_name,
                    log_file=log_file,
                    meta_data_column_list=meta_data_column_list,
                    filter_on_quality=filter_on_quality,
                    keep_highest_read_only=keep_highest_read_only,
                    remove_dual_alpha_and_dual_beta_barcodes=remove_dual_alpha_and_dual_beta_barcodes,
                    remove_dual_beta_barcodes=remove_dual_beta_barcodes,
                    remove_dual_alpha_barcodes=remove_dual_alpha_barcodes,
                    remove_single_chain_barcodes=remove_single_chain_barcodes,
                    remove_only_alpha_or_beta_chain_barcodes=remove_only_alpha_or_beta_chain_barcodes,
                    dual_chain_umi_count_ratio_threshold=dual_chain_umi_count_ratio_threshold,
                    threads=threads,
                )
            )
        else:
            parsed_contig_df_dict[gem_10x_run_name] = parse_10x_from_contig(
                contig_csv=contig_df,
                contig_csv_name=gem_10x_run_name,
                log_file=log_file,
                meta_data_column_list=meta_data_column_list,
                filter_on_quality=filter_on_quality,
                keep_highest_read_only=keep_highest_read_only,
                remove_dual_alpha_and_dual_beta_barcodes=remove_dual_alpha_and_dual_beta_barcodes,
                remove_dual_beta_barcodes=remove_dual_beta_barcodes,
                remove_dual_alpha_barcodes=remove_dual_alpha_barcodes,
                remove_single_chain_barcodes=remove_single_chain_barcodes,
                remove_only_alpha_or_beta_chain_barcodes=remove_only_alpha_or_beta_chain_barcodes,
                dual_chain_umi_count_ratio_threshold=dual_chain_umi_count_ratio_threshold,
                threads=threads,
            )

    for gem_10x_run_name in gem_10x_run_name_list:
        parsed_contig_df_dict[gem_10x_run_name].loc[:, "barcode"] = gem_10x_run_name + "_" + parsed_contig_df_dict[gem_10x_run_name].loc[:, "barcode"]

    combined_contig_df = pd.concat([parsed_contig_df_dict[gem_10x_run_name] for gem_10x_run_name in gem_10x_run_name_list])
    combined_contig_df["TRAV_IMGT"] = combined_contig_df["TRAV"] + "*01"
    combined_contig_df["TRAJ_IMGT"] = combined_contig_df["TRAJ"] + "*01"
    combined_contig_df["TRBV_IMGT"] = combined_contig_df["TRBV"] + "*01"
    combined_contig_df["TRBJ_IMGT"] = combined_contig_df["TRBJ"] + "*01"

    if not remove_dual_alpha_barcodes or not remove_dual_beta_barcodes:
        ax, fig, gs = startfig(7, 5)

        if not remove_dual_alpha_barcodes:
            ax.hist(
                [item for sublist in alpha_dual_chain_umi_count_ratio_list_dict.values() for item in sublist],
                bins=np.arange(0, 1.1, 0.1),
                color="black",
                alpha=0.5,
                label="dual alpha",
                density=True,
            )

        if not remove_dual_beta_barcodes:
            ax.hist(
                [item for sublist in beta_dual_chain_umi_count_ratio_list_dict.values() for item in sublist],
                bins=np.arange(0, 1.1, 0.1),
                color="red",
                alpha=0.5,
                label="dual beta",
                density=True,
            )
        ax.axvline(dual_chain_umi_count_ratio_threshold, label="Filter", linewidth=1)
        ax.set_xlim(0, 1)
        ax.set_xticks(np.arange(0, 1.1, 0.2))
        ax.set_xticklabels([round(xtick, 1) for xtick in np.arange(0, 1.1, 0.2)])
        ax.set_ylabel("Density", fontsize=7)
        ax.set_xlabel("UMI count ratio\n(least abundant / most abundant dual chain)", fontsize=7)
        ax.tick_params("both", labelsize=7)
        ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=7, frameon=False)
        fig.tight_layout()
        fig.savefig(outs_dir / "umi_count_ratio_least_abundant_most_abundant_dual_chain_hist.pdf")
        plt.close()

    return combined_contig_df


def merge_combined_meta_df_with_tcr_contig_dfs(combined_contig_df: pd.DataFrame, combined_meta_df: pd.DataFrame) -> pd.DataFrame:
    """Inner-merge parsed TCR contig data onto cell metadata (e.g. Seurat/Scanpy output), by barcode.

    Cells with no matching barcode in ``combined_contig_df`` (no usable TCR, or filtered out
    during parsing) are dropped, since the merge is an inner join.

    Parameters
    ----------
    combined_contig_df : DataFrame
        Output of :func:`parse_10x_contig_list`. Its ``barcode`` column is prefixed per run as
        ``"{gem_10x_run_name}_{barcode}"``.
    combined_meta_df : DataFrame
        Per-cell metadata dataframe with a ``barcode`` column. Tip: this barcode column must
        carry the same ``"{gem_10x_run_name}_"`` prefix, so if it doesn't already (e.g. raw
        Seurat/Scanpy output), add it yourself before calling this function - otherwise the
        inner join below will silently match nothing and every cell is dropped.

    Returns
    -------
    DataFrame
        ``combined_meta_df`` merged with the matching TCR contig columns, index reset.
    """
    combined_meta_df = combined_meta_df.merge(combined_contig_df, on="barcode", how="inner").copy()

    combined_meta_df.reset_index(inplace=True, drop=True)

    return combined_meta_df


def filter_meta_10x_df_on_tcr_umi_counts(meta_10x_df: pd.DataFrame, sample_col: str = None, umi_count_quantile_threshold: float = 0.15):
    meta_10x_df = meta_10x_df.copy()
    print("Shape combined_meta_df before UMI count quantile filtering:", meta_10x_df.shape[0])

    remove_idx_list = []
    for sample in meta_10x_df.loc[:, sample_col].unique():
        sample_df = meta_10x_df.loc[meta_10x_df.loc[:, sample_col] == sample, :]
        cell_count_sums = sample_df.loc[:, ["umis_a", "umis_b"]].sum(1)
        remove_idx = sample_df.loc[cell_count_sums < cell_count_sums.quantile(umi_count_quantile_threshold), :].index
        remove_idx_list.extend(remove_idx)

    meta_10x_df.drop(remove_idx_list, axis=0, inplace=True)
    print("Shape combined_meta_df after UMI count quantile filtering:", meta_10x_df.shape[0])
    return meta_10x_df
