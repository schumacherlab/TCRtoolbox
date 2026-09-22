# plotting helper functions
import itertools
import os
import re
from pathlib import Path
from typing import Union

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
from numpy.typing import ArrayLike
from scipy.interpolate import interpn


def view_color_dict(color_dict, size=1, label_length_scaler=4, add_hex_code=False, save_path=""):
    palette = list(color_dict.values())
    palette_names = list(color_dict.keys())
    n_patches = len(palette)
    f, ax = plt.subplots(1, 1, figsize=(n_patches * size, size + label_length_scaler))
    heatmap = ax.imshow(np.arange(n_patches).reshape(1, n_patches), cmap=mpl.colors.ListedColormap(list(palette)), interpolation="nearest", aspect="auto")

    if add_hex_code:
        for i, hex_code in enumerate(palette):
            ax.text(
                i, 0, hex_code.upper(), ha="center", va="center", fontsize=8, color="white", bbox=dict(boxstyle="round,pad=0.2", facecolor="black", edgecolor="none", alpha=0.5)
            )

    ax.set_xticks(np.arange(n_patches))
    ax.set_yticks([-0.5, 0.5])
    ax.set_xticklabels([name for name in palette_names], rotation=90)
    ax.yaxis.set_major_locator(mpl.ticker.NullLocator())
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()
    plt.close()


def startfig(w=4, h=2, rows=1, columns=1, wrs=None, hrs=None, frameon=True, return_first_ax=True):
    ratio = 0.393701  # 1 cm in inch
    myfigsize = (w * ratio, h * ratio)
    fig = plt.figure(figsize=(myfigsize))
    gs = mpl.gridspec.GridSpec(rows, columns, width_ratios=wrs, height_ratios=hrs)
    if return_first_ax == True:
        a = fig.add_subplot(gs[0, 0], frameon=frameon)
        return a, fig, gs
    else:
        return fig, gs


def colorize_errorbars(err, colors, ax=None):
    ax = ax or plt.gca()
    scs = []
    for line in err.lines[1]:
        line.remove()
        sc = ax.scatter(*line.get_data(), c=colors, marker=line.get_marker(), s=line.get_markersize() ** 2, zorder=3)
        scs.append(sc)
    err.lines = (err.lines[0], tuple(scs), err.lines[2])
    for col in err.lines[2]:
        col.set_color(colors)


def bootstrap_confidence_interval(preds, true, bootstrap_iterations, metric_formula, beta=None):
    from sklearn.utils import resample

    print(f"starting bootstrap {metric_formula.__name__}")

    random_state = 42  # use for reproducibility
    metric_bootstrap_list = []
    bootstrap_samples = len(preds)
    df_for_bootstrap = pd.DataFrame(data={"preds": preds, "true": true})

    for boot in range(bootstrap_iterations):
        bootstrap_df = resample(df_for_bootstrap, replace=True, n_samples=bootstrap_samples, random_state=random_state)
        bootstrap_preds = bootstrap_df["preds"]
        bootstrap_true = bootstrap_df["true"]

        if beta is not None:
            bootstrap_metric = metric_formula(bootstrap_true, bootstrap_preds, beta=beta)
        else:
            bootstrap_metric = metric_formula(bootstrap_true, bootstrap_preds)
        metric_bootstrap_list.append(bootstrap_metric)
        random_state += 1

    return metric_bootstrap_list


def set_mpl_params(mpl):
    mylines = 0.15 * 2.82
    # Matplotlib plotting params:
    mpl.style.use("classic")
    mpl.rcParams["axes.linewidth"] = mylines
    mpl.rcParams["ytick.direction"] = "out"
    mpl.rcParams["xtick.direction"] = "out"
    mpl.rcParams["xtick.major.size"] = 2
    mpl.rcParams["ytick.major.size"] = 2
    mpl.rcParams["xtick.major.width"] = mylines
    mpl.rcParams["ytick.major.width"] = mylines
    mpl.rcParams["grid.linewidth"] = mylines / 1.5
    mpl.rcParams["grid.color"] = "0.8"
    mpl.rcParams["grid.linestyle"] = "solid"
    mpl.rcParams["legend.frameon"] = False
    mpl.rcParams["font.family"] = "Arial"  # added
    mpl.rcParams["figure.dpi"] = 300
    mpl.rc("savefig", dpi=300)
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    return mpl, mylines


def set_sns_params(sns):
    mylines = 0.15 * 2.82
    # seaborn plotting params:
    sns.set_style(
        "white",
        rc={
            "axes.linewidth": mylines,
            "ytick.direction": "out",
            "xtick.direction": "out",
            "xtick.major.size": 2,
            "ytick.major.size": 2,
            "xtick.major.width": mylines,
            "ytick.major.width": mylines,
            "grid.linewidth": mylines / 1.5,
            "grid.color": "0.8",
            "grid.linestyle": "solid",
            "legend.frameon": False,
            "figure.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 300,
        },
    )
    return sns, mylines


def plot_v_genes_per_well(
    plate_df: pd.DataFrame,
    save_path: Union[str, os.PathLike[str]],
    plate_title: str,
    well_384_or_96: int = 96,
    TRAV_col_name: str = "TRAV",
    TRBV_col_name: str = "TRBV",
    simulation_mode: bool = False,
):
    for column in [TRAV_col_name, TRBV_col_name, "well_coord"]:
        if column not in plate_df.columns:
            raise ValueError(f"column {column} was not in plate_df.columns when plotting v genes per well")

    df = plate_df.copy()
    fig = plt.figure(figsize=(12, 8))
    for i, row in df.iterrows():
        trav = str(df.loc[i, TRAV_col_name])
        trbv = str(df.loc[i, TRBV_col_name])
        well = df.loc[i, "well_coord"]
        if well_384_or_96 == 96:
            plt.text(well[1] - 0.15, well[0], trav + "\n" + trbv, ha="left", va="center", color="b", fontsize=8)
        if well_384_or_96 == 384:
            plt.text(well[1] - 0.15, well[0], trav + "\n" + trbv, ha="left", va="center", color="b", fontsize=3)

    if well_384_or_96 == 96:
        # name x axis and ticks
        plt.xlabel("column")
        plt.xticks(np.arange(1, 13, 1))
        # name y axis and ticks
        plt.ylabel("row")
        plt.yticks(np.arange(1, 10, 1), ["A", "B", "C", "D", "E", "F", "G", "H", ""])
    elif well_384_or_96 == 384:
        # name x axis and ticks
        plt.xlabel("column")
        plt.xticks(np.arange(1, 25, 1))
        # name y axis and ticks
        plt.ylabel("row")
        plt.yticks(np.arange(1, 18, 1), ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", ""])
    else:
        raise ValueError("well_384_or_96 must be either 96 or 384")

    # set title
    plt.title(plate_title)
    #  # reverse y axis labels
    plt.gca().invert_yaxis()
    # remove outline of plot right and top
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["top"].set_visible(False)

    if simulation_mode:
        print("simulation mode: not saving plate map plots")
    elif save_path:
        plt.savefig(save_path, dpi=300)
        plt.close()
    else:
        plt.show()
        plt.close()


def calculate_positional_entropy(aligned_sequences, gap_token):
    # poly_position_dict_all = {}
    poly_position_dict_no_gaps = {}
    all_aa_set = set()

    # check if all sequences are the same length
    first_seq_len = len(aligned_sequences[0])
    for seq in aligned_sequences:
        if len(seq) != first_seq_len:
            raise ValueError("Not all sequences are the same length")

    # create a dictionary of all amino acids at each position
    for position in range(first_seq_len):
        current_pos_aa_list = list()
        for seq in aligned_sequences:
            current_pos_aa_list.append(seq[position] if seq[position] != gap_token else "X")
            all_aa_set.add(seq[position])
        # poly_position_dict_all[position] = current_aa_set
        if len(set(current_pos_aa_list)) == 1 and "-" in current_pos_aa_list:
            pass
        else:
            poly_position_dict_no_gaps[position] = current_pos_aa_list

    # replace gap_token with 'X' in all_aa_set
    if gap_token in all_aa_set:
        if "X" not in all_aa_set:
            all_aa_set = {acid if acid != gap_token else "X" for acid in all_aa_set}
        else:
            raise ValueError("X is used as gap token for plotting but is already in all_aa_set")
    # sort all_aa_set and have gap be the last element
    all_aa_set = sorted(all_aa_set)
    if "X" in all_aa_set:
        all_aa_set.remove("X")
    all_aa_set.append("X")

    columns = {acid: [] for acid in all_aa_set}
    poly_position_df_no_gaps = pd.DataFrame(columns)
    # Calculate the fraction of each amino acid at each position and fill the dataframe
    for position, amino_acid_list in poly_position_dict_no_gaps.items():
        total_length = len(amino_acid_list)
        # replace '-' with X
        if gap_token in amino_acid_list:
            amino_acid_list = [acid if acid != gap_token else "X" for acid in amino_acid_list]

        fractions = {acid: amino_acid_list.count(acid) / total_length for acid in all_aa_set}
        # make all fraction be of type float
        fractions = {acid: float(fraction) for acid, fraction in fractions.items()}
        poly_position_df_no_gaps.loc[position] = fractions

    return poly_position_df_no_gaps


def plot_sequence_logo(aligned_sequences, gap_token, title, out_path, plot_gaps=False, color_scheme="skylign_protein"):
    from logomaker import Logo

    poly_position_df_no_gaps = calculate_positional_entropy(aligned_sequences, gap_token)

    if plot_gaps:
        plot_df = poly_position_df_no_gaps
    else:
        plot_df = poly_position_df_no_gaps.reset_index(drop=True)

    # create a Logo object
    logo = Logo(plot_df, color_scheme=color_scheme, vpad=0.1, width=0.8)
    # style using LogoOptions
    logo.style_glyphs(flip_below=True, flip_above=True)
    logo.style_spines(visible=False)
    logo.style_spines(spines=("left", "bottom"), visible=True)
    logo.style_xticks(rotation=90, fmt="%d", anchor=0)
    # logo.style_yticks(spacing=2)
    logo.ax.set_ylabel("bits")
    logo.ax.set_xlabel("position")
    logo.ax.xaxis.set_ticks_position("none")
    logo.ax.yaxis.set_ticks_position("none")
    logo.ax.set_title(title)
    logo.fig.set_figwidth(len(plot_df) / 2)
    # display as png
    logo.fig.savefig(f"{out_path}", dpi=300, transparent=True, bbox_inches="tight")


def plot_count_histogram(
    count_file: Union[str, os.PathLike[str]],
    count_file_ref_name_col: str,
    count_file_count_col: str,
    bin_min: int,
    bin_max: int,
    bin_size: int,
    xlabel: str,
    title: str,
    ylim_max: int,
    log10: bool = False,
    min_count_for_95th_5th_ratio: float = None,
    save_file_path: Union[str, os.PathLike[str]] = None,
    w: float = 6.25,
    h: float = 6.25,
):
    """Plot a histogram of per-reference counts from `count_file`, and print a 95th/5th percentile ratio as a library-uniformity QC metric.

    Parameters
    ----------
    min_count_for_95th_5th_ratio : float, optional
        Only counts above this threshold are included when computing the 95th/5th percentile
        ratio that gets printed (and drawn as an orange vertical line). This threshold is compared
        directly against the plotted series, so if `log10=True` it must be given in log10 space
        (e.g. `3` means "raw count > 1000"), not in raw counts.
    log10 : bool, default False
        If True, plot and filter on log10(count) instead of the raw count.
    """
    # Only python engine allows regex separators:
    count_df = pd.read_csv(count_file, sep=r"\t|,", engine="python", usecols=[count_file_ref_name_col, count_file_count_col]).copy()

    if count_df[count_file_ref_name_col].duplicated().any():
        raise Exception("There are duplicate reference names in the count file!")

    raw_counts_series = count_df[count_file_count_col].copy()
    counts_series = raw_counts_series.copy()

    if log10:
        n_nonpositive = (raw_counts_series <= 0).sum()
        if n_nonpositive:
            raise Exception(f"{n_nonpositive} row(s) have count <= 0 and cannot be log10-transformed!")
        counts_series = np.log10(raw_counts_series)

    if min_count_for_95th_5th_ratio:
        mask = counts_series > min_count_for_95th_5th_ratio
    else:
        mask = pd.Series(True, index=counts_series.index)

    percentile_counts_series = raw_counts_series[mask]  # always raw, regardless of log10

    p95 = np.percentile(percentile_counts_series, 95)
    p5 = np.percentile(percentile_counts_series, 5)
    ratio = p95 / p5 if p5 != 0 else np.nan
    print(f"95th percentile: {p95:.2f} | 5th percentile: {p5:.2f} | 95th/5th ratio: {ratio:.2f}")

    if bin_max < counts_series.max():
        print("Warning: bin_max is smaller than maximal value in count file:", counts_series.max())

    ax, fig, gs = startfig(w, h)
    ax.hist(counts_series, bins=np.arange(bin_min, bin_max, bin_size), color="grey", edgecolor="none", zorder=4)
    ax.axvline(np.median(counts_series), color="red", linewidth=0.75, zorder=6)
    if min_count_for_95th_5th_ratio:
        ax.axvline(min_count_for_95th_5th_ratio, color="orange", linewidth=0.75, zorder=6)

    ax.set_xlim(0, bin_max)
    ax.set_ylim(0, ylim_max)
    ax.set_xticks(ax.get_xticks())
    ax.set_yticks(ax.get_yticks())
    ax.set_xticklabels([round(x) for x in ax.get_xticks()], fontsize=7, rotation=90)
    ax.set_yticklabels([round(y) for y in ax.get_yticks()], fontsize=7)
    if log10:
        ax.set_xlabel(f"{xlabel} (log10)", fontsize=7)
    else:
        ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel("Frequency", fontsize=7)
    ax.set_title(title, fontsize=7)
    ax.tick_params(top=False, right=False)
    fig.tight_layout()

    if save_file_path:
        fig.savefig(save_file_path)
        plt.close()
    else:
        plt.show()
        plt.close()


def plot_10x_tcr_umi_counts_hist(
    combined_meta_df: pd.DataFrame,
    save_dir: Union[str, os.PathLike[str]],
    fig_name: str,
    bin_max: int = 100,
    bin_size: int = 1,
    ylim_max: float = 0.07,
    umi_count_quantile_threshold: float = None,
):
    combined_meta_df = combined_meta_df.copy()
    combined_meta_df.loc[:, "umis_a_b_sum"] = combined_meta_df.loc[:, ["umis_a", "umis_b"]].sum(1)

    if combined_meta_df.loc[:, "umis_a"].max() > bin_max:
        print("Warning: bin_max is smaller than maximal UMI count alpha:", combined_meta_df.loc[:, "umis_a"].max())

    if combined_meta_df.loc[:, "umis_b"].max() > bin_max:
        print("Warning: bin_max is smaller than maximal UMI count beta:", combined_meta_df.loc[:, "umis_b"].max())

    if combined_meta_df.loc[:, "umis_a_b_sum"].max() > bin_max:
        print("Warning: bin_max is smaller than maximal UMI count alpha + beta chain sums:", combined_meta_df.loc[:, "umis_a_b_sum"].max())

    ax, fig, gs = startfig(11, 7)
    ax.hist(combined_meta_df.loc[:, "umis_a"], color="black", density=True, bins=np.arange(0, bin_max, bin_size), alpha=0.3, label="UMI counts alpha chain")
    ax.hist(combined_meta_df.loc[:, "umis_b"], color="red", density=True, bins=np.arange(0, bin_max, bin_size), alpha=0.3, label="UMI counts beta chain")

    ax.set_ylim(0, ylim_max)
    ax.set_xlim(0, bin_max + bin_size)
    ax.set_ylabel("Density", fontsize=8)
    ax.set_xlabel("UMI counts", fontsize=8)
    ax.legend(fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, fig_name + "_umi_a_umi_b_separate_counts_histogram.pdf"))
    plt.close()

    ax, fig, gs = startfig(11, 7)
    ax.hist(combined_meta_df.loc[:, ["umis_a", "umis_b"]].sum(1), color="black", density=True, bins=np.arange(0, bin_max, bin_size), alpha=0.3)

    ax.set_ylim(0, ylim_max)
    ax.set_xlim(0, bin_max + bin_size)
    ax.set_ylabel("Density", fontsize=8)
    ax.set_xlabel("UMI counts alpha + beta chain", fontsize=8)
    if umi_count_quantile_threshold:
        ax.axvline(np.quantile(combined_meta_df.loc[:, "umis_a_b_sum"], umi_count_quantile_threshold), color="yellow", linewidth=0.5, label="UMI count\nfilter threshold")
        ax.legend(fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    fig.savefig(os.path.join(save_dir, fig_name + "_umi_counts_histogram.pdf"))
    plt.close()


def plot_frac_pair_clonotype_of_most_abundant_pair_clonotype_within_chain_clonotype(
    meta_10x_df: pd.DataFrame,
    save_dir: Union[str, os.PathLike[str]],
    save_prefix: str = None,
    bin_size: float = 0.1,
    ylim_max_alpha: float = None,
    ylim_max_beta: float = None,
    vmax_alpha: float = None,
    vmax_beta: float = None,
):
    meta_10x_df = meta_10x_df.copy()
    alpha_ratio_list = []
    beta_umi_list = []
    beta_umi_non_dupli_alpha_list = []
    beta_ratio_list = []
    alpha_umi_list = []
    alpha_umi_non_dupli_alpha_list = []

    for chain in meta_10x_df.loc[:, "clonotype_alpha_nt"].unique():
        chain_df = meta_10x_df.loc[meta_10x_df.loc[:, "clonotype_alpha_nt"] == chain, :].copy()
        vc = chain_df.loc[:, "clonotype_beta_nt"].value_counts()
        if vc.shape[0] == 1:
            beta_umi_non_dupli_alpha_list.append(chain_df.loc[:, "umis_b"].sum())
            continue
        else:
            for chain_vc, other_chain in zip(vc[1:], vc.index[1:]):
                alpha_ratio_list.append(chain_vc / vc.iloc[0])
                beta_umi_list.append(chain_df.loc[chain_df.loc[:, "clonotype_beta_nt"] == other_chain, "umis_b"].median())

    for chain in meta_10x_df.loc[:, "clonotype_beta_nt"].unique():
        chain_df = meta_10x_df.loc[meta_10x_df.loc[:, "clonotype_beta_nt"] == chain, :].copy()
        vc = chain_df.loc[:, "clonotype_alpha_nt"].value_counts()
        if vc.shape[0] == 1:
            alpha_umi_non_dupli_alpha_list.append(chain_df.loc[:, "umis_a"].sum())
            continue
        else:
            for chain_vc, other_chain in zip(vc[1:], vc.index[1:]):
                beta_ratio_list.append(chain_vc / vc.iloc[0])
                alpha_umi_list.append(chain_df.loc[chain_df.loc[:, "clonotype_alpha_nt"] == other_chain, "umis_a"].median())

    bins = np.arange(0, 1 + bin_size, bin_size)
    alpha_ratio_df = pd.DataFrame({"alpha_ratio": alpha_ratio_list, "beta_umi_counts": beta_umi_list})
    alpha_ratio_df.loc[:, "alpha_ratio_bins"] = pd.cut(alpha_ratio_df.loc[:, "alpha_ratio"], bins=bins)
    average_beta_umi_count_per_bin = alpha_ratio_df.groupby("alpha_ratio_bins")["beta_umi_counts"].median().fillna(0)

    ax, fig, gs = startfig(12, 10)
    counts, plot_bins, patches = ax.hist(alpha_ratio_df.loc[:, "alpha_ratio"], bins=bins)

    if not all(bins == plot_bins):
        raise Exception("Bins do not match between pd.cut and ax.hist!")

    if vmax_alpha:
        if vmax_alpha < average_beta_umi_count_per_bin.max():
            print("Warning: set vmax_alpha is lower than the average_beta_umi_count_per_bin.max()!")
        norm = mpl.colors.Normalize(vmin=0, vmax=vmax_alpha)
    else:
        norm = mpl.colors.Normalize(vmin=0, vmax=average_beta_umi_count_per_bin.max())
    cmap = mpl.colormaps["viridis"]

    for patch, left_bin, right_bin in zip(patches, plot_bins[:-1], plot_bins[1:]):
        interval = pd.Interval(round(left_bin, 6), round(right_bin, 6), closed="right")
        beta_umi_mean = average_beta_umi_count_per_bin.get(interval)
        patch.set_facecolor(cmap(norm(beta_umi_mean)))

    if ylim_max_alpha:
        ax.set_ylim(0, ylim_max_alpha)
    ax.set_ylabel("Frequency")
    ax.set_xlabel("Fraction paired clonotype of most abundant\npaired clonotype within alpha chain clonotype")
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label="Median of median normalized\nUMI counts/beta chain")
    fig.tight_layout()
    if save_prefix:
        fig.savefig(os.path.join(save_dir, save_prefix + "_fraction_pair_clonotype_of_most_abundant_pair_clonotype_within_alpha_chain_clonotype_hist.pdf"))
    else:
        fig.savefig(os.path.join(save_dir, "fraction_pair_clonotype_of_most_abundant_pair_clonotype_within_alpha_chain_clonotype_hist.pdf"))
    plt.close()

    beta_ratio_df = pd.DataFrame({"beta_ratio": beta_ratio_list, "alpha_umi_counts": alpha_umi_list})
    beta_ratio_df.loc[:, "beta_ratio_bins"] = pd.cut(beta_ratio_df.loc[:, "beta_ratio"], bins=bins, right=True)
    average_alpha_umi_count_per_bin = beta_ratio_df.groupby("beta_ratio_bins")["alpha_umi_counts"].median().fillna(0)

    ax, fig, gs = startfig(12, 10)
    counts, plot_bins, patches = ax.hist(beta_ratio_df.loc[:, "beta_ratio"], bins=bins)

    if not all(bins == plot_bins):
        raise Exception("Bins do not match between pd.cut and ax.hist!")

    if vmax_beta:
        if vmax_beta < average_alpha_umi_count_per_bin.max():
            print("Warning: set vmax_beta is lower than the average_alpha_umi_count_per_bin.max()!")
        norm = mpl.colors.Normalize(vmin=0, vmax=vmax_beta)
    else:
        norm = mpl.colors.Normalize(vmin=0, vmax=average_alpha_umi_count_per_bin.max())
    cmap = mpl.colormaps["viridis"]

    for patch, left_bin, right_bin in zip(patches, plot_bins[:-1], plot_bins[1:]):
        interval = pd.Interval(round(left_bin, 6), round(right_bin, 6), closed="right")
        alpha_umi_mean = average_alpha_umi_count_per_bin.get(interval)
        patch.set_facecolor(cmap(norm(alpha_umi_mean)))

    if ylim_max_beta:
        ax.set_ylim(0, ylim_max_beta)
    ax.set_ylabel("Frequency")
    ax.set_xlabel("Fraction paired clonotype of most abundant\npaired clonotype within beta chain clonotype")
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label="Median of median normalized\nUMI counts/alpha chain")
    fig.tight_layout()
    if save_prefix:
        fig.savefig(os.path.join(save_dir, save_prefix + "_fraction_pair_clonotype_of_most_abundant_pair_clonotype_within_beta_chain_clonotype_hist.pdf"))
    else:
        fig.savefig(os.path.join(save_dir, "fraction_pair_clonotype_of_most_abundant_pair_clonotype_within_beta_chain_clonotype_hist.pdf"))
    plt.close()


# From: https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density/53865762#53865762
def plot_density_scatter(
    x: ArrayLike,
    y: ArrayLike,
    save_file_path: Union[str, os.PathLike[str]],
    sort: bool = True,
    bins: int = 20,
    xlim_min: float = 0,
    ylim_min: float = 0,
    xlim_max: float = None,
    ylim_max: float = None,
    xlabel: str = None,
    ylabel: str = None,
    w: int = 10,
    h: int = 10,
    **kwargs,
):
    """
    Scatter plot colored by 2d histogram
    """
    x = np.array(x)
    y = np.array(y)
    ax, fig, gs = startfig(w, h)
    data, x_e, y_e = np.histogram2d(x, y, bins=bins, density=True)
    z = interpn((0.5 * (x_e[1:] + x_e[:-1]), 0.5 * (y_e[1:] + y_e[:-1])), data, np.vstack([x, y]).T, method="splinef2d", bounds_error=False)

    # To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort:
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter(x, y, c=z, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if xlim_min or xlim_max:
        ax.set_xlim(xlim_min, xlim_max)

    if ylim_min or ylim_max:
        ax.set_ylim(ylim_min, ylim_max)

    norm = mpl.colors.Normalize(vmin=np.min(z), vmax=np.max(z))
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm), ax=ax)
    cbar.ax.set_ylabel("Density")
    fig.tight_layout()
    fig.savefig(save_file_path)
    plt.close()


def plot_clonotype_frequency(
    combined_meta_df: pd.DataFrame,
    donor_col: str,
    outs_dir: Union[str, os.PathLike[str]],
    ylim_max: int,
    clonotype_type: str = "aa",
    clonotype_col: str = "unique_clonotype_aa",
    clonotype_freq_col: str = "unique_clonotype_frequency_aa",
    bin_max: float = 0.07,
    bin_size: float = 0.001,
    coreceptor_col: Union[str, None] = "CD4_or_CD8",
    CD8_str: Union[str, None] = "CD8",
    CD4_str: Union[str, None] = "CD4",
):
    save_dir = os.path.join(outs_dir, "clonotype_frequency")
    os.mkdir(save_dir)

    if clonotype_type not in ["aa", "nt"]:
        raise Exception("Clonotype type is only allowed to be 'nt' or 'aa' and not:", clonotype_type)

    for donor in combined_meta_df.loc[:, donor_col].unique():
        ax, fig, gs = startfig(9, 7)
        donor_df = combined_meta_df.loc[combined_meta_df.loc[:, donor_col] == donor, :].copy()

        if donor_df.loc[:, clonotype_col].duplicated().any():
            raise Exception("There are duplicate clonotypes in:", clonotype_col, "for donor:", donor, "you have to run drop_clonotype_duplicates before running this function!")
        if donor_df.loc[:, clonotype_freq_col].max() > bin_max:
            print("Warning: bin_max is smaller than maximal clonotype frequency:", donor_df.loc[:, clonotype_freq_col].max(), ", of donor:", donor)

        if CD4_str:
            ax.hist(
                donor_df.loc[donor_df.loc[:, coreceptor_col] == CD4_str, clonotype_freq_col],
                np.arange(0, bin_max, bin_size),
                color="black",
                alpha=0.3,
                label=CD4_str,
                density=False,
            )

        if CD8_str:
            ax.hist(
                donor_df.loc[donor_df.loc[:, coreceptor_col] == CD8_str, clonotype_freq_col], np.arange(0, bin_max, bin_size), color="red", alpha=0.3, label=CD8_str, density=False
            )

        if not CD8_str and not CD4_str:
            ax.hist(donor_df.loc[:, clonotype_freq_col], np.arange(0, bin_max, bin_size), color="black", alpha=0.3, label="CD8_and_CD4", density=False)

        ax.set_ylim(0, ylim_max)
        ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
        ax.set_ylabel("Frequency", fontsize=7)
        if clonotype_type == "aa":
            ax.set_xlabel("aa clonotype frequency", fontsize=7)
        elif clonotype_type == "nt":
            ax.set_xlabel("nt clonotype frequency", fontsize=7)

        if CD4_str and CD8_str:
            ax.set_title(
                donor
                + "\nCD8 n = "
                + str(donor_df.loc[donor_df.loc[:, coreceptor_col] == CD8_str, :].shape[0])
                + ", CD4 n = "
                + str(donor_df.loc[donor_df.loc[:, coreceptor_col] == CD4_str, :].shape[0]),
                fontsize=7,
            )
        elif CD4_str:
            ax.set_title(donor + "\nCD4 n = " + str(donor_df.loc[donor_df.loc[:, coreceptor_col] == CD4_str, :].shape[0]), fontsize=7)
        elif CD8_str:
            ax.set_title(donor + "\nCD8 n = " + str(donor_df.loc[donor_df.loc[:, coreceptor_col] == CD8_str, :].shape[0]), fontsize=7)
        else:
            ax.set_title(donor + "\nn = " + str(donor_df.shape[0]), fontsize=7)

        ax.tick_params("both", labelsize=7)
        fig.tight_layout()
        if clonotype_type == "aa":
            fig.savefig(os.path.join(save_dir, donor + "_clonotype_aa_freq_hist.pdf"))
        elif clonotype_type == "nt":
            fig.savefig(os.path.join(save_dir, donor + "_clonotype_nt_freq_hist.pdf"))
        plt.close()

    for donor in combined_meta_df.loc[:, donor_col].unique():
        ylim_max = 5
        ax, fig, gs = startfig(9, 7)
        donor_df = combined_meta_df.loc[combined_meta_df.loc[:, donor_col] == donor, :].copy()

        if CD4_str:
            ax.hist(
                donor_df.loc[donor_df.loc[:, coreceptor_col] == CD4_str, clonotype_freq_col],
                np.arange(0, bin_max, bin_size),
                color="black",
                alpha=0.3,
                label=CD4_str,
                density=False,
            )

        if CD8_str:
            ax.hist(
                donor_df.loc[donor_df.loc[:, coreceptor_col] == CD8_str, clonotype_freq_col], np.arange(0, bin_max, bin_size), color="red", alpha=0.3, label=CD8_str, density=False
            )

        if not CD4_str and not CD8_str:
            ax.hist(donor_df.loc[:, clonotype_freq_col], np.arange(0, bin_max, bin_size), color="black", alpha=0.3, label="CD8_and_CD4", density=False)

        ax.set_ylim(0, ylim_max)
        ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
        ax.set_ylabel("Frequency", fontsize=7)
        if clonotype_type == "aa":
            ax.set_xlabel("aa clonotype frequency", fontsize=7)
        elif clonotype_type == "nt":
            ax.set_xlabel("nt clonotype frequency", fontsize=7)

        if CD4_str and CD8_str:
            ax.set_title(
                donor
                + "\nCD8 n = "
                + str(donor_df.loc[donor_df.loc[:, coreceptor_col] == CD8_str, :].shape[0])
                + ", CD4 n = "
                + str(donor_df.loc[donor_df.loc[:, coreceptor_col] == CD4_str, :].shape[0]),
                fontsize=7,
            )
        elif CD4_str:
            ax.set_title(donor + "\nCD4 n = " + str(donor_df.loc[donor_df.loc[:, coreceptor_col] == CD4_str, :].shape[0]), fontsize=7)
        elif CD8_str:
            ax.set_title(donor + "\nCD8 n = " + str(donor_df.loc[donor_df.loc[:, coreceptor_col] == CD8_str, :].shape[0]), fontsize=7)
        else:
            ax.set_title(donor + "\nn = " + str(donor_df.shape[0]), fontsize=7)

        ax.tick_params("both", labelsize=7)
        fig.tight_layout()
        if clonotype_type == "aa":
            fig.savefig(os.path.join(save_dir, donor + "_clonotype_aa_freq_hist_zoom_in.pdf"))
        elif clonotype_type == "nt":
            fig.savefig(os.path.join(save_dir, donor + "_clonotype_nt_freq_hist_zoom_in.pdf"))
        plt.close()


def generate_color_dict(key_list, min_luminance: float = 0.5):
    """
    Generate colors for a list of keys while avoiding colors too dark for black text.

    Parameters
    ----------
    key_list : list
        List of keys to assign colors.
    min_luminance : float, default 0.5
        Minimum relative luminance to keep a color. Lower values may be too dark for black text.
    """

    base_colors = plt.cm.tab20.colors

    def adjust_color(color, factor):
        r, g, b = color
        return (min(r * factor, 1), min(g * factor, 1), min(b * factor, 1))

    def luminance(color):
        r, g, b = color
        return 0.299 * r + 0.587 * g + 0.114 * b  # relative luminance formula

    factors = [1.0, 0.8, 1.2]  # normal, darker, lighter
    expanded_colors = []

    for factor in factors:
        for c in base_colors:
            new_c = adjust_color(c, factor)
            if luminance(new_c) >= min_luminance:
                expanded_colors.append(new_c)

    max_colors = len(expanded_colors)
    if len(key_list) > max_colors:
        raise ValueError(f"At most {max_colors} keys supported with min_luminance={min_luminance}")

    color_cycle = itertools.cycle(expanded_colors)
    color_dict = {s: mcolors.to_hex(next(color_cycle)) for s in key_list}

    return color_dict


def plot_custom_tcr_lib_pool_map(plate_df: pd.DataFrame, lib_col: str, lib_col_dict: dict, plate_name: str, well_col: str = "well", outs_dir: Path = None):
    plate_df = plate_df.copy()
    plate_df.dropna(subset=["cdr3j_alpha_nt_order_primers"], axis=0, inplace=True)

    for lib in plate_df[lib_col].unique():
        if lib not in lib_col_dict:
            raise Exception(f"All libraries in {lib_col}: {plate_df[lib_col].unique()} need a color in lib_col_dict: {list(lib_col_dict.keys())}")

    letter_to_number_map = dict(zip("ABCDEFGHIJKLMNOP", range(16)))

    plot_plate_df = pd.DataFrame(np.nan, index=letter_to_number_map.keys(), columns=np.arange(24))

    for col_idx, lib in enumerate(plate_df[lib_col].unique()):
        lib_df = plate_df.loc[plate_df[lib_col] == lib, :]
        for well in lib_df[well_col]:
            row = re.split(r"\d", well)[0]
            col = int(re.split(r"[A-P]", well)[1]) - 1
            i = letter_to_number_map[row]
            j = col

            if not np.isnan(plot_plate_df.iloc[i, j]):
                raise ValueError(f"Well {well} assigned twice")

            plot_plate_df.iloc[i, j] = col_idx + 1

    libs = list(plate_df[lib_col].unique())
    cmap = ListedColormap([lib_col_dict[lib] for lib in libs])
    bounds = np.arange(0.5, len(libs) + 1.5, 1)
    norm = BoundaryNorm(bounds, cmap.N)

    ax, fig, gs = startfig(w=36, h=55)
    ax.matshow(plot_plate_df, cmap=cmap, norm=norm)

    # Axis ticks
    ax.set_xticks(np.arange(plot_plate_df.shape[1]))
    ax.set_xticklabels([str(i + 1) for i in range(plot_plate_df.shape[1])], fontsize=8)
    ax.set_yticks(np.arange(plot_plate_df.shape[0]))
    ax.set_yticklabels(list(letter_to_number_map.keys()), fontsize=8)
    ax.xaxis.set_ticks_position("top")
    ax.yaxis.set_ticks_position("left")

    # Grid lines around wells
    ax.set_xticks(np.arange(-0.5, plot_plate_df.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, plot_plate_df.shape[0], 1), minor=True)
    ax.grid(which="minor", color="black", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Add well names inside cells
    for i, row_letter in enumerate(letter_to_number_map.keys()):
        for j in range(plot_plate_df.shape[1]):
            well_name = f"{row_letter}{j + 1}"
            ax.text(j, i, well_name, ha="center", va="center", fontsize=6, color="black")

    ax.set_title("Plate: " + str(plate_name), pad=25)

    legend_handles = [Patch(facecolor=lib_col_dict[lib], edgecolor="black", label=lib) for lib in libs]

    ax.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=1, frameon=False, fontsize=8)

    fig.tight_layout()

    if outs_dir:
        fig.savefig(os.path.join(outs_dir, f"plate_{plate_name}_library_pooling_map.pdf"))
        plt.close()
    else:
        plt.show()


# Modified from: https://stackoverflow.com/questions/36153410/how-to-create-a-swarm-plot-with-matplotlib
def simple_beeswarm2(y, nbins=None, width=1.0):
    """
    Returns x coordinates for points in `y` so plotting `x` vs `y`
    creates a deterministic, symmetric beeswarm plot.
    Includes zero values in the first bin.
    """
    y = np.asarray(y)

    if nbins is None:
        nbins = int(np.ceil(len(y) / 6))

    x = np.zeros(len(y))

    nn, ybins = np.histogram(y, bins=nbins)
    nmax = nn.max()

    ibs = []
    for k, (ymin, ymax) in enumerate(zip(ybins[:-1], ybins[1:])):
        if k == 0:
            # include lower bound for first bin
            i = np.nonzero((y >= ymin) & (y <= ymax))[0]
        else:
            i = np.nonzero((y > ymin) & (y <= ymax))[0]
        ibs.append(i)

    dx = width / max(1, nmax // 2)

    for i in ibs:
        if len(i) > 1:
            yy = y[i]
            j = len(i) % 2
            i = i[np.argsort(yy)]

            a = i[j::2]
            b = i[j + 1 :: 2]

            offset = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[a] = -offset
            x[b] = offset

    return x
