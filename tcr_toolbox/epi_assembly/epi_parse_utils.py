import pandas as pd
from Bio.Align import PairwiseAligner


def compute_tile_step_sizes(nt_tile_length: int, overlap_aa: int = 11):
    """Compute the nucleotide step size between consecutive minigene tiles for a given amino-acid overlap.

    Parameters
    ----------
    nt_tile_length : int
        Length, in nucleotides, of each tile; must be a multiple of 3 (one codon).
    overlap_aa : int, default = 11
        Number of amino acids by which consecutive tiles should overlap.

    Returns
    -------
    list[int]
        A 2-element list `[step_nt, step_nt]` (both equal) giving the nucleotide
        offset between the start of one tile and the next; consumed by
        `tile_ag_df`, which indexes it with `step_sizes[i % 2]`.

    Raises
    ------
    ValueError
        If `nt_tile_length` is not a multiple of 3, or if the requested overlap
        (`overlap_aa * 3` nt) is not strictly smaller than `nt_tile_length`.
    """
    if nt_tile_length % 3 != 0:
        raise ValueError("nt_tile_length must be a multiple of 3")

    overlap_nt = overlap_aa * 3

    if overlap_nt >= nt_tile_length:
        raise ValueError("overlap must be smaller than tile length")

    step_nt = nt_tile_length - overlap_nt

    return [step_nt, step_nt]


def tile_ag_df(
    ag_df: pd.DataFrame, nt_tile_length: int = 93, nt_seq_col: str = "minigene_nt", aa_seq_col: str = "minigene_aa", ag_name_col: str = "ag_name", tile_name_col: str = "ag_name"
):
    """Split each antigen's minigene coding sequence in `ag_df` into fixed-length, overlapping tiles.

    Rows whose nucleotide sequence (`nt_seq_col`) is at or under `nt_tile_length`
    are kept as a single tile, named after `ag_name_col`. Longer sequences are
    tiled into consecutive `nt_tile_length`-nt windows advanced by the step size
    from `compute_tile_step_sizes` (built from `nt_tile_length` and a fixed 11-aa
    target overlap, plus one extra tile anchored to the sequence's C-terminus so the
    end of the sequence is always fully covered (this last tile may overlap the
    preceding one by more than the target overlap). Each tile becomes its own row with
    the sliced nt/aa subsequences and a `{ag_name}_{tile_index}` tile name (lower
    dash, since dashes are used to separate TCR and Ag strings in TCR-Ag pair
    definition strings).

    Parameters
    ----------
    ag_df : pd.DataFrame
        One row per antigen, containing at least `nt_seq_col`, `aa_seq_col`, and
        `ag_name_col`.
    nt_tile_length : int, default = 93
        Tile length in nucleotides; must be a multiple of 3.
    nt_seq_col : str, default = "minigene_nt"
        Column holding each antigen's full nucleotide coding sequence.
    aa_seq_col : str, default = "minigene_aa"
        Column holding each antigen's full amino acid sequence (must be in frame
        with `nt_seq_col`).
    ag_name_col : str, default = "ag_name"
        Column holding each antigen's name, used as the prefix for tile names.
    tile_name_col : str, default = "ag_name"
        Column to write the resulting tile name into (`ag_name` for untiled
        rows, `{ag_name}_{tile_index}` for tiled rows).

    Returns
    -------
    pd.DataFrame
        `ag_df` with one row per tile (or per untiled antigen), reindexed from
        0; `nt_seq_col`/`aa_seq_col` are sliced to the tile, `tile_name_col` is
        set to the tile's name, and all other columns are copied from the
        source row.

    Raises
    ------
    Exception
        If `nt_tile_length` is not a multiple of 3.
    ValueError
        Propagated from `compute_tile_step_sizes` if the (fixed 11-aa) overlap
        is not smaller than `nt_tile_length`.
    """
    ag_df = ag_df.copy()
    if nt_tile_length % 3 != 0:
        raise Exception("CDS nt_tile_length needs to be a multiple of 3!")
    step_sizes = compute_tile_step_sizes(nt_tile_length=nt_tile_length)

    updated_rows = []

    for idx, row in ag_df.iterrows():
        nt_seq_len = len(row[nt_seq_col])
        ag_name = row[ag_name_col]

        if nt_seq_len % 3 != 0:
            raise ValueError(f"{ag_name}: nt sequence length ({nt_seq_len}) is not a multiple of 3 - cannot tile in frame.")

        if nt_seq_len <= nt_tile_length:
            row[tile_name_col] = ag_name
            updated_rows.append(row)
        else:
            starts = [0]
            i = 0
            while True:
                next_start = starts[-1] + step_sizes[i % 2]
                if next_start + nt_tile_length > nt_seq_len:
                    break
                starts.append(next_start)
                i += 1

            if starts[-1] + nt_tile_length < nt_seq_len:
                starts.append(nt_seq_len - nt_tile_length)

            for tile_index, start in enumerate(starts):
                nt_tile_seq = row[nt_seq_col][start : start + nt_tile_length]
                aa_tile_seq = row[aa_seq_col][start // 3 : start // 3 + nt_tile_length // 3]
                new_row = row.copy()
                new_row[nt_seq_col] = nt_tile_seq
                new_row[aa_seq_col] = aa_tile_seq
                new_row[tile_name_col] = f"{ag_name}_{tile_index}"  # lower dash because we reserve upper dash for Ag and TCR separation in names
                updated_rows.append(new_row)

    ag_df = pd.DataFrame(updated_rows).reset_index(drop=True)
    return ag_df


def var_start_end_for_kmer_overlap(minigene: str, reference: str, return_var_type: bool = False):
    """Find where a variant sits in `minigene`, and which overlap rule a kmer must satisfy to "contain" it.

    Aligns `minigene` (a variant amino acid sequence, e.g. a minigene
    translation) against `reference` (the wild-type sequence) to work out
    which minigene position(s) the variant affects, and classifies it into
    one of two kinds `run_mixmhcpred` treats differently when deciding whether
    a kmer overlaps the variant:
    - Missense, in-frame insertion, or frameshift: the variant occupies one or
      more minigene positions (the mutated/inserted residues), and a
      kmer overlaps it as soon as it contains at least one of them.
    - In-frame deletion: the deleted residues don't exist in `minigene` at
      all, so instead the two minigene residues immediately flanking the
      deletion — the neojunction it creates — are returned, and a kmer only
      overlaps if it contains both of them.
    A deletion that occurs together with another change nearby is classified
    with that other change instead of as a pure deletion, so it only needs a
    single overlapping residue rather than both flanking ones.

    Parameters
    ----------
    minigene : str
        The variant amino acid sequence (e.g. a minigene's translation) to
        align and locate positions within.
    reference : str
        The wild-type amino acid sequence to align `minigene` against.
    return_var_type : bool, default = False
        If True, also return the variant type classification.

    Returns
    -------
    tuple
        `(start, end)` normally, or `(start, end, var_type)` when
        `return_var_type` is True, where `var_type` is `"inframe_del"` or
        `"missense_or_inframe_ins_or_frameshift"`. If the aligner produces no
        alignment at all, `start`/`end` (and `var_type`, when requested) are
        all `None` — always matching the requested tuple width.

    Raises
    ------
    Exception
        Not raised by this function itself; downstream callers (e.g.
        `run_mixmhcpred`) raise if `var_type` is anything other than the two
        values above.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    alignments = aligner.align(reference, minigene)
    if not alignments:
        return (None, None, None) if return_var_type else (None, None)

    aln = alignments[0]

    fasta = aln.format("fasta").splitlines()
    ref_gapped = fasta[1]
    mini_gapped = fasta[3]

    mini_pos = -1
    start = None
    end = None
    var_type = None
    seen_non_del_change = False

    for r, m in zip(ref_gapped, mini_gapped):
        # print(r, m)
        if m != "-":
            mini_pos += 1
            # print(mini_pos)

        mini_pos = mini_pos if mini_pos >= 0 else 0

        if (r == "-" and m != "-") or (r != m and r != "-" and m != "-"):
            seen_non_del_change = True

        # deletion
        if m == "-" and r != "-":
            # print("m:", m, "r:", r)
            if start is None:
                start = mini_pos
            end = mini_pos + 1
            if not seen_non_del_change:
                var_type = "inframe_del"
            else:
                var_type = "missense_or_inframe_ins_or_frameshift"

        # rest
        elif (r == "-" and m != "-") or r != m:
            # print("m:", m, "r:", r)
            if start is None:
                start = mini_pos
            end = mini_pos + 1
            var_type = "missense_or_inframe_ins_or_frameshift"  # as an inframe_ins can potentially also run until the end like a frameshift and we want to select kmers overlapping the frameshift in the same way as insertions, we return inframe_ins for frameshift

    if return_var_type:
        return start, end, var_type
    else:
        return start, end


def align_to_add_var_start_and_end_aa_seq_ag_df(
    ag_df: pd.DataFrame,
    aa_var_seq_col: str,
    aa_ref_seq_col: str,
    var_start_col: str = "minigene_aa_var_start",
    var_end_col: str = "minigene_aa_var_end",
    var_type_col: str = "var_type_for_kmer_overlap_col",
):
    """Annotate each row of `ag_df` with its variant's amino acid start/end position and type, via alignment.

    For every row, aligns `aa_var_seq_col` (the variant/minigene amino acid
    sequence) against `aa_ref_seq_col` (the wild-type reference) using
    `var_start_end_for_kmer_overlap(..., return_var_type=True)`, and writes the
    resulting `(start, end, var_type)` into `var_start_col`, `var_end_col`, and
    a fixed column named `"var_type_for_kmer_overlap_col"`.

    Parameters
    ----------
    ag_df : pd.DataFrame
        One row per antigen/minigene, containing at least `aa_var_seq_col` and
        `aa_ref_seq_col`.
    aa_var_seq_col : str
        Column holding each row's variant (minigene) amino acid sequence.
    aa_ref_seq_col : str
        Column holding each row's wild-type reference amino acid sequence.
    var_start_col : str, default = "minigene_aa_var_start"
        Column to write each row's variant start position into.
    var_end_col : str, default = "minigene_aa_var_end"
        Column to write each row's variant end position into.

    Returns
    -------
    pd.DataFrame
        `ag_df` with `var_start_col`, `var_end_col`, and
        `"var_type_for_kmer_overlap_col"` added (or overwritten).
    """
    ag_df = ag_df.copy()
    ag_df[[var_start_col, var_end_col, "var_type_for_kmer_overlap_col"]] = ag_df.apply(
        lambda r: var_start_end_for_kmer_overlap(r[aa_var_seq_col], r[aa_ref_seq_col], return_var_type=True), axis=1, result_type="expand"
    )
    return ag_df
