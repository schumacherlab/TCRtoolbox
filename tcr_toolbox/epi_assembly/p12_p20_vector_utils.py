import collections
import json
import os
from itertools import product
from re import search
from typing import Union
import dotenv

import numpy as np
import pysam
from Bio import Restriction, SeqIO
from Bio.Seq import Seq
from numpy.random import default_rng
from pydna.amplify import Anneal
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer

from tcr_toolbox.tcr_assembly.constants import sites_to_check
from tcr_toolbox.utils.codon_optimizer import make_search_strings
from tcr_toolbox.utils.constants import AA_TO_ID
from tcr_toolbox.epi_assembly.constants import p12_cloning_sites, p20_cloning_sites, p12_p20_flanking_nt_dict
from tcr_toolbox.utils.utils import reverse_translate

rng = default_rng(123)
dotenv.load_dotenv()
tcr_toolbox_data_path = os.getenv("tcr_toolbox_data_path")


def generate_barcode_dict(number_barcodes: int, barcode_large_bool: bool):
    """Sample a random, non-overlapping set of barcodes free of BsmBI/other flagged restriction sites.

    Loads the 18-mer, divergence-based error-correcting DNA barcodes (FreeBarcodes) fasta from:
    https://github.com/hawkjo/freebarcodes/blob/master/README.md
    If `barcode_large_bool` is True, a double 18-mer (i.e., 36-mer barcode) fasta is loaded and
    a single otherwise. The function discards any barcode that contains one of the sites in
    `sites_to_check` or a full BsmBI recognition site anywhere in the sequence, or a partial
    (5-of-6 bp) BsmBI site fragment within 5 bp of either edge (which could complete a false
    BsmBI site once flanking sequence is added), and then randomly samples `number_barcodes` of
    the remaining barcodes (seeded via the module-level `rng`).

    Parameters
    ----------
    number_barcodes : int
        Number of barcodes to sample for the library.
    barcode_large_bool : bool
        If True, use the double 18-mer (18-2_18-2) barcode fasta instead of the
        single 18-mer (18-2) barcode fasta.

    Returns
    -------
    collections.defaultdict[str, str]
        Mapping of sampled barcode name to barcode nucleotide sequence.

    Raises
    ------
    Exception
        If fewer barcodes remain after site filtering than `number_barcodes` were
        requested.
    """
    if barcode_large_bool:
        barcodes_fasta_fname = f"{tcr_toolbox_data_path}/tcr_toolbox_datasets/tcr_assembly/minigene_freebarcodes/barcodes_18-2_18-2.fasta"
    else:
        barcodes_fasta_fname = f"{tcr_toolbox_data_path}/tcr_toolbox_datasets/tcr_assembly/minigene_freebarcodes/barcodes18-2.fasta"

    number_barcodes = number_barcodes
    barcodes_dict = collections.defaultdict(str)
    (search_string, restriction_enzyme_sites, restriction_edge_sites) = make_search_strings(enzymes=["BsmBI"], sites_to_check=sites_to_check)

    with pysam.FastxFile(barcodes_fasta_fname) as fasta_in:
        for barcode_entry in fasta_in:
            restriction_edge_bool = (search(restriction_edge_sites, barcode_entry.sequence[:5]) is not None) | (
                search(restriction_edge_sites, barcode_entry.sequence[-5:]) is not None
            )
            if (search(search_string, barcode_entry.sequence) is None) & (restriction_edge_bool == False):
                barcodes_dict[barcode_entry.name] = barcode_entry.sequence
            else:
                continue
    print("Length barcode list after site removal:", len(barcodes_dict))

    barcodes_dict_final = collections.defaultdict(str)
    sample_barcodes_names_list = list(barcodes_dict.keys())
    if len(sample_barcodes_names_list) > number_barcodes:
        selected_indices = rng.choice(len(sample_barcodes_names_list), size=number_barcodes, replace=False)
        sample_barcodes_names_list = [sample_barcodes_names_list[i] for i in selected_indices].copy()

        for barcode_name in sample_barcodes_names_list:
            barcodes_dict_final[barcode_name] = barcodes_dict[barcode_name]

    elif len(sample_barcodes_names_list) < number_barcodes:
        raise Exception("Not enough barcodes for library size! Set barcode_long to True!")
    else:
        barcodes_dict_final = barcodes_dict.copy()

    return barcodes_dict_final


def _print_progress_update(count: int, total: int, all_label: str, label: str, step: int = 2500) -> None:
    if count == total:
        print("All", count, all_label)
    elif total > step and count % step == 0:
        print(count, label)


def write_barcoded_p12_p20_epitope_oligo_order_fasta(
    epitopes_nucl_list: list,
    epitopes_names_list: list,
    primers_fw_list: list,
    primers_rv_list: list,
    barcode_large_bool: bool,
    fasta_out_fname: Union[str, os.PathLike[str]],
    p12_or_p20: str = "p20",
    subpool_indices_dict: dict = None,
    no_uvp_but_only_subpool_primers: bool = False,
    max_oligo_length: int = 300,
):
    """Assemble and write a Ag-encoding minigene or minimal peptide-encoding oligonucleotide
    pool fasta file.

    To assemble an Ag- or minimal-peptide-encoding oligonucleotide, a TAA stop-codon is
    added directly after each codon-optimized epitope coding sequence and followed by
    a unique 18 bp or 36 bp (required for minimal peptide libraries) DNA barcode generated
    via `generate_barcode_dict`. The resulting nucleotide sequence is flanked by BsmBI
    recognition sites (the p12/p20 cloning sites) and, by default, the universal primer
    binding sites (`primers_fw_list[0]`/`primers_rv_list[0]`) needed for PCR amplification
    to generate the dsDNA needed for subsequent Golden Gate assembly into the p12/p20
    lentiviral vector.

    To save costs, epitope libraries that are screened simultaneously (e.g. from different
    patients) can each be added as a subpool to a single oligonucleotide library via
    `subpool_indices_dict`. Each subpool is flanked by its own orthogonal primer pair.
    By default, this orthogonal primer pair is nested inside the universal primer
    binding sites, so that the full pool can also be PCR amplified as a whole using
    only the universal primers. If `no_uvp_but_only_subpool_primers` is set instead,
    the universal primer binding sites are omitted entirely and each subpool is
    flanked only by its own orthogonal primer pair. Either way, each subpool can later
    be PCR amplified out of the full pool individually, to generate the dsDNA that is
    needed for Golden Gate assembly into the p12/p20 lentiviral vector. The resulting
    oligos are written to `fasta_out_fname`, and, when subpools are used,
    `subpool_indices_dict` is also written to a corresponding `.json` file (same
    basename) for later use by `check_p12_p20_oligo_order_fasta_pydna`.

    Parameters
    ----------
    epitopes_nucl_list : list
        Codon-optimized nucleotide coding sequences (length must be a multiple of 3)
        for each epitope, in order.
    epitopes_names_list : list
        Name/identifier for each epitope, matching `epitopes_nucl_list` in length and
        order; used in the fasta header.
    primers_fw_list : list
        Forward primer sequences. If `subpool_indices_dict` is None, only
        `primers_fw_list[0]` (the universal forward primer) is used. If subpools are
        used and `no_uvp_but_only_subpool_primers` is True, `primers_fw_list[i]` is
        the orthogonal forward primer for subpool i (no universal primer is
        included); if False, `primers_fw_list[0]` is the universal forward primer
        and `primers_fw_list[1:]` are the per-subpool orthogonal forward primers
        (index i+1 for subpool i).
    primers_rv_list : list
        Reverse primer sequences, structured analogously to `primers_fw_list`.
    barcode_large_bool : bool
        Passed through to `generate_barcode_dict`; if True, use double-length
        (36 nt) barcodes instead of single-length (18 nt) ones.
    fasta_out_fname : str | os.PathLike[str]
        Output path for the oligo-order fasta.
    p12_or_p20 : str, default = "p20"
        Which vector's BsmBI cloning sites to flank the epitope with; must be "p12"
        or "p20".
    subpool_indices_dict : dict, default = None
        Mapping of zero-indexed, consecutive subpool number to the list of
        `epitopes_*_list` indices belonging to that subpool. Must cover all epitope
        indices exactly once and contain at least 2 subpools. If None, no
        subpool/orthogonal primers or tagging are added and a single flat pool is
        written using only the universal primers.
    no_uvp_but_only_subpool_primers : bool, default = False
        If True, skip the universal primers (UVP) and flank each epitope only with
        its subpool's orthogonal primers. Only relevant when `subpool_indices_dict`
        is given.
    max_oligo_length : int, default = 300
        Maximum allowed length, in bp, for a resulting order oligo; any oligo
        longer than this raises an exception. 300 bp is the current limit for
        most oligonucleotide manufacturers (e.g. Twist); override this if your
        manufacturer supports longer syntheses.

    Returns
    -------
    None
        Writes `fasta_out_fname` (and, for subpools, a matching `.json` file) as a
        side effect and prints progress/completion messages.

    Raises
    ------
    ValueError
        If `p12_or_p20` is not "p12" or "p20".
    Exception
        If `epitopes_names_list` and `epitopes_nucl_list` differ in length; if the
        number of generated barcodes does not match the number of epitopes; if
        `subpool_indices_dict` has fewer than 2 subpools, is not zero-indexed with
        increments of 1, has a primer-count mismatch with the number of subpools,
        does not cover all epitope indices exactly once, or has overlapping indices
        between subpools; if any epitope nucleotide sequence length is not a
        multiple of 3; if a sampled barcode is reused across epitopes; if the stored
        barcode is not the reverse complement of the barcode encoded in the
        resulting oligo; or if any resulting order oligo exceeds `max_oligo_length`.
    """
    if p12_or_p20 not in {"p12", "p20"}:
        raise ValueError(f"p12_or_p20 must be 'p12' or 'p20' and not {p12_or_p20}")

    if p12_or_p20 == "p20":
        end_5_BsmBI = p20_cloning_sites["end_5_BsmBI_p20"]
        end_3_BsmBI = p20_cloning_sites["end_3_BsmBI_p20"]
    elif p12_or_p20 == "p12":
        end_5_BsmBI = p12_cloning_sites["end_5_BsmBI_p12"]
        end_3_BsmBI = p12_cloning_sites["end_3_BsmBI_p12"]

    epitopes_order_dict = collections.defaultdict(str)

    if len(epitopes_names_list) != len(epitopes_nucl_list):
        raise Exception("Length epitope names list is not equal to length epitope nucleotide list!")

    print("Generating barcodes!")
    barcodes_dict_final = generate_barcode_dict(len(epitopes_nucl_list), barcode_large_bool=barcode_large_bool)

    if len(epitopes_names_list) != len(barcodes_dict_final.keys()):
        raise Exception("Length epitopes list is not equal to number of generated barcodes!")

    if subpool_indices_dict:
        print("Adding subpool primers...")

        if len(subpool_indices_dict.keys()) < 2:
            raise Exception("There need to be at least 2 subpools!")

        if (list(subpool_indices_dict.keys())[0] != 0) | (list(subpool_indices_dict.keys())[1] != 1):
            raise Exception("Subpool indices must be zero-indexed list of subpool numbers with increments of 1 (e.g., [0, 1, 2])!")

        if no_uvp_but_only_subpool_primers:
            if (len(subpool_indices_dict.keys()) != len(primers_fw_list)) | (len(subpool_indices_dict.keys()) != len(primers_rv_list)):
                raise Exception("Number of subpools is not equal to the number of additionally added primers!")
        else:
            if (len(subpool_indices_dict.keys()) != len(primers_fw_list) - 1) | (len(subpool_indices_dict.keys()) != len(primers_rv_list) - 1):
                raise Exception("Number of subpools is not equal to the number of additionally added primers!")

        for subpool in subpool_indices_dict.keys():
            # Remove potential np.int64 encoding from the dictionary
            subpool_indices_dict[subpool] = [idx.item() if isinstance(idx, np.int64) else idx for idx in subpool_indices_dict[subpool]].copy()

        if set([idx for subpool_indices in subpool_indices_dict.values() for idx in subpool_indices]) != set([i for i in np.arange(len(epitopes_names_list))]):
            raise Exception("Subpool indices do not cover all indices in the epitopes lists!")

        for subpool_1, subpool_2 in product(subpool_indices_dict.keys(), subpool_indices_dict.keys()):
            if subpool_1 != subpool_2:
                if len(set(subpool_indices_dict[subpool_1]).intersection(set(subpool_indices_dict[subpool_2]))) > 0:
                    raise Exception("Subpool indices overlap between subpools:", subpool_1, subpool_2)

        barcodes_used_list = []
        z = 0
        total_epitopes = len(epitopes_names_list)
        print("Start processing epitopes for fasta writing!")
        for idx, (barcode_name, epitope_name, epitope) in enumerate(zip(barcodes_dict_final.keys(), epitopes_names_list, epitopes_nucl_list)):
            if len(epitope) % 3 != 0:
                raise Exception("Epitope nucleotide input sequence must be a coding-sequence and have a multiple of 3!", "Did you codon optimize it?")

            # This should really be impossible to occur...but sensible check to have:
            if barcodes_dict_final[barcode_name] in barcodes_used_list:
                raise Exception("Barcodes overlap between subpools!")

            subpool_current_idx_tmp_list = [subpool for subpool in subpool_indices_dict.keys() if idx in subpool_indices_dict[subpool]]

            if len(subpool_current_idx_tmp_list) != 1:
                if len(subpool_current_idx_tmp_list) > 1:
                    raise Exception("Current epitope index occurs in multiple subpools!")
                else:
                    raise Exception("Current epitope index occurs in no subpools!")

            if no_uvp_but_only_subpool_primers:
                epitopes_order_dict[barcode_name + "_" + str(epitope_name) + "_" + str(subpool_current_idx_tmp_list[0])] = (
                    primers_fw_list[subpool_current_idx_tmp_list[0]]
                    + end_5_BsmBI
                    + epitope
                    + reverse_translate("*")
                    + barcodes_dict_final[barcode_name]
                    + end_3_BsmBI
                    + str(Seq(primers_rv_list[subpool_current_idx_tmp_list[0]]).reverse_complement())
                )
            else:
                epitopes_order_dict[barcode_name + "_" + str(epitope_name) + "_" + str(subpool_current_idx_tmp_list[0])] = (
                    primers_fw_list[0]
                    + primers_fw_list[subpool_current_idx_tmp_list[0] + 1]
                    + end_5_BsmBI
                    + epitope
                    + reverse_translate("*")
                    + barcodes_dict_final[barcode_name]
                    + end_3_BsmBI
                    + str(Seq(primers_rv_list[subpool_current_idx_tmp_list[0] + 1]).reverse_complement())
                    + str(Seq(primers_rv_list[0]).reverse_complement())
                )

            barcodes_used_list.append(barcodes_dict_final[barcode_name])
            z += 1
            _print_progress_update(count=z, total=total_epitopes, all_label="epitopes processed!", label="epitopes processed!")

        # Check whether barcode names are still the reverse complement of the stored barcodes:
        first_key = list(epitopes_order_dict.keys())[0]
        barcode = first_key.split("_")[0]
        epitope_order_oligo_first_key = epitopes_order_dict[first_key]

        barcode_length = 2 * 18 if barcode_large_bool else 1 * 18
        end_3_len = len(end_3_BsmBI)
        if no_uvp_but_only_subpool_primers:
            start_idx = -(end_3_len + 20 + barcode_length)
            end_idx = -(end_3_len + 20)
        else:
            start_idx = -(end_3_len + 40 + barcode_length)
            end_idx = -(end_3_len + 40)

        if Seq(barcode).reverse_complement() != epitope_order_oligo_first_key[start_idx:end_idx]:
            raise Exception("Barcode names are not the reverse complement of the stored barcodes!")

        # Check whether order oligos are longer than max_oligo_length:
        if any([len(epitope_order_oligo) > max_oligo_length for epitope_order_oligo in epitopes_order_dict.values()]):
            raise Exception(
                f"There are order oligos that are longer than {max_oligo_length} bp! This is not allowed by your oligonucleotide manufacturer! Please remove these and re-run!"
            )

        with open(fasta_out_fname, "w") as fasta_out:
            for name, epitope_order_oligo in epitopes_order_dict.items():
                fasta_out.write(">" + str(name) + "\n")
                fasta_out.write(epitope_order_oligo + "\n")

        with open(os.path.join(os.path.dirname(fasta_out_fname), os.path.basename(fasta_out_fname).split(".")[0] + ".json"), "w") as json_out:
            json.dump(subpool_indices_dict, json_out)

    else:
        print("Start processing epitopes for fasta writing!")
        z = 0
        total_epitopes = len(epitopes_names_list)
        for barcode_name, epitope_name, epitope in zip(barcodes_dict_final.keys(), epitopes_names_list, epitopes_nucl_list):
            if len(epitope) % 3 != 0:
                raise Exception("Epitope nucleotide input sequence must be a coding-sequence and have a multiple of 3!", "Did you codon optimize it?")

            epitopes_order_dict[barcode_name + "_" + str(epitope_name)] = (
                primers_fw_list[0]
                + end_5_BsmBI
                + epitope
                + reverse_translate("*")
                + barcodes_dict_final[barcode_name]
                + end_3_BsmBI
                + str(Seq(primers_rv_list[0]).reverse_complement())
            )

            z += 1
            _print_progress_update(count=z, total=total_epitopes, all_label="epitopes processed!", label="epitopes processed!")

        # Check whether barcode names are still the reverse complement of the stored barcodes:
        first_key = list(epitopes_order_dict.keys())[0]
        barcode = first_key.split("_")[0]
        epitope_order_oligo_first_key = epitopes_order_dict[first_key]

        barcode_length = 2 * 18 if barcode_large_bool else 1 * 18
        end_3_len = len(end_3_BsmBI)
        start_idx = -(end_3_len + 20 + barcode_length)
        end_idx = -(end_3_len + 20)

        if Seq(barcode).reverse_complement() != epitope_order_oligo_first_key[start_idx:end_idx]:
            raise Exception("Barcode names are not the reverse complement of the stored barcodes!")

        # Check whether order oligos are longer than max_oligo_length:
        if any([len(epitope_order_oligo) > max_oligo_length for epitope_order_oligo in epitopes_order_dict.values()]):
            raise Exception(
                f"There are order oligos that are longer than {max_oligo_length} bp! This is not allowed by your oligonucleotide manufacturer! Please remove these and re-run!"
            )

        with open(fasta_out_fname, "w") as fasta_out:
            for name, epitope_order_oligo in epitopes_order_dict.items():
                fasta_out.write(">" + str(name) + "\n")
                fasta_out.write(epitope_order_oligo + "\n")

    return print("Fasta file written!")


def check_p12_p20_oligo_order_fasta_pydna(
    fasta_in_fname: Union[str, os.PathLike[str]],
    primers_fw_list: list,
    primers_rv_list: list,
    aa_check_list: list,
    p12_or_p20: str = "p20",
    source_peptide: str = None,
    fusion_protein_bool_list: list = None,
    also_check_full_pool: bool = False,
    # if epitope_variant_hamming_pool is false,
    # fasta entry.name.split('_')[1] must indicate instead of hamming circle
    # what the fasta entry epitope belongs to (e.g., to what biological sample).
    # (move this descritpion to function documentation when we start writing
    # documentation for all functions):
    epitope_variant_hamming_pool: bool = True,
    no_uvp_but_only_subpool_primers: bool = False,
):
    """Validate a p12/p20 oligo-order fasta by modeling cloning into the vector in silico with pydna.

    To identify design errors before wet lab cloning, cloning of each Ag or minimal
    peptide oligo (as written by `write_barcoded_p12_p20_epitope_oligo_order_fasta`)
    into the p12/p20 lentiviral vector is modeled in silico using pydna, modeling PCR,
    restriction digestion, and DNA ligation. For each entry, the resulting
    insert's vector position, translation, barcode, stop codon and (when
    `epitope_variant_hamming_pool` is True) Hamming distance to `source_peptide` are
    then checked against the fasta header and `aa_check_list` (see Raises for what
    triggers a failure). The Hamming distance is computed over residue 1 and
    residues 3-8 of the peptide only (0-indexed: positions `[0]` and `[2:8]`) —
    the two MHC anchor residues are excluded, and for peptides longer than 9
    residues, positions beyond index 8 are not compared. Proteins entries
    encoded by RNA fusion transcripts are flagged in
    `fusion_protein_bool_list` and skipped. If a subpool `.json` file is present, subpool
    PCRs are checked first (each amplified with its own orthogonal primers); the full pool is
    additionally checked (using only the universal primers) when no subpool file is found or
    `also_check_full_pool` is True. Because the full-pool check always amplifies with
    `primers_fw_list[0]`/`primers_rv_list[0]` under the assumption that these are universal
    primers, `also_check_full_pool` cannot be combined with `no_uvp_but_only_subpool_primers=True`
    (where index 0 is subpool 0's orthogonal primer, not a universal one) — this combination
    raises a `ValueError`.

    Parameters
    ----------
    fasta_in_fname : str | os.PathLike[str]
        Path to the oligo-order fasta to validate (as written by
        `write_barcoded_p12_p20_epitope_oligo_order_fasta`). A corresponding
        `<basename>.json` file, if present, is read as the `subpool_indices_dict`
        used to select subpool-specific primers.
    primers_fw_list : list
        Forward primers, structured as in
        `write_barcoded_p12_p20_epitope_oligo_order_fasta`: if
        `no_uvp_but_only_subpool_primers` is True, `primers_fw_list[i]` is the
        orthogonal forward primer for subpool i; if False, `primers_fw_list[0]` is
        the universal forward primer and `primers_fw_list[i + 1]` is the orthogonal
        forward primer for subpool i.
    primers_rv_list : list
        Reverse primers, structured analogously to `primers_fw_list`.
    aa_check_list : list
        Expected amino acid translation for each cloned fasta entry, in the same order
        as the fasta; also used to size `fusion_protein_bool_list` when the latter is
        not given.
    p12_or_p20 : str, default = "p20"
        Which vector's BsmBI cloning sites/flanking sequences to use for locating
        and excising the insert; must be "p12" or "p20".
    source_peptide : str, default = None
        Wild-type/reference amino acid sequence used as the Hamming-distance
        baseline for non-zero, non-"n" Hamming circles and for negative controls,
        when `epitope_variant_hamming_pool` is True.
    fusion_protein_bool_list : list, default = None
        Per-entry booleans (same order as the fasta) marking entries whose
        translation should be skipped because they encode a protein encoded by a
        RNA fusion transcript. If not given, defaults to all-False.
    also_check_full_pool : bool, default = False
        If True, also run the full-pool check (with universal primers only) even
        when a subpool `.json` file was found. Cannot be True at the same time as
        `no_uvp_but_only_subpool_primers=True` (raises `ValueError`); see above.
    epitope_variant_hamming_pool : bool, default = True
        If True, treat `entry.name.split("_")[1]` as a Hamming-circle number
        ("0"/int/"n") and validate it against the Hamming distance to
        `source_peptide`, and validate `entry.name.split("_")[2]` against
        `aa_check_list`.
    no_uvp_but_only_subpool_primers : bool, default = False
        If True, use only the subpool-specific primers (no universal primers) when
        checking subpools. Only relevant when a subpool `.json` file is found.

    Returns
    -------
    None
        Prints progress and a final "Finished checking all oligos!" message;
        raises on the first detected inconsistency instead of returning a report.

    Raises
    ------
    ValueError
        If `p12_or_p20` is not "p12" or "p20"; if `also_check_full_pool` and
        `no_uvp_but_only_subpool_primers` are both True; or if `aa_check_list` is
        not provided.
    Exception
        If a subpool index appears in zero or multiple subpools, or the subpool
        implied by the entry's index disagrees with the subpool encoded in its
        fasta header; if PCR pre-amplification of an oligo yields other than
        exactly one product, or a product of unexpected length; if BsmBI digestion
        of the amplicon does not yield exactly 3 fragments; if Golden Gate ligation
        into the vector does not yield a circular product; if the insert is not
        found at the expected vector position; if the insert's translation,
        barcode, or stop codon do not match expectations; if a barcode is reused
        across entries; or if the Hamming-circle number (or negative-control
        distance) encoded in the fasta header does not match the actual Hamming
        distance to `source_peptide`.
    """
    if p12_or_p20 not in {"p12", "p20"}:
        raise ValueError(f"p12_or_p20 must be 'p12' or 'p20' and not {p12_or_p20}")

    if not aa_check_list:
        raise ValueError("aa_check_list must be provided and non-empty.")

    if also_check_full_pool and no_uvp_but_only_subpool_primers:
        raise ValueError(
            "also_check_full_pool cannot be combined with no_uvp_but_only_subpool_primers=True: "
            "without universal primers, primers_fw_list[0]/primers_rv_list[0] are subpool 0's "
            "orthogonal primers, not universal primers, so the full-pool check would silently "
            "amplify only subpool 0's entries instead of the full pool."
        )

    p_dna_file = (
        tcr_toolbox_data_path + "/tcr_toolbox_datasets/epi_assembly/plasmids/P20_pTwist_Lenti_SFFV_coNeoUTR2-3_IL2-SP_RT-TCRb1_correct_orient.dna"
        if p12_or_p20 == "p20"
        else tcr_toolbox_data_path + "/tcr_toolbox_datasets/epi_assembly/plasmids/P12_pTwist_Lenti_SFFV_Pfiz_CD74_RT-TCRb1_correct_orient.dna"
    )

    subpool_indices_dict = None

    records = list(SeqIO.parse(p_dna_file, "snapgene"))
    p_record = records[0]
    p = Dseq(str(p_record.seq), circular=True)

    json_in_fname = os.path.dirname(fasta_in_fname) + "/" + os.path.basename(fasta_in_fname).split(".")[0] + ".json"
    if not os.path.isfile(json_in_fname):
        print(
            "No subpool_indices_dict json file detected!",
            "Is it correct that you encoded only a single full pool in this oligo pool?",
            "OR did you accidentally remove the json file that stores the subpool_indices_dict?",
        )
    else:
        with open(json_in_fname, "r") as json_in:
            subpool_indices_dict = json.load(json_in)
        if not subpool_indices_dict:
            raise Exception(f"subpool_indices_dict json file '{json_in_fname}' exists but is empty!")

    if not fusion_protein_bool_list:
        fusion_protein_bool_list = [False] * len(aa_check_list)

    if subpool_indices_dict:
        print("Started checking subpools:", subpool_indices_dict.keys())
        pcr_products_dict = collections.defaultdict(list)
        ligation_products_dict = collections.defaultdict(list)
        barcodes_used_list = []
        total_oligos = len(aa_check_list)
        with pysam.FastxFile(fasta_in_fname) as fasta_in:
            for i, (entry, epitope_variant, fusion_protein_bool) in enumerate(zip(fasta_in, aa_check_list, fusion_protein_bool_list)):
                if fusion_protein_bool:
                    print("Skipping fusion protein check! Check fusion protein fasta idx manually:", i)
                    continue
                _print_progress_update(count=i + 1, total=total_oligos, all_label="oligos are correct!", label="oligos correct!")

                if epitope_variant_hamming_pool:
                    if not epitope_variant == entry.name.split("_")[2]:
                        raise Exception(
                            "Epitope variant translation in fasta header does not equal epitope variant sequence in aa check list!",
                            "Epitope variant in fasta header:",
                            entry.name.split("_")[2],
                            "Epitope variant in aa check list:",
                            epitope_variant,
                        )

                subpool_current_idx_tmp_list = [subpool for subpool in subpool_indices_dict.keys() if i in subpool_indices_dict[subpool]]

                if not len(subpool_current_idx_tmp_list) == 1:
                    if len(subpool_current_idx_tmp_list) > 1:
                        raise Exception("Current epitopes list idx occurs in multiple subpools!")

                    if len(subpool_current_idx_tmp_list) == 0:
                        raise Exception("Current epitopes list idx does not occur in any subpool!")

                if not subpool_current_idx_tmp_list[0] == entry.name.split("_")[-1]:
                    raise Exception(
                        "Idx matches subpool that does not equal subpool in fasta header:",
                        "Idx subpool match:",
                        subpool_current_idx_tmp_list[0],
                        "Subpool in fasta header:",
                        entry.name.split("_")[4],
                    )

                if not no_uvp_but_only_subpool_primers:
                    primer_fw = primers_fw_list[int(subpool_current_idx_tmp_list[0]) + 1]
                    primer_rv = primers_rv_list[int(subpool_current_idx_tmp_list[0]) + 1]

                else:
                    primer_fw = primers_fw_list[int(subpool_current_idx_tmp_list[0])]
                    primer_rv = primers_rv_list[int(subpool_current_idx_tmp_list[0])]

                pcr_products_dict[entry.name] = Anneal((Primer(primer_fw), Primer(primer_rv)), Dseqrecord(entry.sequence)).products

                if not len(pcr_products_dict[entry.name]) == 1:
                    raise Exception("Oligo pre-amplification either has too many products or no products:", pcr_products_dict[entry.name])

                product_len = len(entry.sequence) - len(primer_fw) - len(primer_rv) if not no_uvp_but_only_subpool_primers else len(entry.sequence)

                if not len(pcr_products_dict[entry.name][0].seq) == product_len:
                    raise Exception("Oligo pre-amplification amplicon has wrong size!", pcr_products_dict[entry.name])

                if not len(pcr_products_dict[entry.name][0].seq.cut(Restriction.BsmBI)) == 3:
                    raise Exception("BsmBI digested pre-amplified oligo has more than 3 products:", pcr_products_dict[entry.name][0].seq.cut(Restriction.BsmBI))

                ligation_products_dict[entry.name] = (p.cut(Restriction.BsmBI)[0] + pcr_products_dict[entry.name][0].seq.cut(Restriction.BsmBI)[1]).looped()

                if not ligation_products_dict[entry.name].circular:
                    raise Exception("Golden gate assembly into vector failed!")

                if p12_or_p20 == "p20":
                    # five_end starts at the end of the IL2-Signal Peptide (i.e., start of the epitope):
                    five_end = ligation_products_dict[entry.name].find(p12_p20_flanking_nt_dict["p20_end_5"]) + len(p12_p20_flanking_nt_dict["p20_end_5"])
                    # three_end is at the end of the barcode. So, three end is at the end of epitope + stop codon + barcode:
                    three_end = ligation_products_dict[entry.name].find(p12_p20_flanking_nt_dict["p20_end_3"])

                if p12_or_p20 == "p12":
                    # five_end starts at the end of the CD74-Signal Peptide (i.e., start of the epitope):
                    five_end = ligation_products_dict[entry.name].find(p12_p20_flanking_nt_dict["p12_end_5"]) + len(p12_p20_flanking_nt_dict["p12_end_5"])
                    # three_end is at the end of the barcode. So, three end is at the end of epitope + stop codon + barcode:
                    three_end = ligation_products_dict[entry.name].find(p12_p20_flanking_nt_dict["p12_end_3"])

                # Five_end should be 9059 (P20) or 9079 (P12) and three_end 0:
                if p12_or_p20 == "p20" and not ((five_end == 9059) & (three_end == 0)):
                    raise Exception("Insert inserted at the wrong vector position!")
                if p12_or_p20 == "p12" and not ((five_end == 9079) & (three_end == 0)):
                    raise Exception("Insert inserted at the wrong vector position!")

                insert_region = ligation_products_dict[entry.name][five_end:three_end]
                insert_translation = insert_region.translate(to_stop=True)
                # Extract barcode and stop codon from insert using string slicing.
                # Avoids pydna circular DNA position arithmetic which breaks when three_end=0.
                insert_seq_str = str(insert_region)
                barcode_rc = entry.name.split("_")[0]
                barcode_len = len(barcode_rc)
                barcode_in_insert = insert_seq_str[-barcode_len:]
                stop_in_insert = insert_seq_str[-(barcode_len + 3) : -barcode_len]
                barcode_rc_from_insert = str(Seq(barcode_in_insert).reverse_complement())

                if not (epitope_variant == insert_translation) & (barcode_rc == barcode_rc_from_insert) & ((stop_in_insert == "TAA") | (stop_in_insert == "TAG")):
                    if epitope_variant != insert_translation:
                        raise Exception(
                            "idx:",
                            str(i),
                            entry.name,
                            "translation in ligation_product does not match original AA sequence!",
                            "Original AA sequence:",
                            epitope_variant,
                            "Translation in ligation product:",
                            insert_translation,
                        )
                    elif barcode_rc != barcode_rc_from_insert:
                        raise Exception(
                            "Barcode in ligation_product does not match original Barcode sequence!",
                            "Original Barcode sequence:",
                            barcode_rc,
                            "Barcode in ligation product (RC):",
                            barcode_rc_from_insert,
                        )
                    elif stop_in_insert not in ("TAA", "TAG"):
                        raise Exception("Stop codon not TAA nor TAG!")
                    else:
                        raise Exception("It is not possible to end here!")

                if barcode_rc_from_insert in barcodes_used_list:
                    raise Exception("Barcode in ligation product occurs multiple times!", barcode_rc_from_insert)
                else:
                    barcodes_used_list.append(barcode_rc_from_insert)

                if epitope_variant_hamming_pool:
                    if (entry.name.split("_")[1] != "0") & (entry.name.split("_")[1] != "n"):
                        source_peptide_arr = np.array([AA_TO_ID[res] for res in list(source_peptide)[:1]] + [AA_TO_ID[res] for res in list(source_peptide)[2:8]])
                        epitope_variant_tmp_arr = np.array([AA_TO_ID[res] for res in list(insert_translation)[:1]] + [AA_TO_ID[res] for res in list(insert_translation)[2:8]])
                        # Fast hamming implementation using numpy,
                        # as the ARM tcr-toolbox-environment doesn't allow for scipy installation:
                        if not np.count_nonzero(source_peptide_arr != epitope_variant_tmp_arr) == int(entry.name.split("_")[1]):
                            raise Exception(
                                "Oligo circle name does not equal hamming circle of oligo insert AA translation!",
                                "Hamming circle:",
                                entry.name.split("_")[1],
                                "Hamming distance from source peptide:",
                                np.count_nonzero(source_peptide_arr != epitope_variant_tmp_arr),
                            )

                    elif entry.name.split("_")[1] == "n":
                        source_peptide_arr = np.array([AA_TO_ID[res] for res in list(source_peptide)[:1]] + [AA_TO_ID[res] for res in list(source_peptide)[2:8]])
                        epitope_variant_tmp_arr = np.array([AA_TO_ID[res] for res in list(insert_translation)[:1]] + [AA_TO_ID[res] for res in list(insert_translation)[2:8]])
                        # Fast hamming implementation using numpy,
                        # as the ARM tcr-toolbox-environment doesn't allow for scipy installation:
                        if not np.count_nonzero(source_peptide_arr != epitope_variant_tmp_arr) > 6:
                            raise Exception(
                                "Negative control epitopes are within 6 hamming distance from source epitope!",
                                epitope_variant,
                                "Hamming distance from source peptide:",
                                np.count_nonzero(source_peptide_arr != epitope_variant_tmp_arr),
                            )

    if (subpool_indices_dict is None) | also_check_full_pool:
        print("\nStarted checking full pool:")
        pcr_products_dict = collections.defaultdict(list)
        ligation_products_dict = collections.defaultdict(list)
        barcodes_used_list = []
        total_oligos = len(aa_check_list)
        with pysam.FastxFile(fasta_in_fname) as fasta_in:
            for i, (entry, epitope_variant, fusion_protein_bool) in enumerate(zip(fasta_in, aa_check_list, fusion_protein_bool_list)):
                if fusion_protein_bool:
                    print("Skipping fusion protein check! Check fusion protein fasta idx manually:", i)
                    continue  # hard to check translation
                _print_progress_update(count=i + 1, total=total_oligos, all_label="oligos are correct!", label="oligos correct!")

                if epitope_variant_hamming_pool:
                    if not epitope_variant == entry.name.split("_")[2]:
                        raise Exception(
                            "Epitope variant translation in fasta header does not equal epitope variant sequence in aa check list!",
                            "Epitope variant in fasta header:",
                            entry.name.split("_")[2],
                            "Epitope variant in aa check list:",
                            epitope_variant,
                        )

                primer_fw = primers_fw_list[0]
                primer_rv = primers_rv_list[0]

                pcr_products_dict[entry.name] = Anneal((Primer(primer_fw), Primer(primer_rv)), Dseqrecord(entry.sequence)).products

                if not len(pcr_products_dict[entry.name]) == 1:
                    raise Exception("Oligo pre-amplification either has too many products or no products:", pcr_products_dict[entry.name])

                if not len(pcr_products_dict[entry.name][0].seq) == len(entry.sequence):
                    raise Exception("Oligo pre-amplification amplicon has wrong size!", pcr_products_dict[entry.name])

                if not len(pcr_products_dict[entry.name][0].seq.cut(Restriction.BsmBI)) == 3:
                    raise Exception("BsmBI digested pre-amplified oligo has more than 3 products:", pcr_products_dict[entry.name][0].seq.cut(Restriction.BsmBI))

                ligation_products_dict[entry.name] = (p.cut(Restriction.BsmBI)[0] + pcr_products_dict[entry.name][0].seq.cut(Restriction.BsmBI)[1]).looped()

                if not ligation_products_dict[entry.name].circular:
                    raise Exception("Golden gate assembly into vector failed!")

                if p12_or_p20 == "p20":
                    # five_end starts at the end of the IL2-Signal Peptide (i.e., start of the epitope):
                    five_end = ligation_products_dict[entry.name].find(p12_p20_flanking_nt_dict["p20_end_5"]) + len(p12_p20_flanking_nt_dict["p20_end_5"])
                    # three_end is at the end of the barcode. So, three end is at the end of epitope + stop codon + barcode:
                    three_end = ligation_products_dict[entry.name].find(p12_p20_flanking_nt_dict["p20_end_3"])

                if p12_or_p20 == "p12":
                    # five_end starts at the end of the CD74-Signal Peptide (i.e., start of the epitope):
                    five_end = ligation_products_dict[entry.name].find(p12_p20_flanking_nt_dict["p12_end_5"]) + len(p12_p20_flanking_nt_dict["p12_end_5"])
                    # three_end is at the end of the barcode. So, three end is at the end of epitope + stop codon + barcode:
                    three_end = ligation_products_dict[entry.name].find(p12_p20_flanking_nt_dict["p12_end_3"])

                # Five_end should be 9059 and three_end 0:
                if p12_or_p20 == "p20" and not ((five_end == 9059) & (three_end == 0)):
                    raise Exception("Insert inserted at the wrong vector position!")

                # Five_end should be 9079 and three_end 0:
                if p12_or_p20 == "p12" and not ((five_end == 9079) & (three_end == 0)):
                    raise Exception("Insert inserted at the wrong vector position!")

                insert_region = ligation_products_dict[entry.name][five_end:three_end]
                insert_translation = insert_region.translate(to_stop=True)
                # Extract barcode and stop codon from insert using string slicing.
                # Avoids pydna circular DNA position arithmetic which breaks when three_end=0.
                insert_seq_str = str(insert_region)
                barcode_rc = entry.name.split("_")[0]
                barcode_len = len(barcode_rc)
                barcode_in_insert = insert_seq_str[-barcode_len:]
                stop_in_insert = insert_seq_str[-(barcode_len + 3) : -barcode_len]
                barcode_rc_from_insert = str(Seq(barcode_in_insert).reverse_complement())

                if not (epitope_variant == insert_translation) & (barcode_rc == barcode_rc_from_insert) & ((stop_in_insert == "TAA") | (stop_in_insert == "TAG")):
                    if epitope_variant != insert_translation:
                        raise Exception(
                            "idx:",
                            str(i),
                            entry.name,
                            "translation in ligation_product does not match original AA sequence!",
                            "Original AA sequence:",
                            epitope_variant,
                            "Translation in ligation product:",
                            insert_translation,
                        )
                    elif barcode_rc != barcode_rc_from_insert:
                        raise Exception(
                            "Barcode in ligation_product does not match original Barcode sequence!",
                            "Original Barcode sequence:",
                            barcode_rc,
                            "Barcode in ligation product (RC):",
                            barcode_rc_from_insert,
                        )
                    elif stop_in_insert not in ("TAA", "TAG"):
                        raise Exception("Stop codon not TAA nor TAG!")
                    else:
                        raise Exception("It is not possible to end here!")

                if barcode_rc_from_insert in barcodes_used_list:
                    raise Exception("Barcode in ligation product occurs multiple times!", barcode_rc_from_insert)
                else:
                    barcodes_used_list.append(barcode_rc_from_insert)

                if epitope_variant_hamming_pool:
                    if (entry.name.split("_")[1] != "0") & (entry.name.split("_")[1] != "n"):
                        source_peptide_arr = np.array([AA_TO_ID[res] for res in list(source_peptide)[:1]] + [AA_TO_ID[res] for res in list(source_peptide)[2:8]])
                        epitope_variant_tmp_arr = np.array([AA_TO_ID[res] for res in list(insert_translation)[:1]] + [AA_TO_ID[res] for res in list(insert_translation)[2:8]])
                        # Fast hamming implementation using numpy,
                        # as the ARM tcr-toolbox-environment doesn't allow for scipy installation:
                        if not np.count_nonzero(source_peptide_arr != epitope_variant_tmp_arr) == int(entry.name.split("_")[1]):
                            raise Exception(
                                "Oligo circle name does not equal hamming circle of oligo insert AA translation!",
                                "Hamming circle:",
                                entry.name.split("_")[1],
                                "Hamming distance from source peptide:",
                                np.count_nonzero(source_peptide_arr != epitope_variant_tmp_arr),
                            )

                    elif entry.name.split("_")[1] == "n":
                        source_peptide_arr = np.array([AA_TO_ID[res] for res in list(source_peptide)[:1]] + [AA_TO_ID[res] for res in list(source_peptide)[2:8]])
                        epitope_variant_tmp_arr = np.array([AA_TO_ID[res] for res in list(insert_translation)[:1]] + [AA_TO_ID[res] for res in list(insert_translation)[2:8]])
                        # Fast hamming implementation using numpy,
                        # as the ARM tcr-toolbox-environment doesn't allow for scipy installation:
                        if not np.count_nonzero(source_peptide_arr != epitope_variant_tmp_arr) > 6:
                            raise Exception(
                                "Negative control epitopes are within 6 hamming distance from source epitope!",
                                epitope_variant,
                                "Hamming distance from source peptide:",
                                np.count_nonzero(source_peptide_arr != epitope_variant_tmp_arr),
                            )

    return print("Finished checking all oligos!")


def write_minigene_aa_fasta_from_oligo_order_fasta(
    oligo_order_fasta_out_file: Union[str, os.PathLike[str]], fasta_out_fname: Union[str, os.PathLike[str]], p12_or_p20: str = "p20"
):
    """Extract and translate the epitope insert from each oligo-order entry into an amino-acid fasta.

    For each entry in `oligo_order_fasta_out_file`, splits the oligo sequence on the
    vector's 5' and 3' BsmBI cloning sites (selected via `p12_or_p20`) to isolate the
    nucleotide region between them (epitope coding sequence, stop codon, and
    barcode), translates it up to (and excluding) the first stop codon — which
    trims the translation down to just the epitope — and writes the resulting
    amino acid sequence under the same fasta header to `fasta_out_fname`. This fasta file is useful for verification of the actual amino acid
    sequence that was screened, for instance, the actual amino acid sequence of a recognized
    Ag.


    Parameters
    ----------
    oligo_order_fasta_out_file : str | os.PathLike[str]
        Path to the oligo-order fasta to read (as written by
        `write_barcoded_p12_p20_epitope_oligo_order_fasta`).
    fasta_out_fname : str | os.PathLike[str]
        Output path for the resulting amino-acid fasta.
    p12_or_p20 : str, default = "p20"
        Which vector's BsmBI cloning sites delimit the epitope coding sequence in
        each oligo; must be "p12" or "p20".

    Returns
    -------
    None
        Writes `fasta_out_fname` as a side effect.

    Raises
    ------
    ValueError
        If `p12_or_p20` is not "p12" or "p20".
    """
    if p12_or_p20 not in {"p12", "p20"}:
        raise ValueError(f"p12_or_p20 must be 'p12' or 'p20' and not {p12_or_p20}")

    with pysam.FastxFile(oligo_order_fasta_out_file) as fasta_in, open(fasta_out_fname, "w") as fasta_out:
        for entry in fasta_in:
            print(">" + entry.name, file=fasta_out)
            if p12_or_p20 == "p20":
                aa_seq = Seq(entry.sequence.split(p20_cloning_sites["end_5_BsmBI_p20"])[1].split(p20_cloning_sites["end_3_BsmBI_p20"])[0]).translate(to_stop=True)
            if p12_or_p20 == "p12":
                aa_seq = Seq(entry.sequence.split(p12_cloning_sites["end_5_BsmBI_p12"])[1].split(p12_cloning_sites["end_3_BsmBI_p12"])[0]).translate(to_stop=True)

            print(aa_seq, file=fasta_out)
