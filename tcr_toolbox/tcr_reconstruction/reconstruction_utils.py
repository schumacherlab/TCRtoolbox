"""
This script is imputes the amino acid (aa) and nucleotide (nt) sequences from IMGT.
The constants TRAV, TRAJ, TRBV, TRBV and TRBJ are all types of notations of VDJ segments found in the database.
These entries are added to a translation dict and checked for direct matches in the IMGT fasta files.
Hereafter a series of mutations to the entries are done to see if the database entries will match the IMGT standard format.
The resulting translation dictionaries translate database V, D and J entries to aa and nt sequence.
"""

from collections import defaultdict
from typing import Union
import os
import re
from datetime import datetime
import json
import numpy as np
import pandas as pd
import logging
from collections import Counter

from tcr_toolbox.tcr_reconstruction.constants import CONSTANT_OVERLAP_NTS_HUMAN, CONSTANT_OVERLAP_NTS_MOUSE
from tcr_toolbox.tcr_reconstruction.constants import TRANSLATION_TABLE, VDJ_NAME, CONSTANT_FIRST_AAS
from tcr_toolbox.tcr_reconstruction.constants import START_CODON, STOP_CODONS,AMBIGUOUS_AA_CODE, VDJ_CDR3_ALL_END_MOTIFS
from tcr_toolbox.tcr_reconstruction.constants import VDJ_MAX_CDR3_LEN, VDJ_MIN_CDR3_LEN, TRAV, TRAJ, TRBV, TRBD, TRBJ
from tcr_toolbox.utils.constants import VDJ_CDR3_COMMON_END_MOTIFS, VDJ_CDR3_RARE_END_MOTIFS, AA_STRING, each_codon_to_aa
from tcr_toolbox.utils.logger import init_logger

logger = logging.getLogger(__name__)


def reconstruct_full_tcr(
    df: pd.DataFrame,
    v_column_nt: str,
    v_column_aa: str,
    j_column_nt: str,
    j_column_aa: str,
    cdr3_column_aa: str,
    include_leader=False,
    exclude_C_F_W=False,
    include_constant=False,
    constant_sequence=None,
    human_or_mouse_constant=None,
    verbose: Union[bool, int] = False,
):
    """
    Reconstructs the full TCR sequence from the V, J and CDR3 columns.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe containing the V, J and CDR3 columns.
    v_column_nt : str
        The name of the V column containing the nucleotide sequence.
    v_column_aa : str
        The name of the V column containing the amino acid sequence.
    j_column_nt : str
        The name of the J column containing the nucleotide sequence.
    j_column_aa : str
        The name of the J column containing the amino acid sequence.
    cdr3_column_aa : str
        The name of the CDR3 column containing the amino acid sequence.
    include_leader : bool, optional
        If True, the leader sequence will be included in the full sequence.
        The default is True.
    exclude_C_F_W : bool, optional
        If True, the CDR3 seqeuences without 1st C, and last F and W amino acids will be excluded from the full sequence.
        The default is False.
    include_constant : bool, optional
        If True, the constant sequence will be included in the full sequence.
        The default is False.
    constant_sequence : str, optional
        The constant sequence to include in the full sequence.
        The default is None.
    human_or_mouse_constant : str, optional
        The type of constant sequence (TRAC or TRBC).
        The default is None.
    verbose : Union[bool, int], optional
        If True, the function will print out the progress.
        If 2 the detailed progress will be printed out.
        The default is False.

    Returns
    -------
    full_sequence_column : list
        List of the full TCR sequences.

    """

    # df = df.reset_index(drop=True) # using iteritems, so no need to reset index
    cdr3_column_aa = df[cdr3_column_aa]
    v_column_nt = df[v_column_nt]
    v_column_aa = df[v_column_aa]
    j_column_nt = df[j_column_nt]
    j_column_aa = df[j_column_aa]

    full_sequence_column = []
    na_counter = Counter()
    missing_index_counter = Counter()
    for i, row in df.iterrows():
        cdr3 = cdr3_column_aa[i]

        if pd.isna(cdr3):
            na_counter.update(["cdr3_na"])
            full_sequence_column.append(np.nan)
            continue

        # conditionals to exclude sequences
        if cdr3[0] != "C":
            logger.debug(f"BE AWARE: the cdr3: {cdr3} at index {i} did not start with a 'C'.")
            na_counter.update(["no_C_at_start"])
            if exclude_C_F_W:
                logger.debug(f"The full sequence at index {i} will not be reconstructed!")
                full_sequence_column.append(np.nan)
                continue

        if cdr3[-1] != "F" and cdr3[-1] != "W":
            logger.debug(f"BE AWARE: the cdr3: {cdr3} at index {i} did not end with a 'F' or 'W'.")
            na_counter.update(["did_not_end_with_F_or_W"])
            if exclude_C_F_W:
                full_sequence_column.append(np.nan)
                logger.debug(f"The full sequence at index {i} will not be reconstructed!")
                continue

        # check for None and nan values in the v_column_nt
        if include_leader and pd.isna(v_column_nt[i]):
            logger.debug(f"missing V nt sequence, full sequence at index {i} will not be reconstructed!")
            na_counter.update(["missing_v_nt"])
            full_sequence_column.append(np.nan)
            continue

        if pd.isna(cdr3_column_aa[i]) or pd.isna(j_column_nt[i]) or pd.isna(j_column_aa[i]) or pd.isna(v_column_aa[i]):
            missing_index_counter.update([i])
            na_counter.update(["missing_j_nt_or_cdr3_aa_or_j_aa"])
            full_sequence_column.append(np.nan)
            logger.debug(f"missing V,J or cdr3 sequence, full sequence at index {i} will not be reconstructed!")
            continue

        # start of reconstruction
        if include_leader:
            # find C104 and translate the leader+v
            if not re.match(r"^[ATCG]+$", v_column_nt[i].upper()):
                logger.debug(f"BE AWARE: the v sequence: {v_column_nt[i]} at index {i} contained non ATCG characters.")
                na_counter.update(["non_atcg_in_v"])
                full_sequence_column.append(np.nan)
                continue

            trxv_until_c_104 = trim_until_c104(v_column_nt[i].upper(), nt_or_aa="nt")
            trxv_aa = translate_nt_to_aa(trxv_until_c_104)
        else:
            # the amino acid sequence is not including the leader
            # check for only aa characters
            all_amino_acids = AA_STRING
            if not all([aa in all_amino_acids for aa in v_column_aa[i].upper()]):
                logger.debug(f"BE AWARE: the v sequence: {v_column_aa[i]} at index {i} contained non amino acid characters.")
                na_counter.update(["non_aa_in_v"])
                full_sequence_column.append(np.nan)
                continue

            trxv_aa = trim_until_c104(v_column_aa[i].upper(), nt_or_aa="aa")

        # add v and cdr3 together
        trxv_plus_cdr3 = trxv_aa + cdr3

        # Cannot use the script to match CDR3 end motifs in J segments, because the
        # "TRBJ2-4*01" and "TRBJ1-3*01" have multiple matches in multiple reading frames, both coding for
        # valid amino acid sequences. therefore the aa seq will be used to find the pattern
        trxj_aa = j_column_aa[i]
        index_j = find_cdr3_end_motif(trxj_aa, VDJ_CDR3_ALL_END_MOTIFS)
        full_seq = trxv_plus_cdr3 + trxj_aa[index_j + 1 :]

        if include_constant:
            full_seq = include_constant_region(full_seq=full_seq, j_seq_nt=j_column_nt[i], constant_sequence=constant_sequence, human_or_mouse_constant=human_or_mouse_constant)

        # full_seq = trxv_plus_cdr3
        full_sequence_column.append(full_seq)

    # logger.warning(f"did not reconstruct (reason:count) : {na_counter}")
    if verbose:
        print(f"did not reconstruct (reason:count) : {na_counter}")
        if verbose == 2:
            print(f"missing indices: {missing_index_counter.keys()}")

    return full_sequence_column


def include_constant_region(full_seq, j_seq_nt, constant_sequence, human_or_mouse_constant):
    if constant_sequence is None:
        raise ValueError("No constant sequence provided!, despite include_constant=True")

    constant_sequence = constant_sequence.upper()
    # check if the constant sequence are nucleotides:
    constant_sequence_tokens = [token for token in constant_sequence]
    if set(constant_sequence_tokens) == {"A", "T", "C", "G"}:
        # translate to aa
        constant_sequence = translate_nt_to_aa(constant_sequence)

    # The aminoacid between the j an constant region of the beta TCR
    # is always E (g from J-gene en ag from constant gene both human and mouse)
    if constant_sequence[0:2] == "DL":
        full_seq = full_seq + "E" + constant_sequence
    elif constant_sequence[0:3] == "EDL":
        full_seq = full_seq + constant_sequence
    # The amino acid between the j and constant of the alpha TCR
    # is ambiguous and will be imputed below
    elif constant_sequence[0:3] == "IQN":
        trac_or_trbc = CONSTANT_FIRST_AAS[constant_sequence[0:3]]
        last_nt = j_seq_nt[-1]
        j_constant_overlap_amino_acid = impute_j_constant_amino_acid(last_nt=last_nt, trac_or_trbc=trac_or_trbc, human_or_mouse_constant=human_or_mouse_constant)
        full_seq = full_seq + j_constant_overlap_amino_acid + constant_sequence

    else:
        raise ValueError("The first amino acids of the constant sequence are not (E)DL (TRBC) or IQN (TRAC)")

    return full_seq


def impute_j_constant_amino_acid(last_nt, trac_or_trbc, human_or_mouse_constant):
    """
    Imputes the amino acid of the J and constant overlap, based on the last nucleotide of the J gene and the type of
    constant sequence (TRAC or TRBC)

    Parameters
    ----------
    last_nt : str
        The last nucleotide of the J gene.
    trac_or_trbc : str
        The type of constant sequence (TRAC or TRBC)
    human_or_mouse_constant : str
        The type of constant sequence (human or mouse)

    Returns
    -------
    j_constant_overlap_amino_acid : str
        The imputed amino acid of the J-constant overlap.

    """
    if human_or_mouse_constant == "human":
        constant_overlap_nts = CONSTANT_OVERLAP_NTS_HUMAN
    elif human_or_mouse_constant == "mouse":
        constant_overlap_nts = CONSTANT_OVERLAP_NTS_MOUSE
    else:
        raise ValueError("the type for constant must be specified as either 'human' or 'mouse'")

    j_constant_overlap_nts = last_nt.upper() + constant_overlap_nts[trac_or_trbc]
    j_constant_overlap_amino_acid = translate_nt_to_aa(j_constant_overlap_nts)

    return j_constant_overlap_amino_acid


def translate_nt_to_aa(dna):
    """
    Translates a DNA sequence into an amino acid sequence.

    Parameters
    ----------
    dna : str
        The DNA sequence to translate.

    Returns
    -------
    str
        The translated amino acid sequence.
    """
    protein = []
    end = len(dna) - (len(dna) % 3) - 1
    dna = dna.upper()

    for i in range(0, end, 3):
        codon = dna[i : i + 3]

        aminoacid = TRANSLATION_TABLE.get(codon, AMBIGUOUS_AA_CODE)
        protein.append(aminoacid)

    return "".join(protein)


def trim_until_c104(v_seq, nt_or_aa):
    """
    Finds the first C in the V gene sequence.

    Parameters
    ----------
    v_seq : str
        The V gene sequence.
    nt_or_aa : str
        The sequence is in nucleotides or amino acids ('aa' or 'nt').

    Returns
    -------
    str
        The sequence until the first C in the V gene sequence.
    """
    v_seq = v_seq.upper()

    if nt_or_aa == "nt":
        # assert if the sequence starts with ATG,
        assert v_seq[0:3] == START_CODON
        # assert if only nt characters are present
        all_nucleotides = {"A", "T", "C", "G"}
        if not all([nt in all_nucleotides for nt in v_seq]):
            raise ValueError("The V sequence contains non-nucleotide characters.")

        # the reading frame is determined to remove the last nts
        residual_nts = len(v_seq) % 3
        reading_frame_v = v_seq[: len(v_seq) - residual_nts]

        # find the first C going from the last codon to the left
        for codon_end in range(len(reading_frame_v), 0, -3):
            if TRANSLATION_TABLE[v_seq[codon_end - 3 : codon_end]] == "C":
                # return the sequence without the C, as this is the starting AA from the CDR3
                return v_seq[: codon_end - 3]

    elif nt_or_aa == "aa":
        # assert if only amino acid characters are present
        all_amino_acids = AA_STRING
        if not all([aa in all_amino_acids for aa in v_seq]):
            raise ValueError("The V sequence contains non-amino acid characters.")
        # find the first C going from the last codon to the left
        for index_c104 in range(len(v_seq) - 1, -1, -1):
            if v_seq[index_c104] == "C":
                # return the sequence without the C, as this is the starting AA from the CDR3
                return v_seq[:index_c104]


def find_cdr3_end_motif(j_amino_acids, allowed_end_motifs):
    """
    Search the CDR3 end motif in the amino acid sequence of J

    Parameters
    ----------
    j_amino_acids : str
        The amino acids of the J gene.
    allowed_end_motifs : list
        The list of allowed end motifs.

    Returns
    -------
    int
        The index of the end of the CDR3 region.
    """

    for motif in allowed_end_motifs:
        assert motif in VDJ_CDR3_ALL_END_MOTIFS
        for idx in range(len(j_amino_acids) - len(motif) + 1):
            valid_end = True
            for i, aa in enumerate(motif):
                valid_end = valid_end and (j_amino_acids[idx + i] == aa or aa == "X")
            if valid_end:
                if idx != 0:  # Fix to prevent the first hit to mess up the allignment  because of double end motifs in TRAJ35 and TRAJ16
                    return idx
                else:
                    continue
    return None



def reconstruct_vdj(vdj_column, vdj_type, translation_dict_path_nt, translation_dict_path_aa, verbose=False):
    """
    This function is used in the refinement of the database in the DataParser class (refine_dataset()).
    It matches the raw vdj allele names to the IMGT database using a translation
    dictionary (constructed in vdj_construction_utilt.py). This dictionary should be kept up to date if new data is added.
    :param vdj_column: column containing the v, d and j segments and v and j segments
    :param vdj_type: TRAV, TRAJ, TRBV, TRBD or TRBJ as a string
    :return: IMGT standardized name for aa and nt and
    """

    # RECONSTRUCTION FULL SEQ test script path
    nt_path = translation_dict_path_nt
    aa_path = translation_dict_path_aa

    missing_vdj_counter = {vdj_type + "_nt": Counter(), vdj_type + "_aa": Counter()}

    with open(aa_path) as aa_dict:
        translation_dict_aa = json.load(aa_dict)
    with open(nt_path) as nt_dict:
        translation_dict_nt = json.load(nt_dict)

    parsed_vdj_column_aa = []
    parsed_vdj_column_nt = []
    full_seq_vdj_column_aa = []
    full_seq_vdj_column_nt = []

    for vdj_value in vdj_column[vdj_type]:
        try:
            parsed_vdj_column_aa.append(translation_dict_aa[vdj_type][vdj_value][0])
            full_seq_vdj_column_aa.append(translation_dict_aa[vdj_type][vdj_value][1])
            # logger.debug('MATCHED {} in the amino acid dictionary'.format(vdj_value))
        except (TypeError, KeyError):
            # logger.debug('could not match {} in the amino acid dictionary'.format(vdj_value))
            parsed_vdj_column_aa.append(np.nan)
            full_seq_vdj_column_aa.append(np.nan)
            missing_vdj_counter[vdj_type + "_aa"].update([vdj_value])

        try:
            parsed_vdj_column_nt.append(translation_dict_nt[vdj_type][vdj_value][0])
            full_seq_vdj_column_nt.append(translation_dict_nt[vdj_type][vdj_value][1])
            # logger.debug('MATCHED {} in the nucleotide dictionary'.format(vdj_value))
        except (TypeError, KeyError):
            # logger.debug('could not match {} in the nucleotide dictionary'.format(vdj_value))
            parsed_vdj_column_nt.append(np.nan)
            full_seq_vdj_column_nt.append(np.nan)
            missing_vdj_counter[vdj_type + "_nt"].update([vdj_value])

    if verbose:
        print("Missing {} in the amino acid dictionary:".format(vdj_type))
        print(missing_vdj_counter[vdj_type + "_aa"])
        print("Missing {} in the nucleotide dictionary:".format(vdj_type))
        print(missing_vdj_counter[vdj_type + "_nt"])

    return (parsed_vdj_column_aa, full_seq_vdj_column_aa, parsed_vdj_column_nt, full_seq_vdj_column_nt)


def IMGT_to_vdj_dict(path):
    from pysam import FastaFile

    """
    Initializes the dictionary used for translating to IMGT standard notation (uses downloaded fasta files).
    Returns a dictionary of IMGT-formatted FASTA sequences, where the key is the VDJ allele name (e.g. TRAV1-1*01)
    and the value is a nested dictionary of sequences, where the key is the sequence name (e.g. TRAV1-1*01|TRAJ9*01)
    and the value is the sequence.

    :param path: Path to the directory containing the IMGT-formatted FASTA files.
    :type path: str
    :return: A dictionary of IMGT-formatted FASTA sequences.
    :rtype: dict
    """
    vdj_seqs = defaultdict(dict)

    # loop over all files and fetch sequences for each vdj allele (aa or nt)
    for filename in os.listdir(path):
        if filename.endswith(".fasta"):
            vdj_name_1 = filename[0:4] + "_" + filename[-8:-6]  # 0:4 is the TRAV, -8:-6 is nt or aa
            vdj_seqs[vdj_name_1] = {}

            # add to nested dict (ie: dict[TRAJnt][TRAJ9*01]: sequence)
            with FastaFile(filename) as fasta_file:
                for ref in fasta_file.references:
                    vdj_seqs[vdj_name_1][ref.split("|")[1]] = fasta_file.fetch(reference=ref)  # middle of ref is TRAJ9*01 etc

    # revert to normal dict, to avoid accidentally adding empty entries that are not from IMGT during imputing
    return dict(vdj_seqs)


def make_vdj_segment_translation_dict():
    """
    This function creates a dictionary of dictionaries for the V, D, and J gene segment annotations found in the database.
    :return: A dictionary of dictionaries for the V, D, and J gene segments.
    """
    VDJ_DICT = defaultdict(dict)
    for i, VDJ in enumerate([TRAV, TRAJ, TRBV, TRBD, TRBJ]):
        VDJ_DICT[VDJ_NAME[i]] = defaultdict(list)
        for entry in VDJ:
            entry = entry.strip()  # remove the spaces around the entries
            VDJ_DICT[VDJ_NAME[i]][entry] = None

    return VDJ_DICT


def impute_missing_vdj(vdj_dict, aa_nt, imgt):
    """
    Impute the missing annotations for the V D and J.
    :param vdj_dict: annotations as keys and sequence as values
    :param aa_nt: amino acid or nucleotide dictionary
    :param imgt: imgt gold standard dict with annotations as keys and sequence as values
    :return: vdj translation dict with imputed annotations and corresponding sequences
    """
    # for TRAV, TRAJ, TRBV, TRBD, TRBJ:
    for vdj_ab in VDJ_NAME:
        logger.info("Imputing missing keys for: {0}{1} #########################################\n\n".format(vdj_ab, aa_nt))
        for entry in vdj_dict[vdj_ab]:
            # if not none: the entry is already present in the correct IMGT format
            if vdj_dict[vdj_ab][entry] is not None:
                logger.debug(" entry exists for: {}".format(entry))
                continue

            else:
                if not re.match(vdj_ab, entry):
                    logger.debug("\tentry: {0} did not start with: {1}, adding None to entry in dict.".format(entry, vdj_ab))
                    vdj_dict[vdj_ab][entry] = None

                else:
                    suffix = entry[4:]

                    # if first 2 numbers match to a single key in IMGT, match them:
                    matches = [[key, val] for key, val in imgt[vdj_ab + aa_nt].items() if key.startswith(vdj_ab + suffix[0:2])]

                    if len(matches) == 1:
                        # un-nest NESTED dict because of list comprehension
                        vdj_dict[vdj_ab][entry] = [matches[0][0], matches[0][1]]
                        logger.debug("\t{0} first 2 digits exactly 1 key: {1}".format(entry, matches[0][0]))
                        continue

                    else:
                        # match entry + *01 and if match return key[:-3], which gives the key without *01 info.
                        with_01 = [key[:-3] for key, val in imgt[vdj_ab + aa_nt].items() if key == entry + "*01"]
                        # with_02 = [key[:-3] for key, val in imgt[vdj_ab + aa_nt].items() if key == entry + '*02']
                        # hits_with_01_not_02 = [entry_with_01 for entry_with_01 in with_01 if entry_with_01 not in with_02]
                        if len(with_01) > 0:
                            fix = entry + "*01"
                            vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                            logger.debug("\t added *01 to {0} and found an exact match for: {1}".format(entry, fix))
                            continue

                        # some entries have -0 when this can never be the case (ie should be TRAV1-01 istead of TRAV1-1)
                        if re.search("-0", entry):
                            fix = entry.replace("-0", "-")
                            try:
                                vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                logger.debug("\t replaced -0 for - for entry: {0}, to {1}".format(entry, fix))
                                continue
                            except KeyError:
                                try:
                                    fix = entry.replace("-", "*")
                                    vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                    logger.debug("\t replaced - for * for entry: {0}, to {1}".format(entry, fix))
                                    continue
                                except KeyError:
                                    try:
                                        fix = entry.replace("0", "", 1) + "*01"
                                        vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                        logger.debug("\t replaced -0 for - , and added *01 for entry: {0}, to {1}".format(entry, fix))
                                        continue
                                    except KeyError:
                                        try:
                                            fix = entry.replace("0", "", 2)
                                            vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                            logger.debug("removed first and 2nd 0 for entry: {0}, fix: {1}".format(entry, fix))
                                            continue
                                        except KeyError:
                                            logger.debug("Could not fix -0 for {}".format(entry))
                                            continue

                        # an entry cannot begin with TRAV01, as it should be TRAV1.
                        if re.search(r"[TRAV]0\d", entry):
                            fix = entry.replace("0", "", 1)
                            try:
                                vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                logger.debug("removed 0 for {0}, fix: {1}".format(entry, fix))
                                continue
                            except KeyError:
                                try:
                                    fix = entry.replace("0", "", 1) + "*01"
                                    vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                    logger.debug("removed 0 and added *01 for {0}, fix: {1}".format(entry, fix))
                                    continue
                                except KeyError:
                                    try:
                                        fix = entry.replace("0", "", 2)
                                        vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                        logger.debug("removed first and 2nd 0 {0}, fix: {1}".format(entry, fix))
                                        continue
                                    except KeyError:
                                        try:
                                            fix = entry.replace("0", "", 2) + "*01"
                                            vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                            logger.debug("removed first and 2nd 0 and added *01 {0}, fix: {1}".format(entry, fix))
                                            continue
                                        except KeyError:
                                            logger.debug("adding zeros did not work for {0}".format(entry))

                        if re.search(r":", suffix):  # entries with a ':' are changed to a *
                            suffix = suffix.replace(":", "*")
                            suffix = suffix.strip(",")
                            fix = vdj_ab + suffix
                            try:
                                vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                logger.debug("\treplacing ':' for '*' success for: {0} -> {1}".format(entry, fix))
                                continue
                            except KeyError:
                                logger.debug("\tmissed exception replacing : for * {0} -> {1}".format(entry, fix))

                        # find the /DV exceptions
                        if re.search("/DV", entry):
                            if re.search(" F", entry):
                                try:
                                    fix = entry.strip(" F") + "*01"
                                    vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                    logger.debug('\tstripping " F" successfull for: {0}, {1}'.format(entry, fix))
                                    continue
                                except KeyError as key_error:
                                    logger.debug('\tfailed stripping " F"  for {0}, {1}'.format(entry, key_error))

                        if re.search("DV", entry):
                            if re.search("..", entry):
                                fixing = entry.strip(".")
                                logger.debug("\tthis entry: {} ,has a dot, stripping it now".format(entry))
                            else:
                                fixing = entry

                            try:  # without (also adding) *01 had no successes
                                fix = fixing.replace("DV", "/DV") + "*01"
                                vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                logger.debug("\t replacing DV for /DV succsesfully for : {0}, {1}".format(entry, fix))
                                continue
                            except KeyError as key_error:
                                logger.debug("failed replacing DV for /DV for: {0}, {1}".format(entry, key_error))

                        # replace - for *
                        if bool(re.search(r"-", suffix)):
                            suffix = suffix.replace("-", "*")
                            fix = vdj_ab + suffix
                            try:
                                vdj_dict[vdj_ab][entry] = [fix, imgt[vdj_ab + aa_nt][fix]]
                                logger.debug("\treplacing - for * success for: {0} -> {1}".format(entry, fix))
                                continue
                            except KeyError:
                                logger.debug("\treplacing - for * FAILED for: {0} -> {1}".format(entry, fix))

                    if len(matches) > 1:
                        logger.debug("\t{} first 2 digits matches more than one key".format(entry))
                        # continue
                    if len(matches) == 0:
                        logger.debug("\t{} first 2 digits matches no keys".format(entry))
                        # continue

    return vdj_dict


def match_IMGT_to_translation_dict(IMGT, translation_dict, aa_nt):
    """
    Match all the annotations in the database to IMGT gold standard anntotations and extract corresponding sequences.
    This dictionary will be expanded on in the impute_missing() function)

    :param IMGT: IMGT dictionary with VDJ gene annotations as keys and sequence as value
    :param translation_dict: annotations from database as keys and empty values
    :param aa_nt: amino acid or nucleotide
    :return: translation dict, where exact matches to IMGT contain the corresponding sequence
    """
    logger.info(" Making {} translation dict: ##############################################\n\n".format(aa_nt))
    for vdj in VDJ_NAME:
        logger.info("\t making translation dict for: {}".format(vdj))
        logger.debug("\tlength translation dict {0} before adding new keys: {1}".format(vdj, len(translation_dict[vdj])))
        for entry in IMGT[vdj + aa_nt]:
            if entry != "nan":
                try:
                    translation_dict[vdj][entry] = [entry, IMGT[vdj + aa_nt][entry]]

                except KeyError as key_error:  # should not be needed, because it uses a default dict!
                    logger.debug("KEYERROR: {} adding as a new entry to the translation_dict".format(key_error))
                    translation_dict[vdj][entry] = None
        logger.debug("\tlength translation dict {0} after adding new keys: {1}\n".format(vdj, len(translation_dict[vdj])))

        # return to a normal dict for the next step
        translation_dict[vdj] = dict(translation_dict[vdj])

    return dict(translation_dict)


def codon_to_aa_ord(codon: Union[str, bytes]) -> int:
    """Return amino acid corresponding to a codon.

    If the codon is not in the translation table, a default
    AA is returned (see AMBIGUOUS_AA_CODE in vdj.constants)
    """
    assert len(codon) == 3
    if isinstance(codon, bytes):
        codon = codon.decode("ascii")
    assert all(c in "NACGT" for c in codon)
    code = each_codon_to_aa.get(codon, AMBIGUOUS_AA_CODE)
    assert len(code) == 1
    return ord(code)


# From: https://github.com/10XGenomics/cellranger/blob/7c9ccea73e3b619687fd7b3e426ebc4640a7d653/lib/python/cellranger/vdj/annotations.py#L278
def search_cdr3_signature_no_vj(seq):
    """
    Search the CDR3 signature in a sequence without guides from annotations.

    This could lead to more false positive signature hits than the guided version.

    Return value:
    A tuple (CDR3 DNA seq, CDR3 amino-acid seq, start position in seq, end position)
    """

    min_cdr3_aas = VDJ_MIN_CDR3_LEN // 3

    cys_pos = 0
    for frame in range(3):
        amino_acids = bytes(bytearray(codon_to_aa_ord(seq[i : (i + 3)]) for i in range(frame, len(seq) - 2, 3)))

        for idx in range(min_cdr3_aas, len(amino_acids) - 3):
            # First try to find the end motif
            for motif in VDJ_CDR3_COMMON_END_MOTIFS + VDJ_CDR3_RARE_END_MOTIFS:
                valid_end = all(amino_acids[idx + i] == ord(aa) or aa == "X" for i, aa in enumerate(motif))

                if valid_end:
                    # The CDR3 includes the first F of the signature
                    fgxg_idx = idx
                    fgxg_pos = frame + fgxg_idx * 3

                    cys_idx = None
                    # Find the Cysteine closer to the end but past a minimum
                    # number of amino acids
                    for idx in range(fgxg_idx - min_cdr3_aas, 0, -1):
                        if amino_acids[idx] == ord(b"C"):
                            cys_idx = idx
                            cys_pos = frame + cys_idx * 3
                            break
                        elif chr(amino_acids[idx]) in STOP_CODONS:
                            break

                    if cys_idx:
                        if fgxg_pos - cys_pos < VDJ_MAX_CDR3_LEN:
                            # include the start of the FGXG
                            cdr3_seq = seq[cys_pos : (fgxg_pos + 3)]
                            cdr3_aas = amino_acids[cys_idx : (fgxg_idx + 1)]
                            return (cdr3_seq, cdr3_aas, cys_pos, fgxg_pos)

    return None


# if __name__ == "__main__":
#     """
#     This script will create the translation dictionaries and impute missing annotations.
#     """

#     current_datetime = str(datetime.now())[0:13].replace(" ", "_") + "h_"  # only take date and hour
#     current_datetime = current_datetime.replace(":", "-")

#     with open("imgt_to_py_dict_after_benchmark.json") as imgt_file:
#         imgt = json.load(imgt_file)

#     # make the aa translation dictionary
#     logger = init_logger("vdj_parsing_aa.log")
#     vdj_translation_dict_before_aa = make_vdj_segment_translation_dict()
#     vdj_translation_dict_aa = match_IMGT_to_translation_dict(imgt, vdj_translation_dict_before_aa, "_aa")
#     vdj_translation_dict_aa = impute_missing_vdj(vdj_translation_dict_aa, "_aa", imgt=imgt)
#     with open("vdj_translation_dict_aa_" + current_datetime + ".json", "w", encoding="utf-8") as vdj_transl:
#         json.dump(vdj_translation_dict_aa, vdj_transl, ensure_ascii=False, indent=4)

#     # make the nt translation dictionary
#     logger = init_logger("vdj_parsing_nt.log")
#     vdj_translation_dict_before_nt = make_vdj_segment_translation_dict()
#     vdj_translation_dict_nt = match_IMGT_to_translation_dict(imgt, vdj_translation_dict_before_nt, "_nt")
#     vdj_translation_dict_nt = impute_missing_vdj(vdj_translation_dict_nt, "_nt", imgt=imgt)
#     with open("vdj_translation_dict_nt_" + current_datetime + ".json", "w", encoding="utf-8") as vdj_transl:
#         json.dump(vdj_translation_dict_nt, vdj_transl, ensure_ascii=False, indent=4)


