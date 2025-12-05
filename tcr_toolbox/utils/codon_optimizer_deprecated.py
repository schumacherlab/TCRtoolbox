from typing import List

import numpy as np
from Bio.Restriction import Restriction
from numpy.random import default_rng
from re import search

# import pandas as pd
from Bio.Seq import Seq

# from pydna.dseq import Dseq
# from Bio import SeqIO
# from pydna.amplify import pcr
# from pydna.primer import Primer


from tcr_toolbox.utils.constants import codon_1_to_codon_2, aa_to_codon, codon_to_aa, codon_2_to_codon_1


def make_search_strings(enzymes: List = None, sites_to_check: dict = None):
    if enzymes is None or enzymes == []:
        enzymes = []
        print("No enzymes sites will be searched!!")
    if sites_to_check is None or sites_to_check == {}:
        sites_to_check = {}
        print("No sites_to_check will be searched!")

    search_string = list(sites_to_check.values())[0]
    for site in list(sites_to_check.values())[1:]:
        rev_comp = str(Seq(site).reverse_complement())
        search_string += "|" + site + "|" + rev_comp

    restriction_enzymes = {}
    for enzyme in enzymes:
        restriction_enzymes[enzyme] = getattr(Restriction, enzyme)

    restriction_enzyme_sites = ""
    restriction_edge_sites = ""

    # Loop through the enzymes list
    for enzyme in enzymes:
        # Get the site and its reverse complement from the restriction_enzymes dictionary
        site = restriction_enzymes[enzyme].site
        rev_comp = str(Seq(site).reverse_complement())

        # Add the site and its reverse complement to the string, separated by |
        restriction_enzyme_sites += "|" + site + "|" + rev_comp

        # If 6+ cutter, also add the site and its reverse complement without the first and last base to an edge search string:
        if len(site) == 6:
            restriction_edge_sites += "|" + site[:-1] + "|" + rev_comp[:-1] + "|" + site[1:] + "|" + rev_comp[1:]

    search_string += restriction_enzyme_sites
    restriction_enzyme_sites = restriction_enzyme_sites[1:]
    restriction_edge_sites = restriction_edge_sites[1:]

    return search_string, restriction_enzyme_sites, restriction_edge_sites


# function from DNAchisel: https://github.com/Edinburgh-Genome-Foundry/DnaChisel/
def gc_content(sequence, window_size=None):
    """Compute global or local GC content.
    Parameters
    ----------
    sequence
      An ATGC DNA sequence (upper case!)
    window_size
      If provided, the local GC content for the different sliding windows of
      this size is returned, else the global GC content is returned.
    Returns
    --------
      A number between 0 and 1 indication the proportion
      of GC content. If window_size is provided, returns
      a list of len(sequence)-window_size values indicating
      the local GC contents (sliding-window method). The i-th value
      indicates the GC content in the window [i, i+window_size]
    """
    # The code is a little cryptic as it uses numpy array operations
    # but the speed gain is 300x compared with pure-python string operations

    arr = np.frombuffer((sequence + "").encode(), dtype="uint8")
    arr_GCs = (arr == 71) | (arr == 67)  # 67=C, 71=G

    if window_size is None:
        return 1.0 * arr_GCs.sum() / len(sequence)
    else:
        cs = np.cumsum(arr_GCs)
        a = cs[window_size - 1 :]
        b = np.hstack([[0], cs[:-window_size]])
        return 1.0 * (a - b) / window_size


def codon_optimize(
    sequence: str,
    enzymes: List,
    sites_to_check: dict,
    verbose: int = 0,
    if_fail_return_input_sequence: bool = False,
    if_fail_return_optimized_sequence: bool = False,
    if_fail_return_none: bool = False,
):
    input_sequence = sequence
    codon_dict = codon_1_to_codon_2
    seq_len = len(sequence)

    # print(f'Optimizing sequence... \n{sequence} {Seq(sequence).translate(), round(gc_content(sequence), 2)}')
    if seq_len % 3 != 0:
        raise Exception("Sequence must be a coding-sequence and have a multiple of 3")

    (search_string, restriction_enzyme_sites, restriction_edge_sites) = make_search_strings(enzymes=enzymes, sites_to_check=sites_to_check)

    # for i, sequence in enumerate(sequence_list):
    gc_lower = 0.35
    gc_higher = 0.6

    restriction_edge_bool = (search(restriction_edge_sites, sequence[:5]) is not None) | (search(restriction_edge_sites, sequence[-5:]) is not None)

    if (gc_content(sequence) >= gc_higher) | (gc_content(sequence) <= gc_lower) | (search(search_string, sequence) is not None) | restriction_edge_bool:
        aa_sequence = Seq(sequence).translate()
        sequence_codons = np.array([aa_to_codon[aa] for aa in aa_sequence])
        iterations = 0
        iteration_weights = [2, 4, 6, 8]
        i = 0
        while (gc_content(sequence) >= gc_higher) | (gc_content(sequence) <= gc_lower) | (search(search_string, sequence) is not None) | restriction_edge_bool:
            i += 1

            # less stringent settings:
            if iterations == seq_len * iteration_weights[0]:
                if verbose == 2:
                    print("Cannot optimize sequence in " + str(seq_len * iteration_weights[0]) + " iterations with less stringency:")
                    print(sequence, f"(aa {Seq(sequence).translate()})")

                    if gc_content(sequence) >= gc_higher or gc_content(sequence) <= gc_lower:
                        print(f"gc content outside constraint: ", round(gc_content(sequence), 2))

                    if search(search_string, sequence) is not None:
                        print(f"found search string: {search(search_string, sequence) is not None}")
                        print(f"{search(search_string, sequence)}")

                    if search(restriction_enzyme_sites, sequence) is not None:
                        print(
                            f"found restriction site(s) string: {search(restriction_enzyme_sites, sequence) is not None} (of enzymes: {enzymes}, sites: {restriction_enzyme_sites})"
                        )
                        print(f"{search(restriction_enzyme_sites, sequence)}")

                    if restriction_edge_bool:
                        print(f"found restriction site(s) edges: {restriction_edge_bool}")

                # Less stringent settings:
                gc_higher = 0.65
                gc_lower = 0.25

                # need to adapt the search string to accomodate for the alternative site generation
                # if search(search_string, sequence) is not None:
                #     # search_string = 'GGGGGG|CCCCCC|TTTTTT|AAAAAA|GGGGG|CCCCC|TTTTT|AAAAA|' + restriction_enzyme_sites[:-1]
                #     search_string = 'GGGGGG|CCCCCC|TTTTTT|AAAAAA|GGGGG|CCCCC|TTTTT|AAAAA' + '|' + restriction_enzyme_sites

            if iterations == seq_len * iteration_weights[1]:
                if verbose == 1:
                    print("Cannot optimize sequence in " + str(seq_len * iteration_weights[1]) + " iterations with less stringency:")
                    print(sequence, f"(aa {Seq(sequence).translate()})")
                    if gc_content(sequence) >= gc_higher or gc_content(sequence) <= gc_lower:
                        print(f"gc content outside constraint: ", round(gc_content(sequence), 2))

                    if search(search_string, sequence) is not None:
                        print(f"found search string: {search(search_string, sequence) is not None}")
                        print(f"{search(search_string, sequence)}")

                    if search(restriction_enzyme_sites, sequence) is not None:
                        print(
                            f"found restriction site(s) string: {search(restriction_enzyme_sites, sequence) is not None} (of enzymes: {enzymes}, sites: {restriction_enzyme_sites})"
                        )
                        print(f"{search(restriction_enzyme_sites, sequence)}")

                    if restriction_edge_bool:
                        print(f"found restriction site(s) edges: {restriction_edge_bool}")

                codon_dict = codon_2_to_codon_1

                if gc_content(sequence) >= gc_higher:
                    if seq_len <= 50:
                        if iterations <= seq_len * (iteration_weights[1] + 1):
                            gc_higher = 0.68
                        else:
                            gc_higher = 0.70
                    else:
                        if verbose == 1:
                            print("Input sequence length is longer than 50 bp so max_gc_content_after_iterations cannot be used!")

            if iterations == seq_len * iteration_weights[2]:
                codon_dict = codon_1_to_codon_2

            if iterations == seq_len * iteration_weights[3]:
                # if verbose:
                print("Cannot optimize sequence in " + str(seq_len * iteration_weights[3]) + " iterations with less stringency:")
                print(sequence, f"(aa {Seq(sequence).translate()})")

                if gc_content(sequence) >= gc_higher or gc_content(sequence) <= gc_lower:
                    print(f"gc content outside constraint: ", round(gc_content(sequence), 2))

                if search(search_string, sequence) is not None:
                    print(f"found search string: {search(search_string, sequence) is not None}")
                    print(f"{search(search_string, sequence)}")

                if search(restriction_enzyme_sites, sequence) is not None:
                    print(f"found restriction site(s) string: {search(restriction_enzyme_sites, sequence) is not None} (of enzymes: {enzymes}, sites: {restriction_enzyme_sites})")
                    print(f"{search(restriction_enzyme_sites, sequence)}")

                if restriction_edge_bool:
                    print(f"found restriction site(s) edges: {restriction_edge_bool}")

                if if_fail_return_input_sequence:
                    return input_sequence
                elif if_fail_return_optimized_sequence:
                    return sequence
                elif if_fail_return_none:
                    return None
                else:
                    raise Exception("Cannot optimize sequence in " + str(seq_len * 12) + " iterations with less stringency!")

            rng = default_rng(i)
            codon = rng.choice(sequence_codons.shape[0], size=1)[0]

            iterations += 1
            if sequence_codons[codon] in codon_dict.keys():
                if gc_content(sequence) >= gc_higher:
                    if gc_content(sequence_codons[codon]) >= 1:
                        sequence_codons[codon] = codon_dict[sequence_codons[codon]]
                        sequence = "".join(map(str, [res for res in sequence_codons]))

                elif gc_content(sequence) <= gc_lower:
                    if gc_content(sequence_codons[codon]) <= gc_lower:
                        sequence_codons[codon] = codon_dict[sequence_codons[codon]]
                        sequence = "".join(map(str, [res for res in sequence_codons]))

                elif search(search_string, sequence) is not None:
                    sequence_codons[codon] = codon_dict[sequence_codons[codon]]
                    sequence = "".join(map(str, [res for res in sequence_codons]))

                elif restriction_edge_bool:
                    sequence_codons[codon] = codon_dict[sequence_codons[codon]]
                    sequence = "".join(map(str, [res for res in sequence_codons]))

            restriction_edge_bool = (search(restriction_edge_sites, sequence[:5]) is not None) | (search(restriction_edge_sites, sequence[-5:]) is not None)

        # If this exception is ever raised, implement fix ...
        if search("TTTTT|AAAAA", sequence):
            raise Exception("AT homo-polymer found")

        # Double check for restriction enzyme sites with an external software:
        elif search(search_string, sequence) is not None:
            raise Exception("Enzyme sites still found")

        elif aa_sequence != (Seq(sequence).translate()):
            raise Exception("Codon optmized DNA sequence does not encode for the same AA sequence as the original DNA sequence.")

        else:
            sequence = sequence

    # print('Codon-optimized:', k, 'sequences of', str(i + 1), 'sequences\n')
    # print('Optimized  sequence...  ', sequence, '###################')
    return sequence
