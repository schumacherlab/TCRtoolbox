from typing import List, Literal, Dict
import numpy as np
from Bio.Restriction import Restriction
from numpy.random import default_rng
from re import search
from Bio.Seq import Seq
from tcr_toolbox.utils.constants import aa_to_all_codons_list, codon_usage
from dotenv import load_dotenv
import os
from Bio import SeqIO
from Bio.SeqUtils import CodonAdaptationIndex

# dotenv
load_dotenv()
tcr_toolbox_data_path = os.getenv("tcr_toolbox_data_path")

# CAI
reference_sequences = []
reference_fasta = (
    tcr_toolbox_data_path
    + "/tcr_toolbox_datasets/tcr_assembly/ensembl_human_genome_fasta/reference_highly_expressed.fa"
)
# Load sequences into a list
for record in SeqIO.parse(reference_fasta, "fasta"):
    seq = str(record.seq)
    if len(seq) % 3 == 0:
        if "N" in seq:
            continue
        reference_sequences.append(seq)
cai = CodonAdaptationIndex(sequences=reference_sequences)


def make_search_strings(enzymes: List[str] | None = None, sites_to_check: Dict[str, str] | None = None):
    search_string = ""
    restriction_enzyme_sites = ""
    restriction_edge_sites = ""

    if sites_to_check:
        sites = list(sites_to_check.values())
        for site in sites:
            rc = str(Seq(site).reverse_complement())
            if not search_string:
                search_string = site + "|" + rc
            else:
                search_string += "|" + site + "|" + rc

    # enzyme sites
    if enzymes:
        enzyme_terms = []
        edge_terms = []
        for name in enzymes:
            enzyme = getattr(Restriction, name)
            site = enzyme.site
            rc = str(Seq(site).reverse_complement())

            enzyme_terms.append(site)
            enzyme_terms.append(rc)

            if len(site) == 6:
                edge_terms.extend([site[:-1], rc[:-1], site[1:], rc[1:]])

        # add enzyme sites to search_string (same rule as original)
        if enzyme_terms:
            if search_string:
                search_string += "|" + "|".join(enzyme_terms)
            else:
                search_string = "|".join(enzyme_terms)

        restriction_enzyme_sites = "|".join(enzyme_terms)
        restriction_edge_sites = "|".join(edge_terms)

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


def print_error_message_optimize(sequence, search_string, restriction_enzyme_sites, restriction_edge_bool, gc_higher, gc_lower, iterations, enzymes, cai_lower):
    print("Cannot optimize sequence in " + str(iterations) + " iterations with less stringency:")
    print(sequence, f"(aa {Seq(sequence).translate()})")

    if gc_content(sequence) >= gc_higher or gc_content(sequence) <= gc_lower:
        print(f"gc content outside constraint: {round(gc_content(sequence), 2)}")

    if search(search_string, sequence) is not None:
        print(f"found search string: {search(search_string, sequence) is not None}")
        print(f"{search(search_string, sequence)}")

    if search(restriction_enzyme_sites, sequence) is not None:
        print(f"found restriction site(s) string: {search(restriction_enzyme_sites, sequence) is not None} (of enzymes: {enzymes}, sites: {restriction_enzyme_sites})")
        print(f"{search(restriction_enzyme_sites, sequence)}")

    if restriction_edge_bool:
        print(f"found restriction site(s) edges: {restriction_edge_bool}")
    if cai.calculate(sequence) < cai_lower:
        print(f"cai lower than constraint: {round(cai.calculate(sequence), 2)}")


def make_codon_sampling_probability_distribution(n, steepness):
    """
    Generate n values following a steep decay pattern with adjustable steepness.
    Args:
        n (int): Number of values to generate.
        steepness (float): Factor to adjust the steepness of the decay. Higher values increase steepness.

    Returns:
        list: Decay values normalized to sum to 1.

    Examples:
        n=3, steepness=2: [0.86681333 0.11731043 0.01587624]
        n=3, steepness=1: [0.66524096 0.24472847 0.09003057]
        n=3, steepness=0.5: [0.50648039 0.30719589 0.18632372]

        n=2, steepness=2: [0.88079708 0.11920292]
        n=2, steepness=1: [0.73105858 0.26894142]
        n=2, steepness=0.5: [0.62245933 0.37754067]
    """
    x = np.arange(1, n + 1)
    values = np.exp(-steepness * (x - 1))  # Exponential decay with adjustable steepness
    normalized_values = values / values.sum()  # Normalize to sum to 1

    return normalized_values


def make_codon_sampling_probability_distribution_with_prior(codon_list, steepness):
    if steepness < 1:
        raise ValueError("Steepness must be greater than or equal to 1")

    codon_freq_list = np.array([codon_usage[codon] ** steepness for codon in codon_list])
    codon_freq_list = codon_freq_list / np.sum(codon_freq_list)
    return codon_freq_list


def codon_optimize(
    sequence_aa: str,  # AA sequence
    enzymes: List,
    sites_to_check: dict,
    verbose: int = 0,  # level 0 (none), 1 (low), 2 (medium), 3 (high), 4 (very high)
    if_fail: Literal["raise", "return_optimized_sequence", "return_none"] = "raise",  # removed "return_input_sequence", as you dont want to return AA when you expect NT
    steepness: float = 4,
    idx: int = 0,  # used to determine the seed for the random number generator (optional)
):
    """
    Optimize a protein (amino acid) sequence into a codon-optimized nucleotide DNA sequence
    for expression in Homo sapiens or other organisms, while respecting GC content,
    codon usage bias, restriction enzyme avoidance, and codon adaptation index (CAI).

    The algorithm iteratively selects codons for each amino acid, biased by human
    codon usage frequencies. During optimization, constraints such as GC content,
    CAI, and absence of restriction sites are enforced. If constraints are not met,
    the algorithm gradually relaxes the codon bias (via `steepness`) to attempt a solution.

    Parameters
    ----------
    sequence_aa : str
        Amino acid sequence to be codon-optimized. Must only contain valid amino acid
        single-letter codes. DNA sequences (ATGC) are not allowed.
    enzymes : list
        List of restriction enzyme objects or names to check against in the sequence.
    sites_to_check : dict
        Optional dictionary mapping enzyme names to recognition sequences for additional checks.
    verbose : int, default=0
        Level of logging during optimization:
            0 = none, 1 = low, 2 = medium, 3 = high, 4 = very high.
    if_fail : Literal["raise", "return_optimized_sequence", "return_none"], default="raise"
        Behavior if sequence cannot be optimized within iteration limits:
            - "raise": raise an Exception
            - "return_optimized_sequence": return the best sequence achieved
            - "return_none": return None
    steepness : float, default=4
        Controls codon bias during optimization:
            - Higher values strongly favor frequent optimal codons.
            - Lower values (>=1) give more uniform codon selection.
            - Must be >=1 for `make_codon_sampling_probability_distribution_with_prior`.
            - Can be reduced during optimization if constraints are not met.
    idx : int, default=0
        Optional seed for random number generator to allow reproducible results.

    Returns
    -------
    sequence_nt : str
        Codon-optimized DNA nucleotide sequence encoding the same protein as `sequence_aa`.

    Raises
    ------
    Exception
        If input sequence is invalid, contains non-amino acid characters, or if
        codon optimization fails and `if_fail="raise"`.

    Notes
    -----
    - Optimization ensures:
        - GC content remains between configurable lower and upper bounds.
        - Avoidance of specified restriction enzyme recognition sites.
        - Codon Adaptation Index (CAI) above the threshold.
    - If constraints are not met after multiple iterations, the algorithm relaxes
    constraints (reducing `steepness`, relaxing GC/CAI bounds) progressively.
    - The function uses a randomized iterative strategy, so results may vary unless `idx` is fixed.
    - Final sequence is verified to encode the same amino acid sequence as input.
    - The algorithm is compatible with `aa_to_all_codons_list` and `codon_usage` dictionaries.

    Examples
    --------
    # Codon-optimize a small protein sequence
    optimized_seq = codon_optimize(
        sequence_aa="MKT",
        enzymes=[EcoRI, BamHI],
        sites_to_check={"p5_illumina_seq_primer_1": "GAATTCAGCCTGAG"},
        verbose=2,
        steepness=4,
        idx=42
    )
    """
    if not set(sequence_aa).issubset(set(aa_to_all_codons_list.keys())) or set(sequence_aa) == set("ATGC"):
        raise Exception("Input sequence must be a protein sequence")

    (search_string, restriction_enzyme_sites, restriction_edge_sites) = make_search_strings(enzymes=enzymes, sites_to_check=sites_to_check)

    gc_lower = 0.35
    gc_higher = 0.6
    cai_lower = 0.95

    rng = default_rng(idx)
    sequence_codons = []

    for aa in sequence_aa:
        codon_options = aa_to_all_codons_list[aa]
        probabilities = make_codon_sampling_probability_distribution_with_prior(codon_options, steepness=steepness)
        chosen_codon = rng.choice(codon_options, p=probabilities)
        sequence_codons.append(chosen_codon)

    sequence_codons = np.array(sequence_codons)
    sequence_nt = "".join(sequence_codons)
    seq_len = len(sequence_nt)
    iterations = 0  # used to determine how many iterations have been done and rng seed (1-indexed, so 1 is the first iteration)
    iteration_weights = [2, 4, 6, 8]  # currently we hard coded (number_iterations == seq_len x iteration_weights) to determine stringency
    verbose_threshold = [4, 3, 2, 1]  # verbose level to print error message (0 is the most stringent)
    curr_it_weight = 0  # used to determine which stringency level we are at
    curr_steepness = steepness  # used to determine how steep the probability distribution is
    restriction_edge_bool = (search(restriction_edge_sites, sequence_nt[:5]) is not None) | (search(restriction_edge_sites, sequence_nt[-5:]) is not None)

    while (
        gc_content(sequence_nt) >= gc_higher
        or gc_content(sequence_nt) <= gc_lower
        or search(search_string, sequence_nt) is not None
        or restriction_edge_bool
        or cai.calculate(sequence_nt) < cai_lower
    ):
        if iterations == seq_len * iteration_weights[curr_it_weight]:
            if verbose >= verbose_threshold[curr_it_weight]:
                print_error_message_optimize(
                    sequence=sequence_nt,
                    search_string=search_string,
                    restriction_enzyme_sites=restriction_enzyme_sites,
                    restriction_edge_bool=restriction_edge_bool,
                    gc_higher=gc_higher,
                    gc_lower=gc_lower,
                    iterations=iterations,
                    enzymes=enzymes,
                    cai_lower=cai_lower,
                )

            if curr_it_weight == 0:
                cai_lower = 0.90
                gc_higher = 0.65
                gc_lower = 0.25

            if curr_it_weight == 1:
                cai_lower = 0.85
                if curr_steepness - 1 >= 1:
                    curr_steepness -= 1

                if gc_content(sequence_nt) >= gc_higher:
                    if seq_len <= 50:
                        if iterations <= seq_len * (iteration_weights[curr_it_weight] + 1):
                            gc_higher = 0.68
                        else:
                            gc_higher = 0.70
                    else:
                        if verbose == 1:
                            print("Input sequence length is longer than 50 bp so max_gc_content_after_iterations cannot be used!")

            if curr_it_weight == 2:
                cai_lower = 0.75
                if curr_steepness - 1 >= 1:
                    curr_steepness -= 1

            if curr_it_weight == 3:
                if if_fail == "return_optimized_sequence":
                    return sequence_nt
                elif if_fail == "return_none":
                    return None
                else:
                    raise Exception(
                        f"Cannot optimize sequence in {str(iterations)} iterations with less stringency!" + " To identify why, set codon_optimize(verbose) to 1 or higher."
                    )
            curr_it_weight += 1

        rng = default_rng(iterations + idx)  # seed is iterations + idx to prevent the same starting seed for different sequences
        codon_i = rng.choice(sequence_codons.shape[0], size=1)[0]  # pick a random codon (index)
        current_codon = sequence_codons[codon_i]
        new_codon_options = aa_to_all_codons_list[sequence_aa[codon_i]]
        # probabilities = make_codon_sampling_probability_distribution(len(new_codon_options), steepness=curr_steepness) # without prior
        probabilities = make_codon_sampling_probability_distribution_with_prior(new_codon_options, steepness=curr_steepness)
        new_codon = rng.choice(new_codon_options, p=probabilities)

        if new_codon != current_codon and len(new_codon_options) > 1:
            if gc_content(sequence_nt) >= gc_higher:
                if gc_content(current_codon) >= gc_content(new_codon):
                    sequence_codons[codon_i] = new_codon

            elif gc_content(sequence_nt) <= gc_lower:
                if gc_content(current_codon) <= gc_content(new_codon):
                    sequence_codons[codon_i] = new_codon

            elif search(search_string, sequence_nt) is not None or restriction_edge_bool:
                sequence_codons[codon_i] = new_codon

            sequence_nt = "".join(map(str, [res for res in sequence_codons]))

        restriction_edge_bool = (search(restriction_edge_sites, sequence_nt[:5]) is not None) | (search(restriction_edge_sites, sequence_nt[-5:]) is not None)
        iterations += 1

    if search("TTTTT|AAAAA", sequence_nt):  # If this exception is ever raised, implement fix ... (should be impossible)
        raise Exception("AT homo-polymer found")

    elif search(search_string, sequence_nt) is not None:  # Double check for restriction enzyme sites with an external check:
        raise Exception("Enzyme sites still found")

    elif sequence_aa != (Seq(sequence_nt).translate()):
        raise Exception("Codon optmized DNA sequence does not encode for the same AA sequence as the original DNA sequence.")

    return sequence_nt
