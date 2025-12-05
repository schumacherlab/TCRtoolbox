aa_to_codon = {
    "*": "TAA",
    # changed to TAA because TAA & TAG stop codons ~2.5 less read-through than TGA https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5829715/
    "A": "GCC",  # GCG other website
    "C": "TGC",
    "D": "GAC",  # GAT other website
    "E": "GAG",
    "F": "TTC",
    "G": "GGC",
    "H": "CAC",
    "I": "ATC",
    "K": "AAG",
    "L": "CTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCT",
    # only does not agree on this codon changed from CCC to CCT source: https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs
    "Q": "CAG",
    "R": "AGA",
    "S": "AGC",
    "T": "ACC",
    "V": "GTG",
    "W": "TGG",
    "Y": "TAC",
}

aa_to_codon_2 = {
    "*": "TAG",
    # changed to TAA because TAA & TAG stop codons ~2.5 less read-through than TGA https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5829715/
    "A": "GCT",
    "C": "TGT",
    "D": "GAT",
    "E": "GAA",
    "F": "TTC",  # Kept TTC to not introduce terminator sequences
    "G": "GGA",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",  # Kept AAG to not introduce terminator sequences
    "L": "CTT",
    # Changed from CTC to CTT to not introduce hard to fix BsmBI sites. 17.81 to 14.08 in https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs
    "M": "ATG",
    "N": "AAT",
    "P": "CCT",
    "Q": "CAA",
    "R": "AGG",
    "S": "TCC",
    "T": "ACA",
    "V": "GTT",  # Changed from GTC to GTT to not introduce hard to fix BsmBI sites.
    "W": "TGG",
    "Y": "TAT",
}

codon_1_to_codon_2 = {aa_to_codon[key]: value for key, value in aa_to_codon_2.items()}
codon_2_to_codon_1 = {aa_to_codon_2[key]: value for key, value in aa_to_codon.items()}
codon_to_aa_2 = {v: k for k, v in aa_to_codon_2.items()}
codon_to_aa = {v: k for k, v in aa_to_codon.items()}

# Codons are ordered in these lists on their usage frequency
# using https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs#cptab:
aa_to_all_codons_list = {
    "G": ["GGC", "GGA", "GGG", "GGT"],
    "W": ["TGG"],
    "C": ["TGC", "TGT"],
    "F": ["TTC", "TTT"],
    "L": ["CTG", "CTC", "CTT", "TTG", "TTA", "CTA"],
    "I": ["ATC", "ATT", "ATA"],
    "V": ["GTG", "GTC", "GTT", "GTA"],
    "S": ["AGC", "TCC", "TCT", "TCA", "AGT", "TCG"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "T": ["ACC", "ACA", "ACT", "ACG"],
    "A": ["GCC", "GCT", "GCA", "GCG"],
    "Y": ["TAC", "TAT"],
    "H": ["CAC", "CAT"],
    "Q": ["CAG", "CAA"],
    "N": ["AAT", "AAC"],
    "K": ["AAG", "AAA"],
    "D": ["GAC", "GAT"],
    "E": ["GAG", "GAA"],
    "R": ["AGA", "AGG", "CGG", "CGC", "CGA", "CGT"],
    "M": ["ATG"],
    # Stop codons
    "*": ["TAA", "TAG", "TGA"],
}

AA_STRING = "ACDEFGHIKLMNPQRSTVWY"

AA_TO_ID = {
    "A": 0,
    "C": 1,
    "D": 2,
    "E": 3,
    "F": 4,
    "G": 5,
    "H": 6,
    "I": 7,
    "K": 8,
    "L": 9,
    "M": 10,
    "N": 11,
    "P": 12,
    "Q": 13,
    "R": 14,
    "S": 15,
    "T": 16,
    "V": 17,
    "W": 18,
    "Y": 19,
}

ID_TO_AA = {v: k for k, v in AA_TO_ID.items()}

p20_cloning_sites = {"end_5_BsmBI_p20": "CGATCGTCTCACAGC", "end_3_BsmBI_p20": "CTGATGAGACGGTAG"}

p20_pMX_trim_seqs = {
    "pcr_pMX_rev_muTCRb": "GGGACACTTTTGGAGGTGTGACGTTCCGCAGGTCCTC",
    "pcr_pMX_rev_muTCRa": "GGGGGTCCTTCAGCTGGTATACGGCGGGCTCGGGGTTCTGGATGT",
    #'pcr_pMX_fwd':'GGGGTGGACCATCCTCTAGACTGCCGGATCCGGCGGCCGCC',
    # TGA stop codon must be trimmed to prevent alignment artifacts
    # (probably only a problem w very similar epitope seqs w/o dissimilar barcode seq in between)
    "pcr_p20_TGA": "TCAGTGTTTTTGGAGGTGTGACGTGCTCAGTCA",
    "pcr_p20_TAA": "TCAGTGTTTTTGGAGGTGTGACGTGCTCAGTTA",
    "pcr_p20_barcode": "TCAGTGTTTTTGGAGGTGTGACGTGCTCAG",
    "RT_TCRb_1": "TTTTGGAGGTGTGACGTTCCGCAGGTCCTC",
    # TAA stop codon must be trimmed to prevent alignment artifacts:
    # (probably only a problem w very similar epitope seqs w/o dissimilar barcode seq in between),
    "RT_TCRb_1_p20_barcode": "TTTTGGAGGTGTGACGTGCTCAG",
    "RT_TCRb_1_p20_TGA": "TTTTGGAGGTGTGACGTGCTCAGTCA",
    "RT_TCRb_1_p20_TAA": "TTTTGGAGGTGTGACGTGCTCAGTTA",
    "p20_IL2_SP": "GCTGTTGGTCACCAGGGCCAGGCTCAGGGCGATGCAGCTCAGCAG",  # On purpose partial sequence of p20_IL2_SP. Use in cutadapt with -a and without $ anchor or X internal symbols
}

cigar_string_to_pysam_cigar_number = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8, "B": 9}

pysam_cigar_number_to_cigar_string = {v: k for k, v in cigar_string_to_pysam_cigar_number.items()}


positive_epitopes_list = [
    "ELAGIGILTV",  # MART1-ELA
    "EAAGIGILTV",  # MART1-WT
    "SLLMWITQA",  # NYESO1
    "SLLMWITQC",  # NYESO1 native peptide that should be less well recognized
    "NLVPMVATV",  # NLV, order new C7 TCR that does not have the RT-primer mis-priming site in V gene
    "ARDPHSGHFV",  # CDK4 WT
    "ALDPHSGHFV",  # CDK4mut
    "ACDPHSGHFV",  # CDK4mut2 that should bind HLA-A2 less well
]

known_negative_epitopes_list = [
    "NLNCCSVPV",  # GNL3L
    "LLYDANYFL",  # LLY
    "KLWAQCVQL",  # KLW
    "LLFGYPVYV",  # Tax, used as negative because model TCR has a bad reactivity profile
]

each_codon_to_aa = {
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    "TGG": "W",
    "TGT": "C",
    "TGC": "C",
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "AGT": "S",
    "AGC": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",
    "ATG": "M",
    # Stop codons
    "TAG": "*",
    "TAA": "*",
    "TGA": "*",
}

# https://www.genscript.com/tools/codon-frequency-table
codon_usage = {
    "TTT": 0.45,
    "TTC": 0.55,
    "TTA": 0.07,
    "TTG": 0.13,
    "TAT": 0.43,
    "TAC": 0.57,
    "TAA": 0.28,
    "TAG": 0.2,
    "CTT": 0.13,
    "CTC": 0.2,
    "CTA": 0.07,
    "CTG": 0.41,
    "CAT": 0.41,
    "CAC": 0.59,
    "CAA": 0.25,
    "CAG": 0.75,
    "ATT": 0.36,
    "ATC": 0.48,
    "ATA": 0.16,
    "ATG": 1.0,
    "AAT": 0.46,
    "AAC": 0.54,
    "AAA": 0.42,
    "AAG": 0.58,
    "GTT": 0.18,
    "GTC": 0.24,
    "GTA": 0.11,
    "GTG": 0.47,
    "GAT": 0.46,
    "GAC": 0.54,
    "GAA": 0.42,
    "GAG": 0.58,
    "TCT": 0.18,
    "TCC": 0.22,
    "TCA": 0.15,
    "TCG": 0.06,
    "TGT": 0.45,
    "TGC": 0.55,
    "TGA": 0.52,
    "TGG": 1.0,
    "CCT": 0.28,
    "CCC": 0.33,
    "CCA": 0.27,
    "CCG": 0.11,
    "CGT": 0.08,
    "CGC": 0.19,
    "CGA": 0.11,
    "CGG": 0.21,
    "ACT": 0.24,
    "ACC": 0.36,
    "ACA": 0.28,
    "ACG": 0.12,
    "AGT": 0.15,
    "AGC": 0.24,
    "AGA": 0.2,
    "AGG": 0.2,
    "GCT": 0.26,
    "GCC": 0.4,
    "GCA": 0.23,
    "GCG": 0.11,
    "GGT": 0.16,
    "GGC": 0.34,
    "GGA": 0.25,
    "GGG": 0.25,
}

VDJ_CDR3_COMMON_END_MOTIFS = ["FGXG", "WGXG"]
VDJ_CDR3_RARE_END_MOTIFS = ["XGXG", "FXXG"]
VDJ_CDR3_ALL_END_MOTIFS = VDJ_CDR3_COMMON_END_MOTIFS + VDJ_CDR3_RARE_END_MOTIFS