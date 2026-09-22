model_minigene_aa_dict = {"MEL063_CCSER2_P329L": "GYRMVHPSLLKSSRSLFSGTMTVDGNKNSPA", "MEL063_TNFAIP2_P348A": "DRALELEARRWAEDVAPQRLDGHCHSELAID"}

# Before I removed the 4-bases at the end that were added to increase distance to DNA end but not needed because of orthoprimers, end sequences were:
# end_5_BsmBI_p20 : CGATCGTCTCACAGC and end_3_BsmBI_p20 : CTGATGAGACGGTAG
# Ends consist of BsmBI recognition site + 1 base + 4-base overhang
p20_cloning_sites = {"end_5_BsmBI_p20": "CGTCTCACAGC", "end_3_BsmBI_p20": "CTGATGAGACG"}

p12_cloning_sites = {"end_5_BsmBI_p12": "CGTCTCTGAAG", "end_3_BsmBI_p12": "CTGATGAGACG"}


p12_p20_flanking_nt_dict = {
    "p12_end_5": "AGCTGGAGAACCTGAGAATGAAG",
    "p12_end_3": "CTGAGCACGTCACACCTCCAAAA",
    "p20_end_5": "GCCTGGCCCTGGTGACCAACAGC",
    "p20_end_3": "CTGAGCACGTCACACCTCCAAAA",
}
