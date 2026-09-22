THREE_PRIME_LINKED_ADAPTERS = ["p20_IL2_SP", "p12_CD74"]


def linked_adapter_pairs(five_prime_keys: list[str]) -> list[tuple[str, str]]:
    """Pair each 5' anchor key with every 3' adapter it should be tried against."""
    return [(five, three) for five in five_prime_keys for three in THREE_PRIME_LINKED_ADAPTERS]


trim_adapter_config = {
    "no_epi_barcode_rt_pair_scan_plate": {"-g": ["RT_TCRb_1"], "-a_linked": linked_adapter_pairs(["RT_TCRb_1_p20_TGA", "RT_TCRb_1_p20_TAA"])},
    "epi_barcode_rt_pair_scan_plate": {"-g": ["RT_TCRb_1"], "-a_linked": linked_adapter_pairs(["RT_TCRb_1_p20_barcode"])},
    "no_epi_barcode_pcr_gdna": {"-g": ["pcr_pMX_rev_muTCRb", "pcr_pMX_rev_muTCRa"], "-a_linked": linked_adapter_pairs(["pcr_p20_TGA", "pcr_p20_TAA"])},
    "epi_barcode_pcr_gdna": {"-g": ["pcr_pMX_rev_muTCRb", "pcr_pMX_rev_muTCRa"], "-a_linked": linked_adapter_pairs(["pcr_p20_barcode"])},
}
