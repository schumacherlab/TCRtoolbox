# TCR Assembly CLI Pipeline

This repository provides a command-line interface (CLI) for running the **TCR assembly pipeline**:
```
Input CSV (run_tcr_csv).
        │
Check whether input CSV contains required columns & does not contain NaNs.
        │
Translate V genes to IMGT format.
        │
Filter TCRs based on unique V gene names that are available in stock. 
        │
Generate CDR3-J nucleotide order sequences through codon optimization to
introduce sufficient diversity for distinguishing all TCRs using only
CDR3-Jβ Illumina sequencing.
        │
Assign TCRs to wells in plates (384-well)
        │
Add well-specific orthoprimer combinations to CDR3-J order nt sequences.
        │
Generate CDR3-J oligo pool order sheet. 
        │
Generate V Gene Premix Dispense (idot / echo / Hamilton) instruction files. 
        │
Perform in silico PCR & Golden Gate assembly pydna simulation on written
assembly files (oligo pool order sheet and V gene Premix dispense files). 
        │
Check whether assembly will be successfull by comparing translation of
simulated ligation products to TCR aa sequences reconstructed from the
original input CSV. 
        │
Optional: Generate sequencing reference files (Illumina & Nanopore): alpha, beta,
full-length (nanopore), and an optional combined beta + epitope nt sequences
reference. 
```

The pipeline consists of **three sequential run modes**:

1. **simulation mode** : Check that your TCR input CVS passes all asssembly requirements by simulating the full pipeline.
   Re-design input until you pass all checks. 
2. **oligo_order mode** : Prepare CDR3-J oligo pool synthesis order sheets. Optional: check whether you have enough V gene
   stock tube volume available and identify which V gene stock tubes need replacement. 
4. **run_assembly mode** : Generate robotic dispense instructions files for final TCR assembly wet-lab execution. Optional:
   raise an error if any required V gene stock tubes have not been replaced and update the V gene stock volume database
   based on the latest wet-lab assembly.

*Both `oligo_order` and `run_assembly` mode perform a `pydna` asssembly simulation to check that the assembly will be successful.
The oligo_order mode also generates V gene premix dispense instruction files for the simulation but deletes these afterward.*

> ⚠️ **Important:** Always run the modes in the order: `simulation` → `oligo_order` → `run_assembly`.

---
## Quick Start

Your TCR input CSV must include the following columns:

| Column Name      | Description |
|-----------------|-------------|
| **TRAV_IMGT**    | TRAV gene in IMGT format (10X Genomics format works). |
| **TRAJ_IMGT**    | TRAJ gene in IMGT format (10X Genomics format works). |
| **cdr3_alpha_aa** | CDR3 alpha amino acid sequence. |
| **TRBV_IMGT**    | TRBV gene in IMGT format (10X Genomics format works). |
| **TRBJ_IMGT**    | TRBJ gene in IMGT format (10X Genomics format works). |
| **cdr3_beta_aa** | CDR3 beta amino acid sequence. |
| **custom_name**  | Custom TCR identifier. For example, a TCR with custom name `ywe_151_293` assigned to plate 2, well K14 will have the standardized assembly name: <br>`1_2_K14_253_ywe_151_293 → [oligo subpool]_[plate number]_[well]_[well number]_[custom name]`. |
| **group**        | Assigns TCRs to plates in groups. <br> - Plates cannot contain multiple groups. <br> - If a group does not fully fill its last plate, it must occupy **at least half** of the plate. <br> **Example:** A patient `ywe` has 400 TCRs and a plate holds 384 wells: 1 full plate → 384 wells used, remaining → 16 wells. Half a plate = 192 wells → 16 wells is **not allowed**.|


Start a pipeline run mode as follows:

```bash
tcr_toolbox run-tcr-assembly configs/tcr_assembly/simulation_run_config.json
```

Configuration JSON template files for each run mode are located in:
```bash
TCRtoolbox/configs/
```

Make a copy of the configuration JSON template files and add a run_name suffix. 
For example, when you run_name is "hla_loss" run: 
```bash
cd insert/path/to/your/experiment

# Copy and rename simulation config
cp TCRtoolbox/configs/tcr_assembly/simulation_run_config.json simulation_run_config_hla_loss.json
# Edit simulation_run_config_hla_loss.json as needed and run:
tcr_toolbox run-tcr-assembly simulation_run_config_hla_loss.json

# Copy and rename oligo order config
cp TCRtoolbox/configs/tcr_assembly/oligo_order_run_config.json oligo_order_run_config_hla_loss.json
# Edit oligo_order_run_config_hla_loss.json as needed and run:
tcr_toolbox run-tcr-assembly oligo_order_run_config_hla_loss.json

# Copy and rename assembly config
# Edit run_assembly_run_config_hla_loss.json as needed and run:
tcr_toolbox run-tcr-assembly run_assembly_run_config_hla_loss.json
```

Below we describe the purpose of each run mode and how parameters in their onfiguration JSON template files should 
be adjusted:

## 1. Simulation Mode

**Overview:**  
Simulates the TCR assembly pipeline to check that your input CSV meets all requirements.  
Re-design CSV and simulate again until your CSV passes all checks. The simulation:
- Verifies CDR3-J sequences can be codon optimized and ligated to V genes.  
- Confirms V gene names are in stock and flags low-volume barcoded V gene tubes.  
- Ensures the last plate per TCR group is at least half filled.  
- Writes logs in the directory in which you run the CLI and deletes the simulation run directory afterward.  
- **Only checks; no real oligos or plates are produced for wet-lab execution.**

**CLI & Template:**  
```bash
tcr_toolbox run-tcr-assembly configs/tcr_assembly/simulation_run_config.json
```
- Normally, only edit in `simulation_run_config.json` template:
  - `run_tcr_csv` -> path to your TCR input CSV.
  - `run_name` -> custom run name (replace spaces in name with `_` and avoid special symbols). 
- Other parameters: only change after consulting configuration manual at the end of this tutorial. 

`TCRtoolbox/configs/tcr_assembly/simulation_run_config.json`: 
```bash
{
  "run_mode": "simulation",
  "run_tcr_csv": "insert/path/to/your/run_tcr_df.csv",
  "run_name": "[insert run name here (replace whitespaces with '_' and no '*/@-./&%#!' symbols]",
  "run_path": null,
  "grouping_col": "group",
  "oligo_name_annotation_col_list": null,
  "filter_cannot_be_codon_optimized": true,
  "filter_succeeding_nt_cys_104_beta": false,
  "allow_cdr3j_nt_duplicates": false,
  "steepness": 4,
  "well_plate_size": 384,
  "remove_longer_than_200_nt_oligos": false,
  "v_gene_premix_dispenser": "idot",
  "v_gene_transfer_volume_nl": 150,
  "max_v_gene_source_well_volume_nl": 70000,
  "v_gene_hamilton_ul": 21.12,
  "water_hamilton_ul": 58.88,
  "skip_barcoded_v_gene_stock_tube_vol_update": false, 
  "min_v_gene_stock_tube_vol_ul": 50.00, 
  "generate_illumina_refs": false,
  "epitope_barcode_refs": null,
  "epitope_name_col_name": "name",
  "epitope_order_col_name": "sequence",
  "trimmed_beta_model_tcr_dict": null,
  "tcr_read_length": 150,
  "model_epitope_dict": null,
  "epitope_barcode_length": 18,
  "generate_nanopore_ref": false,
  "add_number_of_negatives_by_mutating_full_refs": null,
  "add_number_of_negatives_by_mutating_cdr3j_in_refs": null,
  "add_number_of_negatives_by_mutating_v_gene_in_refs": null,
  "min_nt_diff_negative_ref_seqs": null,
  "max_nt_diff_negative_ref_seqs": null
}
```

## 2. oligo_order mode 

**Overview:**  
Prepares oligo pools for synthesis and checks that enough V gene stock volume is available.  
This mode must be run **before** `run_assembly` to perform wet-lab assembly. The process:
- Initializes a standardized assembly run directory. Run directories are currently written to the accompanying *dataset* path
'tcr_toolbox_data/tcr_toolbox_tcr_assembly_runs'.
- Assigns TCRs to wells in plates based on the TCR grouping column.  
- Prepares CDR3-J oligo pool order sheets for synthesis.  
- Flags barcoded V gene stock tubes with insufficient volume for replacement during oligo pool manufacturing.  
- Simulates in silico whether assembly yields correct ligation products using `pydna`.

**CLI & Template:**  
```bash
tcr_toolbox run-tcr-assembly configs/tcr_assembly/oligo_order_run_config.json
```
- Normally, only edit in `oligo_order_run_config.json` template:
  - `run_tcr_csv` -> path to your TCR input CSV.
  - `run_name` -> custom run name (replace spaces in name with `_` and avoid special symbols). 
- Other parameters: only change after consulting configuration manual at the end of this tutorial. 

`TCRtoolbox/tcr_toolbox/configs/tcr_assembly/oligo_order_run_config.json`: 
```bash
{
  "run_mode": "oligo_order",
  "run_tcr_csv": "insert/path/to/your/run_tcr_df.csv",
  "run_name": "[insert run name here (replace whitespaces with '_' and no '*/@-./&%#!' symbols]",
  "run_path": null,
  "grouping_col": "group",
  "oligo_name_annotation_col_list": null,
  "filter_cannot_be_codon_optimized": true,
  "filter_succeeding_nt_cys_104_beta": false,
  "allow_cdr3j_nt_duplicates": false,
  "steepness": 4,
  "well_plate_size": 384,
  "remove_longer_than_200_nt_oligos": false,
  "v_gene_premix_dispenser": "idot",
  "v_gene_transfer_volume_nl": 150,
  "max_v_gene_source_well_volume_nl": 70000,
  "v_gene_hamilton_ul": 21.12,
  "water_hamilton_ul": 58.88,
  "skip_barcoded_v_gene_stock_tube_vol_update": false, 
  "min_v_gene_stock_tube_vol_ul": 50.00, 
  "generate_illumina_refs": false,
  "epitope_barcode_refs": null,
  "epitope_name_col_name": "name",
  "epitope_order_col_name": "sequence",
  "trimmed_beta_model_tcr_dict": null,
  "tcr_read_length": 150,
  "model_epitope_dict": null,
  "epitope_barcode_length": 18,
  "generate_nanopore_ref": false,
  "add_number_of_negatives_by_mutating_full_refs": null,
  "add_number_of_negatives_by_mutating_cdr3j_in_refs": null,
  "add_number_of_negatives_by_mutating_v_gene_in_refs": null,
  "min_nt_diff_negative_ref_seqs": null,
  "max_nt_diff_negative_ref_seqs": null
}
```
## 3. run_assembly mode

**Overview:**  
Prepare the robotic assembly steps in wet-lab. 
This mode must be run **after** `oligo_order` mode, once:  
- The CDR3-J oligo pool has been manufactured.  
- Any barcoded V gene stock tubes requiring replacement have been replaced.  

The pipeline:
- Writes V gene premixing instruction files for robotic assembly. 
- Stores all instructions in the standardized assembly run directory initialized in `oligo_order` mode.
- Re-runs the pydna simulation again to check ligation products.  

**CLI & Template:**  
```bash
tcr_toolbox run-tcr-assembly configs/tcr_assembly/run_assembly_run_config.json
```
- Normally, only edit in the `oligo_order_run_config.json` template:  
  - `run_tcr_csv` → path to your TCR input CSV. Must match the CSV used in the corresponding `oligo_order` run.  
  - `run_name` → run name from the corresponding `oligo_order` run. Copy from the `oligo_order` config file.  
  - `run_path` → path to the standardized run directory initialized by the corresponding `oligo_order` run. Run directories are currently written to `tcr_toolbox_data/tcr_toolbox_tcr_assembly_runs`.  
  - Optional: `generate_illumina_refs = True` to generate Illumina reference files. To generate a combined TCR beta + epitope nucleotide reference file for PAIR-scan, include `epitope_barcode_refs` CSV or FASTA file.  
  - Optional: `generate_nanopore_ref = True` to generate Nanopore reference files.
- Other parameters: only change after consulting configuration manual at the end of this tutorial. 

`TCRtoolbox/tcr_toolbox/configs/tcr_assembly/run_assembly_run_config.json`:
```bash
{
  "run_mode": "run_assembly",
  "run_tcr_csv": "insert/path/to/your/run_tcr_df.csv",
  "run_name": "[insert run name that matches numbered_run_name in run_path in the line below this line]",
  "run_path": "insert/path/to/your/numbered_run_name_directory",
  "grouping_col": "group",
  "oligo_name_annotation_col_list": null,
  "filter_cannot_be_codon_optimized": true,
  "filter_succeeding_nt_cys_104_beta": false,
  "allow_cdr3j_nt_duplicates": false,
  "steepness": 4,
  "well_plate_size": 384,
  "remove_longer_than_200_nt_oligos": false,
  "v_gene_premix_dispenser": "idot",
  "v_gene_transfer_volume_nl": 150,
  "max_v_gene_source_well_volume_nl": 70000,
  "v_gene_hamilton_ul": 21.12,
  "water_hamilton_ul": 58.88,
  "skip_barcoded_v_gene_stock_tube_vol_update": false, 
  "min_v_gene_stock_tube_vol_ul": 50.00, 
  "generate_illumina_refs": true,
  "epitope_barcode_refs": null,
  "epitope_name_col_name": "name",
  "epitope_order_col_name": "sequence",
  "trimmed_beta_model_tcr_dict": null,
  "tcr_read_length": 150,
  "model_epitope_dict": null,
  "epitope_barcode_length": 18,
  "generate_nanopore_ref": true,
  "add_number_of_negatives_by_mutating_full_refs": null,
  "add_number_of_negatives_by_mutating_cdr3j_in_refs": null,
  "add_number_of_negatives_by_mutating_v_gene_in_refs": null,
  "min_nt_diff_negative_ref_seqs": null,
  "max_nt_diff_negative_ref_seqs": null
}

```

## Standardized run directory structure (in the toolbox dataset)
```bash
standardized_run_directory/
├── cdr3j_oligo_order_sheets/       # Oligo pool order sheets
├── plate_sheets/                   # .xlsx files per plate
│   └── (TCR names, V genes, CDR3 alpha/beta sequences, J genes, CDR3-J order sequences, well assignments)
├── v_genes_premix_dispense/        # V gene premix dispense instruction files
│   ├── v_gene_plate_maps/          # .pdf showing which V genes are added to which wells
│   └── v_gene_stock/               # V gene usage counts and missing stock information
├── pydna_logs/                     # pydna simulation log files
├── pydna_ligation_dicts/           # All TCRs ligated in pMX_S1_Kana_2 vector
│                                   # Can generate SnapGene files via generate_snapgene_files
└── sequencing_quality_analysis/    # Directory for TCR assembly QC FASTQ files
    └── references/                 # Illumina and nanopore reference files
```
*Do not edit and re-save files manually.*


## Configuration parameter documentation: 
- **`run_mode`** : `str`  
Specifies the mode of the pipeline. Options: `simulation`, `oligo_order`, `run_assembly`. 

- **`run_tcr_csv`** : `str` or `os.PathLike` 
Path to TCR input CSV.
  
- **`run_name`** : `str` or `None` 
Run name for output files and directories.

- **`run_path`** : `str` or `os.PathLike`, optional  
Base directory for outputs. Required argument in `run_assembly` mode. Optional argument in `simulation and `oligo_order` mode. 

- **`grouping_col`** : `str`, default=`'group'`
Column in TCR input CSV used to group TCRs into plates.

- **`oligo_name_annotation_col_list`** : `list[str]` or `None`  
Columns in TCR input CSV used to annotate oligo names. If `None`, defaults to `['custom_name']`. 
For example, a TCR with custom name "ywe_151_293" assigned to plate 2, well K14 will have the standardized assembly name:
1_2_K14_253_ywe_151_293 → [oligo subpool]_[plate number]_[well]_[well number]_[custom name].

- **`filter_cannot_be_codon_optimized`** : `bool`, default=`True`  
Remove TCR sequences that cannot be codon optimized.

- **`filter_succeeding_nt_cys_104_beta`** : `bool`, default=`False`  
  Filter CDR3β sequences based on the expected second amino acid after conserved Cys104.
  
  - The first three nucleotides of the overhang come from the conserved Cys104 codon (TGT). 
  - The fourth nucleotide of the overhang is contributed by the final nucleotide of the TRBV codon immediately before the Cys104 codon. 
  - In human TRBV genes, this nucleotide is highly conserved: * G for most TRBV genes * A for TRBV20-1 and TRBV29-1 (based on sequence analysis of ~89,000 TCRs). 
  - Consequently, the second amino acid in the CDR3β sequence is typically: * Alanine (A) when the nucleotide is G * Serine (S) when the nucleotide is A. 
    Sequences not matching these patterns may represent sequencing artifacts. 
  
  Behavior: 
  - If True: sequences where the second amino acid does not provide the expected nucleotide for the Cys104 overhang are removed.
  - If False: sequences are retained, and the second amino acid is replaced by the V gene compatible amino acid to produce the correct 4-base overhang.

- **`allow_cdr3j_nt_duplicates`** : `bool`, default=`False` 
Whether to allow duplicate CDR3-J nucleotide sequences. Set to `False` if you want to distinguish all TCRs in your library by 
only sequencing the CDR3-J of a single chain (alpha or beta). CDR3-J nucleotide sequence duplicates are made unique by picking
unique codons with the codon optimizer (i.e., codon optimizer induced diversity). 

- **`steepness`** : `float`, default = 4
  Bias toward frequent codons during optimization:
  - Higher values (e.g., 4) favor optimal codons.  
  - Lower values (≥1) make codon choice more uniform and thereby diverse.   
  - Must be ≥1 in `make_codon_sampling_probability_distribution_with_prior`.

- **`well_plate_size`** : `int`, default=`384`  
Plate size.

- **`remove_longer_than_200_nt_oligos`** : `bool`, default=`False`
Remove oligos longer than 200 nucleotides.  
Twist oligo pools have pricing tiers based on oligo length and pool size.

- **`v_gene_premix_dispenser`** : `str`, default=`'idot'` 
Dispenser type for V gene premix: `'idot'` or `'echo'`.
`echo` has been successfully used in wet-lab but is not actively maintained.

- **`v_gene_transfer_volume_nl`** : `int`, default=`150`  
Volume (nL) of V gene to transfer to each assembly well.

- **`max_v_gene_source_well_volume_nl`** : `int`, default=`70000` 
Maximum allowed volume per V gene source well (nL).

- **`v_gene_hamilton_ul`** : `float`, default=`21.12`  
Hamilton V gene stock transfer volume to idot/echo source plate (µL).

- **`water_hamilton_ul`** : `float`, default=`58.88` 
Hamilton water volume for diluting transferred V gene stock (µL).

- **`skip_barcoded_v_gene_stock_tube_vol_update`** : `bool`, default=`False`  
Skip updating V gene stock tube volumes for barcoded tubes. Set to `True` 
if you do not want to use the barcoded V gene stock volume tracking system 
of this package. 

- **`min_v_gene_stock_tube_vol_ul`** : `float`, default=`50.0`  
Minimum volume to maintain in a barcoded V gene stock tube (µL).

- **`generate_illumina_refs`** : `bool`, default=`True`  
Generate Illumina reference FASTA files. Reference files are only generated 
in `run_assembly` mode. 

- **`epitope_barcode_refs`** : `str` or `os.PathLike` or `None`  
CSV or FASTA containing epitope barcode references. Include if you want to write a 
combined TCR beta + epitope nucleotide sequence reference for e.g., PAIR-scan.  

- **`epitope_name_col_name`** : `str`, default=`'name'` 
Column name for epitope identifiers in epitope CSV. 

- **`epitope_order_col_name`** : `str`, default=`'sequence'`  
Column name for epitope nucleotide sequences in epitope CSV. 

- **`trimmed_beta_model_tcr_dict`** : `dict` or `None`  
Optionally include trimmed beta chain nucleotides of common model TCRs in reference files
that are stored in this dictionary.

- **`tcr_read_length`** : `int`, default=`150`  
Illumina sequencing read length.

- **`model_epitope_dict`** : `dict` or `None` 
Optional include epitope minigenes of common model epitopes in reference files
that are stored in this dictionary. 

- **`epitope_barcode_length`** : `int`, default=`18`  
Length of epitope barcode sequences.

- **`generate_nanopore_ref`** : `bool`, default=`True` 
Generate Nanopore reference FASTA files. Reference files are only generated 
in `run_assembly` mode. 

- **Nanopore negative control reference parameters:**
  - `add_number_of_negatives_by_mutating_full_refs` : `int` or `None`  
  - `add_number_of_negatives_by_mutating_cdr3j_in_refs` : `int` or `None`  
  - `add_number_of_negatives_by_mutating_v_gene_in_refs` : `int` or `None`  
  - `min_nt_diff_negative_ref_seqs` : `int` or `None`  
  - `max_nt_diff_negative_ref_seqs` : `int` or `None`  

## Make snapgene files of specific TCR names for further investigation:
```python
from tcr_toolbox.tcr_assembly.order_automation import generate_snapgene_files
from dotenv import load_dotenv

load_dotenv()
tcr_toolbox_data_path = os.getenv('tcr_toolbox_data_path')

# Long format TCR ID name
os.mkdir(os.path.join(tcr_toolbox_data_path, "/tcr_toolbox_tcr_assembly_runs/r7_ywe_t500_nsclc57_no3lam397_1st_ylq_gen/tcr_snapgene_files"))
generate_snapgene_files(
    tcr_list = ["1_2_K14_253_ywe_151_293", "1_3_B16_39_ywe_324_623"] 
    output_path: os.path.join(tcr_toolbox_data_path, "/tcr_toolbox_tcr_assembly_runs/r7_ywe_t500_nsclc57_no3lam397_1st_ylq_gen/tcr_snapgene_files"),
    run_path: os.path.join(tcr_toolbox_data_path, "/tcr_toolbox_tcr_assembly_runs/r7_ywe_t500_nsclc57_no3lam397_1st_ylq_gen/tcr_snapgene_files")
)
# OR 
# Short format TCR ID name: 
os.mkdir(os.path.join(tcr_toolbox_data_path, "/tcr_toolbox_tcr_assembly_runs/r7_ywe_t500_nsclc57_no3lam397_1st_ylq_gen/tcr_snapgene_files"))
generate_snapgene_files(
    tcr_list = ["2_K14", "3_B16"] 
    output_path: os.path.join(tcr_toolbox_data_path, "/tcr_toolbox_tcr_assembly_runs/r7_ywe_t500_nsclc57_no3lam397_1st_ylq_gen/tcr_snapgene_files"),
    run_path: os.path.join(tcr_toolbox_data_path, "/tcr_toolbox_tcr_assembly_runs/r7_ywe_t500_nsclc57_no3lam397_1st_ylq_gen")
)
```