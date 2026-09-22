# PAIR-scan counting tutorial: 

## 0. Generate reference `.fa` files
If you set `generate_illumina_refs = True` in `tcr_toolbox run-tcr-assembly run_assembly_run_config.json` and included barcoded-epitope references, then your reference files can be found in the TCR assembly run directory of your library: `[tcr_toolbox_data_path]/tcr_toolbox_tcr_assembly_runs/[your_run]/sequencing_quality_analysis/references/[your_run]_150bp_beta.fa`. If not, use the function `tcr_toolbox.sequencing_analysis.reference.generate_assembly_nt_refs()` in your counts analysis `.ipynb` notebook to generate reference `.fa` files in the `reference` dir of your TCR assembly `run_dir`:
```python
from pathlib import Path
import pandas as pd
import os
import collections
from dotenv import load_dotenv

load_dotenv()
tcr_toolbox_data_path = os.getenv("tcr_toolbox_data_path")

from tcr_toolbox.sequencing_analysis.reference import generate_assembly_nt_refs
```

```python
plate_sheet_dir = Path(tcr_toolbox_data_path,
    "/tcr_toolbox_tcr_assembly_runs/numbered_run_name/plate_sheets"
)
tcr_df_dict = collections.defaultdict(pd.DataFrame)
for plate_xlsx in plate_sheet_dir.glob("*.xlsx"):
    tcr_df_dict[int(plate_xlsx.stem.split("_")[-1])] = pd.read_excel(plate_xlsx, index_col=0)

tcr_refs_df = pd.concat([tcr_df for tcr_df in tcr_df_dict.values()])
tcr_refs_df.reset_index(drop=True, inplace=True)
print(tcr_refs_df.loc[:, "name"].duplicated().any())
print(tcr_refs_df.loc[:, "name"].isna().any())

Output:
False
False
```

```python
generate_assembly_nt_refs(
     tcr_refs_df=tcr_refs_df,
     epitope_barcode_refs='[/path/to/your/]barcoded_epitope_sequences.fa',
     tcr_name_col_name='name',
     epitope_name_col_name='name',
     fasta_alpha_out_fname='[/path/to/your/]output_alpha.fa',
     fasta_beta_out_fname='[/path/to/your/]output_beta.fa',
     fasta_beta_epitope_plate_seq_out_fname='[/path/to/your/]output_beta_plate.fa',
     fasta_epitope_out_fname='[/path/to/your/]output_epi.fa',
     trimmed_beta_model_tcr_dict=None,
     model_epitope_dict=None,
     alpha_order_col_name='cdr3j_alpha_nt_order_primers',
     beta_order_col_name='cdr3j_beta_nt_order_primers',
     trim_assembly_primers_from_cdr3j=True,
     v_alpha_col='TRAV_IMGT_allele_collapsed',
     v_beta_col='TRBV_IMGT_allele_collapsed',
     add_v_gene_to_duplicate_filter=True,
     epitope_order_col_name='sequence',
     epitope_barcode_length=36, # check your barcode length and if needed adjust to 18 bp 
     read_length=150,
     p12_or_p20 = "p20",
     gtf=True,
     verbose=1
)

Output:
tcr_ref_df shape: 983
Removing TRAV + CDR3 alpha duplicates...
Shape before TRAV + CDR3 alpha duplicate removal: 983
Shape after duplicate removal: 909
Removing TRBV + CDR3 beta duplicates...
Shape before TRBV + CDR3 beta duplicate removal: 983
Shape after duplicate removal: 850
epitope_barcode_refs_df shape: 2778

Translation of first 5 epitope minigene sequences:
['EPPEVGSDCTTIHYNDMCNSSCMGGMNRRPI*', 'EPPEVGSDCTTIHYNDMCNSSCMGGMNRRPI*', ...]

Translation of first 5 reconstructed alpha chains:
['MACPGFLWALVISTCLEFSMAQTVTQSQPEMSVQEAETVTLSCTYDTSESDYYLFWYKQPPSRQMILVIRQEAYKQQNATENRFSVNFQKAAKSFSLKISDSQLGDAAMYFCAYRRGRSGGSEKLVFGKGTKLTVNPYI',
 'MKTFAGFSFLFLWLQLDCMSRGEDVEQSLFLSVREGDSSVINCTYTDSSSTYLYWYKQEPGAGLQLLTYIFSNMDMKQDQRLTVLLNKKDKHLSLRIADTQTGDSAIYFCAEKGSGGGADGLTFGKGTHLIIQPYI', ...]

Translation of first 5 reconstructed beta chains:
['MGTSLLCWMALCLLGADHADTGVSQNPRHKITKRGQNVTFRCDPISEHNRLYWYRQTLGQGPEFLTYFQNEAQLEKSRLLSDRFSAERPKGSFSTLEIQRTEQGDSAMYLCASSPQGSTGELFFGEGSRLTVLE',
 'MGCRLLCCAVLCLLGAVPIDTEVTQTPKHLVMGMTNKKSLKCEQHMGHRAMYWYKQKAKKPPELMFVYSYEKLSINESVPSRFSPECPNSSLLNLHLHALQPEDSALYLCASSQTPGQTGSPLHFGNGTRLTVTE', ...]

Finished writing references!
```

## 1. Initialize Project Directory:
```bash
mkdir pair_scan_project_dir
mkdir pair_scan_project_dir/references
mkdir pair_scan_project_dir/run_logs
mkdir pair_scan_project_dir/bulk_seq_data
mkdir pair_scan_project_dir/bulk_seq_data/epi
mkdir pair_scan_project_dir/bulk_seq_data/epi/run_logs
mkdir pair_scan_project_dir/bulk_seq_data/tcr
mkdir pair_scan_project_dir/bulk_seq_data/tcr/run_logs
```
These commands will create the following directory structure:
```bash
pair_scan_project_dir/
├── references/
├── run_logs/
└── bulk_seq_data/
    ├── epi/
    │   └── run_logs/
    └── tcr/
        └── run_logs/
```
## 2. Transfer Reference and Config Files:
If running on a remote computer, sync reference and configuration files using `rsync`:
```bash
rsync -avzP '[local_project_path]/references/tcr_beta_plate.fa' [your_clustrer]:[cluster_path_to_dir]/pair_scan_project_dir/references/
rsync -avzP '[local_project_path]/references/beta.fa' [your_cluster]:[cluster_path_to_dir]/pair_scan_project_dir/references/
rsync -avzP '[local_project_path]/references/epi.fa' [your_cluster]:[cluster_path_to_dir]/pair_scan_project_dir/references/
```
run_config templates are stored in `configs/sequencing_analysis/`. Choose config files based on the aligner that you will use: 

- **Apple Silicon ARM-based processors** → use **minimap2** (`*_minimap2.json`)  
  `bwa-mem` does **not** run on ARM CPUs only on x86-based processors. 
- **Intel/AMD x86-based processors** → use **bwa-mem** (`*_bwa.json`)  
  Preferrably, run `bwa-mem` as it is slightly more sensitive than minimap2 in internal benchmarking.

  For bulk sequencing, only use the `_udi-tcr-_` and `p20-minigene` config files. 
  
```bash
rsync -avzP '[tcr_toolbox_data_path]/tcr_toolbox_datasets/align_count_run_configs/run_config_count_umi_cell_150bp_pair_scan_bwa.json' [your_cluster]:[cluster_path_to_dir]/pair_scan_project_dir/
rsync -avzP '[tcr_toolbox_data_path]/tcr_toolbox_datasets/align_count_run_configs/run_config_count_reads_bulk_150bp_udi-tcr_bwa.json' [your_cluster]:[cluster_path_to_dir]/pair_scan_project_dir/bulk_seq_data/tcr/
rsync -avzP '[tcr_toolbox_data_path]/tcr_toolbox_datasets/align_count_run_configs/run_config_count_reads_bulk_150bp_p20-minigene_bwa.json' [your_cluster]:[cluster_path_to_dir]/pair_scan_project_dir/bulk_seq_data/epi/
```
## 3. Transfer FASTQ files: 
```bash
rsync -avzP '[local_project_path]/seq_data/pair_scan_project_dir/'*.fastq.gz [your_cluster]:[cluster_path_to_dir]/pair_scan_project_dir/
```

## 4. Adjust parameters in your config files:  
If you run remotely, first rsync your barcodes file from the `tcr_toolbox_dataset` directory on dropbox: 
```bash
rsync -avzP '[tcr_toolbox_data_path]/tcr_toolbox_datasets/pair_scan_luna_plate_seq/barcodes/barcodes_for_counting.tsv' [your_cluster]:~/tcr_toolbox_datasets/pair_scan_luna_plate_seq/barcodes/
```
Then, adjust parameters in your config files. In this example: 
  - `project_dir`: `[cluster_path_to_dir]/pair_scan_project_dir`
  - `reference_file`: `[cluster_path_to_dir]/pair_scan_project_dir/references/tcr_beta_plate.fa`
  - `barcode_file`: `/tcr_toolbox_datasets/pair_scan_luna_plate_seq/barcodes/barcodes_for_counting.tsv`
  - `threads`: number of cpus. For example, `8` if you 8 cpus available. 

**The JSON below is an illustrative example only — do not copy it as-is.** Edit the `project_dir`, `reference_file`, and `barcode_file` values in **your own** rsynced config file to match your setup. The `# adjust` comments are also not valid JSON; they exist here only to highlight which fields to change.
```bash
   {
        "project_dir": "/path/to/your/pairscan_project", # adjust
        "reference_file": "/path/to/reference.fa", # adjust
        "barcode_file": "/path/to/barcodes/barcodes_for_counting.tsv", # adjust
        "bc_pattern": "CCCCCCCCNNNNNNNN",
        "use_minimap2": true,
        "filter_secondary_alignments": true,
        "remove_tmp_dir": true,
        "threads": 8,
        "bwa_index_algo": "is",
        "epitope_barcode": true,
        "minimal_overlap_epitope": 36,
        "minimal_overlap_tcr": 104,
        "max_mismatches_epitope": 0,
        "max_mismatches_tcr": 0,
        "max_soft_5_end_epitope": 76,
        "max_soft_5_end_tcr": 0,
        "max_soft_3_end_epitope": 0,
        "max_soft_3_end_tcr": 0,
        "max_insertion_epitope": 0,
        "max_insertion_tcr": 0,
        "max_deletion_epitope": 0,
        "max_deletion_tcr": 0,
        "min_quality_3_end": 31,
        "minimap2_kmer_length": 19
    }
```

## 5. Start counting run: 
For example, run on a SLURM cluster: 
```bash
rsync -avzP '[tcr_toolbox_data_path]/tcr_toolbox_datasets/align_count_run_configs/run_count_slurm_sh_scripts/run_count_umi_cell.sh' [your_cluster]:[cluster_path_to_dir]/pair_scan_project_dir/
rsync -avzP '[tcr_toolbox_data_path]/tcr_toolbox_datasets/align_count_run_configs/run_count_slurm_sh_scripts/run_count_reads_bulk_slurm.sh' [your_cluster]:[cluster_path_to_dir]/pair_scan_project_dir/

cd [cluster_path_to_dir]/pair_scan_project_dir
sbatch run_count_umi_cell.sh run_config_count_umi_cell_150bp_pair_scan_bwa.json

cd [cluster_path_to_dir]/pair_scan_project_dir/bulk_seq_data/tcr
sbatch run_count_reads_bulk_slurm.sh run_config_count_reads_bulk_150bp_udi-tcr_bwa.json

cd [cluster_path_to_dir]/pair_scan_project_dir/bulk_seq_data/epi
sbatch run_count_reads_bulk_slurm.sh run_config_count_reads_bulk_150bp_p20-minigene_bwa.json
```


## 6. Monitor progress of active counting run: 
Progress of a counting run can be monitored in `run_logs/[file name].err`: 
```bash 
vi run_logs/run_count_umi_cell_315431.err
# OR
tail run_logs/run_count_umi_cell_315431.err
```
Familiarize yourself with CIGAR and MD alignment tags to understand how your alignments are filtered. 
Epi alignments show a 75–78 bp softclip because the epitope reference is trimmed to only 36 bp: the DNA barcode, plus a few epitope-coding bases padded in only if the barcode is shorter than 36 bp. The rest of the read's epitope-coding sequence isn't in the reference, so the aligner soft-clips it.
```bash
Started cutadapt trimming: UDI28_TCCAACGC-AAGTCCAA_S27_R1_001_extracted.fastq.gz
Applying 3' quality trimming: 31
Started bwa mem alignment: UDI28_TCCAACGC-AAGTCCAA_S27_R1_001_extracted_trimmed.fastq.gz
% primary aligned (UDI28_TCCAACGC-AAGTCCAA_S27_R1_001_extracted_trimmed.fastq.gz): 91.00% (3343486 reads)
Alignment UDI28_TCCAACGC-AAGTCCAA_S27_R1_001_extracted_trimmed.fastq.gz completed.
64 secondary alignments detected! Removing secondary alignments!
Top 10 most common cigar strings: [('75S36M', 1322153), ('76S36M', 919481), ('74S36M', 249603), ('72S36M', 233134), ('71S36M', 163714), ('104M', 114465), ('73S36M', 114395), ('105M', 91595), ('70S36M', 69576), ('68S36M', 65370)]
Top 10 most common MD tags: [('36', 4626253), ('104', 99727), ('11G24', 89781), ('105', 71727), ('14C21', 42485), ('32', 41289), ('26G9', 36897), ('8C27', 36160), ('31G4', 33564), ('33', 22067)]
Filtering BAM file...
Top 10 most common cigar strings: [('75S36M', 1219995), ('76S36M', 864257), ('74S36M', 231911), ('72S36M', 214749), ('71S36M', 153979), ('73S36M', 104423), ('104M', 99720), ('105M', 71217), ('70S36M', 65948), ('68S36M', 62399)]
Top 10 most common MD tags: [('36', 4615912), ('104', 99720), ('105', 71217)]
```

## 7. After a counting run is done, you can rsync the data from your remote project_dir to your local project_dir: 
```bash
rsync -avzP --exclude='*fastq.gz' --exclude='reference/' [your_cluster]:[cluster_path_to_dir]/pair_scan_project_dir/* '[local_project_path]/seq_data/pair_scan_project_dir'/
```
