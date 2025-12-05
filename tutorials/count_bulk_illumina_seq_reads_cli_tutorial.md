# Bulk Illumina TCR beta or alpha chain read counting tutorial

## 0. Generate reference `.fa` files:
If you set `generate_illumina_refs = True` in `tcr_toolbox run-tcr-assembly run_assembly_run_config.json` (recommended), then your reference files can be found in the TCR assembly run directory of your library: `/tcr_toolbox_data/tcr_toolbox_tcr_assembly_runs/[your_run]/sequencing_quality_analysis/references/[your_run]_150bp_beta.fa`. If not, use the function `tcr_toolbox.sequencing_analysis.reference.generate_assembly_nt_refs()` in your counts analysis `.ipynb` notebook to generate reference `.fa` files in the `reference` dir of your TCR assembly `run_dir`: 
```python
from pathlib import Path
import pandas as pd

load_dotenv()
tcr_toolbox_data_path = os.getenv("tcr_toolbox_data_path")


plate_sheet_dir = Path(tcr_toolbox_data_path,
    "/tcr_toolbox_tcr_assembly_runs/[your_run]/plate_sheets"
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

generate_assembly_nt_refs(
     tcr_refs_df=tcr_refs_df,
     tcr_name_col_name="name",
     fasta_alpha_out_fname=tcr_toolbox_data_path+'/tcr_toolbox_tcr_assembly_runs/[your_run]/sequencing_quality_analysis/references/r8_alpha.fa',
     fasta_beta_out_fname=tcr_toolbox_data_path+'/tcr_toolbox_tcr_assembly_runs/[your_run]/sequencing_quality_analysis/references/r8_beta.fa',
     alpha_order_col_name="cdr3j_alpha_nt_order_primers",
     beta_order_col_name="cdr3j_beta_nt_order_primers",
     trim_assembly_primers_from_cdr3j=True,
     v_alpha_col="TRAV_IMGT_allele_collapsed",
     v_beta_col="TRBV_IMGT_allele_collapsed",
     add_v_gene_to_duplicate_filter=True,
     read_length=150,
     gtf=True,
     verbose=1,
     trim_constant_seq=True,
 )

    Output:
    tcr_ref_df shape: 2688
    Removing TRAV + CDR3 alpha duplicates...
    Shape before TRAV + CDR3 alpha duplicate removal: 2688
    Shape after duplicate removal: 2135
    Removing TRBV + CDR3 beta duplicates...
    Shape before TRBV + CDR3 beta duplicate removal: 2688
    Shape after duplicate removal: 1537

    Translation of first 5 reconstructed alpha chains:
    ['METLLGLLILWLQLQWVSSKQEVTQIPAALSVPEGENLVLNCSFTDSAIYNLQWFRQDPGKGLTSLLLIQSSQREQTSGRLNASLDKSSGRSTLYIAASQPGDSATYLCAVRPSDKLIFGTGTRLQVFPNI',
    'MKTFAGFSFLFLWLQLDCMSRGEDVEQSLFLSVREGDSSVINCTYTDSSSTYLYWYKQEPGAGLQLLTYIFSNMDMKQDQRLTVLLNKKDKHLSLRIADTQTGDSAIYFCAAAEGGFKTIFGAGTRLFVKANI', ...]

    Translation of first 5 reconstructed beta chains:
    ['MGTSLLCWMALCLLGADHADTGVSQNPRHKITKRGQNVTFRCDPISEHNRLYWYRQTLGQGPEFLTYFQNEAQLEKSRLLSDRFSAERPKGSFSTLEIQRTEQGDSAMYLCASELANEQFFGPGTRLTVLE',
    'MSISLLCCAAFPLLWAGPVNAGVTQTPKFRILKIGQSMTLQCTQDMNHNYMYWYRQDPGMGLKLIYYSVGAGITDKGEVPNGYNVSRSTTEDFPLRLELAAPSQTSVYFCASSTPGQGSYEQYFGPGTRLTVTE', ...]

    Finished writing references!
```

## 1. Initialize Project Directory:
```bash
# Local: 
mkdir /schumi/cr/EPI_CR_1_6_run_1
mkdir /schumi/cr/EPI_CR_1_6_run_1/configs
mkdir /schumi/cr/EPI_CR_1_6_run_1/run_logs

# Remote:
mkdir /[path/to/your/remote/project/dir]
mkdir /[path/to/your/remote/project/dir]/run_logs
mkdir /[path/to/your/remote/project/dir]/references
```

These commands will create the following directory structure:
```bash
# Local:
EPI_CR_1_6_run_1/
    ├── run_logs/
    └── configs/

# Remote:
[path/to/your/remote/project/dir]/
    ├── run_logs/
    └── references/
```

## 2. Copy the config template

Run configuration templates are stored in:  
`[tcr_toolbox_data_path]/tcr_toolbox_datasets/align_count_run_configs`

Choose the config file based on:

### 1. Your DNA library preparation protocol
We provide configuration templates for four bulk sequencing DNA library preparation protocols. Pick the config matching your screen:

- **Ag/epi minigene bulk sequencing:**  
  - `_p20-minigene_` P20 minigene vector.
  - `_mscv-minigene_` MSCV minigene vector.x
  
  *Used for Ag/epi library assembly QC, dropout screens, and PAIR-scan screen Ag baseline.*

- **Single TCR chain bulk sequencing:**  
  - `_udi-tcr_` standard Illumina UDI p5 + p7 primers. 
  - `_custom-tcr_` custom barcoded p5 + p7 primers.  

  *Used for TCR library assembly QC, Moravec pooled TCR screens, and PAIR-scan screen TCR baseline.*

### 2. Your CPU architecture (this determines the aligner)

- **ARM-based processors** → use **minimap2** (`*_minimap2.json`)  
  `bwa-mem` does **not** run on ARM CPUs only on x86-based processors. 

- **x86-based processors** → use **bwa-mem** (`*_bwa.json`)  
  Preferrably, run `bwa-mem` as it is slightly more sensitive than minimap2 in internal benchmarking. 

### Final configuration file format:
`run_config_count_reads_bulk_150bp_{DNA library prep protocol}_{minimap2|bwa}.json`

We will run `_custom-tcr_` counting on a MacBook Pro: `run_config_count_reads_bulk_150bp_custom-tcr_minimap2.json`. Copy the template:
```bash
# Local:
rsync -avzP [tcr_toolbox_data_path]/tcr_toolbox_datasets/align_count_run_configs/run_config_count_reads_bulk_150bp_custom-tcr_minimap2.json \
    /schumi/cr/EPI_CR_1_6_run_1/configs/run_config_count_reads_bulk_150bp_custom-tcr_minimap2.json

# Remote:
rsync -avzP [tcr_toolbox_data_path]/tcr_toolbox_datasets/align_count_run_configs/run_config_count_reads_bulk_150bp_custom-tcr_minimap2.json \
    [your_cluster]:/[path/to/your/remote/project/dir]/run_config_count_reads_bulk_150bp_custom-tcr_minimap2.json
```

## 3. Transfer FASTQ and (if running remote) reference .fa files: 
```bash
# Local
rsync -avzP [path/to/your/]FASTQ/source/dir /schumi/cr/EPI_CR_1_6_run_1/

# Remote:
rsync -avzP [path/to/your/]FASTQ/source/dir [your_cluster]:[/path/to/remote/project/dir]/
rsync -avzP [tcr_toolbox_data_path]/tcr_toolbox_tcr_assembly_runs/[your_run]/sequencing_quality_analysis/references/r8_beta.fa \
    [your_cluster]:[/path/to/remote/project/dir]/references/
```

## 4. Adjust parameters in your config files: 
In your `run_config_count_reads_bulk_150bp_custom-tcr_minimap2.json` copy adjust: 
- `project_dir`: to your project dir. In this example, `/schumi/cr/EPI_CR_1_6_run_1`. 
- `reference_file`: to your reference file. If you set `generate_illumina_refs = True` in `tcr_toolbox run-tcr-assembly run_assembly_run_config.json`, you can find your reference in the TCR assembly run_dir of your library: `[tcr_toolbox_data]/tcr_toolbox_tcr_assembly_runs/[your_run]_tcr_assembly/sequencing_quality_analysis/references/r8_beta.fa`
- `threads`: the number of cpus you want to use if 1 thread = 1 cpu on your computer. 
```bash
/configs/run_config_count_reads_bulk_150bp_custom-tcr_minimap2.json
{
    "project_dir": "/schumi/cr/EPI_CR_1_6_run_1", 
    "reference_file": "[tcr_toolbox_data_path]/tcr_toolbox_tcr_assembly_runs/[your_run]/sequencing_quality_analysis/references/r8_beta.fa",
    "minimal_overlap": 105,
    "max_mismatches": 0,
    "max_soft_5_end": 0,
    "max_soft_3_end": 0,
    "use_minimap2": true,
    "filter_secondary_alignments": true,
    "remove_tmp_dir": true,
    "threads": 8,
    "bwa_index_algo": "is",
    "epitope_barcode": true,
    "max_insertion": 0,
    "max_deletion": 0,
    "min_quality_3_end": 31,
    "minimap2_kmer_length": 19
}
```

## 5. Start counting run: 
On your MacBook in the terminal run:  
```bash
conda activate tcr_toolbox_env
cd /schumi/cr/EPI_CR_1_6_run_1/run_logs
tcr_toolbox count-reads-bulk ../configs/run_config_count_reads_bulk_150bp_custom-tcr_minimap2.json > count_reads_bulk.out 2>count_reads_bulk.err
```

When connected to a remote SLURM cluster first adjust in `run_count_reads_bulk_slurm.sh`:
-  `#SBATCH --cpus-per-task`: to your required number of cpus. In this example: `32`. 
-  `#SBATCH --mem`: to your required amount of memory. In this example: `32GB`. 
-  `#SBATCH --time`: to your required run time. In this example: 10 hours. 
```bash
cd project_dir
vi ./run_count_reads_bulk_slurm.sh
#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32 
#SBATCH --qos=cpu_qos
#SBATCH --mem=32GB
#SBATCH --time=0-10:00:00
#SBATCH --output=./run_logs/run_count_reads_bulk_%A.out
#SBATCH --error=./run_logs/run_count_reads_bulk_%A.err

mkdir -p ./logs

source ~/miniconda3/etc/profile.d/conda.sh
conda activate py312-tcr-toolbox

tcr_toolbox count-reads-bulk $1
```
Then run: 
```bash
conda activate py312-tcr-toolbox
cd [project_dir]
sbatch run_count_reads_bulk_slurm.sh TCRtoolbox/configs/sequencing_analysis/run_config_count_reads_bulk_150bp_custom-tcr_minimap2.json
``` 

## 6. Monitor progress of your active counting run: 
```bash
# Both local and remote:
cd run_logs
# Scroll through whole log file:
vi count_reads_bulk.err

# View tail of log file:
tail count_reads_bulk.err
```
Familiarize yourself with CIGAR and MD alignment tags to understand how your alignments are filtered: 
```bash
Started cutadapt trimming: CD8_lib_GAATCCA_1.fastq.gz
Applying 3' quality trimming: 31
Started minimap2 alignment: CD8_lib_GAATCCA_1_trimmed.fastq.gz
Alignment CD8_lib_GAATCCA_1_trimmed.fastq.gz completed.
filter_secondary_alignments is set to True, but no secondary alignments detected!
Top 10 most common cigar strings: [('106M', 1267758), ('105M', 249385), ('102M', 79346), ('107M', 41876), ('103M', 32986), ('104M', 31249), ('101M', 18709), ('99M', 12215), ('98M', 11213), ('100M', 10923)]
Top 10 most common MD tags: [('106', 1021088), ('105', 207385), ('102', 67747), ('107', 33405), ('103', 29160), ('104', 27111), ('101', 16462), ('99', 10870), ('98', 9864), ('96', 9775)]
Filtering BAM file...
Top 10 most common cigar strings: [('106M', 1021006), ('105M', 206125), ('107M', 33389), ('108M', 3447), ('109M', 635)]
Top 10 most common MD tags: [('106', 1021006), ('105', 206125), ('107', 33389), ('108', 3447), ('109', 635)]
```

## 7. Check whether your sequencing library contains all reference names in your reference library: 
In your count analysis .ipynb run: 
```python
from tcr_toolbox.sequencing_analysis.utils import overlap_between_ref_and_count_names

overlap_between_ref_and_count_names(
        count_file="/counts/CD4_lib_CTCATCT_1_trimmed_sorted_cigarmd_filtered_counts.csv",
        count_file_ref_name_col="reference_name",
        reference_file="[tcr_toolbox_data]/tcr_toolbox_tcr_assembly_runs/[your_run]/sequencing_quality_analysis/references/r8_150bp_beta.fa",
        library_ref_id_list=["P108"],
        print_diff_names=True,
    )

Output: 
"# intersecting: 812"
"# diff: 5"
"# union: 819"
"refs not in count names: {'1_7_I16_207_P108_513', '1_7_J7_222_P108_528', '1_7_G9_152_P108_458', '1_7_J11_226_P108_532', '1_8_C23_70_P108_760'}"
"count names not in ref subset: {'1_6_A10_9_MAP_CD8-Tpexh-IFNG', '1_6_C7_54_MAP_No-RNA'}"
```
