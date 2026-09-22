[![TCRtoolbox DOI](https://img.shields.io/badge/TCRtoolbox-10.5281/zenodo.17806454-blue.svg)](https://doi.org/10.5281/zenodo.17806454)
[![tcr_toolbox_data DOI](https://img.shields.io/badge/tcr_toolbox_data-10.5281/zenodo.17832522-blue.svg)](https://doi.org/10.5281/zenodo.17832521)
[![pre-print DOI](https://img.shields.io/badge/bioRxiv-10.1101/2025.04.28.651095-red.svg)](https://doi.org/10.1101/2025.04.28.651095)


# TCRtoolbox
TCRtoolbox is an end-to-end Python package containing tools to design and prepare assembly of TCR and epitope libraries for functional screening and to analyze the resulting sequencing data.

TCRtoolbox provides command line interface (CLI) tools for:
- reconstructing TCR amino acid or nucleotide sequences from V/J/CDR3 information.
- preparing [T-RAP TCR library assembly](https://doi.org/10.1101/2025.04.28.651095) runs (including TCR reference .fa file generation).
- counting bulk TCR-barcode Illumina sequencing data, for [TCR library reactivity screens](https://doi.org/10.1038/s41587-024-02210-6).
- counting bulk epitope-barcode Illumina sequencing data, for Ag library [HANSolo dropout screens](https://doi.org/10.1038/s41587-022-01547-0).
- counting single-well paired TCR- and epitope-barcode sequencing data, to detect TCR-Ag pairs in library-on-library [PAIR-Scan screens](https://doi.org/10.64898/2026.09.07.749948).

Example notebooks are provided for:
- parsing CellRanger single-cell TCR-sequencing contig files into TCR clonotype tables, and designing TCR libraries by selecting top clonotypes for [T-RAP TCR assembly](https://doi.org/10.1101/2025.04.28.651095).
- designing Ag and minimal peptide libraries, and preparing Ag and minimal peptide library assembly runs (including epitope-barcode reference .fa file generation).
- analyzing PAIR-Scan sequencing data by counting detected TCR–Ag pairs, and calculating PAIR-Scan confidence and TCR–Ag pair enrichment values.
Downstream analysis of TCR library reactivity and Ag library dropout screen read counts is not part of TCRtoolbox itself, and can instead be performed with the external [python DESeq2 package](https://pydeseq2.readthedocs.io/en/stable/) or [R DESeq2 package](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

In addition, TCRtoolbox provides many helper functions for working with TCR and epitope sequences.

## Getting started 

```bash
# Clone the repo, install environment, install TCRtoolbox:
git clone https://github.com/schumacherlab/TCRtoolbox.git
cd TCRtoolbox
conda env create -f tcr-toolbox_ARM_env.yml -n py312-tcr-toolbox 
conda activate py312-tcr-toolbox
pip install -e .

# Download the datasets on which TCRtoolbox depends (see section `TCRtoolbox datasets`)
# Visit https://doi.org/10.5281/zenodo.17832521, download tcr_toolbox_data.zip from the latest version, and move it into the TCRtoolbox folder, then unzip it:
unzip -q tcr_toolbox_data.zip

# Setup TCRtoolbox/.env by changing [path_to] to the path where you downloaded and unzipped these files on your local system:
vi .env
tcr_toolbox_data_path='[path_to]/tcr_toolbox_data'
v_gene_barcode_tracing_path='[path_to]/tcr_toolbox_data/tcr_toolbox_datasets/tcr_assembly/barcode_tracing'

### Tutorials #### 
# For detailed tutorials have a look at the markdown tutorials in: TCRtoolbox/tutorials

# Simulate a TCR assembly run:
tcr_toolbox run-tcr-assembly configs/tcr_assembly/simulation_run_config.json

# Analyse bulk TCR read sequencing data for a DNA sequencing library protocol:
tcr_toolbox count-reads-bulk configs/sequencing_analysis/run_config_count_reads_bulk_150bp_custom-tcr_bwa.json

# Reconstruct TCRs from CDR3 and VDJ input:
tcr_toolbox reconstruct-tcrs-simple configs/tcr_reconstruction/tcr_reconstruction_simple.json
```


## Project layout

- `configs/` - configuration template files for the command line tools. 
- `tutorials/` - tutorials for running command line tools using configuration file templates. 
- `notebooks/` - example notebook tutorials for running TCRtoolbox tool functions. 
- `tcr_toolbox/` - main Python package. 
    - `epi_assembly` - Ag and minimal peptide library design and assembly preparation code.
	- `sequencing_analysis/` - bulk TCR- and epitope-barcode sequencing read counting, single-well paired TCR- and epitope-barcode UMI counting, and calculating PAIR-Scan confidence and TCR–Ag pair enrichment values.
	- `tcr_assembly/` - T-RAP TCR assembly run design and preparation code. 
	- `tcr_parsing/` - parsing CellRanger single-cell TCR-sequencing contig files into TCR clonotype tables, and counting clonotype frequencies. 
	- `tcr_reconstruction/` - full TCR amino acid and nucleotide sequence reconstruction from V/J/CDR3 information. 
	- `utils/` - utilities used across modules. 

## Installation 

Clone the package from github: 
```bash
git clone https://github.com/schumacherlab/TCRtoolbox.git
```
This project is distributed with multiple Conda environment YAML files tailored to specific platforms. Use Conda to create and activate an environment from the YAML file that matches your system.

Example (create and activate):
```bash
cd TCRtoolbox
conda env create -f tcr-toolbox_ARM_env.yml -n py312-tcr-toolbox
conda activate py312-tcr-toolbox
```

Which YAML file to use:
- `tcr-toolbox_x86_env.yml` on x86 machines. Has a bwa dependency that does not (yet) work on ARM machines (e.g. apple M processors).
- `tcr-toolbox_ARM_env.yml` on ARM machines. Replaces bwa with minimap2.

Install the TCRtoolbox package with pip inside your active conda environment (run pip in the same directory as the `pyproject.toml` is located):  
```bash 
cd TCRtoolbox
pip install -e .
```


Next, TCR toolbox depends on several files (see TCRtoolbox datasets for a more detailed explanation).
```bash
# Visit https://doi.org/10.5281/zenodo.17832521, download tcr_toolbox_data.zip from the latest version, and move it into the TCRtoolbox folder, then unzip it:
unzip -q tcr_toolbox_data.zip
```

To be able to use tcr_toolbox_data you need to have a `.env`. You can use `vi`, or any text-editor of your liking to make/adapt the .env file. In the following you need to adapt `[path_to]` so the path points to where these unzipped files are located on your system.
```bash
vi
tcr_toolbox_data_path='[path_to]/tcr_toolbox_data'
v_gene_barcode_tracing_path='[path_to]/tcr_toolbox_data/tcr_toolbox_datasets/tcr_assembly/barcode_tracing'
```

## Tutorials

The repository contains Markdown tutorials that explain how to run our pipelines. These can be found in the `tutorials` directory. Quick summary:

- [tutorials/assembly_tutorial.md](https://github.com/schumacherlab/TCRtoolbox/blob/main/tutorials/assembly_tutorial.md) : [T-RAP](https://doi.org/10.1101/2025.04.28.651095) TCR assembly run design and preparation using a single `run-tcr-assembly` command and a configuration file.
- [tutorials/count_bulk_illumina_seq_reads_cli_tutorial.md](https://github.com/schumacherlab/TCRtoolbox/blob/main/tutorials/count_bulk_illumina_seq_reads_cli_tutorial.md) : unique DNA sequencing library protocols can be counted using a single `count-reads-bulk` command and provided protocol-specific config templates. For example, [Moravec et al., 2024 TCR reactivity screen](https://doi.org/10.1038/s41587-024-02210-6) beta chain bulk sequencing data, or [HANSolo](https://doi.org/10.1038/s41587-022-01547-0) Ag-barcode dropout screen bulk sequencing data, can be counted. Tutorial further explains how to generate a reference .fa file for bulk sequence counting using the TCRtoolbox. 
- [tutorials/pair_scan_counting_tutorial.md](https://github.com/schumacherlab/TCRtoolbox/blob/main/tutorials/pair_scan_counting_tutorial.md) : counting single-well paired TCR- and epitope-barcode sequencing data, for instance for [PAIR-Scan](https://doi.org/10.64898/2026.09.07.749948) screens, using a single `count-umi-cell` command and a configuration file.
- [tutorials/tcr_reconstruction_tutorial.md](https://github.com/schumacherlab/TCRtoolbox/blob/main/tutorials/tcr_reconstruction_tutorial.md) : Minimal pipeline to reconstruct full length TCR sequences (`Leader+Va+CDR3a+Ja+constanta` and `Leader+Vb+CDR3b+Jb+constantb`) from V, J, CDR3 annotations and the constant sequence, using a single `reconstruct-tcrs-simple` command and a configuration file. Leader and constant sequence are optional.

## Notebooks 
The repository contains example .ipynb notebooks that explain library design and analysis workflows. These can be found in the `notebooks` directory. Quick summary:

- [notebooks/tcr_parsing_examples/parse_and_design_tcr_library_example.ipynb](https://github.com/schumacherlab/TCRtoolbox/blob/main/notebooks/tcr_parsing_examples/parse_and_design_tcr_library_example.ipynb) : parsing CellRanger single-cell TCR-sequencing contig files into TCRαβ amino acid clonotype tables, counting clonotype frequency, and designing a TCR library by selecting the top clonotypes for [T-RAP](https://doi.org/10.1101/2025.04.28.651095) TCR assembly command line preparation in [tutorials/assembly_tutorial.md](https://github.com/schumacherlab/TCRtoolbox/blob/main/tutorials/assembly_tutorial.md).
- [notebooks/epi_assembly_examples/design_ag_library_example.ipynb](https://github.com/schumacherlab/TCRtoolbox/blob/main/notebooks/epi_assembly_examples/design_ag_library_example.ipynb) : design and preparation of Ag-encoding minigene library oligonucleotide pools. 
- [notebooks/epi_assembly_examples/design_minimal_peptide_library_example.ipynb](https://github.com/schumacherlab/TCRtoolbox/blob/main/notebooks/epi_assembly_examples/design_minimal_peptide_library_example.ipynb) : design and preparation of minimal peptide-encoding library oligonucleotide pools. 
- [notebooks/sequencing_analysis/pair_scan_analysis_example.ipynb](https://github.com/schumacherlab/TCRtoolbox/blob/main/notebooks/sequencing_analysis/pair_scan_analysis_example.ipynb) : analyzing [PAIR-Scan](https://doi.org/10.64898/2026.09.07.749948) sequencing data by counting detected TCR–Ag pairs, and calculating PAIR-Scan confidence and TCR–Ag pair enrichment values.

## TCRtoolbox datasets

Many tcr toolbox functions and pipelines depend on several files provided outside of this code repository. The tcr_toolbox_data can be downloaded from Zenodo: https://doi.org/10.5281/zenodo.17832521 (this link always resolves to the latest version).
To download:
1. Open https://doi.org/10.5281/zenodo.17832521 in your browser.
2. Download `tcr_toolbox_data.zip` from the Files section.
3. Move it into the `TCRtoolbox` folder and unzip it:
```bash
unzip -q tcr_toolbox_data.zip
```

It contains the following files:
- `tcr_toolbox_datasets/` :
    - `epi_assembly/` : contains files needed for epitope assembly. 
	- `tcr_assembly/` : contains files needed for the tcr assembly and TCR reconstruction pipelines. Some files are required for our robotics TCR assembly platform, and might not be strictly required for your purpose.
    - `tcr_reconstruction/` : contains files needed for tcr reconstruction (such as imgt reference reference sequences)
- `tcr_toolbox_tcr_assembly_runs/` : can be empty, when you use the TCR assembly pipeline output files will be written 
- `test` : minimal dataset and configs to test command line pipelines

## Development

Besides to the provided pipelines, tcr_toolbox can be used as a python package.

- Import package utilities in Python (example):

```python
from tcr_toolbox.tcr_assembly import tcr_assembly_pipeline
# see subpackages for specific functions and helpers
```

## Contributing

Contributions are welcome:
- Open an issue if you find a bug or want a new feature
- Send a pull request with tests that demonstrate fixes/improvements

## Acknowledgements & Contact

This repository is developed and maintained by members of the Schumacher lab. For questions, please open an issue or contact the maintainers listed in the repository metadata.

## License

TCRtoolbox is provided under the Apache 2.0 licence (see LICENSE).
