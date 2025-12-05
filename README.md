[![TCRtoolbox DOI](https://img.shields.io/badge/TCRtoolbox-10.5281/zenodo.17806454-blue.svg)](https://doi.org/10.5281/zenodo.17806454)
[![tcr_toolbox_data DOI](https://img.shields.io/badge/tcr_toolbox_data-10.5281/zenodo.17832522-blue.svg)](https://doi.org/10.5281/zenodo.17832522)
[![pre-print DOI](https://img.shields.io/badge/bioRxiv-10.1101/2025.04.28.651095-red.svg)](https://doi.org/10.1101/2025.04.28.651095)


# TCRtoolbox

TCRtoolbox is a collection of tools related to TCR sequences used in the Schumacher lab. 
TCRtoolbox provides command line interface (CLI) tools for designing and preparing TCR assembly runs (including TCR reference .fa file generation), analysing Illumina bulk TCR and Ag-barcoded sequencing data, and reconstructing TCR amino acid or nucleotide sequences from V/J/CDR3 information. 
In addition, the the tcr_toolbox package provides many helper functions for working with TCR sequences. 

## Getting started 

```bash
# Clone the repo, install environment, install TCRtoolbox:
git clone https://github.com/schumacherlab/TCRtoolbox.git
cd TCRtoolbox
conda env create -f tcr-toolbox_ARM_env.yml -n py312-tcr-toolbox 
conda activate py312-tcr-toolbox
pip install -e .

# Download the datasets on which TCRtoolbox depends (see section `TCRtoolbox datasets`) 
curl -L -o tcr_toolbox_data.zip "https://zenodo.org/record/17832522/files/tcr_toolbox_data.zip?download=1"
unzip -q tcr_toolbox_data.zip

# Setup TCRtoolbox/.env to change [path_to] to the path where you downloaded these files on your local system
vi .env
tcr_toolbox_data_path='[path_to]/tcr_toolbox_data'
v_gene_barcode_tracing_path='[path_to]/tcr_toolbox_data/tcr_toolbox_datasets/tcr_assembly/barcode_tracing'

### Tutorials #### 
# For detailed tutorials have a look at the markdown tutorials in: TCRtoolbox/tutorials

# Simulate a TCR assembly run:
tcr_toolbox run-tcr-assembly configs/tcr_assembly/simulation_run_config.json

# Analyse bulk TCR read sequencing data for a DNA sequencing library protocol:
tcr_toolbox count-reads-bulk configs/tcr_assembly/run_config_count_reads_bulk_150bp_custom-tcr_bwa.json

# Reconstruct TCRs from CDR3 and VDJ input:
tcr_toolbox reconstruct-tcrs-simple configs/tcr_reconstruction/tcr_reconstruction_simple.json
```


## Project layout

- `configs/` - configuration template files for the command line tools. 
- `tutorials/` - tutorials for running command tools using configuration file templates
- `tcr_toolbox/` - main Python package
	- `sequencing_analysis/` - Illumina bulk TCR beta or alpha chain and antigen sequencing read counter scripts.
	- `tcr_assembly/` - TCR assembly run design and preparation code. 
	- `tcr_parsing/` - parsers for different TCR input formats. 
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
curl -L -o tcr_toolbox_data.zip "https://zenodo.org/record/17832522/files/tcr_toolbox_data.zip?download=1"
unzip -q tcr_toolbox_data.zip
```

To be able to use tcr_toolbox_data you need to have a `.env`. You can use `vi`, or any text-editor of your liking to make/adapt the .env file. In the following you need to adapt `[path_to]` so the path point to where these files are located on your system.
```bash
vi
tcr_toolbox_data_path='[path_to]/tcr_toolbox_data'
v_gene_barcode_tracing_path='[path_to]/tcr_toolbox_data/tcr_toolbox_datasets/tcr_assembly/barcode_tracing'
```

## Tutorials

The repository contains a small set of standalone Markdown tutorials that explain how to run our pipelines. These can be found in the `tutorials` directory. Quick summary:

- [tutorials/assembly_tutorial.md](https://github.com/schumacherlab/TCRtoolbox/blob/main/tutorials/assembly_tutorial.md) : TCR assembly run design and preparation using a single command and a configuration file.
- [tutorials/count_bulk_illumina_seq_reads_cli_tutorial.md](https://github.com/schumacherlab/TCRtoolbox/blob/main/tutorials/count_bulk_illumina_seq_reads_cli_tutorial.md) : unique DNA sequencing library protocols can be counted using a single command and provided protocol-specific config templates, reference .fa file generation, and monitoring of bulk read counting jobs. For example, Moravec et al., Nat. Biotech. 2024 TCR reactivity screen beta chain bulk sequencing data can be counted. 
- [tutorials/tcr_reconstruction_tutorial.md](https://github.com/schumacherlab/TCRtoolbox/blob/main/tutorials/tcr_reconstruction_tutorial.md) : Minimal pipeline to reconstruct full length TCR sequences (`Leader+Va+CDR3a+Ja+constanta` and `Leader+Vb+CDR3b+Jb+constantb`) from V, J, CDR3 annotations and the constant sequence. Leader and constant sequence are optional.

## TCRtoolbox datasets

Many tcr toolbox functions and pipelines depend on several files provided outside of this code repository. The tcr_toolbox_data can be downloaded here: https://doi.org/10.5281/zenodo.17832522.

Or accessed through the command line using:
```bash
cd TCRtoolbox
curl -L -o tcr_toolbox_data.zip "https://zenodo.org/record/17832522/files/tcr_toolbox_data.zip?download=1"
unzip -q tcr_toolbox_data.zip
```

It contains the following files:
- `tcr_toolbox_datasets/` :
    - `tcr_assembly/` : contains files needed for the tcr assembly and TCR reconstruction pipelines. Some files are required for our robotics TCR assembly platform, and might not be strictly required for your purpose.
    - `tcr_reconstruction/` : contains files needed for tcr reconstruction (such as imgt reference reference sequences)
- `tcr_toolbox_tcr_assembly_runs/` : can be empty, when you use the TCR assembly pipeline output files will be written 
- `test` : minimal dataset and configs to test command line pipelines

## Development

Besides to the provided pipelines, tcr_toolbox can be used as a python package.

- Import package utilities in Python (example):

```python
from tcr_toolbox import tcr_assembly_pipeline
# see subpackages for specific functions and helpers
```

## Contributing

Contributions are welcome:
- Open an issue if you find a bug or want a new feature
- Send a pull request with tests that demonstrate fixes/improvements

## Acknowledgements & Contact

This repository is developed and maintained by members of the Schumacher lab. For questions, please open an issue or contact the maintainers listed in the repository metadata.

## License

TCRtoolbox is provided under the Apache 2.0 licence (see LICENCE.txt)
