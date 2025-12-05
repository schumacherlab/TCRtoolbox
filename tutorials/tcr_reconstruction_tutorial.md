# TCR_reconstruction
This script allows you to reconstruct full length tcr sequences (`Leader+Va+CDR3a+Ja+constanta` and `Leader+Vb+CDR3b+Jb+constantb`) from V, J, CDR3 annotations and the constant sequence (the leader and constant sequence are optional).

## Input format:

- Your dataset should contain the following components:

| Column name | content |
| --- | --- |
| TRAV | TRAV gene annotation conform to IMGT format |
| TRAJ | TRAJ gene annotation conform to IMGT format |
| TRBV | TRBV gene annotation conform to IMGT format |
| TRBJ | TRBJ gene annotation conform to IMGT format |
| cdr3_alpha_aa | The CDR3 alpha amino acid sequence |
| cdr3_beta_aa  | The CDR3 beta amino acid sequence |

## Running the TCR reconstruction CLI pipeline

Run the TCR reconstruction CLI tool:

```bash
tcr_toolbox reconstruct-tcrs-simple tcr_reconstruction_simple.json
```

The config should contain the following:
```json
{
  "include_leader": false, 
  "include_constant": false,
  "mouse_or_human": "human",
  "constant_beta": null,
  "constant_alpha": null,
  "exclude_c_fw": false,
  "dataframe": null,
  "dataset_file_path": "tcr_reconstruction_example_dataset.csv",
  "output_file_path": "tcr_reconstruction_output.csv",
  "verbose": false
}
```

* `include_leader`: It is possible to include/exclude the leader sequence to the V sequence by adapting:  `include_leader=true` (include) or `include_leader=false` (exclude).

* `include_constant` It is possible to include/exclude a custom constant region behind (behind the J region)
  * Please note that according to IMGT the amino acid at the intersection is made up of 1 nt from the J region, and 2 nts from the constant region. In practice most people regard this amino acid as either part of the variable J region for TRCa ((X)IQN) or part of the constant region for TCRb ((E)DL), eventhough 2/3 nts come from the constant gene. If you choose to add the constant region this amino acid that is shared between the J and constant region is placed between the J and constant region, according to the variable nt that is provided by the J region. If you do not choose to add the constant region, this amino WILL NOT BE INCLUDED at the end of the J sequence!

* `mouse_or_human`: Human or mouse (mouse v-genes are currently not supported)

* `constant_beta` and `constant_alpha`: The constant sequences. If you leave this empty no constant sequence will be added.

* `exclude_c_fw` If true it will remove any TCR of which the CDR3 does not start with a C and or does not end with an F or W

* `dataframe` When importing the function (instead of using the CLI) here you can put the dataframe, leave as null for the CLI.

* `dataset_file_path` If dataframe is empty a dataset file path needs to be provided containing your CDR3 sequences, V and J genes.

* `output_file_path` Path to store the output csv containing the full length TCR sequences.

* `verbose` Enabling this will give verbose output on the progress of the TCR sequence reconstruction.

* It is possible to choose a specific set of V/J reference sequences from IMGT/GENE-DB (https://www.imgt.org/genedb/)
  * 'functional': This set contains all functional V/J sequences as defined by IMGT/GENE-DB (without ORF and pseudogene annotations) 
  * 'after benchmark': this set contains all functional V/J sequences as defined by IMGT/GENE-DB (without ORF and pseudogene annotations). However, here the TRAJ58 (ORF) is added and 2 sequences are changed to reference sequences from Ensemble instead of IMGT.


## Benchmark:

### in vitro
The code served as the basis of the TCR robotics assisted TCR assembly pipeline, which has been used to reconstruct thousands of TCR sequences in vitro.

### in silico
The reconstruction algorithm has been benchmarked on an internal dataset (RootPath), and an publically available dataset (10x). Compared to the reconstruction of RootPath this reconstruction method matched near 100% of their reconstructions (900+ TCRA and TCRB). The few sequences that differed in full sequence are explained by a difference in assumptions (or possible error) between the unknown RootPath script and this method.
The 10x data was biological sample of 10k TCRS (50:50 TCRA:TCRB). The reconstruction accuracy of this dataset was >85% The remaining 15% was explained by biological and/or technical noise and by missing info for allelic differences (alleles are often not annotated).
