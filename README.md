# Condensed BISER

Run BISER on condensed assemblies with masked regions omitted from the assembly. This process is needed
to identify segmental duplications (SDs) without confounding matches produced by mobile elements and simple
repeats.

## Input FASTA

The input FASTA file must be masked. Masking can be soft (masked bases are lower-case) or hard (masked bases
are "N"). Standard references are often soft-masked, but custom masking can be accomplished by running repeat
annotations, generating a merged BED file of them, and using BedTools maskfasta to mask the assembly FASTA
file using the repeat BED.

## Installation

Pull the project from GitHub.

```
git clone --recurse-submodules https://github.com/ChristineRBeck/condensed_biser.git
```

It is recommended to install the pipeline into one directory and run it in another directory. The directory
created when the pipeline is cloned from GitHub is referred to as `PIPELINE_DIR`, and the directory where
analysis takes place is referred to as `PROJECT_DIR`.

The tool depends on:
* BioPython
* Pandas
* Pysam
* Python 3
* Numpy
* BISER


## Running a single sample

To run a single sample, change to a project directory where the tool will run (`PROJECT_DIR`) and link profiles
from the directory where the pipeline was downloaded (`PIPELINE_DIR`). If you are running inside the pipeline
directory (not recommended), then skip the `ln` step.

```
cd PROJECT_DIR
ln -s PIPELINE_DIR/profiles

snakemake -j 1 --profile profiles/local --config sample=SAMPLE assembly=FASTA_PATH
```

Parameters "sample" and "assembly" must be 

## Running multiple samples

To run multiple samples, create a TSV file called "samples.tsv" or "samples.xlsx" in the project directory
(`PROJECT_DIR`) with two columns:

* `SAMPLE`: Name of the sample.
* `ASSEMBLY`: Path to the masked assembly.

By default, the pipeline will run all samples.

```
cd PROJECT_DIR
ln -s PIPELINE_DIR/profiles

snakemake -j 1 --profile profiles/local --config sample=SAMPLE assembly=FASTA_PATH
```

To run mulitple samples simultaneously, adjust the `-j` option.


## Configuration parameters

If the pipeline is being run from an assembly table named `sample.tsv`, no configuration parameters are required.

Optional configuration parameters may be defined on the command-line after `--config` or in a file. Default
configuration file names are "config.yaml" and "config.json" if both are present, "config.yaml" is used. A custom
configuration file name can be defined on the command line with the "config_file" option (e.g.
`--config test_config.yaml`), and the default configuration files are ignored if present.

Parameters:
* assembly: Path to an input assembly for single-sample mode. Use with "sample" option. Cannot be used with
  "sample_table" and no default sample table is loaded if present.
* config_file: Configuration file name. Overrides "config.yaml" and "config.json" default configuration files.
* mask_condense: Shortest unmasked region between masked loci. If the length of an unmasked locus between two masked
  loci is not larger than this value, the masked region is extended to include the short unmaksed locus.
* sample: Sample name for single-sample mode. Use with "assembly" option. Cannot be used with
  "sample_table" and no default sample table is loaded if present.
* sample_table: Override "samples.tsv" and "sample.xlsx" for the default table names. File names must end with ".tsv",
  ".tsv.gz", or ".xlsx". TSV files must be in tab-separated-values format, and XLSX files must be in Excel format.

## Distributing over a cluster

TODO: Configuring profiles
