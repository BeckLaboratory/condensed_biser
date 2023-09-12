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
ln -s PIPELINE_DIR/Snakefile

snakemake -j 1 --profile profiles/default --config sample=SAMPLE assembly=FASTA_PATH
```

Parameters "sample" and "assembly" must be 

## Running multiple samples

To run multiple samples, create a TSV file called "samples.tsv" or "samples.xlsx" in the project directory
(`PROJECT_DIR`) with two columns:

* `SAMPLE`: Name of the sample.
* `ASSEMBLY`: Path to the masked assembly.

By default, the pipeline will run all samples.

The pipeline is distributed with a slurm profile that tells Slurm how to schedule jobs. This profile can be used which
can be used by linking the pipeline's "profiles" directory to the project directory. If you are using a custom profile,
skip linking and point snakemake to your profile. You may use the resource requests in the Slurm profile to develop one
for your environment. See Snakemake profile documentation for additional help.

```
cd PROJECT_DIR
ln -s PIPELINE_DIR/profiles
ln -s PIPELINE_DIR/Snakefile

snakemake -j 1 --profile profiles/default
```

Adjust the `-j` option to the number of simultaneous jobs.


## Configuration parameters

Optional configuration parameters may be defined on the command-line after `--config` or in a file. Default
configuration file names are "config.yaml" and "config.json" if both are present, "config.yaml" is used. A custom
configuration file name can be defined on the command line with the "config_file" option (e.g.
`--config test_config.yaml`), and the default configuration files are ignored if present.

### Parameters

Each parameter is listed with the expected type and default value. 

* assembly [str]: Path to an input assembly for single-sample mode. Use with "sample" option. Cannot be used with
  "sample_table" and no default sample table is loaded if present.
* config_file [str]: Configuration file name. Overrides "config.yaml" and "config.json" default configuration files.
* keep_temp [bool, False]: Keep temp files generated for BISER (location tunable with "tmp_dir" option). Set to "True"
  to tell the pipeline not to clear the temporary directory created for each BISER run.
* mask_condense [int, 1000]: Shortest unmasked region between masked loci. If the length of an unmasked locus between two masked
  loci is not larger than this value, the masked region is extended to include the short unmaksed locus.
* sample [str]: Sample name for single-sample mode. Use with "assembly" option. Cannot be used with
  "sample_table" and no default sample table is loaded if present.
* sample_table [str]: Override "samples.tsv" and "sample.xlsx" for the default table names. File names must end with ".tsv",
  ".tsv.gz", or ".xlsx". TSV files must be in tab-separated-values format, and XLSX files must be in Excel format.
* shell_prefix [str]: Prepend this to shell commands after "set -euo pipefail; " (BASH strict mode). May be used to load
  modules or set the environment needed for commands. For example, if "shell_prefix" is set to "module load biser; ",
  the "biser" command is run as "set -euo pipefail; module load biser; biser". Note that the ";" after commands is
  required or the prefix becomes part of the command itself (useful for setting environment variables, but probably
  not what most users want).
* tmp_dir [str]: Path to BISER runtime temporary files. "temp/tmp" by default. A subdirectory is created for
  each sample in the format "biser_{sample}"


