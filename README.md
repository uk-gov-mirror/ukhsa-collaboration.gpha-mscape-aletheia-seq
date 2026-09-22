# AletheiaSeq

In Greek mythology, Aletheia is the divine personification and spirit of truth, sincerity, and disclosure. Her name translates literally to "unconcealedness" or "state of being unhidden".

AletheiaSeq is a python package intended to be used to parse the output from BLAST and (optionally) Skope according to an input YAML file to determine presence/absence of specific signals in a set of sequence reads. The intention is the loci included in the BLAST database are specific to a species or characteristic and therefore can be used as a signal of a "true" detection.

The config YAML details which (and how many) loci must be present in the read set to positively identify a species or trait (e.g., toxin or resistance loci). It also details any QC thresholds to be applied to the outputs such as blast percent identity thresholds or number of reads the loci should be identified in. The package will take any blastn result as input so could be used on contigs but filtering the assembly to remove any contigs with poor read support is advised as minimum hits would have to be set to `1` in this instance.

## Installation

Clone repo and create environment:

`git clone git@github.com:ukhsa-collaboration/gpha-mscape-aletheia-seq.git`

`conda create -n aletheia_seq python=3.12`

`conda activate aletheia_seq`

Installation for users:

`cd gpha-mscape-aletheia-seq`

`pip install .`

Installation for developers (installs code in editable mode):

`cd gpha-mscape-aletheia-seq`

`pip install --editable '.[dev]'`

## Usage

```
usage: aletheiaseq [-h] --yaml YAML --outfolder OUT_FOLDER {prepare,parse} ...

In Greek mythology, Aletheia is the divine personification and spirit of truth, sincerity, and disclosure. Her name translates literally to 'unconcealedness' or 'state of being unhidden'.
This package is intended to replicate loci-based presence/absence confirmation used in reference laboratories for speciation and characterisation of species of interest as a way of producing a 'ground-truth' from sequence data.

positional arguments:
  {prepare,parse}
    prepare             Use the prepare subcommand to check your input files are in the correct format and to produce an example bash script for running blast and (optionally) skope. AletheiaSeq will
                        not run these commandline programs for you and you do not need to run this validation every time you want to use the programme.
    parse               Use the parse subcommand to process blast and (optionally) skope outputs and produce a report summarising the loci detected against the criteria set out in the input yaml. To
                        produce a json compatible with the onyx analysis table you will need to use the --onyx flag.

options:
  -h, --help            show this help message and exit
  --yaml YAML, -y YAML  Path to YAML file containing loci and QC information required to run AletheiaSeq. See full documentation for YAMl file format.
  --outfolder OUT_FOLDER, -o OUT_FOLDER
                        Folder to store output files. Folder will be created if it does not exist, files will be overwritten if it does.
```

AletheiaSeq has two submethods `prepare` and `parse`.

Prepare should be used to check that the yaml file is valid and the correct sequences are defined. It will generate a bash script to run the blast command.

```
usage: aletheiaseq prepare [-h] --fasta REF_FASTA --blast_db BLAST_DB [--skope_idx SKOPE_IDX]

options:
  -h, --help            show this help message and exit
  --fasta REF_FASTA, -f REF_FASTA
                        Path to FASTA file containing reference sequences for loci of interest. N.B. This must be the fasta file used to make the blast database.
  --blast_db BLAST_DB, -db BLAST_DB
                        Path to blast db (filename prefix only).
  --skope_idx SKOPE_IDX, -s SKOPE_IDX
                        Path to Skope index. Index will not be validated but if provided the skope command will be included in the bash script generated.
```

Parse will process the blast results and produce an onyx analysis table JSON file (with the option to publish) as well as a summary html file. Both methods require the defining yaml file and an output folder so these should be specified before the subcommand in the commandline.

```
usage: aletheiaseq parse [-h] --sample_id SAMPLE --blast_out BLAST_OUT [--skope_out SKOPE_OUT] [--onyx_server {mscape,synthscape,devscape}] [--publish] [--no_dryrun]

options:
  -h, --help            show this help message and exit
  --sample_id SAMPLE, -id SAMPLE
                        Relevant ID for sample being processed. If outputs are being ingested into onyx this should be the CLIMB ID.
  --blast_out BLAST_OUT, -b BLAST_OUT
                        Path to blast output. N.B. outfmt must match what is detailed in the YAML to parse the output correctly.
  --skope_out SKOPE_OUT
                        Path to output of skope classify.
  --onyx_server {mscape,synthscape,devscape}
                        Onyx server for analysis record. If not provided no analysis table result will be generated.
  --publish, -p         Flag to indicate whether onyx analysis object should be pushed to the database
  --no_dryrun           Flag to indicate that onyx upload should not be a test upload (default is to perform test upload to prevent accidental publishing)
```

## Inputs

### YAML format

The YAML file should contain 3 sections:
- Schema details - information about the specific schema
- Blast details - the loci included in the blast database, filters to be applied and the outfmt string to be used in blastn (`prepare` will validate that the column names used in the filters are present in the outfmt string)
- Group details - the species/characteristics to be determined using the loci including details of minimum number of loci and any required loci

An example YAML file is provided in example.yml

YAML files can be validated by using the `prepare` command

## Outputs

AletheiaSeq will produce an HTML report describing the loci identified and optionally a JSON output compatible with the onyx analysis table.
