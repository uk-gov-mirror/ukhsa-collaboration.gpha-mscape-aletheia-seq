# AletheiaSeq

In Greek mythology, Aletheia is the divine personification and spirit of truth, sincerity, and disclosure. Her name translates literally to "unconcealedness" or "state of being unhidden".

AletheiaSeq is a python package intended to be used to parse the output from BLAST and (optionally) Skope according to an input YAML file to determine presence/absence of specific signals in a set of sequence reads. The intention is the loci included in the BLAST database are specific to a species or characteristic and therefore can be used as a signal of a "true" detection.

The config YAML details which (and how many) loci must be present in the read set to positively identify a species or trait (e.g., toxin or resistance loci). It also details any QC thresholds to be applied to the outputs such as blast percent identity thresholds or number of reads the loci should be identified in. The package will take any blastn result as input so could be used on contigs but filtering the assembly to remove any contigs with poor read support is advised as minimum hits would have to be set to `1` in this instance.

## Installation

Clone repo and create environment:

`git clone git@github.com:ukhsa-collaboration/gpha-mscape-aletheia-seq.git`

`conda env create -n aletheia_seq python=3.12`

`conda activate aletheia_seq`

Installation for users:

`cd gpha-mscape-aletheia-seq`

`pip install .`

Installation for developers (installs code in editable mode):

`cd gpha-mscape-aletheia-seq`

`pip install --editable '.[dev]'`

## Usage

```
project-name --input <path> --output <path>
```

## Inputs

### YAML format

The YAML file should contain 3 sections:
- Schema details - information about the specific schema
- Blast details - the loci included in the blast database, filters to be applied and the outfmt string to be used in blastn (`prepare` will validate that the column names used in the filters are present in the outfmt string)
- Group details - the species/characteristics to be determined using the loci including details of minimum number of loci and any required loci

An example YAML file is provided in docs/examples/example.yml

YAML files can be validated by using the `prepare` command

## Outputs

AletheiaSeq will produce an HTML report describing the loci identified and optionally a JSON output compatible with the onyx analysis table.
