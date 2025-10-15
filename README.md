
This is the Nextstrain build for chikunguny virus (CHIKV).

## Input Data
Input metadata and sequences will be made available via <https://data.nextstrain.org>

These data are generously shared by labs around the world and deposited in NCBI genbank by the authors.
Please contact these labs first if you plan to publish using these data.
CHIKV sequences and metadata can be downloaded in the `/ingest` folder using
`nextstrain build --cpus 1 ingest` or `nextstrain build --cpus 1 .` if running directly from the `/ingest` directory.


Running the ingest pipeline produces `ingest/data/metadata.tsv`, `ingest/data/extended_metadata.tsv` and `"ingest/data/sequences.fasta"`.
-- change metadata name (extended)


## Use locally ingested data

Once you have run the ingest pipeline locally you can copy the files into the top-level `data` directory so that the main phylo workflow uses these files rather than downloading from s3:

```sh
mkdir -p data
for i in ingest/data/*/{metadata.tsv,sequences.fasta}; do cp $i ${i#ingest/}; done
```

### `ingest/vendored`

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage copies of ingest scripts in [`ingest/vendored`](./ingest/vendored), from [nextstrain/ingest](https://github.com/nextstrain/ingest). To pull new changes from the central ingest repository, first install `git subrepo`, then run:

See [ingest/vendored/README.md](./ingest/vendored/README.md#vendoring) for instructions on how to update the vendored scripts.


## Run Analysis Pipeline

The workflow produces whole genome and G gene trees for RSV-A and RSV-B.
To run the workflow, use `snakemake -j4 -p --configfile config/configfile.yaml` and `nextstrain view auspice` to visualise results.

## Installation

Follow the standard [installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.

## Data use

We gratefully acknowledge the authors, originating and submitting laboratories of the genetic sequences and metadata for sharing their work. Please note that although data generators have generously shared data in an open fashion, that does not mean there should be free license to publish on this data. Data generators should be cited where possible and collaborations should be sought in some circumstances.



