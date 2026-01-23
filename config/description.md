We gratefully acknowledge the authors, originating and submitting laboratories of the genetic sequences and metadata for sharing their work via INSDC.

#### Analysis

Our bioinformatic processing workflow can be found at [github.com/nextstrain/chikv](https://github.com/nextstrain/chikv) and includes:

- sequence alignment by a combination of [Nextclade](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli.html) and [MAFFT](https://mafft.cbrc.jp/alignment/software/).
- phylogenetic reconstruction using [IQTREE](http://www.iqtree.org/)
- ancestral state reconstruction and temporal inference using [TreeTime](https://github.com/neherlab/treetime)
- clade assignment via clade definitions derived from the literature.

#### Underlying data
These analyses are based on data we sourced from INSDC, and we provide these data here in the hope that it will facilitate further analysis.
While these data have been shared openly by those who generated them, it does not mean there should be free license to publish on this data.
Data generators should be cited where possible and collaborations should be sought in some circumstances. Please try to avoid scooping someone else's work. Reach out if uncertain.

- [sequences](https://data.nextstrain.org/files/workflows/chikv/sequences.fasta.xz)
- [metadata](https://data.nextstrain.org/files/workflows/chikv/metadata.tsv.gz)

