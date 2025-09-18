cd chikungunya
augur index \
  --sequences sequences.fasta \
  --output sequence_index.tsv

augur filter \
  --sequences sequences.fasta \
  --sequence-index sequence_index.tsv \
  --metadata metadata.tsv \
  --metadata-id-columns Accession \
  --output-sequences results/filtered/sequences.fasta \
  --output-metadata results/filtered/metadata.tsv \
