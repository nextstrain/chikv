cd ./data

augur index \
  --sequences full_data/sequences.fasta \
  --output full_data/sequence_index.tsv

augur filter \
  --sequences full_data/sequences.fasta \
  --sequence-index full_data/sequence_index.tsv \
  --metadata full_data/metadata.tsv \
  --metadata-id-columns Accession \
  --exclude-ambiguous-dates-by any \
  --exclude-all \
  --include-where country=Senegal \
  --output-sequences senegal/sequences.fasta \
  --output-metadata senegal/metadata.tsv \
