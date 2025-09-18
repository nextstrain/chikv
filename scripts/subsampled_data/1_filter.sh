cd chikungunya/data

augur index \
  --sequences full_data/sequences.fasta \
  --output full_data/sequence_index.tsv

augur filter \
  --sequences full_data/sequences.fasta \
  --sequence-index full_data/sequence_index.tsv \
  --metadata full_data/metadata.tsv \
  --metadata-id-columns Accession \
  --exclude-ambiguous-dates-by any \
  --output-sequences subsampled_data/sequences.fasta \
  --output-metadata subsampled_data/metadata.tsv \
  --group-by country year month \
  --subsample-max-sequences 100 \
  --probabilistic-sampling
