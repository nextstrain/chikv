cd data/manual

augur index \
  --sequences full_data/sequences.fasta \
  --output full_data/sequence_index.tsv

augur filter \
  --sequences full_data/sequences.fasta \
  --sequence-index full_data/sequence_index.tsv \
  --metadata full_data/metadata.tsv \
  --metadata-id-columns Accession accession\
  --exclude-ambiguous-dates-by any \
  --output-sequences subsampled_data/10/sequences.fasta \
  --output-metadata subsampled_data/10/metadata.tsv \
  --group-by country year month \
  --subsample-max-sequences 10 \
  --probabilistic-sampling
