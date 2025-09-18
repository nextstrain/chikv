cd chikungunya

augur export v2 \
  --tree results/tree.nwk \
  --metadata data/subsampled_data/metadata.tsv \
  --node-data results/branch_lengths.json \
              results/nt_muts.json \
              results/aa_muts.json \
  --output auspice/chikv.json