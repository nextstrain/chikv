cd chikungunya

augur tree \
  --alignment results/10/aligned.fasta \
  --output results/10/tree_raw.nwk

  augur refine \
  --tree results/10/tree_raw.nwk \
  --alignment results/10/aligned.fasta \
  --metadata data/manual/subsampled_data/10/metadata.tsv\
  --output-tree results/10/tree.nwk \
  --output-node-data results/10/branch_lengths.json \
  --timetree \
  --coalescent opt \
  --date-confidence \
  --date-inference marginal \
  --clock-filter-iqd 4

augur traits \
  --tree results/10/tree.nwk \
  --metadata data/manual/subsampled_data/10/metadata.tsv \
  --output-node-data results/10/traits.json \
  --columns region country \
  --confidence


augur ancestral \
  --tree results/10/tree.nwk \
  --alignment results/10/aligned.fasta \
  --output-node-data results/10/nt_muts.json \
  --inference joint


# augur translate \
#   --tree results/10/tree.nwk \
#   --ancestral-sequences results/10/nt_muts.json \
#   --reference-sequence config/chikv_reference.gb \
#   --output-node-data results/10/aa_muts.json