cd chikungunya

augur tree \
  --alignment results/aligned.fasta \
  --output results/tree_raw.nwk

  augur refine \
  --tree results/tree_raw.nwk \
  --alignment results/aligned.fasta \
  --metadata results/filtered/metadata.tsv \
  --output-tree results/tree.nwk \
  --output-node-data results/branch_lengths.json \
  --timetree \
  --coalescent opt \
  --date-confidence \
  --date-inference marginal \
  --clock-filter-iqd 4

augur traits \
  --tree results/tree.nwk \
  --metadata data/metadata.tsv \
  --output-node-data results/traits.json \
  --columns region country \
  --confidence


augur ancestral \
  --tree results/tree.nwk \
  --alignment results/aligned.fasta \
  --output-node-data results/nt_muts.json \
  --inference joint

augur translate \
  --tree results/tree.nwk \
  --ancestral-sequences results/nt_muts.json \
  --reference-sequence config/chikv_reference.gb \
  --output-node-data results/aa_muts.json