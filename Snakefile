sequences_full = "data/full_data/sequences.fasta"
metadata_full = "data/full_data/metadata.tsv"


configfile: "config/config.yaml"


rule all:
    input:
        "data/subsampled_data/sequences.fasta",
        "data/subsampled_data/metadata.tsv",
        "results/aligned.fasta",
        "results/tree.nwk",
        "results/branch_lengths.json",
        "results/traits.json",
        "results/nt_muts.json",
        "auspice/chikv.json",
        expand(
            "data/subsampled_data/{build}/metadata.tsv", build=config.get("builds_to_run")
        ),


rule index:
    input:
        sequences="data/full_data/sequences.fasta",
    output:
        index="data/full_data/sequence_index.tsv",
    shell:
        "augur index \
            --sequences {input.sequences} \
            --output {output.index}"


rule filter:
    input:
        sequences="data/full_data/sequences.fasta",
        index="data/full_data/sequence_index.tsv",
        metadata="data/full_data/metadata.tsv",
    output:
        sequences="data/subsampled_data/sequences.fasta",
        metadata="data/subsampled_data/metadata.tsv",
    # params, pull them from config. can be a function, gets wildcard as parameter, look at rsv repo

    shell:
        "augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.index} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession \
            --exclude-ambiguous-dates-by any \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --group-by country year month \
            --subsample-max-sequences 100 \
            --probabilistic-sampling"


rule filter_country:
    input:
        sequences=sequences_full,
        index="data/full_data/sequence_index.tsv",
        metadata=metadata_full,
    output:
        sequences="data/subsampled_data/{build}/sequences.fasta",
        metadata="data/subsampled_data/{build}/metadata.tsv",
    shell:
        "augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.index} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession \
            --exclude-all \
            --include-where country={wildcards.build} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}"


# should i exclude amgigous dates here?
# should i subsample?


rule align:
    input:
        sequences="data/subsampled_data/sequences.fasta",
        ref_seq="config/chikv_reference.gb",
    output:
        alignment="results/aligned.fasta",
    shell:
        "augur align \
            --sequences {input.sequences}\
            --reference-sequence {input.ref_seq}\
            --output {output.alignment}\
            --fill-gaps"


rule tree:
    input:
        alignment="results/aligned.fasta",
    output:
        tree="results/tree_raw.nwk",
    shell:
        "augur tree \
        --alignment {input.alignment} \
        --output {output.tree}"


rule refine:
    input:
        tree="results/tree_raw.nwk",
        alignment="results/aligned.fasta",
        metadata="data/subsampled_data/metadata.tsv",
    output:
        tree="results/tree.nwk",
        node_data="results/branch_lengths.json",
    shell:
        "augur refine \
        --tree {input.tree} \
        --alignment {input.alignment} \
        --metadata {input.metadata} \
        --output-tree {output.tree} \
        --output-node-data {output.node_data} \
        --timetree \
        --coalescent opt \
        --date-confidence \
        --date-inference marginal \
        --clock-filter-iqd 4"


rule traits:
    input:
        tree="results/tree.nwk",
        metadata="data/subsampled_data/metadata.tsv",
    output:
        node_data="results/traits.json",
    shell:
        "augur traits \
        --tree {input.tree} \
        --metadata {input.metadata} \
        --output-node-data {output.node_data} \
        --columns region country \
        --confidence"


rule ancestral:
    input:
        tree="results/tree.nwk",
        alignment="results/aligned.fasta",
    output:
        node_data="results/nt_muts.json",
    shell:
        "augur ancestral \
        --tree {input.tree} \
        --alignment {input.alignment} \
        --output-node-data {output.node_data} \
        --inference joint"


rule translate:
    input:
        tree="results/tree.nwk",
        ancestral_seq="results/nt_muts.json",
        ref_seq="config/chikv_reference.gb",
    output:
        node_data="results/aa_muts.json",
    shell:
        "augur translate \
        --tree {input.tree} \
        --ancestral-sequences {input.ancestral_seq} \
        --reference-sequence {input.ref_seq} \
        --output-node-data {output.node_data}"


rule export:
    input:
        tree="results/tree.nwk",
        metadata="data/subsampled_data/metadata.tsv",
        branch_lengths="results/branch_lengths.json",
        nt_muts="results/nt_muts.json",
        aa_muts="results/aa_muts.json",
    output:
        auspice="auspice/chikv.json"
    shell:
        "augur export v2 \
        --tree {input.tree} \
        --metadata {input.metadata} \
        --node-data {input.branch_lengths} \
                    {input.nt_muts} \
                    {input.aa_muts} \
        --output {output.auspice}"
