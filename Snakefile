sequences_full = "data/full_data/sequences.fasta"
metadata_full = "data/full_data/metadata.tsv"


configfile: "config/config.yaml"


datadir = "data/subsampled_data"


full_data_dir = "data/{source}/full_data/"
background_data_dir = "data/{source}/subsampled_data/"
country_data_dir = "data/{source}/subsampled_data/country/{build}/"
build_data_dir = "data/{source}/subsampled_data/country_w_background/{build}/"


wildcard_constraints:
    build=r"[^/]+",  #constrain build wildcard to not contain slashes, so i dont get AmbiguousRuleException
    source="manual|ingest",


rule all:
    input:
        "data/manual/subsampled_data/sequences.fasta",
        "data/manual/subsampled_data/metadata.tsv",
        "results/manual/aligned.fasta",
        "results/manual/tree.nwk",
        "results/manual/branch_lengths.json",
        "results/manual/traits.json",
        "results/ingest/traits.json",
        "results/manual/nt_muts.json",
        "auspice/manual/chikv.json",
        "auspice/ingest/chikv.json",
        expand(
            "data/manual/subsampled_data/country/{build}/metadata.tsv",
            build=config.get("builds_to_run"),
        ),
        expand(
            "data/manual/subsampled_data/country_w_background/{build}/sequences.fasta",
            build=config.get("builds_to_run"),
        ),
        expand(
            "data/manual/subsampled_data/country_w_background/{build}/metadata.tsv",
            build=config.get("builds_to_run"),
        ),
        "results/manual/colors.tsv",
        "results/manual/aligned_masked.fasta",


rule index:
    input:
        sequences=full_data_dir + "sequences.fasta",
    output:
        index=full_data_dir + "sequence_index.tsv",
    shell:
        "augur index \
            --sequences {input.sequences} \
            --output {output.index}"


rule filter:
    input:
        sequences=full_data_dir + "sequences.fasta",
        index=full_data_dir + "sequence_index.tsv",
        metadata=full_data_dir + "metadata.tsv",
        exclude="config/outliers.txt",
    output:  # will serve as background
        sequences=background_data_dir + "sequences.fasta",
        metadata=background_data_dir + "metadata.tsv",
    shell:
        "augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.index} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession accession \
            --exclude-ambiguous-dates-by any \
            --exclude-where qc.overallStatus=bad 'abbr_authors=Badar et al.' \
            --exclude {input.exclude} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --group-by country year \
            --subsample-max-sequences 200 \
            --probabilistic-sampling \
            --subsample-seed 1"


rule filter_country:
    input:
        sequences=full_data_dir + "sequences.fasta",
        index=full_data_dir + "sequence_index.tsv",
        metadata=full_data_dir + "metadata.tsv",
    output:
        sequences=country_data_dir + "sequences.fasta",
        metadata=country_data_dir + "metadata.tsv",
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


rule merge_samples:
    input:
        metadata_country=country_data_dir + "metadata.tsv",
        metadata_background=background_data_dir + "metadata.tsv",
        sequences_country=country_data_dir + "sequences.fasta",
        sequences_background=background_data_dir + "sequences.fasta",
    output:
        sequences=build_data_dir + "sequences.fasta",
        metadata=build_data_dir + "metadata.tsv",
    shell:
        "augur merge \
            --metadata country={input.metadata_country} background={input.metadata_background} \
            --metadata-id-columns strain name \
            --sequences {input.sequences_country} {input.sequences_background} \
            --output-metadata {output.metadata} \
            --source-columns source_{{NAME}} \
            --output-sequences {output.sequences}"


rule align:
    input:
        sequences=background_data_dir + "sequences.fasta",
        ref_seq="config/chikv_reference.gb",
    output:
        alignment="results/{source}/aligned.fasta",
    shell:
        "augur align \
            --sequences {input.sequences}\
            --reference-sequence {input.ref_seq}\
            --output {output.alignment}\
            --fill-gaps"


rule mask:
    input:
        alignment="results/{source}/aligned.fasta",
    output:
        alignment_masked="results/{source}/aligned_masked.fasta",
    shell:
        "augur mask \
        --sequences {input.alignment} \
        --mask-from-beginning 76 \
        --mask-from-end 513 \
        --output {output.alignment_masked}"


rule tree:
    input:
        alignment="results/{source}/aligned_masked.fasta",
    output:
        tree="results/{source}/tree_raw.nwk",
    shell:
        "augur tree \
        --alignment {input.alignment} \
        --output {output.tree}"


rule refine:
    input:
        tree="results/{source}/tree_raw.nwk",
        alignment="results/{source}/aligned_masked.fasta",
        metadata=background_data_dir + "metadata.tsv",
    output:
        tree="results/{source}/tree.nwk",
        node_data="results/{source}/branch_lengths.json",
    shell:
        "augur refine \
        --tree {input.tree} \
        --alignment {input.alignment} \
        --metadata {input.metadata} \
        --metadata-id-columns accession Accession \
        --output-tree {output.tree} \
        --output-node-data {output.node_data} \
        --timetree \
        --coalescent opt \
        --date-confidence \
        --date-inference marginal"


rule traits:
    input:
        tree="results/{source}/tree.nwk",
        metadata=background_data_dir + "metadata.tsv",
    output:
        node_data="results/{source}/traits.json",
    shell:
        "augur traits \
        --tree {input.tree} \
        --metadata {input.metadata} \
        --metadata-id-columns accession Accession \
        --output-node-data {output.node_data} \
        --columns region country \
        --confidence"


rule ancestral:
    input:
        tree="results/{source}/tree.nwk",
        alignment="results/{source}/aligned_masked.fasta",
    output:
        node_data="results/{source}/nt_muts.json",
    shell:
        "augur ancestral \
        --tree {input.tree} \
        --alignment {input.alignment} \
        --output-node-data {output.node_data} \
        --inference joint"


rule translate:
    input:
        tree="results/{source}/tree.nwk",
        ancestral_seq="results/{source}/nt_muts.json",
        ref_seq="config/chikv_reference.gb",
    output:
        node_data="results/{source}/aa_muts.json",
    shell:
        "augur translate \
        --tree {input.tree} \
        --ancestral-sequences {input.ancestral_seq} \
        --reference-sequence {input.ref_seq} \
        --output-node-data {output.node_data}"


rule colors:
    input:
        color_schemes="config/color_schemes.tsv",
        color_orderings="config/color_orderings.tsv",
        metadata=background_data_dir + "metadata.tsv",
    output:
        colors="results/{source}/colors.tsv",
    shell:
        """
        python scripts/assign-colors.py \
            --color-schemes {input.color_schemes} \
            --ordering {input.color_orderings} \
            --metadata {input.metadata} \
            --output {output.colors}
        """


rule export:
    input:
        tree="results/{source}/tree.nwk",
        metadata=background_data_dir + "metadata.tsv",
        branch_lengths="results/{source}/branch_lengths.json",
        nt_muts="results/{source}/nt_muts.json",
        aa_muts="results/{source}/aa_muts.json",
        lat_longs="config/lat_longs.tsv",
        colors="results/{source}/colors.tsv",
    output:
        auspice="auspice/{source}/chikv.json",
    params:
        auspice_config="config/auspice_config.json",
        geo_resolutions="country",
        colors="config/colors.tsv",
    shell:
        "augur export v2 \
        --tree {input.tree} \
        --metadata {input.metadata} \
        --metadata-id-columns accession Accession \
        --node-data {input.branch_lengths} \
                    {input.nt_muts} \
                    {input.aa_muts} \
        --geo-resolutions {params.geo_resolutions} \
        --colors {input.colors} \
        --lat-longs {input.lat_longs} \
        --auspice-config {params.auspice_config} \
        --output {output.auspice}"
