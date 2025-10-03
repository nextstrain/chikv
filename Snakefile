## imports
import pandas as pd


## config
configfile: "config/config.yaml"


## important directories
full_data_dir = "data/full_data/"
background_data_dir = "data/subsampled_data/"
country_data_dir = "data/subsampled_data/country/{country_build}/"
region_data_dir = "data/subsampled_data/region/{region_build}/"
country_build_data_dir = "data/subsampled_data/country_w_background/{build}/"
region_build_data_dir = "data/subsampled_data/region_w_background/{build}/"


## 
regions = ["Asia",	"Oceania",	"Africa",	"Europe",	"South America",	"North America"]



## wildcards
wildcard_constraints:
    country_build=r"[^/]+",  # country to filter for
    region_build=r"[^/]+",  # region to filter for
    build=r"[^/]+",  # can be a country or a region
    build_type="country|region",
    #constrain country_build wildcard to not contain slashes, so i dont get AmbiguousRuleException


rule all:
    input:
        # data
        # intermediate results
        # auspice files
        expand(
            "auspice/chikv_{country_build}.json",
            country_build=config.get("country_builds_to_run")
        ),
        expand(
            "auspice/chikv_{region_build}.json",
            region_build=config.get("region_builds_to_run")
        ),
        "auspice/chikv_general.json",




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
            --subsample-max-sequences 100 \
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
            --metadata-id-columns Accession accession \
            --exclude-all \
            --include-where country={wildcards.country_build}\
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}"


# should i exclude amgigous dates here?
# should i subsample?


rule filter_region:
    input:
        sequences=full_data_dir + "sequences.fasta",
        index=full_data_dir + "sequence_index.tsv",
        metadata=full_data_dir + "metadata.tsv",
    output:
        sequences=region_data_dir + "sequences.fasta",
        metadata=region_data_dir + "metadata.tsv",
    shell:
        "augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.index} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession accession \
            --exclude-all \
            --include-where region={wildcards.region_build}\
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}"

rule merge_samples_country:
    input:
        metadata_country=country_data_dir + "metadata.tsv",
        metadata_background=background_data_dir + "metadata.tsv",
        sequences_country=country_data_dir + "sequences.fasta",
        sequences_background=background_data_dir + "sequences.fasta",
    output:
        sequences="data/subsampled_data/country_w_background/{country_build}/"
        + "sequences_merged.fasta",
        metadata="data/subsampled_data/country_w_background/{country_build}/"
        + "metadata_merged.tsv",
    shell:
        "augur merge \
            --metadata country={input.metadata_country} background={input.metadata_background} \
            --metadata-id-columns accession Accession\
            --sequences {input.sequences_country} {input.sequences_background} \
            --output-metadata {output.metadata} \
            --source-columns source_{{NAME}} \
            --output-sequences {output.sequences}"


rule merge_samples_region:
    input:
        metadata_region=region_data_dir + "metadata.tsv",
        metadata_background=background_data_dir + "metadata.tsv",
        sequences_region=region_data_dir + "sequences.fasta",
        sequences_background=background_data_dir + "sequences.fasta",
    output:
        sequences="data/subsampled_data/region_w_background/{region_build}/"
        + "sequences_merged.fasta",
        metadata="data/subsampled_data/region_w_background/{region_build}/"
        + "metadata_merged.tsv",
    shell:
        "augur merge \
            --metadata region={input.metadata_region} background={input.metadata_background} \
            --metadata-id-columns accession Accession\
            --sequences {input.sequences_region} {input.sequences_background} \
            --output-metadata {output.metadata} \
            --source-columns source_{{NAME}} \
            --output-sequences {output.sequences}"


rule filter_out_short_reads:
    input:
        sequences="data/subsampled_data/{build_type}_w_background/{build}/"
        + "sequences_merged.fasta",
        metadata="data/subsampled_data/{build_type}_w_background/{build}/"
        + "metadata_merged.tsv",
    output:
        sequences="data/subsampled_data/{build_type}_w_background/{build}/"
        + "sequences.fasta",
        metadata="data/subsampled_data/{build_type}_w_background/{build}/"
        + "metadata.tsv",
    shell:
        "augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession accession \
            --min-length 1000 \
            --exclude-ambiguous-dates-by year \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}"


def get_sequences(wildcards):
    build_type = "region" if wildcards.build in regions else "country"
    if wildcards.build == "general":
        seq_path = background_data_dir + "sequences.fasta"
    elif wildcards.build in regions:
        seq_path = region_build_data_dir + "sequences.fasta"
    else:
        seq_path = country_build_data_dir + "sequences.fasta"
    return seq_path


def get_metadata(wildcards):

    if wildcards.build == "general":
        met_path = background_data_dir + "metadata.tsv"
    elif wildcards.build in regions:
        met_path = region_build_data_dir + "metadata.tsv"
    else:
        met_path = country_build_data_dir + "metadata.tsv"
    return met_path


rule align:
    input:
        sequences=get_sequences,
        ref_seq="config/chikv_reference.gb",
    output:
        alignment="results/{build}/aligned.fasta",
    shell:
        "augur align \
            --sequences {input.sequences}\
            --reference-sequence {input.ref_seq}\
            --output {output.alignment}\
            --fill-gaps"


rule mask:
    input:
        alignment="results/{build}/aligned.fasta",
    output:
        alignment_masked="results/{build}/aligned_masked.fasta",
    shell:
        "augur mask \
        --sequences {input.alignment} \
        --mask-from-beginning 76 \
        --mask-from-end 513 \
        --output {output.alignment_masked}"


rule tree:
    input:
        alignment="results/{build}/aligned_masked.fasta",
    output:
        tree="results/{build}/tree_raw.nwk",
    shell:
        "augur tree \
        --alignment {input.alignment} \
        --output {output.tree}"


rule refine:
    input:
        tree="results/{build}/tree_raw.nwk",
        alignment="results/{build}/aligned_masked.fasta",
        metadata=get_metadata,
    output:
        tree="results/{build}/tree.nwk",
        node_data="results/{build}/branch_lengths.json",
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
        --date-inference marginal \
        --clock-rate 5e-4"


rule traits:
    input:
        tree="results/{build}/tree.nwk",
        metadata=get_metadata,
    output:
        node_data="results/{build}/traits.json",
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
        tree="results/{build}/tree.nwk",
        alignment="results/{build}/aligned_masked.fasta",
    output:
        node_data="results/{build}/nt_muts.json",
    shell:
        "augur ancestral \
        --tree {input.tree} \
        --alignment {input.alignment} \
        --output-node-data {output.node_data} \
        --inference joint"


rule translate:
    input:
        tree="results/{build}/tree.nwk",
        ancestral_seq="results/{build}/nt_muts.json",
        ref_seq="config/chikv_reference_adjusted.gb",
    output:
        node_data="results/{build}/aa_muts.json",
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
        metadata=get_metadata,
    output:
        colors="results/{build}/colors.tsv",
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
        tree="results/{build}/tree.nwk",
        metadata=get_metadata,
        branch_lengths="results/{build}/branch_lengths.json",
        nt_muts="results/{build}/nt_muts.json",
        aa_muts="results/{build}/aa_muts.json",
        lat_longs="config/lat_longs.tsv",
        colors="results/{build}/colors.tsv",
    output:
        auspice="auspice/chikv_{build}.json",
        
    params:
        auspice_config="config/auspice_config_ingest.json",
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


rule align_all:
    input:
        sequences="data/full_data/sequences.fasta",
        ref_seq="config/chikv_reference.gb",
    output:
        alignment="data/full_data/aligned.fasta",
    shell:
        "augur align \
            --sequences {input.sequences} \
            --reference-name NC_004162 \
            --output {output.alignment} \
            --nthreads auto \
            --fill-gaps"    


rule subsample_E1:
    input:
        alignment="data/full_data/aligned.fasta"
    output:
        sequences="data/subsampled/E1/sequences.fasta"
    shell:
        "python scripts/subsample_E1.py --alignment {input.alignment} --output {output.sequences}"