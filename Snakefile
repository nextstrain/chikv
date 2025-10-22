## imports
import pandas as pd # not sure if this is ever used


## config
configfile: "config/config.yaml"


## important directories
full_data_dir = "data/full_data/"
background_data_dir = "data/subsampled_data/" # random, balanced subsample from the full data
country_data_dir = "data/subsampled_data/country/{country_build}/" # data filtered for just the country
region_data_dir = "data/subsampled_data/region/{region_build}/" # data filtered for just the region
country_build_data_dir = "data/subsampled_data/country_w_background/{build}/" # country + background data
region_build_data_dir = "data/subsampled_data/region_w_background/{build}/" # region + background data


##
regions = ["Asia", "Oceania", "Africa", "Europe", "South America", "North America"] # to know whether we are building a country or data build


## wildcards
wildcard_constraints:
    #constrain build wildcards to not contain slashes, so i dont get AmbiguousRuleException
    country_build=r"[^/]+",  # country to filter for
    region_build=r"[^/]+",  # region to filter for
    build=r"[^/]+",  # can be a country or a region
    build_type="country|region",
    

rule all:
    input:
        # data
        "data/subsampled_data/E1/metadata.tsv",

        "data/full_data/metadata.tsv",
        "data/full_data/sequences.fasta",
        # intermediate results
        "results/E1/aligned.fasta",
        # auspice files
        expand(
            "auspice/chikv_{country_build}.json",
            country_build=config.get("country_builds_to_run"),
        ),
        expand(
            "auspice/chikv_{region_build}.json",
            region_build=config.get("region_builds_to_run"),
        ),
        #"auspice/chikv_E1.json",
        "auspice/chikv_general.json", # general build from background data


# === filtering and subsampling ===

rule index:
    message:
        """
        Creating index of full dataset
        """
    input:
        sequences=full_data_dir + "sequences.fasta",
    output:
        index=full_data_dir + "sequence_index.tsv",
    shell:
        "augur index \
            --sequences {input.sequences} \
            --output {output.index}"


rule filter:
    """
    Creates the global 'background' dataset by performing quality control and
    probabilistic subsampling on the full dataset.
    """
    message: "Creating balanced background dataset"
    input:
        sequences=full_data_dir + "sequences.fasta",
        index=full_data_dir + "sequence_index.tsv",
        metadata=full_data_dir + "metadata.tsv",
        exclude="config/outliers.txt", # accessions of any samples we want to exclude
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
            --exclude-where qc.overallStatus=bad \
            --exclude {input.exclude} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --group-by country year \
            --subsample-max-sequences 100 \
            --probabilistic-sampling \
            --subsample-seed 1"


rule filter_country:
    """Extracts all sequences for a single country to create a 'focal' dataset."""
    message: "Extracting sequences for: {wildcards.country_build}"
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
    """Extracts all sequences for a single region to create a 'focal' dataset."""
    message: "Extracting sequences for region: {wildcards.region_build}"
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
    """Merges a country-specific 'focal' set with the global 'background' set."""
    message: "Merging '{wildcards.country_build}' data with background data..."
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
    """Merges a region-specific 'focal' set with the global 'background' set."""
    message: 
        "Merging '{wildcards.region_build}' data with background data..."
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





rule filter_e1:
    "Use nextclade coverage information to find sequences that fully cover E1"
    message:
        "filtering for E1"
    input:
        sequences=full_data_dir + "sequences.fasta",
        index=full_data_dir + "sequence_index.tsv",
        metadata=full_data_dir + "metadata.tsv",
    output:
        sequences= "data/subsampled_data/E1/sequences_full.fasta",
        metadata= "data/subsampled_data/E1/metadata_full.tsv",
    shell:
        "augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.index} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession accession \
            --query 'E1>0.8' \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}"


rule subsample_e1:
    """Subsamples the E1 gene sequences to create the final analysis set."""
    message:
        "Creating balanced E1 set"
    input:
        sequences= "data/subsampled_data/E1/sequences_full.fasta",
        metadata= "data/subsampled_data/E1/metadata_full.tsv",
    output:
        sequences= "data/subsampled_data/E1/sequences.fasta",
        metadata= "data/subsampled_data/E1/metadata.tsv",
    shell:
        "augur filter \
            --metadata {input.metadata} \
            --sequences {input.sequences} \
            --metadata-id-columns Accession accession \
            --group-by country year \
            --subsample-max-sequences 100 \
            --probabilistic-sampling \
            --subsample-seed 1 \
            --output-metadata {output.metadata} \
            --output-sequences {output.sequences} \
            --output-log data/subsampled_data/E1/filter_log.tsv"

rule remove_ref_e1:
    input:
        sequences="data/subsampled_data/E1/sequences.fasta",
        metadata="data/subsampled_data/E1/metadata.tsv",
    output:
        sequences="data/subsampled_data/E1/sequences_wo_ref.fasta",
        metadata="data/subsampled_data/E1/metadata_wo_ref.tsv",
    shell:
        "augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession accession \
            --exclude config/ref_name.txt \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}"

rule align_e1:
    input:
        sequences= "data/subsampled_data/E1/sequences_wo_ref.fasta",
        ref_seq="config/chikv_reference_E1.gb",
    output:
        alignment="results/E1/aligned.fasta",
    log:
        "logs/align_e1.log",
    shell:
        "augur align \
            --sequences {input.sequences}\
            --reference-sequence {input.ref_seq}\
            --output {output.alignment}\
            --fill-gaps \
            --nthreads auto \
            1> {log}"


rule quality_control:
    """
    Performs final QC on a merged dataset, removing short sequences and those
    with ambiguous dates.
    """
    message: "Final QC on '{wildcards.build}"
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



# === build trees and auspice files ===


def get_sequences(wildcards):
    build_type = "region" if wildcards.build in regions else "country"
    if wildcards.build == "E1":
        return "data/subsampled_data/E1/sequences_wo_ref.fasta"
    if wildcards.build == "general":
        seq_path = background_data_dir + "sequences.fasta"
    elif wildcards.build in regions:
        seq_path = region_build_data_dir + "sequences.fasta"
    else:
        seq_path = country_build_data_dir + "sequences.fasta"
    return seq_path


def get_metadata(wildcards):
    if wildcards.build == "E1":
        return "data/subsampled_data/E1/metadata_wo_ref.tsv"
    if wildcards.build == "general":
        met_path = background_data_dir + "metadata.tsv"
    elif wildcards.build in regions:
        met_path = region_build_data_dir + "metadata.tsv"
    else:
        met_path = country_build_data_dir + "metadata.tsv"
    return met_path


def get_ref(wildcards):
    if wildcards.build == "E1":
        return "config/chikv_reference_E1.gb"
    else:
        return "config/chikv_reference_adjusted.gb"


def get_alignment_for_trees(wildcards): 
    if wildcards.build == "E1": # we don't need to do any masking cause it's just the E1 gene anyway
        return "results/E1/aligned.fasta"
    else:
        print(f"{wildcards.build}")
        return f"results/{wildcards.build}/aligned_masked.fasta"



rule align:
    """align build sequences to reference"""
    message:
        "Aligning {wildcards.build} sequences to reference"
    input:
        sequences=get_sequences,
        ref_seq="config/chikv_reference.gb",
    output:
        alignment="results/{build}/aligned.fasta",
    log:
        "logs/align_{build}.log",
    shell:
        "augur align \
            --sequences {input.sequences}\
            --reference-sequence {input.ref_seq}\
            --output {output.alignment}\
            --fill-gaps \
            --nthreads auto \
            1> {log}"


rule mask:
    input:
        alignment="results/{build}/aligned.fasta",
    output:
        alignment_masked="results/{build}/aligned_masked.fasta",
    wildcard_constraints:
        build="(?!E1$)[^/]+",
    shell:
        "augur mask \
        --sequences {input.alignment} \
        --mask-from-beginning 76 \
        --mask-from-end 513 \
        --output {output.alignment_masked}"


rule tree:
    input:
        alignment=get_alignment_for_trees,
    output:
        tree="results/{build}/tree_raw.nwk",
    log:
        "logs/tree_{build}.log",
    shell:
        "augur tree \
        --alignment {input.alignment} \
        --output {output.tree} \
        1> {log} 2>&1"


rule refine:
    input:
        tree="results/{build}/tree_raw.nwk",
        alignment=get_alignment_for_trees,
        metadata=get_metadata,
    output:
        tree="results/{build}/tree.nwk",
        node_data="results/{build}/branch_lengths.json",
    log:
        "logs/refine_{build}.log",
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
        --clock-rate 5e-4 \
        1> {log} 2>&1"


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
        alignment=get_alignment_for_trees,
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
        ref_seq=get_ref,
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


rule remove_ref:
    input:
        sequences="data/full_data/sequences.fasta",
        metadata="data/full_data/metadata.tsv",
        index="data/full_data/sequence_index.tsv",
    output:
        sequences="data/full_data/sequences_wo_ref.fasta",
        metadata="data/full_data/metadata_wo_ref.tsv",
    shell:
        "augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.index} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession accession \
            --exclude config/ref_name.txt \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}"


rule align_all:
    input:
        sequences="data/full_data/sequences_wo_ref.fasta",
        ref_seq="config/chikv_reference.gb",
    output:
        alignment="data/full_data/aligned.fasta",
    shell:
        "augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.ref_seq} \
            --output {output.alignment} \
            --nthreads auto \
            --fill-gaps"





rule update_example_data_wildcards:
    """This updates the files under example_data/ based on latest available data from data.nextstrain.org.

    The subset of data is generated by an augur filter call which:
    - sets the subsampling size to 50
    - applies the grouping from the config
    """
    message:
        "Update example data"
    input:
        sequences="data/full_data/sequences.fasta",
        metadata="data/full_data/metadata.tsv",
    output:
        sequences="example_data/sequences.fasta",
        metadata="example_data/metadata.tsv",
    params:
        strain_id=config["strain_id_field"],
        group_by=config["filter"]["group_by"],
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --sequences {input.sequences} \
            --group-by {params.group_by} \
            --subsample-max-sequences 50 \
            --subsample-seed 0 \
            --output-metadata {output.metadata} \
            --output-sequences {output.sequences}
        """


rule update_example_data:
    input:
        "example_data/sequences.fasta",
        "example_data/sequences.fasta",


rule clean:
    shell:
        "rm -rf results auspice data/subsampled_data"



rule download:
    message: "downloading sequences and metadata from data.nextstrain.org"
    output:
        metadata =  "data/full_data/metadata.tsv.gz",
        sequences = "data/full_data/sequences.fasta.xz"
    params:
        metadata_url = "http://data.nextstrain.org/files/workflows/chikv/metadata.tsv.gz",
        sequence_url = "http://data.nextstrain.org/files/workflows/chikv/sequences.fasta.xz"
    shell:
        """
        curl -fsSL --compressed {params.metadata_url:q} --output {output.metadata}
        curl -fsSL --compressed {params.sequence_url:q} --output {output.sequences}
        """

rule decompress_sequences:
    message: "decompressing sequences"
    input:
        sequences = "data/full_data/sequences.fasta.xz",
    output:
        sequences = "data/full_data/sequences.fasta"
    shell:
        """
        xz --decompress --keep {input.sequences} -c > {output.sequences}
        """

rule decompress_metadata:
    message: "decompressing sequences"
    input:
        sequences = "data/full_data/metadata.tsv.gz",
    output:
        sequences = "data/full_data/metadata.tsv"
    shell:
        """
        gzip --decompress --keep {input.sequences} -c > {output.sequences}
        """