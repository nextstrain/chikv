## imports
import pandas as pd # not sure if this is ever used


## config
configfile: "config/config.yaml"


## important directories
full_data_dir = "data/full_data/"
background_data_dir = "builds/" # random, balanced subsample from the full data


##
display_names = {"South_America": "South America", "North_America": "North America", "Reunion": "Réunion"}

def display_name(build):
    return display_names.get(build, build)


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
        # intermediate results
        # auspice files
        expand(
            "auspice/chikv_{build}.json",
            build=config.get("country_builds_to_run"),
        ),
        expand(
            "auspice/chikv_{build}.json",
            build=config.get("region_builds_to_run"),
        ),
        "auspice/chikv_E1.json",
        "auspice/chikv_global.json",


# === download and decompress ===



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
        sequences = "data/full_data/sequences_raw.fasta"
    shell:
        """
        xz --decompress --keep {input.sequences} -c > {output.sequences}
        """

rule decompress_metadata:
    message: "decompressing metadata"
    input:
        metadata = "data/full_data/metadata.tsv.gz",
    output:
        metadata = "data/full_data/metadata_raw.tsv"
    shell:
        """
        gzip --decompress --keep {input.metadata} -c > {output.metadata}
        """


# === filter and subsample ===



rule quality_control:
    """
    Performs preliminary QC on full dataset, removing short sequences and those
    with ambiguous dates.
    """
    message: "Preliminary QC on full data"
    input:
        sequences="data/full_data/sequences_raw.fasta",
        metadata="data/full_data/metadata_raw.tsv",
    output:
        sequences="data/full_data/sequences.fasta",
        metadata="data/full_data/metadata.tsv",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession accession \
            --min-length 1000 \
            --exclude-where 'qc.overallStatus=bad' \
            --exclude-where 'qc.overallStatus=' \
            --exclude-ambiguous-dates-by year \
            --query '`qc.overallStatus`.notnull()' \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}
        """



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
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.index}
        """


rule filter_background:
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
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.index} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession accession \
            --exclude-ambiguous-dates-by any \
            --exclude-where 'qc.overallStatus=bad' \
            --query '`qc.overallStatus`.notnull()' \
            --exclude {input.exclude} \
            --min-length 6000 \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --group-by country year \
            --subsample-max-sequences 500 \
            --probabilistic-sampling \
            --subsample-seed 1
        """


rule filter_geo:
    """Extracts all sequences for a single country to create a 'focal' dataset."""
    message: "Extracting sequences for: {wildcards.build}"
    input:
        sequences=full_data_dir + "sequences.fasta",
        index=full_data_dir + "sequence_index.tsv",
        metadata=full_data_dir + "metadata.tsv",
    output:
        sequences="builds/{build}/" + "sequences.fasta",
        metadata="builds/{build}/" + "metadata.tsv",
    params:
        geo_filter = lambda w: f"region='{display_name(w.build)}'" if w.build in config.get("region_builds_to_run") else f"country='{display_name(w.build)}'",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.index} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession accession \
            --exclude-all \
            --include-where {params.geo_filter} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}
        """

rule merge_samples_geo:
    """Merges a geo-specific 'focal' set with the global 'background' set."""
    message: "Merging '{wildcards.build}' data with background data..."
    input:
        metadata_geo="builds/{build}/" + "metadata.tsv",
        metadata_background=background_data_dir + "metadata.tsv",
        sequences_geo="builds/{build}/" + "sequences.fasta",
        sequences_background=background_data_dir + "sequences.fasta",
    output:
        sequences="builds/{build}/" + "all_sequences.fasta",
        metadata="builds/{build}/" + "all_metadata.tsv",
    shell:
        """
        augur merge \
            --metadata geo={input.metadata_geo} background={input.metadata_background} \
            --metadata-id-columns accession\
            --sequences {input.sequences_geo} {input.sequences_background} \
            --output-metadata {output.metadata} \
            --source-columns source_{{NAME}} \
            --output-sequences {output.sequences}
        """



rule filter_e1:
    "Use nextclade coverage information to find sequences that fully cover E1"
    message:
        "filtering for E1"
    input:
        sequences=full_data_dir + "sequences.fasta",
        index=full_data_dir + "sequence_index.tsv",
        metadata=full_data_dir + "metadata.tsv",
    output:
        sequences= "builds/E1/sequences_full.fasta",
        metadata= "builds/E1/metadata_full.tsv",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.index} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession accession \
            --query 'E1>0.8' \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}
        """

rule subsample_e1:
    """Subsamples the E1 gene sequences to create the final analysis set."""
    message:
        "Creating balanced E1 set"
    input:
        sequences= "builds/E1/sequences_full.fasta",
        metadata= "builds/E1/metadata_full.tsv",
    output:
        sequences= "builds/E1/sequences.fasta",
        metadata= "builds/E1/metadata.tsv",
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --sequences {input.sequences} \
            --metadata-id-columns Accession accession \
            --exclude-where authors='Xiao,P.' \
            --exclude-where authors='Feng,G.,Zhang,J.,Zhang,Y.,Li,C.,Zhang,D.,Li,Y.,Zhou,H.,Li,N.,Xiao,P.,Lu,H.' \
            --group-by country year \
            --subsample-max-sequences 2000 \
            --probabilistic-sampling \
            --subsample-seed 1 \
            --output-metadata {output.metadata} \
            --output-sequences {output.sequences} \
            --output-log "builds/E1"/filter_log.tsv"
        """

rule remove_ref_e1:
    input:
        sequences="builds/E1/sequences.fasta",
        metadata="builds/E1/metadata.tsv",
    output:
        sequences="builds/E1/sequences_wo_ref.fasta",
        metadata="builds/E1/metadata_wo_ref.tsv",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession accession \
            --exclude config/ref_name.txt \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}
            """

rule align_e1:
    input:
        sequences= "builds/E1/sequences_wo_ref.fasta",
        ref_seq="config/chikv_reference_E1.gb",
    output:
        alignment="builds/E1/aligned.fasta",
    log:
        "logs/align_e1.log",
    shell:
        """
        augur align \
            --sequences {input.sequences}\
            --reference-sequence {input.ref_seq}\
            --output {output.alignment}\
            --fill-gaps \
            --nthreads auto \
            1> {log}
        """




rule filter_global:
    """
    Subsamples full-length sequences for a global build
    """
    message: "Subsampling for a global build"
    input:
        sequences=full_data_dir + "sequences.fasta",
        index=full_data_dir + "sequence_index.tsv",
        metadata=full_data_dir + "metadata.tsv",
        exclude="config/outliers.txt", # accessions of any samples we want to exclude
    output:  # will serve as background
        sequences= "builds/global/" + "sequences.fasta",
        metadata="builds/global/" + "metadata.tsv",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.index} \
            --metadata {input.metadata} \
            --metadata-id-columns Accession accession \
            --exclude-ambiguous-dates-by any \
            --exclude-where 'qc.overallStatus=bad' \
            --exclude {input.exclude} \
            --min-length 10000 \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --group-by country year \
            --subsample-max-sequences 3000 \
            --probabilistic-sampling \
            --subsample-seed 1
        """


# === build trees and auspice files ===


def get_sequences(wildcards):
    if wildcards.build == "E1":
        return "builds/E1/sequences_wo_ref.fasta"
    else:
        seq_path = f"builds/{wildcards. build}/" + "all_sequences.fasta"
    return seq_path


def get_metadata(wildcards):
    if wildcards.build == "E1":
        return "builds/E1/metadata_wo_ref.tsv"
    else:
        met_path = f"builds/{wildcards.build}/" + "all_metadata.tsv"
    return met_path


def get_ref(wildcards):
    if wildcards.build == "E1":
        return "config/chikv_reference_E1.gb"
    else:
        return "config/chikv_reference_adjusted.gb"

rule colors:
    """generate color mapping file for geographic traits"""
    message:
        "Assigning colors for: {wildcards.build}"
    input:
        color_schemes="config/color_schemes.tsv",
        color_orderings="config/color_orderings.tsv",
        metadata=get_metadata,
    output:
        colors="builds/{build}/colors.tsv",
    shell:
        """
        python scripts/assign-colors.py \
            --color-schemes {input.color_schemes} \
            --ordering {input.color_orderings} \
            --metadata {input.metadata} \
            --output {output.colors}
        """


rule align:
    """align build sequences to reference"""
    message:
        "Aligning {wildcards.build} sequences to reference"
    input:
        sequences=get_sequences,
        ref_seq="config/chikv_reference.gb",
    output:
        alignment="builds/{build}/aligned.fasta",
    log:
        "logs/align_{build}.log",
    shell:
        """
        augur align \
            --sequences {input.sequences}\
            --reference-sequence {input.ref_seq}\
            --output {output.alignment}\
            --fill-gaps \
            --nthreads auto \
            1> {log}
        """


rule mask:
    "mask noncoding regions from analysis"
    message:
        "Masking first 76 and last 513 nucleotides "
    input:
        alignment="builds/{build}/aligned.fasta",
    output:
        alignment_masked="builds/{build}/aligned_masked.fasta",
    wildcard_constraints:
        build="(?!E1$)[^/]+",
    shell:
        """
        augur mask \
        --sequences {input.alignment} \
        --mask-from-beginning 76 \
        --mask-from-end 513 \
        --output {output.alignment_masked}
        """

def get_alignment_for_trees(wildcards):
    if wildcards.build == "E1":
        return "builds/E1/aligned.fasta"
    else:
        return f"builds/{wildcards.build}/aligned_masked.fasta"

rule tree:
    "build a tree usig the IQ-TREE maximum likelihood algorithm"
    message:
        "Building initial maximum likelihood tree for {wildcards.build}"
    input:
        alignment=get_alignment_for_trees,
    output:
        tree="builds/{build}/tree_raw.nwk",
    log:
        "logs/tree_{build}.log",
    shell:
        """
        augur tree \
        --alignment {input.alignment} \
        --nthreads auto \
        --output {output.tree} \
        1> {log}
        """


rule refine:
    "use TreeTime to get a time-resolved tree"
    message:
        "Inferring timetree for {wildcards.build}"
    input:
        tree="builds/{build}/tree_raw.nwk",
        alignment=get_alignment_for_trees,
        metadata="builds/{build}/metadata.tsv",
    output:
        tree="builds/{build}/tree.nwk",
        node_data="builds/{build}/branch_lengths.json",
    log:
        "logs/refine_{build}.log",
    shell:
        """
        augur refine \
        --tree {input.tree} \
        --alignment {input.alignment} \
        --metadata {input.metadata} \
        --metadata-id-columns accession Accession \
        --output-tree {output.tree} \
        --output-node-data {output.node_data} \
        --timetree \
        --coalescent opt \
        --root 'mid_point' \
        --date-confidence \
        --date-inference marginal \
        --clock-rate 5e-4 \
        1> {log}
        """


rule traits:
    "use TreeTime to infer region and country for internal and original nodes"
    message:
        "Reconstructing ancestral traits for {wildcards.build}"
    input:
        tree="builds/{build}/tree.nwk",
        metadata="builds/{build}/metadata.tsv",
    output:
        traits="builds/{build}/traits.json",
    shell:
        """
        augur traits \
        --tree {input.tree} \
        --metadata {input.metadata} \
        --metadata-id-columns accession Accession \
        --output-node-data {output.traits} \
        --columns region country \
        --confidence
        """


rule ancestral:
    """use TreeTime to infer Maximum Likelihood ancestral sequences for internal nodes"""
    message:
        "Inferring ancestral sequences"
    input:
        tree="builds/{build}/tree.nwk",
        alignment=get_alignment_for_trees,
    output:
        node_data="builds/{build}/nt_muts.json",
    log:
        "logs/align_{build}.log",
    shell:
        """
        augur ancestral \
        --tree {input.tree} \
        --alignment {input.alignment} \
        --output-node-data {output.node_data} \
        --inference joint \
        1> {log}
        """


rule translate:
    """translate nucleotide mutations into amino acid mutations for gene regions"""
    message:
        "Translating gene regions from nucleotides to amino acids"
    input:
        tree="builds/{build}/tree.nwk",
        ancestral_seq="builds/{build}/nt_muts.json",
        ref_seq=get_ref,
    output:
        node_data="builds/{build}/aa_muts.json",
    shell:
        """
        augur translate \
        --tree {input.tree} \
        --ancestral-sequences {input.ancestral_seq} \
        --reference-sequence {input.ref_seq} \
        --output-node-data {output.node_data}
        """





rule export:
    """Use pipeline outputs to generate Auspice JSON for visualization"""
    message:
        "Exporting Auspice JSON for {wildcards.build}"
    input:
        tree="builds/{build}/tree.nwk",
        metadata="builds/{build}/metadata.tsv",
        branch_lengths="builds/{build}/branch_lengths.json",
        nt_muts="builds/{build}/nt_muts.json",
        aa_muts="builds/{build}/aa_muts.json",
        traits="builds/{build}/traits.json",
        lat_longs="config/lat_longs.tsv",
        colors="builds/{build}/colors.tsv",
    output:
        auspice="auspice/chikv_{build}.json",
    params:
        auspice_config="config/auspice_config.json",
        geo_resolutions="country",
        colors="config/colors.tsv",
    shell:
        """
        augur export v2 \
        --tree {input.tree} \
        --metadata {input.metadata} \
        --metadata-id-columns accession Accession \
        --node-data {input.branch_lengths} \
                    {input.nt_muts} \
                    {input.aa_muts} \
                    {input.traits} \
        --geo-resolutions {params.geo_resolutions} \
        --colors {input.colors} \
        --lat-longs {input.lat_longs} \
        --auspice-config {params.auspice_config} \
        --output {output.auspice}
        """




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
        """
        rm -rf results auspice data
        """



