cd chikungunya

augur curate normalize-strings --metadata data/raw/metadata_raw.tsv --output-metadata data/tempfiles/metadata_normalized.tsv

augur curate format-dates --metadata data/tempfiles/metadata_normalized.tsv --output-metadata data/tempfiles/metadata_formatted_dates.tsv --date-fields Collection_Date --expected-date-formats %Y %Y-%m

# augur curate parse-genbank-location --metadata data/tempfiles/metadata_formatted_dates.tsv --output-metadata data/tempfiles/metadata_parsed_loc.tsv --location-field Geo_Location -> dont think this worked

augur curate transform-strain-name --metadata data/tempfiles/metadata_formatted_dates.tsv --output-metadata data/tempfiles/metadata_strain_name.tsv --backup-fields Accession
augur curate rename --metadata data/tempfiles/metadata_strain_name.tsv --output-metadata data/tempfiles/metadata_renamed_cols.tsv --field-map Country=country Collection_Date=date

cp data/tempfiles/metadata_renamed_cols.tsv data/full_data/metadata.tsv

cp data/raw/sequences.fasta data/full_data/sequences.fasta