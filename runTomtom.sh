#!/bin/bash

directory_path="EEB498"

# Loop through directories in the input folder
for pri_folder in "$directory_path"/*upstream_regions/; do
  if [ -d "$pri_folder" ]; then
    species_name=$(basename "$pri_folder")
    prefix="$directory_path/${species_name}"

    # Loop through other directories in the same input folder
    for sec_folder in "$directory_path"/*/; do
      # Exclude the primary folder
      if [ "$sec_folder" != "$pri_folder" ]; then
        echo "Processing folder: $sec_folder"
        species_name2=$(basename "$sec_folder")
        prefix2="$directory_path/${species_name2}"

        output_folder="${prefix}${prefix2}TomTom"
        tomtom "$pri_folder"/"streme" "$sec_folder"/"streme" -o "$output_folder"
        echo "TomTom completed for $prefix and $prefix2"
      fi
    done
  fi
done



#streme --p UpstreamSequence/CDOUGupstream_regions.fa --n  UpstreamSequence/CDOUGAllUpstreamRegions.fa