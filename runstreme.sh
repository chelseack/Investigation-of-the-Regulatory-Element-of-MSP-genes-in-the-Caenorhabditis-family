#!/bin/bash

# Define the input folder
input_folder="UpstreamSequence"

# Check if the input folder exists
if [ ! -d "$input_folder" ]; then
  echo "Error: Input folder '$input_folder' not found."
  exit 1
fi

# Loop through files in the input folder
for file_path in "$input_folder"/*upstream_regions.fa; do
  # Check if there are matching files
  if [ -e "$file_path" ]; then
    # Extract file name without extension
    filename=$(basename -- "$file_path")
    prefix="${filename%%upstream_regions.fa}"

    # Create an output folder with the same name as the file
    output_folder="${prefix}Streme"

    # Run streme command
    streme --p "$file_path" --n "UpstreamSequence/${prefix}AllUpstreamRegions.fa" -dna -o "$output_folder"
streme --p UpstreamSequence/CAFRAupstream_regions.fa --n UpstreamSequence/CAFRAAllUpstremeRegions.fa -dna -o CAFRAStreme
    echo "Streme completed for $filename"

  fi
done

