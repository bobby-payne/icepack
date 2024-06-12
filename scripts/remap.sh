#!/bin/bash

# Source directory containing input files
source_dir="./"

# Destination directory for output files
dest_dir="../sic/"

# Loop through each file in the source directory
for file in "$source_dir"/*; do
    # Check if the file is a regular file (not a directory)
    if [ -f "$file" ]; then
        # Extract the filename (without the path)
        filename=$(basename "$file")
        
        # Apply the cdo remapbil command to the file
        cdo remapbil,r360x180 "$file" "$dest_dir/$filename"

    fi
done
