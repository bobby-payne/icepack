#!/bin/bash

# Source directory containing input files
source_dir="./"

# Destination directory for output files
dest_dir="./"

# Loop through each file in the source directory
for file in "$source_dir"/*; do
    # Check if the file is a regular file (not a directory)
    if [ -f "$file" ]; then
        # Extract the filename (without the path)
        filename=$(basename "$file")
        
        # Apply the cdo daymean command to the file, removing the hourly part of the filename
        cdo daymean "$file" "$dest_dir/${filename%???????.*}".nc

    fi
done
