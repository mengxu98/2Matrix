#!/bin/bash

# Pattern
file_pattern="$1"

# Source folder path
source_dir="$2"

# Excluded folder path
exclude_dir="$3"

# Destination folder path
destination_dir="$4"

# Check and copy folders and files
for item in "$source_dir"/"$file_pattern"*; do
    # Get the name of the file or folder
    item_name=$(basename "$item")
    
    # Check if the same folder or file exists in the destination folder
    if [ ! -e "$exclude_dir/$item_name" ]; then
        # Copy folder or file to the destination folder
        cp -R "$item" "$destination_dir"
        echo "Copied $item_name to $destination_dir"
    else
        echo "$item_name exists in $exclude_dir, skipping copy"
    fi
done
