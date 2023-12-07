#!/bin/bash

# The raw script was refence from: https://zhuanlan.zhihu.com/p/627482865?utm_id=0
# This script is used to download GEO data
# Usage:
#   bash GEO_data_download.sh [-S|--silent] [GEO_ID]
# Example:
#   bash GEO_data_download.sh GEO168408
#   bash GEO_data_download.sh -S GEO168408

# Default values
silent=false

# Function to check if a command is available
command_exists() {
    command -v "$1" &>/dev/null
}

# Parse command line options
while [[ "$#" -gt 0 ]]; do
    case "$1" in
    -S | --silent)
        silent=true
        shift
        ;;
    *)

        if [ -z "$geo_id" ]; then
            geo_id="$1"
        else
            echo "Error: Unknown argument '$1'"
            exit 1
        fi
        shift
        ;;
    esac
done

Check if GEO_ID is provided
if [ -z "$geo_id" ]; then
    echo "Error: Please provide a GSE number."
    exit 1
fi

# Create a new directory with the GEO ID
mkdir "$geo_id"

# Define additional wget options for silent or non-silent download
wget_options="-r -nH --cut-dirs=3 --timeout=600 --progress=dot"
if $silent; then
    wget_options="$wget_options -q"
fi

# Use nohup to run wget in the background with output redirected to a log file if in silent mode
if $silent; then
    nohup wget $wget_options "ftp://ftp.ncbi.nlm.nih.gov/geo/series/${geo_id:0:-3}nnn/${geo_id}/" >"${geo_id}/${geo_id}_download.log" 2>&1 &
    echo "The data download for '${geo_id}' is running silently in the background."
    echo "Check '${geo_id}_download.log' for progress."
    echo "The '${geo_id}_download.log' file is in the directory '${geo_id}'."
else
    wget $wget_options "ftp://ftp.ncbi.nlm.nih.gov/geo/series/${geo_id:0:-3}nnn/${geo_id}/"
    echo "The '${geo_id}' successed download."
fi
