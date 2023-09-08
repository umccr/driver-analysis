#!/bin/bash

# Specify the input CSV file
csv_file="./7a28e15d-2038-4584-8caa-a090dc17bcac.csv"

# Read the CSV file and loop through each row
while IFS=',' read -r id file_id subject_id name volume_id volume_name path time_created unique_hash; do
    # Check if the current row is not the header row (assuming the first row is a header)
    if [[ "$id" != "id" ]]; then
        # Extract the "file_id" value (remove double quotes if present)
        file_id=$(echo "$file_id" | sed 's/"//g')
        
        # Download the file using the "ica" command
        ica files download "$file_id" ./out/
    fi
done < "$csv_file"
