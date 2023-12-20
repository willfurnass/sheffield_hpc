#!/bin/bash

# Define the lines to be added
line1=".. include:: /referenceinfo/imports/decommissioned/decom_watermark.rst"

# Find all .rst files recursively and add the lines if not present
find . -type f -name "*.rst" -print | while read -r file; do
    # Check if the line is already present
    if ! grep -q "$line1" "$file"; then
        echo "MISSING LINE 1 in $file"
        # Find the last line starting and ending with ":" to avoid issues with :orphan: etc...
        lastColonLine=$(grep -n '^:.*:$' "$file" | tail -n 1 | cut -d ":" -f 1)
        
        awk -v line1="$line1" -v lastColonLine="$lastColonLine" '
            {
                print;
                if (NR == lastColonLine) {
                    print "\n" line1 "\n";
                }
            }' "$file" > temp_file && mv temp_file "$file"
    fi
done
