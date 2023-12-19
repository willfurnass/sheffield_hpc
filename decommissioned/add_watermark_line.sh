#!/bin/bash

# Define the lines to be added
line1=".. include:: /referenceinfo/imports/decommissioned/decom_watermark.rst"

# Find all .rst files recursively and add the lines if not present
find . -type f -name "*.rst" -exec sh -c "
    # Check if each line is already present
    if ! grep -q \"$line1\" \"\$0\"; then
        # Add the first line to the beginning of the file
        echo \"MISSING LINE 1 in \$0\ \"
        awk -v line1=\"$line1\" 'BEGIN {print line1 \"\\n\"} {print}' \"\$0\" > temp_file && mv temp_file \"\$0\"
    fi
" {} \;
