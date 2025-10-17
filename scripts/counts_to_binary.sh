#!/bin/bash

# Set defaults
threshold=${1}
input_file=${2}
output_file=${3}

awk -v thresh="$threshold" '
NR==1 {
    # Print header as-is
    print $0
    next
}
{
    # Print gene name
    printf "%s", $1
    
    # Convert each value to TRUE/FALSE
    for(i=2; i<=NF; i++) {
        if($i >= thresh)
            printf "\t%s", "TRUE"
        else
            printf "\t%s", "FALSE"
    }
    printf "\n"
}' "$input_file" > "$output_file"

