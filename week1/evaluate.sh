#!/bin/bash

# Add strict error handling
set -euxo pipefail

# Evaluation script for comparing Python and Codon genome assemblers
# Runs both implementations on all datasets and compares runtime and N50

# Function to format time in seconds to MM:SS format
format_time() {
    local seconds=$1
    printf "%d:%02d" $((seconds/60)) $((seconds%60))
}

# Print header
printf "%-10s %-10s %-10s %-10s %-10s %-10s\n" "Dataset" "Language" "Runtime" "N50" "Num_Contigs" "Total_Length"
printf "%s\n" "-------------------------------------------------------------------------------------------------------"

# List of datasets to process
datasets=("data1" "data2" "data3" "data4")

for dataset in "${datasets[@]}"; do
    # Skip if dataset directory doesn't exist
    if [ ! -d "data/${dataset}" ]; then
        echo "Skipping ${dataset}: directory not found"
        continue
    fi
    
    # Run Python implementation
    echo "Running Python implementation on ${dataset}..."
    python_start_time=$(date +%s)
    python3 code/main.py "data/${dataset}" > /dev/null 2>&1
    python_end_time=$(date +%s)
    python_runtime=$((python_end_time - python_start_time))
    python_runtime_formatted=$(format_time $python_runtime)
    
    # Get Python statistics
    python_stats=$(python3 code/n50.py "data/${dataset}/contig.fasta")
    python_n50=$(echo "$python_stats" | grep "N50:" | cut -d' ' -f2)
    python_num_contigs=$(echo "$python_stats" | grep "Number of contigs:" | cut -d' ' -f4)
    python_total_length=$(echo "$python_stats" | grep "Total length:" | cut -d' ' -f3)
    
    # Print Python results
    printf "%-10s %-10s %-10s %-10s %-10s %-10s\n" "${dataset}" "python" "${python_runtime_formatted}" "${python_n50}" "${python_num_contigs}" "${python_total_length}"
    
    # Run Codon implementation
    echo "Running Codon implementation on ${dataset}..."
    codon_start_time=$(date +%s)
    codon run -release code/main_codon_simple.py "${dataset}" > /dev/null 2>&1
    codon_end_time=$(date +%s)
    codon_runtime=$((codon_end_time - codon_start_time))
    codon_runtime_formatted=$(format_time $codon_runtime)
    
    # Get Codon statistics
    codon_stats=$(python3 code/n50.py "data/${dataset}/contig_codon.fasta")
    codon_n50=$(echo "$codon_stats" | grep "N50:" | cut -d' ' -f2)
    codon_num_contigs=$(echo "$codon_stats" | grep "Number of contigs:" | cut -d' ' -f4)
    codon_total_length=$(echo "$codon_stats" | grep "Total length:" | cut -d' ' -f3)
    
    # Print Codon results
    printf "%-10s %-10s %-10s %-10s %-10s %-10s\n" "${dataset}" "codon" "${codon_runtime_formatted}" "${codon_n50}" "${codon_num_contigs}" "${codon_total_length}"
done

# Save runtime comparison to a log file
echo "Runtime comparison completed. Results are displayed above."