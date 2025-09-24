#!/bin/bash

# Evaluation script for comparing Python and Codon genome assemblers
# Runs both implementations on all datasets and compares runtime and N50

# Remove strict error handling for more graceful error handling
set -uo pipefail

# Function to format time in seconds to MM:SS format
format_time() {
    local seconds=$1
    printf "%d:%02d" $((seconds/60)) $((seconds%60))
}

# Function to safely get stats from n50.py output
get_stat() {
    local stats="$1"
    local pattern="$2"
    local field="$3"
    echo "$stats" | grep "$pattern" | cut -d' ' -f"$field" 2>/dev/null || echo "N/A"
}

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to safely run a command and capture output
safe_run() {
    local cmd="$1"
    local logfile="$2"
    
    if eval "$cmd" > "$logfile" 2>&1; then
        return 0
    else
        echo "Command failed: $cmd"
        echo "Error log:"
        cat "$logfile"
        return 1
    fi
}

# Check required commands
echo "Checking required commands..."
if ! command_exists python3; then
    echo "Error: python3 not found"
    exit 1
fi

if ! command_exists codon; then
    echo "Error: codon not found"
    exit 1
fi

echo "All required commands found"

# Print header
printf "%-10s %-10s %-10s %-10s %-10s %-10s\n" "Dataset" "Language" "Runtime" "N50" "Num_Contigs" "Total_Length"
printf "%s\n" "-------------------------------------------------------------------------------------------------------"

# List of datasets to process
datasets=("data1" "data2" "data3" "data4")

for dataset in "${datasets[@]}"; do
    echo "Processing dataset: ${dataset}"
    
    # Skip if dataset directory doesn't exist
    if [ ! -d "data/${dataset}" ]; then
        echo "Skipping ${dataset}: directory not found"
        continue
    fi
    
    # Check if required input files exist
    if [ ! -f "data/${dataset}/short_1.fasta" ] || [ ! -f "data/${dataset}/short_2.fasta" ] || [ ! -f "data/${dataset}/long.fasta" ]; then
        echo "Skipping ${dataset}: missing input files"
        continue
    fi
    
    # Run Python implementation
    echo "Running Python implementation on ${dataset}..."
    python_start_time=$(date +%s)
    
    if safe_run "python3 code/main.py data/${dataset}" "python_${dataset}.log"; then
        python_end_time=$(date +%s)
        python_runtime=$((python_end_time - python_start_time))
        python_runtime_formatted=$(format_time $python_runtime)
        
        # Get Python statistics
        if [ -f "data/${dataset}/contig.fasta" ]; then
            python_stats=$(python3 code/n50.py "data/${dataset}/contig.fasta" 2>/dev/null || echo "Error getting stats")
            python_n50=$(get_stat "$python_stats" "N50:" 2)
            python_num_contigs=$(get_stat "$python_stats" "Number of contigs:" 4)
            python_total_length=$(get_stat "$python_stats" "Total length:" 3)
        else
            echo "Warning: Python contig.fasta not found for ${dataset}"
            python_n50="N/A"
            python_num_contigs="N/A"
            python_total_length="N/A"
        fi
        
        # Print Python results
        printf "%-10s %-10s %-10s %-10s %-10s %-10s\n" "${dataset}" "python" "${python_runtime_formatted}" "${python_n50}" "${python_num_contigs}" "${python_total_length}"
    else
        echo "Python implementation failed for ${dataset}"
        printf "%-10s %-10s %-10s %-10s %-10s %-10s\n" "${dataset}" "python" "FAILED" "N/A" "N/A" "N/A"
    fi
    
    # Run Codon implementation
    echo "Running Codon implementation on ${dataset}..."
    codon_start_time=$(date +%s)
    
    if safe_run "codon run -release code/main_codon_simple.py data/${dataset}" "codon_${dataset}.log"; then
        codon_end_time=$(date +%s)
        codon_runtime=$((codon_end_time - codon_start_time))
        codon_runtime_formatted=$(format_time $codon_runtime)
        
        # Get Codon statistics
        if [ -f "data/${dataset}/contig_codon.fasta" ]; then
            codon_stats=$(python3 code/n50.py "data/${dataset}/contig_codon.fasta" 2>/dev/null || echo "Error getting stats")
            codon_n50=$(get_stat "$codon_stats" "N50:" 2)
            codon_num_contigs=$(get_stat "$codon_stats" "Number of contigs:" 4)
            codon_total_length=$(get_stat "$codon_stats" "Total length:" 3)
        else
            echo "Warning: Codon contig_codon.fasta not found for ${dataset}"
            codon_n50="N/A"
            codon_num_contigs="N/A"
            codon_total_length="N/A"
        fi
        
        # Print Codon results
        printf "%-10s %-10s %-10s %-10s %-10s %-10s\n" "${dataset}" "codon" "${codon_runtime_formatted}" "${codon_n50}" "${codon_num_contigs}" "${codon_total_length}"
    else
        echo "Codon implementation failed for ${dataset}"
        printf "%-10s %-10s %-10s %-10s %-10s %-10s\n" "${dataset}" "codon" "FAILED" "N/A" "N/A" "N/A"
    fi
    
    echo "Completed ${dataset}"
    echo "----------------------------------------"
done

echo ""
echo "Evaluation completed successfully!"
echo "Log files available: python_*.log, codon_*.log"
echo ""

# Show a summary comparison if both implementations worked
echo "Summary:"
echo "Dataset-wise comparison shows runtime and assembly quality metrics."
echo "N50 values indicate assembly contiguity - higher values are generally better."
echo "Number of contigs shows fragmentation - fewer contigs often indicate better assembly."