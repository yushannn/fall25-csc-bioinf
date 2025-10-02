#!/usr/bin/env python3

def calculate_n50(fasta_file):
    """Calculate N50 for a FASTA file of contigs"""
    lengths = []
    
    with open(fasta_file, 'r') as f:
        current_length = 0
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                if current_length > 0:
                    lengths.append(current_length)
                    current_length = 0
            else:
                current_length += len(line)
        
        # Add the last contig
        if current_length > 0:
            lengths.append(current_length)
    
    if not lengths:
        return 0
    
    # Count number of contigs
    num_contigs = len(lengths)
    
    # Calculate total length
    total_length = sum(lengths)
    
    # Sort lengths in descending order
    lengths.sort(reverse=True)
    
    # Find N50
    cumulative_length = 0
    for length in lengths:
        cumulative_length += length
        if cumulative_length >= total_length / 2:
            return length
    
    return 0

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <fasta_file>")
        sys.exit(1)
    
    n50 = calculate_n50(sys.argv[1])
    print(f"N50: {n50}")
    
    # Additional statistics
    with open(sys.argv[1], 'r') as f:
        content = f.read()
        
    # Count contigs (number of header lines)
    num_contigs = content.count('>')
    
    # Calculate total length (excluding headers and newlines)
    lines = content.split('\n')
    total_length = 0
    for line in lines:
        if not line.startswith('>') and line:
            total_length += len(line)
            
    print(f"Number of contigs: {num_contigs}")
    print(f"Total length: {total_length} bp")