import sys
import os

def join_path(path, name):
    if path.endswith('/'):
        return path + name
    else:
        return path + '/' + name

def read_fasta(path, name):
    data = []
    full_path = join_path(path, name)
    
    with open(full_path, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) > 0 and line[0] != '>':
                data.append(line)
                
    print(f"{name}: {len(data)} sequences, length {len(data[0]) if data else 0}")
    return data

def read_data(path):
    short1 = read_fasta(path, "short_1.fasta")
    short2 = read_fasta(path, "short_2.fasta")
    long1 = read_fasta(path, "long.fasta")
    return short1, short2, long1

def reverse_complement(seq):
    """Calculate reverse complement of a DNA sequence"""
    result = ""
    for i in range(len(seq) - 1, -1, -1):
        base = seq[i]
        if base == 'A':
            result += 'T'
        elif base == 'T':
            result += 'A'
        elif base == 'G':
            result += 'C'
        elif base == 'C':
            result += 'G'
        else:
            result += 'N'
    return result

def main():
    if len(sys.argv) < 2:
        print("Usage: codon run main_codon_simple.py <data_directory>")
        return
    
    # Read data - fix the path issue
    data_path = join_path('.', sys.argv[1])
    short1, short2, long1 = read_data(data_path)
    
    print(f"Read {len(short1)} short1 sequences")
    print(f"Read {len(short2)} short2 sequences")  
    print(f"Read {len(long1)} long sequences")
    
    # Implement a simple De Bruijn Graph algorithm
    k = 25
    
    # Build k-mer dictionary
    kmer_dict = {}
    kmer_count = 0
    
    # Process all sequences
    all_sequences = short1 + short2 + long1
    
    for seq in all_sequences:
        # Add forward k-mers
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if kmer not in kmer_dict:
                kmer_dict[kmer] = 0
            kmer_dict[kmer] += 1
            
        # Add reverse complement k-mers
        rc_seq = reverse_complement(seq)
        for i in range(len(rc_seq) - k + 1):
            kmer = rc_seq[i:i+k]
            if kmer not in kmer_dict:
                kmer_dict[kmer] = 0
            kmer_dict[kmer] += 1
    
    print(f"Built k-mer dictionary with {len(kmer_dict)} unique {k}-mers")
    
    # Find most frequent k-mers and try to extend them
    sorted_kmers = sorted(kmer_dict.keys(), key=lambda x: kmer_dict[x], reverse=True)
    
    contigs = []
    used_kmers = set()
    
    for start_kmer in sorted_kmers[:100]:  # Try top 100 k-mers as starts
        if start_kmer in used_kmers:
            continue
            
        contig = start_kmer
        current_kmer = start_kmer
        used_kmers.add(start_kmer)
        
        # Extend to the right
        while True:
            found_extension = False
            best_next = None
            best_count = 0
            
            # Try all possible next k-mers
            prefix = current_kmer[1:]  # Remove first base
            for base in ['A', 'C', 'G', 'T']:
                next_kmer = prefix + base
                if next_kmer in kmer_dict and next_kmer not in used_kmers:
                    if kmer_dict[next_kmer] > best_count:
                        best_count = kmer_dict[next_kmer]
                        best_next = next_kmer
                        found_extension = True
            
            if found_extension and best_count >= 2:  # Minimum coverage threshold
                contig += best_next[-1]  # Add last base
                used_kmers.add(best_next)
                current_kmer = best_next
            else:
                break
        
        # Only keep contigs longer than k+10
        if len(contig) > k + 10:
            contigs.append(contig)
    
    # Sort contigs by length
    contigs.sort(key=len, reverse=True)
    
    # Output contigs
    output_path = join_path(data_path, 'contig_codon.fasta')
    with open(output_path, 'w') as f:
        for i, contig in enumerate(contigs[:20]):  # Output top 20 contigs
            f.write(f'>contig_{i}\n')
            f.write(f'{contig}\n')
            print(f"Contig {i}: length {len(contig)}")
    
    print(f"Wrote {len(contigs[:20])} contigs to {output_path}")

if __name__ == "__main__":
    main()