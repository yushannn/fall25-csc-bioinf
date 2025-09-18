import sys
import os
from python import sys as pysys
from python import os as pyos
from python import copy
from python import dbg

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

def main():
    if len(sys.argv) < 2:
        print("Usage: codon run main_codon.py <data_directory>")
        return
    
    # Read data
    data_path = join_path('./', sys.argv[1])
    short1, short2, long1 = read_data(data_path)
    
    print(f"Read {len(short1)} short1 sequences")
    print(f"Read {len(short2)} short2 sequences")  
    print(f"Read {len(long1)} long sequences")
    
    # Build DBG
    k = 25
    print(f"Building DBG with k={k}")
    graph = dbg.DBG(k=k, data_list=[short1, short2, long1])
    
    # Output contigs
    output_path = join_path(data_path, 'contig_codon.fasta')
    with open(output_path, 'w') as f:
        for i in range(20):
            contig = graph.get_longest_contig()
            if contig is None:
                print(f"No more contigs found after {i} iterations")
                break
            print(f"Contig {i}: length {len(contig)}")
            f.write(f'>contig_{i}\n')
            f.write(f'{contig}\n')
    
    print(f"Results written to {output_path}")

if __name__ == "__main__":
    main()