import os
from typing import List, Tuple

def join_path(path: str, name: str) -> str:
    """Join path components"""
    if path.endswith('/'):
        return path + name
    else:
        return path + '/' + name

def read_fasta(path: str, name: str) -> List[str]:
    """Read sequences from a FASTA file"""
    data: List[str] = []
    full_path = join_path(path, name)
    
    with open(full_path, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) > 0 and line[0] != '>':
                data.append(line)
                
    print(f"{name}: {len(data)} sequences, length {len(data[0]) if data else 0}")
    return data

def read_data(path: str) -> Tuple[List[str], List[str], List[str]]:
    """Read all data files from the specified path"""
    short1 = read_fasta(path, "short_1.fasta")
    short2 = read_fasta(path, "short_2.fasta")
    long1 = read_fasta(path, "long.fasta")
    return short1, short2, long1