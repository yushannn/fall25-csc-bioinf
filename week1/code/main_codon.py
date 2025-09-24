import sys
import os

def reverse_complement(key):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    key_list = list(key[::-1])
    for i in range(len(key_list)):
        key_list[i] = complement[key_list[i]]
    return ''.join(key_list)

class Node:
    def __init__(self, kmer):
        self._children = set[int]()
        self._count = 0
        self.kmer = kmer
        self.visited = False
        self.depth = 0
        self.max_depth_child = -1

    def add_child(self, kmer_idx):
        self._children.add(kmer_idx)

    def increase(self):
        self._count += 1

    def reset(self):
        self.visited = False
        self.depth = 0
        self.max_depth_child = -1

    def get_count(self):
        return self._count

    def get_children(self):
        return list(self._children)

    def remove_children(self, target):
        self._children = self._children - target

class DBG:
    def __init__(self, k, data_list):
        self.k = k
        self.nodes = dict[int, Node]()
        self.kmer2idx = dict[str, int]()
        self.kmer_count = 0
        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list):
        assert len(data_list) > 0
        assert self.k <= len(data_list[0][0])

    def _build(self, data_list):
        for data in data_list:
            for original in data:
                rc = reverse_complement(original)
                for i in range(len(original) - self.k):
                    self._add_arc(original[i: i + self.k], original[i + 1: i + 1 + self.k])
                    self._add_arc(rc[i: i + self.k], rc[i + 1: i + 1 + self.k])

    def _add_node(self, kmer):
        if kmer not in self.kmer2idx:
            self.kmer2idx[kmer] = self.kmer_count
            self.nodes[self.kmer_count] = Node(kmer)
            self.kmer_count += 1
        idx = self.kmer2idx[kmer]
        self.nodes[idx].increase()
        return idx

    def _add_arc(self, kmer1, kmer2):
        idx1 = self._add_node(kmer1)
        idx2 = self._add_node(kmer2)
        self.nodes[idx1].add_child(idx2)

    def _get_count(self, child):
        return self.nodes[child].get_count()

    def _get_sorted_children(self, idx):
        children = self.nodes[idx].get_children()
        children.sort(key=self._get_count, reverse=True)
        return children

    def _get_depth(self, idx):
        if not self.nodes[idx].visited:
            self.nodes[idx].visited = True
            children = self._get_sorted_children(idx)
            max_depth, max_child = 0, -1
            for child in children:
                depth = self._get_depth(child)
                if depth > max_depth:
                    max_depth, max_child = depth, child
            self.nodes[idx].depth = max_depth + 1
            self.nodes[idx].max_depth_child = max_child
        return self.nodes[idx].depth

    def _reset(self):
        for idx in self.nodes.keys():
            self.nodes[idx].reset()

    def _get_longest_path(self):
        max_depth, max_idx = 0, -1
        for idx in self.nodes.keys():
            depth = self._get_depth(idx)
            if depth > max_depth:
                max_depth, max_idx = depth, idx

        path = []
        current_idx = max_idx
        while current_idx >= 0:
            path.append(current_idx)
            current_idx = self.nodes[current_idx].max_depth_child
        return path

    def _delete_path(self, path):
        path_set = set(path)
        for idx in path:
            if idx in self.nodes:
                del self.nodes[idx]
        for idx in self.nodes.keys():
            self.nodes[idx].remove_children(path_set)

    def _concat_path(self, path):
        if len(path) < 1:
            return ""
        concat = self.nodes[path[0]].kmer
        for i in range(1, len(path)):
            concat += self.nodes[path[i]].kmer[-1]
        return concat

    def get_longest_contig(self):
        self._reset()
        path = self._get_longest_path()
        contig = self._concat_path(path)
        self._delete_path(path)
        return contig

def read_fasta(path, name):
    data = []
    file_path = path + '/' + name
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) > 0 and line[0] != '>':
                data.append(line)
    print(f"{name} {len(data)} {len(data[0]) if data else 0}")
    return data

def read_data(path):
    short1 = read_fasta(path, "short_1.fasta")
    short2 = read_fasta(path, "short_2.fasta")
    long1 = read_fasta(path, "long.fasta")
    return short1, short2, long1

def main():
    argv = sys.argv
    if len(argv) < 2:
        print("Usage: codon run main_codon.py <dataset_path>")
        sys.exit(1)
    
    dataset_path = './' + argv[1]
    short1, short2, long1 = read_data(dataset_path)

    k = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])
    
    output_path = dataset_path + '/contig_codon.fasta'
    with open(output_path, 'w') as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c == "" or len(c) == 0:
                break
            print(f"{i} {len(c)}")
            f.write(f'>contig_{i}\n')
            f.write(c + '\n')

if __name__ == "__main__":
    main()