def reverse_complement(key: str) -> str:
    """Calculate reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    result = []
    for i in range(len(key) - 1, -1, -1):
        result.append(complement.get(key[i], 'N'))
    return ''.join(result)

class DBG:
    """De Bruijn Graph implementation"""
    def __init__(self, k: int, data_list):
        self.k = k
        self.nodes = {}
        
        # Validate input
        assert len(data_list) > 0
        if data_list[0]:
            assert self.k <= len(data_list[0][0])
            
        # Build the graph
        for data in data_list:
            for original in data:
                if len(original) < self.k + 1:
                    continue
                    
                rc = reverse_complement(original)
                
                # Process original sequence
                orig_len = len(original)
                for i in range(orig_len - self.k):
                    kmer1 = original[i:i + self.k]
                    kmer2 = original[i + 1:i + 1 + self.k]
                    self._add_arc(kmer1, kmer2)
                    
                    # Process reverse complement
                    rc_kmer1 = rc[i:i + self.k]
                    rc_kmer2 = rc[i + 1:i + 1 + self.k]
                    self._add_arc(rc_kmer1, rc_kmer2)

    def _add_node(self, kmer):
        if kmer not in self.nodes:
            self.nodes[kmer] = {
                'children': [],
                'count': 0,
                'visited': False,
                'depth': 0,
                'max_child': None
            }
        self.nodes[kmer]['count'] += 1

    def _add_arc(self, kmer1, kmer2):
        self._add_node(kmer1)
        self._add_node(kmer2)
        if kmer2 not in self.nodes[kmer1]['children']:
            self.nodes[kmer1]['children'].append(kmer2)

    def _get_count(self, child):
        return self.nodes[child]['count']

    def _get_sorted_children(self, kmer):
        if kmer not in self.nodes:
            return []
        children = self.nodes[kmer]['children']
        children.sort(key=self._get_count, reverse=True)
        return children

    def _get_depth(self, kmer):
        if not self.nodes[kmer]['visited']:
            self.nodes[kmer]['visited'] = True
            children = self._get_sorted_children(kmer)
            max_depth = 0
            max_child = None
            
            for child in children:
                depth = self._get_depth(child)
                if depth > max_depth:
                    max_depth = depth
                    max_child = child
                    
            self.nodes[kmer]['depth'] = max_depth + 1
            self.nodes[kmer]['max_child'] = max_child
        
        return self.nodes[kmer]['depth']

    def _reset(self):
        for kmer in self.nodes:
            self.nodes[kmer]['visited'] = False
            self.nodes[kmer]['depth'] = 0
            self.nodes[kmer]['max_child'] = None

    def _get_longest_path(self):
        max_depth = 0
        max_kmer = None
        
        for kmer in self.nodes:
            depth = self._get_depth(kmer)
            if depth > max_depth:
                max_depth = depth
                max_kmer = kmer

        path = []
        current_kmer = max_kmer
        while current_kmer is not None:
            path.append(current_kmer)
            current_kmer = self.nodes[current_kmer]['max_child']
            
        return path

    def _delete_path(self, path):
        # Delete nodes
        for kmer in path:
            if kmer in self.nodes:
                del self.nodes[kmer]
        
        # Remove edges pointing to deleted nodes
        path_set = set(path)
        for kmer in list(self.nodes.keys()):
            self.nodes[kmer]['children'] = [child for child in self.nodes[kmer]['children'] if child not in path_set]

    def _concat_path(self, path):
        if len(path) < 1:
            return None
            
        # Get first kmer
        concat = path[0]
        
        # Append last character of subsequent kmers
        for i in range(1, len(path)):
            if path[i]:
                concat += path[i][-1]
                    
        return concat

    def get_longest_contig(self):
        if len(self.nodes) == 0:
            return None
            
        self._reset()
        path = self._get_longest_path()
        contig = self._concat_path(path)
        self._delete_path(path)
        return contig