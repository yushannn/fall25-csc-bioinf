def reverse_complement(seq: str) -> str:
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

def build_graph(sequences: list[list[str]], k: int) -> tuple[dict[str, list[str]], dict[str, int]]:
    """Build a de Bruijn graph from sequences"""
    graph = {}  # kmer -> list of next kmers
    counts = {}  # kmer -> count
    
    for seq_list in sequences:
        for seq in seq_list:
            if len(seq) < k + 1:
                continue
                
            # Process forward sequence
            for i in range(len(seq) - k):
                kmer1 = seq[i:i+k]
                kmer2 = seq[i+1:i+1+k]
                
                # Add nodes
                if kmer1 not in graph:
                    graph[kmer1] = []
                    counts[kmer1] = 0
                if kmer2 not in graph:
                    graph[kmer2] = []
                    counts[kmer2] = 0
                    
                # Add edge
                if kmer2 not in graph[kmer1]:
                    graph[kmer1].append(kmer2)
                    
                # Increment counts
                counts[kmer1] += 1
                counts[kmer2] += 1
                
            # Process reverse complement
            rc_seq = reverse_complement(seq)
            for i in range(len(rc_seq) - k):
                kmer1 = rc_seq[i:i+k]
                kmer2 = rc_seq[i+1:i+1+k]
                
                # Add nodes
                if kmer1 not in graph:
                    graph[kmer1] = []
                    counts[kmer1] = 0
                if kmer2 not in graph:
                    graph[kmer2] = []
                    counts[kmer2] = 0
                    
                # Add edge
                if kmer2 not in graph[kmer1]:
                    graph[kmer1].append(kmer2)
                    
                # Increment counts
                counts[kmer1] += 1
                counts[kmer2] += 1
                
    return graph, counts

def find_longest_path(graph: dict[str, list[str]], counts: dict[str, int]) -> list[str]:
    """Find the longest path in the graph"""
    # Track visited nodes and their depths
    visited = {}
    depths = {}
    max_children = {}
    
    # DFS to find longest path
    def dfs(node: str) -> int:
        if node in visited:
            return depths[node]
            
        visited[node] = True
        max_depth = 0
        max_child = None
        
        # Sort children by count
        if node in graph:
            children = sorted(graph[node], key=lambda x: counts[x] if x in counts else 0, reverse=True)
        else:
            children = []
        
        for child in children:
            depth = dfs(child)
            if depth > max_depth:
                max_depth = depth
                max_child = child
                
        depths[node] = max_depth + 1
        max_children[node] = max_child
        return depths[node]
    
    # Find node with maximum depth
    max_depth = 0
    max_node = None
    
    for node in graph:
        visited = {}  # Reset visited for each starting node
        depth = dfs(node)
        if depth > max_depth:
            max_depth = depth
            max_node = node
    
    # Reconstruct path
    path = []
    current = max_node
    while current is not None:
        path.append(current)
        if current in max_children:
            current = max_children[current]
        else:
            current = None
        
    return path

def concat_path(path: list[str]) -> str:
    """Concatenate kmers in path to form contig"""
    if not path:
        return ""
        
    result = path[0]
    for i in range(1, len(path)):
        result += path[i][-1]
        
    return result

def remove_path(graph: dict[str, list[str]], counts: dict[str, int], path: list[str]) -> None:
    """Remove a path from the graph"""
    path_set = set(path)
    
    # Remove nodes in path
    for node in path:
        if node in graph:
            del graph[node]
        if node in counts:
            del counts[node]
    
    # Remove edges to deleted nodes
    for node in list(graph.keys()):
        graph[node] = [child for child in graph[node] if child not in path_set]

def get_longest_contig(sequences: list[list[str]], k: int) -> str:
    """Get the longest contig from sequences using de Bruijn graph"""
    graph, counts = build_graph(sequences, k)
    path = find_longest_path(graph, counts)
    contig = concat_path(path)
    remove_path(graph, counts, path)
    return contig

def assemble_contigs(sequences: list[list[str]], k: int, max_contigs: int = 20) -> list[str]:
    """Assemble contigs from sequences"""
    graph, counts = build_graph(sequences, k)
    contigs = []
    
    for i in range(max_contigs):
        if len(graph) == 0:  # Graph is empty
            break
            
        path = find_longest_path(graph, counts)
        if len(path) == 0:  # No path found
            break
            
        contig = concat_path(path)
        contigs.append(contig)
        remove_path(graph, counts, path)
        
    return contigs