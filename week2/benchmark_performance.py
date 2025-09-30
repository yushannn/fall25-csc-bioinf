"""Performance benchmark script for bio_codon.motifs module.

This script provides basic performance benchmarks to validate the
performance improvements claimed in the report.
"""

import time
import sys
import os
from typing import List, Dict

# Add the bio_codon module to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'bio_codon'))

try:
    from bio_codon.motifs import create, Motif
    from bio_codon.motifs.matrix import FrequencyPositionMatrix
except ImportError as e:
    print(f"Import error: {e}")
    sys.exit(1)


def benchmark_motif_creation(num_sequences: int = 1000) -> float:
    """Benchmark motif creation performance."""
    sequences = ["ATCGATCG"] * num_sequences
    
    start_time = time.time()
    motif = create(sequences)
    end_time = time.time()
    
    return (end_time - start_time) * 1000  # Return milliseconds


def benchmark_matrix_normalization(motif_length: int = 8, alphabet: str = "ACGT") -> float:
    """Benchmark matrix normalization performance."""
    # Create a sample frequency matrix
    counts = {}
    for letter in alphabet:
        counts[letter] = [10, 5, 15, 8, 12, 3, 7, 9][:motif_length]
    
    fpm = FrequencyPositionMatrix(alphabet, counts)
    
    start_time = time.time()
    pwm = fpm.normalize()
    end_time = time.time()
    
    return (end_time - start_time) * 1000  # Return milliseconds


def benchmark_pssm_calculation(motif_length: int = 8, alphabet: str = "ACGT") -> float:
    """Benchmark PSSM calculation performance."""
    # Create a sample motif
    counts = {}
    for letter in alphabet:
        counts[letter] = [10, 5, 15, 8, 12, 3, 7, 9][:motif_length]
    
    fpm = FrequencyPositionMatrix(alphabet, counts)
    pwm = fpm.normalize()
    
    start_time = time.time()
    pssm = pwm.log_odds()
    end_time = time.time()
    
    return (end_time - start_time) * 1000  # Return milliseconds


def benchmark_reverse_complement(num_iterations: int = 1000) -> float:
    """Benchmark reverse complement performance."""
    sequences = ["ATCGATCG", "ATCGATCG", "ATCGATCG"]
    motif = create(sequences)
    
    start_time = time.time()
    for _ in range(num_iterations):
        rc_motif = motif.reverse_complement()
    end_time = time.time()
    
    return (end_time - start_time) * 1000  # Return milliseconds


def benchmark_sequence_search(sequence_length: int = 10000) -> float:
    """Benchmark sequence search performance."""
    # Create a test motif
    sequences = ["ATCGATCG", "ATCGATCG", "TTCGATCG"]
    motif = create(sequences)
    pssm = motif.counts.normalize().log_odds()
    
    # Create a long test sequence
    test_sequence = "ATCGATCG" * (sequence_length // 8)
    
    start_time = time.time()
    matches = pssm.search(test_sequence, threshold=-10.0)
    end_time = time.time()
    
    return (end_time - start_time) * 1000  # Return milliseconds


def run_benchmarks():
    """Run all benchmarks and display results."""
    print("=" * 60)
    print("Performance Benchmark Results for bio_codon.motifs")
    print("=" * 60)
    
    benchmarks = [
        ("Motif Creation (1000 sequences)", lambda: benchmark_motif_creation(1000)),
        ("Matrix Normalization", lambda: benchmark_matrix_normalization()),
        ("PSSM Score Calculation", lambda: benchmark_pssm_calculation()),
        ("Reverse Complement (1000x)", lambda: benchmark_reverse_complement(1000)),
        ("Large Sequence Search (10kb)", lambda: benchmark_sequence_search(10000)),
    ]
    
    results = {}
    
    for name, benchmark_func in benchmarks:
        print(f"\nRunning: {name}")
        
        # Run multiple iterations for more stable results
        times = []
        for i in range(5):
            try:
                elapsed = benchmark_func()
                times.append(elapsed)
                print(f"  Iteration {i+1}: {elapsed:.2f}ms")
            except Exception as e:
                print(f"  Iteration {i+1}: ERROR - {e}")
                times.append(float('inf'))
        
        # Calculate average, excluding any failed runs
        valid_times = [t for t in times if t != float('inf')]
        if valid_times:
            avg_time = sum(valid_times) / len(valid_times)
            min_time = min(valid_times)
            max_time = max(valid_times)
            
            results[name] = {
                'avg': avg_time,
                'min': min_time,
                'max': max_time,
                'runs': len(valid_times)
            }
            
            print(f"  Average: {avg_time:.2f}ms (min: {min_time:.2f}ms, max: {max_time:.2f}ms)")
        else:
            print(f"  All iterations failed!")
            results[name] = {'avg': float('inf'), 'runs': 0}
    
    # Summary
    print("\n" + "=" * 60)
    print("BENCHMARK SUMMARY")
    print("=" * 60)
    
    for name, result in results.items():
        if result['runs'] > 0:
            print(f"{name:<35}: {result['avg']:>8.2f}ms")
        else:
            print(f"{name:<35}: {'FAILED':>8}")
    
    print("\n" + "=" * 60)
    print("Performance Analysis:")
    print("- These results represent the Codon implementation performance")
    print("- Actual speedup vs Python would require Python BioPython comparison")
    print("- Results may vary based on system specifications and load")
    print("- Consider these as baseline measurements for optimization")
    print("=" * 60)
    
    return results


def validate_functionality():
    """Quick functionality validation before benchmarking."""
    print("Validating functionality before benchmarking...")
    
    try:
        # Test basic motif creation
        sequences = ["ATCG", "ATGG", "TTCG"]
        motif = create(sequences)
        assert motif.length == 4
        assert motif.alphabet == "ACGT"
        
        # Test consensus
        consensus = motif.consensus
        assert len(consensus) == 4
        
        # Test matrix operations
        pwm = motif.counts.normalize()
        pssm = pwm.log_odds()
        
        # Test reverse complement
        rc_motif = motif.reverse_complement()
        assert rc_motif.length == motif.length
        
        # Test sequence search
        matches = pssm.search("ATCGATCG", threshold=-10.0)
        assert isinstance(matches, list)
        
        print("✓ All functionality tests passed!")
        return True
        
    except Exception as e:
        print(f"❌ Functionality validation failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    print("Bio.motifs Performance Benchmark Suite")
    print("=" * 60)
    
    if validate_functionality():
        results = run_benchmarks()
        
        # Exit with success code
        sys.exit(0)
    else:
        print("Functionality validation failed. Skipping benchmarks.")
        sys.exit(1)