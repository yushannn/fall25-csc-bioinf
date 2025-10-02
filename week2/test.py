# File test.py
try:
    __codon__
    from bio_codon import motifs
except NameError:
    from Bio import motifs

def test(func):
    func.is_test = True
    return func

def assertEqual(a, b, msg=""):
    if a != b:
        raise AssertionError(f"Expected {a} == {b}. {msg}")

def assertTrue(condition, msg=""):
    if not condition:
        raise AssertionError(f"Expected True. {msg}")

@test
def test_motif_creation():
    instances = ["TACAA", "TACGC", "TACAC", "TACCC"]
    motif = motifs.create(instances)
    assertEqual(motif.length, 5)
    assertTrue("A" in motif.alphabet)

@test
def test_frequency_matrix():
    instances = ["TACAA", "TACGC", "TACAC", "TACCC"]
    motif = motifs.create(instances)
    fpm = motif.counts
    assertEqual(fpm.length, 5)

@test
def test_weight_matrix():
    instances = ["TACAA", "TACGC", "TACAC", "TACCC"]
    motif = motifs.create(instances)
    pwm = motif.counts.normalize()
    assertEqual(pwm.length, 5)

@test
def test_scoring_matrix():
    instances = ["TACAA", "TACGC", "TACAC", "TACCC"]
    motif = motifs.create(instances)
    pssm = motif.counts.normalize().log_odds()
    assertEqual(pssm.length, 5)

@test
def test_consensus():
    instances = ["TACAA", "TACGC", "TACAC", "TACCC"]
    motif = motifs.create(instances)
    consensus = motif.consensus
    assertEqual(len(consensus), 5)

@test
def test_reverse_complement():
    instances = ["ATCG", "ATGG"]
    motif = motifs.create(instances)
    rc = motif.reverse_complement()
    assertEqual(rc.length, 4)

@test
def test_pssm_scoring():
    instances = ["TACAA", "TACGC"]
    motif = motifs.create(instances)
    pssm = motif.pssm
    score = pssm.calculate("TACAA")
    assertTrue(score > 0)

@test
def test_search_functionality():
    instances = ["ATG"] * 3
    motif = motifs.create(instances)
    matches = list(motif.pssm.search("ATGATGATG"))
    assertTrue(len(matches) > 0)

@test
def test_single_sequence():
    instances = ["ATCG"]
    motif = motifs.create(instances)
    assertEqual(motif.length, 4)
    assertTrue(hasattr(motif, 'counts'))

@test
def test_short_sequences():
    instances = ["A", "T", "G"]
    motif = motifs.create(instances)
    assertEqual(motif.length, 1)

@test
def test_long_sequences():
    long_seq = "ATCGATCGATCGATCG"
    instances = [long_seq] * 2
    motif = motifs.create(instances)
    assertEqual(motif.length, 16)

def test_empty_input():
    try:
        motifs.create([])
        assertTrue(False, "Should have raised an error")
    except (ValueError, IndexError):
        pass

def test_invalid_characters():
    # Test behavior with invalid characters - implementation may vary
    try:
        motif = motifs.create(["ATCGX"])
        # If no error is raised, verify motif was still created
        assertTrue(hasattr(motif, 'length'))
        assertTrue(motif.length == 5)
    except (ValueError, KeyError):
        # Error handling is also acceptable
        pass

@test
def test_mixed_length_sequences():
    instances = ["ATCG", "AT"]
    try:
        motif = motifs.create(instances)
        assertTrue(False, "Should have raised an error for different lengths")
    except ValueError:
        pass

@test
def test_type_errors():
    # Test string input behavior - implementation may vary
    try:
        result = motifs.create("ATCG")
        # Some implementations may handle this differently
        assertTrue(hasattr(result, 'length') or result is None)
    except (TypeError, AttributeError):
        # Error handling is acceptable
        pass
    
    # Test numeric input
    try:
        motifs.create([123])
        assertTrue(False, "Should reject non-string elements")
    except (TypeError, AttributeError):
        pass

@test 
def test_performance_large_dataset():
    # Test with larger dataset to verify basic performance
    many_instances = ["ATCGATCG"] * 50
    motif = motifs.create(many_instances)
    assertEqual(motif.length, 8)
    pssm = motif.pssm
    assertTrue(hasattr(pssm, 'calculate'))

@test
def test_integration_motif_to_matrices():
    instances = ["TACAA", "TACGC", "TACAC"]
    motif = motifs.create(instances)
    
    # Test the flow from motif -> counts -> pwm -> pssm
    counts = motif.counts
    assertEqual(counts.length, 5)
    
    pwm = counts.normalize()
    assertEqual(pwm.length, 5)
    
    pssm = pwm.log_odds()
    assertEqual(pssm.length, 5)
    
    # Test scoring with original sequence
    score = pssm.calculate("TACAA")
    # Score should be a numeric value (may be different types in different implementations)
    assertTrue(score is not None)

# Run all tests
tests = [
    test_motif_creation,
    test_frequency_matrix,
    test_weight_matrix,
    test_scoring_matrix,
    test_consensus,
    test_reverse_complement,
    test_pssm_scoring,
    test_search_functionality,
    test_single_sequence,
    test_short_sequences,
    test_long_sequences,
    test_empty_input,
    test_invalid_characters,
    test_mixed_length_sequences,
    test_type_errors,
    test_performance_large_dataset,
    test_integration_motif_to_matrices
]

for test_func in tests:
    test_func()

print(f"All {len(tests)} tests passed")
