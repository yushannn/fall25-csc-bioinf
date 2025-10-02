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

test_motif_creation()
test_frequency_matrix()
test_weight_matrix()
test_scoring_matrix()
test_consensus()

print("All tests passed")
