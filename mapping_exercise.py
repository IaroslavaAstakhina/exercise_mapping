
#____________      TASK 1      ______________________
import random
from Bio import SeqIO
from pathlib import Path
import time

def brute_force(text, p):
    positions = []

    for i in range(len(text) - len(p) + 1):
        if text[i:i + len(p)] == p:
            positions.append(i)

    return positions


def suffix_array(text):
    return sorted(range(len(text)), key=lambda i: text[i:])


def find_block(text, sa, p):
    n = len(sa)
    m = len(p)
    lo = 0
    hi = n

    while lo < hi:
        mid = (lo + hi) // 2
        prefix = text[sa[mid]:sa[mid] + m]

        if prefix < p:
            lo = mid + 1
        else:
            hi = mid

    a = lo
    lo = 0
    hi = n

    while lo < hi:
        mid = (lo + hi) // 2
        prefix = text[sa[mid]:sa[mid] + m]

        if prefix <= p:
            lo = mid + 1
        else:
            hi = mid

    b = lo

    return a, b

def bwt(text, sa):
    result = ""

    for i in sa:
        result += text[i - 1]

    return result


def first_occurrence(letters):
    first = {}
    sorted_letters = sorted(letters)

    for i in range(len(sorted_letters)):
        letter = sorted_letters[i]

        if letter not in first:
            first[letter] = i

    return first


def inverse_bwt(bwt):
    first = first_occurrence(bwt)

    # Count occurrences of each letter in L
    counts = {}
    k = []

    for letter in bwt:
        counts[letter] = counts.get(letter, 0) + 1
        k.append(counts[letter])

    # Reconstruct the text backwards
    row = 0
    result = ""

    for i in range(len(bwt) - 1):
        letter = bwt[row]
        result += letter

        row = first[letter] + k[row] - 1

    return result[::-1] + "$"


#____________      TASK 3 - FMIndex      ______________________

class FMIndex:

    def __init__(self, text, s):
        self.s = s

        # Build suffix array and BWT
        sa = suffix_array(text)
        self.bwt = bwt(text, sa)

        # First occurrence of each letter
        self.first_occurrence = first_occurrence(self.bwt)

        # Count table for each letter
        self.counts = {}

        for letter in set(self.bwt):
            self.counts[letter] = [0] * (len(text) + 1)

        for i in range(len(self.bwt)):
            for letter in self.counts:
                self.counts[letter][i + 1] = self.counts[letter][i]

            current = self.bwt[i]
            self.counts[current][i + 1] += 1

        # Store only sampled SA values
        self.sampled = {}

        for row in range(len(sa)):
            if sa[row] % s == 0:
                self.sampled[row] = sa[row]


    def lf(self, row):
        letter = self.bwt[row]

        return (self.first_occurrence[letter]
                + self.counts[letter][row])


    def count(self, p):
        top = 0
        bottom = len(self.bwt)

        # Backward search
        for letter in reversed(p):

            if letter not in self.first_occurrence:
                return (0, 0)

            top_new = (self.first_occurrence[letter]
                       + self.counts[letter][top])

            bottom_new = (self.first_occurrence[letter]
                          + self.counts[letter][bottom])

            top = top_new
            bottom = bottom_new

            if top == bottom:
                break

        return (top, bottom)


    def locate(self, p):
        top, bottom = self.count(p)
        positions = []

        for row in range(top, bottom):
            current_row = row
            steps = 0

            # Follow LF until a sampled SA value is found
            while current_row not in self.sampled:
                current_row = self.lf(current_row)
                steps += 1

            position = self.sampled[current_row] + steps
            positions.append(position)

        return positions


# ============================================================
# TASK 4 - MAPPING READS
# ============================================================

def reverse_complement(seq):
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    result = ""

    for letter in reversed(seq):
        result += complement[letter]

    return result


def hamming(a, b):
    assert len(a) == len(b)

    mismatches = 0

    for x, y in zip(a, b):
        if x != y:
            mismatches += 1

    return mismatches


def seeds(read, n):
    result = []

    for i in range(n):
        start = i * len(read) // n
        end = (i + 1) * len(read) // n

        result.append((start, read[start:end]))

    return result


def map_read(read, fm, ref, d=2, n_seeds=None):
    if n_seeds is None:
        n_seeds = d + 1

    hits = []
    total_candidates = 0

    for strand, sequence in [
        ("+", read),
        ("-", reverse_complement(read))
    ]:

        candidates = set()

        # Find candidates from exact seed matches
        for offset, seed in seeds(sequence, n_seeds):

            for p in fm.locate(seed):
                start = p - offset

                if 0 <= start <= len(ref) - len(sequence):
                    candidates.add(start)

        # Verify each unique candidate
        total_candidates += len(candidates)

        for start in candidates:
            fragment = ref[start:start + len(sequence)]
            mismatches = hamming(sequence, fragment)

            if mismatches <= d:
                hits.append((strand, start, mismatches))

    return sorted(hits), total_candidates


def brute_map(read, ref, d=2):
    hits = []

    for strand, sequence in [
        ("+", read),
        ("-", reverse_complement(read))
    ]:

        for start in range(len(ref) - len(sequence) + 1):
            fragment = ref[start:start + len(sequence)]

            mismatches = hamming(sequence, fragment)

            if mismatches <= d:
                hits.append((strand, start, mismatches))

    return sorted(hits)


# ============================================================
# TASK 1 - TESTS
# ============================================================

text = "AAGGCGGCTTCA$"
p = "GGC"

print("________ TASK 1 _______")

positions = brute_force(text, p)
sa = suffix_array(text)
a, b = find_block(text, sa, p)

print("Text:", text)
print("Pattern:", p)
print("Brute force:", positions)
print("Suffix array:", sa)
print("Block:", (a, b))
print("Positions from SA:", sorted(sa[a:b]))

# Self-checks
assert sorted(sa) == list(range(len(text)))

for i in range(len(sa) - 1):
    assert text[sa[i]:] < text[sa[i + 1]:]

assert sorted(sa[a:b]) == positions


# ============================================================
# TASK 2 - TESTS
# ============================================================

print("\n________ TASK 2 _______")

transformed = bwt(text, sa)
first = first_occurrence(transformed)
reconstructed = inverse_bwt(transformed)

print("Original text:", text)
print("BWT:", transformed)
print("First occurrence:", first)
print("Inverse BWT:", reconstructed)

# Self-checks
assert sorted(transformed) == sorted(text)
assert transformed.count("$") == 1
assert first["$"] == 0
assert reconstructed == text


# ============================================================
# TASK 3 - TESTS
# ============================================================

print("\n________ TASK 3 _______")

# Example from assignment: CAT$
fm_cat = FMIndex("CAT$", 2)

print("CAT$ example:")
print("BWT:", fm_cat.bwt)
print("Counts T:", fm_cat.counts["T"])
print("Sampled SA:", fm_cat.sampled)
print("LF(1):", fm_cat.lf(1))
print("Count AT:", fm_cat.count("AT"))
print("Locate AT:", fm_cat.locate("AT"))

assert fm_cat.bwt == "TC$A"
assert fm_cat.counts["T"] == [0, 1, 1, 1, 1]
assert fm_cat.sampled == {2: 0, 3: 2}
assert fm_cat.lf(1) == 2
assert fm_cat.count("AT") == (1, 2)
assert fm_cat.locate("AT") == [1]

# Handout example: GGC
fm = FMIndex(text, 4)

print("\nHandout example:")
print("BWT:", fm.bwt)
print("Count GGC:", fm.count("GGC"))
print("Locate GGC:", fm.locate("GGC"))

assert fm.count("GGC") == (9, 11)
assert sorted(fm.locate("GGC")) == [2, 5]


# ============================================================
# RANDOM TESTS - TASKS 1, 2 AND 3
# ============================================================

print("\n________ RANDOM TESTS _______")

for test in range(300):

    length = random.randint(1, 30)
    random_text = "".join(
        random.choices("ACGT", k=length)
    ) + "$"

    random_sa = suffix_array(random_text)

    # --------------------------------------------------------
    # Task 1 - Suffix Array and Binary Search
    # --------------------------------------------------------

    assert sorted(random_sa) == list(range(len(random_text)))

    for i in range(len(random_sa) - 1):
        assert random_text[random_sa[i]:] < random_text[random_sa[i + 1]:]

    # --------------------------------------------------------
    # Task 2 - BWT and Inverse BWT
    # --------------------------------------------------------

    random_bwt = bwt(random_text, random_sa)

    assert sorted(random_bwt) == sorted(random_text)
    assert random_bwt.count("$") == 1
    assert inverse_bwt(random_bwt) == random_text

    # --------------------------------------------------------
    # Task 3 - FMIndex (s = 1, 4, 32)
    # --------------------------------------------------------

    for s in [1, 4, 32]:

        index = FMIndex(random_text, s)

        # Check sampled SA values
        assert sorted(index.sampled.values()) == list(
            range(0, len(random_text), s)
        )

        for j in range(5):

            pattern = "".join(
                random.choices("ACGTN", k=random.randint(1, 5))
            )

            # Task 1 - Brute force and binary search
            expected = brute_force(random_text, pattern)
            a, b = find_block(random_text, random_sa, pattern)

            assert sorted(random_sa[a:b]) == expected

            # Task 3 - Backward search and locate
            top, bottom = index.count(pattern)
            positions = index.locate(pattern)

            assert bottom - top == len(positions)
            assert len(positions) == len(expected)
            assert sorted(positions) == expected

            # Compare intervals from Task 1 and Task 3
            if expected:
                assert (top, bottom) == (a, b)
            else:
                assert top == bottom
                assert a == b

print("random texts - done")


# ============================================================
# TASK 4 - FUNCTION TESTS
# ============================================================

print("\n________ TASK 4 - FUNCTION TESTS _______")

assert reverse_complement("AACG") == "CGTT"
assert reverse_complement(reverse_complement("AACG")) == "AACG"

assert hamming("ACGT", "AGGA") == 2

assert seeds("ACGTACGTAC", 3) == [
    (0, "ACG"),
    (3, "TAC"),
    (6, "GTAC")
]

test_ref = "AAGGCGGCTTCA"
test_fm = FMIndex(test_ref + "$", 4)

hits, candidates = map_read("GGC", test_fm, test_ref, d=0)

assert hits == brute_map("GGC", test_ref, d=0)
assert ("+", 2, 0) in hits
assert ("+", 5, 0) in hits

print("Reverse complement - OK")
print("Hamming - OK")
print("Seeds - OK")
print("Mapping - OK")


# ============================================================
# TASK 4 - MAP ALL READS
# ============================================================

folder = Path(__file__).resolve().parent

reference_record = SeqIO.read(folder / "reference.fasta", "fasta")
ref = str(reference_record.seq).upper()

read_records = list(SeqIO.parse(folder / "reads.fasta", "fasta"))

print("\n________ TASK 4 - MAPPING READS _______")
print("Reference length:", len(ref))
print("Number of reads:", len(read_records))

# Build the FMIndex only once
fm = FMIndex(ref + "$", 32)

results = []
total_candidates = 0

start_time = time.perf_counter()

for record in read_records:
    read = str(record.seq).upper()

    hits, candidates = map_read(read, fm, ref, d=2)

    results.append((record.id, read, hits))
    total_candidates += candidates

elapsed = time.perf_counter() - start_time

mapped = sum(1 for name, read, hits in results if hits)
unmapped = sum(1 for name, read, hits in results if not hits)
multiple = sum(1 for name, read, hits in results if len(hits) > 1)

print("\nMapping results:")
print("Mapped reads:", mapped)
print("Unmapped reads:", unmapped)
print("Reads with multiple hits:", multiple)
print("Total candidates:", total_candidates)
print("Average candidates per read:", total_candidates / len(results))
print("Mapping time:", round(elapsed, 3), "seconds")

print("\nFirst five reads:")

for name, read, hits in results[:5]:
    print("\nRead:", name)

    if hits:
        for strand, position, mismatches in hits:
            print("Strand:", strand,
                  "Position:", position,
                  "Mismatches:", mismatches)
    else:
        print("Unmapped")


# ============================================================
# TASK 4 - SELF-CHECKS
# ============================================================

print("\n________ TASK 4 - SELF-CHECKS _______")

# 1. Verify every reported hit
for name, read, hits in results:

    for strand, position, mismatches in hits:
        assert 0 <= position <= len(ref) - len(read)

        if strand == "+":
            sequence = read
        else:
            sequence = reverse_complement(read)

        fragment = ref[position:position + len(read)]

        assert hamming(sequence, fragment) == mismatches
        assert mismatches <= 2

print("All reported hits verified - OK")


# 2. Compare 30 random reads with brute force
sample = random.sample(results, min(30, len(results)))

for name, read, hits in sample:
    expected = brute_map(read, ref, d=2)
    assert hits == expected

print("30 random reads vs brute force - OK")


# 3. Check every unmapped read with brute force
for name, read, hits in results:
    if not hits:
        assert brute_map(read, ref, d=2) == []

print("All unmapped reads verified - OK")


# 4. Check the reference ends
first_read = ref[:100]
last_read_rc = reverse_complement(ref[-100:])

first_hits, _ = map_read(first_read, fm, ref)
last_hits, _ = map_read(last_read_rc, fm, ref)

assert ("+", 0, 0) in first_hits
assert ("-", len(ref) - 100, 0) in last_hits

print("Reference ends - OK")
print("All Task 4 self-checks passed!")


# ============================================================
# TASK 4 - THINK ABOUT IT: SEED EXPERIMENT
# ============================================================

print("\n________ TASK 4 - SEED EXPERIMENT _______")

for n_seeds in [3, 10, 12, 16]:

    experiment_candidates = 0
    all_hits = []

    start_time = time.perf_counter()

    for record in read_records:
        read = str(record.seq).upper()

        hits, candidates = map_read(
            read, fm, ref, d=2, n_seeds=n_seeds
        )

        all_hits.append(hits)
        experiment_candidates += candidates

    elapsed = time.perf_counter() - start_time

    same_hits = all(
        all_hits[i] == results[i][2]
        for i in range(len(results))
    )

    seed_lengths = sorted(set(
        len(seed) for offset, seed in seeds("A" * 100, n_seeds)
    ))

    print("\nNumber of seeds:", n_seeds)
    print("Seed lengths:", seed_lengths)
    print("Average candidates:",
          round(experiment_candidates / len(read_records), 2))
    print("Time:", round(elapsed, 3), "seconds")
    print("Same hits:", same_hits)

    assert same_hits

print("\nSeed experiment completed!")
