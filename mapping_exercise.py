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

