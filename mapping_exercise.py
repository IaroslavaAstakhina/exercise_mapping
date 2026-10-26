#____________      TASK 1      ______________________
import random
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


# TESTY - Handout Example

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

# Task 1 self-checks
assert sorted(sa) == list(range(len(text)))

for i in range(len(sa) - 1):
    assert text[sa[i]:] < text[sa[i + 1]:]

assert sorted(sa[a:b]) == positions

#____________      TASK 2 TEST      ______________________

print("\n________ TASK 2 _______")

transformed = bwt(text, sa)
first = first_occurrence(transformed)
reconstructed = inverse_bwt(transformed)

print("Original text:", text)
print("BWT:", transformed)
print("First occurrence:", first)
print("Inverse BWT:", reconstructed)

# Task 2 self-checks
assert sorted(transformed) == sorted(text)
assert transformed.count("$") == 1
assert first["$"] == 0
assert reconstructed == text


# -------------------RANDOM TESTS--------------------------------------------------------------


print("\n________ RANDOM TESTS _______")

for test in range(300):
    length = random.randint(1, 30)
    random_text = "".join(
        random.choices("ACGT", k=length)
    ) + "$"

    random_sa = suffix_array(random_text)

    # Test Task 1
    assert sorted(random_sa) == list(range(len(random_text)))

    for i in range(len(random_sa) - 1):
        assert random_text[random_sa[i]:] < random_text[random_sa[i + 1]:]

    for j in range(5):
        pattern = "".join(
            random.choices("ACGT", k=random.randint(1, 5))
        )

        a, b = find_block(random_text, random_sa, pattern)

        assert sorted(random_sa[a:b]) == brute_force(
            random_text, pattern
        )

    # Test Task 2
    random_bwt = bwt(random_text, random_sa)

    assert sorted(random_bwt) == sorted(random_text)
    assert random_bwt.count("$") == 1
    assert inverse_bwt(random_bwt) == random_text

print(" random tests - done")
