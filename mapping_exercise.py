#____________      TASK 1      ______________________

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


# Test TASK1
text = "AAGGCGGCTTCA$"
p = "GGC"

positions = brute_force(text, p)
sa = suffix_array(text)
a, b = find_block(text, sa, p)

print("_________ TEST TASK 1 ____________")
print("Brute force:", positions)
print("Suffix array:", sa)
print("Block:", (a, b))
print("Positions from SA:", sorted(sa[a:b]))

assert sorted(sa) == list(range(len(text)))

for i in range(len(sa) - 1):
    assert text[sa[i]:] < text[sa[i + 1]:]

assert sorted(sa[a:b]) == positions
