
def associated(x, y):
    return abs(x) == abs(y)


def divisor_cascade1(S):
    # Step 1
    k = 0
    t = [1]
    d = {t[0]}
    for si in S:
        if not any(associated(si, sj) for sj in d):
            k = k + 1
            t.append(si)
            d.add(si)

    # Skip to step 8?
    if len(t) > 2:

        # Step 2
        A = {(i, j) for j in range(1, k+1) for i in range(1, j)}
        B = set(A)

        while A or B:
            # Step 3 (go to step 7?)
            s = None
            if A:
                i, j = A.pop()

                # Step 4 (skip to step 6?)
                if t[i] % t[j] == 0:
                    s = t[i] // t[j]
                    # Go to step 5
                # Step 6:
                elif t[j] % t[i] == 0:
                    s = t[j] // t[i]
                    # Go to step 5
                else:
                    # Go to step 3
                    continue

            # Step 7 (go to step 8?)
            elif B:
                i, j = B.pop()
                s = gcd(t[i], t[j])
                # Go to step 5
            else:
                # Go to step 8
                break

            # Step 5
            if any(associated(s, sj) for sj in d):
                # Go to step 3
                continue
            k = k + 1
            t.append(s)
            d.add(s)
            A = A | {(i, k) for i in range(1, k-1)}
            B = B | {(i, k) for i in range(1, k-1)}
            # Go to step 3
            continue

    # Step 8
    d.update(t)
    return d


def divisor_cascade2(S):
    d = {1}
    queue = list(S)
    for si in queue:
        if any(associated(si, sj) for sj in d):
            continue
        for sj in d:
            if si % sj == 0:
                s = si // sj
            elif sj % si == 0:
                s = sj // si
            else:
                s = gcd(si, sj)
            if not associated(s, 1):
                queue.append(s)
        d.add(si)
    return d


from collections import defaultdict, Counter

def factorial_basis(S):
    d = set()
    factormap = {}
    queue = list(S)
    for si in queue:
        if any(associated(si, sj) for sj in d):
            continue
        for sj in d:
            if si % sj == 0:
                s = si // sj
                factormap[si] = [s, sj]
            elif sj % si == 0:
                s = sj // si
                factormap[sj] = [s, si]
            else:
                s = igcd(si, sj)
            if not associated(s, 1):
                queue.append(s)
        d.add(si)

    # The basis is the subset of the divisors that were not factored.
    basis = set(d) - factormap.keys()

    return basis, factormap


def factorial_basis2(S):
    basis = set()
    seen = set()
    factormap = {}
    queue = list(S)
    for si in queue:
        if associated(si, 1) or si in basis or si in factormap:
            continue
        for sj in basis:
            g = gcd(si, sj)
            if g != 1:
                if g == si:
                    ci = sj // si
                    factormap[sj] = [ci, si]
                    queue.extend([ci, si])
                    # Swap si for sj in the basis
                    basis.remove(sj)
                    basis.add(si)
                elif g == sj:
                    cj = si // sj
                    factormap[si] = [cj, sj]
                    queue.append(cj)
                else:
                    ci = si // g
                    cj = sj // g
                    basis.remove(sj)
                    factormap.update({si: [ci, g], sj: [cj, g]})
                    queue.extend([ci, si, g])
                # Exit the loop without adding si If it was found to divide
                # something then it has already been added
                break
        else:
            basis.add(si)

    return basis, factormap


assert divisor_cascade1([14700, 5040]) == {1, 12, 35, 144, 420, 1225, 5040, 14700}
assert divisor_cascade2([14700, 5040]) == {1, 12, 35, 144, 420, 1225, 5040, 14700}
