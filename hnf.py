from sympy import *
from sympy.utilities.iterables import iterable


def iround(num, den):
    # round away from 0 on tie
    if den < 0:
        num = -num
        den = -den
    neg = (num < 0)
    w, r = divmod(-num if neg else num, den)
    w = int(w)
    if 2*r >= den:
        w += 1
    return -w if neg else w


def idiv(a, b):
    # safe a/b asserting that b divides a evenly
    # if a%b != 0: print('non-int division of ',(a, b))
    return a//b


def _swap(k, cols, A, B, la, D):
    # swap rows k and k - 1 and calculate new la and D
    j = k - 1
    for M in (A, B):
        M[k], M[j] = M[j], M[k]
    for i in range(j):
        la[k][i], la[j][i] = la[j][i], la[k][i]
    for i in range(k + 1, cols):
        t = la[i][j]*D[k] - la[i][k]*la[k][j]
        la[i][j] = idiv(
            la[i][j]*la[k][j] + la[i][k]*D[j - 1],
            D[j])
        la[i][k] = idiv(t, D[j])
    D[j] = idiv(D[j - 1]*D[k] + la[k][j]**2, D[j])


def hnf(G, _n=3, _d=4):
    """return ``n``, ``H``, ``P`` where ``n`` is the number of
    iterations needed to find ``H`` and ``P`` such that ``H*G = P``
    with ``H`` being in Hermite normal form and ``P`` being unimodular.
    Input ``G`` should be a list of lists where the lists of ``G``
    are the rows of literal integers in G.
    The method used is that of Algorithm 4 in the paper by
    Havas, Majewski and Matthews [1].
    This use of this routine to solve a wide variety of problems
    is described briefly in [2]. XXX how to cite properly?
    See Also
    ========
    igcdLLL - finding extended Euclidean gcd of several integers
    References
    ==========
    .. [1] https://projecteuclid.org/euclid.em/1048515660/
    .. [2] https://rosettacode.org/wiki/Diophantine_linear_system_solving
    """
    def pivot(i):
        # return the index of the first non-zero column
        # else C (one more than last column index)
        for j in range(C):
            if A[i][j]:
                return j
        return C

    def minus(i):
        # negate la in row or col i
        for r in range(1, R):
            for c in range(r - 1):
                if i in (r, c):
                    la[r][c] = -la[r][c]

    def reduce2(k, i):
        # update A, B and la
        col1 = pivot(i)
        if col1 < C:
            if A[i][col1] < 0:
                minus(i)
                A[i] = [-_ for _ in A[i]]
                B[i] = [-_ for _ in B[i]]
            q = A[k][col1]//A[i][col1] # floor
        elif 2*abs(la[k][i]) > D[i]:
            q = iround(la[k][i], D[i])
        else:
            q = 0
        if q != 0:
            for M in (A, B):
                M[k] = [mk - q*mi for mk, mi in zip(M[k], M[i])]
            la[k][i] -= q*D[i]
            for j in range(i):
                la[k][j] -= q*la[i][j]
        col2 = pivot(k)
        return col1, col2

    G = G.tolist()

    if not iterable(G) and G and all(iterable(i) for i in G):
        raise ValueError('G should be list of lists')
    # alpha = _n/_d
    ok = (_d > 0) and (_n > 0) and (_n <= _d)
    if not ok:
        raise ValueError('alpha must be in (0, 1], but is %s/%s' % (_n, _d))

    # get dimensions
    R = len(G)
    C = set([len(i) for i in G])
    if len(C) != 1:
        raise ValueError('all rows should have the same length')
    C = C.pop()

    # set constants and working space
    m = R - 1
    n = C - 1
    B = [[0]*R for i in range(R)]
    for i in range(R):
        B[i][i] = 1
    la = [[0]*i for i in range(R)]
    A = [i[:] for i in G]
    D = [1]*(R+1)  # -1 references last element

    # check for special case where first non-zero column has, as the
    # only entry, a negative value in the last row
    for j in range(C):
        if any(A[i][j] for i in range(R)):
            if A[m][j] < 0 and not any(A[i][j] for i in range(m)):
                A[m] = [-i for i in A[m]]
                B[m][m] = -1
            break  # this was the first non-zero column

    # begin iterations
    k = 1
    loop = 0
    while k <= m:
        loop += 1
        j = k - 1
        col1, col2 = reduce2(k, j)
        if col1 <= min(col2, n) or (col1 == col2 == C and
                _d*(D[j - 1]*D[k] + la[k][j]**2) < _n*D[j]**2):
            _swap(k, R, A, B, la, D)
            if k > 1:
                k = j
        else:
            for i in range(j -1, -1, -1):
                reduce2(k, i)
            k += 1
    H = A[::-1]  # Hermitian norm form
    P = B[::-1]  # unimodular matrix
    # P*G = H
    return loop, Matrix(H), Matrix(P)


def test():

    G = Matrix([[12, 19, 28, 34], [19, 30, 44, 53]])

    n, H, P = hnf(G)
    H = Matrix(H)
    P = Matrix(P)

    assert H == Matrix([[1, 0, -4, -13], [0, 1,  4,  10]])
    assert P == Matrix([[-30,  19], [ 19, -12]])
    assert P*G == H

    M = Matrix([
        [   3,    7,   13,    21,    31,    43,    57,    73,    91,    111],
        [  11,   36,   77,   134,   207,   296,   401,   522,   659,    812],
        [  31,  113,  249,   439,   683,   981,  1333,  1739,  2199,   2713],
        [  69,  262,  583,  1032,  1609,  2314,  3147,  4108,  5197,   6414],
        [ 131,  507, 1133,  2009,  3135,  4511,  6137,  8013, 10139,  12515],
        [ 223,  872, 1953,  3466,  5411,  7788, 10597, 13838, 17511,  21616],
        [ 351, 1381, 3097,  5499,  8587, 12361, 16821, 21967, 27799,  34317],
        [ 521, 2058, 4619,  8204, 12813, 18446, 25103, 32784, 41489,  51218],
        [ 739, 2927, 6573, 11677, 18239, 26259, 35737, 46673, 59067,  72919],
        [1011, 4012, 9013, 16014, 25015, 36016, 49017, 64018, 81019, 100020]])

    H = Matrix([
        [1, 0,  7, 22, 45,  76, 115, 162, 217, 280],
        [0, 1,  4,  9, 16,  25,  36,  49,  64,  81],
        [0, 0, 12, 36, 72, 120, 180, 252, 336, 432],
        [0, 0,  0,  0,  0,   0,   0,   0,   0,   0],
        [0, 0,  0,  0,  0,   0,   0,   0,   0,   0],
        [0, 0,  0,  0,  0,   0,   0,   0,   0,   0],
        [0, 0,  0,  0,  0,   0,   0,   0,   0,   0],
        [0, 0,  0,  0,  0,   0,   0,   0,   0,   0],
        [0, 0,  0,  0,  0,   0,   0,   0,   0,   0],
        [0, 0,  0,  0,  0,   0,   0,   0,   0,   0]])

    # This is different from the paper:
    P = Matrix([
        [-11,  -8, -4,  1,  3,  4,  4,  1,  0, -3],
        [ -2,  -1,  0,  1, -1,  0,  1,  0,  1, -1],
        [-15, -11, -4,  0,  4,  5,  4,  3,  1, -5],
        [ -1,   2, -1, -1,  2,  0, -2,  1,  0,  0],
        [ -1,   1,  0,  1, -1,  0,  0,  1, -2,  1],
        [ -1,   0,  2, -1,  1, -1,  1, -1, -1,  1],
        [  0,  -1,  2,  0, -1, -1,  0,  2, -1,  0],
        [ -1,   0,  1,  1,  1, -2,  0, -1,  1,  0],
        [  0,  -1,  1,  1, -1,  1, -2,  1,  0,  0],
        [ -1,   1,  1,  0, -2,  1,  0,  0,  0,  0]])

    assert P*M == H
    _, H1, P1 = hnf(M)
    assert P1 == P
    assert H1 == H


#test()

M2 = Matrix(100, 100, lambda i, j: (i + 1) ** 3 * (j + 1) ** 2 + i + j + 2)
