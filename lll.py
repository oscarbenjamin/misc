#!/usr/bin/env python


from sympy import sympify
from mpmath import mp, mpf


def dot(a, b):
    return sum(ai*bi for ai, bi in zip(a, b))

def sub(a, b):
    return [ai - bi for ai, bi in zip(a, b)]

def mul(a, b):
    return [a*bi for bi in b]


def gram_schmidt(b):
    n = len(b)
    mu = [[None]*n for _ in range(n)]
    bstar = [[None]*n for _ in range(n)]

    for i in range(n):
        bi = b[i]
        for j in range(i):
            mu[i][j] = dot(b[i], bstar[j]) / dot(bstar[j], bstar[j])
            bi = sub(bi, mul(mu[i][j], bstar[j]))
        bstar[i] = bi

    return bstar, mu


def lll(b, delta=0.75):
    n = len(b)
    bstar, mu = gram_schmidt(b)
    k = 1
    while k < n:
        for j in reversed(range(k)):
            if abs(mu[k][j]) > 0.5:
                b[k] = sub(b[k], mul(round(mu[k][j]), b[j]))
                bstar, mu = gram_schmidt(b)
        if dot(bstar[k], bstar[k]) >= (delta - mu[k][k-1]**2) * dot(bstar[k-1], bstar[k-1]):
            k = k + 1
        else:
            b[k], b[k-1] = b[k-1], b[k]
            bstar, mu = gram_schmidt(b)
            k = max(k-1, 1)
    return b


from sympy.core.evalf import dps_to_prec
from sympy.polys.matrices import DomainMatrix


class PrecisionExhausted(Exception):
    pass


def pslq(x, dps=15):
    gamma = 2

    n = len(x)
    x = [sympify(xi).evalf(dps) for xi in x]
    A = eye(n)
    B = eye(n)

    srev = []
    ssq = RR(0)
    for xj in reversed(x):
        ssq += xj ** 2
        srev.append(ssq ** 0.5)

    s = srev[::-1]
    t = 1/s[0]
    y = [t * xi for xi in x]
    s = [t * si for si in s]

    def Hij(i, j):
        if i < j:
            return 0
        elif i == j:
            return s[j+1]/s[j]
        else:
            return -y[i]*y[j]/(s[j]*s[j+1])

    H = Matrix(n, n-1, Hij)

    for i in range(1, n):
        for j in reversed(range(i)):
            t = round(H[i,j]/H[j,j])
            y[j] = y[j] + t*y[i]
            for k in range(j+1):
                H[i,k] = H[i,k] - t*H[j,k]
            for k in range(n):
                A[i,k] = A[i,k] - t*A[j,k]
                B[k,j] = B[k,j] + t*B[k,i]

    maxy = max(abs(yi) for yi in y)

    while True:
        m = max(range(n-1), key=lambda i: gamma**i * abs(H[i,i]))
        y[m], y[m+1] = y[m+1], y[m]
        A[m,:], A[m+1,:] = A[m+1,:], A[m,:]
        H[m,:], H[m+1,:] = H[m+1,:], H[m,:]
        B[:,m], B[:,m+1] = B[:,m+1], B[:,m]

        if m < n-2:
            t0 = sqrt(H[m,m]**2 + H[m,m+1]**2)
            t1 = H[m,m] / t0
            t2 = H[m,m+1] / t0
            for i in range(m, n):
                t3 = H[i,m]
                t4 = H[i,m+1]
                H[i,m] = t1*t3 + t2*t4
                H[i,m+1] = -t2*t3 + t1*t4

        for i in range(m+1, n):
            for j in reversed(range(min(i-1,m+1)+1)):
                t = round(H[i,j]/H[j,j])
                y[j] = y[j] + t*y[i]
                for k in range(j+1):
                    H[i,k] = H[i,k] - t*H[j,k]
                for k in range(n):
                    A[i,k] = A[i,k] - t*A[j,k]
                    B[k,j] = B[k,j] + t*B[k,i]

        imin = None
        for i in range(n):
            if y[i] == 0:
                break
        else:
            i = min(range(n), key=lambda i: abs(y[i]))
        #if abs(y[i]) < 10**-dps or abs(y[i])/maxy < 10**-10:
        if abs(y[i])/maxy < 10**-10:
            return B[:,i].flat()
        maxy = abs(y[i])

        if max(map(abs, A.flat())) > 10**dps:
            raise PrecisionExhausted

        #print([H[i,i] for i in range(n-1)]           )
        maxabs = max(abs(H[j,j]) for j in range(n-1))
        M = 1/maxabs
        print(M, min(map(abs, y)))
        #print(B.T)


def pslq_eq(x, dps=15):
    coeffs = pslq(x, dps)
    terms = (Mul(ci, xi, evaluate=False) for ci, xi in zip(coeffs, x))
    return Eq(Add(*terms, evaluate=False), 0, evaluate=False)


def minpoly_pslq(a, deg, dps=15):
    powers = [a**i for i in range(deg+1)]
    coeffs = pslq(powers, dps)
    while coeffs[0] == 0:
        coeffs.pop(0)
    p = Poly(coeffs[::-1], z)
    if p.LC() < 0:
        p = -p
    return p


def integer_relation(num, order, dps=10, N=None):
    if N is None:
        N = 10 ** int(dps)
    mp.dps = dps
    num = mpf(sympify(num).evalf(dps)._mpf_)

    b = [[0] * (order+2) for _ in range(order+1)]
    for i in range(order+1):
        b[i][i] = 1
        b[i][-1] = N*num**(order-i)

    #b = [[QQ(i) for i in row] for row in b]

    b_reduced = lll(b)

    b_best = min(b_reduced, key=lambda bi: abs(bi[-1]))[:-1]
    p = Poly(b_best, x, domain=QQ)

    return p, Matrix(b_reduced)


#b = [[1, 1, 1], [-1, 0, 2], [3, 5, 6]]
#
#for bi in lll(b):
#    print(bi)
