
def dot(a, b):
    return sum(ai*bi for ai, bi in zip(a, b))

def sub(a, b):
    return [ai - bi for ai, bi in zip(a, b)]

def mul(a, b):
    return [a*bi for bi in b]

def add(a, b):
    return [ai + bi for ai, bi in zip(a, b)]


def lll(b):
    # initialise
    n = len(b)
    bstar = [None] * n
    B = [None] * n
    mu = [[None] * n for _ in range(n)]

    k = 1
    kmax = 0
    bstar[0] = b[0]
    B[0] = dot(b[0], b[0])

    def red(k, l):
        if abs(mu[k][l]) <= 0.5:
            return
        q = round(mu[k][l])
        b[k] = sub(b[k], mul(q, b[l]))
        mu[k][l] = mu[k][l] - q
        for i in range(l):
            mu[k][i] = mu[k][i] - q*mu[l][i]

    def swap(k):
        b[k], b[k-1] = b[k-1], b[k]
        if k > 1:
            for j in range(k-1):
                mu[k][j], mu[k-1][j] = mu[k-1][j], mu[k][j]

        mu_ = mu[k][k-1]
        B_ = B[k] + mu_**2*B[k-1]
        mu[k][k-1] = mu_*B[k-1]/B_
        b_ = bstar[k-1]
        bstar[k-1] = add(bstar[k], mul(mu_, b_))
        bstar[k] = add(mul(-mu[k][k-1], bstar[k]), mul(B[k]/B_, b_))
        B[k] = B[k-1]*B[k]/B_
        B[k-1] = B_

        for i in range(k+1, kmax+1):
            t = mu[i][k]
            mu[i][k] = mu[i][k-1] - mu_*t
            mu[i][k-1] = t + mu[k][k-1]*mu[i][k]

    while k < n:
        # 2. incremental Gram-Schmidt
        if k > kmax:
            kmax = k
            bstar[k] = b[k]
            for j in range(k):
                mu[k][j] = dot(b[k], bstar[j]) / B[j]
                bstar[k] = sub(b[k], mul(mu[k][j], bstar[j]))
            B[k] = dot(bstar[k], bstar[k])

            if B[k] == 0:
                raise ValueError("Not a valid basis")


        # 3. Test LLL condition
        red(k, k-1)

        while B[k] < (0.75 - mu[k][k-1]**2) * B[k-1]:
            swap(k)
            k = max(1, k-1)

        for l in reversed(range(k-1)):
            red(k, l)

        k = k + 1

    # Finished
    return b

