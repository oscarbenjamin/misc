#
# Parse and format big integers in different bases.
#
# Richard P. Brent and Paul Zimmermann, Modern Computer Arithmetic
#
# https://members.loria.fr/PZimmermann/mca/mca-cup-0.5.9.pdf
#


# Algorithm 1.25:

def parse_int(S: str, B: int = 10) -> int:
    """Parse string S as an integer in base B"""
    m = len(S)
    l = list(map(int, S[::-1]))
    b, k = B, m    # invariant len(l) == k
    while k > 1:
        last = [l[-1]] if k % 2 == 1 else []
        l = [l1 + b*l2 for l1, l2 in zip(l[::2], l[1::2])]
        l.extend(last)
        b, k = b**2, (k + 1) >> 1
    [l0] = l
    return l0


from functools import lru_cache

@lru_cache
def integer_log(y: int, x: int):
    r = e = 0
    while y >= x:
        d = x
        m = 1
        while y >= d:
            y, rem = divmod(y, d)
            r = r or rem
            e += m
            if y > d:
                d *= d
                m *= 2
    return e

# Algorithm 1.26:

#@profile
def format_int(A: int, B: int = 10) -> str:
    """Format integer A as a string in base B"""
    if A.bit_length() < 1000:
        # Here we use str for the base case but really this should be something
        # that works for arbitrary bases using the standard O(n^2) algorithm.
        # The threshold should be tuned.
        return str(A)
    else:
        # find k so that B**(2*k-2) <= A < B**(2*k)
        e = integer_log(A, B)
        if e % 2 == 1:
            k = (e + 1) // 2
        else:
            k = e // 2 + 1
        #assert B**(2*k - 2) <= A < B**(2*k)

        Bk = B**k
        Q, R = divmod(A, Bk)
        r = format_int(R, B)
        q = format_int(Q, B)
        pad = '0'*(k - len(r))
        return ''.join([q, pad, r])


#n = parse_int('1'*10**6)
n = parse_int('1' * 10**6)
print(len(format_int(n)))
