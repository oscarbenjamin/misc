#@profile
def shift(x, m, b):
    # Calculate (x << b) % (2**m + 1)
    b %= 2 * m

    if b >= m:
        x = -x
        b -= m

    upper = x >> (m - b)
    x -= upper << (m - b)
    return (x << b) - upper

# FFT using Cooley-Tukey, a divide and conquer algorithm
# NOTE that Cooley-Tukey requires n to be a power of two
#@profile
def NTT(A, m, root = 1):
    n = len(A)
    if n == 1:
        return A

    root2 = root * 2 % (2 * m)
    even = NTT(A[::2], m, root2)
    odd  = NTT(A[1::2], m, root2)

    # Multiply by roots of unity
    odd = [shift(odd[k], m, root * k) for k in range(n//2)]

    return [even[k] + odd[k] for k in range(n//2)] +\
           [even[k] - odd[k] for k in range(n//2)]

def INTT(A, m, root = 1):
    A2 = NTT(A, m, -root)
    inverse_shift = -(len(A).bit_length() - 1)
    return [shift(a, m, inverse_shift) for a in A2]

def int_to_limbs(x, limb, n):
    # Convert big integer x to base 'limb' (a list of length n)
    if n == 1:
        return [x]
    n2 = n//2
    upper = x >> n2 * limb
    lower = x - (upper << n2 * limb)
    return int_to_limbs(lower, limb, n2) + int_to_limbs(upper, limb, n - n2)

def limbs_to_int(A, limb):
    # Convert a number given in base 'limb' into a big integer
    n = len(A)
    if n == 1:
        return A[0]
    n2 = n//2
    return limbs_to_int(A[:n2], limb) + (limbs_to_int(A[n2:], limb) << (n2 * limb))

#@profile
def negacyclic_convolution(A,B,m):
    n = len(A)
    assert len(A) == len(B) <= m
    assert m % n == 0

    root = m // n

    # Weight A and B
    A = [shift(a, m, root * i) for i,a in enumerate(A)]
    B = [shift(b, m, root * i) for i,b in enumerate(B)]

    A = NTT(A, m, 2 * root)
    B = NTT(B, m, 2 * root)

    A = [shift(a, m, 0) for a in A]
    B = [shift(b, m, 0) for b in B]

    if m >= 1e5: # Use ssa if num of bits in a and b >= 1e5
        C = [SSA(a, b, m) for a,b in zip(A,B)]
    else:
        C = [shift(a * b, m, 0) for a,b in zip(A,B)]

    C = INTT(C, m, 2 * root)
    C = [shift(c, m, -root * i) for i,c in enumerate(C)]
    C = [shift(c, m, 0) for c in C]
    return C

#@profile
def SSA(x1, x2, M = None):
    """
        Calculates x1 * x2 if M is None
        Otherwise calculates x1 * x2 % (2**M + 1)
    """

    #n = 256 # Sṕlit factor, configurable, needs to be a factor of 2

    if M is None:
        n1 = x1.bit_length()
        n2 = x2.bit_length()
        limb = (n1 + n2 + n - 1) // n
        M = limb * n
    else:
        assert M % n == 0
        limb = M // n

    lower_lim_m = 2 * limb + n.bit_length() + 1
    m = ((lower_lim_m + n - 1)//n) * n

    A = int_to_limbs(x1, limb, n)
    B = int_to_limbs(x2, limb, n)

    C = negacyclic_convolution(A, B, m)

    # Detect underflow in output of the negacyclic convolution
    C = [c - ((1 << m) + 1) if c >> (m - 1) else c for c in C]

    return limbs_to_int(C, limb)


#a = 12345**20000; b = 13245**20000
a = 12345**500000; b = 13245**500000
#a = 12345**2000000; b = 13245**2000000
#a = 12345**4000000; b = 13245**4000000
a = b = int('1' * 100000000, 2)
#n = 2**13
n = 2**12

print('Multiplying',a.bit_length(),'numbers with n=', n)
import time

for c in range(1):
    L = time.perf_counter()
    y1 = SSA(a,b)
    R = time.perf_counter()
    print(c, 'SSA took', R - L, flush=True)

#L = time.perf_counter()
#y2 = a * b
#R = time.perf_counter()
#print('Built in took', R - L)

#assert y1 == y2
