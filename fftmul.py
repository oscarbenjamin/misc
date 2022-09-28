def check(A, omega, K, modulo):
    """Check that omega, K and modulo are suitable for fft of A"""
    assert len(A) == K
    assert K % 2 == 0
    assert pow(omega, K, modulo) == 1
    assert pow(omega, K // 2, modulo) == -1 % modulo


def fft_slow(A, omega, K, modulo):
    """Basic quadratic DFT algoritm"""
    check(A, omega, K, modulo)
    return [sum(omega**(i*j)*A[j] for j in range(K)) % modulo for i in range(K)]


def ifft_slow(A, omega, K, modulo):
    """Basic quadratic inverse DFT algoritm"""
    check(A, omega, K, modulo)
    return [sum(omega**((-i*j) % K)*A[j] for j in range(K)) % modulo for i in range(K)]


def fft(A, omega, K, modulo):
    """Fast N*log(N) DFT algorithm"""
    #check(A, omega, K, modulo)
    if K == 2:
        a0, a1 = A
        #return [(a0 + a1) % modulo, (a0 - a1) % modulo]
        return [a0 + a1, a0 - a1]
    else:
        omega2 = omega ** 2
        K2 = K // 2
        E = fft(A[::2], omega2, K2, modulo)
        O = fft(A[1::2], omega2, K2, modulo)

        A2 = [None] * K
        for j in range(K2):
            p, q = E[j], omega**j * O[j]
            #A2[j], A2[K2+j] = (p + q) % modulo, (p - q) % modulo
            A2[j], A2[K2+j] = (p + q), (p - q)
        return A2


def ifft(A, omega, K, modulo):
    """Fast N*log(N) inverse DFT algorithm"""
    #check(A, omega, K, modulo)
    #omegainv = pow(omega, K-1, modulo)
    omegainv = omega**(K-1)
    return fft(A, omegainv, K, modulo)


def bitsplit(A, M, K=None):
    """Split integer A into K M-bit integers, padding with zeros."""
    As = bin(A)[2:]
    B = len(As)
    num, start = divmod(B, M)
    Adigits = [int(As[:start], 2)] if start else []
    for i in range(num):
        segment = As[i*M+start:(i+1)*M+start]
        Adigits.append(int(segment, 2))
    assert As == ''.join(bin(a)[2:].rjust(M, '0') for a in Adigits).lstrip('0')
    Adigits = Adigits[::-1]
    if K is not None:
        Adigits += [0] * (K - len(Adigits))
    return Adigits


def bitjoin(Adigits, M):
    """Inverse of bitsplit"""
    assert all(a.bit_length() <= M for a in Adigits)
    As = ''.join(bin(a)[2:].rjust(M, '0') for a in Adigits[::-1])
    return int(As, 2)


#@profile
def fftmulmod(A, B, n, K, k):
    """Compute A*B using FFT multiplication.

    The algorithm generally computes A*B mod 2**n + 1 but we require n to be
    large enough that this is just A*B. Also we require n = M*K and K = 2**k.

    Basic steps:
    - Split the bits of A and B and pad with zeros so that both are
      represented by K (=2**k) digits in base 2**M.
    - Compute n' ~ sqrt(n) so that we can work modulo 2**n' + 1.
    - Premultiply both digit sequences by theta**j
    - Fourier transform both digit sequences in Z / <2**n' + 1>
    - Compute the convolution with pointwise multiplication.
    - Post divide by K*theta**j (note: ifft(fft(x)) = K*x).
    - Inverse transform to get digits of C = A*B in base 2**M
    - Assemble the digits to get C as a single int.
    """
    assert A.bit_length() + B.bit_length() < n
    assert 2 ** k == K
    M = n // K
    assert M * K == n
    Adigits = bitsplit(A, M, K)
    Bdigits = bitsplit(B, M, K)
    assert len(Adigits) == K
    assert len(Bdigits) == K

    bound = (2*n) // K + k
    nprime = bound + K - bound % K
    assert nprime % K == 0
    theta = 2 ** (nprime // K)
    omega = theta ** 2

    modulus = 2**nprime + 1

    for j in range(K):
        Adigits[j] = (theta**j * Adigits[j]) % modulus
        Bdigits[j] = (theta**j * Bdigits[j]) % modulus

    Afreq = fft(Adigits, omega, K, modulus)
    Bfreq = fft(Bdigits, omega, K, modulus)

    Cfreq = [None] * K
    for j in range(K):
        Cfreq[j] = (Afreq[j]*Bfreq[j])
        Cfreq[j] %= modulus

    Cdigits = ifft(Cfreq, omega, K, modulus)

    for j in range(K):
        Cdigits[j] = (Cdigits[j] * pow(K * theta**j, -1, modulus))
        Cdigits[j] %= modulus
        if Cdigits[j] >= (j + 1) * 2**(2*M):
            Cdigits[j] = Cdigits[j] - modulus

    # We can't use bitjoin here because the digits are larger than the base.
    # (the bitstrings overlap). It should be possible to do this more
    # efficiently than below though by iteratively normalising the digits.
    C = sum(Cdigits[j]*2**(j*M) for j in range(K))
    return C


def fft_mont(A, omega, K):
    """Fast N*log(N) DFT algorithm"""
    #check(A, omega, K, modulo)
    if K == 2:
        a0, a1 = A
        #return [(a0 + a1) % modulo, (a0 - a1) % modulo]
        return [a0 + a1, a0 - a1]
    else:
        omega2 = omega ** 2
        K2 = K // 2
        E = fft_mont(A[::2], omega2, K2, modulo)
        O = fft_mont(A[1::2], omega2, K2, modulo)

        A2 = [None] * K
        for j in range(K2):
            p, q = E[j], omega**j * O[j]
            #A2[j], A2[K2+j] = (p + q) % modulo, (p - q) % modulo
            A2[j], A2[K2+j] = (p + q), (p - q)
        return A2


def ifft_mont(A, omega, K):
    """Fast N*log(N) inverse DFT algorithm"""
    #check(A, omega, K, modulo)
    omegainv = pow(omega, K-1)
    return fft_mont(A, omegainv, K)


def fftmulmod_mont(A, B, n, K, k):
    """Compute A*B using FFT multiplication.

    The algorithm generally computes A*B mod 2**n + 1 but we require n to be
    large enough that this is just A*B. Also we require n = M*K and K = 2**k.

    Basic steps:
    - Split the bits of A and B and pad with zeros so that both are
      represented by K (=2**k) digits in base 2**M.
    - Compute n' ~ sqrt(n) so that we can work modulo 2**n' + 1.
    - Premultiply both digit sequences by theta**j
    - Fourier transform both digit sequences in Z / <2**n' + 1>
    - Compute the convolution with pointwise multiplication.
    - Post divide by K*theta**j (note: ifft(fft(x)) = K*x).
    - Inverse transform to get digits of C = A*B in base 2**M
    - Assemble the digits to get C as a single int.
    """
    assert A.bit_length() + B.bit_length() < n
    assert 2 ** k == K
    M = n // K
    assert M * K == n
    Adigits = bitsplit(A, M, K)
    Bdigits = bitsplit(B, M, K)
    assert len(Adigits) == K
    assert len(Bdigits) == K

    bound = (2*n) // K + k
    nprime = bound + K - bound % K
    assert nprime % K == 0
    theta = 2 ** (nprime // K)
    omega = theta ** 2

    mont = montgomery_gen(nprime)

    modulus = 2**nprime + 1

    theta = mont.from_int(theta)
    omega = mont.from_int(omega)
    Kmont = mont.from_int(K)
    Adigits = [mont.from_int(ai) for ai in Adigits]
    Bdigits = [mont.from_int(bi) for bi in Bdigits]

    for j in range(K):
        Adigits[j] = theta**j * Adigits[j]
        Bdigits[j] = theta**j * Bdigits[j]

    #Adigits = [ai.to_int() for ai in Adigits]
    #Bdigits = [bi.to_int() for bi in Bdigits]

    Afreq = fft(Adigits, omega, K, modulus)
    Bfreq = fft(Bdigits, omega, K, modulus)

    Cfreq = [None] * K
    for j in range(K):
        Cfreq[j] = (Afreq[j]*Bfreq[j])
        #Cfreq[j] %= modulus

    Cdigits = ifft(Cfreq, omega, K, modulus)

    for j in range(K):
        Cdigits[j] = (Cdigits[j] * pow(Kmont * theta**j, -1, modulus))
        Cdigits[j] %= modulus
        if Cdigits[j] >= (j + 1) * 2**(2*M):
            Cdigits[j] = Cdigits[j] - modulus

    # We can't use bitjoin here because the digits are larger than the base.
    # (the bitstrings overlap). It should be possible to do this more
    # efficiently than below though by iteratively normalising the digits.
    C = sum(Cdigits[j]*2**(j*M) for j in range(K))
    return C


def montgomery_gen(n):
    """Set up a Montgomery form for arithmetic modulo 2**n + 1

    Returns 4 functions:
    - convert_to: Convert an integer into Montgomery form
    - convert_from: Convert from Montgomery form back ordinary integer
    - add: Efficiently add in Montgomery form
    - multiply: Efficiently mulyiply in Montgomery form
    """
    N = 2**n + 1
    R = 2**n
    Ninv = 2**n - 1

    def to_mont(e):
        return (e*R) % N

    def from_mont(e):
        return multiply(e, 1)

    def add(a, b):
        c = a + b
        if c > N:
            c -= N
        return c

    def sub(a, b):
        c = a - b
        if c < 0:
            c += N
        return c

    def mul(c, d):
        x = c*d
        y = x + N*((x*Ninv) & (R - 1))
        #assert y % R == 0
        z = y // R
        if z >= N:
            z -= N
        return z

    def pow(a, n):
        if n == 0:
            return R
        elif n % 2 == 0:
            return pow(x * x, n // 2)
        else:
            return x * pow(x * x, (n - 1) // 2)


    class MontgomeryForm:
        def __init__(self, value):
            self.value = value

        @classmethod
        def from_int(cls, integer):
            return cls((integer*R) % N)

        def to_int(self):
            return mul(self.value, 1)

        def __add__(self, other):
            assert isinstance(other, MontgomeryForm)
            return MontgomeryForm(add(self.value, other.value))

        def __sub__(self, other):
            assert isinstance(other, MontgomeryForm)
            return MontgomeryForm(sub(self.value, other.value))

        def __mul__(self, other):
            assert isinstance(other, MontgomeryForm)
            return MontgomeryForm(mul(self.value, other.value))

        def __pow__(self, other):
            assert isinstance(other, int) and other >= 0
            return MontgomeryForm(pow(self.value, other))

    return MontgomeryForm



#to_mont, from_mont, add, multiply = montgomery_gen(10)


a = int('1'*100000)
b = fftmulmod(a, a, 2**20, 2**3, 3)
#b = fftmulmod_mont(a, a, 2**20, 2**3, 3)
assert b == a*a
