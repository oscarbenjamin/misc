def bitrev(j, K):
    k = 1
    while 2**k < K:
        k += 1
    return int(bin(j)[2:].rjust(k, '0')[::-1], 2)


def forward_fft(A, omega, K, modulo):
    assert len(A) == K
    assert 1 == (omega ** K) % modulo < (omega ** (K // 2)) % modulo
    if K == 2:
        a0, a1 = A
        return [(a0 + a1) % modulo, (a0 - a1) % modulo]
    else:
        omega2 = omega ** 2
        K2 = K // 2
        A2 = [None] * K
        A2[::2] = forward_fft(A[::2], omega2, K2, modulo)
        A2[1::2] = forward_fft(A[1::2], omega2, K2, modulo)
        for j in range(K2):
            a2j, a2j1 = A2[2*j], A2[2*j+1]
            omega_a2j1 = omega ** (K-j) * a2j1
            A2[2*j], A[2*j+1] = (a2j + omega_a2j1) % modulo, (a2j - omega_a2j1) % modulo
        return A2


#def forward_fft(A, omega, K, modulo):
    #if N = 1

def check(A, omega, K, modulo):
    assert len(A) == K
    assert 1 == (omega ** K) % modulo < (omega ** (K // 2)) % modulo


def fft_slow(A, omega, K, modulo):
    check(A, omega, K, modulo)
    return [sum(omega**(i*j)*A[j] for j in range(K)) % modulo for i in range(K)]


def ifft_slow(A, omega, K, modulo):
    check(A, omega, K, modulo)
    return [sum(omega**((-i*j) % K)*A[j] for j in range(K)) % modulo for i in range(K)]


def fft(A, omega, K, modulo):
    check(A, omega, K, modulo)
    if K == 2:
        a0, a1 = A
        return [(a0 + a1) % modulo, (a0 - a1) % modulo]
    else:
        omega2 = omega ** 2
        K2 = K // 2
        E = fft(A[::2], omega2, K2, modulo)
        O = fft(A[1::2], omega2, K2, modulo)

        A2 = [None] * K
        for j in range(K2):
            p, q = E[j], omega**j * O[j]
            A2[j], A2[K2+j] = (p + q) % modulo, (p - q) % modulo
        return A2


def ifft(A, omega, K, modulo):
    check(A, omega, K, modulo)
    omegaprime = (omega ** (K-1)) % modulo
    return fft(A, omegaprime, K, modulo)


def bitsplit(A, M, K=None):
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
    assert all(a.bit_length() <= M for a in Adigits)
    As = ''.join(bin(a)[2:].rjust(M, '0') for a in Adigits[::-1])
    return int(As, 2)


def fftmulmod(A, B, n, K, k):
    assert 2 ** k == K
    M = n // K
    assert M * K == n
    Adigits = bitsplit(A, M, K)
    Bdigits = bitsplit(B, M, K)

    bound = (2*n) // K + k
    nprime = bound + K - bound % K
    assert nprime % K == 0
    theta = 2 ** (nprime // K)
    omega = theta ** 2

    modulus = 2**nprime + 1

    print(modulus)

    for j in range(K):
        Adigits[j] = (theta**j * Adigits[j]) % modulus
        Bdigits[j] = (theta**j * Bdigits[j]) % modulus

    Afreq = fft(Adigits, omega, K, modulus)
    Bfreq = fft(Bdigits, omega, K, modulus)

    Cfreq = [None] * K
    for j in range(K):
        Cfreq[j] = (Afreq[j]*Bfreq[j]) % modulus

    Cdigits = ifft(Cfreq, omega, K, modulus)

    print(Adigits)
    print(Bdigits)
    print(Cdigits)

    for j in range(K):
        Cdigits[j] = (Cdigits[j] // (K * theta**j)) % modulus
        if Cdigits[j] >= (j + 1) * 2**(2*M):
            Cdigits[j] = Cdigits[j] - modulus

    C = sum(Cdigits[j]*2**(j*M) for j in range(K))
    return C



#def forward_fft(A, omega, K, modulo):
#    assert len(A) == K
#    assert 1 == (omega ** K) % modulo < (omega ** (K // 2)) % modulo
#    if K == 2:
#        a0, a1 = A
#        return [(a0 + a1) % modulo, (a0 - a1) % modulo]
#    else:
#        omega2 = opmega**2


def backward_fft(A, omega, K, modulo):
    assert len(A) == K
    assert 1 == (omega ** K) % modulo < (omega ** (K // 2)) % modulo
    if K == 2:
        a0, a1 = A
        return [(a0 + a1) % modulo, (a0 - a1) % modulo]
    else:
        omega2 = omega ** 2
        K2 = K // 2
        A2 = [None] * K
        A2[:K2] = backward_fft(A[:K2], omega2, K2, modulo)
        A2[K2:] = backward_fft(A[K2:], omega2, K2, modulo)
        for j in range(K2):
            aj, aK2j = A2[j], A2[K2+j]
            omega_aK2j = omega ** (K-j) * aK2j
            A2[j], A2[K2+j] = (aj + omega_aK2j) % modulo, (aj - omega_aK2j) % modulo
        return A2


def primitive_root(omega, K, modulo):
    if not (1 == (omega ** K) % modulo < (omega ** (K // 2)) % modulo):
        return False
    else:
        return all(sum(omega**(i*j) for j in range(K)) for i in range(K))

