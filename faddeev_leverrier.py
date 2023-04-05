
def faddeev_leverrier(A):
    assert A.is_square

    n = A.shape[0]
    Mk = zeros(n, n)
    In = eye(n)
    c = [S.One]

    for k in range(1, n+1):
        Mk = A*Mk + c[-1]*In
        ck = -(A*Mk).trace()/k
        c.append(ck)

    detA = c[-1]
    adjA = Mk
    pA = c

    if n % 2 == 1:
        detA = -detA
    else:
        adjA = -adjA

    #assert detA == A.det()
    #assert adjA == A.adjugate()
    #assert pA == A.charpoly().all_coeffs()

    return detA, adjA, pA


def inv(A):
    detA, adjA, pA = faddeev_leverrier(A)
    if detA == 0:
        raise ZeroDivisionError
    return adjA / detA
