from sympy import *
from sympy.polys.matrices import DomainMatrix


def faddeev_leverrier(A):
    Adm = DomainMatrix.from_Matrix(A)
    K = Adm.domain
    domain = Adm.domain.unify(QQ)
    Adm = Adm.convert_to(domain)
    return _fl(Adm)

def trace(A):
    n = A.shape[0]
    return sum((A[i,i].element for i in range(n)))

def _fl(A):
    K = A.domain
    n = A.shape[0]
    ck = K.one
    cs = [ck]
    I = DomainMatrix.eye(n, K)
    Mk = I
    for k in range(1, n):
        AMk = A*Mk
        ck = (-K.one/k)*trace(AMk)
        Mk = AMk + ck*I
        cs.append(ck)
    s = (-1)**(n-1)
    #detA = s*(-S.One/n)*(A*Mk).trace()
    #cs.append(detA)
    adjA = s*Mk
    return adjA
