from sympy import igcd
from collections import deque

def associated(x, y):
    return abs(x) == abs(y)

def factorial_basis(d):
    d = (di for di in d if not associated(di, 1))
    basis = set()
    queue = deque(d)
    for si in queue:
        for bi in basis:
            pass
