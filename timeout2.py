import sympy as sp
from sympy import *
from sympy import Eq

def loop_forever():
    x = symbols('x')
    e = S('''-8*A - x**5*(A + B + C + E) - x**4*(4*A - 2**(3/4)*C + 4*C + D + 2**(3/4)*E + 4*E + F) - x**3*(-4*2**(3/4)*C + sqrt(2)*C - 2**(3/4)*D + 4*D + sqrt(2)*E + 4*2**(3/4)*E + 2**(3/4)*F + 4*F) - x**2*(4*sqrt(2)*C - 4*2**(3/4)*D + sqrt(2)*D + 4*sqrt(2)*E + sqrt(2)*F + 4*2**(3/4)*F) - x*(2*A + 2*B + 4*sqrt(2)*D + 4*sqrt(2)*F) + 5''')
    eqn = e.subs(exp(1), Symbol("E"))
    solve_undetermined_coeffs(eqn,var('A:F'),x)

import signal

def handler(signum, frame):
    print("Forever is over!")
    raise Exception("end of time")

signal.signal(signal.SIGALRM, handler)
signal.alarm(5)

try:
    loop_forever()
except Exception:
    pass
else:
    assert False
