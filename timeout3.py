import signal
from sympy import *


class Timeout(Exception):
    pass


def loop_forever():
    F, G, H, a, b, c, d, e, f, g, h, i, k, m, n, p, r, s, t, x=symbols('F G H a b c d e f g h i k m n p r s t x')
    integrand=((-x*exp(exp(5))*ln(x)+(x**2+13*x+42)*exp(exp(5)))*exp(ln(x)/(6+x))+((-2*x**4-24*x**3-73*x**2-12*x-36)*exp(x**2)+(-x**2-12*x-36)*exp(2))*exp(exp(5)))/(x**2+12*x+36)
    integrate(integrand,x)

def handler(signum, frame):
    print("Forever is over!")
    raise Timeout("end of time")

signal.signal(signal.SIGALRM, handler)
signal.alarm(5)

try:
    loop_forever()
except Timeout:
    1/0
except Exception:
    pass
else:
    pass
