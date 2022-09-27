from functools import total_ordering
from sympy import *
from sympy.core.evalf import dps_to_prec
from mpmath import workprec, mpf, mp, mpi, MPIntervalContext


def coeffs(polyexpr, x):
    """Convert Expr to a poly

    >>> x = symbols('x')
    >>> coeffs(x**2 - 1, x)
    [mpz(1), mpz(0), mpz(-1)]
    """
    return Poly(polyexpr, x).rep.rep


def diff(poly):
    """Differentiate a poly

    >>> diff([1, 0, -1])
    [2, 0]
    """
    exponents = range(len(poly))[::-1]
    return [n*c for c, n in zip(poly[:-1], exponents)]


def horner(poly, x0):
    """Evaluate a poly

    >>> horner([1, 0, -1], 3)
    8
    """
    val = poly[0]
    for c in poly[1:]:
        val = c + val*x0
    return val


def float_poly(poly):
    """Convert exact poly to float

    >>> float_poly([2, 3])
    [2.0, 3.0]
    """
    return [float(c) for c in poly]


def newton_float(poly, polydiff, a, b):
    """Newton method with floats

    Finds an approximate root of poly between a and b

    >>> res, conv = newton_float([1.0, 0.0, -2.0], [2.0, 0.0], 1, 2)
    >>> res
    1.414213562373095
    >>> conv
    True
    """
    xi = (a + b)/2
    while True:
        f = horner(poly, xi)
        fp = horner(polydiff, xi)
        xi1 = xi - f/fp
        if not (a < xi1 < b):
            return xi, False
        if abs(xi1 - xi) < 1e-15*(abs(xi1) + abs(xi)):
            return xi1, True
        xi = xi1


def newton_bounded(poly, polydiff, a, b, eps=1e-15):
    """Newton method with bounds

    Finds an approximate root of poly between a an b using eps as an error
    tolerance for the root.

    Alo returns a boolean conv indicating whether Newton iteration diverged
    from the interval (a, b).

    >>> res, conv = newton_bounded([1.0, 0.0, -2.0], [2.0, 0.0], 1, 2)
    >>> res
    1.414213562373095
    >>> conv
    True
    """
    xi = (a + b)/2
    while True:
        f = horner(poly, xi)
        fp = horner(polydiff, xi)
        xi1 = xi - f/fp
        if not (a < xi1 < b):
            return xi, False
        if abs(xi1 - xi) < eps*(abs(xi1) + abs(xi)):
            return xi1, True
        xi = xi1


def newton_mpf(poly, polydiff, a, b, prec=None):
    """Newton method with multiprecision

    Finds an approximate root of poly between a and b using a working precision
    of prec binary digits. If prec=None (the default) then the calculation is
    done with ordinary 53 bit precision floats.

    >>> prec = 100
    >>> res, conv = newton_mpf([1, 0, -2], [2, 0], 1, 2, prec)
    >>> Float(res, precision=prec)
    1.4142135623730950488016887242
    >>> conv
    True
    """
    if prec is None:
        return newton_bounded(poly, polydiff, a, b, eps=1e-15)
    else:
        with workprec(prec):
            if QQ.of_type(a):
                a = mpf(a.numerator) / mpf(a.denominator)
            else:
                a = mpf(a)
            if QQ.of_type(b):
                b = mpf(b.numerator) / mpf(b.denominator)
            else:
                b = mpf(b)
            eps = 100*mp.eps()
            return newton_bounded(poly, polydiff, a, b, eps)


def is_proper_subset(X, Y):
    """Check if interval X is a proper subset of interval Y

    >>> is_proper_subset(mpi(1, 2), mpi(0, 3))
    True
    """
    return X.a > Y.a and X.b < Y.b


def intersection(X, Y):
    """Compute the intersection of two intervals X and Y

    Returns None if the intersection is the empty set:

    >>> intersection(mpi(0, 2), mpi(1, 3))
    mpi('1.0', '2.0')
    >>> print(intersection(mpi(0, 1), mpi(2, 3)))
    None
    """
    a = max(X.a, Y.a)
    b = min(X.b, Y.b)
    if a <= b:
        return mpi(a, b)
    else:
        return None


def interval_newton(poly, polydiff, X1):
    """Single step of interval-Newton method

    Given a poly and its derivative and an interval X1 perform a single step of
    the interval Newton method and return a new interval X2 that is not larger
    than X1 (usually a lot smaller) and certainly contains any roots of poly
    that are in X1.

    Returns a tuple (X2, status). There status variable can have one of three
    possible values:

    - 'empty' means that the interval X1 does not contain any roots. In this
      case X2 will be None

    - 'unique' means that the interval X2 (and likewise X1) contains exactly
      one root.

    - 'continue' means that some progress was made but more iterations are
      needed.

    - 'stalled' means that no progress was made and the new interval is no
      smaller. In this case the interval could be split or more precision might
      be needed.

    Examples
    ========

    >>> interval_newton([1.0, 0.0, -2.0], [2.0, 0.0], mpi(1, 2))
    (mpi('1.375', '1.4375'), 'unique')
    >>> interval_newton([1.0, 0.0, -2.0], [2.0, 0.0], mpi(2, 3))
    (None, 'empty')
    >>> interval_newton([1.0, 0.0, -2.0], [2.0, 0.0], mpi(-1, 2))
    (mpi('-1.0', '2.0'), 'stalled')
    >>> interval_newton([1.0, 0.0, -2.0], [2.0, 0.0], mpi(0, 2))
    (mpi('1.25', '2.0'), 'iterate')

    When the status is 'iterate' it is potentially possible to refine the
    interval with another step:

    >>> X2, status = interval_newton([1.0, 0.0, -2.0], [2.0, 0.0], mpi(0, 2))
    >>> X2
    mpi('1.25', '2.0')
    >>> status
    'iterate'
    >>> interval_newton([1.0, 0.0, -2.0], [2.0, 0.0], X2)
    (mpi('1.3687499999999999', '1.46484375'), 'unique')

    References
    ==========

    https://hal.inria.fr/inria-00072253/document
    https://en.wikipedia.org/wiki/Newton%27s_method#Interval_Newton's_method
    https://en.wikipedia.org/wiki/Interval_arithmetic#Interval_Newton_method
    """

    x1 = X1.mid
    f = horner(poly, x1)
    Fp = horner(polydiff, X1)
    Z2 = x1 - f/Fp
    X2 = intersection(Z2, X1)

    if X2 is None:
        status = 'empty'
    elif is_proper_subset(Z2, X1):
        status = 'unique'
    elif X2 == X1:
        status = 'stalled'
    else:
        status = 'iterate'
    return X2, status


class Infinite:

    def mul_q(self, other):
        if other == 0:
            raise ValueError('0 * oo')
        elif other > 0:
            return self
        else:
            return -self

    def __add__(self, other):
        if QQ.of_type(other):
            return self
        else:
            return NotImplemented

    def __mul__(self, other):
        if QQ.of_type(other):
            return self.mul_q(other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if QQ.of_type(other):
            return self.mul_q(other)
        else:
            return NotImplemented


@total_ordering
class Infinity(Infinite):

    def __repr__(self):
        return 'inf'

    def __eq__(self, other):
        if isinstance(other, Infinity):
            return True
        else:
            return NotImplemented

    def __neg__(self):
        return neginf

    def __gt__(self, other):
        if isinstance(other, Infinity):
            return False
        elif QQ.of_type(other):
            return True
        else:
            return NotImplemented


@total_ordering
class NegativeInfinity(Infinite):

    def __repr__(self):
        return '-inf'

    def __eq__(self, other):
        if isinstance(other, NegativeInfinity):
            return True
        else:
            return NotImplemented

    def __neg__(self):
        return inf

    def __lt__(self, other):
        if isinstance(other, Infinity):
            return False
        elif QQ.of_type(other):
            return True
        else:
            return NotImplemented



inf = Infinity()
neginf = NegativeInfinity()


class RationalInterval:

    def __init__(self, a, b):
        if not self._is_valid(a) and self._is_valid(b):
            raise TypeError('QQ or infinite expected')
        if not a <= b:
            raise ValueError('start and end out of order')
        self.a = a
        self.b = b

    def __repr__(self):
        return f'RationalInterval({repr(self.a)}, {repr(self.b)})'

    def _new(self, a, b):
        return RationalInterval(a, b)

    def _is_valid(self, value):
        return QQ.of_type(value) or isinstance(value, Infinite)

    @property
    def mid(self):
        return (self.a + self.b) / 2

    @property
    def width(self):
        return self.b - self.a

    def n(self, dps=15):
        return [S(self.a).n(dps), S(self.b).n(dps)]

    def is_proper_subset(self, other):
        return self.a > other.a and self.b < other.b

    def intersect(self, other):
        a = max(self.a, other.a)
        b = min(self.b, other.b)
        if a <= b:
            return self._new(a, b)
        else:
            return None

    def contains(self, other):
        return self.a < other < self.b

    def split_zero(self):
        zero = QQ.zero
        if self.contains(0):
            return [self._new(self.a, zero), self._new(zero, self.b)]
        else:
            return [self]

    def add(self, other):
        return self._new(self.a + other.a, self.b + other.b)

    def mul(self, other):
        a1, a2 = self.a, other.a
        b1, b2 = self.b, other.b
        products = [a1*a2, a1*b2, b1*a2, b1*b2]
        a = min(products)
        b = max(products)
        return self._new(a, b)

    def mul_scalar(self, other):
        if other > 0:
            return self._new(other*self.a, other*self.b)
        elif other < 0:
            return self._new(other*self.b, other*self.a)
        elif self.a == neginf and self.b == inf:
            return self._new(neginf, inf)
        elif self.a == neginf:
            return self._new(neginf, QQ.zero)
        elif self.b == inf:
            return self._new(QQ.zero, inf)
        else:
            return self._new(QQ.zero, QQ.zero)

    def add_scalar(self, other):
        return self._new(self.a + other, self.b + other)

    def invert(self):
        a, b = self.a, self.b
        if self.contains(0):
            raise ZeroDivisionError
        elif a != 0 and b != 0:
            return self._new(1/b, 1/a)
        elif a == 0:
            return self._new(1/b, inf)
        elif b == 0:
            return self._new(-inf, 1/a)
        else:
            raise RuntimeError

    def neg(self):
        return self._new(-self.b, -self.a)

    def div(self, other):
        return self.mul(other.invert())

    def div_rscalar(self, other):
        return self.invert().mul_scalar(other)

    def contains(self, other):
        return self.a < other < self.b

    def __add__(self, other):
        if QQ.of_type(other):
            return self.add_scalar(other)
        elif isinstance(self, type(other)):
            return self.add(other)
        else:
            return NotImplemented

    def __radd__(self, other):
        if QQ.of_type(other):
            return self.add_scalar(other)

    def __sub__(self, other):
        if QQ.of_type(other):
            return self.add_scalar(-other)
        elif isinstance(self, type(other)):
            return self.add(other.neg())
        else:
            return NotImplemented

    def __rsub__(self, other):
        if QQ.of_type(other):
            return self.neg().add_scalar(other)
        else:
            return NotImplemented

    def __mul__(self, other):
        if QQ.of_type(other):
            return self.mul_scalar(other)
        elif isinstance(self, type(other)):
            return self.mul(other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if QQ.of_type(other):
            return self.mul_scalar(other)
        return NotImplemented

    def __truediv__(self, other):
        if QQ.of_type(other):
            return self.div_scalar(other)
        elif isinstance(self, type(other)):
            return self.div(other)
        else:
            return NotImplemented

    def __rtruediv__(self, other):
        if QQ.of_type(other):
            return self.div_rscalar(other)
        else:
            return NotImplemented

    def __contains__(self, element):
        if QQ.of_type(element):
            return self.contains(element)
        return NotImplemented


def interval_newton_rat_1step(poly, polydiff, X1):
    x1 = X1.mid
    f = horner(poly, x1)
    Fp = horner(polydiff, X1)

    if QQ.zero not in Fp:
        Z2 = x1 - f/Fp
        X2 = Z2.intersect(X1)

        if X2 is None:
            Xs = []
            status = 'empty'
        else:
            Xs = [X2]
            if Z2.is_proper_subset(X1):
                status = 'unique'
            elif X2 == X1:
                status = 'stalled'
            else:
                status = 'iterate'

    else:
        # The derivative is possibly zero in X1. We split the Fp into a
        # positive interval and a negative interval and proceed with each
        # separately.
        Fp1, Fp2 = Fp.split_zero()
        X21 = (x1 - f/Fp1).intersect(X1)
        X22 = (x1 - f/Fp2).intersect(X1)
        Xs = [Xi for Xi in [X21, X22] if Xi is not None]
        status = 'split'

    return Xs, status


def interval_newton_rat(poly, polydiff, X, prec=53):
    a, b = X
    X = RationalInterval(QQ(a), QQ(b))
    # Iterate over the queue as it extends:
    queue = [X]
    roots = []
    while queue:
        nextqueue = []
        for Xi in queue:
            print(Xi.n())
            Xs, status = interval_newton_rat_1step(poly, polydiff, Xi)
            if status == 'unique':
                [Xj] = Xs
                a, b = Xj.a, Xj.b
                roots.append(Xj)
                continue
                if abs(a - b) < (abs(a)+abs(b)) / 2**prec:
                    roots.append(Xj)
                else:
                    nextqueue.append(Xj)
            else:
                nextqueue.extend(Xs)
        queue = nextqueue
    starts = [X.a for X in queue]
    ends = [X.b for X in queue]
    return roots


def qq_to_iv(num, iv):
    """Convert rational QQ to real MPF interval"""
    return iv.mpf(num.numerator) / iv.mpf(num.denominator)


def poly_to_iv(poly, iv):
    """Convert poly with QQ coefficients to MPF intervals"""
    return [qq_to_iv(coeff, iv) for coeff in poly]


def intersect_iv(X, Y, iv):
    a = max(X.a, Y.a)
    b = min(X.b, Y.b)
    if a <= b:
        return iv.mpf((a, b))
    else:
        return None


def is_proper_subset(X, Y):
    return X.a > Y.a and X.b < Y.b


def split_iv(X, iv):
    return iv.mpf((X.a, X.mid)), iv.mpf((X.mid, X.b))


def split_iv_zero(X, iv):
    assert 0 in X
    return iv.mpf((X.a, 0)), iv.mpf((0, X.b))


def interval_newton(poly, polydiff, X, prec=53):
    iv = MPIntervalContext()
    iv.dps = prec + 20

    poly_iv = poly_to_iv(poly, iv)
    polydiff_iv = poly_to_iv(polydiff, iv)
    X_iv = iv.mpf(X)

    queue = [X_iv]
    while queue:
        if len(queue) == 1:
            [Xj] = queue
            a, b = Xj.a, Xj.b
            if abs(a - b) * 2**(prec+3) < (abs(a) + abs(b)):
                return Float(Xj.mid, precision=prec)

        nextqueue = []

        for X1 in queue:
            x1 = X1.mid
            f = horner(poly_iv, x1)
            Fp = horner(polydiff_iv, X1)

            if 0 not in Fp:
                Z2 = x1 - f/Fp
                X2 = intersect_iv(Z2, X1, iv)

                if X2 is None:
                    pass # no roots
                elif is_proper_subset(Z2, X1):
                    # this interval is the one -> discard all others
                    nextqueue = [X2]
                    break
                elif 4*X2.delta > 3*X1.delta:
                    # converging -> continue
                    nextqueue.append(X2)
                else:
                    # insufficient progress -> split
                    nextqueue.extend(split_iv(X2, iv))

            else:
                # Split at the zero of the denominator
                Fp1, Fp2 = split_iv_zero(Fp, iv)
                X21 = intersect_iv(x1 - f/Fp1, X1, iv)
                X22 = intersect_iv(x1 - f/Fp2, X1, iv)
                X2s = [Xi for Xi in [X21, X22] if Xi is not None]
                nextqueue.extend(X2s)

        queue = nextqueue


def refine_root(p, a, b, dps=15):
    poly = Poly(p).rep.rep
    polydiff = Poly(p.diff()).rep.rep
    prec = dps_to_prec(dps)
    return interval_newton(poly, polydiff, (a, b), prec)



x = Symbol('x')
#p = (x-1)*(x-2)*(x-S(3)/2)*(x-4)*(x-7)
p = x**2 - 2
roots = refine_root(p, 0, 1.5)
#poly = Poly(p, domain=QQ).rep.rep
#polydiff = Poly(p.diff(), domain=QQ).rep.rep
#roots = interval_newton(poly, polydiff, [0, 1.5])
#print(len(roots), 'roots')
#for r in roots:
    #print(r.n())
print(roots)


def evalf_crootof(root, dps=None):
    prec = dps_to_prec(dps)

    if not root.is_real:
        raise NotImplementedError('non-real roots')

    # QQ arithmetic:
    interval = root._get_interval()
    a = interval.a
    b = interval.b
    if 2**prec * abs(a - b) < abs(a) + abs(b):
        r = (a + b) / 2
        return Float(r.numerator, precision=prec) / r.denominator

    poly = root.poly.rep.rep
    polydiff = diff(poly)
    #a = mpf(interval.a.numerator) / mpf(interval.a.denominator)
    #b = mpf(interval.b.numerator) / mpf(interval.b.denominator)
    rapprox, converged = newton_mpf(poly, polydiff, a, b, prec+10)
    sign, man, exp, bc = rapprox._mpf_
    if exp > 0:
        num, den = man*2**exp, 1
    else:
        num, den = man, 2**-exp
    sign = (-1)**sign
    rapprox_q = QQ(sign*num, den)
    #if not (a < rapprox_q < b):
    #    raise ConvergenceError
    return Float(rapprox, precision=prec)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
