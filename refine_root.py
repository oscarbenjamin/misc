from sympy import *
from sympy.core.evalf import dps_to_prec
from mpmath import MPIntervalContext


def qq_to_iv(num, iv):
    """Convert rational QQ to real MPF interval"""
    return iv.mpf(num.numerator) / iv.mpf(num.denominator)


def poly_to_iv(poly, iv):
    """Convert poly with QQ coefficients to MPF intervals"""
    return [qq_to_iv(coeff, iv) for coeff in poly]


def horner(poly, x0):
    """Evaluate a poly

    >>> horner([1, 0, -1], 3)
    8
    """
    val = poly[0]
    for c in poly[1:]:
        val = c + val*x0
    return val


def eval_poly_iv(poly, x0, iv):
    return horner(poly_to_iv(poly, iv), x0)


def intersect_iv(X, Y, iv):
    a = max(X.a, Y.a)
    b = min(X.b, Y.b)
    if a <= b:
        return iv.mpf((a, b))
    else:
        return None


def subtract_iv(X, Y, iv):
    if X.b <= Y.a or X.a >= Y.b:
        # no overlap
        return [X]
    elif Y.a <= X.a and X.b <= Y.b:
        # X is a subset of Y (or equal)
        return []
    elif X.a <= Y.a and X.b <= Y.b:
        # Y overlaps right of X
        return [iv.mpf((X.a, Y.a))]
    elif Y.a <= X.a and Y.b <= X.b:
        # Y overlaps left of X
        return [iv.mpf((X.a, Y.a))]
    else:
        # Y is a strict subset of X
        return [iv.mpf((X.a, Y.a)), iv.mpf((Y.b, X.b))]


def is_proper_subset(X, Y):
    return X.a > Y.a and X.b < Y.b


def split_iv(X, iv):
    return iv.mpf((X.a, X.mid)), iv.mpf((X.mid, X.b))


def split_iv_zero(X, iv):
    assert 0 in X
    return iv.mpf((X.a, 0)), iv.mpf((0, X.b))

def interval_newton(poly, polydiff, X, prec=53):
    iv = MPIntervalContext()
    iv.prec = prec + 20

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
            F = horner(poly_iv, X1)
            if not 0 in F:
                continue

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
                elif 4*X2.delta < 3*X1.delta:
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
                for X2 in X21, X22:
                    if X2 is None:
                        continue
                    elif 4*X2.delta < 3*X1.delta:
                        nextqueue.append(X2)
                    else:
                        nextqueue.extend(split_iv(X2, iv))

        queue = nextqueue


def interval_newton_1step(f, F, X1, iv):
    x1 = X1.mid
    f1 = f(x1)
    FX1 = F(X1)

    # Split at the zero of the denominator
    if FX1.a < 0 < FX1.b:
        FX1s = [iv.mpf((FX1.a, 0)), iv.mpf((0, FX1.b))]
    else:
        FX1s = [FX1]

    X2s = []
    for FXi in FX1s:
        Z2i = x1 - f1/FX1
        X2i = intersect_iv(Z2i, X1, iv)
        if X2i is not None:
            X2s.append(X2i)

    return X2s


class FPIntervalContext:
    def mpf(self, ab):
        return fpi(ab)


class fpi:

    inf = float('inf')
    ninf = float('-inf')

    def __new__(cls, ab):
        if isinstance(ab, tuple):
            ar, br = ab
            if isinstance(ar, fpi):
                a = ar.a
            else:
                a = float(ar)
            if isinstance(br, fpi):
                b = br.b
            else:
                b = float(br)
        else:
            a = b = float(ab)
        return cls._new(a, b)

    @classmethod
    def _new(cls, a, b):
        assert isinstance(a, float)
        assert isinstance(b, float)
        obj = super().__new__(cls)
        obj.a = a
        obj.b = b
        return obj

    def __repr__(self):
        return f'fp({self.a, self.b})'

    def __eq__(self, other):
        if isinstance(other, fpi):
            return (self.a, self.b) == (other.a, other.b)
        else:
            return NotImplemented

    @property
    def mid(self):
        m = (self.a + self.b) / 2
        return self._new(m, m)

    @property
    def delta(self):
        d = self.b - self.a
        return self._new(d, d)

    def __abs__(self):
        b = max(abs(self.a), abs(self.b))
        return self._new(0.0, b)

    def __add__(self, other):
        if isinstance(other, fpi):
            return self.add(other)
        elif isinstance(other, float):
            return self.add(self._new(other, other))
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, float):
            return self.add(self._new(other, other))
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, fpi):
            return self.sub(other)
        elif isinstance(other, float):
            return self.sub(self._new(other, other))
        else:
            return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, float):
            return self.sub(self._new(other, other))
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, fpi):
            return self.mul(other)
        elif isinstance(other, float):
            return self.mul(self._new(other, other))
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, float):
            return self.mul(self._new(other, other))
        else:
            return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, fpi):
            return self.div(other)
        elif isinstance(other, (float, int)):
            return self.div(fpi(other))
        else:
            return NotImplemented

    def __gt__(self, other):
        if isinstance(other, fpi):
            return self.gt(other)
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, fpi):
            return self.ge(other)
        else:
            return NotImplemented

    def __contains__(self, value):
        return self.a <= value <= self.b

    def add(self, other):
        return self._new(self.a + other.a, self.b + other.b)

    def sub(self, other):
        return self._new(self.a - other.b, self.b - other.a)

    def mul(self, other):
        x1, x2, y1, y2 = self.a, self.b, other.a, other.b
        if x1 == x2 == 0 and y1 == self.ninf and y2 == self.inf:
            return self._new(self.ninf, self.inf)
        if y1 == y2 == 0 and x1 == self.ninf and x2 == self.inf:
            return self._new(self.ninf, self.inf)
        combinations = [x1*y1, x1*y2, x2*y1, x2*y2]
        return self._new(min(combinations), max(combinations))

    def inverse(self):
        a, b = self.a, self.b
        if not (a <= 0 <= b):
            return self._new(1/b, 1/a)
        elif a == 0 and b != 0:
            return self._new(1/b, self.inf)
        elif b == 0 and a != 0:
            return self._new(self.ninf, 1/a)
        else:
            return self._new(self.ninf, self.inf)

    def div(self, other):
        return self.mul(other.inverse())

    def gt(self, other):
        return self.a > other.b

    def ge(self, other):
        return self.a >= other.b


class PolyEvaluator:

    def __init__(self, poly):
        self.poly = Poly(p).rep.rep
        self.polydiff = Poly(p.diff(), Poly(p).gens).rep.rep
        self._precs = {}

    def _prec(self, prec, usefp):
        if usefp and prec == 53:
            iv = FPIntervalContext()
        else:
            iv = MPIntervalContext()
            iv.prec = prec
        polycoeffs = poly_to_iv(self.poly, iv)
        polydiffcoeffs = poly_to_iv(self.polydiff, iv)

        def f(X):
            return horner(polycoeffs, X)

        def fp(X):
            return horner(polydiffcoeffs, X)

        return iv, f, fp

    def get(self, prec, use=True):
        X, prec = prec
        key = (use, prec)
        try:
            iv, f, fp = self._precs[key]
        except KeyError:
            iv, f, fp = self._precs[key] = self._prec(prec, usefp=use)

        if isinstance(X, fpi) and isinstance(iv, MPIntervalContext):
            X = iv.mpf((X.a, X.b))

        if isinstance(X, tuple):
            X = iv.mpf(X)
        return iv, f, fp, X, prec


# https://hal.inria.fr/inria-00072253/document
def isolate_real_roots(poly, X0, maxprec, use=True):

    def f(X):
        iv.prec, old_prec = prec, iv.prec
        result = eval_poly_iv(poly, X, iv)
        iv.prec = old_prec
        return result

    def fp(X):
        iv.prec, old_prec = prec, iv.prec
        result = eval_poly_iv(polydiff, X, iv)
        iv.prec = old_prec
        return result

    def prec_insufficient(Xi, X):
        if Xi is None:
            return False
        return Xi.delta >= X.delta or f(Xi).delta >= f(X).delta

    def add_root(Xj, prec):
        # Add root to the output list
        roots.append((Xj, prec))
        # Remove this interval from all intervals in the stack
        newL = []
        for Xi, prec in L:
            Xis = subtract_iv(Xi, Xj, iv)
            for Xip in Xis:
                newL.append((Xip, prec))
        L[:] = newL

    def precision(Xi):
        a, b = Xi.a, Xi.b
        delta = abs(a - b)
        if delta == 0:
            return float('inf')
        return int((abs(a)+abs(b)/delta).mid).bit_length()

    #iv = MPIntervalContext()
    pe = PolyEvaluator(poly)
    iv, f, fe, X0, prec = pe.get((X0, 53), use=use)
    alpha = 0.75
    L = [(X0, 64)]
    roots = []

    while L:
        iv, f, fp, X, prec = pe.get(L.pop(), use=use)

        #fXa = abs(f(X.a))
        #fXb = abs(f(X.b))
        #if fXa <= 1e-10*fXb:
        #    X1 = iv.mpf((X.a-X.delta/999, X.a+X.delta/999))
        #    X2 = iv.mpf((X.a+X.delta/999, X.b))
        #    L.append((X1, prec))
        #    X = X2
        #elif fXb <= 1e-10*fXa:
        #    X1 = iv.mpf((X.a, X.b-X.delta/999))
        #    X2 = iv.mpf((X.b-X.delta/999, X.b+X.delta/999))
        #    L.append((X1, prec))
        #    X = X2

        Xs = interval_newton_1step(f, fp, X, iv)

        if any(Xi.delta > alpha*X.delta for Xi in Xs):
            Xs = iv.mpf((X.a, X.mid)), iv.mpf((X.mid, X.b))

        if any(prec_insufficient(Xi, X) for Xi in Xs):
            newprec = 2*prec
        else:
            newprec = prec

        for Xi in Xs:
            if 0 in f(Xi):
                if is_proper_subset(Xi, X) and precision(Xi) >= maxprec:
                    add_root(Xi, prec)
                else:
                    L.append((Xi, newprec))

    return [Float(r.mid, precision=maxprec) for r, _ in roots][::-1]


def refine_root(p, a, b, dps=15, use=False):
    poly = Poly(p).rep.rep
    polydiff = Poly(p.diff(), Poly(p).gens).rep.rep
    prec = dps_to_prec(dps)
    return isolate_real_roots(p, (a, b), prec, use=use)


#x = Symbol('x')
#p = x**2 - 2
#root = refine_root(p, 0, 1.5)
#print(root)
