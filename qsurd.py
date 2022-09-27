#
#
#
#
#
#
# {{(1, 2), (4, 5)}
class SMonomial(frozenset):

    @classmethod
    def _new(cls, powermap):
        return super().__new__(cls, powermap)

    def __repr__(self):
        return f'SMonomial({set(self)})'

    def as_dict(self):
        return dict(self)

    def __mul__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.mul(other)

    def __pow__(self, exponent):
        if not isinstance(exponent, int):
            return NotImplemented
        else:
            return self.pow(exponent)

    def mul(self, other):
        powermap = dict(self)
        for g, n in other:
            other_n = powermap.get(g)
            if other_n is None:
                powermap[g] = n
            else:
                powermap_n = other_n + n
                if powermap_n:
                    powermap[g] = powermap_n
                else:
                    powermap.pop(g)
        return self._new(powermap.items())

    def pow(self, exponent):
        if exponent < 0:
            raise ValueError("Negative power")
        elif exponent == 0:
            return self.one
        else:
            powermap = [(g, m*exponent) for g, m in self]
            return self._new(powermap)

    def is_one(self):
        return self == self.one


class SPoly(dict):

    def __init__(self, arg, ring):
        dict.__init__(self, arg)
        self.ring = ring

    @classmethod
    def _new(cls, arg, domain):
        return cls(arg, domain)

    def __repr__(self):
        if not self:
            return '0'

        terms = []
        symbols = self.ring.symbols
        for m, c in sorted(self.items()):
            if m.is_one():
                term = str(c)
            else:
                if c != 1:
                    factors = [str(c)]
                else:
                    factors = []
                for s, e in m.as_dict().items():
                    if e != 1:
                        factors.append(f'{symbols[s]}**{e}')
                    else:
                        factors.append(f'{symbols[s]}')
                term = '*'.join(factors)
            terms.append(term)
        return ' + '.join(terms)

    def __add__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        elif self.ring != other.ring:
            raise ValueError("Domain mismatch")
        return self.add(other)

    def __sub__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        elif self.ring != other.ring:
            raise ValueError("Domain mismatch")
        return self.sub(other)

    def __neg__(self):
        return self.neg()

    def __mul__(self, other):
        if isinstance(other, type(self)):
            if self.ring != other.ring:
                raise ValueError("Domain mismatch")
            return self.pmul(other)
        elif self.ring.domain.of_type(other) or isinstance(other, int):
            return self.cmul(other)
        elif isinstance(other, SMonomial):
            return self.mmul(other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if self.ring.domain.of_type(other) or isinstance(other, int):
            return self.cmul(other)
        elif isinstance(other, SMonomial):
            return self.mmul(other)
        else:
            return NotImplemented

    def __pow__(self, exponent):
        if not isinstance(exponent, int):
            return NotImplemented
        else:
            return self.pow(exponent)

    def add(self, other):
        coeffmap = dict(self)
        for m, c in other.items():
            self_c = coeffmap.get(m)
            if self_c is None:
                coeffmap[m] = c
            else:
                new_c = self_c + c
                if new_c:
                    coeffmap[m] = new_c
                else:
                    coeffmap.pop(m)
        return self._new(coeffmap, self.ring)

    def sub(self, other):
        return self.add(other.neg())

    def neg(self):
        return self._new({m: -c for m, c in self.items()}, self.ring)

    def cmul(self, coeff):
        return self._new({m: c*coeff for m, c in self.items()}, self.ring)

    def mmul(self, monomial):
        return self._new({m*monomial: c for m, c in self.items()}, self.ring)

    def pmul(self, other):
        coeffmap = {}
        for m1, c1 in self.items():
            for m2, c2 in other.items():
                m3 = m1.mul(m2)
                c3 = c1 * c2
                coeff3 = coeffmap.get(m3)
                if coeff3 is None:
                    coeffmap[m3] = c3
                else:
                    c3 += coeff3
                    if c3:
                        coeffmap[m3] = c3
                    else:
                        coeffmap.pop(m3)
        return self._new(coeffmap, self.ring)

    def pow(self, exponent):
        # Iterated multiplication is faster than iterated squaring:
        #
        # Journal of Computer and System Sciences
        # Volume 6, Issue 1, February 1972, Pages 1-8
        # Journal of Computer and System Sciences
        # Computation of powers of multivariate polynomials over the integers
        # Lee E.Heindel*
        if exponent == 0:
            return self.one
        elif exponent < 0:
            raise ValueError("Negative power")
        else:
            prod = self
            for n in range(1, exponent):
                prod = prod.pmul(self)
            return prod


class SDomain:

    def __init__(self, domain, symbols):
        self.domain = domain
        self.symbols = symbols
        self.zero = SPoly({}, self)
        self.one = SPoly({SMonomial.one: domain.one}, self)
        self.ngens = len(symbols)
        self.monomial_basis = [SMonomial([(n, 1)]) for n in range(self.ngens)]
        self.monomial_map = {s:m for s, m in zip(symbols, self.monomial_basis)}

    def __repr__(self):
        gens = ','.join(map(str, self.symbols))
        return f'{self.domain}[{gens}]'

    def __call__(self, expr):
        return self.from_sympy(expr)

    def from_sympy(self, expr):
        if expr in self.monomial_map:
            monom = self.monomial_map[expr]
            return SPoly({monom: self.domain.one}, self)
        elif isinstance(expr, Mul):
            return prod([self.from_sympy(arg) for arg in expr.args])
        elif isinstance(expr, Add):
            return sum([self.from_sympy(arg) for arg in expr.args], self.zero)
        elif isinstance(expr, Pow) and isinstance(expr.exp, Integer) and expr.exp > 0:
            return self.from_sympy(expr.base) ** expr.exp.p
        else:
            k = self.domain.from_sympy(expr)
            return k * self.one

    def to_sympy(self, poly):
        terms = []
        for m, c in poly.items():
            powers = [Pow(self.symbols[i], e) for i, e in monom]
            coeff = self.domain.to_sympy(c)
            term = Mul(*([coeff] + powers))
            terms.append(term)
        return Add(*terms)



SMonomial.one = SMonomial()
SPoly.one = SPoly({SMonomial.one:1}, ZZ)


x0, x1, x2, x3 = [SMonomial([(n, 1)]) for n in range(4)]
ZZxy = SDomain(ZZ, symbols('x:4'))

m1 = x1**2*x3**4
m2 = x0
p1 = SPoly({m1:2,m2:3}, ZZxy)
assert str(m1) == repr(m1) == 'x1**2*x3**4'
assert str(m2) == repr(m2) == 'x0'
assert str(p1) == repr(p1) == '2*x1**2*x3**4 + 3*x0'

m3 = m1 * m2
assert str(m3) == repr(m3) == 'x0*x1**2*x3**4'

assert (x0 != x1) is True
assert (x0 == x1) is False
assert (x0 != x0) is False
assert (x0 == x0) is True

assert m1*m2 == SMonomial([(0, 1), (1, 2), (3, 4)])
assert m2**2 == SMonomial([(0, 2)])
assert m1**2 == SMonomial([(1, 4), (3, 8)])
assert m1**0 == SMonomial([]) == SMonomial.one

p2 = SPoly({m3:4}, ZZxy)
assert -p1 == SPoly({m1:-2,m2:-3}, ZZxy)
assert p1 + p1 == SPoly({m1:4,m2:6}, ZZxy)
assert p1 + p2 == SPoly({m1:2,m2:3, m3:4}, ZZxy)
assert p1 - p2 == SPoly({m1:2,m2:3, m3:-4}, ZZxy)
assert p1 - p1 == SPoly({}, ZZxy)

assert SPoly({m2:2}, ZZxy) * SPoly({m2:4}, ZZxy) == SPoly({m2**2: 8}, ZZxy)
