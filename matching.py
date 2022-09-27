from sympy import sympify, Add, Pow, Mul, Symbol, Integer, parse_expr
from matchpy import match, Operation, Arity, Pattern, Wildcard
from matchpy.utils import Multiset
from matchpy import Symbol as mSymbol

mAdd = Operation.new('Add', Arity.polyadic, associative=True, commutative=True)
mMul = Operation.new('Mul', Arity.polyadic, associative=True, commutative=True)
mPow = Operation.new('Pow', Arity.binary)

s2m_map = {
    Add: mAdd,
    Mul: mMul,
    Pow: mPow,
}
m2s_map = {v:k for k, v in s2m_map.items()}

def s2m(expr, sequence=False):
    """sympy expression -> matchpy expression"""
    expr = sympify(expr)
    t = type(expr)
    if t in s2m_map:
        is_sequence = t in (Add, Mul)
        args = [s2m(arg, is_sequence) for arg in expr.args]
        return s2m_map[t](*args)
    elif isinstance(expr, Wild):
        name = expr.name
        if sequence:
            return Wildcard.star(name)
        else:
            return Wildcard.dot(name)
    else:
        return mSymbol(str(expr))

def m2s(mexpr):
    """matchpy expression -> sympy expression"""
    t = type(mexpr)
    if t in m2s_map:
        args = [m2s(arg) for arg in mexpr.operands]
        return m2s_map[t](*args)
    elif isinstance(mexpr, Wildcard):
        name = mexpr.variable_name
        return Wild(name)
    elif isinstance(mexpr, Multiset):
        results = [m2s(res) for res in mexpr]
        if len(results) == 1:
            return results[0]
        else:
            return results
    else:
        return parse_expr(mexpr.name)

def mpy_match(subject, pattern):
    """match sympy expression using matchpy internally"""
    msubject = s2m(subject)
    mpattern = Pattern(s2m(pattern))
    msubstitution = next(match(msubject, mpattern), None)
    if msubstitution is None:
        return
    else:
        substitution = {Wild(k):m2s(v) for k, v in msubstitution.items()}
        return substitution
