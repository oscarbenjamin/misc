from __future__ import annotations

from collections import deque
from typing import overload, Iterable, Optional, Dict

from sympy.core import Basic, Expr, Mul, Pow
from sympy.multipledispatch import dispatch

from sympy.physics.quantum import BraBase, KetBase, InnerProduct, OuterProduct


_QOpt = Dict[str, bool]


def qapply(q_expr: Basic, *, dagger=False, ip_doit=True) -> Basic:
    """Apply operators to states.

    >>> from sympy.physics.quantum import Bra, Ket
    >>> b1, b2, k1, k2 = Bra(1), Bra(2), Ket(1), Ket(2)
    >>> expr = b1*b2*k1*k2
    >>> expr
    <1|*<2|*|1>*|2>
    >>> qapply(expr)
    <1|2>*<2|1>
    >>> expr2 = k1*k2*b1*b2
    >>> expr2
    |1>*|2>*<1|*<2|
    >>> qapply(expr2)
    |1>*|2><1|*<2|
    """
    options: _QOpt = {'dagger': dagger, 'ip_doit': ip_doit}
    return _qapply(q_expr, options)


@overload
def _qapply(q_expr: Expr, options: _QOpt) -> Expr:
    ...


@overload
def _qapply(q_expr: Basic, options: _QOpt) -> Basic:
    ...


def _qapply(q_expr: Basic, options: _QOpt) -> Basic:
    """Inner _qapply for recursion."""
    if q_expr.is_Atom:
        return q_expr
    elif isinstance(q_expr, Mul):
        return _qapply_mul(q_expr, options)
    elif isinstance(q_expr, Pow):
        return _qapply_pow(q_expr, options)
    else:
        newargs = [_qapply(arg, options) for arg in q_expr.args]
        return q_expr.func(*newargs)


def _qapply_pow(q_expr: Pow, options: _QOpt) -> Expr:
    """qapply a Pow"""
    return _qapply_mul_args([q_expr], options)


def _qapply_mul(q_expr: Mul, options: _QOpt) -> Expr:
    """qapply a Mul"""
    return _qapply_mul_args(list(q_expr.args), options)


def _qapply_mul_args(args: list[Expr], options: _QOpt) -> Expr:
    """qapply(Mul(*args))

    It is assumed that lhs has already been processed.

    >>> from sympy.physics.quantum import Bra, Ket
    >>> b1, b2, k1, k2 = Bra(1), Bra(2), Ket(1), Ket(2)
    >>> options = {}
    >>> _qapply_mul_args([b1, b2, b1, b2], options)
    <1|*<2|*<1|*<2|
    >>> _qapply_mul_args([b1, b2, k1, k2], options)
    <1|2>*<2|1>
    """
    cpart: list[Expr] = []
    ncpart: list[Expr] = []

    for arg in args:
        cpart_i, ncpart_i = _qapply_args_cnc(arg, options)
        cpart.extend(cpart_i)
        for arg_i in ncpart_i:
            _qapply_mul_nc_append(cpart, ncpart, arg_i, options)

    return Mul(*(cpart + ncpart))


def _qapply_mul_nc_append(cpart: list[Expr],
                          ncpart: list[Expr],
                          arg: Expr,
                          options: _QOpt,
                          ) -> None:
    """qapply(Mul(*(cpart + ncpart)) * rhs) in-place.

    >>> from sympy.physics.quantum import Bra, Ket
    >>> b1, b2, k1, k2 = Bra(1), Bra(2), Ket(1), Ket(2)
    >>> options = {}
    >>> cpart, ncpart = [], [b1, b2]
    >>> cpart, ncpart
    ([], [<1|, <2|])
    >>> _qapply_mul_nc_append(cpart, ncpart, b1, options)
    >>> cpart, ncpart
    ([], [<1|, <2|, <1|])
    >>> _qapply_mul_nc_append(cpart, ncpart, k1, options)
    >>> cpart, ncpart
    ([<1|1>], [<1|, <2|])
    """
    while ncpart:
        result = _qapply_mul_nc_2args(ncpart[-1], arg, options)
        if result is None:
            break
        else:
            ncpart.pop()
        if result.is_commutative:
            cpart.append(result)
            return
        arg = result
    ncpart.append(arg)


def _qapply_mul_nc_2args(lhs: Expr, rhs: Expr, options: _QOpt
                         ) -> Optional[Expr]:
    """qapply(lhs * rhs) -> result, interacted

    >>> from sympy.physics.quantum import Bra, Ket
    >>> b1, b2, k1, k2 = Bra(1), Bra(2), Ket(1), Ket(2)
    >>> options = {}
    >>> print(_qapply_mul_nc_2args(b1, b2, options))
    None
    >>> print(_qapply_mul_nc_2args(b1, k2, options))
    <1|2>
    """
    try:
        return qmultiply(lhs, rhs)
    except NotImplementedError:
        return None


def _qapply_args_cnc(q_expr: Expr, options: _QOpt
                     ) -> tuple[list[Expr], list[Expr]]:
    """split commutative and non-commutative args and qapply"""
    #
    # Needs expansion...
    #
    return q_expr.args_cnc()


@dispatch(BraBase, KetBase)
def qmultiply(b, k):
    return InnerProduct(b, k)


@dispatch(KetBase, BraBase)
def qmultiply(k, b):
    return OuterProduct(k, b)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
