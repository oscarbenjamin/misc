# Complexity Results for Fourier-Motzkin Elimination
# Delaram Talaashrafi, The University of Western Ontario
#
# https://ir.lib.uwo.ca/cgi/viewcontent.cgi?article=8167&context=etd

from sympy import Le, Ge, linear_eq_to_matrix, Matrix

def fme(inequalities, syms, elim):
    """
    Fourier-Motzkin Elimination

    >>> from sympy import symbols
    >>> x, y, z = symbols('x, y, z')
    >>> inequalities = [x + 2*y + 3*z >= 1, 2*x + 3*y - 4*z >= 5]
    >>> fme(inequalities, [x, y, z], [z])
    [10*x + 17*y <= 19]
    """
    inequalities = list(inequalities)
    syms = list(syms)
    elim = list(elim)
    if not (inequalities and syms and elim):
        return inequalities

    syms = [s for s in syms if s not in elim] + elim

    exprs = []
    for i in inequalities:
        if isinstance(i, Le):
            exprs.append(i.lhs - i.rhs)
        elif isinstance(i, Ge):
            exprs.append(i.rhs - i.lhs)
        else:
            raise ValueError("inequalities must be a list of Le or Ge")

    A, b = linear_eq_to_matrix(exprs, syms)

    elim_indices = range(len(syms) - len(elim), len(syms))

    for i in elim_indices:
        A, b = fme1(A, b, i)

    lhs_sides = A*Matrix(syms)
    inequalities = [lhs_sides[i] <= b[i] for i in range(len(b))]
    # How to handle False?
    inequalities = [i.canonical for i in inequalities if i not in [True, False]]

    return inequalities


def fme1(A, b, i):
    """
    Fourier-Motzkin Elimination for one variable
    """
    S_plus = []
    S_minus = []
    S_zero = []
    for j in range(A.rows):
        if A[j, i] > 0:
            S_plus.append(j)
        elif A[j, i] < 0:
            S_minus.append(j)
        else:
            S_zero.append(j)

    S_nonzero = []
    for j in S_plus:
        for k in S_minus:
            S_nonzero.append(combine(A, b, i, j, k))

    A_nonzero = Matrix([s[0] for s in S_nonzero])
    b_nonzero = Matrix([s[1] for s in S_nonzero])

    A_zero = A[S_zero, :]
    b_zero = b[S_zero, :]

    A_new = Matrix.vstack(A_nonzero, A_zero)
    b_new = Matrix.vstack(b_nonzero, b_zero)

    return A_new, b_new


def combine(A, b, i, j, k):
    """
    Combine two rows of A and b
    """
    Aijk = -A[k, i]*A[j, :] + A[j, i]*A[k, :]
    bij = -A[k, i]*b[j, 0] + A[j, i]*b[k, 0]
    return Aijk, bij


if __name__ == "__main__":
    import doctest
    doctest.testmod()
