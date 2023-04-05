from sympy import *


def linsolve_cond(eqs, unknowns, unique=False):
    if not eqs:
        return S.Complexes**len(unknowns)

    # Preprocessing
    A, b = linear_eq_to_matrix(eqs, unknowns)
    Aaug = Matrix.hstack(A, b).tolist()

    # Main workhorse:
    sols_conds = _linsolve_cond(Aaug)

    # sols_conds is a list of 3-tuples:
    #   [(solset, pivot_conds, consistency_conds),...]
    #
    # solset: solution set as a FiniteSet or ImageSet
    # pivot_conds: list of conditions (e.g. a!=0) assumed in pivoting
    # consistency_conds: list of conditions needed for existence of solutions

    # Build all the separate cases into a Piecewise:
    sets_conds = []
    for solset, pivot_conds, consistency_conds in sols_conds:
        pivot_cond = And(*pivot_conds)
        consistency_cond = And(*consistency_conds)
        if consistency_cond is not S.false:
            sets_conds.append((solset, pivot_cond & consistency_cond))
        if consistency_cond is not S.true:
            sets_conds.append((S.EmptySet, pivot_cond & Not(consistency_cond)))

    sets_conds_d = {}
    for ss, conds in sets_conds:
        if ss not in sets_conds_d:
            sets_conds_d[ss] = conds
        else:
            sets_conds_d[ss] = Or(sets_conds_d[ss], conds)

    if unique:
        sets_conds_d = {s: c for s, c in sets_conds_d.items() if isinstance(s, FiniteSet)}

    return Piecewise(*sets_conds_d.items())


def _linsolve_cond(Aaug, _recurse=None):
    Nr, Nc = len(Aaug), len(Aaug[0])

    Aorig = Matrix(Aaug)

    if _recurse is None:
        row, col, pivots, pivot_conds = 0, 0, [], []
    else:
        row, col, pivots, pivot_conds = _recurse

    if pivots:
        row, col = pivots[-1]
        row += 1
        col += 1
    else:
        row, col = 0, 0

    sols_conds = []

    # Call self recursively for alternate pivots
    def recurse_zero_pivot(r, c):
        pivot = Aaug[r][c]
        Aaugr = [[Arc.subs(pivot, 0) for Arc in Arow] for Arow in Aaug]
        pivot_condsr = pivot_conds[:] + [Eq(pivot, 0)]
        _recurse = (r, c, pivots[:], pivot_condsr)
        sols_conds.extend(_linsolve_cond(Aaugr, _recurse=_recurse))

    while row < Nr and col < Nc-1:
        # Find pivot row and swap into position
        for r in range(row, Nr):
            is_zero = Aaug[r][col].is_zero
            if not is_zero:
                if is_zero is None:
                    # Recurse for the case that the pivot is zero
                    recurse_zero_pivot(r, col)
                    pivot_conds.append(Ne(Aaug[r][col], 0))
                if r != row:
                    Aaug[r], Aaug[row] = Aaug[row], Aaug[r]
                break
        else:
            # All zeros, next column
            col += 1
            continue

        if pivots:
            assert pivots[-1][0] != row

        pivots.append((row, col))
        pivot_row = Aaug[row]
        pivot_div = Aaug[row][col]
        for r in range(row+1, Nr):
            pivot_mul = Aaug[r][col]
            if pivot_mul.is_zero:
                continue
            Aaug[r][col] = S.Zero
            for c in range(col+1, Nc):
                Aaug[r][c] = Aaug[r][c]*pivot_div - pivot_row[c]*pivot_mul

        # Next row/column...
        row += 1
        col += 1

    # Back substitute and list of possibilities
    sol_set, consistency_conds = _back_substitute(Aaug, pivots)

    sols_conds.append((sol_set, pivot_conds, consistency_conds))

    return sols_conds


def _back_substitute(Aaug, pivots):
    Nc = len(Aaug[0])

    # Check conditions for existence of solutions then find solutions by
    # back-substitution below
    consistency_conds = []
    for row in reversed(range(len(Aaug))):
        is_zero = [e.is_zero for e in Aaug[row]]
        if not all(x is True for x in is_zero[:-1]):
            break
        elif is_zero[-1] is False:
            consistency_conds.append(S.false)
        elif is_zero[-1] is None:
            consistency_conds.append(Eq(Aaug[row][-1], 0))

    assert (row == 0 and not pivots) or row == pivots[-1][0]

    # Matrix of all zeros?
    if not pivots:
        solset = S.Complexes**(Nc-1)
        return solset, consistency_conds

    gen = numbered_symbols('tau')
    params = []
    sol = [None] * (Nc-1)

    pivots_cols = {c:r for r, c in pivots}
    for col in reversed(range(Nc-1)):
        if col in pivots_cols:
            r = pivots_cols[col]
            lhsterms = (Aaug[r][c]*sol[c] for c in range(col+1, Nc-1))
            sol[col] = (Aaug[r][-1] - Add(*lhsterms)) / Aaug[r][col]
        else:
            # Non-pivot gets a free symbol
            sym = next(gen)
            params.append(sym)
            sol[col] = sym

    if params:
        solset = ImageSet(Lambda(tuple(params), tuple(sol)), *[S.Complexes]*len(params))
    else:
        solset = FiniteSet(tuple(sol))

    return solset, consistency_conds


x, y, z = symbols('x, y, z')
a, b, c, d, e = symbols('a, b, c, d, e', finite=True)

unknowns = (x, y, z)
eqs = [sqrt(3)*x+y, sqrt(2)*z]

sol = linsolve_cond(eqs, unknowns)
pprint(sol)


M = Matrix(symbols('M:9')).reshape(3, 3)
xs = Matrix(symbols('x:3'))
b = Matrix(symbols('b:3'))
sol3 = linsolve_cond(list(M*xs - b), list(xs), unique=True)
