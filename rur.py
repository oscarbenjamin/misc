from sympy import *


def solve_poly_numeric(polys, syms, exact=False, prec=None):
    """Solve a system of polynomials having rational coefficients."""
    _, info = parallel_poly_from_expr(polys)
    domain = info['domain']
    if domain not in (ZZ, QQ):
        raise ValueError("Poly should have rational coefficients")

    # Compute a preliminary Groebner basis
    gb = groebner(polys, syms)

    # Handle inconsistent or infinite cases
    if 1 in gb:
        return []
    elif not gb.is_zero_dimensional:
        raise ValueError("Infinitely many solutions")

    # Split the system by factorising the final polynomial
    c, fms = factor_list(gb[-1])
    gbs = []
    for factor, m in fms:
        gb_new = groebner(gb[:-1] + [factor], syms)
        gbs.append(gb_new)

    # Now solve each subsystem
    solutions = []
    for gbi in gbs:
        solutions.extend(solve_separating(gbi))

    # Make the solutions approximate (this is because otherwise you'll see
    # complicated RootOf expressions).
    if not exact:
        solutions = [[s.evalf(prec) for s in sol] for sol in solutions]

    return solutions


def solve_separating(gb):
    syms = gb.gens
    N = len(syms)
    s = Dummy('s')
    i = 0
    while True:
        eq_s = s - sum(j**i*syms[j] for j in range(N))
        gb = groebner(list(gb) + [eq_s], syms + (s,))
        if is_separated(gb):
            return solve_rur(gb)
        i += 1


def is_separated(gb):
    """Test if a Groebner basis is separated"""
    for p in gb.polys[:-1]:
        if sum(p.degree_list()[:-1]) != 1:
            return False
    return sum(gb.polys[-1].degree_list()[:-1]) == 0


def solve_rur(gb):
    [sol] = linsolve(gb[:-1], gb.gens[:-1])
    s = gb.gens[-1]
    s_sols = set(gb.polys[-1].as_poly(s).real_roots())
    return [sol.subs(s, s_sol) for s_sol in s_sols]



def rand_poly_system(syms, degree, sparseness=1):
    monomials = list(itermonomials(syms, degree))
    M = len(monomials)
    N = len(syms)
    coeffs = randMatrix(N, M, -10, 10, percent=sparseness*100)
    return list(coeffs * Matrix(monomials))
