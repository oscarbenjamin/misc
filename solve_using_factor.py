from itertools import product

def factors(expr):
    return [f for f, m in factor_list(expr)[1]]


def solve_factor(eqs, syms):
    eq_factors = [factors(eq) for eq in eqs]
    systems = list(product(*eq_factors))
    solutions = set()
    for system in systems:
        sols = nonlinsolve(system, syms)
        solutions.update(sols)
        print(sols)
    return list(ordered(solutions))


x, y, z, w = syms = symbols('x, y, z, w')

eqs = [2*x*(w + 1), 2*y*(w - 1), 2*z*(w - 1), x**2 + y**2 + z**2 - 1]

print(solve_factor(eqs, syms))
