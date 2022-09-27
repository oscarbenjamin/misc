operations = {
    Add: np.add,
    Mul: np.multiply,
    Pow: np.power,
    exp: np.exp,
}

def nsubs(expr, replacement):
    rep = replacement.get(expr, None)
    if rep is not None:
        return rep
    args = expr.args
    if args:
        args_n = [nsubs(a, replacement) for a in args]
        return operations[type(expr)](*args_n)
    elif expr.is_Number:
        return np.array(float(expr))
    elif expr.is_ImaginaryUnit:
        return np.array(complex(expr))
    else:
        raise NotImplementedError

from scipy.integrate import quad

def quad_lam(expr, sym, start, end):
    f = lambdify(sym, expr)
    return quad(f, start, end)

def quad_sub(expr, sym, start, end):
    f = lambda v: nsubs(expr, {sym:v})
    return quad(f, start, end)
