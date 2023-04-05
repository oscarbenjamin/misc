from sympy import EXRAW, symbols, Matrix, nan, simplify

def is_nonzero(e, K):
    if K == EXRAW:
        # Evaluate numerically with random values
        e = e._random()
    return bool(e)

def exquo(num, den, K):
    """Exact quotient for a division that should cancel."""
    if K == EXRAW:
        numquo, numrem = _exquo(num, den)
        if numrem is not S.Zero and is_nonzero(numrem, K):
            breakpoint()
            numquo = num / den
        return numquo
    else:
        return K.exquo(num, den)

def _exquo(num, den):
    if num == den:
        result = S.One, S.Zero
    elif num.is_Pow:
        result = _exquo_pow(num, den)
    elif num.is_Add:
        result = _exquo_add(num, den)
    elif num.is_Mul:
        result = _exquo_mul(num, den)
    else:
        result = S.Zero, num
    quo, rem = result
    # assert (quo*den + rem - num).cancel() == 0
    return result

def _exquo_add(num, den):
    termsquo = []
    termsrem = []
    for arg in num.args:
        argquo, argrem = _exquo(arg, den)
        termsquo.append(argquo)
        termsrem.append(argrem)
    return Add(*termsquo), Add(*termsrem)

def _exquo_mul(num, den):
    for n, arg in enumerate(num.args):
        argquo, argrem = _exquo(arg, den)
        if argrem is S.Zero:
            numquo = Mul(*num.args[:n], argquo, *num.args[n+1:])
            return numquo, S.Zero
        if argquo is not S.Zero:
            prev = Mul(*num.args[:n])
            rest = Mul(*num.args[n+1:])
            restquo, restrem = _exquo(rest, den)
            numquo = argquo*restquo*den + argquo*restrem + argrem*restquo
            numrem = argrem * restrem
            return prev*numquo, prev*numrem
    return S.Zero, num

def _exquo_pow(num, den):
    base, exp = num.args
    #assert exp.is_Integer and exp > 1
    basequo, baserem = _exquo(base, den)
    if baserem is S.Zero:
        return base**(exp-1) * basequo, S.Zero
    else:
        return S.Zero, num

def fflu_solve_mat(LU, b):
    LU = LU.tolist()
    b = b.tolist()
    swaps = ddm_ifflu(LU, EXRAW)
    y = ddm_fffs(LU, b, swaps, EXRAW)
    x, d = ddm_ffbs(LU, y, EXRAW)
    return Matrix(x), d

def ddm_ifflu(a, K):
    """a  <--  LU(a)"""
    if not a:
        return a
    m, n = len(a), len(a[0])

    swaps = []

    for i in range(min(m, n)):
        if not is_nonzero(a[i][i], K):
            for ip in range(i+1, m):
                if is_nonzero(a[ip][i], K):
                    swaps.append((i, ip))
                    a[i], a[ip] = a[ip], a[i]
                    break
            else:
                # M = Matrix([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 1, 2]])
                continue
        for j in range(i+1, m):
            for k in range(i+1, n):
                a[j][k] = a[i][i]*a[j][k] - a[j][i]*a[i][k]
                if i > 1:
                    a[j][k] = exquo(a[j][k], a[i-1][i-1], K)

    return swaps

def ddm_fffs(LU, b, swaps, K):
    """Solve L y = b"""
    m, n, o = len(LU), len(LU[0]), len(b[0])

    y = [row[:]  for row in b]
    for i1, i2 in swaps:
        y[i1], y[i2] = y[i2], y[i1]

    for i in range(m-1):
        for j in range(i+1, m):
            for k in range(o):
                y[j][k] = LU[i][i]*y[j][k] - LU[j][i]*y[i][k]
                if i > 1:
                    y[j][k] = exquo(y[j][k], LU[i-1][i-1], K)

    if m > n:
        for i in range(n, m):
            for j in range(o):
                if y[i][j]:
                    raise NonInvertibleMatrixError

    return y

def ddm_ffbs(LU, y, K):
    """Solve Ux = y"""
    m, n, o = len(LU), len(LU[0]), len(y[0])

    x = [row[:]  for row in y]

    d = LU[n-1][n-1]  # determinant
    for k in range(o):
        for i in reversed(range(n)):
            if not LU[i][i]:
                raise NonInvertibleMatrixError
            x[i][k] = d * x[i][k]
            for j in range(i+1, n):
                x[i][k] = x[i][k] - LU[i][j] * x[j][k]
            x[i][k] = exquo(x[i][k], LU[i][i], K)
    return x, d


x = symbols('x')
M = Matrix([[(x + 1)**2 - (x**2 + 2*x + 1), x], [x, 0]])
b = Matrix([1, 1])

sol_num, sol_den = fflu_solve_mat(M, b)

sol = sol_num / sol_den

assert sol.subs(x, 1).has(nan) is False
assert sol_num.subs(x, 1).has(nan) is False
assert simplify(M*sol - b).is_zero_matrix is True

