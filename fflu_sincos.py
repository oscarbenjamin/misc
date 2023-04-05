from x import t

A, b = linear_eq_to_matrix(equation, q_list)


def matrix_solve(A, b):

    Ab = Matrix.hstack(A, b)

    replacements = {}
    reverse = {}
    expi = {}

    num = 0

    def new_sym(name):
        nonlocal num
        num += 1
        return Dummy(f'{name}{num}')

    polys, opts = parallel_poly_from_expr(Ab)

    gens = opts['gens']
    transcendentals = []
    sincos = []
    newgens = []

    for g in gens:
        if g.func == exp:
            transcendentals.append(g)
        elif g.func in (sin, cos):
            arg = g.args[0]
            if sin(arg) in gens and cos(arg) in gens:
                sincos.append((sin(arg), cos(arg)))
            else:
                transcendentals.append(g)
        else:
            newgens.append(g)

    for g in transcendentals:
        sym = new_sym('t')
        replacements[g] = sym
        newgens.append(sym)

    denom = 1

    for s, c in sincos:
        sym = new_sym('z')
        replacements[s] = (sym - 1/sym)/I
        replacements[c] = (sym + 1/sym)
        reverse[sym] = (c + I*s)
        denom *= 2*sym
        newgens.append(sym)

    # This is multiplying both sides of the equation
    f = [expand_mul(denom * p.as_expr().subs(replacements)) for p in polys]

    Ab = Matrix(f).reshape(*Ab.shape)

    A = Ab[:, :A.shape[1]]
    b = Ab[:, A.shape[1]:]

    sol_num, sol_den = fflu_solve(A, b)

    sol_num = factor_terms(sol_num)
    sol_den = factor_terms(sol_den)

    powers = sol_den.as_powers_dict()
    for num in sol_num:
        new_powers = num.as_powers_dict()
        powers = {expr: min(powers[expr], new_powers.get(expr, 0)) for expr in powers}

    gcd_terms = Mul(*(expr ** p for expr, p in powers.items()))

    sol_num /= gcd_terms
    sol_den /= gcd_terms

    return sol_num, sol_den, reverse

    for s in Ab.atoms(exp):
        sym = new_sym()
        replacements[s] = sym
        reverse[sym] = s

    sincos_atoms = Ab.atoms(sin, cos)

    denom = 1

    for s in sincos_atoms:
        if s in replacements:
            continue
        arg = s.args[0]
        sym = new_sym()
        if s.func == cos:
            other = sin(arg)
        else:
            other = cos(arg)
        if other not in sincos_atoms:
            replacements[s] = sym
            reverse[sym] = s
        else:
            replacements[sin(arg)] = (sym - 1/sym)/(2*I)
            replacements[cos(arg)] = (sym + 1/sym)/2
            reverse[sym] = cos(arg) + I*sin(arg)

    A = A.subs(replacements)
    b = b.subs(replacements)

    breakpoint()

    A = DomainMatrix.from_Matrix(A.subs(replacements))
    b = DomainMatrix.from_Matrix(b.subs(replacements))
    A, b = A.unify(b)
    print(A.domain, '<--- hopefully not EX')

    sol, den = charpoly_solve(A, b, denom)

    return sol.subs(reverse) / den.subs(reverse)


def charpoly_solve(A, b, denom):
    pA = A.charpoly()

    pA = [A.domain.to_sympy(pAi) for pAi in pA]
    A = A.to_Matrix()
    b = b.to_Matrix()

    # pA(A) = pA[0]*A**n + PA[1]*A**(n-1) + ... + pA[n-1]*I = 0
    # Therefore
    # A**-1 = (pA[0]*A**(n-1) + pA[1]*A**(n-2) + ... + pA[n-2]) / (-pA[n-1])
    # Hence
    # x = A**-1*b
    # x = (pA[0]*A + pA[1]*(A + pA[2]*(A + ... pA[n-2]*I)*b)) / (-pA[n-1])
    result = b
    for i in range(1, A.shape[0]):
        result = pA[i]*b + A*result

    f = - denom * pA[-1]

    return result, f


def f():

    sol_num, sol_den = fflu_solve(A.subs(replacements), b.subs(replacements))

    _expand = lambda e: expand_mul(expand_multinomial(e))

    denominator = (denom * sol_den).subs(reverse)
    numerator = sol_num.subs(reverse)

    return numerator, denominator

def fflu_solve(M, b):
    m, n = M.shape
    m, o = b.shape
    LU = DomainMatrix.from_Matrix(M)
    b = DomainMatrix.from_Matrix(b)
    LU, b = LU.unify(b)
    LU = LU.to_dense().rep
    b = b.to_dense().rep
    K = LU.domain
    print(K)
    swaps = ddm_ifflu(LU, K)
    y = ddm_fffs(LU, b, swaps, K)
    x, d = ddm_ffbs(LU, y, K)
    sol_num = DomainMatrix(x, (n, o), K).to_Matrix()
    sol_den = K.to_sympy(d)
    return sol_num, sol_den

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
        if not a[i][i]._random():
            for ip in range(i+1, m):
                if a[ip][i]._random():
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
                    a[j][k] = K.exquo(a[j][k], a[i-1][i-1])

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
                    y[j][k] = K.exquo(y[j][k], LU[i-1][i-1])

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
            x[i][k] = K.exquo(x[i][k], LU[i][i])
    return x, d
