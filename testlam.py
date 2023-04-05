from random import normalvariate
from sympy import *
from flint import arb

# Everything from sympy.functions.__all__

functions = [
    factorial, factorial2, rf, ff, binomial, RisingFactorial, FallingFactorial,
    subfactorial, carmichael, fibonacci, lucas,
    # XXX: motzkin not a symbolic function
    #motzkin,
    tribonacci, harmonic,
    bernoulli, bell, euler, catalan, genocchi, andre, partition, sqrt, root,
    Min, Max, Id, real_root, cbrt, Rem, re, im, sign, Abs, conjugate, arg,
    # XXX: Not sure any of these make sense for mpmath/arb:
    #polar_lift, periodic_argument, unbranched_argument, principal_branch,
    #transpose, adjoint, polarify, unpolarify,
    sin, cos, tan, sec, csc, cot,
    sinc, asin, acos, atan, asec, acsc, acot, atan2,
    # exp_polar,
    exp, ln, log,
    LambertW, sinh, cosh, tanh, coth, sech, csch, asinh, acosh, atanh, acoth,
    asech, acsch, floor, ceiling, frac,
    # Probably Piecewise should be tested, needs special consideration for Arb
    # Piecewise,
    # Non-symbolic functions
    # piecewise_fold, piecewise_exclusive,
    erf, erfc, erfi, erf2, erfinv, erfcinv, erf2inv, Ei,
    expint, E1, li, Li, Si, Ci, Shi, Chi, fresnels, fresnelc, gamma,
    lowergamma, uppergamma, polygamma, loggamma, digamma, trigamma, multigamma,
    dirichlet_eta, zeta, lerchphi, polylog, stieltjes, riemann_xi,
    # These possibly don't make sense for mpmath/arb
    #Eijk, LeviCivita, KroneckerDelta, SingularityFunction, DiracDelta,
    Heaviside,
    # These possibly don't make sense for mpmath/arb
    # bspline_basis, bspline_basis_set, interpolating_spline,
    besselj, bessely,
    besseli, besselk,
    hankel1, hankel2, jn, yn,
    # Non-symbolic
    # jn_zeros,
    hn1, hn2, airyai,
    airybi, airyaiprime, airybiprime, marcumq,
    # These have different signatures so need bespoke testing:
    # hyper, meijerg, appellf1,
    legendre, assoc_legendre, hermite, hermite_prob, chebyshevt, chebyshevu,
    # Non-symbolic
    #chebyshevu_root, chebyshevt_root,
    laguerre, assoc_laguerre, gegenbauer,
    jacobi,
    # Non-symbolic
    # jacobi_normalized,
    Ynm,
    # Non-symbolic
    # Ynm_c,
    Znm, elliptic_k, elliptic_f,
    elliptic_e, elliptic_pi, beta, betainc, betainc_regularized, mathieus,
    mathieuc, mathieusprime, mathieucprime
]

values = [0, 1, 2]
values += [normalvariate(0, 1) for _ in range(10)]
#values = [v for v in values if v >= 0]

x, y, z = symbols('x, y, z')

def pairs1(v):
    return [(v, vp) for vp in values]

def pairs2(v):
    return [(vp, v) for vp in values]

def allpairs():
    return [(v1,v2) for v1 in values for v2 in values]

def alltriples():
    return [(v1,v2,v3) for v1 in values for v2 in values for v3 in values]

failures = {
    # DomainError
    sqrt: [v for v in values if v < 0],
    cbrt: [v for v in values if v < 0],
    root: [(v1,v2) for v1,v2 in allpairs() if v1 <= 0 or v2 == 0 or v2 != int(v2)] + [
        (v1,v2,v3) for v1,v2,v3 in alltriples() if v1 <= 0 or v2 == 0 or v2 != int(v2) or v3 != int(v3)],
    # nan etc
    beta: [0] + pairs1(0) + pairs2(0),
    # TypeError: unsupported operand type(s) for %
    factorial2: values,
    frac: values,
    # ValueError: gamma function pole (probably due to rewrite)
    rf: pairs1(0) + [(x,y) for (x,y) in allpairs() if x == y or x < 0 or y < 0 or y - x == 1],
    ff: pairs1(0) + [(x,y) for (x,y) in allpairs() if x == y or x < 0 or y < 0 or y - x == 1],
    binomial: [(0,1), (0,2), (1,2)],
    # e.g. NameError: name 'uppergamma' is not defined
    # TypeError: f() takes 1 positional argument but 2 were given
    subfactorial: values,
    carmichael: values + allpairs(),
    fibonacci: values + allpairs(),
    lucas: values,
    tribonacci: values + allpairs(),
    harmonic: values + allpairs(),
    bernoulli: values + allpairs(),
    bell: values + allpairs() + alltriples(),
    euler: values + allpairs(),
    genocchi: values + allpairs(),
    andre: values + allpairs(),
    partition: values,
    real_root: allpairs(),
    Rem: allpairs(),
    re: values,
    im: values,
    conjugate: values,
    arg: values,
    LambertW: values + allpairs(),
    erfi: values,
    erfinv: values,
    erfcinv: values,
    erf2inv: values + allpairs(),
    Ei: values,
    expint: allpairs(),
    E1: values,
    li: values,
    Li: values,
    Si: values,
    Ci: values,
    Shi: values,
    Chi: values,
    lowergamma: values + allpairs(),
    uppergamma: values + allpairs(),
    polygamma: [(0,0)] + allpairs(),
    dirichlet_eta: values + allpairs(),
    zeta: values + allpairs(),
    polylog: allpairs(),
    stieltjes: values + allpairs(),
    riemann_xi: values,
    hankel1: allpairs(),
    hankel2: allpairs(),
    jn: allpairs(),
    yn: allpairs(),
    hn1: allpairs(),
    hn2: allpairs(),
    airyai: values,
    airybi: values,
    airyaiprime: values,
    airybiprime: values,
    carmichael: values + allpairs() + alltriples(),
    marcumq: alltriples(),
    legendre: allpairs(),
    assoc_legendre: alltriples(),
    hermite: allpairs(),
    hermite_prob: allpairs(),
    chebyshevt: allpairs(),
    chebyshevu: allpairs(),
    assoc_laguerre: alltriples(),
    gegenbauer: alltriples(),
    elliptic_k: values,
    elliptic_f: allpairs(),
    elliptic_e: values + allpairs(),
    elliptic_pi: allpairs() + alltriples(),
    mathieus: alltriples(),
    mathieuc: alltriples(),
    mathieusprime: alltriples(),
    mathieucprime: alltriples(),
    # Zero should be skipped but these also fail with NameError:
    digamma: [0] + values,
    trigamma: [0] + values,
    # Fails in mpmath as well:
    multigamma: allpairs(),
    lerchphi: values + allpairs() + alltriples(),
    # ZeroDivisionError (probably due to rewrites)
    csc: [0],
    cot: [0],
    asec: [v for v in values if not (v <= -1) or (v >= 1)],
    acsc: [v for v in values if not (v <= -1) or (v >= 1)],
    acot: [0],
    log: [v for v in values if v <= 0] + [(v1,v2) for v1, v2 in allpairs() if v1 <= 0 or v2 <= 0 or v2 == 1],
    coth: [0],
    csch: [0],
    acoth: [v for v in values if not (v < -1) or (v > 1)],
    asech: [0] + [v for v in values if v > 1 or v < 0],
    acsch: [0],
    acosh: [v for v in values if v < 1],
    acos: [v for v in values if not (-1 < v < 1)],
    asin: [v for v in values if not (-1 < v < 1)],
    atanh: [v for v in values if not (-1 < v < 1)],
    # gamma pole:
    gamma: [0],
    # TypeError: laguerre() missing 1 required positional argument: 'z'
    # fails in mpmath...
    laguerre: allpairs(),
}

def nargs(f):
    if f in [sqrt, cbrt, E1]:
        return {1}
    elif f is root:
        return {2,3}
    elif f is real_root:
        return {1,2}
    elif f is loggamma:
        # KeyError: 'loggamma' (during lambdify)
        return set()
    else:
        return f.nargs


skipped = 0
passed = 0
failed = 0


def update_count(fstr, val1, val2, eps=1e-11):
    global skipped, passed, failed
    diff = complex(val1) - complex(val2)
    if abs(diff) < eps:
        passed += 1
    else:
        print(fstr, '->', diff, val1, val2)
        failed += 1


for func in functions:
    fail_values = failures.get(func, ())

    if 1 in nargs(func):
        f_mpmath = lambdify(x, func(x), modules='mpmath')
        f_arb = lambdify(x, func(x), modules='arb')
        for value in values:
            fstr = f'{func}({value})'
            #print(fstr)
            if value not in fail_values:
                mpval = f_mpmath(value)
                arbval = f_arb(arb(value))
                update_count(fstr, mpval, arbval)
            else:
                #print('skipping', fstr)
                skipped += 1

    if 2 in nargs(func):
        f_mpmath = lambdify((x, y), func(x, y), modules='mpmath')
        f_arb = lambdify((x, y), func(x, y), modules='arb')
        for value1 in values:
            for value2 in values:
                fstr = f'{func}({value1}, {value2})'
                #print(fstr)
                if (value1, value2) not in fail_values:
                    mpval = f_mpmath(value1, value2)
                    arbval = f_arb(arb(value1), arb(value2))
                    update_count(fstr, mpval, arbval)
                else:
                    #print('skipping', fstr)
                    skipped += 1

    if 3 in nargs(func):
        f_mpmath = lambdify((x, y, z), func(x, y, z), modules='mpmath')
        f_arb = lambdify((x, y, z), func(x, y, z), modules='arb')
        for value1 in values:
            for value2 in values:
                for value3 in values:
                    fstr = f'{func}({value1}, {value2}, {value3})'
                    #print(fstr)
                    if (value1, value2, value3) not in fail_values:
                        mpval = f_mpmath(value1, value2, value3)
                        arbval = f_arb(arb(value1), arb(value2), arb(value3))
                        update_count(fstr, mpval, arbval)
                    else:
                        #print('skipping', fstr)
                        skipped += 1

    if nargs(func) - {1, 2, 3} and func not in [
            carmichael, bell, Min, Max, lerchphi, # <-- any number of args
            jacobi, Ynm, Znm, betainc, betainc_regularized # <-- these take 4 args
            ]:
        assert False


print()
print('total:', skipped + passed + failed)
print('skipped:', skipped)
print('passed:', passed)
print('failed:', failed)
