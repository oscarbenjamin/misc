import numpy as np
from sympy import *

x, y = symbols('x, y')

funcs = [jn, yn, besseli, besselj, besselk, bessely]

special_points = [
    0, S.Half, -S.Half, S(3)/2, -S(3)/2, 1, 2, 3, -1, -2, -3,
    pi, pi/2,
    I, -I, 1+I, 1-I, -1+I, -1-I,
    #oo, -oo, I*oo, -I*oo, zoo,
]
generic_points = (list(np.random.randn(5))
                + list(1j*np.random.randn(5))
                + list(np.random.randn(5)+1j*np.random.randn(5)))

points = special_points
points = generic_points

bad = 0
total = 0
for f1 in funcs:
    for f2 in funcs:
        total += 1
        expr1 = f1(x, y)
        try:
            expr2 = expr1.rewrite(f2)
        except Exception as e:
            print('rewrite failed:', f1, '->', f2)
            1/0
        fine = True
        for p1 in points:
            for p2 in points:
                try:
                    rep = {x:p1, y:p2}
                    if (expr1.subs(rep).n().n(4, chop=True)
                            != expr2.subs(rep).n().n(4, chop=True)
                            and not expr2.subs(rep).n().has(jn, yn)):
                        print()
                        print(f1, f2, rep)
                        print(expr1, ' --> ', expr2)
                        print(expr1.subs(rep), expr2.subs(rep))
                        print(expr1.subs(rep).n(5), expr2.subs(rep).n(5))
                        fine = False
                        break
                    expr1.subs(rep).rewrite(f2)
                except (Exception, RecursionError) as e:
                    print(e)
                    print(f1, f2, rep)
        if not fine:
            bad += 1

print(bad, 'bad rewrites out of', total)
