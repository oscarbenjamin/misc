import numpy as np
from sympy import *

x, y = symbols('x, y')

funcs = [besseli, besselj, besselk, bessely, jn, yn]

points = [
    0, 1, -1, I, -I, 1+I, 1-I, -1+I, -1-I,
    #oo, -oo, I*oo, -I*oo, zoo,
    pi, pi/2,
]
points = list(np.random.randn(10)) + list(np.random.randn(10)+1j*np.random.randn(10))
points = list(range(-5, 5))

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
            continue
        fine = True
        for p1 in points:
            for p2 in points:
                try:
                    rep = {x:p1, y:p2}
                    if expr1.subs(rep).n() != expr2.subs(rep).n():
                        print()
                        print(f1, f2, rep)
                        print(expr1, ' --> ', expr2)
                        print(expr1.subs(rep), expr2.subs(rep))
                        print(expr1.subs(rep).n(5), expr2.subs(rep).n(5))
                        fine = False
                        break
                except (Exception, RecursionError) as e:
                    print(e)
                    print(f1, f2, rep)
            if not fine:
                break
        if not fine:
            bad += 1

print(bad, 'bad rewrites out of', total)
