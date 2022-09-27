import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def latexd(expr):
    return rf'${latex(expr)}$'

def show_zeros(expr, x, y, subs={}, xlim=(-2, 2), ylim=(-2, 2), npoints=1000):
    if isinstance(expr, list):
        multiple = True
        exprs = expr
    else:
        multiple = False
        exprs = [expr]

    xvals = np.linspace(*xlim, 1000)
    yvals = np.linspace(*ylim, 1000)
    X, Y = np.meshgrid(xvals, yvals)
    fs1 = [lambdify((x, y), expr.subs(subs)) for expr in exprs]
    fs2 = [lambdify((x, y), sqf_part(expr.subs(subs))) for expr in exprs]
    Fs1 = [f(X, Y) for f in fs1]
    Fs2 = [f(X, Y) for f in fs2]
    Fcs = [np.sign(F)*np.log(1 + np.abs(F)) for F in Fs1]

    plt.figure()
    if not multiple:
        [Fc] = Fcs
        maxabs = np.max(np.abs(Fc))
        plt.imshow(Fc[::-1,:], extent=xlim+ylim, vmin=-maxabs, vmax=maxabs, aspect='auto')
        for F in Fs2:
            obj = plt.contour(X, Y, F, levels=[0])
    else:
        for F, col in zip(Fs2, ['red', 'blue']):
            obj = plt.contour(X, Y, F, levels=[0], colors=[col])
    title = ',    '.join([latexd(e) for e in exprs])
    if subs:
        title += '\t' + latexd(subs)
        title += '\n' + ',    '.join([latexd(e.subs(subs)) for e in exprs])
    plt.title(title)
    plt.xlabel(latexd(x), size='x-large')
    plt.ylabel(latexd(y), size='x-large', rotation=0)
