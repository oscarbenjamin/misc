from sympy import Symbol, Matrix, I

gamma_c = Symbol("gamma_c")
n_c = Symbol("n_c")
gamma_h = Symbol("gamma_h")
n_h = Symbol("n_h")
epsilon_c = Symbol("epsilon_c")
g = Symbol("g")
epsilon_h = Symbol("epsilon_h")

M = Matrix([[-gamma_c*(n_c + 1) - gamma_h*(n_h + 1), gamma_c*n_c, 0, gamma_h*n_h, 0], [gamma_c*(n_c + 1) - gamma_h*n_h, -gamma_c*n_c - gamma_h*n_h - gamma_h*(n_h + 1), I*g, -gamma_h*n_h, -I*g], [0, I*g, -gamma_c*n_c/2 - gamma_c*(n_c + 1)/2 - gamma_h*n_h/2 - gamma_h*(n_h + 1)/2 - I*(-epsilon_c + epsilon_h), -I*g, 0], [-gamma_c*n_c + gamma_h*(n_h + 1), -gamma_c*n_c, -I*g, -gamma_c*n_c - gamma_c*(n_c + 1) - gamma_h*n_h, I*g], [0, -I*g, 0, I*g, -gamma_c*n_c/2 - gamma_c*(n_c + 1)/2 - gamma_h*n_h/2 - gamma_h*(n_h + 1)/2 - I*(epsilon_c - epsilon_h)]])
