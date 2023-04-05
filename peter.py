import sympy as sm
import random
import time


random.seed(105)

N = 3     # dimension of the problem

q_list = [sm.symbols('q' + str(i)) for i in range(N)]     # variables to be solved for
param_list = [sm.symbols('a' + str(i)) for i in range(5)] # list of parameters, '5' is arbitary
pL_vals    = [i for i in range(5)]

# set up arbitrary linear equations
gleichung = []
for _ in range(N):                 # number of equations
    summe1 = 0
    for q in q_list:
        summe = random.choice(param_list)
        for _ in range(len(param_list)):
            wert = random.choice(param_list)
            c1 = random.choice(param_list)
            c2 = random.choice([sm.sin(c1), sm.cos(c1), sm.exp(c1)])
            summe += wert * c2
        summe1 += summe * q + random.choice(param_list)
    gleichung.append(summe1)

# convert the list to an sm.Matrix, so operation may be done
equation = sm.Matrix(gleichung) # equation = 0 to be solved for q_list
# lambdification to check the accuracy of the solutions later
equation_lam = sm.lambdify(q_list + param_list, equation, cse=True) 

# LU decomposition
start = time.time()
matrix_A = equation.jacobian(q_list)
vector_b = equation.subs({i: 0 for i in q_list})
solution = matrix_A.LUsolve(-vector_b)
print('it took {:.5f} sec to solve with LU.solve'.format(time.time() - start))

op_zahl = sum([solution[l].count_ops(visual=False) for l in range(N)])
print('solution contains {} operations'.format(op_zahl))

solution_lam = sm.lambdify(param_list, solution, cse=True )
loesung = [solution_lam(*pL_vals)[i][0] for i in range(N)]
print('the solution with ULsolve is:', loesung)
print('The closer these values are to zero, the better the solution', equation_lam(*loesung, *pL_vals),'\n')


# using sm.linsolve
if N < 3:
    start = time.time() 
    solution, = sm.linsolve(equation, q_list)
    print('it took {:.5f} sec to solve with sm.linsolve'.format(time.time() - start))
    
    op_zahl = sum([solution[l].count_ops(visual=False) for l in range(N)])
    print('solution contains {} operations'.format(op_zahl))
       
    solution_lam = sm.lambdify(param_list, solution, cse=True)
    loesung = [solution_lam(*pL_vals)[i] for i in range(N)]
    print('the solution with linsolve is:', loesung, '\n')

    
# using sm.solve
if N < 3:
    start = time.time() 
    solution = sm.solve(equation, q_list)
    print('it took {:.5f} sec to solve with sm.solve'.format(time.time() - start))
    
    op_zahl = sum([solution[l].count_ops(visual=False) for l in q_list])
    print('solution contains {} operations'.format(op_zahl))
    
    solution = [solution[i] for i in q_list]
    solution_lam = sm.lambdify(param_list, solution, cse=True)
    loesung = solution_lam(*pL_vals)
    print('the solution with sm.sove is:', loesung)
