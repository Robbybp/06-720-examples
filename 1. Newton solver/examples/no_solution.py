import os
import pyomo.environ as pyo

m = pyo.ConcreteModel(name='06-720 A1.3')

FILE_DIR = os.path.dirname(__file__)

executable = os.path.join(FILE_DIR, os.pardir, 'newton')
solver = pyo.SolverFactory(executable)

"""
A simple system that has no solution.
(A convex quadratic with a line below it.)
"""

def c1_rule(m):
    return m.x1**2/8 - m.x2 + 2 == 0

def c2_rule(m):
    return 4*m.x1 - 5*m.x2 - 5 == 0

m.x1 = pyo.Var()
m.x2 = pyo.Var()

m.c1 = pyo.Constraint(rule=c1_rule)
m.c2 = pyo.Constraint(rule=c2_rule)

m.x1 = 0.
m.x2 = 0.

solver.solve(m, tee=True)

print()
print('Solution as recorded by Pyomo:')
m.x1.display()
m.x2.display()

"""
If this solver was doing something "more advanced," like minimizing
infeasibility, we would expect it to converge somewhere near
(x1, x2) = (3.5, 3.0).
Instead it fails to converge and ends up at a seemingly random point.
This is a strong motivator for reframing "solving a system of nonlinear
equations" as "minimizing infeasibility of a system of equations."
"""
