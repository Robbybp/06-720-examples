import os
import pyomo.environ as pyo
from pyutilib.common._exceptions import ApplicationError

m = pyo.ConcreteModel(name='Redundant equations')

FILE_DIR = os.path.dirname(__file__)

executable = os.path.join(FILE_DIR, os.pardir, 'newton')
solver = pyo.SolverFactory(executable)

"""
A simple system with no true solution, but a solution
within any arbitrarily small positive tolerance.
"""

def c1_rule(m):
    return m.x2 - pyo.cos(m.x1) + 1 == 0.

def c2_rule(m):
    return m.x2 - pyo.exp(-m.x1) == 0.

m.x1 = pyo.Var()
m.x2 = pyo.Var()

m.c1 = pyo.Constraint(rule=c1_rule)
m.c2 = pyo.Constraint(rule=c2_rule)

m.x1 = -10.
m.x2 = -10.

try:
    solver.solve(m, tee=True)
except ApplicationError:
    pass

print()
print('Solution as recorded by Pyomo:')
m.x1.display()
m.x2.display()

"""
Interesting that we can converge this system from such a poor starting point
"""
