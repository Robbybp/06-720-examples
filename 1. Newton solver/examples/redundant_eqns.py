import os
import pyomo.environ as pyo
from pyutilib.common._exceptions import ApplicationError

m = pyo.ConcreteModel(name='Redundant equations')

FILE_DIR = os.path.dirname(__file__)

executable = os.path.join(FILE_DIR, os.pardir, 'newton')
solver = pyo.SolverFactory(executable)

"""
A simple system with redundant equations
"""

def c1_rule(m):
    return m.x1**2 + m.x2**2 == 1.

def c2_rule(m):
    return m.x3 == 3.

def c3_rule(m):
    return m.x1**2 + m.x2**2 + m.x3**2 == 10.

m.x1 = pyo.Var()
m.x2 = pyo.Var()
m.x3 = pyo.Var()

m.c1 = pyo.Constraint(rule=c1_rule)
m.c2 = pyo.Constraint(rule=c2_rule)
m.c3 = pyo.Constraint(rule=c3_rule)

m.x1 = 10.
m.x2 = 10.
m.x3 = 10.

try:
    solver.solve(m, tee=True)
except ApplicationError:
    pass

print()
print('Solution as recorded by Pyomo:')
m.x1.display()
m.x2.display()
m.x3.display()

"""
We fail to factorize in the first iteration because our Jacobian is singular.
I believe it is an open question what a solver should do if given equations
like this.
"""
