import os
import pyomo.environ as pyo

m = pyo.ConcreteModel(name='06-720 A1.3')

FILE_DIR = os.path.dirname(__file__)

executable = os.path.join(FILE_DIR, os.pardir, 'newton')
solver = pyo.SolverFactory(executable)

"""
This is a system that converges slowly because it starts close to
a point of singular Jacobian.
"""

def c1_rule(m):
    return m.x2 - (m.x1-2)**3/8 + 3*m.x1 - 15 == 0

def c2_rule(m):
    return m.x2 - m.x1 + 2 == 0

m.x1 = pyo.Var()
m.x2 = pyo.Var()

m.c1 = pyo.Constraint(rule=c1_rule)
m.c2 = pyo.Constraint(rule=c2_rule)

m.x1 = 5.
m.x2 = 5.

solver.solve(m, tee=True)

print()
print('Solution as recorded by Pyomo:')
m.x1.display()
m.x2.display()
