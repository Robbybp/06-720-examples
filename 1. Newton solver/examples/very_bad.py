import os
import pyomo.environ as pyo

import math

m = pyo.ConcreteModel(name='Pathological case')

FILE_DIR = os.path.dirname(__file__)

executable = os.path.join(FILE_DIR, os.pardir, 'newton')
solver = pyo.SolverFactory(executable)

"""
This is about as bad of a situation as I can imagine.
We have points of very nearly singular Jacobian at 2\pi*i.
Starting near one of these points pushes us to an
essentially random point in space, from which it is essentially
impossible to converge due to the nearly indistinuishable
nature of our surface.
"""

def c1_rule(m):
    return m.x2 + pyo.cos(m.x1) - 2 == 0

def c2_rule(m):
    return m.x2 - pyo.exp(-10*m.x1**2)/2 - 1/2 == 0
#    return m.x2 + m.x1 == 0
# A more reasonable function like this line converges quickly,
# even from a bad starting point.

m.x1 = pyo.Var()
m.x2 = pyo.Var()

m.c1 = pyo.Constraint(rule=c1_rule)
m.c2 = pyo.Constraint(rule=c2_rule)

m.x1 = 2*math.pi - 1e-8
m.x2 = 1.

solver.solve(m, tee=True)

print()
print('Solution as recorded by Pyomo:')
m.x1.display()
m.x2.display()
