import os
import pyomo.environ as pyo

import math

FILE_DIR = os.path.dirname(__file__)

executable = os.path.join(FILE_DIR, os.pardir, 'newton')
solver = pyo.SolverFactory(executable)

m = pyo.ConcreteModel(name='shrinking-core')

m.comp = pyo.Set(initialize=['A', 'B'])

m.x = pyo.Var(m.comp, initialize=0.5)

m.comp_dens_pure = pyo.Param(m.comp, initialize={'A': 8., 'B': 6.})

m.dens = pyo.Var(initialize=3.5)

m.d_particle = pyo.Param(mutable=True, initialize=1.)
m.V_particle = pyo.Expression(expr=1/6*math.pi*m.d_particle**3)
m.interface_radius = pyo.Var(initialize=0.5)

def rad_rule(m):
    return (
            m.interface_radius**3 ==
            (m.d_particle/2)**3 *
            m.x['A']*m.comp_dens['A']/m.comp_dens_pure['A']
            )

m.comp_holdup = pyo.Var(m.comp, initialize={'A': 4., 'B': 3.})
m.comp_holdup.fix()
def holdup_rule(m, j):
    return m.comp_holdup[j] == m.V_particle * m.dens * m.x[j]
m.holdup_eqn = pyo.Constraint(m.comp, rule=holdup_rule)

def sum_rule(m):
    return 1. == sum(m.x[j] for j in m.comp)
m.sum_eqn = pyo.Constraint(rule=sum_rule)

solver.solve(m, tee=True)

import pdb; pdb.set_trace()
