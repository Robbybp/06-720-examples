import os
import pyomo.environ as pyo

from pyutilib.common._exceptions import ApplicationError

import math

FILE_DIR = os.path.dirname(__file__)

executable = os.path.join(FILE_DIR, os.pardir, 'newton')
solver = pyo.SolverFactory(executable)

"""
One time-step of a dynamic shrinking core reaction model.
"""

m = pyo.ConcreteModel(name='shrinking-core')

m.comp = pyo.Set(initialize=['A', 'B'])

m.x = pyo.Var(m.comp, initialize=0.5)

m.comp_dens_pure = pyo.Param(m.comp, mutable=True, initialize={'A': 8., 'B': 6.})

m.dens = pyo.Var(initialize=3.5)

m.d_particle = pyo.Param(mutable=True, initialize=1.)
m.V_particle = pyo.Expression(expr=1/6*math.pi*m.d_particle**3)
m.interface_radius = pyo.Var(initialize=0.5)

def rad_rule(m):
    return (
            m.interface_radius**3 ==
            (m.d_particle/2)**3 *
            (m.x['A']*m.dens/m.comp_dens_pure['A'])
            #m.interface_radius ==
            #(m.d_particle/2) *
            #(m.x['A']*m.dens/m.comp_dens_pure['A'])**(1/3)
            # "Cube-root formulation" leads to ASL evaluation errors
            # as we evaluate 1/3 as a decimal float.
            )
m.radius_eqn = pyo.Constraint(rule=rad_rule)

m.comp_holdup = pyo.Var(m.comp, initialize={'A': 4., 'B': 3.})
m.comp_holdup.fix()
def holdup_rule(m, j):
    return m.comp_holdup[j] == m.V_particle * m.dens * m.x[j]
m.holdup_eqn = pyo.Constraint(m.comp, rule=holdup_rule)

def sum_rule(m):
    #return 1. == sum(m.x[j] for j in m.comp)
    return m.V_particle * m.dens == sum(m.comp_holdup[j] for j in m.comp)
m.sum_eqn = pyo.Constraint(rule=sum_rule)

m.k0 = pyo.Var(initialize=1.)
m.k0.fix()
m.temperature = pyo.Var(initialize=1275.)
m.temperature.fix()
m.act_energy = pyo.Var(initialize=5.)
m.act_energy.fix()
m.gas_const = pyo.Param(initialize=0.00831)
m.rate_coef = pyo.Var(initialize=100.)
def rate_coef_rule(m):
    return (
            m.rate_coef ==
            m.k0 * pyo.exp(-m.act_energy/(m.gas_const*m.temperature))
            )
m.rate_coef_eqn = pyo.Constraint(rule=rate_coef_rule)

m.reaction_rate = pyo.Var(initialize=1.)
def reaction_rate_rule(m):
    return (
            m.reaction_rate ==
            m.rate_coef * m.interface_radius**2 * 2/m.d_particle**3
            )
m.reaction_rate_eqn = pyo.Constraint(rule=reaction_rate_rule)

m.comp_accum = pyo.Var(m.comp, initialize={'A': -1., 'B': 1.})
m.stoich = pyo.Param(m.comp, initialize={'A': -1, 'B': 1.})
def diff_rule(m, j):
    return m.comp_accum[j] == m.V_particle * m.stoich[j] * m.reaction_rate
m.diff_eqn = pyo.Constraint(m.comp, rule=diff_rule)

m.comp_holdup_prev = pyo.Var(m.comp, initialize={'A': 0., 'B': 4.})
m.comp_holdup_prev.fix()
m.comp_holdup.unfix()
m.delta_t = pyo.Param(initialize=1.)
def disc_rule(m, j):
    # Implicit Euler
    return m.comp_holdup[j] == m.comp_holdup_prev[j] + m.delta_t*m.comp_accum[j]
m.disc_eqn = pyo.Constraint(m.comp, rule=disc_rule)

try:
    solver.solve(m, tee=True)
except ApplicationError:
    pass

import pdb; pdb.set_trace()
