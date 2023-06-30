import openmc
import json
from cylinder import main
import openmc.deplete as od
from pathlib import Path

r = 1.71069 # ib
#r = 0.762 # ia
base_power = 2.25e9
#linear_power_dens = base_power / (449.58) / 1436
linear_power = base_power / 1436
model = main(r)
n_steps = 10 
timesteps = n_steps * [3]
#power = n_steps * [linear_power_dens]
power = n_steps * [linear_power]

with open('../serpent_fissq.json') as f:
    fission_q = json.load(f)

openmc.config['cross_sections'] = '/home/oleg/projects/cross-section-libraries/endfb71_h5/cross_sections.xml'
cycle_time_dict = {3: (['Pa'], 'd'),
                   50: (['Y', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Gd', 'Sm'], 'd'),
                   500: (['Eu'], 'd'),
                   20: (['Se', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Sb', 'Te'], 's'),
                   200: (['Zr', 'Cd', 'In', 'Sn'], 'd'),
                   20: (['Kr', 'Xe', 'H'], 's'),
                   60: (['Br', 'I'], 'd'),
                   3435: (['Sr', 'Ba', 'Rb', 'Cs'], 'd')
                   }
                   #3435: (['Th', 'Li', 'Be', 'F'], 'd'),
                   #16: (['Np', 'Pu'], 'a')
#}

#res = od.Results('depletion_results.h5')
#operator = od.CoupledOperator(model, "../chain_endfb71_ace.xml", fission_q=fission_q)#, prev_results=res)
#integrator = od.PredictorIntegrator(operator,  timesteps[0:7], power=power[0:7], timestep_units='d')

#for c_time, (nucs, units) in cycle_time_dict.items():
#    integrator.add_transfer_rate('fuel', nucs, 1/c_time, transfer_rate_units=f'1/{units}')

#integrator.integrate()
paths = ['depletion_results.h5', 'summary.h5', 'tallies.out']
for path in paths:
    p = Path(path)
    p.rename(f'ref_{path}')
f = 'openmc_simulation_n'
for n in range(0,8):
    p = Path(f'{f}{n}.h5')
    p.rename(f'ref_{f}{n}.h5')

res = od.Results('ref_depletion_results.h5')

operator = od.CoupledOperator(model, "../chain_endfb71_ace.xml", fission_q=fission_q, prev_results=res)
integrator = od.PredictorIntegrator(operator,  timesteps, power=power, timestep_units='d')

for c_time, (nucs, units) in cycle_time_dict.items():
    integrator.add_transfer_rate('fuel', nucs, 1/c_time, transfer_rate_units=f'1/{units}')

integrator.add_batchwise('refuel',
        mats_id_or_name=['fuel'],
        mat_vector={'U233': 0.45, 'Th232': 0.55}, 
        bracket=[0, 100],
        bracket_limit=[0, 6500],
        tol=10)

integrator.integrate()
