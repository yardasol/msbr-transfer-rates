import openmc
import openmc.deplete as od
import json

model = openmc.Model.from_xml()
n_steps = 33
timesteps = n_steps * [3]
power = n_steps * [2.25e9]

with open('../serpent_fissq.json') as f:
    fission_q = json.load(f)

openmc.config['cross_sections'] = '/home/oleg/projects/cross-section-libraries/endfb71_h5/cross_sections.xml'
operator = od.CoupledOperator(model, "../chain_endfb71_ace.xml", fission_q=fission_q)
integrator = od.PredictorIntegrator(operator,  timesteps, power=power, timestep_units='d')

cycle_time_dict = {3: (['Pa'], 'd'),
                   50: (['Y', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Gd', 'Sm'], 'd'),
                   500: (['Eu'], 'd'),
                   20: (['Se', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Sb', 'Te'], 's'),
                   200: (['Zr', 'Cd', 'In', 'Sn'], 'd'),
                   20: (['Kr', 'Xe', 'H'], 's'),
                   60: (['Br', 'I'], 'd'),
                   3435: (['Sr', 'Ba', 'Rb', 'Cs'], 'd')
}
for c_time, (nucs, units) in cycle_time_dict.items():
    integrator.add_transfer_rate('fuel', nucs, 1/c_time, transfer_rate_units=f'1/{units}')

integrator.integrate()
