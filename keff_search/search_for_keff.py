import openmc
from openmc_msbr_model import main

crit_cr_height, guesses, keffs = openmc.search_for_keff(main, bracket=[0.1, 449.57], tol=1e-2, print_iterations=True, run_args={'output': False})
print('Critical CR Height: {:4.0f} cm'.format(crit_cr_height))
