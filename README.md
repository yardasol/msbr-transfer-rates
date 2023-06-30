# msbr-transfer-rates
This repository holds various python scripts for simulating a CSG model of the MSBR.
The MSBR geometry is extremely complicated, so it is strongly reccomended to use the default
particle and batch settings, as they will give convergence in both the neutron source distribution
and in Keff.

To reproduce the results, use the cross section libary created by running the `download_endbf71.bash`
and `process_endfb71_to_openmc.bash` scripts found [here](https://github.com/arfc/saltproc/tree/master/scripts)
. The models also require using fission-q values from serpent, but these
are automatically provided in the `serpent_fissq.json` file.

Some of the scripts require downloading the [2022-yardas-ms](https://github.com/arfc/2022-yardas-ms/tree/master) repository in the same root directory.

Model Parameters:
- Material temperature: 900K
- Cross section library: ENDF B/VII.1

## Directory structure

### `olek-work-summary.ipynb`
Summary of the results and analysis from the various cases.

### `Th232-U233-cycle-times.ipynb`
Calculation of cycle times from SaltProc results.

### `compare_with_saltproc.ipynb`
Comparison of the results fromn `with_feeding` and `no_feeding` with SaltProc results.

### `with_feeding`
SaltProc MSBR model with reprocessing and feeding. Run the `run_with_transfer_rates.py`
script to reproduce the results.

### `no_feeds`
SaltProc MSBR model with reproecssing and no feeding. Run the `run_with_transfer_rates.py`
script to reproduce the results.

### `li6`
SaltProc MSBR model with 99.995\% enriched Li. Run the `openmc_msbr_model.py` script
to generate the model.

### `cr_holes`
Supereset of `li6` but with channels in the control rods. Run the `openmc_msbr_model.py` script
to generate the model. You can use the `openmc-plotter` tool to inspect the geometry. Adjust the 
arugment of the `main` function in `openmc_msbr_model.py` to adjust the control rod height.

### `absorber_rods`
Superset of `cr_holes` with the addition boron absorbing rods.

### `neutron_mesh`
Superset of `absorber_rods`, but we add a mesh to tally flux and absorption reactions and visualize it.
Use the `msbr_flux_plot` notebook to visiualize the flux and neutron absorptions in the core.

### `salt_cylinder`
MSBR pincell model. Used to investigate reactivity control by varying the amount of feed material.
Requires using the `batchwise` branch of the `openmsr` fork of OpenMC.
