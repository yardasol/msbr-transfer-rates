## Directory structure

### `with_feeds`
SaltProc MSBR model with reprocessing and feeding

### `no_feeds`
SaltProc MSBR model with reproecssing and no feeding

### `li6`
SaltProc MSBR model with 99.995\% enriched Li

### `cr_holes`
Supereset of `li6` but with channels in the control rods

### `absorber_rods`
Superset of `cr_holes` but with boron absorbing rods 

### `neutron_mesh`
Superset of `absorber_rods`, but we add a mesh to tally flux and absorption reactions and visualize it.

### `thorium_search`
Reactivity control by controlling thorium concentration in batchwise refueling.
