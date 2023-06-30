import openmc
import math
import numpy as np

def main(r):
    cyl = openmc.ZCylinder(r=r)
    gr_sq_d = 4.953 * 2
    gr_sq = openmc.rectangular_prism(gr_sq_d,
                                     gr_sq_d)
    bound = openmc.rectangular_prism(10.16, 10.16, boundary_type='reflective')
    h = 449.58 #cm
    core_base = openmc.ZPlane(z0=-h/2, boundary_type='vacuum')
    core_top = openmc.ZPlane(z0=h/2, boundary_type='vacuum')

    comps = np.array([1, 1, 1, 2, 1, 4, 1, 4])

    # Ref values
    #vals = np.array([71.7,71.7,16.0,16.0, 12.0,12.0, 0.3,0.3])
    vals = np.array([71.75,71.75,16.0,16.0, 12.0,12.0, 0.25,0.25])
    #thf4_conc = 12.09
    #diff = thf4_conc - 12.0
    #uf4_conc = 0.3 - diff
    #vals = np.array([71.7,71.7,16.0,16.0,thf4_conc,thf4_conc,uf4_conc,uf4_conc])

    nucs = (['Li', 'F', 'Be', 'F', 'Th', 'F', 'U', 'F'])
    vals = comps*vals
    tots = ([vals[0], vals[1] + vals[3] + vals[5] + vals[7], vals[2], vals[4], vals[6]])
    tots = tots / np.sum(tots) * 100
    nucs = [nucs[0], nucs[1], nucs[2], nucs[4], nucs[6]]

    # reference
    components = {'Li': {'percent': tots[0]/100,
                         'enrichment': 99.995,
                         'enrichment_target': 'Li7',
                         'enrichment_type': 'wo'},
                  'F19': tots[1]/100,
                  'Be9': tots[2]/100,
                  'Th232': tots[3]/100,
                  'U233': tots[4]/100}
    fuel = openmc.Material(name='fuel')
    fuel.add_components(components, percent_type='ao')

    fuel.set_density('g/cm3', density=3.35)
    fuel.depletable = True
    fuel.volume = h * (math.pi * r**2 + (10.16**2 - gr_sq_d**2))

    moder = openmc.Material(name='graphite')
    moder.set_density('g/cm3', density=1.84)
    moder.add_nuclide('C0', 1.000, percent_type='wo')
    moder.add_s_alpha_beta('c_Graphite')

    mats = openmc.Materials(materials=[fuel, moder])


    c1 = openmc.Cell(fill=fuel, region=-cyl & -core_top & +core_base)
    c2 = openmc.Cell(fill=moder, region=+cyl & gr_sq & -core_top & +core_base)
    c3 = openmc.Cell(fill=fuel, region=~gr_sq & bound & -core_top & +core_base)
    u = openmc.Universe(cells=[c1, c2, c3])
    geo = openmc.Geometry()
    geo.root_universe = u

    settings = openmc.Settings()
    settings.particles = 2000
    settings.batches = 50 
    settings.inactive = 10 
    settings.temperature = {'default': 900,
                            'method': 'interpolation',
                            'range': (800, 1000)}

    return openmc.Model(geo, mats, settings)



