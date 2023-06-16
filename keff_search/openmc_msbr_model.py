import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import openmc
from openmc.deplete import CoupledOperator, PredictorIntegrator, CELIIntegrator

import core_elements as ce
import control_rods as cr
import root_geometry as rg

def main(height):
    # Materials
    fuel = openmc.Material(name='fuel')
    fuel.set_density('g/cm3', density=3.35)
    fuel.add_components({'Li7': 0.0787474673879085,
                         'Be9': 0.0225566879138321,
                         'F19': 0.454003012179284,
                         'Th232': 0.435579130482336,
                         'U233': 0.00911370203663893},
                        percent_type='wo')
    fuel.depletable = True
    fuel.volume = 48710000.0

    moder = openmc.Material(name='graphite')
    moder.set_density('g/cm3', density=1.84)
    moder.add_nuclide('C0', 1.000, percent_type='wo')
    moder.add_s_alpha_beta('c_Graphite')

    hast = openmc.Material(name='hastelloyN')
    hast.set_density('g/cm3', density=8.671)
    hast.add_components({'Al27': 0.003,
                         'Ni': 0.677,
                         'W': 0.250,
                         'Cr': 0.070},
                        percent_type='wo')

    mat = openmc.Materials(materials=[fuel, moder, hast])

    colormap = {moder: 'purple',
                hast: 'blue',
                fuel: 'yellow'}

    def parse_arguments():
        """Parses arguments from command line.

        Returns
        -------
        deplete : bool
            Flag indicated whether or not to run a depletion simulation.
        volume_calculation : bool
            Flag indicating whether or not to run a stochastic volume calcuation.
        entropy : bool
            Flag indicating whether or not to calculate Shannon entropy.
        """
        parser = argparse.ArgumentParser()
        parser.add_argument('--deplete',
                            type=bool,
                            default=False,
                            help='flag for running depletion')
        parser.add_argument('--volume_calculation',
                            type=bool,
                            default=False,
                            help='flag for running stochastic volume calculation')
        parser.add_argument('--entropy',
                            type=bool,
                            default=False,
                            help='flag for including entropy mesh')



        args = parser.parse_args()
        return bool(args.deplete), bool(args.volume_calculation), bool(args.entropy)

    def shared_elem_geometry(elem_type='core',
                             gr_sq_d=4.953,
                             gr_sq_r=0.46,
                             r_rib=0.66802,
                             l1=4.28498,
                             l2=4.53898,
                             l3=5.62102,
                             r_es=2.2225,
                             es_name='gr_round_4'):
        """Creates surfaces and regions for lattice elements.

        Parameters
        ----------
        elem_type : 'core', 'cr'
            Indicates the type of element. 'core' inidcates for Zones IB and IIA.
            'cr' indicates control rod.
        gr_sq_d : float
            Half-width of graphite square element in cm.
        gr_sq_r : float
            Radius of graphite square rounded corners in cm.
        r_rib : float
            Radius of graphite element rib section.
        l1 : float
            Coordinate used to position graphite element ribs and rib tips.
        l2 : float
            Coordinate used to position graphite element ribs.
        l3 : float
            Coordinate used to position graphite element rib tips.
        r_es : float
            Radius of extra cylindrical surface used for element
        es_name : str
            Name of extra cylindrical surface.

        Returns
        -------
        gr_sq_neg : openmc.Intersection
            The region bounding the outer surface of the graphite
            element.
        gr_extra_regions : list of (openmc.Region, str)
            'Add-on' regions and their names for the graphite element.
            Includes ribs, rib tips, and gap-filling regions.
        inter_elem_channel : openmc.Region, list of (openmc.Region, str)
            Inter-element channel region(s)
        extra_surf : openmc.ZCylinder
            Extra cylindrical surface used in the element.
        """
        # Square graphite element
        gr_sq_neg = openmc.rectangular_prism(gr_sq_d*2,
                                             gr_sq_d*2,
                                             corner_radius=gr_sq_r)

        # Rib tip surfaces
        ul_t = openmc.ZCylinder(-l1, -l3, r_rib, name='rib_ul_tip')
        br_t = openmc.ZCylinder(l1, l3, r_rib, name='rib_br_tip')
        ru_t = openmc.ZCylinder(-l3, l1, r_rib, name='rib_ru_tip')
        lb_t = openmc.ZCylinder(l3, -l1, r_rib, name='rib_lb_tip')

        # Graphite element rib tip regions
        rib_ul_t = -ul_t
        rib_br_t = -br_t
        rib_ru_t = -ru_t
        rib_lb_t = -lb_t

        if elem_type == 'core':
            # Graphite element ribs for zones I-B and II-A
            ul = openmc.ZCylinder(-l1, l2, r_rib, name='rib_ul')
            br = openmc.ZCylinder(l1, -l2, r_rib, name='rib_br')
            lb = openmc.ZCylinder(-l2, -l1, r_rib, name='rib_lb')
            ru = openmc.ZCylinder(l2, l1, r_rib, name='rib_ru')

            # Graphite element rib regions.
            rib_ul = -ul
            rib_br = -br
            rib_lb = -lb
            rib_ru = -ru

            # inter-element fuel channel region
            inter_elem_channel = +ul & +br & +lb & +ru

        elif elem_type == 'cr':
            # Parameters for control rod element
            r_d = 1.16
            e_d = 2 * r_d / np.sqrt(3)
            r_c = 0.18

            # Base rib region
            ul = openmc.model.hexagonal_prism(origin=(-l1, l2), edge_length=e_d,
                                              orientation='x', corner_radius=r_c)
            br = openmc.model.hexagonal_prism(origin=(l1, -l2), edge_length=e_d,
                                              orientation='x',corner_radius=r_c)
            lb = openmc.model.hexagonal_prism(origin=(-l2, -l1), edge_length=e_d,
                                              orientation='y',corner_radius=r_c)
            ru = openmc.model.hexagonal_prism(origin=(l2, l1), edge_length=e_d,
                                              orientation='y',corner_radius=r_c)

            rib_ul = ul
            rib_lb = br
            rib_br = lb
            rib_ru = ru

            inter_elem_channel = ~ul & ~br & ~lb & ~ru

        ribs = [[rib_ul, 'rib_ul'],
                [rib_br, 'rib_br'],
                [rib_ru, 'rib_ru'],
                [rib_lb, 'rib_lb'],
                [rib_ul_t, 'rib_ul_t'],
                [rib_br_t, 'rib_br_t'],
                [rib_ru_t, 'rib_ru_t'],
                [rib_lb_t, 'rib_lb_t']]

        gr_extra_regions = ribs
        inter_elem_channel = inter_elem_channel & +ul_t & +br_t & +ru_t & +lb_t

        extra_surf = openmc.ZCylinder(r=r_es, name=es_name)

        return gr_sq_neg, gr_extra_regions, inter_elem_channel, extra_surf

    def cr_lattice(cr_boundary, core_base, core_top, height=0.0):
        """Creates the control rod lattice.

        Parameters
        ----------
        cr_boundary : openmc.Intersection
            Outer bound of the lattice in the xy-plane.
        core_base : openmc.ZPlane
            Core bottom bounding surface.
        core_top : openmc.ZPlane
            Core top bounding surface.

        Returns
        -------
        c1 : openmc.Cell
            Cell containing the control rod lattice.

        """
        (gr_sq_neg,
         gr_extra_regions,
         inter_elem_channel,
         fuel_hole) = shared_elem_geometry(elem_type='cr',
                                           gr_sq_d=7.23646,
                                           gr_sq_r=0.99,
                                           r_rib=0.8,
                                           l1=5.8801,
                                           l2=6.505,
                                           l3=8.03646,
                                           r_es=5.08,
                                           es_name='cr_fuel_hole')
        f = cr.control_rod(gr_sq_neg,
                           gr_extra_regions,
                           inter_elem_channel,
                           fuel_hole,
                           fuel,
                           moder,
                           height=height)
        e = cr.control_rod_channel(gr_sq_neg,
                                   gr_extra_regions,
                                   inter_elem_channel,
                                   fuel_hole,
                                   fuel,
                                   moder)

        cl = openmc.RectLattice()
        cl.pitch = np.array([15.24, 15.24])
        N = 2 / 2
        cl.lower_left = -1 * cl.pitch * N
        cl.universes = [[f, e],
                        [e, f]]
        c1 = openmc.Cell(fill=cl, region=(+core_base & -core_top & cr_boundary),
                         name='cr_lattice')

        return c1

    def main_lattice(zone_i_boundary, cr_boundary, core_base, core_top):
        """Creates the core lattice.

        Parameters
        ----------
        zone_i_boundary : 3-tuple of openmc.model.IsogonalOctagon
            Zone I bounding surfaces in the xy-plane
        cr_boundary : openmc.Intersection
            Outer bound of the lattice in the xy-plane.
        core_base : openmc.ZPlane
            Core bottom bounding surface.
        core_top : openmc.ZPlane
            Core top bounding surface.

        Returns
        -------
        main_cells : list of openmc.Cell
            Cells containing the main lattice.

        """
        (gr_sq_neg,
         gr_extra_regions,
         inter_elem_channel,
         gr_round_4) = shared_elem_geometry()

        l = ce.zoneIB(gr_sq_neg,
                      gr_extra_regions,
                      inter_elem_channel,
                      gr_round_4,
                      moder,
                      fuel,
                      hast)

        (gr_sq_neg,
         gr_extra_regions,
         inter_elem_channel,
         gr_round_4) = shared_elem_geometry()

        z = ce.zoneIIA(gr_sq_neg,
                       gr_extra_regions,
                       inter_elem_channel,
                       gr_round_4,
                       moder,
                       fuel)
        v = ce.void_cell()
        # tres, uno, dos, quatro
        t, u, d, q = ce.graphite_triangles(fuel, moder)

        s1, s2, s3 = zone_i_boundary

        main = openmc.RectLattice()
        main.pitch = np.array([10.16, 10.16])
        N = 45 / 2
        main.lower_left = -1 * main.pitch * N
        main.universes = [[v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, d, z, z, z, z, z, z, z, z, z, u, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v],
                          [v, v, v, v, v, v, v, v, v, v, v, v, v, d, z, z, z, z, l, l, l, l, l, l, l, l, l, z, z, z, z, u, v, v, v, v, v, v, v, v, v, v, v, v, v],
                          [v, v, v, v, v, v, v, v, v, v, v, d, z, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, z, u, v, v, v, v, v, v, v, v, v, v, v],
                          [v, v, v, v, v, v, v, v, v, v, d, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, u, v, v, v, v, v, v, v, v, v, v],
                          [v, v, v, v, v, v, v, v, v, d, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, u, v, v, v, v, v, v, v, v, v],
                          [v, v, v, v, v, v, v, v, d, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, u, v, v, v, v, v, v, v, v],
                          [v, v, v, v, v, v, v, d, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, u, v, v, v, v, v, v, v],
                          [v, v, v, v, v, v, d, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, u, v, v, v, v, v, v],
                          [v, v, v, v, v, d, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, u, v, v, v, v, v],
                          [v, v, v, v, d, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, u, v, v, v, v],
                          [v, v, v, d, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, u, v, v, v],
                          [v, v, d, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, u, v, v],
                          [v, v, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, v, v],
                          [v, d, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, u, v],
                          [v, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, v],
                          [v, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, v],
                          [v, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, v],
                          [d, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, u],
                          [z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z],
                          [z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z],
                          [z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z],
                          [z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, v, v, v, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z],
                          [z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, v, v, v, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z],
                          [z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, v, v, v, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z],
                          [z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z],
                          [z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z],
                          [z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z],
                          [t, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, q],
                          [v, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, v],
                          [v, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, v],
                          [v, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, v],
                          [v, t, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, q, v],
                          [v, v, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, v, v],
                          [v, v, t, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, q, v, v],
                          [v, v, v, t, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, q, v, v, v],
                          [v, v, v, v, t, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, q, v, v, v, v],
                          [v, v, v, v, v, t, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, q, v, v, v, v, v],
                          [v, v, v, v, v, v, t, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, q, v, v, v, v, v, v],
                          [v, v, v, v, v, v, v, t, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, q, v, v, v, v, v, v, v],
                          [v, v, v, v, v, v, v, v, t, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, q, v, v, v, v, v, v, v, v],
                          [v, v, v, v, v, v, v, v, v, t, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, q, v, v, v, v, v, v, v, v, v],
                          [v, v, v, v, v, v, v, v, v, v, t, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, q, v, v, v, v, v, v, v, v, v, v],
                          [v, v, v, v, v, v, v, v, v, v, v, t, z, z, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, z, z, q, v, v, v, v, v, v, v, v, v, v, v],
                          [v, v, v, v, v, v, v, v, v, v, v, v, v, t, z, z, z, z, l, l, l, l, l, l, l, l, l, z, z, z, z, q, v, v, v, v, v, v, v, v, v, v, v, v, v],
                          [v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, t, z, z, z, z, z, z, z, z, z, q, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v]]

        c1 = openmc.Cell(fill=main, region=(+core_base &
                                            -core_top &
                                            +zone_i_boundary[2] &
                                            -zone_i_boundary[1] &
                                            ~cr_boundary), name='main_lattice_smaller_octader')
        c2 = openmc.Cell(fill=main, region=(+core_base &
                                            -core_top &
                                            -zone_i_boundary[2] &
                                            ~cr_boundary),
                         name='main_lattice_insite_smallest_octader')
        c3 = openmc.Cell(fill=main, region=(+core_base &
                                            -core_top &
                                            -zone_i_boundary[0] &
                                            +zone_i_boundary[1] &
                                            +zone_i_boundary[2] &
                                            ~cr_boundary),
                         name=('main_lattice_inside_base_octader'
                              '_deducted_smaller_smallest'))
        main_cells = [c1, c2, c3]
        return main_cells

    deplete, volume_calculation, entropy = parse_arguments()

    (zone_bounds,
     core_bounds,
     reflector_bounds,
     vessel_bounds) = rg.shared_root_geometry()

    (cr_boundary,
     zone_i_boundary,
     zone_ii_boundary) = zone_bounds

    (annulus_boundary,
     lower_plenum_boundary,
     core_base,
     core_top) = core_bounds

    (radial_reflector_boundary,
     bottom_reflector_boundary,
     top_reflector_boundary) = reflector_bounds

    (radial_vessel_boundary,
     bottom_vessel_boundary,
     top_vessel_boundary) = vessel_bounds

    main = main_lattice(zone_i_boundary,
                        cr_boundary,
                        core_base,
                        core_top)

    cr_lat = cr_lattice(cr_boundary,
                    core_base,
                    core_top,
                    height=height)

    iib = rg.zoneIIB(zone_i_boundary,
                     zone_ii_boundary,
                     core_base,
                     core_top,
                     fuel,
                     moder)

    a = rg.annulus(zone_ii_boundary,
                   annulus_boundary,
                   core_base,
                   core_top,
                   fuel)

    lp = rg.lower_plenum(core_base,
                         lower_plenum_boundary,
                         annulus_boundary,
                         fuel)

    rr, rb, rt = rg.reflectors(annulus_boundary,
                               radial_reflector_boundary,
                               lower_plenum_boundary,
                               bottom_reflector_boundary,
                               core_top,
                               top_reflector_boundary,
                               moder)

    vr, vb, vt = rg.vessel(radial_reflector_boundary,
                           radial_vessel_boundary,
                           bottom_vessel_boundary,
                           top_vessel_boundary,
                           top_reflector_boundary,
                           bottom_reflector_boundary,
                           hast)

    geo = openmc.Geometry()
    univ = openmc.Universe(cells=[cr_lat, lp, a, rr, rb, rt, vr, vb, vt])
    univ.add_cells(main)
    univ.add_cells(iib)

    geo.root_universe = univ
    geo.remove_redundant_surfaces()

# Settings
    settings = openmc.Settings()
    settings.particles = 60000
    settings.batches = 200
    settings.inactive = 80
    settings.temperature = {'default': 900,
                            'method': 'interpolation',
                            'range': (800, 1000)}

    ll, ur = geo.root_universe.bounding_box
    if volume_calculation:
        msbr_volume_calc = openmc.VolumeCalculation([fuel, moder], int(1e10), ll, ur)
        #msbr_volume_calc.set_trigger(1e-03, 'rel_err')
        settings.volume_calculations = [msbr_volume_calc]
        settings.run_mode = 'volume'
    if entropy:
        entropy_mesh = openmc.RegularMesh()
        entropy_mesh.lower_left = ll
        entropy_mesh.upper_right = ur
        entropy_mesh.dimension = (20, 20, 20)
        settings.entropy_mesh = entropy_mesh

    return openmc.Model(geo, mat, settings)
