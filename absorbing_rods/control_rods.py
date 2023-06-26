import openmc
import numpy as np

def control_rod(gr_sq_neg,
                gr_extra_regions,
                inter_elem_channel,
                fuel_hole,
                fuel,
                moder):
    """Create universe for control rod element with control rod fully inserted.
    Based on specification in Roberton, 1971.

    Parameters
    ----------
    gr_sq_neg : openmc.Intersection
        The region bounding the outer surface of the 6 in. x 6 in. graphite
        element.
    gr_extra_regions : list of (openmc.Region, str)
        'Add-on' regions and their names for the graphite element.
        Includes ribs, rib tips, and gap-filling regions.
    inter_elem_channel : openmc.Region, list of (openmc.Region, str)
        Inter-element channel region(s).
    fuel_hole : openmc.ZCylinder
        Central fuel hole in graphite element.
    fuel : openmc.Material
        Fuel salt material
    moder : openmc.Material
        Graphite material

    Returns
    -------
    cr : openmc.Universe
        Univerese for control rod element with control rod fully insterted.
    """

    height = 0.01
    lower_bound = openmc.ZPlane(z0=height, name='control_rod_lower_bound')
    s1 = openmc.ZCylinder(r=4.7625, name='control_rod')
    s1b = openmc.ZCylinder(r=1.71069, name='control_rod_hole')
    c1 = openmc.Cell(fill=moder, region=-s1 & +s1b & +lower_bound, name='control_rod')
    c1b = openmc.Cell(fill=fuel, region=-s1b & +lower_bound, name='control_rod_fuel')
    c2 = openmc.Cell(fill=fuel, region=(+s1 & -fuel_hole & +lower_bound), name='cr_fuel_inner')
    c2b = openmc.Cell(fill=fuel, region=(-fuel_hole & -lower_bound))

    #c1 = openmc.Cell(fill=moder, region=-s1 & +s1b, name='control_rod')
    #c1b = openmc.Cell(fill=fuel, region=-s1b, name='control_rod_fuel')
    #c2 = openmc.Cell(fill=fuel, region=(+s1 & -fuel_hole), name='cr_fuel_inner')
    c3 = openmc.Cell(fill=moder, region=(+fuel_hole &
                                         gr_sq_neg &
                                         inter_elem_channel),
                     name='cr_moderator')
    c4 = openmc.Cell(fill=fuel, region= (~gr_sq_neg & inter_elem_channel),
                     name='cr_fuel_outer')
    #universe_id=3
    cr = openmc.Universe(name='control_rod', cells=[c1, c1b, c2, c2b, c3, c4])
    #cr = openmc.Universe(name='control_rod', cells=[c1, c1b, c2, c3, c4])

    for (reg, name) in gr_extra_regions:
            cr.add_cell(openmc.Cell(fill=moder, region=reg,
                                    name=f'cr_moderator_{name}'))

    return cr

def absorbing_rod(gr_sq_neg,
                gr_extra_regions,
                inter_elem_channel,
                fuel_hole,
                fuel,
                moder,
                b4c,
                height=0.0):
    """Create universe for absorbing rod element

    Parameters
    ----------
    gr_sq_neg : openmc.Intersection
        The region bounding the outer surface of the 6 in. x 6 in. graphite
        element.
    gr_extra_regions : list of (openmc.Region, str)
        'Add-on' regions and their names for the graphite element.
        Includes ribs, rib tips, and gap-filling regions.
    inter_elem_channel : openmc.Region, list of (openmc.Region, str)
        Inter-element channel region(s).
    fuel_hole : openmc.ZCylinder
        Central fuel hole in graphite element.
    fuel : openmc.Material
        Fuel salt material
    b4c : openmc.Material
        Absorbing material

    Returns
    -------
    ar : openmc.Universe
        Univerese for absorbing rod element.
    """

    lower_bound = openmc.ZPlane(z0=height, name='control_rod_lower_bound')
    s1 = openmc.ZCylinder(r=4.7625, name='absorbing_rod')
    c1 = openmc.Cell(fill=b4c, region=-s1 & +lower_bound, name='absorbing_rod')
    c2 = openmc.Cell(fill=fuel, region=(+s1 & -fuel_hole & +lower_bound), name='ar_fuel_inner')
    c2b = openmc.Cell(fill=fuel, region=(-fuel_hole & -lower_bound))
    c3 = openmc.Cell(fill=moder, region=(+fuel_hole &
                                         gr_sq_neg &
                                         inter_elem_channel),
                     name='ar_moderator')
    c4 = openmc.Cell(fill=fuel, region= (~gr_sq_neg & inter_elem_channel),
                     name='ar_fuel_outer')
    #universe_id=3
    ar = openmc.Universe(name='absorbing_rod', cells=[c1, c2, c2b, c3, c4])

    for (reg, name) in gr_extra_regions:
            ar.add_cell(openmc.Cell(fill=moder, region=reg,
                                    name=f'ar_moderator_{name}'))

    return ar



def control_rod_channel(gr_sq_neg,
                        gr_extra_regions,
                        inter_elem_channel,
                        fuel_hole,
                        fuel,
                        moder):
    """Create universe for control rod element with control rod fully withdrawn.
    Based on specification in Roberton, 1971.

    Parameters
    ----------
    gr_sq_neg : openmc.Intersection
        The region bounding the outer surface of the 6 in. x 6 in. graphite
        element.
    gr_extra_regions : list of (openmc.Region, str)
        'Add-on' regions and their names for the graphite element.
        Includes ribs, rib tips, and gap-filling regions.
    inter_elem_channel : openmc.Region, list of (openmc.Region, str)
        Inter-element channel region(s).
    fuel_hole : openmc.ZCylinder
        Central fuel hole in graphite element.
    fuel : openmc.Material
        Fuel salt material
    moder : openmc.Material
        Graphite material

    Returns
    -------
    crc : openmc.Universe
        Universe for control rod element with control rod fully withdrawn.
    """

    c1 = openmc.Cell(fill=fuel, region=(-fuel_hole), name='crc_fuel_inner')

    c2 = openmc.Cell(fill=moder, region=(+fuel_hole &
                                         gr_sq_neg &
                                         inter_elem_channel),
                     name='crc_moderator')
    c3 = openmc.Cell(fill=fuel, region=(~gr_sq_neg & inter_elem_channel),
                     name='crc_fuel_outer')

    # universe_id=4
    crc = openmc.Universe(name='control_rod_channel', cells=[c1, c2, c3])

    for (reg, name) in gr_extra_regions:
        crc.add_cell(openmc.Cell(fill=moder, region=reg,
                                 name=f'crc_moderator_{name}'))

    return crc
