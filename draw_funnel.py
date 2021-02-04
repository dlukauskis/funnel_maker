# Script for drawing a funnel in PyMol from the origin point along a vector parallel to the input one.
#
# TIP: load it automatically by adding this line to your ~/.pymolrc file
# load {path_to_user_scripts}/draw_funnel.py
#
# USAGE:
#
# draw_funnel p1 [, p2 [, s_cent [, beta_cent [, wall_width,
#     [, wall_buffer [, lower_wall [, upper_wall [, vec_step [, angle_sample ]]]]]]]]]
#
# ARGUMENTS:
#
# p1 = string, origin point of a funnel vector
# p2 = string, terminal point of a funnel vector
# s_cent = float, inflection point of the sigmoid function
# beta_cent = float, steepness of the sigmoid function
# wall_width = float, the radius of the funnel at the widest point (excluding the buffer)
# wall_buffer = float, the radius of the funnel at the narrowest point
# lower_wall = float, the starting point of the funnel from the origin along the selected vector (in A)
# upper_wall = float, the ending point of the funnel from the origin along the selected vector (in A)
# vec_step = float, determines how often the funnel points are drawn along the selected vector (in A)
# angle_sample = float, sets the angle between adjacent points at each funnel ring

# by Antonija (based on Giulio's script)
# edited by Vladas for Python 3 compatibility.
#
# edited again by Dom, to only use 2 points to define the vector, just as used in PLUMED and openMM scripts

import numpy as np
from pymol import cmd

def draw_funnel(p1, p2, s_cent=16, beta_cent=0.3, wall_width=15.5, wall_buffer=1.5, lower_wall=0, upper_wall=32, vec_step=2.5, angle_sample=18):
    s_cent = float(s_cent)
    beta_cent = float(beta_cent)
    wall_width = float(wall_width)
    wall_buffer = float(wall_buffer)
    lower_wall = float(lower_wall)
    upper_wall = float(upper_wall)
    vec_step = float(vec_step)
    angle_sample = float(angle_sample)
    
    # get coords of the origin and vector points
    #origin = cmd.get_coords(selection, 1)[0]
    origin = cmd.get_coords(p1, 1)[0]
    print('Origin:', origin)
    
    v1 = cmd.get_coords(p1, 1)[0]
    v2 = cmd.get_coords(p2, 1)[0]
    # calculate the vector defined by points p1 and p2
    vec = np.array(v2, dtype=float) - np.array(v1, dtype=float)
    # BEWARE: inconsistency with linalg, if vec is a list and not an array!!!
#    print(np.linalg.norm(vec), np.linalg.norm(v2 - v1))
    # make it a unit vector
    unit_vec = vec/np.linalg.norm(vec)
#    print(np.linalg.norm(vec), vec, np.linalg.norm(unit_vec), unit_vec) 
    # how to get orthogonal vectors
    # https://math.stackexchange.com/questions/133177/finding-a-unit-vector-perpendicular-to-another-vector
    # determine 1st orthogonal vector
    a0 = np.random.randint(1,10)
    a1 = np.random.randint(1,10)
    a2 = -(a0*vec[0] + a1*vec[1])/vec[2]
    a = np.asarray([a0, a1, a2])
    unit_a = a/np.linalg.norm(a)
    # determine 2nd orthogonal vector
    unit_b = np.cross(unit_a, unit_vec)
#    print(unit_vec, unit_a, unit_b)
#    print(np.linalg.norm(unit_vec), np.linalg.norm(unit_a), np.linalg.norm(unit_b))
    # iterate along the selected vector
    for step in np.arange(lower_wall, upper_wall, vec_step):
        # iterate around a circle with its radius defined by the sigmoid function
        radius = (wall_width / (1 + np.exp(beta_cent * (step - s_cent)))) + wall_buffer
        for angle in np.arange(-np.pi, np.pi, 2 * np.pi / angle_sample):
            # calculate parametric functions for this specific case
            # https://math.stackexchange.com/questions/73237/parametric-equation-of-a-circle-in-3d-space
            # generate pseudoatoms along the axis
            pos = origin + unit_vec*step + radius*(np.cos(angle)*unit_a + np.sin(angle)*unit_b)
            cmd.pseudoatom('funnel', pos=pos.tolist())

    cmd.color('orange', selection='funnel')
    cmd.show_as('nonbonded', 'funnel')
    #cmd.show('spheres', 'not funnel and ' + selection)

cmd.extend("draw_funnel", draw_funnel)
