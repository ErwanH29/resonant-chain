#!/usr/bin/env python

# from amuse.community.rebound.interface import Rebound
from amuse.community.ph4.interface import ph4
from amuse.units import units, constants, nbody_system, quantities
from amuse.lab import Particles, Particle
from amuse.couple import bridge
from amuse.ext.orbital_elements import new_binary_from_orbital_elements, orbital_elements_from_binary

from tqdm import tqdm
import numpy as np
from amuse.io import write_set_to_file, read_set_from_file

import matplotlib.pyplot as plt
from fractions import Fraction

from generate_resonant_chain import bring_planet_pair_in_resonance

def semi_to_orbital_period(a, Mtot) :
    return 2*np.pi * (a**3/(constants.G*Mtot)).sqrt()

def orbital_period_to_semi(P, Mtot) :
    return ((constants.G*Mtot) * (P/(2*np.pi))**2)**(1./3.)

def resonant_chain_planetary_system(bodies):

    star = bodies[bodies.type == "star"][0]    
    planets = bodies[bodies.type == "planet"]
    for pi in range(len(planets)-1):
        bodies = resonant_pair_planetary_system(bodies, inner_planet_id=pi)
    return bodies

def resonant_pair_planetary_system(bodies, inner_planet_id=0, outer_planet_id=1):

    star = bodies[bodies.type == "star"][0]    
    planets = bodies[bodies.type == "planet"]
    
    for pi in planets:
        orbital_elements = orbital_elements_from_binary(star + pi)
        pi.semimajor_axis = orbital_elements[2]
        pi.eccentricity = orbital_elements[3]
        pi.inclination = orbital_elements[5]

    P1 = planets[inner_planet_id]
    P2 = planets[outer_planet_id]
    print(f"Resonance for {P1.name} and {P2.name}")
    Porb_a = semi_to_orbital_period(P1.semimajor_axis, star.mass + P1.mass)
    Porb_b = semi_to_orbital_period(P2.semimajor_axis, star.mass + P2.mass)
    fraction = Fraction(Porb_a/Porb_b).limit_denominator(10)
    print(f"{P1.name}, {P2.name} F={fraction}, {fraction.numerator/fraction.denominator},  {Porb_a/Porb_b}")
    print(star.name)
    bring_planet_pair_in_resonance(bodies, P1, P2)

    return bodies

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--inner", dest="inner_planet_id",
                      type="int",
                      default = 0,
                      help="inner planet id [%default]")
    result.add_option("--outer", dest="outer_planet_id",
                      type="int",
                      default = 1,
                      help="outer planet id [%default]")
    result.add_option("-f", dest="infilename", 
                      default = "input_filename.amuse",
                      help="input infilename [%default]")
    result.add_option("-F", dest="outfilename", 
                      default = "output_filename.amuse",
                      help="output infilename [%default]")
    return result
    
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    bodies = read_set_from_file(o.infilename)
    
    if o.inner_planet_id<0:
        bodies = resonant_chain_planetary_system(bodies)
    else:
        bodies = resonant_pair_planetary_system(bodies, o.inner_planet_id, o.outer_planet_id)
    write_set_to_file(bodies,
                      o.outfilename,
                      overwrite_file=True,
                      append_to_file=False,
                      close_file=True)
