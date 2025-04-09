"""
   Minimalistic routine for running a gravity code
"""
from amuse.lab import *
from amuse.couple import bridge
from amuse.units import nbody_system
from amuse.ic.kingmodel import new_king_model
from amuse.community.fractalcluster.interface import new_fractal_cluster_model

from amuse.community.symple.interface import symple	# symplectic
from amuse.community.huayno.interface import Huayno	# symplectic

from amuse.ext.orbital_elements import get_orbital_elements_from_arrays
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ic import make_planets_oligarch
from amuse.ext import solarsystem 

from matplotlib import pyplot as plt
import numpy as np
import math

from fractions import Fraction
from collections import Counter

def orbital_period(a, Mtot) :
    return 2*np.pi*(a**3/(constants.G*Mtot)).sqrt()

def nearest_integer_ratio(float_num, max_numerator=20, max_denominator=20):
    # Convert the float to a Fraction object
    original_frac = Fraction(float_num).limit_denominator(max_denominator)
    
    # Initialize best fraction and error
    best_frac = Fraction(1, 1)
    best_error = abs(float_num - 1)

    # Iterate through all possible fractions within the limits
    for denominator in range(1, max_denominator + 1):
        for numerator in range(1, max_numerator + 1):
            frac = Fraction(numerator, denominator)
            error = abs(float_num - frac)
            
            if error < best_error:
                best_frac = frac
                best_error = error

    return best_frac.numerator, best_frac.denominator

def plot_numerator_denominator_pairs(numbers, max_numden=20):
    pairs = [nearest_integer_ratio(num,
                                   max_numerator=max_numden,
                                   max_denominator=max_numden) for num in numbers]
    pair_counts = Counter(pairs)
    n_inner = 0
    n_outer = 0
    for pj, pi in pairs:
        print(pj, pi)
        if pi>pj:
            n_inner += 1
        else:
            n_outer += 1
    print(f"Inner resonances: {n_inner}, outer resonances: {n_outer} ratio inner/outer= {n_inner/n_outer}")
    #print(pair_counts)
    print("N=", sum(pair_counts.values()), pair_counts[(max_numden, 1)], pair_counts[(1, max_numden)])
    n_rgtn = pair_counts[(max_numden, 1)]
    n_rltn = pair_counts[(1, max_numden)]
    #for pi in pair_counts:
    #    print("x=", pi)

    del pair_counts[(max_numden, 1)]
    del pair_counts[(1, max_numden)]
    #print("N=", len(pair_counts))
    #print(pair_counts)
    
    fig, ax = plt.subplots(figsize=(6, 6))

    max_count = max(pair_counts.values())

    for (num, den), count in pair_counts.items():
        size = 50*(count / max_count) * max_numden  # Adjust 500 to change the maximum square size
        color = [(num-1)/(max_numden-1), 0, (den-1)/(max_numden-1)]  # Red for numerator, Blue for denominator
        ax.scatter(num, den, s=size, c=[color], alpha=0.6, edgecolors='black')

    x = [1,max_numden]
    plt.plot(x,x, ls=":", lw=1, c='k')
    x = [1,max_numden]
    y = [1,max_numden/3]
    plt.plot(x,y, ls=":", lw=1, c='r')
    plt.plot(y,x, ls=":", lw=1, c='b')

    n_tot = n_inner + n_outer + n_rgtn + n_rltn
    plt.text(max_numden-3.5, max_numden-1.0, f"$N_{{i>j}}=${n_inner/n_tot:.4f}")
    plt.text(max_numden-3.5, max_numden-2.0, f"$N_{{i<j}}=${n_outer/n_tot:.4f}")
    plt.text(max_numden-3.5, max_numden-3.0, f"$N_{{inner}}=${n_rgtn/n_tot:.4f}")
    plt.text(max_numden-3.5, max_numden-4.0, f"$N_{{outer}}=${n_rltn/n_tot:.4f}")
    
    ax.set_xlim(0.5, max_numden+0.5)
    ax.set_ylim(0.5, max_numden+0.5)
    ax.set_xticks(range(1, max_numden+1))
    ax.set_yticks(range(1, max_numden+1))
    ax.set_xlabel('i (outer)')
    ax.set_ylabel('j (inner)')
    #ax.set_title('resonance i/j Pairs Distribution')
    ax.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.show()


def print_reonances(bodies, resonant_body, max_numden=20):

    bodies = bodies - resonant_body
    sun = bodies[bodies.type=="star"]
    planets = bodies - sun
    
    Mtot = sun.mass.sum()
    Pplanet = orbital_period(planets.semimajor_axis, Mtot)
    a_planet = resonant_body.semimajor_axis
    Pplanets = orbital_period(planets.semimajor_axis, Mtot)
    Pplanets = Pplanets/orbital_period(resonant_body.semimajor_axis, Mtot)
    
    pairs = [nearest_integer_ratio(num,
                                   max_numerator=max_numden,
                                   max_denominator=max_numden) for num in Pplanets]
    print(f"{resonant_body.name} forms resonant pairs with {pairs}")
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="fcluster", default = "",
                      help="input filename [%default]")
    return result
    
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    bodies = read_set_from_file(o.fcluster, close_file=True)
    print(bodies)
    bodies.semimajor_axis = 0 | units.au
    #for pi in reversed(list(bodies.iter_history())):
    #    print("time=", pi[0].age.in_(units.Myr))
    #pi = list(bodies.iter_history())[-1]
    #print(f"Time={pi[0].age.in_(units.Myr)}")

    sun = bodies[bodies.type=="star"]
    planets = bodies - sun
    for pi in planets:
        orbit = orbital_elements_from_binary(sun+pi, G=constants.G)
        pi.semimajor_axis = orbit[2]
        pi.eccentricity = orbit[3]
        pi.true_anomal = orbit[4]
        pi.inclination = orbit[5]
        pi.long_of_ascending_node = orbit[6]
        pi.arg_of_periapse = orbit[7]
    Pplanet = orbital_period(planets.semimajor_axis, sun.mass+planets.mass)

    for pi in range(len(planets)):
        Porb = Pplanet/Pplanet[pi]
        print(f"Number reonances with {planets[pi].name}: {len(Porb)} of which inner {len(Porb[Porb<=1])} and outer {len(Porb[Porb>1])}.")

    max_numden = 5
    for pi in planets:
        print_reonances(bodies, pi, max_numden=max_numden)    
    #plot_numerator_denominator_pairs(Pplanet, max_numden=max_numden)

    
