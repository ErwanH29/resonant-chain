#!/usr/bin/env python

# from amuse.community.rebound.interface import Rebound
import sys

from amuse.community.ph4.interface import ph4
from amuse.units import units, constants, nbody_system, quantities
from amuse.lab import Particles, Particle
from amuse.couple import bridge
from amuse.ext.orbital_elements import new_binary_from_orbital_elements, orbital_elements_from_binary

import glob
from tqdm import tqdm
import numpy as np
from amuse.io import write_set_to_file, read_set_from_file

from fractions import Fraction

from generate_resonant_chain import bring_planet_pair_in_resonance

def semi_to_orbital_period(a, Mtot) :
    a = abs(a)
    Mtot = abs(Mtot)
    return 2*np.pi * (a**3/(constants.G*Mtot)).sqrt()

def orbital_period_to_semi(P, Mtot) :
    return ((constants.G*Mtot) * (P/(2*np.pi))**2)**(1./3.)

def resonant_chain_planetary_system(bodies,
                                    tau_a_factor,
                                    t_integration, n_steps):

    star = bodies[bodies.type=="HOST"]
    star = star[star.mass.argmax()]
    planets = bodies[bodies.type=="PLANET"]
    print("resonant_chain_planetary_system")
    for pi in range(len(planets)-1):
        bodies = resonant_pair_planetary_system(bodies,
                                                inner_planet_id=pi,
                                                tau_a_factor=tau_a_factor,
                                                t_integration=t_integration,
                                                n_steps=n_steps,
                                                plot_results=False)
    print("resonant_chain_planetary_system: returning")

    return bodies

def resonant_pair_planetary_system(bodies, inner_planet_id=0, outer_planet_id=1,
                                   tau_a_factor=-1e5,
                                   t_integration=100, n_steps=100,
                                   plot_results=True):
    star = bodies[bodies.type=="HOST"]
    star = star[star.mass.argmax()]
    planets = bodies[bodies.type=="PLANET"]
    print("resonant_pair_planetary_system: for pi in planets")
    for pi in planets:
        orbital_elements = orbital_elements_from_binary(star + pi)
        pi.semimajor_axis = orbital_elements[2]
        pi.eccentricity = orbital_elements[3]
        pi.inclination = orbital_elements[5]

    print("resonant_pair_planetary_system: getting periods")
    P1 = planets[inner_planet_id]
    P2 = planets[outer_planet_id]
    Porb_a = semi_to_orbital_period(P1.semimajor_axis, star.mass + P1.mass)[0]
    Porb_b = semi_to_orbital_period(P2.semimajor_axis, star.mass + P2.mass)[0]
    if np.isnan(Porb_a.value_in(units.s)) or np.isnan(Porb_b.value_in(units.s)):
        return bodies

    print("resonant_pair_planetary_system: getting fractions")
    fraction = Fraction(Porb_a/Porb_b).limit_denominator(10)
    print(f"{P1.key}, {P2.key}, F={fraction}, {fraction.numerator/fraction.denominator},  {Porb_a/Porb_b}")
    bring_planet_pair_in_resonance(bodies, P1, P2,
                                   tau_a_factor,
                                   t_integration, n_steps,
                                   plot_results)
    print("resonant_pair_planetary_system: returning")
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
    result.add_option("--n_steps", dest="n_steps", 
                      default = 100, type="int",
                      help="number of steps [%default]")
    result.add_option("--t_integration", dest="t_integration", 
                      default = 1000, type="float",
                      help="integration time in units of the outer orbital period [%default]")
    result.add_option("--tau", 
                      dest="tau_a_factor", type="float", 
                      default = -1e5,
                      help="migration parameter (in terms of outer orbital period) [%default]")
    return result
    
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    from amuse.ext.orbital_elements import orbital_elements

    bodies = read_set_from_file("initial_particles/run_2.hdf5")
    bodies -= bodies[bodies.syst_id < 0]
    print(len(bodies), end="   ")
    bodies -= bodies[bodies.mass == 0 | units.MSun]
    print(len(bodies))

    start = int(sys.argv[1])
    end = int(sys.argv[2])
    stalled_df = [ ]
    #with open("to_ignore.txt", 'r') as f:
    #    f = f.readlines()
    #    f = f[0].strip().split(", ")
    #    for data in f:
    #        stalled_df.append(int(data))

    stable_cluster = Particles()
    for syst_id in np.unique(bodies.syst_id)[start:end]:
        data = glob.glob("systems/*.hdf5")
        if f"systems/resonant_system{syst_id}.hdf5" in data:
            print(f"System {syst_id} already simulated")
            continue
        elif syst_id in stalled_df:
            print(f"System {syst_id} stalls")
            continue
        
        if syst_id > 0:
            print(f"!!!!!!!!!!!!!!!!!!!! SYSTEM {syst_id}/{bodies.syst_id.max()} !!!!!!!!!!!")
            system = bodies[bodies.syst_id == syst_id]
            keep = np.ones(len(system), dtype=bool)
            print(f"System size: {len(system)}")
            if 0:#for i in range(len(system)):
                if system[i].mass > 0.08 | units.MSun:
                    continue

                if not keep[i]:
                    continue
                for j in range(i+1, len(system)):
                    if not keep[j]:
                        continue
                    if max(system[i].sma/system[j].sma, system[j].sma/system[i].sma) < 1.1 and system[j].sma < 10 | units.au:
                        system[i].mass += system[j].mass
                        keep[j] = False
                    elif system[j].sma < 1 | units.au and system[i].sma < 1 | units.au:
                        system[i].mass += system[j].mass
                        keep[j] = False


            system = system[keep]
            if len(system) > 2:
                original_pos = system[system.mass.argmax()].position
                original_vel = system[system.mass.argmax()].velocity
                copied_system = system.copy()

                #system = system.sorted_by_attribute("sma")
                system.position -= original_pos
                system.velocity -= original_vel
                host = system[system.mass.argmax()]
                p = system - host
                
                print(f"PLANET ORIGINAL SYSTEM")
                print(f"Masses: {system.mass}")
                for pl in p:
                    ke = orbital_elements(pl + host, G=constants.G)
                    print(ke[3], ke[2].in_(units.au), ke[5].in_(units.deg))
            
                system = resonant_chain_planetary_system(system, o.tau_a_factor, o.t_integration, o.n_steps)
                host = system[system.mass.argmax()]
                p = system - host
                print(f"AFTER SHAKING")
                for pl in p:
                    ke = orbital_elements(pl + host, G=constants.G)
                    print(ke[3], ke[2].in_(units.au), ke[5].in_(units.deg))
            
                print(f"BACK INTO CLUSTER")
                system.position -= host.position
                system.velocity -= host.velocity
                system.position += original_pos
                system.velocity += original_vel
                host = system[system.mass.argmax()]
                p = system - host
                for pl in p:
                    ke = orbital_elements(pl + host, G=constants.G)
                    print(ke[3], ke[2].in_(units.au), ke[5].in_(units.deg))

                print(f"PLANET ORIGINAL SYSTEM")
                print((original_pos-host.position).in_(units.pc))
                print((original_vel-host.velocity).in_(units.kms))
                print(":"*50)
                if np.isnan(host.x.value_in(units.pc)):
                    print("!!! NAN ERROR: Restoring OG !!!")
                    system = copied_system
                    
                    host = system[system.mass.argmax()]
                    p = system - host
                    for pl in p:
                        ke = orbital_elements(pl + host, G=constants.G)
                        print(ke[3], ke[2].in_(units.au), ke[5].in_(units.deg))
            
            write_set_to_file(system,
                      f"systems/resonant_system{syst_id}.hdf5",
                      overwrite_file=True,
                      append_to_file=False,
                      close_file=True)

