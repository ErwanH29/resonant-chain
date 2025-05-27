#!/usr/bin/env python


from fractions import Fraction

from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.io import write_set_to_file, read_set_from_file

from generate_resonant_chain import bring_planet_pair_in_resonance, semi_to_orbital_period


def resonant_chain_planetary_system(bodies, tau_a_factor, t_integration, n_steps, coll_check=False):
    """
    Create a resonant chain of planetary bodies in a system.
    Args:
        bodies (Particles): The set of bodies in the system.
        tau_a_factor (float): Migration parameter in terms of outer orbital period.
        t_integration (float): Integration time in units of the outer orbital period.
        n_steps (int): Number of integration steps.
        coll_check (bool): Whether to check for collisions or not.
    Returns:
        bodies (Particles): The updated set of bodies in the system.
    """  
    planets = bodies[bodies.type == "planet"]
    for pi in range(len(planets)-1):
        bodies = resonant_pair_planetary_system(
                        bodies,
                        inner_planet_id=pi,
                        tau_a_factor=tau_a_factor,
                        t_integration=t_integration,
                        n_steps=n_steps,
                        plot_results=False,
                        coll_check=coll_check
                        )
    return bodies

def resonant_pair_planetary_system(bodies, inner_planet_id=0, outer_planet_id=1, tau_a_factor=-1e5,
                                   t_integration=100, n_steps=100, plot_results=True, coll_check=False):
    """
    Create a resonant pair of planetary bodies in a system.
    Args:
        bodies (Particles): The set of bodies in the system.
        inner_planet_id (int): Index of the inner planet.
        outer_planet_id (int): Index of the outer planet.
        tau_a_factor (float): Migration parameter in terms of outer orbital period.
        t_integration (float): Integration time in units of the outer orbital period.
        n_steps (int): Number of integration steps.
        plot_results (bool): Whether to plot the results or not.
        coll_check (bool): Whether to check for collisions or not.
    Returns:
        bodies (Particles): The updated set of bodies in the system.
    """

    star = bodies[bodies.type == "star"][0]    
    planets = bodies[bodies.type == "planet"]
    print(f"Number of planets is {len(planets)}")
    
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
    print(f"{P1.name}, {P2.name} Orbital Period Ratio={fraction}, {fraction.numerator/fraction.denominator},  {Porb_a/Porb_b}")
    
    bring_planet_pair_in_resonance(
        planetary_system=bodies, 
        inner_planet=P1, 
        outer_planet=P2,
        tau_a_factor=tau_a_factor,
        t_integration=t_integration, 
        n_steps=n_steps,
        plot_results=plot_results,
        coll_check=coll_check
    )

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

    bodies = read_set_from_file(o.infilename)

    if o.inner_planet_id<0:
        bodies = resonant_chain_planetary_system(bodies,
                                                 o.tau_a_factor,
                                                 o.t_integration, 
                                                 o.n_steps)

    else:
        bodies = resonant_pair_planetary_system(bodies, 
                                                o.inner_planet_id, 
                                                o.outer_planet_id,
                                                o.tau_a_factor,
                                                o.t_integration, 
                                                o.n_steps)
    write_set_to_file(
        bodies,
        o.outfilename,
        overwrite_file=True,
        append_to_file=False,
        close_file=True
    )
