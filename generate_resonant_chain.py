#!/usr/bin/env python

from amuse.community.ph4.interface import ph4
from amuse.ext.basicgraph import UnionFind
from amuse.units import units, constants, nbody_system, quantities
from amuse.lab import Particles, Particle
from amuse.couple import bridge
from amuse.ext.orbital_elements import new_binary_from_orbital_elements, orbital_elements_from_binary

from tqdm import tqdm
import numpy as np
from amuse.io import write_set_to_file, read_set_from_file

import matplotlib.pyplot as plt


STAR_MASS = 0.08 | units.MSun

def semi_to_orbital_period(a, Mtot) :
    """
    Calculate the orbital period from the semi-major axis and total mass.
    Args:
        a (units.length): Semi-major axis in units of length.
        Mtot (units.mass): Total mass of the system in units of mass.
    Returns:
        units.time: Orbital period in units of time.
    """
    return 2*np.pi * (a**3/(constants.G*Mtot)).sqrt()

def orbital_period_to_semi(P, Mtot) :
    """
    Calculate the semi-major axis from the orbital period and total mass.
    Args:
        P (units.time): Orbital period in units of time.
        Mtot (units.mass): Total mass of the system in units of mass.
    Returns:
        units.length: Semi-major axis in units of length.
    """
    return ((constants.G*Mtot) * (P/(2*np.pi))**2)**(1./3.)

def resolve_merger(coll_set):
    """
    Resolve a merger using sticky-sphere approximation
    Args:
        coll_set (UnionFind): Colliding particle pairs.
    Returns:
        Particle: A new particle representing the remnant of the merger.
    """
    remnant = Particle()
    remnant.mass = coll_set.mass.sum()
    remnant.position = coll_set.center_of_mass()
    remnant.velocity = coll_set.center_of_mass_velocity()
    if remnant.mass < STAR_MASS:
        remnant.type = "planet"
    else:
        remnant.type = "star"
    
    return remnant

def ZAMS_radius(star_mass) -> units.RSun:
    """
    Define stellar radius at ZAMS.
    Args:
        star_mass (units.mass): Mass of the star
    Returns:
        units.length: The ZAMS radius of the star
    """
    mass_in_sun = star_mass.value_in(units.MSun)
    mass_sq = (mass_in_sun)**2.

    numerator = mass_in_sun**1.25 * (0.1148 + 0.8604 * mass_sq)
    denominator = (0.04651 + mass_sq)
    r_zams = numerator / denominator

    return r_zams | units.RSun

def planet_radius(planet_mass) -> units.REarth:
    """
    Compute planet radius (arXiv:2311.12593).
    Args:
        planet_mass (units.mass):  Mass of the planet
    Returns:
        units.length:  Planet radius
    """
    mass_in_earth = planet_mass.value_in(units.MEarth)

    if mass_in_earth < 7.8:
        return (1. | units.REarth)*(mass_in_earth)**0.41
    elif mass_in_earth < 125:
        return (0.55 | units.REarth)*(mass_in_earth)**0.65
    return (14.3 | units.REarth)*(mass_in_earth)**(-0.02) 
    

class CodeWithMigration():
    # not elegant: 1. need the additional input of particles. 2. need additional asignment of timestep
    def __init__(self, code, particles, do_sync=True, verbose=False):
        self.code = code
        if hasattr(self.code, 'model_time'):
            self.time = self.code.model_time
        else:
            self.time = quantities.zero
        self.do_sync=do_sync
        self.verbose=verbose
        self.timestep=None
        required_attributes = [
            'mass', 'x', 'y', 'z', 
            'vx', 'vy', 'vz', 
            'tau_a', 'tau_e'
            ]
        self.required_attributes = lambda p, x : x in required_attributes
        self.particles = particles
        
    def kick_with_field_code(self, particles, dt):
        """
        Apply the kick to the particles using the field code.
        
        Args:
            particles (Particles): The particles to update.
            dt (units.time): The time step for the kick.
        """
        star = particles[particles.mass.argmax()]
        if len(particles.mass > STAR_MASS) > 1:
            raise Exception("More than one star in the planetary system, cannot continue.")
        
        planets = particles - star
        
        dvx = planets.vx - star.vx
        dvy = planets.vy - star.vy
        dvz = planets.vz - star.vz
        dx  = planets.x  - star.x
        dy  = planets.y  - star.y
        dz  = planets.z  - star.z
        r2  = dx*dx + dy*dy + dz*dz
        
        ax  = dvx * (0 | units.s)**-1
        ay  = dvy * (0 | units.s)**-1
        az  = dvz * (0 | units.s)**-1
        
        # Migration
        if hasattr(planets, 'tau_a'):
            ax += dvx/(2. * planets.tau_a)
            ay += dvy/(2. * planets.tau_a)
            az += dvz/(2. * planets.tau_a)
        
        # Eccentricity damping
        if hasattr(planets, 'tau_e'):
            vdotr  = dx*dvx + dy*dvy + dz*dvz
            prefac = (2. * vdotr/r2) / planets.tau_e
            ax += prefac * dx
            ay += prefac * dy
            az += prefac * dz
                
        self.update_velocities(planets, dt, ax, ay, az)
            
    def update_velocities(self, particle, dt,  ax, ay, az):
        """
        Update the velocities of the particles.
        Args:
            particle (Particles): The particles to update.
            dt (units.time): The time step for the update.
            ax (units.length/units.time**2): Acceleration in x direction.
            ay (units.length/units.time**2): Acceleration in y direction.
            az (units.length/units.time**2): Acceleration in z direction.
        """
        particle.vx += dt * ax
        particle.vy += dt * ay
        particle.vz += dt * az
    
    def evolve_model(self, tend):
        """
        Evolve the model to a specified end time.
        Args:
            tend (units.time): The end time to evolve the model to.
        """
        ### Require new channels every step?
        timestep = self.timestep
        while self.time < tend:
            dt = min(timestep, tend-self.time)
            codetopart = self.code.particles.new_channel_to(self.particles)
            codetopart.copy()
            parts = self.particles.copy(filter_attributes = self.required_attributes)
            self.kick_with_field_code(parts, dt)
            copytopart = parts.new_channel_to(
                            self.particles, 
                            attributes=["vx", "vy", "vz"], 
                            target_names=["vx", "vy", "vz"]
                            )
            copytocode = parts.new_channel_to(
                            self.code.particles,
                            attributes=["vx", "vy", "vz"], 
                            target_names=["vx", "vy", "vz"]
                            )
            copytocode.copy()
            copytopart.copy()
            self.time+=timestep
        self.time = tend

def bring_planet_pair_in_resonance(planetary_system, 
                                   inner_planet, 
                                   outer_planet,
                                   tau_a_factor = -1.e5,
                                   t_integration=100, 
                                   n_steps=100,
                                   plot_results=False,
                                   coll_check=False):
    
    star = planetary_system[planetary_system.mass >= STAR_MASS][0]    
    planets = planetary_system[planetary_system.mass < STAR_MASS]
    first_planet = planets[0]
    last_planet = planets[1]
    orbital_elements = orbital_elements_from_binary(star+planets[-1])
    Porbit = semi_to_orbital_period(orbital_elements[2], np.sum(orbital_elements[:2]))

    # set migration timescale
    planetary_system.tau_a = -np.inf | units.yr
    planetary_system.tau_e = -np.inf | units.yr
    outer_planet.tau_a = tau_a_factor * Porbit
    outer_planet.tau_e = outer_planet.tau_a/n_steps

    converter = nbody_system.nbody_to_si(planetary_system.mass.sum(), Porbit)
    nbody = ph4(convert_nbody=converter)
    nbody.particles.add_particles(planetary_system)
    if coll_check:
        particles = nbody.particles
        star_mask = particles.mass > STAR_MASS
        
        particles[star_mask].radius = ZAMS_radius(particles[star_mask].mass)
        particles[~star_mask].radius = planet_radius(particles[~star_mask].mass)
        
        nbody.parameters.epsilon_squared = (0.01*Porbit)**2
        coll_sc = nbody.stopping_conditions.collision_detection
        coll_sc.enable()

    #setup nbody
    migration_code = CodeWithMigration(nbody, planetary_system, do_sync=True, verbose=False)
    planet_migration = bridge.Bridge(use_threading=False)
    planet_migration.add_system(nbody)
    planet_migration.add_code(migration_code)
    migration_code.timestep = 0.1 * Porbit

    # Look to convert [] --> np.zeros((n_planets, N)) to optimise performance
    N = n_steps
    Porb = []
    sma = []
    ecc = []
    inc = []
    phi_a = []
    phi_b = []
    a1 = np.zeros(N)|units.au
    a2 = np.zeros(N)|units.au
    e1 = np.zeros(N)
    e2 = np.zeros(N)
    ele1 = []
    ele2 = []
    ts = np.linspace(1, t_integration, N) * Porbit

    channel_from_system_to_framework = nbody.particles.new_channel_to(planetary_system)
    for i,t in enumerate(tqdm(ts)):
        planet_migration.evolve_model(t)
        if coll_check:
            if coll_sc.is_set():
                print("!!! Collision Detected !!!")
                coll_set = UnionFind()
                for p, q in zip(coll_sc.particles(0), coll_sc.particles(1)):
                    coll_set.union(p, q)
                coll_set = coll_set.sets()
                
                remnant = resolve_merger(coll_set)
                nbody.particles.remove_particles(coll_sc.particles(0))
                nbody.particles.remove_particles(coll_sc.particles(1))
                nbody.particles.add_particle(remnant)
                nbody.particles.synchronize_to(planetary_system)
                
        channel_from_system_to_framework.copy()
        
        orbit_a = orbital_elements_from_binary(star + first_planet)
        orbit_b = orbital_elements_from_binary(star + last_planet)
        
        a1[i], e1[i]=orbit_a[2], orbit_a[3]
        a2[i], e2[i]=orbit_b[2], orbit_b[3]
        ele1.append(orbit_a)
        ele2.append(orbit_b)
        
        P = [] | units.yr
        a = [] | units.au
        e = []
        i = [] | units.deg
        p_a = [] 
        p_b = []
        phi = []
        inner_orbit = None
        outer_orbit = None
        for pi in planets:
            if outer_orbit is not None:
                inner_orbit = outer_orbit
            outer_orbit = orbital_elements_from_binary(star+pi)
            P.append(semi_to_orbital_period(outer_orbit[2], star.mass))
            a.append(outer_orbit[2])
            e.append(outer_orbit[3])
            i.append(outer_orbit[5]|units.deg)

            if inner_orbit is not None:
                ta_a = inner_orbit[4]
                ta_b = outer_orbit[4]
                aop_a = inner_orbit[7]
                aop_b = outer_orbit[7]
                phi = (ta_a+aop_a)-2*(ta_b+aop_b)
                p_a.append(phi + aop_a)
                p_b.append(phi + aop_b)
            
        Porb.append(P.value_in(units.yr))
        sma.append(a.value_in(units.au))
        ecc.append(e)
        inc.append(i.value_in(units.deg))
        phi_a.append(p_a)
        phi_b.append(p_b)

    if plot_results:
        phi_a = np.array(phi_a).T
        phi_b = np.array(phi_b).T
        
        Porb = np.array(Porb).T
        sma = np.array(sma).T
        ecc = np.array(ecc).T
        inc = np.array(inc).T
    
        for ai in sma[:]:
            plt.plot(ts.value_in(units.yr), ai, lw=3)
        #plt.axhline(y=20.8, linestyle='-', lw=1)
        plt.xlabel('Time[yr]')
        plt.ylabel('a [au]')
        plt.show()

        for Pi in Porb[:]:
            plt.plot(ts.value_in(units.yr), Pi, lw=3)
        #plt.axhline(y=20.8, linestyle='-', lw=1)
        plt.xlabel('Time[yr]')
        plt.ylabel('P [yr]')
        plt.show()
    
        for ei in ecc[:]:
            plt.plot(ts.value_in(units.yr), ei)
        plt.xlabel('Time[yr]')
        plt.ylabel('e')
        plt.show()

        for ii in inc[:]:
            plt.plot(ts.value_in(units.yr), ii)
        plt.xlabel('Time[yr]')
        plt.ylabel('i [deg]')
        plt.show()

        for pi in phi_a[:]:
            plt.scatter(ts.value_in(units.yr), pi%(360))
        for pi in phi_b[:]:
            plt.scatter(ts.value_in(units.yr), pi%(360))
        plt.xlabel('Time[yr]')
        plt.ylabel('phi [deg]')
        plt.show()

        for pi in range(len(phi_a[:])):
            plt.scatter(phi_a[pi]%(360), phi_b[pi]%(360))
        plt.xlabel('phi a [deg]')
        plt.ylabel('phi b [deg]')
        plt.show()

    return planetary_system
    
def add_planet_in_resonant_chain(bodies, name_star, semimajor_axis,
                                 eccentricity,
                                 inclination, mplanet, Pratio=0, name="Aa",
                                 tau_a_factor=-1e5,
                                 t_integration=100,
                                 n_steps=100):

    star = bodies[bodies.mass >= STAR_MASS][0]    
    planets = bodies[bodies.mass < STAR_MASS]
    if len(planets) == 0:
        print(f"add planet to star")
        bodies = new_binary_from_orbital_elements(star.mass,
                                                  mplanet, semimajor_axis, eccentricity,
                                                  inclination=inclination)
        bodies[0].type = "star"
        bodies[0].name = name_star
        bodies[1].type = "planet"
        bodies[1].name = name
        orbital_elements = orbital_elements_from_binary(bodies)
        return bodies

    star = bodies[bodies.mass >= STAR_MASS][0]    
    planets = bodies[bodies.mass < STAR_MASS]
    last_planet = planets[-1]
    orbital_elements = orbital_elements_from_binary(star+last_planet)
    Porbit = semi_to_orbital_period(orbital_elements[2], np.sum(orbital_elements[:2]))
    print(f"outer orbital period: {Porbit.in_(units.yr)}")

    if Pratio>0:
        Pouter = Pratio*Porbit
        semimajor_axis = orbital_period_to_semi(Pouter, bodies.mass.sum())
    
    second_planet = new_binary_from_orbital_elements(bodies.mass.sum(),
                                                     mplanet, semimajor_axis, 0,
                                                     inclination=inclination)

    bodies.add_particle(second_planet[1])
    bodies[-1].type = "planet"
    bodies[-1].name = name
    bodies.move_to_center()

    bodies = bring_planet_pair_in_resonance(bodies, bodies[-2], bodies[-1],
                                            tau_a_factor=tau_a_factor,
                                            t_integration=t_integration,
                                            n_steps=n_steps)

    return bodies

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n_planets", type="int",
                      default = 2,
                      help="number of planets [%default]")
    result.add_option("-f", dest="infilename", 
                      default = "input_filename.amuse",
                      help="input infilename [%default]")
    result.add_option("-a", dest="semimajor_axis", type="float", unit=units.au,
                      default = 1|units.au,
                      help="semi-major axis of the first planet [%default]")
    result.add_option("-e", dest="eccentricity", type="float", 
                      default = 0.0,
                      help="eccenticity of the first planet [%default]")
    result.add_option("-i", dest="inclination", type="float", unit=units.deg,
                      default = 1|units.deg,
                      help="semi-major axis of the first planet [%default]")
    result.add_option("-P", "--Period_ratio",
                      dest="Pratio", type="float", 
                      default = 0,
                      help="Resonant period ratio [%default]")
    result.add_option("--tau", 
                      dest="tau_a_factor", type="float", 
                      default = -1e5,
                      help="migration parameter (in terms of outer orbital period) [%default]")
    result.add_option("-t", "--t_integration", 
                      dest="t_integration", type="float", 
                      default = 100,
                      help="migration time scale (in terms of outer orbital period) [%default]")
    result.add_option("--nsteps", 
                      dest="n_steps", type="int", 
                      default = 100,
                      help="number of migration steps [%default]")
    result.add_option("-m", dest="mplanet", type="float", unit=units.MEarth,
                      default = 1|units.MEarth,
                      help="mass of the planet planet [%default]")
    result.add_option("--name", dest="name", 
                      default = "planet",
                      help="planet name [%default]")
    result.add_option("--name_star", dest="name_star", 
                      default = "Sun",
                      help="stellar name [%default]")
    result.add_option("-M", dest="Mstar", type="float", unit=units.MSun,
                      default = -1|units.MSun,
                      help="mass of the star [%default]")
    return result
    
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    if o.Mstar<=0|units.MSun:
        bodies = read_set_from_file(o.infilename)
    else:
        bodies = Particles(1)
        bodies.type = "star"
        bodies.mass = o.Mstar
        bodies.name = o.name_star

    bodies = add_planet_in_resonant_chain(bodies,
                                          o.name_star,
                                          o.semimajor_axis,
                                          o.eccentricity,
                                          o.inclination, o.mplanet, o.Pratio, o.name,
                                          o.tau_a_factor,
                                          o.t_integration,
                                          o.n_steps)
    write_set_to_file(bodies,
                      o.infilename,
                      overwrite_file=True,
                      append_to_file=False,
                      close_file=True)
