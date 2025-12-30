#!/usr/bin/env python

# from amuse.community.rebound.interface import Rebound
from amuse.community.ph4.interface import ph4
from amuse.community.huayno.interface import Huayno
from amuse.units import units, constants, nbody_system, quantities
from amuse.lab import Particles, Particle
from amuse.couple import bridge
from amuse.ext.orbital_elements import new_binary_from_orbital_elements, orbital_elements_from_binary

from tqdm import tqdm
import numpy as np
from amuse.io import write_set_to_file, read_set_from_file


# # 1. Initialize a resonance

def semi_to_orbital_period(a, Mtot) :
    return 2*np.pi * (a**3/(constants.G*Mtot)).sqrt()

def orbital_period_to_semi(P, Mtot) :
    return ((constants.G*Mtot) * (P/(2*np.pi))**2)**(1./3.)

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
        required_attributes = ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'tau_a', 'tau_e']
        self.required_attributes = lambda p, x : x in required_attributes
        self.particles = particles
        
    def kick_with_field_code(self, particles, dt):
        #!!!need a for loop for different particles?
        for i in range(1,len(particles)):
            particle = particles[i]
            
            dvx = particle.vx-particles[0].vx
            dvy = particle.vy-particles[0].vy
            dvz = particle.vz-particles[0].vz
            dx  = particle.x-particles[0].x
            dy  = particle.y-particles[0].y
            dz  = particle.z-particles[0].z
            r2  = dx*dx+dy*dy+dz*dz
            
            ax = 0|units.m/units.s**2
            ay = 0|units.m/units.s**2
            az = 0|units.m/units.s**2
            #migration
            if hasattr(particle,'tau_a'):
                ax += dvx/(2.*particle.tau_a)
                ay += dvy/(2.*particle.tau_a)
                az += dvz/(2.*particle.tau_a)

            #ecc-damping
            if hasattr(particle,'tau_e'):
                vdotr  = dx*dvx+dy*dvy+dz*dvz
                prefac = 2*vdotr/r2/particle.tau_e
                ax += prefac*dx
                ay += prefac*dy
                az += prefac*dz
            self.update_velocities(particle, dt, ax, ay, az)
            
    def update_velocities(self, particle, dt,  ax, ay, az):
        particle.vx += dt * ax
        particle.vy += dt * ay
        particle.vz += dt * az
    
    def evolve_model(self, tend):
        timestep = self.timestep
        while self.time < tend:
            dt = min(timestep, tend-self.time)
            codetopart = self.code.particles.new_channel_to(self.particles)
            codetopart.copy()
            parts = self.particles.copy(filter_attributes = self.required_attributes)
            self.kick_with_field_code(parts, dt) # kick
            copytopart = parts.new_channel_to(self.particles)
            copytocode = parts.new_channel_to(self.code.particles)
            copytocode.copy()
            copytopart.copy()
            self.time+=timestep
        self.time = tend

def bring_planet_pair_in_resonance(planetary_system, inner_planet, outer_planet,
                                   tau_a_factor = -1.e5,
                                   t_integration=100, n_steps=100,
                                   plot_results=False):
    
    star = planetary_system[planetary_system.mass.argmax()]
    planets = planetary_system - star
    planets = planets[planets.mass > 0. | units.MSun]
    first_planet = planets[0]
    last_planet = planets[1]
    orbital_elements = orbital_elements_from_binary(star+planets[-1])
    Porbit = semi_to_orbital_period(orbital_elements[2], np.sum(orbital_elements[:2]))

    # set migration timescale
    #planetary_system.tau_a = -np.inf | units.yr
    #planetary_system.tau_e = -np.inf | units.yr
    #outer_planet.tau_a  = tau_a_factor * Porbit
    #outer_planet.tau_e  = outer_planet.tau_a/500
    #planetary_system.tau_a = 0.01*tau_a_factor * Porbit
    #planetary_system.tau_e = 0.01*outer_planet.tau_a/n_steps #500
    planetary_system.tau_a = -np.inf | units.yr
    planetary_system.tau_e = -np.inf | units.yr
    outer_planet.tau_a  = tau_a_factor * Porbit
    outer_planet.tau_e  = outer_planet.tau_a/n_steps #500
    ######
    #outer_planet.tau_a  = -(1e5*t_integration) | units.yr
    #outer_planet.tau_e  = -(2e4*t_integration) | units.yr
    #outer_planet.tau_a  = -1e5 | units.yr
    #outer_planet.tau_e  = -2e4 | units.yr
    #planetary_system[-1].tau_a = -1e6 | units.yr
    #planetary_system.tau_e = -2e5 | units.yr

    converter = nbody_system.nbody_to_si(planetary_system.mass.sum(), Porbit)
    nbody = ph4(convert_nbody=converter)
    #nbody = Huayno(convert_nbody=converter)
    nbody.particles.add_particles(planetary_system)
    #nbody.particles.add_particles(star.as_set())
    #nbody.particles.add_particles(planets)

    #setup nbody
    migration_code = CodeWithMigration(nbody, planetary_system, do_sync=True, verbose=False)
    planet_migration = bridge.Bridge(use_threading=False)
    planet_migration.add_system(nbody)
    planet_migration.add_code(migration_code)
    migration_code.timestep = 0.1 * Porbit

    #orbit_a = orbital_elements_from_binary(planetary_system[0:2], G=constants.G)
    # hostmass, planetmass, semimajor axis, eccentricity, true_anomaly, inclination, loan, aop
    #print(planetary_system)

    N=n_steps
    Porb = []
    sma = []
    ecc = []
    inc = []
    phi_a = []
    phi_b = []
    a1=np.zeros(N)|units.au
    a2=np.zeros(N)|units.au
    e1=np.zeros(N)
    e2=np.zeros(N)
    ele1=[]
    ele2=[]
    ts=np.linspace(1,t_integration,N) * Porbit

    channel_from_system_to_framework = nbody.particles.new_channel_to(planetary_system)
    # nbody.get_time_step()
    for i,t in enumerate(tqdm(ts)):
        planet_migration.evolve_model(t)
        channel_from_system_to_framework.copy()
        orbit_a = orbital_elements_from_binary(star + first_planet)
        a1[i], e1[i]=orbit_a[2], orbit_a[3]
        orbit_b = orbital_elements_from_binary(star + last_planet)
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
        #aop_a = []
        #aop_b = []
        inner_orbit = None
        outer_orbit = None
        #star = bodies[0]
        for pi in planets:
            if outer_orbit is not None:
                inner_orbit = outer_orbit
            outer_orbit = orbital_elements_from_binary(star+pi)
            P.append(semi_to_orbital_period(outer_orbit[2], star.mass))
            a.append(outer_orbit[2])
            e.append(outer_orbit[3])
            i.append(outer_orbit[5]|units.deg)

            # hostmass, planetmass, semimajor axis, eccentricity, true_anomaly, inclination, loan, aop

            if inner_orbit is not None:
                ta_a = inner_orbit[4]
                ta_b = outer_orbit[4]
                aop_a = inner_orbit[7]
                aop_b = outer_orbit[7]
                phi = (ta_a+aop_a)-2*(ta_b+aop_b)
                p_a.append(phi + aop_a)
                p_b.append(phi + aop_b)
            #p_a=[(ele1[i][4]+ele1[i][7])-2*(ele2[i][4]+ele2[i][7])+ele1[i][7] for i in range(N)]
            #p_b=[(ele1[i][4]+ele1[i][7])-2*(ele2[i][4]+ele2[i][7])+ele2[i][7] for i in range(N)]
            
        Porb.append(P.value_in(units.yr))
        sma.append(a.value_in(units.au))
        ecc.append(e)
        inc.append(i.value_in(units.deg))
        phi_a.append(p_a)
        phi_b.append(p_b)

    #print("phi_a:", phi_a)

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

    star = bodies[bodies.type == "star"][0]    
    planets = bodies[bodies.type == "planet"]
    if len(planets) == 0:
        print(f"add planet to star")
        bodies = new_binary_from_orbital_elements(star.mass,
                                                  mplanet, semimajor_axis, eccentricity,
                                                  inclination=inclination)
        bodies[0].type = "star"
        bodies[0].name = name_star
        bodies[1].type = "planet"
        bodies[1].name = name
        #bodies.add_particle(planet)
        #print(bodies)
        orbital_elements = orbital_elements_from_binary(bodies)
        # hostmass, planetmass, semimajor axis, eccentricity, true_anomaly, inclination, loan, aop
        #print(orbital_elements)
        return bodies

    star = bodies[bodies.type == "star"][0]    
    planets = bodies[bodies.type == "planet"]
    last_planet = planets[-1]
    orbital_elements = orbital_elements_from_binary(star+last_planet)
    #print(orbital_elements)
    Porbit = semi_to_orbital_period(orbital_elements[2], np.sum(orbital_elements[:2]))
    print(f"outer orbital period: {Porbit.in_(units.yr)}")

    if Pratio>0:
        Pouter = Pratio*Porbit
        semimajor_axis = orbital_period_to_semi(Pouter, bodies.mass.sum())
    
    #planetary_system = Particles()
    first_planet = last_planet
    second_planet = new_binary_from_orbital_elements(bodies.mass.sum(),
                                                     mplanet, semimajor_axis, 0,
                                                     inclination=inclination)

    bodies.add_particle(second_planet[1])
    bodies[-1].type = "planet"
    bodies[-1].name = name
    bodies.move_to_center()
    #print(bodies)

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

