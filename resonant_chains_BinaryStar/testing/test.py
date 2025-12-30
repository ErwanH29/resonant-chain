from amuse.lab import *
import glob


original = read_set_from_file("initial_conds.hdf5", "hdf5")

others = Particles()
others.add_particles(original[original.syst_id == -1])

shaken_systems = glob.glob("systems/*")
for s in shaken_systems:
    p = read_set_from_file(s)
    print(len(p), p)
    STOP
    others.add_particles(p)
    
print(len(original), len(others))