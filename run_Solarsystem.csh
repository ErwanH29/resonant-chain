#Period ratios observed:
#JS= 5/2
#SU= 3/1
#UN= 2/1

#Solarsystem_giants_R_JS13_SU13_UN12.amuse
python generate_resonant_chain.py -M 1.0 -a 5.204 -e 0.049 -m 317.83 -i 1.3 --name Jupiter --name_star Sun --tau -1e5 -t 10 --nsteps 10
python generate_resonant_chain.py -a 9.573  -e 0.052 -m 95.16 -i 2.5 --name Saturn --tau -1e7 -t 10000 --nsteps 1000
python generate_resonant_chain.py -a 19.165 -e 0.047 -m 14.54 -i 0.8 --name Uranus --tau -1e7 -t 10000 --nsteps 1000
python generate_resonant_chain.py -a 30.178 -e 0.010 -m 17.15 -i 1.8 --name Neptune --tau -1e7 -t 10000 --nsteps 1000
mv input_filename.amuse Solarsystem_giants.amuse

#Solarsystem_giants_R_JS13_SU13_UN12.amuse
#python generate_resonant_chain.py -M 1.0 -a 5.204 -e 0.049 -m 317.83 -i 1.3 --name Jupiter --name_star Sun --tau -1e8 -t 1000 --nsteps 10000
#python generate_resonant_chain.py -a 9.573  -e 0.052 -m 95.16 -i 2.5 --name Saturn --tau -1e8 -t 1000 --nsteps 10000
#python generate_resonant_chain.py -a 19.165 -e 0.047 -m 14.54 -i 0.8 --name Uranus --tau -1e8 -t 1000 --nsteps 10000
#python generate_resonant_chain.py -a 30.178 -e 0.010 -m 17.15 -i 1.8 --name Neptune --tau -1e8 -t 1000 --nsteps 10000
#mv input_filename.amuse Solarsystem_giants.amuse

#Solarsystem_giants_R_JS25_SU13_UN12.amuse
#python generate_resonant_chain.py -M 1.0 -a 5.204 -e 0.049 -m 317.83 -i 1.3 --name Jupiter --name_star Sun --tau -1e8 -t 100 --nsteps 2000
#python generate_resonant_chain.py -P 2.5 -e 0.052 -m 95.16 -i 2.5 --name Saturn --tau -1e8 -t 100 --nsteps 2000
#python generate_resonant_chain.py -P 3.0 -e 0.047 -m 14.54 -i 0.8 --name Uranus --tau -1e8 -t 100 --nsteps 2000
#python generate_resonant_chain.py -P 2.0 -e 0.010 -m 17.15 -i 1.8 --name Neptune --tau -1e8 -t 100 --nsteps 2000
#mv input_filename.amuse Solarsystem_giants.amuse   


