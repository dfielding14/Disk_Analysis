"""
ANALYZER

This program is meant to read ORION I & II data outputs and 
analyze them in the context of protostellar disk misalignment.

The pertinent information in each data output are the star 
particles' indices, masses, postions, and angular momenta 
directions, as well as the angular momentum and structure of
the gas within ~100-150 AU of each star particle.

Created by Drummond Fielding 
on June 1, 2013 
in Berkeley, California
while eating a peach
"""
#==============================================================================#
#==============================================================================#
"""
The function STAR_GATHERER takes all the data from an output
and returns:

nstars 		= 	number of stars
indices 	= 	array containing the index of each star particle
masses 		= 	array containing the mass of each star particle
positions 	= 	array containing the (x,y,z) position of each star 
		  		particle
L_star 		= 	array containing the unit angular momentum vector of
         		each star particle
"""
def STAR_GATHERER(pf,data):
	indices = data['particle_id']
	masses = data['particle_mass']
	nstars = len(masses)
	positions = np.zeros((nstars, 3))
	L_star = np.zeros((nstars, 3))
	for i in xrange(nstars):
		positions[i] = [data['particle_position_' + direction][i] for direction in ['x', 'y', 'z']]
		L_star_vec = np.array([data['particle_angmomen_' + direction][i] for direction in ['x', 'y', 'z']])
		L_star[i] = L_star_vec/np.sqrt(np.square(L_star_vec).sum())
	return nstars, indices, masses, positions, L_star

#==============================================================================#
"""
The function STAR_CLEANER takes the results of STAR_GATHERER
and removes the stars that do not match the given criterion.

M_min = minimum mass in units of M_sun
D_min = minumum distance between stars that have mass above
		M_min but maybe be too close to be distinguished from
		a nearby companion star. Returns only the most massive
		of the two. D_min should be in AU

TO DO: Distance cutoff
"""
def STAR_CLEANER(M_min, D_min, nstars, indices, masses, positions, L_star):
	# Mass Cleaning
	clean_nstars	= 0
	clean_indices	= np.array([])
	clean_masses	= np.array([])
	for i in xrange(nstars):
		if masses[i]/M_sun >= M_min:
			clean_nstars +=1
			clean_indices 	= np.append(clean_indices, indices[i])
			clean_masses 	= np.append(clean_masses, masses[i])
	clean_positions	= np.zeros((clean_nstars,3))
	clean_L_star	= np.zeros((clean_nstars,3))
	for i in xrange(clean_nstars):
		for j in xrange(nstars):
			if clean_indices[i] == indices[j]:			
				clean_positions[i] 	= positions[j]
				clean_L_star[i] 	= L_star[j]
	# # Distance Cleaning, tc stands for totally clean
	# tc_nstars 	= 0
	# tc_indices 	= np.array([])
	# tc_masses	= np.array([])
	# for i in xrange(clean_nstars):
	# 	for j in xrange(clean_nstars):
	# 		if i != j and np.sqrt(np.sum( (clean_positions[i]-clean_positions[j])**2 )) <= D_min*1.5e13:
	# 			if clean_masses[i] > clean_masses[j]:
	# 				tc_nstars +=1
	# 				tc_indices = np.append()
	# 		else:
	# 			break

	# tc_positions 	= np.zeros((tc_nstars,3))
	# tc_L_star 		= np.zeros((tc_nstars,3))

	return clean_nstars, clean_indices, clean_masses, clean_positions, \
		clean_L_star
#==============================================================================#
"""
The function DISK_HUNTER looks at the circumstellar material
and determines the misalignment angle as a function of radius, 
as well as the mass profile as a function of R and PHI (the 
polar azimuthal angle). The two mass profiles allow drastic 
changes in misalignment to be distinguished from over/under
densities near the star that may be transient inflowing or
outflowing material.

TO DO: Angular mass profile
for now the angular mass profile will be left out because
defining cylindrical-wedge shaped objects is non trivial.
"""
def DISK_HUNTER(pf, position, L_star, radii):
	mass_profile = np.array([])
	angle_profile = np.array([])
	for radius in radii:
		sp = pf.h.sphere(position, radius)
		L_disk  = -1.*sp.quantities['AngularMomentumVector']()
		angle_profile = np.append(angle_profile,np.arccos(np.dot(L_disk, L_star)))
		mass_profile = np.append(mass_profile,sp.quantities["TotalQuantity"]("CellMass"))
	return angle_profile, mass_profile
#==============================================================================#
"""
The function FILE_COUNTER returns the number of orion outputs are in the 
directory supplied.
"""
def FILE_COUNTER(directory):
	filelist = np.array([])
	for file in os.listdir(directory):
		if fnmatch.fnmatch(file, '*.orion') or fnmatch.fnmatch(file, '*.hdf5'):
			filelist = np.append(filelist, file)
	nfiles = int(len(filelist))
	
	if nfiles == 0:
		print "no files found, make sure they end with .orion or .hdf5 \
		and are in the directory given"
	return nfiles
#==============================================================================#

from yt.pmods import *
from physical_constants import *
import glob
import fnmatch

# This is for use on the NERSC machines
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
# import pylab as plt
#from matplotlib import rc
#rc('text', usetex=True)
#rc('font', family='serif')

#########################################################################
num_procs = 2 # make sure to change this when using different computers #
#########################################################################

nradii = 20
min_radii = 5.
max_radii = 100.
radii = np.logspace(np.log10(min_radii*1.5e13), np.log10(max_radii*1.5e13),nradii)

