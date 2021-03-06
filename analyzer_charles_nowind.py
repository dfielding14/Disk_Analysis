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
The function CAN_OPENER finds, then opens the output files in the 
supplied directory and outputs the parameterfile (pf) and all the data (dd)

This might be better done using yt's built in TimeSeriesData
"""
def CAN_OPENER(directory):
	filelist = np.array([])
	for file in os.listdir(loc):
		if fnmatch.fnmatch(file, '*.orion') or fnmatch.fnmatch(file, '*.hdf5'):
			filelist = np.append(filelist, file)
	nfiles = int(len(filelist))
	
	if nfiles == 0:
		print "no files found, make sure they end with .orion or .hdf5 \
		and are in the directory given"

	pfs = np.array([])
	all_data = np.array([])
	for i in xrange(nfiles):
		pf = load(directory+filelist)
		data = pf.h.all_data()
		pfs = np.append(pfs,pf)
		all_data = np.append(all_data,data)
	return pfs, all_data
#==============================================================================#

from yt.pmods import *
import glob
import fnmatch
import matplotlib.pyplot as plt
day        = 8.64e4                         # seconds
year       = 365.2425 * day                 # seconds
M_sun 	   = 1.9891e33        				# gm

##########################################################################
#num_procs = 3 # make sure to change this when using different computers # ~*~*~TEST~*~*~
num_procs = 12 # make sure to change this when using different computers #
##########################################################################

t0=time.time()
ts = TimeSeriesData.from_filenames("/clusterfs/henyey/dfielding/charles/charles/nowind/pltm*") #charles w/ wind
#ts = TimeSeriesData.from_filenames("/clusterfs/henyey/dfielding/charles/charles/wind_test/pltm*") #charles w/ wind test
nfiles = len(ts)

nradii = 48
min_radii = 40.
max_radii = 400.
radii = np.logspace(np.log10(min_radii*1.5e13), np.log10(max_radii*1.5e13),nradii)

my_rank = ytcfg.getint("yt", "__topcomm_parallel_rank")

if my_rank == 0:
	print 'the number of files in the time series is ' + str(nfiles)


my_storage = {}
for sto, pf in ts.piter(storage = my_storage):
	print 'working on', pf.parameter_filename, 'which is at time:', pf.current_time/year
	data = pf.h.all_data()
	nstars, indices, masses, positions, L_star = STAR_GATHERER(pf,data)
	nstars, indices, masses, positions, L_star = STAR_CLEANER(0.05, 0., nstars, indices, masses, positions, L_star)
	angle_profiles = np.zeros((nstars,nradii))
	mass_profiles  = np.zeros((nstars,nradii))
	for i in xrange(int(nstars)):
		print "I am processor "+str( my_rank )+" and I am working on " + str(i+1)+' out of '+str(nstars)
		angle_profiles[i], mass_profiles[i] = DISK_HUNTER(pf,positions[i], L_star[i], radii)
	print 'processor '+str(my_rank)+' is done'
	sto.result = nstars, indices, masses, positions, L_star, angle_profiles, mass_profiles, pf.current_time
t1=time.time()


if my_rank == 0:
	print 'the time it took to gather and clean all stars, and hunt their disks was:',t1-t0, 'seconds'

max_nstar = 0
uniq_indices = np.array([])
for i in range(nfiles):
	uniq_indices = np.append(uniq_indices, my_storage[i][1])
	if my_storage[i][0] > max_nstar:
		max_nstar = my_storage[i][0]
uniq_indices = np.unique(uniq_indices)

if my_rank == 0:
	print 'most number of stars in an output: ', max_nstar
	print 'the unique indices are:', uniq_indices
t2 = time.time()

stars = {}	
for i in xrange(len(uniq_indices)):
	ntimes=0
	matches1 = np.array([])
	matches2 = np.array([])
	for j in xrange(nfiles):
		for k in xrange(int(my_storage[j][0])):
			if uniq_indices[i] == my_storage[j][1][k]:
				ntimes += 1
				matches1 = np.append(matches1, j)
				matches2 = np.append(matches2, k)
	index = uniq_indices[i]
	mass_hist     = np.zeros(ntimes)
	position_hist = np.zeros((ntimes, 3))
	L_star_hist   = np.zeros((ntimes, 3))
	angle_profile_hist = np.zeros((ntimes, nradii))
	mass_profile_hist  = np.zeros((ntimes, nradii))
	age = np.zeros(ntimes)
	for k in xrange(ntimes):
		mass_hist[k] 			= my_storage[matches1[k]][2][matches2[k]]
		position_hist[k] 		= my_storage[matches1[k]][3][matches2[k]]
		L_star_hist[k] 			= my_storage[matches1[k]][4][matches2[k]]
		angle_profile_hist[k] 	= my_storage[matches1[k]][5][matches2[k]]
		mass_profile_hist[k]  	= my_storage[matches1[k]][6][matches2[k]]
		age[k] 					= my_storage[matches1[k]][7]
	stars[i] = (ntimes, index, mass_hist, position_hist, L_star_hist,angle_profile_hist,mass_profile_hist,age)
t3 = time.time()

if my_rank == 1:
	print 'the time it took arrange the stars in their dictionary was:',t3-t2, 'seconds'

for i in range(0+my_rank,len(stars),num_procs):
	ntimes = int(stars[i][0])
	star_masses=np.zeros(nradii)
	star_ages = np.zeros(nradii)
	misalignment_angle_profiles = np.zeros((ntimes, nradii))
	mass_profiles = np.zeros((ntimes, nradii))
	for j in range(ntimes):
		star_masses[j] = stars[i][2][j]
		misalignment_angle_profiles[j] = stars[i][5][j]
		mass_profiles[j] = stars[i][6][j]
		star_ages[j] = stars[i][7][j] / year
	filename='nowind_star_'+str(stars[i][1])+'_misalignment_mass_profile_farout.txt'
	np.savetxt(filename,np.c_[star_ages, star_masses, np.transpose(misalignment_angle_profiles),np.transpose(mass_profiles),radii/1.5e13], header = 'disk star misalignment analysis of myers data. column 0: age(years), column 1: masses(g), column 2-2+ntimes: misalignment profiles, column 3+ntimes - 3+2ntimes: mass profiles, column 4+2ntimes: radii(AU)') 
