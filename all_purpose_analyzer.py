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

TO DO: Distance cutoff
"""
def STAR_CLEANER(M_min, nstars, indices, masses, positions, L_star):
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

UPDATE 7/12/13:
no longer am I going to use the whole sphere within a radius
I am now going to use the shell ie the difference between
the two spheres
"""
def DISK_HUNTER(pf, position, L_star, radii):
	nradii = len(radii)

	L_sphere = np.zeros((nradii,3))
	M_sphere = np.zeros(nradii)
	for i in xrange(nradii):
		sp = pf.h.sphere(position,radii[i])
		L_sphere[i] = np.array([np.sum(-1.*sp['AngularMomentumX']),np.sum(-1.*sp['AngularMomentumY']),np.sum(-1.*sp['AngularMomentumZ'])])
		M_sphere[i] = sp.quantities["TotalQuantity"]("CellMass")[0]


	L_shell_vec = np.diff(L_sphere, axis=0)
	L_shell = np.zeros((nradii,3))
	for i in xrange(nradii-1):
		L_shell[i+1] = L_shell_vec[i]
	L_shell_unit_vec = np.array([shell_vec/LA.norm(shell_vec) for shell_vec in L_shell_vec])
	M_shell = np.diff(M_sphere)

	angle_profile = np.array([])
	for i in xrange(nradii-1):
		angle_profile = np.append(angle_profile,np.arccos(np.dot(L_shell_unit_vec[i], L_star)))

	angle_profile = np.append(0., angle_profile)
	M_shell = np.append(0., M_shell)
	return angle_profile, L_shell, L_sphere, M_shell, M_sphere
#==============================================================================#
"""
MAX_RESOLVER finds the maximum spatial resolution of the given parameter files

It returns the max resolution in AU
"""
def MAX_RESOLVER(pf):
	length = abs(pf.domain_right_edge[0] - pf.domain_left_edge[0])
	coursest_ncells = pf.domain_dimensions[0]
	max_refinement = 0.
	for i in xrange(len(pf.h.grids)):
		if pf.h.grids[i].Level > max_refinement:
			max_refinement = pf.h.grids[i].Level
	Highest_Resolution = length / coursest_ncells / 2.**max_refinement # in cm
	return Highest_Resolution/1.5e13 # in AU
#==============================================================================#

from yt.pmods import *
import glob
import fnmatch
import matplotlib.pyplot as plt
from mpi4py import MPI
from numpy import linalg as LA

day        = 8.64e4                         # seconds
year       = 365.2425 * day                 # seconds
M_sun 	   = 1.9891e33        				# gm

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
print my_rank
num_procs = comm.size

def radiusizer(ts, min_rad, max_rad, nrad):	
	nfiles = len(ts)
	for i in xrange(nfiles):
		max_res = MAX_RESOLVER(ts[i])
		if min_rad + 1.0 < max_res:
			print 'the minimum radius you supplied was too small and was increased from ' + str(min_rad) + 'AU to '+ str(max_res+1.0) + 'AU, which is 1 AU more than the highest res.'
			min_rad = max_res+1.0
	radii = np.logspace(np.log10(min_rad*1.5e13), np.log10(max_rad*1.5e13),nrad)
	return radii, nrad

def work_horse(ts, radii, nradii, fn_prefix):
	nfiles = len(ts)
	"""
	Making the array of radii to be used in the coming calculations. There is a built in check so that 
	the minumum radius is not smaller than the actual highest resolution of the data outputs
	"""
	t0=time.time()

	if my_rank == 0:
		print 'the number of files in the time series is ' + str(nfiles)
	
	my_storage = {}
	for sto, pf in ts.piter(storage = my_storage):
		print 'working on', pf.parameter_filename, 'which is at time:', pf.current_time/year
		data = pf.h.all_data()
		nstars, indices, masses, positions, L_star = STAR_GATHERER(pf,data)
		nstars, indices, masses, positions, L_star = STAR_CLEANER(0.01, nstars, indices, masses, positions, L_star)
		angle_profiles = np.zeros((nstars,nradii))
		L_shell_vecs = np.zeros((nstars, nradii,3))
		L_spheres= np.zeros((nstars,nradii,3))
		M_shells = np.zeros((nstars,nradii))
		M_spheres= np.zeros((nstars,nradii))
		for i in xrange(int(nstars)):
			print "I am processor "+str( my_rank )+" and I am working on " + str(i+1)+' out of '+str(nstars)
			angle_profiles[i], L_shell_vecs[i], L_spheres[i], M_shells[i], M_spheres[i]= DISK_HUNTER(pf,positions[i], L_star[i], radii)
		# timing
		time_processor_finished = time.time()
		time_processor = time_processor_finished - t0
		#storage of results
		sto.result = (nstars, indices, masses, positions, L_star, angle_profiles, L_shell_vecs, L_spheres, M_shells, M_spheres, pf.current_time)
		print 'processor '+str(my_rank)+' is done and it took ' + str(time_processor) + ' seconds = '+ str(time_processor/60.) + ' minutes'
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
	
	for i in range(0+my_rank, len(uniq_indices), num_procs):
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
		angle_profile_hist 	= np.zeros((ntimes, nradii))
		L_shell_vec_hist	= np.zeros((ntimes,nradii,3))
		L_sphere_hist		= np.zeros((ntimes,nradii,3))
		M_shell_hist 		= np.zeros((ntimes,nradii))
		M_sphere_hist		= np.zeros((ntimes,nradii))
		age = np.zeros(ntimes)
		
		for k in xrange(ntimes):
			mass_hist[k] 			= my_storage[matches1[k]][2][matches2[k]]
			position_hist[k] 		= my_storage[matches1[k]][3][matches2[k]]
			L_star_hist[k] 			= my_storage[matches1[k]][4][matches2[k]]
			angle_profile_hist[k] 	= my_storage[matches1[k]][5][matches2[k]]
			L_shell_vec_hist[k]  	= my_storage[matches1[k]][6][matches2[k]]
			L_sphere_hist[k]  		= my_storage[matches1[k]][7][matches2[k]]
			M_shell_hist[k]  		= my_storage[matches1[k]][8][matches2[k]]
			M_sphere_hist[k]  		= my_storage[matches1[k]][9][matches2[k]]
			age[k] 					= my_storage[matches1[k]][10]
	
	
		radii_time = np.zeros((ntimes,nradii))
		for i in xrange(ntimes):
			radii_time[i] = radii
		filename= fn_prefix + str(index) + '_shells.txt'
		
		L_shell_for_saving = np.zeros((ntimes, 3*nradii))
		L_sphere_for_saving = np.zeros((ntimes, 3*nradii))
		for j in xrange(ntimes):
			for i in xrange(nradii):
				L_shell_for_saving[j,3*i] = L_shell_vec_hist[j,i,0]
				L_shell_for_saving[j,3*i+1] = L_shell_vec_hist[j,i,1]
				L_shell_for_saving[j,3*i+2] = L_shell_vec_hist[j,i,2]
				L_sphere_for_saving[j,3*i] = L_sphere_hist[j,i,0]
				L_sphere_for_saving[j,3*i+1] = L_sphere_hist[j,i,1]
				L_sphere_for_saving[j,3*i+2] = L_sphere_hist[j,i,2]
		np.savetxt(filename,np.c_[age, mass_hist, L_star_hist, angle_profile_hist, L_shell_for_saving, L_sphere_for_saving, M_shell_hist, M_sphere_hist, radii_time/1.5e13], header = 'disk star misalignment analysis of myers data. column 0: age(years) | column 1: masses(g) | column 2,3,4:L_star x,y,z |  angle profile hist |  L_shell_vec |  L_sphere_hist | M_shell_hist | M_sphere_hist | radii') 
	t3 = time.time()

# if my_rank < 3:
# 	print "hi from ", my_rank
# 	ts = TimeSeriesData.from_filenames("/clusterfs/henyey/dfielding/stella/pltrt2704*") #stella 1
# 	radii , nradii = radiusizer(ts, 5., 300., 32)
# 	work_horse(ts, radii, nradii, 'stella_04_star_')
# elif my_rank < 6 and my_rank > 2:
# 	print "hi from ", my_rank
# 	ts = TimeSeriesData.from_filenames("/clusterfs/henyey/dfielding/stella/pltrt2705*") #stella 2
# 	radii , nradii = radiusizer(ts, 5., 300., 32)
# 	work_horse(ts, radii, nradii, 'stella_05_star_')
# elif my_rank < 9 and my_rank > 5:
# 	print "hi from ", my_rank
# 	ts = TimeSeriesData.from_filenames("/clusterfs/henyey/dfielding/stella/pltrt2708*") #stella 3
# 	radii , nradii = radiusizer(ts, 5., 300., 32)
# 	work_horse(ts, radii, nradii, 'stella_08_star_')
# elif my_rank < 12 and my_rank > 8:
# 	print "hi from ", my_rank
# 	ts = TimeSeriesData.from_filenames("/clusterfs/henyey/dfielding/stella/pltrt2713*") #stella 4
# 	radii , nradii = radiusizer(ts, 5., 300., 32)
# 	work_horse(ts, radii, nradii, 'stella_13_star_')


# ts = TimeSeriesData.from_filenames("/clusterfs/henyey/dfielding/stella/pltrt2704*") #stella 1
# radii , nradii = radiusizer(ts, 5., 300., 16)
# work_horse(ts, radii, nradii, 'stella_04_star_')

# ts = TimeSeriesData.from_filenames("/clusterfs/henyey/dfielding/stella/pltrt2705*") #stella 2
# radii , nradii = radiusizer(ts, 5., 300., 16)
# work_horse(ts, radii, nradii, 'stella_05_star_')

# ts = TimeSeriesData.from_filenames("/clusterfs/henyey/dfielding/stella/pltrt2708*") #stella 3
# radii , nradii = radiusizer(ts, 5., 300., 16)
# work_horse(ts, radii, nradii, 'stella_08_star_')

# ts = TimeSeriesData.from_filenames("/clusterfs/henyey/dfielding/stella/pltrt2713*") #stella 4
# radii , nradii = radiusizer(ts, 5., 300., 16)
# work_horse(ts, radii, nradii, 'stella_13_star_')


ts = TimeSeriesData.from_filenames("/clusterfs/henyey/dfielding/charles/charles/wind/pltm*") # charles wind
radii , nradii = radiusizer(ts, 40., 500., 20)
work_horse(ts, radii, nradii, 'charles_wind_star')

ts = TimeSeriesData.from_filenames("/clusterfs/henyey/dfielding/charles/charles/nowind/pltm*") # charles wind
radii , nradii = radiusizer(ts, 40., 500., 20)
work_horse(ts, radii, nradii, 'charles_no_wind_star')


