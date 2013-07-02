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
def L_DISK(pf, position):
	sp = pf.h.sphere(position, radius)
	L_disk = -1.*sp.quantities['AngularMomentumVector']()
	L_disk = L_disk/np.sqrt(np.square(L_disk).sum())
	return L_disk
#==============================================================================#

from yt.mods import *
import glob
import fnmatch
import matplotlib.pyplot as plt
from mpi4py import MPI
day        = 8.64e4                         # seconds
year       = 365.2425 * day                 # seconds
M_sun 	   = 1.9891e33        				# gm

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
num_procs = comm.size


"""
Here I am using the built in time series function to open all the 
parameter files and load then into an array: ts
"""
t0=time.time()
ts = TimeSeriesData.from_filenames("../data/stella/*.orion")
nfiles = len(ts)

"""
Find the number of stars in each pf, then get their:
index, mass, position, L_star, L_disk(@50AU), time
"""
radius = 1.5e13 * 50. #the radius for disk angular momentum calculation
my_storage = {}
for sto, pf in ts.piter(storage = my_storage):
	data = pf.h.all_data()
	nstars, indices, masses, positions, L_star = STAR_GATHERER(pf,data)
	nstars, indices, masses, positions, L_star = STAR_CLEANER(0.0001, 0., nstars, indices, masses, positions, L_star)
	L_disk = np.zeros((nstars,3))
	for i in xrange(int(nstars)):
		L_disk[i] = L_DISK(pf,positions[i])
	sto.result = (nstars, indices, masses, positions, L_star, L_disk, pf.current_time)

"""
Create an array of unique indicies and find how many stars there are
"""
max_nstar = 0
uniq_indices = np.array([])
for i in range(nfiles):
	uniq_indices = np.append(uniq_indices, my_storage[i][1])
	if my_storage[i][0] > max_nstar:
		max_nstar = my_storage[i][0]
uniq_indices = np.unique(uniq_indices)


"""
repackage everything so that each dictionary element contains
all that we need to know about 1 star only:
how many times we see the star, star ages, star angular momenta,
disk angular momenta, star masses, index, & the star positions
"""
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
	L_disk_hist	  = np.zeros((ntimes, 3))
	ages = np.zeros(ntimes)
	for l in xrange(ntimes):
		mass_hist[l] 			= my_storage[matches1[l]][2][matches2[l]]
		position_hist[l] 		= my_storage[matches1[l]][3][matches2[l]]
		L_star_hist[l] 			= my_storage[matches1[l]][4][matches2[l]]
		L_disk_hist[l] 			= my_storage[matches1[l]][5][matches2[l]]
		ages[l] 				= my_storage[matches1[l]][6]
	stars[i] = (ntimes, ages, L_star_hist, L_disk_hist, mass_hist, index, position_hist)

"""
For each star the angle between the initial angular momentum and the subsequent
angular momenta is calculated, and so is the angle between the disk and the 
initial angular momentum of the star. 
"""
for i in xrange(len(uniq_indices)):
	ntimes 		= stars[i][0]
	ages        = stars[i][1]
	L_star_hist = stars[i][2]
	L_disk_hist = stars[i][3]
	mass_hist	= stars[i][4]/M_sun
	star_angles = np.zeros(ntimes)
	disk_angles = np.zeros(ntimes)
	for j in xrange(int(ntimes)):
		star_angles[j] = np.arccos(np.dot(L_star_hist[0],L_star_hist[j]))
		disk_angles[j] = np.arccos(np.dot(L_star_hist[0],L_disk_hist[j]))
	star_angles[0] = 0.0 # to avoid the NaN at time 0
	filename = 'star_'+str(stars[i][5])+'_single_angles.txt'
	np.savetxt(filename,np.c_[ages, star_angles, disk_angles, mass_hist],header = 'ages, star angle, disk angle, masses')

"""
For each star set up a coordinate system where the star's initial ang. mom.
is the z axis and then calculate the *theta* *phi* angles of the star's & disk's
angular momentum at each time.
"""
for i in xrange(len(uniq_indices)):
	ntimes 		= stars[i][0]
	ages        = stars[i][1]
	L_star_hist = stars[i][2]
	L_disk_hist = stars[i][3]
	mass_hist	= stars[i][4]/M_sun
	star_angles = np.zeros((ntimes,2)) # the first angle is theta, the second is phi
	disk_angles = np.zeros((ntimes,2))

	zhat = L_star_hist[0]
	xhat = np.cross(zhat, np.array([0.,0.,1.]) )
	yhat = np.cross(zhat, xhat)

	for j in xrange(int(ntimes)):
		theta 	= np.arccos( np.dot(L_star_hist[j],zhat))
		phi 	= np.arccos( np.dot(L_star_hist[j],xhat)/np.sin(theta))
		star_angles[j,0] = (180./np.pi)*theta
		star_angles[j,1] = (180./np.pi)*phi
	star_angles[0,0]=0.0 # to avoid the NaN at time 0
	star_angles[0,1]=0.0 # to avoid the NaN at time 0

	for k in xrange(int(ntimes)):
		theta 	= np.arccos( np.dot(L_disk_hist[k],zhat))
		phi 	= np.arccos( np.dot(L_disk_hist[k],xhat)/np.sin(theta))
		disk_angles[k,0] = (180./np.pi)*theta
		disk_angles[k,1] = (180./np.pi)*phi

	filename = 'star_'+str(stars[i][5])+'_double_angles.txt'
	np.savetxt(filename,np.c_[ages, star_angles, disk_angles, mass_hist], header = 'ages, star theta, star phi, disk theta, disk phi, masses')
















