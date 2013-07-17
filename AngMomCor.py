"""
This Program is meant to measure 2 lengths for each star particle in the given
ORION outputs. 

One is Zeta_*, which is the radial distance away from the star
that encloses the same amount of angular momentum as the star. A few nuaces 
must be addressed. 1) Orion does not accurately handle the stars' angular mom.
so it may be totally worthless to compare the stars' angular momenta to the 
realistically handled circumstellar gas' angular momentum. Solution would be to
look in the ORION code and use the details of how it deals with the stars' 
angular momentum to better inform the comparison. 2) A possible solution to 1)
and an interesting question on its own is whether I should use the 
gravitationally bound gas in the comparison or not. I could use the grav bound
gas to represent the ang mom of the star and look further out to the unbound gas
to represent gas that could/will be accreted.

Two is Ksi, which is the shell angular momentum correlation length. I am going
to define it such that L_shell(r)~L_shell(r+Ksi) and if you were to consider
any thicker shells then the angular momenta of the neighboring shells would be 
substantially different from each other. 

I am hoping to show that if Zeta_star >> Ksi then the misalignment angle will be
randomly distributed, where as if Ksi<<Zeta_star then the angle will be small. 
"""
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
def findzeta(zd, CSM_SPHERE_L):
	izd = find_nearest_index(radii,zd)
	LSTAR = CSM_SPHERE_L[izd]
	LSTARMAG = LA.norm(LSTAR)
	find_zeta = np.array([LA.norm(sphere-LSTAR) for sphere in CSM_SPHERE_L])
	izeta=izd+find_nearest_index(find_zeta[izd:], LSTARMAG )
	ZETA = radii[izeta] - radii[izd] 
	return ZETA
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
def find_ksi(percent, CSM_L_SHELL,kd):
	ikd = find_nearest_index(radii,kd)
	KSIs = np.zeros(ikd)
	for i in range(ikd):
		for j in range(i+1,nrad-1):
			difference = CSM_L_SHELL[j]-CSM_L_SHELL[i]
			mag_of_diff = np.sqrt(np.sum((difference)**2))
			mag_of_shell = np.sqrt(np.sum((CSM_L_SHELL[i])**2))
			if mag_of_diff < percent*mag_of_shell:
				continue
			else:
				KSIs[i]=radii[j]-radii[i]
				break
	return KSIs
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
def LENGTHS(pf, position, radii): 
	"""
	This function will calculate ZETA and KSI for the given star particle. 
	In order to do so it must be supplied with the star's ang mom and an
	array containing the radii of the shells that are going to be looked at. 
	"""
	nrad = len(radii)
	CSM_SPHERE_L = np.zeros((nrad,3))
	for i in xrange(nrad):
		sp = pf.h.sphere(position,radii[i]*1.5e13)
		CSM_SPHERE_L[i] = np.array([np.sum(-1.*sp['AngularMomentumX']),np.sum(-1.*sp['AngularMomentumY']),np.sum(-1.*sp['AngularMomentumZ'])])
	CSM_L_SHELL = np.diff(CSM_SPHERE_L, axis=0)
	
	ZETAs=np.zeros(len(range(20,200)))
	for zd in range(20,200):
		ZETAs[zd-20] = findzeta(zd, CSM_SPHERE_L)
		print 'zeta when the disk angular momentum is calculated at ',zd, ' AU is', ZETAs[zd-20]
	
	all_ksis = np.zeros((20, find_nearest_index(radii,200)))
	percents = np.linspace(0.05, 1, num=20)
	for i in xrange(len(percents)):
		percent = percents[i]
		KSIs = find_ksi(percent, CSM_L_SHELL,200)
		all_ksis[i] = KSIs
		#print KSIs
		print 'the average ksi when taken at', percent, 'percent is', np.mean(KSIs), 'with deviation of', np.std(KSIs)
	return CSM_L_SHELL
#==============================================================================#
def ClEAN_STARS(M_min,pf,data):

	"""
	The function ClEAN_STARS takes all the data from an output
	and returns:
	
	nstars 		= 	number of stars
	indices 	= 	array containing the index of each star particle
	masses 		= 	array containing the mass of each star particle
	positions 	= 	array containing the (x,y,z) position of each star 
			  		particle
	L_star 		= 	array containing the unit angular momentum vector of
	         		each star particle
	
	iff the star particles are above the M_min supplied
	"""
	
	indices = data['particle_id']
	masses = data['particle_mass']
	nstars = len(masses)
	positions = np.zeros((nstars, 3))
	L_star = np.zeros((nstars, 3))
	for i in xrange(nstars):
		positions[i] = [data['particle_position_' + direction][i] for direction in ['x', 'y', 'z']]
		L_star_vec = np.array([data['particle_angmomen_' + direction][i] for direction in ['x', 'y', 'z']])
		L_star[i] = L_star_vec/np.sqrt(np.square(L_star_vec).sum())

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

	return clean_nstars, clean_indices, clean_masses, clean_positions, clean_L_star
#==============================================================================#
def DISK_HUNTER(pf, position, L_star, radii):
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
	mass_profile = np.array([])
	angle_profile = np.array([])
	for radius in radii:
		sp = pf.h.sphere(position, radius)
		L_disk  = -1.*sp.quantities['AngularMomentumVector']()
		angle_profile = np.append(angle_profile,np.arccos(np.dot(L_disk, L_star)))
		mass_profile = np.append(mass_profile,sp.quantities["TotalQuantity"]("CellMass"))
	return angle_profile, mass_profile
#==============================================================================#
def MAX_RESOLVER(pf):
	"""
	MAX_RESOLVER finds the maximum spatial resolution of the given parameter files

	It returns the max resolution in AU
	"""
	length = abs(pf.domain_right_edge[0] - pf.domain_left_edge[0])
	coursest_ncells = pf.domain_dimensions[0]
	max_refinement = 0.
	for i in xrange(len(pf.h.grids)):
		if pf.h.grids[i].Level > max_refinement:
			max_refinement = pf.h.grids[i].Level
	Highest_Resolution = length / coursest_ncells / 2.**max_refinement # in cm
	return Highest_Resolution/1.5e13 # in AU
#==============================================================================#
def CAN_OPENER(directory):
	"""
	The function CAN_OPENER finds, then opens the output files in the 
	supplied directory and outputs the parameterfile (pf) and all the data (dd)
	
	This might be better done using yt's built in TimeSeriesData
	"""
	filelist = np.array([])
	for file in os.listdir(directory):
		if fnmatch.fnmatch(file, '*.hdf5'):
			filelist = np.append(filelist, file)
	nfiles = int(len(filelist))
	
	if nfiles == 0:
		print "no files found, make sure they end with .hdf5 \
		and are in" + directory

	pfs = np.array([])
	all_data = np.array([])
	for i in xrange(nfiles):
		pf = load(directory+filelist[i])
		data = pf.h.all_data()
		pfs = np.append(pfs,pf)
		all_data = np.append(all_data,data)
	return pfs, all_data
#==============================================================================#

from yt.mods import *
import glob
import fnmatch
import matplotlib.pyplot as plt
from mpi4py import MPI
from numpy import linalg as LA

day        = 8.64e4                         # seconds
year       = 365.2425 * day                 # seconds
M_sun 	   = 1.9891e33        				# gm

comm 		= MPI.COMM_WORLD
my_rank = comm.Get_rank()
num_procs = comm.size

ts = TimeSeriesData.from_filenames('/clusterfs/henyey/dfielding/andrew/data*.hdf5')


final_pf   = ts[-1]
final_data = final_pf.h.all_data()
final_masses = final_data['particle_mass']
i_max_mass = final_masses.argmax()
# MM_index will be the index of the most massive star
MM_index  = final_data['particle_id'][i_max_mass]

radii = np.logspace(np.log10(5.0), np.log10(350.), num=250)

SHELL_STORAGE = {}
for sto, pf in ts.piter(storage = SHELL_STORAGE):
	data = pf.h.all_data()
	indices = data['particle_id']
	""" getting only the most massive star and if it hasn't formed yet move on """
	index = 0
	for i in range(len(indices)):
		if indices[i] == MM_index:
			index = i
			print 'the out put at time '+ str(pf.current_time/year) + ' has the targeted star'
	if index == 0:
		print 'the out put at time '+ str(pf.current_time/year) + ' does not have the targeted star'
		break

	mass = data['particle_mass'][index]
	position = [data['particle_position_' + direction][index] for direction in ['x', 'y', 'z']]

	CSM_L_SHELL = LENGTHS(pf,position,radii)
	CSM_L_SHELL = np.append(0,CSM_L_SHELL)
	sto.result = (CSM_L_SHELL,pf.current_time/year)

fnlist = np.array(['MM_L_SHELL_HR_'+str(i+1) for i in xrange(len(SHELL_STORAGE))])
for i in range(len(SHELL_STORAGE)):
	np.savetxt(fnlist[i],np.c_[SHELL_STORAGE[i][0],radii], header = 'current time = '+str(SHELL_STORAGE[i][1])+'column 0: L_shell ||| column 1: radii')




























