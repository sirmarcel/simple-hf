# Must be run with 'ipython --gui=wx --pylab=wx' and then 'run plot.py'
import hf

from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from mayavi import mlab

import numpy as np
import os


"""
# Run fortran code
os.system("gfortran hf.F90 -o hf && ./hf")

# Read results from fortran code
data = np.genfromtxt("fort.15")
grid = data[:,0]
density = data[:,1]

data = np.genfromtxt("fort.16")
potential = data[:,1]
plt.plot(grid, density)
"""
a = -408.75
b =  1093.0

param = {
    'N': 100,
    'max_iterations': 1000,
    'dr' : 0.1,              # grid spacing
    'initial_parameter': 1.4,
    'damping': 0.0001,

    # HF Potential parameters
    'potential_f': lambda density: 2*a*density + 3*b*density**2
}

# Solve Hartree-Fock equation
grid, wavefunction, density, potential = hf.hartreeFock(**param)

# MATPLOTLIB

# show 2d plots
plt.figure(1)

plt.subplot(311)
plt.title('wavefunction')
plt.ylabel(r'$\psi$')
plt.plot(grid, wavefunction)

plt.subplot(312)
plt.title('density')
plt.ylabel(r'$\rho$')
plt.plot(grid, density)

plt.subplot(313)
plt.title('potential')
plt.ylabel(r'$V$')
plt.xlabel(r'$r$')
plt.plot(grid, potential)

plt.tight_layout()
plt.show()


# MAYAVI

# Make a 2D grid out of a radial 1D function
def Make2DGrid(x, y):
	f = interp1d(x, y, fill_value=0.)
	g_max = max(x)
	n = len(x)

	# evaluate result on 2D grid
	s = np.zeros((n, n))
	for x in xrange(0, n):
		for y in xrange(0, n):
			s[x,y] = f( np.sqrt( (x - n/2)**2 + (y - n/2)**2 ) * param['dr'] )
	return s



# Make grid to shift the objects later
g_max = max(grid)
x, y = np.mgrid[ -g_max:g_max:param['dr'], -g_max:g_max:param['dr'] ]

# Plot objects shifted on x-axis
mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0))
shift = g_max+1


scale = lambda y: 5./max(abs(y))

mlab.surf(x - shift,y, Make2DGrid(grid, potential),    warp_scale = scale(potential))
mlab.surf(x,y,         Make2DGrid(grid, density),      warp_scale = scale(density))
mlab.surf(x + shift,y, Make2DGrid(grid, wavefunction), warp_scale = scale(wavefunction))

# Add description
t_z = 8
ts_scale = 0.5
mlab.text3d(-shift,-g_max, t_z, 'potential')
mlab.text3d(-shift,-g_max, t_z-1, 'scale={:.1f}'.format(scale(potential)), scale=ts_scale)

mlab.text3d(     0,-g_max, t_z, 'density')
mlab.text3d(     0,-g_max, t_z-1, 'scale={:.1f}'.format(scale(density)), scale=ts_scale)

mlab.text3d( shift,-g_max, t_z, 'wavefunction')
mlab.text3d( shift,-g_max, t_z-1, 'scale={:.1f}'.format(scale(wavefunction)), scale=ts_scale)
