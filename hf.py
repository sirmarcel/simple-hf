#!/usr/bin/env python
#     (the first line allows execution directly from the linux shell) 
#
# --- simple hartree-fock simulation, basically a port of http://personal.ph.surrey.ac.uk/~phs3ps/simple-hf.html
# Author:        M. Langer, T. Horstschaefer
# dependencies:  PYTHON v2.7, numpy
# last modified: 

# coding=utf-8

import math
import numpy as np

def gaussian(x,parameter):
  return (x/parameter)*np.exp(-(x/parameter)**2/2)

hbar = 41.437

def hartreeFock(N, max_iterations, dr, initial_parameter, damping, potential_f):
	grid = np.arange(0, N*dr, dr)
	density = np.zeros(N)

	# init wf with gaussian
	wavefunction = gaussian(grid, initial_parameter)

	for i in range(max_iterations):  
	  wavefunction /= np.sqrt(dr*sum(wavefunction**2))
	  
	  # define density with a new, "virtual" wf (wf / grid) that is nomalised in spherical coordinates
	  density[1:] = 4 * (wavefunction[1:]/grid[1:])**2 / (4*np.pi)
	  density[0] = 4*( (4*wavefunction[1]/3-wavefunction[2]/6) /dr)**2/ (4*np.pi)
	  
	  potential = potential_f(density)

	  # calculate laplacian on the wave function with three-point approach, treating the edges with two-point
	  kinetic = np.empty(N)
	  kinetic[0] = hbar*wavefunction[0]/(dr**2)
	  kinetic[N-1] = -(hbar/(2*dr**2))*(wavefunction[N-2]-wavefunction[N-1])
	  for j in xrange(1,N-1):
	    kinetic[j] = -(hbar/(2*dr**2))*(wavefunction[j-1]-2*wavefunction[j]+ wavefunction[j+1])
	  
	  transformed = kinetic + potential*wavefunction
	  energy = dr * sum(transformed*wavefunction)
	  wavefunction -= damping * (transformed - energy * wavefunction)
	
	return np.array([ grid, wavefunction, density, potential ])
