import argparse
import numpy as np
import pandas as pd

from ase.io import read as aread
from ase.geometry import get_distances
from math import log
import sys
from ase import *

class UnitCellBounds:
	def __init__(self):
		self.normals = np.zeros((6,3))
		self.plane_points = np.zeros((6,3))
		self.areas = np.zeros((6,))

	

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def cuboid_data(o, size=(1,1,1)):
	# code taken from
	# https://stackoverflow.com/a/35978146/4124317
	# suppose axis direction: x: to left; y: to inside; z: to upper
	# get the length, width, and height
	l, w, h = size
	x = [[o[0], o[0] + l, o[0] + l, o[0], o[0]],  
		 [o[0], o[0] + l, o[0] + l, o[0], o[0]],  
		 [o[0], o[0] + l, o[0] + l, o[0], o[0]],  
		 [o[0], o[0] + l, o[0] + l, o[0], o[0]]]  
	y = [[o[1], o[1], o[1] + w, o[1] + w, o[1]],  
		 [o[1], o[1], o[1] + w, o[1] + w, o[1]],  
		 [o[1], o[1], o[1], o[1], o[1]],          
		 [o[1] + w, o[1] + w, o[1] + w, o[1] + w, o[1] + w]]   
	z = [[o[2], o[2], o[2], o[2], o[2]],                       
		 [o[2] + h, o[2] + h, o[2] + h, o[2] + h, o[2] + h],   
		 [o[2], o[2], o[2] + h, o[2] + h, o[2]],               
		 [o[2], o[2], o[2] + h, o[2] + h, o[2]]]                  
	return np.array(x), np.array(y), np.array(z)

def plotCubeAt(pos=(0,0,0), size=(1,1,1), ax=None,**kwargs):
	# Plotting a cube element at position pos
	if ax !=None:
		X, Y, Z = cuboid_data( pos, size )
		ax.plot_surface(X, Y, Z, rstride=1, cstride=1, **kwargs)


if __name__=="__main__":
	atoms = Atoms("SrTiO3",

				  cell=[[4.00, 0.00, 0.00],
						[0.00, 4.00, 0.00],
						[0.00, 0.00, 4.00]],

				  positions=[[0, 0, 0],
							 [2, 2, 2],
							 [0, 2, 2],
							 [2, 0, 2],
							 # [1.5, .5, 2], # put Os too close
							 [2, 2, 0]],
				  pbc=True)

	vects = atoms.get_cell()
	points1 = np.array([[0,0,0],
			[9,9,10]])
	points2 = np.array([[3,0,3],
			[4,9,4]])
	points3 = np.array([[1,1,1],
			[-1,-1,3]])

	ucb = UnitCellBounds()
	ucb.set_face_normals(atoms.get_cell())
	
	move = points3
	ion_point = move[0]
	move_vector = move[1]-move[0]

	# move = [[5,5,5],
	# 		[1,1,1],
	# 		[0,0,0],
	# 		[0,0,0],
	# 		[0,0,0]]

	# new_positions = ucb.move_ions(np.array(atoms.positions), 
	# 	np.array(move), len(atoms.positions))
	# print(atoms.positions)
	# print(new_positions)


	fig = plt.figure()
	ax = fig.gca(projection='3d')

	plt.plot(points1[:,0],points1[:,1],points1[:,2],label="line1")
	plt.plot(points2[:,0],points2[:,1],points2[:,2],label="line2")
	plt.plot(points3[:,0],points3[:,1],points3[:,2],label="line3")

	ax.legend()

	plotCubeAt(ax=ax, size=(4,4,4), alpha=0.35)

	plt.show()