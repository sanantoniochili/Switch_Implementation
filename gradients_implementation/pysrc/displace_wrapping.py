
################## NEEDS MORE DEBUGGING ############################

import argparse
import numpy as np
import pandas as pd

from ase.io import read as aread
# from cysrc.potential import *
from ase.geometry import get_distances
from math import log
import sys
from ase import *

class UnitCellWrap:
	def __init__(self):
		self.normals = np.zeros((6,3))
		self.plane_points = np.zeros((6,3))
		self.areas = np.zeros((6,))

	def set_face_normals(self, vects):
		"""Sets the normal vectors for each
		of the 6 unit cell faces in the following order:
		ab-plane, bc-plane, ac-plane. Parallel faces
		are stored in adjacent rows of the returned matrix.
		"""
		
		for i in range(0,6,2):
			self.plane_points[i+1,] = vects[(i//2+2)%3,]
			if i <= 4 :
				v1 = i//2%3
				v2 = (i//2+1)%3
			else:
				v2 = i//2%3
				v1 = (i//2+1)%3

			self.normals[i,] = np.cross(vects[v1,],vects[v2,])
			self.normals[i+1,] = -np.cross(vects[v1,],vects[v2,])
			
			self.areas[i] = np.linalg.norm(self.normals[i,])
			self.areas[i+1] = -np.linalg.norm(self.normals[i+1,])
	 
	def get_intersection(self, move_vector, 
		ion_point, faceno, vects, epsilon=1e-6):
		"""Returns the point of intersection between move_vector and 
		a given face of the unit cell parallilepiped. If the intersection
		point lies outside the face area then it returns None.
		"""
		opp_flag = -1
		if faceno%2 == 0:
			opp_flag = 1

		intersect_point = ion_point

		ndotu = np.dot(self.normals[faceno],move_vector)
		if abs(ndotu) < epsilon:
			si = -1
		else:
			ion_plane_vector = ion_point - self.plane_points[faceno]
			si = -np.dot(self.normals[faceno],ion_plane_vector) / ndotu
			intersect_point = ion_plane_vector + si * move_vector + self.plane_points[faceno]

			ends_test = intersect_point-ion_point
			if (ends_test < epsilon).all():
				return None
			ends_test = intersect_point-ion_point+move_vector
			if (ends_test < epsilon).all():
				return None

		if faceno <= 4 :
			v1 = faceno//2%3
			v2 = (faceno//2+1)%3

			qcrossi = -np.cross(
			intersect_point-self.plane_points[faceno,], vects[v1,]) \
			 / self.areas[faceno]
			coef1 = np.linalg.norm(qcrossi)
			icrossr = np.cross(
				intersect_point-self.plane_points[faceno,], vects[v2,]) \
				 / self.areas[faceno]
			coef2 = np.linalg.norm(icrossr)
		else:
			v2 = faceno//2%3
			v1 = (faceno//2+1)%3

			qcrossi = np.cross(
				intersect_point-self.plane_points[faceno,], vects[v1,]) \
				 / self.areas[faceno]
			coef1 = np.linalg.norm(qcrossi)
			icrossr = -np.cross(
				intersect_point-self.plane_points[faceno,], vects[v2,]) \
				 / self.areas[faceno]
			coef2 = np.linalg.norm(icrossr)

		# print("coef1",coef1,"coef2",coef2,"I",intersect_point,
		# 	"plane point",self.plane_points[faceno,],
		# 	"normal",self.normals[faceno])

		if ((np.dot(qcrossi,self.normals[faceno])>=0) & \
			(0 <= si) & (si <= 1) & (0 <= coef1) & \
			(coef1 <= 1) & (0 <= coef2) & (coef2 <= 1)): # check if intersection point is inside face borders

			return intersect_point
		else:
			return None

	def get_wrap_intersection(self, faceno, vects):
		v3 = (faceno//2+2)%3
		if faceno%2 == 1:
			return -vects[v3]
		else:
			return vects[v3]

	def move_ion(self, ion_point, move_vector, vects):
		while(True):
			# print("\nIon movement {} : {} {}".format(
				# moveno,ion_point,ion_point+move_vector))

			dx = np.array([0.,0.,0.])
			intersect_point = None

			for faceno in range(6):
				# print("\nFACE {}".format(faceno))
				intersect_point_temp = self.get_intersection(
					move_vector, ion_point, faceno, vects)
				# print("Intersection point:",intersect_point_temp)
				if intersect_point_temp is not None:
					if intersect_point is None:
						intersect_point = intersect_point_temp
					dx += self.get_wrap_intersection(faceno, vects)
					# print("Moving ion by",dx)
			
			if intersect_point is None: 
				ion_point = ion_point+move_vector
				break
			
			move_vector = move_vector-(intersect_point-ion_point)
			ion_point = intersect_point+dx
			
			# print("\nWrapped intersection point:",ion_point)
			# print("Remaining displacement:",move_vector)
		return ion_point

	def move_ions(self, ion_positions, displacement_vects, vects, N):
		assert(len(ion_positions)==N)
		assert(ion_positions.shape==displacement_vects.shape)
		new_positions = np.zeros((N,3))

		for ioni in range(N):
			new_positions[ioni] = self.move_ion(ion_positions[ioni], 
									displacement_vects[ioni], vects)
			print("OLD:",ion_positions[ioni],"NEW:",new_positions[ioni])
		return new_positions


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

	ucb = UnitCellWrap()
	ucb.set_face_normals(atoms.get_cell())
	
	move = points3
	ion_point = move[0]
	move_vector = move[1]-move[0]
	
	# ion_point = ucb.move_ion(ion_point, move_vector)
	# print("\nFINAL POINT:",ion_point)

	move = [[5,5,5],
			[1,1,1],
			[0,0,0],
			[0,0,0],
			[0,0,0]]

	new_positions = ucb.move_ions(np.array(atoms.positions), 
		np.array(move), len(atoms.positions))
	print(atoms.positions)
	print(new_positions)


	fig = plt.figure()
	ax = fig.gca(projection='3d')

	plt.plot(points1[:,0],points1[:,1],points1[:,2],label="line1")
	plt.plot(points2[:,0],points2[:,1],points2[:,2],label="line2")
	plt.plot(points3[:,0],points3[:,1],points3[:,2],label="line3")

	ax.legend()

	plotCubeAt(ax=ax, size=(4,4,4), alpha=0.35)

	plt.show()

# first working commit 805d98ce4756973683edd8a2ade9fe2c3ac1dbcc
# wrap point working commit fdddfdbca77f12d7d6d87bb783e052f799eefef6
# wrap move working commit 0b1aed3396eada6394784d596518c3dc72f1e18f
# wrap correct displacement dff817323f747e42e3f7203ef44e5d162fcd5d79
