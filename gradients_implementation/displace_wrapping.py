import argparse
import numpy as np
import pandas as pd

from ase.io import read as aread
from cysrc.potential import *
from ase.geometry import get_distances
from math import log
import sys
from ase import *

class UnitCellBounds:
	def __init__(self):
		self.normals = np.zeros((6,3))
		self.points = np.zeros((6,3))
		self.areas = np.zeros((6,))


	def set_face_normals(self, vects):
		"""Sets the normal vectors for each
		of the 6 unit cell faces in the following order:
		ab-plane, bc-plane, ac-plane. Parallel faces
		are stored in adjacent rows of the returned matrix.
		"""
		
		for i in range(0,6,2):
			self.points[i+1,] = vects[(i//2+2)%3,]
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

			# print("vects[{},] x vects[{},]".format(v1,v2))
			# print("(vects[{},]+vects[{},])-vects[{},] x (vects[{},]+vects[{},])-vects[{},])\n".format(
			# 	v1,(v+2)%3,(v+2)%3,v2,(v+2)%3,(v+2)%3))

	 
	def get_intersection(
		self, plane_p, dis_vector, 
		ion_p, faceno, vects, epsilon=1e-6):
		"""Returns the point of intersection between dis_vector and 
		a given face of the unit cell parallilepiped. If the intersection
		point lies outside the face area then it returns None.
		"""
		ndotu = np.dot(self.normals[faceno],dis_vector)
		if abs(ndotu) < epsilon:
			inter_p = ion_p+dis_vector
			si = -1
		else:
			ion_plane_vector = ion_p - plane_p
			si = -np.dot(self.normals[faceno],ion_plane_vector) / ndotu
			inter_p = ion_plane_vector + si * dis_vector + plane_p

		if faceno <= 4 :
			v1 = faceno//2%3
			v2 = (faceno//2+1)%3

			qcrossi = -np.cross(
			inter_p-self.points[faceno,], vects[v1,]) \
			 / self.areas[faceno]
			coef1 = np.linalg.norm(qcrossi)
			icrossr = np.cross(
				inter_p-self.points[faceno,], vects[v2,]) \
				 / self.areas[faceno]
			coef2 = np.linalg.norm(icrossr)
		else:
			v2 = faceno//2%3
			v1 = (faceno//2+1)%3

			qcrossi = np.cross(
				inter_p-self.points[faceno,], vects[v1,]) \
				 / self.areas[faceno]
			coef1 = np.linalg.norm(qcrossi)
			icrossr = -np.cross(
				inter_p-self.points[faceno,], vects[v2,]) \
				 / self.areas[faceno]
			coef2 = np.linalg.norm(icrossr)

		print("coef1",coef1,"coef2",coef2,"I",inter_p,
			"source",self.points[faceno,],
			"normal",self.normals[faceno],
			"normal direction",np.dot(qcrossi,self.normals[faceno]))
		# print(qcrossi)

		if ((np.dot(qcrossi,self.normals[faceno])>=0) & \
			(0 <= si) & (si <= 1) & (0 <= coef1) & \
			(coef1 <= 1) & (0 <= coef2) & (coef2 <= 1)): # check if intersection point
			return inter_p					# is inside face
		else:
			return None

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def cuboid_data(o, size=(1,1,1)):
	# code taken from
	# https://stackoverflow.com/a/35978146/4124317
	# suppose axis direction: x: to left; y: to inside; z: to upper
	# get the length, width, and height
	l, w, h = size
	# x = [[o[0], o[0] + l],  
	# 	[o[0], o[0] + l]]
	x = [[o[0]+l, o[0]+l],  
		[o[0]+w, o[0]+w]]  
	y =	 [[o[1]+w, o[1]+w],          
		 [o[1], o[1]]]   
	z = [[o[2], o[2]+h],                       
		[o[2], o[2]+h]]   
		 # [o[2], o[2]],               
		 # [o[2], o[2]]]               
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
			[5,5,5]])
	points2 = np.array([[3,-2,3],
			[3,3,3]])
	points3 = np.array([[0,0,0],
			[3,3,0]])

	ucb = UnitCellBounds()
	ucb.set_face_normals(atoms.get_cell())

	faceno = 5

	print("line1",ucb.get_intersection(
		ucb.points[faceno], points1[1,]-points1[0,], 
		points1[0], faceno, vects),"\n")
	print("line2",ucb.get_intersection(
		ucb.points[faceno], points2[1,]-points2[0,], 
		points2[0], faceno, vects),"\n")
	print("line3",ucb.get_intersection(
		ucb.points[faceno], points3[1,]-points3[0,], 
		points3[0], faceno, vects),"\n")

	# print(ucb.points[faceno],ucb.normals[faceno])

	fig = plt.figure()
	ax = fig.gca(projection='3d')

	plt.plot(points1[:,0],points1[:,1],points1[:,2],label="line1")
	plt.plot(points2[:,0],points2[:,1],points2[:,2],label="line2")
	plt.plot(points3[:,0],points3[:,1],points3[:,2],label="line3")

	ax.legend()

	plotCubeAt(ax=ax, size=(4,4,4))

	plt.show()