# -*- coding: utf-8 -*-

#Try to read the fits file

import numpy as np
from numpy import unravel_index

#import astropy.units as u
from astropy.io import fits

import matplotlib.cm as cm
import matplotlib.pyplot as plt

import operator

import astropy.units as u

def read_fits_bz(fits_name):
#Read fits file

	opened_fits=fits.open(fits_name)

	data= opened_fits[0].data

#	print np.size(data)
#	fig, ax = plt.subplots()
#	CS = ax.contour(data)
#	plt.show()

	return data[345:475,370:500] #data[370:450,370:500]

def flat_fits(data,nt):
#flatten fits data so only strong flux tubes are present

#	Set a threshold for B_z
#	flat_dat=np.where((np.abs(data) < bmax), 0.0,data)
#	flat_dat=np.where((data > -1250.0), 0.0,data)
	flat_dat=np.where((np.abs(data) < 100.0), 0.0,data)

	nrows=len(flat_dat)
	ncols=len(flat_dat[0])

	print nrows,ncols

	flat_dat_grad=np.gradient(flat_dat)
	flat_dat_grad_i=flat_dat_grad[0]
	flat_dat_grad_j=flat_dat_grad[1]

	tube_loc_i=[]
	tube_loc_j=[]
	tube_loc_B=[]

	for ii in range(2,nrows-1):
		for jj in range(2,ncols-1):
			if (flat_dat_grad_i[ii-1,jj]*flat_dat_grad_i[ii+1,jj] < 0) and (flat_dat_grad_j[ii,jj-1]*flat_dat_grad_j[ii,jj+1] < 0) and (np.abs(flat_dat[ii,jj]) > 100):

				tube_loc_i.append([ii])
				tube_loc_j.append([jj])
				tube_loc_B.append([flat_dat[ii,jj]]) 
#				print ii, jj, flat_dat[ii,jj]
			
#	print tube_loc_i
#	print tube_loc_j

	fti=[]
	ftj=[]
	ftB=[]
	
	for ntube in range(0,nt):
	#	Find the maximum value
#		loc=unravel_index(np.abs(tube_loc_B).argmax(), tube_loc_B.shape)

#		print tube_loc_B
		loc, value = max(enumerate(np.abs(tube_loc_B)), key=operator.itemgetter(1))
#		print loc
		#print tube_loc_i[loc]
		#print tube_loc_j[loc]
			
		fti.append(tube_loc_i[loc])
		ftj.append(tube_loc_j[loc])
		ftB.append(tube_loc_B[loc])	

		tube_loc_B[loc]=[0.0]
	"""#	Remove said flux tube and find next maximum
	#	flat_dat[loc[0]-5:loc[0]+5,loc[1]-5:loc[1]+5]=0.0
		i=loc[0]
		i2=loc[0]
		j=loc[1]
		j2=loc[1]
		while flat_dat[i-1,loc[1]] > flat_dat[i,loc[1]]:
			i=i-1
		while flat_dat[i2+1,loc[1]] > flat_dat[i2,loc[1]]:
			i2=i2+1
		while flat_dat[loc[0],j-1] > flat_dat[loc[0],j]:
			j=j-1
		while flat_dat[loc[0],j2+1] > flat_dat[loc[0],j2]:
			j2=j2+1
		flat_dat[i:i2,j:j2]=0.0
	"""
	
#	loc=unravel_index(np.abs(flat_dat).argmax(), flat_dat.shape)

#	print loc
#	print flat_dat[loc]

#	flat_dat[loc[0]-5:loc[0]+5,loc[1]-5:loc[1]+5]=0.0

#	print np.size(flat_dat)
#	fig, ax = plt.subplots()
#	CS = ax.contour(flat_dat)
#	cbar = fig.colorbar(CS)
#	plt.show()


###############################################################
#	This doesnt work yet....
###############################################################
#	sort out x,y locations scaled between 0 and 1
#	fti *= 0.00769230769  #(nrows-1.0) 
#	ftj *= 0.00769230769  #(ncols-1.0)

#	Scale x,y locations to Mm
#	fti *= 10.0*u.Mm
#	fti -= 5.0*u.Mm
#	ftj *= 10.0*u.Mm
#	ftj -= 5.0*u.Mm
	
#	Flux tubes units to Tesla
#	ftB *= 0.0001*u.T

	return fti,ftj,ftB

