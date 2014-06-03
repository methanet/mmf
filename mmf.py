from __future__ import division
from itertools import combinations
import numpy as np, networkx as nx, scipy
from scipy import identity, diag, linalg, dot
from classic_jacobi import classic_jacobi_rotate

from matplotlib import pylab, pyplot as plt
import numpy.testing as npt

'''Performs orthogonal rotation of A and P by given 
rotation angle tuple phi=(sin,cos) and a list of active indices.
'''

'''Performs orthogonal rotation of A and P by given 
rotation angle tuple phi=(sin,cos) and a list of active indices.
'''
#phi is (2,) array
def rotate_active_by (A, P, p, q, phi, active):
	c = phi[1]#1.0/np.sqrt(t*t+1.0)
	s = phi[0]#c*t
	R = np.array([[c,-s],[s,c]],dtype=np.float64)
	
	rows = np.dot(R, A[[p,q],:][:,active])
	A [p,active] = rows[0,:]
	A [q,active] = rows[1,:]
	
	cols = np.dot(A [:,[p,q]][active,:],R.T)
	A[active,p] = cols[:,0]
	A[active,q] = cols[:,1]
	
	P [:,[p,q]] = np.dot(P [:,[p,q]],R.T)

#def calc_parameters(aaT,active,C,m,):	
#return cost, phis, temp remove
	
'''The binary parallel MMF algorithm '''
def mmf_approx (a):
	n = a.shape[0]
	p = identity(n,dtype=np.float64) #corresponds to the the wavelets	

	active01 = np.ones(n)
	active = np.where(active01)[0]

	activelist = range(n)
	m = identity(2,dtype=np.float64)
	C = np.empty((2,2),dtype=np.float64)
	cost = np.empty((n,n),dtype=np.float64)
	
	while active.size > 1:
		cost[:,:] = 0.0
		aTemp = a[:,active][active,:]
		aaT = np.dot(aTemp,aTemp)
		print a.nbytes,aaT.nbytes
		active_size = active.size
		masterdict = {} # {(i,j): [nd.array of dim(2,),ind_to_remove]}
		
		#for k in xrange(active.size-1):
		#	for l in xrange(k+1,active.size):
		for k,l in combinations(xrange(active.size),2): 
			i = active[k]; j = active[l]
			
			#this would contain [[cos,sin],[-sin,cos]] after rotation
			m[(0,1),(0,1)]= 1.0 
			m[(0,1),(1,0)]= 0.0
			C[:,:]= 0.0 #contains the 'mass' moved to the diagonal

			C[0,0] = aaT[k,k]
			C[0,1] = aaT[k,l]
			C[1,0] = aaT[k,l]
			C[1,1] = aaT[l,l]
			
			eigval,eigvec = linalg.eigh(C)
			classic_jacobi_rotate(C,m,0,1)
			
			if eigval[0]<eigval[1]: 
				cost[i,j] = eigval[0]
				masterdict[(i,j)] = [eigvec[:2,1],i] #second column view instead of (eigvec[0,1],eigvec[1,1])
			else:
				cost[i,j] = eigval[1]
				masterdict[(i,j)] = [eigvec[:2,1],j]
				
		cost[:,:] = cost[:,:] + cost.T
		G = nx.from_numpy_matrix(cost*(-1)) #since the diffusion matrix is used
	
		if len(G.edges())==0:
			break
			
		matching = nx.max_weight_matching(G,maxcardinality=True)
		#remove repeated vertices from the matching (G is undirected)
		s = set([])
		for key in matching.iterkeys():
			if key < matching[key]:
				s.add((key,matching[key]))
		#print "matching: ", s

		to_remove=[]
		for edge in s:
			rotate_active_by(a,p,edge[0],edge[1],masterdict[edge][0],active)
			to_remove.append(masterdict[edge][1])
			
		active01[to_remove] = 0
		active = np.where(active01)[0]
		activelist = list(active)
	npt.assert_almost_equal(np.dot(p,p.T),identity(n),decimal=10) #check orthogonality
	return diag(a), p
