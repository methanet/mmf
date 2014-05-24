from __future__ import division
import numpy as np
from numpy import random as random, diag
from scipy import linalg as linalg

from mmf_approx import mmf_approx

'''Use MMF for matrix completion
Input: Matrix with observed entries, 2xm matrix of observed entries
Output: Recovered matrix
'''
def mmfcomplete(A, observed):
	#from classic_jacobi import classic_jacobi
	m,n = A.shape
	b = A[observed[:,0],observed[:,1]] #flattened observed entries in A
	AtA = np.dot(A.T,A)

	
	eigval, eigvec = mmf(AtA) #destroys A
	
	S = eigval #similary to the singular value
	largest = np.argsort(S)[::-1]
	mindim = min(A.shape)
	S = np.sqrt(S[:,largest[:mindim]])
	
	V = eigvec
	for i in xrange(eigvec.shape[1]):
		V[i,:] = V[i,:]/linalg.norm(V[i,:])
	U = np.dot(A,V)
	for i in xrange(U.shape[1]):
		U[:,i] = U[:,i]/linalg.norm(U[:,i])

	V = V[:,largest]
	U = U [:,largest[:mindim]]
	print U.shape, S.shape, V.shape
	
	P = np.empty((observed.shape[0],mindim))
	for i in xrange(mindim):
		P[:,i] = (np.outer(U[:,observed[i,0]],V.T[:,observed[i,1]]))[observed[:,0],observed[:,1]].flat
		
	#least square problem Px=b
	x,residuals, rank, s = linalg.lstsq(P,b)

	#recover the matrix 
	sigmas = np.zeros((m,n))
	sigmas[xrange(m),xrange(m)]= x
	A_complete = np.dot(np.dot(U,sigmas),V.T)
	return A_complete