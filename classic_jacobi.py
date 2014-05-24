from __future__ import division
import numpy as np
from scipy import identity, diag, linalg
from numpy import random
import numpy.testing as npt
import math
'''
The original Jacobi's method for finding eigendecomposition of 
a matrix. 
'''

def classic_jacobi_rotate (A, P, p, q):
	h = A [q,q] - A[p,p]
	if np.abs(A[p,q]) < np.abs(h)*1.0e-36: # approximate t if theta is very big
		t = A[p,q]/np.float64(h)
	else: 
		theta = 0.5*h/A[p,q]
		t = 1.0/(np.abs(theta) + np.sqrt(theta*theta + 1.0)) 
		if theta < 0.0: t = -t
		
	c = 1.0/np.sqrt(t*t+1.0)
	s = c*t
	R = np.array([[c,-s],[s,c]],dtype=np.float64)
	A [[p,q],:] = np.dot(R, A[[p,q],:])
	A [:,[p,q]] = np.dot(A [:,[p,q]],R.T)
	P [:,[p,q]] = np.dot(P [:,[p,q]],R.T)
	#print A
	
def classic_jacobi (A):
	n = A.shape[0]
	P = identity(n,dtype=np.float)
	
	while True:
		p = np.argmax(np.abs(A-diag(diag(A))), axis=0) #max(indices) along columns
		m1 = np.amax(np.abs(A-diag(diag(A))), axis=0) #max(values) along columns
		q = np.argmax(m1)
		p = p[q]

		if np.abs(np.amax(m1)) < 1.0e-15:
			npt.assert_almost_equal(np.dot(P,P.T),identity(n),decimal=10)
			return P, diag(A)	
		classic_jacobi_rotate(A,P,p,q)


'''Some tests'''	
def test1():
	a = np.array([[8,-1,3,-1],[-1,6,2,0],[3,2,9,1],[-1,0,1,7]],dtype=np.float)
	#a = np.array([[1,2],[2,1]], dtype =np.float)
	acopy = a.copy()

	v,a = classic_jacobi (a)
	npt.assert_almost_equal(a[1]*v[:,1],np.dot(acopy,v[:,1]),decimal=10) 
	npt.assert_almost_equal(a[2]*v[:,2],np.dot(acopy,v[:,2]),decimal=10)#sanity check
	print 'eigendecomposition', a
	print v
	
def test2():
	a = np.array([[ 1.55209707,  0.03023631],[ 0.03023631,1.55209707]])
	p = np.identity(2)
	classic_jacobi_rotate(a,p,0,1)
	print a

#test SVD from eigendecomposition
def test3():
	a = random.randn(3,5)
	aTa = np.dot(a.T,a)

	eigvec, eigval = classic_jacobi (aTa)
	S = eigval
	largest = np.argsort(S)[::-1] #inds from largest to smallest
	mindim = min(a.shape)
	S = np.sqrt(S[:,largest[:mindim]])

	V = eigvec
	for i in xrange(eigvec.shape[0]):
		V[i,:] = V[i,:]/linalg.norm(V[i,:])	
	U = np.dot(a,V)
	for i in xrange(U.shape[1]):
		U[:,i] = U[:,i]/linalg.norm(U[:,i])

	V = V[:,largest]#permute columns to match singular values
	U = U [:,largest[:mindim]]
	
	#make the sigma matrix of appropriate dimensions
	sigm = np.zeros(a.shape)
	sigm [:mindim,:mindim] = diag(S)
	
	if not (np.allclose(a,np.dot(U,np.dot(sigm,V.T)))):
		print "WARNING: svd problem"
		return
test3()