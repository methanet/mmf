from __future__ import division
from IPython.parallel import Client
from itertools import combinations
import numpy as np, networkx as nx, scipy, time
from scipy import identity, diag, linalg, dot
from classic_jacobi import classic_jacobi_rotate
from numpy import nbytes
from matplotlib import pylab, pyplot as plt
import numpy.testing as npt

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

#returns the index of nxn matrix m in the flattened upper triangle of m
#NB: i<=j, i.e. an upper triangle index
def ind_matrix_to_vec(i, j, n):
	return i*n - (i-1)*i/2 + j-i
	
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
	
	while active.size > 1:
		masterdict = {} # {(i,j): [nd.array of dim(2,),ind_to_remove]}
		
		aTemp = a[:,active][active,:]
		aaT = np.dot(aTemp,aTemp)
		aaTshape = aaT.shape[0]
		
		aaT = aaT[np.triu_indices(aaTshape)] #flatten just the upper triangle of aaT
		active_size = active.size
		
		for k,l in combinations(xrange(active.size),2): 
			i = active[k]; j = active[l]
			
			#this would contain [[cos,sin],[-sin,cos]] after rotation
			m[(0,1),(0,1)]= 1.0 
			m[(0,1),(1,0)]= 0.0
			C[:,:]= 0.0 #contains the 'mass' moved to the diagonal

			C[0,0] = aaT[ind_matrix_to_vec(k,k,aaTshape)]#aaT[k,k]
			C[0,1] = aaT[ind_matrix_to_vec(k,l,aaTshape)]#aaT[k,l]
			C[1,0] = aaT[ind_matrix_to_vec(k,l,aaTshape)]#aaT[k,l]
			C[1,1] = aaT[ind_matrix_to_vec(l,l,aaTshape)]#aaT[l,l]

			#eigval,eigvec = linalg.eigh(C)
			#print eigval
			#print eigvec
			
			classic_jacobi_rotate(C,m,0,1)
			
			print m
			print C
			#if eigval[0]<eigval[1]:
			if C[0,0] <C[1,1]:
				masterdict[(i,j)] = [m[:2,1],i, C[0,0]] 
				#masterdict[(i,j)] = [eigvec[:2,1],i, eigval[0]] #take view of 2nd column instead of (eigvec[0,1],eigvec[1,1])
			else:
				masterdict[(i,j)] = [m[:2,1],i, C[1,1]] 
				#masterdict[(i,j)] = [eigvec[:2,1],j,eigval[1]]
		
		H=nx.Graph()
		
		#TODO use map here or delete each item after use?
		for pair, info in masterdict.iteritems():
			H.add_edge(pair[0],pair[1],weight=-info[2]) #TODO MAKE THIS (i.e. -1) A PARAMETER

		if len(H.edges())==0:
			break
				
		matching = nx.max_weight_matching(H,maxcardinality=True)
		#TODO MAKE THIS MORE EFFICIENT
		#remove repeated vertices from the matching (G is undirected)
		s = set([])
		for key in matching.iterkeys():
			if key < matching[key]:
				s.add((key,matching[key]))
		print "matching: ", s

		to_remove=[]
		for edge in s:
			rotate_active_by(a,p,edge[0],edge[1],masterdict[edge][0],active)
			to_remove.append(masterdict[edge][1])
			
		active01[to_remove] = 0
		active = np.where(active01)[0]
		activelist = list(active)
	npt.assert_almost_equal(np.dot(p,p.T),identity(n),decimal=10) #check orthogonality
	
	print diag(a)
	return diag(a), p

def test1(number_of_nodes):
	G = nx.cycle_graph(number_of_nodes)
	adj_matrix = np.array(nx.adjacency_matrix(G))
	degree = np.identity(number_of_nodes)*2
	laplacian = degree + adj_matrix 
	laplacian= np.array(laplacian,dtype=np.float64)
	
	t = time.time()
	a,wavelets = mmf_approx (linalg.expm3(laplacian*0.1))
	print time.time() - t
	for i in xrange(number_of_nodes):
		plt.plot(np.arange(0, number_of_nodes), wavelets[:,i]/linalg.norm(wavelets[:,i]))
		plt.show()
		
test1(16)

'''
rc = Client()
print rc.ids
dview = rc[:]
dview.execute('import numpy as np')
dview['P'] = np.random.random((3,5))
print dview['P']
#parallel_result = dview.map_sync(lambda x: x*2, pup)
#print parallel_result

def mymap(inds):
	return P[inds[0],inds[1]]

result = dview.map_sync(mymap,[(1,2),(0,2)])
#result = dview.apply(mymap,np.zeros((3,5)))
print result
'''