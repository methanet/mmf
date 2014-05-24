from __future__ import division
import numpy as np
from numpy import random as random,diag
from scipy import linalg as linalg

import matplotlib.pyplot as plt

from mmfcomplete import mmfcomplete
def svdcomplete(A): 
	m,n = A.shape
	U,s,V = linalg.svd(A,full_matrices=False)
	S = np.diag(s)
	
	err = np.empty(11)
	for i in xrange(11):
 		App = np.dot(U[:,:i], np.dot(S[:i,:][:,:i], V[:i,:]))
		print np.allclose(A,App)
		err[i] = linalg.norm(A-App)
	print "er"
	print err
	plt.plot(err)
	plt.show()
	
A = random.rand(10,20)
svdcomplete(A)
missingx = [1,2,3]
missingy = [5,6,2]
A [missingx,missingy] = 0.0

observed = np.vstack(np.nonzero(A)).T
mmfcomplete(A,observed)