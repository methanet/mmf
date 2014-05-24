import numpy as np 

cost = np.ones((4,4))*8
a = np.array([[8,-1,3,-1],[-1,6,2,0],[3,2,9,1],[-1,0,1,7]],dtype=np.float)

def fun(i,c):
	return (a[i,i],a[i,i]*10)
	
result = map (fun, cost,range(3))
print result

	