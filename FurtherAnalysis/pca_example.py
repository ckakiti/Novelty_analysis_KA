import numpy as np
from scipy import linalg as LA
#mock data #imagine rows are syllables and columns are mice
#change values and/or size of array to see how it effects the eigenvectors/eigenvalues 
mockdata = np.array([
        [5.1,5, 13],
        [5,9,28],
        [4,6,2],
        [5,8,16],])
#centering the data
mockdata -= np.mean(x, axis = 0)  
#calculate covariance matrix
cov = np.cov(mockdata, rowvar = False)
#get eigenvectors and eigenvalues
evals , evecs = LA.eigh(cov)
#formatting and SORT by eigenvalues (aka sort by how much variance each eigen vector explains)
idx = np.argsort(evals)[::-1]
evecs = evecs[:,idx]
evals = evals[idx]
#this is your data represented as PCS 
mockdata_represented_as_pcs = np.dot(x, evecs)