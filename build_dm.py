from sklearn.neighbors import kneighbors_graph
import numpy as np
from scipy.sparse import csgraph
from sklearn import preprocessing
import pandas as pd

# takes data as nsamples x nfeature
def diffusion_net(in_csv, num_nbhrs, out_prefix):
    d = np.genfromtxt(in_csv, delimiter=',') # getting the data
    
    mms = preprocessing.MinMaxScaler() # normalize the data
    dn = mms.fit_transform(d) # normalized
    
    G = kneighbors_graph(dn, 
                         n_neighbors = num_nbhrs, 
                         mode = 'connectivity', 
                         metric = 'minkowski', 
                         p = 2, 
                         metric_params = None, 
                         include_self = False) # create graph from the Euclidean distances 
    G = G.toarray()
    L = csgraph.laplacian(G, normed=True) # Getting the normalized graph laplacian matrix 
    eigval, eigvec = np.linalg.eig(L) # getting the eigenvalues and eigenvectors
    #eigval_abs = abs(eigval) # getting the absolute values of eigenvals (these are complex numbers)
    eigval_abs = eigval.real
    
    
    minvals = np.argsort(eigval_abs)
    mineig1, mineig2, mineig3 = 125,100,111 # initialize the min eigenvalue indices to nonzero eigenvals 
    for minval in minvals:
        if eigval_abs[minval] == 0.0:
            continue
        if eigval_abs[mineig1] > eigval_abs[minval]:
            mineig1 = minval
    for minval in minvals:
        if eigval_abs[minval] == 0.0:
            continue
        if eigval_abs[mineig2] > eigval_abs[minval] and eigval_abs[mineig1] != eigval_abs[minval]:
            mineig2 = minval
    for minval in minvals:
        if eigval_abs[minval] == 0.0:
            continue
        if eigval_abs[mineig3] > eigval_abs[minval] and eigval_abs[mineig1] != eigval_abs[minval] and and eigval_abs[mineig2] != eigval_abs[minval]:
            mineig2 = minval
        
    # the code below finds the two smallest non zero eigenvalues: 
    minvals = np.argsort(eigval_abs)
    mineig2 = 100 # randomly initialize the minimum eigenvalue indices to index values which are known to correspond to nonzero eigenvals 
    mineig1 = np.min(np.nonzero(eigval_abs))
    mineig3 = 125
    for minval in minvals:
        if eigval_abs[minval] == 0.0:
            continue
        if eigval_abs[mineig2] > eigval_abs[minval] and eigval_abs[mineig1] != eigval_abs[minval]:
            mineig2 = minval
        if eigval_abs[mineig3] > eigval_abs[minval] and eigval_abs[mineig2] != eigval_abs[minval]:
            mineig3 = minval
    
    DC1 = eigvec[:, mineig1] # based off of the eigenvalues found, get the corresponding eigenvectors
    DC2 = eigvec[:, mineig2]
    DC3 = eigvec[:, mineig3]
    
    f = open(in_csv)
    labels = f.readline().split(',')[1:] # read in the labels # read in the labels
    df = pd.DataFrame({'DC1':np.real(DC1), 'DC2':np.real(dim2), 'DC3':np.real(DC3),})
    df.index = labels
    df.to_csv(out_prefix + str(num_nbhrs) + ".csv")