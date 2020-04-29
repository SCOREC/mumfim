import numpy as np
from scipy.sparse.linalg import eigs, eigsh

def check_validity(spmat):
    print "Number of nonzeros: " + str(spmat.nnz)
    symmetric = np.abs((spmat-spmat.transpose()) >= 1E-15).nnz == 0
    # note this is not the best method for determining if it is symmetric
    # due to the fact that we are using an absolute tolerance...
    print "Symmetric: " + str(symmetric)
    print "Diagonal Positive: " + str(np.all(spmat.diagonal() > -1E-15))
    D = np.abs(spmat.diagonal())
    S = np.sum(np.abs(spmat), axis=1) - D
    diagonal_dominant = np.all(D > S)
    print "Diagonal Dominant: ", str(diagonal_dominant)
    print "Finding Eigenvalues"
    try:
        eigvals_max = eigsh(spmat, k=10, which='LA', return_eigenvectors=False)
        eigvals_min = eigsh(spmat, k=10, which='SA', return_eigenvectors=False)
        positive_definite = np.all(eigvals_min > -1E-15)
        print "Matrix is Positive Definite: " + str(positive_definite)
        # we take the condition number ignoring the last six eigenvalues because
        # these will be rigid body modes in the unconstrained matrix
        print "Approximage Condition Number: %e"%(np.abs(np.max(eigvals_max))/
                                                  np.abs(eigvals_min[-7]))
        print eigvals_max
        print eigvals_min
    except:
        print "There was an error in the eigenvalue computation"
        print "This may mean there are issues with the input file."
        print "Try checking the validity of the input file."
