import numpy as np
import scipy
from scipy.sparse import linalg
from scipy import linalg
from scipy.sparse import csr_matrix
from scipy import io

#mat_file = "/fasttmp/deoges/Delaunay/DELscript_automate/Del_networks/set_4/Del_2925L10_alpha0_3/Del_3_stiffness_STIF1.mtx"
#mat_file = "./Load-Only_STIF2.mtx"
mat_file = "./abaqus/Job-1_STIF2.mtx"
abaqus_mat_data = np.loadtxt(mat_file,
                            dtype={'names':('n1','d1','n2','d2','data'),'formats':('i8','i8','i8','i8','f8')},
                            delimiter=',')
row_idx = (abaqus_mat_data['n1']-1)*3+(abaqus_mat_data['d1']-1)
col_idx = (abaqus_mat_data['n2']-1)*3+(abaqus_mat_data['d2']-1)
mat_size = np.max(row_idx+1)
#assert(mat_size == np.max(col_idx)+1)
aba_mat = csr_matrix((abaqus_mat_data['data'], (row_idx, col_idx)),
            shape=(mat_size,mat_size))
#aba_mat = aba_mat + scipy.sparse.triu(aba_mat, 1).transpose()
aba_mat = aba_mat + scipy.sparse.tril(aba_mat, -1).transpose()

symmetric = np.abs((aba_mat-aba_mat.transpose())>1E-10).nnz == 0
print "Diagonal all positive: " + str(np.all(np.diag(aba_mat.toarray())>=-1E-15))
# directly compare data to avoid memory issues
#print symmetric
print "Symmetric: " + str(symmetric)
if symmetric:
	eigvals_max = scipy.sparse.linalg.eigsh(aba_mat, k=10, which='LA', 
									   return_eigenvectors=False)
	eigvals_min = scipy.sparse.linalg.eigsh(aba_mat, k=10, which='SA', 
									   return_eigenvectors=False)
positive_definite = False if not symmetric else np.all(eigvals_min > 0);
print "Abaqus Matrix is Positive Definite: " + str(positive_definite)
print "Approximate Condition Number: %e"%(np.abs(eigvals_max[0])/np.abs(eigvals_min[0]))
print eigvals_max
print eigvals_min

#scipy.io.savemat("aba_mat.mat", mdict={'aba_mat': aba_mat})
