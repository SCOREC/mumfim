#! /lore/mersoj/anaconda2/bin/python
import numpy as np
import scipy
from scipy.sparse import csr_matrix
import stiffness_validity

mat_file = "./abaqus/Job-1_STIF2.mtx"
abaqus_mat_data = np.loadtxt(mat_file,
                            dtype={'names':('n1','d1','n2','d2','data'),'formats':('i8','i8','i8','i8','f8')},
                            delimiter=',')
row_idx = (abaqus_mat_data['n1']-1)*3+(abaqus_mat_data['d1']-1)
col_idx = (abaqus_mat_data['n2']-1)*3+(abaqus_mat_data['d2']-1)
mat_size = np.max(row_idx+1)
aba_mat = csr_matrix((abaqus_mat_data['data'], (row_idx, col_idx)),
            shape=(mat_size,mat_size))
aba_mat = aba_mat + scipy.sparse.tril(aba_mat, -1).transpose()
stiffness_validity.check_validity(aba_mat)
