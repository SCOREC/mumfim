#! /usr/bin/env python
import numpy as np
from scipy.sparse import csr_matrix
import scipy.io as scio
import scipy.sparse
import scipy.sparse.linalg
import sys

class Node:
    def __init__(self, x, y, z, dof1, dof2, dof3):
        self.coords = [x,y,z]
        # dof numbers
        self.dof_num = [dof1, dof2, dof3]
class Element:
    def __init__(self, n1, n2, youngs_modulus, fiber_area):
        self.node1 = n1
        self.node2 = n2
        self.fiber_area = fiber_area
        self.youngs_modulus = youngs_modulus
    def get_nodes(self):
        return [self.node1, self.node2]
    def get_length(self):
        return np.linalg.norm(np.array(self.node2.coords)-np.array(self.node1.coords))
    def get_unit_vector(self):
        return (np.array(self.node2.coords)-np.array(self.node1.coords))/self.get_length()
    def get_k(self):
        L = self.get_length()
        L0 = L # for this case we just want the initial stiffness
        stretch = compute_stretch(L,L0)
        assert(np.isclose(stretch, 1))
        P = compute_PK1_stress(self.youngs_modulus,self.fiber_area,stretch)
        # Microscale stiffness Matrix
        n = self.get_unit_vector()
        LVec = L*n
        LNorm = np.linalg.norm(LVec) # note this should be L
        I = np.identity(3)
        dpdl = (self.youngs_modulus*self.fiber_area/2*(3*L**2/L0**3-1/L0))*n
        dfdu = np.outer(dpdl,n)+(1/LNorm*I- np.outer((1/LNorm**3)*LVec,LVec))*P
        ke = np.zeros((6,6))
        #ke[0:3,0:3] = dfdu
        #ke[0:3,3:6] = -dfdu
        #ke[3:6,0:3] = -dfdu
        #ke[3:6,3:6] = dfdu
        ke[0:3,0:3] = -dfdu
        ke[0:3,3:6] = dfdu
        ke[3:6,0:3] = dfdu
        ke[3:6,3:6] = -dfdu
        return ke
         
class FiberNetworkParams:
    def __init__(self, E, R):
       self.youngs_modulus = E
       self.fiber_radius = R
class FiberNetwork:
    def __init__(self, filename, param_filename):
        self.params, self.param_assignment = read_network_params(param_filename)
        self.nodes, self.elements = read_network(filename, self.params,
                                                 self.param_assignment)
        assert(len(self.elements) == len(self.param_assignment))
        # assume that the last dof added to the system is the max
        self.dof_max = self.nodes[-1].dof_num[-1]+1
        self.K = np.zeros((self.dof_max, self.dof_max))
    def constuct_global_stiffness(self):
        # sum element contributions into global stiffness matrix
        #elmts = [self.elements[0], self.elements[10]]
        for element in self.elements:
        #for element in elmts:
            ke = element.get_k()
            assert(np.allclose(ke, ke.T))
            global_dofs = element.node1.dof_num+element.node2.dof_num
            for i in range(6):
                for j in range(6):
                    P = global_dofs[i]
                    Q = global_dofs[j]
                    self.K[P,Q] += ke[i,j]
                    # forcing K to be symmetric because floating point
                    # ops messed this up
        self.K = csr_matrix(self.K)
        self.K = 0.5*(self.K+self.K.transpose())

def read_network_params(filename):
    """
    reads network and returns
    """
    param_assignment = np.loadtxt(filename,skiprows=2, dtype=np.int)
    with open(filename,'r') as f:
        line1 = f.readline()
        num_params, num_fibers = [int(x) for x in line1.split()]
        params = [None]*num_params
        for i in range(num_params):
            line = f.readline()
            # ignore the param number...e.g. assume params are in order
            _,R,E = [float(x) for x in line.split()]
            params[i] = FiberNetworkParams(E,R)
    return (params, param_assignment)

def read_network(filename, params, param_assignment):
    with open(filename, 'r') as f:
        line1 = f.readline()
        num_nodes, num_fibers, _ = [int(x) for x in line1.split()]
        nodes = [None]*num_nodes
        fibers = [None]*num_fibers
        dof = 0;
        for i in range(num_nodes):
            x,y,z = [float(x) for x in f.readline().split()]
            nodes[i] = Node(x, y, z, dof, dof+1, dof+2)
            dof = dof + 3
        for i in range(num_fibers):
            n1, n2 = [int(x) for x in f.readline().split()]
            E = params[param_assignment[0]].youngs_modulus
            A0 = np.pi * np.power(params[param_assignment[0]].fiber_radius,2)
            print A0, E
            fibers[i] = Element(nodes[n1], nodes[n2], E, A0)
        return (nodes, fibers)

def compute_stretch(L,L0):
    return L/L0

def compute_PK1_stress(E, A0, stretch):
    return E*1/2*(stretch**2-1)*A0*stretch

def main():
    if len(sys.argv) != 2:
        print "Usage: "+sys.argv[0]+" mesh_file"
        sys.exit()
    print sys.argv
    mesh_file_name = sys.argv[1]
    mesh_params_file_name = mesh_file_name+".params"
    fn = FiberNetwork(mesh_file_name, mesh_params_file_name)
    #fn = FiberNetwork("./del_rho120_alpha45_new_3.txt",
    #                  "./del_rho120_alpha45_new_3.txt.params")

    print "Constructing global stiffness"
    fn.constuct_global_stiffness()
    print "Converting to a sparse matrix"
    K = fn.K
    print "Number of nonzeros: " + str(K.nnz)

    symmetric = np.abs((K-K.transpose()) >=1E-15).nnz == 0
    # note this is not the best method for determining if it is symmetric
    # due to the fact that we are using an absolute tolerance...
    print "Symmetric: " + str(symmetric)
    print "Diagonal Positive: " + str(np.all(K.diagonal() > -1E-15))
    print "Finding Eigenvalues"
    #eigvals_max = scipy.sparse.linalg.eigs(K, k=2, which='LR',
    #                                        return_eigenvectors=False)
    #eigvals_min = scipy.sparse.linalg.eigs(K, k=2, which='SR',
    #                                        return_eigenvectors=False)
    #positive_definite = np.all(eigvals_min > -1E-15)
    #print "Matrix is Positive Definite: " + str(positive_definite)
    #print "Approximage Condition Number: %e"%(np.abs(np.max(eigvals_max))/np.abs(np.min(eigvals_min)))
    #print eigvals_max
    #print eigvals_min

    print "Writing to matrix market file"
    scio.mmio.mmwrite("globalK", K)
    scio.savemat("globalK.mat", mdict={'K': K})

if __name__ == "__main__":
    main()
