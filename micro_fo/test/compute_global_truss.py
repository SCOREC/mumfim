#! /lore/mersoj/anaconda2/bin/python
"""
This module reads a network mesh generated for the biotissue code and presents
some information for debugging.
"""
import sys
import numpy as np
from scipy.sparse import csr_matrix
import scipy.io as scio
import matplotlib.pyplot as plt
#import scipy.sparse.linalg
from scipy.sparse.linalg import eigsh, eigs
import os
import stiffness_validity

class Node(object):
    """
    A simple mesh node.
    """
    def __init__(self, num, x, y, z, dof1, dof2, dof3):
        self.num = num
        self.coords = [x, y, z]
        # dof numbers
        self.dof_num = [dof1, dof2, dof3]
class Element(object):
    """
    A truss element
    """
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
        stretch = compute_stretch(L, L0)
        assert np.isclose(stretch, 1)
        P = compute_PK1_stress(self.youngs_modulus, self.fiber_area, stretch)
        # Microscale stiffness Matrix
        n = self.get_unit_vector()
        LVec = L*n
        LNorm = np.linalg.norm(LVec) # note this should be L
        I = np.identity(3)
        dpdl = (self.youngs_modulus*self.fiber_area/2*(3*L**2/L0**3-1/L0))*n
        dfdu = np.outer(dpdl, n)+(1/LNorm*I- np.outer((1/LNorm**3)*LVec, LVec))*P
        ke = np.zeros((6, 6))
        #ke[0:3, 0:3] = -dfdu
        #ke[0:3, 3:6] = dfdu
        #ke[3:6, 0:3] = dfdu
        #ke[3:6, 3:6] = -dfdu

        ke[0:3, 0:3] = dfdu
        ke[0:3, 3:6] = -dfdu
        ke[3:6, 0:3] = -dfdu
        ke[3:6, 3:6] = dfdu
        return ke

class FiberNetworkParams(object):
    def __init__(self, E, R):
        self.youngs_modulus = E
        self.fiber_radius = R
class FiberNetwork(object):
    def __init__(self, filename, param_filename):
        self.params, self.param_assignment = read_network_params(param_filename)
        self.nodes, self.elements = read_network(filename, self.params,
                                                 self.param_assignment)
        assert len(self.elements) == len(self.param_assignment)
        # assume that the last dof added to the system is the max
        self.dof_max = self.nodes[-1].dof_num[-1]+1
        # consider constructing the sparsity pattern and directly creating csr here
        self.K = np.zeros((self.dof_max, self.dof_max))
    def constuct_global_stiffness(self):
        # sum element contributions into global stiffness matrix
        #elmts = [self.elements[0], self.elements[10]]
        for element in self.elements:
        #for element in elmts:
            ke = element.get_k()
            assert np.allclose(ke, ke.T)
            global_dofs = element.node1.dof_num+element.node2.dof_num
            for i in range(6):
                for j in range(6):
                    P = global_dofs[i]
                    Q = global_dofs[j]
                    self.K[P, Q] += ke[i, j]
                    # forcing K to be symmetric because floating point
                    # ops messed this up
        self.K = csr_matrix(self.K)
        self.K = 0.5*(self.K+self.K.transpose())

def read_network_params(filename):
    """
    reads network parameters
    """
    param_assignment = np.loadtxt(filename, skiprows=2, dtype=np.int)
    with open(filename, 'r') as f:
        line1 = f.readline()
        num_params, _ = [int(x) for x in line1.split()]
        params = [None]*num_params
        for i in range(num_params):
            line = f.readline()
            # ignore the param number...e.g. assume params are in order
            _, R, E = [float(x) for x in line.split()]
            params[i] = FiberNetworkParams(E, R)
    return (params, param_assignment)

def read_network(filename, params, param_assignment):
    """
    reads the network
    """
    with open(filename, 'r') as f:
        line1 = f.readline()
        num_nodes, num_fibers, _ = [int(x) for x in line1.split()]
        nodes = [None]*num_nodes
        fibers = [None]*num_fibers
        dof = 0
        for i in range(num_nodes):
            x, y, z = [float(v) for v in f.readline().split()]
            nodes[i] = Node(i, x, y, z, dof, dof+1, dof+2)
            dof = dof + 3
        for i in range(num_fibers):
            n1, n2 = [int(v) for v in f.readline().split()]
            E = params[param_assignment[0]].youngs_modulus
            A0 = np.pi * np.power(params[param_assignment[0]].fiber_radius, 2)
            fibers[i] = Element(nodes[n1], nodes[n2], E, A0)
        return (nodes, fibers)

def compute_stretch(L, L0):
    """
    compute the stretch ratio
    """
    return L/L0

def compute_PK1_stress(E, A0, stretch):
    """
    compute the 1st piola kirchoff stress
    """
    return E*1/2*(stretch**2-1)*A0*stretch

def get_graph_from_network(elements):
    import networkx as nx
    fiberGraph = nx.Graph()
    for i,element in enumerate(elements):
        n1 = element.node1.num
        n2 = element.node2.num
        fiberGraph.add_edge(n1,n2, weight=i)
    return fiberGraph
    

def main():
    """
    Main function reads the command line arguments and drives the mesh loader
    to compute important statistics
    """
    if len(sys.argv) != 2:
        print "Usage: "+sys.argv[0]+" mesh_file/yaml file"
        sys.exit()
    mesh_file_names = []
    mesh_file_name = sys.argv[1]
    if os.path.splitext(mesh_file_name)[1] == ".yaml":
        try:
            import yaml
        except ImportError:
            print "You must have yaml module installed to read yaml file for mesh"
            sys.exit()

        with open(mesh_file_name, 'r') as stream:
            doc = yaml.load(stream)
            for case in doc:
                mesh_file_names.append(case["Problem Definition"]["Mesh"])
    else:
        mesh_file_names.append(mesh_file_name)

    for mesh_file_name in mesh_file_names:
        mesh_file_name = os.path.abspath(mesh_file_name)
        print "Loading Mesh: " + mesh_file_name
        mesh_params_file_name = mesh_file_name+".params"
        fn = FiberNetwork(mesh_file_name, mesh_params_file_name)


        print "Constructing global stiffness"
        fn.constuct_global_stiffness()
        name = os.path.splitext(os.path.split(mesh_file_name)[1])[0]
        output_name = os.path.abspath("./"+name)+"_globalK"
        print "Writing to matrix market file to: " + output_name + ".mtx"
        scio.mmio.mmwrite(output_name, fn.K)
        #print "Computing Basic Stiffness matrix information"
        #stiffness_validity.check_validity(fn.K)

if __name__ == "__main__":
    main()
