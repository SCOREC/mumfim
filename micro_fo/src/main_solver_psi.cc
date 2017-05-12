#include "RepresentVolElem.h"
#include "globals.h"
#include "RVE_Util.h"
#include "Sparskit_Externs.h"
#include "Util.h"
#include <iostream>
#include <cstring>
namespace bio
{
  int MicroFO::main_solver(double * coords_loc,double * fedisp)
  {
    int failed = 0;
    // Give initial guess for the position of the micronodes,
    // Calculate the position of the RVE boundary and
    //double vol_p = vol;
    //double dvol_p =
    create_element_shape_vars(vol,&dvol[0],&fedisp[0]);
    // Jacobian matrix of the microscopic problem
    matrix.resize(sparse_structure->numNonzeros());
    matrix_axial.resize(sparse_structure->numNonzeros());
    // Fill vector with x,y,z,(rx,ry,rz) coords of nodes
    update_coordinate_vector();
    // Solution of the microscopic problem
    failed = Solver();
/*
    if(failed)
    {
      vol = vol_p;
    }
*/
    return failed;
  }
  void MicroFO::post_processing(double * sigma, double * dSdx, double * Q, double & fem_res_norm)
  {
    /* Sum up forces of fiber nodes that lie on the boundary of the RVE. This is the summation part of Eq. (7) in
       T. Stylianopoulos, V. H. Barocas, Comput. Methods Appl. Mech. Engrg. 196 (2007) 2981-2990:
       \sum_{boundary cross-links} x_i F_j.
    */
    double stress[6];
    calc_stress(stress);
    // Calculate the components of the Q term of the macroscopic stress balance, loc_vastrx, loc_vastry, loc_vastrz
    avg_vol_stress(stress,Q[0],Q[1],Q[2],vol,fem_res_norm);
    /* calc_tdydxr calculates tdydxr, which is change of fiber node positions as a function of RVE vertex positions.
       calc_femjacob_method calculates dSdx, which is the change of macroscale stress as a function of the macroscale
       FE node positions.
    */
    calc_tdydxr();
    calc_femjacob_newmethod(dSdx,vol,dvol,stress);
    /* Calculate the volume averaged (macroscale) stresses. This is the volume averaging part of Eq. (7) in
       T. Stylianopoulos, V. H. Barocas, Comput. Methods Appl. Mech. Engrg. 196 (2007) 2981-2990:
       S_ij = 1/V \sum_{boundary cross-links} x_i F_j
       Volume averaged stresses are dimensionalized by multiplying by scale_conversion.
    */
    for(int ii = 0; ii < 6; ++ii)
      sigma[ii] = (stress[ii] / vol) * scale_conversion;
    // todo (m) : hacky, change/fix this
    /*
    double orientation_tensor[9]={};
    calcFiberOrientationTensor(*fiber_network,orientation_tensor);
    for (int ii=0; ii<9; ii++)
      rve_info[4*3*6+9+ii] = orientation_tensor[ii];
    */
    double P2 = 0.0;
    calcP2(*fiber_network,P2);
    rve_info[4 * 3 * 6 + 9 + 0] = P2;
  }
  // For size effect testing. Hardcoded
  double MicroFO::calc_stiffness()
  {
    double traction = 0.0;
    int num_boundary_nodes = fiber_network->numBoundaryNodes(FiberNetwork::RIGHT);
    for(int ii = 0; ii < num_boundary_nodes; ii++)
    {
      int n = fiber_network->boundaryNode(FiberNetwork::RIGHT,ii);
      traction += force_vector[n*3];
    }
    // define a displacement (10% strain)
    double displace = fiber_network->sideCoord(FiberNetwork::RIGHT) / 5.0;
    double poisson_ratio = 0.3;
    double lateral_displace = -displace*poisson_ratio;
    double lateral_length = 2 * fiber_network->sideCoord(FiberNetwork::RIGHT) + lateral_displace;
    double stiffness = traction/pow(lateral_length,2.0);
    return stiffness;
  }
  void MicroFO::calc_stress(double * stress)
  {
    double xfx = 0.0, xfy = 0.0, xfz = 0.0;
    double yfx = 0.0, yfy = 0.0, yfz = 0.0;
    double zfx = 0.0, zfy = 0.0, zfz = 0.0;
    /* for decomposing forces on x faces. "p" and "n" denote positive and negative, respectively */
    //double fxp = 0.0, fxn = 0.0;
    int num_dofs_per_node = fiber_network->dofsPerNode();
    int num_boundary_nodes = fiber_network->numBoundaryNodes(FiberNetwork::ALL);
    for(int ii = 0; ii < num_boundary_nodes; ii++)
    {
      int n = fiber_network->boundaryNode(FiberNetwork::ALL,ii);
      int dof = n * num_dofs_per_node;
      const Node & node = fiber_network->node(n);
      if (!TRUSS) //beam
      {
        xfx += force_vector[dof] * node.x;
        yfx += force_vector[dof] * node.y;
        zfx += force_vector[dof] * node.z;
        xfy += force_vector[dof + 1] * node.x;
        yfy += force_vector[dof + 1] * node.y;
        zfy += force_vector[dof + 1] * node.z;
        xfz += force_vector[dof + 2] * node.x;
        yfz += force_vector[dof + 2] * node.y;
        zfz += force_vector[dof + 2] * node.z;
      }
      else
      {
        xfx += force_vector[dof] * node.x;
        yfx += force_vector[dof] * node.y;
        zfx += force_vector[dof] * node.z;
        xfy += force_vector[dof + 1] * node.x;
        yfy += force_vector[dof + 1] * node.y;
        zfy += force_vector[dof + 1] * node.z;
        xfz += force_vector[dof + 2] * node.x;
        yfz += force_vector[dof + 2] * node.y;
        zfz += force_vector[dof + 2] * node.z;
      }
    }
    stress[0] = xfx;
    stress[1] = 0.5 * (xfy + yfx);
    stress[2] = 0.5 * (zfx + xfz);
    stress[3] = yfy;
    stress[4] = 0.5 * (yfz + zfy);
    stress[5] = zfz;

    /*
  std::cout << "Microscale stress: " << stress[0] << " " << stress[1] << " "
  << stress[2] << " " << stress[3] << " "
  << stress[4] << " " << stress[5] << std::endl;
    */
  }
  void MicroFO::calc_tdydxr()
  {
    buffers->zero();
    int num_dofs = fiber_network->numDofs();
    double drop_tol = 0.0;
    int ierr = 0;
    int len = buffers->matrixLength();
    ilut_(&num_dofs,
          &matrix[0],
          sparse_structure->getCols(),
          sparse_structure->getRows(),
          &num_dofs,
          &drop_tol,
          buffers->matrixBuffer(),
          buffers->colsBuffer(),
          buffers->rowsBuffer(),
          &len,
          buffers->doubleWorkBuffer(),
          buffers->intWorkBuffer(),
          &ierr);
    if(ierr != 0)
      std::cerr << "Error: error value return from ilut_ is " << ierr << std::endl;
    // calculation of dridx: derivative of the residuals of the microscopic problem with respect to the edges of the RVE
    double dridx[num_dofs * 24];
    calc_dridx(&dridx[0]);
    double solution[num_dofs];
    memset(&solution[0],0,sizeof(double)*num_dofs);
    zero(force_vector);
    for(int ii = 0; ii < 24; ii++)
    {
      for(int jj = 0; jj < num_dofs; jj++)
        force_vector[jj] = dridx[jj * 24 + ii];
      lusol_(&num_dofs,
             &force_vector[0],
             &solution[0],
             buffers->matrixBuffer(),
             buffers->colsBuffer(),
             buffers->rowsBuffer());
      for(int jj = 0; jj < num_dofs; jj++)
        tdydxr[jj * 24 + ii] = solution[jj];
      memset(&solution[0],0,sizeof(double)*num_dofs);
    }
  }
  void MicroFO::calc_femjacob_newmethod(double * dSdx,
                                        double vol,
                                        double * dvol,
                                        double * stress)
  {
    int num_dofs = fiber_network->numDofs();
    buffers->zero();
    double dSdxr[144] = {};
    double dsdxr[144] = {};
    double * dxrdx = new double[24 * num_rve_doubles]();
//      memset(dSdx,0,sizeof(double)*6*num_rve_doubles);
//      memset(dSdxr,0,sizeof(double)*144);
//      memset(dsdxr,0,sizeof(double)*144);
    // The purpose of this function is to calculate the dSdx. In order to do this we have first
    // to calculate some other derivatives and particularly the dridx, dydx, dsdxr, dSdxr, dxrdx
    /* Multiply ttdSdy and tdydxr to obtain dsdxr.
       - ttdSdy is the change of microscale stresses on the boundary of RVE as function of fiber node positions.
       Calculated with calc_pre_cond function in solver_gmres_df.cc.
       - tdydxr is the change of fiber node positions as function of RVE vertex positions.
       Calculated with calc_tdydxr function in main_solver_psi.cc.
       - dsdxr is the change of microscale stresses on the boundary of RVE as function of RVE vertex positions.
    */
    matrix_multiply(ttdSdy.data(),6,num_dofs,tdydxr.data(),num_dofs,24,dsdxr);
    /* dSdxr is the macroscopic (volume-averaged) stress as a function of the RVE vertex positions.
       The expression to determine dSdxr is obtained (I believe) by taking the directional derivative
       of the volume averaged stresses with respect to the positions of the RVE vertex positions.
    */
    for(int jj = 0; jj < 24; jj++)
    {
      dSdxr[0*24 + jj] = ((dsdxr[0*24 + jj]/vol) - (dvol[jj]*stress[0]) / (vol*vol)) * scale_conversion;
      dSdxr[1*24 + jj] = ((dsdxr[1*24 + jj]/vol) - (dvol[jj]*stress[1]) / (vol*vol)) * scale_conversion;
      dSdxr[2*24 + jj] = ((dsdxr[2*24 + jj]/vol) - (dvol[jj]*stress[2]) / (vol*vol)) * scale_conversion;
      dSdxr[3*24 + jj] = ((dsdxr[3*24 + jj]/vol) - (dvol[jj]*stress[3]) / (vol*vol)) * scale_conversion;
      dSdxr[4*24 + jj] = ((dsdxr[4*24 + jj]/vol) - (dvol[jj]*stress[4]) / (vol*vol)) * scale_conversion;
      dSdxr[5*24 + jj] = ((dsdxr[5*24 + jj]/vol) - (dvol[jj]*stress[5]) / (vol*vol)) * scale_conversion;
    }
    // calculate the derivative dxrdx: derivative of the position of the edges of the RVE with respect to the position of the FE nodes
//      make_dRVEdFE(dxrdx,coords);
    make_dRVEdFE(dxrdx,&initial_coords[0]);
    // Multiply dSdxr with dxrdx to give dSdx
    matrix_multiply(dSdxr,6,24,dxrdx,24,num_rve_doubles,dSdx);
    delete [] dxrdx;
  }
  void MicroFO::calc_quantities(int side,
                                double * dridx,
                                int node,
                                int local_node1,
                                int local_node2,
                                int local_node3,
                                int local_node4)
  {
    // Calculation of dridx
    const Node & initial = init_fiber_network->node(node);
    double x = initial.x;
    double y = initial.y;
    double z = initial.z;
    int dofs_per_node = fiber_network->dofsPerNode();
    double bottom = fiber_network->sideCoord(FiberNetwork::BOTTOM);
    double top    = fiber_network->sideCoord(FiberNetwork::TOP);
    double left   = fiber_network->sideCoord(FiberNetwork::LEFT);
    double right  = fiber_network->sideCoord(FiberNetwork::RIGHT);
    double back   = fiber_network->sideCoord(FiberNetwork::BACK);
    double front  = fiber_network->sideCoord(FiberNetwork::FRONT);
    double area0 = (front-back)
      * (top-bottom)
      * (right-left); // area of the whole RVE ... volume
    double area1,area2,area3,area4;
    // initialize area variables.
    area1=0.0; area2=0.0; area3=0.0; area4=0.0;
    switch(side)
    {
    case FiberNetwork::BOTTOM:
      area1 = (right - x) * (top - y) * (front - z);
      area2 = (-left + x) * (top - y) * (front - z);
      area3 = (-left + x) * (top - y) * (-back + z);
      area4 = (right - x) * (top - y) * (-back + z);
      break;
    case FiberNetwork::RIGHT:
      area1= (top - y)     * (front - z) * (-left + x);
      area2= (top - y)     * (-back + z) * (-left + x);
      area3= (-bottom + y) * (-back + z) * (-left + x);
      area4= (-bottom + y) * (front - z) * (-left + x);
      break;
    case FiberNetwork::TOP:
      area1=(right - x) * (-bottom + y) * (front - z);
      area2=(-left + x) * (-bottom + y) * (front - z);
      area3=(-left + x) * (-bottom + y) * (-back + z);
      area4=(right - x) * (-bottom + y) * (-back + z);
      break;
    case FiberNetwork::LEFT:
      area1=(top - y)     * (front - z) * (right - x);
      area2=(top - y)     * (-back + z) * (right - x);
      area3=(-bottom + y) * (-back + z) * (right - x);
      area4=(-bottom + y) * (front - z) * (right - x);
      break;
    case FiberNetwork::BACK:
      area1=(right - x) * (top - y)     * (front - z);
      area2=(-left + x) * (top - y)     * (front - z);
      area3=(-left + x) * (-bottom + y) * (front - z);
      area4=(right - x) * (-bottom + y) * (front - z);
      break;
    case FiberNetwork::FRONT:
      area1=(right - x) * (top - y)     * (-back + z);
      area2=(-left + x) * (top - y)     * (-back + z);
      area3=(-left + x) * (-bottom + y) * (-back + z);
      area4=(right - x) * (-bottom + y) * (-back + z);
      break;
    }
    dridx[(node * dofs_per_node)     * 24 + (local_node1 * dofs_per_node)    ] = area1 / area0;
    dridx[(node * dofs_per_node)     * 24 + (local_node2 * dofs_per_node)    ] = area2 / area0;
    dridx[(node * dofs_per_node)     * 24 + (local_node3 * dofs_per_node)    ] = area3 / area0;
    dridx[(node * dofs_per_node)     * 24 + (local_node4 * dofs_per_node)    ] = area4 / area0;
    dridx[(node * dofs_per_node + 1) * 24 + (local_node1 * dofs_per_node + 1)] = area1 / area0;
    dridx[(node * dofs_per_node + 1) * 24 + (local_node2 * dofs_per_node + 1)] = area2 / area0;
    dridx[(node * dofs_per_node + 1) * 24 + (local_node3 * dofs_per_node + 1)] = area3 / area0;
    dridx[(node * dofs_per_node + 1) * 24 + (local_node4 * dofs_per_node + 1)] = area4 / area0;
    dridx[(node * dofs_per_node + 2) * 24 + (local_node1 * dofs_per_node + 2)] = area1 / area0;
    dridx[(node * dofs_per_node + 2) * 24 + (local_node2 * dofs_per_node + 2)] = area2 / area0;
    dridx[(node * dofs_per_node + 2) * 24 + (local_node3 * dofs_per_node + 2)] = area3 / area0;
    dridx[(node * dofs_per_node + 2) * 24 + (local_node4 * dofs_per_node + 2)] = area4 / area0;
  }
  void MicroFO::calc_dridx(double * dridx)
  {
    // need to assert that the allocated memory for dridx is enough to hold num_dofs*24
    int num_dofs = fiber_network->numDofs();
    memset(dridx,0,num_dofs*24*sizeof(double));
    int local_node[8];
    for(int ii = 0; ii < 8; ii++)
      local_node[ii] = ii;
    //int node;
    int verts[4];
    //int num;
    //int * bounds;
    for(int side = FiberNetwork::TOP;
        side != FiberNetwork::ALL;
        side++)
    {
      switch(side)
      {
      case FiberNetwork::BOTTOM: // y = 0
        verts[0] = 2;
        verts[1] = 3;
        verts[2] = 1;
        verts[3] = 0;
        break;
      case FiberNetwork::RIGHT: // x = x0
        verts[0] = 3;
        verts[1] = 1;
        verts[2] = 5;
        verts[3] = 7;
        break;
      case FiberNetwork::TOP: // y = y0
        verts[0] = 6;
        verts[1] = 7;
        verts[2] = 5;
        verts[3] = 4;
        break;
      case FiberNetwork::LEFT: // x = 0
        verts[0] = 2;
        verts[1] = 0;
        verts[2] = 4;
        verts[3] = 6;
        break;
      case FiberNetwork::BACK: // z = 0
        verts[0] = 2;
        verts[1] = 3;
        verts[2] = 7;
        verts[3] = 6;
        break;
      case FiberNetwork::FRONT: // z = z0
        verts[0] = 0;
        verts[1] = 1;
        verts[2] = 5;
        verts[3] = 4;
        break;
      }
      FiberNetwork::Side s = (FiberNetwork::Side)side;
      int num = fiber_network->numBoundaryNodes(s);
      for(int jj = 0; jj < num; jj++)
      {
        int n = fiber_network->boundaryNode(s,jj);
        calc_quantities(s,
                        dridx,
                        n,
                        local_node[verts[0]],
                        local_node[verts[1]],
                        local_node[verts[2]],
                        local_node[verts[3]]);
      }
    }
  }
}
