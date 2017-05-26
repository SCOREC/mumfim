#include "RepresentVolElem.h"
#include "globals.h"
#include "RVE_Util.h"
#include "LagrangeMapping.h"// scorecutil stuff to be taken out eventually
#include <cstring>
#include <vector>
/*
  This file contains the functions
  create_element_shape_vars
  make_dRVEdFE
*/
namespace bio
{
  using namespace std;
  void MicroFO::create_element_shape_vars(double & vol,
                                          double * dvol,
                                          double * F)
    {
      double ddett[24];
      double dett;
      //double ttky[8],ttkx[8],ttkz[8],dkdy[8][24],dkdx[8][24],dkdz[8][24];
      double phi[8],phic[8],phie[8],phis[8];
      double dxrdx[24 * num_rve_doubles];
      // Calculation of dxrdx
      make_dRVEdFE(dxrdx,&initial_coords[0]);
//      make_dRVEdFE(dxrdx,coords);
      // Multiply dxrdx with fedisp to get rvedisp, the displacement of the RVE boundary
      double rvedisp[24];
      matrix_multiply(dxrdx,24,num_rve_doubles,fedisp,num_rve_doubles,1,rvedisp);
/*
      AMSI_DEBUG(std::cout <<"fedisp:"<<"("<<fedisp[0]<<","<<fedisp[1]<<","<<fedisp[2]<<")"
                                     <<", ("<<fedisp[3]<<","<<fedisp[4]<<","<<fedisp[5]<<")"
                                     <<", ("<<fedisp[6]<<","<<fedisp[7]<<","<<fedisp[8]<<")"
                                     <<", ("<<fedisp[9]<<","<<fedisp[10]<<","<<fedisp[11]<<")"
                                     <<std::endl);
      AMSI_DEBUG(std::cout <<"rvedisp:"<<"("<<rvedisp[0]<<","<<rvedisp[1]<<","<<rvedisp[2]<<")"
                                     <<", ("<<rvedisp[3]<<","<<rvedisp[4]<<","<<rvedisp[5]<<")"
                                     <<", ("<<rvedisp[6]<<","<<rvedisp[7]<<","<<rvedisp[8]<<")"
                                     <<", ("<<rvedisp[9]<<","<<rvedisp[10]<<","<<rvedisp[11]<<")"
                                     <<", ("<<rvedisp[12]<<","<<rvedisp[13]<<","<<rvedisp[14]<<")"
                                     <<", ("<<rvedisp[15]<<","<<rvedisp[16]<<","<<rvedisp[17]<<")"
                                     <<", ("<<rvedisp[18]<<","<<rvedisp[19]<<","<<rvedisp[20]<<")"
                                     <<", ("<<rvedisp[21]<<","<<rvedisp[22]<<","<<rvedisp[23]<<")"
                 <<std::endl);
*/
      // New implementation using deformation gradient
      double rvedisp[24] = {};
      getRVECornerDisp(F,rvedisp);

      // To test fiber only RVE (boundary conditions)
      if(FIBER_ONLY_SIZE_EFFECT_TEST)
      {
        for(int ii=0;ii<24;ii++)
          rvedisp[ii] = 0.0;
        double displace = fiber_network->sideCoord(FiberNetwork::RIGHT) / 5.0;
        double poisson_ratio = 0.3;
        double lateral_displace = -displace*poisson_ratio;
        // x
        rvedisp[3] = displace;
        rvedisp[9] = displace;
        rvedisp[15] = displace;
        rvedisp[21] = displace;
        // y
        rvedisp[13] = lateral_displace;
        rvedisp[16] = lateral_displace;
        rvedisp[19] = lateral_displace;
        rvedisp[22] = lateral_displace;
        // z
        rvedisp[2] = lateral_displace;
        rvedisp[5] = lateral_displace;
        rvedisp[14] = lateral_displace;
        rvedisp[17] = lateral_displace;
      }
      double bottom = fiber_network->sideCoord(FiberNetwork::BOTTOM);
      double top    = fiber_network->sideCoord(FiberNetwork::TOP);
      double left   = fiber_network->sideCoord(FiberNetwork::LEFT);
      double right  = fiber_network->sideCoord(FiberNetwork::RIGHT);
      double back   = fiber_network->sideCoord(FiberNetwork::BACK);
      double front  = fiber_network->sideCoord(FiberNetwork::FRONT);
      // Initial guess for the position of the micronodes
      // Apply linear interpolation initial guess to all nodes
      if(firstTimeThrough)
      {
        int num_nodes = fiber_network->numNodes();
        for(int ii = 0; ii < num_nodes; ii++)
        {
          const Node & ni = init_fiber_network->node(ii);
          Node n = fiber_network->node(ii);
          double x = ni.x;
          double y = ni.y;
          double z = ni.z;
          // weight the contribution of each corner's displacement by the sub-cube formed with the opposite corner (this gives the correct weights for linear interpolation)
          n.x +=
            ((right - x) * (top - y)     * (-back + z) * rvedisp[0]  +
             (-left + x) * (top - y)     * (-back + z) * rvedisp[3]  +
             (right - x) * (top - y)     * (front - z) * rvedisp[6]  +
             (-left + x) * (top - y)     * (front - z) * rvedisp[9]  +
             (right - x) * (-bottom + y) * (-back + z) * rvedisp[12] +
             (-left + x) * (-bottom + y) * (-back + z) * rvedisp[15] +
             (right - x) * (-bottom + y) * (front - z) * rvedisp[18] +
             (-left + x) * (-bottom + y) * (front - z) * rvedisp[21]) /
            ((right-left)*(top-bottom)*(front-back));
          n.y +=
            ((right - x) * (top - y)     * (-back + z) * rvedisp[1]  +
             (-left + x) * (top - y)     * (-back + z) * rvedisp[4]  +
             (right - x) * (top - y)     * (front - z) * rvedisp[7]  +
             (-left + x) * (top - y)     * (front - z) * rvedisp[10] +
             (right - x) * (-bottom + y) * (-back + z) * rvedisp[13] +
             (-left + x) * (-bottom + y) * (-back + z) * rvedisp[16] +
             (right - x) * (-bottom + y) * (front - z) * rvedisp[19] +
             (-left + x) * (-bottom + y) * (front - z) * rvedisp[22]) /
            ((right-left)*(top-bottom)*(front-back));
          n.z +=
            ((right - x) * (top - y)     * (-back + z) * rvedisp[2]  +
             (-left + x) * (top - y)     * (-back + z) * rvedisp[5]  +
             (right - x) * (top - y)     * (front - z) * rvedisp[8]  +
             (-left + x) * (top - y)     * (front - z) * rvedisp[11] +
             (right - x) * (-bottom + y) * (-back + z) * rvedisp[14] +
             (-left + x) * (-bottom + y) * (-back + z) * rvedisp[17] +
             (right - x) * (-bottom + y) * (front - z) * rvedisp[20] +
             (-left + x) * (-bottom + y) * (front - z) * rvedisp[23]) /
            ((right-left)*(top-bottom)*(front-back));
          fiber_network->setNode(ii,n);
        }
      }
      else
      {
        // After first time, initial guess is given by first order continuation
        int num_dofs = fiber_network->numDofs();
        double update[num_dofs];
        matrix_multiply(tdydxr.data(),num_dofs,24,rvedisp,24,1,update);
        // Linear interpolation is applied to boundary nodes
        int num_boundary_nodes = fiber_network->numBoundaryNodes(FiberNetwork::ALL);
        for(int ii = 0; ii < num_boundary_nodes; ii++)
        {
          int node = fiber_network->boundaryNode(FiberNetwork::ALL,ii);
          update[3*node  ] = 0.0;
          update[3*node+1] = 0.0;
          update[3*node+2] = 0.0;
          const Node & ni = init_fiber_network->node(node);
          Node n = fiber_network->node(node);
          double x = ni.x;
          double y = ni.y;
          double z = ni.z;
          // weight the contribution of each corner's displacement by the sub-cube formed with the opposite corner (this gives the correct weights for linear interpolation)
          n.x +=
            ((right - x) * (top - y)     * (-back + z) * rvedisp[0]  +
             (-left + x) * (top - y)     * (-back + z) * rvedisp[3]  +
             (right - x) * (top - y)     * (front - z) * rvedisp[6]  +
             (-left + x) * (top - y)     * (front - z) * rvedisp[9]  +
             (right - x) * (-bottom + y) * (-back + z) * rvedisp[12] +
             (-left + x) * (-bottom + y) * (-back + z) * rvedisp[15] +
             (right - x) * (-bottom + y) * (front - z) * rvedisp[18] +
             (-left + x) * (-bottom + y) * (front - z) * rvedisp[21]) /
            ((right-left)*(top-bottom)*(front-back));
          n.y +=
            ((right - x) * (top - y)     * (-back + z) * rvedisp[1]  +
             (-left + x) * (top - y)     * (-back + z) * rvedisp[4]  +
             (right - x) * (top - y)     * (front - z) * rvedisp[7]  +
             (-left + x) * (top - y)     * (front - z) * rvedisp[10] +
             (right - x) * (-bottom + y) * (-back + z) * rvedisp[13] +
             (-left + x) * (-bottom + y) * (-back + z) * rvedisp[16] +
             (right - x) * (-bottom + y) * (front - z) * rvedisp[19] +
             (-left + x) * (-bottom + y) * (front - z) * rvedisp[22]) /
            ((right-left)*(top-bottom)*(front-back));
          n.z +=
            ((right - x) * (top - y)     * (-back + z) * rvedisp[2]  +
             (-left + x) * (top - y)     * (-back + z) * rvedisp[5]  +
             (right - x) * (top - y)     * (front - z) * rvedisp[8]  +
             (-left + x) * (top - y)     * (front - z) * rvedisp[11] +
             (right - x) * (-bottom + y) * (-back + z) * rvedisp[14] +
             (-left + x) * (-bottom + y) * (-back + z) * rvedisp[17] +
             (right - x) * (-bottom + y) * (front - z) * rvedisp[20] +
             (-left + x) * (-bottom + y) * (front - z) * rvedisp[23]) /
            ((right-left)*(top-bottom)*(front-back));
          fiber_network->setNode(node,n);
        }
        // Linear interpolation applied to support nodes
//      if (SPECIFY_FIBER_TYPE)
        if (false)
        {
        // Extract nodes of support fibers.
          std::vector<int> supportNodes;
          int num_elements = fiber_network->numElements();
          for (int ii = 0; ii < num_elements; ii++)
          {
            const Element & e = fiber_network->element(ii);
            if (e.fiber_type != 0)
              supportNodes.push_back(e.node2_id);
          }
          int num_support_nodes = supportNodes.size();
          for (int ii = 0; ii < num_support_nodes; ii++)
          {
            int node = supportNodes[ii];
            update[3*node  ] = 0.0;
            update[3*node+1] = 0.0;
            update[3*node+2] = 0.0;
            const Node & ni = init_fiber_network->node(node);
            Node n = fiber_network->node(node);
            double x = ni.x;
            double y = ni.y;
            double z = ni.z;
            // weight the contribution of each corner's displacement by the sub-cube formed with the opposite corner (this gives the correct weights for linear interpolation)
            n.x +=
              ((right - x) * (top - y)     * (-back + z) * rvedisp[0]  +
               (-left + x) * (top - y)     * (-back + z) * rvedisp[3]  +
               (right - x) * (top - y)     * (front - z) * rvedisp[6]  +
               (-left + x) * (top - y)     * (front - z) * rvedisp[9]  +
               (right - x) * (-bottom + y) * (-back + z) * rvedisp[12] +
               (-left + x) * (-bottom + y) * (-back + z) * rvedisp[15] +
               (right - x) * (-bottom + y) * (front - z) * rvedisp[18] +
               (-left + x) * (-bottom + y) * (front - z) * rvedisp[21]) /
              ((right-left)*(top-bottom)*(front-back));
            n.y +=
              ((right - x) * (top - y)     * (-back + z) * rvedisp[1]  +
               (-left + x) * (top - y)     * (-back + z) * rvedisp[4]  +
               (right - x) * (top - y)     * (front - z) * rvedisp[7]  +
               (-left + x) * (top - y)     * (front - z) * rvedisp[10] +
               (right - x) * (-bottom + y) * (-back + z) * rvedisp[13] +
               (-left + x) * (-bottom + y) * (-back + z) * rvedisp[16] +
               (right - x) * (-bottom + y) * (front - z) * rvedisp[19] +
               (-left + x) * (-bottom + y) * (front - z) * rvedisp[22]) /
              ((right-left)*(top-bottom)*(front-back));
            n.z +=
              ((right - x) * (top - y)     * (-back + z) * rvedisp[2]  +
               (-left + x) * (top - y)     * (-back + z) * rvedisp[5]  +
               (right - x) * (top - y)     * (front - z) * rvedisp[8]  +
               (-left + x) * (top - y)     * (front - z) * rvedisp[11] +
               (right - x) * (-bottom + y) * (-back + z) * rvedisp[14] +
               (-left + x) * (-bottom + y) * (-back + z) * rvedisp[17] +
               (right - x) * (-bottom + y) * (front - z) * rvedisp[20] +
               (-left + x) * (-bottom + y) * (front - z) * rvedisp[23]) /
              ((right-left)*(top-bottom)*(front-back));
            fiber_network->setNode(node,n);
          } // End for loop through supportNodes array.
        }// End interpolation of support nodes.
        int num_nodes = fiber_network->numNodes();
        for(int ii = 0; ii < num_nodes; ii++)
        {
          Node n = fiber_network->node(ii);
          n.x += update[3*ii    ];
          n.y += update[3*ii + 1];
          n.z += update[3*ii + 2];
          fiber_network->setNode(ii,n);
        }
      }
      // Calculation of the deformed position of the RVE boundary
      for(int ii = 0; ii < 8; ii++)
      {
        rve[ii][0] += rvedisp[3*ii];
        rve[ii][1] += rvedisp[3*ii+1];
        rve[ii][2] += rvedisp[3*ii+2];
      }
      // Calculation of the RVE volume
      double gp[2], gw[2];
      get_gauss_points(gp,gw);
      // make sure volume is set to 0.0
      vol = 0.0;
      memset(dvol,0,24*sizeof(double));
      apf::Vector3 t[8];
      for(int ii = 0; ii < 8; ii++)
      {
        t[ii][0] = rve[ii][0];
        t[ii][1] = rve[ii][1];
        t[ii][2] = rve[ii][2];
      }
      // j = ii
      for(int ii = 0; ii < 2; ii++)
      {
        // k = jj
        for(int jj = 0; jj < 2; jj++)
        {
          // kj = kk
          for (int kk = 0; kk < 2; kk++)
          {
            tsfun(phi,phic,phie,phis,gp[ii],gp[jj],gp[kk]);
            //TODO: currently gp and gw are hard coded in RVE_Util.h.
            //gp[0] = -0.5774502....
            //gp[1] = 0.5774502 ....
            //gw[0] = 1.0;
            //gw[1] = 1.0;
            calc_dos(t,phic,phie,phis,&dett,ddett);
            vol += gw[ii] * gw[jj] * gw[kk] * dett;
            for(int nn = 0; nn < 24; nn++)
              dvol[nn] += gw[ii] * gw[jj] * gw[kk] * ddett[nn];
          }
        }
      }
    }
  // This is called with init_coords by create_element_shape_vars
  // and with coords by calc_femjacob_newmethod
  void MicroFO::make_dRVEdFE(double * dRVEdFE, double * lcoords)
  {
    /* Calculates the derivative dxrdx:
       derivative of the position of RVE edges (nodes?) with respect to the position of the FE nodes
    */
    int rve_dof = 24; //total number of DOF at each RVE = number of vertices (8) x DOF of each vertex (3)
    double xgp=0.0;
    double ygp=0.0;
    double zgp=0.0;
    /* Size of RVE in terms of physical domain in the
       x, y, and z directions.
    */
    //todo: pull this out of here
    for(int ii = 0; ii < num_rve_doubles * rve_dof; ii++)
      dRVEdFE[ii] = 0;
    /* find the Gauss (Integration) point where RVE is defined.
       Gauss point is in terms of Barycentric coordinates (natural coordinates) of
       macroscale tetrahedral finite element.
     */
    double gpx = gauss_pt[0];
    double gpy = gauss_pt[1];
    double gpz = gauss_pt[2];
    /* Set up the SCOREC lagrange mapping to
       1. determine xgp, ygp, and zgp, which defines the Gauss point
          in terms of the physical domain.
       2. determine u.x, u.y, and u.z, which define the coordinates of the RVE vertices with respect to
          the natural coordinates (isoparametric representation).
     */
    std::vector<mPoint> knots;
    for(int ii = 0; ii < num_element_nodes; ii++)
    {
      mPoint pt(lcoords[3*ii],
                      lcoords[3*ii+1],
                      lcoords[3*ii+2]);
      knots.push_back(pt);
    }
    LagrangeMapping lmp(0,&knots);
    // determine Gauss points in terms of physical domain.
    if(firstTimeThrough)
    {
      /*
      apf::Vector3 local(gpx,gpy,gpz);
      apf::Vector3 global;
      apf::mapLocalToGlobal(macro_element,local,global);
      xgp = global[0];
      ygp = global[1];
      zgp = global[2];
      */
      lmp.eval(gpx,gpy,gpz,xgp,ygp,zgp);
      /*
      //std::cout<<"xgp, ygp, zgp = "<<xgp<<", "<<ygp<<", "<<zgp<<std::endl;
      // Equivalent to:
      double xgp2 = 0.0; double ygp2 = 0.0; double zgp2 = 0.0;
      for (int ii = 0; ii < num_element_nodes; ii++)
      {
        xgp += lcoords[3*ii] * PHI(ii,gpx,gpy,gpz);
        ygp += lcoords[3*ii+1] * PHI(ii,gpx,gpy,gpz);
        zgp += lcoords[3*ii+2] * PHI(ii,gpx,gpy,gpz);
      }
      //std::cout<<"xgp2, ygp2, zgp2 = "<<xgp2<<", "<<ygp2<<", "<<zgp2<<std::endl;
      // where PHI is linear tetrahedral shape function defined at gpx, gpy, gpz.
      // lcoords (input to function) contains the physical position of the vertices
      // of the macroscale tetrahedral element.
      */
    }
    int ex_sgn =  1;
    int ey_sgn =  1;
    int ez_sgn = -1;
    double ex = 0; double ey = 0; double ez = 0;
    int count = 0;
    double e, n, l;
    // this looping results in the same corner ordering as the hardcoded rve.xyz value assignment during construction
    for(int ii = 0; ii < 2; ii++)
    {
      ey_sgn *= -1;
      ey = ey_sgn * half_rve_dim;
      /*
      if (ey_sgn < 0)
	ey = fiber_network->sideCoord(FiberNetwork::BOTTOM)/0.5 * half_rve_dim;
      else
	ey = fiber_network->sideCoord(FiberNetwork::TOP)/0.5 * half_rve_dim;
      */
      for(int jj = 0; jj < 2; jj++)
      {
	ez_sgn *= -1;
	ez = ez_sgn * half_rve_dim;
	/*
	if (ez_sgn < 0)
	  ez = fiber_network->sideCoord(FiberNetwork::BACK)/0.5 * half_rve_dim;
	else
	  ez = fiber_network->sideCoord(FiberNetwork::FRONT)/0.5 * half_rve_dim;
	*/
        for(int kk = 0; kk < 2; kk++)
        {
	  ex_sgn *= -1;
	  ex = ex_sgn * half_rve_dim;
	  /*
	  if (ex_sgn < 0)
	    ex = fiber_network->sideCoord(FiberNetwork::LEFT)/0.5 * half_rve_dim;
	  else
	    ex = fiber_network->sideCoord(FiberNetwork::RIGHT)/0.5 * half_rve_dim;
	  */
          /* calculate the coordinates of the vertices of the RVE with respect to physical domain.
            (RVEs are centered at the Gauss point in terms of physical domain (xgp,ygp,zgp)).
            Therefore, the vertices of the RVE is plus/minus half the RVE dimensions in each direction).
           */
          if(firstTimeThrough)
          {
            // Position of the RVE vertices with respect to the physical domain
            double xr = xgp + ex;
            double yr = ygp + ey;
            double zr = zgp + ez;
            // Position of the RVE vertices with respect to natural coordiantes of macroscale tetrahedral element.
            //
            double tmpu, tmpv, tmpw;
            lmp.invert(xr, yr, zr, tmpu, tmpv, tmpw);
            u[count][0] = e = tmpu;
            u[count][1] = n = tmpv;
            u[count][2] = l = tmpw;
          }
          else
          {
            e = u[count][0];
            n = u[count][1];
            l = u[count][2];
          }
          for(int ii = 0; ii < num_rve_doubles / 3; ii++)
          {
            /* For each vertex of RVE, calculate difference of linear tetrahedral shape functions
               defined on vertex positions (e,n,l) and Gauss point (gpx,gpy,gpz).
               Both are defined with respect to natural coordinates of macroscale tetrahedral element.
             */
            double diff = PHI(ii, e, n, l) - PHI(ii, gpx, gpy, gpz);
            double v = diff / rve_dim;
            dRVEdFE[(3*count)   * num_rve_doubles + ii*3    ] = v;
            dRVEdFE[(3*count+1) * num_rve_doubles + ii*3 + 1] = v;
            dRVEdFE[(3*count+2) * num_rve_doubles + ii*3 + 2] = v;
          }
          count++;
        }
      }
    }
  }
  void MicroFO::getRVECornerDisp(const double F[], double rvedisp[])
  {
    // Reference coordinates of RVE. X[node index][component index]:
    apf::Vector3 X[8];
    X[0][0] = -0.5; X[0][1] = -0.5; X[0][2] = 0.5;
    X[1][0] = 0.5;  X[1][1] = -0.5; X[1][2] = 0.5;
    X[2][0] = -0.5; X[2][1] = -0.5; X[2][2] = -0.5;
    X[3][0] = 0.5;  X[3][1] = -0.5; X[3][2] = -0.5;
    X[4][0] = -0.5; X[4][1] = 0.5;  X[4][2] = 0.5;
    X[5][0] = 0.5;  X[5][1] = 0.5;  X[5][2] = 0.5;
    X[6][0] = -0.5; X[6][1] = 0.5;  X[6][2] = -0.5;
    X[7][0] = 0.5;  X[7][1] = 0.5;  X[7][2] = -0.5;

    // incremental displacement, u_inc.
    double u_inc[3] = {};
    // total displacement, u.
    double u[3] = {}; double I = 0;
    for (int RVE_nd = 0; RVE_nd < 8; ++RVE_nd)
    {
      for (int i = 0; i < 3; ++i)
      {
	for (int j = 0; j < 3; ++j)
	{
	  if (i == j)
	    I = 1.0;
	  u[i] += (F[i*3 + j] - I)*X[RVE_nd][j];
	  I = 0.0;
	}
	// rve - X = previous displacement.
	u_inc[i] = u[i] - (rve[RVE_nd][i] - X[RVE_nd][i]);
      }
      rvedisp[RVE_nd * 3] = u_inc[0];
      rvedisp[RVE_nd * 3 + 1] = u_inc[1];
      rvedisp[RVE_nd * 3 + 2] = u_inc[2];
      u[0] = 0.0; u[1] = 0.0; u[2] = 0.0;
      u_inc[0] = 0.0; u_inc[1] = 0.0; u_inc[2] = 0.0;
    }
  }
  void getdRVEdFE(double* dRVEdFE, const double* FE_disp, const double* RVE_disp)
  {
    /* Dimensions of variables
       output:
              dRVEdFE: 24 x 12 = 288. 24 = 8 RVE corners x 3 dofs per corner. 12 = 4 FE nodes x 3 dofs per node.
       intput:
              FE_disp:  4 FE nodes x 3 dofs per node = 12.
	      RVE_disp: 8 corner nodex x 3 dofs per node = 24.
    */
    /*1. Find displacement at Gauss point. The displacement at this point is identical for both micro and macroscale. */
    
  }
}
