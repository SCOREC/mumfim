namespace bio
{
  /*
  void MicroFO::avg_vol_stress(double * stress,
                               double & loc_vastrx,
                               double & loc_vastry,
                               double & loc_vastrz,
                               double vol,
                               double fem_res_norm)
  {
    double gp[num_gauss_pts];
    double gw[num_gauss_pts];
    int nrve = 1;
    // This function is used for the calculation of: loc_vastrx, loc_vastry, loc_vastrz
    // These variables are the x,y,z components of the Q term in the macroscopic stress balance
    loc_vastrx = 0.0;
    loc_vastry = 0.0;
    loc_vastrz = 0.0;
    double gauss_pt[3];
    if(false)
    {
      double fi = 0;
      double fj = 0;
      double fk = 0;
      double dundx = 0;
      double dundy = 0;
      double dundz = 0;
      // one for each vertex on the RVE cube...
      double ux[8], uy[8], uz[8];
      double Fx[8], Fy[8], Fz[8];
      for(int ii = 0; ii < 8; ii++)
      {
        // FE displacements
        ux[ii] = coords[ii * 3]     - initial_coords[ii * 3];
        uy[ii] = coords[ii * 3 + 1] - initial_coords[ii * 3 + 1];
        uz[ii] = coords[ii * 3 + 2] - initial_coords[ii * 3 + 2];
        // FE positions
        Fx[ii] = coords[ii * 3];
        Fy[ii] = coords[ii * 3 + 1];
        Fz[ii] = coords[ii * 3 + 2];
      }
      get_gauss_points(gp,gw);
      // Find to which gauss point the current RVE corresponds
      int ll = 0;
      for(int ii = 0; ii < 2; ii++)
      {
        gauss_pt[1] = gp[ii];
        for(int jj = 0; jj < 2; jj++)
        {
          gauss_pt[2] = gp[jj];
          for(int kk = 0; kk < 2; kk++)
          {
            gauss_pt[0] = gp[kk];
            ll++;
            if(ll == nrve)
              break;
          }
        }
      }
      // CALCULATE THE Q TERM IN  MOMENTUM BALANCE
      double xyz[3];
      for(int side = FiberNetwork::TOP;
          side != FiberNetwork::ALL;
          side++)
      {
        xyz[0] = gauss_pt[0];
        xyz[1] = gauss_pt[1];
        xyz[2] = gauss_pt[2];
        switch(side)
        {
        case FiberNetwork::BOTTOM:
        {
          xyz[1] -= half_rve_dim;
          break;
        }
        case FiberNetwork::RIGHT:
        {
          xyz[0] += half_rve_dim;
          break;
        }
        case FiberNetwork::TOP:
        {
          xyz[1] += half_rve_dim;
          break;
        }
        case FiberNetwork::LEFT:
        {
          xyz[0] -= half_rve_dim;
          break;
        }
        case FiberNetwork::FRONT:
        {
          xyz[2] += half_rve_dim;
          break;
        }
        case FiberNetwork::BACK:
        {
          xyz[2] -= half_rve_dim;
          break;
        }
        }
        // CALCULATE THE MICRO-SCOPIC TERM OF THE Q
        double phi[8],phic[8],phie[8],phis[8];
        double dudx,dudy,dudz;
        double dvdx,dvdy,dvdz;
        double dgdx,dgdy,dgdz;
        double dett;
        double norm = 0.0;
        double normal[3];
        int num_local_boundary = fiber_network->numBoundaryNodes((FiberNetwork::Side) side);
        for(int ii = 0; ii < num_local_boundary; ii++)
        {
          double loc[3];
          int node = fiber_network->boundaryNode((FiberNetwork::Side)side,
                                                 ii);
          const Node & n = fiber_network->node(node);
          loc[0] = 2 * (n.x) * half_rve_dim + xyz[0];
          loc[1] = 2 * (n.y) * half_rve_dim + xyz[1];
          loc[2] = 2 * (n.z) * half_rve_dim + xyz[2];
          tsfun(phi,phic,phie,phis,loc[0],loc[1],loc[2]);
          calc_deriv(Fx,Fy,Fz,
                     ux,uy,uz,
                     phic,phie,phis,
                     &dudx,&dudy,&dudz,
                     &dvdx,&dvdy,&dvdz,
                     &dgdx,&dgdy,&dgdz,
                     &dett);
          if((fabs(force_vector[3*node  ])<1e-15) &&
             (fabs(force_vector[3*node+1])<1e-15) &&
             (fabs(force_vector[3*node+2])<1e-15))
            continue;
          if(side == FiberNetwork::TOP)
          {
            norm = sqrt(dvdx*dvdx + (dvdy-1)*(dvdy-1) + dvdz*dvdz);
            normal[0] = dvdx/norm;
            normal[1] = (dvdy-1)/norm;
            normal[2] = dvdz/norm;
          }
          else if(side == FiberNetwork::LEFT)
          {
            norm = sqrt(dvdx*dvdx + (dvdy-1)*(dvdy-1) + dvdz*dvdz);
            normal[0] = dvdx/norm;
            normal[1] = (1-dvdy)/norm;
            normal[2] = dvdz/norm;
          }
          else if(side == FiberNetwork::BOTTOM)
          {
            norm = sqrt(dudy*dudy + (dudx-1)*(dudx-1) + dudz*dudz);
            normal[0] = (1-dudx)/norm;
            normal[1] = dudy/norm;
            normal[2] = dudz/norm;
          }
          else if(side == FiberNetwork::RIGHT)
          {
            norm = sqrt(dudy*dudy + (dudx-1)*(dudx-1) + dudz*dudz);
            normal[0] = (dudx-1)/norm;
            normal[1] = dudy/norm;
            normal[2] = dudz/norm;
          }
          else if(side == FiberNetwork::FRONT)
          {
            norm = sqrt(dgdy*dgdy + (dgdz-1)*(dgdz-1) + dgdx*dgdx);
            normal[0] = dgdx/norm;
            normal[1] = dgdy/norm;
            normal[2] = (1-dgdz)/norm;
          }
          else if(side == FiberNetwork::BACK)
          {
            norm = sqrt(dgdy*dgdy + (dgdz-1)*(dgdz-1) + dgdx*dgdx);
            normal[0] = dgdx/norm;
            normal[1] = dgdy/norm;
            normal[2] = (dgdz-1)/norm;
          }
          double ndudx = normal[0]*dudx + normal[1]*dvdx + normal[2]*dgdx;
          double ndudy = normal[0]*dudy + normal[1]*dvdy + normal[2]*dgdy;
          double ndudz = normal[0]*dudz + normal[1]*dvdz + normal[2]*dgdz;
          fi += normal[0]*force_vector[3*node]*ndudx +
            normal[1]*force_vector[3*node]*ndudy +
            normal[2]*force_vector[3*node]*ndudz;
          fj += normal[0]*force_vector[3*node+1]*ndudx +
            normal[1]*force_vector[3*node+1]*ndudy +
            normal[2]*force_vector[3*node+1]*ndudz;
          fk += normal[0]*force_vector[3*node+2]*ndudx +
            normal[1]*force_vector[3*node+2]*ndudy +
            normal[2]*force_vector[3*node+2]*ndudz;
        }
        // CALCULATE THE MACRO-SCOPIC TERM OF THE Q
        tsfun(phi,phic,phie,phis,xyz[0],xyz[1],xyz[2]);
        calc_deriv(Fx,Fy,Fz,
                   ux,uy,uz,
                   phic,phie,phis,
                   &dudx,&dudy,&dudz,
                   &dvdx,&dvdy,&dvdz,
                   &dgdx,&dgdy,&dgdz,
                   &dett);
        double area = 0.0;
        if(side == FiberNetwork::TOP)
        {
          norm = sqrt(dvdx*dvdx + (dvdy-1)*(dvdy-1) + dvdz*dvdz);
          normal[0] = dvdx/norm;
          normal[1] = (dvdy-1)/norm;
          normal[2] = dvdz/norm;
          area = 0.5*((rve[1][0]-rve[2][0])*(rve[0][2]-rve[3][2])-
                      (rve[1][2]-rve[2][2])*(rve[0][0]-rve[3][0]));
        }
        else if(side == FiberNetwork::LEFT)
        {
          norm = sqrt(dvdx*dvdx + (dvdy-1)*(dvdy-1)+dvdz*dvdz);
          normal[0] = dvdx/norm;
          normal[1] = (1-dvdy)/norm;
          normal[2] = dvdz/norm;
          area = 0.5*((rve[5][0]-rve[6][0])*(rve[4][2]-rve[7][2])-
                      (rve[5][2]-rve[6][2])*(rve[0][4]-rve[7][0]));
        }
        else if(side == FiberNetwork::BOTTOM)
        {
          norm = sqrt(dudy*dudy + (dudx-1)*(dudx-1) + dudz*dudz);
          normal[0] = (1-dudx)/norm;
          normal[1] = dudy/norm;
          normal[2] = dudz/norm;
          area = 0.5*((rve[5][2]-rve[3][2])*(rve[7][1]-rve[1][1])-
                      (rve[5][1]-rve[3][1])*(rve[7][2]-rve[1][2]));
        }
        else if(side == FiberNetwork::RIGHT)
        {
          norm = sqrt(dudy*dudy + (dudx-1)*(dudx-1) + dudz*dudz);
          normal[0] = (dudx-1)/norm;
          normal[1] = dudy/norm;
          normal[2] = dudz/norm;
          area = 0.5*((rve[4][2]-rve[2][2])*(rve[6][1]-rve[0][1])-
                      (rve[4][1]-rve[2][1])*(rve[6][2]-rve[0][2]));
        }
        else if(side == FiberNetwork::FRONT)
        {
          norm = sqrt(dgdy*dgdy + (dgdz-1)*(dgdz-1) + dgdx*dgdx);
          normal[0] = dgdx/norm;
          normal[1] = dgdy/norm;
          normal[2] = (1-dgdz)/norm;
          area = 0.5*((rve[5][0]-rve[0][0])*(rve[4][1]-rve[1][1])-
                      (rve[5][1]-rve[0][1])*(rve[4][0]-rve[1][0]));
        }
        else if(side == FiberNetwork::BACK)
        {
          norm = sqrt(dgdy*dgdy + (dgdz-1)*(dgdz-1) + dgdx*dgdx);
          normal[0] = dgdx/norm;
          normal[1] = dgdy/norm;
          normal[2] = (dgdz-1)/norm;
          area = 0.5*((rve[7][0]-rve[2][0])*(rve[6][1]-rve[3][1])-
                      (rve[7][1]-rve[2][1])*(rve[6][0]-rve[3][0]));
        }
        dundx += (normal[0]*dudx + normal[1]*dvdx + normal[2]*dgdx) * area;
        dundy += (normal[0]*dudy + normal[1]*dvdy + normal[2]*dgdy) * area;
        dundz += (normal[0]*dudz + normal[1]*dvdz + normal[2]*dgdz) * area;
      }
      loc_vastrx = (fi - (stress[0]*dundx +
                          stress[1]*dundy +
                          stress[2]*dundz) / vol) / vol;
      loc_vastry = (fj - (stress[1]*dundx +
                          stress[3]*dundy +
                          stress[4]*dundz) / vol) / vol;
      loc_vastrz = (fk - (stress[2]*dundx +
                          stress[4]*dundy +
                          stress[5]*dundz) / vol) / vol;
    }
  */
}
