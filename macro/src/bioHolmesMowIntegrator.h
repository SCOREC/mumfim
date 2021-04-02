#ifndef BIO_HOLMESMOW_INTEGRATOR_H_
#define BIO_HOLMESMOW_INTEGRATOR_H_
#include "bioNonlinearTissue.h"
#include <ElementalSystem.h>
#include <apfShape.h>
#include <cstring>
#include <math.h> //natural log
// Holmes-Mow Constitutive relationship.
// Equations are based on http://help.mrl.sci.utah.edu/help/index.jsp (FEBIO documentation)
namespace bio
{
  class HolmesMowIntegrator : public amsi::ElementalSystem
  {
  public:
  HolmesMowIntegrator(NonlinearTissue * n,
                      apf::Field * field,
                      double shear_modulus,
                      double poisson_ratio,
                      int o)
    : ElementalSystem(field,o)
      , current_integration_point(0)
      , analysis(n)
      , dim(0)
    {
      ShearModulus  = shear_modulus;
      PoissonsRatio = poisson_ratio;
    }
    void inElement(apf::MeshElement * me)
    {
      ElementalSystem::inElement(me);
      fs = apf::getShape(f);
      es = fs->getEntityShape(apf::getMesh(f)->getType(apf::getMeshEntity(me)));
      dim = apf::getDimension(me);
    }
    void outElement()
    {
      current_integration_point = 0;
      ElementalSystem::outElement();
    }
    bool includesBodyForces() { return true; }
    void atPoint(apf::Vector3 const &p, double w, double)
    {
      int & nen = nenodes; // = 4 (tets)
      int & nedof = nedofs; // = 12 (tets)
      // Calculate current coordinates using accumulated displacements up to this point.
      // 1. Get coordinates on underlying mesh
      apf::Mesh * mesh = apf::getMesh(f);
      apf::Field * apf_coord_field = mesh->getCoordinateField();
      apf::Element * mesh_coord_elem = apf::createElement(apf_coord_field,me);
      apf::NewArray<apf::Vector3> mesh_xyz;
      apf::getVectorNodes(mesh_coord_elem,mesh_xyz);
      // 2. Get coordinates from apf_primary_field (passed in), which contains the accumulated displacement
      apf::NewArray<apf::Vector3> primary_field_xyz;
      apf::getVectorNodes(e,primary_field_xyz);
      // 3. Calculate current coordinates
      apf::DynamicMatrix xyz(nen,dim); xyz.zero();
      for (int ii = 0; ii < nen; ii++)
        for (int jj = 0; jj < dim; jj++)
          xyz(ii,jj) = mesh_xyz[ii][jj] + primary_field_xyz[ii][jj];
      // For Updated Lagrangian, the Jacobian of the current coordinate is used
      // Note: that entires of Jacobian is hard coded for Linear tetrahedra elements.
      // TO DO: Generalize Jacobian for current configuration.
      // Note: Form of Jacobian is
      // [dx/dr dx/ds dx/dt]    (x,y,z): physical domain coordinates.
      // [dy/dr dy/ds dy/dt]    (r,s,t): parent domain coordiantes.
      // [dz/dr dz/ds dz/dt]
      // This is equivalent to how Jacobian is calculated in apf.
      apf::Matrix<3,3> J, Jinv;
      J[0][0] = xyz(1,0) - xyz(0,0); // x2-x1
      J[1][0] = xyz(2,0) - xyz(0,0); // x3-x1
      J[2][0] = xyz(3,0) - xyz(0,0); // x4-x1
      J[0][1] = xyz(1,1) - xyz(0,1); // y2-y1
      J[1][1] = xyz(2,1) - xyz(0,1); // y3-y1
      J[2][1] = xyz(3,1) - xyz(0,1); // y4-y1
      J[0][2] = xyz(1,2) - xyz(0,2); // z2-z1
      J[1][2] = xyz(2,2) - xyz(0,2); // z3-z1
      J[2][2] = xyz(3,2) - xyz(0,2); // z4-z1
      double detJ = getDeterminant(J); // For Gauss-Quadrature Integration.
      Jinv = apf::invert(J);
      // Calculate derivative of shape fxns wrt current coordinates.
      apf::DynamicMatrix grads(nen,dim);
      grads.zero();
      apf::NewArray<apf::Vector3> local_grads;
      es->getLocalGradients(mesh, apf::getMeshEntity(me), p, local_grads);
      for (int ii = 0; ii < nen; ii++)
        for (int jj = 0; jj < dim; jj++)
          for (int kk = 0; kk < dim; kk++)
            grads(ii,jj) += Jinv[jj][kk] * local_grads[ii][kk];
      /*=======================================================
        Constitutive Component of Tangent Stiffness Matrix: K0
        - See Bathe PP.557, Table 6.6 (Updated Lagrangian)
        =======================================================
      */
      // hard-coded for 3d, make a general function... to produce this
      apf::DynamicMatrix BL(6,nedof); // linear strain disp
      BL.zero();
      for(int ii = 0; ii < nen; ii++)
      {
        BL(0,dim*ii  ) = grads(ii,0); // N_(ii,1)
        BL(1,dim*ii+1) = grads(ii,1); // N_(ii,2)
        BL(2,dim*ii+2) = grads(ii,2); // N_(ii,3)
        BL(3,dim*ii  ) = grads(ii,1); // N_(ii,2)
        BL(3,dim*ii+1) = grads(ii,0); // N_(ii,1)
        BL(4,dim*ii+1) = grads(ii,2); // N_(ii,3)
        BL(4,dim*ii+2) = grads(ii,1); // N_(ii,2)
        BL(5,dim*ii  ) = grads(ii,2); // N_(ii,3)
        BL(5,dim*ii+2) = grads(ii,0); // N_(ii,1)
      }
      /*
        Determine stress-strain matrix, D, for Holmes-Mow material.
        See http://help.mrl.sci.utah.edu/help/index.jsp
      */
      // Calculate Cauchy stress for Holmes-Mow material.
      apf::NewArray<apf::Vector3> grads0;
      apf::getShapeGrads(e,p,grads0); // derivative of shape function with respect to initial configuration.
      apf::Matrix<3,3> F; // Deformation gradient tensor
      apf::DynamicVector u(nedofs);
      getDisplacements(u);
      for (int ii = 0; ii < 3; ii++)
        for(int jj = 0; jj < 3; jj++)
        {
          F[ii][jj] = (ii == jj);
          for (int kk = 0; kk < nen; kk++)
            F[ii][jj] += grads0[kk][jj] * u(kk*3+ii);
//          F[ii][jj] += xyz(kk,ii) * grads0[kk][jj];
        }
      double detF = getDeterminant(F);
      apf::DynamicMatrix leftCauchyGreen(3,3); leftCauchyGreen.zero();
      apf::DynamicMatrix rightCauchyGreen(3,3); rightCauchyGreen.zero();
      apf::DynamicMatrix FT(3,3); FT.zero();
      apf::transpose(fromMatrix(F),FT);
      apf::multiply(fromMatrix(F),FT,leftCauchyGreen);
      apf::multiply(FT,fromMatrix(F),rightCauchyGreen);
      /** Determine Invariants I1 and I2 */
      double I1 = 0.0;
      double trAsq = 0.0;
      for (int ii = 0; ii < 3; ii++)
      {
        I1 += rightCauchyGreen(ii,ii);
        trAsq += ( rightCauchyGreen(ii,ii) * rightCauchyGreen(ii,ii) );
      }
      double I2 = 0.5 * ( I1 * I1 - trAsq );
      /** Determine Cauchy stress tensor
       * lambda and mu are Lame constants. see http://scienceworld.wolfram.com/physics/LameConstants.html.
       * lambda = 2*PoissonsRatio*ShearModulus/(1-2*PoissonsRatio)
       * mu = lambda * (1-2*PoissonsRatio)/(2*PoissonsRatio)
       */
      double lambda = ( 2.0 * PoissonsRatio * ShearModulus )/(1.0 - 2.0 * PoissonsRatio);
      double mu = lambda * (1.0 - 2.0 * PoissonsRatio)/(2.0 * PoissonsRatio);
      apf::DynamicMatrix Cauchy(3,3);
      double beta = 4.0;
      double Q = beta/(lambda + 2.0 * mu) * (
        ( 2.0 * mu - lambda ) * ( I1 - 3.0 )
        + lambda * ( I2 - 3.0 )
        - ( lambda + 2.0 * mu ) * std::log(detF*detF)
        );
      for (int ii = 0; ii < 3; ii++)
        for (int jj = 0; jj < 3; jj++)
        {
          if (ii == jj)
            Cauchy(ii,jj) = 0.5/detF * std::exp(Q) * (
              ( 2.0 * mu + lambda * (I1 - 1.0) ) * leftCauchyGreen(ii,jj)
              - lambda * leftCauchyGreen(ii,jj) * leftCauchyGreen(ii,jj)
              - ( lambda + 2 * mu )
              );
          else
            Cauchy(ii,jj) = 0.5/detF * std::exp(Q) * (
              ( 2.0 * mu + lambda * (I1 - 1.0) ) * leftCauchyGreen(ii,jj)
              - lambda * leftCauchyGreen(ii,jj) * leftCauchyGreen(ii,jj)
              );
        }
      /** Elasticity Tensor */
      apf::DynamicMatrix D(6,6); D.zero();
      double coeff1 = (4.0 * beta)/(lambda + 2.0 * mu) * detF/exp(Q);
      double coeff2 = exp(Q)/detF;
      apf::DynamicMatrix bxb(6,6); bxb.zero();
      apf::DynamicMatrix b_b(6,6); b_b.zero();
      apf::DynamicMatrix sxs(6,6); sxs.zero();
      apf::DynamicMatrix I_I(6,6); I_I.zero();
      apf::Matrix<3,3> IDmat;
      for (int ii = 0; ii < 3; ii++)
        for (int jj = 0; jj < 3; jj++)
        {
          if (ii==jj)
            IDmat[ii][jj] = 1.0;
          else
            IDmat[ii][jj] = 0.0;
        }
      for (int ii = 0; ii < 3; ii++)
        for (int jj = 0; jj < 3; jj++)
          for (int kk = 0; kk < 3; kk++)
            for (int ll = 0; ll < 3; ll++)
            {
              int a = ToVoigt(ii,jj);
              int b = ToVoigt(kk,ll);
              sxs(a,b) = Cauchy(ii,jj) * Cauchy(kk,ll);
              bxb(a,b) = leftCauchyGreen(ii,jj) * leftCauchyGreen(kk,ll);
              b_b(a,b) = 0.5 * (leftCauchyGreen(ii,kk) * leftCauchyGreen(jj,ll) + leftCauchyGreen(ii,ll) * leftCauchyGreen(jj,kk) );
              I_I(a,b) = 0.5 * (IDmat[ii][kk] * IDmat[jj][ll] + IDmat[ii][ll] * IDmat[jj][kk] );
            }
      apf::DynamicMatrix T1(6,6); T1.zero();
      apf::DynamicMatrix T2(6,6); T2.zero();
      apf::DynamicMatrix T3(6,6); T3.zero();
      /** Form D matrix:
       *  D = coeff1 * sxs + coeff2 * (lambda * (bxb - b_b) + (lambda + 2.0 * mu) * I_I) */
      T1 = sxs;
      T1 *= coeff1;
      T2 = bxb;
      T2 -= b_b;
      T2 *= lambda;
      T2 *= coeff2;
      T3 = I_I;
      T3 *= (lambda + 2.0 * mu);
      T3 *= coeff2;
      D += T1;
      D += T2;
      D += T3;
      apf::DynamicMatrix K0(nedof,nedof);
      apf::DynamicMatrix BLT(nedof,6);
      apf::DynamicMatrix DBL(6,nedof);
      apf::transpose(BL,BLT);
      apf::multiply(D,BL,DBL);
      apf::multiply(BLT,DBL,K0);
      /*=======================================================================
        Initial Stress (or Geometric) component of Tangent Stiffness Matrix: K1
        - See Bathe PP.557, Table 6.6 (Updated Lagrangian)
        ========================================================================*/
      apf::DynamicMatrix BNL(9,nedof); //nonlinear strain disp
      BNL.zero();
      double Bp[3][10] = {};
      for(int ii = 0; ii < 3; ii++)
        for(int jj = 0; jj < 4; jj++)
          Bp[ii][jj*3] = grads(jj,ii);
      for(int ii = 0; ii < 3; ii++)
        for(int jj = 0; jj < 10; jj++)
        {
          BNL(ii,jj) = Bp[ii][jj];
          BNL(3+ii,jj+1) = Bp[ii][jj];
          BNL(6+ii,jj+2) = Bp[ii][jj];
        }
      apf::DynamicMatrix tau(9,9);
      tau.zero();
      for (int ii = 0; ii < 3; ii++)
        for (int jj = 0; jj < 3; jj++)
          for (int kk = 0; kk < 3; kk++)
            tau(jj+ii*3,kk+ii*3) = Cauchy(jj,kk);
      apf::DynamicMatrix K1(nedof,nedof); K1.zero();
      apf::DynamicMatrix BNLT(nedof,9); BNLT.zero();
      apf::DynamicMatrix BNLTxtau(9,nedof); BNLTxtau.zero();
      apf::transpose(BNL,BNLT);
      apf::multiply(BNLT,tau,BNLTxtau);
      apf::multiply(BNLTxtau,BNL,K1);
      /*=======================================================================
        Terms for force vector (RHS)
        - See Bathe PP.557, Table 6.6 (Updated Lagrangian)
        ========================================================================*/
      double CauchyVoigt[6] = {};
      CauchyVoigt[0] = Cauchy(0,0);
      CauchyVoigt[1] = Cauchy(1,1);
      CauchyVoigt[2] = Cauchy(2,2);
      CauchyVoigt[3] = Cauchy(0,1);
      CauchyVoigt[4] = Cauchy(1,2);
      CauchyVoigt[5] = Cauchy(0,2);
      apf::DynamicMatrix BLTxCauchyVoigt(nedof,1); BLTxCauchyVoigt.zero();
      for (int ii = 0; ii < nedof; ii++)
        for (int jj = 0; jj < 6; jj++)
          BLTxCauchyVoigt(ii,0) += BLT(ii,jj) * CauchyVoigt[jj];
      for (int ii = 0; ii < nedof; ii++)
      {
        fe[ii] -= w * detJ * BLTxCauchyVoigt(ii,0);
        for (int jj = 0; jj < nedof; jj++)
          Ke(ii,jj) += w * detJ * (K0(ii,jj) + K1(ii,jj));
      }
      // store stress and strain values for post processing.
//      apf::DynamicVector u(nedofs);
//      getDisplacements(u);
/*
      std::cout<<"Ke:"<<std::endl;
      std::cout<<Ke<<std::endl;
      std::cout<<"fe:"<<std::endl;
      std::cout<<fe<<std::endl;
*/
      apf::DynamicVector strain(6);
      apf::multiply(BL,u,strain);
      analysis->storeStrain(me,strain.begin());
      analysis->storeStress(me,CauchyVoigt);
      current_integration_point++;
    }
    void getDisplacements(apf::DynamicVector & u)
    {
      apf::Element * e_disp = apf::createElement(f,me);
      apf::NewArray<apf::Vector3> disp;
      apf::getVectorNodes(e_disp,disp);
      for (int ii=0;ii<nedofs;ii++)
        u(ii) = disp[ii/3][ii%3];
    }
    int current_integration_point;
  private:
    NonlinearTissue * analysis;
    int dim;
    apf::FieldShape * fs;
    apf::EntityShape * es;
    double ShearModulus;
    double PoissonsRatio;
    int ToVoigt(int i,int j)
    {
      if (i==0 && j==0)
        return 0;
      else if (i==1 && j==1)
        return 1;
      else if (i==2 && j==2)
        return 2;
      else if (i==0 && j==1)
        return 3;
      else if (i==1 && j==0)
        return 3;
      else if (i==1 && j==2)
        return 4;
      else if (i==2 && j==1)
        return 4;
      else if (i==0 && j==2)
        return 5;
      else if (i==2 && j==0)
        return 5;
      else
      {
        std::cout<<"invalid indices in ToVoigt!"<<std::endl;
        return -1;
      }
    }
  };
}
#endif
