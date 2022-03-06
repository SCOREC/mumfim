#ifndef BIO_NEOHOOKEAN_INTEGRATOR_FAD_H_
#define BIO_NEOHOOKEAN_INTEGRATOR_FAD_H_
#include "bioNonlinearTissue.h"
#include <ElementalSystem.h>
#include <apfShape.h>
#include <apfSIM.h>
#include <cstring>
#include <math.h> //natural log
#include <mthAD.h> //Auto Differentiation Library
#include <mthVector.h>
#include <mthMatrix.h>
namespace bio
{
  class NeoHookeanIntegratorFAD : public amsi::ElementalSystem
  {
  public:
  NeoHookeanIntegratorFAD(NonlinearTissue * n,
                          apf::Field * field,
                          apf::Field * strain_ip_field,
                          apf::Field * stress_ip_field,
                          double shear_modulus,
                          double poisson_ratio,
                          apf::Field * type_field,
                          int o)
    : ElementalSystem(field,o)
      , current_integration_point(0)
      , analysis(n)
      , dim(0)
      , strain_field(strain_ip_field)
      , stress_field(stress_ip_field)
      , micro_type_field(type_field)
    {
    // Parameters from V. Lai et al. Journal of Biomechanical Engineering, Vol 135, 071007 (2013)
      ShearModulus  = shear_modulus;// kPa
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
    void atPoint(apf::Vector3 const &p, double w, double dV)
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
      // Note: that entrees of Jacobian is hard coded for Linear tetrahedra elements.
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
      /// Calculate derivative of shape fxns wrt current coordinates.
      apf::DynamicMatrix grads(nen,dim); //<for Linear Tetrahedral Elements, nen = 4.
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
        Determine stress-strain matrix, D, for NeoHookean material.
        See Bonet and Wood 2nd Edition, PP.250.
        - lambda and mu are effective lame parameters.
        - lambda = ShearModulus/detJ, mu = ( PoissonsRatio - ShearModulus * ln(detJ) )/detJ
       */
      // Calculate Cauchy stress from compressible NeoHookean equation. See Bonet and Wood 2nd Ed. PP.163.
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
      /// Test Auto Differentiation
      mth::Vector<mth::AD<double>> uAD;
      uAD.resize(nedofs);
      for (int ii=0; ii < nedofs; ii++)
      {
        uAD[ii] = u[ii];
        uAD[ii].diff(ii,nedofs);
        std::cout<<"val of dof "<<ii<<" is "<<uAD[ii].val()<<std::endl;
        std::cout<<"disp of dof "<<ii<<" is "<<u[ii]<<std::endl;
        for (int jj=0; jj < nedofs;jj++)
          std::cout<<jj<<"th dx of dof "<<ii<<" is "<<uAD[ii].dx(jj)<<std::endl;
      }
      mth::Matrix3x3<mth::AD<double>> FAD;
      mth::Matrix3x3<mth::AD<double>> leftCauchyGreenFAD;
      FAD.zero(); leftCauchyGreenFAD.zero();
      for (int ii = 0; ii < 3; ii++)
        for(int jj = 0; jj < 3; jj++)
        {
          FAD(ii,jj) = (ii == jj);
          for (int kk = 0; kk < nen; kk++)
            FAD(ii,jj) += grads0[kk][jj] * uAD(kk*3+ii);
          leftCauchyGreenFAD(ii,jj) = FAD(ii,jj) * FAD(jj,ii);
        }
      for (int ii = 0; ii < 3; ii++)
        for (int jj = 0; jj < 3; jj++)
        {
          std::cout<<"value of FAD "<<ii<<","<<jj<<"="<<FAD(ii,jj).val()<<std::endl;
          std::cout<<"value of leftCauchyGreen "<<ii<<","<<jj<<"="<<leftCauchyGreenFAD(ii,jj).val()<<std::endl;
          for (int kk = 0; kk < nedofs; kk++)
          {
            std::cout<<kk<<"th derivative of FAD"<<ii<<","<<jj<<"="<<FAD(ii,jj).dx(kk)<<std::endl;
            std::cout<<kk<<"th derivative of leftCauchyGreen"<<ii<<","<<jj<<"="<<leftCauchyGreenFAD(ii,jj).dx(kk)<<std::endl;
          }
        }
      apf::DynamicMatrix leftCauchyGreen(3,3); leftCauchyGreen.zero();
      apf::DynamicMatrix FT(3,3); FT.zero();
      apf::transpose(fromMatrix(F),FT);
      apf::multiply(fromMatrix(F),FT,leftCauchyGreen);
      apf::DynamicMatrix Cauchy(3,3);
      double mu = ShearModulus;
      double lambda = ( 2.0 * mu * PoissonsRatio )/( 1.0 - 2.0 * PoissonsRatio );
      for (int ii = 0; ii < 3; ii++)
        for (int jj = 0; jj < 3; jj++)
        {
          if (ii == jj)
            Cauchy(ii,jj) = 1.0/detF * ( mu * (leftCauchyGreen(ii,jj) - 1.0) + lambda  * log(detF) );
          else
            Cauchy(ii,jj) = mu/detF * leftCauchyGreen(ii,jj);
        }
      mth::Matrix3x3<mth::AD<double>> CauchyFAD; CauchyFAD.zero();
      for (int ii = 0; ii < 3; ii++)
        for (int jj = 0; jj < 3; jj++)
        {
          if (ii == jj)
            CauchyFAD(ii,jj) = 1.0/detF * ( mu * (leftCauchyGreenFAD(ii,jj) - 1.0) + lambda  * log(detF) );
          else
            CauchyFAD(ii,jj) = mu/detF * leftCauchyGreenFAD(ii,jj);
        }
      for (int ii = 0; ii < 3; ii++)
        for (int jj = 0; jj < 3; jj++)
        {
          std::cout<<"value of Cauchy "<<ii<<","<<jj<<"="<<CauchyFAD(ii,jj).val()<<std::endl;
          for (int kk = 0; kk < nedofs; kk++)
            std::cout<<kk<<"th derivative of "<<ii<<","<<jj<<"="<<CauchyFAD(ii,jj).dx(kk)<<std::endl;
        }
      double lambda_prime = lambda/detF;
      double mu_prime = ( mu - lambda * log(detF) )/detF;
      apf::DynamicMatrix D(6,6);
      D.zero();
      D(0,0) = lambda_prime + (2.0 * mu_prime);
      D(0,1) = lambda_prime;
      D(0,2) = lambda_prime;
      D(1,0) = lambda_prime;
      D(1,1) = lambda_prime + (2.0 * mu_prime);
      D(1,2) = lambda_prime;
      D(2,0) = lambda_prime;
      D(2,1) = lambda_prime;
      D(2,2) = lambda_prime + (2.0 * mu_prime);
      D(3,3) = mu_prime;
      D(4,4) = mu_prime;
      D(5,5) = mu_prime;
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
      // Calculate internal elemental forces via Auto Differentiation
      mth::Vector<mth::AD<double>> feAD;
      feAD.resize(nedofs);
      for (int ii=0; ii<nedofs; ii++)
        feAD(ii) = feAD(ii) * 0.0; //<zero out vector.
      std::cout<<"zero out feAD:"<<feAD<<std::endl;
      for (int a=0; a<nen; a++)
      {
        feAD(a * 3 + 0) += w * detJ * ( CauchyFAD(0,0) * grads(a,0) + CauchyFAD(0,1) * grads(a,1) + CauchyFAD(0,2) * grads(a,2) );
        feAD(a * 3 + 1) += w * detJ * ( CauchyFAD(1,0) * grads(a,0) + CauchyFAD(1,1) * grads(a,1) + CauchyFAD(1,2) * grads(a,2) );
        feAD(a * 3 + 2) += w * detJ * ( CauchyFAD(2,0) * grads(a,0) + CauchyFAD(2,1) * grads(a,1) + CauchyFAD(2,2) * grads(a,2) );
      }
      apf::DynamicMatrix KefromfeAD(nedofs,nedofs);
      KefromfeAD.zero();
      for (int ii=0; ii<nedofs; ii++)
        for (int jj=0; jj<nedofs; jj++)
          KefromfeAD(ii,jj) = feAD(ii).dx(jj);
      std::cout<<"feAD ="<<std::endl;
      std::cout<<feAD<<std::endl;
      apf::DynamicMatrix BLTxCauchyVoigt(nedof,1); BLTxCauchyVoigt.zero();
      for (int ii = 0; ii < nedof; ii++)
        for (int jj = 0; jj < 6; jj++)
          BLTxCauchyVoigt(ii,0) += BLT(ii,jj) * CauchyVoigt[jj];
      for (int ii = 0; ii < nedof; ii++)
      {
        fe[ii] -= w * detJ * BLTxCauchyVoigt(ii,0);
//      fe[ii] = -1.0 * feAD(ii).val();
        for (int jj = 0; jj < nedof; jj++)
        {
          Ke(ii,jj) = w * detJ * (K0(ii,jj) + K1(ii,jj));
//        Ke(ii,jj) = feAD(ii).dx(jj);
        }
      }
      std::cout<<"fe = "<<std::endl;
      std::cout<<fe<<std::endl;
      std::cout<<"Ke = "<<std::endl;
      std::cout<<Ke<<std::endl;
      std::cout<<"KefromfeAD = "<<std::endl;
      std::cout<<KefromfeAD<<std::endl;
      // store stress and strain values for post processing.
//      apf::DynamicVector u(nedofs);
//      getDisplacements(u);
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
    apf::Field * strain_field;
    apf::Field * stress_field;
    apf::Field * micro_type_field;
    double ShearModulus;
    double PoissonsRatio;
  };
}
#endif