#ifndef BIO_MULTISCALE_INTEGRATOR_H_
#define BIO_MULTISCALE_INTEGRATOR_H_
#include "bioMultiscaleTissue.h"
#include <ElementalSystem.h>
#include <apfShape.h>
#include <apfSIM.h>
#include <cstring>
namespace bio
{
  class MultiscaleIntegrator : public amsi::ElementalSystem
  {
  public:
    MultiscaleIntegrator(MultiscaleTissue * n,
                         apf::Field * field,
                         apf::Field * type_field,
                         int o)
    : ElementalSystem(field,o)
      , current_integration_point(0)
      , analysis(n)
      , micro_type_field(type_field)
  {
    matrixShearModulus  = 0.0;
    matrixPoissonsRatio = 0.45;
    matrixBulkModulus = 2.0*matrixShearModulus*(1+matrixPoissonsRatio)/(3.0*(1.0-2.0*matrixPoissonsRatio));
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
      int microscaleType = apf::getScalar(micro_type_field,apf::getMeshEntity(me),current_integration_point);
      // Choose appropraite atPoint function based on microscale type
      if(microscaleType == MultiscaleTissue::NONE)
        atPointLinearElastic(p,w,dV);
      else if(microscaleType == MultiscaleTissue::FIBER_ONLY)
      {
        atPointUseMicro(p,w,dV);
        if(matrixShearModulus!=0.0)  // matrixShearModulus = 0 corresponds to "fiber only"
          atPointMatrix(p,w,dV);
      }
      else if(microscaleType == MultiscaleTissue::FIBER_MATRIX)
        atPointUseMicro(p,w,dV);
    }
    void atPointUseMicro(apf::Vector3 const &p, double w, double dV)
    {
        RVE_Info * rve_info = analysis->getRVEResult(apf::getMeshEntity(me),
                                                              current_integration_point);
        // Note: Macro and micro solvers index stress tensor differently
        //
        //   micro    ----->    macro
        //
        // [ 0 1 2 ]          [ 0 3 5 ]
        // [   3 4 ]  ----->  [   1 4 ]
        // [     5 ]          [     2 ]
        //
        apf::Matrix3x3 Jac;
        apf::getJacobian(me,p,Jac);
        // Note: Form of Jacobian in apf is
        // [dx/dr dx/ds dx/dt]    (x,y,z): physical domain coordinates.
        // [dy/dr dy/ds dy/dt]    (r,s,t): parent domain coordiantes.
        // [dz/dr dz/ds dz/dt]
        double wxdetjac = w * apf::getJacobianDeterminant(Jac,dim);
        int & nen = nenodes; // = 4 (tets)
        int & nedof = nedofs; // = 12 (tets)
        int offset = 9;
        double * stress_deriv[6];
        stress_deriv[0] = &(rve_info->derivS[offset]);
        stress_deriv[1] = &(rve_info->derivS[offset + 3*nedof]);
        stress_deriv[2] = &(rve_info->derivS[offset + 5*nedof]);
        stress_deriv[3] = &(rve_info->derivS[offset +   nedof]);
        stress_deriv[4] = &(rve_info->derivS[offset + 4*nedof]);
        stress_deriv[5] = &(rve_info->derivS[offset + 2*nedof]);
        apf::NewArray<apf::Vector3> grads;
        apf::getShapeGrads(e,p,grads);
        // hard-coded for 3d, make a general function... to produce this
        apf::DynamicMatrix BL(6,nedof); // linear strain disp
        BL.zero();
        for(int ii = 0; ii < nen; ii++)
        {
          BL(0,dim*ii  ) = grads[ii][0]; // N_(ii,1)
          BL(1,dim*ii+1) = grads[ii][1]; // N_(ii,2)
          BL(2,dim*ii+2) = grads[ii][2]; // N_(ii,3)
          BL(3,dim*ii  ) = grads[ii][1]; // N_(ii,2)
          BL(3,dim*ii+1) = grads[ii][0]; // N_(ii,1)
          BL(4,dim*ii+1) = grads[ii][2]; // N_(ii,3)
          BL(4,dim*ii+2) = grads[ii][1]; // N_(ii,2)
          BL(5,dim*ii  ) = grads[ii][2]; // N_(ii,3)
          BL(5,dim*ii+2) = grads[ii][0]; // N_(ii,1)
        }
        apf::DynamicMatrix K0(nedof,nedof);
        K0.zero();
        // K0 = BL^T * stress_deriv
        for(int ii = 0; ii < nedof; ii++)
          for(int jj = 0; jj < nedof; jj++)
            for(int kk = 0; kk < 6; kk++)
              K0(ii,jj) += BL(kk,ii) * stress_deriv[kk][jj];
        apf::DynamicMatrix BNL(9,nedof); //nonlinear strain disp
        BNL.zero();
        double Bp[3][10] = {};
        for(int ii = 0; ii < 3; ii++)
          for(int jj = 0; jj < 4; jj++)
            Bp[ii][jj*3] = grads[jj][ii];
        for(int ii = 0; ii < 3; ii++)
          for(int jj = 0; jj < 10; jj++)
          {
            BNL(ii,jj) = Bp[ii][jj];
            BNL(3+ii,jj+1) = Bp[ii][jj];
            BNL(6+ii,jj+2) = Bp[ii][jj];
          }
        // Fill S matrix - in block notation contains stress tensor along diagonal, zero elsewhere
        double S[9][9] = {};
        double SV[6];
        double factor = 1;
        SV[0] = factor*(rve_info->derivS[0]);
        SV[1] = factor*(rve_info->derivS[3]);
        SV[2] = factor*(rve_info->derivS[5]);
        SV[3] = factor*(rve_info->derivS[1]);
        SV[4] = factor*(rve_info->derivS[4]);
        SV[5] = factor*(rve_info->derivS[2]);
        S[0][0] = S[0+3][0+3] = S[0+6][0+6] = rve_info->derivS[0];
        S[0][1] = S[0+3][1+3] = S[0+6][1+6] = rve_info->derivS[1];
        S[0][2] = S[0+3][2+3] = S[0+6][2+6] = rve_info->derivS[2];
        S[1][0] = S[1+3][0+3] = S[1+6][0+6] = rve_info->derivS[1];
        S[1][1] = S[1+3][1+3] = S[1+6][1+6] = rve_info->derivS[3];
        S[1][2] = S[1+3][2+3] = S[1+6][2+6] = rve_info->derivS[4];
        S[2][0] = S[2+3][0+3] = S[2+6][0+6] = rve_info->derivS[2];
        S[2][1] = S[2+3][1+3] = S[2+6][1+6] = rve_info->derivS[4];
        S[2][2] = S[2+3][2+3] = S[2+6][2+6] = rve_info->derivS[5];
        apf::DynamicMatrix BNLTxS(nedof,9);
        BNLTxS.zero();
        // BNLTxS = BNL^T * S
        for(int ii = 0; ii < nedof; ii++)
          for(int jj = 0; jj < 9; jj++)
            for(int kk = 0; kk < 9; kk++)
              BNLTxS(ii,jj) += BNL(kk,ii) * S[kk][jj];
        // for force vector calculation
        apf::DynamicMatrix BLTxSV(nedof,1);
        BLTxSV.zero();
        for(int ii = 0; ii < nedof; ii++)
          for(int jj = 0; jj < 6; jj++)
            BLTxSV(ii,0) += BL(jj,ii) * SV[jj];
        // retrieve virtual strain/stress for force vector calc
        double Q[3];
        Q[0] = rve_info->derivS[6];
        Q[1] = rve_info->derivS[7];
        Q[2] = rve_info->derivS[8];
        apf::DynamicMatrix N(num_field_components,nedof);
        N.zero();
        for(int ii = 0; ii < num_field_components; ii++)
          for(int jj = 0; jj < nen; jj++)
            N(ii,ii+(jj*num_field_components)) = p[ii];
            // assumption for linear tet... should be all shape function values...
        apf::DynamicMatrix NTxQ(nedof,1);
        NTxQ.zero();
        for(int ii = 0; ii < nedof; ii++)
          for(int jj = 0; jj < num_field_components; jj++)
            NTxQ(ii,0) += N(jj,ii) * Q[jj];
        apf::DynamicMatrix K1(nedof,nedof);
        K1.zero();
        apf::multiply(BNLTxS,BNL,K1);
        for(int ii = 0; ii < nedof; ii++)
        {
          fe[ii] += wxdetjac * (NTxQ(ii,0) - BLTxSV(ii,0)); // P - F
          for(int jj = 0; jj < nedof; jj++)
            Ke(ii,jj) += wxdetjac * (K0(ii,jj) + K1(ii,jj));
        }
      // store stress and strain values for post processing.
      // Get displacements
      apf::DynamicVector u(nedofs);
      getDisplacements(u);
      // Compute strains and stresses (for force vector)
      apf::DynamicVector strain(6);
      apf::multiply(BL,u,strain);
      analysis->storeStrain(me,strain.begin());
      analysis->storeStress(me,SV);
      current_integration_point++;
    }
    void atPointMatrix(apf::Vector3 const &p, double w, double dV)
    {
      apf::Matrix3x3 Jac;
      apf::getJacobian(me,p,Jac);
      int dim = apf::getDimension(me);
      double wxdetjac = w * apf::getJacobianDeterminant(Jac,dim);
      apf::NewArray<apf::Vector3> grads;
      apf::getShapeGrads(e,p,grads);
      // Get displacements at integration point
      apf::DynamicVector u(nedofs);
      getDisplacements(u);
      // Calculate deformation gradient tensor, F = d_ij + u_i,j
      apf::Matrix3x3 F;
      for(int ii=0;ii<3;ii++)
        for(int jj=0;jj<3;jj++)
        {
          F[ii][jj] = (ii==jj);
          for(int kk=0;kk<nenodes;kk++)
            F[ii][jj] += grads[kk][jj] * u(kk*3 + ii);
        }
      // Calculate J, the determinant of F
      double J = apf::getDeterminant(F);
      // Calculate left Cauchy-Green deformation tensor, B = F * F^T
      apf::DynamicMatrix B(3,3);
      for(int ii = 0; ii < 3; ii++)
        for(int jj = 0; jj < 3; jj++)
        {
          B(ii,jj) = 0.0;
          for(int kk=0;kk<3;kk++)
            B(ii,jj) += F[ii][kk] * F[jj][kk];
        }
      double Bqq = B(0,0)+B(1,1)+B(2,2);
      double oneThird = 1.0/3.0;
      // Calculate stresses from neo-Hookean law
      double GJ = matrixShearModulus/pow(J,5.0*oneThird);
      apf::DynamicMatrix cauchy(3,3);
      cauchy(0,0) = GJ * B(0,0);
      cauchy(0,1) = GJ * B(0,1);
      cauchy(0,2) = GJ * B(0,2);
      cauchy(1,0) = cauchy(0,1);
      cauchy(1,1) = GJ * B(1,1);
      cauchy(1,2) = GJ * B(1,2);
      cauchy(2,0) = cauchy(0,2);
      cauchy(2,1) = cauchy(1,2);
      cauchy(2,2) = GJ * B(2,2);
      for(int ii=0;ii<3;ii++)
      {
        //maxtrixStress(ii) -= matrixShearModulus*pow(J,-2.0*(matrixPoissonsRatio/(1.0-2.0*matrixPoissonsRatio))-1.0);
        //maxtrixStress(ii) += 2.0*GJ*matrixPoissonsRatio/(1.0-2.0*matrixPoissonsRatio) * log(J) - GJ;
        cauchy(ii,ii) += matrixBulkModulus*(J-1) - GJ*oneThird * (Bqq);
      }
      // Store matrix stress
      /*
      double stress[6];
      stress[0] = cauchy(0,0);
      stress[1] = cauchy(1,1);
      stress[2] = cauchy(2,2);
      stress[3] = cauchy(0,1);
      stress[4] = cauchy(1,2);
      stress[5] = cauchy(0,2);
      //analysis->storeMatrixStress(me,stress);
      */
      double dl[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
      // Calculate C tensor
      double C[3][3][3][3];
      for(int ii=0;ii<3;ii++)
        for(int jj=0;jj<3;jj++)
          for(int kk=0;kk<3;kk++)
            for(int ll=0;ll<3;ll++)
              C[ii][jj][kk][ll] = ( dl[ii][kk]*B(jj,ll) + B(ii,ll)*dl[jj][kk] -
                                    2.0*oneThird * (B(ii,jj)*dl[kk][ll]+B(kk,ll)*dl[ii][jj]) +
                                    2.0*oneThird*oneThird * dl[ii][jj]*dl[kk][ll] * (Bqq)
                ) * matrixShearModulus/pow(J,2.0*oneThird) +
                matrixBulkModulus*(2.0*J-1.0)*J*dl[ii][jj]*dl[kk][ll];
      // Fill force vector
      for(int aa=0;aa<nenodes;aa++)
        for(int ii=0;ii<3;ii++)
        {
          int row = 3*aa + ii;
          for(int jj=0;jj<3;jj++)
            fe[row] -= wxdetjac * cauchy(ii,jj)*grads[aa][jj];
        }
      // Fill stiffness matrix
      for(int aa=0;aa<nenodes;aa++)
        for(int ii=0;ii<3;ii++)
          for(int bb=0;bb<nenodes;bb++)
            for(int kk=0;kk<3;kk++)
            {
              int row = 3*aa + ii;
              int col = 3*bb + kk;
              for(int jj=0;jj<3;jj++)
              {
                for(int ll=0;ll<3;ll++)
                  Ke(row,col) += wxdetjac * C[ii][jj][kk][ll] * grads[bb][ll] * grads[aa][jj];
                Ke(row,col) -= wxdetjac * cauchy(ii,jj) * grads[aa][kk] * grads[bb][jj];
              }
            }
    }
    void atPointLinearElastic(apf::Vector3 const &p, double w, double dV)
    {
      apf::Matrix3x3 Jac;
      apf::getJacobian(me,p,Jac);
      double wxdetjac = w * apf::getJacobianDeterminant(Jac,apf::getDimension(me));
      int & nen = nenodes; // = 4 (tets)
      int & nedof = nedofs; // = 12 (tets)
      apf::NewArray<apf::Vector3> grads;
      apf::getShapeGrads(e,p,grads);
      // Linear strain-displacement matrix see Bathe pgs 555-556
      apf::DynamicMatrix BL(6,nedof); // linear strain disp
      for(int ii = 0; ii < nen; ii++)
      {
        BL(0,3*ii) = grads[ii][0]; // N_(ii,1)
        BL(0,3*ii+1) = BL(0,3*ii+2) = 0.0;
        BL(1,3*ii+1) = grads[ii][1]; // N_(ii,2)
        BL(1,3*ii) = BL(1,3*ii+2) = 0.0;
        BL(2,3*ii+2) = grads[ii][2]; // N_(ii,3)
        BL(2,3*ii) = BL(2,3*ii+1) = 0.0;
        BL(3,3*ii) = grads[ii][1]; // N_(ii,2)
        BL(3,3*ii+1) = grads[ii][0]; // N_(ii,1)
        BL(3,3*ii+2) = 0.0;
        BL(4,3*ii) = 0.0;
        BL(4,3*ii+1) = grads[ii][2]; // N_(ii,3)
        BL(4,3*ii+2) = grads[ii][1]; // N_(ii,2)
        BL(5,3*ii) = grads[ii][2]; // N_(ii,3)
        BL(5,3*ii+1) = 0.0;
        BL(5,3*ii+2) = grads[ii][0];  // N_(ii,1)
      }
      /*
        Determine stress-strain matrix, D. See Bathe pg 195.
        This implementation is identical to GetIsotropicStressStrainTensor function in
        LinearElasticIntegrator.cc
       */
      apf::DynamicMatrix D(6,6);
      double E = 0.005868; // Young's modulus, not right AT ALL
      double v = 0.236249; // Poisson's ratio
      double lambda = ( v * E ) / ( ( 1 + v ) * ( 1 - 2 * v ) );
      double mu = E / ( 2 * ( 1 + v ) );
      D(0,0) = lambda + (2 * mu);
      D(0,1) = lambda;
      D(0,2) = lambda;
      D(0,3) = D(0,4) = D(0,5) = 0.0;
      D(1,0) = lambda;
      D(1,1) = lambda + (2 * mu);
      D(1,2) = lambda;
      D(1,3) = D(1,4) = D(1,5) = 0.0;
      D(2,0) = lambda;
      D(2,1) = lambda;
      D(2,2) = lambda + (2 * mu);
      D(2,3) = D(2,4) = D(2,5) = 0.0;
      D(3,0) = D(3,1) = D(3,2) = D(3,4) = D(3,5) = 0.0;
      D(3,3) = mu;
      D(4,0) = D(4,1) = D(4,2) = D(4,3) = D(4,5) = 0.0;
      D(4,4) = mu;
      D(5,0) = D(5,1) = D(5,2) = D(5,3) = D(5,4) = 0.0;
      D(5,5) = mu;
      /*
        Numerical integration. Identical format to LinearElasticIntegrator::atPoint function in
        LinearElasticIntegrator.cc
       */
      apf::DynamicMatrix K0(nedof,nedof);
      apf::DynamicMatrix DBL(6,nedofs);
      apf::multiply(D,BL,DBL);
      apf::DynamicMatrix BLT(nedofs,6);
      apf::transpose(BL,BLT);
      apf::multiply(BLT,DBL,K0);
      /*
        Determine body forces
       */
      // Get displacements
      apf::DynamicVector u(nedofs);
      getDisplacements(u);
      // Compute strains and stresses (for force vector)
      apf::DynamicVector strain(6);
      apf::multiply(BL,u,strain);
      apf::DynamicVector stress(6);
      apf::multiply(D,strain,stress);
      // store stress and strain values for post processing.
      analysis->storeStrain(me,strain.begin());
      analysis->storeStress(me,stress.begin());
      // for force vector calculation
      apf::DynamicVector BLTxSV(nedof);
      apf::multiply(BLT,stress,BLTxSV);
      for(int ii = 0; ii < nedof; ii++)
      {
        fe[ii] += wxdetjac * (-BLTxSV(ii)); // P - F
          for(int jj = 0; jj < nedof; jj++)
            Ke(ii,jj) += wxdetjac * K0(ii,jj);
      }
      current_integration_point++;
    }
    void getDisplacements(apf::DynamicVector & u)
    {
      apf::Element * e_disp = apf::createElement(f,me);
      apf::NewArray<apf::Vector3> disp;
      apf::getVectorNodes(e_disp,disp);
      for(int ii=0;ii<nedofs;ii++)
        u(ii) = disp[ii/3][ii%3];
    }
    int current_integration_point;
  private:
    MultiscaleTissue * analysis;
    int dim;
    apf::FieldShape * fs;
    apf::EntityShape * es;
    apf::Field * micro_type_field;
    double linearYoungsModulus;
    double linearPoissonsRatio;
    double matrixShearModulus;
    double matrixPoissonsRatio;
    double matrixBulkModulus;
  };
}
#endif
