#ifndef BIO_TRANSVERSELY_ISOTROPIC_NEOHOOKEAN_INTEGRATOR_H_
#define BIO_TRANSVERSELY_ISOTROPIC_NEOHOOKEAN_INTEGRATOR_H_
#include "bioNonlinearTissue.h"
#include <ElementalSystem.h>
#include <apfShape.h>
#include <apfSIM.h>
#include <cstring>
#include <math.h> //natural log
namespace bio
{
  // 01/23/17: Transversely Isotropic NeoHookean Constitutive relation based on J. Bonet and A. J. Burton Comput. Methods Appl. Mech. Engrg. 162 (1998) 151-164.
  class TrnsIsoNeoHookeanIntegrator : public amsi::ElementalSystem
  {
  public:
  TrnsIsoNeoHookeanIntegrator(NonlinearTissue * n,
			      apf::Field * field,
			      apf::Field * stf_vrtn,
			      apf::Field * det_dfm_grd,
			      apf::Field * axl_yngs_mod,
			      double youngs_modulus,
			      double poisson_ratio,
			      double * axis,
			      double axial_shear_modulus,
			      double axial_youngs_modulus,
			      int o)

    : ElementalSystem(field,o)
      , current_integration_point(0)
      , stf_vrtn_fld(stf_vrtn)
      , EA_fld(axl_yngs_mod)
      , dfm_grd_fld(det_dfm_grd)
      , analysis(n)
      , dim(0)
      , ShearModulus(0.0)
      , PoissonsRatio(poisson_ratio)
      , AxialShearModulus(axial_shear_modulus)
      , AxialYoungsModulus(axial_youngs_modulus)
    {
      ShearModulus = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
      Axis[0] = axis[0];
      Axis[1] = axis[1];
      Axis[2] = axis[2];
    }
    void inElement(apf::MeshElement * me)
    {
      ElementalSystem::inElement(me);
      fs = apf::getShape(f);
      msh = apf::getMesh(f);
      es = fs->getEntityShape(msh->getType(apf::getMeshEntity(me)));
      dim = apf::getDimension(me);
      ce = apf::createElement(msh->getCoordinateField(),me);
      apf::getVectorNodes(ce,msh_xyz);
      apf::getVectorNodes(e,fld_xyz);
      current_integration_point = 0;
    }
    void outElement()
    {
      apf::destroyElement(ce);
      ElementalSystem::outElement();
    }
    bool includesBodyForces() { return true; }
    void atPoint(apf::Vector3 const &p, double w, double dV)
    {
      int & nen = nenodes; // = 4 (tets)
      int & nedof = nedofs; // = 12 (tets)
      // 3. Calculate current coordinates
      apf::DynamicMatrix xyz(nen,dim); xyz.zero();
      for (int ii = 0; ii < nen; ii++)
        for (int jj = 0; jj < dim; jj++)
          xyz(ii,jj) = msh_xyz[ii][jj] + fld_xyz[ii][jj];
      // For Updated Lagrangian, the Jacobian of the current coordinate is used
      // Note: entires of Jacobian is hard coded for Linear tetrahedra elements.
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
      es->getLocalGradients(msh, apf::getMeshEntity(me), p, local_grads);
      for (int ii = 0; ii < nen; ii++)
        for (int jj = 0; jj < dim; jj++)
          for (int kk = 0; kk < dim; kk++)
            grads(ii,jj) += Jinv[jj][kk] * local_grads[ii][kk];
      /*=======================================================
        Constitutive Component of Tangent Stiffness Matrix: K0
        - See Bathe PP.557, Table 6.6 (Updated Lagrangian)
        - Based on following Voigt format: [11,22,33,12,23,13]
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
        Calculate Deformation Tensors
        - F = deformation gradient tensor
        - detF = determinant of F
        - leftCauchyGreen = F*F^T
        - rightCauchyGreen = F^T*F
       */
      apf::NewArray<apf::Vector3> grads0;
      apf::getShapeGrads(e,p,grads0); // derivative of shape function with respect to initial configuration.
      apf::Matrix3x3 F(0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0); // Deformation gradient tensor
      apf::DynamicVector u(nedofs);
      getDisplacements(u);
      for (int ii = 0; ii < 3; ++ii)
        for(int jj = 0; jj < 3; ++jj)
        {
          F[ii][jj] = (ii == jj);
          for (int kk = 0; kk < nen; kk++)
            F[ii][jj] += grads0[kk][jj] * u(kk*3+ii);
//          F[ii][jj] += xyz(kk,ii) * grads0[kk][jj];
        }
      double detF = getDeterminant(F);
      apf::setScalar(dfm_grd_fld, apf::getMeshEntity(me), current_integration_point, detF);
      apf::Matrix3x3 leftCauchyGreen(0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0);
      apf::Matrix3x3 rightCauchyGreen(0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0);
      apf::Matrix3x3 FT(0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0);
      FT = transpose(F);
      leftCauchyGreen = F*FT;
      rightCauchyGreen = FT*F;
      // Specify Axial direction
      apf::Vector3 A(&Axis[0]);
      apf::Vector3 a(0.0,0.0,0.0);
      a = F*A;
      /*
        Calculate Cauchy stress from compressible NeoHookean equation. See Bonet and Wood 2nd Ed. PP.163.
         - sigma = sigma_nh + sigma_trns (see J. Bonet and A.J. Burton Comput. Methods Appl. Mech. Engrg. 162 (1998) 151-164)
         - sigma_nh_ij   = mu/J*(b_ij - I_ij) + lambda*(J-1)*I_ij
         - sigma_trns_ij = 2*beta/J*(a_r*a_r - 1.0)*I_ij + 2/J[alpha + beta*(C_rs*I_rs-3.0) - 2*gamma*(a_r*a_r - 1.0)]a_i*a_j
                         - alpha/J*(b_is*a_s*a_j+a_i*b_jr*a_r)
         - C_rs: elements of right Cauchy-Green deformation tensor (rightCauchyGreen).
         - b_is: elements of left Cauchy-Green deformation tensor (leftCauchyGreen).
         - J   : determinant of deformation gradient (detF).
         - Material Constants:
           + mu = Shear Modulus, nu = Poisson Ratio, G_A = axial shear modulus, E_A = axial Young's modulus.
           + lambda = 2*mu*(nu+n*nu^2)/m
           + alpha = mu - G_A
           + beta = mu*nu^2*(1-n)/(2*m)
           + gamma = E_A*(1-nu)/(8*m) - (lambda+2*mu)/8 + alpha/2 - beta
           + m = 1 - nu - 2*n*nu^2
           + n = E_A/(2*mu*(1+nu))
       */
      // Material Constants:
      double mu = ShearModulus;
      double nu = PoissonsRatio;
      double GA = AxialShearModulus;
      double ET = ShearModulus * (2.0 * (1.0 + nu));

      /* stf_vrtn_coeff = 1 - t/x */
      double EA = 0.0;
      double stf_vrtn_coeff = apf::getScalar(stf_vrtn_fld, apf::getMeshEntity(me), current_integration_point);
/*
      if (stf_vrtn_coeff > 0.0)
        EA = stf_vrtn_coeff * ET + (1.0 - stf_vrtn_coeff) * AxialYoungsModulus;
      else
        EA = AxialYoungsModulus;
*/
      if (stf_vrtn_coeff > 0.0)
	EA = stf_vrtn_coeff * AxialYoungsModulus;
      else
	EA = AxialYoungsModulus;
      apf::setScalar(EA_fld, apf::getMeshEntity(me), current_integration_point, EA);
      double n = ( 2.0 * mu * (1.0 + nu) )/EA; // Note typo in paper.
      double m = 1 - nu - 2.0 * n * nu * nu;
      double lambda = ( 2.0 * mu * ( nu + n * nu * nu ) )/m;
      double alpha = mu - GA;
      double beta = ( mu * nu * nu * (1.0 - n) )/( 2.0 * m );
      double gamma = EA * (1.0 - nu )/( 8.0 * m ) - ( lambda + 2.0 * mu )/8.0 + alpha/2.0 - beta;
      // Cauchy Stress Tensor
      apf::Matrix3x3 Cauchy(0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0);
      apf::Matrix3x3 Cauchy_nh(0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0);
      apf::Matrix3x3 Cauchy_trns(0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0);
      if (detF < 0.0){
        std::cout<<"error: detF < 0"<<std::endl;
        exit (EXIT_FAILURE);
      }
      for (int ii = 0; ii < 3; ++ii)
        for (int jj = 0; jj < 3; ++jj)
        {
          // From Bonet and Burton Comput. Methods Appl. Mech. Engrg. 162 (1998) 151-164 Eq.(37)
//        Cauchy_nh[ii][jj] = mu/detF * ( leftCauchyGreen[ii][jj] - kronDel(ii,jj) ) + lambda * ( detF - 1.0 ) * kronDel(ii,jj);
          // From Bonet and Wood PP.163
          Cauchy_nh[ii][jj] = mu/detF * ( leftCauchyGreen[ii][jj] - kronDel(ii,jj) ) + lambda/detF * log(detF) * kronDel(ii,jj);
          // From Bonet and Burton Comput. Methods Appl. Mech. Engrg. 162 (1998) 151-164 Eq.(68)
          Cauchy_trns[ii][jj] = 2.0 * beta/detF * ( a * a - 1.0 ) * kronDel(ii,jj)
            + 2.0/detF * ( alpha + 2.0 * beta * log(detF) + 2 * gamma * ( a * a - 1.0 ) ) * a[ii] * a[jj]
            - alpha/detF * ( (leftCauchyGreen*a)[ii] * a[jj] + a[ii] * (leftCauchyGreen*a)[jj] );
        }
      Cauchy = Cauchy_nh + Cauchy_trns;
      /*
        Voigt Format for spatial elasticity tensor (c++ index begins with 0)
          [ c_0000 c_0011 c_0022 c_0001 c_0012 c_0002 ]   [ D_00 D_01 D_02 D_03 D_04 D_05 ]
          [        c_1111 c_1122 c_1101 c_1112 c_1102 ]   [      D_11 D_12 D_13 D_14 D_15 ]
          [               c_2222 c_2201 c_2212 c_2202 ] = [           D_22 D_23 D_24 D_25 ]
          [                      c_0101 c_0112 c_0102 ]   [                D_33 D_34 D_35 ]
          [                             c_1212 c_1202 ]   [                     D_44 D_45 ]
          [                                    c_0202 ]   [                          D_55 ]
       */
      apf::DynamicMatrix D(6,6); D.zero();
      D(0,0) = SpatialElasticityTensor(0,0,0,0,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma);
      D(0,1) = SpatialElasticityTensor(0,0,1,1,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(1,0) = D(0,1);
      D(0,2) = SpatialElasticityTensor(0,0,2,2,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(2,0) = D(0,2);
      D(0,3) = SpatialElasticityTensor(0,0,0,1,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(3,0) = D(0,3);
      D(0,4) = SpatialElasticityTensor(0,0,1,2,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(4,0) = D(0,4);
      D(0,5) = SpatialElasticityTensor(0,0,0,2,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(5,0) = D(0,5);

      D(1,1) = SpatialElasticityTensor(1,1,1,1,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma);
      D(1,2) = SpatialElasticityTensor(1,1,2,2,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(2,1) = D(1,2);
      D(1,3) = SpatialElasticityTensor(1,1,0,1,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(3,1) = D(1,3);
      D(1,4) = SpatialElasticityTensor(1,1,1,2,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(4,1) = D(1,4);
      D(1,5) = SpatialElasticityTensor(1,1,0,2,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(5,1) = D(1,5);

      D(2,2) = SpatialElasticityTensor(2,2,2,2,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma);
      D(2,3) = SpatialElasticityTensor(2,2,0,1,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(3,2) = D(2,3);
      D(2,4) = SpatialElasticityTensor(2,2,1,2,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(4,2) = D(2,4);
      D(2,5) = SpatialElasticityTensor(2,2,0,2,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(5,2) = D(2,5);

      D(3,3) = SpatialElasticityTensor(0,1,0,1,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma);
      D(3,4) = SpatialElasticityTensor(0,1,1,2,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(4,3) = D(3,4);
      D(3,5) = SpatialElasticityTensor(0,1,0,2,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(5,3) = D(3,5);

      D(4,4) = SpatialElasticityTensor(1,2,1,2,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma);
      D(4,5) = SpatialElasticityTensor(1,2,0,2,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma); D(5,4) = D(4,5);

      D(5,5) = SpatialElasticityTensor(0,2,0,2,a,leftCauchyGreen,detF,lambda,alpha,beta,gamma);
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
            tau(jj+ii*3,kk+ii*3) = Cauchy[jj][kk];
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
      CauchyVoigt[0] = Cauchy[0][0];
      CauchyVoigt[1] = Cauchy[1][1];
      CauchyVoigt[2] = Cauchy[2][2];
      CauchyVoigt[3] = Cauchy[0][1];
      CauchyVoigt[4] = Cauchy[1][2];
      CauchyVoigt[5] = Cauchy[0][2];
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
      apf::destroyElement(e_disp);
    }
    int current_integration_point;
  private:
    apf::Field * stf_vrtn_fld;
    apf::Field * EA_fld;
    apf::Field * dfm_grd_fld;
    NonlinearTissue * analysis;
    int dim;
    apf::FieldShape * fs;
    apf::EntityShape * es;
    apf::Element * ce;
    apf::Mesh * msh;
    apf::NewArray<apf::Vector3> msh_xyz;
    apf::NewArray<apf::Vector3> fld_xyz;
    double kronDel(int const ii, int const jj)
    {
      double r = 0.0;
      if (ii==jj)
        r = 1.0;
      return r;
    }
    double SpatialElasticityTensor(int const ii, int const jj, int const kk, int const ll, apf::Vector3 const& a, apf::Matrix3x3 const& b,
                                   double detF, double lambda, double alpha, double beta, double gamma)
    {
      /*
        Spatial Elasticity Tensor
        - c = c_nh + c_trns (see J. Bonet and A.J. Burton Comput. Methods Appl. Mech. Engrg. 162 (1998) 151-164)
        - c_nh_ijkl = lambda*(2*J-1)*I_ij*I_kl + 2/J( mu - lambda*J*(2*J-1) )delta_ik*delta_jl
        - c_trns_ijkl = 8*gamma/J*a_i*a_j*a_k*a_l
                      + 4*beta/J*(a_i*a_j*delta_kl + delta_ij*a_k*a_l)
                      - alpha/J*(a_i*a_l*b_jk + b_ik*a_j*a_l)
                      - 4*beta/J(a_r*a_r - 1.0)*delta_ik*delta_jl

       - C_rs: elements of right Cauchy-Green deformation tensor (rightCauchyGreen).
       - b_is: elements of left Cauchy-Green deformation tensor (leftCauchyGreen).
       - J   : determinant of deformation gradient (detF).
       - Material Constants:
         + mu = Shear Modulus, nu = Poisson Ratio, G_A = axial shear modulus, E_A = axial Young's modulus.
         + lambda = 2*mu*(nu+n*nu^2)/m
         + alpha = mu - G_A
         + beta = mu*nu^2*(1-n)/(2*m)
         + gamma = E_A*(1-nu)/(8*m) - (lambda+2*mu)/8 + alpha/2 - beta
         + m = 1 - nu - 2*n*nu^2
         + n = E_A/(2*mu*(1+nu))
      */
      if (ii > 2 || jj > 2 || kk > 2 || ll > 2)
      {
        std::cout<<"index error in spatial elasticity tensor"<<std::endl;
        exit (EXIT_FAILURE);
      }
      double mu = ShearModulus;
      double c = 0.0; double c_nh = 0.0; double c_trns = 0.0;
      // From Bonet and Burton Comput. Methods Appl. Mech. Engrg. 162 (1998) 151-164 Eq.(39)
//      c_nh = lambda * ( 2.0 * detF - 1.0 ) * kronDel(ii,jj) * kronDel(kk,ll) + 2.0/detF * ( mu - lambda * detF * ( 2.0 * detF - 1.0 ) ) * kronDel(ii,kk) * kronDel(jj,ll);
      // From Bonet and Wood PP.164
      c_nh = lambda/detF * kronDel(ii,jj) * kronDel(kk,ll) + 2.0/detF * ( mu - lambda * log(detF) ) * 0.5 * ( kronDel(ii,kk) * kronDel(jj,ll) + kronDel(ii,ll) * kronDel(jj,kk) );
      // From Bonet and Burton Comput. Methods Appl. Mech. Engrg. 162 (1998) 151-164 Eq.(70)
      c_trns = 8.0 * gamma/detF * a[ii] * a[jj] * a[kk] * a[ll]
        + 4.0 * beta/detF * ( a[ii] * a[jj] * kronDel(kk,ll) + kronDel(ii,jj) * a[kk] * a[ll] )
        - alpha/detF * ( a[ii] * a[ll] * b[jj][kk] + b[ii][kk] * a[jj] * a[ll] )
        - 4.0 * beta/detF * ( a * a - 1.0 ) * kronDel(ii,kk) * kronDel(jj,ll);
      c = c_nh + c_trns;
      return c;
    }
    double ShearModulus;
    double PoissonsRatio;
    double Axis[3];
    double AxialShearModulus;
    double AxialYoungsModulus;
  };
}
#endif
