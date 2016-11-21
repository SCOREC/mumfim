#include "NonLinFibMtx.h"
#include "NonFiniteElement.h"
#include <ConvenienceFunctions.h>
#include <apfSIM.h>
#include <cassert>
#include <memory.h>
#include <vector>
#include <list>
#include <map>
#include <iostream>
namespace Biotissue
{
  // get averaged stress
  // S_ij = 1/V * integral(x_j * S_ik * n_k)dA
  // All coordinates expressed in current coorindates (current equilibrium equation)
  void NonLinFibMtx::AverageStress(apf::Matrix3x3 & stress)
    {
      AverageMatrixStress(matrix_stress);
      AverageFiberStress(fiber_stress);
      stress = matrix_stress + fiber_stress;
    }
    void NonLinFibMtx::UnbalancedForce(apf::Vector3 & term)
    {
      UnbalancedMatrixForce(matrix_term);
      UnbalancedFiberForce(fiber_term);
      term = matrix_term + fiber_term;
    }
    void NonLinFibMtx::AverageMatrixStress(apf::Matrix3x3 & S)
    {
      apf::Matrix3x3 stress;
      GFIter gfiter = GM_faceIter(model);
      while(pGFace gface = GFIter_next(gfiter))
      {
        apf::Vector3 normal;
        FaceNormByGrad(gface,normal);
        std::list<pEntity> faces;
        Model_GetClassifiedEntities(part,gface,2,faces);
        for(std::list<pEntity>::iterator iter = faces.begin(),
              iterend = faces.end(); iter != iterend; iter++)
        {
          pFace face = (pFace)*iter;
          pRegion adjacent_region;
          Mesh_BoundaryFace_GetRegion(face, adjacent_region);
          CauchyStress(adjacent_region, stress);
          // compute integral(x_j * S_ik * n_k)dA at this boundary face(not sure if it needs to be deformed surface or not)
          //compute the traction term S_ik * n_k
          apf::Vector3 traction;
          for(int ii = 0; ii < 3; ii++)
            for(int jj = 0; jj < 3; jj++)
              traction[ii] += stress[ii][jj] * normal[jj]; // S_ij * n_j
          apf::MeshElement * me = apf::createMeshElement(apf_mesh,apf::castEntity(face));
          int num_pts = apf::countIntPoints(me,1);
          for(int gauss_index = 0; gauss_index < num_pts; gauss_index++)
          {
            apf::Vector3 gauss_point;
            apf::getIntPoint(me,1,gauss_index,gauss_point);
            double w = apf::getIntWeight (me,1,gauss_index);
            apf::Matrix3x3 J;
            apf::getJacobian(me,gauss_point,J);
            double detjac = apf::getJacobianDeterminant(J,getDimension(me));
            apf::Vector3 x;
            apf::mapLocalToGlobal(me,gauss_point,x);
            for(int ii = 0; ii < 3; ii++) // really should be number of field compone
              for(int jj = 0; jj < 3; jj++)
                stress[ii][jj] += x[jj] * traction[ii] * w * detjac;
          }
        }
      }
      GFIter_delete(gfiter);
      /*
      S[0] = stress[0][0];
      S[1] = stress[1][1];
      S[2] = stress[2][2];
      S[3] = 0.5 * (stress[0][1] + stress[1][0]);
      S[4] = 0.5 * (stress[1][2] + stress[2][1]);
      S[5] = 0.5 * (stress[2][0] + stress[0][2]);
      */
      double volume = apf::getDeterminant(deformation_gradient) * standard_volume;
      S = stress / volume;
    }
  // ??????????
    void NonLinFibMtx::Mesh_BoundaryFace_GetRegion(pFace face,
                                                   pRegion & region)
    {
      pPList regions = PList_new();
      regions = F_regions(face);
      assert(PList_size(regions) == 1);
      void * iter = 0;
      region = (pRegion)PList_next(regions, &iter);
      PList_delete(regions);
    }
    // todo: move into seperate .cc for utility functions
    void Stress_VectorToMatrix(const apf::Vector<6> & s1, apf::Matrix3x3 & s2)
    {
      s2[0][0] = s1[0];
      s2[1][1] = s1[1];
      s2[2][2] = s1[2];
      s2[1][0] = s2[0][1] = s1[3];
      s2[2][1] = s2[1][2] = s1[4];
      s2[0][2] = s2[2][0] = s1[5];
    }
    void NonLinFibMtx::FaceNormByGrad(pGFace gface, apf::Vector3 & normal)
    {
      //double par[2] = {0,0};
      //double xyz[3] = {0,0,0};
      //GF_normal(gface,par,xyz); // the unit normal in undeformed configuration
      //apf::Vector3 N(xyz);
      // Nanson's formula
      //double J = apf::getDeterminant(deformation_gradient);
      //apf::Matrix3x3 grad_inv = apf::invert(deformation_gradient);
      // TODO: fix below
      //N = grad_inv * N;
      //N = J * N;
      //normal = N.normalize();
    }
// compute the unbalanced term from the matrix
    void NonLinFibMtx::UnbalancedMatrixForce(apf::Vector3 & term)
    {
      matrix_term = matrix_term * 0.0; // reset the term
      GFIter gfiter = GM_faceIter(model);
      pGFace gface;
      apf::Vector3 normal;
      while((gface = GFIter_next(gfiter)))
      {
        double norm[3] = {0.0,0.0,0.0};
        double p[3] = {0.33,0.33,0.33}; // assumes barycentric coordinates
        GF_normal(gface,&p[0],&norm[0]);
        normal.fromArray(norm);
        std::list<pEntity> faces;
        // get all mesh faces on the model face
        Model_GetClassifiedEntities(part,gface,2,faces);
        for(std::list<pEntity>::iterator iter = faces.begin(),
              iterend = faces.end(); iter != iterend; iter++)
        {
          pFace face = (pFace)*iter;
          pRegion adjacent_region;
          Mesh_BoundaryFace_GetRegion(face, adjacent_region);
          apf::Matrix3x3 cauchy_stress;
          CauchyStress(adjacent_region,cauchy_stress);
          apf::Matrix3x3 coefficients;
          coefficients = cauchy_stress - total_stress;
          apf::MeshElement * me = apf::createMeshElement(apf_mesh,apf::castEntity((pEntity)face));
          // compute a surface integral for the unbalanced term
          int num_gauss_points = apf::countIntPoints(me,1);
          for(int gauss_index = 0; gauss_index < num_gauss_points; gauss_index++)
          {
            apf::Vector3 gauss_point;
            apf::getIntPoint(me,1,gauss_index,gauss_point);
            double w = apf::getIntWeight(me,1,gauss_index);
            apf::Matrix3x3 J;
            apf::getJacobian(me,gauss_point,J);
            double detjac = apf::getJacobianDeterminant(J,apf::getDimension(me));
            apf::Matrix3x3 disp_grad;
            SpatialDispGrad(disp_grad);
            apf::Vector3 derivUNormal;
            for(int ii = 0; ii < 3; ii++)
              for(int jj = 0; jj < 3; jj++)
                derivUNormal[ii] += disp_grad[jj][ii] * normal[jj];
            for(int ii = 0; ii < 3; ii++)
              for(int jj = 0; jj < 3; jj++)
                term[ii] += w * detjac * coefficients[jj][ii] * derivUNormal[jj];
          }
        }
      }
    }
    // calculates I - F^{-1} for the specified location in the given region, where F is the deformation gradient at the same point
    void NonLinFibMtx::LocalSpatialDispGrad(pRegion region,
                                            const apf::Vector3 & pos,
                                            apf::Matrix3x3 & disp_grad)
    {
      apf::MeshEntity * me = apf::castEntity((pEntity)region);
      apf::Element * e = apf::createElement(apf_primary_field,me);
      apf::NewArray<apf::Vector3> region_nodes;
      apf::getVectorNodes(e,region_nodes);
      apf::Matrix3x3 gradient;
      apf::getVectorGrad(e,pos,gradient);
      // disp_deriv holds def_grad_inverse here;
      disp_grad = apf::invert(gradient);
      disp_grad = Identity3x3 - disp_grad;
    }
    void NonLinFibMtx::CauchyStress(pEntity entity,
                                    apf::Matrix3x3 & cauchy_stress)
    {
      // should probably implement as an integrator...
      apf::MeshElement * me = apf::createMeshElement(apf_mesh,apf::castEntity(entity));
      apf::Element * e = apf::createElement(apf_primary_field,me);
      apf::NewArray<apf::Vector3> nodal_displacements;
      apf::getVectorNodes(e,nodal_displacements);
      int num_gauss_pts = apf::countIntPoints(me,1);
      for(int gauss_index = 0; gauss_index < num_gauss_pts; gauss_index++)
      {
        apf::Vector3 gauss_pt;
        apf::getIntPoint(me,1,gauss_index,gauss_pt);
        double w = apf::getIntWeight(me,1,gauss_index);
        apf::Matrix3x3 Jac;
        apf::getJacobian(me,gauss_pt,Jac);
        double detjac = apf::getJacobianDeterminant(Jac,apf::getDimension(me));
        apf::Matrix3x3 def_grad;
        apf::getVectorGrad(e,gauss_pt,def_grad);
        double det_def_grad = apf::getDeterminant(def_grad);
        apf::Matrix3x3 right_cauchy;
        RightCauchy(def_grad,right_cauchy);
        apf::Matrix3x3 pk2_stress;
        PK2Stress(right_cauchy,
                  det_def_grad,
                  poisson_ratio,
                  shear_modulus,
                  pk2_stress);
        apf::Matrix3x3 local_cauchy;
        PK2ToCauchy(def_grad,pk2_stress,local_cauchy);
        cauchy_stress = cauchy_stress + (local_cauchy * w * detjac);
      }
    }
    void NonLinFibMtx::PK2ToCauchy(const apf::Matrix3x3 & def_grad,
                                   const apf::Matrix3x3 & pk2_stress,
                                   apf::Matrix3x3 & cauchy_stress)
    {
      double inv_det_def_grad = 1.0 / apf::getDeterminant(def_grad);
      // sigma = 1/det(F) * F * S * F^T
      // sigma - cauchy stress
      // F - deformation gradient
      // S - pk2 stress
      cauchy_stress = pk2_stress * apf::transpose(def_grad);
      cauchy_stress = def_grad * cauchy_stress;
      cauchy_stress = cauchy_stress * inv_det_def_grad;
      //cauchy_stress = inv_det_def_grad * def_grad * pk2_stress * apf::transpose(def_grad);
    }
  }
