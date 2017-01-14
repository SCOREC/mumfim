#ifndef BIO_VOLUME_CONSTRAINT_SURFACE_H_
#define BIO_VOLUME_CONSTRAINT_SURFACE_H_
#include "bioNonlinearTissue.h"
#include <apfFunctions.h> //<faceNormal and triangleArea.
#include <amsiLAS.h>
#include <ElementalSystem.h>
#include <simAnalysis.h>
#include <apfShape.h>
#include <apfSIM.h>
#include <cstring>
/** VolumeConstraintSurface class that contains implementation of volume constraint based on surface mesh.
 * Calculations are based on Hong et al., "Fast Volume Preservation for a Mass-Spring System," IEEE Computer Graphics and Applications (2006) */
namespace bio
{
  class VolumeConstraintSurface : public amsi::ElementalSystem
  {
  public:
  VolumeConstraintSurface(pGFace rg, int tag, pMesh part, apf::Field * field, int o)
    : ElementalSystem(field,o)
      , local(false)
      , updateG_flag(true)
      , updateH_flag(false)
      , dim(0)
      , dof(-2)
      , rgn(rg)
      , rgn_tag(tag)
      , d2Vdu2(0,0)
      , dVdu(0,0)
      , Lambda_d2Vdu2(0,0)
      , Lambda_dVdu(0,0)
      , DeltaV(0)
      , Lambda(0.0)
      , Beta(-2.5)
      , Vol_glb(0)
      , initVol_glb(0)
      , prevVol_glb(0)
      , elem_num(0)
    {
      //  right now all dofs associated with this constraint are just assigned to the "last" rank in the scale
      int sz = 0;
      int rnk = -1;
      MPI_Comm_rank(AMSI_COMM_SCALE,&rnk);
      MPI_Comm_size(AMSI_COMM_SCALE,&sz);
      if(rnk == sz - 1)
        local = true;
    }
    // TODO:
    // passing in the mesh partition here is a stopgap to avoid having to develop an intermediary at the moment,
    // really some abstract mesh iterator range should be constructed and passed to the constraint constructor, and then used
    // by apply to iterate over the applicable mesh entities
    void apply(amsi::LAS * las, apf::Mesh * msh, pMesh prt, apf::Numbering * nm)
    {
      int dm = 2; ///< For surface mesh entities.
      std::vector<pEntity> msh_rgns;
      amsi::getClassifiedEnts(prt,rgn,dm,std::back_inserter(msh_rgns));
      elem_num = 0; //< reset element number counter.
      for(std::vector<pEntity>::iterator mrgn = msh_rgns.begin(); mrgn != msh_rgns.end(); mrgn++)
      {
        apf::MeshEntity * mnt = apf::castEntity(*mrgn);
        apf::MeshElement * mlm = apf::createMeshElement(msh,mnt);
        process(mlm);
        apf::NewArray<int> ids;
        apf::getElementNumbers(nm,mnt,ids);
        apf::DynamicMatrix & lG = getLambdaFirstVolDeriv();
        // Multiply appropriate terms by -1
        lG *= -1.0;
        amsi::assembleVector(las,nedofs,&ids[0],&lG(0,0));
        elem_num++;
      }
      // Turn off updateG_flag and updateH_flag.
      setGflag(false);
      setHflag(false);
    }
    void inElement(apf::MeshElement * me)
    {
      ElementalSystem::inElement(me);
      fs = apf::getShape(f);
      es = fs->getEntityShape(apf::getMesh(f)->getType(apf::getMeshEntity(me)));
      d2Vdu2.setSize(nedofs,nedofs);
      dVdu.setSize(1,nedofs);
      // hard-coded, should be derived from me
      // bill: can't we just add one to the dim of the meshelement, since this reduces to an area constraint calculated from edges in 2d, and a length constraint from vertices in 1d... the dim value should work correctly in all cases
      dim = 3;
    }
    bool includesBodyForces() { return true; }
    void atPoint(apf::Vector3 const &p, double w, double dV)
    {
      int & nen = nenodes; // = 3 (triangle)
      apf::Mesh * mesh = apf::getMesh(f);
      apf::Field * apf_coord_field = mesh->getCoordinateField();
      apf::Element * mesh_coord_elem = apf::createElement(apf_coord_field,me);
      apf::NewArray<apf::Vector3> mesh_xyz;
      apf::getVectorNodes(mesh_coord_elem,mesh_xyz);
      apf::NewArray<apf::Vector3> primary_field_xyz;
      apf::getVectorNodes(e,primary_field_xyz);
      apf::DynamicMatrix xyz(nen,dim);
      xyz.zero();
      for (int ii = 0; ii < nen; ii++)
        for (int jj = 0; jj < dim; jj++)
          xyz(ii,jj) = mesh_xyz[ii][jj] + primary_field_xyz[ii][jj];
      apf::Vector3 pt0;
      apf::Vector3 pt1;
      apf::Vector3 pt2;
      for (int jj = 0; jj < dim; jj++)
      {
        pt0[jj] = xyz(0,jj);
        pt1[jj] = xyz(1,jj);
        pt2[jj] = xyz(2,jj);
      }
      calcG(pt0, pt1, pt2, dVdu);
      if (updateH_flag)
        calcH(pt0, pt1, pt2, d2Vdu2, nen);
      // Determine normal direction (based on measureVol_pGFace in apfsimWrapper.cc)
      // todo: make this a standalone function... should be pretty straightforward
      pRegion mshRgnIn = F_region((pFace)apf::getMeshEntity(me),0);
      pRegion mshRgnOut = F_region((pFace)apf::getMeshEntity(me),1);
      int normal_dir = 1;
      int GRgnIn_tag = mshRgnIn == NULL ? -2 : GEN_tag(R_whatIn(mshRgnIn));
      int GRgnOut_tag = mshRgnOut == NULL ? -2 : GEN_tag(R_whatIn(mshRgnOut));
      if (GRgnIn_tag == rgn_tag)
        normal_dir = 1;
      else if (GRgnOut_tag == rgn_tag)
        normal_dir = -1;
      else // have to deal with partitioning..
      {
        if (GRgnIn_tag == -2 && GRgnOut_tag != rgn_tag)
          normal_dir = 1;
        else if (GRgnOut_tag == -2 && GRgnIn_tag != rgn_tag)
          normal_dir = -1;
      }
      DeltaV = Vol_glb - prevVol_glb;
      dVdu *= normal_dir;
      d2Vdu2 *= normal_dir;
      Lambda_d2Vdu2 = d2Vdu2;
      Lambda_d2Vdu2 *= Lambda;
      // Modifications to force-vector.
      Lambda_dVdu = dVdu;
      Lambda_dVdu *= Lambda;
      // Augmented Lagrangian method:
      // A   = Beta * (GG^T + DeltaV * H)
      // VH  = DeltaV * H
      // BVG = Beta * DeltaV * G
      // Lambda * H modified to Lambda * H + Beta * (GG^T + DeltaV * H)
      apf::DynamicMatrix A(nedofs,nedofs);
      apf::DynamicMatrix GT(nedofs,1);
      apf::transpose(dVdu,GT);
      apf::multiply(GT,dVdu,A);
      apf::DynamicMatrix VH(nedofs,nedofs);
      VH = d2Vdu2;
      VH *= DeltaV;
      A += VH;
      A *= Beta;
      A /= initVol_glb; //scale entire A by V0^2.
      Lambda_d2Vdu2 += A;
      apf::DynamicMatrix BVG(1,nedofs);
      BVG = dVdu;
      BVG *= DeltaV;
      BVG *= Beta;
      BVG /= prevVol_glb;
      Lambda_dVdu += BVG;
    }
    void calcG(apf::Vector3 const &pt0, apf::Vector3 const &pt1, apf::Vector3 const &pt2, apf::DynamicMatrix& dVdu)
    {
      dVdu.zero();
      // volume derivatives
      // dVdu is of size 1 x nedofs, the below shouldn't be hardcoded for triangles in that case..
      dVdu(0,0) = 0.5 * (pt2[1] * pt1[2] - pt1[1] * pt2[2]);  ///<dVdx1
      dVdu(0,1) = 0.5 * (-pt2[0] * pt1[2] + pt1[0] * pt2[2]); ///<dVdy1
      dVdu(0,2) = 0.5 * (pt2[0] * pt1[1] - pt1[0] * pt2[1]);  ///<dVdz1
      dVdu(0,3) = 0.5 * (-pt2[1] * pt0[2] + pt0[1] * pt2[2]); ///<dVdx2
      dVdu(0,4) = 0.5 * (pt2[0] * pt0[2] - pt0[0] * pt2[2]);  ///<dVdy2
      dVdu(0,5) = 0.5 * (-pt2[0] * pt0[1] + pt0[0] * pt2[1]); ///<dVdz2
      dVdu(0,6) = 0.5 * (pt1[1] * pt0[2] - pt0[1] * pt1[2]);  ///<dVdx3
      dVdu(0,7) = 0.5 * (-pt1[0] * pt0[2] + pt0[0] * pt1[2]); ///<dVdy3
      dVdu(0,8) = 0.5 * (pt1[0] * pt0[1] - pt0[0] * pt1[1]);  ///<dVdz3
    } // end of calcG
    void calcH(apf::Vector3 const & pt0,
               apf::Vector3 const & pt1,
               apf::Vector3 const & pt2,
               apf::DynamicMatrix & dVdu,
               int nen)
    {
      d2Vdu2.zero();
      /*
        submat1 = [ 0  -z3  y3]
                  [ z3   0 -x3]
                  [-y3  x3   0]
        submat2 = [  0  z2 -y2]
                  [-z2   0  x2]
                  [ y2 -x2   0]
        submat3 = [  0 -z1  y1]
                  [ z1   0 -x1]
                  [-y1  x1   0]
      */
      apf::Matrix<3,3> submat1;
      apf::Matrix<3,3> submat2;
      apf::Matrix<3,3> submat3;
      submat1[0][0] = 0.0;     submat1[0][1] = -pt2[2]; submat1[0][2] = pt2[1];
      submat1[1][0] = pt2[2];  submat1[1][1] = 0.0;     submat1[1][2] = -pt2[0];
      submat1[2][0] = -pt2[1]; submat1[2][1] = pt2[0];  submat1[2][2] = 0.0;
      submat2[0][0] = 0.0;     submat2[0][1] = pt1[2];  submat2[0][2] = -pt1[1];
      submat2[1][0] = -pt1[2]; submat2[1][1] = 0.0;     submat2[1][2] = pt1[0];
      submat2[2][0] = pt1[1];  submat2[2][1] = -pt1[0]; submat2[2][2] = 0.0;
      submat3[0][0] = 0.0;     submat3[0][1] = -pt0[2]; submat3[0][2] = pt0[1];
      submat3[1][0] = pt0[2];  submat3[1][1] = 0.0;     submat3[1][2] = -pt0[0];
      submat3[2][0] = -pt0[1]; submat3[2][1] = pt0[0];  submat3[2][2] = 0.0;
      for (int ii = 0; ii < nen; ii++)
        for (int jj = 0; jj < nen; jj++)
        {
          // Diagonal submatrices
          d2Vdu2(ii,jj) = 0.0;
          d2Vdu2(ii+dim,jj+dim) = 0.0;
          d2Vdu2(ii+2*dim,jj+2*dim) = 0.0;
          // off-diagonal submatrices
          d2Vdu2(ii,jj+dim) = 0.0;
          d2Vdu2(ii+dim,jj) = 0.0;
          d2Vdu2(ii,jj+2*dim) = 0.0;
          d2Vdu2(ii+2*dim,jj) = 0.0;
          d2Vdu2(ii+dim,jj+2*dim) = 0.0;
          d2Vdu2(ii+2*dim,jj+dim) = 0.0;
        }
    } // end of calcH
    apf::DynamicMatrix& getFirstVolDeriv(){return dVdu;}
    apf::DynamicMatrix& getSecondVolDeriv(){return d2Vdu2;}
    apf::DynamicMatrix& getLambdaFirstVolDeriv(){return Lambda_dVdu;}
    apf::DynamicMatrix& getLambdaSecondVolDeriv(){return Lambda_d2Vdu2;}
    double getDeltaV(){return DeltaV;}
    void setDof(int d) { dof = d; }
    int getDof() { return dof; }
    bool isLocal() { return local; }
    double getLambda(){return Lambda;}
    double getBeta(){return Beta;}
    double getVol(){return Vol_glb;}
    double getinitVol(){return initVol_glb;}
    void setLambda(double Delta_Lambda) { Lambda = Delta_Lambda; }
    void setBeta(double factor) { Beta = factor; }
    void setGflag(bool flag){updateG_flag = flag;}
    void setHflag(bool flag){updateH_flag = flag;}
    void setVol(double vol){Vol_glb = vol;}
    void setInitVol(double vol){initVol_glb = vol;}
    void setPrevVol(double vol){prevVol_glb = vol;}
    pGFace getFace() { return rgn; }
    int getRegionTag() {return rgn_tag;}
  private:
    bool local;
    bool updateG_flag;
    bool updateH_flag;
    int dim;
    int dof;
    pGFace rgn;
    int rgn_tag;
    apf::FieldShape * fs;
    apf::EntityShape * es;
    apf::DynamicMatrix d2Vdu2;
    apf::DynamicMatrix dVdu;
    apf::DynamicMatrix Lambda_d2Vdu2;
    apf::DynamicMatrix Lambda_dVdu;
    double DeltaV;
    double Lambda; // Lagrange multiplier for Volume Preservation.
    double Beta;   // Penalty parameter for Augmented Lagrangian Method.
    double Vol_glb;     // global current volume of entire structure being constrained.
    double initVol_glb; // global initial volume of entire structure being constrained.
    double prevVol_glb; // global previous volume of entire structure being constrained.
    int elem_num; // keeps count of element number within region.
  };
}
#endif
