#ifndef BIO_VOLUME_CONSTRAINT_H_
#define BIO_VOLUME_CONSTRAINT_H_
#include "bioNonlinearTissue.h"
#include <amsiLAS.h>
#include <ElementalSystem.h>
#include <simAnalysis.h>
#include <apfShape.h>
#include <apfSIM.h>
#include <cstring>
/*
  Contains two classes VolumeConstraint and VolumeConstraintADMM which are used for imposing Volume Constraint on the system.
  TO DO: VolumeConstraint and VolumeConstraintADMM classes have ~90% of the system code. Need to refactor to reduce duplication of code.
*/
namespace bio
{
  class VolumeConstraint : public amsi::ElementalSystem
  {
  public:
  VolumeConstraint(pGEntity rg, apf::Field * field,int o)
    : ElementalSystem(field,o)
      , local(false)
      , dim(0)
      , dof(-2)
      , rgn(rg)
      , d2Vdu2(0,0)
      , dVdu(0,0)
      , Lambda_d2Vdu2(0,0)
      , Lambda_dVdu(0,0)
      , DeltaV(0)
      , Lambda(0)
      , Beta(0)
      , Vol(0)
      , initVol(0)
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
      dim = msh->getDimension();
      std::vector<pEntity> rgns;
      amsi::getClassifiedEnts(prt,rgn,dim,std::back_inserter(rgns));
      for(std::vector<pEntity>::iterator mrgn = rgns.begin(); mrgn != rgns.end(); mrgn++)
      {
        apf::MeshEntity * mnt = apf::castEntity(*mrgn);
        apf::MeshElement * mlm = apf::createMeshElement(msh,mnt);
        process(mlm);
        apf::NewArray<int> ids;
        apf::getElementNumbers(nm,mnt,ids);
        apf::DynamicMatrix & lG = getLambdaFirstVolDeriv();
        apf::DynamicMatrix & lH = getLambdaSecondVolDeriv();
        apf::DynamicMatrix & G = getFirstVolDeriv();
        // Multiply appropriate terms by -1
        lG *= -1.0;
        double dvl = getDeltaV() * -1.0;
        amsi::assembleMatrix(las,nedofs,&ids[0],nedofs,&ids[0],&lH(0,0));
        amsi::assembleMatrix(las,nedofs,&ids[0],1,&dof,&G(0,0));
        amsi::assembleMatrix(las,1,&dof,nedofs,&ids[0],&G(0,0)); // Transpose of G.
        amsi::assembleVector(las,nedofs,&ids[0],&lG(0,0));
        amsi::assembleVector(las,1,&dof,&dvl);
      }
    }
    void inElement(apf::MeshElement * me)
    {
      ElementalSystem::inElement(me);
      fs = apf::getShape(f);
      es = fs->getEntityShape(apf::getMesh(f)->getType(apf::getMeshEntity(me)));
      d2Vdu2.setSize(nedofs,nedofs);
      dVdu.setSize(1,nedofs);
      d2Vdu2.zero();
      dVdu.zero();
      Beta = 1e6;
    }
    bool includesBodyForces() { return true; }
    void atPoint(apf::Vector3 const &p, double w, double dV)
    {
      int & nen = nenodes; // = 4 (tets)
      // Jacobian of underlying mesh. (Initial configuration).
      apf::Matrix3x3 Jac;
      apf::getJacobian(me,p,Jac);
      // Note: Form of Jacobian in apf is
      // [dx/dr dx/ds dx/dt]    (x,y,z): physical domain coordinates.
      // [dy/dr dy/ds dy/dt]    (r,s,t): parent domain coordiantes.
      // [dz/dr dz/ds dz/dt]
      double wxdetjac = w * apf::getJacobianDeterminant(Jac,dim);
      // Jacobian of updated configuration
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
      // For Updated Lagrangian, the Jacobian of the updated coordinates are used
      // Note: that entires of Jacobian is hard coded for Linear tetrahedra elements.
      // TO DO: Generalize Jacobian for current configuration.
      apf::Matrix<3,3> J;
      J[0][0] = xyz(1,0) - xyz(0,0); // x2-x1
      J[0][1] = xyz(2,0) - xyz(0,0); // x3-x1
      J[0][2] = xyz(3,0) - xyz(0,0); // x4-x1
      J[1][0] = xyz(1,1) - xyz(0,1); // y2-y1
      J[1][1] = xyz(2,1) - xyz(0,1); // y3-y1
      J[1][2] = xyz(3,1) - xyz(0,1); // y4-y1
      J[2][0] = xyz(1,2) - xyz(0,2); // z2-z1
      J[2][1] = xyz(2,2) - xyz(0,2); // z3-z1
      J[2][2] = xyz(3,2) - xyz(0,2); // z4-z1
      double detJ = getDeterminant(J);
      // Determine derivatives of volume wrt to coordinates of updated configuration.
      apf::Matrix<3,3> dx(J), dy(J), dz(J);
      apf::Matrix<3,3> dxdx(J), dydy(J), dzdz(J);
      apf::Matrix<3,3> dxdy(J), dydz(J), dxdz(J);
      /*
        Derivatives of Shape function wrt parent domain coordinates (xi1,xi2,xi3)
        - local_grads: derivative of shape functions wrt to parent domain coordinates.
        - p          : integration point of mesh element wrt to parent domain coordinates.
        - mesh       : mesh on which the field f is defined. In this case, f is the apf_primary_field.
        Note: that for Linear Tetrahedral elements, derivatives of shape function wrt parent domain coodinates are constant.
      */
      apf::NewArray<apf::Vector3> local_grads;
      es->getLocalGradients(mesh, apf::getMeshEntity(me), p, local_grads);
      for (int a = 0; a < nen; a++) // loop through node a
      {
        for (int dxi = 0; dxi < dim; dxi++)
        {
          dx[0][dxi] = local_grads[a][dxi];
          dy[1][dxi] = local_grads[a][dxi];
          dz[2][dxi] = local_grads[a][dxi];
        }
        // Fill "G" vector
        dVdu(0,a*dim) = apf::getDeterminant(dx);
        dVdu(0,a*dim+1) = apf::getDeterminant(dy);
        dVdu(0,a*dim+2) = apf::getDeterminant(dz);
        for (int b = 0; b < nen; b++) // loop through node b
        {
          for (int dxi = 0; dxi < dim; dxi++)
          {
            // diagonal terms
            dxdx[0][dxi] = local_grads[a][dxi] + local_grads[b][dxi];
            dydy[1][dxi] = local_grads[a][dxi] + local_grads[b][dxi];
            dzdz[2][dxi] = local_grads[a][dxi] + local_grads[b][dxi];
            // off-diagonal terms
            dxdy[0][dxi] = local_grads[a][dxi];
            dxdy[1][dxi] = local_grads[b][dxi];
            dydz[1][dxi] = local_grads[a][dxi];
            dydz[2][dxi] = local_grads[b][dxi];
            dxdz[0][dxi] = local_grads[a][dxi];
            dxdz[2][dxi] = local_grads[b][dxi];
          }
          // Fill "H" matrix
          // diagonal terms
          d2Vdu2(a*dim,b*dim) = apf::getDeterminant(dxdx);
          d2Vdu2(a*dim+1,b*dim+1) = apf::getDeterminant(dydy);
          d2Vdu2(a*dim+2,b*dim+2) = apf::getDeterminant(dzdz);
          // off-diagonal terms
          d2Vdu2(a*dim,b*dim+1) = apf::getDeterminant(dxdy);
          d2Vdu2(a*dim,b*dim+2) = apf::getDeterminant(dxdz);
          d2Vdu2(a*dim+1,b*dim+2) = apf::getDeterminant(dydz);
          d2Vdu2(a*dim+1,b*dim) = apf::getDeterminant(dxdy) * -1.0;
          d2Vdu2(a*dim+2,b*dim) = apf::getDeterminant(dxdz) * -1.0;
          d2Vdu2(a*dim+2,b*dim+1) = apf::getDeterminant(dydz) * -1.0;
        }
      }
      // Calculate Delta V
      Vol = w * detJ;
      initVol = wxdetjac;
      DeltaV = Vol - initVol;
//      DeltaV = w * detJ - wxdetjac;
      // Modifications to tangent-stiffness matrix
      d2Vdu2 *= w;
      Lambda_d2Vdu2 = d2Vdu2;
      Lambda_d2Vdu2 *= Lambda;
      // Modifications to force-vector.
      dVdu *= w;
      Lambda_dVdu = dVdu;
      Lambda_dVdu *= Lambda;
      // Augmented Lagrangian method
      apf::DynamicMatrix A(nedofs,nedofs);
      apf::DynamicMatrix GT(nedofs,1);
      apf::transpose(dVdu,GT);
      apf::multiply(GT,dVdu,A);
      apf::DynamicMatrix VH(nedofs,nedofs);
      VH = d2Vdu2;
      VH *= DeltaV;
      A += VH;
      A *= Beta;
      Lambda_d2Vdu2 += A;
      apf::DynamicMatrix BVG(1,nedofs);
      BVG = dVdu;
      BVG *= DeltaV;
      BVG *= Beta;
      Lambda_dVdu += BVG;
    }
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
    double getVol(){return Vol;}
    double getinitVol(){return initVol;}
    void setLambda(double Delta_Lambda){Lambda += Delta_Lambda;}
    pGEntity getRegion() { return rgn; }
  private:
    bool local;
    int dim;
    int dof;
    pGEntity rgn;
    apf::FieldShape * fs;
    apf::EntityShape * es;
    apf::DynamicMatrix d2Vdu2;
    apf::DynamicMatrix dVdu;
    apf::DynamicMatrix Lambda_d2Vdu2;
    apf::DynamicMatrix Lambda_dVdu;
    double DeltaV;
    double Lambda; // Lagrange multiplier for Volume Preservation.
    double Beta;   // Penalty parameter for Augmented Lagrangian Method.
    double Vol;
    double initVol;
  };
/*
  Start of VolumeConstraintADMM class
*/
  class VolumeConstraintADMM : public amsi::ElementalSystem
  {
  public:
  VolumeConstraintADMM(pGEntity rg, pMesh part, apf::Field * field,int o)
    : ElementalSystem(field,o)
      , local(false)
      , dim(0)
      , dof(-2)
      , rgn(rg)
      , d2Vdu2(0,0)
      , dVdu(0,0)
      , Lambda_d2Vdu2(0,0)
      , Lambda_dVdu(0,0)
      , DeltaV(0)
      , Lambda(0)
      , Beta(0.0)
      , Vol(0)
      , initVol(0)
      , elem_num(0)
    {
      //  right now all dofs associated with this constraint are just assigned to the "last" rank in the scale
      int sz = 0;
      int rnk = -1;
      MPI_Comm_rank(AMSI_COMM_SCALE,&rnk);
      MPI_Comm_size(AMSI_COMM_SCALE,&sz);
      if(rnk == sz - 1)
        local = true;
      /** initialize prevVol vector 
       ** f is member variable of ElementalSystem and stores accumulated displacement (apf_primary_field)*/
      apf::Mesh * mesh = apf::getMesh(f);
      dim = mesh -> getDimension();
      std::vector<pEntity> rgns;
      amsi::getClassifiedEnts(part,rgn,dim,std::back_inserter(rgns));
      prevVol.resize(rgns.size());
    }
    // TODO:
    // passing in the mesh partition here is a stopgap to avoid having to develop an intermediary at the moment,
    // really some abstract mesh iterator range should be constructed and passed to the constraint constructor, and then used
    // by apply to iterate over the applicable mesh entities
    void apply(amsi::LAS * las, apf::Mesh * msh, pMesh prt, apf::Numbering * nm)
    {
      dim = msh->getDimension();
      std::vector<pEntity> rgns;
      amsi::getClassifiedEnts(prt,rgn,dim,std::back_inserter(rgns));
      elem_num = 0; //< reset element number counter.
      for(std::vector<pEntity>::iterator mrgn = rgns.begin(); mrgn != rgns.end(); mrgn++)
      {
        apf::MeshEntity * mnt = apf::castEntity(*mrgn);
        apf::MeshElement * mlm = apf::createMeshElement(msh,mnt);
        process(mlm);
        apf::NewArray<int> ids;
        apf::getElementNumbers(nm,mnt,ids);
        apf::DynamicMatrix & lG = getLambdaFirstVolDeriv();
        apf::DynamicMatrix & lH = getLambdaSecondVolDeriv();
        // Multiply appropriate terms by -1
        lG *= -1.0;
        amsi::assembleMatrix(las,nedofs,&ids[0],nedofs,&ids[0],&lH(0,0));
        amsi::assembleVector(las,nedofs,&ids[0],&lG(0,0));
	elem_num++;
      }
    }
    void inElement(apf::MeshElement * me)
    {
      ElementalSystem::inElement(me);
      fs = apf::getShape(f);
      es = fs->getEntityShape(apf::getMesh(f)->getType(apf::getMeshEntity(me)));
      d2Vdu2.setSize(nedofs,nedofs);
      dVdu.setSize(1,nedofs);
      d2Vdu2.zero();
      dVdu.zero();
    }
    bool includesBodyForces() { return true; }
    void atPoint(apf::Vector3 const &p, double w, double dV)
    {
      int & nen = nenodes; // = 4 (tets)
      // Jacobian of underlying mesh. (Initial configuration).
      apf::Matrix3x3 Jac;
      apf::getJacobian(me,p,Jac);
      // Note: Form of Jacobian in apf is
      // [dx/dr dx/ds dx/dt]    (x,y,z): physical domain coordinates.
      // [dy/dr dy/ds dy/dt]    (r,s,t): parent domain coordiantes.
      // [dz/dr dz/ds dz/dt]
      double wxdetjac = w * apf::getJacobianDeterminant(Jac,dim);
      // Jacobian of updated configuration
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
      // For Updated Lagrangian, the Jacobian of the updated coordinates are used
      // Note: that entires of Jacobian is hard coded for Linear tetrahedra elements.
      // TO DO: Generalize Jacobian for current configuration.
      apf::Matrix<3,3> J;
      J[0][0] = xyz(1,0) - xyz(0,0); // x2-x1
      J[0][1] = xyz(2,0) - xyz(0,0); // x3-x1
      J[0][2] = xyz(3,0) - xyz(0,0); // x4-x1
      J[1][0] = xyz(1,1) - xyz(0,1); // y2-y1
      J[1][1] = xyz(2,1) - xyz(0,1); // y3-y1
      J[1][2] = xyz(3,1) - xyz(0,1); // y4-y1
      J[2][0] = xyz(1,2) - xyz(0,2); // z2-z1
      J[2][1] = xyz(2,2) - xyz(0,2); // z3-z1
      J[2][2] = xyz(3,2) - xyz(0,2); // z4-z1
      double detJ = getDeterminant(J);
      // Determine derivatives of volume wrt to coordinates of updated configuration.
      apf::Matrix<3,3> dx(J), dy(J), dz(J);
      apf::Matrix<3,3> dxdx(J), dydy(J), dzdz(J);
      apf::Matrix<3,3> dxdy(J), dydz(J), dxdz(J);
      /*
        Derivatives of Shape function wrt parent domain coordinates (xi1,xi2,xi3)
        - local_grads: derivative of shape functions wrt to parent domain coordinates.
        - p          : integration point of mesh element wrt to parent domain coordinates.
        - mesh       : mesh on which the field f is defined. In this case, f is the apf_primary_field.
        Note: that for Linear Tetrahedral elements, derivatives of shape function wrt parent domain coodinates are constant.
      */
      apf::NewArray<apf::Vector3> local_grads;
      es->getLocalGradients(mesh, apf::getMeshEntity(me), p, local_grads);
      for (int a = 0; a < nen; a++) // loop through node a
      {
        for (int dxi = 0; dxi < dim; dxi++)
        {
          dx[0][dxi] = local_grads[a][dxi];
          dy[1][dxi] = local_grads[a][dxi];
          dz[2][dxi] = local_grads[a][dxi];
        }
        // Fill "G" vector
        dVdu(0,a*dim) = apf::getDeterminant(dx);
        dVdu(0,a*dim+1) = apf::getDeterminant(dy);
        dVdu(0,a*dim+2) = apf::getDeterminant(dz);
        for (int b = 0; b < nen; b++) // loop through node b
        {
          for (int dxi = 0; dxi < dim; dxi++)
          {
            // diagonal terms
            dxdx[0][dxi] = local_grads[a][dxi] + local_grads[b][dxi];
            dydy[1][dxi] = local_grads[a][dxi] + local_grads[b][dxi];
            dzdz[2][dxi] = local_grads[a][dxi] + local_grads[b][dxi];
            // off-diagonal terms
            dxdy[0][dxi] = local_grads[a][dxi];
            dxdy[1][dxi] = local_grads[b][dxi];
            dydz[1][dxi] = local_grads[a][dxi];
            dydz[2][dxi] = local_grads[b][dxi];
            dxdz[0][dxi] = local_grads[a][dxi];
            dxdz[2][dxi] = local_grads[b][dxi];
          }
          // Fill "H" matrix
          // diagonal terms
          d2Vdu2(a*dim,b*dim) = apf::getDeterminant(dxdx);
          d2Vdu2(a*dim+1,b*dim+1) = apf::getDeterminant(dydy);
          d2Vdu2(a*dim+2,b*dim+2) = apf::getDeterminant(dzdz);
          // off-diagonal terms
          d2Vdu2(a*dim,b*dim+1) = apf::getDeterminant(dxdy);
          d2Vdu2(a*dim,b*dim+2) = apf::getDeterminant(dxdz);
          d2Vdu2(a*dim+1,b*dim+2) = apf::getDeterminant(dydz);
          d2Vdu2(a*dim+1,b*dim) = apf::getDeterminant(dxdy) * -1.0;
          d2Vdu2(a*dim+2,b*dim) = apf::getDeterminant(dxdz) * -1.0;
          d2Vdu2(a*dim+2,b*dim+1) = apf::getDeterminant(dydz) * -1.0;
        }
      }
      /** Calculate Delta V */
      Vol = w * detJ;
      initVol = wxdetjac;
      /** Modifications to tangent-stiffness matrix */
      d2Vdu2 *= w;
      Lambda_d2Vdu2 = d2Vdu2;
      Lambda_d2Vdu2 *= Lambda;
      /** Modifications to force-vector */
      dVdu *= w;
      Lambda_dVdu = dVdu;
      Lambda_dVdu *= Lambda;
/*
      std::cout<<"ELEMENT="<<elem_num<<std::endl;
      std::cout<<"d2Vdu2:"<<std::endl;
      std::cout<<d2Vdu2<<std::endl;

      std::cout<<"dVdu:"<<std::endl;
      std::cout<<dVdu<<std::endl;
*/      
      /** Augmented Lagrangian method with volume difference constraint.*/
      DeltaV = Vol - initVol;
      apf::DynamicMatrix A(nedofs,nedofs);
      apf::DynamicMatrix GT(nedofs,1);
      apf::transpose(dVdu,GT);
      apf::multiply(GT,dVdu,A);
//      std::cout<<"GG^T"<<std::endl;
//      std::cout<<A<<std::endl;
      apf::DynamicMatrix VH(nedofs,nedofs);
      VH = d2Vdu2;
      VH *= DeltaV;
//      std::cout<<"DeltaV = "<<DeltaV<<std::endl;
//      std::cout<<"DeltaV H"<<std::endl;
//      std::cout<<VH<<std::endl;
      A += VH;
      A *= Beta;
      Lambda_d2Vdu2 += A;
//      std::cout<<"LHS mod term"<<std::endl;
//      std::cout<<Lambda_d2Vdu2<<std::endl;
      apf::DynamicMatrix BVG(1,nedofs);
      BVG = dVdu;
      BVG *= DeltaV;
      BVG *= Beta;
      Lambda_dVdu += BVG;
	
      /** Augmented Lagrangian method with relative volume constraint (V-Vp)/Vp 
      // prevVol is updated subsequently in NonLinTissue::updatePrevVolumes(); 
      double prevV = 0.0;
      if (std::abs(prevVol[elem_num]) < 1e-10)
      {
	DeltaV = Vol - initVol;
	prevV = initVol;
      }
      else
      {
	DeltaV = Vol - prevVol[elem_num];
	prevV = prevVol[elem_num];
      }
      prevV = 1.0; ///< This will remove relative difference implementation.
      Lambda_dVdu /= prevV; ///<relative volume constraint modification on lambda G
      Lambda_d2Vdu2 /= prevV; ///<relative volume constraint modification on lambda H
      apf::DynamicMatrix A(nedofs,nedofs);
      apf::DynamicMatrix GT(nedofs,1);
      apf::transpose(dVdu,GT);
      apf::multiply(GT,dVdu,A);
//      std::cout<<"GG^T"<<std::endl;
//      std::cout<<A<<std::endl;
      A /= prevV; ///<relative volume constraint modification on GG^T term.
//      std::cout<<"prevV="<<prevV<<std::endl;
//      std::cout<<"GG^T/preV"<<std::endl;
//      std::cout<<A<<std::endl;
      apf::DynamicMatrix VH(nedofs,nedofs);
      VH = d2Vdu2;
      VH *= DeltaV;
//      std::cout<<"DeltaV = "<<DeltaV<<std::endl;
//      std::cout<<"DeltaV H"<<std::endl;
//      std::cout<<VH<<std::endl;
      VH /= prevV; ///<relative volume constraint modification on VH term.
      A += VH;
      A *= Beta;
      A /= prevV; ///<relative volume constraint modification on entire A term.
      Lambda_d2Vdu2 += A;
//      std::cout<<"LHS mod term"<<std::endl;
//      std::cout<<Lambda_d2Vdu2<<std::endl;
      apf::DynamicMatrix BVG(1,nedofs);
      BVG = dVdu;
      BVG *= DeltaV;
      BVG *= Beta;
      BVG /= std::pow(prevV,2); ///<relative volume constraint modification on RHS term.
      Lambda_dVdu += BVG;
      */
	
      /** Augmented Lagrangian method with relative volume constraint (V-V0)/V0 
      DeltaV = Vol - initVol;
      Lambda_dVdu /= initVol; ///<relative volume constraint modification on lambda G
      Lambda_d2Vdu2 /= initVol; ///<relative volume constraint modification on lambda H
      apf::DynamicMatrix A(nedofs,nedofs);
      apf::DynamicMatrix GT(nedofs,1);
      apf::transpose(dVdu,GT);
      apf::multiply(GT,dVdu,A);
      A /= initVol; ///<relative volume constraint modification on GG^T term.
      apf::DynamicMatrix VH(nedofs,nedofs);
      VH = d2Vdu2;
      VH *= DeltaV;
      VH /= initVol; ///<relative volume constraint modification on VH term.
      A += VH;
      A *= Beta;
      A /= initVol; ///<relative volume constraint modification on entire A term.
      Lambda_d2Vdu2 += A;
      apf::DynamicMatrix BVG(1,nedofs);
      BVG = dVdu;
      BVG *= DeltaV;
      BVG *= Beta;
      BVG /= std::pow(initVol,2); ///<relative volume constraint modification on RHS term.
      Lambda_dVdu += BVG;
      */
/*
      std::cout<<"Lambda d2Vdu2:"<<std::endl;
      std::cout<<Lambda_d2Vdu2<<std::endl;

      std::cout<<"Lamdba dVdu:"<<std::endl;
      std::cout<<Lambda_dVdu<<std::endl;
*/

    }
    apf::DynamicMatrix& getFirstVolDeriv(){return dVdu;}
    apf::DynamicMatrix& getSecondVolDeriv(){return d2Vdu2;}
    apf::DynamicMatrix& getLambdaFirstVolDeriv(){return Lambda_dVdu;}
    apf::DynamicMatrix& getLambdaSecondVolDeriv(){return Lambda_d2Vdu2;}
    double getDeltaV(){return DeltaV;}
    double getLambda(){return Lambda;}
    double getBeta(){return Beta;}
    double getVol(){return Vol;}
    double getinitVol(){return initVol;}
    void setLambda(double Delta_Lambda) { Lambda = Delta_Lambda; }
    void setBeta(double factor) { Beta = factor; }
    void setPrevVol(int idx, double vol){prevVol[idx] = vol;}
    void getPrevVols(std::vector<double> & v){v = prevVol;}
    pGEntity getRegion() { return rgn; }
  private:
    bool local;
    int dim;
    int dof;
    pGEntity rgn;
    apf::FieldShape * fs;
    apf::EntityShape * es;
    apf::DynamicMatrix d2Vdu2;
    apf::DynamicMatrix dVdu;
    apf::DynamicMatrix Lambda_d2Vdu2;
    apf::DynamicMatrix Lambda_dVdu;
    double DeltaV;
    double Lambda; // Lagrange multiplier for Volume Preservation.
    double Beta;   // Penalty parameter for Augmented Lagrangian Method.
    double Vol;
    double initVol;
    std::vector<double> prevVol;
    int elem_num; // keeps count of element number within region.
  }; // End of VolumeConstraintADMM
}
#endif
