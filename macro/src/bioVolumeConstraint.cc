#include "bioVolumeConstraint.h"
#include <amsiCasters.h>
#include <amsiFEA.h>
#include <apfMeasure.h>
#include <gmi.h>
#include <simClassified.h>
namespace bio
{
  VolumeConstraint * buildVolumeConstraint(pACase cs, pANode nd, apf::Numbering * un)
  {
    VolumeConstraint * cnst = NULL;
    std::vector<apf::ModelEntity*> mdl_ents;
    amsi::getAssociatedModelItems(cs,nd,std::back_inserter(mdl_ents));
    std::cout << "Volume constraint discovered, effects model entities : ";
    gmi_model * mdl = apf::getMesh(un)->getModel();
    for(auto mdl_ent = mdl_ents.begin(); mdl_ent != mdl_ents.end(); ++mdl_ent)
      std::cout << mdl->ops->tag(mdl,(gmi_ent*)*mdl_ent) << " ";
    std::cout << std::endl;
    pANode mthd = AttNode_childByType(nd,"constraint type");
    pANode vrsn = AttNode_childByType(nd,"version");
    assert(!AttInfoInt_isExpression((pAttInfoInt)vrsn));
    int version = AttInfoInt_value((pAttInfoInt)vrsn);
    char * tp = AttNode_imageClass(mthd);
    if(std::string("penalty method").compare(tp) == 0)
    {
      pANode bt = AttNode_childByType(mthd,"beta");
      assert(!AttInfoDouble_isExpression((pAttInfoDouble)bt));
      double beta = AttInfoDouble_value((pAttInfoDouble)bt);
      switch(version)
      {
      case 0: // region
        //cnst = new PenaltyConstraint_Volume(mdl_ents.begin(),mdl_ents.end(),un,beta);
        std::cerr << "WARNING: unimplemented constraint: Penalty method volume constraint." << std::endl;
        break;
      case 1: // surface
        cnst = new PenaltyConstraint_VolumeSurface(mdl_ents.begin(),mdl_ents.end(),un,beta);
        break;
      }
    }
    else if (std::string("lagrange multiplier").compare(tp) == 0)
    {
      pANode bt = AttNode_childByType(mthd,"beta");
      assert(!AttInfoDouble_isExpression((pAttInfoDouble)bt));
      double beta = AttInfoDouble_value((pAttInfoDouble)bt);
      switch(version)
      {
      case 0: // region
        cnst = new LagrangeConstraint_Volume(mdl_ents.begin(),mdl_ents.end(),un,beta);
        break;
      case 1: // surface
        cnst = new LagrangeConstraint_VolumeSurface(mdl_ents.begin(),mdl_ents.end(),un,0.0,beta);
        break;
      }
    }
    std::cout << "\tInitial volume : " << cnst->getVolume() << std::endl;
    Sim_deleteString(tp);
    return cnst;
  }
  double VolumeConstraint::calcVolume()
  {
    return amsi::measureDisplacedModelEntities(mdl_ents.begin(),mdl_ents.end(),apf::getField(nm));
  }
  void VolumeConstraint::iter()
  {
    prev_vol = vol;
    vol = calcVolume();
    _iter();
  }
  void VolumeConstraint::step()
  {
    load_vol = calcVolume();
    prev_vol = vol = load_vol;
    _step();
  }
  void VolumeConstraint::inElement(apf::MeshEntity * m)
  {
    apf::Field * fld = apf::getField(nm);
    me = apf::createMeshElement(apf::getMesh(fld),m);
    e = apf::createElement(fld,me);
    nen = apf::countNodes(e);
    nedofs = nen * apf::countComponents(fld);
    _inElement(m);
  }
  void VolumeConstraint::outElement()
  {
    apf::destroyMeshElement(me);
    apf::destroyElement(e);
  }
  void LagrangeConstraint_Volume::_iter()
  {
    //lambda += beta * ((vol-prev_vol)/prev_vol);
    //amsi::Log(lg) << "-1, -1," << lambda << ", " << beta << std::endl;
  }
  void LagrangeConstraint_Volume::_step()
  {
    lambda += beta * ((vol-prev_vol)/prev_vol);
  }
  void LagrangeConstraint_Volume::_inElement(apf::MeshEntity * me)
  {
    d2Vdu2.setSize(nedofs,nedofs);
    dVdu.setSize(1,nedofs);
    d2Vdu2.zero();
    dVdu.zero();
  }
  void LagrangeConstraint_Volume::apply(amsi::LAS * las)
  {
    apf::Mesh * msh = apf::getMesh(apf::getField(nm));
    int dim = msh->getDimension();
    for(auto mdl_ent = mdl_ents.begin(); mdl_ent != mdl_ents.end(); ++mdl_ent)
    {
      for(auto ent = amsi::beginClassified(msh,*mdl_ent,dim); ent != amsi::endClassified(ent); ++ent)
      {
        inElement(*ent);
        apf::MeshElement * mlm = apf::createMeshElement(msh,*ent);
        process(mlm);
        apf::NewArray<int> ids;
        apf::getElementNumbers(nm,*ent,ids);
        apf::DynamicMatrix ldVdu = dVdu;
        ldVdu *= -1.0 * lambda;
        apf::DynamicMatrix ld2Vdu2 = d2Vdu2;
        ld2Vdu2 *= -1.0 * lambda;
        //double ndv = -1.0 * (vol - prev_vol);
        amsi::assembleMatrix(las,nedofs,&ids[0],nedofs,&ids[0],&ld2Vdu2(0,0));
        amsi::assembleMatrix(las,nedofs,&ids[0],1,&dof,&dVdu(0,0));
        //amsi::assembleMatrix(las,1,&dof,nedofs,&ids[0],&dVdu(0,0)); // Transpose of G.??
        amsi::assembleVector(las,nedofs,&ids[0],&ldVdu(0,0));
        //amsi::assembleVector(las,1,&dof,&ndv);
        apf::destroyMeshElement(mlm);
      }
    }
  }
  void LagrangeConstraint_Volume::atPoint(apf::Vector3 const & p, double w, double dV)
  {
    // Jacobian of underlying mesh. (Initial configuration).
    apf::Matrix3x3 Jac;
    apf::getJacobian(me,p,Jac);
    // Note: Form of Jacobian in apf is
    // [dx/dr dx/ds dx/dt]    (x,y,z): physical domain coordinates.
    // [dy/dr dy/ds dy/dt]    (r,s,t): parent domain coordiantes.
    // [dz/dr dz/ds dz/dt]
    // Jacobian of updated configuration
    // 1. Get coordinates on underlying mesh
    apf::Field * u = apf::getField(nm);
    apf::Mesh * msh = apf::getMesh(u);
    int dim = msh->getDimension();
    //double wxdetjac = w * apf::getJacobianDeterminant(Jac,dim);
    apf::Field * xyz = msh->getCoordinateField();
    apf::Element * xyz_elem = apf::createElement(xyz,me);
    apf::NewArray<apf::Vector3> xyz_nds;
    apf::getVectorNodes(xyz_elem,xyz_nds);
    // 2. Get coordinates from apf_primary_field (passed in), which contains the accumulated displacement
    apf::NewArray<apf::Vector3> u_nds;
    apf::getVectorNodes(e,u_nds);
    // 3. Calculate current coordinates
    apf::DynamicMatrix xyz_u(nen,dim);
    for (int ii = 0; ii < nen; ii++)
      for (int jj = 0; jj < dim; jj++)
        xyz_u(ii,jj) = xyz_nds[ii][jj] + u_nds[ii][jj];
    // For Updated Lagrangian, the Jacobian of the updated coordinates are used
    // Note: that entires of Jacobian is hard coded for Linear tetrahedra elements.
    // TO DO: Generalize Jacobian for current configuration.
    apf::Matrix<3,3> J;
    J[0][0] = xyz_u(1,0) - xyz_u(0,0); // x2-x1
    J[0][1] = xyz_u(2,0) - xyz_u(0,0); // x3-x1
    J[0][2] = xyz_u(3,0) - xyz_u(0,0); // x4-x1
    J[1][0] = xyz_u(1,1) - xyz_u(0,1); // y2-y1
    J[1][1] = xyz_u(2,1) - xyz_u(0,1); // y3-y1
    J[1][2] = xyz_u(3,1) - xyz_u(0,1); // y4-y1
    J[2][0] = xyz_u(1,2) - xyz_u(0,2); // z2-z1
    J[2][1] = xyz_u(2,2) - xyz_u(0,2); // z3-z1
    J[2][2] = xyz_u(3,2) - xyz_u(0,2); // z4-z1
    //double detJ = getDeterminant(J);
    // Determine derivatives of volume wrt to coordinates of updated configuration.
    apf::Matrix<3,3> dx(J), dy(J), dz(J);
    apf::Matrix<3,3> dxdx(J), dydy(J), dzdz(J);
    apf::Matrix<3,3> dxdy(J), dydz(J), dxdz(J);
    // Derivatives of Shape function wrt parent domain coordinates (xi1,xi2,xi3)
    //   - local_grads: derivative of shape functions wrt to parent domain coordinates.
    //   - p          : integration point of mesh element wrt to parent domain coordinates.
    //   - mesh       : mesh on which the field f is defined. In this case, f is the apf_primary_field.
    //   Note: that for Linear Tetrahedral elements, derivatives of shape function
    //         wrt parent domain coodinates are constant.
    apf::NewArray<apf::Vector3> grads;
    apf::getShapeGrads(e,p,grads);
    for (int a = 0; a < nen; a++) // loop through node a
    {
      for (int dxi = 0; dxi < dim; dxi++)
      {
        dx[0][dxi] = grads[a][dxi];
        dy[1][dxi] = grads[a][dxi];
        dz[2][dxi] = grads[a][dxi];
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
          dxdx[0][dxi] = grads[a][dxi] + grads[b][dxi];
          dydy[1][dxi] = grads[a][dxi] + grads[b][dxi];
          dzdz[2][dxi] = grads[a][dxi] + grads[b][dxi];
          // off-diagonal terms
          dxdy[0][dxi] = grads[a][dxi];
          dxdy[1][dxi] = grads[b][dxi];
          dydz[1][dxi] = grads[a][dxi];
          dydz[2][dxi] = grads[b][dxi];
          dxdz[0][dxi] = grads[a][dxi];
          dxdz[2][dxi] = grads[b][dxi];
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
    //vol = w * detJ;
    //init_vol = wxdetjac;
    //double delta_v = vol - init_vol;
    // Modifications to tangent-stiffness matrix
    d2Vdu2 *= w;
    // Modifications to force-vector.
    dVdu *= w;
    // Augmented Lagrangian method
    /*
    apf::DynamicMatrix A(nedofs,nedofs);
    apf::DynamicMatrix GT(nedofs,1);
    apf::transpose(dVdu,GT);
    apf::multiply(GT,dVdu,A);
    apf::DynamicMatrix VH(nedofs,nedofs);
    VH = d2Vdu2;
    VH *= delta_v;
    A += VH;
    A *= beta;
    lambda_d2Vdu2 += A;
    apf::DynamicMatrix BVG(1,nedofs);
    BVG = dVdu;
    BVG *= delta_v;
    BVG *= beta;
    lambda_dVdu += BVG;
    */
  }
  void calcdVdu(apf::DynamicMatrix & dVdu, apf::Vector3 const & pt0, apf::Vector3 const & pt1, apf::Vector3 const & pt2)
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
  }
  void LagrangeConstraint_VolumeSurface::_inElement(apf::MeshEntity * me)
  {
    dVdu.setSize(1,nedofs);
    dVdu.zero();
  }
  void LagrangeConstraint_VolumeSurface::_iter()
  {
    lambda += beta * ((vol - prev_vol)/prev_vol);
  }
  void LagrangeConstraint_VolumeSurface::apply(amsi::LAS * las)
  {
    int dm = msh->getDimension() - 1;
    for(auto rgn = mdl_ents.begin(); rgn != mdl_ents.end(); ++rgn)
    {
      crt_rgn = *rgn;
      gmi_model * mdl = msh->getModel();
      gmi_set * fcs = mdl->ops->adjacent(mdl,(gmi_ent*)crt_rgn,2);
      for(int ii = 0; ii < fcs->n; ++ii)
      {
        crt_fc = reinterpret_cast<apf::ModelEntity*>(fcs->e[ii]);
        for(auto ent = amsi::beginClassified(msh,crt_fc,dm); ent != amsi::endClassified(ent); ++ent)
        {
          inElement(*ent);
          apf::MeshElement * mlm = apf::createMeshElement(msh,*ent);
          process(mlm);
          apf::NewArray<int> ids;
          apf::getElementNumbers(nm,*ent,ids);
          apf::DynamicMatrix ndVdu = dVdu;
          ndVdu *= -1.0;
          amsi::assembleVector(las,nedofs,&ids[0],&ndVdu(0,0));
          apf::destroyMeshElement(mlm);
        }
      }
      gmi_free_set(fcs);
    }
  }
  void LagrangeConstraint_VolumeSurface::atPoint(apf::Vector3 const & p, double w, double dV)
  {
    int dim = msh->getDimension();
    apf::Field * xyz = msh->getCoordinateField();
    apf::Element * xyz_elem = apf::createElement(xyz,me);
    apf::NewArray<apf::Vector3> xyz_nds;
    apf::getVectorNodes(xyz_elem,xyz_nds);
    apf::NewArray<apf::Vector3> u_nds;
    apf::getVectorNodes(e,u_nds);
    apf::Vector3 xyz_u[3];
    //apf::DynamicMatrix xyz_u(nen,dim);
    for (int ii = 0; ii < nen; ii++)
      for (int jj = 0; jj < dim; jj++)
        xyz_u[ii][jj] = xyz_nds[ii][jj] + u_nds[ii][jj];
    calcdVdu(dVdu,xyz_u[0], xyz_u[1], xyz_u[2]);
    int sd = amsi::side(crt_rgn,msh,apf::getMeshEntity(me));
    assert(sd != 0); // if both adjacent regions are in the same model region, we haven't implemented a way to handle that, (should basically ignore)
    double dVol = vol - load_vol; // should match volume convergence operator related to this constraint ??
    dVdu *= (lambda * sd) + (dVol * beta / prev_vol);
    /*
    apf::DynamicMatrix BVG = dVdu;
    dVdu *= lambda; // <-- division by volume is contained in lambda via update of constraint.
    // Modifications to force-vector.
    // Augmented Lagrangian method:
    // A   = beta * (GG^T + dVol * H)
    // VH  = dVol * H
    // BVG = beta * delta_v * G
    // lambda * H modified to lambda * H + beta * (GG^T + dVol * H)
    // Second Derivative is not used (see N2P106 for details).
    BVG *= dVol;
    BVG *= beta;
    BVG /= prev_vol;
    dVdu += BVG; // lambda * n * dVdu + (dVol * beta / pv) * dVdu;
    */
  }
  void PenaltyConstraint_VolumeSurface::_inElement(apf::MeshEntity * me)
  {
    dVdu.setSize(1,nedofs);
    dVdu.zero();
  }
  void PenaltyConstraint_VolumeSurface::apply(amsi::LAS * las)
  {
    int dm = msh->getDimension() - 1;
    for(auto rgn = mdl_ents.begin(); rgn != mdl_ents.end(); ++rgn)
    {
      crt_rgn = *rgn;
      gmi_model * mdl = msh->getModel();
      gmi_set * fcs = mdl->ops->adjacent(mdl,(gmi_ent*)crt_rgn,2);
      for(int ii = 0; ii < fcs->n; ++ii)
      {
        crt_fc = reinterpret_cast<apf::ModelEntity*>(fcs->e[ii]);
        for(auto ent = amsi::beginClassified(msh,crt_fc,dm); ent != amsi::endClassified(ent); ++ent)
        {
          inElement(*ent);
          apf::MeshElement * mlm = apf::createMeshElement(msh,*ent);
          process(mlm);
          apf::NewArray<int> ids;
          apf::getElementNumbers(nm,*ent,ids);
          apf::DynamicMatrix ndVdu = dVdu;
          ndVdu *= -1.0;
          amsi::assembleVector(las,nedofs,&ids[0],&ndVdu(0,0));
          apf::destroyMeshElement(mlm);
        }
      }
      gmi_free_set(fcs);
    }
  }
  // this is really only safe for mesh faces with a single integration point atm
  void PenaltyConstraint_VolumeSurface::atPoint(apf::Vector3 const &p, double w, double dV)
  {
    int dim = msh->getDimension();
    apf::Field * xyz = msh->getCoordinateField();
    apf::Element * xyz_elem = apf::createElement(xyz,me);
    apf::NewArray<apf::Vector3> xyz_nds;
    apf::getVectorNodes(xyz_elem,xyz_nds);
    apf::NewArray<apf::Vector3> u_nds;
    apf::getVectorNodes(e,u_nds);
    apf::Vector3 xyz_u[3];
    //apf::DynamicMatrix xyz_u(nen,dim);
    for (int ii = 0; ii < nen; ii++)
      for (int jj = 0; jj < dim; jj++)
        xyz_u[ii][jj] = xyz_nds[ii][jj] + u_nds[ii][jj];
    calcdVdu(dVdu,xyz_u[0], xyz_u[1], xyz_u[2]);
    int sd = amsi::side(crt_rgn,msh,apf::getMeshEntity(me));
    assert(sd != 0); // if both adjacent regions are in the same model region, we haven't implemented a way to handle that, (should basically ignore)
    double dVol = vol - load_vol;
    dVdu *= (sd * dVol * beta / prev_vol);
    /*
    apf::DynamicMatrix BVG = dVdu;
    dVdu *= lambda; // <-- division by volume is contained in lambda via update of constraint.
    // Modifications to force-vector.
    // Augmented Lagrangian method:
    // A   = beta * (GG^T + dVol * H)
    // VH  = dVol * H
    // BVG = beta * delta_v * G
    // lambda * H modified to lambda * H + beta * (GG^T + dVol * H)
    // Second Derivative is not used (see N2P106 for details).
    BVG *= dVol;
    BVG *= beta;
    BVG /= prev_vol;
    dVdu += BVG; // lambda * n * dVdu + (dVol * beta / pv) * dVdu;
    */
  }
  /*
  void PenaltyConstraint_VolumeSurface::calcd2Vdu2(apf::Vector3 const & pt0, apf::Vector3 const & pt1, apf::Vector3 const & pt2)
  {
    d2Vdu2.zero();
    // submat1 = [ 0  -z3  y3]
    // [ z3   0 -x3]
    // [-y3  x3   0]
    // submat2 = [  0  z2 -y2]
    // [-z2   0  x2]
    // [ y2 -x2   0]
    // submat3 = [  0 -z1  y1]
    // [ z1   0 -x1]
    // [-y1  x1   0]
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
    for (int ii = 0; ii < nenodes; ii++)
      for (int jj = 0; jj < nenodes; jj++)
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
  }
  */
  /*
  VolumeConstraintADMM::VolumeConstraintADMM(apf::ModelEntity * me, apf::Field * fld, int o)
    : VolumeConstraint(me,fld,o)
    , prev_vols()
    , elem_num(0)
  { }
  void VolumeConstraintADMM::apply(amsi::LAS * las, apf::Mesh * msh, pMesh prt, apf::Numbering * nm)
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
  void VolumeConstraintADMM::atPoint(apf::Vector3 const &p, double w, double dV)
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
    // Derivatives of Shape function wrt parent domain coordinates (xi1,xi2,xi3)
    //  - local_grads: derivative of shape functions wrt to parent domain coordinates.
    //  - p          : integration point of mesh element wrt to parent domain coordinates.
    //  - mesh       : mesh on which the field f is defined. In this case, f is the apf_primary_field.
    //  Note: that for Linear Tetrahedral elements, derivatives of shape function wrt parent domain coodinates are constant.
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
    init_Vol = wxdetjac;
    // Modifications to tangent-stiffness matrix 
    d2Vdu2 *= w;
    lambda_d2Vdu2 = d2Vdu2;
    lambda_d2Vdu2 *= lambda;
    // Modifications to force-vector 
    dVdu *= w;
    lambda_dVdu = dVdu;
    lambda_dVdu *= lambda;
    // Augmented Lagrangian method with volume difference constraint.
    delta_v = Vol - init_vol;
    apf::DynamicMatrix A(nedofs,nedofs);
    apf::DynamicMatrix GT(nedofs,1);
    apf::transpose(dVdu,GT);
    apf::multiply(GT,dVdu,A);
    apf::DynamicMatrix VH(nedofs,nedofs);
    VH = d2Vdu2;
    VH *= delta_v;
    A += VH;
    A *= beta;
    lambda_d2Vdu2 += A;
    apf::DynamicMatrix BVG(1,nedofs);
    BVG = dVdu;
    BVG *= delta_v;
    BVG *= beta;
    lambda_dVdu += BVG;
  }
 */
}
