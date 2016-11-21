#include "NonLinFibMtx.h"
#include "NonFiniteElement.h"

#include <ConvenienceFunctions.h> // Identity3x3

#include <memory.h>
#include <vector>
#include <list>
#include <map>
#include <iostream>

namespace Biotissue {

// this likely needs an intensive rewrite... 
    void NonLinFibMtx::AverageFiberStress(apf::Matrix3x3 & stress)
    {
      /*
      apf::Matrix3x3 volume_stress;

      int global, local, offset;
      GetDOFInfo(global, local, offset);
      ////////////////////////////////////////
      ComputeFiberForceVector();
      ////////////////////////////////////////
      // get the force vector ...
  
      std::list<pVertex> bdryMVtx;
      GFIter gfiter = GM_faceIter(model);
      pGFace gf;
      while (gf=GFIter_next(gfiter)) {
	int ftag = GEN_tag(gf);
	getMeshVtxOnGeomFace(ftag,bdryMVtx);
	list<pVertex>::iterator Iter = bdryMVtx.begin();
	list<pVertex>::iterator IterEnd = bdryMVtx.end();
	for(; Iter != IterEnd; Iter++) {
	  pVertex vtx = (pVertex)*Iter;
	  double xyz[3], val[3], frc[3];
	  V_coord (vtx, xyz);
	  for(int icomp=0; icomp<3; icomp++){
	    SCOREC::Field::DofKey key1((pEntity)vtx, 1, GetField()->getTag(), icomp);
	    SCOREC::Field::Dof *dof1 = DofManager::Instance()->getDof(key1);
	    assert(dof1);
	    int index = dof1->getGlobalNumbering(); 
	    frc[icomp] = RHS[index-1];   // check start number and the fiber force vector
	    val[icomp] = dof1->getValue(1);
	    xyz[icomp] += val[icomp];
	  }
      
	  for(int i=0; i<3; i++)
	    for(int j=0; j<3; j++)
	      volStress[i][j] += xyz[j]*frc[i];
	}
	bdryMVtx.clear();
      }
  
      computeUnbalancedTermFromFiber(RHS, numEq);
  
      delete []RHS;	
      GFIter_delete(gfiter);
      VIter_delete(viter);
  
      sigmaF_[0] = volStress[0][0]; sigmaF_[1] = volStress[1][1]; sigmaF_[2] = volStress[2][2];
      sigmaF_[3] = 0.5*(volStress[0][1]+volStress[1][0]);
      sigmaF_[4] = 0.5*(volStress[1][2]+volStress[2][1]);
      sigmaF_[5] = 0.5*(volStress[2][0]+volStress[0][2]);
  
      //compute the current volume of RVE
      double det = SCOREC::Util::det(dfmGrad_);
      //original size of RVE is 1 unit
      double undeformed_vol=1;
      double vol = det * undeformed_vol;
  
      for(int i=0; i<6; i++)
	sigmaF_[i] = sigmaF_[i]/vol;
      */  
    }
  
  
  void NonLinFibMtx::UnbalancedFiberForce(apf::Vector3 & term)
  {
    /*
      GFIter gfiter = GM_faceIter(getModel());
      pGFace gf;
      while(gf=GFIter_next(gfiter)) {
	int ftag = GEN_tag(gf);
	mVector normal;
	computeRVEFaceNormalFromDfmGrad(gf, normal);
	list<pVertex> bdryMVtx;
	getMeshVtxOnGeomFace(ftag,bdryMVtx);
	list<pVertex>::iterator iter = bdryMVtx.begin();
	list<pVertex>::iterator iterEnd = bdryMVtx.end();		
	int numBdryV = bdryMVtx.size();
	for(;iter != iterEnd; iter++) {
	  pVertex bdryV = *iter;
	  int numAdjEdge = V_numEdges(bdryV);
	  for(int j=0; j<numAdjEdge; j++) {
	    pEdge edge = V_edge(bdryV, j);
	    gType entType = E_whatInType(edge);
	    int adjFibers = 0;
	    if(entType == Gedge) {
	      adjFibers++;
	      double frc[3];
	      for(int icomp=0; icomp<3; icomp++){
		SCOREC::Field::DofKey key1((pEntity)bdryV, 1, GetField()->getTag(), icomp);
		SCOREC::Field::Dof *dof1 = DofManager::Instance()->getDof(key1);
		assert(dof1);
		int index = dof1->getGlobalNumbering();
		frc[icomp] = RHS[index-1];   // check start number and the fiber force vector
	      }
	      mTensor2 s;
	      for(int icomp=0; icomp<3; icomp++)
		for(int jcomp=0; jcomp<3; jcomp++)
		  s(icomp,jcomp) = frc[icomp]*normal(jcomp);		
					

	      mVector UxN = SCOREC::Util::operator*(normal, dispGrad_);
	      mVector Q = SCOREC::Util::operator*(UxN, s);
			
	      for(int jcomp=0; jcomp<3; jcomp++)
		QF_[jcomp] += Q(jcomp);				
	    }
	  }
	}
	bdryMVtx.clear();
      }

      //compute the current volume of RVE
      double det = SCOREC::Util::det(dfmGrad_);
      //original size of RVE is 1 unit
      double undeformed_vol=1;
      double vol = det * undeformed_vol;

      for(int i=0; i<3; i++)
	QF_[i] = QF_[i]/vol;


      GFIter_delete(gfiter);
      */
    }


/*
    void NonLinFibMtx::getVertexDisplacementOnEdge(pEdge edge, double u[6])
    {
      pPList vlist = PList_new();

      pVertex v0 = E_vertex(edge,0);
      pVertex v1 = E_vertex(edge,1);	

      for (int icomp=0; icomp<3; icomp++) {
	SCOREC::Field::DofKey key((pEntity)v0, 1, GetField()->getTag(), icomp);
	SCOREC::Field::Dof *dof = DofManager::Instance()->getDof(key);
	assert(dof);
	// Get the accumulated displacement
	u[icomp] = dof->getValue(1);

	SCOREC::Field::DofKey key1((pEntity)v1, 1, GetField()->getTag(), icomp);
	SCOREC::Field::Dof *dof1 = DofManager::Instance()->getDof(key1);
	assert(dof1);
	// Get the accumulated displacement
	u[3+icomp] = dof1->getValue(1);

      }

    }
*/
// calculates the spatial displacement tensor $\nabla_{\mathbf x} \mathbf U = \mathbf I - \mathbf F^{-1}$ where $\mathbf F$ is the deformation gradient tensor
    void NonLinFibMtx::SpatialDispGrad(apf::Matrix3x3 & disp_grad)
    {
      apf::Matrix3x3 def_grad_inverse = apf::invert(deformation_gradient);
      disp_grad = Identity3x3 - def_grad_inverse;
    }

  } // end of namespace











