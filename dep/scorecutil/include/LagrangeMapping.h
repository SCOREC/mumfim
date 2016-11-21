/*
   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of Trellis written and maintained by the
   Scientific Computation Research Center (SCOREC) at Rensselaer Polytechnic
   Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.

   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA
*/

#ifndef H_LagrangeMapping
#define H_LagrangeMapping

#include "Mapping.h"

  class LagrangeMapping : public Mapping { 
    //std::vector<mPoint> knots; - move to the base class
    public:
      // TRELLIS_MTYPE elementType; - this is the reduncency of the elementType of the base class Mapping.h
      virtual ~LagrangeMapping(){};
      LagrangeMapping(pEntity me);
      LagrangeMapping(pEntity me, LevelSetFunction _levelSetFunc, int condition);
      LagrangeMapping(pEntity me, std::vector<mPoint> *input_knots);
      inline std::vector<mPoint> * getKnots(){return &knots;}
      virtual void normalVector (pEntity border , 
				 double u, double v, double w, 
				 mVector &n) const;
      virtual void eval(double u, double v, double w,
			double &x, double &y, double &z) const;
      
      virtual void deval(double u, double v, double w,
			 double &dxdu, double &dydu, double &dzdu,
			 double &dxdv, double &dydv, double &dzdv,
			 double &dxdw, double &dydw, double &dzdw) const;
      virtual void boundingBox (  mPoint &min, mPoint &max ) const;
      
      double detJacDeriv(double &u, double &v, double &w, double *jDeriv, mVector *dus);
      int geomOrder() const;
      int order() const;
    protected:
      void   GeomShapeFunction          (double u, double v, double w, double * ) const;
      void   GradGeomShapeFunction      (double u, double v ,double w , mVector*) const;

      // need to move all of the topo shape function here
      void GeomShapeFunctionLine (double u, double v, double w, double *result) const;
      void GeomShapeFunctionTri ( double u, double v, double dum, double *x) const;
      void GeomShapeFunctionQuad (double u, double v, double dum, double *x) const;
      void GeomShapeFunctionTet (double u, double v, double w, double *x) const;
      void GeomShapeFunctionHex (double u, double v, double w, double *x) const;
      void GeomShapeFunctionPrism (double u, double v, double w, double *x) const;
      void GeomShapeFunctionPyr (double u, double v, double w, double *x) const;
      
      void GradGeomShapeFunctionLine (double u, double v, double w, mVector *grad) const;
      void GradGeomShapeFunctionTri (double u, double v, double dum , mVector *result) const;
      void GradGeomShapeFunctionQuad (double u, double v, double dum , mVector *result) const;
      void GradGeomShapeFunctionTet (double r, double s, double t, mVector *result) const;
      void GradGeomShapeFunctionHex (double u, double v, double w , mVector *grad) const;
      void GradGeomShapeFunctionPrism (double u, double v, double w , mVector *grad) const;
      void GradGeomShapeFunctionPyr (double u, double v, double w , mVector *result) const;
      
  };

  class ConstantLagrangeMapping : public LagrangeMapping
    {
      double det;
      mTensor2 jac;
    public:
      virtual ~ConstantLagrangeMapping(){}
      ConstantLagrangeMapping(pEntity me);
      virtual double jacInverse(double x, double y, double z, 
				mTensor2&) const;
      virtual double detJac(double u, double v, double w) const;
    };

  class CylindricalCoordinatesLagrangeMapping : public LagrangeMapping
    {
    public:
      virtual ~CylindricalCoordinatesLagrangeMapping(){}
      CylindricalCoordinatesLagrangeMapping(pEntity me);
      virtual double jacInverse(double x, double y, double z, 
				mTensor2&) const;
      virtual double detJac(double u, double v, double w) const;
      virtual int order() const;
    };
  
  class RegularCubeLagrangeMapping : public LagrangeMapping
    {
      double x0,y0,z0;
      double dx,dy,dz;
    public:
      virtual ~RegularCubeLagrangeMapping(){}
      RegularCubeLagrangeMapping(pEntity me);
      virtual void normalVector (pEntity border , 
				 double u, double v, double w, 
				 mVector &n) const;
      virtual void eval(double u, double v, double w,
			double &x, double &y, double &z) const;
      
      virtual void deval(double u, double v, double w,
			 double &dxdu, double &dydu, double &dzdu,
			 double &dxdv, double &dydv, double &dzdv,
			 double &dxdw, double &dydw, double &dzdw) const;
      virtual bool invert(double x, double y, double z,
			  double &u, double &v, double &w) const;
      virtual double jacInverse(double x, double y, double z, 
				mTensor2&) const;
      virtual double detJac(double u, double v, double w) const;
      virtual double PushBack (double u, double v, double w, int vsize, std::vector<mVector> &vec) const;
    protected:
    };

#endif
