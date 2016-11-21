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
#ifndef H_BezierMapping
#define H_BezierMapping
#include "Mapping.h"

  class BezierMapping : public Mapping
  {
  
    public:
      virtual ~BezierMapping(){}
      BezierMapping(pEntity me, int order);
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
      virtual int geomOrder() const;
      virtual int order() const;
    protected:
      std::vector<mPoint> knots;
      int Bezierorder;
      void entGetKnots();
      void getEdgeKnots();
      void getTriFaceKnots() ;
      void getTetKnots();
      double BezierBase(int, double, double, double) const;
      double edgeBezierBase(int, double, double, double) const;
      double triFaceBezierBase(int, double, double, double) const;
      double tetRegionBezierBase(int, double, double, double) const;
      void GradBezierShapeFunction(int, double, double, double,mVector&) const;
      
      void GradBezierShapeFunctionLine(int, double, double, double,mVector&) const;
      void GradBezierShapeFunctionTri(int, double, double, double,mVector&) const;
      void GradBezierShapeFunctionTet(int, double, double, double,mVector&) const;

      void triNormalVector(pEntity,double, double, double, mVector&) const;
      void tetEdgeNormalVector(pEntity,double, double, double, mVector&) const;
      void tetFaceNormalVector(pEntity,double, double, double, mVector&) const;

    };
#endif
