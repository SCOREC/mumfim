/******************************************************************************

  (c) 2004-2010 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE-SCOREC file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/

/**
 * @brief     SCORECUtil
 * @file      file name (relative path should be included e.g. /cint/FMDB.h)
 *
 * @brief     Gauss Quadrature functions for pyramids.
 *
 */

#include "GaussQuadrature.h"
#include "IntPt.h"

/* these are the rule 2 int points and weights */

#define mpt455 -0.455341801261479548
#define ppt455  0.455341801261479548
#define mpt577 -0.577350269189625764
#define ppt577  0.577350269189625764
#define mpt122 -0.122008467928146216
#define ppt122  0.122008467928146216
#define ppt622  0.622008467928146212
#define ppt044  0.0446581987385204510

/* these are the rule 3 int points and weights */

#define zero   0.0000000000000000
#define mpt687 -0.6872983346207417
#define ppt687  0.6872983346207417
#define mpt774 -0.7745966692414834
#define ppt774  0.7745966692414834
#define mpt387 -0.3872983346207416
#define ppt387  0.3872983346207416
#define mpt087 -0.08729833462074170
#define ppt087  0.08729833462074170
#define ppt134 0.1349962850858610
#define ppt215 0.2159940561373775
#define ppt345 0.3455904898198040
#define ppt068 0.06858710562414265
#define ppt109 0.1097393689986282
#define ppt175 0.1755829903978052
#define ppt002 0.002177926162424265
#define ppt003 0.003484681859878818
#define ppt005 0.005575490975806120


/* Rule 1*/
IntPt3d GQPyramid1[1] = {
	{{ 0.0, 0.0, 0.0}, 8.0 }
};

/* Rule 2*/
IntPt3d GQPyramid2[8] = {
	{{mpt455,mpt455,mpt577},ppt622},
	{{ppt455,mpt455,mpt577},ppt622},
	{{mpt455,ppt455,mpt577},ppt622},
	{{ppt455,ppt455,mpt577},ppt622},
	{{mpt122,mpt122,ppt577},ppt044},
	{{ppt122,mpt122,ppt577},ppt044},
	{{mpt122,ppt122,ppt577},ppt044},
	{{ppt122,ppt122,ppt577},ppt044}
};

IntPt3d *GQPyramid[3] = {GQPyramid1, GQPyramid2, GQPyramid2};
int GQPyramidnPt[3] = {1, 8, 8};

IntPt3d *getGQPyramidPts(int order)
{
	return GQPyramid[order];
}

int getNGQPyramidPts(int order)
{
	return GQPyramidnPt[order];
}

