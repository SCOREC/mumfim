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
 * @brief     Gauss Quadrature functions for tetrahedrons.
 *
 */

#include "GaussQuadrature.h"
#include "IntPt.h"
#include "GaussLegendreSimplex.h"
//#include "Global.h"

/* constants for 4-point rule */
#define  a4  0.5854101966249685
#define  b4  0.138196601125015

/* constants for 5-point rule */
#define  a5  0.500000000000000
#define  b5  0.166666666666667

/* constants for 16-point rule */
#define  a16  0.0503737941001228
#define  b16  0.0665420686332923
#define  c16  0.7716429020672371
#define  d16  0.0761190326442543
#define  e16  0.1197005277978019
#define  f16  0.0718316452676693
#define  g16  0.4042339134672644

/* constants for 17-point rule */
#define  a17  0.1884185567365411
#define  b17  0.0670385837260428
#define  c17  0.0452855923632739
#define  p17  0.7316369079576180
#define  q17  0.0894543640141273
#define  e17  0.1325810999384657
#define  f17  0.0245400397290300
#define  g17  0.4214394310662522

/* constants for 29-point rule */
#define  a29  0.0904012904601475
#define  b29  0.0191198342789912
#define  c29  0.0436149384066657
#define  d29  0.0258116759619916
#define  p29  0.8277192480479295
#define  q29  0.0574269173173568
#define  e29  0.0513518841255634
#define  f29  0.4860510285706072
#define  g29  0.2312985436519147
#define  h29  0.2967538129690260
#define  i29  0.6081079894015281
#define  j29  0.0475690988147229

IntPt3d GQTet1[1] = {
	{{.25,.25,.25},1.0}
};

IntPt3d GQTet2[4] = {
	{ {a4,b4,b4},.25 },
	{ {b4,a4,b4},.25 },
	{ {b4,b4,a4},.25 },
	{ {b4,b4,b4},.25 }
};

IntPt3d GQTet3[5] = {
	{{0.25,0.25,0.25},-0.8},
	{{a5,b5,b5},0.45},
	{{b5,a5,b5},0.45},
	{{b5,b5,a5},0.45},
	{{b5,b5,b5},0.45}
};

IntPt3d GQTet4[16] = {
	{{c16,d16,d16},a16},
	{{d16,c16,d16},a16},
	{{d16,d16,c16},a16},
	{{d16,d16,d16},a16},
	{{e16,f16,g16},b16},
	{{f16,e16,g16},b16},
	{{e16,g16,g16},b16},
	{{f16,g16,g16},b16},
	{{g16,g16,e16},b16},
	{{g16,g16,f16},b16},
	{{g16,e16,f16},b16},
	{{g16,f16,e16},b16},
	{{e16,g16,f16},b16},
	{{f16,g16,e16},b16},
	{{g16,e16,g16},b16},
	{{g16,f16,g16},b16}
};

IntPt3d GQTet5[17] = {
	{{0.25,0.25,0.25},a17},
	{{p17,q17,q17},b17},
	{{q17,p17,q17},b17},
	{{q17,q17,p17},b17},
	{{q17,q17,q17},b17},
	{{e17,f17,g17},c17},
	{{f17,e17,g17},c17},
	{{e17,g17,g17},c17},
	{{f17,g17,g17},c17},
	{{g17,g17,e17},c17},
	{{g17,g17,f17},c17},
	{{g17,e17,f17},c17},
	{{g17,f17,e17},c17},
	{{e17,g17,f17},c17},
	{{f17,g17,e17},c17},
	{{g17,e17,g17},c17},
	{{g17,f17,g17},c17}
};

IntPt3d GQTet6[29] = {
	{{0.25,0.25,0.25},a29},

	{{p29,q29,q29},b29},
	{{q29,p29,q29},b29},
	{{q29,q29,p29},b29},
	{{q29,q29,q29},b29},

	{{e29,f29,g29},c29},
	{{f29,e29,g29},c29},
	{{e29,g29,g29},c29},
	{{f29,g29,g29},c29},
	{{g29,g29,e29},c29},
	{{g29,g29,f29},c29},
	{{g29,e29,f29},c29},
	{{g29,f29,e29},c29},
	{{e29,g29,f29},c29},
	{{f29,g29,e29},c29},
	{{g29,e29,g29},c29},
	{{g29,f29,g29},c29},

	{{h29,i29,j29},d29},
	{{i29,h29,j29},d29},
	{{h29,j29,j29},d29},
	{{i29,j29,j29},d29},
	{{j29,j29,h29},d29},
	{{j29,j29,i29},d29},
	{{j29,h29,i29},d29},
	{{j29,i29,h29},d29},
	{{h29,j29,i29},d29},
	{{i29,j29,h29},d29},
	{{j29,h29,j29},d29},
	{{j29,i29,j29},d29}
};

IntPt3d * GQTet[7] = {GQTet1,GQTet1,GQTet2,GQTet3,GQTet4,GQTet5,GQTet6};
IntPt3d * GQTetDegen[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int GQTetnPt[7] = {1,1,4,5,16,17,29};

IntPt3d *getGQTetPts(int order)
{
	if(order < 7) {
		return GQTet[order];
	}

	int n = (order + 4) / 2;
	int index = n - 5;
	if (!GQTetDegen[index]) {
		int npts = n * n * n;
		GQTetDegen[index] = new IntPt3d[npts];
		GaussLegendreTet(n,n,n,GQTetDegen[index]);
	}
	return GQTetDegen[index];
}


int getNGQTetPts(int order)
{
	if (order < 7) {
		return GQTetnPt[order];
	} else {
		int n = (order + 4) / 2;
		return n*n*n;
	}
}

