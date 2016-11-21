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
 * @brief     Gauss Quadrature functions for hexahedrons.
 *
 */

#include "GaussQuadrature.h"
#include "IntPt.h"
#include "GaussLegendreSimplex.h"
//#include <stdio.h>

int GaussLegendre1D(int, double **, double **);

const double a1  = 0.40824826;
const double ma1 = -0.40824826;
const double a2  = 0.81649658;
const double ma2 = -0.81649658;
const double b1  = 0.70710678;
const double mb1 = -0.70710678;
const double c1  = 0.57735027;
const double mc1 = -0.57735027;
const double w1  = 1.3333333333;
const double xh6[6] = { a1,  a1, ma1, ma1, ma2,  a2};
const double yh6[6] = { b1, mb1,  b1, mb1,  0.,  0.};
const double zh6[6] = {mc1, mc1,  c1,  c1, mc1,  c1};
const double ph6[6] = { w1,  w1,  w1,  w1,  w1,  w1};

IntPt3d GQH1[1] = {
	{ {0.0, 0.0, 0.0}, 8.0 }
};

IntPt3d GQH6[6] = {
	{ {xh6[0], yh6[0], zh6[0]}, ph6[0] },
	{ {xh6[1], yh6[1], zh6[1]}, ph6[1] },
	{ {xh6[2], yh6[2], zh6[2]}, ph6[2] },
	{ {xh6[3], yh6[3], zh6[3]}, ph6[3] },
	{ {xh6[4], yh6[4], zh6[4]}, ph6[4] },
	{ {xh6[5], yh6[5], zh6[5]}, ph6[5] }
};

const double xh8[8] = {0.577350269189626, -0.577350269189626, 0.577350269189626, -0.577350269189626,
                       0.577350269189626, -0.577350269189626, 0.577350269189626, -0.577350269189626
                      };
const double yh8[8] = {0.577350269189626, 0.577350269189626, -0.577350269189626, -0.577350269189626,
                       0.577350269189626, 0.577350269189626, -0.577350269189626, -0.577350269189626
                      };
const double zh8[8] = { -0.577350269189626, -0.577350269189626, -0.577350269189626, -0.577350269189626,
                        0.577350269189626, 0.577350269189626, 0.577350269189626, 0.577350269189626
                      };
const double ph8[8] = {1., 1., 1., 1., 1., 1., 1., 1.};

IntPt3d GQH8[8] = {
	{ {xh8[0], yh8[0], zh8[0]}, ph8[0] },
	{ {xh8[1], yh8[1], zh8[1]}, ph8[1] },
	{ {xh8[2], yh8[2], zh8[2]}, ph8[2] },
	{ {xh8[3], yh8[3], zh8[3]}, ph8[3] },
	{ {xh8[4], yh8[4], zh8[4]}, ph8[4] },
	{ {xh8[5], yh8[5], zh8[5]}, ph8[5] },
	{ {xh8[6], yh8[6], zh8[6]}, ph8[6] },
	{ {xh8[7], yh8[7], zh8[7]}, ph8[7] }
};

IntPt3d *GQH[17] = {GQH1, GQH1, GQH6, GQH8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//int GaussLegendre1D(int, double **, double **);
int GQHnPt[4] = {1, 4, 5, 15};

IntPt3d *getGQHPts(int order)
{
	if (order < 2) {
		return GQH[order];
	}

	if (2 == order) {
		return GQH[3];
	}

	if (3 == order) {
		return GQH[3];
	}

	int n = (order + 3) / 2;
	int index = n - 2 + 4;

	if (!GQH[index]) {
		double *pt, *wt;
		//      printf("o = %d n = %d i = %d\n",order,n*n*n,index);
		GaussLegendre1D(n, &pt, &wt);
		GQH[index] = new IntPt3d[n * n * n];
		int l = 0;

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				for (int k = 0; k < n; ++k) {
					GQH[index][l].pt[0] = pt[i];
					GQH[index][l].pt[1] = pt[j];
					GQH[index][l].pt[2] = pt[k];
					GQH[index][l++].weight = wt[i] * wt[j] * wt[k];
					//      printf ("%f %f %f %f\n",pt[i],pt[j],pt[k],wt[i]*wt[j]*wt[k]);
				}
			}
		}
	}

	return GQH[index];
}

int getNGQHPts(int order)
{
	if (3 == order) {
		return 8;
	}

	if (2 == order) {
		return 8;
	}

	if (order < 2) {
		return GQHnPt[order];
	}

	return ((order + 3) / 2) * ((order + 3) / 2) * ((order + 3) / 2);
}

