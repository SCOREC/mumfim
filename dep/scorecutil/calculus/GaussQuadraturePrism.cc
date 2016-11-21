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
 * @brief     Gauss Quadrature functions for prisms.
 *
 */

#include "GaussQuadrature.h"
#include "IntPt.h"

/* constants for 3D - 2pt rule */

#define Qp21 -0.577350269189626
#define Qp22  0.577350269189626
#define Pp23  0.166666666666667
#define Pp24  0.666666666666667
#define Ww32  0.666666666666667 /* Pw22 * Qw2 */

/* constants for 3D - rule 3 (18 point) */

#define Qp31  -0.774596669241483
#define Qp32   0.000000000000000
#define Qp33   0.774596669241483
#define Qw31   0.555555555555556
#define Qw32   0.888888888888889
#define Qw33   0.555555555555556

#define Pp34   0.816847572980459
#define Pp35   0.091576213509771
#define Pp36   0.108103018168070
#define Pp37   0.445948490915965
#define Pw33   0.219903487310644
#define Pw34   0.446763179356022

#define Ww3331 0.122168604061469  /* Pw33 * Qw31 */
#define Ww3332 0.195469766498350  /* Pw33 * Qw32 */
#define Ww3333 0.122168604061469  /* Pw33 * Qw33 */
#define Ww3431 0.248201766308901  /* Pw34 * Qw31 */
#define Ww3432 0.397122826094242  /* Pw34 * Qw32 */
#define Ww3433 0.248201766308901  /* Pw34 * Qw33 */

/* constants for 3D rule 4 (48 point) */

#define Qp41 -0.861136311594053
#define Qp42 -0.339981043584856
#define Qp43  0.339981043584856
#define Qp44  0.861136311594053
#define Qw41  0.347854845137454
#define Qw42  0.652145154862544
#define Qw43  0.652145154862544
#define Qw44  0.347854845137544

#define Pp41  0.873821971016996
#define Pp42  0.063089014491502
#define Pp43  0.501426509658179
#define Pp44  0.249286745170910
#define Pp45  0.636502499121399
#define Pp46  0.310352451033785
#define Pp47  0.053145049844816
#define Pw41  0.101689812740414      /* S in symtri.c times 2 */
#define Pw42  0.233572551452758      /* T in symtri.c times 2 */
#define Pw43  0.165702151236748      /* U in symtri.c times 2 */

/* the weights should go */
/*(Pw41,Pw41,Pw41,Pw42,Pw42,Pw42,Pw43,Pw43,Pw43,Pw43,Pw43,Pw43)*Qw41 */
/*                   >>                                   *Qw42,Qw43,Qw44 */

#define Ww4141 0.0353732940628734 /* Pw41 * Qw41 */
#define Ww4142 0.0663165186775404 /* Pw41 * Qw42 */
#define Ww4143 0.0663165186775404 /* Pw41 * Qw43 */
#define Ww4144 0.0353732940628734 /* Pw41 * Qw44 */
#define Ww4241 0.0812493437139591 /* Pw42 * Qw41 */
#define Ww4242 0.152323207738798  /* Pw42 * Qw42 */
#define Ww4243 0.152323207738798  /* Pw42 * Qw43 */
#define Ww4244 0.0812493437139591 /* Pw42 * Qw44 */
#define Ww4341 0.057640296157402  /* Pw43 * Qw41 */
#define Ww4342 0.108061855079346  /* Pw43 * Qw42 */
#define Ww4343 0.108061855079346  /* Pw43 * Qw43 */
#define Ww4344 0.057640296157402  /* Pw43 * Qw44 */

IntPt3d GQPrism1[1] = {
	{{ 0.333333333333, 0.333333333333, 0.0}, 4.0}
};


IntPt3d GQPrism2[6] = {
	{{Pp23,Pp23,Qp21},Ww32},
	{{Pp24,Pp23,Qp21},Ww32},
	{{Pp23,Pp24,Qp21},Ww32},
	{{Pp23,Pp23,Qp22},Ww32},
	{{Pp24,Pp23,Qp22},Ww32},
	{{Pp23,Pp24,Qp22},Ww32}
};


IntPt3d *GQPrism[3] = {GQPrism1, GQPrism2, GQPrism2};
int GQPrismnPt[3] = {1, 6, 6};

IntPt3d *getGQPrismPts(int order)
{
	return GQPrism[order];
}

int getNGQPrismPts(int order)
{
	return GQPrismnPt[order];
}

