#include "bioRVE2.h"
#include <apfFunctions.h>
#include <apfMeasure.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
namespace bio
{
  void alignFiberNetwork( RVE * rve, FiberNetwork * fn, const double algn_vec[3] )
  {
    /// Populate disp array based on direction of alignment vector.
    double disp[6] = {};
    if (algn_vec[0] > 0)
    {
      double d = -1.0 + std::sqrt(1.0/(1.0 + algn_vec[0]));
      disp[0] = algn_vec[0];
      disp[1] = 0.0;
      disp[2] = d/2.0; disp[3] = -d/2.0; disp[4] = d/2.0; disp[5] = -d/2.0;
    }
    else if(algn_vec[1] > 0)
    {
      double d = -1.0 + std::sqrt(1.0/(1.0 + algn_vec[1]));
      disp[2] = algn_vec[1];
      disp[3] = 0.0;
      disp[0] = d/2.0; disp[1] = -d/2.0; disp[4] = d/2.0; disp[5] = -d/2.0;
    }
    else if(algn_vec[2] > 0)
    {
      double d = -1.0 + std::sqrt(1.0/(1.0 + algn_vec[2]));
      disp[4] = algn_vec[2];
      disp[5] = 0.0;
      disp[0] = d/2.0; disp[1] = -d/2.0; disp[2] = d/2.0; disp[3] = -d/2.0;
    }
    double init_dens = calcFiberDensity(rve,fn);
    affineDeformation(rve, fn, disp); ///< Align fibers via affine deformation.
    // apply disp to rve nodes
    amsi::AccumOp acc;
    amsi::ApplyVector(rve->getNumbering(),rve->getUField(),disp,0,&acc).run();
    //updateRVEBounds(rve, fn, disp);   ///< update RVE boundaries
    /// calculate size of hydrostatic expansion based on density
    double dens = calcFiberDensity(rve,fn);
    double r = -(dens - init_dens)/dens;
    double a = -3.995; double b = -3.995; double c=0.00127;
    double d = (-b - std::sqrt(b * b - 4 * a * (c - r)) )/(2*a);
    disp[0] = d; disp[1] = -d;
    disp[2] = d; disp[3] = -d;
    disp[4] = d; disp[5] = -d;
    affineDeformation(rve, fn, disp);
    updateRVEBounds(rve, fn, disp);
  } /// End AlignFiberNetwork function.
  /// Function to affinely deform RVE box.
  //  - Input variabl disp = [positive x, negative x, positive y, negative y, positive z, negative z]
  void affineDeformation( RVE * rve, FiberNetwork * fn, const double disp[6] )
  {
    (void)fn;
    double rvedisp[24] = {};
    // positive x face of RVE
    rvedisp[1 * 3] = disp[0]; rvedisp[3 * 3] = disp[0]; rvedisp[5 * 3] = disp[0]; rvedisp[7 * 3] = disp[0];
    // negative x face of RVE
    rvedisp[0 * 3] = disp[1]; rvedisp[2 * 3] = disp[1]; rvedisp[4 * 3] = disp[1]; rvedisp[6 * 3] = disp[1];
    // positive y face of RVE
    rvedisp[4 * 3 + 1] = disp[2]; rvedisp[5 * 3 + 1] = disp[2]; rvedisp[6 * 3 + 1] = disp[2]; rvedisp[7 * 3 + 1] = disp[2];
    // negative y face of RVE
    rvedisp[0 * 3 + 1] = disp[3]; rvedisp[1 * 3 + 1] = disp[3]; rvedisp[2 * 3 + 1] = disp[3]; rvedisp[3 * 3 + 1] = disp[3];
    // positive z face of RVE
    rvedisp[2 * 3 + 2] = disp[4]; rvedisp[3 * 3 + 2] = disp[4]; rvedisp[6 * 3 + 2] = disp[4]; rvedisp[7 * 3 + 2] = disp[4];
    // negative z face of RVE
    rvedisp[0 * 3 + 2] = disp[5]; rvedisp[1 * 3 + 2] = disp[5]; rvedisp[4 * 3 + 2] = disp[5]; rvedisp[5 * 3 + 2] = disp[5];
    /// length of RVE box in x, y, z directions.
    double xlen = std::abs( rve->sideCoord(RVE::side::rgt) - rve->sideCoord(RVE::side::lft) );
    double ylen = std::abs( rve->sideCoord(RVE::side::top) - rve->sideCoord(RVE::side::bot) );
    double zlen = std::abs( rve->sideCoord(RVE::side::frt) - rve->sideCoord(RVE::side::bck) );
    (void)xlen;
    (void)ylen;
    (void)zlen;
    /*
    int num_nodes = fn.numNodes();
    for (int ii = 0; ii < num_nodes; ii++)
    {
      Node n = fn.node(ii);
      /// Shift coordinates from -0.5 to 0.5 --> 0 to 1.
      double xx = ( n.x + xlen/2.0 )/xlen;
      double yy = ( n.y + xlen/2.0 )/ylen;
      double zz = ( n.z + xlen/2.0 )/zlen;
      /// Affinely deform nodes according to rvedisp. (trilinterpelation?)
      n.x += (1.0 - zz ) * ( (rvedisp[0 * 3] * (1.0 - xx) + rvedisp[1 * 3] * xx) * (1.0 - yy) +
                             (rvedisp[4 * 3] * (1.0 - xx) + rvedisp[5 * 3] * xx) * yy )
           + zz *          ( (rvedisp[2 * 3] * (1.0 - xx) + rvedisp[3 * 3] * xx) * (1.0 - yy) +
                             (rvedisp[6 * 3] * (1.0 - xx) + rvedisp[7 * 3] * xx) * yy );
      n.y += (1.0 - zz ) * ( (rvedisp[0 * 3 + 1] * (1.0 - xx) + rvedisp[1 * 3 + 1] * xx) * (1.0 - yy) +
                             (rvedisp[4 * 3 + 1] * (1.0 - xx) + rvedisp[5 * 3 + 1] * xx) * yy )
           + zz *          ( (rvedisp[2 * 3 + 1] * (1.0 - xx) + rvedisp[3 * 3 + 1] * xx) * (1.0 - yy) +
                             (rvedisp[6 * 3 + 1] * (1.0 - xx) + rvedisp[7 * 3 + 1] * xx) * yy );
      n.z += (1.0 - zz ) * ( (rvedisp[0 * 3 + 2] * (1.0 - xx) + rvedisp[1 * 3 + 2] * xx) * (1.0 - yy) +
                             (rvedisp[4 * 3 + 2] * (1.0 - xx) + rvedisp[5 * 3 + 2] * xx) * yy )
           + zz *          ( (rvedisp[2 * 3 + 2] * (1.0 - xx) + rvedisp[3 * 3 + 2] * xx) * (1.0 - yy) +
                             (rvedisp[6 * 3 + 2] * (1.0 - xx) + rvedisp[7 * 3 + 2] * xx) * yy );
      fn.setNode(ii,n);
    }
    */
  }
  /*
  void updateRVEBounds(RVE * rve, FiberNetwork * fn, const double disp[6])
  {
    fn.updateSideCoord(FiberNetwork::RIGHT, disp[0]);
    fn.updateSideCoord(FiberNetwork::LEFT, disp[1]);
    fn.updateSideCoord(FiberNetwork::TOP, disp[2]);
    fn.updateSideCoord(FiberNetwork::BOTTOM, disp[3]);
    fn.updateSideCoord(FiberNetwork::FRONT, disp[4]);
    fn.updateSideCoord(FiberNetwork::BACK, disp[5]);
  }
  */
  double calcFiberDensity(RVE * rve, FiberNetwork * fn)
  {
    std::vector<double> lngths;
    calcDimMeasures(fn->getNetworkMesh(),1,std::back_inserter(lngths)); // calc lengths
    double ttl = std::accumulate(lngths.begin(),lngths.end(),0.0);
    double vol = amsi::measureDisplacedMeshEntity(rve->getMeshEnt(),rve->getUField());
    return ttl / vol;
  }
} /// end bio namespace
