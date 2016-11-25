#include "RepresentVolElem.h"
#include "SparseMatrix.h"
#include "RVE_Util.h"
#include <iostream>
#include <list>
#include <sstream>
#include <cstring>
#include <getopt.h>
#include <string>
#include <fenv.h>
#include <globals.h>
int main(int argc, char **argv)
{
  using namespace bio;
  FiberNetwork * network = NULL;
  if (SPECIFY_FIBER_TYPE)
    network = new SupportFiberNetwork();
  else
    network = new FiberNetwork();
  network->readFromFile("fiber_networks/clipped_del_8/txt");
  SparseMatrix * matrix_struct = (Make_Structure(network));
  SparskitBuffers buffers(network->numDofs());
  int num_displacements = 20;
  double displacements[] = { 0.05, 0.1, 0.15, 0.2, 0.25,
                             0.3, 0.35, 0.4, 0.45, 0.5,
                             0.55, 0.6, 0.65, 0.7, 0.75,
                             0.8, 0.85, 0.9, 0.95, 1};
  for(int ii = 0; ii < 20; ii++)
    displacements[ii] *= -1;
  std::string sigma_names[] = {"sigma11.dat",
                               "sigma12.dat",
                               "sigma13.dat",
                               "sigma22.dat",
                               "sigma23.dat",
                               "sigma33.dat"};
  std::string d_sigma_names[] = {"dsigma11_du1.dat",
                                 "dsigma12_du1.dat",
                                 "dsigma13_du1.dat",
                                 "dsigma22_du1.dat",
                                 "dsigma23_du1.dat",
                                 "dsigma33_du1.dat"};
  for(int ii = 0; ii < 6; ii++)
  {
    std::ofstream clear_file1(sigma_names[ii],std::ofstream::trunc);
    std::ofstream clear_file2(d_sigma_names[ii],std::ofstream::trunc);
  }
  for(int ii = 0; ii < num_displacements; ii++)
  {
    // Setup test data (Not used for size effect test)
    double pt[3] = {0.25, 0.25, 0.25};
    double theta = -2 * atan(sqrt(2) - sqrt(3));
    apf::Matrix3x3 roty( cos(theta), 0,sin(theta),
                         0, 1,         0,
                         -sin(theta), 0, cos(theta));
    apf::Vector3 coord1(-1.0,  0.0, -1.0/sqrt(2));
    apf::Vector3 coord2( 1.0,  0.0, -1.0/sqrt(2));
    apf::Vector3 coord3( 0.0,  1.0,  1.0/sqrt(2));
    apf::Vector3 coord4( 0.0, -1.0,  1.0/sqrt(2));
    apf::Vector3 ncoord1 = roty * coord1;
    apf::Vector3 ncoord2 = roty * coord2;
    apf::Vector3 ncoord3 = roty * coord3;
    apf::Vector3 ncoord4 = roty * coord4;
    double init_coords[12] = {};
    ncoord1.toArray(init_coords);
    ncoord2.toArray(init_coords+3);
    ncoord3.toArray(init_coords+6);
    ncoord4.toArray(init_coords+9);
    apf::Vector3 disp(displacements[ii],0.0,0.0);
    double data[24] = {};
    apf::Vector3 coord_disp = ncoord1 + disp;
    coord_disp.toArray(data);
    disp.toArray(data+12);
    double results[81];
    micro_fo_header hdr;
    hdr.data[ELEMENT_TYPE] = 4; // tet
    hdr.data[GAUSS_ID] = 0;
    hdr.data[FIBER_REACTION] = NONLINEAR;
    micro_fo_params prms;
    prms.data[FIBER_RADIUS] = 3.5e-8;
    prms.data[VOLUME_FRACTION] = 0.003;
    prms.data[YOUNGS_MODULUS] = 43200000;
    prms.data[NONLINEAR_PARAM] = 1.2;
    MicroFO rve(&hdr.data[0],&pt[0],0,network,matrix_struct,&buffers,&prms.data[0],init_coords,4);
    P_computeFiberOnlyRVE(&rve,data,results);
    for(int jj = 0; jj < 6; jj++)
    {
      std::ofstream sigma_file(sigma_names[jj],std::ofstream::app);
      sigma_file << displacements[ii] << "\t" << results[jj] << std::endl;
      sigma_file.close();
      std::ofstream dsigma_file(d_sigma_names[jj], std::ofstream::app);
      dsigma_file << displacements[ii] << "\t" << results[9 + 12*jj] << std::endl;
      dsigma_file.close();
    }
  }
  return 0;
}
