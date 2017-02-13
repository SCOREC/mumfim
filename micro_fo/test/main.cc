#include "RepresentVolElem.h"
#include "SparseMatrix.h"
#include "Sparskit_Externs.h"
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
  // Clear results file
  std::ofstream out_file;
  out_file.open("results.dat");
  out_file << "cube_length stiffness11" << std::endl;
  out_file.close();
  int num_networks = atoi(argv[1]);
  std::vector<FiberNetwork*> networks;
  std::vector<SparseMatrix*> sparse_structs;
  bio::SparskitBuffers * buffers = NULL;
  networks.reserve(num_networks);
  std::stringstream current_file;
  int max_dofs = -1;
  for(int ii=0;ii<num_networks;ii++)
  {
    current_file << "test_networks/perpendicular_trusses.txt";
    // current_file <<"/fasttmp/chanv3/Develop/biotissue_amsi/macro/test/SimpleCube/fiber_networks/clipped_del_8.txt";
    bio::FiberNetwork * fn = NULL;
    if (SPECIFY_FIBER_TYPE)
      fn = new SupportFiberNetwork();
    else
      fn = new FiberNetwork();
    fn->readFromFile(current_file.str());
    current_file.str(std::string());
    sparse_structs.push_back(Make_Structure(networks.back()));
    networks.push_back(fn);
    int num_dofs = networks.back()->numDofs();
    max_dofs = num_dofs > max_dofs ? num_dofs : max_dofs;
  }
  buffers = new bio::SparskitBuffers(max_dofs);
  // Number of RVEs
  int n = num_networks;
  std::list<bio::MicroFO*> rves;
  rves.resize(n,NULL);
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
  std::cout << "coordinates : "
            << ncoord1 << " "
            << ncoord2 << " "
            << ncoord3 << " "
            << ncoord4 << std::endl;
  double init_coords[12] = {};
  ncoord1.toArray(init_coords);
  ncoord2.toArray(init_coords+3);
  ncoord3.toArray(init_coords+6);
  ncoord4.toArray(init_coords+9);
  // using face formed by 2-3-4
  apf::Vector3 v1 = ncoord2 - ncoord3;
  apf::Vector3 v2 = ncoord2 - ncoord4;
  apf::Vector3 norm = apf::cross(v1,v2).normalize();
  double normal[3] = {};
  norm.toArray(normal);
  std::cout << "normal : " << norm << std::endl;
  double disp_mag = 1.0;
//  apf::Vector3 disp = norm * disp_mag;
  apf::Vector3 disp2 = norm * disp_mag;
  apf::Vector3 disp3 = norm * 0.0;
  apf::Vector3 disp4 = norm * 0.0;
  double data[24] = {};
  // displace coordinates 2-3-4
  apf::Vector3 coord2disp = ncoord2 + disp2;
  apf::Vector3 coord3disp = ncoord3 + disp3;
  apf::Vector3 coord4disp = ncoord4 + disp4;
  // set displaced coords
  coord2disp.toArray(data+3);
  coord3disp.toArray(data+6);
  coord4disp.toArray(data+9);
  std::cout << "displaced coordinates :"
            << ncoord1 << " "
            << coord2disp << " "
            << coord3disp << " "
            << coord4disp << std::endl;
  // set displacements
  disp2.toArray(data+15);
  disp3.toArray(data+18);
  disp4.toArray(data+21);
  std::cout << "displacements :"
            << "(" << 0.0 << ", " << 0.0 << ", " << 0.0 << ")"
            << disp2 << " "
            << disp3 << " "
            << disp4 << std::endl;
  double results[81];
  // Initialize RVEs
  std::cout<<"Initializing RVEs..."<<std::endl;
  int nn=0;
  micro_fo_header hdr;
  hdr.data[RVE_TYPE] = 0;
  hdr.data[ELEMENT_TYPE] = 4; // tet
  hdr.data[GAUSS_ID] = 0;
  hdr.data[FIBER_REACTION] = NONLINEAR;
  micro_fo_params prms;
  prms.data[FIBER_RADIUS] = 3.5e-8;
  prms.data[VOLUME_FRACTION] = 0.003;
  prms.data[YOUNGS_MODULUS] = 43200000;
  prms.data[NONLINEAR_PARAM] = 1.2;
  for(std::list<bio::MicroFO*>::iterator currentRVE = rves.begin(); currentRVE != rves.end(); ++currentRVE)
  {
    std::cout<<"Initializing RVEs ("<< nn <<")..."<<std::endl;
    (*currentRVE) = new bio::MicroFO(&hdr.data[0],
                                           &pt[0],
                                           0,
                                           networks[nn],
                                           sparse_structs[nn],
                                           buffers,
                                           &prms.data[0],
                                           init_coords,
                                           4);
    nn++;
  }
  // Compute RVEs
  std::cout<<"Running RVEs..."<<std::endl;
  nn = 0;
  for(std::list<bio::MicroFO*>::iterator currentRVE = rves.begin(); currentRVE != rves.end(); ++currentRVE)
  {
    std::cout<<"Running RVEs ("<< nn <<")..."<<std::endl;
    P_computeFiberOnlyRVE((*currentRVE),data,results);
    nn++;
  }
  std::cout << "Dimensionalized stress: "
            << results[0] << " " << results[1] << " "
            << results[2] << " " << results[3] << " "
            << results[4] << " " << results[5] << std::endl;
  std::cout << "Unbalanced force term: " << results[6] << " " << results[7] << " " << results[8] << std::endl;
  std::cout << "Dimensionalized stress derivatives: " << std::endl;
  for(int ii = 0; ii < 6; ii++)
  {
    for(int jj = 0; jj < 12; jj++)
      std::cout << results[9 + ii*12 + jj] << " ";
    std::cout << std::endl;
  }
  std::cout << std::endl;
  return 0;
}
