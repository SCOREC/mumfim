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
// Test RVEs by themselves
int main(int argc, char **argv)
{
  using namespace bio;
  // std::ofstream sigma11_file("multi_disp_sigma11.dat",std::ofstream::trunc);
  // std::ofstream sigma22_file("multi_disp_sigma22.dat",std::ofstream::trunc);
  // std::ofstream sigma33_file("multi_disp_sigma33.dat",std::ofstream::trunc);
  int num_networks = atoi(argv[1]);
  std::vector<FiberNetwork*> networks;
  std::vector<SparseMatrix*> sparse_structs;
  networks.reserve(num_networks);
  std::stringstream current_file;
  int max_dofs = -1;
  for(int ii=0;ii<num_networks;ii++)
  {
    current_file << "fiber_networks/del_8200seedL6_1.txt";
//    current_file << "fiber_networks/del_6150seedL5.5_1.txt";
//    current_file << "fiber_networks/del_4450seedL5_1.txt";
    std::cout<<"attempting to read "<<current_file.str()<<std::endl;
//    current_file << "fiber_networks/del_rho100_1.txt";
    bio::FiberNetwork * fn = SPECIFY_FIBER_TYPE ? new SupportFiberNetwork() : new FiberNetwork();
    fn->readFromFile(current_file.str());
    current_file.str(std::string());
    sparse_structs.push_back(Make_Structure(fn));
    networks.push_back(fn);
    int dofs = fn->numDofs();
    max_dofs = dofs > max_dofs ? dofs : max_dofs;
  }
  SparskitBuffers buffers(max_dofs);
  int num_displacements = 1;
  double displacements[] = {0.07};
/*
  int num_displacements = 7;
  double displacements[] = {0.001, 0.005, 0.009, 0.01, 0.05, 0.09, 0.10};
*/
  for(int ii = 0; ii < num_displacements; ii++)
  {
    std::cout<<"displacement="<<displacements[ii]<<std::endl;
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
    apf::Vector3 disp(displacements[ii],0.0,0.0);
    double data[24] = {};
    // displace coordinates 2-3-4
    apf::Vector3 coord2disp = ncoord2 + disp;
    apf::Vector3 coord3disp = ncoord3 + disp;
    apf::Vector3 coord4disp = ncoord4 + disp;
    // set displaced coords
    ncoord1.toArray(data);
    coord2disp.toArray(data+3);
    coord3disp.toArray(data+6);
    coord4disp.toArray(data+9);
    // set displacements
    disp.toArray(data+15);
    disp.toArray(data+18);
    disp.toArray(data+21);
    double results[81];
    // Initialize RVEs
    std::cout<<"Initializing RVEs..."<<std::endl;
    int nn=0;
    micro_fo_header hdr;
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
      (*currentRVE) = new bio::MicroFO(&hdr.data[0],&pt[0],0,networks[nn],sparse_structs[nn],&buffers,&prms.data[0],init_coords,4);
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
//    sigma11_file << displacements[ii] << "\t" << results[0] << std::endl;
//    sigma22_file << displacements[ii] << "\t" << results[1] << std::endl;
//    sigma33_file << displacements[ii] << "\t" << results[2] << std::endl;
//    dsigma11_df_file << displacements[ii] << "\t" << results[9 + 3] << "\t" << results[9 + 6] << "\t" << results[9 + 9] << std::endl;
  }
  return 0;
}
