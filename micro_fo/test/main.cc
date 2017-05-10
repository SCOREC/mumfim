#include "bioFiberRVEAnalysis.cc"
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
  int num_fns = atoi(argv[1]);
  std::string fn_prfx(argv[2]);
  std::vector<FiberNetwork*> fns;
  std::vector<las::CSR*> csrs;
  int dof_max = -1;
  NetworkLoader ldr;
  for(int ii = 0; ii < num_fns; ++ii)
  {
    std::stringstream fl;
    fl << fn_prfx.c_str() << ii+1 << ".txt";
    std::ifstream strm(fl.str());
    FiberNetwork * fn = ldr.fromStream(strm);
    int dofs = fn->getDofCount();
    fns.push_back(fn);
    csrs.push_back(las::createCSR(fn->getUNumbering(),dofs));
    dof_max = dofs > dof_max ? dofs : dof_max;
  }
  las::SparskitBuffers bfrs(dof_max);
  // Number of RVEs
  std::vector<FiberRVEAnalysis*> rves;
  for(int ii = 0; ii < num_fns; ++ii)
    rves.push_back(makeAnalysis(fns[ii],&bfrs));
  // Compute RVEs
  std::cout << "Running RVEs..." << std::endl;
  for(auto rve = rves.begin(); rve != rves.end(); ++rve)
  {
    FiberRVEIteration itr(*rve);
    FiberRVEConvergence cnvrg(*rve);
    amsi::numericalSolve(&itr,&cnvrg);
  }
  /*
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
  */
  return 0;
}
