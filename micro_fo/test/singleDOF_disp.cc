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
  FiberNetwork * network = NULL;
  if (SPECIFY_FIBER_TYPE)
    network = new SupportFiberNetwork();
  else
    network = new FiberNetwork();
  network->readFromFile("test_networks/perpendicular_trusses.txt");
  SparseMatrix * matrix_struct = (Make_Structure(network));
  SparskitBuffers buffers(network->numDofs());
/*
  int num_displacements = 20;
  double displacements[] = { 0.05, 0.1, 0.15, 0.2, 0.25,
                             0.3, 0.35, 0.4, 0.45, 0.5,
                             0.55, 0.6, 0.65, 0.7, 0.75,
                             0.8, 0.85, 0.9, 0.95, 1};
*/
  int num_displacements = 24;
  double displacements[] = {-0.03, -0.0275, -0.025, -0.0225, -0.02, -0.0175,
                            -0.015, -0.0125, -0.01, -0.0075, -0.005, -0.0025,
                            0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015,
                            0.0175, 0.02, 0.0225, 0.025, 0.0275, 0.03};
  // Test sent to Ehsan
/*
  int num_displacements = 1;
  double displacements[] = {1};
*/
   /* Specify node number that will be displaced
     Options are nodenum = 1, 2, 3, or 4 for linear tetrahedral element.
     nodenum = 1 is not chosen because node 1 is held fixed.
  */
  int nodenum = 2;
  /* Specify direction of displacement
     Options are disp_dir = x, y, or z.
   */
  std::string disp_dir = "x";
  std::stringstream ss; ss << nodenum; std::string nodenum_string = ss.str();
  std::string sigma_names[] = {"sigma11_node"+nodenum_string+"_"+disp_dir+".dat",
                               "sigma12_node"+nodenum_string+"_"+disp_dir+".dat",
                               "sigma13_node"+nodenum_string+"_"+disp_dir+".dat",
                               "sigma22_node"+nodenum_string+"_"+disp_dir+".dat",
                               "sigma23_node"+nodenum_string+"_"+disp_dir+".dat",
                               "sigma33_node"+nodenum_string+"_"+disp_dir+".dat"};
  std::string d_sigma_names[] = {"dsigma11_node"+nodenum_string+"_"+disp_dir+".dat",
                                 "dsigma12_node"+nodenum_string+"_"+disp_dir+".dat",
                                 "dsigma13_node"+nodenum_string+"_"+disp_dir+".dat",
                                 "dsigma22_node"+nodenum_string+"_"+disp_dir+".dat",
                                 "dsigma23_node"+nodenum_string+"_"+disp_dir+".dat",
                                 "dsigma33_node"+nodenum_string+"_"+disp_dir+".dat"};
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
    int i = nodenum;
    int j;
    if (disp_dir == "x") j=0;
    if (disp_dir == "y") j=1;
    if (disp_dir == "z") j=2;
/* disp_array is        x  y  z
                 node 1[       ]
                 node 2[       ]
                 node 3[       ]
                 node 4[       ]
   Since C++ begins indexing from 0, we use nodenum-1 to
   index disp_array.
   disp_array[0][all] = 0 because node 1 is never moved.
*/
    double disp_array[4][3] = {};
    disp_array[i-1][j] = displacements[ii];
    apf::Vector3 disp2(disp_array[1][0],disp_array[1][1],disp_array[1][2]);
    apf::Vector3 disp3(disp_array[2][0],disp_array[2][1],disp_array[2][2]);
    apf::Vector3 disp4(disp_array[3][0],disp_array[3][1],disp_array[3][2]);
    /*
    apf::Vector3 disp2(displacements[ii],0.0,0.0);
    apf::Vector3 disp3(0.0,0.0,0.0);
    apf::Vector3 disp4(0.0,0.0,0.0);
    */
    double data[24] = {};
    // displace coordinates 2-3-4
    apf::Vector3 coord2disp = ncoord2 + disp2;
    apf::Vector3 coord3disp = ncoord3 + disp3;
    apf::Vector3 coord4disp = ncoord4 + disp4;
    // set displaced coords
    coord2disp.toArray(data+3);
    coord3disp.toArray(data+6);
    coord4disp.toArray(data+9);
    // set displacements (i.e., fedisp)
    disp2.toArray(data+15);
    disp3.toArray(data+18);
    disp4.toArray(data+21);
  std::cout << "displacements :"
            << "(" << 0.0 << ", " << 0.0 << ", " << 0.0 << ")"
            << disp2 << " "
            << disp3 << " "
            << disp4 << std::endl;
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
    // Initialize RVE
    MicroFO rve(&hdr.data[0],&pt[0],0,network,matrix_struct,&buffers,&prms.data[0],init_coords,4);
    // Compute RVE
    std::cout<<"Running RVE"<<std::endl;
    P_computeFiberOnlyRVE(&rve,data,results);
/*
    sigma11_file << displacements[ii] << "\t" << results[0] << std::endl;
    dsigma11_df_file << displacements[ii] << "\t" << results[9 + 3] << "\t" << results[9 + 6] << "\t" << results[9 + 9] << std::endl;
*/
    /* Print out results*/
    for(int jj = 0; jj < 6; jj++)
    {
      std::ofstream sigma_file(sigma_names[jj],std::ofstream::app);
      sigma_file << displacements[ii] << "\t" << results[jj] << std::endl;
      sigma_file.close();
      /* Print out the stress derivatives with respect to x, y, and z of node nodenum */
      std::ofstream dsigma_file(d_sigma_names[jj], std::ofstream::app);
      dsigma_file << displacements[ii]
                  << "\t" << results[9 + 12*jj + (nodenum-1)*3]
                  << "\t" << results[9 + 12*jj + (nodenum-1)*3 + 1]
                  << "\t" << results[9 + 12*jj + (nodenum-1)*3 + 2]
                  << std::endl;
      dsigma_file.close();
    }
  }
  return 0;
}
