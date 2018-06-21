//#include "RepresentVolElem.h"
//#include "RVE_Util.h"
//#include "SparseMatrix.h"
//#include <amsiReporter.h>
int main(int argc, char ** argv)
{
  (void)argc;
  (void)argv;
  /*
  using namespace bio;
  int num_rves = atoi(argv[1]);
  amsi::Log l = amsi::activateLog("serialization");
  FiberNetwork * network = new FiberNetwork;
  network->readFromFile("fiber_networks/clipped_del_1.txt");
  FiberNetwork ** fbr_ntwrks = & network;
  SparseMatrix * sparse_struct(Make_Structure(network));
  SparseMatrix ** sprs_strcts = & sparse_struct;
  SparskitBuffers buffers(network->numDofs());
  std::vector<MicroFO*> rve(num_rves);
  std::vector< std::vector<char> > rve_data;
  // determine initial coordinates of the tetrahedron
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
  double pt[3] = {0.25, 0.25, 0.25};
  rve_data.resize(num_rves);
  micro_fo_header hdr;
  hdr.data[ELEMENT_TYPE] = 4; // tet
  hdr.data[GAUSS_ID] = 0;
  hdr.data[FIBER_REACTION] = NONLINEAR;
  micro_fo_params prms;
  prms.data[FIBER_RADIUS] = 3.5e-8;
  prms.data[VOLUME_FRACTION] = 0.003;
  prms.data[YOUNGS_MODULUS] = 43200000;
  prms.data[NONLINEAR_PARAM] = 1.2;
  for(int ii = 0; ii < num_rves; ii++)
  {
    rve[ii] = new MicroFO(&hdr.data[0],
                          &pt[0],
                          0,
                          network,
                          sparse_struct,
                          &buffers,
                          &prms.data[0],
                          init_coords,
                          4);
  }
  double t_0 = amsi::getElapsedTime(l);
  amsi::log(l) << num_rves << " SERIALIZATION_START " << t_0 << std::endl;
  for(int ii = 0; ii < num_rves; ii++)
  {
    rve[ii]->collectMigrationData();
    rve[ii]->getMigrationData(rve_data[ii]);
  }
  double t_1 = amsi::getElapsedTime(l);
  amsi::log(l) << num_rves << " SERIALIZATION_STOP " << t_1 << std::endl;
  for(int ii = 0; ii < num_rves; ii++)
    delete rve[ii];
  rve.reserve(num_rves);
  double t_2 = amsi::getElapsedTime(l);
  amsi::log(l) << num_rves << " DESERIALIZATION_START " << t_2 << std::endl;
  for(int ii = 0; ii < num_rves; ii++)
  {
    rve[ii] = new MicroFO();
    rve[ii]->setMigrationData(rve_data[ii]);
    rve[ii]->constructRVEFromMigrationData(&fbr_ntwrks,
                                           &sprs_strcts,
                                           &buffers);
  }
  double t_3 = amsi::getElapsedTime(l);
  amsi::log(l) << num_rves << " DESERIALIZATION_STOP " << t_3 << std::endl;
  std::string fnm(argv[1]);
  fnm += "_rve_serialization.log";
  std::fstream strm(fnm);
  amsi::flushToStream(l,strm);
*/
  return 0;
}