//#include <amsiMultiscale.h>
//#include <amsiUtil.h>
/*
void planMigration(std::vector<int> & to_serialize,
                   std::vector<int> & send_to,
                   MPI_Comm comm,
                   int local_size,
                   double * local_weights)
{
  int rank = -1;
  MPI_Comm_rank(comm,&rank);
  if(!rank)
  {
    to_serialize.push_back(0);
    send_to.push_back(1);
  }
}
*/
int main(int argc, char ** argv)
{
  /*
  using namespace bio;
  amsi::initMultiscale(argc,argv);
  amsi::Migration migrator(AMSI_COMM_SCALE,
                           amsi::Migration::USER_ALGO,
                           &planMigration);
  int rank = -1;
  MPI_Comm_rank(AMSI_COMM_SCALE,&rank);
  FiberNetwork * network = new FiberNetwork;
  network->readFromFile("fiber_networks/clipped_del_1.txt");
  SparseMatrix * sparse_struct(Make_Structure(network));
  SparskitBuffers buffers(network->numDofs());
  std::vector<MicroFO*> rve(1);
  std::vector< std::vector<char> > rve_data;
  if(rank == 0)
  {
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
    micro_fo_header hdr;
    hdr.data[ELEMENT_TYPE] = 4; // tet
    hdr.data[GAUSS_ID] = 0;
    hdr.data[FIBER_REACTION] = NONLINEAR;
    micro_fo_params prms;
    prms.data[FIBER_RADIUS] = 3.5e-8;
    prms.data[VOLUME_FRACTION] = 0.003;
    prms.data[YOUNGS_MODULUS] = 43200000;
    prms.data[NONLINEAR_PARAM] = 1.2;
    rve[0] = new MicroFO(&hdr.data[0],
                         &pt[0],
                         0,
                         network,
                         sparse_struct,
                         &buffers,
                         &prms.data[0],
                         init_coords,
                         4);
    rve_data.resize(1);
    rve[0]->collectMigrationData();
    rve[0]->getMigrationData(rve_data[0]);
  }
  std::vector<int> to_serialize;
  migrator.plan(to_serialize,1);
  migrator.execute(rve_data);
  if(rve_data.size()) // only rank 1 should have the RVE now...
  {
    rve.reserve(1);
    rve[0] = new MicroFO();
    rve[0]->setMigrationData(rve_data[0]);
    FiberNetwork ** ntwrk = & network;
    SparseMatrix ** sprs = & sparse_struct;
    rve[0]->constructRVEFromMigrationData(&ntwrk,
                                          &sprs,
                                          &buffers);
  }
  amsi::freeMultiscale();
  */
  return 0;
}

