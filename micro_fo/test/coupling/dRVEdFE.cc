#include "bioMultiscaleRVE.h"
#include "bioFiberNetworkIO.h"
#include "io.h"
#include <PCU.h>
#include <fstream>
std::string fn_dx_rve_dx_fe1("dx_rve_dx_fe1.test");
std::string fn_dx_rve_dx_fe2("dx_rve_dx_fe2.test");
std::string fn_fn("fn.test");
bool test1()
{
  bio::RVE rve;
  bio::FiberNetwork fn(bio::loadFromFile(fn_fn));
  // fiber reaction is obsolete now
  bio::micro_fo_header hdr;
  hdr.data[bio::RVE_TYPE] = 0;
  hdr.data[bio::FIELD_ORDER] = 1;
  hdr.data[bio::ELEMENT_TYPE] = 4;
  hdr.data[bio::GAUSS_ID] = 0;
  // many of these are obsolete now
  bio::micro_fo_params prm;
  prm.data[bio::FIBER_RADIUS] = 0.0025;
  prm.data[bio::VOLUME_FRACTION] = 1;
  prm.data[bio::YOUNGS_MODULUS] = 60000;
  prm.data[bio::NONLINEAR_PARAM] = 0;
  prm.data[bio::ORIENTATION_AXIS_X] = 0;
  prm.data[bio::ORIENTATION_AXIS_Y] = 0;
  prm.data[bio::ORIENTATION_AXIS_Z] = 0;
  bio::micro_fo_init_data ini;
  ini.init_data[0] = 0;
  ini.init_data[1] = -1;
  ini.init_data[2] = 1;
  ini.init_data[3] = 1;
  ini.init_data[4] = -1;
  ini.init_data[5] = 0;
  ini.init_data[6] = -1;
  ini.init_data[7] = -1;
  ini.init_data[8] = 0;
  ini.init_data[9] = 0;
  ini.init_data[10] = 0;
  ini.init_data[11] = 0;
  bio::MultiscaleRVE mrve(&rve,&fn,hdr,prm,ini);
  apf::DynamicMatrix dRVEdFE;
  mrve.calcdRVEdFE(dRVEdFE);
  std::ifstream fin_dx_rve_dx_fe(fn_dx_rve_dx_fe1.c_str());
  apf::DynamicMatrix dx_rve_dx_fe;
  fin_dx_rve_dx_fe >> dx_rve_dx_fe;
  int pmt[] = {3,4,5, 12,13,14, 0,1,2, 6,7,8, 15,16,17, 21,22,23, 9,10,11, 18,19,20};
  apf::DynamicMatrix dx_rve_dx_fe_pmt;
  rowPermute(dx_rve_dx_fe,pmt,dx_rve_dx_fe_pmt);
  return (dRVEdFE == dx_rve_dx_fe_pmt);
}
bool test2()
{
  bio::RVE rve;
  bio::FiberNetwork fn(bio::loadFromFile(fn_fn));
  // fiber reaction is obsolete now
  bio::micro_fo_header hdr;
  hdr.data[bio::RVE_TYPE] = 0;
  hdr.data[bio::FIELD_ORDER] = 1;
  hdr.data[bio::ELEMENT_TYPE] = 4;
  hdr.data[bio::GAUSS_ID] = 0;
  // many of these are obsolete now
  bio::micro_fo_params prm;
  prm.data[bio::FIBER_RADIUS] = 0.0025;
  prm.data[bio::VOLUME_FRACTION] = 1;
  prm.data[bio::YOUNGS_MODULUS] = 60000;
  prm.data[bio::NONLINEAR_PARAM] = 0;
  prm.data[bio::ORIENTATION_AXIS_X] = 0;
  prm.data[bio::ORIENTATION_AXIS_Y] = 0;
  prm.data[bio::ORIENTATION_AXIS_Z] = 0;
  bio::micro_fo_init_data ini;
  ini.init_data[0] = -1;
  ini.init_data[1] = 0;
  ini.init_data[2] = 1;
  ini.init_data[3] = -1;
  ini.init_data[4] = -1;
  ini.init_data[5] = 2;
  ini.init_data[6] = -1;
  ini.init_data[7] = -1;
  ini.init_data[8] = 0;
  ini.init_data[9] = 0;
  ini.init_data[10] = -1;
  ini.init_data[11] = 1;
  bio::MultiscaleRVE mrve(&rve,&fn,hdr,prm,ini);
  apf::DynamicMatrix dRVEdFE;
  mrve.calcdRVEdFE(dRVEdFE);
  std::ifstream fin_dx_rve_dx_fe(fn_dx_rve_dx_fe2.c_str());
  apf::DynamicMatrix dx_rve_dx_fe;
  fin_dx_rve_dx_fe >> dx_rve_dx_fe;
  int pmt[] = {3,4,5, 12,13,14, 0,1,2, 6,7,8, 15,16,17, 21,22,23, 9,10,11, 18,19,20};
  apf::DynamicMatrix dx_rve_dx_fe_pmt;
  rowPermute(dx_rve_dx_fe,pmt,dx_rve_dx_fe_pmt);
  return (dRVEdFE == dx_rve_dx_fe_pmt);
}
int main(int ac, char * av[])
{
  MPI_Init(&ac,&av);
  PCU_Comm_Init();
  bool failed = false;
  if(test1())
    std::cout << "dRVEdFE calculation (1) successful!" << std::endl;
  else
  {
    std::cout << "dRVEdFE calculation (1) failed!" << std::endl;
    failed = true;
  }
  if(test2())
    std::cout << "dRVEdFE calculation (2) successful!" << std::endl;
  else
  {
    std::cout << "dRVEdFE calculation (2) failed!" << std::endl;
    failed = true;
  }
  PCU_Comm_Free();
  MPI_Finalize();
  return failed;
}
