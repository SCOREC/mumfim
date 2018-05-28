#include "bioMultiscaleRVEAnalysis_impl.h"
#include "io.h"
#include <mpi.h>
#include <fstream>
#include <string>
// in
std::string fn_ds_dx_rve("ds_dx_rve.test");
std::string fn_strs("strs.test");
std::string fn_dV_dx_rve("dV_dx_rve.test");
std::string fn_vol_cnv("vol_cnv.test");
std::string fn_dS_dx_rve("dS_dx_rve.test");
inline void notify(std::ostream & out,const std::string & fn)
{
  out << "Attempting to read " << fn << std::endl;
}
int main(int ac, char * av[])
{
  MPI_Init(&ac,&av);
  // read ds_dx_rve
  notify(std::cout,fn_ds_dx_rve);
  std::ifstream fin_ds_dx_rve(fn_ds_dx_rve.c_str());
  apf::DynamicMatrix ds_dx_rve;
  fin_ds_dx_rve >> ds_dx_rve;
  fin_ds_dx_rve.close();
  // read strs
  notify(std::cout,fn_strs);
  std::ifstream fin_strs(fn_strs.c_str());
  apf::DynamicVector strs;
  fin_strs >> strs;
  fin_strs.close();
  // read dV_dx_rve
  notify(std::cout,fn_dV_dx_rve);
  std::ifstream fin_dV_dx_rve(fn_dV_dx_rve.c_str());
  apf::DynamicVector dV_dx_rve;
  fin_dV_dx_rve >> dV_dx_rve;
  fin_dV_dx_rve.close();
  notify(std::cout,fn_vol_cnv);
  std::ifstream fin_vol_cnv(fn_vol_cnv.c_str());
  double vol = 0.0;
  double cnv = 0.0;
  fin_vol_cnv >> vol >> cnv;
  // convert
  std::cout << "Converting..." << std::endl;
  apf::DynamicMatrix dS_dx_rve;
  bio::convertStressDiv(ds_dx_rve,
                        strs,
                        dV_dx_rve,
                        vol,
                        cnv,
                        dS_dx_rve);
  // read comparison dS_dx_rve
  notify(std::cout,fn_dS_dx_rve);
  std::ifstream fin_dS_dx_rve(fn_dS_dx_rve.c_str());
  apf::DynamicMatrix dS_dx_rve_in;
  fin_dS_dx_rve >> dS_dx_rve_in;
  fin_dS_dx_rve.close();
  bool eq = (dS_dx_rve_in == dS_dx_rve);
  if(eq)
    std::cout << "Stress derivative conversion successful!" << std::endl;
  else
    std::cout << "Stress derivative conversion failed!" << std::endl;
  MPI_Finalize();
  return !eq;
}
