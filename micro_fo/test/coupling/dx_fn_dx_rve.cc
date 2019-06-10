#include <PCU.h>
#include <lasCSRCore.h>
#include <lasConfig.h>
#include <mpi.h>
#include <fstream>
#include <string>
#include "bioFiberNetworkIO.h"
#include "bioMultiscaleCoupling.h"
#include "bioMultiscaleRVEAnalysis.h"
#include "io.h"
std::string fn_fn("fn.test");
std::string fn_dR_dx_rve("dR_dx_rve.test");
std::string fn_k("K_dR_dx_rve.test");
std::string fn_dx_fn_dx_rve("dx_fn_dx_rve.test");
inline void notify(std::ostream & out,const std::string & fn)
{
  out << "Attempting to read " << fn << std::endl;
}
int main(int ac, char * av[])
{
  MPI_Init(&ac,&av);
  PCU_Comm_Init();
  bio::FiberNetwork fn(bio::loadFromFile(fn_fn));
  las::Sparsity* csr =
      (las::Sparsity*)las::createCSR(fn.getUNumbering(), fn.getDofCount());
  las::SparskitBuffers* bfrs = new las::SparskitBuffers(fn.getDofCount());
  bio::LinearStructs<las::MICRO_BACKEND> * vecs =
      bio::createLinearStructs(fn.getDofCount(), 1E-6, csr, bfrs);
  // we can leave the solver params uninitialized since we don't solve anything
  // in this test
  bio::micro_fo_solver slvr;
  bio::micro_fo_int_solver slvr_int;
  slvr_int.data[bio::MICRO_SOLVER_TYPE] = 0;
  bio::FiberRVEAnalysis* ans =
      bio::createFiberRVEAnalysis(&fn, vecs, slvr, slvr_int);
  // load k and overwrite ans->k
  apf::DynamicVector kv;
  notify(std::cout,fn_k);
  std::ifstream fin_k(fn_k.c_str());
  fin_k >> kv;
  fin_k.close();
  double * ka = &(*reinterpret_cast<las::csrMat*>(ans->getK()))(0,0);
  memcpy(ka,&kv[0],sizeof(double)*kv.size());
  apf::DynamicMatrix dR_dx_rve;
  notify(std::cout,fn_dR_dx_rve);
  std::ifstream fin_dR_dx_rve(fn_dR_dx_rve.c_str());
  fin_dR_dx_rve >> dR_dx_rve;
  fin_dR_dx_rve.close();
  apf::DynamicMatrix dx_fn_dx_rve;
  calcdx_fn_dx_rve(dx_fn_dx_rve,ans,dR_dx_rve);
  notify(std::cout,fn_dx_fn_dx_rve);
  std::ifstream fin_dx_fn_dx_rve(fn_dx_fn_dx_rve);
  apf::DynamicMatrix dx_fn_dx_rve_tst;
  fin_dx_fn_dx_rve >> dx_fn_dx_rve_tst;
  bool eq = (dx_fn_dx_rve == dx_fn_dx_rve_tst);
  std::cout << "Calculation of dx_fn_dx_rve " << ( eq ? "succeeded!" : "failed!") << std::endl;
  delete vecs;
  delete bfrs;
  las::destroySparsity<las::CSR*>(csr);
  PCU_Comm_Free();
  MPI_Finalize();
  return !eq;
}
