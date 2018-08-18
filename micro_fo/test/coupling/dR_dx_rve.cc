#include <PCU.h>
#include <apfDynamicMatrix.h>
#include <lasCSRCore.h>
#include <lasConfig.h>
#include <fstream>
#include <string>
#include "bioFiberNetworkIO.h"
#include "bioMultiscaleCoupling.h"
#include "bioMultiscaleRVEAnalysis.h"
#include "io.h"
std::string fn_dR_dx_rve("dR_dx_rve.test");
std::string fn_fn("fn.test");
inline void notify(std::ostream & out,const std::string & fn)
{
  out << "Attempting to read " << fn << std::endl;
}
int main(int ac, char*av[])
{
  MPI_Init(&ac,&av);
  PCU_Comm_Init();
  bio::FiberNetwork fn(bio::loadFromFile(fn_fn));
  las::Sparsity* csr =
      (las::Sparsity*)las::createCSR(fn.getUNumbering(), fn.getDofCount());
  las::SparskitBuffers* bfrs = new las::SparskitBuffers(fn.getDofCount());
  bio::LinearStructs* vecs =
      bio::createLinearStructs(fn.getDofCount(),1E-6,csr, bfrs);
  // we can leave the solver params uninitialized since we don't solve anything
  // in this test
  bio::micro_fo_solver slvr;
  bio::micro_fo_int_solver slvr_int;
  bio::FiberRVEAnalysis* ans = bio::createFiberRVEAnalysis(&fn, vecs, slvr, slvr_int);
  apf::DynamicMatrix dRdx_rve;
  bio::calcdR_dx_rve(dRdx_rve,ans);
  notify(std::cout,fn_dR_dx_rve);
  std::ifstream fin_dR_dx_rve(fn_dR_dx_rve.c_str());
  apf::DynamicMatrix dRdx_rve_tst;
  fin_dR_dx_rve >> dRdx_rve_tst;
  int pmt[] = {3,4,5, 15,16,17, 0,1,2, 9,10,11, 12,13,14, 21,22,23, 6,7,8, 18,19,20};
  //int pmt[] = {3,4,5, 12,13,14, 0,1,2, 6,7,8, 15,16,17, 21,22,23, 9,10,11, 18,19,20};
  // this is a different permutation than the one in the dRVEdFE test case
  apf::DynamicMatrix dRdx_rve_tst_pmt;
  colPermute(dRdx_rve_tst,&pmt[0],dRdx_rve_tst_pmt);
  bool eq = (dRdx_rve == dRdx_rve_tst_pmt);
  std::cout << "dRdx_rve calculation " << (eq ? "successful!" : "failed!") << std::endl;
  delete vecs;
  delete bfrs;
  las::destroySparsity<las::CSR*>(csr);
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}
