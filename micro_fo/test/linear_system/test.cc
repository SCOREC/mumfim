#include "SparskitLinearSystem.h"
#include "FiberNetwork.h"
#include <PCU.h>
#include <mpi.h>
#include <cassert>
int main(int argc, char * argv[])
{
  assert(argv[1]);
  int result = 0;
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  bio::FiberNetwork * fn = bio::loadFromFile(std::string(argv[1]));
  apf::Numbering * num = fn->getNumbering();
  int ndofs = apf::NaiveOrder(num);
  bio::CSR * csr = bio::createCSR(num,ndofs);
  bio::skMat k(csr);
  std::cout << "ndofs : " << csr->getNumEqs() << std::endl
	    << "nnz   : " << csr->getNumNonzero() << std::endl;
  bio::skVec f = bio::makeVec(csr->getNumEqs());
  double vls[7] = {8.0, 6.0, 7.0, 5.0, 3.0, 0.0, 9.0};
  int rws[7] = {0, 2, 4, 6, 7, 8, 9};
  bio::setVecValues(&f,vls,rws,7,false); // set vaules
  result += f[0] == 8.0 ? 0 : 1;
  bio::setVecValues(&f,vls,rws,7,true);  // add values
  result += f[0] == 16.0 ? 0 : 1;
  bio::setVecValues(&f,vls,rws,7,false); // set values
  result += f[0] == 8.0 ? 0 : 1;
  bio::destroyVec(f);
  PCU_Comm_Free();
  MPI_Finalize();
  return result;
}
