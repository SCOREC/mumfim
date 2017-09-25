#include "lasSparskit.h"
#include <mpi.h>
#include <cassert>
int main(int argc, char * argv[])
{
  MPI_Init(&argc,&argv);
  double eye[] = { 1.0 , 0.0 , 0.0 ,
                   0.0 , 1.0 , 0.0 ,
                   0.0 , 0.0 , 1.0 };
  las::CSR * eye_csr = las::csrFromFull(&eye[0],3,3);
  assert(eye_csr);
  las::Mat * mat_csr = las::createSparskitMatrix(eye_csr);
  assert(mat_csr);
  las::LasOps * ops = las::initSparskitOps();
  assert(ops);
  int rwcls[] = {0, 1, 2};
  double vl = 1.0;
  ops->set(mat_csr,1,&rwcls[0],1,&rwcls[0],&vl);
  ops->set(mat_csr,1,&rwcls[1],1,&rwcls[1],&vl);
  ops->set(mat_csr,1,&rwcls[2],1,&rwcls[2],&vl);
  MPI_Finalize();
  return 0;
}
