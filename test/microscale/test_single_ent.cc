#include <apfMeshUtil.h>
#include <vector>
#include <iostream>
#include <PCU.h>
#include <mpi.h>

template <typename O>
void originCenterCube(O cbe_crnrs, double crd)
{
  // APF orders differently than canonical
  // iosparametric hexs
  // front
  *cbe_crnrs++ = apf::Vector3(crd * -1.0,crd * -1.0,crd * -1.0);
  *cbe_crnrs++ = apf::Vector3(crd *  1.0,crd * -1.0,crd * -1.0);
  *cbe_crnrs++ = apf::Vector3(crd *  1.0,crd *  1.0,crd * -1.0);
  *cbe_crnrs++ = apf::Vector3(crd * -1.0,crd *  1.0,crd * -1.0);
  // back
  *cbe_crnrs++ = apf::Vector3(crd * -1.0,crd * -1.0,crd *  1.0);
  *cbe_crnrs++ = apf::Vector3(crd *  1.0,crd * -1.0,crd *  1.0);
  *cbe_crnrs++ = apf::Vector3(crd *  1.0,crd *  1.0,crd *  1.0);
  *cbe_crnrs++ = apf::Vector3(crd * -1.0,crd *  1.0,crd *  1.0);
}

int main(int argc, char ** argv)
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  std::cout<<"Making box"<<std::endl;
  std::vector<apf::Vector3> cbe_crnrs;
  double crd = 0.5;
  originCenterCube(std::back_inserter(cbe_crnrs),crd);
  apf::Mesh::Type tp = apf::Mesh::HEX;
  apf::Mesh2 * msh = amsi::makeSingleEntityMesh(tp,&cbe_crnrs[0]);
  
  std::cout<<"Cleaning up"<<std::endl;
  msh->destroyNative();
  apf::destroyMesh(msh);
  PCU_Comm_Free();
  MPI_Finalize();
}
