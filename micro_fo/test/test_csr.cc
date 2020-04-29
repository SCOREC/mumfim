#include "lasSparskit.h"
#include "bioFiberNetworkIO.h"
#include "bioFiberRVEAnalysis.h"
#include "bioRVE.h"
#include <lasCSRCore.h>
#include <amsiAnalysis.h>
#include <apfMeshIterator.h>
#include <mpi.h>
#include <cassert>
#include <sstream>
int main(int argc, char * argv[])
{
  amsi::initAnalysis(argc,argv,MPI_COMM_WORLD);
  double eye[] = { 1.0 , 0.0 , 0.0 ,
                   0.0 , 1.0 , 0.0 ,
                   0.0 , 0.0 , 1.0 };
  las::Sparsity * eye_csr = las::csrFromFull(&eye[0],3,3);
  assert(eye_csr);
  auto mb = las::getMatBuilder<las::sparskit>(0);
  las::Mat * mat_csr = mb->create(LAS_IGNORE,
                                  LAS_IGNORE,
                                  eye_csr,
                                  MPI_COMM_SELF);
  assert(mat_csr);
  auto ops = las::getLASOps<las::sparskit>();
  assert(ops);
  int rwcls[] = {0, 1, 2};
  double vl = 1.0;
  ops->set(mat_csr,1,&rwcls[0],1,&rwcls[0],&vl);
  ops->set(mat_csr,1,&rwcls[1],1,&rwcls[1],&vl);
  ops->set(mat_csr,1,&rwcls[2],1,&rwcls[2],&vl);
  std::string fbr_ntwrk_str("7 6 0\n0.0  0.0  0.0\n-0.5  0.0  0.0\n0.5  0.0  0.0\n0.0 -0.5  0.0\n0.0  0.5  0.0\n0.0  0.0 -0.5\n0.0  0.0  0.5\n0 1\n0 2\n0 3\n0 4\n0 5\n0 6\n0 7");
  std::stringstream strm(fbr_ntwrk_str);
  apf::Mesh2 * fn_msh = bio::loadFromStream(strm);
  bio::RVE rve;
  bio::FiberNetwork fn(fn_msh);
  apf::Numbering * nm = fn.getUNumbering();
  std::vector<apf::MeshEntity*> bnds;
  auto bgn = amsi::apfMeshIterator(fn_msh,0);
  decltype(bgn) end = amsi::apfEndIterator(fn_msh);
  bio::getBoundaryVerts(&rve,fn_msh,bgn,end,bio::RVE::side::all,std::back_inserter(bnds));
  //bio::applyRVEBC(bnds.begin(),bnds.end(),nm,ops,k,f);
  int ndofs = apf::NaiveOrder(nm);
  //las::Vec * f = las::createSparskitVector(ndofs);
  las::Sparsity * fn_csr = las::createCSR(nm,ndofs);
  las::Mat * k = mb->create(LAS_IGNORE,
                            LAS_IGNORE,
                            fn_csr,
                            MPI_COMM_SELF);
  apf::MeshEntity * vrt = NULL;
  apf::MeshIterator * it = NULL;
  double vls[9] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
  for(it = fn_msh->begin(0); (vrt = fn_msh->iterate(it));)
  {
    apf::NewArray<int> dofs(3);
    apf::getElementNumbers(nm,vrt,dofs);
    ops->set(k,3,&dofs[0],3,&dofs[0],&vls[0]);
  }
  fn_msh->end(it);
  las::printSparskitMat(std::cout,k);
  // more complicated fiber network and csr
  std::string fbr_ntwrk_str2("13 12 0\n0.0 0.0 0.0\n-0.5 0.0 0.0\n-0.25 0.0 0.0\n0.25 0.0 0.0\n0.5 0.0 0.0\n0.0 -0.5 0.0\n0.0 -0.25 0.0\n0.0 0.25 0.0\n0.0 0.5 0.0\n0.0 0.0 -0.5\n0.0 0.0 -0.25\n0.0 0.0 0.25\n0.0 0.0 0.5\n1 2\n2 0\n0 3\n3 4\n5 6\n6 0\n0 7\n7 8\n9 10\n10 0\n0 11\n11 12");
  std::stringstream strm2(fbr_ntwrk_str2);
  apf::Mesh2 * fn_msh2 = bio::loadFromStream(strm2);
  bio::RVE rve2;
  bio::FiberNetwork fn2(fn_msh2);
  apf::Numbering * nm2 = fn2.getUNumbering();
  std::vector<apf::MeshEntity*> bnds2;
  auto bgn2 = amsi::apfMeshIterator(fn_msh2,0);
  decltype(bgn2) end2 = amsi::apfEndIterator(fn_msh2);
  bio::getBoundaryVerts(&rve2,fn_msh2,bgn2,end2,bio::RVE::side::all,std::back_inserter(bnds2));
  //bio::applyRVEBC(bnds2.begin(),bnds2.end(),nm2,ops,k,f);
  int ndofs2 = apf::NaiveOrder(nm2);
  las::Sparsity * fn_csr2 = las::createCSR(nm2,ndofs2);
  las::Mat * k2 = mb->create(LAS_IGNORE,
                             LAS_IGNORE,
                             fn_csr2,
                             MPI_COMM_SELF);
  apf::MeshEntity * edg = NULL;
  apf::MeshIterator * it2 = NULL;
  double vls2[36] { };
  for(int ii = 0; ii < 36; ++ii)
    vls2[ii] = 1.0;
  std::stringstream k2_strm;
  for(it2 = fn_msh2->begin(1); (edg = fn_msh2->iterate(it2));)
  {
    apf::NewArray<int> dofs;
    int nedofs = apf::getElementNumbers(nm2,edg,dofs);
    ops->assemble(k2,nedofs,&dofs[0],nedofs,&dofs[0],&vls2[0]);
  }
  fn_msh2->end(it2);
  // todo : refactor this because it is terrible
  // hacky way to get the matrix out, but it works!
  /*
    las::printSparskitMat(k2_strm,k2);
    double k2_bare[21][21] {};
    std::string ln;
    int ln_cnt = 0;
    while(ln_cnt < 21 && std::getline(k2_strm,ln))
    {
    std::stringstream ln_strm(ln);
    std::copy(std::istream_iterator<double>(ln_strm),
    std::istream_iterator<double>(),
    &k2_bare[ln_cnt][0]);
    ln_cnt++;
    }
    // compare against expected values
    double insidence[49] = { 6, 1, 1, 1, 1, 1, 1,
    1, 2, 0, 0, 0, 0, 0,
    1, 0, 2, 0, 0, 0, 0,
    1, 0, 0, 2, 0, 0, 0,
    1, 0, 0, 0, 2, 0, 0,
    1, 0, 0, 0, 0, 2, 0,
    1, 0, 0, 0, 0, 0, 2 };
    int idx = 0;
    for(int rr = 1; rr < 20; rr+=3)
    for(int cc = 1; cc < 20; cc+=3)
    {
    assert(insidence[idx] == k2_bare[rr][cc]);
    idx++;
    }
  */
  amsi::freeAnalysis();
  return 0;
}
