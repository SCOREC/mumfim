#include <amsiAnalysis.h>
#include <amsiDetectOscillation.h>
#include <apf.h>
#include <lasCSRCore.h>
#include <mpi.h>
#include <iostream>
#include "bioFiberNetworkIO.h"
#include "bioFiberRVEAnalysis.h"
#include "bioMicroFOParams.h"
#include "bioMultiscaleRVEAnalysis.h"
#include "bioVerbosity.h"
#include <cmath>
// see python implementation...
bool isclose(double a, double b, double rel_tol=1E-5, double abs_tol=1E-8) {
  return fabs(a-b) <= fmax(rel_tol * fmax(fabs(a), fabs(b)), abs_tol);
}
class TrussIntegratorElemStiff : public bio::TrussIntegrator
{
  public:
  TrussIntegratorElemStiff(bio::FiberNetwork * fn, las::Mat * k, las::Vec * f)
      : bio::TrussIntegrator(fn->getUNumbering(),
                             fn->getUField(),
                             fn->getXpUField(),
                             &(fn->getFiberReactions()[0]),
                             k,
                             f,
                             1)
  {
  }
  virtual void process(apf::MeshElement * e) override
  {
    this->inElement(e);
    int np = countIntPoints(e, this->order);
    for (int p = 0; p < np; ++p)
    {
      ipnode = p;
      apf::Vector3 point;
      getIntPoint(e, this->order, p, point);
      double w = getIntWeight(e, this->order, p);
      double dV = getDV(e, point);
      this->atPoint(point, w, dV);
    }
    std::cout<<"Stiffness Element: "<<id<<"\n";
    for(int i=0;i<6;++i){
      for(int j=0; j<6;++j) {
        std::cout<<es->ke(i, j)<<" ";
        // gaurantee that the elemental stiffness matrix is diagonal
        assert(isclose(es->ke(i,j), es->ke(j,i)));
      }
      std::cout<<"\n";
    }
    std::cout<<"\n\n";
    this->outElement();
  }
};
int main(int argc, char * argv[])
{
  amsi::initAnalysis(argc, argv, MPI_COMM_WORLD);
  std::vector<bio::MicroCase> cases;
  bio::loadMicroFOFromYamlFile(
      "/fasttmp/mersoj/develop/biotissue/micro_fo/test/fiber_only.yaml", cases);
  bio::printMicroFOCase(cases[0]);
  std::string file_name = cases[0].pd.meshFile;
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::stringstream prm_name_ss;
  prm_name_ss << file_name << ".params";
  bio::FiberNetworkReactions rctns;
  apf::Mesh2 * fn_msh = bio::loadFromFile(file_name);
  bio::loadParamsFromFile(fn_msh, prm_name_ss.str(),
                          std::back_inserter(rctns.rctns));
  apf::Field * u = apf::createLagrangeField(fn_msh, "u", apf::VECTOR, 1);
  apf::Numbering * n = apf::createNumbering(u);
  int ndofs = apf::NaiveOrder(n);
  // do we need to zero field? if this assert fails we need to zero the field.
  std::cout << "Problem has " << ndofs << " degrees of freedom" << std::endl;
  assert(ndofs > 0);
  las::Sparsity * sprs = las::createCSR(n, ndofs);
  // clean up the un-needed field and numbering
  apf::destroyField(u);
  apf::destroyNumbering(n);
  las::SparskitBuffers * bfrs = new las::SparskitBuffers(ndofs);
  bio::FiberNetwork * fn = new bio::FiberNetwork(fn_msh);
  fn->setFiberReactions(rctns.rctns);
  bio::LinearStructs * vecs =
      bio::createLinearStructs(ndofs, cases[0].ss.slvrTolerance, sprs, bfrs);
  // get the stiffness matrix
  auto ops = las::getLASOps<las::sparskit>();
  ops->zero(vecs->getK());
  ops->zero(vecs->getU());
  ops->zero(vecs->getF());
  // apf::Mesh * fn = an.getFn()->getNetworkMesh();
  //apf::Integrator * truss_es =
  //    bio::createMicroElementalSystem(fn, vecs->getK(), vecs->getF());
  apf::Integrator * truss_es = new TrussIntegratorElemStiff(fn, vecs->getK(), vecs->getF());
  //truss_es->process(fn_msh);
  apf::MeshEntity * me = NULL;
  apf::MeshIterator * itr = fn_msh->begin(1);
  int ii = 0;
  while ((me = fn_msh->iterate(itr)))
  {
    apf::MeshElement * mlm = apf::createMeshElement(fn_msh, me);
    truss_es->process(mlm);
    apf::destroyMeshElement(mlm);
    ++ii;
  }
  fn_msh->end(itr);
  //las::printSparskitMat(std::cout, vecs->getK());
  delete vecs;
  delete bfrs;
  las::destroySparsity<las::CSR *>(sprs);
  amsi::freeAnalysis();
  return 0;
}
