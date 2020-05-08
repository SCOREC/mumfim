#include <amsiAnalysis.h>
#include <amsiDetectOscillation.h>
#include <apf.h>
#include <lasCSRCore.h>
#include <mpi.h>
#include <cmath>
#include <iostream>
#include "bioFiberNetworkIO.h"
#include "bioFiberRVEAnalysis.h"
#include "bioMicroFOParams.h"
#include "bioMultiscaleRVEAnalysis.h"
#include "bioVerbosity.h"
// see python implementation...
bool isclose(double a, double b, double rel_tol = 1E-5, double abs_tol = 1E-8)
{
  return fabs(a - b) <= fmax(rel_tol * fmax(fabs(a), fabs(b)), abs_tol);
}
class TrussIntegratorElemStiff : public bio::TrussIntegrator
{
  public:
  TrussIntegratorElemStiff(bio::FiberNetwork * fn, las::Mat * k, las::Vec * f)
      : bio::TrussIntegrator(fn->getUNumbering(),
                             fn->getUField(),
                             fn->getXpUField(),
                             fn->getFiberReactions(),
                             k,
                             f,
                             1)
  {
  }
  virtual void outElement()
  {
    std::cout << "Stiffness Element: " << id << "\n";
    for (int i = 0; i < 6; ++i)
    {
      for (int j = 0; j < 6; ++j)
      {
        std::cout << es->ke(i, j) << " ";
        // gaurantee that the elemental stiffness matrix is diagonal
        assert(isclose(es->ke(i, j), es->ke(j, i)));
      }
      std::cout << "\n";
    }
    std::cout << "\n\n";
    bio::TrussIntegrator::outElement();
  }
};
int main(int argc, char * argv[])
{
  amsi::initAnalysis(argc, argv, MPI_COMM_WORLD);
  std::vector<bio::MicroCase> cases;
  bio::loadMicroFOFromYamlFile(argv[1], cases);
  bio::printMicroFOCase(cases[0]);
  std::string file_name = cases[0].pd.meshFile;
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bio::FiberNetworkLibrary network_library;
  network_library.load(file_name,file_name+".params",0,0);
  auto fiber_network = network_library.getOriginalNetwork(0,0);
  // I'm not confident that the move thing here works as intended
  auto an = bio::createFiberRVEAnalysis(std::move(fiber_network), std::move(cases[0].ss));
  apf::Integrator * truss_es =
      new TrussIntegratorElemStiff(an->getFn(), an->getK(), an->getF());
  auto fn_msh = an->getFn()->getNetworkMesh();
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
  amsi::freeAnalysis();
  return 0;
}
