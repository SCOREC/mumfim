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
int main(int argc, char * argv[])
{
  amsi::initAnalysis(argc, argv, MPI_COMM_WORLD);
  std::vector<bio::MicroCase> cases;
  // bio::loadMicroFOFromYamlFile(
  //    "/fasttmp/mersoj/develop/biotissue/micro_fo/test/fiber_only.yaml",
  //    cases);
  bio::loadMicroFOFromYamlFile(
      "./test_global_stiffness_data/global_stiffness.yaml", cases);
  for (std::size_t i = 0; i < cases.size(); ++i)
  //for (std::size_t i = 0; i < 1; ++i)
  {
    bio::printMicroFOCase(cases[i]);
    std::string file_name = cases[i].pd.meshFile;
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout<<"Reading network"<<std::endl;
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
    std::cout<<"Constructing sparsity pattern"<<std::endl;
    las::Sparsity * sprs = las::createCSR(n, ndofs);
    int nnz = ((las::CSR*)sprs)->getNumNonzero();
    std::cout<<"NNZ: "<<nnz<<std::endl;
    // clean up the un-needed field and numbering
    apf::destroyField(u);
    apf::destroyNumbering(n);
    std::cout<<"Constructing buffers"<<std::endl;
    las::SparskitBuffers * bfrs = new las::SparskitBuffers(ndofs, nnz);
    std::cout<<"createing fiber network class"<<std::endl;
    bio::FiberNetwork * fn = new bio::FiberNetwork(fn_msh);
    fn->setFiberReactions(rctns.rctns);
    std::cout<<"creating sparskit linear structures"<<std::endl;
    bio::LinearStructs<las::MICRO_BACKEND> * vecs =
        bio::createLinearStructs(ndofs, cases[i].ss.slvrTolerance, sprs, bfrs);
    std::cout<<"Zeroing sparkit data"<<std::endl;
    // get the stiffness matrix
    auto ops = las::getLASOps<las::sparskit>();
    ops->zero(vecs->getK());
    ops->zero(vecs->getU());
    ops->zero(vecs->getF());
    // apf::Mesh * fn = an.getFn()->getNetworkMesh();
    std::cout<<"Constructing the elemental stiffness integrator"<<std::endl;
    apf::Integrator * truss_es =
        bio::createMicroElementalSystem(fn, vecs->getK(), vecs->getF());
    std::cout<<"Computing the  global stiffness matrix"<<std::endl;
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
    std::ofstream out("GlobalKMatrix.mtx");
    if (!out.is_open())
    {
      std::cerr << "Could not open file for writing!" << std::endl;
      std::abort();
    }
    // we write the matrix to output to apply output filtering
    std::cout << "Writing matrix" << std::endl;
    las::printSparskitMat(out, vecs->getK(), las::PrintType::mmarket, true);
    out.close();

    std::cout<<"Reading biotissue matrix"<<std::endl;
    std::ifstream in("GlobalKMatrix.mtx");
    if (!in.is_open())
    {
      std::cerr << "Could not open file for reading!" << std::endl;
      std::abort();
    }
    las::Mat * readMat = las::readSparskitMat(in, las::PrintType::mmarket);
    in.close();
    std::size_t found = file_name.find_last_of("/");
    std::size_t found_dot = file_name.find_last_of(".");
    std::stringstream ss;
    ss << "./test_global_stiffness_data/" << file_name.substr(found+1, found_dot-found-1) << "_globalK.mtx";
    std::cout << "Reading Python Matrix: "<<ss.str() << std::endl;
    std::ifstream in2(ss.str());
    if (!in2.is_open())
    {
      std::cerr << "Could not open "<< ss.str() << " for reading!" << std::endl;
      std::abort();
    }
    las::Mat * readMat2 = las::readSparskitMat(in2, las::PrintType::mmarket);
    in2.close();
    std::cout << "Comparing matrix" << std::endl;
    assert(readMat);
    bool close = las::sparskitMatClose(readMat, readMat2, 1E-10, 1E-15);
    std::string comp =
         close ? "True" : "False";
    std::cout << "Matrix was close " << comp << std::endl;
    assert(close);
    las::LasCreateMat* mb = las::getMatBuilder<las::sparskit>(0);
    mb->destroy(readMat);
    mb->destroy(readMat2);
    delete vecs;
    vecs = NULL;
    delete bfrs;
    bfrs = NULL;
    las::destroySparsity<las::CSR *>(sprs);
    delete fn;
    delete truss_es;
  }
  las::LasCreateMat* mb = las::getMatBuilder<las::sparskit>(0);
  delete mb;
  mb = NULL;
  amsi::freeAnalysis();
  return 0;
}
