#include "bioTissueAnalysis.h"
#include "bioTissueMultiscaleAnalysis.h" // only for convergence ops, delete after
#include <amsiAnalysis.h>
#include <amsiUtil.h>
#include <SimError.h>
#include <iostream>
#include <sstream>
#include <cstring>
#include <getopt.h>
#include <string>
#include <fenv.h>
#include <mpi.h>
void display_help_string()
{
  std::cout << "Usage: singlescale [OPTIONS]\n"
            << "  [-h, --help]                              display this help text\n"
            << "  [-g, --model model_file]                  the model file (.smd)\n"
            << "  [-m, --mesh mesh_file]                    the mesh file (.sms)\n";
}
std::string model_filename;
std::string mesh_filename;
bool parse_options(int & argc, char ** & argv)
{
  bool result = true;
  bool quit_loop = false;
  opterr = 0;
  while(!quit_loop)
  {
    static struct option long_options[] =
      {
        {"help",        no_argument,        0, 'h'},
        {"model",       required_argument,  0, 'g'},
        {"mesh",        required_argument,  0, 'm'},
      };
    int option_index = 0;
    int option = getopt_long(argc, argv, "hg:m:", long_options, &option_index);
    switch (option)
    {
    case 'h':
      display_help_string();
      result = false;
      break;
    case 'g':
      model_filename = optarg;
      break;
    case 'm':
      mesh_filename = optarg;
      break;
    case -1:
      //end of options
      quit_loop = true;
      break;
    default:
      break;
    }
  }
  optind = 0; // reset global variable used by getopt_long in case of later use in the program
  opterr = 1;
  return result;
}
int main(int argc, char ** argv)
{
  int result = 0;
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
  if(parse_options(argc,argv))
  {
    amsi::use_simmetrix = true;
    amsi::use_petsc = true;
    amsi::initAnalysis(argc,argv);
    int rnk = -1;
    MPI_Comm_rank(AMSI_COMM_WORLD,&rnk);
    if(rnk > 0)
      amsi::suppressOutput(std::cout);
    int sz = 0;
    MPI_Comm_size(MPI_COMM_WORLD,&sz);
    AMSI_DEBUG(Sim_logOn("simmetrix_log"));
    pGModel mdl = GM_load(model_filename.c_str(),NULL,NULL);
    pParMesh msh = PM_load(mesh_filename.c_str(), sthreadNone, mdl, NULL);
    std::vector<pACase> css;
    pAManager att_mn = SModel_attManager(mdl);
    amsi::getTypeCases(att_mn,"analysis",std::back_inserter(css));
    auto css_nd = css.end();
    for(auto cs = css.begin(); cs != css_nd; ++cs)
    {
      amsi::initCase(mdl,*cs);
      pACase pd = (pACase)AttNode_childByType((pANode)*cs,amsi::getSimCaseAttributeDesc(amsi::PROBLEM_DEFINITION));
      pACase ss = (pACase)AttNode_childByType((pANode)*cs,amsi::getSimCaseAttributeDesc(amsi::SOLUTION_STRATEGY));
      int nm_stps = AttInfoInt_value((pAttInfoInt)AttNode_childByType((pANode)ss,"num timesteps"));
      bio::NonlinearTissue tssu(mdl,msh,pd,AMSI_COMM_SCALE);
      amsi::PetscLAS las;
      bio::TissueIteration itr(&tssu,&las);
      amsi::RelativeResidualConvergence rs_cnvrg(&las,1e-8);
      bio::VolumeConvergenceAccm_Incrmt dv_cnvrg(&tssu,1e-3); // %dv
      amsi::MultiConvergence cnvrg(&rs_cnvrg,&dv_cnvrg);
      for(int stp = 1; stp <= nm_stps; ++stp)
      {
        tssu.setSimulationTime((double)stp/nm_stps);
        if(amsi::numericalSolve(&itr,&cnvrg))
        {
          tssu.step();
          las.step();
          std::stringstream fnm;
          fnm << amsi::fs->getResultsDir() << "/msh_stp_" << stp;
          apf::writeVtkFiles(fnm.str().c_str(),tssu.getMesh());
          std::cout << "Load step " << stp << " completed successfully, continuing..." << std::endl;
        }
        else
        {
          std::cerr << "Load step " << stp << " failed! Exiting..." << std::endl;
          break;
        }
      }
      std::cout << "Analysis case exited, continuing..." << std::endl;
      amsi::freeCase(*cs);
    }
    if(rnk > 0)
      amsi::expressOutput(std::cout);
    M_release(msh);
    //AMAN_release(att_mn);
    GM_release(mdl);
    amsi::freeAnalysis();
  }
  else result++;
  AMSI_DEBUG(Sim_logOff());
  return result;
}
