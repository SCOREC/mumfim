#include "bioTissueAnalysis.h"
#include <amsiInterface.h>
#include <SimFEA.h>
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
    amsi::amsiInit(argc,argv);
    int sz = 0;
    MPI_Comm_size(MPI_COMM_WORLD,&sz);
    AMSI_DEBUG(Sim_logOn("simmetrix_log"));
    pGModel mdl = GM_load(model_filename.c_str(),NULL,NULL);
    pParMesh msh = PM_load(mesh_filename.c_str(), sthreadNone, mdl, NULL);
    std::vector<pACase> css;
    amsi::getTypeCases(SModel_attManager(mdl),"analysis",std::back_inserter(css));
    auto css_nd = css.end();
    for(auto cs = css.begin(); cs != css_nd; ++cs)
    {
      amsi::initCase(mdl,*cs);
      pACase pd = (pACase)AttNode_childByType((pANode)*cs,amsi::getSimCaseAttributeDesc(amsi::PROBLEM_DEFINITION));
      bio::NonLinTissue tssu(mdl,msh,pd,AMSI_COMM_SCALE);
      amsi::PetscLAS las;
      bio::TissueIteration itr(&tssu,&las);
      amsi::RelativeResidualConvergence cnvrg(&las,1e-8);
      int stp = 0;
      int nm_stps = 10;
      do
      {
        stp++;
        tssu.step();
        tssu.setSimulationTime((double)stp/nm_stps);
      } while (amsi::numericalSolve(&itr,&cnvrg));
        std::cout << "Analysis case completed successfully, continuing..." << std::endl;
      amsi::freeCase(*cs);
    }
    amsi::amsiFree();
  }
  else result++;
  AMSI_DEBUG(Sim_logOff());
  return result;
}
