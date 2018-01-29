#include "bioAnalysisIO.h"
#include "bioTissueAnalysis.h"
#include "bioTissueMultiscaleAnalysis.h" // only for convergence ops, delete after
#include <amsiAnalysis.h>
#include <amsiCasters.h>
#include <amsiUtil.h>
#include <SimError.h>
#include <iostream>
#include <limits>
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
            << "  [-c, --case string]                       a string specifying the analysis case to run\n"
            << "  [-m, --mesh mesh_file]                    the mesh file (.sms)\n";
}
std::string model_filename("");
std::string mesh_filename("");
std::string analysis_case("");
std::string vol_log("volume");
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
        {"case",        required_argument,  0, 'c'},
        {"mesh",        required_argument,  0, 'm'}
      };
    int option_index = 0;
    int option = getopt_long(argc, argv, "hg:m:c:", long_options, &option_index);
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
    case 'c':
      analysis_case = optarg;
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
    amsi::initAnalysis(argc,argv);
    int rnk = -1;
    MPI_Comm_rank(AMSI_COMM_WORLD,&rnk);
    if(rnk > 0)
      amsi::suppressOutput(std::cout);
    int sz = 0;
    MPI_Comm_size(AMSI_COMM_WORLD,&sz);
    AMSI_DEBUG(Sim_logOn("simmetrix_log"));
    pGModel mdl = GM_load(model_filename.c_str(),NULL,NULL);
    pParMesh msh = PM_load(mesh_filename.c_str(),mdl,NULL);
    for(auto cs = amsi::getNextAnalysisCase(mdl,analysis_case); cs != NULL; cs = amsi::getNextAnalysisCase(mdl,analysis_case))
    {
      amsi::initCase(mdl,cs);
      bio::TissueAnalysis an(mdl,msh,cs,AMSI_COMM_WORLD);
      an.init();
      an.run();
      amsi::freeCase(cs);
    }
    if(rnk > 0)
      amsi::expressOutput(std::cout);
    M_release(msh);
    GM_release(mdl);
    amsi::freeAnalysis();
  }
  return result;
}
