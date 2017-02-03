#include "bioTissueMultiscaleAnalysis.h"
#include "bioMultiscaleTissue.h"
#include "NonLinFibMtx.h"
#include "RVE_Util.h"
#include <amsiMultiscale.h>
#include <amsiAnalysis.h>
#include <amsiUtil.h>
#include <SimError.h>
#include <iostream>
#include <sstream>
#include <cstring>
#include <getopt.h>
#include <string>
#include <fenv.h>
void display_help_string()
{
  std::cout << "Usage: multiscale [OPTIONS]\n"
            << "  [-f, --fibernetwork fiber_network_file]   filename for fiber network (.vect)\n"
            << "  [-h, --help]                              display this help text\n"
            << "  [-s, --model model_file]                  the model file (.smd)\n"
            << "  [-m, --mesh mesh_file]                    the mesh file (.sms)\n"
            << "  [-c, --case string]                       a string specifying the analysis case to run"
            << "  [-b, --balancing]                         specify if load balancing of RVEs is desired";
}
std::string model_filename("");
std::string mesh_filename("");
std::string fiber_network_filename("");
std::string analysis_case("");
bool parse_options(int & argc, char ** & argv)
{
  bool result = true;
  bool quit_loop = false;
  opterr = 0;
  while(!quit_loop)
  {
    static struct option long_options[] =
      {
        {"fibernetwork", required_argument, 0, 'f'},
        {"help",        no_argument,        0, 'h'},
        {"model",       required_argument,  0, 's'},
        {"mesh",        required_argument,  0, 'm'},
        {"balancing",   required_argument,  0, 'b'},
        {"case",        required_argument,  0, 'c'}
      };
    int option_index = 0;
    int option = getopt_long(argc, argv, "hf:l:m:s:b:c:", long_options, &option_index);
    switch (option)
    {
    case 'h':
      display_help_string();
      result = false;
      break;
    case 'f':
      fiber_network_filename = optarg;
      break;
    case 's':
      model_filename = optarg;
      break;
    case 'm':
      mesh_filename = optarg;
      break;
    case 'b':
      bio::rve_load_balancing = atoi(optarg);
      break;
    case 'c':
      analysis_case = optarg;
      break;
    case -1:
// end of options
      quit_loop = true;
      break;
    default:
      break;
    }
  }
  optind = 0; // reset global variable used by getopt_long
  opterr = 1;
  return result;
}
int run_micro_fo(int & argc, char ** & argv, MPI_Comm comm)
{
  int rnk = MPI_Comm_rank(comm,&rnk);
  srand(8675309+rnk);
  bio::P_computeRVEs(fiber_network_filename,30);
  std::cout << "Microscale successfully exited." << std::endl;
  return 0;
}
int run_micro_fm(int & argc, char ** & argv, MPI_Comm comm)
{
  return 0;
}
int run_macro(int & argc, char ** & argv, MPI_Comm cm)
{
  AMSI_DEBUG(Sim_logOn("simmetrix_log"));
  int result = 0;
  amsi::getLocal()->createDD("micro_fo_data"); // coupling data sent every newton iteration
  try
  {
    pGModel mdl = GM_load(model_filename.c_str(),NULL,NULL);
    pParMesh msh = PM_load(mesh_filename.c_str(), sthreadNone, mdl, NULL);
    for(auto cs = amsi::getNextAnalysisCase(mdl,analysis_case); cs != NULL;)
    {
      amsi::initCase(mdl,cs);
      bio::TissueMultiScaleAnalysis TMSA(mdl,msh,cs,cm);
      result += TMSA.run();
      amsi::freeCase(cs);
    }
  } catch (pSimError err) {
    std::cout << "Simmetrix error caught: " << std::endl
              << "  Code  : " << SimError_code(err) << std::endl
              << "  String: " << SimError_toString(err) << std::endl;
    SimError_delete(err);
    return -1;
  }
# ifdef LOGRUN
  amsi::Log macro_stress = amsi::activateLog("macro_stresses");
  int rank = -1;
  MPI_Comm_rank(AMSI_COMM_SCALE,&rank);
  std::stringstream fname;
  fname << amsi::fs->getResultsDir() << "/macro_stress." << rank << ".log";
  {
    std::fstream strss_fs(fname.str().c_str(), std::ios::out | std::ios::app);
    amsi::flushToStream(macro_stress,strss_fs);
  }
# endif
  AMSI_DEBUG(Sim_logOff());
  return result;
}
int main(int argc, char **argv)
{
  int result = 0;
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
  if(parse_options(argc,argv))
  {
    amsi::use_simmetrix = true;
    amsi::use_petsc = true;
    amsi::initMultiscale(argc,argv);
    amsi::initAnalysis(argc,argv);
#   ifdef LOGRUN
    amsi::Log execution_time = amsi::activateLog("execution_time");
#   endif
    amsi::ControlService * control = amsi::ControlService::Instance();
//    control->setSuppressOutput(true);
    control->setScaleMain("macro",&run_macro);
    control->setScaleMain("micro_fo",&run_micro_fo);
    result = control->Execute(argc,argv);
#   ifdef LOGRUN
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    // write the execution time log to an output file
    amsi::log(execution_time) << amsi::getElapsedTime(execution_time) << std::endl;
    std::stringstream filename;
    filename << amsi::fs->getResultsDir() << "execution_time." << rank << ".log";
    {
      std::fstream tm_fs(filename.str().c_str(), std::ios::out | std::ios::app);
      amsi::flushToStream(execution_time, tm_fs);
    }
    amsi::deleteLog(execution_time);
#   endif
    amsi::freeAnalysis();
    amsi::freeMultiscale();
  }
  else
    result++;
  return result;
  // speculative work on new amsi interface under here, mostly just here to see if things compile and link
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
  if(parse_options(argc,argv))
  {
    amsi::Scale * scls[2];
    scls[0] = new amsi::Scale;
    scls[1] = new amsi::Scale;
    double rts[2] = {1.0/128.0, 127.0/128.0};
    amsi::DefaultRankSet rnks;
    amsi::getGlobalRanks(&rnks);
    amsi::assignRatios(&rnks,2,&rts[0],&scls[0]);
    scls[0]->synchronize();
    scls[1]->synchronize();
    if(amsi::assignedTo(scls[0]))
      run_macro(argc,argv,scls[0]->getComm());
    else if(amsi::assignedTo(scls[1]))
      run_micro_fo(argc,argv,scls[1]->getComm());
  }
}
