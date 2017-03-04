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
            << "  [-c, --case string]                       a string specifying the analysis case to run"
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
struct eps_updt
{
  bio::TissueIteration * itr;
  double & eps;
  double operator()()
  {
    int it = itr->iteration();
    bool inc = it > 20;
    double e = inc ? (pow(10,(it-20)/5))*eps : eps;
    std::cout << "epsilon update (" << it << "): " << e << std::endl;
    return e;
  }
  eps_updt(bio::TissueIteration * i, double & e)
    : itr(i)
    , eps(e)
  {}
};
struct dv_updt
{
  int & dv_its;
  double & eps;
  double operator()()
  {
    bool inc = dv_its > 5;
    double r = inc ? eps + (dv_its-5)*2.5e-4 : eps;
    std::cout << "dv epsilon update (" << dv_its << "): " << r << std::endl;
    ++dv_its;
    if(dv_its > 20)
      r = std::numeric_limits<double>::max();
    return r;
  }
  dv_updt(int & di, double & e)
    : dv_its(di)
    , eps(e)
  {}
};
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
    MPI_Comm_size(MPI_COMM_WORLD,&sz);
    AMSI_DEBUG(Sim_logOn("simmetrix_log"));
    pGModel mdl = GM_load(model_filename.c_str(),NULL,NULL);
    pParMesh msh = PM_load(mesh_filename.c_str(),mdl,NULL);
    for(auto cs = amsi::getNextAnalysisCase(mdl,analysis_case); cs != NULL; )
    {
      amsi::initCase(mdl,cs);
      // initialize logs... should really only happen if they are used....?
      amsi::Log vols;
      amsi::Log dsps;
      if(rnk == 0)
      {
        vols = amsi::activateLog("volume");
        amsi::log(vols) << "LOADSTEP, ENT_TAG, VOL" << std::endl;
        dsps = amsi::activateLog("displacement");
        amsi::log(dsps) << "LOADSTEP, ENT_TAG, X, Y, Z" << std::endl;
      }
      pACase pd = (pACase)AttNode_childByType((pANode)cs,amsi::getSimCaseAttributeDesc(amsi::PROBLEM_DEFINITION));
      pACase ss = (pACase)AttNode_childByType((pANode)cs,amsi::getSimCaseAttributeDesc(amsi::SOLUTION_STRATEGY));
      // discover output policies
      std::vector<pModelItem> dp_itms;
      std::vector<apf::ModelEntity*> dsp_itms;
      amsi::getTrackedModelItems(cs,"output displacement",std::back_inserter(dp_itms));
      std::transform(dp_itms.begin(),dp_itms.end(),std::back_inserter(dsp_itms),amsi::reinterpret_caster<pModelItem,apf::ModelEntity*>());
      std::vector<pModelItem> vl_itms;
      std::vector<apf::ModelEntity*> vol_itms;
      amsi::getTrackedModelItems(cs,"output volume",std::back_inserter(vl_itms));
      std::transform(vl_itms.begin(),vl_itms.end(),std::back_inserter(vol_itms),amsi::reinterpret_caster<pModelItem,apf::ModelEntity*>());
      int nm_stps = AttInfoInt_value((pAttInfoInt)AttNode_childByType((pANode)ss,"num timesteps"));
      bio::NonlinearTissue tssu(mdl,msh,pd,AMSI_COMM_SCALE);
      tssu.recoverSecondaryVariables(0);
      amsi::PetscLAS las;
      bio::logVolumes(vol_itms.begin(),vol_itms.end(),vols,0,tssu.getUField());
      bio::logDisps(dsp_itms.begin(),dsp_itms.end(),dsps,0,tssu.getUField());
      for(int stp = 1; stp <= nm_stps; ++stp)
      {
        bio::TissueIteration itr(&tssu,&las);
        double eps = 1e-8;
        eps_updt eps_scheme(&itr,eps);
        amsi::RelativeResidualConvergence<decltype(eps_scheme)> rs_cnvrg(&las,eps_scheme);
        double dv_eps = 1e-3;
        int dv_its = 0;
        //dv_updt dv_eps_scheme(dv_its,dv_eps);
        //bio::VolumeConvergence<decltype(dv_eps_scheme)> dv_cnvrg(&tssu,dv_eps_scheme); // %dv
        //amsi::MultiConvergence cnvrg(&rs_cnvrg,&dv_cnvrg);
        tssu.setSimulationTime((double)stp/nm_stps);
        if(amsi::numericalSolve(&itr,&rs_cnvrg))
        {
          tssu.recoverSecondaryVariables(stp);
          tssu.step();
          las.step();
          std::stringstream fnm;
          std::string msh_prfx("msh_stp_");
          fnm << amsi::fs->getResultsDir() << "/" << msh_prfx << stp;
          apf::writeVtkFiles(fnm.str().c_str(),tssu.getMesh());
          std::string pvd("/out.pvd");
          amsi::writePvdFile(pvd,msh_prfx,stp);
          bio::logDisps(dsp_itms.begin(),dsp_itms.end(),dsps,stp,tssu.getUField());
          bio::logVolumes(vol_itms.begin(),vol_itms.end(),vols,stp,tssu.getUField());
          if(rnk == 0)
          {
            std::string vl_fl(amsi::fs->getResultsDir() + std::string("/vols.log"));
            std::fstream vl_st(vl_fl.c_str(),std::ios::out | std::ios::app);
            amsi::flushToStream(vols,vl_st);
            std::string ds_fl(amsi::fs->getResultsDir() + std::string("/dsps.log"));
            std::fstream ds_st(ds_fl.c_str(),std::ios::out | std::ios::app);
            amsi::flushToStream(dsps,ds_st);
          }
          std::cout << "Load step " << stp << " completed successfully, continuing..." << std::endl;
        }
        else
        {
          std::cerr << "Load step " << stp << " failed! Exiting..." << std::endl;
          break;
        }
      }
      //retrieve and print logs
      amsi::Log vol = amsi::activateLog(vol_log.c_str());
      amsi::deleteLog(vol);
      std::cout << "Analysis case exited, continuing..." << std::endl;
      if(rnk == 0)
      {
        amsi::deleteLog(dsps);
        amsi::deleteLog(vols);
      }
      amsi::freeCase(cs);
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
