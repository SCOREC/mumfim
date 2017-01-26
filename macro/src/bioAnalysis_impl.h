#include <apfsimWrapper.h>
#include <apfWrapper.h>
#include <simWrapper.h>
#include <simAnalysis.h>
namespace bio
{
  template <typename I>
    void logVolumes(I bgn_mdl_itm, I nd_mdl_itm, amsi::Log log, int stp, pMesh msh, apf::Field * U)
  {
    int rnk = -1;
    MPI_Comm_rank(AMSI_COMM_SCALE,&rnk);
    for(auto mdl_itm = bgn_mdl_itm; mdl_itm != nd_mdl_itm; ++mdl_itm)
    {
      double v = amsi::measureDisplacedEntity((pGEntity)*mdl_itm,msh,U);
      if(rnk == 0)
        amsi::log(log) << stp << ", "
                       << GEN_tag((pGEntity)*mdl_itm) << ", "
                       << v << std::endl;
    }
  }
  template <typename I>
    void logDisps(I bgn_mdl_itm, I nd_mdl_itm, amsi::Log log, int stp, pMesh msh, apf::Field * U)
  {
    int rnk = -1;
    MPI_Comm_rank(AMSI_COMM_SCALE,&rnk);
    for(auto mdl_itm = bgn_mdl_itm; mdl_itm != nd_mdl_itm; ++mdl_itm)
    {
      std::list<pEntity> ents;
      int dm = amsi::modelItemTypeDim(GEN_type((pGEntity)*mdl_itm));
      amsi::getClassifiedDimEnts(msh,(pGEntity)*mdl_itm,0,dm,std::back_inserter(ents));
      double dsp[3] = {};
      amsi::getAvgFieldValue(U,ents.begin(),ents.end(),&dsp[0]);
      if(rnk == 0)
      {
        amsi::log(log) << stp << ", " << GEN_tag((pGEntity)*mdl_itm)
                       << ", " << dsp[0]
                       << ", " << dsp[1]
                       << ", " << dsp[2] << std::endl;
      }
    }
  }
}
