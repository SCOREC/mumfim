#include <apfsimWrapper.h>
#include <apfMeasure.h>
#include <apfWrapper.h>
#include <simClassified.h>
#include <simWrapper.h>
#include <simAnalysis.h>
namespace bio
{
  template <typename I>
    void logVolumes(I bgn_mdl_itm, I nd_mdl_itm, amsi::Log log, int stp, apf::Field * U)
  {
    apf::Mesh * msh = apf::getMesh(U);
    int rnk = -1;
    MPI_Comm_rank(AMSI_COMM_SCALE,&rnk);
    for(auto mdl_itm = bgn_mdl_itm; mdl_itm != nd_mdl_itm; ++mdl_itm)
    {
      double v = amsi::measureDisplacedModelEntity(*mdl_itm,U);
      if(rnk == 0)
        amsi::log(log) << stp << ", "
                       << msh->getModelTag(*mdl_itm) << ", "
                       << v << std::endl;
    }
  }
  // only gets vertex nod
  template <typename I>
    void logDisps(I bgn_mdl_itm, I nd_mdl_itm, amsi::Log log, int stp, apf::Field * U)
  {
    apf::Mesh * msh = apf::getMesh(U);
    int rnk = -1;
    MPI_Comm_rank(AMSI_COMM_SCALE,&rnk);
    for(auto mdl_ent = bgn_mdl_itm; mdl_ent != nd_mdl_itm; ++mdl_ent)
    {
      double dsp[3] = {};
      auto bgn = amsi::beginClassified(msh,*mdl_ent,0);
      auto end = amsi::endClassified(bgn);
      amsi::getAvgFieldValue(U,bgn,end,&dsp[0]);
      if(rnk == 0)
      {
        amsi::log(log) << stp << ", " << msh->getModelTag(*mdl_ent)
                       << ", " << dsp[0]
                       << ", " << dsp[1]
                       << ", " << dsp[2] << std::endl;
      }
    }
  }
}
