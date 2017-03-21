#include "bioVolumeConvergence.h"
#include <simNonlinearAnalysis.h>
namespace bio
{
  template <typename O>
    void buildLASConvergenceOperators(pACase ss, amsi::Iteration * it, amsi::LAS * las, O out)
  {
    std::vector<pAttribute> cnvrg_atts;
    amsi::cutPaste<pAttribute>(AttCase_attributes(ss,"convergence operator"),std::back_inserter(cnvrg_atts));
    for(auto cnvrg_att = cnvrg_atts.begin(); cnvrg_att != cnvrg_atts.end(); ++cnvrg_att)
    {
      auto cvg = amsi::buildSimConvergenceOperator(ss,*cnvrg_att,it,las);
      if(cvg != NULL)
        *out++ = cvg;
    }
  }
  template <typename I, typename O>
    void buildVolConvergenceOperators(pACase ss, amsi::Iteration * it, I vl_tks, apf::FIeld * u, O out)
  {
    std::vector<pAttribute> cnvrg_atts;
    amsi::cutPaste<pAttribute>(AttCase_attributes(ss,"convergence operator"),std::back_inserter(cnvrg_atts));
    for(auto cnvrg_att = cnvrg_atts.begin(); cnvrg_att != cnvrg_atts.end(); ++cnvrg_att)
    {
      char * tp = Attribute_imageClass(*cnvrg_att);
      if(std::string("volume convergence").compare(tp) == 0)
        *out++ = buildVolConvergenceOperator(ss,*cnvrg_att,it,u);
      Sim_deleteString(tp);
    }
  }
}
