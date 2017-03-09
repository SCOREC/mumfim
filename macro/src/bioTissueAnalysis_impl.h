#include "bioVolumeConvergence.h"
#include <simNonlinearAnalysis.h>
namespace bio
{
  template <typename O>
    void buildConvergenceOperators(pACase ss, amsi::Iteration * it, amsi::LAS * las, apf::Field * u, O out)
  {
    std::vector<pANode> cnvrg_nds;
    amsi::cutPaste<pANode>(AttNode_childrenByType((pANode)ss,"convergence operator"),std::back_inserter(cnvrg_nds));
    for(auto cnvrg_nd = cnvrg_nds.begin(); cnvrg_nd != cnvrg_nds.end(); ++cnvrg_nd)
    {
      char * tp = AttNode_imageClass(*cnvrg_nd);
      if(std::string("volume convergence").compare(tp) == 0)
        *out++ = buildBioConvergenceOperator(ss,*cnvrg_nd,it,u);
      else
        *out++ = amsi::buildSimConvergenceOperator(ss,*cnvrg_nd,it,las);
      Sim_deleteString(tp);
    }
  }
}
