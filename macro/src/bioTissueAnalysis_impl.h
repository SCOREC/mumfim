#include "bioVolumeConvergence.h"
#include <simNonlinearAnalysis.h>
namespace bio
{
  /**
   * Extracts the volume convergence regions from a simmetrix case
   * @param ss simmetrix case
   * @param las the amsi las system used for convergence
   * @param out typically std::back_inserter to some sort of list type.
   */
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
  /**
   * Extracts the volume convergence regions from a simmetrix case
   * @param ss simmetrix case
   * @param it amsi iterator
   * @param u displacement
   * @param vl_tks list of volumes tracked for volume convergence
   * @param out typically std::back_inserter to some sort of list type.
   */
  template <typename I, typename O>
    void buildVolConvergenceOperators(pACase ss, amsi::Iteration * it, apf::Field * u, I vl_tks, O out)
  {
    std::vector<pAttribute> cnvrg_atts;
    amsi::cutPaste<pAttribute>(AttCase_attributes(ss,"convergence operator"),std::back_inserter(cnvrg_atts));
    for(auto cnvrg_att = cnvrg_atts.begin(); cnvrg_att != cnvrg_atts.end(); ++cnvrg_att)
    {
      char * tp = Attribute_imageClass(*cnvrg_att);
      if(std::string("volume convergence").compare(tp) == 0)
      {
        pAttribute rgn_nd_att = Attribute_childByType(*cnvrg_att,"regions");
        pANode rgn_nd = AttributeRefNode_value((pAttributeRefNode)rgn_nd_att);
        *out++ = buildVolConvergenceOperator(ss,*cnvrg_att,it,vl_tks[rgn_nd],u);
      }
      Sim_deleteString(tp);
    }
  }
}
