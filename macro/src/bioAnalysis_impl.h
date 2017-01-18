#include <apfsimWrapper.h>
namespace bio
{
  template <typename I>
    void logVolumes(I bgn_mdl_itm, I nd_mdl_itm, amsi::Log log, int stp, pMesh msh, apf::Field * U)
  {
    for(auto mdl_itm = bgn_mdl_itm; mdl_itm != nd_mdl_itm; ++mdl_itm)
    {
      double v = amsi::measureDisplacedEntity((pGEntity)*mdl_itm,msh,U);
      int rnk = -1;
      MPI_Comm_rank(AMSI_COMM_SCALE,&rnk);
      if(rnk == 0)
        amsi::log(log) << stp << ", "
                       << GEN_tag((pGEntity)*mdl_itm) << ", "
                       << v << std::endl;
    }
  }
}
