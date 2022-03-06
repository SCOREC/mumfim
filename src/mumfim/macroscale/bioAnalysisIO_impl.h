#include <apfMeasure.h>
#include <apfWrapper.h>
#include "bioNonlinearTissue.h"
namespace mumfim
{
  template <typename I>
  void logVolumes(I bgn_mdl_itm,
                  I nd_mdl_itm,
                  amsi::Log log,
                  int stp,
                  apf::Field * U)
  {
    apf::Mesh * msh = apf::getMesh(U);
    int rnk = -1;
    MPI_Comm_rank(AMSI_COMM_SCALE, &rnk);
    for (auto mdl_itm = bgn_mdl_itm; mdl_itm != nd_mdl_itm; ++mdl_itm)
    {
      double v = amsi::measureDisplacedModelEntity(*mdl_itm, U);
      if (rnk == 0)
        amsi::log(log) << stp << ", " << msh->getModelTag(*mdl_itm) << ", " << v
                       << std::endl;
    }
  }
  // only gets vertex nod
  template <typename I>
  void logDisps(I bgn_mdl_itm,
                I nd_mdl_itm,
                amsi::Log log,
                int stp,
                apf::Field * U)
  {
    apf::Mesh * msh = apf::getMesh(U);
    int rnk = -1;
    MPI_Comm_rank(AMSI_COMM_SCALE, &rnk);
    for (auto mdl_ent = bgn_mdl_itm; mdl_ent != nd_mdl_itm; ++mdl_ent)
    {
      std::set<apf::MeshEntity *> classified_verts;
      // ugly, but we don't have the reverse classification in pumi...
      for (int dim = 0; dim <= msh->getModelType(*mdl_ent); ++dim)
      {
        apf::MeshEntity * ent;
        auto * it = msh->begin(dim);
        while ((ent = msh->iterate(it)))
        {
          bool is_classified = (msh->toModel(ent) == *mdl_ent);
          if (is_classified)
          {
            if(dim == 0) {
              classified_verts.insert(ent);
            }
            else {
              apf::Adjacent adjacent_verts;
              msh->getAdjacent(ent,0,adjacent_verts);
              for(const auto& vert: adjacent_verts) {
                classified_verts.insert(vert);
              }
            }
          }
        }
        msh->end(it);
      }
      double dsp[3] = {};
      amsi::getAvgFieldValue(U, classified_verts.begin(),
                             classified_verts.end(), &dsp[0]);
      if (rnk == 0)
      {
        amsi::log(log) << stp << ", " << msh->getModelTag(*mdl_ent) << ", "
                       << dsp[0] << ", " << dsp[1] << ", " << dsp[2]
                       << std::endl;
      }
    }
  }
  template <typename I>
  void logForces(I bgn, I end, amsi::Log lg, int stp, NonlinearTissue * tssu)
  {
    int rnk = -1;
    MPI_Comm_rank(AMSI_COMM_SCALE, &rnk);
    for (auto mdl_ent = bgn; mdl_ent != end; ++mdl_ent)
    {
      double frc[3] = {};
      tssu->getLoadOn(*mdl_ent, &frc[0]);
      if (rnk == 0)
      {
        amsi::log(lg) << stp << ", " << tssu->getMesh()->getModelTag(*mdl_ent)
                      << ", " << frc[0] << ", " << frc[1] << ", " << frc[2]
                      << std::endl;
      }
    }
  }
}  // namespace mumfim
