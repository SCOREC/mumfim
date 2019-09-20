#include <apfMesh2.h>
#include <cassert>
#include <cmath>
#include <fstream>
#include "bioFiberNetwork.h"
namespace bio
{
  class ParamLoader
  {
    protected:
    apf::Mesh2 * msh;
    apf::MeshTag * id_tg;
    apf::MeshTag * rct_tg;
    std::vector<FiberReaction *> rctns;
    std::map<int, int> fbr_2_rctn;
    template <typename O>
    void processReactionLine(std::istream & strm, O out)
    {
      int tp = -1;
      strm >> tp;
      if (tp == FiberConstitutive::linear)
      {
        double fiber_radius;
        LinearReaction * rct = new LinearReaction;
        strm >> fiber_radius >> rct->E >> rct->fiber_density;
        rct->fiber_area = fiber_radius * fiber_radius * M_PI;
        *out++ = rct;
      }
      else if (tp == FiberConstitutive::nonlinear)
      {
        NonlinearReaction * rct = new NonlinearReaction;
        double fiber_radius;
        strm >> fiber_radius >> rct->E >> rct->B >> rct->length_ratio_trns >>
            rct->fiber_density;
        rct->cmp_ratio_trns = 1.0;
        rct->fiber_area = fiber_radius * fiber_radius * M_PI;
        *out++ = rct;
      }
    }
    int processEdgeReactionLine(std::istream & strm)
    {
      int r = -1;
      strm >> r;
      return r;
    }
    void applyReactionLabels()
    {
      apf::MeshEntity * me = NULL;
      for (apf::MeshIterator * it = msh->begin(1); (me = msh->iterate(it));)
      {
        long id = -1;
        msh->getLongTag(me, id_tg, &id);
        long rctn = fbr_2_rctn[id];
        msh->setLongTag(me, rct_tg, &rctn);
      }
    }
    public:
    ParamLoader(apf::Mesh2 * m)
        : msh(m)
        , id_tg(msh->findTag("id"))
        , rct_tg(msh->createIntTag("fiber_reaction", 1))
        , rctns()
        , fbr_2_rctn()
    {
    }
    template <typename O>
    void fromStream(std::istream & strm, O out)
    {
      int nr = -1;
      int ne = -1;
      strm >> nr >> ne;
      assert(ne == apf::countEntitiesOfType(msh, apf::Mesh::EDGE) &&
             "Must have the same number of edges in the mesh as specified in "
             "the parameter file");
      auto rctns_out = std::back_inserter(rctns);
      // TODO these process functions should verify the number of items
      // in the line as to not cause errors in strange places if the input files
      // are not correct
      for (int ii = 0; ii < nr; ++ii)
        processReactionLine(strm, rctns_out);
      for (int ii = 0; ii < ne; ++ii)
        fbr_2_rctn[ii] = processEdgeReactionLine(strm);
      applyReactionLabels();
      std::copy(rctns.begin(), rctns.end(), out);
    }
  };
  template <typename O>
  void loadParamsFromFile(apf::Mesh2 * msh, const std::string & fnm, O rctns)
  {
    std::ifstream strm(fnm.c_str());
    if (!strm.is_open())
    {
      std::cerr << "Could not open parameter file " << fnm << " for reading.\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    ParamLoader ldr(msh);
    ldr.fromStream(strm, rctns);
  }
  template <typename O>
  void loadParamsFromStream(apf::Mesh2 * msh, std::istream & strm, O rctns)
  {
    ParamLoader ldr(msh);
    ldr.fromStream(strm, rctns);
  }
}  // namespace bio
