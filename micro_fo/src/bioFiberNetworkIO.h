#ifndef BIO_FIBER_NETWORK_IO_H_
#define BIO_FIBER_NETWORK_IO_H_
#include "bioFiberNetwork.h"
#include "bioFiberReactions.h"
#include <iostream>
#include <apf.h>
#include <apfMesh2.h>
namespace bio
{
  apf::Mesh2 * loadFromStream(std::istream & strm);
  apf::Mesh2 * loadFromFile(const std::string & fnm);
  template <typename O>
    void loadParamsFromFile(apf::Mesh2 * msh, const std::string & fnm, O rctns);
  /**
   * Only here because templates
   */
  class ParamLoader
  {
  protected:
    apf::Mesh2 * msh;
    apf::MeshTag * id_tg;
    apf::MeshTag * rct_tg;
    std::vector<FiberReaction*> rctns;
    std::map<int,int> fbr_2_rctn;
    template <typename O>
      void processReactionLine(std::istream & strm, O out)
    {
      int tp = -1;
      strm >> tp;
      if(tp == FiberConstitutive::linear)
      {
        LinearReaction * rct = new LinearReaction;
        strm >> rct->fiber_area >> rct->E;
        *out++ = rct;
      }
      else if(tp == FiberConstitutive::nonlinear)
      {
        NonlinearReaction * rct = new NonlinearReaction;
        strm >> rct->fiber_area >> rct->E >> rct->B >> rct->length_ratio_trns;
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
      for(apf::MeshIterator * it = msh->begin(1); (me = msh->iterate(it));)
      {
        long id = -1;
        msh->getLongTag(me,id_tg,&id);
        long rctn = fbr_2_rctn[id];
        msh->setLongTag(me,rct_tg,&rctn);
      }
    }
  public:
    ParamLoader(apf::Mesh2 * m)
      : msh(m)
      , id_tg(msh->findTag("id"))
      , rct_tg(msh->createIntTag("fiber_reaction",1))
      , rctns()
      , fbr_2_rctn()
    { }
    template <typename O>
    void fromStream(std::istream & strm, O out)
    {
      int nr = -1;
      int ne = -1;
      strm >> nr >> ne;
      assert(ne == apf::countEntitiesOfType(msh,apf::Mesh::EDGE) && "Must have the same number of edges in the mesh as specified in the parameter file");
      auto rctns_out = std::back_inserter(rctns);
      for(int ii = 0; ii < nr; ++ii)
        processReactionLine(strm,rctns_out);
      for(int ii = 0; ii < ne; ++ii)
        fbr_2_rctn[ii] = processEdgeReactionLine(strm);
      applyReactionLabels();
      std::copy(rctns.begin(),rctns.end(),out);
    }
  };
}
#include "bioFiberNetworkIO_impl.h"
#endif
