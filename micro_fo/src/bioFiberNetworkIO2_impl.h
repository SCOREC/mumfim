#include <fstream>
namespace bio
{
  template <typename O>
    void NetworkLoader::paramsFromStream(std::istream & strm, O out)
  {
    if(!msh)
    {
      std::cerr << "ERROR: cannot load fiber network parameters before the network is loaded!" << std::endl;
      return;
    }
    rct_tg = msh->createIntTag("fiber_reaction",1);
    int nr = -1;
    int ne = -1;
    strm >> nr >> ne;
    for(int ii = 0; ii < nr; ++ii)
      processReactionLine(strm,out);
    for(int ii = 0; ii < ne; ++ii)
      processEdgeReactionLine(strm,ii);
  }
  template <typename O>
    void NetworkLoader::processReactionLine(std::istream & is, O out)
  {
    int tp = -1;
    is >> tp;
    if(tp == FiberConstitutive::linear)
    {
      LinearReaction * rct = new LinearReaction;
      is >> rct->fiber_area >> rct->E;
      *out++ = rct;
    }
    else if(tp == FiberConstitutive::nonlinear)
    {
      NonlinearReaction * rct = new NonlinearReaction;
      is >> rct->fiber_area >> rct->E >> rct->B >> rct->length_ratio_trns;
      *out++ = rct;
    }
  }

}
