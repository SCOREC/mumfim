#include <cassert>
namespace bio
{
  template <typename I>
    FiberNetwork::FiberNetwork(apf::Mesh * f,
                               FiberMember t,
                               I rctn_bgn,
                               I rctn_end)
    : fn(f)
    , u(NULL)
    , du(NULL)
    , dw(NULL)
    , udof(NULL)
    , wdof(NULL)
    , ucnt(0)
    , wcnt(0)
    , tp(t)
    , dim(f->getDimension())
    , rctns()
  {
    assert(f);
    du = apf::createLagrangeField(fn,"du",apf::VECTOR,1);
    u  = apf::createLagrangeField(fn,"u",apf::VECTOR,1);
    if(tp == FiberMember::euler_bernoulli || tp == FiberMember::timoshenko)
      apf::createLagrangeField(fn,"dw",apf::VECTOR,1);
    udof = apf::createNumbering(du);
    wdof = apf::createNumbering(dw);
    ucnt = apf::AdjReorder(udof);
    wcnt = apf::AdjReorder(wdof);
    apf::SetNumberingOffset(wdof,ucnt);
    std::copy(rctn_bgn,rctn_end,std::back_inserter(rctns));
  }

}
