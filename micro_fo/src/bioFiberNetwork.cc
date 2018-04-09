#include "bioFiberNetwork.h"
#include <iostream>
#include <cassert> // assert
#include <numeric> // accumulate
#include <string>
namespace bio
{
  FiberNetwork::FiberNetwork(apf::Mesh * f)
    : fn(f)
    , u(NULL)
    , xpufnc(NULL)
    , xpu(NULL)
    , du(NULL)
    , dw(NULL)
    , udof(NULL)
    , wdof(NULL)
    , ucnt(0)
    , wcnt(0)
    , tp(FiberMember::truss)
    , dim(f->getDimension())
  {
    assert(f);
    du = apf::createLagrangeField(fn,"du",apf::VECTOR,1);
    u  = apf::createLagrangeField(fn,"u",apf::VECTOR,1);
    xpufnc = new amsi::XpYFunc(fn->getCoordinateField(),u);
    xpu = apf::createUserField(fn,"xpu",apf::VECTOR,apf::getShape(u),xpufnc);
    apf::zeroField(du);
    apf::zeroField(u);
    udof = apf::createNumbering(du);
    //wdof = apf::createNumbering(dw);
    ucnt = apf::NaiveOrder(udof);
    //wcnt = apf::AdjReorder(wdof);
    //apf::SetNumberingOffset(wdof,ucnt);
  }
}
