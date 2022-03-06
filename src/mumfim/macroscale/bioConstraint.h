#ifndef MUMFIM_CONSTRAINT_H_
#define MUMFIM_CONSTRAINT_H_
#include <apf.h>
#include <amsiLAS.h>
namespace mumfim
{
  class Constraint
  {
  public:
    virtual void update() = 0;
    virtual void apply(amsi::LAS * ls) = 0;
  };
}
#endif
