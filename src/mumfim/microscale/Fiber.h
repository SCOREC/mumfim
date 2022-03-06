#ifndef MUMFIM_FIBER_H_
#define MUMFIM_FIBER_H_
#include <amsiEnumOps.h>
namespace mumfim
{
  #define FBR_MEMBERS(OP) OP(truss), OP(euler_bernoulli), OP(timoshenko), OP(num_fiber_members)
  enum FiberMember { FBR_MEMBERS(MAKE_ENUM_OP) };
}
#endif
