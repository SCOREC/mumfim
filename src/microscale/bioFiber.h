#ifndef BIO_FIBER_H_
#define BIO_FIBER_H_
#include <amsiEnumOps.h>
namespace bio
{
  #define FBR_MEMBERS(OP) OP(truss), OP(euler_bernoulli), OP(timoshenko), OP(num_fiber_members)
  enum FiberMember { FBR_MEMBERS(MAKE_ENUM_OP) };
  const char * getFiberMemberString(int ii);
}
#endif
