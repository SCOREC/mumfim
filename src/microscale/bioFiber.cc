#include "bioFiber.h"
#include <cassert>
namespace bio
{
  const char * getFiberMemberString(int ii)
  {
    static const char * FiberMemberStrings[] = {FBR_MEMBERS(MAKE_STRING_OP)};
    assert(ii < FiberMember::num_fiber_members && ii >= 0);
    return FiberMemberStrings[ii];
  }
}
