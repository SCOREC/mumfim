#include "bioFiberReactions.h"
namespace bio
{
  const char * getFiberConstitutiveString(int ii)
  {
    const char * const strings[] = {FBR_RCT_TYPES(MAKE_STRING_OP)};
    assert(tp < FiberConstitutive::num_fbr_constitutive);
    return strings[ii];
  }
}
