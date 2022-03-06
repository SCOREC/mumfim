#include "Analysis.h"
#include <cassert>
#include <cstring>
namespace mumfim
{
  static const char * const ConstitutiveTypes[] = {CONSTITUTIVE_TYPES(MAKE_STRING_OP)};
  const char * getConstitutiveTypeString(int ii)
  {
    assert(ii < num_constitutive_types);
    return ConstitutiveTypes[ii];
  }
  int getConstitutiveTypeFromString(const char * tp)
  {
    for(int ii = 0; ii < num_constitutive_types; ++ii)
      if(strcmp(tp,getConstitutiveTypeString(ii)) == 0)
        return ii;
    return -1;
  }
}
