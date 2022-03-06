#ifndef MUMFIM_ANALYSIS_H_
#define MUMFIM_ANALYSIS_H_
#include <amsiEnumOps.h>
namespace mumfim
{
  #define CONSTITUTIVE_TYPES(OP) OP(isotropic_neohookean), OP(transverse_isotropic), OP(num_constitutive_types)
  enum ConstitutiveType{CONSTITUTIVE_TYPES(MAKE_ENUM_OP)};
  const char * getConstitutiveTypeString(int ii);
  int getConstitutiveTypeFromString(const char * tp);
}
#endif
