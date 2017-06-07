#ifndef BIO_ANALYSIS_H_
#define BIO_ANALYSIS_H_
#include <amsiEnumOps.h>
namespace bio
{
  #define CONSTITUTIVE_TYPES(OP) OP(isotropic_neohookian), OP(transverse_isotropic), OP(num_constitutive_types)
  enum ConstitutiveType{CONSTITUTIVE_TYPES(MAKE_ENUM_OP)};
  const char * getConstitutiveTypeString(int ii);
  int getConstitutiveTypeFromString(const char * tp);
}
#endif
