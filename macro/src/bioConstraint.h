#ifndef BIO_CONSTRAINT_H_
#define BIO_CONSTRAINT_H_
#include <apf.h>
#include <amsiLAS.h>
namespace bio
{
  class Constraint
  {
  public:
    virtual void update() = 0;
    virtual void apply(amsi::LAS * ls) = 0;
  };
}
#endif
