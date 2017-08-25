#include "bioRVE2.h"
#include <cassert>
int main(int agrc, char * argv[])
{
  bio::RVE rve(3);
  assert(rve.getMeshEnt());
  assert(rve.getElement());
  assert(rve.getNumbering());
  assert(rve.getField());
}
