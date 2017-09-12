#include "bioRVE2.h"
#include <cassert>
int main(int argc, char * argv[])
{
  (void)argc;
  (void)argv;
  bio::RVE rve(3);
  assert(rve.getMeshEnt());
  assert(rve.getElement());
  assert(rve.getNumbering());
  assert(rve.getField());
}
