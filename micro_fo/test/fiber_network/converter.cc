#include "FiberNetworkIO.h"

#include <ifstream>

int main(int argc, char * argv[])
{
  assert(argv[1] && argv[2]);
  FiberNetwork * fn = bio::parseFromStream_oldFormat(std::ifstream(argv[1]));
  bio::writeToStream(fn,std::ifstream(argv[2]));
  return 0;
}
