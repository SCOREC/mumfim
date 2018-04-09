#include <fstream>
namespace bio
{
  template <typename O>
    void loadParamsFromFile(apf::Mesh2 * msh, const std::string & fnm, O rctns)
  {
    std::fstream strm(fnm);
    ParamLoader ldr(msh);
    ldr.fromStream(strm,rctns);
  }
}
