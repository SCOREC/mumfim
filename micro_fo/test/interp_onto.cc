#include "bioRVE.h"
#include <amsiAnalysis.h>
#include <algorithm>
int main(int argc, char * argv[])
{
  int result = 0;
  amsi::initAnalysis(argc,argv);
  bio::RVE rve;
  bio::RVE inr(0.25);
  std::vector<double> u {-1.0, -1.0, -1.0,
                        1.0, -1.0, -1.0,
                        1.0, -1.0,  1.0,
                        -1.0, -1.0, 1.0,
                        -1.0, 1.0, -1.0,
                        1.0, 1.0, -1.0,
                        1.0, 1.0, 1.0,
                        -1.0, 1.0, 1.0};
  double sclr = 0.125;
  std::transform(u.begin(),u.end(),u.begin(),std::bind2nd(std::multiplies<double>(),sclr));
  apf::DynamicVector du(24);
  std::copy(u.begin(),u.end(),du.begin());
  bio::displaceRVE(&rve,du);
  //bio::InterpOnto(rve.getElement(),inr.getUField()).apply(inr.getUField());
  apf::Element * inr_elmt = inr.getElement();
  apf::NewArray<apf::Vector3> inr_du;
  apf::getVectorNodes(inr_elmt,inr_du);
  std::cout << inr_du[0]  << std::endl;
  //assert(fabs(inr_du[0].x() - (sclr*0.5) < 1e-8));
  amsi::freeAnalysis();
  return result;
}
