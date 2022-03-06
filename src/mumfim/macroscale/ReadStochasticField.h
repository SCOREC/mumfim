#ifndef MUMFIM_ANALYSIS_H_
#define MUMFIM_ANALYSIS_H_
#include <vector>

namespace mumfim
{
using vec3d = std::vector<std::vector<std::vector<double>>>;
//typedef std::vector<std::vector<std::vector<double>>> vec3d;
enum class InterpType { Nearest };

class GridData
{
  public:
    GridData(const char * path);
    double interpolate(double px, double py, double pz, InterpType tp = InterpType::Nearest);
    double getMaxVal() {return maxVal;}
    double getMinVal() {return minVal;}
    ~GridData() {}

  private:
    double x0;
    double y0;
    double z0;
    double dx;
    double dy;
    double dz;
    int Nx;
    int Ny;
    int Nz;
    vec3d grid;
    double maxVal;
    double minVal;
};
// gives a string output with the appropriate indices replacing the
// "$alignment and $orientation variables"
std::string getNetworkSuffix(GridData & alignment_grid, GridData & orientation_grid, 
                             const std::string & prefix, double px, double py, double pz,
                             int num_alignment_bins, int num_orientation_bins);
}
#endif
