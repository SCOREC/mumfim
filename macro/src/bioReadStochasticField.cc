#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <SimUtil.h>
#include <cassert>
#include <limits>
#include <bioReadStochasticField.h>
#include <cmath>


namespace bio {
//using vec3d = std::vector<std::vector<std::vector<double>>>;

static int readRFGrid(const char * RFGFileName, double & x0, double & y0, double & z0,
                double & dx, double & dy, double & dz, int & Nx, int & Ny, int & Nz,
                vec3d & grid, double & maxVal, double & minVal)
{
  std::ifstream RFGFile(RFGFileName);
  if(!RFGFile.is_open())
  {
    std::cerr<<"Could not open "<<RFGFileName<<std::endl;
    std::abort();
    return 1;
  }
  std::string line,tmp;
  getline(RFGFile, tmp, ',');
  x0 = std::atof(tmp.c_str());
  getline(RFGFile, tmp, ',');
  y0 = std::atof(tmp.c_str());
  getline(RFGFile, tmp, ',');
  z0 = std::atof(tmp.c_str());
  getline(RFGFile, tmp, ',');
  dx = std::atof(tmp.c_str());
  getline(RFGFile, tmp, ',');
  dy = std::atof(tmp.c_str());
  getline(RFGFile, tmp, ',');
  dz = std::atof(tmp.c_str());
  getline(RFGFile, tmp, ',');
  Nx = std::atoi(tmp.c_str());
  getline(RFGFile, tmp, ',');
  Ny = std::atoi(tmp.c_str());
  getline(RFGFile, tmp, ',');
  Nz = std::atoi(tmp.c_str());
  grid = vec3d(Nx,std::vector<std::vector<double>>(Ny,std::vector<double>(Nz,0)));
  maxVal = std::numeric_limits<double>::min();
  minVal = std::numeric_limits<double>::max();
  
  for(int i=0; i<Nz; ++i)
  {
    for(int j=0; j<Ny; ++j)
    {
      for(int k=0; k<Nz; ++k)
      {
        getline(RFGFile, tmp, ',');
        grid[k][j][i] = std::atof(tmp.c_str());
        minVal = grid[k][j][i] < minVal ? grid[k][j][i] : minVal;
        maxVal = grid[k][j][i] > maxVal ? grid[k][j][i] : maxVal;
      }
    }
  }
  RFGFile.close();
  return 0;
}

//enum class InterpType { Nearest };
// take a grid and 
static int interpolateGridData(double x0, double y0, double z0, double dx, double dy, double dz,
                        int Nx, int Ny, int Nz,const vec3d & grid,
                        double px, double py, double pz, double & val,
                        InterpType tp)
{
  if ((px < x0 || px > x0+Nx*dx) ||
      (py < y0 || py > y0+Ny*dy) ||
      (pz < z0 || pz > z0+Nz*dz))
  {
    std::cerr << "Point ("<<px<<","<<py<<","<<pz<<") is outside of bounds";
    std::cerr << "("<<x0<<","<<y0<<","<<z0<<")--"<<"("<<(x0+dx*Nx)<<","<<(y0+dy*Ny)<<","<<(z0+dz*Nz)<<")."<<std::endl;
  }
  //assert(px >= 0 && py >= 0 && pz >= 0);
  int ix = std::floor((px-x0)/(dx));
  int iy = std::floor((py-y0)/(dy));
  int iz = std::floor((pz-z0)/(dz));
  // Put a halo around the grid
  if (ix >= Nx) ix = Nx-1;
  else if (ix < 0) ix = 0;
  if (iy >= Ny) iy = Ny-1;
  else if (iy < 0) iy = 0;
  if (iz >= Nz) iz = Nz-1;
  else if (iz < 0) iz = 0;
  assert(ix >= 0 && ix < Nx);
  assert(iy >= 0 && iy < Ny);
  assert(iz >= 0 && iz < Nz);
  val = grid[ix][iy][iz];
  return 0;
}

GridData::GridData(const char* path)
{
  readRFGrid(path,x0,y0,z0,dx,dy,dz,Nx,Ny,Nz,grid,maxVal,minVal);
}
double GridData::interpolate(double px, double py, double pz, InterpType tp)
{
  double interpVal;
  interpolateGridData(x0,y0,z0,dx,dy,dz,Nx,Ny,Nz,grid,px,py,pz,interpVal,tp);
  return interpVal;
}
static bool replace(std::string& str, const std::string& from, const std::string & to)
{
  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
    return false;
  str.replace(start_pos,from.length(),to);
  return true;
}
static int digitize(GridData & grid, double px, double py, double pz, int num_bins)
{
  double interpVal = grid.interpolate(px,py,pz);
  double l = grid.getMaxVal()-grid.getMinVal();
  // we add a small tolerance to the length such that if the max value is selected
  // it will still be less than the number of bins
  double delta = (l+1E-15)/num_bins;
  int digVal=std::floor((interpVal-grid.getMinVal())/delta);
  if(digVal == num_bins)
    digVal = num_bins-1;
  else if ((digVal < 0) || (digVal > num_bins))
  {
    std::cerr<<"An incorrect digitized value was generated on the stochastic grid"<<std::endl;
    std::abort();
  }
  assert(digVal >= 0 && digVal<num_bins);
  return digVal;
}
std::string getNetworkSuffix(GridData & alignment_grid, GridData & orientation_grid, 
                             const std::string & prefix, double px, double py, double pz,
                             int num_alignment_bins, int num_orientation_bins)
{
  std::string output(prefix);

  int alignment_bin = digitize(alignment_grid,px,py,pz,num_alignment_bins);
  int orientation_bin = digitize(orientation_grid,px,py,pz,num_orientation_bins);
  replace(output,"$alignment",std::to_string(alignment_bin));
  replace(output,"$orientation",std::to_string(orientation_bin));
  return output;
}

}
