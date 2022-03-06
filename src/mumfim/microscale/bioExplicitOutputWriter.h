#ifndef MUMFIM_OUTPUT_WRITER_H
#define MUMFIM_OUTPUT_WRITER_H
#include <apf.h>
#include <string>

namespace mumfim
{
// stolen from amsi
// FIXME remove and use amsi version (this is here from ETFEM transition)
struct PvdData
{
  PvdData(std::string filename, double timestep, int part = -1)
      : filename(filename), timestep(timestep), part(part)
  {
    if (part < 0) part = 0;
  }
  std::string filename;
  double timestep;
  int part;
};
class ExplicitOutputWriter
{
  private:
  unsigned long outputFrame;
  std::vector<PvdData> pvdData;
  std::string folder;
  std::string pvdName;
  int header_freq;
  apf::Mesh * mesh;
  std::string fname;
  bool initd;

  public:
  ExplicitOutputWriter(apf::Mesh * mesh, std::string folder,
                       std::string pvdName);
  void writeHistoryData(unsigned long iteration, double W_int, double W_ext, double W_damp, double W_kin, double time);
  void writeFieldData(unsigned long iteration, double time);
  void init();
};
/// Write a paraview collection file for meshes with the format msh_prfx(ii)
/// where ii ranges from 1 to sz. TODO USE AMSI VERSION
void writePvdFile(const std::string & col_fnm,
                  const std::vector<PvdData> & pvd_data);
}
#endif
