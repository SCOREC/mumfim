#ifndef __OUTPUT_WRITER_H__
#define __OUTPUT_WRITER_H__
#include <apf.h>
#include <string>

namespace bio
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
  double & W_int;
  double & W_ext;
  double & W_damp;
  double & W_kin;
  double & time;
  apf::Mesh * mesh;

  public:
  ExplicitOutputWriter(apf::Mesh * mesh, std::string folder,
                       std::string pvdName, double & W_int, double & W_ext,
                       double & W_damp, double & W_kin, double & time);
  void writeHistoryData(unsigned long iteration);
  void writeFieldData(unsigned long iteration);
};
/// Write a paraview collection file for meshes with the format msh_prfx(ii)
/// where ii ranges from 1 to sz.
void writePvdFile(const std::string & col_fnm,
                  const std::vector<PvdData> & pvd_data);
}
#endif
