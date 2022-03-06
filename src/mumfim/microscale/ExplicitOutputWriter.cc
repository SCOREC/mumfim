#include "ExplicitOutputWriter.h"
#include <iomanip>
#include <ostream>
//#include <unistd.h>
#include <sys/stat.h>
#include <fstream>
#include <cassert>
#include <sstream>
#include <iostream>
#include <stdio.h>

namespace mumfim
{
  /*
ExplicitOutputWriter::ExplicitOutputWriter(apf::Mesh * mesh, std::string folder, std::string pvdName,
                     double & W_int, double & W_ext, double & W_damp,
                     double & W_kin, double & time)
    : outputFrame(0)
    , pvdData(std::vector<PvdData>())
    , folder(folder)
    , pvdName(pvdName)
    , header_freq(0)
    , W_int(W_int)
    , W_ext(W_ext)
    , W_damp(W_damp)
    , W_kin(W_kin)
    , time(time)
    , mesh(mesh)
{
  // if the folder doesn't exist create it
  std::string fname = folder + "/";
  struct stat sb;
  if (stat(fname.c_str(), &sb) != 0)
  {
    mkdir(fname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }
  // this clears the micro csv file.
  fname = folder + "/" + "micro_results.csv";
  std::ofstream strm(fname);
}
*/
ExplicitOutputWriter::ExplicitOutputWriter(apf::Mesh * mesh, std::string folder, std::string pvdName)
    : outputFrame(0)
    , pvdData(std::vector<PvdData>())
    , folder(folder)
    , pvdName(pvdName)
    , header_freq(0)
    , mesh(mesh)
    , initd(false)
{
}
void ExplicitOutputWriter::init()
{
  // if the folder doesn't exist create it
  fname = folder + "/";
  struct stat sb;
  if (stat(fname.c_str(), &sb) != 0)
  {
    mkdir(fname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }
  // this clears the micro csv file.
  fname = folder + "/" + "energies.csv";
  std::ofstream strm(fname);
  initd = true;
}
void ExplicitOutputWriter::writeHistoryData(unsigned long iteration, double W_int, double W_ext, double W_damp, double W_kin, double time)
{
  if(!initd) init();
  std::ofstream strm;
  strm.open(fname, std::ios::out | std::ios::app);
  assert(strm.is_open());
  if ((header_freq != 0) &&
      (iteration % header_freq) == 0)
  {
    strm << std::right << std::setw(2) << "# " << std::setw(10) << "Itr,"
         << std::setw(25) << "t," << std::setw(25) << "W_int," << std::setw(25)
         << "W_ext," << std::setw(25) << "W_kin," << std::setw(25) << "W_damp,"
         << std::setw(24) << "W_total"
         << "\n";
  }
  else if(header_freq == 0 && iteration == 0)
  {
    strm << std::right << std::setw(2) << "# " << std::setw(10) << "Itr,"
         << std::setw(25) << "t," << std::setw(25) << "W_int," << std::setw(25)
         << "W_ext," << std::setw(25) << "W_kin," << std::setw(25) << "W_damp,"
         << std::setw(24) << "W_total"
         << "\n";
  }
    double W_total = fabs(W_kin + W_int + W_damp - W_ext);
    strm << std::right << std::setw(11) << iteration << ", " << std::scientific
         << std::setprecision(17) << std::setw(18) << time << ", " << std::setw(18)
         << W_int << ", " << std::setw(18) << W_ext << ", " << std::setw(18)
         << W_kin << ", " << std::setw(18) << fabs(W_damp) << ", " << std::setw(18)
         << W_total << "\n";
  strm.close();
}
void ExplicitOutputWriter::writeFieldData(unsigned long iteration, double time)
{
  if(!initd) init();
  std::stringstream sout;
  sout << "frame_" << outputFrame++;
  apf::writeVtkFiles((folder + "/" + sout.str()).c_str(), mesh, 1);
  pvdData.push_back(PvdData(sout.str(), time, 0));
  writePvdFile(folder + "/" + pvdName, pvdData);
}
void writePvdFile(const std::string & col_fnm,
                  const std::vector<PvdData> & pvd_data)
{
  std::stringstream pvd;
  pvd << col_fnm;
  // std::string pvd(fs->getResultsDir() + col_fnm);
  std::fstream pvdf(pvd.str().c_str(), std::ios::out);
  pvdf << "<VTKFile type=\"Collection\" version=\"0.1\">" << std::endl;
  pvdf << "  <Collection>" << std::endl;
  for (std::size_t ii = 0; ii < pvd_data.size(); ii++)
  {
    pvdf << "    <DataSet timestep=\"" << ii
         << "\" group=\"\" ";
    pvdf << "part=\"" << pvd_data[ii].part << "\" file=\""
         << pvd_data[ii].filename << "/" << pvd_data[ii].filename;
    pvdf << ".pvtu\"/>" << std::endl;
  }
  pvdf << "  </Collection>" << std::endl;
  pvdf << "</VTKFile>" << std::endl;
}
}
