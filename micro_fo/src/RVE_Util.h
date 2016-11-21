#ifndef RVEUTIL_H_
#define RVEUTIL_H_
#include "RepresentVolElem.h"
#include "SparseMatrix.h"
#ifdef UNIX
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#ifdef _HPUX_SOURCE
#include <sys/unistd.h>
#else
#include <unistd.h>
#endif
#endif
#include <string>
namespace bio
{
  extern int rve_load_balancing;
  extern bool lb_per_iteration;
  int P_computeRVEs(const std::string & fiber_network_filename, int num_fiber_files);
  int P_computeFiberOnlyRVE(MicroFO*,double*,double*);
  void initCoupling(size_t &,size_t &);
  int NumElementNodes(int element_type);
  // Old fiber-only solution stuff
  inline void get_gauss_points(double gp[], double gw[])
  {
    gp[0] = -0.577350269189626;
    gp[1] = 0.577350269189626;
    gw[0] = 1.0;
    gw[1] = 1.0;
  }
  void tsfun(double phi[], double phic[], double phie[], double phis[], double x, double y, double z);
  void tsfun_2d(double phi[], double phic[], double phie[], double x, double y);
  SparseMatrix * Make_Structure(FiberNetwork *);
  void fill_extra_periodic_connections(const std::vector< PBCRelation > & bcs,
                                       FiberNetwork * fiber_network,
                                       SparseMatrix * sparse_struct);
  void store_periodic_locations(const std::vector<PBCRelation> & bcs,
                                FiberNetwork * fiber_network,
                                SparseMatrix * sparse_struct);
  void matrix_multiply(double *amat, int arow, int acol, double *bmat, int brow, int bcol, double *cmat);
  void matrix_multiply_ATranspose(double *amat, int arow, int acol, double *bmat, int brow, int bcol, double *cmat);
  void solve_matrix_system(double *A, double *B, double *soln, int order);
  void calc_ksi_eta_zeta(double xr, double yr, double zr, double u[], double init_coords_loc[], double gpx, double gpy, double gpz);
  void calc_dos(apf::Vector3 rve[8],
                double phic[],
                double phie[],
                double phis[],
                double *dett,
                double ddett[]);
  void calc_deriv(double x[], double y[], double z[], double ux[], double uy[], double uz[], double phic[], double phie[], double phis[], double *dudx, double *dudy, double *dudz, double *dvdx, double *dvdy, double *dvdz, double *dgdx, double *dgdy, double *dgdz, double *dett);
  double calc_norm(double fvec[], int gsize);
} // end of namespace Biotissue
#endif
