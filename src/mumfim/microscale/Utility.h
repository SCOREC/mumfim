#ifndef MUMFIM_UTIL_H_
#define MUMFIM_UTIL_H_
#include "apfFieldOp.h" //amsi
#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <type_traits>
#include <Kokkos_Core.hpp> // for KOKKOS_INLINE_FUNCTION
namespace mumfim
{
  /**
   * Calculate the primary measure (lenght, area, volume) for every entity
   * of the specified dimension in the mesh.
   * @param msh The mesh.
   * @param dim The dimensionality of the entities to measure.
   * @param msrs The resulting measures of all entities.
   */
  template <typename O>
  void calcDimMeasures(apf::Mesh * msh, int dim, O msrs);
  /**
   * Make a 3x3 identity matrix.
   */
  apf::Matrix3x3 eye();
  /**
   * Make a 3x3 matrix containing all 1's
   */
  apf::Matrix3x3 ones();
  /**
   * Apply an incremental deformaiton gradient to the given fields
   */
  class ApplyIncrementalDeformationGradient : public amsi::FieldOp
  {
    protected:
    apf::Field * xyz;
    apf::Field * du;
    apf::Field * u;
    apf::MeshEntity * ent;
    apf::Matrix3x3 FmI;
    std::vector<apf::MeshEntity *>::iterator ent_bgn;
    std::vector<apf::MeshEntity *>::iterator ent_end;
    apf::Vector3 nd_xyz;
    bool hasEntList; // bool is true if the iterators to a list are passed in
    void init(apf::Mesh * msh, const apf::Matrix3x3 & F)
    {
      int d = msh->getDimension();
      for (int ii = 0; ii < d; ++ii)
        for (int jj = 0; jj < d; ++jj)
          FmI[ii][jj] = F[ii][jj] - (ii == jj ? 1.0 : 0.0);
    }

    public:
    ApplyIncrementalDeformationGradient(apf::Matrix3x3 F,
                             apf::Mesh * msh,
                             apf::Field * du_,
                             apf::Field * u_)
        : xyz(msh->getCoordinateField())
        , du(du_)
        , u(u_)
        , ent(NULL)
        , hasEntList(false)
    {
      init(msh, F);
    }
    ApplyIncrementalDeformationGradient(std::vector<apf::MeshEntity *>::iterator ent_bgn,
                             std::vector<apf::MeshEntity *>::iterator ent_end,
                             apf::Matrix3x3 F,
                             apf::Mesh * msh,
                             apf::Field * du_,
                             apf::Field * u_)
        : xyz(msh->getCoordinateField())
        , du(du_)
        , u(u_)
        , ent(NULL)
        , ent_bgn(ent_bgn)
        , ent_end(ent_end)
        , hasEntList(true)
    {
      init(msh, F);
    }
    virtual bool inEntity(apf::MeshEntity * m)
    {
      ent = m;
      return true;
    }
    virtual void outEntity() {}
    virtual void atNode(int nd)
    {
      apf::Vector3 nd_u_old;
      // the xyz coordinates in the original frame
      apf::getVector(xyz, ent, nd, nd_xyz);
      apf::getVector(u, ent, nd, nd_u_old);
      // the increment in displacemnt is the increment in displacement gradient
      // multiplied by the current coordinate
      apf::Vector3 nd_du = FmI * (nd_xyz+nd_u_old);
      // the current displacement is the original displacement
      // plus the increment in affine displacement
      apf::Vector3 nd_u = nd_du+nd_u_old;
      apf::setVector(u, ent, nd, nd_u);
      apf::setVector(du, ent, nd, nd_du);
    }
    void run()
    {
      if (hasEntList)
        apply(ent_bgn, ent_end, u);
      else
        apply(u);
    }
  };

  /**
   * templated function to check if two floating point values are close
   * \param a first number to check
   * \param b second number to check
   * \param rtol relative tolerance
   * \param atol absolute tolerance
   */
  template<typename T, typename std::enable_if<std::is_floating_point<T>::value,int>::type=0>
  KOKKOS_INLINE_FUNCTION
  bool isClose(T a, T b, T rtol = 1E-8, T atol = 1E-10)
  {
    return fabs(a - b) <= fmax(rtol * fmax(fabs(a), fabs(b)), atol);
  }
}
#include "Utility_impl.h"
#endif
