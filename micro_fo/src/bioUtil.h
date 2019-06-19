#ifndef BIO_UTIL_H_
#define BIO_UTIL_H_
#include "apfFieldOp.h" //amsi
#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <cassert>
namespace bio
{
  /**
   * Convert a global coordinate to a parametric coordinate on the specified entity.
   * @note This basically only works for linear entities (maybe only tets?)
   * @param msh The mesh containing the mesh entity.
   * @param ent The mesh entity to perform the conversion on.
   * @param glb The global coordinate to convert.
   * @return The barycentric coordinates of the global coordinate.
   */
  void mapGlobalToLocal(apf::Mesh * msh,
                        apf::MeshEntity * e,
                        apf::Vector3 const& global,
                        apf::Vector3 & local);
  /**
   * Convert many global coordinates to local coordinates for a single mesh entity.
   * This function avoids multiply-calculating the inversion matrix (for linear elements)
   * which results from calling the calcLocalCoord() function multiple times.
   * @param lcl_crds A dynamic array containing the local coordinates converted from
   *                 global coordinates.
   * @param macro_msh The mesh containing the entity to perform the conversion on.
   * @param macro_ent The mesh entity the local coordinates are in respect to.
   * @param gbl_crds The global coordinates to convert.
   */
  template <typename I, typename O>
    void mapGlobalsToLocals(apf::Mesh * msh,
                            apf::MeshEntity * e,
                            I gbl_crds_bgn,
                            I gbl_crds_end,
                            O lcl_crds);
  /**
   * Calculate the primary measure (lenght, area, volume) for every entity
   * of the specified dimension in the mesh.
   * @param msh The mesh.
   * @param dim The dimensionality of the entities to measure.
   * @param msrs The resulting measures of all entities.
   */
  template <typename O>
    void calcDimMeasures(apf::Mesh * msh, int dim, O msrs);
  void getCoords(apf::Mesh * msh,
                 apf::MeshEntity ** vrts,
                 apf::Vector3 * crds,
                 int nm);
  double calcDeformedLength(apf::Field * f, apf::MeshEntity * e);
  double calcDistance(const apf::Vector3 & a, const apf::Vector3 & b);
  /**
   * Make a 3x3 identity matrix.
   */
  apf::Matrix3x3 eye();
  /**
   * Make a 3x3 matrix containing all 1's
   */
  apf::Matrix3x3 ones();
  /**
   * Takes a numbering and array of doubles with an initial dof offset and either
   *  sets or accumulates the array values at the index corresponding to the
   *  numbering minus the initial dof offset.
   */
  class ApplySolution : public amsi::FieldOp
  {
  protected:
    apf::Numbering * num;
    apf::Field * fld;
    int fldcmp;
    const double * sol;
    int dof0;
    apf::MeshEntity * me;
    bool add;
  public:
    ApplySolution(apf::Field * fld, apf::Numbering * nm, const double * s, int d = 0, bool a = true)
      : amsi::FieldOp()
      , num(nm)
      , fld(fld)
      , fldcmp()
      , sol(s)
      , dof0(d)
      , me(NULL)
      , add(a)
    {
      assert(nm);
      fldcmp = apf::countComponents(fld);
    }
    void run()
    {
     apply(fld); 
    }
    bool inEntity(apf::MeshEntity * m)
    {
      me = m;
      return true;
    }
    void outEntity() {}
    void atNode(int)
    {
      double cmps[fldcmp];
      apf::getComponents(fld,me,0,&cmps[0]);
      if(add)
        for(int ii = 0; ii < fldcmp; ii++)
        {
          int dof = apf::getNumber(num,me,0,ii);
          cmps[ii] += sol[dof - dof0];
        }
      else
        for(int ii = 0; ii < fldcmp; ii++)
        {
          int dof = apf::getNumber(num,me,0,ii);
          cmps[ii] = sol[dof - dof0];
        }
      apf::setComponents(fld,me,0,&cmps[0]);
    }
  };
  //void fromArray(apf::DynamicVector & to, const double * from, int sz);
  //void fromArray(apf::DynamicMatrix & to, const double * from, int nr, int nc);
  class ApplyDeformationGradient : public amsi::FieldOp
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
    ApplyDeformationGradient(apf::Matrix3x3 F,
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
    ApplyDeformationGradient(std::vector<apf::MeshEntity *>::iterator ent_bgn,
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
      apf::Vector3 nd_u = FmI * nd_xyz;
      // FIXME this is a place where floating point math could kill us
      apf::Vector3 nd_du = nd_u - nd_u_old;
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
}
#include "bioUtil_impl.h"
#endif
