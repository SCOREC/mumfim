#ifndef BIO_APF_UTIL_H_
#define BIO_APF_UTIL_H_
#include <apf.h>
#include <apfElement.h>
#include <apfField.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <apfNumbering.h>
#include <gmi.h>
#include <gmi_null.h>
#include <maMap.h>
#include <cassert>
namespace bio
{
  /**
   * Create an empty mesh associated with a 'null' model;
   * @return Pointer to an empty mesh.
   */
  apf::Mesh2 * makeNullMdlEmptyMesh();
  /**
   * Create a mesh containing a single entity related to a 'null' model.
   * @param t The entity type of the single mesh entity.
   * @param vs An array of correctly ordered coordinates defining the vertices of
   *           mesh entity to be created.
   * @return Pointer to a mesh with a single mesh entity.
   */
  apf::Mesh2 * makeSingleEntityMesh(apf::Mesh::Type t, const apf::Vector3 * vs);
  /**
   * Convert a global coordinate to a parametric coordinate on the specified entity.
   * @note This basically only works for linear entities (maybe only tets?)
   * @param msh The mesh containing the mesh entity.
   * @param ent The mesh entity to perform the conversion on.
   * @param glb The global coordinate to convert.
   * @return The barycentric coordinates of the global coordinate.
   */
  apf::Vector3 calcLocalCoord(apf::Mesh * msh,
			      apf::MeshEntity * ent,
			      const apf::Vector3 & glb);
  /**
   * Convert many global coordinates to local coordinates for a single mesh entity.
   * This function avoid multiply-calculating the inversion matrix (for linear elements)
   * which results from calling the calcLocalCoord() function multiple times.
   * @param lcl_crds A dynamic array containing the local coordinates converted from
   *                 global coordinates.
   * @param macro_msh The mesh containing the entity to perform the conversion on.
   * @param macro_ent The mesh entity the local coordinates are in respect to.
   * @param gbl_crds The global coordinates to convert.
   */
  void calcLocalCoords(apf::DynamicArray<apf::Vector3> & lcl_crds,
		       apf::Mesh * macro_msh,
		       apf::MeshEntity * macro_ent,
		       const apf::DynamicArray<apf::Vector3> & gbl_crds);
  /**
   * Calculate the primary measure (lenght, area, volume) for every entity
   * of the specified dimension in the mesh.
   * @param msh The mesh.
   * @param dim The dimensionality of the entities to measure.
   * @param msrs The resulting measures of all entities.
   */
  void calcDimMeasures(apf::Mesh * msh, int dim, std::vector<double> & msrs);
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
  class ApplySolution : public apf::FieldOp
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
  ApplySolution(apf::Numbering * nm, const double * s, int d = 0, bool a = true)
      : apf::FieldOp()
      , num(nm)
      , fld(NULL)
      , fldcmp()
      , sol(s)
      , dof0(d)
      , me(NULL)
      , add(a)
    {
      assert(nm);
      fld = apf::getField(num);
      fldcmp = apf::countComponents(fld);
    }
    bool inEntity(apf::MeshEntity * m) { me = m; }
    void outEntity() {}
    void atNode(int nde)
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
  void fromArray(apf::DynamicVector & to,
		 const double * from,
		 int sz)
  {
    to.resize(sz);
    to.zero();
    for(int ii = 0; ii < sz; ii++)
      to[ii] = from[ii];
  }
  void fromArray(apf::DynamicMatrix & to,
		 const double * from,
		 int nr, int nc)
  {
    to.resize(nr,nc);
    to.zero();
    for(int ii = 0; ii < nr; ii++)
      for(int jj = 0; jj < nc; jj++)
	to(ii,jj) = from[ii*nc + jj];
  }
}
#endif
