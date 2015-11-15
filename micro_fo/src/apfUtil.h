#ifndef BIO_APF_UTIL_H_
#define BIO_APF_UTIL_H_

#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <gmi.h>
#include <gmi_null.h>
#include <maMap.h>

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

  /**
   * Calculate the difference between all coordinates of the downwardly-adjacent
   *  vertices of the specified edge.
   * @param msh The mesh containing the edge.
   * @param me The edge.
   * @param diffs The difference between the coordinates of the vertices associated with
   *              the edge.
   */
  void calcEdgeVertDiffs(apf::Mesh * msh, apf::MeshEntity * me, apf::Vector3 & diffs);

  /**
   * Make a 3x3 identity matrix.
   */
  apf::Matrix3x3 eye();

  /**
   * Make a 3x3 matrix containing all 1's
   */
  apf::Matrix3x3 ones();
}

#endif
