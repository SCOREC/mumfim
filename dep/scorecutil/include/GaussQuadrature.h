/******************************************************************************

  (c) 2004-2013 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE-SCOREC file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/

/**
 * @brief     SCORECUtil
 * @file      file name (relative path should be included e.g. /cint/FMDB.h)
 *
 * @brief     Functions to find Gauss integration points.
 *
 * @remark    The functions here can find the quantity and locations
 *  of Gauss integration points on a variety of different element types.
 * They also work for a number of different orders of shape functions but
 *  do have limits on how high the order can be.
 *
 */

#ifndef H_GAUSS_QUADRATURE
#define H_GAUSS_QUADRATURE

#include "IntPt.h"
#include "Mapping.h"

// Hexahedron
/**
 * @brief     Get the locations of Gauss integration points on hexahedral element.
 * @remark    extended description (reserve space)
 *
 * @param     order  (In) order of the shape functions
 * @return    array of Gauss points
 */
IntPt3d *getGQHPts(int order);

/**
 * @brief     Get the number of Gauss integration points on hexahedral element.
 *
 * @param     order  (In) order of the shape functions
 * @return    the number of Gauss points
 */
int getNGQHPts(int order);

// Prism
/**
 * @brief     Get the locations of Gauss integration points on prismatic element.
 * @remark    extended description (reserve space)
 *
 * @param     order  (In) order of the shape functions
 * @return    array of Gauss points
 */
IntPt3d *getGQPrismPts(int order);

/**
 * @brief     Get the number of Gauss integration points on prismatic element.
 *
 * @param     order  (In) order of the shape functions
 * @return    the number of Gauss points
 */
int getNGQPrismPts(int order);

// Pyramid
/**
 * @brief     Get the locations of Gauss integration points on pyramidal element.
 * @remark    extended description (reserve space)
 *
 * @param     order  (In) order of the shape functions
 * @return    array of Gauss points
 */
IntPt3d *getGQPyramidPts(int order);

/**
 * @brief     Get the number of Gauss integration points on pyramidal element.
 *
 * @param     order  (In) order of the shape functions
 * @return    the number of Gauss points
 */
int getNGQPyramidPts(int order);

// Quadralateral
/**
 * @brief     Get the locations of Gauss integration points on quadralateral element.
 * @remark    extended description (reserve space)
 *
 * @param     order  (In) order of the shape functions
 * @return    array of Gauss points
 */
IntPt2d *getGQQPts(int order);

/**
 * @brief     Get the number of Gauss integration points on quadralateral element.
 *
 * @param     order  (In) order of the shape functions
 * @return    the number of Gauss points
 */
int getNGQQPts(int order);

// Tetrahedron
/**
 * @brief     Get the locations of Gauss integration points on tetrahedral element.
 * @remark    extended description (reserve space)
 *
 * @param     order  (In) order of the shape functions
 * @return    array of Gauss points
 */
IntPt3d *getGQTetPts(int order);

/**
 * @brief     Get the number of Gauss integration points on tetrahedral element.
 *
 * @param     order  (In) order of the shape functions
 * @return    the number of Gauss points
 */
int getNGQTetPts(int order);

// Triangle
/**
 * @brief     Get the locations of Gauss integration points on triangular element.
 * @remark    extended description (reserve space)
 *
 * @param     order  (In) order of the shape functions
 * @return    array of Gauss points
 */
IntPt2d *getGQTPts(int order);

/**
 * @brief     Get the number of Gauss integration points on triangular element.
 *
 * @param     order  (In) order of the shape functions
 * @return    the number of Gauss points
 */
int getNGQTPts(int order);

/**
 * @brief     Get the number of Gauss integration points on the passed mesh entity.
 *
 * @param     ent    (In) A mesh entity used as an element
 * @param     order  (In) The order of the shape functions during integration
 *
 * @return    the number of Gauss points
 */
int numIntegrationPoints(int elementType, int order);

#endif
