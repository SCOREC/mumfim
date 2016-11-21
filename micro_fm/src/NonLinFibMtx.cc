#include "NonLinFibMtx.h"
#include <Solvers.h>
#include <NonLinearElasticIntegrator.h>
#include <TrussIntegrator.h>
#include <ConvenienceFunctions.h>
#include <apfFunctions.h>
#include <SimMeshTools.h> // E_length
#include <sim.h>
#include <apfSIM.h>
#include <apfShape.h> // countNodesIn
#include <apfNumbering.h> // adjreorder and naivereorder
#include <amsiControlService.h>
namespace Biotissue
{
  using namespace amsi;
  NonLinFibMtx::NonLinFibMtx(pGModel in_model,
                             pParMesh in_mesh,
                             pACase ipd,
                             double in_A, double in_B,
                             double in_E, double in_V,
                             MPI_Comm cm)
    : FEA(cm)
    , apfSimFEA(in_model,in_mesh,ipd,cm)
    , integration_order(1)
    , A(in_A)
    , B(in_B)
    , shear_modulus(in_E)
    , poisson_ratio(in_V)
  {
    apf_primary_field = apf::createFieldOn(apf_mesh,"displacement",apf::VECTOR);
    apf_primary_numbering = apf::createNumbering(apf_primary_field);
    disp_field_shape = apf::getShape(apf_primary_field);
    accumulated_disp = apf::createFieldOn(apf_mesh,"accumulated_displacement",apf::VECTOR);
    Initialize();
  }
  void NonLinFibMtx::Initialize()
  {
    // determine which geometric edges are embedded fibers
    // should only be one geometric region
    GRIter griter = GM_regionIter(model);
    while((pGRegion region = GRIter_next(griter)))
    {
      GEIter eiter = GM_edgeIter(model);
      while((pGEdge edge = GEIter_next(eiter)))
      {
        if(!GEN_inClosure(region, edge))
        {
          embedded_fibers.push_back(edge);
          /*
          // not working at the moment for some reason... but also not used
          double low, high;
          GE_parRange(edge,&low,&high);
          total_embedded_length += GE_length(edge,low,high);
          */
        }
      }
    }
    SetBoundaryGeomFaces();
  }
  void NonLinFibMtx::SetBoundaryGeomFaces()
  {
    pGRegion geom_region = GRIter_next(GM_regionIter(model));
    assert(geom_region);
    pPList geom_faces = GR_faces(geom_region);
    pGFace face;
    void *iter = 0;
    while((face = (pGFace)PList_next(geom_faces,&iter)))
      boundary_geom_faces.push_back(face);
    PList_delete(geom_faces);
  }
  /**
   *@brief Give the initial guess for the very first step. When step zero,
   * use the linear interpolation of displacement along the x direction
   * while zero initial guess for the x and y directions. For subsequent
   * steps, use zero initial guess for all directions.
   * This seems non-generic... should be a more general method
   */
  void NonLinFibMtx::InitialGuess()
  {
    double min[3];
    double max[3];
    GM_bounds(model, min, max);
    std::cout << std::setprecision(16) <<
      "Model max bound: " <<
      std::scientific << max[0] << ", " <<
      std::scientific << max[1] << ", " <<
      std::scientific << max[2] << std::endl;
    std::cout << std::setprecision(16) <<
      "Model min bound: " <<
      std::scientific << min[0] << ", " <<
      std::scientific << min[1] << ", " <<
      std::scientific << min[2] << std::endl;
    double xyz[3];
    VIter viter = M_vertexIter(part);
    pVertex pv;
    while((pv = VIter_next(viter)))
    {
      apf::MeshEntity * vertex = apf::castEntity(static_cast<pEntity>(pv));
      V_coord(pv, xyz);
      apf::Vector3 disp;
      apf::getVector(apf_primary_field,vertex,0,disp);
      if(!apf::isFixed(apf_primary_numbering,vertex,0,0))
      {
        if(load_step == 0)
          disp[0] = 0.00005542834123703 / 2 / (max[0]-min[0]) * xyz[0]
            - 0.00005542834123703 / 2 / (max[0]-min[0]) * min[0]; // 689fibers
        else
          disp[0] = 0.0;
        disp[2] = disp[1] = 0.0;
        apf::setVector(apf_primary_field,vertex,0,disp);
      }
    }
    // when step 0, need to put the displacement from both the initial guess and boundary condition into the history 1 of the displacement field
    if(load_step == 0)
      AccumulateDisplacement();
  }
  void NonLinFibMtx::DirectionAndStretch(apf::Vector3 & net_direction,
                                         double & mean_stretch)
  {
    double xyz0[3], xyz1[3];
    double total_length = 0.0;
    mean_stretch = net_direction[2] = net_direction[1] = net_direction[0] = 0.0;
    int total_num_mesh_edges = 0;
    for(std::list<pGEntity>::iterator iter = embedded_fibers.begin(),
          iterend = embedded_fibers.end();
        iter != iterend; ++iter)
    {
      std::list<pEntity> mesh_edges;
      Model_GetClassifiedEntities(part,*iter,1,mesh_edges);
      total_num_mesh_edges += mesh_edges.size();
      for(std::list<pEntity>::iterator mesh_edge_iter = mesh_edges.begin(),
            mesh_edge_iterend = mesh_edges.end();
          mesh_edge_iter != mesh_edge_iterend; mesh_edge_iter++)
      {
        pEdge edge = static_cast<pEdge>(*mesh_edge_iter);
        pVertex vtx0, vtx1;
        vtx0 = E_vertex(edge,0);
        vtx1 = E_vertex(edge,1);
        V_coord(vtx0, xyz0);
        V_coord(vtx1, xyz1);
        apf::Vector3 pos0(xyz0[0],xyz0[1],xyz0[2]);
        apf::Vector3 pos1(xyz1[0],xyz1[1],xyz1[2]);
        apf::Vector3 disp0, disp1;
        apf::getVector(apf_primary_field,apf::castEntity(static_cast<pEntity>(vtx0)),0,disp0);
        apf::getVector(apf_primary_field,apf::castEntity(static_cast<pEntity>(vtx1)),0,disp1);
        pos0 = pos0 + disp0;
        pos1 = pos1 + disp1;
        double length = sqrt((pos0[0]-pos1[0])*(pos0[0]-pos1[0]) +
                             (pos0[1]-pos1[1])*(pos0[1]-pos1[1]) +
                             (pos0[2]-pos1[2])*(pos0[2]-pos1[2]));
        double orig_length = E_length(edge);
        mean_stretch += length / orig_length;
        total_length += length;
        net_direction[0] += ((pos0[0]-pos1[0])*(pos0[0]-pos1[0])) / length;
        net_direction[1] += ((pos0[1]-pos1[1])*(pos0[1]-pos1[1])) / length;
        net_direction[2] += ((pos0[2]-pos1[2])*(pos0[2]-pos1[2])) / length;
      }
    }
    mean_stretch /= total_num_mesh_edges;
    net_direction[0] /= total_length;
    net_direction[1] /= total_length;
    net_direction[2] /= total_length;
  }
  void NonLinFibMtx::BoundaryFiberEdges(std::list<pEdge> & boundary_edges)
  {
    std::list<pGEntity> embedded_fibers;
    std::list<pGVertex> boundary_geom_vertices;
    std::list<pGEdge> boundary_geom_edge;
    std::list<pVertex> boundary_vertices;
    for(std::list<pGEntity>::iterator fiber_iter = embedded_fibers.begin(),
          fiber_iterend = embedded_fibers.end();
        fiber_iter != fiber_iterend; fiber_iter++)
    {
      pGEdge geom_edge = static_cast<pGEdge>(*fiber_iter);
      pGVertex geom_vtx0 = GE_vertex(geom_edge,0);
      pGVertex geom_vtx1 = GE_vertex(geom_edge,1);
      GRIter griter = GM_regionIter(model);
      while((pGRegion geom_region = GRIter_next(griter)))
      {
        bool vtx0_in_closure, vtx1_in_closure;
        if(vtx0_in_closure = GR_inClosure(geom_region,geom_vtx0))
        {
          boundary_vertices.push_back(M_classifiedVertex(part,geom_vtx0));
          boundary_geom_vertices.push_back(geom_vtx0);
        }
        if(vtx1_in_closure = GR_inClosure(geom_region,geom_vtx1))
        {
          boundary_vertices.push_back(M_classifiedVertex(part,geom_vtx1));
          boundary_geom_vertices.push_back(geom_vtx1);
        }
        if(vtx0_in_closure || vtx1_in_closure)
          boundary_geom_edge.push_back(geom_edge);
      }
    }
    std::list<pEntity> mesh_edges;
    for(std::list<pGEdge>::iterator geom_edge_iter = boundary_geom_edge.begin(),
          geom_edge_iter_end = boundary_geom_edge.end();
        geom_edge_iter != geom_edge_iter_end; geom_edge_iter++)
    {
      Model_GetClassifiedEntities(part,*geom_edge_iter,1,mesh_edges);
      for(std::list<pEntity>::iterator edge_iter = mesh_edges.begin(),
            edge_iter_end = mesh_edges.end();
          edge_iter != edge_iter_end; edge_iter++)
      {
        pEdge edge = static_cast<pEdge>(*edge_iter);
        pVertex vtx0, vtx1;
        vtx0 = E_vertex(edge,0);
        vtx1 = E_vertex(edge,1);
        bool vtx0_in = std::find(boundary_vertices.begin(),
                                 boundary_vertices.end(),
                                 vtx0) != boundary_vertices.end();
        bool vtx1_in = std::find(boundary_vertices.begin(),
                                 boundary_vertices.end(),
                                 vtx1) != boundary_vertices.end();
        bool edge_in = std::find(boundary_edges.begin(),
                                 boundary_edges.end(),
                                 edge) != boundary_edges.end();
        if((vtx0_in || vtx1_in) && !edge_in)
          boundary_edges.push_back(edge);
      }
    }
  }
  void NonLinFibMtx::GetAverageFiberStress(LAS * las,
                                           apf::Matrix3x3 & fiber_stress)
  {
    int global_size, local_size, offset;
    GetDOFInfo(global_size,local_size,offset);
    las->Reinitialize(local_size, global_size, offset);
    las->ZeroVector();
    //todo: either uncomment this or calculate the average fiber stress without this 'hack'
//      FiberForceVector(las);
    // old code only got the 'right' model face???
    std::list<pEntity> boundary_vertices;
    GFIter gfiter = GM_faceIter(model);
    pGFace gf;
    while ((gf = GFIter_next(gfiter)))
      Model_GetClassifiedEntities(part,static_cast<pGEntity>(gf),0,boundary_vertices);
    fiber_stress = fiber_stress * 0.0;
    double * RHS;
    las->GetVector(RHS);
    apf::Vector3 total_tractions;
    std::list<pEntity>::iterator iter = boundary_vertices.begin();
    std::list<pEntity>::iterator iter_end = boundary_vertices.end();
    for(; iter != iter_end; iter++)
    {
      pVertex vtx = static_cast<pVertex>(*iter);
      double xyz[3];
      V_coord(vtx, xyz);
      // assuming 3d analysis
      for(int ii = 0; ii < 3; ii++)
      {
        int global_number = apf::getNumber(apf_primary_numbering, apf::castEntity(*iter),0,ii);
        total_tractions[ii] += RHS[global_number];
        for(int jj = 0; jj < 3; jj++)
          fiber_stress[ii][jj] += RHS[global_number] * xyz[ii];
      }
      std::cout.precision(16);
      std::cout << "Total traction components: "
                << std::scientific << total_tractions[0] << ", "
                << std::scientific << total_tractions[1] << ", "
                << std::scientific << total_tractions[2] << std::endl;
    }
  }
  void NonLinFibMtx::GetAverageMatrixStress(apf::Matrix3x3 & matrix_stress)
  {
    matrix_stress = matrix_stress * 0.0;
    apf::Vector3 total_traction;
    std::list<pEntity> boundary_mesh_faces;
    std::list<pRegion> boundary_mesh_regions;
    int geom_face_index = 0;
    GFIter gfit = GM_faceIter(model);
    for(std::vector<pGEntity>::iterator geom_face_iter = boundary_geom_faces.begin(),
          geom_face_iterend = boundary_geom_faces.end();
        geom_face_iter != geom_face_iterend; geom_face_iter++)
    {
      // get mesh faces on this geometric face
      Model_GetClassifiedEntities(part,*geom_face_iter,2,boundary_mesh_faces);
      std::list<pEntity>::iterator face_iter = boundary_mesh_faces.begin();
      std::list<pEntity>::iterator face_iterend = boundary_mesh_faces.end();
      for(; face_iter != face_iterend; face_iter++)
      {
        pFace face = static_cast<pFace>(*face_iter);
        double area = F_area(face);
        pPList region_list = F_regions(face);
        //assumes triangle... true for linear tets but not general
        double verts[3][3];
        F_vertexCoords(face,&verts[0],1);
        apf::Vector3 a(verts[0]);
        apf::Vector3 b(verts[1]);
        apf::Vector3 c(verts[2]);
        apf::Vector3 V(b-a);
        apf::Vector3 W(c-a);
        apf::Vector3 face_norm;
        face_norm[0] = (V.y() * W.z()) - (V.z() * W.y());
        face_norm[1] = (V.z() * W.x()) - (V.x() * W.z());
        face_norm[2] = (V.x() * W.y()) - (V.y() * W.x());
        void * it = 0;
        while((pRegion region = static_cast<pRegion>(PList_next(region_list, &it))))
        {
          apf::MeshEntity * me = apf::castEntity(static_cast<pEntity>(region));
          apf::MeshElement * melm = apf::createMeshElement(apf_mesh,me);
          apf::Element * e = apf::createElement(apf_primary_field,me);
          int num_gauss_pts = apf::countIntPoints(melm,1);
          apf::NewArray<apf::Vector3> displacements;
          apf::getVectorNodes(e,displacements);
          for(int gauss_index = 0; gauss_index < num_gauss_pts; gauss_index++)
          {
            apf::Vector3 gauss_pt;
            apf::getIntPoint(melm,1,gauss_index,gauss_pt);
            double w = apf::getIntWeight(melm,1,gauss_index);
            apf::Matrix3x3 J;
            apf::getJacobian(melm,gauss_pt,J);
            double det_jac = apf::getJacobianDeterminant(J,apf::getDimension(melm));
            apf::NewArray<apf::Vector3> shape_derivs;
            apf::getShapeGrads(e,gauss_pt,shape_derivs);
            apf::Matrix3x3 def_grad;
            double det_def_grad;
            DeformationGradient(e,
                                gauss_pt,
                                def_grad,
                                det_def_grad);
            apf::Matrix3x3 C;
            RightCauchy(def_grad,C);
            apf::Matrix3x3 pk2;
            PK2Stress(C,
                      det_def_grad,
                      poisson_ratio,
                      shear_modulus,
                      pk2);
            // assuming 3d analysis
            apf::Vector3 traction;
            for(int ii = 0; ii < field_components; ii++)
            {
              traction[ii] = pk2 * face_norm * Identity3x3[ii];
              total_traction[ii] += traction[ii] * area;
            }
            apf::Vector3 face_values;
            SurfItg(face, face_values);
            for(int ii = 0; ii < 3; ii++)
              for(int jj = 0; jj < 3; jj++)
                matrix_stress[ii][jj] += traction[jj] * face_values[ii];
          }
        }
        PList_delete(region_list);
      }
    }
    // assuming 1 geometric region
    double vol;
    RIter riter = M_regionIter(part);
    pRegion rgn = RIter_next(riter);
    vol = R_volume(rgn);
    std::cout << std::setprecision(16) <<
      "Total traction components: " <<
      std::scientific << total_traction[0] << ", " <<
      std::scientific << total_traction[1] << ", " <<
      std::scientific << total_traction[2] << std::endl;
  }
  // surface integration of icomp coordinate value on face
  void NonLinFibMtx::SurfItg(pFace face, apf::Vector3 & values)
  {
    values = values * 0.0;
    apf::MeshEntity * me = apf::castEntity(static_cast<pEntity>(face));
    apf::MeshElement * melm = apf::createMeshElement(apf_mesh,me);
    apf::Element * e = apf::createElement(apf_primary_field,me);
    int num_gauss_pts = apf::countIntPoints(melm,1);
    for(int gauss_index = 0; gauss_index < num_gauss_pts; gauss_index++)
    {
      apf::Vector3 gauss_pt;
      apf::getIntPoint(melm,1,gauss_index,gauss_pt);
      double w = apf::getIntWeight(melm,1,gauss_index);
      apf::Matrix3x3 J;
      apf::getJacobian(melm,gauss_pt,J);
      double det_jac = apf::getJacobianDeterminant(J,apf::getDimension(melm));
      //apf::NewArray<apf::Vector3> shape_vals;
      //apf::getShapeValues(e,gauss_pt,shape_vals);
      apf::Vector3 global_coords;
      apf::mapLocalToGlobal(melm,gauss_pt,global_coords);
      //assumes 3d analysis
      for(int ii = 0; ii < 3; ii++) // num_field_components...
        values[ii] += w * det_jac * global_coords[ii];
    }
    // translate to global coords?
    //mapping -> eval(pos[0], pos[1], pos[2], x[0], x[1], x[2]);
  }
// accumulate iterative displacements into total displacement field
  void NonLinFibMtx::AccumulateDisplacement()
  {
    // iterate over dimensions
    for(int ii = 0; ii <= 3; ii++)
    {
      if(!disp_field_shape->hasNodesIn(ii))
        break;
      apf::MeshIterator * iter = apf_mesh->begin(ii);
      apf::MeshEntity * current_ent;
      while((current_ent = apf_mesh->iterate(iter)))
      {
        int ent_type = apf_mesh->getType(current_ent);
        int num_nodes = disp_field_shape->countNodesOn(ent_type);
        for(int jj = 0; jj < num_nodes; jj++)
        {
          apf::Vector3 disp;
          apf::getVector(apf_primary_field,current_ent,jj,disp);
          apf::Vector3 accum_disp;
          apf::getVector(accumulated_disp,current_ent,jj,accum_disp);
          // to prevent accumulating the essential bc's.. this wasn't in the orginal
          for(int kk = 0; kk < 3; kk++)
          {
            if(!apf::isFixed(apf_primary_numbering, current_ent, jj, kk))
              accum_disp[kk] += disp[kk];
          }
          apf::setVector(accumulated_disp,current_ent,jj,accum_disp);
        }
      }
    }
  }
  // todo: pull out into support function file
  void NonLinFibMtx::Config_CheckRegionValidity()
  {
    RIter riter = M_regionIter(part);
    pRegion rgn;
    while((rgn = RIter_next(riter)))
    {
      int is_good = R_isValidElement(rgn);
      if(is_good == 0)
      {
        std::cout << "INVALID REGION!" << std::endl;
        double vol = R_volume(rgn);
        std::cout.precision(16);
        std::cout << "Volume is " << std::scientific << vol << std::endl;
        int region_id = EN_id(rgn);
        std::cout << "Region ID is " << region_id << std::endl;
        pPList vlist = PList_new();
        vlist = R_vertices(rgn, 1);
        void *iter = 0;
        pEntity ent;
        std::cout << "Loop over the vertices to see the vertex classification: " << std::endl;
        while((ent = static_cast<pEntity>(PList_next(vlist, &iter))))
        {
          int vtxID = EN_id(ent);
          std::cout << "Vertex ID is " << vtxID << std::endl;
          pGEntity gent = EN_whatIn(ent);
          int gtype = GEN_type (gent);
          std::cout << "Geometric type is " << gtype << std::endl;
          int genTag = GEN_tag(gent);
          std::cout << "Geometric tag is " << genTag << std::endl;
          if(gtype == 1)
          {
            pGVertex vtx0 = GE_vertex ((pGEdge)gent,0);
            pGVertex vtx1 = GE_vertex ((pGEdge)gent,1);
            std::cout << "Edge tag " << genTag << " has vertices tags " << GEN_tag(vtx0) << " " << GEN_tag(vtx1) << std::endl;
          }
        }
      }
    }
  }
  void NonLinFibMtx::Run(LAS * las,
                         int max_load_step)
  {
    for(; load_step < max_load_step; load_step++)
    {
      std::cout << "Begin solving load_step " << load_step << ":" << std::endl;
      double residual_norm = 0.0;
      NewtonSolver(this,las,100,1.e-14,1.0,residual_norm);
      displaceMesh(apf_primary_field);
    }
    AverageStress(total_stress);
    UnbalancedForce(unbalanced_force);
  }
  void NonLinFibMtx::ComputeMatrix(LAS * las)
  {
    las->ZeroMatrix();
    NonLinearElasticIntegrator matrix_elemental_system(apf_primary_field,
                                                       integration_order,
                                                       poisson_ratio,
                                                       shear_modulus);
    TrussIntegrator fiber_elemental_system(apf_primary_field,
                                           integration_order,
                                           A,B);
    apf::DynamicMatrix & mKe = matrix_elemental_system.Ke();
    apf::DynamicVector & mFe = matrix_elemental_system.Fe();
    apf::DynamicMatrix & fKe = fiber_elemental_system.Ke();
    apf::DynamicVector & fFe = fiber_elemental_system.Fe();
    GRIter griter = GM_regionIter(model);
    while((pGEntity gregion = GRIter_next(griter)))
    {
      std::list<pEntity> mesh_regions;
      Model_GetClassifiedEntities(part,gregion,3,mesh_regions);
      int num_local_regions = mesh_regions.size();
      for(std::list<pEntity>::iterator iter = mesh_regions.begin(),
            iterend = mesh_regions.end(); iter != iterend; iter++)
      {
        apf::MeshEntity * me = apf::castEntity(*iter);
        apf::MeshElement * melm = apf::createMeshElement(apf_mesh,
                                                         me);
        matrix_elemental_system.process(melm);
        apf::NewArray<int> dof_numbers;
        apf::getElementNumbers(apf_primary_numbering,
                               me,
                               dof_numbers);
        FEA::AssembleDOFs(las,
                          (int)mKe.getColumns(),
                          &dof_numbers[0],
                          &matrix_elemental_system.NodeValues()[0],
                          &mKe(0,0),
                          &mFe(0),
                          true);
      }
    }
    // Iterate over the list of embedded geometric edges
    for(std::list<pGEntity>::iterator iter = embedded_fibers.begin(),
          iterend = embedded_fibers.end(); iter != iterend; iter++)
    {
      std::list<pEntity> edges;
      Model_GetClassifiedEntities(part,*iter,1,edges);
      for(std::list<pEntity>::iterator eiter = edges.begin(),
            eiterend = edges.end(); eiter != eiterend; eiter++)
      {
        apf::MeshEntity * me = apf::castEntity(*eiter);
        apf::MeshElement * melm = apf::createMeshElement(apf_mesh,me);
        fiber_elemental_system.process(melm);
        apf::NewArray<int> dof_numbers;
        apf::getElementNumbers(apf_primary_numbering,
                               me,
                               dof_numbers);
        AssembleDOFs(las,
                     (int)fKe.getColumns(),
                     &dof_numbers[0],
                     &fiber_elemental_system.NodeValues()[0],
                     &fKe(0,0),
                     &fFe(0),
                     true);
      }
    }
    if(nonlinear_iteration == 0)
      AccumulateDisplacement();
    nonlinear_iteration++;
  }
  void NonLinFibMtx::ComputeVector(LAS * las, double & residual_norm)
  {
    if(nonlinear_iteration == 0)
      InitialGuess();
    las->GetVectorNorm(residual_norm);
  }
/*
  void NonLinFibMtx::MatrixStiffnessMatrix(LAS * las)
  {
  GRIter griter = GM_regionIter(model);
  while (pGEntity gregion = GRIter_next(griter))
  {
  std::list<pEntity> mesh_regions;
  Model_GetClassifiedEntities(part,gregion,3,mesh_regions);
  int num_local_regions = mesh_regions.size();
  // Assembly begins
  for(std::list<pEntity>::iterator iter = mesh_regions.begin(),
  iterend = mesh_regions.end(); iter != iterend; iter++)
  {
  apf::MeshEntity * me = apf::castEntity(*iter);
  apf::MeshElement * melm = apf::createMeshElement(apf_mesh,me);
  apf::Element * e = apf::createElement(apf_mesh->getCoordinateField(),melm);
  // Get the element stiffness matrix for the region
  apf::DynamicMatrix Ke;
  MatrixElementStiffnessMatrix(me,Ke);
  apf::NewArray<int> dof_numbers;
  apf::getElementNumbers(apf_primary_numbering,me,dof_numbers);
  apf::NewArray<apf::Vector3> node_values;
  apf::getVectorNodes(e,node_values);
  AssembleDOFs(las,
  Ke.getColumns(),
  &dof_numbers[0],
  &node_values[0],
  &Ke(0,0),
  NULL);
  }
  }
  }
*/
/*
  void NonLinFibMtx::MatrixForceVector(LAS * las)
  {
  las->ZeroVector();
  GRIter griter = GM_regionIter(model);
  while(pGEntity gregion = GRIter_next(griter))
  {
  std::list<pEntity> mesh_regions;
  Model_GetClassifiedEntities(part,gregion,3,mesh_regions);
  // Iterate over all region entities
  for(std::list<pEntity>::iterator iter = mesh_regions.begin(),
  iterend = mesh_regions.end(); iter != iterend; iter++)
  {
  apf::MeshEntity * me = apf::castEntity(*iter);
  apf::DynamicVector fe;
  MatrixElementForceVector(me,fe);
  // assemble the result into the global linear system
  }
  }
  }
  void NonLinFibMtx::FiberStiffnessMatrix(LAS * las)
  {
  // Iterate over the list of embedded geometric edges
  for(std::list<pGEntity>::iterator iter = embedded_fibers.begin(),
  iterend = embedded_fibers.end(); iter != iterend; iter++)
  {
  std::list<pEntity> edges;
  Model_GetClassifiedEntities(part,*iter,1,edges);
  for(std::list<pEntity>::iterator eiter = edges.begin(),
  eiterend = edges.end(); eiter != eiterend; eiter++)
  {
  apf::MeshEntity * me = apf::castEntity(*eiter);
  apf::DynamicMatrix Ke;
  TrussElementStiffness(me,Ke);
  // assemble the result into the global linear system
  }
  }
  }
  void NonLinFibMtx::FiberForceVector(LAS * las)
  {
  for(std::list<pGEntity>::iterator geom_iter = embedded_fibers.begin(),
  geom_iterend = embedded_fibers.end(); geom_iter != geom_iterend; geom_iter++)
  {
  std::list<pEntity> edges;
  Model_GetClassifiedEntities(part,*geom_iter,1,edges);
  for(std::list<pEntity>::iterator iter = edges.begin(),
  iterend = edges.end(); iter != iterend; iter++)
  {
  apf::MeshEntity * me = apf::castEntity(*iter);
  apf::DynamicVector fe;
  TrussElementForce(me,fe);
  // assemble result into the global linear system
  }
  }
  }
  void NonLinFibMtx::AssembleFiberNetwork(LAS * las,
  pGEntity gentity)
  {
  las->Zero();
  for(std::list<pGEntity>::iterator geom_iter = embedded_fibers.begin(),
  geom_iterend = embedded_fibers.end(); geom_iter != geom_iterend; geom_iter++)
  {
  std::list<pEntity> edges;
  Model_GetClassifiedEntities(part,*geom_iter,1,edges);
  for(std::list<pEntity>::iterator iter = edges.begin(),
  iterend = edges.end(); iter != iterend; iter++)
  {
  apf::MeshEntity * me = apf::castEntity(*iter);
  //TrussElementForce (Ent,  d_ElmtFrcVec, dfmFibLen, FibFrc, DfrcDlen);
  //TrussElementStiffness (Ent, dfmFibLen,  FibFrc,  DfrcDlen, d_ElmtTgtStfs);
  // retrieve dof numbers
  // assemble force vector and stiffness matrix
  }
  }
  }
*/
/*
void NonLinFibMtx::MatrixElementStiffnessMatrix(apf::MeshEntity * me,
apf::DynamicMatrix & Ke)
{
apf::MeshElement * melm = apf::createMeshElement(apf_mesh,me);
apf::Element * e = apf::createElement(apf_mesh->getCoordinateField(),melm);
int num_nodes = apf::countNodes(e);
int num_dofs = num_nodes * apf::countComponents(apf_primary_field);
Ke.setSize(num_dofs,num_dofs);
int num_gauss_pts = apf::countIntPoints(melm,1);
apf::NewArray<apf::Vector3> displacements;
apf::getVectorNodes(e,displacements);
for(int gauss_index = 0; gauss_index < num_gauss_pts; gauss_index++)
{
apf::Vector3 gauss_pt;
apf::getIntPoint(melm,1,gauss_index,gauss_pt);
double w = apf::getIntWeight(melm,1,gauss_index);
apf::Matrix3x3 J;
apf::getJacobian(melm,gauss_pt,J);
double det_jac = apf::getJacobianDeterminant(J,apf::getDimension(melm));
apf::NewArray<apf::Vector3> shape_derivs;
apf::getShapeGrads(e,gauss_pt,shape_derivs);
apf::Matrix3x3 grad;
double grad_det;
DeformationGradient(me,gauss_pt,grad,grad_det);
apf::Matrix3x3 C;
RightCauchy(grad,C);
// todo: pull version in macro up to femanalysis
apf::Matrix<6,6> hookes;
LinearizedConstitutive(C,det_jac,hookes);
apf::Matrix3x3 pk2;
PK2Stress(C,det_jac,pk2);
apf::DynamicMatrix BL;
LinearStrainDisp(num_nodes,
shape_derivs,
displacements,
BL);
apf::DynamicMatrix BNL;
NonLinearStrainDisp(num_nodes,
shape_derivs,
BNL);
apf::DynamicMatrix K0;
LinearStiffnessK0(num_nodes,
BL,
hookes,
K0);
apf::DynamicMatrix KNL;
NonLinearStiffness(num_nodes,
BNL,
pk2,
KNL);
K0 *= w * det_jac;
KNL *= w * det_jac;
Ke += K0;
Ke += KNL;
}
}
*/
 /*
   void NonLinFibMtx::MatrixElementForceVector(apf::MeshEntity * me,
   apf::DynamicVector & fe)
   {
   apf::MeshElement * melm = apf::createMeshElement(apf_mesh,me);
   apf::Element * e = apf::createElement(apf_primary_field,me);
   int num_nodes = apf::countNodes(e);
   int num_dofs = num_nodes * apf::countComponents(apf_primary_field);
   fe.setSize(num_dofs);
   apf::NewArray<apf::Vector3> displacements;
   apf::getVectorNodes(e,displacements);
   int num_gauss_pts = apf::countIntPoints(melm,1);
   for(int gauss_id = 0; gauss_id < num_gauss_pts; gauss_id++)
   {
   apf::Vector3 gauss_pt;
   apf::getIntPoint(melm,1,gauss_id,gauss_pt);
   double w = apf::getIntWeight(melm,1,gauss_id);
   apf::Matrix3x3 J;
   apf::getJacobian(melm,gauss_pt,J);
   double det_jac = apf::getJacobianDeterminant(J,apf::getDimension(melm));
   double wxdetjac = w * det_jac;
   apf::NewArray<apf::Vector3> shape_derivs;
   apf::getShapeGrads(e,gauss_pt,shape_derivs);
   apf::Matrix3x3 grad;
   double grad_det;
   DeformationGradient(e,
   gauss_pt,
   grad,
   grad_det);
   apf::Matrix3x3 C;
   RightCauchy(grad,C);
   apf::Matrix3x3 pk2;
   PK2Stress(C,
   det_jac,
   poisson_ratio,
   shear_modulus,
   pk2);
   apf::Vector<6> pk2_vec;
   SymmMatrixToVector(pk2,pk2_vec);
   apf::DynamicMatrix BL;
   LinearStrainDisp(num_nodes,
   shape_derivs,
   displacements,
   BL);
   for(int ii = 0; ii < num_dofs; ii++)
   for(int jj = 0; jj < 6; jj++)
   fe[ii] += BL(jj,ii) * pk2_vec[jj];
   fe *= (-wxdetjac);
   }
   }
 */
void NonLinFibMtx::TrussElementForce(apf::MeshEntity * edge,
                                     apf::DynamicVector & fe)
  {
    //todo: need access to original coords and length...
    pEdge entity = static_cast<pEdge>(reinterpret_cast<pEntity>(edge));
    pVertex vtx0 = E_vertex(entity, 0);
    pVertex vtx1 = E_vertex(entity, 1);
    // Get coordinates for the two mesh vertices
    double xyz0[3], xyz1[3];
    V_coord(vtx0, xyz0);
    V_coord(vtx1, xyz1);
    // Get undeformed length of the mesh edge
    double init_length = E_length(entity);
    apf::MeshEntity * me = apf::castEntity(entity);
    apf::Element * e = apf::createElement(apf_primary_field,me);
    apf::NewArray<apf::Vector3> disps;
    apf::getVectorNodes(e,disps);
    for(int ii = 0; ii < 3; ii++)
    {
      xyz0[ii] += disps[0][ii];
      xyz1[ii] += disps[1][ii];
    }
    // compute the deformed length of the mesh edge
    double current_length = sqrt(pow((xyz1[0]-xyz0[0]),2) +
                                 pow((xyz1[1]-xyz0[1]),2) +
                                 pow((xyz1[2]-xyz0[2]),2) );
    // Compute the stretch ratio
    double stretch_ratio = current_length / init_length;
    // Compute fiber axial force and related force derivative
    double fiber_force, dforce_dlength;
    FiberAxialForce(A,B,
                    stretch_ratio,
                    init_length,
                    dforce_dlength,
                    fiber_force);
    //[x1,x2,y1,y2,z1,z2]
    fe[0] = -fiber_force * (xyz0[0] - xyz1[0]) / current_length;
    fe[1] = -fe[0];
    fe[2] = -fiber_force * (xyz0[1] - xyz1[1]) / current_length;
    fe[3] = -fe[2];
    fe[4] = -fiber_force * (xyz0[2] - xyz1[2]) / current_length;
    fe[5] = -fe[4];
  }
  /*
    void NonLinFibMtx::TrussElementStiffness(apf::MeshEntity * me,
    apf::DynamicMatrix & Ke)
    {
    // todo: need access to original coords and length
    pEdge entity = static_cast<pEdge>(reinterpret_cast<pEntity>(me));
    double xyz0[3], xyz1[3];
    // get initial coordinates
    pVertex vtx0 = E_vertex(entity, 0);
    pVertex vtx1 = E_vertex(entity, 1);
    V_coord(vtx0, xyz0);
    V_coord(vtx1, xyz1);
    // Get initial/undeformed fiber length
    double init_len = E_length(entity);
    apf::Element * e = apf::createElement(apf_primary_field,me);
    apf::NewArray<apf::Vector3> disps;
    apf::getVectorNodes(e,disps);
    for(int ii = 0; ii < 3; ii++)
    {
    xyz0[ii] += disps[0][ii];
    xyz1[ii] += disps[1][ii];
    }
    int num_nodes = apf::countNodes(e);
    int field_components = apf::countComponents(apf_primary_field);
    int num_dofs = num_nodes * field_components;
    Ke.setSize(num_dofs,num_dofs);
    // Get deformed fiber length
    double current_len = sqrt(pow((xyz1[0]-xyz0[0]),2) +
    pow((xyz1[1]-xyz0[1]),2) +
    pow((xyz1[2]-xyz0[2]),2) );
    // Compute fiber stretch ratio
    double stretch_ratio = current_len / init_len;
    // Compute fiber axial force and related force derivative
    double fiber_force, dforce_dlength;
    FiberAxialForce(A,B,
    stretch_ratio,
    init_len,
    dforce_dlength,
    fiber_force);
    double force_per_length = fiber_force / current_len;
    double difference = dforce_dlength - force_per_length;
  */
/*
  [
  Ke(0,0)  Ke(0,1)  Ke(0,2) -Ke(0,0) -Ke(0,1) -Ke(0,2)
  Ke(1,0)  Ke(1,1)  Ke(1,2) -Ke(1,0) -Ke(1,1) -Ke(1,2)
  Ke(2,0)  Ke(2,1)  Ke(2,2) -Ke(2,0) -Ke(2,1) -Ke(2,2)
  -Ke(0,0) -Ke(0,1) -Ke(0,2)  Ke(0,0)  Ke(0,1)  Ke(0,2)
  -Ke(1,0) -Ke(1,1) -Ke(1,2)  Ke(1,0)  Ke(1,1)  Ke(1,2)
  -Ke(2,0) -Ke(2,1) -Ke(2,2)  Ke(2,0)  Ke(2,1)  Ke(2,2)
  ]
*/
  /*
  // i=0..2 j=0..2
  for(int ii = 0; ii < field_components; ii++)
  for(int jj = 0; jj < field_components; jj++)
  {
  double to_add = ii == jj ? force_per_length : 0.0;
  Ke(ii,jj) = difference * (xyz0[ii] - xyz1[ii]) / current_len
  * (xyz0[jj] - xyz1[jj]) / current_len + to_add;
  int & n = field_components;
  Ke(ii+n,jj) = -Ke(ii,jj);
  Ke(ii,jj+n) = -Ke(ii,jj);
  Ke(ii+n,jj+n) = Ke(ii,jj);
  }
  }
  */
  /*****************************************************************************************
The fiber material law is:
F = A * L^2 * ( exp( B * GS ) - 1);
where:
A: fiber force parameter -- A =  E*Af/B
B: fiber unitless parameter
GS: fiber green strain  -- not exactly green strain -- missing 1/2 factor
L: fiber stretch ratio
  ****************************************************************************************/
// todo: pull out of the class and into a seperate utility function, since it doesn't depend on the class data members at all (except Acoeff and BCoeff which can be passed)
/*    void NonLinFibMtx::FiberAxialForce(double A, double B,
      double stretch_ratio,
      double init_length,
      double & dforce_dlength,
      double & fiber_force)
      {
      double ratio_limit = 1.4;
      double ratio = stretch_ratio > ratio_limit ? ratio_limit : stretch_ratio;
      double exp_term = exp(B * (ratio * ratio - 1));
      double ratio_sqrd = ratio * ratio;
      // = A^2 s^2 exp(B (s^2 - 1))
      fiber_force = (A * ratio_sqrd) * exp_term;
      // = 2 A s/l ( exp(B (s^2 - 1)) - 1 + B s^2 exp(B (s^2 - 1)) )
      dforce_dlength = 2 * A * (ratio / init_length) *
      ( exp_term - 1 + B * ratio_sqrd * exp_term );
      }
*/
/*
  void NonLinFibMtx::DeformationGradient(apf::MeshEntity * me,
  const apf::Vector3 & p,
  apf::Matrix3x3 & grad,
  double & grad_det)
  {
  apf::NewArray<apf::Vector3> u;
  apf::MeshElement * melm = apf::createMeshElement(apf_mesh,me);
  apf::Element * e = apf::createElement(apf_mesh->getCoordinateField(),melm);
  apf::getVectorNodes(e,u);
  apf::NewArray<apf::Vector3> shape_derivs;
  apf::getShapeGrads(e,p,shape_derivs);
  grad = Identity;
  int num_nodes = apf::countNodes(e);
  int field_components = apf::countComponents(apf_primary_field);
  for(int ii = 0; ii < field_components; ii++)
  for(int jj = 0; jj < field_components; jj++)
  for(int kk = 0; kk < num_nodes; kk++)
  grad[ii][jj] += u[kk][ii] * shape_derivs[kk][jj];
  grad_det = apf::getDeterminant(grad);
  }
  // todo: pull out into seperate function file
  void NonLinFibMtx::LeftCauchy(const apf::Matrix3x3 & F,
  apf::Matrix3x3 & B)
  {
  for(int ii = 0; ii < 3; ii++)
  for(int jj = 0; jj < 3; jj++)
  {
  B[ii][jj] = 0.0;
  for(int kk = 0; kk < 3; kk++)
  B[ii][jj] += F[ii][kk] * F[jj][kk];
  }
  }
  // todo: pull into seperate function file
  void NonLinFibMtx::RightCauchy(const apf::Matrix3x3 & F,
  apf::Matrix3x3 & C)
  {
  for(int ii = 0; ii < 3; ii++)
  for(int jj = 0; jj < 3; jj++)
  {
  C[ii][jj] = 0.0;
  for(int kk = 0; kk < 3; kk++)
  C[ii][jj] += F[kk][ii] * F[kk][jj];
  }
  }
  void NonLinFibMtx::LinearizedConstitutive(const apf::Matrix3x3 & C,
  double J,
  apf::Matrix<6,6> & strain)
  {
  // use shear modulus and poisson's ratio as material constant
  double beta = poisson_ratio / ( 1 - 2 * poisson_ratio);
  double c1 = shear_modulus * 0.5;
  double J_2beta = pow(J, -2 * beta);
  double coeff1 = 4 * beta * c1 * J_2beta;
  double coeff2 = 4 * c1 * J_2beta;
  apf::Matrix3x3 C_inv = apf::invert(C);
  for(int ii = 0; ii < 3; ii++)
  for(int jj = 0; jj < 3; jj++)
  for(int kk = 0; kk < 3; kk++)
  for(int ll = 0; ll < 3; ll++)
  {
  // use the voigt matrix to convert indices from 4th-order tensor to 6x6 matrix
  int x = Voigt[ii][jj];
  int y = Voigt[kk][ll];
  strain[x][y] = 0.0;
  strain[x][y] += coeff1 * C_inv[ii][jj] * C_inv[kk][ll];
  strain[x][y] += coeff2 * 0.5 *
  (C_inv[ii][kk] * C_inv[jj][ll] + C_inv[ii][ll] * C_inv[jj][kk]);
  }
  }
  void NonLinFibMtx::PK2Stress(const apf::Matrix3x3 & C,
  double J,
  apf::Matrix3x3 & pk2)
  {
  double beta = poisson_ratio / ( 1 - 2 * poisson_ratio );
  double gamma1 = shear_modulus;
  double gamma2 = -shear_modulus * pow(J, -2 * beta);
  apf::Matrix3x3 C_inv = apf::invert(C);
  for(int ii = 0; ii < 3; ii++)
  for(int jj = 0; jj < 3; jj++)
  pk2[ii][jj] = gamma1 * Identity[ii][jj] + gamma2 * C_inv[ii][jj];
  }
*/
}
