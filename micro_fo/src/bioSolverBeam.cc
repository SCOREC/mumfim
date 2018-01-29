namespace bio
{
  /*
  // Solver for beam network
  // Uses corotational formulation, K_e = T_e^T k_e T_e,
  // where k_e is the stiffness matrix for an element e
  // evaluated in the reference configuration, and T_e is
  // the transformation matrix, taking k_e to global coordinates.
  // This formulation is valid for small strain, but large
  // displacement and rotation.
  / *
  int MicroFO::Solver_Beam()
  {
    int result = 0;
    int num_dofs = fiber_network->numDofs();
    //int num_elements = fiber_network->numElements();
    //int num_nodes = fiber_network->numNodes();
    int step = 0;
    //int cont;
    double tol = 0.0;
    double tolerance = 1e-8; // Convergence seems pretty tough for this problem
    double norm = 1.0;
    int ierr = 0;
    double solution[num_dofs];
    //double fib_str[num_elements];
    //double dfdE[num_elements];
    std::vector<double> len;
    while(norm > tolerance)
    {
      calcFiberLengths(*fiber_network,len);
      // Calculate initial displacements
      calc_current_disps_beam(&solution[0]);
      // Get tangent stiffness matrix
      calc_matrix_beam(&matrix[0],&len[0]);
      //std::vector<double> fm;
      //sparse_structure->debugFullMatrix(&matrix[0],fm,num_dofs);
      // Calculate force vector
      calc_force_vector_beam(&solution[0]);
      // Apply boundary conditions
      matrix_bcs_beam();
      force_vector_bcs_beam();
      //std::vector<double> fmbc;
      //sparse_structure->debugFullMatrix(&matrix[0],fmbc,num_dofs);
      buffers->zero();
      // Preconditioner
      int buffer_length = buffers->matrixLength();
      ilut_(&num_dofs,
            &matrix[0],
            sparse_structure->getCols(),
            sparse_structure->getRows(),
            &num_dofs,
            &tol,
            buffers->matrixBuffer(),
            buffers->colsBuffer(),
            buffers->rowsBuffer(),
            &buffer_length,
            buffers->doubleWorkBuffer(),
            buffers->intWorkBuffer(),
            &ierr);
      if(ierr != 0)
      {
        std::cerr << "ERROR: ilut_ returned error code " << ierr << std::endl;
        result++;
      }
      // Solver
      lusol_(&num_dofs,
             &force_vector[0],
             &solution[0],
             buffers->matrixBuffer(),
             buffers->colsBuffer(),
             buffers->rowsBuffer());
      norm = calc_norm_beam(&solution[0]);
      for(int ii = 0; ii < num_dofs; ii++){
        coordinate_vector[ii] += solution[ii];
      }
      update_nodes();
      // Getting this from truss solver
      // Limit iterations by easing tolerance condition
      step++;
      if(step == 20)
        tolerance *= 10;
      else if(step == 30)
        tolerance *= 10;
      else if(step == 40)
        tolerance *= 10;
      else if(step == 50)
        tolerance *= 10;
      else if(step > 60)
      {
        std::cerr << "Warning: unusual number of newton iterations in micro_fo!" << std::endl;
        norm = tolerance / 2;
      }
    }
    // After final iteration get forces once more without BCs applied
    calcFiberLengths(*fiber_network,len);
    calc_matrix_beam(&matrix[0],&matrix_axial[0],&len[0]);
    calc_current_disps_beam(&solution[0]);
    calc_force_vector_beam(&solution[0]);
    calc_dsdy_beam();
    // Newton iteration ended. The network got a new equilibrium position and the fiber length and fiber force are updated
    //calc_fiber_str(&len[0],
//                 &fib_str[0],
//                 &dfdE[0]);
//    calc_force_vector(&len[0],
//                    &fib_str[0]);
    //calc_mean_fiber_stretch_omega();
    //update the values of arrnode, arrelmt with the new equilibrium position of the micronodes
    // is this call necessary ?
    //update_nodes();
    return result;
  }
  */
  /*
void MicroFO::calc_current_disps_beam(double * solution)
{
  int num_nodes = fiber_network->numNodes();
  for(int ii=0;ii<num_nodes;ii++)
  {
    Node n = fiber_network->node(ii);
    Node in = init_fiber_network->node(ii);
    solution[ii*6]   = n.x - in.x;
    solution[ii*6+1] = n.y - in.y;
    solution[ii*6+2] = n.z - in.z;
    solution[ii*6+3] = n.rx - in.rx;
    solution[ii*6+4] = n.ry - in.ry;
    solution[ii*6+5] = n.rz - in.rz;
  }
}
*/
/*
double MicroFO::calc_norm_beam(double * x)
{
  double total = 0.0;
  for(int ii=0;ii<fiber_network->numNodes();ii++)
  {
    const Node & n = fiber_network->node(ii);
    total += pow( x[ii*6  ] - (n.x  - coordinate_vector[ii*6    ]) ,2.0) +
             pow( x[ii*6+1] - (n.y  - coordinate_vector[ii*6 + 1]) ,2.0) +
             pow( x[ii*6+2] - (n.z  - coordinate_vector[ii*6 + 2]) ,2.0) +
             pow( x[ii*6+3] - (n.rx - coordinate_vector[ii*6 + 3]) ,2.0) +
             pow( x[ii*6+4] - (n.ry - coordinate_vector[ii*6 + 4]) ,2.0) +
             pow( x[ii*6+5] - (n.rz - coordinate_vector[ii*6 + 5]) ,2.0) ;
  }
  return sqrt(total);
}
*/
 /*
void MicroFO::force_vector_bcs_beam()
{
  // Boundary conditions on force vector for the microscopic problem
  if(PERIODIC)
  {
    // For each periodic pair add the equations together and set
    // the other equation to equality
    // Except for displacement corresponding to RVE side, which is
    // zero displacement dirichlet
    const std::vector<PBCRelation> & rl_bcs = fiber_network->getRightLeftBCs();
    const std::vector<PBCRelation> & tb_bcs = fiber_network->getTopBottomBCs();
    const std::vector<PBCRelation> & fb_bcs = fiber_network->getFrontBackBCs();
    // Right / left
    int pdofs_rl[5] = {1,2,3,4,5};
    force_vector_periodic_bcs_helper(rl_bcs,0,pdofs_rl);
    // Top / bottom
    int pdofs_tb[5] = {0,2,3,4,5};
    force_vector_periodic_bcs_helper(tb_bcs,1,pdofs_tb);
    // Front / back
    int pdofs_fb[5] = {0,1,3,4,5};
    force_vector_periodic_bcs_helper(fb_bcs,2,pdofs_fb);
  }
  else
  {
    // Dirichlet - all boundary dofs have zero displacement
    int num_boundary_nodes = fiber_network->numBoundaryNodes(FiberNetwork::ALL);
    for(int ii = 0; ii < num_boundary_nodes; ii++)
    {
      int node = fiber_network->boundaryNode(FiberNetwork::ALL,ii);
      force_vector[node*6] = 0.0;
      force_vector[node*6+1] = 0.0;
      force_vector[node*6+2] = 0.0;
      force_vector[node*6+3] = 0.0;
      force_vector[node*6+4] = 0.0;
      force_vector[node*6+5] = 0.0;
    }
  }
}
*/
/*
void MicroFO::force_vector_periodic_bcs_helper(const std::vector<PBCRelation> & bcs, int ddof, int * pdofs)
{
  int num_bcs = bcs.size();
  for(int ii = 0; ii < num_bcs; ii++)
  {
    // Dirichlet
    int node1 = bcs[ii].node1_id;
    int node2 = bcs[ii].node2_id;
    force_vector[node1*6 + ddof] = 0.0;
    force_vector[node2*6 + ddof] = 0.0;
    // Periodic
    for(int dd=0;dd<5;dd++)
    {
      int jj = node1*6 + pdofs[dd];
      int kk = node2*6 + pdofs[dd];
      force_vector[kk] += force_vector[jj];
      force_vector[jj] = 0.0;
    }
  }
}
*/

 /*
void MicroFO::calc_matrix_beam(double * matrix,double * lengths)
{
  // Construct global stiffness matrix for beam frame
  int num_elements = fiber_network->numElements();
  memset(matrix,0,sparse_structure->numNonzeros()*sizeof(double));
  double t[12][12];
  double material_matrix[144];
  double temp_matrix[144];
  memset(t,0,144*sizeof(double));
  double cxx,cyx,czx,cxy,cyy,czy,cxz,cyz,czz;
  for(int ii = 0; ii < num_elements; ii++)
  {
    const Element & e = fiber_network->element(ii);
    const Node & node1 = fiber_network->node(e.node1_id);
    const Node & node2 = fiber_network->node(e.node2_id);
    double L = lengths[ii];
    // Calculate transformation matrix
    cxx = (node2.x - node1.x) / L;
    cyx = (node2.y - node1.y) / L;
    czx = (node2.z - node1.z) / L;
    if( close(cxx,0.0) && close(cyx,0.0) )
    {
      if(czx>0)
      {
        cxx =  0.0; cyx = 0.0; czx = 1.0;
        cxy =  0.0; cyy = 1.0; czy = 0.0;
        cxz = -1.0; cyz = 0.0; czz = 0.0;
      }
      else
      {
        cxx = 0.0; cyx = 0.0; czx = -1.0;
        cxy = 0.0; cyy = 1.0; czy = 0.0;
        cxz = 1.0; cyz = 0.0; czz = 0.0;
      }
    }
    else
    {
        double D = sqrt(cxx*cxx + cyx*cyx);
        cxy = -cyx/D;
        cyy = cxx/D;
        czy = 0.0;
        cxz = -cxx*czx/D;
        cyz = -cyx*czx/D;
        czz = D;
    }
    for(int jj=0;jj<12;jj+=3)
    {
      t[jj][jj] = cxx;
      t[jj][jj+1] = cyx;
      t[jj][jj+2] = czx;
      t[jj+1][jj] = cxy;
      t[jj+1][jj+1] = cyy;
      t[jj+1][jj+2] = czy;
      t[jj+2][jj] = cxz;
      t[jj+2][jj+1] = cyz;
      t[jj+2][jj+2] = czz;
    }
    // Get material stiffness matrix
    make_material_stiffness_matrix(&material_matrix[0],ii,L);
    // Get geometric stiffness matrix
    if(GEOMETRIC_STIFFNESS)
    {
      double geometry_matrix[144];
      make_geometry_stiffness_matrix(&geometry_matrix[0],ii,L);
      for(int ii=0;ii<144;ii++)
        material_matrix[ii] += geometry_matrix[ii];
    }
    // Transform stiffness matrix from local to global coordinates
    matrix_multiply_ATranspose(&t[0][0],12,12,&material_matrix[0],12,12,temp_matrix);
    matrix_multiply(temp_matrix,12,12,&t[0][0],12,12,&material_matrix[0]);
    // Fill global stiffness matrix
    double * material_matrix_access = &material_matrix[0];
    for(int jj=0;jj<144;jj++)
      matrix[e.Ke[jj]] -= material_matrix_access[jj]; //calculating -K (negative of tangent stiffness matrix).
  }
}
*/
  /*
// overloaded function
  void MicroFO::calc_matrix_beam(double * matrix,double * matrix_axial,double * lengths)
{
  // Construct global stiffness matrix for beam frame
  int num_elements = fiber_network->numElements();
  memset(matrix,0,sparse_structure->numNonzeros()*sizeof(double));
  memset(matrix_axial,0,sparse_structure->numNonzeros()*sizeof(double));
  double t[12][12];
  double material_matrix[144];
  double material_matrix_axial[144];
  double temp_matrix[144];
  double temp_matrix_axial[144];
  memset(t,0,144*sizeof(double));
  double cxx,cyx,czx,cxy,cyy,czy,cxz,cyz,czz;
  for(int ii = 0; ii < num_elements; ii++)
  {
    const Element & e = fiber_network->element(ii);
    const Node & node1 = fiber_network->node(e.node1_id);
    const Node & node2 = fiber_network->node(e.node2_id);
    double L = lengths[ii];
    // Calculate transformation matrix
    cxx = (node2.x - node1.x) / L;
    cyx = (node2.y - node1.y) / L;
    czx = (node2.z - node1.z) / L;
    if( close(cxx,0.0) && close(cyx,0.0) )
    {
      if(czx>0)
      {
        cxx =  0.0; cyx = 0.0; czx = 1.0;
        cxy =  0.0; cyy = 1.0; czy = 0.0;
        cxz = -1.0; cyz = 0.0; czz = 0.0;
      }
      else
      {
        cxx = 0.0; cyx = 0.0; czx = -1.0;
        cxy = 0.0; cyy = 1.0; czy = 0.0;
        cxz = 1.0; cyz = 0.0; czz = 0.0;
      }
    }
    else
    {
        double D = sqrt(cxx*cxx + cyx*cyx);
        cxy = -cyx/D;
        cyy = cxx/D;
        czy = 0.0;
        cxz = -cxx*czx/D;
        cyz = -cyx*czx/D;
        czz = D;
    }
    for(int jj=0;jj<12;jj+=3)
    {
      t[jj][jj] = cxx;
      t[jj][jj+1] = cyx;
      t[jj][jj+2] = czx;
      t[jj+1][jj] = cxy;
      t[jj+1][jj+1] = cyy;
      t[jj+1][jj+2] = czy;
      t[jj+2][jj] = cxz;
      t[jj+2][jj+1] = cyz;
      t[jj+2][jj+2] = czz;
    }
    // Get axial material stiffness matrix
    memset(material_matrix_axial,0,144*sizeof(double));
    double dL = 1.0/L;
    double E = 0; //E(ii);
    double A = fiber_area;
    // Create local stiffness matrix (only diagonal and above)
    // Its symmetric
    // Axial
    material_matrix_axial[0*12 + 0] =  E*A*dL;
    material_matrix_axial[0*12 + 6] = -E*A*dL;
    material_matrix_axial[6*12 + 6] =  E*A*dL;
    // Flip over diagonal
    for(int jj=0;jj<12;jj++)
      for(int kk=0;kk<jj;kk++)
        material_matrix_axial[jj*12 + kk] = material_matrix_axial[kk*12 + jj];
    // Get material stiffness matrix
    make_material_stiffness_matrix(&material_matrix[0],ii,L);
    // Get geometric stiffness matrix
    if(GEOMETRIC_STIFFNESS)
    {
      double geometry_matrix[144];
      make_geometry_stiffness_matrix(&geometry_matrix[0],ii,L);
      for(int ii=0;ii<144;ii++)
        material_matrix[ii] += geometry_matrix[ii];
    }
    // Transform stiffness matrix from local to global coordinates
    matrix_multiply_ATranspose(&t[0][0],12,12,&material_matrix_axial[0],12,12,temp_matrix_axial);
    matrix_multiply(temp_matrix_axial,12,12,&t[0][0],12,12,&material_matrix_axial[0]);
    matrix_multiply_ATranspose(&t[0][0],12,12,&material_matrix[0],12,12,temp_matrix);
    matrix_multiply(temp_matrix,12,12,&t[0][0],12,12,&material_matrix[0]);
    // Fill global stiffness matrix
    double * material_matrix_access = &material_matrix[0];
    double * material_matrix_axial_access = &material_matrix_axial[0];
    for(int jj=0;jj<144;jj++)
    {
      matrix_axial[e.Ke[jj]] -= material_matrix_axial_access[jj]; // calculating -K (negative of tangent stiffness matrix)
      matrix[e.Ke[jj]] -= material_matrix_access[jj]; // calculating -K (negative of tangent stiffness matrix)
    }
  }
}
*/
/*
void MicroFO::make_material_stiffness_matrix(double * matrix, int ielem, double L)
{
    memset(matrix,0,144*sizeof(double));
    double dL = 1.0/L;
    double dL2 = 1.0/(L*L);
    double dL3 = 1.0/(L*L*L);
    double E = 0; //fiber_get_E(ielem);
    double A = fiber_area; //fiber_network->fiberArea();
    double G = 0; //fiber_get_G(ielem);
    double J = 0; //fiber_get_J(ielem);
    double I1 = 0; //fiber_get_I1(ielem);
    double I2 = 0; //fiber_get_I2(ielem);
    // Create local stiffness matrix (only diagonal and above)
    // Its symmetric
    // Axial
    matrix[0*12 + 0] =  E*A*dL;
    matrix[0*12 + 6] = -E*A*dL;
    matrix[6*12 + 6] =  E*A*dL;
    // Torsion
    matrix[3*12 + 3] =  G*J*dL;
    matrix[3*12 + 9] = -G*J*dL;
    matrix[9*12 + 9] =  G*J*dL;
    // Bending
    matrix[1*12 + 1]   =  12*E*I2*dL3;
    matrix[2*12 + 2]   =  12*E*I1*dL3;
    matrix[4*12 + 4]   =   4*E*I1*dL;
    matrix[5*12 + 5]   =   4*E*I2*dL;
    matrix[7*12 + 7]   =  12*E*I2*dL3;
    matrix[8*12 + 8]   =  12*E*I1*dL3;
    matrix[10*12 + 10] =   4*E*I1*dL;
    matrix[11*12 + 11] =   4*E*I2*dL;
    matrix[1*12 + 5]   =   6*E*I2*dL2;
    matrix[2*12 + 4]   =  -6*E*I1*dL2;
    matrix[4*12 + 8]   =   6*E*I1*dL2;
    matrix[5*12 + 7]   =  -6*E*I2*dL2;
    matrix[7*12 + 11]  =  -6*E*I2*dL2;
    matrix[8*12 + 10]  =   6*E*I1*dL2;
    matrix[1*12 + 7]   = -12*E*I2*dL3;
    matrix[2*12 + 8]   = -12*E*I1*dL3;
    matrix[4*12 + 10]  =   2*E*I1*dL;
    matrix[5*12 + 11]  =   2*E*I2*dL;
    matrix[1*12 + 11]  =   6*E*I2*dL2;
    matrix[2*12 + 10]  =  -6*E*I1*dL2;
    // Add shearing
    if(TIMOSHENKO)
    {
    }
    // Flip over diagonal
    for(int jj=0;jj<12;jj++)
      for(int kk=0;kk<jj;kk++)
        matrix[jj*12 + kk] = matrix[kk*12 + jj];
}
*/
 /*
void MicroFO::make_geometry_stiffness_matrix(double * matrix, int ielem, double L)
{
    memset(matrix,0,144*sizeof(double));
    / *
    double dL = 1.0/L;
    double dL2 = 1.0/(L*L);
    double dL3 = 1.0/(L*L*L);
    double E = 0; //fiber_get_E(ielem);
    double A = fiber_area; //fiber_network->fiberArea();
    double G = 0; //fiber_G;
    double J = 0;  //fiber_J;
    double I1 = 0; //fiber_I1;
    double I2 = 0; //fiber_I2;
    * /
    // Add shearing
    if(TIMOSHENKO)
    {
    }
    // Flip over diagonal
    for(int jj=0;jj<12;jj++)
      for(int kk=0;kk<jj;kk++)
        matrix[jj*12 + kk] = matrix[kk*12 + jj];
}
 */
/*
void MicroFO::calc_force_vector_beam(double * x)
{
  int num_dofs = fiber_network->numDofs();
  amux_(&num_dofs,
        x,
<<<<<<< HEAD:micro_fo/src/bioSolverBeam.cc
        &force_vector[0] / *force_vector.data()* /,
=======
        &force_vector[0] //force_vector.data(),
>>>>>>> 584b07b234754b7455adc47f885b19e4b0eaadff:micro_fo/src/bioSolverBeam.cc
        &matrix[0],
        sparse_structure->getCols(),
        sparse_structure->getRows());
  amux_(&num_dofs,
        x,
        &force_vector_axial[0],
        &matrix_axial[0],
        sparse_structure->getCols(),
        sparse_structure->getRows());
  for (int ii = 0; ii < num_dofs; ii++){
    force_vector[ii] = -1.0*force_vector[ii];
    force_vector_axial[ii] = -1.0*force_vector_axial[ii];
  } // to get correct sign for force vector because matrix and matrix_axial are -K.
}
*/
 /*
void MicroFO::matrix_bcs_beam()
{
  // Boundary conditions on stiffness matrix for the microscopic problem
  if(PERIODIC)
  {
    // For each periodic pair add the equations together and set
    // the other equation to equality
    // Except for displacement corresponding to RVE side, which is
    // zero displacement dirichlet
    // Note:
    // Assuming that we won't have multiple elements connecting to the
    // same boundary node. This is what the 'con' vector checks for
    // in the old dirichlet code below. It almost certainly won't
    // happen for the networks we're creating, unless the network is
    // purposefully constructed like that, at which point this would
    // need to be updated.
    const std::vector<PBCRelation> & rl_bcs = fiber_network->getRightLeftBCs();
    const std::vector<PBCRelation> & tb_bcs = fiber_network->getTopBottomBCs();
    const std::vector<PBCRelation> & fb_bcs = fiber_network->getFrontBackBCs();
    // Right / left
    int pdofs_rl[5] = {1,2,3,4,5};
    matrix_periodic_bcs_helper(rl_bcs,0,pdofs_rl);
    // Top / bottom
    int pdofs_tb[5] = {0,2,3,4,5};
    matrix_periodic_bcs_helper(tb_bcs,1,pdofs_tb);
    // Front / back
    int pdofs_fb[5] = {0,1,3,4,5};
    matrix_periodic_bcs_helper(fb_bcs,2,pdofs_fb);
  }
  else
  {
    // Dirichlet - all boundary dofs have zero displacement
    int num_elements = fiber_network->numElements();
    std::vector<int> con;
    int num_boundary_nodes = fiber_network->numBoundaryNodes(FiberNetwork::ALL);
    for(int ii = 0; ii < num_boundary_nodes; ii++)
    {
      int node = fiber_network->boundaryNode(FiberNetwork::ALL,ii);
      con.push_back(node);
      for(int jj = 0; jj < num_elements; jj++)
      {
        const Element & e = fiber_network->element(jj);
        if(e.node1_id == node)
          con.push_back(e.node2_id);
        if(e.node2_id == node)
          con.push_back(e.node1_id);
      }
      for(std::vector<int>::iterator iter = con.begin(), iterend = con.end();
          iter != iterend; iter++)
      {
        for(int nn = 0; nn < 6; nn++)
        {
          int jj = 6 * (*iter) + nn;
          matrix[sparse_structure->sparseLocation(jj,node*6  )] = 0.0;
          matrix[sparse_structure->sparseLocation(jj,node*6+1)] = 0.0;
          matrix[sparse_structure->sparseLocation(jj,node*6+2)] = 0.0;
          matrix[sparse_structure->sparseLocation(jj,node*6+3)] = 0.0;
          matrix[sparse_structure->sparseLocation(jj,node*6+4)] = 0.0;
          matrix[sparse_structure->sparseLocation(jj,node*6+5)] = 0.0;
        }
      }
      matrix[sparse_structure->sparseLocation(node*6  ,node*6  )] = 1.0;
      matrix[sparse_structure->sparseLocation(node*6+1,node*6+1)] = 1.0;
      matrix[sparse_structure->sparseLocation(node*6+2,node*6+2)] = 1.0;
      matrix[sparse_structure->sparseLocation(node*6+3,node*6+3)] = 1.0;
      matrix[sparse_structure->sparseLocation(node*6+4,node*6+4)] = 1.0;
      matrix[sparse_structure->sparseLocation(node*6+5,node*6+5)] = 1.0;
      con.erase(con.begin(),con.end());
    }
  }
}
 */
/*
void MicroFO::matrix_periodic_bcs_helper(const std::vector<PBCRelation> & bcs, int ddof, int * pdofs)
{
  // This is a bit of a mess, but trying to avoid using sparseLocation()
  int num_dofs_per_node = fiber_network->dofsPerNode();
  int num_dofs_per_element = 2 * num_dofs_per_node;
  int num_rl_bcs = bcs.size();
  for(int ii = 0; ii < num_rl_bcs; ii++)
  {
    //int pnode1 = bcs[ii].node1_id;
    //int pnode2 = bcs[ii].node2_id;
    int elem1 = bcs[ii].elem1;
    int elem2 = bcs[ii].elem2;
    int offset1,offset2;
    // Set offsets on whether boundary nodes are first or second node in element
    if(bcs[ii].elem1first)
      offset1 = 0;
    else
      offset1 = num_dofs_per_node*num_dofs_per_element;
    if(bcs[ii].elem2first)
      offset2 = 0;
    else
      offset2 = num_dofs_per_node*num_dofs_per_element;
    // Periodic - add equations from pnode1 to equations of pnode2
    for(int jj = 0; jj < num_dofs_per_node-1; jj++)
      for(int kk = 0; kk < num_dofs_per_element; kk++)
      {
        int locationTo = bcs[ii].Ke1[pdofs[jj]*num_dofs_per_element + kk];
        const Element & e = fiber_network->element(elem1);
        int locationFrom = e.Ke[offset1 + pdofs[jj]*num_dofs_per_element + kk];
        matrix[locationTo] += matrix[locationFrom];
      }
    // Periodic - set equations for pnode1 to equality conditions
    // First zero out what was just moved to pnode2
    for(int jj = 0; jj < num_dofs_per_node; jj++)
      for(int kk = 0; kk < num_dofs_per_element; kk++)
      {
        const Element & e = fiber_network->element(elem1);
        int location = e.Ke[offset1 + jj*num_dofs_per_element + kk];
        matrix[location] = 0.0;
      }
    // Set equality conditions
    for(int jj = 0; jj < num_dofs_per_node-1; jj++)
    {
      matrix[bcs[ii].Ke2[2*pdofs[jj]    ]] =  1.0;
      matrix[bcs[ii].Ke2[2*pdofs[jj] + 1]] = -1.0;
    }
    const Element & e = fiber_network->element(elem2);
    // Dirichlet - pnode2
    for(int jj = 0; jj < num_dofs_per_element; jj++)
      matrix[e.Ke[offset2 + ddof*num_dofs_per_element + jj]] = 0.0;
    int KeIndex = offset2 + ddof*num_dofs_per_element + ddof + num_dofs_per_node*(!bcs[ii].elem2first);
    matrix[e.Ke[KeIndex]] = 1.0;
    // Dirichlet - pnode1 (zeroing out should already be done)
    matrix[bcs[ii].Ke2[2*ddof]] = 1.0;
  }
}
<<<<<<< HEAD:micro_fo/src/bioSolverBeam.cc
// Calculates dSdy from beam values
=======
>>>>>>> 584b07b234754b7455adc47f885b19e4b0eaadff:micro_fo/src/bioSolverBeam.cc
*/
/*
void MicroFO::calc_dsdy_beam()
{
  int num_dofs = fiber_network->numDofs();
  int num_elements = fiber_network->numElements();
  ttdSdy.assign(num_dofs*6,0);
  std::vector<int> con;
  // Calculation of the derivative dSdy which is used later for the calculation
  // of the derivative dSdx in the function calc_femjacob_newmethod
  int num_boundary_nodes = fiber_network->numBoundaryNodes(FiberNetwork::ALL);
  for(int ii = 0; ii < num_boundary_nodes; ii++)
  {
    int node = fiber_network->boundaryNode(FiberNetwork::ALL,ii);
    Node n = fiber_network->node(node);
    double xyz[3];
    xyz[0] = n.x;
    xyz[1] = n.y;
    xyz[2] = n.z;
    con.push_back(node);
    for(int jj = 0; jj < num_elements; jj++)
    {
      const Element & e = fiber_network->element(jj);
      if(e.node1_id == node)
        con.push_back(e.node2_id);
      if(e.node2_id == node)
        con.push_back(e.node1_id);
    }
    for(std::vector<int>::iterator iter = con.begin(), iterend = con.end();
      iter != iterend; iter++)
    {
        for(int nn = 0; nn < 6; nn++)
        {
          int kk = 6 * (*iter) + nn;
          int location = sparse_structure->sparseLocation(kk,node * 6);
          ttdSdy[               kk] -= matrix[location] * xyz[0];
          ttdSdy[    num_dofs + kk] -= matrix[location] * xyz[1] / 2;
          ttdSdy[2 * num_dofs + kk] -= matrix[location] * xyz[2] / 2;
          if(kk == node * 6)
            ttdSdy[               kk] += force_vector[node * 6];
          if(kk == node * 6 + 1)
            ttdSdy[    num_dofs + kk] += force_vector[node * 6] / 2;
          if(kk == node * 6 + 2)
            ttdSdy[2 * num_dofs + kk] += force_vector[node * 6] / 2;
          matrix[location] = 0.0;
          location = sparse_structure->sparseLocation(kk,node * 6 + 1);
          ttdSdy[    num_dofs + kk] -= matrix[location] * xyz[0] / 2;
          ttdSdy[3 * num_dofs + kk] -= matrix[location] * xyz[1];
          ttdSdy[4 * num_dofs + kk] -= matrix[location] * xyz[2] / 2;
          if(kk == node * 6)
            ttdSdy[    num_dofs + kk] += force_vector[node * 6 + 1] / 2;
          if(kk == node * 6 + 1)
            ttdSdy[3 * num_dofs + kk] += force_vector[node * 6 + 1];
          if(kk == node * 6 + 2)
            ttdSdy[4 * num_dofs + kk] += force_vector[node * 6 + 1] / 2;
          matrix[location] = 0.0;
          location = sparse_structure->sparseLocation(kk,node * 6 + 2);
          ttdSdy[2 * num_dofs + kk] -= matrix[location] * xyz[0] / 2;
          ttdSdy[4 * num_dofs + kk] -= matrix[location] * xyz[1] / 2;
          ttdSdy[5 * num_dofs + kk] -= matrix[location] * xyz[2];
          if(kk == node * 6)
            ttdSdy[2 * num_dofs + kk] += force_vector[node * 6 + 2] / 2;
          if(kk == node * 6 + 1)
            ttdSdy[4 * num_dofs + kk] += force_vector[node * 6 + 2] / 2;
          if(kk == node * 6 + 2)
            ttdSdy[5 * num_dofs + kk] += force_vector[node * 6 + 2];
          matrix[location] = 0.0;
        }
      }
      matrix[sparse_structure->sparseLocation(node * 6    ,node * 6    )] = 1.0;
      matrix[sparse_structure->sparseLocation(node * 6 + 1,node * 6 + 1)] = 1.0;
      matrix[sparse_structure->sparseLocation(node * 6 + 2,node * 6 + 2)] = 1.0;
      con.erase(con.begin(),con.end());
    }
}
*/
}
