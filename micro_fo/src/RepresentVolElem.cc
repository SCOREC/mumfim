#include "RepresentVolElem.h"

#include "globals.h"
#include "RVE_Util.h"
#include "FiberNetwork.h"

// scorecutil stuff to phase out
#include "LagrangeMapping.h"

#include <amsiInterface.h>

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <gmi_null.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <iostream>
#include <valarray>

namespace Biotissue {


  MicroFO::MicroFO() :
    num_field_comp(3)
  { }
  
  MicroFO::MicroFO(int macro_type,
		   double * gauss,
		   int id,
		   FiberNetwork * fn,
		   SparseMatrix * st,
		   SparskitBuffers * b,
		   double * init_coords,
		   int num_nodes) :
    buffers(b),
    fiber_reactions(),
    firstTimeThrough(true),
    post_migration(false),
    element_type(macro_type),
    num_field_comp(3), // HACKY hard-coded
    num_element_nodes(num_nodes),
    network_id(id),
    num_rve_doubles(num_field_comp*num_element_nodes),
    init_fiber_network(fn),
    sparse_structure(st),
    num_gauss_pts(2) // HACKY hard-coded
  {
    fiber_G = 2500;
    fiber_J = 2.51327e-07;
    fiber_I1 = 1.25664e-07;
    fiber_I2 = 1.25664e-07;
    HALF_RVE_DIM = 0.011;
    fiber_B = 1.2;
    lexp = 1.5;

    /*
    else // beam
    {
      fiber_radius = 3.49911271e-05;
      fiber_E = 6500;
    }
    */

    gauss_pt[0] = gauss[0];
    gauss_pt[1] = gauss[1];
    gauss_pt[2] = gauss[2];

    fiber_network = fn->clone();

    initial_coords.assign(init_coords,init_coords+num_rve_doubles);
    // build a single-element mesh using the passed information about the corresponding element in macroscale
    /*
    apf::Vector3 verts[num_element_nodes];
    for(int ii = 0; ii < num_element_nodes; ii++)
      verts[ii] = apf::Vector3(&init_coords[ii*num_field_comp]);
    
    gmi_model * null_model = gmi_load(".null");
    int type_dim = apf::Mesh::typeDimension[macro_type];
    element_mesh = apf::makeEmptyMdsMesh(null_model,type_dim,false);
    apf::MeshEntity * e = apf::buildOneElement(element_mesh,
					       element_mesh->findModelEntity(type_dim,0),
					       macro_type,
					       verts);
    element_mesh->acceptChanges();
    // assumption about shape functions at macroscale
    apf::changeMeshShape(element_mesh,apf::getLagrange(1),true);

    // get the only mesh entity
    apf::MeshIterator * it = NULL;
    it = element_mesh->begin(type_dim);
    apf::MeshEntity * macro_entity = element_mesh->iterate(it);
    element_mesh->end(it);
    macro_element = apf::createMeshElement(element_mesh,
					 macro_entity);
    */
					 
    // todo: implement copy-constructor
//    fiber_network = new FiberNetwork(*fn);
    //fiber_network->copy(fn);

    
    SetRveCorners();
    
    int num_dofs = fiber_network->numDofs();
    force_vector.resize(num_dofs);
    force_vector_axial.resize(num_dofs);
    coordinate_vector.resize(num_dofs);
    tdydxr.resize(24*num_dofs);
    ttdSdy.resize(6*num_dofs);

    // fiber length using the dimensionless coordinate system of the RVE
    total_fiber_length = compute_total_fiber_length();
    /* volume fraction defined in Eq. (9) of T. Stylianopoulos and V.H. Barocas, Comput. Methods Appl. Mech. Engrg. 196 (2007) 2981-2990.*/
    fiber_volume_fraction = (total_fiber_length * fiber_network->fiberArea()) / ((2*HALF_RVE_DIM)*(2*HALF_RVE_DIM));
    scale_conversion = fiber_volume_fraction / (total_fiber_length * fiber_network->fiberArea());
  }
  
  MicroFO::~MicroFO()
  {
    //element_mesh->destroyNative();
    //apf::destroyMesh(element_mesh);
    
    // point to outside memory, so just point to null
    coords = NULL;
    displacement = NULL;
    rve_info = NULL;
  }

  void MicroFO::SetRveCorners()
  {
    rve[0][0] = fiber_network->sideCoord(FiberNetwork::LEFT);
    rve[0][1] = fiber_network->sideCoord(FiberNetwork::BOTTOM);
    rve[0][2] = fiber_network->sideCoord(FiberNetwork::FRONT);
  
    rve[1][0] = fiber_network->sideCoord(FiberNetwork::RIGHT);
    rve[1][1] = fiber_network->sideCoord(FiberNetwork::BOTTOM);
    rve[1][2] = fiber_network->sideCoord(FiberNetwork::FRONT);
  
    rve[2][0] = fiber_network->sideCoord(FiberNetwork::LEFT);
    rve[2][1] = fiber_network->sideCoord(FiberNetwork::BOTTOM);
    rve[2][2] = fiber_network->sideCoord(FiberNetwork::BACK);
  
    rve[3][0] = fiber_network->sideCoord(FiberNetwork::RIGHT);
    rve[3][1] = fiber_network->sideCoord(FiberNetwork::BOTTOM);
    rve[3][2] = fiber_network->sideCoord(FiberNetwork::BACK);
  
    rve[4][0] = fiber_network->sideCoord(FiberNetwork::LEFT);
    rve[4][1] = fiber_network->sideCoord(FiberNetwork::TOP);
    rve[4][2] = fiber_network->sideCoord(FiberNetwork::FRONT);
  
    rve[5][0] = fiber_network->sideCoord(FiberNetwork::RIGHT);
    rve[5][1] = fiber_network->sideCoord(FiberNetwork::TOP);
    rve[5][2] = fiber_network->sideCoord(FiberNetwork::FRONT);
  
    rve[6][0] = fiber_network->sideCoord(FiberNetwork::LEFT);
    rve[6][1] = fiber_network->sideCoord(FiberNetwork::TOP);
    rve[6][2] = fiber_network->sideCoord(FiberNetwork::BACK);
  
    rve[7][0] = fiber_network->sideCoord(FiberNetwork::RIGHT);
    rve[7][1] = fiber_network->sideCoord(FiberNetwork::TOP);
    rve[7][2] = fiber_network->sideCoord(FiberNetwork::BACK); 
  }
   
void MicroFO::SetDisplacement(double * input_coords)
{
  // Point coords and displacements to correct portion of input array
  coords = input_coords;
  displacement = &input_coords[num_rve_doubles]; // num_rve_doubles = size of input_coords / 2
}


  // these are just the shape functinos...
double MicroFO::PHI(int i, double u, double v, double w)
{
  double x[8], t, k, ui, vi;
  
  switch(element_type) {
  case apf::Mesh::HEX: // Hex;
      x[0] = 0.125 * (1.-u) * (1.-v) * (1.-w);
      x[1] = 0.125 * (1.+u) * (1.-v) * (1.-w);
      x[2] = 0.125 * (1.+u) * (1.+v) * (1.-w);
      x[3] = 0.125 * (1.-u) * (1.+v) * (1.-w);
      x[4] = 0.125 * (1.-u) * (1.-v) * (1.+w);
      x[5] = 0.125 * (1.+u) * (1.-v) * (1.+w);
      x[6] = 0.125 * (1.+u) * (1.+v) * (1.+w);
      x[7] = 0.125 * (1.-u) * (1.+v) * (1.+w);
    break;
  case apf::Mesh::PRISM: // Prism;
    k = 1.-u-v;
    x[0] =  (k) * (0.5*(1.-w));
    x[1] =  (u) * (0.5*(1.-w));
    x[2] =  (v) * (0.5*(1.-w));
    x[3] =  (k) * (0.5*(1.+w));
    x[4] =  (u) * (0.5*(1.+w));
    x[5] =  (v) * (0.5*(1.+w));
    break;
  case apf::Mesh::TET://Tet
    t = 1.-u-v-w;
    x[0] = t;
    x[1] = u;
    x[2] = v;
    x[3] = w;
    break;
  case apf::Mesh::PYRAMID://Pyramid
    ui = 2.*u/(1.-w); vi = 2.*v/(1.-w);      
    x[0] = 0.125*(1.-ui)*(1.-vi)*(1.-w);
    x[1] = 0.125*(1.+ui)*(1.-vi)*(1.-w);
    x[2] = 0.125*(1.+ui)*(1.+vi)*(1.-w);
    x[3] = 0.125*(1.-ui)*(1.+vi)*(1.-w);
    x[4] = 0.5*(1.+w);
    break;
  }
  return x[i];
}

  int MicroFO::run()
  {
    int result = 0;
    double fem_res_norm = 1.0;

    if(post_migration)
    {
      int num_elements = fiber_network->numElements();
      double fib_str[num_elements];
      double dfdE[num_elements];
      std::vector<double> len(num_elements);
      
      calcFiberLengths(*fiber_network,len);
      calc_fiber_str(&len[0],
		     &fib_str[0],
		     &dfdE[0]);
      
      // Construct "matrix" for following function
      calc_precond(&matrix[0],
		   &len[0],
		   &fib_str[0],
		   &dfdE[0],
		   false);
      
      // Reconstruct tdydx which is needed at the beginning
      // of the solve for the first order continuation guess
      calc_tdydxr();
      post_migration = false;
    }
    
    
    // wrap on the main_solver
    main_solver(coords,
		displacement,
		rve_info[0], rve_info[1], rve_info[2],
		rve_info[3], rve_info[4],
		rve_info[5], // average stress
		&rve_info[9],                          // derivative of stress
		rve_info[6], rve_info[7], rve_info[8], // volume balance stress
		fem_res_norm);
    
    return result; 
  }

void MicroFO::output(const std::string & filename)
{
  FILE * fp = fopen(filename.c_str(), "w");

  fprintf(fp, "%.15lf %.15lf %.15lf\n", rve_info[0], rve_info[1], rve_info[2]);
  fprintf(fp, "%.15lf %.15lf\n",rve_info[3], rve_info[4]);
  fprintf(fp, "%.15lf\n", rve_info[5]);
    
  fprintf(fp, "%.15lf %.15lf %.15lf\n", rve_info[6], rve_info[7], rve_info[8]);   
    
  // rotate the output order for debugging  
  for(int i = 0; i < 6; i++)
  {
    fprintf(fp, "%.15lf %.15lf %.15lf\n", rve_info[i*24+9],  rve_info[i*24+10], rve_info[i*24+11]);
    fprintf(fp, "%.15lf %.15lf %.15lf\n", rve_info[i*24+12], rve_info[i*24+13], rve_info[i*24+14]);
    fprintf(fp, "%.15lf %.15lf %.15lf\n", rve_info[i*24+15], rve_info[i*24+16], rve_info[i*24+17]);
    fprintf(fp, "%.15lf %.15lf %.15lf\n", rve_info[i*24+18], rve_info[i*24+19], rve_info[i*24+20]);
    fprintf(fp, "%.15lf %.15lf %.15lf\n", rve_info[i*24+21], rve_info[i*24+22], rve_info[i*24+23]);
    fprintf(fp, "%.15lf %.15lf %.15lf\n", rve_info[i*24+24], rve_info[i*24+25], rve_info[i*24+26]);
    fprintf(fp, "%.15lf %.15lf %.15lf\n", rve_info[i*24+27], rve_info[i*24+28], rve_info[i*24+29]);
    fprintf(fp, "%.15lf %.15lf %.15lf\n", rve_info[i*24+30], rve_info[i*24+31], rve_info[i*24+32]);
  }
   
  fclose(fp);  
}

void MicroFO::output(std::ofstream & outFile)
{
  //output the init coords first
  for(int i = 0; i < 8; i++)
    outFile << std::fixed << std::setprecision(16) << initial_coords[3*i] << " " << initial_coords[3*i+1] << " " << initial_coords[3*i+2] << std::endl;

  outFile << " init coords *******" << std::endl;

  for(int i = 0; i < 8; i++)
    outFile << std::fixed << std::setprecision(16) << coords[3*i] << " " << coords[3*i+1] << " " << coords[3*i+2] << std::endl;
  
  outFile << " coords *******" << std::endl;
  
  for(int i = 0; i < 8; i++)
    outFile << std::fixed << std::setprecision(16) << displacement[3*i] << " " << displacement[3*i+1] << " " << displacement[3*i+2] << std::endl;
  
  outFile << " disp *******" << std::endl;

  outFile << std::fixed << std::setprecision(16) << rve_info[0] << " " << rve_info[1] << " " << rve_info[2] << std::endl;

  outFile << std::fixed << std::setprecision(16) << rve_info[3] << " " << rve_info[4] << std::endl;
  
  outFile << std::fixed << std::setprecision(16) << rve_info[5] << std::endl;
  outFile << " stress *******" << std::endl;
  
  outFile << std::fixed << std::setprecision(16) << rve_info[6] << " " << rve_info[7] << " " << rve_info[8] << std::endl;

  outFile << " volstress *******" << std::endl;
    
  for(int i = 0; i < 6; i++)
  {
    outFile << std::fixed << std::setprecision(16) << rve_info[i*24+21] << " " << rve_info[i*24+22] << " " << rve_info[i*24+23] << std::endl;
    outFile << std::fixed << std::setprecision(16) << rve_info[i*24+9]  << " " << rve_info[i*24+10] << " " << rve_info[i*24+11] << std::endl;
    outFile << std::fixed << std::setprecision(16) << rve_info[i*24+27] << " " << rve_info[i*24+28] << " " << rve_info[i*24+29] << std::endl;
    outFile << std::fixed << std::setprecision(16) << rve_info[i*24+15] << " " << rve_info[i*24+16] << " " << rve_info[i*24+17] << std::endl;
    outFile << std::fixed << std::setprecision(16) << rve_info[i*24+24] << " " << rve_info[i*24+25] << " " << rve_info[i*24+26] << std::endl;
    outFile << std::fixed << std::setprecision(16) << rve_info[i*24+12] << " " << rve_info[i*24+13] << " " << rve_info[i*24+14] << std::endl;
    outFile << std::fixed << std::setprecision(16) << rve_info[i*24+30] << " " << rve_info[i*24+31] << " " << rve_info[i*24+32] << std::endl;
    outFile << std::fixed << std::setprecision(16) << rve_info[i*24+18] << " " << rve_info[i*24+19] << " " << rve_info[i*24+20] << std::endl;
  }
 
  outFile << " derivative stress *******" << std::endl;
}

void MicroFO::outputFiber(const std::string & filename) const
{
    FILE * fp = fopen(filename.c_str(), "w");

    int num_nodes = fiber_network->numNodes();
    int num_elements = fiber_network->numElements();
    
    fprintf(fp, "%d %d %d\n", num_nodes, 3*num_elements, num_elements);
    
    for(int ii = 0; ii < num_elements; ii++) 
    {
      int node1 = fiber_network->element(ii).node1_id;
      int node2 = fiber_network->element(ii).node2_id;
/*
      std::cout<<ii+1<<","<<node1<<","<<node2<<","
	       <<fiber_network->node(node1).x<<","<<fiber_network->node(node1).y<<","<<fiber_network->node(node1).z<<","
	       <<fiber_network->node(node2).x<<","<<fiber_network->node(node2).y<<","<<fiber_network->node(node2).z<<","
	       <<fiber_network->element(ii).fiber_type<<std::endl;
*/
	
      fprintf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf %i\n",
	      ii+1,
	      node1,
	      node2,
	      fiber_network->node(node1).x,
	      fiber_network->node(node1).y,
	      fiber_network->node(node1).z,
	      fiber_network->node(node2).x,
	      fiber_network->node(node2).y,
	      fiber_network->node(node2).z,
	      fiber_network->element(ii).fiber_type);
    }
    
    fclose(fp);
}
  
void MicroFO::printD(char * file, double * x, int num)
{
  FILE * fp = fopen(file, "w");
  fprintf(fp, "%d\n", num);

  for(int i = 0; i < num; i++)
    fprintf(fp, "%.15lf\n", x[i]);

  fclose(fp);
}

void MicroFO::printI(char * file, int * x, int num)
{
  FILE * fp = fopen(file, "w");
  fprintf(fp, "%d\n", num);

  for(int i = 0; i < num; i++)
    fprintf(fp, "%d\n", x[i]);

  fclose(fp);
}

void MicroFO::getRVEs(double * ss)
{
  memcpy(ss, rve_info, (num_rve_doubles * 6 + 9)*sizeof(double));
}

// Avoid memcpy in above function
void MicroFO::setResults(double * ss)
{
  rve_info = ss;
}
   
/*
void MicroFO::eval_stress(SCOREC::Util::mTensor2 & T)
{
  T[0][0] = rve_info[0]; T[0][1] = rve_info[1]; T[0][2] = rve_info[2];
  T[1][0] = rve_info[1]; T[1][1] = rve_info[3]; T[1][2] = rve_info[4];
  T[2][0] = rve_info[2]; T[2][1] = rve_info[4]; T[2][2] = rve_info[5];
}

void MicroFO::eval_balancestress(SCOREC::Util::mVector & V)
{
  V[0] = rve_info[6]; V[1] = rve_info[7]; V[2] = rve_info[8];
}
*/

void MicroFO::eval_derivstress(double * dSdx)
{
  memcpy(dSdx, rve_info + 9, 6 * num_rve_doubles * sizeof(double));
}

void MicroFO::SetDeformationGradient(int guass_pt, 
				     double * grad)
{
  for (int ii = 0; ii < 3; ii++) // ii = j
    for (int jj = 0; jj < 3; jj++) // jj = k
      F[ii][jj] = grad[guass_pt*9 + ii*3 + jj];
}

void MicroFO::SetDeformationGradients(double * grads)
{
  for(int ii = 0; ii < 8; ii++)  //loop over 8 vertices
    for(int jj = 0; jj < 3; jj++)
      for(int kk = 0; kk < 3; kk++)
	FItp[ii][jj][kk] = grads[ii*9 + jj*3 + kk];
}


void MicroFO::collectMigrationData()
{

  // Make sure lists are cleared
  clearMigrationData();

  // First int is how many ints there are (set at end)
  // Need this because ints and doubles are collected into
  // the same buffer to send
  intMigrationData.push_back(0);

  // Element type
  intMigrationData.push_back(element_type);
  intMigrationData.push_back(num_element_nodes);

  // Parameters
  doubleMigrationData.push_back(fiber_volume_fraction);
  doubleMigrationData.push_back(HALF_RVE_DIM);
  doubleMigrationData.push_back(fiber_E);
  doubleMigrationData.push_back(fiber_B);
  doubleMigrationData.push_back(lexp);
  doubleMigrationData.push_back(fiber_G);
  doubleMigrationData.push_back(fiber_J);
  doubleMigrationData.push_back(fiber_I1);
  doubleMigrationData.push_back(fiber_I2);

  // Gauss point
  doubleMigrationData.push_back(gauss_pt[0]);
  doubleMigrationData.push_back(gauss_pt[1]);
  doubleMigrationData.push_back(gauss_pt[2]);

  // Initial macro coordinates
  for(int ii=0;ii<num_rve_doubles;ii++)
    doubleMigrationData.push_back(initial_coords[ii]);

  // Network id for getting info from file
  intMigrationData.push_back(network_id);

  // Current node x,y,z coords
  fiber_network->getNodeCoordinates(doubleMigrationData);

  // RVE edges info
  for(int ii=0;ii<8;ii++)
  {
    doubleMigrationData.push_back(u[ii][0]);
    doubleMigrationData.push_back(u[ii][1]);
    doubleMigrationData.push_back(u[ii][2]);
  }

  // RVE edges info
  for(int ii=0;ii<8;ii++)
  {
    doubleMigrationData.push_back(rve[ii][0]);
    doubleMigrationData.push_back(rve[ii][1]);
    doubleMigrationData.push_back(rve[ii][2]);
  }

  // Set first int to number of ints
  intMigrationData[0] = intMigrationData.size();

}

void MicroFO::constructRVEFromMigrationData(FiberNetwork ** networks,
					    SparseMatrix ** sparse_structs,
					    SparskitBuffers * b)
{
  buffers = b;
  
  // For now just assume migration data has been placed in intMigrationData & doubleMigrationData
  int i_int  = 1; // Skip first integer
  int i_double = 0;

  // Use first order continuation not linear interpolation 
  firstTimeThrough = false;
  post_migration = true;

  // Get element type/gauss point and setup initial parameters
  element_type = intMigrationData[i_int++];
  num_field_comp = 3;
  num_gauss_pts = 2;
  num_element_nodes = intMigrationData[i_int++];
  num_rve_doubles = num_field_comp * num_element_nodes;

  fiber_volume_fraction = doubleMigrationData[i_double++];
  HALF_RVE_DIM = doubleMigrationData[i_double++];
  fiber_E      = doubleMigrationData[i_double++];
  fiber_B      = doubleMigrationData[i_double++];
  lexp         = doubleMigrationData[i_double++];
  fiber_G      = doubleMigrationData[i_double++];
  fiber_J      = doubleMigrationData[i_double++];
  fiber_I1     = doubleMigrationData[i_double++];
  fiber_I2     = doubleMigrationData[i_double++];

  gauss_pt[0] = doubleMigrationData[i_double++];
  gauss_pt[1] = doubleMigrationData[i_double++];
  gauss_pt[2] = doubleMigrationData[i_double++];

  // Get initial coordinates
  initial_coords.assign(&doubleMigrationData[i_double],
                        &doubleMigrationData[i_double]+num_rve_doubles);
  i_double += num_rve_doubles;

  // Reconstruct element info, initial nodes coords, and boundary nodes from networks data
  network_id = intMigrationData[i_int++];

  //fiber_network = new FiberNetwork(*networks[network_id]);
  fiber_network = networks[network_id]->clone();

  init_fiber_network = networks[network_id];

  sparse_structure = sparse_structs[network_id];

  SetRveCorners();

  int num_dofs = fiber_network->numDofs();
  force_vector.resize(num_dofs);
  coordinate_vector.resize(num_dofs);
  tdydxr.resize(24*num_dofs);
  ttdSdy.resize(6*num_dofs);

  total_fiber_length = compute_total_fiber_length();
  scale_conversion = fiber_volume_fraction / (total_fiber_length * fiber_network->fiberArea());

  // Set current node coord data
  fiber_network->setNodeCoordinates(&doubleMigrationData[i_double]);
  i_double += fiber_network->numNodes() * 3;

  // Get RVE boundary info
  for(int ii=0;ii<8;ii++)
  {
    u[ii][0] = doubleMigrationData[i_double++];
    u[ii][1] = doubleMigrationData[i_double++];
    u[ii][2] = doubleMigrationData[i_double++];
  }

  // Get RVE boundary info
  for(int ii=0;ii<8;ii++)
  {
    rve[ii][0] = doubleMigrationData[i_double++];
    rve[ii][1] = doubleMigrationData[i_double++];
    rve[ii][2] = doubleMigrationData[i_double++];
  }

  // Set up remaining data
  matrix.resize(sparse_structure->numNonzeros());
  matrix_axial.resize(sparse_structure->numNonzeros());

  clearMigrationData();
}

void MicroFO::clearMigrationData()
{
  // This just amounts to clearing the migration data vectors
  intMigrationData.clear();
  doubleMigrationData.clear();
}

// Fill a char vector with migration data
void MicroFO::getMigrationData(std::vector<char> & rveData)
{
  int rveDataSize = doubleMigrationData.size()*sizeof(double)
                    + intMigrationData.size()*sizeof(int);
  rveData.resize(rveDataSize);

  memcpy(rveData.data(),intMigrationData.data(),intMigrationData.size()*sizeof(int));
  memcpy(rveData.data()+intMigrationData.size()*sizeof(int),
         doubleMigrationData.data(),
         doubleMigrationData.size()*sizeof(double));

}

void MicroFO::setMigrationData(std::vector<char> & rveData)
{
  int numInts = ((int*)rveData.data())[0];
  int numDoubles = (rveData.size()-numInts*sizeof(int))/sizeof(double);

  intMigrationData.resize(numInts);
  doubleMigrationData.resize(numDoubles);

  memcpy(intMigrationData.data(),rveData.data(),numInts*sizeof(int));
  memcpy(doubleMigrationData.data(),rveData.data()+numInts*sizeof(int),numDoubles*sizeof(double));
}

  // Return rve_iterations per load step
double MicroFO::getWeight()
{
  double weight = 0.0;
  for(std::vector<double>::iterator iter = rve_iterations.begin(); iter != rve_iterations.end(); ++iter)
    weight += *iter;
  return weight;
}
  

  // Return rve_timing per load step
double MicroFO::getTiming()
{
  double timing = 0.0;
  for(std::vector<double>::iterator iter = rve_timing.begin(); iter != rve_timing.end(); ++iter)
    timing += *iter;
  return timing;
}

 }
