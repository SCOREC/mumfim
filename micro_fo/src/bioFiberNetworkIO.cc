#include "bioFiberNetwork.h"
#include <amsiMultiscale.h>
#include <amsiUtil.h>
#include <cassert>
#include <iostream>
namespace bio
{
  void FiberNetwork::readFromFile(const std::string & filename)
  {
    AMSI_DEBUG(std::cout << "Attempting to access file: " << filename.c_str() << " for reading." << std::endl);
    FILE * file = fopen(filename.c_str(),"r");
    assert(file);
    // read header info
    fseek(file,0,SEEK_SET);
    // Assume no periodic connections
    int num_new_rl_bcs = 0;
    int num_new_tb_bcs = 0;
    int num_new_fb_bcs = 0;
    int num_new_nodes = 0;
    int num_new_dofs = 0;
    int num_new_elements = 0;
    // Periodic has extra info in header (total numbers of periodic pairs)
    fscanf(file,"%d %d %d\n",
           &num_new_nodes,
           &num_new_dofs,
           &num_new_elements);
    num_dofs += dofs_per_node * num_new_nodes;
    elements.resize(num_elements+num_new_elements);
    nodes.resize(num_nodes+num_new_nodes);
    rl_bcs.resize(num_rl_bcs+num_new_rl_bcs);
    tb_bcs.resize(num_tb_bcs+num_new_tb_bcs);
    fb_bcs.resize(num_fb_bcs+num_new_fb_bcs);
    // Fill element array
    for(int ii = num_elements; ii < num_elements+num_new_elements; ii++)
    {
      double node1[3], node2[3];
      int temp;
      fscanf(file,"%d %d %d %lf %lf %lf %lf %lf %lf\n",
             &temp,
             &elements[ii].node1_id,
             &elements[ii].node2_id,
             &node1[0],
             &node1[1],
             &node1[2],
             &node2[0],
             &node2[1],
             &node2[2]);
      elements[ii].fiber_type = 0;
      elements[ii].node1_id--; // File counts from 1 not 0, so subtract 1
      elements[ii].node2_id--;
      nodes[elements[ii].node1_id].x = node1[0];
      nodes[elements[ii].node1_id].y = node1[1];
      nodes[elements[ii].node1_id].z = node1[2];
      nodes[elements[ii].node2_id].x = node2[0];
      nodes[elements[ii].node2_id].y = node2[1];
      nodes[elements[ii].node2_id].z = node2[2];
      elements[ii].orig_len = sqrt(pow((node2[0] - node1[0]),2) +
                                   pow((node2[1] - node1[1]),2) +
                                   pow((node2[2] - node1[2]),2));
    }
    // Get periodic connections
    // Or could redetermine these connections instead of reading from file
    for(int ii = num_rl_bcs; ii < num_rl_bcs+num_new_rl_bcs; ii++)
    {
      fscanf(file,"%d %d\n",
             &rl_bcs[ii].node1_id,
             &rl_bcs[ii].node2_id);
      rl_bcs[ii].node1_id--;
      rl_bcs[ii].node2_id--;
    }
    for(int ii = num_tb_bcs; ii < num_tb_bcs+num_new_tb_bcs; ii++)
    {
      fscanf(file,"%d %d\n",
             &tb_bcs[ii].node1_id,
             &tb_bcs[ii].node2_id);
      tb_bcs[ii].node1_id--;
      tb_bcs[ii].node2_id--;
    }
    for(int ii = num_fb_bcs; ii < num_fb_bcs+num_new_fb_bcs; ii++)
    {
      fscanf(file,"%d %d\n",
             &fb_bcs[ii].node1_id,
             &fb_bcs[ii].node2_id);
      fb_bcs[ii].node1_id--;
      fb_bcs[ii].node2_id--;
    }
    collectPeriodicConnectionInfo(rl_bcs);
    collectPeriodicConnectionInfo(tb_bcs);
    collectPeriodicConnectionInfo(fb_bcs);
    // Check end of fiber network file for rve boundary location
    //  if not, then use default values from globals
    side[BOTTOM] = -0.5;
    side[TOP]    =  0.5;
    side[LEFT]   = -0.5;
    side[RIGHT]  =  0.5;
    side[BACK]   = -0.5;
    side[FRONT]  =  0.5;
    // Initial rotation for each node is 0 (only used for beams)
    for(int ii = num_nodes; ii < num_nodes+num_new_nodes; ii++)
    {
      nodes[ii].rx = 0.0;
      nodes[ii].ry = 0.0;
      nodes[ii].rz = 0.0;
    }
    fclose(file);
    bound.resize(BOUNDARY_TYPES); // For each boundary type
    for(int ii = num_nodes; ii < num_nodes+num_new_nodes; ii++)
    {
      bool top    = fabs(nodes[ii].y - side[TOP])    < 1e-8;
      bool bottom = fabs(nodes[ii].y - side[BOTTOM]) < 1e-8;
      bool right  = fabs(nodes[ii].x - side[RIGHT])  < 1e-8;
      bool left   = fabs(nodes[ii].x - side[LEFT])   < 1e-8;
      bool front  = fabs(nodes[ii].z - side[FRONT])  < 1e-8;
      bool back   = fabs(nodes[ii].z - side[BACK])   < 1e-8;
      if( top || bottom || right || left || front || back )
        bound[ALL].push_back(ii);
      if( top )
        bound[TOP].push_back(ii);
      if( bottom )
        bound[BOTTOM].push_back(ii);
      if( left )
        bound[LEFT].push_back(ii);
      if( right )
        bound[RIGHT].push_back(ii);
      if( back )
        bound[BACK].push_back(ii);
      if( front )
        bound[FRONT].push_back(ii);
    }
    num_nodes += num_new_nodes;
    num_elements += num_new_elements;
    num_rl_bcs += num_new_rl_bcs;
    num_tb_bcs += num_new_tb_bcs;
    num_fb_bcs += num_new_fb_bcs;
  }
  /* Redefinition of readFromFile for SupportFiberNetwork class */
  void SupportFiberNetwork::readFromFile(const std::string & filename)
  {
    AMSI_DEBUG
      (
        std::cout << "Attempting to access file: " << filename.c_str() << " for reading." << std::endl;
        )
      FILE * file = fopen(filename.c_str(),"r");
    assert(file);
    // read header info
    fseek(file,0,SEEK_SET);
    // Assume no periodic connections
    int num_new_rl_bcs = 0;
    int num_new_tb_bcs = 0;
    int num_new_fb_bcs = 0;
    int num_new_nodes = 0;
    int num_new_dofs = 0;
    int num_new_elements = 0;
    // Periodic has extra info in header (total numbers of periodic pairs)
    fscanf(file,"%d %d %d\n",
           &num_new_nodes,
           &num_new_dofs,
           &num_new_elements);
    num_dofs += dofs_per_node*num_new_nodes;
    elements.resize(num_elements+num_new_elements);
    nodes.resize(num_nodes+num_new_nodes);
    rl_bcs.resize(num_rl_bcs+num_new_rl_bcs);
    tb_bcs.resize(num_tb_bcs+num_new_tb_bcs);
    fb_bcs.resize(num_fb_bcs+num_new_fb_bcs);
    // Fill element array
    for(int ii = num_elements; ii < num_elements+num_new_elements; ii++)
    {
      double node1[3], node2[3];
      int temp;
      fscanf(file,"%d %d %d %lf %lf %lf %lf %lf %lf %i\n",
             &temp,
             &elements[ii].node1_id,
             &elements[ii].node2_id,
             &node1[0],
             &node1[1],
             &node1[2],
             &node2[0],
             &node2[1],
             &node2[2],
             &elements[ii].fiber_type);
      elements[ii].node1_id--; // File counts from 1 not 0, so subtract 1
      elements[ii].node2_id--;
      nodes[elements[ii].node1_id].x = node1[0];
      nodes[elements[ii].node1_id].y = node1[1];
      nodes[elements[ii].node1_id].z = node1[2];
      nodes[elements[ii].node2_id].x = node2[0];
      nodes[elements[ii].node2_id].y = node2[1];
      nodes[elements[ii].node2_id].z = node2[2];
      elements[ii].orig_len = sqrt(pow((node2[0] - node1[0]),2) +
                                   pow((node2[1] - node1[1]),2) +
                                   pow((node2[2] - node1[2]),2));
    }
    // Get periodic connections
    // Or could redetermine these connections instead of reading from file
    for(int ii = num_rl_bcs; ii < num_rl_bcs+num_new_rl_bcs; ii++)
    {
      fscanf(file,"%d %d\n",
             &rl_bcs[ii].node1_id,
             &rl_bcs[ii].node2_id);
      rl_bcs[ii].node1_id--;
      rl_bcs[ii].node2_id--;
    }
    for(int ii = num_tb_bcs; ii < num_tb_bcs+num_new_tb_bcs; ii++)
    {
      fscanf(file,"%d %d\n",
             &tb_bcs[ii].node1_id,
             &tb_bcs[ii].node2_id);
      tb_bcs[ii].node1_id--;
      tb_bcs[ii].node2_id--;
    }
    for(int ii = num_fb_bcs; ii < num_fb_bcs+num_new_fb_bcs; ii++)
    {
      fscanf(file,"%d %d\n",
             &fb_bcs[ii].node1_id,
             &fb_bcs[ii].node2_id);
      fb_bcs[ii].node1_id--;
      fb_bcs[ii].node2_id--;
    }
    collectPeriodicConnectionInfo(rl_bcs);
    collectPeriodicConnectionInfo(tb_bcs);
    collectPeriodicConnectionInfo(fb_bcs);
    side[BOTTOM] = -0.5;
    side[TOP]    =  0.5;
    side[LEFT]   = -0.5;
    side[RIGHT]  =  0.5;
    side[BACK]   = -0.5;
    side[FRONT]  =  0.5;
    // Initial rotation for each node is 0 (only used for beams)
    for(int ii = num_nodes; ii < num_nodes+num_new_nodes; ii++)
    {
      nodes[ii].rx = 0.0;
      nodes[ii].ry = 0.0;
      nodes[ii].rz = 0.0;
    }
    fclose(file);
    bound.resize(BOUNDARY_TYPES); // For each boundary type
    for(int ii = num_nodes; ii < num_nodes+num_new_nodes; ii++)
    {
      bool top    = fabs(nodes[ii].y - side[TOP])    < 1e-8;
      bool bottom = fabs(nodes[ii].y - side[BOTTOM]) < 1e-8;
      bool right  = fabs(nodes[ii].x - side[RIGHT])  < 1e-8;
      bool left   = fabs(nodes[ii].x - side[LEFT])   < 1e-8;
      bool front  = fabs(nodes[ii].z - side[FRONT])  < 1e-8;
      bool back   = fabs(nodes[ii].z - side[BACK])   < 1e-8;
      if( top || bottom || right || left || front || back )
        bound[ALL].push_back(ii);
      if( top )
        bound[TOP].push_back(ii);
      if( bottom )
        bound[BOTTOM].push_back(ii);
      if( left )
        bound[LEFT].push_back(ii);
      if( right )
        bound[RIGHT].push_back(ii);
      if( back )
        bound[BACK].push_back(ii);
      if( front )
        bound[FRONT].push_back(ii);
    }
    num_nodes += num_new_nodes;
    num_elements += num_new_elements;
    num_rl_bcs += num_new_rl_bcs;
    num_tb_bcs += num_new_tb_bcs;
    num_fb_bcs += num_new_fb_bcs;
    // Populate support vector with nodes corresponding to support fibers
    for (int ii = 0; ii < num_elements; ii++)
    {
      if (elements[ii].fiber_type != 0)
        support.push_back(elements[ii].node2_id);
    }
    num_support_nodes = support.size();
  }
  void FiberNetwork::output(const std::string & filename) const
  {
    FILE * fp = fopen(filename.c_str(), "w");
    fprintf(fp, "%d %d %d\n", num_nodes, 3*num_elements, num_elements);
    for(int ii = 0; ii < num_elements; ii++)
    {
      int node1 = elements[ii].node1_id;
      int node2 = elements[ii].node2_id;
      fprintf(fp, "%d %d %d %.16e %.16e %.16e %.16e %.16e %.16e %i\n",
              ii+1,
              node1,
              node2,
              nodes[node1].x,
              nodes[node1].y,
              nodes[node1].z,
              nodes[node2].x,
              nodes[node2].y,
              nodes[node2].z,
              elements[ii].fiber_type);
    }
    fclose(fp);
  }
  void FiberNetwork::outputVTK(const std::string & filename) const
  {
    FILE * fp = fopen(filename.c_str(), "w");
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Fibers \n");
    fprintf(fp, "ASCII\n\n");
    // write out the points
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp, "POINTS %d float\n", num_nodes);
    for(int ii = 0; ii < num_nodes; ii++)
      fprintf(fp, "%.16e %.16e %.16e\n",
              nodes[ii].x,
              nodes[ii].y,
              nodes[ii].z);
    fprintf(fp, "\n");
    fprintf(fp, "CELLS %d %d\n", num_elements, 3*num_elements);
    for(int ii = 0; ii < num_elements; ii++)
      fprintf(fp, "%d %d %d\n",
              2,
              elements[ii].node1_id,
              elements[ii].node2_id);
    fprintf(fp, "\n");
    fprintf(fp, "CELL_TYPES %d\n", num_elements);
    for(int ii = 0; ii < num_elements; ii++)
      fprintf(fp, "%d\n", 3);
    fprintf(fp, "\n");
    fprintf(fp, "POINT_DATA %d\n", num_nodes);
    fprintf(fp, "SCALARS fibers float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for(int ii = 0; ii < num_nodes; ii++)
      fprintf(fp, "%.9f \n", (double)(ii%5)+1); // ????
    fclose(fp);
  }
}
