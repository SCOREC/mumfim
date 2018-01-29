#include <string.h>
#include <iostream>
namespace bio
{
  /*
  void fill_dof_ids(int * dof_ids, int n1, int n2)
  {
    if(TRUSS)
    {
      dof_ids[0] = n1 * 3;
      dof_ids[1] = n1 * 3 + 1;
      dof_ids[2] = n1 * 3 + 2;
      dof_ids[3] = n2 * 3;
      dof_ids[4] = n2 * 3 + 1;
      dof_ids[5] = n2 * 3 + 2;
    }
    else // Beam
    {
      dof_ids[0] = n1 * 6;
      dof_ids[1] = n1 * 6 + 1;
      dof_ids[2] = n1 * 6 + 2;
      dof_ids[3] = n1 * 6 + 3;
      dof_ids[4] = n1 * 6 + 4;
      dof_ids[5] = n1 * 6 + 5;
      dof_ids[6]  = n2 * 6;
      dof_ids[7]  = n2 * 6 + 1;
      dof_ids[8]  = n2 * 6 + 2;
      dof_ids[9]  = n2 * 6 + 3;
      dof_ids[10] = n2 * 6 + 4;
      dof_ids[11] = n2 * 6 + 5;
    }
  }
  void fill_dof_ids_1(int * dof_ids, int n)
  {
    if(TRUSS)
    {
      dof_ids[0] = n * 3;
      dof_ids[1] = n * 3 + 1;
      dof_ids[2] = n * 3 + 2;
    }
    else // Beam
    {
      dof_ids[0] = n * 6;
      dof_ids[1] = n * 6 + 1;
      dof_ids[2] = n * 6 + 2;
      dof_ids[3] = n * 6 + 3;
      dof_ids[4] = n * 6 + 4;
      dof_ids[5] = n * 6 + 5;
    }
  }
  SparseMatrix * Make_Structure(FiberNetwork * fiber_network)
  {
    int num_dofs = fiber_network->numDofs();
    int num_dofs_per_element = 2 * fiber_network->dofsPerNode();
    SparseMatrix * sparse_struct = new SparseMatrix(num_dofs);
    int dof_ids[num_dofs_per_element];
    int num_elements = fiber_network->numElements();
    for(int ii = 0; ii < num_elements; ii++)
    {
      const Element & e = fiber_network->element(ii);
      fill_dof_ids(dof_ids,
                   e.node1_id,
                   e.node2_id);
      for(int jj = 0; jj < num_dofs_per_element; jj++)
        for(int kk = 0; kk < num_dofs_per_element; kk++)
        {
          if(dof_ids[kk] > num_dofs || dof_ids[jj] > num_dofs)
            std::cerr << "ERROR: dof id out of range!" << std::endl;
          if(sparse_struct->needMtxElem(dof_ids[kk],dof_ids[jj]))
            sparse_struct->addMtxElem(dof_ids[kk],dof_ids[jj],num_dofs);
        }
    }
    if(PERIODIC)
    {
      const std::vector<PBCRelation> & rl_bcs = fiber_network->getRightLeftBCs();
      const std::vector<PBCRelation> & tb_bcs = fiber_network->getTopBottomBCs();
      const std::vector<PBCRelation> & fb_bcs = fiber_network->getFrontBackBCs();
      fill_extra_periodic_connections(rl_bcs,fiber_network,sparse_struct);
      fill_extra_periodic_connections(tb_bcs,fiber_network,sparse_struct);
      fill_extra_periodic_connections(fb_bcs,fiber_network,sparse_struct);
    }
    // sparse structure should be finalized here..
    sparse_struct->finalize();
    // this actually shouldn't really go here since it modifies the fibernetwork... should be a better place to put this
    for(int ii = 0; ii < num_elements; ii++)
    {
      // this is very slow and we should not do this, find a better way
      Element e = fiber_network->element(ii);
      fill_dof_ids(dof_ids,
                   e.node1_id,
                   e.node2_id);
      e.Ke.resize(num_dofs_per_element*num_dofs_per_element);
      int count = 0;
      for(int jj = 0; jj < num_dofs_per_element; jj++)
        for(int kk = 0; kk < num_dofs_per_element; kk++)
          e.Ke[count++] = sparse_struct->sparseLocation(dof_ids[kk],dof_ids[jj]);
      fiber_network->setElement(ii,e);
    }
    if(PERIODIC)
    {
      const std::vector<PBCRelation> & rl_bcs = fiber_network->getRightLeftBCs();
      const std::vector<PBCRelation> & tb_bcs = fiber_network->getTopBottomBCs();
      const std::vector<PBCRelation> & fb_bcs = fiber_network->getFrontBackBCs();
      store_periodic_locations(rl_bcs,fiber_network,sparse_struct);
      store_periodic_locations(tb_bcs,fiber_network,sparse_struct);
      store_periodic_locations(fb_bcs,fiber_network,sparse_struct);
    }
    return sparse_struct;
  }
  void fill_extra_periodic_connections(const std::vector<PBCRelation> & bcs,
                                       FiberNetwork * fiber_network,
                                       SparseMatrix * sparse_struct)
  {
    int num_dofs = fiber_network->numDofs();
    int num_dofs_per_node = fiber_network->dofsPerNode();
    int num_dofs_per_element = 2 * num_dofs_per_node;
    int dof_ids[num_dofs_per_element];
    int dof_ids_1[fiber_network->dofsPerNode()];
    int dof_ids_2[fiber_network->dofsPerNode()];
    // If using periodic BCs then need a few more locations in the sparse matrix
    int num_bcs = bcs.size();
    for(int ii=0;ii<num_bcs;ii++)
    {
      int pnode1 = bcs[ii].node1_id;
      int pnode2 = bcs[ii].node2_id;
      int node1;
      if(bcs[ii].elem1first)
      {
        node1 = fiber_network->element(bcs[ii].elem1).node2_id;
        fill_dof_ids(dof_ids, pnode1, node1);
      }
      else
      {
        node1 = fiber_network->element(bcs[ii].elem1).node1_id;
        fill_dof_ids(dof_ids, node1, pnode1);
      }
      fill_dof_ids_1(dof_ids_1,pnode1);
      fill_dof_ids_1(dof_ids_2,pnode2);
      // Equations for pnode1 will be added to equations for pnode2
      for(int jj = 0; jj < num_dofs_per_node; jj++)
        for(int kk = 0; kk < num_dofs_per_element; kk++)
        {
          if(dof_ids[kk] > num_dofs || dof_ids_2[jj] > num_dofs)
            std::cerr << "ERROR: dof id out of range!" << std::endl;
          if(sparse_struct->needMtxElem(dof_ids[kk],dof_ids_2[jj]))
            sparse_struct->addMtxElem(dof_ids[kk],dof_ids_2[jj],num_dofs);
        }
      // Also equations for pnode1 will be set to equality conditions
      for(int jj = 0; jj < num_dofs_per_node; jj++)
      {
        if(dof_ids_2[jj] > num_dofs || dof_ids_1[jj] > num_dofs)
          std::cerr << "ERROR: dof id out of range!" << std::endl;
        if(sparse_struct->needMtxElem(dof_ids_1[jj],dof_ids_1[jj]))
          sparse_struct->addMtxElem(dof_ids_1[jj],dof_ids_1[jj],num_dofs);
        if(sparse_struct->needMtxElem(dof_ids_2[jj],dof_ids_1[jj]))
          sparse_struct->addMtxElem(dof_ids_2[jj],dof_ids_1[jj],num_dofs);
      }
    }
  }
  void store_periodic_locations(std::vector<PBCRelation> & bcs,FiberNetwork * fiber_network,SparseMatrix * sparse_struct)
  {
    int num_dofs_per_node = fiber_network->dofsPerNode();
    int num_dofs_per_element = 2 * num_dofs_per_node;
    int dof_ids[num_dofs_per_element];
    int dof_ids_1[fiber_network->dofsPerNode()];
    int dof_ids_2[fiber_network->dofsPerNode()];
    int num_bcs = bcs.size();
    for(int ii=0;ii<num_bcs;ii++)
    {
      int pnode1 = bcs[ii].node1_id;
      int pnode2 = bcs[ii].node2_id;
      int node1;
      if(bcs[ii].elem1first)
      {
        node1 = fiber_network->element(bcs[ii].elem1).node2_id;
        fill_dof_ids(dof_ids, pnode1, node1);
      }
      else
      {
        node1 = fiber_network->element(bcs[ii].elem1).node1_id;
        fill_dof_ids(dof_ids, node1, pnode1);
      }
      // Get relevant dofs
      fill_dof_ids_1(dof_ids_1,pnode1);
      fill_dof_ids_1(dof_ids_2,pnode2);
      // Ke1 is used for referencing sparse matrix locations when adding equations
      // from pnode1 to equations of pnode2
      bcs[ii].Ke1.resize(num_dofs_per_node*num_dofs_per_element);
      int count = 0;
      for(int jj = 0; jj < num_dofs_per_node; jj++)
        for(int kk = 0; kk < num_dofs_per_element; kk++)
          bcs[ii].Ke1[count++] = sparse_struct->sparseLocation(dof_ids[kk],dof_ids_2[jj]);
      // Ke2 is used for referencing sparse matrix locations when setting equality
      // conditions for equations associated with pnode1
      bcs[ii].Ke2.resize(2*num_dofs_per_node);
      count = 0;
      for(int jj = 0; jj < num_dofs_per_node; jj++)
      {
        bcs[ii].Ke2[count++] = sparse_struct->sparseLocation(dof_ids_1[jj],dof_ids_1[jj]);
        bcs[ii].Ke2[count++] = sparse_struct->sparseLocation(dof_ids_2[jj],dof_ids_1[jj]);
      }
    }
  }
  */
}
