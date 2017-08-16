#include "bioFiberNetwork.h"
#include <iostream>
#include <cassert> // assert
#include <numeric> // accumulate
namespace bio
{
  FiberNetwork::FiberNetwork(apf::Mesh * f)
    : fn(f)
    , u(NULL)
    , du(NULL)
    , dw(NULL)
    , udof(NULL)
    , wdof(NULL)
    , ucnt(0)
    , wcnt(0)
    , tp(FiberMember::truss)
    , dim(f->getDimension())
  {
    assert(f);
    du = apf::createLagrangeField(fn,"du",apf::VECTOR,1);
    u  = apf::createLagrangeField(fn,"u",apf::VECTOR,1);
    udof = apf::createNumbering(du);
    //wdof = apf::createNumbering(dw);
    ucnt = apf::NaiveOrder(udof);
    //wcnt = apf::AdjReorder(wdof);
    //apf::SetNumberingOffset(wdof,ucnt);
  }
  /*
  FiberNetwork * FiberNetwork::clone()
  {
    return new FiberNetwork(*this);
  }
  FiberNetwork::FiberNetwork(const FiberNetwork & fn)
  {
    num_nodes = fn.numNodes();
    num_dofs = fn.numDofs();
    num_elements = fn.numElements();
    dofs_per_node = fn.dofsPerNode();
    for(int ii = 0; ii < ALL; ii++)
      side[ii] = fn.side[ii];
    nodes.resize(num_nodes);
    elements.resize(num_elements);
    std::copy(fn.nodes.begin(),fn.nodes.end(),nodes.begin());
    std::copy(fn.elements.begin(),fn.elements.end(),elements.begin());
    bound.resize(BOUNDARY_TYPES); // For each boundary type
    for(int ii = 0; ii < BOUNDARY_TYPES; ii++)
    {
      bound[ii].resize(fn.numBoundaryNodes((Side)ii));
      std::copy(fn.bound[ii].begin(),fn.bound[ii].end(),bound[ii].begin());
    }
    // Periodic Connections
    std::copy(fn.rl_bcs.begin(),fn.rl_bcs.end(),rl_bcs.begin());
    std::copy(fn.tb_bcs.begin(),fn.tb_bcs.end(),tb_bcs.begin());
    std::copy(fn.fb_bcs.begin(),fn.fb_bcs.end(),fb_bcs.begin());
    num_rl_bcs = rl_bcs.size();
    num_tb_bcs = tb_bcs.size();
    num_fb_bcs = fb_bcs.size();
    side[TOP] = fn.sideCoord(TOP);
    side[BOTTOM] = fn.sideCoord(BOTTOM);
    side[LEFT] = fn.sideCoord(LEFT);
    side[RIGHT] = fn.sideCoord(RIGHT);
    side[BACK] = fn.sideCoord(BACK);
    side[FRONT] = fn.sideCoord(FRONT);
  }
  FiberNetwork::FiberNetwork()
    : bound()
    , elements()
    , nodes()
    , num_nodes(0)
    , num_dofs(0)
    , num_elements(0)
    , dofs_per_node(0)
    , side()
    , num_rl_bcs(0)
    , num_tb_bcs(0)
    , num_fb_bcs(0)
    , rl_bcs()
    , tb_bcs()
    , fb_bcs()
  {
    dofs_per_node = 3;
  }
  void FiberNetwork::collectPeriodicConnectionInfo(std::vector<PBCRelation> & bcs)
  {
    for(int ii = 0; ii < (int)bcs.size(); ii++)
    {
      for(int jj = 0; jj < num_elements; jj++)
      {
        const Element & e = element(jj);
        if(e.node1_id == bcs[ii].node1_id)
        {
          bcs[ii].elem1 = jj;
          bcs[ii].elem1first = true;
          break;
        }
        else if(e.node2_id == bcs[ii].node1_id)
        {
          bcs[ii].elem1 = jj;
          bcs[ii].elem1first = false;
          break;
        }
      }
      for(int jj=0;jj<num_elements;jj++)
      {
        const Element & e = element(jj);
        if(e.node1_id == bcs[ii].node2_id)
        {
          bcs[ii].elem2 = jj;
          bcs[ii].elem2first = true;
          break;
        }
        else if(e.node2_id == bcs[ii].node2_id)
        {
          bcs[ii].elem2 = jj;
          bcs[ii].elem2first = false;
          break;
        }
      }
    }
  }
  int FiberNetwork::numBoundaryNodes(Side side) const
  {
    return bound[side].size();
  }
  int FiberNetwork::boundaryNode(Side side, size_t index) const
  {
    assert(index < bound[side].size());
    return bound[side][index];
  }
  void FiberNetwork::setNode(size_t index, const Node & n)
  {
    assert(index < nodes.size());
    nodes[index] = n;
  }
  void FiberNetwork::setElement(size_t index, const Element & e)
  {
    assert(index < elements.size());
    elements[index] = e;
  }
  void FiberNetwork::setNodeCoordinates(double * coords)
  {
    int ii = 0;
    for(std::vector<Node>::iterator it = nodes.begin(); it != nodes.end(); it++)
    {
      it->x = coords[ii+0];
      it->y = coords[ii+1];
      it->z = coords[ii+2];
      ii += 3;
    }
  }
  void FiberNetwork::getNodeCoordinates(std::vector<double> & ns)
  {
    for(std::vector<Node>::iterator it = nodes.begin(); it != nodes.end(); it++)
    {
      ns.push_back(it->x);
      ns.push_back(it->y);
      ns.push_back(it->z);
    }
  }
  void FiberNetwork::updateOrigLen()
  {
    double node1[3], node2[3];
    double num_elements = elements.size();
    for (int ii = 0; ii < num_elements; ii++)
    {
      node1[0] = nodes[elements[ii].node1_id].x;
      node1[1] = nodes[elements[ii].node1_id].y;
      node1[2] = nodes[elements[ii].node1_id].z;
      node2[0] = nodes[elements[ii].node2_id].x;
      node2[1] = nodes[elements[ii].node2_id].y;
      node2[2] = nodes[elements[ii].node2_id].z;
      elements[ii].orig_len = sqrt(pow((node2[0] - node1[0]),2) +
                                   pow((node2[1] - node1[1]),2) +
                                   pow((node2[2] - node1[2]),2));
    }
  }
  SupportFiberNetwork::SupportFiberNetwork() : FiberNetwork()
  { }
  SupportFiberNetwork::SupportFiberNetwork(const SupportFiberNetwork & spfn) : FiberNetwork(spfn)
  {
    num_support_nodes = spfn.num_support_nodes;
    support = spfn.support;
  }
  FiberNetwork * SupportFiberNetwork::clone()
  {
    return static_cast<FiberNetwork*>(new SupportFiberNetwork(*this));
  }
  void calcFiberLengths(const FiberNetwork & fn, std::vector<double> & lengths)
  {
    int num_elements = fn.numElements();
    lengths.resize(num_elements);
    for(int ii = 0; ii < num_elements; ii++)
    {
      const Element & e = fn.element(ii);
      const Node & n1 = fn.node(e.node1_id);
      const Node & n2 = fn.node(e.node2_id);
      lengths[ii] = calcFiberLength(n1,n2);
    }
  }
  double calcNetworkOrientation(const FiberNetwork & fn)
  {
    double result = 0.0;
    int num_elements = fn.numElements();
    double ortn[3] = {};
    for(int ii = 0; ii < num_elements; ii++)
    {
      const Element & e = fn.element(ii);
      calcFiberOrientation(fn.node(e.node1_id),
                           fn.node(e.node2_id),
                           ortn);
      result += ortn[0]*ortn[0];
      result += ortn[1]*ortn[1];
      result += ortn[2]*ortn[2];
    }
    result /= num_elements;
    result = (3*result - 1)/2;
    return result;
  }
  void calcFiberOrientation(const Node & n1,
                            const Node & n2,
                            double (&rslt)[3])
  {
    static double axis[3] = {1.0,0.0,0.0};
    double length = calcFiberLength(n1,n2);
    rslt[0] = (n2.x - n1.x / length) * axis[0];
    rslt[1] = (n2.y - n1.y / length) * axis[1];
    rslt[2] = (n2.z - n1.z / length) * axis[2];
  }

  void calcFiberOrientationTensor(const FiberNetwork & fn,
                                  double (&rslt)[9])
  {
    int num_elements = fn.numElements();
    std::vector<double> len(num_elements);
    calcFiberLengths(fn,len);
    double total_len = 0.0;
    for (int ii = 0; ii < num_elements; ii++)
    {
      const Element & e = fn.element(ii);
      const Node & n1 = fn.node(e.node1_id);
      const Node & n2 = fn.node(e.node2_id);
      double cosalpha = (n2.x - n1.x) / len[ii];
      double cosbeta  = (n2.y - n1.y) / len[ii];
      double cosgamma = (n2.z - n1.z) / len[ii];
   
      rslt[0] += len[ii] * cosalpha * cosalpha;
      rslt[1] += len[ii] * cosalpha * cosbeta;
      rslt[2] += len[ii] * cosalpha * cosgamma;

      rslt[3] += len[ii] * cosbeta * cosalpha;
      rslt[4] += len[ii] * cosbeta * cosbeta;
      rslt[5] += len[ii] * cosbeta * cosgamma;

      rslt[6] += len[ii] * cosgamma * cosalpha;
      rslt[7] += len[ii] * cosgamma * cosbeta;
      rslt[8] += len[ii] * cosgamma * cosgamma;

      total_len += len[ii];
    }
    rslt[0] /= total_len;
    rslt[1] /= total_len;
    rslt[2] /= total_len;

    rslt[3] /= total_len;
    rslt[4] /= total_len;
    rslt[5] /= total_len;

    rslt[6] /= total_len;
    rslt[7] /= total_len;
    rslt[8] /= total_len;
  }

  void calcP2(const FiberNetwork & fn, double & rslt)
  {
    int num_elements = fn.numElements();
    std::vector<double> len(num_elements);
    std::vector<double> cos2(num_elements);
    calcFiberLengths(fn,len);
    // Specify axis of alignment.
    static double axis[3] = {1.0,0.0,0.0};
    for (int ii = 0; ii < num_elements; ii++)
    {
      const Element & e = fn.element(ii);
      const Node & n1 = fn.node(e.node1_id);
      const Node & n2 = fn.node(e.node2_id);
      // direction cosines
      double alpha = (n2.x - n1.x) / len[ii];
      double beta  = (n2.y - n1.y) / len[ii];
      double gamma = (n2.z - n1.z) / len[ii];
      cos2[ii] += ( alpha * axis[0] + beta * axis[1] + gamma * axis[2] ) * ( alpha * axis[0] + beta * axis[1] + gamma * axis[2] );
    }
    double cos2sum = std::accumulate(cos2.begin(), cos2.end(), 0.0);
    double cos2avg = cos2sum / cos2.size();
    rslt = ( 3.0 * cos2avg - 1.0 ) / 2.0;
  }
  
  void calcAvgFiberDirection(const FiberNetwork & fn,
                             double (&rslt)[3])
  {
    int num_elements = fn.numElements();
    double dir[3] = {};
    for(int ii = 0; ii < num_elements; ii++)
    {
      const Element & e = fn.element(ii);
      calcFiberDirection(fn.node(e.node1_id),
                         fn.node(e.node2_id),
                         dir);
      rslt[0] += dir[0];
      rslt[1] += dir[1];
      rslt[2] += dir[2];
    }
    rslt[0] /= num_elements;
    rslt[1] /= num_elements;
    rslt[2] /= num_elements;
  }
  void calcFiberDirection(const Node & n1,
                          const Node & n2,
                          double (&rslt)[3])
  {
    rslt[0] = n2.x - n1.x;
    rslt[1] = n2.y - n1.y;
    rslt[2] = n2.z - n1.z;
    if(rslt[0] < 0)
    {
      rslt[0] *= -1.0;
      rslt[1] *= -1.0;
      rslt[2] *= -1.0;
    }
  }
  */
}
