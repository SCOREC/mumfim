#include "FiberNetwork.h"

#include "globals.h" // CUBE_*, SPECIFY_FIBER_TYPE
#include "Util.h"

#include <amsiConfig.h>
#include <iostream> 
#include <cassert> // assert

namespace bio
{

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
    fiber_types.resize(fn.fiber_types.size());

    fiber_radius = fn.fiber_radius;
    fiber_area = fn.fiber_area;
    
    std::copy(fn.fiber_types.begin(),fn.fiber_types.end(),fiber_types.begin());
    std::copy(fn.nodes.begin(),fn.nodes.end(),nodes.begin());
    std::copy(fn.elements.begin(),fn.elements.end(),elements.begin());
  }

  FiberNetwork::FiberNetwork(int nnds, int nfbrs) :
    elements(nfbrs),
    nodes(nnds),
    num_nodes(nnds),
    num_dofs(0),
    num_elements(nfbrs),
    side(),
    fiber_radius(3.49911271e-08),
    fiber_area(M_PI * fiber_radius * fiber_radius)
  {
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

  void FiberNetwork::setAllCoordinates(double * coords)
  {
    assert(coords);
    int ii = 0;
    double * head = coords;
    for(std::vector<Node>::iterator it = nodes.begin(); it != nodes.end(); it++)
    {
      int num_dofs = it->numDofs();
      it->setDofs(head);
      head += num_dofs;
    }
  }

  void FiberNetwork::getAllCoordinates(std::vector<double> & ns)
  {
    ns.erase(ns.begin(),ns.end());
    std::back_insert_iterator< std::vector<double> > bit(ns);
    for(std::vector<Node>::iterator it = nodes.begin(); it != nodes.end(); it++)
    {
      int num_dofs = it->numDofs();
      double * dofs = it->getDofs();
      std::copy(dofs,dofs+num_dofs,bit);
    }
  }
   
  SupportFiberNetwork::SupportFiberNetwork() : FiberNetwork()
  {
    NonlinearReaction * r1 = new NonlinearReaction;
    fiber_types[0] = r1;
    r1->fiber_area = fiber_area;
    r1->B = 1.2;
    r1->E = 43200000; // 43.2 MPa
    r1->lexp = 3.0; // limit of the length ratio

    NonlinearReaction * r2 = new NonlinearReaction;
    fiber_types[1] = r2;
    r2->fiber_area = fiber_area;
    r2->B = 1.2;
    r2->E = 432000;
    r2->lexp = 3.0;
    
  }

  SupportFiberNetwork::SupportFiberNetwork(const SupportFiberNetwork & spfn) : FiberNetwork(spfn)
  {
    num_support_nodes = spfn.num_support_nodes;
    support = spfn.support;
  }

  FiberNetwork * SupportFiberNetwork::clone()
  {
    return static_cast<FiberNetwork*>(new SupportFiberNetwork(*this));
  }

  bool onBoundary(FiberNetwork & fn, const Node & n, FiberNetwork::Side s)
  {
    int idx[] = {0,0,1,1,2,2}; // left-right, top-bottom, front-back
    bool result = false;
    double coord = 0.0;
    coord = n[idx[side]];
    result = close(coord,fn.sideCoord(s));
    return result;
  }
  
  void gatherBoundaryNodes(FiberNetwork & fn,
			   FiberNetwork::Side s,
			   std::vector<int> & nds)
  {
    nds.erase();
    int num_nodes = fn.numNodes();
    for(int ii = 0; ii < num_nodes; ii++)
    {
      const Node & n = fn.node(ii);
      if(onBoundary(fn,n,s))
	nds.push_back(n);
    }
  }
  
  void calcFiberLengths(const FiberNetwork & fn, std::vector<double> & lengths)
  {
    int num_elements = fn.numElements();
    lengths.resize(num_elements);
    for(int ii = 0; ii < num_elements; ii++)
      lengths[ii] = calcFiberLength(fn.element(ii));
  }

  double calcNetworkOrientation(const FiberNetwork & fn)
  {
    double result = 0.0;
    int num_elements = fn.numElements();
    double ortn[3] = {};
    for(int ii = 0; ii < num_elements; ii++)
    {
      const Element & e = fn.element(ii);
      calcFiberOrientation(fn.node(e.getNode(0)),
			   fn.node(e.getNode(1)),
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
    rslt[0] = (n2[0] - n1[0] / length) * axis[0];
    rslt[1] = (n2[1] - n1[1] / length) * axis[1];
    rslt[2] = (n2[2] - n1[2] / length) * axis[2];
  }

  void calcAvgFiberDirection(const FiberNetwork & fn,
			     double (&rslt)[3])
  {
    int num_elements = fn.numElements();
    double dir[3] = {};
    for(int ii = 0; ii < num_elements; ii++)
    {
      const Element & e = fn.element(ii);
      calcFiberDirection(fn.node(e.getNode(0)),
			 fn.node(e.getNode(1)),
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
    rslt[0] = n2[0] - n1[0];
    rslt[1] = n2[1] - n2[1];
    rslt[2] = n2[2] - n2[2];
    if(rslt[0] < 0)
    {
      rslt[0] *= -1.0;
      rslt[1] *= -1.0;
      rslt[2] *= -1.0;
    }
  }
}
