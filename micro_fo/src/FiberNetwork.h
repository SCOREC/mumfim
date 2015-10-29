#ifndef FIBERNETWORK_H_
#define FIBERNETWORK_H_

#include "FiberReactions.h"

#include <string>
#include <vector>

namespace Biotissue {

struct Element
{  
  int node1_id;    // index of first fiber node
  int node2_id;    // index of second fiber node
  
  double orig_len; // initial (undeformed) length of the fiber
  std::vector<int> Ke; // Will be either size 36 = 6^2 for truss or 144 = 12^2 for beam
  int fiber_type;

  Element(const Element & e)
  {
    node1_id = e.node1_id;
    node2_id = e.node2_id;
    orig_len = e.orig_len;
    Ke.resize(e.Ke.size());
    std::copy(e.Ke.begin(), e.Ke.end(), Ke.begin());
    fiber_type = e.fiber_type;
  }
};

class Node
{
protected:
  int ndofs;
  double * dofs;
  virtual init() : dofs() { ndofs = 3; dofs = new double[ndofs];}
public:
  Node() : ndofs(), dofs() { init(); }
  virtual int numDofs() {return 3;}
  double * getDofs() {return dofs;}
  const double & dof(int idx)
  {
    assert(idx < ndofs);
    return dofs[idx];
  }
};

struct PBCRelation
{
  int node1_id; // Node on right, top, or front
  int node2_id; // Node on left, bottom, or back

  int elem1;  // Element which node1_id is in
  int elem2;  // Element which node2_id is in
  bool elem1first; // True - node1_id is first node in elem1 , False - node1_id is second node in elem1
  bool elem2first;

  std::vector<int> Ke1; // Sparse matrix reference for node 1
  std::vector<int> Ke2; // Sparse matrix reference for node 2
};

class FiberNetwork
{
public:
  enum Side
  {
    TOP = 0,
    BOTTOM = 1,
    LEFT = 2,
    RIGHT = 3,
    FRONT = 4,
    BACK = 5,
    ALL = 6
  };

  // CONSTRUCTORS ============================================ //
  
  FiberNetwork();
  virtual FiberNetwork * clone();

  // IO ====================================================== //
  
  virtual void readFromFile(const std::string & filename);  
  void output(const std::string & filename) const;
  void outputVTK(const std::string & filename) const;

  // GET/SET-ERS ============================================ //
  
  int numDofs() const {return num_dofs;}
  int numNodes() const {return num_nodes;}
  int numElements() const {return num_elements;}
  int dofsPerNode() const {return dofs_per_node;}

  double fiberRadius() const {return fiber_radius;}
  double fiberArea() const {return fiber_area;}
  
  const Node & node(int index) const {return nodes[index];}
  const Element & element(int index) const {return elements[index];}
  const FiberReaction * getFiberParams(int index) const { return fiber_types[elements[index].fiber_type];}

  double sideCoord(Side s) const {return side[s];}  

  void setNode(size_t index, const Node & n);
  void setElement(size_t index, const Element & e);
  
  void setAllCoordinates(double * coords);
  void getAllCoordinates(std::vector<double> & ns);

  // Periodic BCs
  int numRLBCs() {return num_rl_bcs;} // Right-left
  int numTBBCs() {return num_tb_bcs;} // Top-bottom
  int numFBBCs() {return num_fb_bcs;} // Front-back

  const std::vector<PBCRelation> & getRightLeftBCs() { return rl_bcs; }
  const std::vector<PBCRelation> & getTopBottomBCs() { return tb_bcs; }
  const std::vector<PBCRelation> & getFrontBackBCs() { return fb_bcs; }
  
protected:

  // FUNCTIONS ============================================== //

  FiberNetwork(const FiberNetwork &);
  void collectPeriodicConnectionInfo(std::vector<PBCRelation> & bcs);

  // MEMBER VARIABLES ======================================= //
  
  std::vector<FiberReaction*> fiber_types;

  std::vector<Element> elements;
  std::vector<Node> nodes;

  int num_nodes;
  int num_dofs;
  int num_elements;
  int dofs_per_node;

  double side[ALL];

  // these may get pushed into the FiberReaction structures if they vary per fiber
  double fiber_radius;
  double fiber_area;

  int num_rl_bcs;
  int num_tb_bcs;
  int num_fb_bcs;

  std::vector<PBCRelation> rl_bcs;
  std::vector<PBCRelation> tb_bcs;
  std::vector<PBCRelation> fb_bcs;
};

class SupportFiberNetwork : public FiberNetwork
{
private:

public:
  SupportFiberNetwork();
  virtual FiberNetwork * clone();
  
  virtual void readFromFile(const std::string & filename);
  
  int numSupportNodes() const{return num_support_nodes;}
  int supportNode(size_t index) const{return support[index];}
  
protected:
  SupportFiberNetwork(const SupportFiberNetwork &); //Copy constructor.

  //vector to store support nodes. Analagous to bound in FiberNetwork superclass.
  std::vector<int> support; 
  int num_support_nodes;
};

bool onBoundary(FiberNetwork & nm, const Node & n, FiberNetwork::Side s);

void gatherBoundaryNodes(FiberNetwork & fn,
			 FiberNetwork::Side s,
			 std::vector<int> & nds)

// todo (h) : rather than taking nodes, these should take elements themselves

void calcFiberLengths(const FiberNetwork & fn,
		      std::vector<double> & lengths);
inline double calcFiberLength(const Node & n1,
			      const Node & n2)
{
  return sqrt((n2.x-n1.x)*(n2.x-n1.x) + (n2.y-n1.y)*(n2.y-n1.y) + (n2.z-n1.z)*(n2.z-n1.z));
}

double calcNetworkOrientation(const FiberNetwork & fn);
void calcFiberOrientation(const Node & n1,
			  const Node & n2,
			  double (&rslt)[3]);

void calcAvgFiberDirection(const FiberNetwork & fn,
			   double (&rslt)[3]);

void calcFiberDirection(const Node & n1,
			const Node & n2,
			double (&rslt)[3]);
}
#endif
