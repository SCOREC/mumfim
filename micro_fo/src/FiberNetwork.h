#ifndef FIBERNETWORK_H_
#define FIBERNETWORK_H_
#include "FiberReactions.h"
#include <string>
#include <vector>
namespace bio
{
  struct Element
  {
    int node1_id;    // index of first fiber node
    int node2_id;    // index of second fiber node
    double orig_len; // initial (undeformed) length of the fiber
    std::vector<int> Ke; // Will be either size 36 = 6^2 for truss or 144 = 12^2 for beam
    int fiber_type;
    Element() {};
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
// Nodes are equivalent to vertices in this formulation, as 2 node truss and
// beam elements have the two nodes at the two ends of the element
  struct Node
  {
    double x; // x-coordinate of the node
    double y; // y-coordinate of the node
    double z; // z-coordinate of the node
    double rx;
    double ry;
    double rz;
  Node() : x(0.0), y(0.0), z(0.0), rx(0.0), ry(0.0), rz(0.0) {};
    Node(const Node & n)
      {
        x = n.x;
        y = n.y;
        z = n.z;
        rx = n.rx;
        ry = n.ry;
        rz = n.rz;
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
      ALL = 6,
      BOUNDARY_TYPES = 7
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
    // Periodic BCs
    int numRLBCs() {return num_rl_bcs;} // Right-left
    int numTBBCs() {return num_tb_bcs;} // Top-bottom
    int numFBBCs() {return num_fb_bcs;} // Front-back
    const std::vector<PBCRelation> & getRightLeftBCs() { return rl_bcs; }
    const std::vector<PBCRelation> & getTopBottomBCs() { return tb_bcs; }
    const std::vector<PBCRelation> & getFrontBackBCs() { return fb_bcs; }
    const Node & node(int index) const {return nodes[index];}
    const Element & element(int index) const {return elements[index];}
    double sideCoord(Side s) const {return side[s];}
    void updateSideCoord(Side s, double val) {side[s] += val;}
    void setNode(size_t index, const Node & n);
    void setElement(size_t index, const Element & e);
    void setNodeCoordinates(double * coords);
    void getNodeCoordinates(std::vector<double> & ns);
    int numBoundaryNodes(Side side) const;
    int boundaryNode(Side side, size_t index) const;
    void updateOrigLen();
  protected:
    // FUNCTIONS ============================================== //
    FiberNetwork(const FiberNetwork &);
    void collectPeriodicConnectionInfo(std::vector<PBCRelation> & bcs);
    // MEMBER VARIABLES ======================================= //
    std::vector< std::vector<int> > bound;
    std::vector<Element> elements;
    std::vector<Node> nodes;
    int num_nodes;
    int num_dofs;
    int num_elements;
    int dofs_per_node;
    double side[ALL];
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
  void calcFiberOrientationTensor(const FiberNetwork & fn,
				  double (&rslt)[9]);
  void calcP2(const FiberNetwork & fn, double & rslt);
  void calcAvgFiberDirection(const FiberNetwork & fn,
                             double (&rslt)[3]);
  void calcFiberDirection(const Node & n1,
                          const Node & n2,
                          double (&rslt)[3]);
  void AlignFiberNetwork(FiberNetwork & fn, const double align_vec[3]);
  void AffineDeformation(FiberNetwork & fn, const double disp[6]);
  void updateRVEBounds(FiberNetwork & fn, const double disp[6]);
  double calcFiberDensity(const FiberNetwork & fn);
}
#endif
