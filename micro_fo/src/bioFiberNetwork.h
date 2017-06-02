#ifndef BIO_FIBER_NETWORK_H_
#define BIO_FIBER_NETWORK_H_
#include "bioFiberReactions.h"
#include "bioFiber.h"
#include <apf.h>
#include <apfNumbering.h>
#include <string>
#include <vector>
namespace bio
{
  /**
   * Responsible for managing the internal state of a single fiber-network
   *  quasistatics simulation.
   */
  class FiberNetwork
  {
  private:
    apf::Mesh * fn;
    apf::Field * u;
    apf::Field * du;
    apf::Field * dw;
    apf::Numbering * udof;
    apf::Numbering * wdof;
    int ucnt;
    int wcnt;
    FiberMember tp;
    int dim;
    std::vector<FiberReaction*> rctns;
  public:
    /**
     * Construct a FiberNetwork object.
     * @param f A pointer to a fiber network mesh (contains only vertices and edges)
     *          typically loaded using the NetworkLoader classes
     */
    template <typename I>
      FiberNetwork(apf::Mesh * f, FiberMember t, I rcnt_bgn, I rctn_end);
    /**
     *  Gives the dimensionality of the managed fiber network
     *  @return the dimensionality of the fiber network (2 or 3)
     */
    int getDim() const { return fn->getDimension(); }
    int getDofCount() const { return ucnt + wcnt; }
    FiberMember getFiberMember()      { return tp;   }
    apf::Mesh * getNetworkMesh()      { return fn;   }
    apf::Field * getUField()          { return u;    }
    apf::Field * getdUField()         { return du;   }
    apf::Field * getdWField()         { return dw;   }
    apf::Numbering * getUNumbering()  { return udof; }
    apf::Numbering * getdWNumbering() { return wdof; }
    FiberReaction ** getFiberReactions() { return &rctns[0]; }
  };
  /**
   * TODO (m) Bill : move to FiberNetworkIO files
   */
  FiberNetwork * loadFromFile(const std::string & fnm);
  /*
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
  */
  /*
  void calcFiberLengths(const FiberNetwork & fn,
                        std::vector<double> & lengths);
  inline double calcFiberLength(const Node & n1,
                                const Node & n2)
  {
    return sqrt((n2.x-n1.x)*(n2.x-n1.x) + (n2.y-n1.y)*(n2.y-n1.y) + (n2.z-n1.z)*(n2.z-n1.z));
  }
  */
  double calcNetworkOrientation(const FiberNetwork & fn);
  /*
  void calcFiberOrientation(const Node & n1,
                            const Node & n2,
                            double (&rslt)[3]);
  */
  void calcFiberOrientationTensor(const FiberNetwork & fn,
				  double (&rslt)[9]);
  void calcP2(const FiberNetwork & fn, double & rslt);
  void calcAvgFiberDirection(const FiberNetwork & fn,
                             double (&rslt)[3]);
  /*
  void calcFiberDirection(const Node & n1,
                          const Node & n2,
                          double (&rslt)[3]);
  */
}
#include "bioFiberNetwork_impl.h"
#endif
