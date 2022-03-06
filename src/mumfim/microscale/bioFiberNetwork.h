#ifndef MUMFIM_FIBER_NETWORK_H_
#define MUMFIM_FIBER_NETWORK_H_
#include <apf.h>
#include <apfFunctions.h>  // amsi
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <memory>
#include <string>
#include <vector>
#include "bioApfPointers.h"
#include "bioFiber.h"
#include "bioFiberReactions.h"
namespace mumfim
{
  class FiberNetworkReactions
  {
    public:
    FiberNetworkReactions(apf::Mesh2 * msh, std::istream & strm);
    ~FiberNetworkReactions();
    FiberReaction & operator[](size_t idx);
    //const FiberReaction & operator[](size_t idx) const;

    private:
    std::vector<FiberReaction *> mReactionsList;
  };
  // The fiber network base is an initial attempt to pull out the
  // fundamental data that the fiber network will have. The FiberNetwork class
  // should be split into FibetNetworkExplicit and FiberNetworkImplicit
  class FiberNetworkBase
  {
    public:
    using mesh_ptr_type = mumfim::mesh_unique_ptr_type;
    using reaction_ptr_type = std::shared_ptr<FiberNetworkReactions>;
    FiberNetworkBase(mesh_ptr_type mesh, reaction_ptr_type reactions);
    FiberNetworkBase(const FiberNetworkBase & other);
    // FiberNetworkBase() : fn(nullptr) {};
    virtual int getDofCount() const {return 0;}
    apf::Mesh * getNetworkMesh() const { return mMesh.get(); }
    virtual ~FiberNetworkBase();
    // TODO this should be made const, but it will take some refactor work
    // in the Truss integrator
    [[nodiscard]] const FiberReaction & getFiberReaction(size_t idx) const
    {
      return (*mReactions)[idx];
    }
    reaction_ptr_type getFiberReactions() { return mReactions; }
    /*
     * returns the type of rve (e.g. filename as an integer)
     * the mapping between this integer and the filename can
     * be found in the rve_tp log
     */
    int getRVEType() const { return mRVEType; }
    void setRVEType(int rve_t) { mRVEType = rve_t; }
    protected:
    mesh_ptr_type mMesh;
    reaction_ptr_type mReactions;
    int mRVEType;
  };
  /**
   * Responsible for managing the internal state of a single fiber-network
   * FIXME this class has a lot of state which is not needed ...
   * We don't necessarily need all of these fields in every type of analysis
   * especially for the explicit analysis I'd like to get away from having any
   * apf mesh stored since we never really want to deal with the data in that
   * format since it is incompatible with the GPU, and we should probably be
   * using an alternative output format which doesn't require the conversion.
   * Possibly Adios2 which allows us to write a vtk file as part of its binary
   * files.
   */
  class FiberNetwork : public FiberNetworkBase
  {
    protected:
    apf::Field * u;  // displacement field
    amsi::XpYFunc * xpufnc;
    apf::Field * xpu;
    apf::Field * du;
    apf::Field * v;  // velocity field
    apf::Field * a;  // acceleration field
    apf::Field * f;  // force field
    apf::Numbering * udof;
    apf::Numbering * vdof;
    apf::Numbering * adof;
    apf::Numbering * fdof;
    int ucnt;
    FiberMember tp;
    // the rve filename map
    double scale_factor;

    public:
    /**
     * Construct a FiberNetwork object.
     * @param f A pointer to a fiber network mesh (contains only vertices and
     * edges) typically loaded using the NetworkLoader classes
     */
    FiberNetwork(mesh_ptr_type mesh, reaction_ptr_type reactions);
    FiberNetwork(const FiberNetwork & fn);
    virtual ~FiberNetwork();
    /*
     * returns the scale factor for this network...
     * \note This is only relevant for a multiscale analysis, and should
     * possibly be moved to a multiscale fiber network which sublcasses from
     * here?
     */
    [[nodiscard]] double getScaleConversion() const noexcept
    {
      return scale_factor;
    }
    /*
     * sets the scale conversion factor (default 1). Typically this is
     * done in the multiscale initialization phase
     * \param sf scale factor
     */
    void setScaleConversion(double sf) noexcept { scale_factor = sf; }
    /**
     *  Gives the dimensionality of the managed fiber network
     *  @return the dimensionality of the fiber network (2 or 3)
     */
    int getDim() const { return mMesh->getDimension(); }
    /**
     * @return the number of degrees of freedom in the fiber
     * network
     */
    virtual int getDofCount() const override { return ucnt; }
    FiberMember getFiberMember() const { return tp; }
    apf::Field * getUField() const { return u; }
    apf::Field * getdUField() const { return du; }
    apf::Field * getXpUField() const { return xpu; }
    apf::Field * getVField() const { return v; }
    apf::Field * getAField() const { return a; }
    apf::Field * getFField() const { return f; }
    apf::Numbering * getUNumbering() const { return udof; }
    // velocity numbering
    apf::Numbering * getVNumbering() const { return vdof; }
    // acceleration numbering
    apf::Numbering * getANumbering() const { return adof; }
    apf::Numbering * getFNumbering() const { return fdof; }
  };
  /**
   * get the mean orientation vector of the rve
   * \param network pointer to the the fiber network
   * \param \out omega the orientation tensor
   */
  void get3DOrientationTensor(mumfim::FiberNetwork * network, double omega[9]);
  /**
   * get the mean orientation vector of the rve projected into the plane defined
   * by the normal \param network pointer to the the fiber network \param normal
   * the normal vector to the plane to project the vectors into \param \out
   * omega orientation tensor
   */
  void get2DOrientationTensor(mumfim::FiberNetwork * network,
                              double const normal[3],
                              double omega[9]);
}  // namespace mumfim
#endif
