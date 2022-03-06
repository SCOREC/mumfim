#ifndef MUMFIM_FIBER_NETWORK_LIBRARY_H
#define MUMFIM_FIBER_NETWORK_LIBRARY_H
#include <mumfim/microscale/bioFiberNetwork.h>
#include <istream>
#include <map>
#include <memory>
#include <utility>
namespace mumfim
{
  class FiberNetwork;
  class FiberNetworkLibrary
  {
    public:
    std::shared_ptr<FiberNetwork> load(const std::string & network_name,
                                       const std::string & params_name,
                                       std::size_t net_type,
                                       std::size_t net_id);
    std::shared_ptr<FiberNetwork> load(std::istream & mesh_istream,
                                       std::istream & params_istream,
                                       std::size_t net_type,
                                       std::size_t net_id);
    [[nodiscard]] std::unique_ptr<FiberNetwork> getUniqueCopy(
        std::size_t net_type,
        std::size_t net_id) const;
    [[nodiscard]] std::shared_ptr<FiberNetwork> getSharedCopy(
        std::size_t net_type,
        std::size_t net_id) const;
    // this extracts the original network from the library
    [[nodiscard]] std::shared_ptr<FiberNetwork> extractNetwork(
        std::size_t net_type,
        std::size_t net_id);
    FiberNetworkLibrary();
    ~FiberNetworkLibrary();
    private:
    // stores a map with the key as a pair of net type, and net_id
    // and returns a pair of a unique pointer to a mesh, and reaction for
    // that fiber network
    // TODO this map type can be replaced with a unordered map type to bring
    // the complexity from logarithmic to constant for the search. Note,
    // typically we will have a small container, so th cost shouldn't be
    // terrible. transitioning to an unordered map requires writing a hash
    // function for the pair.
    using LibType = typename std::map<std::pair<std::size_t, std::size_t>,
                                      std::shared_ptr<FiberNetwork>>;
    LibType mLibrary;
  };
}  // namespace mumfim
#endif
