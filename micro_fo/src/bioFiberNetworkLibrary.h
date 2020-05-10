#ifndef BIO_FIBER_NETWORK_LIBRARY_H__
#define BIO_FIBER_NETWORK_LIBRARY_H__
#include <bioFiberNetwork.h>
#include <istream>
#include <map>
#include <memory>
#include <utility>
namespace bio
{
  class FiberNetwork;
  class FiberNetworkLibrary
  {
    public:
    void load(const std::string & network_name, const std::string& params_name, 
              size_t net_type, size_t net_id);
    void load(std::istream & mesh_istream,
              std::istream & params_istream,
              size_t net_type,
              size_t net_id);
    std::unique_ptr<FiberNetwork> getCopy(size_t net_type, size_t net_id);
    // this extracts the original network from the library
    std::unique_ptr<FiberNetwork> getOriginalNetwork(size_t net_type, size_t net_id);
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
    using LibType = typename std::map<
        std::pair<size_t, size_t>, std::unique_ptr<FiberNetwork>>;
    LibType mLibrary;
  };
}  // namespace bio
#endif
