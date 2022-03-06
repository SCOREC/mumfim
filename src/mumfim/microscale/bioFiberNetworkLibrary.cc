#include "bioFiberNetworkLibrary.h"
#include <PCU.h>
#include <mpi.h>
#include <algorithm>
#include <istream>
#include <utility>
#include "bioFiberNetwork.h"
#include "bioApfPointers.h"
#include "bioFiberNetworkIO.h"
namespace mumfim
{
  FiberNetworkLibrary::FiberNetworkLibrary() : mLibrary(LibType())
  {
  }
  FiberNetworkLibrary::~FiberNetworkLibrary(){};
  std::shared_ptr<FiberNetwork> FiberNetworkLibrary::load(
      const std::string & network_name,
      const std::string & params_name,
      std::size_t net_type,
      std::size_t net_id)
  {
    std::ifstream network_stream(network_name);
    if (!network_stream.is_open())
    {
      std::cerr << "Cannot open the microscale network file " << network_name
                << " aborting analysis" << std::endl;
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    std::ifstream params_stream(params_name);
    if (!network_stream.is_open())
    {
      std::cerr << "Cannot open the microscale parameters file " << params_name
                << " aborting analysis" << std::endl;
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    return load(network_stream, params_stream, net_type, net_id);
  }
  std::shared_ptr<FiberNetwork> FiberNetworkLibrary::load(
      std::istream & network_stream,
      std::istream & params_stream,
      std::size_t net_type,
      std::size_t net_id)
  {
    auto network_key =  std::make_pair(net_type, net_id);
    auto libIt = mLibrary.find(network_key);
    // the Fiber network isn't in the library yet, so add it
    if (libIt == std::end(mLibrary))
    {
      auto mesh = mumfim::make_unique(loadFromStream(network_stream));
      FiberNetwork::reaction_ptr_type reactions(
          new FiberNetworkReactions(static_cast<apf::Mesh2*>(mesh.get()), params_stream));
      auto network = std::shared_ptr<FiberNetwork>{
          new FiberNetwork(std::move(mesh), std::move(reactions))};
      mLibrary[network_key] = network;
      return network;
    }
    return libIt->second;
  }
  std::unique_ptr<FiberNetwork> FiberNetworkLibrary::getUniqueCopy(
      std::size_t net_type,
      std::size_t net_id) const
  {
    auto libIt = mLibrary.find(std::make_pair(net_type, net_id));
    if (libIt != std::end(mLibrary))
    {
      auto fn = *((libIt->second).get());
      return std::unique_ptr<FiberNetwork>{new FiberNetwork(fn)};
    }
    // return a nullptr if the fiber network doesn't exist in the library
    return std::unique_ptr<FiberNetwork>(nullptr);
  }
  std::shared_ptr<FiberNetwork> FiberNetworkLibrary::getSharedCopy(
      std::size_t net_type,
      std::size_t net_id) const
  {
    auto libIt = mLibrary.find(std::make_pair(net_type, net_id));
    if (libIt != std::end(mLibrary))
    {
      return libIt->second;
    }
    // return a nullptr if the fiber network doesn't exist in the library
    return std::unique_ptr<FiberNetwork>(nullptr);
  }
  std::shared_ptr<FiberNetwork> FiberNetworkLibrary::extractNetwork(
      std::size_t net_type,
      std::size_t net_id)
  {
    auto libIt = mLibrary.find(std::make_pair(net_type, net_id));
    if (libIt != std::end(mLibrary))
    {
      auto network = std::move(libIt->second);
      // clear out the item in the library to
      // make sure nobody uses the moved from network
      mLibrary.erase(libIt);
      return network;
    }
    // return a nullptr if the fiber network doesn't exist in the library
    return std::unique_ptr<FiberNetwork>(nullptr);
  }
}  // namespace mumfim
