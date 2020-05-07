#include "bioFiberNetworkLibrary.h"
#include <PCU.h>
#include <mpi.h>
#include <algorithm>
#include <istream>
#include <utility>
#include "bioFiberNetwork.h"
#include "bioApfPointers.h"
#include "bioFiberNetworkIO.h"
namespace bio
{
  FiberNetworkLibrary::FiberNetworkLibrary() : mLibrary(LibType()), mNonZeroMax(0)
  {
  }
  void FiberNetworkLibrary::load(const std::string & network_name,
                                const std::string & params_name,
                                size_t net_type,
                                size_t net_id)
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
    load(network_stream, params_stream, net_type, net_id);
  }
  void FiberNetworkLibrary::load(std::istream & network_stream,
                                 std::istream & params_stream,
                                 size_t net_type,
                                 size_t net_id)
  {
    auto network_key =  std::make_pair(net_type, net_id);
    auto libIt = mLibrary.find(network_key);
    // the Fiber network isn't in the library yet, so add it
    if (libIt == std::end(mLibrary))
    {
      auto mesh = bio::make_unique(loadFromStream(network_stream));
      FiberNetworkBase::reaction_ptr_type reactions(
          new FiberNetworkReactions(static_cast<apf::Mesh2*>(mesh.get()), params_stream));
      auto network = std::unique_ptr<FiberNetworkBase>{new FiberNetworkBase(std::move(mesh),std::move(reactions))};
      int nnz = network->getNumNonZero();
      mNonZeroMax = mNonZeroMax < nnz ? nnz : mNonZeroMax;
      mLibrary[network_key] = std::move(network);
    }
  }
  std::unique_ptr<FiberNetworkBase> FiberNetworkLibrary::getCopy(size_t net_type, size_t net_id)
  {
    auto libIt = mLibrary.find(std::make_pair(net_type, net_id));
    if (libIt != std::end(mLibrary))
    {
      auto fn = *((libIt->second).get());
      return std::unique_ptr<FiberNetworkBase>{new FiberNetworkBase(fn)};
    }
    // return a nullptr if the fiber network doesn't exist in the library
    return std::unique_ptr<FiberNetworkBase>(nullptr);
  }
  std::unique_ptr<FiberNetworkBase> FiberNetworkLibrary::getOriginalNetwork(size_t net_type, size_t net_id)
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
    return std::unique_ptr<FiberNetworkBase>(nullptr);
  }
}  // namespace bio
