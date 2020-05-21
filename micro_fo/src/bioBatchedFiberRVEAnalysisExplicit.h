#ifndef BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_H__
#define BIO_BATCHED_FIBER_RVE_ANALYSIS_EXPLICIT_H__
#include "bioBatchedRVEAnalysis.h"
#include "bioMicroTypeDefinitions.h"


namespace bio
{
template <typename Scalar, typename LocalOrdinal, typename ExeSpace=Kokkos::DefaultExecutionSpace>
class BatchedFiberRVEAnalysisExplicit : public BatchedRVEAnalysis<Scalar>
{
  private:
    using LO = LocalOrdinal;
  public:
    virtual bool run(const std::vector<DeformationGradient> & dfmGrds,
                     std::vector<Scalar[6]> sigma, bool update_coords=true) final;
    virtual void computeMaterialStiffness(std::vector<Scalar[36]> C) final;
};

// the actual explicit instantionations can be found in the associated cc file
extern template class BatchedFiberRVEAnalysisExplicit<Scalar, LocalOrdinal>;

// implementation file for templated free functions
#include "bioBatchedFiberRVEAnalysisExplicit_impl.h"

// templated function definitions
template <typename Scalar, typename LocalOrdinal, typename ExeSpace>
bool BatchedFiberRVEAnalysisExplicit<Scalar, LocalOrdinal, ExeSpace>::run(const std::vector<DeformationGradient> & dfmGrds,std::vector<Scalar[6]> sigma, bool update_coords)
{
  impl::getLinearReactionForce<double>(1,1.5,4,10);
  return false;
};
template <typename Scalar, typename LocalOrdinal, typename ExeSpace>
void BatchedFiberRVEAnalysisExplicit<Scalar, LocalOrdinal, ExeSpace>::computeMaterialStiffness(std::vector<Scalar[36]> C){};

}
#endif
