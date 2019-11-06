#include "bioRVECoupling.h"
namespace bio{
MicroscaleType getMicroscaleType(pAttribute multiscale_model)
{
  if(multiscale_model)
  {
    char * multiscale_model_type = Attribute_imageClass(multiscale_model); 
    if(strcmp(multiscale_model_type, "isotropic_neohookean") == 0)
    {
        delete [] multiscale_model_type;
        return MicroscaleType::ISOTROPIC_NEOHOOKEAN;
    }
    else if(strcmp(multiscale_model_type, "fiber only") == 0)
    {
        delete [] multiscale_model_type;
        return MicroscaleType::FIBER_ONLY;
    }
    else
    {
        std::cerr<<multiscale_model_type<<" is not a valid multiscale model type."<<std::endl;
        delete []  multiscale_model_type;
        MPI_Abort(AMSI_COMM_WORLD, 1);
    }
  }
  return MicroscaleType::NONE;
}
}
