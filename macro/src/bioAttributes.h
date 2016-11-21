#ifndef BIO_ATTRIBUTES_H_
#define BIO_ATTRIBUTES_H_
#include <simAttributes.h>
namespace bio
{
  // the structure should be defined in micro_fo probably
  struct FiberNetworkInfo
  {
    int elmt_tp;
    int fbr_rctn;
    double fbr_rds;
    double vol_frc;
  };
  void getMultiscaleInfo(pGEntity ent, FiberNetworkInfo & fn_tp)
  {
    pAttribute mltscl = GEN_attrib(ent,"multiscale");
    if(mltscl != NULL)
    {
      fn_tp.elmt_tp = 0; // truss
      fn_tp.fbr_rctn = AttributeInt_value((pAttributeInt)Attribute_childByType(mltscl,"force_reaction"));
      fn_tp.fbr_rds = AttributeTensor0_value((pAttributeTensor0)Attribute_childByType(mltscl,"radius"));
      fn_tp.vol_frc = AttributeTensor0_value((pAttributeTensor0)Attribute_childByType(mltscl,"force reaction"));
    }
  }
  template <typename O>
    void updateMultiscaleAssignment(apf::Mesh * msh, apf::Field * micro_type, O add, O delete)
  {
    apf::MeshEntity * me = NULL;
    for(apf::MeshIterator * it = msh->begin(msh->getDimension()), (me = msh->iterate(it));)
    {
      int nw_tp = Element_UpdateMicroscaleType(me);
      int ld_tp = apf::getScalar(micro_type,me,0);
      if(ld_tp == NONE && nw_tp != NONE)
        *add++ = me;
      else if(ld_tp != NONE && nw_tp == NONE)
        *delete++ = me;
      apf::setScalar(micro_type,me,0,nw_tp);
    }
  }
  // during assembly
  template <typename I, typename O>
    void buildAddTypes(I begin, I end, O out)
  {
    for(I it = begin; it != end; ++it)
    {
      apf::MeshEntity * ent = *it;
      pEntity sim_ent = (pEntity)ent;
      pGEntity gent = EN_whatIn(sim_ent);
      *out++ = getMultiscaleInfo(gent);
    }
  }
  //assert(EN_whatInType(msh_ent) == Gregion);
  //pGEntity ent = EN_whatIn(msh_ent);
  //FiberNetworkInfo * fn_tp = makeMultiscaleInfo(ent);
}
#endif
