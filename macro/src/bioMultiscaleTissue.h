#ifndef BIO_MULTISCALE_TISSUE_H_
#define BIO_MULTISCALE_TISSUE_H_
#include "bioNonlinearTissue.h"
namespace bio
{
  // refactor so this is only dealing with a single type of RVE
  class MultiscaleTissue : public NonlinearTissue
  {
  public:
    enum MicroscaleType{NONE = 0,
                        FIBER_ONLY = 1,
                        FIBER_MATRIX = 2,
                        MICROSCALE_TYPE_COUNT = 3};
  private:
    amsi::ElementalSystem * mltscl;
    apf::Field * crt_rve;
    apf::Field * prv_rve;
    apf::Field * fbr_ornt;
    int nm_rves;
    int nm_rve_rgns;
    std::list<apf::MeshEntity *> rve_ents;
    std::vector<micro_fo_result> rslts;
    std::map<apf::MeshEntity*,std::vector<RVE_Info> > rslt_mp;
    // multiscale coupling communication stuff
    enum PATTERN
    {
      SEND = 0,
      RECV = 1,
      SEND_INIT = 2,
      NUM_PATTERNS = 3
    };
    size_t rve_ptrns[NUM_PATTERNS];
    size_t snd_ptrns[MICROSCALE_TYPE_COUNT];
    size_t ini_ptrns[MICROSCALE_TYPE_COUNT];
    size_t rcv_ptrns[MICROSCALE_TYPE_COUNT];
    size_t M2m_id;
    size_t m2M_id;
    std::map<std::string,int> rve_tps;
    MicroFOMultiscaleDataTypes mtd;
  public:
    MultiscaleTissue(pGModel g, pParMesh m, pACase pd, MPI_Comm cm);
    ~MultiscaleTissue();
    virtual void Assemble(amsi::LAS * las);
    void computeRVEs();
    void initMicro();
    int countRVEsOn(apf::MeshEntity * me);
    void updateMicro();
    void updateRVETypes();
    int updateRVEType(apf::MeshEntity * me);
    void updateRVEExistence();
    // see if we can replace RVE_Info with fiber_only_result or something
    RVE_Info * getRVEResult(apf::MeshEntity * me, int ip);
    template <typename O>
      void updateRVEDeletion(O o, bool all = false);
    template <typename O1, typename O2, typename O3, typename O4>
      void updateRVEAddition(O1 nw_ents,
                             O2 nw_hdrs,
                             O3 nw_prms,
                             O4 nw_data,
                             bool all = false);
  private:
    // private member functions
    amsi::ElementalSystem * getIntegrator(apf::MeshEntity * me, int ii);
    template <typename O>
      void serializeRVEData(O o);
    void computeRVETypeInfo();
    int getRVEType(apf::ModelEntity *);
  };
}
#include "bioMultiscaleTissue_impl.h"
#endif
