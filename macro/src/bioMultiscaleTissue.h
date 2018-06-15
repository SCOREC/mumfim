#ifndef BIO_MULTISCALE_TISSUE_H_
#define BIO_MULTISCALE_TISSUE_H_
#include "bioNonlinearTissue.h"
#include "bioRVECoupling.h"
namespace bio
{
  // refactor so this is only dealing with a single type of RVE
  class MultiscaleTissue : public NonlinearTissue
  {
  public:
    MultiscaleTissue(pGModel g, pParMesh m, pACase pd, MPI_Comm cm);
    ~MultiscaleTissue();
    virtual void Assemble(amsi::LAS * las);
    void computeRVEs();
    void initMicro();
    void updateMicro();
    void updateRVETypes();
    int updateRVEType(apf::MeshEntity * me);
    void updateRVEExistence();
    virtual void recoverSecondaryVariables(int step);
  private:
    amsi::ElementalSystem * mltscl;
    apf::Field * crt_rve;
    apf::Field * prv_rve;
    bool compute_ornt_3D;
    bool compute_ornt_2D;
    double ornt_2D_axis[3];
    apf::Field * ornt_3D;
    apf::Field * ornt_2D;
    RVECoupling fo_cplg;
    int nm_rves;
    //int nm_rve_rgns;
    // multiscale coupling communication stuff
    enum PATTERN
    {
      SEND = 0,
      RECV = 1,
      SEND_INIT = 2,
      NUM_PATTERNS = 3
    };
    //size_t rve_ptrns[NUM_PATTERNS];
    //size_t snd_ptrns[MICROSCALE_TYPE_COUNT];
    //size_t ini_ptrns[MICROSCALE_TYPE_COUNT];
    //size_t rcv_ptrns[MICROSCALE_TYPE_COUNT];
    size_t M2m_id;
    size_t m2M_id;
    std::vector<std::string> rve_dirs;
    std::vector<int> rve_dir_cnts;
    amsi::ElementalSystem * getIntegrator(apf::MeshEntity * me, int ii);
    template <typename O>
      void serializeRVEData(O o);
    template <typename O1, typename O2, typename O3>
      void serializeNewRVEData(O1 nw_hdrs, O2 nw_prms, O3 nw_data, bool all = false);
    void getExternalRVEData(apf::MeshEntity * ent,
                            micro_fo_header & hdr,
                            micro_fo_params & prms);
    void getInternalRVEData(apf::MeshEntity * ent,
                            micro_fo_header & hdr,
                            micro_fo_params & prms,
                            micro_fo_init_data & dat);
    void loadRVELibraryInfo();
  };
  /**
   * Get the index in [tp_bgn, tp_end] for the RVE directory specified in the simmetrix
   *  simmodeler data for the ModelEntity ent.
   */
  template <typename I>
    int getRVEDirectoryIndex(I tp_bgn, I tp_end, apf::ModelEntity * ent);

}
#include "bioMultiscaleTissue_impl.h"
#endif
