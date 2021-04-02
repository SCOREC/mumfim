#ifndef BIO_MULTISCALE_TISSUE_H_
#define BIO_MULTISCALE_TISSUE_H_
#include <memory>
#include <unordered_map>
#include "bioMultiscaleTissue.h"
#include "bioNonlinearTissue.h"
#include "bioRVECoupling.h"
#include "bioReadStochasticField.h"
namespace bio
{
  using StochasticFieldMap = std::map<std::string, std::shared_ptr<GridData> >;
  // refactor so this is only dealing with a single type of RVE
  class MultiscaleTissue : public NonlinearTissue
  {
    public:
    MultiscaleTissue(apf::Mesh * mesh,
                     const mt::CategoryNode & analysis_case,
                     MPI_Comm cm);
    ~MultiscaleTissue();
    virtual void Assemble(amsi::LAS * las) override;
    void computeRVEs();
    void initMicro();
    void updateMicro();
    void updateRVETypes();
    MicroscaleType updateRVEType(apf::MeshEntity * me);
    void updateRVEExistence();
    virtual void recoverSecondaryVariables(int step) override;
    virtual void preRun() override { updateMicro(); }

    private:
    amsi::ElementalSystem * mltscl;
    apf::Field * crt_rve;
    apf::Field * prv_rve;
    bool compute_ornt_3D;
    bool compute_ornt_2D;
    double ornt_2D_axis[3];
    apf::Field * ornt_3D;
    apf::Field * ornt_2D;
    // DEBUG
    apf::Field * test_inc_dfm;
    // END DEBUG
    RVECoupling fo_cplg;
    int nm_rves;
    StochasticFieldMap stochastic_field_map;
    // int nm_rve_rgns;
    // multiscale coupling communication stuff
    enum PATTERN
    {
      SEND = 0,
      RECV = 1,
      SEND_INIT = 2,
      NUM_PATTERNS = 3
    };
    // size_t rve_ptrns[NUM_PATTERNS];
    // size_t snd_ptrns[MICROSCALE_TYPE_COUNT];
    // size_t ini_ptrns[MICROSCALE_TYPE_COUNT];
    // size_t rcv_ptrns[MICROSCALE_TYPE_COUNT];
    size_t M2m_id;
    size_t m2M_id;
    std::vector<std::string> rve_dirs;
    std::vector<int> rve_dir_cnts;
    amsi::ElementalSystem * getIntegrator(apf::MeshEntity * me, int ii);
    template <typename O>
    void serializeRVEData(O o);
    template <typename O1, typename O2, typename O3, typename O4, typename O5>
    void serializeNewRVEData(O1 nw_hdrs,
                             O2 nw_prms,
                             O3 nw_data,
                             O4 nw_slvr_prms,
                             O5 nw_int_slvr_prms,
                             bool all = false);
    void getExternalRVEData(apf::MeshEntity * ent,
                            micro_fo_header & hdr,
                            micro_fo_params & prms,
                            micro_fo_solver & slvr,
                            micro_fo_int_solver & int_slvr);
    void getInternalRVEData(apf::MeshEntity * ent,
                            micro_fo_header & hdr,
                            micro_fo_params & prms,
                            micro_fo_init_data & dat);
    void loadRVELibraryInfo();
    /**
     * Get the index in the RVE directory specified in the
     * simmetrix simmodeler data for the ModelEntity ent.
     */
    int getRVEDirectoryIndex(apf::MeshEntity* rgn);
  };
}  // namespace bio
#include "bioMultiscaleTissue_impl.h"
#endif
