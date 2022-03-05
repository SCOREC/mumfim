#ifndef BIO_RVE_COUPLING_H_
#define BIO_RVE_COUPLING_H_
#include <bioMultiscaleMicroFOParams.h>
#include <apf.h>
#include <amsiMultiscale.h>
#include <apfMesh.h>
#include <model_traits/AssociatedModelTraits.h>
namespace bio
{
  MicroscaleType getMicroscaleType(
      const mt::AssociatedCategoryNode * category_node);

  class RVECoupling
  {
  private:
    apf::Mesh * msh;
    apf::Field * rst_fst;
    apf::Field * crt_rve;
    apf::Field * prv_rve;
    int fld_ord;
    int rve_cnt;
    std::vector<micro_fo_result> rsts;
    std::vector<micro_fo_step_result> stp_rslt;
    size_t snd_ptrn;
    size_t rcv_ptrn;
  public:
  RVECoupling(apf::Mesh * m, apf::Field * crt, apf::Field * old, int o)
      : msh(m)
      , rst_fst(NULL)
      , crt_rve(crt)
      , prv_rve(old)
      , fld_ord(o)
      , rve_cnt(0)
      , rsts()
      , stp_rslt()
      , snd_ptrn()
      , rcv_ptrn()
    {
      rst_fst = apf::createIPField(msh,"micro_fo_rve_offset",apf::SCALAR,1);
      int ip_id = 0;
      apf::MeshIterator * it = m->begin(m->getDimension());
      // assuming 1 integration point per element
      for(apf::MeshEntity * me = NULL; (me = m->iterate(it));)
      {
        apf::MeshElement * mlmt = apf::createMeshElement(m,me);
        int nip = apf::countIntPoints(mlmt,0);
        apf::destroyMeshElement(mlmt);
        for(int ip = 0; ip < nip; ++ip)
          apf::setScalar(rst_fst,me,ip,ip_id++);
      }
      m->end(it);
    }
    ~RVECoupling()
    {
      apf::destroyField(rst_fst);
    }
    void initCoupling()
    {
      amsi::ControlService * cs = amsi::ControlService::Instance();
      amsi::Task * macro = amsi::getLocal();
      amsi::DataDistribution * dd = amsi::createDataDistribution(macro,"micro_fo_data");
      (*dd) = 0;
      amsi::Assemble(dd,macro->comm());
      snd_ptrn = cs->CreateCommPattern("micro_fo_data","macro","micro_fo");
      cs->CommPattern_Reconcile(snd_ptrn);
      rcv_ptrn = cs->RecvCommPattern("macro_fo_data","micro_fo","micro_fo_results","macro");
      cs->CommPattern_Reconcile(rcv_ptrn);
    }
    template <typename I>
      void deleteRVEs(I dlt_bgn, I dlt_end)
    {
      amsi::ControlService * cs = amsi::ControlService::Instance();
      std::vector<int> to_dlt;
      std::copy(dlt_bgn,dlt_end,std::back_inserter(to_dlt));
      cs->RemoveData(snd_ptrn,to_dlt);
    }
    template <typename I1, typename I2, typename I3, typename I4, typename I5>
      void sendNewRVEs(size_t ptrn, I1 hdr, I2 prm, I3 dat, I4 micro_slvr, I5 micro_int_slvr)
    {
      amsi::ControlService * cs = amsi::ControlService::Instance();
      cs->Communicate(ptrn,hdr,amsi::mpi_type<bio::micro_fo_header>());
      cs->Communicate(ptrn,prm,amsi::mpi_type<bio::micro_fo_params>());
      cs->Communicate(ptrn,dat,amsi::mpi_type<bio::micro_fo_init_data>());
      cs->Communicate(ptrn,micro_slvr,amsi::mpi_type<bio::micro_fo_solver>());
      cs->Communicate(ptrn,micro_int_slvr,amsi::mpi_type<bio::micro_fo_int_solver>());
    }
    size_t addRVEs(int sz)
    {
      amsi::ControlService * cs = amsi::ControlService::Instance();
      size_t cnt = sz;
      size_t add_id = cs->addData(snd_ptrn,cnt);
      rve_cnt += sz;
      return add_id;
    }
    void updateRecv()
    {
      amsi::ControlService::Instance()->CommPattern_Reconcile(rcv_ptrn);
    }
    // todo (h) : template or pointer
    void sendRVEData(std::vector<micro_fo_data> & bfr)
    {
      amsi::ControlService * cs = amsi::ControlService::Instance();
      cs->Communicate(snd_ptrn,bfr,amsi::mpi_type<bio::micro_fo_data>());
    }
    void recvRVEData()
    {
      amsi::ControlService * cs = amsi::ControlService::Instance();
      cs->Communicate(rcv_ptrn,rsts, amsi::mpi_type<bio::micro_fo_result>());
    }
    micro_fo_result * getRVEResult(apf::MeshEntity * me, int ip)
    {
      int fst = apf::getScalar(rst_fst,me,ip);
      return &rsts[fst];
    }
    void recvRVEStepData()
    {
      amsi::ControlService * cs = amsi::ControlService::Instance();
      cs->Communicate(rcv_ptrn,stp_rslt,amsi::mpi_type<bio::micro_fo_step_result>());

    }
    micro_fo_step_result * getRVEStepResult(apf::MeshEntity * me, int ip) {
      int fst = apf::getScalar(rst_fst, me, ip);
      return &stp_rslt[fst];
    }
    int countRVEsOn(apf::MeshEntity * me)
    {
      int cnt = 0;
      apf::MeshElement * mlm = apf::createMeshElement(msh,me);
      int ng = apf::countIntPoints(mlm,fld_ord);
      for(int ip = 0; ip < ng; ++ip)
        if(static_cast<MicroscaleType>(apf::getScalar(crt_rve,me,ip)) == 
           MicroscaleType::FIBER_ONLY)
          ++cnt;
      apf::destroyMeshElement(mlm);
      return cnt;
    }
    template <typename O>
      void updateRVEDeletion(O out, bool all = false)
    {
      int iid = 0;
      apf::MeshEntity * rgn = NULL;
      for(apf::MeshIterator * it = msh->begin(3); (rgn = msh->iterate(it)); )
      {
        apf::MeshElement * mlm = apf::createMeshElement(msh,rgn);
        int ng = apf::countIntPoints(mlm,getOrder(mlm));
        for(int ip = 0; ip < ng; ++ip)
        {
          MicroscaleType crt = static_cast<MicroscaleType>(apf::getScalar(crt_rve,rgn,ip));
          MicroscaleType prv = static_cast<MicroscaleType>(apf::getScalar(prv_rve,rgn,ip));
          if((crt == MicroscaleType::NONE && prv != MicroscaleType::NONE) || all)
            out++ = iid; // hacky id method, why?
          iid++;
        }
        apf::destroyMeshElement(mlm);
      }
    }
  };
}
#endif
