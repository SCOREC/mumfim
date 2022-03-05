#ifndef BIO_MULTISCALE_CONVERGENCE_H_
#define BIO_MULTISCALE_CONVERGENCE_H_
namespace bio
{
  // a multi-convergence wrapper that broadcasts to the
  // micro-scale if the rest of the convergence operations
  // are true. could implement as a standalone convergence
  // and provide it as the final operation to a normal
  // multi-convergence insteadr
  class MultiscaleConvergence : public amsi::MultiConvergence
  {
  private:
    amsi::ControlService * cs;
    size_t cplg;
  public:
    template <typename I>
      MultiscaleConvergence(I bgn, I end, size_t c)
        : amsi::MultiConvergence(bgn,end)
      , cs(amsi::ControlService::Instance())
      , cplg(c)
    { }
    virtual bool converged()
    {
      bool rslt = MultiConvergence::converged();
      cs->scaleBroadcast(cplg,&rslt);
      return rslt;
    }
  };
}
#endif
