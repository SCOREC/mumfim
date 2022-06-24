#ifndef MUMFIM_SRC_MUMFIM_MACROSCALE_PETSCSNES_H
#define MUMFIM_SRC_MUMFIM_MACROSCALE_PETSCSNES_H
#include <amsiLAS.h>
#include <amsiPETScLAS.h>
#include <petscsnes.h>
#include <memory>
#include "mumfim/exceptions.h"
#define MumfimPetscCall(petsc_error_code)        \
  if (petsc_error_code) [[unlikely]]             \
  {                                              \
    throw mumfim::petsc_error(petsc_error_code); \
  }
namespace mumfim
{
  struct SnesSolverResult
  {
    size_t num_iterations;
    SNESConvergedReason converged;
  };
  class SNES
  {
    public:
    explicit SNES(MPI_Comm cm) { MumfimPetscCall(SNESCreate(cm, &snes_)); }
    ~SNES() { SNESDestroy(&snes_); }
    operator ::SNES() const { return snes_; }

    private:
    ::SNES snes_;
  };
  /*
  class PetscSNES
  {
    public:
    PetscSNES(MPI_Comm cm, amsi::PetscLAS * las) : snes_(SNES(cm))
    {
      if (las == nullptr)
      {
        throw mumfim_error("Las cannot be a null ptr");
      }
      SNESGetKSP(snes_, &ksp_);
      // MumfimPetscCall(SNESSetFunction(snes, f, ComputeResidual,
      // (void*)this)); MumfimPetscCall(SNESSetJacobian(snes,JMat,JMat,
      // ComputeJacobian, (void*)this));
    }
    template <typename Residual, typename Jacobian>
    SnesSolverResult Solve(const Residual & compute_residual,
                           const Jacobian & compute_jacobian)
    {
      MumfimPetscCall(SNESSolve(snes_, nullptr, solution_vector_));
      SnesSolverResult result;
      result.converged = std::invoke(
          [this]()
          {
            SNESConvergedReason converged;
            MumfimPetscCall(SNESGetConvergedReason(snes_, &converged));
            return converged;
          });
      result.num_iterations = std::invoke(
          [this]()
          {
            PetscInt iteration;
            MumfimPetscCall(SNESGetIterationNumber(snes_, &iteration));
            return iteration;
          });
      return result;
    }
    private:
    SNES snes_;
    ::Vec solution_vector_;
    ::KSP ksp_;
  };
  */
}  // namespace mumfim
#endif  // MUMFIM_SRC_MUMFIM_MACROSCALE_PETSCSNES_H
