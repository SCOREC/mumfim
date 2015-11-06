c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c         Basic Iterative Solvers with Reverse Communication           c
c----------------------------------------------------------------------c
c     This file currently has several basic iterative linear system    c
c     solvers. They are:                                               c
c     CG       -- Conjugate Gradient Method                            c
c     CGNR     -- Conjugate Gradient Method (Normal Residual equation) c
c     BCG      -- Bi-Conjugate Gradient Method                         c
c     DBCG     -- BCG with partial pivoting                            c
c     BCGSTAB  -- BCG stabilized                                        c
c     TFQMR    -- Transpose-Free Quasi-Minimum Residual method         c
c     GMRES    -- Generalized Minimum RESidual method                  c
c     FGMRES   -- Flexible version of Generalized Minimum              c
c                 RESidual method                                      c
c     DQGMRES  -- Direct versions of Quasi Generalize Minimum           c
c                 Residual method                                      c
c----------------------------------------------------------------------c
c     They all have the following calling sequence:
c      subroutine solver(n, rhs, sol, ipar, fpar, w)
c      integer n, ipar(16)
c      real*8 rhs(n), sol(n), fpar(16), w(*)
c     Where
c     (1) 'n' is the size of the linear system,
c     (2) 'rhs' is the right-hand side of the linear system,
c     (3) 'sol' is the solution to the linear system,
c     (4) 'ipar' is an integer parameter array for the reverse
c     communication protocol,
c     (5) 'fpar' is an floating-point parameter array storing
c     information to and from the iterative solvers.
c     (6) 'w' is the work space (size is specified in ipar)
c
c     They are preconditioned iterative solvers with reverse
c     communication. The preconditioners can be applied from either
c     from left or right or both (specified by ipar(2), see below).
c
c     NOTE:
c
c     (1) Work space required by each of the iterative solver
c     routines is as follows:
c       CG      == 5 * n
c       CGNR    == 5 * n
c       BCG     == 7 * n
c       DBCG    == 11 * n
c       BCGSTAB == 9 * n
c       TFQMR   == 11 * n
c       GMRES   == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
c       FGMRES  == n*(2m+1) + (m+1)*m/2 + 3*m + 2 (m = ipar(5),
c                  default m=15)
c       DQGMRES == n + lb * (2*n+4) (lb=ipar(5)+1, default lb = 16)
c
c     (2) ALL iterative solvers require a user-supplied DOT-product
c     routine named DISTDOT. The prototype of DISTDOT is
c
c     real*8 function distdot(n,x,ix,y,iy)
c     integer n, ix, iy
c     real*8 x(1+(n-1)*ix), y(1+(n-1)*iy)
c
c     This interface of DISTDOT is exactly the same as that of
c     DDOT (or SDOT if real == real*8) from BLAS-1. It should have
c     same functionality as DDOT on a single processor machine. On a
c     parallel/distributed environment, each processor can perform
c     DDOT on the data it has, then perform a summation on all the
c     partial results.
c
c     (3) To use this set of routines under SPMD/MIMD program paradigm,
c     several things are to be noted: (a) 'n' should be the number of
c     vector elements of 'rhs' that is present on the local processor.
c     (b) if RHS(i) is on processor j, it is expected that SOL(i)
c     will be on the same processor, i.e. the vectors are distributed
c     to each processor in the same way. (c) the preconditioning and
c     stopping criteria specifications have to be the same on all
c     processor involved, ipar and fpar have to be the same on each
c     processor. (d) DISTDOT should be replaced by a distributed
c     dot-product function.
c
c     ..................................................................
c     Reverse Communication Protocols
c
c     When a reverse-communication routine returns, it could be either
c     that the routine has terminated or it simply requires the caller
c     to perform one matrix-vector multiplication. The possible matrices
c     that involve in the matrix-vector multiplications are:
c     A       (the matrix of the linear system),
c     A^T     (A transposed),
c     Ml^{-1} (inverse of the left preconditioner),
c     Ml^{-T} (inverse of the left preconditioner transposed),
c     Mr^{-1} (inverse of the right preconditioner),
c     Mr^{-T} (inverse of the right preconditioner transposed).
c     For all the matrix vector multiplication, v = A u. The input and
c     output vectors are supposed to be part of the work space 'w', and
c     the starting positions of them are stored in ipar(8:9), see below.
c
c     The array 'ipar' is used to store the information about the solver.
c     Here is the list of what each element represents:
c
c     ipar(1) -- status of the call/return.
c     A call to the solver with ipar(1) == 0 will initialize the
c     iterative solver. On return from the iterative solver, ipar(1)
c     carries the status flag which indicates the condition of the
c     return. The status information is divided into two categories,
c     (1) a positive value indicates the solver requires a matrix-vector
c     multiplication,
c     (2) a non-positive value indicates termination of the solver.
c     Here is the current definition:
c       1 == request a matvec with A,
c       2 == request a matvec with A^T,
c       3 == request a left preconditioner solve (Ml^{-1}),
c       4 == request a left preconditioner transposed solve (Ml^{-T}),
c       5 == request a right preconditioner solve (Mr^{-1}),
c       6 == request a right preconditioner transposed solve (Mr^{-T}),
c      10 == request the caller to perform stopping test,
c       0 == normal termination of the solver, satisfied the stopping
c            criteria,
c      -1 == termination because iteration number is greater than the
c            preset limit,
c      -2 == return due to insufficient work space,
c      -3 == return due to anticipated break-down / divide by zero,
c     -10 == return due to some non-numberical reasons, e.g. invalid
c            floating-point numbers etc.
c
c     ipar(2) -- status of the preconditioning:
c       0 == no preconditioning
c       1 == left preconditioning only
c       2 == right preconditioning only
c       3 == both left and right preconditioning
c
c     ipar(3) -- stopping criteria (details of this will be
c     discussed later).
c
c     ipar(4) -- number of elements in the array 'w'. if this is less
c     than the desired size, it will be over-written with the minimum
c     requirement. In which case the status flag ipar(1) = -2.
c
c     ipar(5) -- size of the Krylov subspace (used by GMRES and its
c     variants), e.g. GMRES(ipar(5)), FGMRES(ipar(5)),
c     DQGMRES(ipar(5)).
c
c     ipar(6) -- maximum number of matrix-vector multiplies, if not a
c     positive number the iterative solver will run till convergence
c     test is satisfied.
c
c     ipar(7) -- current number of matrix-vector multiplies. It is
c     incremented after each matrix-vector multiplication. If there
c     is preconditioning, the counter is incremented after the
c     preconditioning associated with each matrix-vector multiplication.
c
c     ipar(8) -- pointer to the input vector to the requested matrix-
c     vector multiplication.
c
c     ipar(9) -- pointer to the output vector of the requested matrix-
c     vector multiplication.
c
c     To perform v = A * u, it is assumed that u is w(ipar(8):ipar(8)+n-1)
c     and v is stored as w(ipar(9):ipar(9)+n-1).
c
c     ipar(10) -- the return address (used to determine where to go to
c     inside the iterative solvers after the caller has performed the
c     requested services).
c
c     ipar(11) -- the result of the external convergence test
c     On final return from the iterative solvers, this value
c     will be reflected by ipar(1) = 0 (details discussed later)
c
c     ipar(12) to ipar(16) are NOT defined, they are NOT USED by
c     any iterative solver at this time.
c
c     Information about the error and tolerance are stored in the array
c     FPAR. So are some internal variables that need to be saved from
c     one iteration to the next one. Since the internal variables are
c     not the same for each routine, we only define the common ones.
c
c     The first two are input parameters:
c     fpar(1) -- the relative tolerance,
c     fpar(2) -- the absolute tolerance (details discussed later),
c
c     When the iterative solver terminates,
c     fpar(3) -- initial residual/error norm,
c     fpar(4) -- target residual/error norm,
c     fpar(5) -- current residual norm (if available),
c     fpar(6) -- current residual/error norm,
c     fpar(7) -- convergence rate,
c
c     fpar(8:10) are used by some of the iterative solvers to save some
c     internal information.
c
c     fpar(11) -- number of floating-point operations. The iterative
c     solvers will add the number of FLOPS they used to this variable,
c     but they do NOT initialize it, nor add the number of FLOPS due to
c     matrix-vector multiplications (since matvec is outside of the
c     iterative solvers). To insure the correct FLOPS count, the
c     caller should set fpar(11) = 0 before invoking the iterative
c     solvers and account for the number of FLOPS from matrix-vector
c     multiplications and preconditioners.
c
c     fpar(12:16) are not used in current implementation.
c
c     Whether the content of fpar(3), fpar(4) and fpar(6) are residual
c     norms or error norms depends on ipar(3). If the requested
c     convergence test is based on the residual norm, they will be
c     residual norms. If the caller want to test convergence based the
c     error norms (estimated by the norm of the modifications applied
c     to the approximate solution), they will be error norms.
c     Convergence rate is defined by (Fortran 77 statement)
c     fpar(7) = log10(fpar(3) / fpar(6)) / ipar(7)
c     If fpar(7) = 0.5, it means that approximately every 2 (= 1/0.5)
c     steps the residual/error norm decrease by a factor of 10.
c
c     ..................................................................
c     Stopping criteria,
c
c     An iterative solver may be terminated due to (1) satisfying
c     convergence test; (2) exceeding iteration limit; (3) insufficient
c     work space; (4) break-down. Checking of the work space is
c     only done in the initialization stage, i.e. when it is called with
c     ipar(1) == 0. A complete convergence test is done after each
c     update of the solutions. Other conditions are monitored
c     continuously.
c
c     With regard to the number of iteration, when ipar(6) is positive,
c     the current iteration number will be checked against it. If
c     current iteration number is greater the ipar(6) than the solver
c     will return with status -1. If ipar(6) is not positive, the
c     iteration will continue until convergence test is satisfied.
c
c     Two things may be used in the convergence tests, one is the
c     residual 2-norm, the other one is 2-norm of the change in the
c     approximate solution. The residual and the change in approximate
c     solution are from the preconditioned system (if preconditioning
c     is applied). The DQGMRES and TFQMR use two estimates for the
c     residual norms. The estimates are not accurate, but they are
c     acceptable in most of the cases. Generally speaking, the error
c     of the TFQMR's estimate is less accurate.
c
c     The convergence test type is indicated by ipar(3). There are four
c     type convergence tests: (1) tests based on the residual norm;
c     (2) tests based on change in approximate solution; (3) caller
c     does not care, the solver choose one from above two on its own;
c     (4) caller will perform the test, the solver should simply continue.
c     Here is the complete definition:
c      -2 == || dx(i) || <= rtol * || rhs || + atol
c      -1 == || dx(i) || <= rtol * || dx(1) || + atol
c       0 == solver choosing the test
c       1 == || residual || <= rtol * || initial residual || + atol
c       2 == || residual || <= rtol * || rhs || + atol
c     999 == caller will perform the test
c     where dx(i) denote the change in the solution at the ith update.
c     ||.|| denotes 2-norm. rtol = fpar(1) and atol = fpar(2).
c
c     If the caller is to perform the convergence test, the outcome
c     should be stored in ipar(11).
c     ipar(11) = 0 -- failed the convergence test, iterative solver
c     should continue
c     ipar(11) = 1 -- satisfied convergence test, iterative solver
c     should perform the clean up job and stop.
c
c     NOTE: the caller should allow the iterative solver to perform
c     clean up job after the external convergence test is satisfied,
c     since some of the iterative solvers do not directly
c     update the 'sol' array. A typical clean-up stage includes
c     performing the final update of the approximate solution and
c     computing the convergence information (e.g. values of fpar(3:7)).
c
c     ..................................................................
c     Usage:
c
c     To start solving a linear system, the user needs to specify
c     first 6 elements of the ipar, and first 2 elements of fpar.
c     The user may optionally set fpar(11) = 0 if one wants to count
c     the number of floating-point operations. (Note: the iterative
c     solvers will only add the floating-point operations inside
c     themselves, the caller will have to add the FLOPS from the
c     matrix-vector multiplication routines and the preconditioning
c     routines in order to account for all the arithmetic operations.)
c
c     Here is an example:
c     ipar(1) = 0	! always 0 to start an iterative solver
c     ipar(2) = 2	! right preconditioning
c     ipar(3) = 1	! use convergence test scheme 1
c     ipar(4) = 10000	! the 'w' has 10,000 elements
c     ipar(5) = 10	! use *GMRES(10) (e.g. FGMRES(10))
c     ipar(6) = 100	! use at most 100 matvec's
c     fpar(1) = 1.0E-6	! relative tolerance 1.0E-6
c     fpar(2) = 1.0E-10 ! absolute tolerance 1.0E-10
c     fpar(11) = 0.0	! clearing the FLOPS counter
c
c     After the above specifications, one can start to call an iterative
c     solver, say BCG. Here is a piece of pseudo-code showing how it can
c     be done,
c
c 10   call bcg(n,rhs,sol,ipar,fpar,w)
c      if (ipar(1).eq.1) then
c         call amux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
c         goto 10
c      else if (ipar(1).eq.2) then
c         call atmux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
c         goto 10
c      else if (ipar(1).eq.3) then
c         left preconditioner solver
c         goto 10
c      else if (ipar(1).eq.4) then
c         left preconditioner transposed solve
c         goto 10
c      else if (ipar(1).eq.5) then
c         right preconditioner solve
c         goto 10
c      else if (ipar(1).eq.6) then
c         right preconditioner transposed solve
c         goto 10
c      else if (ipar(1).eq.10) then
c         call my own stopping test routine
c         goto 10
c      else if (ipar(1).gt.0) then
c         ipar(1) is an unspecified code
c      else
c         the iterative solver terminated with code = ipar(1)
c      end if
c
c     This segment of pseudo-code assumes the matrix is in CSR format,
c     AMUX and ATMUX are two routines from the SPARSKIT MATVEC module.
c     They performs matrix-vector multiplications for CSR matrices. The
c     input vectors are the second arguments, and the output vectors are
c     the third arguments in the calling sequences. For simplicity, we
c     did not show the name of the routine that performs the
c     preconditioning solves or the convergence tests.
c-----------------------------------------------------------------------
      subroutine cg(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,*)
c-----------------------------------------------------------------------
c     This is a implementation of the Conjugate Gradient (CG) method
c     for solving linear system.
c
c     NOTE: This is not the PCG algorithm. It is a regular CG algorithm.
c     To be consistent with the other solvers, the preconditioners are
c     applied by performing Ml^{-1} A Mr^{-1} P in place of A P in the
c     CG algorithm. The PCG uses its preconditioners very differently.
c
c     fpar(7) is used here internally to store <r, r>.
c     w(:,1) -- residual vector
c     w(:,2) -- P, the conjugate direction
c     w(:,3) -- A P, matrix multiply the conjugate direction
c     w(:,4) -- temporary storage for results of preconditioning
c     w(:,5) -- change in the solution (sol) is stored here until
c               termination of this solver
c-----------------------------------------------------------------------
c     external functions used
c
      real*8 distdot
      logical stopbis, brkdn
      external distdot, stopbis, brkdn, bisinit
c
c     local variables
c
      integer i
      real*8 alpha
      logical lp,rp
      save
c
c     check the status of the call
c
      if (ipar(1).eq.0) ipar(10) = 0

c      print *, 'ipar(10) is', ipar(10)
c      print *, 'lp is', lp
c      print *, 'rp is', rp

      goto (10, 20, 40, 50, 60, 70, 80), ipar(10)
c
c     initialization
c
c      print *, 'Hello before bisinit'
      call bisinit(n,ipar,5*n,1,lp,rp,w(1,5))

      if (ipar(1).lt.0) return
      ipar(1) = 1
      ipar(8) = n+1
      ipar(9) = ipar(8) + n
      ipar(10) = 1
      do i = 1, n
         w(i,2) = sol(i)
      end do
c
c     request for matrix vector multiplication A*x in the initialization
c
      fpar(11) = fpar(11) + 2*n
      return
 10   do i = 1, n
         w(i,2) = rhs(i) - w(i,3)
      end do
c
c     if left preconditioned
c
      if (lp) then
         ipar(1) = 3
         ipar(9) = 1
         ipar(10) = 2
         return
      endif
c
 20   if (lp) then
         do i = 1, n
            w(i,2) = w(i,1)
         end do
c         call dcopy1(n, w(1,1),1, w(1,2),1)
      else
         do i = 1, n
            w(i,1) = w(i,2)
         end do
c         call dcopy1(n, w(1,2),1, w(1,1),1)
      endif
c
      fpar(7) = distdot(n,w,1,w,1)
      fpar(11) = fpar(11) + 2 * n
      fpar(3) = sqrt(fpar(7))
      fpar(5) = fpar(3)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + 2 * n
      else
         fpar(4) = fpar(1) * fpar(3) + fpar(2)
      end if
c
c     before iteration can continue, we need to compute A * p, which
c     includes the preconditioning operations
c
 30   if (rp) then
         ipar(1) = 5
         ipar(8) = n + 1
         if (lp) then
            ipar(9) = ipar(8) + n
         else
            ipar(9) = 3*n + 1
         end if
         ipar(10) = 3
         return
      end if
c
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = n + 1
      endif
      if (lp) then
         ipar(9) = 3*n+1
      else
         ipar(9) = n+n+1
      end if
      ipar(10) = 4
      return
c
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = n+n+1
         ipar(10) = 5
         return
      end if
c
c     continuing with the iterations
c
 60   alpha = distdot(n,w(1,2),1,w(1,3),1)
      fpar(11) = fpar(11) + 2*n
      if (brkdn(alpha,ipar)) goto 900
      alpha = fpar(7) / alpha
      do i = 1, n
         w(i,5) = w(i,5) + alpha * w(i,2)
         w(i,1) = w(i,1) - alpha * w(i,3)
      end do
      ipar(7) = ipar(7) + 1
      fpar(11) = fpar(11) + 4*n
c
c     are we ready to terminate ?
c
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(10) = 6
         return
      end if
 70   if (ipar(3).eq.999.and.ipar(11).eq.1) then
         goto 900
      else if (stopbis(n,ipar,1,fpar,w,w(1,2),alpha)) then
         goto 900
      else
         alpha = fpar(5)*fpar(5) / fpar(7)
         fpar(7) = fpar(5)*fpar(5)
         do i = 1, n
            w(i,2) = w(i,1) + alpha * w(i,2)
         end do
         fpar(11) = fpar(11) + 2*n
         goto 30
      end if
c
c     tide up -- necessary to accommodate the right-preconditioning
c
 900  if (rp) then
         ipar(1) = 5
         ipar(8) = 4*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 7
         return
      end if
 80   if (rp) then
         call tidecg(n,ipar,fpar,sol,w(1,4))
      else
         call tidecg(n,ipar,fpar,sol,w(1,5))
      end if
c
      return
      end
c-----end-of-cg
c-----------------------------------------------------------------------
      subroutine cgnr(n,rhs,sol,ipar,fpar,wk)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n),sol(n),fpar(16),wk(n,*)
c-----------------------------------------------------------------------
c     CGNR -- Using CG algorithm solving A x = b by solving
c     Normal Residual equation: A^T A x = A^T b
c     As long as the matrix is not singular, A^T A is symmetric
c     positive definite, therefore CG (CGNR) will converge.
c
c     Usage of the work space:
c     wk(:,1) == residual vectore R
c     wk(:,2) == the conjugate direction vector P
c     wk(:,3) == a temp holds A P, or A^T R
c     wk(:,4) == a temp holds intermediate results of the preconditionings
c     wk(:,5) == a place to hold the moldification to SOL
c
c     size of the work space WK is required = 5*n
c-----------------------------------------------------------------------
c     external functions used
c
      real*8 distdot
      logical stopbis, brkdn
      external distdot, stopbis, brkdn, bisinit
c
c     local variables
c
      integer i
      real*8 alpha, zz, zzm1
      logical lp, rp
      save
c
c     check the status of the call
c
      if (ipar(1).eq.0) ipar(10) = 0
      goto (10, 20, 40, 50, 60, 70, 80, 90, 100, 110), ipar(10)
c
c     initialization
c
      call bisinit(n,ipar,5*n,1,lp,rp,wk(1,5))
      if (ipar(1).lt.0) return
      ipar(1) = 1
      ipar(8) = 1
      ipar(9) = 1 + n
      ipar(10) = 1
      do i = 1, n
         wk(i,1) = sol(i)
      end do
c      call dcopy1(n,sol,1,wk,1)
c
c     request for matrix vector multiplication A*x in the initialization
c
      return
 10   do i = 1, n
         wk(i,1) = rhs(i) - wk(i,2)
      end do
      fpar(11) = fpar(11) + n
c
c     if left preconditioned
c
      if (lp) then
         ipar(1) = 3
         ipar(10) = 2
         return
      endif
c
 20   if (lp) then
         do i = 1, n
            wk(i,1) = wk(i,2)
         end do
c         call dcopy1(n, wk(1,2),1, wk(1,1),1)
      endif
c
      zz = distdot(n,wk,1,wk,1)
      fpar(11) = fpar(11) + 2 * n
      fpar(3) = sqrt(zz)
      fpar(5) = fpar(3)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + 2 * n
      else
         fpar(4) = fpar(1) * fpar(3) + fpar(2)
      end if
c
c     normal iteration begins here, first half of the iteration
c     computes the conjugate direction
c
 30   continue
c
c     request the caller to perform a A^T r --> wk(:,3)
c
      if (lp) then
         ipar(1) = 4
         ipar(8) = 1
         if (rp) then
            ipar(9) = n + n + 1
         else
            ipar(9) = 3*n + 1
         end if
         ipar(10) = 3
         return
      end if
c
 40   ipar(1) = 2
      if (lp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 1
      end if
      if (rp) then
         ipar(9) = 3*n + 1
      else
         ipar(9) = n + n + 1
      end if
      ipar(10) = 4
      return
c
 50   if (rp) then
         ipar(1) = 6
         ipar(8) = ipar(9)
         ipar(9) = n + n + 1
         ipar(10) = 5
         return
      end if
c
 60   zzm1 = zz
      zz = distdot(n,wk(1,3),1,wk(1,3),1)
      fpar(11) = fpar(11) + 2 * n
      if (brkdn(zz,ipar)) goto 900
      ipar(7) = ipar(7) + 1
      if (ipar(7).gt.1) then
         alpha = zz / zzm1
         do i = 1, n
            wk(i,2) = wk(i,3) + alpha * wk(i,2)
         end do
         fpar(11) = fpar(11) + 2 * n
      else
         do i = 1, n
            wk(i,2) = wk(i,3)
         end do
c         call dcopy1(n,wk(1,3),1,wk(1,2),1)
      end if
c
c     before iteration can continue, we need to compute A * p
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = n + 1
         if (lp) then
            ipar(9) = ipar(8) + n
         else
            ipar(9) = 3*n + 1
         end if
         ipar(10) = 6
         return
      end if
c
 70   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = n + 1
      endif
      if (lp) then
         ipar(9) = 3*n+1
      else
         ipar(9) = n+n+1
      end if
      ipar(10) = 7
      return
c
 80   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = n+n+1
         ipar(10) = 8
         return
      end if
c
c     update the solution -- acclumate the changes in w(:,5)
c
 90   alpha = distdot(n,wk(1,3),1,wk(1,3),1)
      fpar(11) = fpar(11) + 2 * n
      if (brkdn(alpha,ipar)) goto 900
      alpha = zz / alpha
      do i = 1, n
         wk(i,5) = wk(i,5) + alpha * wk(i,2)
         wk(i,1) = wk(i,1) - alpha * wk(i,3)
      end do
      ipar(7) = ipar(7) + 1
      fpar(11) = fpar(11) + 4 * n
c
c     are we ready to terminate ?
c
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(10) = 9
         return
      end if
 100  if (ipar(3).eq.999.and.ipar(11).eq.1) then
         goto 900
      else if (stopbis(n,ipar,1,fpar,wk,wk(1,2),alpha)) then
         goto 900
      else
         goto 30
      end if
c
c     tide up -- necessary to accommodate the right-preconditioning
c
 900  if (rp) then
         ipar(1) = 5
         ipar(8) = 4*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 10
         return
      end if
 110  if (rp) then
         call tidecg(n,ipar,fpar,sol,wk(1,4))
      else
         call tidecg(n,ipar,fpar,sol,wk(1,5))
      end if
      return
      end
c-----end-of-cgnr
c-----------------------------------------------------------------------
      subroutine bcg(n,rhs,sol,ipar,fpar,w)
      implicit none
      integer n, ipar(16)
      real*8  fpar(16), rhs(n), sol(n), w(n,*)
c-----------------------------------------------------------------------
c     BCG: Bi Conjugate Gradient method. Programmed with reverse
c     communication, see the header for detailed specifications
c     of the protocol.
c
c     in this routine, before successful return, the fpar's are
c     fpar(3) == initial residual norm
c     fpar(4) == target residual norm
c     fpar(5) == current residual norm
c     fpar(7) == current rho (rhok = <r, s>)
c     fpar(8) == previous rho (rhokm1)
c
c     w(:,1) -- r, the residual
c     w(:,2) -- s, the dual of the 'r'
c     w(:,3) -- p, the projection direction
c     w(:,4) -- q, the dual of the 'p'
c     w(:,5) -- v, a temp to store A*p, or A*q.
c     w(:,6) -- a temp to store intermediate results
c     w(:,7) -- changes in the solution
c-----------------------------------------------------------------------
c     external routines used
c
      real*8 distdot
      logical stopbis,brkdn
      external distdot, stopbis, brkdn
c
c     local variables
c
      integer i
      real*8 alpha
      logical rp, lp
c
c     status of the program
c
      if (ipar(1).eq.0) ipar(10) = 0
      goto (10, 20, 40, 50, 60, 70, 80, 90, 100, 110), ipar(10)
c
c     initialization, initial residual
c
      call bisinit(n,ipar,7*n,1,lp,rp,w(1,7))
      if (ipar(1).lt.0) return
      ipar(1) = 1
      ipar(8) = 3*n+1
      ipar(9) = ipar(8) + n
      do i = 1, n
         w(i,4) = sol(i)
      end do
c      call dcopy1(n,sol,1,w(1,4),1)
      ipar(10) = 1
      return
 10   do i = 1, n
         w(i,1) = rhs(i) - w(i,5)
      end do
      if (lp) then
         ipar(1) = 3
         ipar(8) = 1
         ipar(9) = n+1
         ipar(10) = 2
         return
      end if
c
 20   if (lp) then
         do i = 1, n
            w(i,1) = w(i,2)
            w(i,3) = w(i,2)
            w(i,4) = w(i,2)
         end do
      else
         do i = 1, n
            w(i,2) = w(i,1)
            w(i,3) = w(i,1)
            w(i,4) = w(i,1)
         end do
      end if
c
      fpar(11) = fpar(11) + 2 * n
      fpar(7) = distdot(n,w,1,w,1)
      fpar(3) = sqrt(fpar(7))
      fpar(5) = fpar(3)
      fpar(8) = 1.0D0
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + 2 * n
      else
         fpar(4) = fpar(1) * fpar(3) + fpar(2)
      end if
      if (ipar(3).ge.0.and.fpar(5).le.fpar(4)) then
         fpar(6) = fpar(5)
         goto 900
      end if
c
c     end of initialization, begin iteration, v = A p
c
 30   if (rp) then
         ipar(1) = 5
         ipar(8) = n + n + 1
         if (lp) then
            ipar(9) = 4*n + 1
         else
            ipar(9) = 5*n + 1
         end if
         ipar(10) = 3
         return
      end if
c
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = n + n + 1
      end if
      if (lp) then
         ipar(9) = 5*n + 1
      else
         ipar(9) = 4*n + 1
      end if
      ipar(10) = 4
      return
c
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 4*n + 1
         ipar(10) = 5
         return
      end if
c
 60   alpha = distdot(n,w(1,4),1,w(1,5),1)
      fpar(11) = fpar(11) + 2 * n
      if (brkdn(alpha,ipar)) goto 900
      alpha = fpar(7) / alpha
      do i = 1, n
         w(i,7) = w(i,7) + alpha * w(i,3)
         w(i,1) = w(i,1) - alpha * w(i,5)
      end do
      ipar(7) = ipar(7) + 1
      fpar(11) = fpar(11) + 4 * n
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(10) = 6
         return
      end if
 70   if (ipar(3).eq.999.and.ipar(11).eq.1) then
         goto 900
      else if (stopbis(n,ipar,1,fpar,w,w(1,3),alpha)) then
         goto 900
      end if
c
c     A^t * x
c
      if (lp) then
         ipar(1) = 4
         ipar(8) = 3*n + 1
         if (rp) then
            ipar(9) = 4*n + 1
         else
            ipar(9) = 5*n + 1
         end if
         ipar(10) = 7
         return
      end if
c
 80   ipar(1) = 2
      if (lp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 3*n + 1
      end if
      if (rp) then
         ipar(9) = 5*n + 1
      else
         ipar(9) = 4*n + 1
      end if
      ipar(10) = 8
      return
c
 90   if (rp) then
         ipar(1) = 6
         ipar(8) = ipar(9)
         ipar(9) = 4*n + 1
         ipar(10) = 9
         return
      end if
c
 100  continue
      do i = 1, n
         w(i,2) = w(i,2) - alpha * w(i,5)
      end do
c      call daxpy1(n,-alpha,w(1,5),1,w(1,2),1)
      fpar(8) = fpar(7)
      fpar(7) = distdot(n,w,1,w(1,2),1)
      fpar(11) = fpar(11) + 4 * n
      if (brkdn(fpar(7), ipar)) return
      alpha = fpar(7) / fpar(8)
      do i = 1, n
         w(i,3) = w(i,1) + alpha * w(i,3)
         w(i,4) = w(i,2) + alpha * w(i,4)
      end do
      fpar(11) = fpar(11) + 4 * n
c
c     end of the iterations
c
      ipar(7) = ipar(7) + 1
      goto 30
c
c     some clean up job to do
c
 900  if (rp) then
         ipar(1) = 5
         ipar(8) = 6*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 10
         return
      end if
 110  if (rp) then
         call tidecg(n,ipar,fpar,sol,w(1,6))
      else
         call tidecg(n,ipar,fpar,sol,w(1,7))
      end if
      return
      end
c-----end-of-bcg
c-----------------------------------------------------------------------
      subroutine bcgstab(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,*)
c-----------------------------------------------------------------------
c     BCGSTAB --- Bi Conjugate Gradient stablized (BCGSTAB)
c     This is an improved BCG routine. (1) no matrix transpose is
c     involved. (2) the convergence is smoother.
c
c     in this routine, before successful return, the fpar's are
c     fpar(3) == initial residual norm squared
c     fpar(4) == target residual norm squared
c     fpar(5) == current residual norm squared
c     fpar(7) == current rho (rhok = <r, s>)
c     fpar(8) == previous rho (rhokm1)
c     fpar(9) == beta (= rhok / rhokm1) and alpha (= rhok / <pk, v>)
c     fpar(10) == omega
c-----------------------------------------------------------------------
c     external routines used
c
      real*8 distdot
      logical stopbis, brkdn
      external distdot, stopbis, brkdn
c
c     local variables
c
      integer i
      real*8 beta
      logical lp, rp
      save
c
c     where to go
c
      if (ipar(1).eq.0) ipar(10) = 0
      goto (10, 20, 40, 50, 60, 70, 80, 90, 100, 110) ipar(10)
c
c     compute initial residual
c
      call bisinit(n,ipar,9*n,1,lp,rp,w(1,9))
      if (ipar(1).lt.0) return
c
      ipar(1) = 1
      ipar(8) = 3*n+1
      ipar(9) = ipar(8) + n
      do i = 1, n
         w(i,4) = sol(i)
      end do
c      call dcopy1(n,sol,1,w(1,4),1)
      ipar(10) = 1
      return
c
c     A * x
c
 10   do i = 1, n
         w(i,4) = rhs(i) - w(i,5)
      end do
      fpar(11) = fpar(11) + n
      if (lp) then
         ipar(1) = 3
         ipar(9) = 1
         ipar(10) = 2
         return
      end if
c
 20   if (lp) then
         do i = 1, n
            w(i,4) = w(i,1)
            w(i,6) = w(i,1)
         end do
      else
         do i = 1, n
            w(i,1) = w(i,4)
            w(i,6) = w(i,4)
         end do
      end if
c
      fpar(3) = sqrt(distdot(n,w,1,w,1))
      fpar(11) = fpar(11) + 2 * n
      fpar(5) = fpar(3)
      fpar(7) = fpar(3)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + 2 * n
      else
         fpar(4) = fpar(1) * fpar(3) + fpar(2)
      end if
      if (ipar(3).gt.0 .and. fpar(5).le.fpar(4)) then
         fpar(6) = fpar(5)
         goto 100
      end if
c
c     beginning of the iterations, A*x including preconditioning
c
 30   if (rp) then
         ipar(1) = 5
         ipar(8) = 3*n+1
         if (lp) then
            ipar(9) = ipar(8)+n
         else
            ipar(9) = 7*n + 1
         end if
         ipar(10) = 3
         return
      end if
c
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 3*n+1
      end if
      if (lp) then
         ipar(9) = 7*n + 1
      else
         ipar(9) = 4*n + 1
      end if
      ipar(10) = 4
      return
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 4*n + 1
         ipar(10) = 5
         return
      end if
c
 60   fpar(9) = distdot(n,w(1,6),1,w(1,5),1)
      ipar(7) = ipar(7) + 1
      fpar(11) = fpar(11) + 2 * n
      if (brkdn(fpar(9), ipar)) goto 900
      fpar(9) = fpar(7) / fpar(9)
      do i = 1, n
         w(i,2) = w(i,1) - fpar(9) * w(i,5)
      end do
      fpar(11) = fpar(11) + 2 * n
c
c     the second A * x
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = n+1
         if (lp) then
            ipar(9) = ipar(8)+n
         else
            ipar(9) = 7*n + 1
         end if
         ipar(10) = 6
         return
      end if
c
 70   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = n+1
      end if
      if (lp) then
         ipar(9) = 7*n + 1
      else
         ipar(9) = n + n + 1
      end if
      ipar(10) = 7
      return
 80   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = n + n + 1
         ipar(10) = 8
         return
      end if
 90   fpar(10) = distdot(n,w(1,3),1,w(1,3),1)
      fpar(11) = fpar(11) + 2 * n
      if (brkdn(fpar(10),ipar)) goto 900
      fpar(10) = distdot(n,w(1,2),1,w(1,3),1) / fpar(10)
      do i = 1, n
         w(i,7) = fpar(9) * w(i,4) + fpar(10) * w(i,2)
         w(i,9) = w(i,9) + w(i,7)
         w(i,1) = w(i,2) - fpar(10) * w(i,3)
      end do
      fpar(11) = fpar(11) + 8 * n
c
      if (ipar(7).eq.1) then
         fpar(3) = sqrt(distdot(n,w(1,7),1,w(1,7),1))
         fpar(11) = fpar(11) + 2 * n
         if (ipar(3).eq.-1) then
            fpar(4) = fpar(1) * sqrt(fpar(3)) + fpar(2)
         end if
      end if
      ipar(7) = ipar(7) + 1
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(10) = 9
         return
      end if
 100  if (ipar(3).eq.999.and.ipar(11).eq.1) then
         goto 900
      else if (stopbis(n,ipar,2,fpar,w,w(1,7),1.0D0)) then
         goto 900
      end if
c
      fpar(8) = fpar(7)
      fpar(7) = distdot(n,w(1,6),1,w(1,1),1)
      beta = fpar(7) * fpar(9) / (fpar(8) * fpar(10))
      do i = 1, n
         w(i,4) = w(i,1) + beta * (w(i,4) - fpar(10) * w(i,5))
      end do
      fpar(11) = fpar(11) + 6 * n
c
c     end of an iteration
c
      goto 30
c
c     some clean up job to do
c
 900  if (rp) then
         ipar(1) = 5
         ipar(8) = 8*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 10
         return
      end if
 110  if (rp) then
         call tidecg(n,ipar,fpar,sol,w(1,8))
      else
         call tidecg(n,ipar,fpar,sol,w(1,9))
      end if
c
      return
      end
c-----end-of-bcgstab
c-----------------------------------------------------------------------
      subroutine tfqmr(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,*)
c-----------------------------------------------------------------------
c     TFQMR --- transpose-free Quasi-Minimum Residual method
c     This is developed from BCG based on the principle of Quasi-Munimum
c     Residual, and it is transpose-free.
c
c     It uses approximate resudual norm.
c
c     Internally, the fpar's are used as following:
c     fpar(3) --- initial residual norm squared
c     fpar(4) --- target residual norm squared
c     fpar(5) --- current residual norm squared
c     fpar(7) --- alpha (= rho / sigma) and beta (= rhok / rhokm1)
c     fpar(8) --- rho (= <r_0, r>)
c     fpar(9) --- tao (= tao * theta * c)
c     fpar(10) --- theta * eta
c
c     w(:,1) -- R, residual
c     w(:,2) -- R0, the initial residual
c     w(:,3) -- W
c     w(:,4) -- Y
c     w(:,5) -- Z
c     w(:,6) -- A * Y
c     w(:,7) -- A * Z
c     w(:,8) -- V
c     w(:,9) -- D
c     w(:,10) -- intermediate results of preconditioning
c     w(:,11) -- changes in the solution
c-----------------------------------------------------------------------
c     external functions
c
      real*8 distdot
      logical stopbis, brkdn
      external stopbis, brkdn, distdot
c
c     local variables, variables that are not to be saved across calls
c
      integer i
      logical lp, rp
      real*8 eta,sigma,theta,te,alpha,rho,tao
      save
c
c     status of the call (where to go)
c
      if (ipar(1).eq.0) ipar(10) = 0
      goto (10,20,40,50,60,70,80,90,100,110), ipar(10)
c
c     initializations
c
      call bisinit(n,ipar,11*n,2,lp,rp,w(1,11))
      if (ipar(1).lt.0) return
      ipar(1) = 1
      ipar(8) = 1
      ipar(9) = 1 + 6*n
      do i = 1, n
         w(i,1) = sol(i)
      end do
c      call dcopy1(n,sol,1,w,1)
      ipar(10) = 1
      return
 10   do i = 1, n
         w(i,1) = rhs(i) - w(i,7)
         w(i,9) = 0.0D0
      end do
      fpar(11) = fpar(11) + n
c
      if (lp) then
         ipar(1) = 3
         ipar(9) = n+1
         ipar(10) = 2
         return
      end if
 20   if (lp) then
         do i = 1, n
            w(i,1) = w(i,2)
            w(i,3) = w(i,2)
         end do
      else
         do i = 1, n
            w(i,2) = w(i,1)
            w(i,3) = w(i,1)
         end do
      end if
c
      fpar(5) = sqrt(distdot(n,w,1,w,1))
      fpar(3) = fpar(5)
      tao = fpar(5)
      fpar(11) = fpar(11) + n + n
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + n + n
      else
         fpar(4) = fpar(1) * tao + fpar(2)
      end if
      te = 0.0D0
c
c     begin iteration
c
 30   sigma = rho
      rho = distdot(n,w(1,2),1,w(1,3),1)
      fpar(11) = fpar(11) + n + n
      if (brkdn(rho,ipar,fpar)) goto 900
      if (ipar(7).eq.0) then
         alpha = 0.0D0
      else
         alpha = rho / sigma
      end if
      do i = 1, n
         w(i,4) = w(i,3) + alpha * w(i,5)
      end do
      fpar(11) = fpar(11) + n + n
c
c     A * x -- with preconditioning
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = 3*n + 1
         if (lp) then
            ipar(9) = 5*n + 1
         else
            ipar(9) = 9*n + 1
         end if
         ipar(10) = 3
         return
      end if
c
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 3*n + 1
      end if
      if (lp) then
         ipar(9) = 9*n + 1
      else
         ipar(9) = 5*n + 1
      end if
      ipar(10) = 4
      return
c
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 5*n + 1
         ipar(10) = 5
         return
      end if
 60   do i = 1, n
         w(i,8) = w(i,6) + alpha * (w(i,7) + alpha * w(i,8))
      end do
      sigma = distdot(n,w(1,2),1,w(1,8),1)
      fpar(11) = fpar(11) + 4 * n
      if (brkdn(sigma,ipar,fpar)) goto 900
      alpha = rho / sigma
      do i = 1, n
         w(i,5) = w(i,4) - alpha * w(i,8)
      end do
      ipar(7) = ipar(7) + 1
      fpar(11) = fpar(11) + 2*n
c
c     the second A * x
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = 4*n + 1
         if (lp) then
            ipar(9) = 6*n + 1
         else
            ipar(9) = 9*n + 1
         end if
         ipar(10) = 6
         return
      end if
c
 70   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 4*n + 1
      end if
      if (lp) then
         ipar(9) = 9*n + 1
      else
         ipar(9) = 6*n + 1
      end if
      ipar(10) = 7
      return
c
 80   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 6*n + 1
         ipar(10) = 8
         return
      end if
 90   continue
      do i = 1, n
         w(i,3) = w(i,3) - alpha * w(i,6)
      end do
c     call daxpy1(n, -alpha, w(1,6), 1, w(1,3), 1)
      ipar(7) = ipar(7) + 1
c
c     update I
c
      theta = distdot(n,w(1,3),1,w(1,3),1) / (tao*tao)
      sigma = 1.0D0 / (1.0D0 + theta)
      tao = tao * sqrt(sigma * theta)
      fpar(11) = fpar(11) + 4*n
      if (brkdn(tao,ipar,fpar)) goto 900
      eta = sigma * alpha
      sigma = te / alpha
      te = theta * eta
      do i = 1, n
         w(i,9) = w(i,4) + sigma * w(i,9)
         w(i,11) = w(i,11) + eta * w(i,9)
         w(i,3) = w(i,3) - alpha * w(i,7)
      end do
      if (ipar(7).eq.0) then
         if (ipar(3).eq.-1) then
            fpar(3) = eta * sqrt(distdot(n,w(1,9),1,w(1,9),1))
            fpar(4) = fpar(1)*fpar(3) + fpar(2)
            fpar(11) = fpar(11) + n + n
         end if
      end if
c
c     update II
c
      theta = distdot(n,w(1,3),1,w(1,3),1) / (tao*tao)
      sigma = 1.0D0 / (1.0D0 + theta)
      tao = tao * sqrt(sigma * theta)
      fpar(11) = fpar(11) + 8*n
      if (brkdn(tao,ipar,fpar)) goto 900
      eta = sigma * alpha
      sigma = te / alpha
      te = theta * eta
      do i = 1, n
         w(i,9) = w(i,5) + sigma * w(i,9)
         w(i,11) = w(i,11) + eta * w(i,9)
      end do
c
c     error/resdual norm
c
      fpar(11) = fpar(11) + 4*n
c
c     this is the correct over-estimate
c      fpar(5) = sqrt(real(ipar(7)+1)) * tao
c     this is an approximation
      fpar(5) = tao
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(10) = 9
         return
      end if
      if (ipar(3).lt.0) then
         fpar(6) = eta * sqrt(distdot(n,w(1,9),1,w(1,9),1))
         fpar(11) = fpar(11) + n + n
      else
         fpar(6) = fpar(5)
      end if
      if (fpar(6).gt.fpar(4) .and. (ipar(7).lt.ipar(6)
     +     .or. ipar(6).le.0)) goto 30
 100  if (ipar(3).eq.999.and.ipar(11).eq.0) goto 30
c
c     clean up
c
 900  if (rp) then
         ipar(1) = 5
         ipar(8) = 10*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 10
         return
      end if
 110  if (rp) then
         call tidecg(n,ipar,fpar,sol,w(1,10))
      else
         call tidecg(n,ipar,fpar,sol,w(1,11))
      end if
c
      return
      end
c-----end-of-tfqmr
c-----------------------------------------------------------------------
      subroutine gmres(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
c-----------------------------------------------------------------------
c     This a version of GMRES implemented with reverse communication.
c     It is a simple restart version of the GMRES algorithm.
c
c     ipar(5) == the dimension of the Krylov subspace
c     after every ipar(5) iterations, the GMRES will restart with
c     the updated solution and recomputed residual vector.
c
c     the space of the `w' is used as follows:
c     (1) the basis for the Krylov subspace, size n*(m+1);
c     (2) the Henssenberg matrix, only the upper triangular
c     portion of the matrix is stored, size (m+1)*m/2 + 1
c     (3) three vectors, all are of size m, they are
c     the cosine and sine of the Given rotations, the third one helds
c     the residuals, it is of size m+1.
c
c     TOTAL SIZE REQUIRED == (n+3)*(m+2) + (m+1)*m/2
c     Note: m == ipar(5). The default value for this is 15 if
c     ipar(5) <= 1.
c-----------------------------------------------------------------------
c     external functions used
c
      real*8 distdot
      external distdot
c
c     local variables, ptr and p2 are temporary pointers,
c     hens points to the Hessenberg matrix,
c     vc, vs point to the cosines and sines of the Givens rotations
c     vrn points to the vectors of residual norms, more precisely
c     the right hand side of the least square problem solved.
c
      integer i,ii,idx,j,k,m,ptr,p2,hens,vc,vs,vrn
      real*8 alpha, c, s, deps
      logical lp, rp
      parameter (deps=1.0D-33)
      save
c
c     check the status of the call
c
      if (ipar(1).eq.0) ipar(10) = 0
      goto (10, 20, 30, 40, 50, 60, 70) ipar(10)
c
c     initialization
c
      if (ipar(5).le.1) then
         m = 15
      else
         m = ipar(5)
      end if
      idx = n * (m+1)
      hens = idx + n
      vc = hens + (m+1) * m / 2 + 1
      vs = vc + m
      vrn = vs + m
      i = vrn + m + 1
      call bisinit(n,ipar,i,1,lp,rp,w(idx+1))
      if (ipar(1).lt.0) return
c
c     request for matrix vector multiplication A*x in the initialization
c
 100  ipar(1) = 1
      ipar(8) = n+1
      ipar(9) = 1
      ipar(10) = 1
      k = 0
      do i = 1, n
         w(n+i) = sol(i)
      end do
c      call dcopy1(n,sol,1,w(1+n),1)
      return
 10   if (lp) then
         do i = 1, n
            w(n+i) = rhs(i) - w(i)
         end do
         ipar(1) = 3
         ipar(10) = 2
         return
      else
         do i = 1, n
            w(i) = rhs(i) - w(i)
         end do
      end if
      fpar(11) = fpar(11) + n
c
 20   alpha = sqrt(distdot(n,w,1,w,1))
      fpar(11) = fpar(11) + 2*n
      if (ipar(7).eq.0) then
         if (abs(ipar(3)).eq.2) then
            fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
            fpar(11) = fpar(11) + 2*n
         else
            fpar(4) = fpar(1) * alpha + fpar(2)
         end if
         fpar(3) = alpha
      end if
      fpar(5) = alpha
      w(vrn+1) = alpha
      if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999) then
         ipar(1) = 0
         fpar(6) = alpha
         goto 300
      end if
      alpha = 1.0D0 / alpha
      do ii = 1, n
         w(ii) = alpha * w(ii)
      end do
c      call dscal1(n,alpha,w,1)
      fpar(11) = fpar(11) + n
c
c     request for a matrix vector multiplication
c
 110  k = k + 1
      if (rp) then
         ipar(1) = 5
         ipar(8) = k*n - n + 1
         if (lp) then
            ipar(9) = k*n + 1
         else
            ipar(9) = idx + 1
         end if
         ipar(10) = 3
         return
      end if
c
 30   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = (k-1)*n + 1
      end if
      if (lp) then
         ipar(9) = idx + 1
      else
         ipar(9) = 1 + k*n
      end if
      ipar(10) = 4
      return
c
 40   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = k*n + 1
         ipar(10) = 5
         return
      end if
c
c     Modified Gram-Schmidt orthogonalization procedure
c     temporary pointer 'ptr' is pointing to the current column of the
c     Henssenberg matrix. 'p2' points to the new basis vector
c
 50   ptr = k * (k - 1) / 2 + hens
      p2 = ipar(9)
      ipar(7) = ipar(7) + 1
      do i = 1, k
         j = (i-1)*n + 1
         alpha = distdot(n,w(p2),1,w(j),1)
         do ii = 0, n-1
            w(p2+ii) = w(p2+ii) - alpha * w(j+ii)
         end do
c         call daxpy1(n,-alpha,w(j),1,w(p2),1)
         w(ptr+i) = alpha
      end do
      alpha = sqrt(distdot(n,w(p2),1,w(p2),1))
      if (alpha.le.deps) then
         k = k - 1
         goto 200
      end if
      w(ptr+i) = alpha
      alpha = 1.0D0 / alpha
      do ii = 0, n-1
         w(p2+ii) = w(p2+ii) * alpha
      end do
c      call dscal1(n,alpha,w(p2),1)
c
c     apply previous Givens rotations and generate a new one to eliminate
c     the subdiagonal element.
c
      p2 = ptr + 1
      do i = 1, k-1
         ptr = p2
         p2 = p2 + 1
         alpha = w(ptr)
         c = w(vc+i)
         s = w(vs+i)
         w(ptr) = c * alpha + s * w(p2)
         w(p2) = c * w(p2) - s * alpha
      end do
      call givens(w(p2), w(p2+1), c, s)
      w(vc+k) = c
      w(vs+k) = s
      p2 = vrn + k
      alpha = - s * w(p2)
      w(p2) = c * w(p2)
      w(p2+1) = alpha
c
c     end of one Arnoldi iteration, alpha will store the estimated
c     residual norm at current stage
c
      fpar(11) = fpar(11) + 4*k*n + 3*(n + k + k)
      alpha = abs(alpha)
      fpar(5) = alpha
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(10) = 6
         return
      end if
      if ((ipar(3).lt.0 .or. alpha.gt.fpar(4)) .and.
     +     k.lt.m) goto 110
 60   if (ipar(3).eq.999 .and. ipar(11).eq.0) goto 110
c
c     update the approximate solution, first solve the upper triangular
c     system, temporary pointer ptr points to the Henssenberg matrix,
c     p2 points to the right-hand-side (also the solution) of the system.
c
 200  ptr = hens + k * (k + 1) / 2
      w(p2) = w(p2) / w(ptr)
      do i = k-1, 1, -1
         ptr = ptr - i
         do ii = 1, i
            w(vrn+ii) = w(vrn+ii) - w(p2) * w(ptr-1+ii)
         end do
c         call daxpy1(i, -w(p2), w(ptr), 1, w(vrn+1), 1)
         p2 = p2 - 1
         ptr = ptr - 1
         w(p2) = w(p2) / w(ptr)
      end do
c
      do ii = 1, n
         w(ii) = w(ii) * w(p2)
      end do
c      call dscal1(n, w(p2), w, 1)
      do i = 1, k-1
         p2 = p2 + 1
         do ii = 1, n
            w(ii) = w(ii) + w(p2) * w(i*n+ii)
         end do
c         call daxpy1(n, w(p2), w(i*n+1), 1, w, 1)
      end do
      fpar(11) = fpar(11) + 4*(k-1)*n + n
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = 1
         ipar(9) = idx + 1
         ipar(10) = 7
         return
      end if
c
 70   if (rp) then
         do i = 1, n
            sol(i) = sol(i) + w(idx+i)
         end do
      else
         do i = 1, n
            sol(i) = sol(i) + w(i)
         end do
      endif
c
c     process the complete stopping criteria
c
      if (ipar(3).lt.0) then
         if (ipar(7).le.m) then
            fpar(3) = abs(w(vrn+1))
            if (m.eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
         end if
         alpha = abs(w(vrn+k))
      end if
      fpar(6) = alpha
      fpar(11) = fpar(11) + n
c
c     do we need to restart ?
c
      if ((ipar(7).lt.ipar(6) .or. ipar(6).le.0).and.
     $     ((ipar(3).eq.999.and.ipar(11).eq.0) .or.
     $     (ipar(3).ne.999.and.alpha.gt.fpar(4)))) goto 100
c
c     termination, set error code, compute convergence rate
c
      if (ipar(1).gt.0) then
         if (ipar(3).eq.999 .and. ipar(11).eq.1) then
            ipar(1) = 0
         else if (ipar(3).ne.999 .and. alpha.le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else
            ipar(1) = -10
         end if
      end if
 300  if (alpha.ne.0.0D0) then
         fpar(7) = log10(fpar(3) / alpha) / ipar(7)
      else
         fpar(7) = 0.0D0
      end if
      return
      end
c-----end-of-gmres
c-----------------------------------------------------------------------
      subroutine dqgmres(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
c-----------------------------------------------------------------------
c     DQGMRES -- Flexible Direct version of Quasi-General Minimum
c     Residual method. The right preconditioning can be varied from
c     step to step.
c
c     Work space used = n + lb * (2*n+4)
c     where lb = ipar(5) + 1 (default 16 if ipar(5) <= 1)
c-----------------------------------------------------------------------
c     local variables
c
      integer i,ii,j,jp1,j0,k,ptrw,ptrv,iv,iw,ic,is,ihm,ihd,lb
      real*8 alpha,beta,psi,c,s,deps,distdot
      logical lp,rp,full
      external distdot,bisinit
      parameter (deps = 1.0D-16)
      save
c
c     where to go
c
      if (ipar(1).eq.0) ipar(10) = 0
      goto (10, 20, 40, 50, 60, 70) ipar(10)
c
c     locations of the work arrays. The arrangement is:
c     w(1:n) -- temporary storage for the results of the preconditioning
c     w(iv+1:iw) -- the V's
c     w(iw+1:ic) -- the W's
c     w(ic+1:is) -- the COSINEs of the Givens rotations
c     w(is+1:ihm) -- the SINEs of the Givens rotations
c     w(ihm+1:ihd) -- the last column of the Hessenberg matrix
c     w(ihd+1:i) -- the inverse of the diagnals of the Hessenberg matrix
c
      if (ipar(5).le.1) then
         lb = 16
      else
         lb = ipar(5) + 1
      end if
      iv = n
      iw = iv + lb * n
      ic = iw + lb * n
      is = ic + lb
      ihm = is + lb
      ihd = ihm + lb
      i = ihd + lb
c
c     parameter check, initializations
c
      full = .false.
      call bisinit(n,ipar,i,1,lp,rp,w(1))
      if (ipar(1).lt.0) return
      ipar(1) = 1
      if (lp) then
         do ii = 1, n
            w(iv+ii) = sol(ii)
         end do
c         call dcopy1(n, sol, 1, w(iv+1), 1)
         ipar(8) = iv+1
         ipar(9) = 1
      else
         do ii = 1, n
            w(ii) = sol(ii)
         end do
c         call dcopy1(n, sol, 1, w(1), 1)
         ipar(8) = 1
         ipar(9) = iv+1
      endif
      ipar(10) = 1
      return
c
 10   if (lp) then
         do i = 1, n
            w(i) = rhs(i) - w(i)
         end do
         ipar(1) = 3
         ipar(8) = 1
         ipar(9) = iv+1
         ipar(10) = 2
         return
      else
         do i = 1, n
            w(iv+i) = rhs(i) - w(iv+i)
         end do
      end if
      fpar(11) = fpar(11) + n
c
 20   alpha = sqrt(distdot(n, w(iv+1), 1, w(iv+1), 1))
      fpar(11) = fpar(11) + (n + n)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + 2*n
      else
         fpar(4) = fpar(1) * alpha + fpar(2)
      end if
      fpar(3) = alpha
      fpar(5) = alpha
      psi = alpha
      if (alpha.le.fpar(4)) then
         ipar(1) = 0
         fpar(6) = alpha
         goto 80
      end if
      alpha = 1.0D0 / alpha
      do i = 1, n
         w(iv+i) = w(iv+i) * alpha
      end do
c      call dscal1(n,alpha,w(iv+1),1)
      fpar(11) = fpar(11) + n
      j = 0
c
c     iterations start here
c
 30   j = j + 1
      if (j.gt.lb) j = j - lb
      jp1 = j + 1
      if (jp1.gt.lb) jp1 = jp1 - lb
      ptrv = iv + (j-1)*n + 1
      ptrw = iv + (jp1-1)*n + 1
      if (.not.full) then
         if (j.gt.jp1) full = .true.
      end if
      if (full) then
         j0 = jp1+1
         if (j0.gt.lb) j0 = j0 - lb
      else
         j0 = 1
      end if
c
c     request the caller to perform matrix-vector multiplication and
c     preconditioning
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = ptrv
         ipar(9) = ptrv + iw - iv
         ipar(10) = 3
         return
      else
         do i = 0, n-1
            w(ptrv+iw-iv+i) = w(ptrv+i)
         end do
c         call dcopy1(n,w(ptrv),1,w(ptrv+iw-iv),1)
      end if
c
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = ptrv
      end if
      if (lp) then
         ipar(9) = 1
      else
         ipar(9) = ptrw
      end if
      ipar(10) = 4
      return
c
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = ptrw
         ipar(10) = 5
         return
      end if
c
c     compute the last column of the Hessenberg matrix
c     modified Gram-schmidt procedure, orthogonalize against (lb-1)
c     previous vectors
c
 60   i = j0
      do while (i.ne.jp1)
         k = iv+(i-1)*n+1
         alpha = distdot(n, w(k), 1, w(ptrw), 1)
         do ii = 0, n-1
            w(ptrw+ii) = w(ptrw+ii) - alpha * w(k+ii)
         end do
c         call daxpy1(n, -alpha, w(k), 1, w(ptrw), 1)
         w(ihm+i) = alpha
         i = i + 1
         if (i.gt.lb) i = 1
      end do
      beta = sqrt(distdot(n, w(ptrw), 1, w(ptrw), 1))
      if (beta .gt. deps) then
         alpha = 1.0D0 / beta
         do ii = 0, n-1
            w(ptrw+ii) = w(ptrw+ii) * alpha
         end do
c         call dscal1(n, alpha, w(ptrw), 1)
      end if
c
c     imcomplete factorization (QR factorization through Givens rotations)
c     (1) apply previous rotations [(lb-1) of them]
c     (2) generate a new rotation
c
      if (full) then
         w(ihm+jp1) = w(ihm+j0) * w(is+jp1)
         w(ihm+j0) = w(ihm+j0) * w(ic+jp1)
      end if
      i = j0
      do while (i.ne.j)
         k = i+1
         if (k.gt.lb) k = k - lb
         c = w(ic+i)
         s = w(is+i)
         alpha = w(ihm+i)
         w(ihm+i) = c * alpha + s * w(ihm+k)
         w(ihm+k) = c * w(ihm+k) - s * alpha
         i = k
      end do
      call givens(w(ihm+j), beta, c, s)
c
c     detect whether diagonal element of this column is zero
c
      if (abs(w(ihm+j)).lt.deps) then
         ipar(1) = -3
         goto 80
      end if
      w(ihd+j) = 1.0D0 / w(ihm+j)
      w(ic+j) = c
      w(is+j) = s
c
c     update the W's (the conjugate directions) -- essentially this is one
c     step triangular solve. Accumulate the influence from the previous
c     W's in w(1:n).
c
      ptrw = iw+(j-1)*n + 1
      if (full) then
         do i = j+1, lb
            alpha = -w(ihm+i)*w(ihd+i)
            do ii = 0, n-1
               w(ptrw+ii) = w(ptrw+ii) + alpha * w(iw+(i-1)*n+1+ii)
            end do
c            call daxpy1(n,alpha,w(iw+(i-1)*n+1),1,w(ptrw),1)
         end do
      end if
      do i = 1, j-1
         alpha = -w(ihm+i)*w(ihd+i)
         do ii = 0, n-1
            w(ptrw+ii) = w(ptrw+ii) + alpha * w(iw+(i-1)*n+1+ii)
         end do
c         call daxpy1(n,alpha,w(iw+(i-1)*n+1),1,w(ptrw),1)
      end do
c
c     update the solution to the linear system
c
      alpha = psi * c * w(ihd+j)
      psi = - s * psi
      do i = 1, n
         sol(i) = sol(i) + alpha * w(ptrw-1+i)
      end do
c      call daxpy1(n, alpha, w(ptrw), 1, sol, 1)
      fpar(11) = fpar(11) + min(ipar(7),lb) * 6 * (n+1)
c
c     determine whether to continue,
c     compute the desired error/residual norm
c
      ipar(7) = ipar(7) + 1
      fpar(5) = abs(psi)
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(10) = 6
         return
      end if
      if (ipar(3).lt.0) then
         alpha = abs(alpha)
         fpar(11) = fpar(11) + 2 * n
         if (ipar(7).eq.1 .and. ipar(3).eq.-1) then
            fpar(3) = alpha*sqrt(distdot(n, w(ptrw), 1, w(ptrw), 1))
            fpar(4) = fpar(1) * fpar(3) + fpar(2)
            fpar(6) = fpar(3)
            fpar(11) = fpar(11) + (n + n)
         else
            fpar(6) = alpha*sqrt(distdot(n, w(ptrw), 1, w(ptrw), 1))
         end if
      else
         fpar(6) = fpar(5)
      end if
      if (fpar(6).gt.fpar(4) .and. (ipar(6).le.0 .or.
     $     ipar(7).lt.ipar(6))) goto 30
 70   if (ipar(3).eq.999 .and. ipar(11).eq.0) goto 30
c
c     clean up the iterative solver
c
 80   fpar(7) = 0.0D0
      if (ipar(7).ne.0) fpar(7) = log10(fpar(3) / fpar(6)) / ipar(7)
      if (ipar(1).gt.0) then
         if (ipar(3).eq.999 .and. ipar(11).ne.0) then
            ipar(1) = 0
         else if (fpar(6).le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(6).gt.0 .and. ipar(7).ge.ipar(6)) then
            ipar(1) = -1
         else
            ipar(1) = -10
         end if
      end if
      return
      end
c-----end-of-dqgmres
c-----------------------------------------------------------------------
      subroutine fgmres(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
c-----------------------------------------------------------------------
c     This a version of FGMRES implemented with reverse communication.
c
c     ipar(5) == the dimension of the Krylov subspace
c
c     the space of the `w' is used as follows:
c     >> the bases for the Krylov subspace, size n*(m+1);
c     >> the above bases after (left-)multiplying with
c     the right-preconditioner invers, size m*n;
c     >> a temporary vector of size n;
c     >> the Henssenberg matrix, only the upper triangular
c     portion of the matrix is stored, size (m+1)*m/2 + 1
c     >> three vectors, first are of size m, they are
c     the cosine and sine of the Given rotations, the third one helds
c     the residuals, it is of size m+1.
c
c     TOTAL SIZE REQUIRED == n*(2m+1) + (m+1)*m/2 + 3*m + 2
c     Note: m == ipar(5). The default value for this is 15 if
c     ipar(5) <= 1.
c-----------------------------------------------------------------------
c     external functions used
c
      real*8 distdot
      external distdot
c
c     local variables, ptr and p2 are temporary pointers,
c     hens points to the Hessenberg matrix,
c     vc, vs point to the cosines and sines of the Givens rotations
c     vrn points to the vectors of residual norms, more precisely
c     the right hand side of the least square problem solved.
c
      integer i,ii,idx,iz,j,k,m,ptr,p2,hens,vc,vs,vrn
      real*8 alpha, c, s, deps
      logical lp, rp
      parameter (deps=1.0D-33)
      save
c
c     check the status of the call
c
      if (ipar(1).eq.0) ipar(10) = 0
      goto (10, 20, 30, 40, 50, 60) ipar(10)
c
c     initialization
c
      if (ipar(5).le.1) then
         m = 15
      else
         m = ipar(5)
      end if
      idx = n * (m+1)
      iz = idx + n
      hens = iz + n*m
      vc = hens + (m+1) * m / 2 + 1
      vs = vc + m
      vrn = vs + m
      i = vrn + m + 1
      call bisinit(n,ipar,i,1,lp,rp,w(idx+1))
      if (ipar(1).lt.0) return
c
c     request for matrix vector multiplication A*x in the initialization
c
 100  ipar(1) = 1
      ipar(8) = n+1
      ipar(9) = 1
      ipar(10) = 1
      k = 0
      do ii = 1, n
         w(ii+n) = sol(ii)
      end do
c      call dcopy1(n,sol,1,w(1+n),1)
      fpar(11) = fpar(11) + n
      return
 10   if (lp) then
         do i = 1, n
            w(n+i) = rhs(i) - w(i)
         end do
         ipar(1) = 3
         ipar(10) = 2
         return
      else
         do i = 1, n
            w(i) = rhs(i) - w(i)
         end do
      end if
      fpar(11) = fpar(11) + n
c
 20   alpha = sqrt(distdot(n,w,1,w,1))
      if (ipar(7).eq.0) then
         if (abs(ipar(3)).eq.2) then
            fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
            fpar(11) = fpar(11) + 2*n
         else
            fpar(4) = fpar(1) * alpha + fpar(2)
         end if
         fpar(3) = alpha
      end if
      fpar(5) = alpha
      w(vrn+1) = alpha
      if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999) then
         ipar(1) = 0
         fpar(6) = alpha
         goto 300
      end if
      alpha = 1.0D0 / alpha
      do ii = 1, n
         w(ii) = w(ii) * alpha
      end do
c      call dscal1(n,alpha,w,1)
      fpar(11) = fpar(11) + n
c
c     request for a matrix vector multiplication -- regular iterations
c
 110  k = k + 1
      if (rp) then
         ipar(1) = 5
         ipar(8) = k*n - n + 1
         ipar(9) = iz + ipar(8)
         ipar(10) = 3
         return
      else
         do ii = 0, n-1
            w(iz+k*n-ii) = w(k*n-ii)
         end do
c         call dcopy1(n,w(k*n-n+1),1,w(iz+k*n-n+1),1)
      end if
c
 30   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = (k-1)*n + 1
      end if
      if (lp) then
         ipar(9) = idx + 1
      else
         ipar(9) = 1 + k*n
      end if
      ipar(10) = 4
      return
c
 40   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = k*n + 1
         ipar(10) = 5
         return
      end if
c
c     Modified Gram-Schmidt orthogonalization procedure
c     temporary pointer 'ptr' is pointing to the current column of the
c     Henssenberg matrix. 'p2' points to the new basis vector
c
 50   ptr = k * (k - 1) / 2 + hens
      p2 = ipar(9)
      ipar(7) = ipar(7) + 1
      do i = 1, k
         j = (i-1)*n + 1
         alpha = distdot(n,w(p2),1,w(j),1)
         do ii = 0, n-1
            w(p2+ii) = w(p2+ii) - alpha * w(j+ii)
         end do
c         call daxpy1(n,-alpha,w(j),1,w(p2),1)
         w(ptr+i) = alpha
      end do
      alpha = sqrt(distdot(n,w(p2),1,w(p2),1))
      if (alpha.le.deps) then
         k = k - 1
         goto 200
      end if
      w(ptr+i) = alpha
      alpha = 1.0D0 / alpha
      do ii = 0, n-1
         w(p2+ii) = w(p2+ii) * alpha
      end do
c      call dscal1(n,alpha,w(p2),1)
c
c     apply previous Givens rotations and generate a new one to eliminate
c     the subdiagonal element.
c
      p2 = ptr + 1
      do i = 1, k-1
         ptr = p2
         p2 = p2 + 1
         alpha = w(ptr)
         c = w(vc+i)
         s = w(vs+i)
         w(ptr) = c * alpha + s * w(p2)
         w(p2) = c * w(p2) - s * alpha
      end do
      call givens(w(p2), w(p2+1), c, s)
      w(vc+k) = c
      w(vs+k) = s
      p2 = vrn + k
      alpha = - s * w(p2)
      w(p2) = c * w(p2)
      w(p2+1) = alpha
      fpar(11) = fpar(11) + (4 * k + 3) * n + 6 * k
c
c     end of one Arnoldi iteration, alpha will store the estimated
c     residual norm at current stage
c
      alpha = abs(alpha)
      fpar(5) = alpha
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(10) = 6
         return
      end if
      if ((ipar(3).lt.0 .or. alpha.gt.fpar(4)) .and.
     +     k.lt.m) goto 110
 60   if (ipar(3).eq.999 .and. ipar(11).eq.0) goto 110
c
c     update the approximate solution, first solve the upper triangular
c     system, temporary pointer ptr points to the Henssenberg matrix,
c     p2 points to the right-hand-side (also the solution) of the system.
c
 200  ptr = hens + k * (k + 1 ) / 2
      w(p2) = w(p2) / w(ptr)
      do i = k-1, 1, -1
         ptr = ptr - i - 1
         do ii = 1, i
            w(vrn+ii) = w(vrn+ii) - w(p2) * w(ptr+ii)
         end do
c         call daxpy1(i, -w(p2), w(ptr), 1, w(vrn+1), 1)
         p2 = p2 - 1
         w(p2) = w(p2) / w(ptr)
      end do
c
      do i = 0, k-1
         do ii = 1, n
            sol(ii) = sol(ii) + w(p2)*w(iz+i*n+ii)
         end do
c         call daxpy1(n, w(p2), w(iz+i*n+1), 1, sol, 1)
         p2 = p2 + 1
      end do
      fpar(11) = fpar(11) + 2*(k+k-1)*n
c
c     process the complete stopping criteria
c
      if (ipar(3).lt.0) then
         if (ipar(7).le.m) then
            fpar(3) = abs(w(vrn+1))
            if (ipar(3).eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
         end if
         alpha = abs(w(vrn+k))
      end if
      fpar(6) = alpha
c
c     do we need to restart ?
c
      if ((ipar(7).lt.ipar(6) .or. ipar(6).le.0).and.
     $     ((ipar(3).eq.999.and.ipar(11).eq.0) .or.
     $     (ipar(3).ne.999.and.alpha.gt.fpar(4)))) goto 100
c
c     termination, set error code, compute convergence rate
c
      if (ipar(1).gt.0) then
         if (ipar(3).eq.999 .and. ipar(11).eq.1) then
            ipar(1) = 0
         else if (ipar(3).ne.999 .and. alpha.le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else
            ipar(1) = -10
         end if
      end if
 300  if (alpha.ne.0.0D0) then
         fpar(7) = log10(fpar(3) / alpha) / ipar(7)
      else
         fpar(7) = 0.0D0
      end if
      return
      end
c-----end-of-fgmres
c-----------------------------------------------------------------------
      subroutine dbcg (n,rhs,sol,ipar,fpar,w)
      implicit none
      integer n,ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,*)
c-----------------------------------------------------------------------
c Quasi GMRES method for spolving a linear
c system of equations a * sol = y.  double precision version.
c this version is without restarting and without preconditioning.
c paramaters :
c -----------
c n     = dimension of the problem
c
c y     = w(:,1) a temporary storage used for various operations
c z     = w(:,2) a work vector of length n.
c v     = w(:,3:4)  size n x 2
c w     = w(:,5:6) size n x 2
c p     = w(:,7:9) work array of dimension n x 3
c
c sol   = the solution of the problem . at input sol must contain an
c         initial guess to the solution.
c    ***  note:   y is destroyed on return.
c
c-----------------------------------------------------------------------
c subroutines and functions called:
c 1) matrix vector multiplication and preconditioning through reverse
c     communication
c
c 2) implu, uppdir, distdot (blas)
c-----------------------------------------------------------------------
c aug. 1983  version.    author youcef saad. yale university computer
c science dept. some  changes made july 3, 1986.
c references: siam j. sci. stat. comp., vol. 5, pp. 203-228 (1984)
c-----------------------------------------------------------------------
c     local variables
c
      real*8 t,sqrt,distdot,ss,res,beta,ss1,delta,x,zeta,umm
      integer k,j,i,i2,ip2,ju,lb,lbm1,np,indp
      logical lp,rp,full, perm(3)
      real*8 ypiv(3),u(3),usav(3)
      external tidecg
      save
c
c     where to go
c
      if (ipar(1).eq.0) ipar(10) = 0
      goto (110, 120, 130, 140, 150, 160, 170, 180, 190, 200) ipar(10)
      call bisinit(n,ipar,11*n,1,lp,rp,w(1,10))
      if (ipar(1).lt.0) return
c-----------------------------------------------------------------------
c     initialize constants for outer loop :
c-----------------------------------------------------------------------
      lb = 3
      lbm1 = 2
c
c     get initial residual vector and norm
c
      ipar(1) = 1
      ipar(8) = 1
      ipar(9) = 1 + n
      do i = 1, n
         w(i,1) = sol(i)
      end do
c      call dcopy1(n,sol,1,w(1,1),1)
      ipar(10) = 1
      return
 110  if (lp) then
         do i = 1, n
            w(i,1) = rhs(i) - w(i,2)
         end do
         ipar(1) = 3
         ipar(8) = 1
         ipar(9) = n+n+1
         ipar(10) = 2
         return
      else
         do i = 1, n
            w(i,3) = rhs(i) - w(i,2)
         end do
      end if
      fpar(11) = fpar(11) + n
c
 120  fpar(3) = sqrt(distdot(n,w(1,3),1,w(1,3),1))
      fpar(11) = fpar(11) + n + n
      fpar(5) = fpar(3)
      fpar(7) = fpar(3)
      zeta = fpar(3)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + 2*n
      else
         fpar(4) = fpar(1) * zeta + fpar(2)
      end if
      if (ipar(3).ge.0.and.fpar(5).le.fpar(4)) then
         fpar(6) = fpar(5)
         goto 900
      end if
c
c     normalize first arnoldi vector
c
      t = 1.0d0/zeta
      do 22 k=1,n
         w(k,3) = w(k,3)*t
         w(k,5) = w(k,3)
 22   continue
      fpar(11) = fpar(11) + n
c
c     initialize constants for main loop
c
      beta = 0.0D0
      delta = 0.0D0
      i2 = 1
      indp = 0
      i = 0
c
c     main loop: i = index of the loop.
c
c-----------------------------------------------------------------------
 30   i = i + 1
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = (1+i2)*n+1
         if (lp) then
            ipar(9) = 1
         else
            ipar(9) = 10*n + 1
         end if
         ipar(10) = 3
         return
      end if
c
 130  ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = (1+i2)*n + 1
      end if
      if (lp) then
         ipar(9) = 10*n + 1
      else
         ipar(9) = 1
      end if
      ipar(10) = 4
      return
c
 140  if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 1
         ipar(10) = 5
         return
      end if
c
c     A^t * x
c
 150  ipar(7) = ipar(7) + 1
      if (lp) then
         ipar(1) = 4
         ipar(8) = (3+i2)*n + 1
         if (rp) then
            ipar(9) = n + 1
         else
            ipar(9) = 10*n + 1
         end if
         ipar(10) = 6
         return
      end if
c
 160  ipar(1) = 2
      if (lp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = (3+i2)*n + 1
      end if
      if (rp) then
         ipar(9) = 10*n + 1
      else
         ipar(9) = n + 1
      end if
      ipar(10) = 7
      return
c
 170  if (rp) then
         ipar(1) = 6
         ipar(8) = ipar(9)
         ipar(9) = n + 1
         ipar(10) = 8
         return
      end if
c-----------------------------------------------------------------------
c     orthogonalize current v against previous v's and
c     determine relevant part of i-th column of u(.,.) the
c     upper triangular matrix --
c-----------------------------------------------------------------------
 180  ipar(7) = ipar(7) + 1
      u(1) = 0.0d0
      ju = 1
      k = i2
      if (i .le. lbm1) ju = 0
      if (i .lt. lb) k = 0
 31   if (k .eq. lbm1) k=0
      k=k+1
c
      if (k .ne. i2) then
         ss  = delta
         ss1 = beta
         ju = ju + 1
         u(ju) = ss
      else
         ss = distdot(n,w(1,1),1,w(1,4+k),1)
         fpar(11) = fpar(11) + 2*n
         ss1= ss
         ju = ju + 1
         u(ju) = ss
      endif
c
      do 32  j=1,n
         w(j,1) = w(j,1) - ss*w(j,k+2)
         w(j,2) = w(j,2) - ss1*w(j,k+4)
 32   continue
      fpar(11) = fpar(11) + 4*n
c
      if (k .ne. i2) goto 31
c
c     end of Mod. Gram. Schmidt loop
c
      t = distdot(n,w(1,2),1,w(1,1),1)
c
      beta   = sqrt(abs(t))
      delta  = t/beta
c
      ss = 1.0d0/beta
      ss1 = 1.0d0/ delta
c
c     normalize and insert new vectors
c
      ip2 = i2
      if (i2 .eq. lbm1) i2=0
      i2=i2+1
c
      do 315 j=1,n
         w(j,i2+2)=w(j,1)*ss
         w(j,i2+4)=w(j,2)*ss1
 315  continue
      fpar(11) = fpar(11) + 4*n
c-----------------------------------------------------------------------
c     end of orthogonalization.
c     now compute the coefficients u(k) of the last
c     column of the  l . u  factorization of h .
c-----------------------------------------------------------------------
      np = min0(i,lb)
      full = (i .ge. lb)
      call implu(np, umm, beta, ypiv, u, perm, full)
c-----------------------------------------------------------------------
c     update conjugate directions and solution
c-----------------------------------------------------------------------
      do 33 k=1,n
         w(k,1) = w(k,ip2+2)
 33   continue
      call uppdir(n, w(1,7), np, lb, indp, w, u, usav, fpar(11))
c-----------------------------------------------------------------------
      if (i .eq. 1) goto 34
      j = np - 1
      if (full) j = j-1
      if (.not.perm(j)) zeta = -zeta*ypiv(j)
 34   x = zeta/u(np)
      if (perm(np))goto 36
      do 35 k=1,n
         w(k,10) = w(k,10) + x*w(k,1)
 35   continue
      fpar(11) = fpar(11) + 2 * n
c-----------------------------------------------------------------------
 36   if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(10) = 9
         return
      end if
      res = abs(beta*zeta/umm)
      fpar(5) = res * sqrt(distdot(n, w(1,i2+2), 1, w(1,i2+2), 1))
      fpar(11) = fpar(11) + 2 * n
      if (ipar(3).lt.0) then
         fpar(6) = x * sqrt(distdot(n,w,1,w,1))
         fpar(11) = fpar(11) + 2 * n
         if (ipar(7).le.2) then
            fpar(3) = fpar(6)
            if (ipar(3).eq.-1) then
               fpar(4) = fpar(1) * sqrt(fpar(3)) + fpar(2)
            end if
         end if
      else
         fpar(6) = fpar(5)
      end if
c---- convergence test -----------------------------------------------
 190  if (ipar(3).eq.999.and.ipar(11).eq.0) then
         goto 30
      else if (fpar(6).gt.fpar(4) .and. (ipar(6).gt.ipar(7) .or.
     +        ipar(6).le.0)) then
         goto 30
      end if
c-----------------------------------------------------------------------
c     here the fact that the last step is different is accounted for.
c-----------------------------------------------------------------------
      if (.not. perm(np)) goto 900
      x = zeta/umm
      do 40 k = 1,n
         w(k,10) = w(k,10) + x*w(k,1)
 40   continue
      fpar(11) = fpar(11) + 2 * n
c
c     right preconditioning and clean-up jobs
c
 900  if (rp) then
         ipar(1) = 5
         ipar(8) = 9*n + 1
         ipar(9) = ipar(8) + n
         ipar(10) = 10
         return
      end if
 200  if (rp) then
         call tidecg(n,ipar,fpar,sol,w(1,11))
      else
         call tidecg(n,ipar,fpar,sol,w(1,10))
      end if
      return
c-----------------------------------------------------------------------
c-------end-of-dbcg-----------------------------------------------------
      end
c
      subroutine implu(np,umm,beta,ypiv,u,permut,full)
      real*8 umm,beta,ypiv(*),u(*),x, xpiv
      logical full, perm, permut(*)
      integer np,k,npm1
c-----------------------------------------------------------------------
c     performs implicitly one step of the lu factorization of a
c     banded hessenberg matrix.
c-----------------------------------------------------------------------
      if (np .le. 1) goto 12
      npm1 = np - 1
c
c     -- perform  previous step of the factorization-
c
      do 6 k=1,npm1
         if(.not. permut(k))goto 5
         x=u(k)
         u(k) = u(k+1)
         u(k+1) = x
 5       u(k+1) = u(k+1) - ypiv(k)*u(k)
 6    continue
c-----------------------------------------------------------------------
c     now determine pivotal information to be used in the next call
c-----------------------------------------------------------------------
 12   umm = u(np)
      perm = (beta .gt. abs(umm))
      if (.not. perm) goto 4
      xpiv = umm / beta
      u(np) = beta
      goto 8
 4    xpiv = beta/umm
 8    permut(np) = perm
      ypiv(np) = xpiv
      if (.not. full) return
c     shift everything up if full...
      do 7 k=1,npm1
         ypiv(k) = ypiv(k+1)
         permut(k) = permut(k+1)
 7    continue
      return
c-----end-of-implu
      end
c-----------------------------------------------------------------------
      subroutine uppdir(n,p,np,lbp,indp,y,u,usav,flops)
      real*8 p(n,lbp), y(*), u(*), usav(*), x, flops
      integer k,np,n,npm1,j,ju,indp,lbp
c-----------------------------------------------------------------------
c     updates the conjugate directions p given the upper part of the
c     banded uppertriangular matrix u.  u contains the non zero
c     elements of the column of the triangular matrix..
c-----------------------------------------------------------------------
      npm1=np-1
      if (np .le. 1) goto 12
      j=indp
      ju = npm1
 10   if (j .le. 0) j=lbp
      x = u(ju) /usav(j)
      if (x .eq. 0.0d0) goto 115
      do 11 k=1,n
         y(k) = y(k) - x*p(k,j)
 11   continue
      flops = flops + 2*n
 115  j = j-1
      ju = ju -1
      if (ju .ge. 1) goto 10
 12   indp = indp + 1
      if (indp .gt. lbp) indp = 1
      usav(indp) = u(np)
      do 13 k=1,n
         p(k,indp) = y(k)
 13   continue
 208  return
c-----------------------------------------------------------------------
c-------end-of-uppdir---------------------------------------------------
      end
      subroutine givens(x,y,c,s)
      real*8 x,y,c,s
c-----------------------------------------------------------------------
c     Given x and y, this subroutine generates a Givens' rotation c, s.
c     And apply the rotation on (x,y) ==> (sqrt(x**2 + y**2), 0).
c     (See P 202 of "matrix computation" by Golub and van Loan.)
c-----------------------------------------------------------------------
      real*8 t
c
      if (x.eq.0.0D0 .and. y.eq.0.0D0) then
         c = 1.0D0
         s = 0.0D0
      else if (abs(y).gt.abs(x)) then
         t = x / y
         x = sqrt(1.0d0+t*t)
         s = sign(1.0D0 / x, y)
         c = t*s
      else
         t = y / x
         y = sqrt(1.0d0+t*t)
         c = sign(1.0D0 / y, x)
         s = t*c
      end if
      x = abs(x*y)
c
c     end of givens
c
      return
      end
c-----end-of-givens
c-----------------------------------------------------------------------
      logical function stopbis(n,ipar,mvpi,fpar,r,delx,sx)
      implicit none
      integer n,mvpi,ipar(16)
      real*8 fpar(16), r(n), delx(n), sx, distdot
      external distdot
c-----------------------------------------------------------------------
c     function for determining the stopping criteria. return value of
c     true if the stopbis criteria is satisfied.
c-----------------------------------------------------------------------
      if (ipar(11) .eq. 1) then
         stopbis = .true.
      else
         stopbis = .false.
      end if
      if (ipar(6).gt.0 .and. ipar(7).ge.ipar(6)) then
         ipar(1) = -1
         stopbis = .true.
      end if
      if (stopbis) return
c
c     computes errors
c
      fpar(5) = sqrt(distdot(n,r,1,r,1))
      fpar(11) = fpar(11) + 2 * n
      if (ipar(3).lt.0) then
         fpar(6) = sx * sqrt(distdot(n,delx,1,delx,1))
         fpar(11) = fpar(11) + 2 * n
         if (ipar(7).eq.mvpi) then
            fpar(3) = fpar(6)
            if (ipar(3).eq.-1) then
               fpar(4) = fpar(1) * sqrt(fpar(3)) + fpar(2)
            end if
         end if
      else
         fpar(6) = fpar(5)
      end if
c
      if (fpar(6).le.fpar(4)) then
         stopbis = .true.
         ipar(11) = 1
      end if
c
      return
      end
c-----end-of-stopbis
c-----------------------------------------------------------------------
      subroutine tidecg(n,ipar,fpar,sol,delx)
      implicit none
      integer i,n,ipar(16)
      real*8 fpar(16),sol(n),delx(n)
c-----------------------------------------------------------------------
c     Some common operations required before terminating the CG routines
c-----------------------------------------------------------------------
      if (ipar(1).gt.0) then
         if ((ipar(3).eq.999 .and. ipar(11).eq.1) .or.
     +        fpar(6).le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else
            ipar(1) = -10
         end if
      end if
      if (fpar(6).lt.0.0D0) then
         fpar(7) = 0.0D0
      else
         fpar(7) = log10(fpar(3) / fpar(6)) / ipar(7)
      end if
      do i = 1, n
         sol(i) = sol(i) + delx(i)
      end do
      return
      end
c-----end-of-tidecg
c-----------------------------------------------------------------------
      logical function brkdn(alpha, ipar)
      implicit none
      integer ipar(16)
      real*8 alpha, eps
      parameter (eps=1.0D-33)
c-----------------------------------------------------------------------
c     test whether alpha is smaller than eps in magnitude. If yes,
c     this routine will return .true.
c-----------------------------------------------------------------------
      if (abs(alpha).lt.eps) then
         brkdn = .true.
         ipar(1) = -3
      else
         brkdn = .false.
      end if
      return
      end
c-----end-of-brkdn
c-----------------------------------------------------------------------
      subroutine bisinit(n,ipar,wksize,dsc,lp,rp,delx)
      implicit none
      integer i,n,ipar(16),wksize,dsc
      logical lp,rp
      real*8 delx(n)
c-----------------------------------------------------------------------
c     some common initializations for the BIS (basic iterative solvers)
c-----------------------------------------------------------------------
      if (ipar(4).lt.wksize) then
         ipar(1) = -2
         ipar(4) = wksize
         return
      endif
c
      if (ipar(2).gt.2) then
         lp = .true.
         rp = .true.
      else if (ipar(2).eq.2) then
         lp = .false.
         rp = .true.
      else if (ipar(2).eq.1) then
         lp = .true.
         rp = .false.
      else
         lp = .false.
         rp = .false.
      endif
      if (ipar(3).eq.0) ipar(3) = dsc
      ipar(7) = 0
      ipar(11) = 0
      do i = 1, n
         delx(i) = 0.0D0
      end do
c
      return
      end
c-----end-of-bisinit
c-----------------------------------------------------------------------
