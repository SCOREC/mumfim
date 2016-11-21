c----------------------------------------------------------------------c
       subroutine pgmres (n, im, rhs, sol, vv, eps, maxits, iout,
     *                    aa, ja, ia, alu, jlu, ju, ierr)
c-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       integer n, im, maxits, iout, ierr, ja(*), ia(n+1), jlu(*), ju(n)
       real*8 vv(n,*), rhs(n), sol(n), aa(*), alu(*), eps
	common /statu/its
c----------------------------------------------------------------------*
c                                                                      *
c                 *** Preconditioned GMRES ***                         *
c                                                                      *
c----------------------------------------------------------------------*
c----------------------------------------------------------------------*
c Coded by Y. Saad - This version dated May, 7, 1990.                  *
c----------------------------------------------------------------------*
c parameters                                                           *
c-----------                                                           *
c on entry:                                                            *
c==========                                                            *
c                                                                      *
c n     == integer. The dimension of the matrix.                       *
c im    == size of krylov subspace:  should not exceed 50 in this      *
c          version (can be reset by changing parameter command for     *
c          kmax below)                                                 *
c rhs   == real vector of length n containing the right hand side.     *
c          Destroyed on return.                                        *
c sol   == real vector of length n containing an initial guess to the  *
c          solution on input. approximate solution on output           *
c eps   == tolerance for stopping criterion. process is stopped        *
c          as soon as ( ||.|| is the euclidean norm):                  *
c          || current residual||/||initial residual|| <= eps           *
c maxits== maximum number of iterations allowed                        *
c iout  == output unit number number for printing intermediate results *
c          if (iout .le. 0) nothing is printed out.                    *
c                                                                      *
c aa, ja,                                                              *
c ia    == the input matrix in compressed sparse row format:           *
c          aa(1:nnz)  = nonzero elements of A stored row-wise in order *
c          ja(1:nnz) = corresponding column indices.                   *
c          ia(1:n+1) = pointer to beginning of each row in aa and ja.  *
c          here nnz = number of nonzero elements in A = ia(n+1)-ia(1)  *
c                                                                      *
c alu,jlu== A matrix stored in Modified Sparse Row format containing   *
c           the L and U factors, as computed by subroutine ilut.       *
c                                                                      *
c ju     == integer array of length n containing the pointers to       *
c           the beginning of each row of U in alu, jlu as computed     *
c           by subroutine ILUT.                                        *
c                                                                      *
c on return:                                                           *
c==========                                                            *
c sol   == contains an approximate solution (upon successful return).  *
c ierr  == integer. Error message with the following meaning.          *
c          ierr = 0 --> successful return.                             *
c          ierr = 1 --> convergence not achieved in itmax iterations.  *
c          ierr =-1 --> the initial guess seems to be the exact        *
c                       solution (initial residual computed was zero)  *
c                                                                      *
c----------------------------------------------------------------------*
c                                                                      *
c work arrays:                                                         *
c=============                                                         *
c vv    == work array of length  n x (im+1) (used to store the Arnoli  *
c          basis)                                                      *
c----------------------------------------------------------------------*
c subroutines called :                                                 *
c amux   : SPARSKIT routine to do the matrix by vector multiplication  *
c          delivers y=Ax, given x  -- see SPARSKIT/BLASSM/amux         *
c lusol : combined forward and backward solves (Preconditioning ope.) *
c BLAS1  routines.                                                     *
c----------------------------------------------------------------------*
       parameter (kmax=1100)
       real*8 hh(kmax+1,kmax), c(kmax), s(kmax), rs(kmax+1),t
c-------------------------------------------------------------
c arnoldi size should not exceed kmax=50 in this version..
c to reset modify paramter kmax accordingly.
c-------------------------------------------------------------
       data epsmac/1.d-16/
       n1 = n + 1
       its = 0
c-------------------------------------------------------------
c outer loop starts here..
c-------------- compute initial residual vector --------------
       call amux (n, sol, vv, aa, ja, ia)
       do 21 j=1,n
          vv(j,1) = rhs(j) - vv(j,1)
 21    continue
c-------------------------------------------------------------
 20    ro = dnrm3(n, vv, 1)
       if (iout .gt. 0 .and. its .eq. 0)
     *      write(iout, 199) its, ro
       if (ro .eq. 0.0d0) goto 999
       t = 1.0d0/ ro
       do 210 j=1, n
          vv(j,1) = vv(j,1)*t
 210   continue
       if (its .eq. 0) eps1=eps*ro
c     ** initialize 1-st term  of rhs of hessenberg system..
       rs(1) = ro
       i = 0
 4     i=i+1
       its = its + 1
       i1 = i + 1
       call lusol (n, vv(1,i), rhs, alu, jlu, ju)
       call amux (n, rhs, vv(1,i1), aa, ja, ia)
c-----------------------------------------
c     modified gram - schmidt...
c-----------------------------------------
       do 55 j=1, i
          t = ddot1(n, vv(1,j),1,vv(1,i1),1)
          hh(j,i) = t
          call daxpy1(n, -t, vv(1,j), 1, vv(1,i1), 1)
 55    continue
       t = dnrm3(n, vv(1,i1), 1)
       hh(i1,i) = t
       if ( t .eq. 0.0d0) goto 58
       t = 1.0d0/t
       do 57  k=1,n
          vv(k,i1) = vv(k,i1)*t
 57    continue
c
c     done with modified gram schimd and arnoldi step..
c     now  update factorization of hh
c
 58    if (i .eq. 1) goto 121
c--------perfrom previous transformations  on i-th column of h
       do 66 k=2,i
          k1 = k-1
          t = hh(k1,i)
          hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
          hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
 66    continue
 121   gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)
c
c     if gamma is zero then any small value will do...
c     will affect only residual estimate
c
       if (gam .eq. 0.0d0) gam = epsmac
c
c     get  next plane rotation
c
       c(i) = hh(i,i)/gam
       s(i) = hh(i1,i)/gam
       rs(i1) = -s(i)*rs(i)
       rs(i) =  c(i)*rs(i)
c
c     detrermine residual norm and test for convergence-
c
       hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
       ro = abs(rs(i1))
 131   format(1h ,2e14.4)
       if (iout .gt. 0)
     *      write(iout, 199) its, ro
       if (i .lt. im .and. (ro .gt. eps1))  goto 4
c
c     now compute solution. first solve upper triangular system.
c
       rs(i) = rs(i)/hh(i,i)
       do 30 ii=2,i
          k=i-ii+1
          k1 = k+1
          t=rs(k)
          do 40 j=k1,i
             t = t-hh(k,j)*rs(j)
 40       continue
          rs(k) = t/hh(k,k)
 30    continue
c
c     form linear combination of v(*,i)'s to get solution
c
       t = rs(1)
       do 15 k=1, n
          rhs(k) = vv(k,1)*t
 15    continue
       do 16 j=2, i
          t = rs(j)
          do 161 k=1, n
             rhs(k) = rhs(k)+t*vv(k,j)
 161      continue
 16    continue
c
c     call preconditioner.
c
       call lusol (n, rhs, rhs, alu, jlu, ju)
       do 17 k=1, n
          sol(k) = sol(k) + rhs(k)
 17    continue
c
c     restart outer loop  when necessary
c
       if (ro .le. eps1) goto 990
       if (its .gt. maxits) goto 991
c
c     else compute residual vector and continue..
c
       do 24 j=1,i
          jj = i1-j+1
          rs(jj-1) = -s(jj-1)*rs(jj)
          rs(jj) = c(jj-1)*rs(jj)
 24    continue
       do 25  j=1,i1
          t = rs(j)
          if (j .eq. 1)  t = t-1.0d0
          call daxpy1 (n, t, vv(1,j), 1,  vv, 1)
 25    continue
 199   format(i4, e20.6)
c old output format: format('   its =', i4, ' res. norm =', d20.6)
c     restart outer loop.
       goto 20
 990   ierr = 0
       return
 991   ierr = 1
       return
 999   continue
       ierr = -1
       return
c-----------------end of pgmres ---------------------------------------
c-----------------------------------------------------------------------
       end
