      subroutine runrc(n,rhs,sol,ipar,fpar,wk,guess,a,ja,ia,
     +     au,jau,ju,solver)
      implicit none
      integer n,ipar(16),ia(n+1),ja(*),ju(*),jau(*)
      real*8 fpar(16),rhs(n),sol(n),guess(n),wk(*),a(*),au(*)
c-----------------------------------------------------------------------
c     the actual tester. It starts the iterative linear system solvers
c     with a initial guess suppied by the user.
c
c     The structure {au, jau, ju} is assumed to have the output from
c     the ILU* routines in ilut.f.
c
c-----------------------------------------------------------------------
c     local variables
c
      integer solver
      integer i, iou, its
      real*8 res, dnrm3
      real dt(2), time
      external cg,bcg,dbcg,bcgstab,tfqmr,gmres,fgmres,dqgmres
      external cgnr
      external dnrm3
      save its,res
c
c     ipar(2) can be 0, 1, 2, please don't use 3
c
c      print *, 'ipar(2) is', ipar(2) 
      if (ipar(2).gt.2) then
         print *, 'I can not do both left and right preconditioning.'
         return
      end if
c
c     normal execution
c
      its = 0
      res = 0.0D0
c
      do i = 1, n
         sol(i) = guess(i)
      end do
c
      iou = 6
      ipar(1) = 0
 10   if(solver.eq.1) then
      call cg(n,rhs,sol,ipar,fpar,wk)
c      NICO
c      write(iou, 99) fpar(6)
c 99   format(e20.6)
c      end NICO      
      else if(solver.eq.2) then 
      call bcg(n,rhs,sol,ipar,fpar,wk)
      else if(solver.eq.3)then
      call dbcg(n,rhs,sol,ipar,fpar,wk)
      else if(solver.eq.4)then
      call bcgstab(n,rhs,sol,ipar,fpar,wk)
      else if(solver.eq.5)then
      call tfqmr(n,rhs,sol,ipar,fpar,wk)
      else if(solver.eq.6)then
      call gmres(n,rhs,sol,ipar,fpar,wk)
      else if(solver.eq.7)then
      call fgmres(n,rhs,sol,ipar,fpar,wk)
      else if(solver.eq.8)then
      call dqgmres(n,rhs,sol,ipar,fpar,wk)
      else if(solver.eq.9)then
      call cgnr(n,rhs,sol,ipar,fpar,wk)
      endif
c
c     output the residuals
c
      if (ipar(7).ne.its) then
c NICO
c      write (iou, *) its, real(res)

      its = ipar(7)
      end if
      res = fpar(5)
c
c NICO
c      print *,'ipar(1) is', ipar(1)

      if (ipar(1).eq.1) then
         call amux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
         goto 10
      else if (ipar(1).eq.2) then
         call atmux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
         goto 10
      else if (ipar(1).eq.3 .or. ipar(1).eq.5) then
         call lusol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
         goto 10
      else if (ipar(1).eq.4 .or. ipar(1).eq.6) then
         call lutsol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
         goto 10
      end if
c      time = dtime(dt)
c  nico
c      write (iou, *) ipar(7), real(fpar(6))
c      write (iou, *) '# retrun code =', ipar(1),
c     +     '	convergence rate =', fpar(7)
c      write (iou, *) '# total execution time (sec)', time
c
c     check the error
c
      call amux(n,sol,wk,a,ja,ia)
      do i = 1, n
         wk(n+i) = sol(i) -1.0D0
         wk(i) = wk(i) - rhs(i)
      end do
c      write (iou, *) '# the actual residual norm is', dnrm3(n,wk,1)
c      write (iou, *) '# the error norm is', dnrm3(n,wk(1+n),1)
c
      if (iou.ne.6) close(iou)
      return
      end
c-----end-of-runrc
c-----------------------------------------------------------------------
      real*8 function distdot(n,x,ix,y,iy)
      integer n, ix, iy
      real*8 x(1+(n-1)*ix), y(1+(n-1)*iy), ddot1
      external ddot1
      distdot = ddot1(n,x,ix,y,iy)
      return
      end
c-----end-of-distdot
c-----------------------------------------------------------------------
c
      function afun (x,y,z)
      real*8 afun, x,y, z 
      afun = -1.0D0
      return 
      end
      
      function bfun (x,y,z)
      real*8 bfun, x,y, z 
      bfun = -1.0D0
      return 
      end
      
      function cfun (x,y,z)
      real*8 cfun, x,y, z 
      cfun = -1.0D0
      return 
      end
      
      function dfun (x,y,z)
      real*8 dfun, x,y, z, gammax, gammay, alpha
      common /func/ gammax, gammay, alpha
      dfun = gammax*exp(x*y)
      return 
      end
      
      function efun (x,y,z)
      real*8 efun, x,y, z, gammax, gammay, alpha
      common /func/ gammax, gammay, alpha
      efun = gammay*exp(-x*y) 
      return 
      end
      
      function ffun (x,y,z)
      real*8 ffun, x,y, z 
      ffun = 0.0D0
      return 
      end
      
      function gfun (x,y,z)
      real*8 gfun, x,y, z, gammax, gammay, alpha
      common /func/ gammax, gammay, alpha
      gfun = alpha 
      return 
      end
      
      function hfun (x,y,z)
      real*8 hfun, x,y, z, gammax, gammay, alpha
      common /func/ gammax, gammay, alpha
      hfun = alpha * sin(gammax*x+gammay*y-z)
      return 
      end
      
      function betfun(x,y,z)
      real*8 betfun,x,y,z
      betfun = 1.0D0
      return
      end
      
      function gamfun(x,y,z)
      real*8 gamfun,x,y,z
      gamfun = 0.0D0
      return
      end
c-----------------------------------------------------------------------
c     functions for the block PDE's 
c-----------------------------------------------------------------------
      subroutine afunbl (nfree,x,y,z,coeff)
      return
      end
c     
      subroutine bfunbl (nfree,x,y,z,coeff)
      return 
      end
      
      subroutine cfunbl (nfree,x,y,z,coeff)
c     
      return 
      end
      
      subroutine dfunbl (nfree,x,y,z,coeff)
      
      return
      end
c     
      subroutine efunbl (nfree,x,y,z,coeff)
      return 
      end
c     
      subroutine ffunbl (nfree,x,y,z,coeff)
      return 
      end
c     
      subroutine gfunbl (nfree,x,y,z,coeff)
      return 
      end
      
