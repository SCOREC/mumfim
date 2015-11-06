      subroutine myplot (ncol,ja,ia,iunt,mode)

      integer ja(*),ia(*),iunt,ncol,id,in,mode

      call pspltm ( ncol , ncol , mode , ja , ia , ' ' ,
     &     0,5.0,in,0,id,iunt)
      return
      end

c-----------------------------------------------------------------------
      subroutine pspltm(nrow,ncol,mode,ja,ia,title,ptitle,size,munt,
     *     nlines,lines,iunt)
c-----------------------------------------------------------------------
      integer nrow,ncol,ptitle,mode,iunt, ja(*), ia(*), lines(nlines) 
      real size
      character title*(*), munt*2 
c----------------------------------------------------------------------- 
c PSPLTM - PostScript PLoTer of a (sparse) Matrix
c This version by loris renggli (renggli@masg1.epfl.ch), Dec 1991
c and Youcef Saad 
c------
c Loris RENGGLI, Swiss Federal Institute of Technology, Math. Dept
c CH-1015 Lausanne (Switzerland)  -- e-mail:  renggli@masg1.epfl.ch
c Modified by Youcef Saad -- June 24, 1992 to add a few features:
c separation lines + acceptance of MSR format.
c-----------------------------------------------------------------------
c input arguments description :
c
c nrow   = number of rows in matrix
c
c ncol   = number of columns in matrix 
c
c mode   = integer indicating whether the matrix is stored in 
c           CSR mode (mode=0) or CSC mode (mode=1) or MSR mode (mode=2) 
c
c ja     = column indices of nonzero elements when matrix is
c          stored rowise. Row indices if stores column-wise.
c ia     = integer array of containing the pointers to the 
c          beginning of the columns in arrays a, ja.
c
c title  = character*(*). a title of arbitrary length to be printed 
c          as a caption to the figure. Can be a blank character if no
c          caption is desired.
c
c ptitle = position of title; 0 under the drawing, else above
c
c size   = size of the drawing  
c
c munt   = units used for size : 'cm' or 'in'
c
c nlines = number of separation lines to draw for showing a partionning
c          of the matrix. enter zero if no partition lines are wanted.
c
c lines  = integer array of length nlines containing the coordinates of 
c          the desired partition lines . The partitioning is symmetric: 
c          a horizontal line across the matrix will be drawn in 
c          between rows lines(i) and lines(i)+1 for i=1, 2, ..., nlines
c          an a vertical line will be similarly drawn between columns
c          lines(i) and lines(i)+1 for i=1,2,...,nlines 
c
c iunt   = logical unit number where to write the matrix into.
c----------------------------------------------------------------------- 
c additional note: use of 'cm' assumes european format for paper size
c (21cm wide) and use of 'in' assumes american format (8.5in wide).
c The correct centering of the figure depends on the proper choice. Y.S.
c-----------------------------------------------------------------------
c external 
      integer LENSTR
      external LENSTR
c local variables ---------------------------------------------------
      integer n,nr,nc,maxdim,istart,ilast,ii,k,ltit
      real lrmrgn,botmrgn,xtit,ytit,ytitof,fnstit,siz
      real xl,xr, yb,yt, scfct,u2dot,frlw,delt,paperx,conv,xx,yy
      logical square 
c change square to .true. if you prefer a square frame around
c a rectangular matrix
      data haf /0.5/, zero/0.0/, conv/2.54/,square/.false./
c-----------------------------------------------------------------------
      siz = size
      nr = nrow
      nc = ncol
      n = nc
      if (mode .eq. 0) n = nr
c      nnz = ia(n+1) - ia(1) 
      maxdim = max(nrow, ncol)
      m = 1 + maxdim
      nc = nc+1
      nr = nr+1
c
c units (cm or in) to dot conversion factor and paper size
c 
      if (munt.eq.'cm' .or. munt.eq.'CM') then
         u2dot = 72.0/conv
        paperx = 21.0
      else
        u2dot = 72.0
        paperx = 8.5*conv
        siz = siz*conv
      end if
c
c left and right margins (drawing is centered)
c 
      lrmrgn = (paperx-siz)/2.0
c
c bottom margin : 2 cm
c
      botmrgn = 2.0
c scaling factor
      scfct = siz*u2dot/m
c matrix frame line witdh
      frlw = 0.25
c font size for title (cm)
      fnstit = 0.5
      ltit = LENSTR(title)
c position of title : centered horizontally
c                     at 1.0 cm vertically over the drawing
      ytitof = 1.0
      xtit = paperx/2.0
      ytit = botmrgn+siz*nr/m + ytitof
c almost exact bounding box
      xl = lrmrgn*u2dot - scfct*frlw/2
      xr = (lrmrgn+siz)*u2dot + scfct*frlw/2
      yb = botmrgn*u2dot - scfct*frlw/2
      yt = (botmrgn+siz*nr/m)*u2dot + scfct*frlw/2
      if (ltit.gt.0) then
        yt = yt + (ytitof+fnstit*0.70)*u2dot
      end if
c add some room to bounding box
      delt = 10.0
      xl = xl-delt
      xr = xr+delt
      yb = yb-delt
      yt = yt+delt
c
c correction for title under the drawing
      if (ptitle.eq.0 .and. ltit.gt.0) then
      ytit = botmrgn + fnstit*0.3
      botmrgn = botmrgn + ytitof + fnstit*0.7
      end if
c begin of output
c
      write(iunt,10) '%!'
      write(iunt,10) '%%Creator: PSPLTM routine'
      write(iunt,12) '%%BoundingBox:',xl,yb,xr,yt
      write(iunt,10) '%%EndComments'
      write(iunt,10) '/cm {72 mul 2.54 div} def'
      write(iunt,10) '/mc {72 div 2.54 mul} def'
      write(iunt,10) '/pnum { 72 div 2.54 mul 20 string'
      write(iunt,10) 'cvs print ( ) print} def'
      write(iunt,10)
     1  '/Cshow {dup stringwidth pop -2 div 0 rmoveto show} def'
c
c we leave margins etc. in cm so it is easy to modify them if
c needed by editing the output file
      write(iunt,10) 'gsave'
      if (ltit.gt.0) then
      write(iunt,*) '/Helvetica findfont ',fnstit,
     &             ' cm scalefont setfont '
      write(iunt,*) xtit,' cm ',ytit,' cm moveto '
      write(iunt,'(3A)') '(',title(1:ltit),') Cshow'
      end if
      write(iunt,*) lrmrgn,' cm ',botmrgn,' cm translate'
      write(iunt,*) siz,' cm ',m,' div dup scale '
c------- 
c draw a frame around the matrix
      write(iunt,*) frlw,' setlinewidth'
      write(iunt,10) 'newpath'
      write(iunt,11) 0, 0, ' moveto'
      if (square) then
      write(iunt,11) m,0,' lineto'
      write(iunt,11) m, m, ' lineto'
      write(iunt,11) 0,m,' lineto'
      else
      write(iunt,11) nc,0,' lineto'
      write(iunt,11) nc,nr,' lineto'
      write(iunt,11) 0,nr,' lineto'
      end if
      write(iunt,10) 'closepath stroke'
c
c     drawing the separation lines 
c 
      write(iunt,*)  ' 0.2 setlinewidth'
      do 22 kol=1, nlines 
         isep = lines(kol) 
c
c     horizontal lines 
c
         yy =  real(nrow-isep) + haf 
         xx = real(ncol+1) 
         write(iunt,13) zero, yy, ' moveto '
         write(iunt,13)  xx, yy, ' lineto stroke '
c
c vertical lines 
c
         xx = real(isep) + haf 
         yy = real(nrow+1)  
         write(iunt,13) xx, zero,' moveto '
         write(iunt,13) xx, yy, ' lineto stroke '             
 22     continue
c 
c----------- plotting loop ---------------------------------------------
c
      write(iunt,10) '1 1 translate'
      write(iunt,10) '0.8 setlinewidth'
      write(iunt,10) '/p {moveto 0 -.40 rmoveto '
      write(iunt,10) '           0  .80 rlineto stroke} def'
c     
      do 1 ii=1, n
        istart = ia(ii)
        ilast  = ia(ii+1)-1 
        if (mode .eq. 1) then
          do 2 k=istart, ilast
            write(iunt,11) ii-1, nrow-ja(k), ' p'
 2        continue 
        else
          do 3 k=istart, ilast
            write(iunt,11) ja(k)-1, nrow-ii, ' p'
 3        continue          
c add diagonal element if MSR mode.
          if (mode .eq. 2) 
     *         write(iunt,11) ii-1, nrow-ii, ' p' 
c
        endif
 1    continue
c-----------------------------------------------------------------------
      write(iunt,10) 'showpage'
      return
c
 10   format (A)
 11   format (2(I6,1x),A)
 12   format (A,4(1x,F9.2))
 13   format (2(F9.2,1x),A)
c-----------------------------------------------------------------------
      end
c
      integer function lenstr(s)
c-----------------------------------------------------------------------
c return length of the string S
c-----------------------------------------------------------------------
      character*(*) s
      integer len
      intrinsic len
      integer n
c----------------------------------------------------------------------- 
      n = len(s)
10    continue
        if (s(n:n).eq.' ') then
          n = n-1
          if (n.gt.0) go to 10
        end if
      lenstr = n
c
      return
c--------end-of-pspltm--------------------------------------------------
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
c----------------------------------------------------------------------- 
      subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),job
      real*8 a(*),ao(*) 
c-----------------------------------------------------------------------
c this subroutine permutes the rows of a matrix in CSR format. 
c rperm  computes B = P A  where P is a permutation matrix.  
c the permutation P is defined through the array perm: for each j, 
c perm(j) represents the destination row number of row number j. 
c Youcef Saad -- recoded Jan 28, 1991.
c-----------------------------------------------------------------------
c on entry:
c----------
c n 	= dimension of the matrix
c a, ja, ia = input matrix in csr format
c perm 	= integer array of length nrow containing the permutation arrays
c	  for the rows: perm(i) is the destination of row i in the
c         permuted matrix. 
c         ---> a(i,j) in the original matrix becomes a(perm(i),j) 
c         in the output  matrix.
c
c job	= integer indicating the work to be done:
c 		job = 1	permute a, ja, ia into ao, jao, iao 
c                       (including the copying of real values ao and
c                       the array iao).
c 		job .ne. 1 :  ignore real values.
c                     (in which case arrays a and ao are not needed nor
c                      used).
c
c------------
c on return: 
c------------ 
c ao, jao, iao = input matrix in a, ja, ia format
c note : 
c        if (job.ne.1)  then the arrays a and ao are not used.
c----------------------------------------------------------------------c
c           Y. Saad, May  2, 1990                                      c
c----------------------------------------------------------------------c
      logical values
      values = (job .eq. 1) 
c     
c     determine pointers for output matix. 
c     
      do 50 j=1,nrow
         i = perm(j)
         iao(i+1) = ia(j+1) - ia(j)
 50   continue
c
c get pointers from lengths
c
      iao(1) = 1
      do 51 j=1,nrow
         iao(j+1)=iao(j+1)+iao(j)
 51   continue
c
c copying 
c
      do 100 ii=1,nrow
c
c old row = ii  -- new row = iperm(ii) -- ko = new pointer
c        
         ko = iao(perm(ii)) 
         do 60 k=ia(ii), ia(ii+1)-1 
            jao(ko) = ja(k) 
            if (values) ao(ko) = a(k)
            ko = ko+1
 60      continue
 100  continue
c
      return
c---------end-of-rperm ------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine cperm (nrow,a,ja,ia,ao,jao,iao,perm,job) 
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*), job
      real*8 a(*), ao(*) 
c-----------------------------------------------------------------------
c this subroutine permutes the columns of a matrix a, ja, ia.
c the result id written in the output matrix  ao, jao, iao.
c cperm computes B = A P, where  P is a permutation matrix
c that maps column j into column perm(j), i.e., on return 
c      a(i,j) becomes a(i,perm(j)) in new matrix 
c Y. Saad, May 2, 1990 / modified Jan. 28, 1991. 
c-----------------------------------------------------------------------
c on entry:
c----------
c nrow 	= row dimension of the matrix
c
c a, ja, ia = input matrix in csr format. 
c
c perm	= integer array of length ncol (number of columns of A
c         containing the permutation array  the columns: 
c         a(i,j) in the original matrix becomes a(i,perm(j))
c         in the output matrix.
c
c job	= integer indicating the work to be done:
c 		job = 1	permute a, ja, ia into ao, jao, iao 
c                       (including the copying of real values ao and
c                       the array iao).
c 		job .ne. 1 :  ignore real values ao and ignore iao.
c
c------------
c on return: 
c------------ 
c ao, jao, iao = input matrix in a, ja, ia format (array ao not needed)
c
c Notes:
c------- 
c 1. if job=1 then ao, iao are not used.
c 2. This routine is in place: ja, jao can be the same. 
c 3. If the matrix is initially sorted (by increasing column number) 
c    then ao,jao,iao  may not be on return. 
c 
c----------------------------------------------------------------------c
c local parameters:
      integer k, i, nnz
c
      nnz = ia(nrow+1)-1
      do 100 k=1,nnz
         jao(k) = perm(ja(k)) 
 100  continue
c
c     done with ja array. return if no need to touch values.
c
      if (job .ne. 1) return
c
c else get new pointers -- and copy values too.
c 
      do 1 i=1, nrow+1
         iao(i) = ia(i)
 1    continue
c
      do 2 k=1, nnz
         ao(k) = a(k)
 2    continue
c
      return
c---------end-of-cperm-------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine dperm (nrow,a,ja,ia,ao,jao,iao,perm,qperm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),
     +        qperm(*),job
      real*8 a(*),ao(*) 
c-----------------------------------------------------------------------
c This routine permutes the rows and columns of a matrix stored in CSR
c format. i.e., it computes P A Q, where P, Q are permutation matrices. 
c P maps row i into row perm(i) and Q maps column j into column qperm(j): 
c      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix
c In the particular case where Q is the transpose of P (symmetric 
c permutation of A) then qperm is not needed. 
c note that qperm should be of length ncol (number of columns) but this
c is not checked. 
c-----------------------------------------------------------------------
c Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991. 
c-----------------------------------------------------------------------
c on entry: 
c---------- 
c n 	= dimension of the matrix
c a, ja, 
c    ia = input matrix in a, ja, ia format
c perm 	= integer array of length n containing the permutation arrays
c	  for the rows: perm(i) is the destination of row i in the
c         permuted matrix -- also the destination of column i in case
c         permutation is symmetric (job .le. 2) 
c
c qperm	= same thing for the columns. This should be provided only
c         if job=3 or job=4, i.e., only in the case of a nonsymmetric
c	  permutation of rows and columns. Otherwise qperm is a dummy
c
c job	= integer indicating the work to be done:
c * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P)
c 		job = 1	permute a, ja, ia into ao, jao, iao 
c 		job = 2 permute matrix ignoring real values.
c * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q 
c 		job = 3	permute a, ja, ia into ao, jao, iao 
c 		job = 4 permute matrix ignoring real values.
c		
c on return: 
c-----------
c ao, jao, iao = input matrix in a, ja, ia format
c
c in case job .eq. 2 or job .eq. 4, a and ao are never referred to 
c and can be dummy arguments. 
c Notes:
c------- 
c  1) algorithm is in place 
c  2) column indices may not be sorted on return even  though they may be 
c     on entry.
c----------------------------------------------------------------------c
c local variables 
      integer locjob, mod
c
c     locjob indicates whether or not real values must be copied. 
c     
      locjob = mod(job,2) 
c
c permute rows first 
c 
      call rperm (nrow,a,ja,ia,ao,jao,iao,perm,locjob)
c
c then permute columns
c
      locjob = 0
c
      if (job .le. 2) then
         call cperm (nrow,ao,jao,iao,ao,jao,iao,perm,locjob) 
      else 
         call cperm (nrow,ao,jao,iao,ao,jao,iao,qperm,locjob) 
      endif 
c     
      return
c-------end-of-dperm----------------------------------------------------
c-----------------------------------------------------------------------
      end
