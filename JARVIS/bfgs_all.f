      subroutine bfgs_all(xyz,natoms,input_names,stress_in,usenrm,
     $                maxcyc,cell,afrc,cfrc,arms,crms,fixed_frame,
     $                 scale_factor)

      integer          nmax, mmax,iwrite,natoms
      parameter        (nmax=15000, mmax=17)
      character*60     task, csave
      logical          lsave(4),fixed_frame,usenrm
      integer          n, m, iprint,ii,jj,kk,mm,
     +                 nbd(nmax), iwa(3*nmax), isave(44)
      double precision f, factr, pgtol, force,xyz(3,natoms)
     +               ,x(3*natoms+9), l(nmax), u(nmax), g(3*natoms+9),
     $     dsave(29), denom,stress_in(3,3),
     +    wa(2*mmax*nmax+4*nmax+12*mmax*mmax+12*mmax),rnorma,rnormc
      integer          i,k,NCYCLE,ispot,maxcyc,ibestpt
      double precision g2(3*natoms+9),dummy(3,3),x2(3*natoms+9),rmsnrm
     $,volu,gmaxc,rmsmax,beststr(3,3),bestxyz(3,natoms),besten,rlow
     $,sout(3,3),cd,cell(3,3),cell_copy(3,3),stress(3,3),gmaxa
     $,enew,volume,forces(3,natoms),cinv(3,3),smat(3,3),cell_un(3,3)
     $,flowa,flowc,afrc,cfrc,arms,crms,cell_rot(3,3),rxyz_rot(3,natoms)
     $,stress_save(3,3),cell_copy2(3,3),stress_copy(3,3),cell_deriv(3,3)
     $,scale_factor(3,3) 
      character*4 input_names(natoms)
      character*4 w1,w2,w3,w4
      integer jopt,irow,icount,indx(3),j,icstart
      iwrite=2
      iprint = 1
      pgtol=0.0D0
      factr=0.
      n=3*natoms+9
      icstart=3*natoms+1
      m=20
      rlow=1.0d10
      flowa=1d10
      flowc=1d10

c
c     input coordinates in bohr
c
c      do i=1,natoms
c      write(*,*)input_names(i),xyz(1,i),xyz(2,i),xyz(3,i)
c      end do
 


      icount=0
      do 100 i=1,natoms
      do 101 j=1,3
      icount=icount+1   
      x(icount)=xyz(j,i)
 101  continue
 100  continue



        do i=1,3
        do j=1,3
        icount=icount+1
        x(icount)=cell(i,j)
        end do
        end do







c     set bounds on variables (none for optimization)
      do 10 i=1,n
         nbd(i)=0
  10  continue

c      do  i=3*natoms+1,n
c         nbd(i)=1
c         l(i)=0.5d0
c         end do



c 
      task = 'START'
      NCYCLE=0
      open(unit=761,file='OUTPUT.DTPOLY',access='append')
      write(761,*)''
      write(761,*)'    /////////////////////////////////////////////'
      write(761,*)'              Geometry+Cell optimization      '
      write(761,30)'  Convergence tolerance on max atomic force:  ',afrc
      write(761,30)'   Convergence tolerance on atomic RMS norm  ',arms
      write(761,30)'  Convergence tolerance on max cell force:  ',cfrc
      write(761,30)'   Convergence tolerance on cell RMS norm  ',crms
      write(761,*)'               Maximum Number of Steps:',maxcyc
      write(761,*)'    ///////////////////////////////////////'
 30   format(A,f6.5)
      
      call space(3)
          write(761,*)''
          write(761,*)'Target Stress Tensor (GPa)'
          write(761,23)stress_in(1,1),stress_in(1,2),stress_in(1,3)
          write(761,23)stress_in(2,1),stress_in(2,2),stress_in(2,3)
          write(761,23)stress_in(3,1),stress_in(3,2),stress_in(3,3)

         write(761,*)''
          write(761,*)'Tensor Scaling'
      write(761,23)scale_factor(1,1),scale_factor(1,2),scale_factor(1,3)
      write(761,23)scale_factor(2,1),scale_factor(2,2),scale_factor(2,3)
      write(761,23)scale_factor(3,1),scale_factor(3,2),scale_factor(3,3)


          call space(3)

      open(unit=15,file='bfgs_history.xyz')
c      write(761,*)' Cycle          Energy(a.u.)            ||F||   RMS'
c      write(761,*)'---------------------------------------------'
c        ------- the beginning of the loop ----------


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         do i=1,3
         do j=1,3
         cell_copy(j,i)=cell(i,j)
         end do
         end do

          if(fixed_frame)then
          write(761,*)''
          write(761,*)''
          write(761,*)''
          write(761,*)''



          write(761,*)'Rotating to internal orientation'
           call rotate_cell(cell_copy,xyz,natoms,cell_rot,rxyz_rot)
           call cell_volume(cell_rot,volume)


            do i=1,3
           do j=1,3
           cell(i,j)=cell_rot(j,i)
           end do
           end do

                 icount=0
      do i=1,natoms
      do j=1,3
      icount=icount+1
      x(icount)=rxyz_rot(j,i)
      end do
      end do


        do i=1,3
        do j=1,3
        icount=icount+1
        x(icount)=cell(i,j)
        end do
        end do
        end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





















 111  continue

c     This is the call to the L-BFGS-B code.



      call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,
     +            csave,lsave,isave,dsave)

      if (task(1:2) .eq. 'FG') then
c        the minimization routine has returned to request the
c        function f and gradient g values at the current x.
c
c
c
         NCYCLE=NCYCLE+1
         close(761)
         if(ncycle.gt.maxcyc)then
            stop
         end if


         icount=0
         do i=1,natoms
         do j=1,3
         icount=icount+1
         xyz(j,i)=x(icount)
         end do   
         end do
        
        do i=1,3
        do j=1,3
        icount=icount+1
        cell(i,j)=x(icount)
        end do
        end do

         do i=1,3
         do j=1,3
         cell_copy(j,i)=cell(i,j)
         end do
         end do

         open(unit=761,file='OUTPUT.DTPOLY',access='append')
              write(761,*)''
              write(761,*)''
              write(761,*)''
              write(761,*)'------------------------------------------'
               write(761,*)'   Optimization Step =',ncycle
              write(761,*)''
              write(761,*)''
         call cell_volume(cell_copy,volume)
         close(761)

         call get_stress_forces(cell_copy,xyz,stress,enew,volume,forces)
         f=enew


             open(unit=761,file='OUTPUT.DTPOLY',access='append')
              write(761,*)''
              write(761,*)''
              write(761,*)''
           write(761,*)'                     Stress Tensor (GPa)'
          write(761,23)stress(1,1),stress(1,2),stress(1,3)
          write(761,23)stress(2,1),stress(2,2),stress(2,3)
          write(761,23)stress(3,1),stress(3,2),stress(3,3)
 23                            format(f25.15,f25.15,f25.15)

         do i=1,3
         do j=1,3
         stress(j,i)=stress(j,i)*scale_factor(j,i)
         end do
         end do

          write(761,*)'Stress after scaling...'
          write(761,23)stress(1,1),stress(1,2),stress(1,3)
          write(761,23)stress(2,1),stress(2,2),stress(2,3)
          write(761,23)stress(3,1),stress(3,2),stress(3,3)






c copy stress and cell for use again
          do i=1,3
          do j=1,3
          stress_save(j,i)=stress(j,i)
          end do
          end do
         do i=1,3
         do j=1,3
         cell_copy2(j,i)=cell(j,i)
         end do
         end do

          if(fixed_frame)then
          stress(1,2)=0.0d0
          stress(1,3)=0.0d0
          stress(2,3)=0.0d0

          do i=1,3
          do j=1,3
          stress(j,i)=-stress(j,i)+stress_in(j,i)
          stress_copy(j,i)=stress(j,i)
          end do
          end do
c convert cp2k stress tensor to cell derivatives
        do i=1,3
        do j=1,3
           irow=irow+1
        stress(i,j)=1.0d0*stress(i,j)/29421.0107637093d0
        end do
        end do
        call migs(cell,3,cinv,indx)
        call mmult(stress,cinv,smat,3)

        do i=1,3
        do j=1,3
        cell_deriv(j,i)=smat(j,i)*volume
        end do
        end do

        write(761,*)''
          write(761,*)'Current Cell Derivatives - Projected Rotations'
          write(761,23)cell_deriv(1,1),cell_deriv(1,2),cell_deriv(1,3)
          write(761,23)cell_deriv(2,1),cell_deriv(2,2),cell_deriv(2,3)
          write(761,23)cell_deriv(3,1),cell_deriv(3,2),cell_deriv(3,3)
          end if



c now get derivatives without projections

          do i=1,3
          do j=1,3
          stress(j,i)=stress_save(j,i)
          end do
          end do
         do i=1,3
         do j=1,3
         cell(j,i)=cell_copy2(j,i)
         end do
         end do

          do i=1,3
          do j=1,3
          stress(j,i)=-stress(j,i)+stress_in(j,i)
          end do
          end do
c convert cp2k stress tensor to cell derivatives
        do i=1,3
        do j=1,3
        stress(i,j)=1.0d0*stress(i,j)/29421.0107637093d0
        end do
        end do
        call migs(cell,3,cinv,indx)
        call mmult(stress,cinv,smat,3)

        do i=1,3
        do j=1,3
        cell_un(j,i)=smat(j,i)*volume
        end do
        end do

        write(761,*)''
          write(761,*)'Current Cell Derivatives - Unprojected'
          write(761,23)cell_un(1,1),cell_un(1,2),cell_un(1,3)
          write(761,23)cell_un(2,1),cell_un(2,2),cell_un(2,3)
          write(761,23)cell_un(3,1),cell_un(3,2),cell_un(3,3)


c copy gradient to g vector

      icount=0
      do i=1,natoms
      do j=1,3
      icount=icount+1
      g(icount)=-forces(j,i)
      end do
      end do


        do i=1,3
        do j=1,3
         icount=icount+1
        if(fixed_frame)then
        g(icount)=cell_deriv(i,j)
        else
        g(icount)=cell_un(i,j)
        end if
        end do
        end do



         rnorma=0.0D0
         rnormc=0.0d0
         gmaxc=1.0d-12
         gmaxa=1.0d-12
         do 80 i=1,3*natoms
            if( dabs(g(i)) .gt. gmaxa)then
            gmaxa=dabs(g(i))
            end if
            rnorma=rnorma+g(i)**2
 80               continue

         do  i=icstart,n
            if( dabs(g(i)) .gt. gmaxc)then
            gmaxc=dabs(g(i))
            end if
            rnormc=rnormc+g(i)**2
            end do

            if(usenrm)f=dsqrt(rnorma+rnormc)

               denom=9.0d0
          if(fixed_frame)denom=6.0d0



            rnorma=dsqrt(rnorma/float(3*natoms))
            rnormc=dsqrt(rnormc/denom)




c      if(gmaxa .lt. flowa .and. gmaxc .lt. flowc)then
        if( rnormc .lt. flowc)then
           flowa=rnorma
           flowc=rnormc
           ibestpt=ncycle
           open(unit=34,file='best_point')
           write(34,*)flowc,flowa
           write(34,*)cell_copy
           write(34,*)xyz
           close(34)
        end if




c         write(761,*)''
c         write(761,*)'****************************************************'
c         write(761,*)'*       AT THE END OF CYCLE NUMBER ',NCYCLE,':      '
c         WRITE(761,*)'* GRADIENT NORM OF OPTIMIZED COORDINATES = ',RNORM
c         WRITE(761,*)'*            ENERGY = ',F,' Hartrees                 '
c         write(761,*)'*****************************************************'
         w1='NOPE'
         w2='NOPE'
         w3='NOPE'
         w4='NOPE'
         if(rnormc.lt.crms)w1='YEP'
         if(rnorma .lt. arms)w2='YEP'
         if(gmaxc.lt. cfrc)w3='YEP'
         if(gmaxa .lt. afrc)w4='YEP'

c          open(unit=761,file='OUTPUT.DTPOLY',access='append')
              write(761,*)''
              write(761,*)''
          write(761,*)'            Energy=',f
          write(761,12)'         Max atomic gradient',gmaxa,w4
          write(761,12)'           Max cell gradient',gmaxc,w3
          write(761,12)'                Atomic ||F||',rnorma,w2
          write(761,12)'                  Cell ||F||',rnormc,w1
           write(761,*)''
          write(761,13)'Best Point',ibestpt,'Cell=',flowc,'Atomic',flowa
         close(761)
c 12      format(A,I3,A,f20.10,A,f12.8,f12.8,A5,f12.8,A5)
 12            format(A,f20.10,A5)
 13            format(A,i5,A6,f12.8,A7,f12.8)




         rfac=0.5291772108d0
c write info to history file
         write(15,*)natoms
         write(15,*)'optimization step number',ncycle
         do kk=1,natoms
         jj=3*(kk-1)+1
         write(15,50)input_names(kk),x(jj)*rfac,
     $        x(jj+1)*rfac,x(jj+2)*rfac
         end do
 50    format(A,f20.15,f20.15,f20.15)

c          if(rnorm.lt.force)then
         if(rnorma .lt. arms  .and. rnormc .lt. crms .and. gmaxa 
     $  .lt. afrc .and. gmaxc .lt. cfrc)then
            call space(1)
            open(unit=761,file='OUTPUT.DTPOLY',access='append')
           write(761,*)''
           write(761,*)''
           write(761,*)''
       write(761,*)'Success! A stationary point has been found! :) '
           write(761,*)''
      write(761,*)'Final coordinates in ./bfgs_history.xyz '
           write(761,*)''
           write(761,*)''
           write(761,*)''

      close(761)

       open(unit=123,file='bfgs.xyz')
       rewind 123
         do kk=1,natoms
         jj=3*(kk-1)+1
         write(123,59)input_names(kk),x(jj)*rfac,
     $        x(jj+1)*rfac,x(jj+2)*rfac

         xyz(1,kk)=x(jj)
         xyz(2,kk)=x(jj+1)
         xyz(3,kk)=x(jj+2)

         end do
         close(123)

           open(unit=34,file='best_point')
           write(34,*)rnormc,rnorma
           write(34,*)cell_copy
           write(34,*)xyz
           close(34)






         return
         end if

 59      format(A,f30.18,f30.18,f30.18)
 57      format(A,f20.15,f20.15,f20.15,i2,i2,i2,i2,i2,i2)


c          go back to the minimization routine.
         goto 111
      endif
c
      if (task(1:5) .eq. 'NEW_X') then
         goto 111
       end if

c        the minimization routine has returned with a new iterate,
c         and we have opted to continue the iteration.

c           ---------- the end of the loop -------------

c     If task is neither FG nor NEW_X we terminate execution.
c
c
c    if we get here then optimization failed
        open(unit=761,file='OUTPUT.DTPOLY',access='append')
         write(761,*)'*********OPTIMIZATION FAILED********'
         write(761,*)'Returning data for point ',ibestpt
         do kk=1,natoms
         do jj=1,3
         xyz(jj,kk)=bestxyz(jj,kk)
         end do
         end do

         do kk=1,3
         do jj=1,3
         sout(jj,kk)=beststr(jj,kk)
         end do
         end do
         f = besten
         close(761)

      return

      end




 
         
      




c======================= The end of driver1 ============================

c     --------------------------------------------------------------
c             DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     n is an INTEGER variable that must be set by the user to the
c       number of variables.  It is not altered by the routine.
c
c     m is an INTEGER variable that must be set by the user to the
c       number of corrections used in the limited memory matrix.
c       It is not altered by the routine.  Values of m < 3  are
c       not recommended, and large values of m can result in excessive
c       computing time. The range  3 <= m <= 20 is recommended. 
c
c     x is a DOUBLE PRECISION array of length n.  On initial entry
c       it must be set by the user to the values of the initial
c       estimate of the solution vector.  Upon successful exit, it
c       contains the values of the variables at the best point
c       found (usually an approximate solution).
c
c     l is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the lower bounds on the variables. If
c       the i-th variable has no lower bound, l(i) need not be defined.
c
c     u is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the upper bounds on the variables. If
c       the i-th variable has no upper bound, u(i) need not be defined.
c
c     nbd is an INTEGER array of dimension n that must be set by the
c       user to the type of bounds imposed on the variables:
c       nbd(i)=0 if x(i) is unbounded,
c              1 if x(i) has only a lower bound,
c              2 if x(i) has both lower and upper bounds, 
c              3 if x(i) has only an upper bound.
c
c     f is a DOUBLE PRECISION variable.  If the routine setulb returns
c       with task(1:2)= 'FG', then f must be set by the user to
c       contain the value of the function at the point x.
c
c     g is a DOUBLE PRECISION array of length n.  If the routine setulb
c       returns with taskb(1:2)= 'FG', then g must be set by the user to
c       contain the components of the gradient at the point x.
c
c     factr is a DOUBLE PRECISION variable that must be set by the user.
c       It is a tolerance in the termination test for the algorithm.
c       The iteration will stop when
c
c        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c       where epsmch is the machine precision which is automatically
c       generated by the code. Typical values for factr on a computer
c       with 15 digits of accuracy in double precision are:
c       factr=1.d+12 for low accuracy;
c             1.d+7  for moderate accuracy; 
c             1.d+1  for extremely high accuracy.
c       The user can suppress this termination test by setting factr=0.
c
c     pgtol is a double precision variable.
c       On entry pgtol >= 0 is specified by the user.  The iteration
c         will stop when
c
c                 max{|proj g_i | i = 1, ..., n} <= pgtol
c
c         where pg_i is the ith component of the projected gradient.
c       The user can suppress this termination test by setting pgtol=0.
c
c     wa is a DOUBLE PRECISION  array of length 
c       (2mmax + 4)nmax + 12mmax^2 + 12mmax used as workspace.
c       This array must not be altered by the user.
c
c     iwa is an INTEGER  array of length 3nmax used as
c       workspace. This array must not be altered by the user.
c
c     task is a CHARACTER string of length 60.
c       On first entry, it must be set to 'START'.
c       On a return with task(1:2)='FG', the user must evaluate the
c         function f and gradient g at the returned value of x.
c       On a return with task(1:5)='NEW_X', an iteration of the
c         algorithm has concluded, and f and g contain f(x) and g(x)
c         respectively.  The user can decide whether to continue or stop
c         the iteration. 
c       When
c         task(1:4)='CONV', the termination test in L-BFGS-B has been 
c           satisfied;
c         task(1:4)='ABNO', the routine has terminated abnormally
c           without being able to satisfy the termination conditions,
c           x contains the best approximation found,
c           f and g contain f(x) and g(x) respectively;
c         task(1:5)='ERROR', the routine has detected an error in the
c           input parameters;
c       On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task
c         contains additional information that the user can print.
c       This array should not be altered unless the user wants to
c          stop the run for some reason.  See driver2 or driver3
c          for a detailed explanation on how to stop the run 
c          by assigning task(1:4)='STOP' in the driver.
c
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     csave  is a CHARACTER working array of length 60.
c
c     lsave is a LOGICAL working array of dimension 4.
c       On exit with task = 'NEW_X', the following information is
c         available:
c       lsave(1) = .true.  the initial x did not satisfy the bounds;
c       lsave(2) = .true.  the problem contains bounds;
c       lsave(3) = .true.  each variable has upper and lower bounds.
c
c     isave is an INTEGER working array of dimension 44.
c       On exit with task = 'NEW_X', it contains information that
c       the user may want to access:
c         isave(30) = the current iteration number;
c         isave(34) = the total number of function and gradient
c                         evaluations;
c         isave(36) = the number of function value or gradient
c                                  evaluations in the current iteration;
c         isave(38) = the number of free variables in the current
c                         iteration;
c         isave(39) = the number of active constraints at the current
c                         iteration;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in isave
c
c     dsave is a DOUBLE PRECISION working array of dimension 29.
c       On exit with task = 'NEW_X', it contains information that
c         the user may want to access:
c         dsave(2) = the value of f at the previous iteration;
c         dsave(5) = the machine precision epsmch generated by the code;
c         dsave(13) = the infinity norm of the projected gradient;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in dsave
c
c     --------------------------------------------------------------
c           END OF THE DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     << An example of subroutine 'timer' for AIX Version 3.2 >>
c
c     subroutine timer(ttime)
c     double precision ttime
c     integer itemp, integer mclock
c     itemp = mclock()
c     ttime = dble(itemp)*1.0d-2
c     return
c     end
c-----------------------------------------------------------------------
