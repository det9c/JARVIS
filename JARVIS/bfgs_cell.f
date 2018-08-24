      subroutine bfgs_cell(xyz,natoms,input_names,maxcyc,cell
     $,stress_in,refcell,usenrm,optalgo,fixed_frame,cfrc,crms
     $,scale_factor) 
      integer          nmax, mmax,iwrite
      parameter        (nmax=1024, mmax=17)
c        nmax is the dimension of the largest problem to be solved.
c        mmax is the maximum number of limited memory corrections.
c     Declare the variables needed by the code.
c       A description of all these variables is given at the end of
c       the driver.
      character*60     task, csave
      logical          lsave(4),restart,refcell,fixed_frame
      integer          n, m, iprint,MPTS,lk,
     +                 nbd(nmax), iwa(3*nmax), isave(44)
      double precision f, factr, pgtol,a,b,c,d,e,
     +                 x(nmax), l(nmax), u(nmax), g(nmax), dsave(29),
     +                 wa(2*mmax*nmax+4*nmax+12*mmax*mmax+12*mmax),denom
c     Declare a few additional variables for this sample problem.
      double precision t1, t2
      integer          i,ntime,j,k,ibestpt

C***************************************************
C***************************************************
C***************************************************
      integer nparams,ndata,iopt,natoms,irow
      parameter(nparams=9,ndata=1)
      DOUBLE PRECISION cell(3,3),factor,xyz(3,natoms)
     $,fracs(3,natoms),sinv(3,3),cell_copy(3,3),stress(3,3)
     $,cinv(3,3),scr(3),smat(3,3),cell_deriv(3,3),volume
     $,cell_tran(3,3),rnorm,stress_in(3,3),force,press,
     $     stress_copy(3,3),eold,enew,estrain,estrold
     $,delstr,flowest,xss,yss,zss,ssx,ssy,ssz,xyz_magic(3,natoms)     
     $,cell_rot(3,3),rxyz_rot(3,natoms),stress_save(3,3),
     $     cell_un(3,3),cell_copy2(3,3),cfrc,crms,gmaxc,rnormc,
     $scale_factor(3,3)
      integer optx(3),opty(3),optz(3),indx(3)
      integer itimes,iunit,i1,i2,i3,maxcyc,optalgo
      character(5) name(natoms)
      character*4 input_names(natoms),w1,w3
      logical usenrm
C***************************************************
C***************************************************
C***************************************************

c     We wish to have output at every iteration.
      iwrite=2
      iprint = 1

c     We specify the tolerances in the stopping criteria.

      delta=1.0D-3
      factr=0.0D0
      pgtol=0.0D0
      restart=.false.
      itimes=0
      factor=.5291772108d0
      iunit=6
      enew=0.0
      estrain=0.0
      delstr=.001d0
      flowest=1d10
c     We specify the dimension n of the sample problem and the number
c        m of limited memory corrections stored.  (n and m should not
c        exceed the limits nmax and mmax respectively.)
C

        n=nparams
        m=3
C
C     DEFINE STARTING POINT AND THE VARIABLE BOUNDS

c        open(unit=1,file='cell.input')
c        read(1,*)cell(1,1),cell(1,2),cell(1,3)
        optx(1)=1
        opty(1)=1
        optz(1)=1
c        read(1,*)cell(2,1),cell(2,2),cell(2,3)
        optx(2)=1
        opty(2)=1
        optz(2)=1
c        read(1,*)cell(3,1),cell(3,2),cell(3,3)
        optx(3)=1
        opty(3)=1
        optz(3)=1
c        close(1)

	

         do i=1,3
         do j=1,3
c         cell(j,i)=cell(j,i)*factor
         cell_copy(j,i)=cell(j,i)
         end do
         end do


c         call matprt(stress_in,3,3,3,3)

         call migs(cell_copy,3,sinv,indx)

          do i=1,natoms
c             print*,input_names(i),xyz(1,i),xyz(2,i),xyz(3,i)
          ssx=sinv(1,1)*xyz(1,i)+sinv(2,1)*xyz(2,i)+sinv(3,1)*xyz(3,i)
          ssy=sinv(1,2)*xyz(1,i)+sinv(2,2)*xyz(2,i)+sinv(3,2)*xyz(3,i)
          ssz=sinv(1,3)*xyz(1,i)+sinv(2,3)*xyz(2,i)+sinv(3,3)*xyz(3,i)
          xss=ssx
          yss=ssy
          zss=ssz
          fracs(1,i)=xss
          fracs(2,i)=yss
          fracs(3,i)=zss
            end do

          do i=1,natoms
          xyz_magic(1,i)=xyz(1,i)
          xyz_magic(2,i)=xyz(2,i)
          xyz_magic(3,i)=xyz(3,i)
          end do



c        read(1,*)stress_in(1,1),stress_in(1,2),stress_in(1,3)
c        read(1,*)stress_in(2,1),stress_in(2,2),stress_in(2,3)
c        read(1,*)stress_in(3,1),stress_in(3,2),stress_in(3,3)

c        close(1)

        irow=0
        do i=1,3
        do j=1,3
           irow=irow+1
        x(irow)=cell(i,j)
        end do
        end do

      do 10 i=1,n
         nbd(i)=0.
 10        continue

          
           open(unit=761,file='OUTPUT.DTPOLY',access='append')
      write(761,*)'-----------------------------------------------'
      write(761,*)'|                 CELL OPT                    |'
      write(761,*)'|        WRITTEN BY DeCARLOS TAYLOR (2010)    |'
      write(761,*)'-----------------------------------------------'

       write(761,*)'Initial Cell'


         do i=1,3
         do j=1,3
         cell_copy(j,i)=cell(i,j)
         end do
         end do

       
      call cell_volume(cell_copy,volume)
      

       write(761,*)''
       write(761,*)'OPT A:',optx(1),opty(1),optz(1)
       write(761,*)'OPT B:',optx(2),opty(2),optz(2)
       write(761,*)'OPT C:',optx(3),opty(3),optz(3)
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







ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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

           irow=0
           do i=1,3
           do j=1,3
           irow=irow+1
           x(irow)=cell(i,j)
           end do
           end do

         do i=1,natoms
          xyz_magic(1,i)=rxyz_rot(1,i)
          xyz_magic(2,i)=rxyz_rot(2,i)
          xyz_magic(3,i)=rxyz_rot(3,i)
          end do




           end if











          if(usenrm)then
          write(761,*)'Following gradient norm'
          else
          write(761,*)'Following energy'
          end if
      write(761,30)'  Convergence tolerance on max cell force:  ',cfrc
      write(761,30)'   Convergence tolerance on cell RMS norm  ',crms
      write(761,*)'               Maximum Number of Steps:',maxcyc
      write(761,*)'    ///////////////////////////////////////'
 30      format(A,f6.5)
         write(761,*)''
         write(761,*)''
         write(761,*)''


      if(optalgo.eq.0)write(761,*)'Using fracs. Fixed Algorithm'
      if(optalgo.eq.1)write(761,*)'Using fracs. Update Algorithm'
      if(optalgo.eq.2)write(761,*)'Using XYZ Algorithm'
      write(761,*)''
      write(761,*)''
      write(761,*)''







          write(761,*)'******Starting optimization*********'
          close(761)

c
      task = 'START'

c        ------- the beginning of the loop ----------

 111    continue

c     This is the call to the L-BFGS-B code.

      call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,
     +            csave,lsave,isave,dsave)

      open(unit=761,file='OUTPUT.DTPOLY',access='append')
      write(761,*)'Task is',task
      close(761)

      if (task(1:2) .eq. 'FG') then
c        the minimization routine has returned to request the
c        function f and gradient g values at the current x.
         itimes=itimes+1
         open(unit=761,file='OUTPUT.DTPOLY',access='append')
         write(761,*)''
         write(761,*)''
         write(761,*)''
         write(761,*)''
         write(761,*)''
         write(761,*)' -------------------------------'
         write(761,*)'|  Iteration Number',itimes,'   |'
         write(761,*)' -------------------------------'
         close(761)
         f=0.0D0

         if(itimes.gt.maxcyc)then
            open(unit=761,file='OUTPUT.DTPOLY',access='append')
         write(761,*)'Optimization has exceeded ',maxcyc,' steps'
         close(761)
         return
         end if




         open(unit=133,file='HISTORY.CELL')
         write(133,*)'---------------------------------'
          write(133,*)'optimization step ',itimes
          write(133,*)'cell is:'
         write(133,33)'A',x(1)*factor,x(2)*factor,x(3)*factor
         write(133,33)'B',x(4)*factor,x(5)*factor,x(6)*factor
         write(133,33)'C',x(7)*factor,x(8)*factor,x(9)*factor
         close(133)

        irow=0
        do i=1,3
        do j=1,3
           irow=irow+1
        cell(i,j)=x(irow)
        end do
        end do

         do i=1,3
         do j=1,3
         cell_copy(j,i)=cell(i,j)
         end do
         end do


         open(unit=761,file='OUTPUT.DTPOLY',access='append')
        call cell_volume(cell_copy,volume)
        close(761)
c        open(unit=761,file='OUTPUT.DTPOLY',access='append')

 33           format(A2,f20.15,f20.15,f20.15)




          do i=1,natoms
          xss=fracs(1,i)
          yss=fracs(2,i)
          zss=fracs(3,i)
          xyz(1,i)=x(1)*xss+x(4)*yss+x(7)*zss
          xyz(2,i)=x(2)*xss+x(5)*yss+x(8)*zss
          xyz(3,i)=x(3)*xss+x(6)*yss+x(9)*zss
          end do


          if(optalgo.eq.2)then
          open(unit=761,file='OUTPUT.DTPOLY',access='append')
           write(761,*)''
           write(761,*)''
          write(761,*)'Using Coordinates From Previous Step'
          write(761,*)''
          write(761,*)''
          close(761)   
          do i=1,natoms
          xyz(1,i)=xyz_magic(1,i)
          xyz(2,i)=xyz_magic(2,i)
          xyz(3,i)=xyz_magic(3,i)
          end do
          end if





          eold=enew
          call get_stress(cell_copy,xyz, stress, enew,volume)

          
          do i=1,natoms
          xyz_magic(1,i)=xyz(1,i)
          xyz_magic(2,i)=xyz(2,i)
          xyz_magic(3,i)=xyz(3,i)
          end do





         open(unit=133,file='HISTORY.CELL',access='append')
          do i=1,natoms
          write(133,33)input_names(i),xyz(1,i)*factor,xyz(2,i)*factor,
     $         xyz(3,i)*factor
          end do
          close(133)





          open(unit=761,file='OUTPUT.DTPOLY',access='append')
          write(761,*)''
          write(761,*)'energy is',enew
c          write(761,*)'delta E is',enew-eold
          write(761,*)''
          write(761,*)'Current Stress Tensor (GPa)'
          write(761,23)stress(1,1),stress(1,2),stress(1,3)
          write(761,23)stress(2,1),stress(2,2),stress(2,3)
          write(761,23)stress(3,1),stress(3,2),stress(3,3)
 23              format(f25.15,f25.15,f25.15)

          press=(stress(1,1)+stress(2,2)+stress(3,3))/3.0d0       
          write(761,*)''
         write(761,*)'Current internal pressure',press

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
          stress_copy(j,i)=stress(j,i)
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
          rnorm=0.
         irow=0
        do i=1,3
        do j=1,3
           irow=irow+1
        if(fixed_frame)then
        g(irow)=cell_deriv(i,j)
        else
        g(irow)=cell_un(i,j)
        end if
        rnorm=rnorm+g(irow)**2
        end do
        end do

        denom=9.0d0
        if(fixed_frame)denom=6.0d0
        
        rnormc=dsqrt(rnorm/denom)
        write(761,*)''

            gmaxc=0.0d0
            do  i=1,9
            if( dabs(g(i)) .gt. gmaxc)then
            gmaxc=dabs(g(i))
            end if
            end do






        if(rnormc.lt.flowest)then
           flowest=rnormc
           ibestpt=itimes
           open(unit=34,file='best_point')
           write(34,*)flowest,flowest
           write(34,*)cell_copy
           write(34,*)xyz
           close(34)
        end if


         w1='NOPE'
         w3='NOPE'

         if(rnormc.lt.crms)w1='YEP'
         if(gmaxc.lt. cfrc)w3='YEP'


          write(761,13)'          Max Cell Gradient',gmaxc,w3
          write(761,13)'                 Cell ||F||',rnormc,w1
           write(761,13)''
          write(761,*)'Best so far',ibestpt,'Max Grad =',flowest

 13       format(A,f20.10,A5)
















c        f=enew+estrain
        f=estrain
        f=enew
        if(usenrm)f=rnormc
        write(761,*)'f returned =',f
        close(761)


        if(rnormc.lt.crms .and. gmaxc.lt.cfrc)then



                      open(unit=34,file='best_point')
           write(34,*)flowest,flowest
           write(34,*)cell_copy
           write(34,*)xyz
           close(34)




           open(unit=761,file='OUTPUT.DTPOLY',access='append')
           write(761,*)''
           write(761,*)''
           write(761,*)''
           write(761,*)''
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
           write(761,*)'    CELL OPTIMIZATION HAS CONVERGED  '
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
        call dump_final(cell_copy,volume,stress_copy)

           write(761,*)''
           write(761,*)''
           write(761,*)''
           write(761,*)''
        close(761)
        return


        end if






c compute new fractionals for new cell above

        if(optalgo.eq.1)then
           open(unit=761,file='OUTPUT.DTPOLY',access='append')
           write(761,*)'Updating Fractionals'
           close(761)
          do i=1,natoms
          ssx=cinv(1,1)*xyz(1,i)+cinv(2,1)*xyz(2,i)+cinv(3,1)*xyz(3,i)
          ssy=cinv(1,2)*xyz(1,i)+cinv(2,2)*xyz(2,i)+cinv(3,2)*xyz(3,i)
          ssz=cinv(1,3)*xyz(1,i)+cinv(2,3)*xyz(2,i)+cinv(3,3)*xyz(3,i)
          xss=ssx
          yss=ssy
          zss=ssz
          fracs(1,i)=xss
          fracs(2,i)=yss
          fracs(3,i)=zss
          end do
        end if   




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
c      stop
            open(unit=761,file='OUTPUT.DTPOLY',access='append')
            write(761,*)'quitting',task
            close(761)


            return
       end

              SUBROUTINE mmult(r,t,z,NBASIS)
c
c   s and b are matrices to be multiplied
c   n is dimension of matrices
c
c
c
c
c
      DOUBLE PRECISION r(NBASIS,NBASIS),t(NBASIS,NBASIS),
     $z(NBASIS,NBASIS),total
      do 61 i=1,NBASIS
         do 62 j=1,NBASIS
            total=0.0D0
            do 63 k=1,NBASIS
               total=total+r(i,k)*t(k,j)
 63                                 continue
               z(i,j)=total
 62                                 continue
 61                                                            continue
         return
         end

