c
c       --------------------------------------------------------------------
c       Main program for running the conjugate gradient methods described in 
c       the paper:
c
c       Gilbert, J.C. and Nocedal, J. (1992). "Global Convergence Properties 
c       of Conjugate Gradient Methods", SIAM Journal on Optimization, Vol. 2,
c       pp. 21-42.
c
c       A web-based Server which solves unconstrained nonlinear optimization
c       problems using this Conjugate Gradient code can be found at:
c
c       http://www-neos.mcs.anl.gov/neos/solvers/UCO:CGPLUS/
c
c       Written by G. Liu, J. Nocedal and R. Waltz
c       October 1998
c       --------------------------------------------------------------------
c 
c	subroutine conjugate(natoms,nopt,force,q,internal,icount2)
        subroutine cg_xyz(xyz,natoms,force,rmsmax,input_names,
     $                 maxcyc,f,volu,sout)



c
c Change the maximum size of the problem dimension here
c
	parameter        (ndim=10000)
	double precision x(3*natoms),g(3*natoms),d(3*natoms),gold(3*natoms),
     $w(3*natoms),x2(3*natoms)       
	double precision f,eps,tlev
        double precision time1,time2,tottime,q(natoms,3)
        logical          finish,refcell,usenrm
	integer          iprint(2),iflag,icall,n,method,mp,lp,i
        integer          iter,nfun,ii,jj,kk
	integer icount

c*******************************
c md stuff 
	integer natoms,maxcyc,irow,itimes,j,mm,ibestpt
	double precision xyz(3,natoms),force,cell(3,3),stress_in(3,3)
     $,fracs(3,natoms),rmsmax,volu,sout(3,3),dummy(3,3),rfac,cd,
     $rnorm,gmaxc,rmsnrm,rlow,beststr(3,3),bestxyz(3,natoms),besten       
	character*4 input_names(natoms),w1,w2
	
	


c end md stuff
c*********************************


	common /cgdd/    mp,lp
	common /runinf/  iter,nfun
        data one/1.0D+0/

        FINISH= .FALSE.
c
c Read problem input information
c
	itimes=0
	n =         3*natoms
	method =    1
	irest =     1   
	iprint(1) = 1 
	iprint(2) = 0  
         rfac=0.5291772108d0
         rlow=1.0d10
c
c Check for correct dimension value n
c
	if (n .lt. 0) then
	   iflag = -3
	   write(*,850)
	   go to 50
	end if
	if (n .gt. ndim) then
	   iflag = -3
	   write(*,860)
	   go to 50
	end if
	
	open(unit=761,file='OUTPUT.DTPOLY',access='append')
       write(761,*)'///////////////////////////////////////////////////'
       write(761,*)'     Conjugate Gradient Geometry optimization      '
      write(761,31)'        Convergence tolerance on max force:  ',force
      write(761,31)'        Convergence tolerance on RMS norm  ',rmsmax
      write(761,*)'               Maximum Number of Steps:',maxcyc
      write(761,*)'    ///////////////////////////////////////'
 31   format(A,f6.5)

      call space(3)
      write(761,*)' Cycle          Energy(a.u.)            ||F||   RMS'
      write(761,*)'---------------------------------------------'
       write(761,*)'///////////////////////////////////////////////////'
	close(761)
c
c Get the initial vector x
c
        irow=0
       do  i=1,natoms
       do  j=1,3
       irow=irow+1
       x(irow)=xyz(j,i)
	end do
	end do


c
c Print parameters
c
	if (iprint(1) .ge. 0) then
	   write (*,820)
	   write (*,840) n, method, irest
	end if

        ICALL=0
c
c This is the convergence constant 
c
        EPS= 1.0D-5

c IFLAG=0 indicates an initial entry to program

        IFLAG=0
c
c Begin counting CPU time. 
c (Note: This function ma y not work on all operating systems.)
c
c        call timer(time1)

  20    CONTINUE
c
c Calculate the function and gradient values here
c
c Rosenbrock test function
	itimes=itimes+1
         if(itimes.gt.maxcyc)then
            open(unit=761,file='OUTPUT.DTPOLY',access='append')
         write(761,*)'Optimization has exceeded ',maxcyc,' steps'
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

         cd=29421.0107637093d0/volu
          write(761,*)'Stress tensor for point ',ibestpt
          write(761,*)sout(1,1)*cd,sout(1,2)*cd,sout(1,3)*cd
          write(761,*)sout(2,1)*cd,sout(2,2)*cd,sout(2,3)*cd
          write(761,*)sout(3,1)*cd,sout(3,2)*cd,sout(3,3)*cd
         close(761)
         return
          end if


         do i=1,3*natoms
            x2(i)=x(i)
         end do
         call get_forces(x2, g, f, dummy)

         do ii=1,3
         do  jj=1,3
            dummy(jj,ii)=dummy(jj,ii)*29421.0107637093d0/volu
            end do
            end do



          if(itimes.eq.1)then
             open(unit=761,file='OUTPUT.DTPOLY',access='append')
                       write(761,*)'Initial Stress Tensor (GPa)'
          write(761,23)dummy(1,1),dummy(1,2),dummy(1,3)
          write(761,23)dummy(2,1),dummy(2,2),dummy(2,3)
          write(761,23)dummy(3,1),dummy(3,2),dummy(3,3)
 23                                   format(f25.15,f25.15,f25.15)
                        close(761)
          end if




          do ii=1,3
          do mm=1,3

             if(dabs(dummy(mm,ii)).gt.100.0)then

          open(unit=761,file='OUTPUT.DTPOLY',access='append')
          write(761,*)'Stress tensor is awful. Taking new step.'
          write(761,*)dummy(1,1),dummy(1,2),dummy(1,3)
          write(761,*)dummy(2,1),dummy(2,2),dummy(2,3)
          write(761,*)dummy(3,1),dummy(3,2),dummy(3,3)
          f=dabs(dummy(1,1)*10.0)

         do kk=1,3
         do jj=1,3
         sout(jj,kk)=dummy(jj,kk)*volu/29421.0107637093d0
         end do
         end do
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
         return
         end if


          end do
          end do



         rnorm=0.0D0
         gmaxc=1.0d-12
         do 80 i=1,3*natoms
            g(i)=-g(i)
            if( dabs(g(i)) .gt. gmaxc)then
            gmaxc=dabs(g(i))
            end if
         rnorm=rnorm+g(i)**2
 80            continue
         rmsnrm=dsqrt(rnorm/float(3*natoms))
         rnorm=dsqrt(rnorm)
c         write(761,*)''
c         write(761,*)'****************************************************'
c         write(761,*)'*       AT THE END OF CYCLE NUMBER ',NCYCLE,':      '
c         WRITE(761,*)'* GRADIENT NORM OF OPTIMIZED COORDINATES = ',RNORM
c         WRITE(761,*)'*            ENERGY = ',F,' Hartrees                 '
c         write(761,*)'*****************************************************'
         w1='NOPE'
         w2='NOPE'
         if(rmsnrm.lt.rmsmax)w1='YEP'
         if(gmaxc .lt. force)w2='YEP'

          open(unit=761,file='OUTPUT.DTPOLY',access='append')
         write(761,12)' ',itimes,'   ',f,'    ',rnorm,rmsnrm,w1,gmaxc,w2
         close(761)
 12            format(A,I3,A,f20.10,A,f12.8,f12.8,A5,f12.8,A5)


         if(gmaxc .lt. rlow)then
            rlow=gmaxc
            ibestpt = itimes
         do kk=1,natoms
         jj=3*(kk-1)+1
         bestxyz(1,kk)=x(jj)
         bestxyz(2,kk)=x(jj+1)
         bestxyz(3,kk)=x(jj+2)
         end do

         do kk=1,3
         do jj=1,3
         beststr(jj,kk)=dummy(jj,kk)*volu/29421.0107637093d0
         end do
         end do

         besten = f
         end if


         if(rmsnrm.lt.rmsmax  .and. gmaxc .lt. force)then
            call space(1)
            open(unit=761,file='OUTPUT.DTPOLY',access='append')
       write(761,*)'Success! A stationary point has been found! :) '
       write(761,*)'RMS norm =',rmsnrm
       write(761,*)'||F||  =',rnorm
       write(761,*)'Max gradient component =',gmaxc
       call space(1)
          write(761,*)' Stress tensor at convergence  (GPa)'
          write(761,23)dummy(1,1),dummy(1,2),dummy(1,3)
          write(761,23)dummy(2,1),dummy(2,2),dummy(2,3)
          write(761,23)dummy(3,1),dummy(3,2),dummy(3,3)
          call space(2)

      write(761,*)'Final coordinates in ./bfgs_history.xyz '
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

         do kk=1,3
         do jj=1,3
         sout(jj,kk)=dummy(jj,kk)*volu/29421.0107637093d0
         end do
         end do







         return
         end if









 59            format(A,f30.18,f30.18,f30.18)
 57                  format(A,f20.15,f20.15,f20.15,i2,i2,i2,i2,i2,i2)







 30	CONTINUE
c
c Call the main optimization code
c
        CALL CGFAM(N,X,F,G,D,GOLD,IPRINT,EPS,W,
     *            IFLAG,IREST,METHOD,FINISH )
c
c       IFLAG=
c              0 : successful termination
c              1 : return to evaluate F and G
c              2 : return with a new iterate, try termination test
c             -i : error
c
	IF(IFLAG.LE.0.OR.ICALL.GT.10000) GO TO 50
        IF(IFLAG.EQ.1) THEN
	   ICALL=ICALL + 1	   
	   GO TO 20
        ENDIF	 
        IF(IFLAG.EQ.2) THEN
c
c Termination Test.  The user may replace it by some other test. However, 
c the parameter 'FINISH' must be set to 'TRUE' when the test is satisfied.
c
	   TLEV= EPS*(ONE + DABS(F))
	   I=0
  40	   I=I+1
	   IF(I.GT.N) THEN
              FINISH = .FALSE.
	      GO TO 30
	   ENDIF
	   IF(DABS(G(I)).GT.TLEV) THEN
	      GO TO 30
	   ELSE
	      GO TO 40
	   ENDIF

        ENDIF

  50	continue
c
c End CPU counting
c
        call timer(time2)
c
c Calculate the elapsed time
c
        tottime = time2-time1
c
c Code has terminated; print final results
c
	if (iprint(1).ge.0.and.iflag.ge.0) then
	   write (*,890) f
	   write (*,900) tottime
	end if
	return
	
c
c Formatting
c
 800	format (12x, i3)
 820	format (//, ' Conjugate Gradient Minimization Routine', /)
 840	format (/, ' n	  =', i6, /,
     *             ' method   =', i6,/,
     *             ' irest    =', i6,/)
 850	format (/,'  Error: negative N value'/)
 860	format (/,'  Error: N too large, increase parameter ndim'/)
 890	format (/,' f(x*) =', 1pd16.8)
 900	format (' It took ',1pd12.6,' CPU seconds'/)

	end
























