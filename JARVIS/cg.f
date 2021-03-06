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
	subroutine conjugate(xyz,natoms,force,input_names,maxcyc,cell
     $,stress_in,refcell,usenrm)



c
c Change the maximum size of the problem dimension here
c
	parameter        (ndim=10000)
	double precision x(ndim),g(ndim),d(ndim),gold(ndim),w(ndim)
	double precision f,eps,tlev,rnorm
        double precision time1,time2,tottime,q(natoms,3)
        logical          finish,refcell,usenrm
	integer          iprint(2),iflag,icall,n,method,mp,lp,i
        integer          iter,nfun
	integer icount

c*******************************
c md stuff 
	integer natoms,maxcyc,irow,itimes
	double precision xyz(3,natoms),force,cell(3,3),stress_in(3,3)
     $,fracs(3,natoms)   
	character*4 input_names(natoms)
	
	


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
	n =         9
	method =    1
	irest =     1   
	iprint(1) = 1 
	iprint(2) = 0  
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
       write(761,*)'            Conjugate Gradient Cell optimization '
       write(761,*)'///////////////////////////////////////////////////'
	close(761)
c
c Get the initial vector x
c
        irow=0
        do i=1,3
        do j=1,3
           irow=irow+1
        x(irow)=cell(i,j)
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
	call fcn(n,x,f,g,natoms,xyz,input_names
     1    ,force,stress_in,itimes,fracs,refcell,usenrm)

c	  if(rnorm.lt.force)then
c         write(*,*)'OPTIMIZATION CONVERGED AFTER',ncycle,' CYCLES'
c        do 500 i=1,n
c        q(  internal(1,i), internal(2,i)  ) = x(i)
c500     continue
c	 stop
c         end if


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
























