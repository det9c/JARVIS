c*********************************************************************
      function urand()
c=====================================================================
c     Return the next pseudo-random deviate from a sequence which is
c     uniformly distributed in the interval [0,1]
c
c     Uses the function ran0, the "minimal standard" random number
c     generator of Park and Miller (Comm. ACM 31, 1192-1201, Oct 1988;
c     Comm. ACM 36 No. 7, 105-110, July 1993).
c=====================================================================
      implicit none
c
c     Input - none
c
c     Output
      double precision     urand
c
c     Local
      integer  iseed
      double precision     ran0
      external ran0
c
c     Common block to make iseed visible to rninit (and to save
c     it between calls)
      common /rnseed/ iseed
c
      urand = ran0( iseed )
      return
      end
c*********************************************************************
      subroutine rninit( seed )
c=====================================================================
c     Initialize random number generator urand with given seed
c=====================================================================
      implicit none
c
c     Input
      integer seed
c
c     Output - none
c
c     Local
      integer iseed
c
c     Common block to communicate with urand
      common /rnseed/ iseed
c
c     Set the seed value
      iseed = seed
      if(iseed.le.0) iseed=123456
      return
      end
c*********************************************************************
      function ran0( seed )
c=====================================================================
c     "Minimal standard" pseudo-random number generator of Park and
c     Miller.  Returns a uniform random deviate r s.t. 0 < r < 1.0.
c     Set seed to any non-zero integer value to initialize a sequence,
c     then do not change seed between calls for successive deviates
c     in the sequence.
c
c     References:
c        Park, S. and Miller, K., "Random Number Generators: Good Ones
c           are Hard to Find", Comm. ACM 31, 1192-1201 (Oct. 1988)
c        Park, S. and Miller, K., in "Remarks on Choosing and Imple-
c           menting Random Number Generators", Comm. ACM 36 No. 7,
c           105-110 (July 1993)
c=====================================================================
c *** Declaration section ***
c
      implicit none
c
c     Input/Output:
      integer seed
c
c     Output:
      double precision ran0
c
c     Constants:
      integer A,M,Q,R
      parameter (A=48271,M=2147483647,Q=44488,R=3399)
      double precision SCALE,EPS,RNMX
      parameter (SCALE=1./M,EPS=1.2e-7,RNMX=1.-EPS)
c
c     Local:
      integer j
c
c *** Executable section ***
c
      j = seed/Q
      seed = A*(seed-j*Q)-R*j
      if (seed .lt. 0) seed = seed+M
      ran0 = min(seed*SCALE,RNMX)
      return
      end

