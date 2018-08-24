 subroutine morse_pot(ax,ay,az,atmsym,nat,rxyz0, fxyz, ener)
use control_module
use morse_module 
use constants 
use atom_types_module
  implicit double precision (a-h,o-z)


  double precision,dimension(3,nat)::rxyz0,fxyz
  integer::nat
  character(4),dimension(nat)::atmsym
  double precision::ener,ax,ay,az
   integer,dimension(nat)::ilabel
  integer,save:: icounter=10


rcut=cutoff



       fxyz=0.0d0
       ener=0.0d0

       do i=1,nat-1
       do j=i+1,nat
         rxij=rxyz0(1,i)-rxyz0(1,j)
         ryij=rxyz0(2,i)-rxyz0(2,j)
         rzij=rxyz0(3,i)-rxyz0(3,j)
         if(pbc)then
         rxij=rxij-ax*anint(rxij/ax)
         ryij=ryij-ay*anint(ryij/ay)
         rzij=rzij-az*anint(rzij/az)
         end if
         rijsq=rxij*rxij + ryij*ryij + rzij*rzij
         rij=dsqrt(rijsq)
         if(rij.gt.rcut)goto 123
! load parameters for easier coding
         indi=type_atom(i)
         indj=type_atom(j)
         ii=max(indi,indj)
         jj=min(indi,indj)
         ij=jj+offset_morse(ii)
         eo=morse_paramval(1,ij)
         alpha=morse_paramval(2,ij)
         req=morse_paramval(3,ij)

         term1=  exp(-alpha*(rij-req)) 
         ener=ener+eo * ( (1.0d0 -term1)**2 -1.0d0)
         dedr=2.0d0*eo*(1.0d0-term1)*alpha*term1/rij

       f2x=dedr*rxij
       f2y=dedr*ryij
       f2z=dedr*rzij

       fxyz(1,i)=fxyz(1,i)+f2x
       fxyz(2,i)=fxyz(2,i)+f2y
       fxyz(3,i)=fxyz(3,i)+f2z
       fxyz(1,j)=fxyz(1,j)-f2x
       fxyz(2,j)=fxyz(2,j)-f2y
       fxyz(3,j)=fxyz(3,j)-f2z






 123   continue


       end do
       end do
       fxyz=-fxyz




end subroutine morse_pot


