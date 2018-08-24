 subroutine exp6_pot(ax,ay,az,atmsym,nat,rxyz0, fxyz, ener, virial)
use control_module
use exp6_module 
use constants 
use atom_types_module
use cell_module
  implicit double precision (a-h,o-z)


  double precision,dimension(3,nat)::rxyz0,fxyz
  integer::nat
  character(4),dimension(nat)::atmsym
  double precision::ener,ax,ay,az
   integer,dimension(nat)::ilabel
   double precision,dimension(3,3)::virial
  integer,save:: icounter=10


rcut=cutoff



       fxyz=0.0d0
       ener=0.0d0
       virial=0.
       do i=1,nat-1
       do j=i+1,nat
         rxij=rxyz0(1,i)-rxyz0(1,j)
         ryij=rxyz0(2,i)-rxyz0(2,j)
         rzij=rxyz0(3,i)-rxyz0(3,j)
         if(pbc)then
            if(cubic)then
             rxij=rxij-ax*anint(rxij/ax)
             ryij=ryij-ay*anint(ryij/ay)
             rzij=rzij-az*anint(rzij/az)
            else
             sxij=fracs(1,i)-fracs(1,j)
             syij=fracs(2,i)-fracs(2,j)
             szij=fracs(3,i)-fracs(3,j)
             sxij=sxij-anint(sxij)
             syij=syij-anint(syij)
             szij=szij-anint(szij)
             rxij=cell(1,1)*sxij+cell(1,2)*syij+cell(1,3)*szij
             ryij=cell(2,1)*sxij+cell(2,2)*syij+cell(2,3)*szij
             rzij=cell(3,1)*sxij+cell(3,2)*syij+cell(3,3)*szij
            end if
          end if

         rijsq=rxij*rxij + ryij*ryij + rzij*rzij
         if(rijsq.gt.cutoff_sq)goto 123
         rij=dsqrt(rijsq)
! load parameters for easier coding
         indi=type_atom(i)
         indj=type_atom(j)
         ii=max(indi,indj)
         jj=min(indi,indj)
         ij=jj+offset_exp6(ii)
         Aij=exp6_paramval(1,ij)
         Bij=exp6_paramval(2,ij)
         Cij=exp6_paramval(3,ij)
         qi=atom_charge(indi)
         qj=atom_charge(indj)

         ener=ener+Aij*dexp(-1.0D0*Bij*rij)-Cij/rijsq**3
         ener=ener+qi*qj/rij
         dedr=-Bij*Aij*dexp(-1.0D0*Bij*rij)/rij + 6.0D0*Cij/rij**8 -qi*qj/rijsq



       f2x=dedr*rxij
       f2y=dedr*ryij
       f2z=dedr*rzij

       fxyz(1,i)=fxyz(1,i)+f2x
       fxyz(2,i)=fxyz(2,i)+f2y
       fxyz(3,i)=fxyz(3,i)+f2z
       fxyz(1,j)=fxyz(1,j)-f2x
       fxyz(2,j)=fxyz(2,j)-f2y
       fxyz(3,j)=fxyz(3,j)-f2z


         virial(1,1)=virial(1,1)-f2x*rxij
         virial(2,1)=virial(2,1)-f2x*ryij
         virial(3,1)=virial(3,1)-f2x*rzij
         virial(2,2)=virial(2,2)-f2y*ryij
         virial(3,2)=virial(3,2)-f2y*rzij
         virial(3,3)=virial(3,3)-f2z*rzij




 123   continue


       end do
       end do
       fxyz=-fxyz
fac=27.21/.52917d0
!do i=1,natoms
!print*,'lj',fxyz(1,i)*fac,fxyz(2,i)*fac,fxyz(3,i)*fac
!end do
print*,'energy',ener*27.21d0
!stop
end subroutine exp6_pot


