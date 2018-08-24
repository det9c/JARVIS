 subroutine lj_pot(ax,ay,az,atmsym,nat,rxyz0, fxyz, ener, virial)
!  use parameters_bks !module with parameters read from supp. input file
  use control_module
  use lj_module
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
        virial=0.0d0
       do i=1,nat-1
       do j=i+1,nat
         if(i.eq.j)goto 123
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
         rij=dsqrt(rijsq)
!         if(rij.gt.rcut)goto 123
         if(rij.lt.4.724d0 .or. rij.gt.18.897d0)goto 123
! load parameters for easier coding
         indi=type_atom(i)
         indj=type_atom(j)
         ii=max(indi,indj)
         jj=min(indi,indj)
         ij=jj+offset_lj(ii)
         eps=lj_paramval(1,ij)
         sig=lj_paramval(2,ij)
!         expon1=lj_paramval(3,ij)
!         expon2=lj_paramval(4,ij)
         expon1=12.
         expon2=6.
         sigdivr=sig/rij
         ener=ener+4.0d0 * eps *(sigdivr**expon1  - sigdivr**expon2)
         dedr=-4.0d0*eps*sig*(expon1*sigdivr**(expon1-1.0d0)-expon2*sigdivr**(expon2-1.0d0))/(rijsq*rij)

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
!print*,'lj energy',ener
!do i=1,natoms
!print*,'lj',fxyz(1,i),fxyz(2,i),fxyz(3,i)
!end do

end subroutine lj_pot

