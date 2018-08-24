subroutine bks_pot(ax,ay,az,atmsym,nat,rxyz0, fxyz, ener,virial)
  use parameters_bks !module with parameters read from supp. input file 
  use control_module 
    use cell_module
  implicit double precision (a-h,o-z)
  double precision,dimension(3,nat)::rxyz0,fxyz,vir2
  integer::nat
  character(4),dimension(nat)::atmsym
  double precision::ener,ax,ay,az!,virial
  double precision,parameter::convert=14.39d0
  double precision,parameter::oxchg=-1.2D0,sichg=2.4D0
  integer,dimension(nat)::ilabel
  integer,save:: icounter=10
  double precision,dimension(3,3)::cell_temp,virial

  rcut=cutoff*0.5291772108d0
! print*,rcut,'cut'
  virial=0.
  cell_temp=cell
    do i=1,nat
      if(atmsym(i).eq.'SI')then
         ilabel(i)=14
      else
         ilabel(i)=8
      end if
      end do

       fxyz=0.0d0
       ener=0.0d0

       do i=1,nat-1
       do j=i+1,nat
          isum=ilabel(i)+ilabel(j)
      if(isum.eq.16)then
         a=aoo
         b=boo
         c=coo
         qi=qo
         qj=qo
         eps=.0010511d0
         sig=1.779239d0
!     O-SI
      else if(isum.eq.22)then
         a=asio
         b=bsio
         c=csio
         qi=qo
         qj=qsi
         eps=.00309795d0
         sig=1.313635d0
!     SI-SI
      else if(isum.eq.28)then
         a=asisi
         b=bsisi
         c=csisi
         qi=qsi
         qj=qsi
         eps=0.0d0
         sig=0.0d0
         end if

         rx=rxyz0(1,i)-rxyz0(1,j)
         ry=rxyz0(2,i)-rxyz0(2,j)
         rz=rxyz0(3,i)-rxyz0(3,j)

!         if(pbc)then
!         rx=rx-ax*anint(rx/ax)
!         ry=ry-ay*anint(ry/ay)
!         rz=rz-az*anint(rz/az)
!         end if

         if(pbc)then
            if(cubic)then
             rx=rx-ax*anint(rx/ax)
             ry=ry-ay*anint(ry/ay)
             rz=rz-az*anint(rz/az)
            else
             sxij=fracs(1,i)-fracs(1,j)
             syij=fracs(2,i)-fracs(2,j)
             szij=fracs(3,i)-fracs(3,j)
             sxij=sxij-anint(sxij)
             syij=syij-anint(syij)
             szij=szij-anint(szij)
             rx=cell_temp(1,1)*sxij+cell_temp(1,2)*syij+cell_temp(1,3)*szij
             ry=cell_temp(2,1)*sxij+cell_temp(2,2)*syij+cell_temp(2,3)*szij
             rz=cell_temp(3,1)*sxij+cell_temp(3,2)*syij+cell_temp(3,3)*szij
            end if
          end if




         rijsq=rx*rx + ry*ry + rz*rz
         rij=dsqrt(rijsq)
eps=0.
sig=0.


     if(rij.gt.rcut)goto 123
       U=convert*qi*qj/rij
       U=U+a*dexp(-1.0D0*b*rij)
       sigdivr=sig/rij
       U=U-(c/rij**6) + 4.0d0*eps*( (sigdivr)**30 - (sigdivr)**6 )
        ener=ener+U

       dedr=(-1.0D0*convert*qi*qj/rij**3) + 6.0D0*c/rij**8 - &
     b*a*dexp(-1.0D0*b*rij)/rij - 4.0d0*eps*( (30.0d0*sigdivr**29)*sig/(rij**3) &
     - (6.0d0*sigdivr**5)*sig/(rij**3) )

       f2x=dedr*rx
       f2y=dedr*ry
       f2z=dedr*rz
       fxyz(1,i)=fxyz(1,i)+f2x
       fxyz(2,i)=fxyz(2,i)+f2y
       fxyz(3,i)=fxyz(3,i)+f2z
       fxyz(1,j)=fxyz(1,j)-f2x
       fxyz(2,j)=fxyz(2,j)-f2y
       fxyz(3,j)=fxyz(3,j)-f2z

       
         virial(1,1)=virial(1,1)-f2x*rx
         virial(2,1)=virial(2,1)-f2x*ry
         virial(3,1)=virial(3,1)-f2x*rz
         virial(2,2)=virial(2,2)-f2y*ry
         virial(3,2)=virial(3,2)-f2y*rz
         virial(3,3)=virial(3,3)-f2z*rz


 123   continue

       end do
       end do

       fxyz=-1.0d0*fxyz
! convert to hartree/bohr
       fxyz=fxyz*.529177d0/27.21d0
       ener=ener/27.21138386d0

!      do i=1,nat
!       print*,'ters',fxyz(1,i),fxyz(2,i),fxyz(3,i)
!       end do


       virial=virial/27.21138386d0 !convert to au
!	print*,'bks e and virial',ener,virial
!        stop
!        open(unit=1,file='virial')
 !       write(1,*)virial
 !       close(1)







end subroutine bks_pot

subroutine space(n)
implicit integer (a-z)
do i=1,n
write(761,*)''
end do
end subroutine space
