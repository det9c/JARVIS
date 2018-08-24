subroutine get_hessian
use atom_types_module
use control_module
use constants
use cell_module
use virial_mod
use elastic_module
implicit double precision (a-h,o-z)
!interface
! subroutine get_forces(rxyz0, fxyz, ener, virial) ! input coordinates are BOHR and output force/ener should be Hartree
!  double precision,dimension(:,:),intent(inout)::rxyz0,fxyz
!  double precision,intent(inout)::ener! ,virial
!  double precision,dimension(3,3),intent(inout)::virial
! end subroutine get_forces
!end interface 

double precision,dimension(3*natoms)::amass,scrvec9,vec29
double precision,dimension(3,natoms)::gplus,gminus
double precision::kinetic
double precision,dimension(3,natoms)::rtdt2
double precision,dimension(3,3)::dummy
double precision,dimension(3*natoms,3*natoms)::hesscopy
double precision,dimension(3*natoms,3*natoms)::vecs

print*,'here'
i3n=3*natoms


hessmat=0.0d0

delta=delta_hessian

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)''
write(761,*)'Delta =',delta
write(761,*)''
close(761)



two=2.0d0
one=1.0d0

icol=0
ifirstcol=1
ilastcol=natoms
if(parallel_elastic)then
ifirstcol=ifirst_hess_atom
ilastcol=ilast_hess_atom
icol=3*(ifirst_hess_atom-1)
end if
ifirstrow=icol+1

do i=ifirstcol,ilastcol

open(unit=19,file='current_atom')
write(19,*)'atom =',i,' out of',ilastcol
close(19)


icol=icol+1
amass(icol)=atom_mass(type_atom(i))
xold=xyz(1,i)
xyz(1,i)=xyz(1,i)+delta
rtdt2=xyz


call get_forces(rtdt2, gplus, ener, dummy)
xyz(1,i)=xyz(1,i)-2.0d0*delta
rtdt2=xyz
call get_forces(rtdt2, gminus, ener, dummy)
gplus=(gplus-gminus)/two/delta

xyz(1,i)=xold

icount=0
do k=1,natoms
do j=1,3
icount=icount+1
hessmat(icount,icol)=gplus(j,k)
!print*,icount,icol,gplus(j,k)
end do
end do
!stop   

icol=icol+1
amass(icol)=atom_mass(type_atom(i))
yold=xyz(2,i)
xyz(2,i)=xyz(2,i)+delta
rtdt2=xyz
call get_forces(rtdt2, gplus, ener, dummy)
xyz(2,i)=xyz(2,i)-2.0d0*delta
rtdt2=xyz
call get_forces(rtdt2, gminus, ener, dummy)
gplus=(gplus-gminus)/two/delta
xyz(2,i)=yold
icount=0
do k=1,natoms
do j=1,3
icount=icount+1
hessmat(icount,icol)=gplus(j,k)
end do
end do



icol=icol+1
amass(icol)=atom_mass(type_atom(i))
zold=xyz(3,i)
xyz(3,i)=xyz(3,i)+delta
rtdt2=xyz
call get_forces(rtdt2, gplus, ener, dummy)
xyz(3,i)=xyz(3,i)-2.0d0*delta
rtdt2=xyz
call get_forces(rtdt2, gminus, ener, dummy)
gplus=(gplus-gminus)/two/delta
xyz(3,i)=zold

icount=0
do k=1,natoms
do j=1,3
icount=icount+1
hessmat(icount,icol)=gplus(j,k)
end do
end do


end do

ilastrow=icol

!print*,hessmat

hessmat=-hessmat

! symmetrize hessian
if( .not. parallel_elastic)then
hessmat=hessmat+transpose(hessmat)
hessmat=hessmat/two
end if

tol=1d-12
do i=1,i3n
do j=1,i3n
if(dabs (hessmat(j,i)) .lt.  tol)hessmat(j,i)=0.0d0
end do
end do


open(unit=12,file='HESSIAN')
write(12,*)ifirstrow,ilastrow,hessmat
close(12)

if(parallel_elastic)return




hesscopy=hessmat

amass=1.0d0/dsqrt(amass)
do i=1,i3n
do j=1,i3n
hesscopy(j,i)=hesscopy(j,i)*amass(j)*amass(i)
end do
end do

!print*,hesscopy
!stop
call tred3(i3n,i3n,hesscopy,scrvec9,vec29,vecs)
call tql3(i3n,i3n,scrvec9,vec29,vecs,iout) ! eigenvals are in scrvec







!call eig(hesscopy,vecs,i3n,i3n,0)

factor=5140.36636949d0
do i=1,i3n
scrvec9(i)=dsqrt(abs(scrvec9(i)))*factor*dsign(one,scrvec9(i))
!write(761,*)scrvec9(i)
end do

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'Frequencies cm-1'

write(761,'(6f8.2)')(scrvec9(i),i=1,i3n)
close(761)












end subroutine get_hessian


