subroutine elascon_hess
use atom_types_module
use virial_mod
use cell_module
use elastic_module
use constants
use control_module
implicit double precision (a-h,o-z)

interface
!subroutine get_stress(cell_bfgs,xyz_bfgs, stress_bfgs, energy_bfgs,vol)
!double precision,dimension(:,:),intent(inout)::xyz_bfgs
!double precision,dimension(3,3),intent(inout)::cell_bfgs,stress_bfgs
!double precision,intent(inout)::energy_bfgs
!double precision,intent(in)::vol
!end subroutine get_stress

   subroutine cell_volume(cell_in,volume)
   double precision,dimension(3,3),intent(in)::cell_in
   double precision,intent(inout)::volume
   end subroutine cell_volume




end interface



double precision,dimension(3,3)::strain,voigt,cell2,cell_zero,stress_plus,stress_minus,scratch3,cell_temp
double precision,dimension(:,:),allocatable::fracs_zero
character(4)::ch
integer,dimension(6):: indx
integer,dimension(3*natoms)::indx2
double precision,dimension(6,6)::cij,tensorav,vecs,copy,cinverse
double precision,dimension(3*natoms,6)::dei
double precision,dimension(6,3*natoms)::deixhess
double precision,dimension(3*natoms)::grad_plus,grad_minus
double precision,dimension(3*natoms,3*natoms)::hinverse
elas_is_on=.true.
cij=0.

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'          -----------------------------------------------'
write(761,*)'          |            ELASTIC CONSTANTS                |'
write(761,*)'          |          HESSIAN ALGORITHM                  |'
write(761,*)'          |        WRITTEN BY DeCARLOS TAYLOR (2011)    |'
write(761,*)'          -----------------------------------------------'


if(parallel_elastic)then
write(761,*)'Parallel hessian evaluation'
write(761,*)'Evaluating derivatives for atoms',ifirst_hess_atom,' to',ilast_hess_atom
else
write(761,*)'Evaluating full hessian'
end if


write(761,*)'Stress corrections to Cij tensor:'
!call matprt(stress_in,3,3,3,3)
write(761,17)stress_in(1,1),stress_in(1,2),stress_in(1,3)
write(761,17)stress_in(2,1),stress_in(2,2),stress_in(2,3)
write(761,17)stress_in(3,1),stress_in(3,2),stress_in(3,3)
call space(2)
close(761)
17 format(f20.6,f20.6,f20.6)


stress_in=-stress_in
maxcyc=1
!zero=0.
!eps_elas=.001
cell_zero=cell_elastic
allocate(fracs_zero(3,natoms))
fracs_zero=fracs

call cell_volume(cell_zero,vol_zero)

if(.not. skipelas)then
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'Zero strain'
close(761)
end if





do i=1,natoms

          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)

          xyz_elastic(1,i)=cell_zero(1,1)*xss+cell_zero(1,2)*yss+cell_zero(1,3)*zss
          xyz_elastic(2,i)=cell_zero(2,1)*xss+cell_zero(2,2)*yss+cell_zero(2,3)*zss
          xyz_elastic(3,i)=cell_zero(3,1)*xss+cell_zero(3,2)*yss+cell_zero(3,3)*zss
end do
cell2=cell_zero


if(skipelas)goto 641

call get_stress(cell2,xyz_elastic, stress_plus,f,vol_zero)



! do e1
voigt=zero
voigt(1,1)=eps_elas
strain=unit_matrix+voigt
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'e1 Voigt strain'
call matprt(strain,3,3,3,3)
close(761)
!cell_temp=cell_zero*au_to_ang
!call matprt(cell_temp,3,3,3,3)
cell2=matmul(strain,cell_zero)
!cell2=cell2/au_to_ang
!call matprt(cell2*au_to_ang,3,3,3,3)

do i=1,natoms

          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          
          xyz_elastic(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz_elastic(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz_elastic(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do

cell_stats=.false.
just_volume=.true.
call cell_volume(cell2,vol2)
just_volume=.false.
cell_stats=.true.

call get_stress(cell2,xyz_elastic, stress_plus,f,vol2)

stress_plus=stress_elastic
icount=0
do i=1,natoms
do j=1,3
icount=icount+1
grad_plus(icount)=fxyz_global(j,i)
end do
end do
!call matprt(auforce_to_bar*1.0d-4*stress_plus/vol_zero,3,3,3,3)


voigt=zero
voigt(1,1)=-1.0*eps_elas
strain=unit_matrix+voigt
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'-e1 Voigt strain'
!cell_temp=cell_zero*au_to_ang
cell2=matmul(strain,cell_zero)
!cell2=cell2/au_to_ang
call matprt(strain,3,3,3,3)
close(761)
!call matprt(cell2*au_to_ang,3,3,3,3)
do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz_elastic(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz_elastic(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz_elastic(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss

end do
cell_stats=.false.
just_volume=.true.
call cell_volume(cell2,vol2)
just_volume=.false.
cell_stats=.true.

call get_stress(cell2,xyz_elastic, stress_minus,f,vol2)


stress_minus=stress_elastic
scratch3=auforce_to_bar*1.0d-4*(stress_plus-stress_minus)/(2.0*eps_elas*vol_zero)
!print*,scratch3(1,1),auforce_to_bar*1.0d-4*stress_plus(1,1)/vol_zero,auforce_to_bar*1.0d-4*stress_minus(1,1)/vol_zero
!stop
scratch3=-scratch3
icount=0
do i=1,natoms
do j=1,3
icount=icount+1
grad_minus(icount)=fxyz_global(j,i)
end do
end do
do i=1,3*natoms
dei(i,1)=grad_plus(i)-grad_minus(i)
end do

!call matprt(scratch3,3,3,3,3)
cij(1,1)=scratch3(1,1)
cij(1,2)=scratch3(2,2)
cij(1,3)=scratch3(3,3)
cij(1,4)= scratch3(2,3)
cij(1,5)=(scratch3(1,3)) 
cij(1,6)=(scratch3(1,2) )

! do e2
voigt=zero
voigt(2,2)=eps_elas
strain=unit_matrix+voigt
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'e2 Voigt strain'
call matprt(strain,3,3,3,3)
close(761)
cell2=matmul(strain,cell_zero)
!call matprt(cell2,3,3,3,3)

do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz_elastic(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz_elastic(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz_elastic(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do

cell_stats=.false.
just_volume=.true.
call cell_volume(cell2,vol2)
just_volume=.false.
cell_stats=.true.

call get_stress(cell2,xyz_elastic, stress_plus,f,vol2)

stress_plus=stress_elastic
icount=0
do i=1,natoms
do j=1,3
icount=icount+1
grad_plus(icount)=fxyz_global(j,i)
end do
end do


voigt=zero
voigt(2,2)=-1.0*eps_elas
strain=unit_matrix+voigt
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'-e2 Voigt strain'
cell2=matmul(strain,cell_zero)
call matprt(strain,3,3,3,3)
close(761)
do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz_elastic(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz_elastic(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz_elastic(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do

cell_stats=.false.
just_volume=.true.
call cell_volume(cell2,vol2)
just_volume=.false.
cell_stats=.true.

call get_stress(cell2,xyz_elastic,stress_minus,f,vol2)

stress_minus=stress_elastic
scratch3=auforce_to_bar*1.0d-4*(stress_plus-stress_minus)/(2.0*eps_elas*vol_zero)
scratch3=-scratch3
cij(2,1)=scratch3(1,1)
cij(2,2)=scratch3(2,2)
cij(2,3)=scratch3(3,3)
cij(2,4)=scratch3(2,3)
cij(2,5)=scratch3(1,3)
cij(2,6)=scratch3(1,2)
icount=0
do i=1,natoms
do j=1,3
icount=icount+1
grad_minus(icount)=fxyz_global(j,i)
end do
end do
do i=1,3*natoms
dei(i,2)=grad_plus(i)-grad_minus(i)
end do






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!e3
! do e3
voigt=zero
voigt(3,3)=eps_elas
strain=unit_matrix+voigt
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'e3 Voigt strain'

call matprt(strain,3,3,3,3)
close(761)
cell2=matmul(strain,cell_zero)
!call matprt(cell2,3,3,3,3)

do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz_elastic(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz_elastic(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz_elastic(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do

cell_stats=.false.
just_volume=.true.
call cell_volume(cell2,vol2)
just_volume=.false.
cell_stats=.true.

call get_stress(cell2,xyz_elastic, stress_plus,f,vol2)


stress_plus=stress_elastic
icount=0
do i=1,natoms
do j=1,3
icount=icount+1
grad_plus(icount)=fxyz_global(j,i)
end do
end do

voigt=zero
voigt(3,3)=-1.0*eps_elas
strain=unit_matrix+voigt
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'-e3 Voigt strain'

cell2=matmul(strain,cell_zero)
call matprt(strain,3,3,3,3)
close(761)
do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz_elastic(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz_elastic(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz_elastic(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do

cell_stats=.false.
just_volume=.true.
call cell_volume(cell2,vol2)
just_volume=.false.
cell_stats=.true.

call get_stress(cell2,xyz_elastic,stress_minus,f,vol2)


stress_minus=stress_elastic
scratch3=auforce_to_bar*1.0d-4*(stress_plus-stress_minus)/(2.0*eps_elas*vol_zero)
scratch3=-scratch3
cij(3,1)=scratch3(1,1)
cij(3,2)=scratch3(2,2)
cij(3,3)=scratch3(3,3)
cij(3,4)=scratch3(2,3)
cij(3,5)=scratch3(1,3)
cij(3,6)=scratch3(1,2)
icount=0
do i=1,natoms
do j=1,3
icount=icount+1
grad_minus(icount)=fxyz_global(j,i)
end do
end do
do i=1,3*natoms
dei(i,3)=grad_plus(i)-grad_minus(i)
end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!e4
! do e4
voigt=zero
voigt(2,3)=eps_elas
voigt(3,2)=eps_elas
strain=unit_matrix+voigt
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'e4 Voigt strain'

call matprt(strain,3,3,3,3)
close(761)
cell2=matmul(strain,cell_zero)
!call matprt(cell2,3,3,3,3)

do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz_elastic(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz_elastic(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz_elastic(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do

cell_stats=.false.
just_volume=.true.
call cell_volume(cell2,vol2)
just_volume=.false.
cell_stats=.true.

call get_stress(cell2,xyz_elastic, stress_plus,f,vol2)

stress_plus=stress_elastic
icount=0
do i=1,natoms
do j=1,3
icount=icount+1
grad_plus(icount)=fxyz_global(j,i)
end do
end do

voigt=zero
voigt(2,3)=-eps_elas
voigt(3,2)=-eps_elas
strain=unit_matrix+voigt
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'-e4 Voigt strain'

cell2=matmul(strain,cell_zero)
call matprt(strain,3,3,3,3)
close(761)
do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz_elastic(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz_elastic(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz_elastic(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do
cell_stats=.false.
just_volume=.true.
call cell_volume(cell2,vol2)
just_volume=.false.
cell_stats=.true.

call get_stress(cell2,xyz_elastic, stress_minus,f,vol2)


stress_minus=stress_elastic
scratch3=auforce_to_bar*1.0d-4*(stress_plus-stress_minus)/(4.0*eps_elas*vol_zero)
scratch3=-scratch3
cij(4,1)=scratch3(1,1)
cij(4,2)=scratch3(2,2)
cij(4,3)=scratch3(3,3)
cij(4,4)=scratch3(2,3)
cij(4,5)=scratch3(1,3)
cij(4,6)=scratch3(1,2)
icount=0
do i=1,natoms
do j=1,3
icount=icount+1
grad_minus(icount)=fxyz_global(j,i)
end do
end do
do i=1,3*natoms
dei(i,4)=(grad_plus(i)-grad_minus(i))/2.0d0
end do





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!e5
! do e5
voigt=zero
voigt(1,3)=eps_elas
voigt(3,1)=eps_elas
strain=unit_matrix+voigt
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'e5 Voigt strain'

call matprt(strain,3,3,3,3)
close(761)
cell2=matmul(strain,cell_zero)
!call matprt(cell2,3,3,3,3)

do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz_elastic(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz_elastic(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz_elastic(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do

cell_stats=.false.
just_volume=.true.
call cell_volume(cell2,vol2)
just_volume=.false.
cell_stats=.true.

call get_stress(cell2,xyz_elastic, stress_plus,f,vol2)

stress_plus=stress_elastic
icount=0
do i=1,natoms
do j=1,3
icount=icount+1
grad_plus(icount)=fxyz_global(j,i)
end do
end do

voigt=zero
voigt(1,3)=-eps_elas
voigt(3,1)=-eps_elas
strain=unit_matrix+voigt
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'-e5 Voigt strain'

cell2=matmul(strain,cell_zero)
call matprt(strain,3,3,3,3)
close(761)
do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz_elastic(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz_elastic(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz_elastic(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do

cell_stats=.false.
just_volume=.true.
call cell_volume(cell2,vol2)
just_volume=.false.
cell_stats=.true.

call get_stress(cell2,xyz_elastic, stress_minus,f,vol2)


stress_minus=stress_elastic
scratch3=auforce_to_bar*1.0d-4*(stress_plus-stress_minus)/(4.0*eps_elas*vol_zero)
scratch3=-scratch3
cij(5,1)=scratch3(1,1)
cij(5,2)=scratch3(2,2)
cij(5,3)=scratch3(3,3)
cij(5,4)=scratch3(2,3)
cij(5,5)=scratch3(1,3)
cij(5,6)=scratch3(1,2)
icount=0
do i=1,natoms
do j=1,3
icount=icount+1
grad_minus(icount)=fxyz_global(j,i)
end do
end do
do i=1,3*natoms
dei(i,5)=(grad_plus(i)-grad_minus(i))/2.0d0
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!e6
! do e6
voigt=zero
voigt(1,2)=eps_elas
voigt(2,1)=eps_elas
strain=unit_matrix+voigt
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'e6 Voigt strain'

call matprt(strain,3,3,3,3)
close(761)
cell2=matmul(strain,cell_zero)
!call matprt(cell2,3,3,3,3)

do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz_elastic(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz_elastic(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz_elastic(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do

cell_stats=.false.
just_volume=.true.
call cell_volume(cell2,vol2)
just_volume=.false.
cell_stats=.true.

call get_stress(cell2,xyz_elastic, stress_plus,f,vol2)

stress_plus=stress_elastic
icount=0
do i=1,natoms
do j=1,3
icount=icount+1
grad_plus(icount)=fxyz_global(j,i)
end do
end do

voigt=zero
voigt(1,2)=-eps_elas
voigt(2,1)=-eps_elas
strain=unit_matrix+voigt
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'-e6 Voigt strain'

cell2=matmul(strain,cell_zero)
call matprt(strain,3,3,3,3)
close(761)
do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz_elastic(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz_elastic(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz_elastic(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do
cell_stats=.false.
just_volume=.true.
call cell_volume(cell2,vol2)
just_volume=.false.
cell_stats=.true.

call get_stress(cell2,xyz_elastic, stress_minus,f,vol2)


stress_minus=stress_elastic
scratch3=auforce_to_bar*1.0d-4*(stress_plus-stress_minus)/(4.0*eps_elas*vol_zero)
scratch3=-scratch3
cij(6,1)=scratch3(1,1)
cij(6,2)=scratch3(2,2)
cij(6,3)=scratch3(3,3)
cij(6,4)=scratch3(2,3)
cij(6,5)=scratch3(1,3)
cij(6,6)=scratch3(1,2)
icount=0
do i=1,natoms
do j=1,3
icount=icount+1
grad_minus(icount)=fxyz_global(j,i)
end do
end do
do i=1,3*natoms
dei(i,6)=(grad_plus(i)-grad_minus(i))/2.0d0
end do

dei=-dei/(2.0d0*eps_elas)
open(unit=99,file='dei')
write(99,*)dei
close(99)
open(unit=99,file='cij-unrelax')
write(99,*)cij
close(99)



open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)''
write(761,*)''
write(761,*)''
write(761,*)'               Unrelaxed Elastic Constant Tensor (Voigt Ordering/GPa): '
write(761,*)'            -----------------------------------------------------------'
call matprt(cij,6,6,6,6)

641 continue
write(761,*)'Building hessian'
do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)

          xyz_elastic(1,i)=cell_zero(1,1)*xss+cell_zero(1,2)*yss+cell_zero(1,3)*zss
          xyz_elastic(2,i)=cell_zero(2,1)*xss+cell_zero(2,2)*yss+cell_zero(2,3)*zss
          xyz_elastic(3,i)=cell_zero(3,1)*xss+cell_zero(3,2)*yss+cell_zero(3,3)*zss
end do
cell=cell_zero
xyz=xyz_elastic
close(761)
print*,'calling'
call get_hessian


open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'Hessian complete'
   if(parallel_elastic .or. skipelas)then
      close(761)
   stop
   end if
!!!cij=-cij
write(761,*)''
write(761,*)''
write(761,*)''
write(761,*)'               Unrelaxed Elastic Constant Tensor (Voigt Ordering/GPa): '
write(761,*)'            -----------------------------------------------------------'
call matprt(cij,6,6,6,6)


hessold=hessmat

call migs(hessmat,3*natoms,hinverse,indx2)

!hessmat=matmul(hessold,hinverse)




deixhess=matmul(transpose(dei),hinverse)
copy=matmul(deixhess,dei)
copy=-copy*auforce_to_bar*1.0d-4/vol_zero

write(761,*)''
write(761,*)''
write(761,*)'                   Relaxation Contributions (Voigt Ordering/GPa):      '
write(761,*)'            -----------------------------------------------------------'
call matprt(copy,6,6,6,6)






cij=cij+copy
write(761,*)''
write(761,*)''
write(761,*)'                 Relaxed Elastic Constant Tensor (Voigt Ordering/GPa):       '
write(761,*)'            -----------------------------------------------------------'
call matprt(cij,6,6,6,6)







tensorav=cij+transpose(cij)
tensorav=tensorav/2.0d0

write(761,*)''
write(761,*)''
write(761,*)'              Symmetrized Elastic constant tensor (Voigt Ordering/GPa):'
write(761,*)'              ---------------------------------------------------------'
call matprt(tensorav,6,6,6,6)

B_v=tensorav(1,1)+tensorav(2,2)+tensorav(3,3)+2.0d0*(tensorav(1,2)+tensorav(1,3)+tensorav(2,3))
B_v=B_v/9.0d0
shear_v=tensorav(1,1)+tensorav(2,2)+tensorav(3,3)+3.0d0*(tensorav(4,4)+tensorav(5,5)+tensorav(6,6))-tensorav(1,2)-tensorav(1,3)-tensorav(2,3)
shear_v=shear_v/15.0d0

copy=tensorav
call eig(tensorav,vecs,6,6,6,6)

write(761,*)''
write(761,*)''
write(761,*)'                                      Eigenvalues '
write(761,*)'                                      -----------'
call matprt(tensorav,6,6,6,6)


write(761,*)''
write(761,*)''
write(761,*)'                                      Eigenmodes     '
write(761,*)'                                      ----------'
call matprt(vecs,6,6,6,6)


call migs(copy,6,cinverse,indx)
write(761,*)''
write(761,*)''
write(761,*)'                                Compliance Tensor (1/GPa)       '
write(761,*)'                                -------------------------'
call matprt(cinverse,6,6,6,6)

B_r=cinverse(1,1)+cinverse(2,2)+cinverse(3,3)+2.0d0*(cinverse(1,2)+cinverse(1,3)+cinverse(2,3))
B_r=1.0d0/B_r
shear_r=4.0d0*(cinverse(1,1)+cinverse(2,2)+cinverse(3,3))+3.0d0*(cinverse(4,4)+cinverse(5,5)+cinverse(6,6))-4.0d0*(cinverse(1,2)+cinverse(1,3)+cinverse(2,3))
shear_r=15.0d0/shear_r

write(761,*)''
write(761,*)''
write(761,*)'                               Mechanical Properties (GPa)'
write(761,*)'                               ---------------------------'
write(761,*)'   Convention:                Voigt                  Reuss                Hill'
write(761,*)'   -----------                -----                  -----                ----'
write(761,40)'   Bulk Modulus             ',B_v,'            ',B_r,'            ',(B_v+B_r)/2.0d0
write(761,40)'   Shear Modulus            ',shear_v,'            ',shear_r,'            ',(shear_v+shear_r)/2.0d0
write(761,*)'                               ----------------------'
write(761,50)'           Youngs Moduli(x/y/z): ',1.0d0/cinverse(1,1),' / ',1.0d0/cinverse(2,2),' / ',1.0d0/cinverse(3,3)
40 format(A,f10.4,A,f10.4,A,f10.4)
50 format(A,f10.4,A,f10.4,A,f10.4)



close(761)
!print*,'modes'
!call matprt(vecs,6,6,6,6)







stop



return
end subroutine elascon_hess
