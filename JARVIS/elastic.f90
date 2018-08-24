subroutine elascon
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
integer:: indx(6)
double precision,dimension(6,6)::cij,tensorav,vecs,copy,cinverse

elas_is_on=.true.
cij=0.

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'          -----------------------------------------------'
write(761,*)'          |            ELASTIC CONSTANTS                |'
write(761,*)'          |        WRITTEN BY DeCARLOS TAYLOR (2011)    |'
write(761,*)'          -----------------------------------------------'


write(761,*)'Stress corrections to Cij tensor:'
!call matprt(stress_in,3,3,3,3)
write(761,17)stress_in(1,1),stress_in(1,2),stress_in(1,3)
write(761,17)stress_in(2,1),stress_in(2,2),stress_in(2,3)
write(761,17)stress_in(3,1),stress_in(3,2),stress_in(3,3)
call space(2)
close(761)
17 format(f20.6,f20.6,f20.6)


stress_in=-stress_in

!zero=0.
!eps_elas=.001
cell_zero=cell_elastic
allocate(fracs_zero(3,natoms))
fracs_zero=fracs

call cell_volume(cell_zero,vol_zero)
open(unit=761,file='OUTPUT.DTPOLY',access='append')
if(parallel_elastic)then
write(761,*)'Parallel evaluation of elastic constant tensor'
if(e1_only)write(761,*)'Current job will evaluate 1st row of tensor'
if(e2_only)write(761,*)'Current job will evaluate 2nd row of tensor'
if(e3_only)write(761,*)'Current job will evaluate 3rd row of tensor'
if(e4_only)write(761,*)'Current job will evaluate 4th row of tensor'
if(e5_only)write(761,*)'Current job will evaluate 5th row of tensor'
if(e6_only)write(761,*)'Current job will evaluate 6th row of tensor'
else
write(761,*)'Evaluating full tensor'
end if


write(761,*)'Zero strain'
close(761)

do i=1,natoms

          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)

          xyz_elastic(1,i)=cell_zero(1,1)*xss+cell_zero(1,2)*yss+cell_zero(1,3)*zss
          xyz_elastic(2,i)=cell_zero(2,1)*xss+cell_zero(2,2)*yss+cell_zero(2,3)*zss
          xyz_elastic(3,i)=cell_zero(3,1)*xss+cell_zero(3,2)*yss+cell_zero(3,3)*zss
end do
cell2=cell_zero


call get_stress(cell2,xyz_elastic, stress_plus,f,vol_zero)



if(e1_only)goto 801
if(e2_only)goto 802
if(e3_only)goto 803
if(e4_only)goto 804
if(e5_only)goto 805
if(e6_only)goto 806






801 continue
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



!call matprt(scratch3,3,3,3,3)
cij(1,1)=scratch3(1,1)+1.0d0*stress_in(1,1)
cij(1,2)=scratch3(2,2)-stress_in(1,1)
cij(1,3)=scratch3(3,3)-1.0*stress_in(1,1)
cij(1,4)= scratch3(2,3)
cij(1,5)=(scratch3(1,3)) +1.0d0*stress_in(1,3)
cij(1,6)=(scratch3(1,2) )+1.0d0*stress_in(1,2)
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'C11=',cij(1,1)
write(761,*)'C12=',cij(1,2)
write(761,*)'C13=',cij(1,3)
write(761,*)'C14=',cij(1,4)
write(761,*)'C15=',cij(1,5)
write(761,*)'C16=',cij(1,6)
close(761)
if(e1_only)stop
! call matprt(stress_minus,3,3,3,3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!e2
! do e2
802 continue
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

call get_stress(cell2,xyz_elastic, stress_minus,f,vol2)

stress_minus=stress_elastic
scratch3=auforce_to_bar*1.0d-4*(stress_plus-stress_minus)/(2.0*eps_elas*vol_zero)
scratch3=-scratch3
cij(2,1)=scratch3(1,1)-1.0*stress_in(2,2)
cij(2,2)=scratch3(2,2)+1.0*stress_in(2,2)
cij(2,3)=scratch3(3,3)-1.0*stress_in(2,2)
cij(2,4)=scratch3(2,3)+1.0*stress_in(2,3)
cij(2,5)=scratch3(1,3)
cij(2,6)=scratch3(1,2)+1.0*stress_in(2,1)

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'C21=',cij(2,1)
write(761,*)'C22=',cij(2,2)
write(761,*)'C23=',cij(2,3)
write(761,*)'C24=',cij(2,4)
write(761,*)'C25=',cij(2,5)
write(761,*)'C26=',cij(2,6)
close(761)
if(e2_only)stop



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!e3
! do e3
803 continue
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

call get_stress(cell2,xyz_elastic, stress_minus,f,vol2)


stress_minus=stress_elastic
scratch3=auforce_to_bar*1.0d-4*(stress_plus-stress_minus)/(2.0*eps_elas*vol_zero)
scratch3=-scratch3
cij(3,1)=scratch3(1,1)-1.0*stress_in(3,3)
cij(3,2)=scratch3(2,2)-1.0*stress_in(3,3)
cij(3,3)=scratch3(3,3)+1.0*stress_in(3,3)
cij(3,4)=scratch3(2,3)+1.0*stress_in(3,2)
cij(3,5)=scratch3(1,3)+1.0*stress_in(3,1)
cij(3,6)=scratch3(1,2)

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'C31=',cij(3,1)
write(761,*)'C32=',cij(3,2)
write(761,*)'C33=',cij(3,3)
write(761,*)'C34=',cij(3,4)
write(761,*)'C35=',cij(3,5)
write(761,*)'C36=',cij(3,6)
close(761)
if(e3_only)stop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!e4
! do e4
804 continue
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
cij(4,1)=scratch3(1,1)-1.0*stress_in(2,3)
cij(4,2)=scratch3(2,2) +1.0*stress_in(3,2)-1.0*stress_in(2,3)
cij(4,3)=scratch3(3,3) !!!!+2.0*stress_in(3,3)
cij(4,4)=scratch3(2,3)+.5*stress_in(3,3)+.5*stress_in(2,2)
cij(4,5)=scratch3(1,3)+.5*stress_in(2,1)
cij(4,6)=scratch3(1,2)+.5*stress_in(3,1)
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'C41=',cij(4,1)
write(761,*)'C42=',cij(4,2)
write(761,*)'C43=',cij(4,3)
write(761,*)'C44=',cij(4,4)
write(761,*)'C45=',cij(4,5)
write(761,*)'C46=',cij(4,6)
close(761)
if(e4_only)stop




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!e5
! do e5
805 continue
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
cij(5,1)=scratch3(1,1)+1.0*stress_in(3,1)-1.0*stress_in(1,3)
cij(5,2)=scratch3(2,2)-1.0*stress_in(1,3)
cij(5,3)=scratch3(3,3)
cij(5,4)=scratch3(2,3)+.5*stress_in(1,2)
cij(5,5)=scratch3(1,3)+.5*stress_in(3,3)+.5*stress_in(1,1)
cij(5,6)=scratch3(1,2)+.5*stress_in(3,2)

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'C51=',cij(5,1)
write(761,*)'C52=',cij(5,2)
write(761,*)'C53=',cij(5,3)
write(761,*)'C54=',cij(5,4)
write(761,*)'C55=',cij(5,5)
write(761,*)'C56=',cij(5,6)
close(761)
if(e5_only)stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!e6
! do e6
806 continue
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
cij(6,1)=scratch3(1,1)+1.0*stress_in(2,1)-1.0*stress_in(1,2)
cij(6,2)=scratch3(2,2)
cij(6,3)=scratch3(3,3)-1.0*stress_in(1,2)
cij(6,4)=scratch3(2,3)+.5*stress_in(1,3)
cij(6,5)=scratch3(1,3)+.5*stress_in(2,3)
cij(6,6)=scratch3(1,2)+.5*stress_in(2,2)+.5*stress_in(1,1)
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'C61=',cij(6,1)
write(761,*)'C62=',cij(6,2)
write(761,*)'C63=',cij(6,3)
write(761,*)'C64=',cij(6,4)
write(761,*)'C65=',cij(6,5)
write(761,*)'C66=',cij(6,6)
close(761)
if(e6_only)stop

!!!cij=-cij
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)''
write(761,*)'            Unsymmetrized Elastic constant tensor (Voigt Ordering/GPa):'
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
end subroutine elascon
