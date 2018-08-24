subroutine test_stress
use atom_types_module
use virial_mod
use cell_module
use elastic_module
use constants
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



double precision,dimension(3,3)::strain,voigt,cell2,cell_zero,stress_plus,stress_minus,scratch3,cell_temp,tensor,virial,vref
double precision,dimension(:,:),allocatable::fracs_zero
character(4)::ch
double precision,dimension(6,6)::cij,tensorav,vecs
double precision,dimension(3,natoms)::rxyz0,fxyz,forces


!elas_is_on=.true.


stress_in=-stress_in
xyz=xyz/0.5291772108d0
call cell_volume(cell,vol_zero)
call get_forces(xyz, fxyz, ener, vref)
cell_zero=cell
allocate(fracs_zero(3,natoms))
fracs_zero=fracs
vref=vref*29421.0107637093d0/vol_zero

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'Analytic stress tensor'
close(761)

17 format(f20.6,f20.6,f20.6)

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'                 |-------------------------------------|'
write(761,*)'                 |  Analytic Tensor From External Code |'
write(761,*)'                 |-------------------------------------|'
write(761,17)vref(1,1),vref(1,2),vref(1,3)
write(761,17)vref(2,1),vref(2,2),vref(2,3)
write(761,17)vref(3,1),vref(3,2),vref(3,3)
call space(2)
close(761)








! do e1
voigt=zero
voigt(1,1)=eps_elas
strain=unit_matrix+voigt
cell2=matmul(strain,cell_zero)
do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          
          xyz(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do
cell=cell2
call get_forces(xyz, fxyz, ener, virial)
e_plus=ener
!call matprt(auforce_to_bar*1.0d-4*stress_plus/vol_zero,3,3,3,3)


voigt=zero
voigt(1,1)=-1.0*eps_elas
strain=unit_matrix+voigt
cell2=matmul(strain,cell_zero)
!cell2=cell2/au_to_ang

!call matprt(cell2*au_to_ang,3,3,3,3)
do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss

end do
cell=cell2
call get_forces(xyz, fxyz, ener, virial)

e_minus=ener
tensor(1,1)=(e_plus-e_minus)/(2.0d0*eps_elas)

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'Numerical Pxx is',tensor(1,1)*29421.0107637093d0/vol_zero
close(761)


! do e2
voigt=zero
voigt(2,2)=eps_elas
strain=unit_matrix+voigt
cell2=matmul(strain,cell_zero)
!call matprt(cell2,3,3,3,3)

do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do
cell=cell2
call get_forces(xyz, fxyz, ener, virial)
e_plus=ener


voigt=zero
voigt(2,2)=-1.0*eps_elas
strain=unit_matrix+voigt
cell2=matmul(strain,cell_zero)
do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do
cell=cell2
call get_forces(xyz, fxyz, ener, virial)
e_minus=ener
tensor(2,2)=(e_plus-e_minus)/(2.0d0*eps_elas)

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'Numerical Pyy is',tensor(2,2)*29421.0107637093d0/vol_zero
close(761)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!e3
! do e3
voigt=zero
voigt(3,3)=eps_elas
strain=unit_matrix+voigt
cell2=matmul(strain,cell_zero)
!call matprt(cell2,3,3,3,3)

do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do
cell=cell2
call get_forces(xyz, fxyz, ener, virial)
e_plus=ener



voigt=zero
voigt(3,3)=-1.0*eps_elas
strain=unit_matrix+voigt
cell2=matmul(strain,cell_zero)
do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do
cell=cell2
call get_forces(xyz, fxyz, ener, virial)
e_minus=ener
tensor(3,3)=(e_plus-e_minus)/(2.0d0*eps_elas)
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'Numerical Pzz is',tensor(3,3)*29421.0107637093d0/vol_zero
close(761)





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!e4
! do e4
voigt=zero
voigt(2,3)=eps_elas
voigt(3,2)=eps_elas
strain=unit_matrix+voigt
cell2=matmul(strain,cell_zero)
!call matprt(cell2,3,3,3,3)

do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do
cell=cell2
call get_forces(xyz, fxyz, ener, virial)
e_plus=ener



voigt=zero
voigt(2,3)=-eps_elas
voigt(3,2)=-eps_elas
strain=unit_matrix+voigt
cell2=matmul(strain,cell_zero)
do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do
cell=cell2
call get_forces(xyz, fxyz, ener, virial)
e_minus=ener
tensor(2,3)=(e_plus-e_minus)/(4.0d0*eps_elas)
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'Numerical Pyz is',tensor(2,3)*29421.0107637093d0/vol_zero
close(761)






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!e5
! do e5
voigt=zero
voigt(1,3)=eps_elas
voigt(3,1)=eps_elas
strain=unit_matrix+voigt


cell2=matmul(strain,cell_zero)
!call matprt(cell2,3,3,3,3)

do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do
cell=cell2
call get_forces(xyz, fxyz, ener, virial)
e_plus=ener



voigt=zero
voigt(1,3)=-eps_elas
voigt(3,1)=-eps_elas
strain=unit_matrix+voigt

cell2=matmul(strain,cell_zero)

do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do
cell=cell2
call get_forces(xyz, fxyz, ener, virial)
e_minus=ener
tensor(1,3)=(e_plus-e_minus)/(4.0d0*eps_elas)
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'Numerical Pxz is',tensor(1,3)*29421.0107637093d0/vol_zero
close(761)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!e6
! do e6
voigt=zero
voigt(1,2)=eps_elas
voigt(2,1)=eps_elas
strain=unit_matrix+voigt


cell2=matmul(strain,cell_zero)
!call matprt(cell2,3,3,3,3)

do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do
cell=cell2
call get_forces(xyz, fxyz, ener, virial)
e_plus=ener


voigt=zero
voigt(1,2)=-eps_elas
voigt(2,1)=-eps_elas
strain=unit_matrix+voigt
cell2=matmul(strain,cell_zero)

do i=1,natoms
          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          xyz(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do
cell=cell2
call get_forces(xyz, fxyz, ener, virial)
e_minus=ener
tensor(1,2)=(e_plus-e_minus)/(4.0d0*eps_elas)
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'Numerical Pxy is',tensor(1,2)*29421.0107637093d0/vol_zero
close(761)
stop


!!!cij=-cij

print*,'Unsymmetrized tensor:'
call matprt(cij,6,6,6,6)

tensorav=cij+transpose(cij)
tensorav=tensorav/2.0d0

B=tensorav(1,1)+tensorav(2,2)+tensorav(3,3)+2.0d0*(tensorav(1,2)+tensorav(1,3)+tensorav(2,3))
B=B/9.0d0
print*,''
print*,'Bulk Modulus =',B



shear=tensorav(1,1)+tensorav(2,2)+tensorav(3,3)+3.0d0*(tensorav(4,4)+tensorav(5,5)+tensorav(6,6))-tensorav(1,2)-tensorav(1,3)-tensorav(2,3)
shear=shear/15.0d0
print*,'Shear Modulus =',shear





open(unit=45,file='cijs')
write(45,*)1
write(45,*)tensorav(1,1)
write(45,*)tensorav(1,2)
write(45,*)tensorav(1,3)
write(45,*)tensorav(1,4)
write(45,*)tensorav(3,3)
write(45,*)tensorav(4,4)
write(45,*)B
close(45)

call eig(tensorav,vecs,6,6,6,6)

print*,'eigenvalues'
call matprt(tensorav,6,6,6,6)

!print*,'modes'
!call matprt(vecs,6,6,6,6)







stop



return
end subroutine test_stress
