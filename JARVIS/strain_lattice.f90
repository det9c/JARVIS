subroutine strain_lattice
use atom_types_module
use virial_mod
use cell_module
use constants
use hugoniot_module
implicit double precision (a-h,o-z)

interface
   subroutine cell_volume(cell_in,volume)
   double precision,dimension(3,3),intent(in)::cell_in
   double precision,intent(inout)::volume
   end subroutine cell_volume
end interface



double precision,dimension(3,3)::strain,voigt,cell2,cell_zero,sinv!,stress_plus,stress_minus,scratch3,cell_temp
double precision,dimension(:,:),allocatable::fracs_zero
integer,dimension(3)::indx
character(4)::ch

write(761,*)'Initial strain to be applied to lattice:'
call space(1)
!write(761,17)strain_in(1,1),strain_in(1,2),strain_in(1,3)
!write(761,17)strain_in(2,1),strain_in(2,2),strain_in(2,3)
!write(761,17)strain_in(3,1),strain_in(3,2),strain_in(3,3)
17 format(f20.10,f20.10,f20.10)

cell_zero=cell
allocate(fracs_zero(3,natoms))

!fracs_zero=fracs
         cell_copy=cell
         call migs(cell_copy,3,sinv,indx)


         do i=1,natoms
          ssx=sinv(1,1)*xyz(1,i)+sinv(1,2)*xyz(2,i)+sinv(1,3)*xyz(3,i)
          ssy=sinv(2,1)*xyz(1,i)+sinv(2,2)*xyz(2,i)+sinv(2,3)*xyz(3,i)
          ssz=sinv(3,1)*xyz(1,i)+sinv(3,2)*xyz(2,i)+sinv(3,3)*xyz(3,i)
          fracs_zero(1,i)=ssx
          fracs_zero(2,i)=ssy
          fracs_zero(3,i)=ssz
          end do
















! apply strain
strain=unit_matrix+strain_in

write(761,17)strain(1,1),strain(1,2),strain(1,3)
write(761,17)strain(2,1),strain(2,2),strain(2,3)
write(761,17)strain(3,1),strain(3,2),strain(3,3)
call space(2)
cell2=matmul(strain,cell_zero)
write(761,*)'                 |----------------------------------|'
write(761,*)'                 |     Strained Cell Parameters     |'
write(761,*)'                 |----------------------------------|'
call cell_volume(cell2,volume_strain)
hugoniot_volume=volume_strain
call space(1)


do i=1,natoms

          xss=fracs_zero(1,i)
          yss=fracs_zero(2,i)
          zss=fracs_zero(3,i)
          
          xyz(1,i)=cell2(1,1)*xss+cell2(1,2)*yss+cell2(1,3)*zss
          xyz(2,i)=cell2(2,1)*xss+cell2(2,2)*yss+cell2(2,3)*zss
          xyz(3,i)=cell2(3,1)*xss+cell2(3,2)*yss+cell2(3,3)*zss
end do

cell=cell2


! set up scaling matrix for berendsen subroutine
!scale_factor=1.0d0
!do i=1,3
!do j=1,3
!if(strain_in(j,i) .ne. 0.0d0)scale_factor(j,i)=0.0d0
!end do
!end do

!if(hugoniot_stat)scale_factor=1.0d0



call space(2)
!write(761,*)'                 |-------------------------------------|'
!write(761,*)'                 |  Stress Tensor Constraint Matrix    |'
!write(761,*)'                 |-------------------------------------|'
!write(761,17)scale_factor(1,1) ,scale_factor(1,2) ,scale_factor(1,3) 
!write(761,17)scale_factor(2,1) ,scale_factor(2,2) ,scale_factor(2,3) 
!write(761,17)scale_factor(3,1) ,scale_factor(3,2) ,scale_factor(3,3) 
call space(2)



return
end subroutine strain_lattice
