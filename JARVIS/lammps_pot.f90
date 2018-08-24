subroutine lammps_pot(atmsym,nat,rxyz0, fxyz, ener, virial)
  use control_module
  use constants
  use cell_module
  use atom_types_module
implicit double precision (a-h,o-z)
interface
subroutine rotate_cell(cell,rxyz0,nat,cell_rot,rxyz_rot)
double precision,dimension(3,nat)::rxyz0,rxyz_rot
integer::nat
double precision,dimension(3,3)::cell,cell_rot
end subroutine rotate_cell

   subroutine cell_volume(cell_in,volume)
   double precision,dimension(3,3),intent(in)::cell_in
   double precision,intent(inout)::volume
   end subroutine cell_volume

subroutine rotate_cell_stress(cell,rxyz0,nat,cell_rot,virial)
double precision,dimension(3,nat)::rxyz0
integer::nat
double precision,dimension(3,3)::cell,cell_rot,virial
end subroutine rotate_cell_stress

end interface

  double precision,dimension(3,nat)::rxyz0,fxyz,rxyz_rot
  integer::nat
  character(4),dimension(nat)::atmsym
  double precision::ener
   integer,dimension(nat)::ilabel
   double precision,dimension(3,3)::virial,cell_lammps,cell_rot
 character(4)::ch

open(unit=14,file='lammps.conf')
write(14,*)nat,'  ARL MD input'
write(14,*)nat,'  atoms'
write(14,*)''
write(14,*)ntypes,'  atom types'
write(14,*)''

! here we will have to rotate cell properly, but for now,

if(cell_opt)call rotate_cell(cell,rxyz0,nat,cell_rot,rxyz_rot)

!call cell_volume(cell_rot,volume)
cell_lammps=transpose(cell_rot)

cell_lammps=cell_lammps*au_to_ang
if(cell_opt)then
write(14,*)0.0,cell_lammps(1,1),'xlo xhi'
write(14,*)0.0,cell_lammps(2,2),'ylo yhi'
write(14,*)0.0,cell_lammps(3,3),'zlo zhi'
write(14,11)cell_lammps(2,1),cell_lammps(3,1),cell_lammps(3,2),' xy xz yz '
else
write(14,*)-100.,100.,'xlo xhi'
write(14,*)-100.,100.,'ylo yhi'
write(14,*)-100.,100.,'zlo zhi'
end if

write(14,*)''
write(14,*)'Masses'
write(14,*)''
do i=1,ntypes
write(14,*)i,atom_mass(i)
end do

write(14,*)''
write(14,*)'Atoms'
write(14,*)''

if(.not. cell_opt)rxyz_rot=rxyz0

do i=1,nat
!write(14,10)i,type_atom(i),atom_charge(type_atom(i)),rxyz0(1,i)*au_to_ang,rxyz0(2,i)*au_to_ang,rxyz0(3,i)*au_to_ang
write(14,10)i,type_atom(i),atom_charge(type_atom(i)),rxyz_rot(1,i)*au_to_ang,rxyz_rot(2,i)*au_to_ang,rxyz_rot(3,i)*au_to_ang
end do

10 format(i6,i4,f20.10,f25.15,f25.15,f25.15)
11 format(3f25.15,A)

close(14)

call system('./run_lammps')

open(unit=1,file='lammps.output')
do i=1,10000
read(1,*)ch
!print*,ch
if(ch == 'Step')exit
end do
! will change this line to read different data
read(1,*)junk,vol,ener,virial(1,1),virial(2,2),virial(3,3),virial(1,2),virial(1,3),virial(2,3)

close(1)

virial(2,1)=virial(1,2)
virial(3,1)=virial(1,3)
virial(3,2)=virial(2,3)

!call matprt(virial,3,3,3,3)

if(lreals)then
factor=.98692d4  !    atm per gpa
factor2=au_to_ang/627.51d0
factor3=627.51d0
!print*,'converting from REAL lammps units'
elseif(metal)then
 factor=10000.0d0  !    bars per gpa
factor2=au_to_ang/27.21138386d0
factor3=27.21138386d0
!print*,'converting from METALS lammps units'
else
print*,'need to set metals or real units for lammps'
stop
end if
virial=virial/factor ! pressure in gpa
virial=virial/29421.0107637093d0 ! in a.u.
vol=vol/au_to_ang**3
virial=virial*vol
ener=ener/factor3

open(unit=1,file='forces')
read(1,*)ch
read(1,*)ch
read(1,*)ch
read(1,*)ch
read(1,*)ch
read(1,*)ch
read(1,*)ch
read(1,*)ch
read(1,*)ch
do i=1,nat
read(1,*)fxyz(1,i),fxyz(2,i),fxyz(3,i)
!print*,fxyz(1,i),fxyz(2,i),fxyz(3,i)
end do
fxyz=fxyz*factor2

close(1)

if(cell_opt)call rotate_cell_stress(cell,fxyz,nat,cell_rot,virial)

!do i=1,nat
!print*,fxyz(1,i),fxyz(2,i),fxyz(3,i)
!end do

!do i=1,3
!do j=1,3
!if(j.ne.i)virial(j,i)=0.
!end do
!end do
!call matprt(virial,3,3,3,3)

return
end

