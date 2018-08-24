subroutine ke_and_t(velocity,kinetic,temperature,scale,scaled_t,pbcflag,kin_vir,temp_tol)
use atom_types_module
use constants
use hugoniot_module
!use control_module
implicit double precision (a-h,o-z)
double precision,dimension(:,:)::velocity
double precision,intent(inout)::kinetic,temperature
double precision,intent(in)::scaled_t,temp_tol
double precision,intent(inout),dimension(3,3)::kin_vir
logical,intent(inout)::scale,pbcflag
logical::doit

freedom=6.0d0
if(pbcflag)freedom=3.0d0

temperature=zero
vx=zero
vy=zero
vz=zero
kin_vir=zero

do i=1,natoms
vx=vx+atom_mass(type_atom(i))*velocity(1,i)**2
vy=vy+atom_mass(type_atom(i))*velocity(2,i)**2
vz=vz+atom_mass(type_atom(i))*velocity(3,i)**2

         kin_vir(1,1)=kin_vir(1,1)+atom_mass(type_atom(i))*velocity(1,i)*velocity(1,i)
         kin_vir(2,1)=kin_vir(2,1)+atom_mass(type_atom(i))*velocity(1,i)*velocity(2,i)
         kin_vir(3,1)=kin_vir(3,1)+atom_mass(type_atom(i))*velocity(1,i)*velocity(3,i)
         kin_vir(2,2)=kin_vir(2,2)+atom_mass(type_atom(i))*velocity(2,i)*velocity(2,i)
         kin_vir(3,2)=kin_vir(3,2)+atom_mass(type_atom(i))*velocity(2,i)*velocity(3,i)
         kin_vir(3,3)=kin_vir(3,3)+atom_mass(type_atom(i))*velocity(3,i)*velocity(3,i)


end do

kin_vir=kin_vir*half
kin_vir(1,2)=kin_vir(2,1)
kin_vir(1,3)=kin_vir(3,1)
kin_vir(2,3)=kin_vir(3,2)


vx=vx*half
vy=vy*half
vz=vz*half

kinetic=vx+vy+vz !kinetic energy in Eo

!temperature=twothirds*kinetic/boltzman
temperature=2.0d0*kinetic/boltzman/(3*natoms-freedom)
temp_unscaled=temperature
if(hugoniot_stat)return


doit=.false.
differ=abs(temperature-scaled_t)
if(differ.ge.temp_tol)doit=.true.
!open(unit=761,file='OUTPUT.DTPOLY',access='append')
!write(761,*)'scale',differ,temperature,scale,doit,temp_tol
!close(761)

if(scale  .or.  doit )then
scale=.false.
factor=dsqrt(scaled_t/temperature)
velocity=velocity*factor
temperature=zero
vx=zero
vy=zero
vz=zero
kin_vir=zero
do i=1,natoms
vx=vx+atom_mass(type_atom(i))*velocity(1,i)**2
vy=vy+atom_mass(type_atom(i))*velocity(2,i)**2
vz=vz+atom_mass(type_atom(i))*velocity(3,i)**2

         kin_vir(1,1)=kin_vir(1,1)+atom_mass(type_atom(i))*velocity(1,i)*velocity(1,i)
         kin_vir(2,1)=kin_vir(2,1)+atom_mass(type_atom(i))*velocity(1,i)*velocity(2,i)
         kin_vir(3,1)=kin_vir(3,1)+atom_mass(type_atom(i))*velocity(1,i)*velocity(3,i)
         kin_vir(2,2)=kin_vir(2,2)+atom_mass(type_atom(i))*velocity(2,i)*velocity(2,i)
         kin_vir(3,2)=kin_vir(3,2)+atom_mass(type_atom(i))*velocity(2,i)*velocity(3,i)
         kin_vir(3,3)=kin_vir(3,3)+atom_mass(type_atom(i))*velocity(3,i)*velocity(3,i)



end do
vx=vx*half
vy=vy*half
vz=vz*half
kinetic=vx+vy+vz
kin_vir=kin_vir*half
kin_vir(1,2)=kin_vir(2,1)
kin_vir(1,3)=kin_vir(3,1)
kin_vir(2,3)=kin_vir(3,2)

!temperature=twothirds*kinetic/boltzman
temperature=2.0d0*kinetic/boltzman/(3*natoms-freedom)
end if 






return
end







