subroutine com_momentum(velocity,total_mass,rxyz)
use atom_types_module
use control_module
use constants
implicit double precision (a-h,o-z)
double precision,dimension(3)::xmom,indx
double precision,dimension(:,:),intent(inout)::velocity,rxyz
double precision,intent(in)::total_mass
double precision,dimension(3,3)::tensor,teninv



xmom=zero

do i=1,natoms
do j=1,3
xmom(j)=xmom(j)+atom_mass(type_atom(i))*velocity(j,i)
end do
end do

xmom=xmom/total_mass

do i=1,natoms
do j=1,3
velocity(j,i)=velocity(j,i)-xmom(j)
end do
end do

if(pbc)return


!remove rotation
! compute center of mass and translate to origin
comx=0.
comy=0.
comz=0.
do i=1,natoms
comx=comx+atom_mass(type_atom(i))*rxyz(1,i)
comy=comy+atom_mass(type_atom(i))*rxyz(2,i)
comz=comz+atom_mass(type_atom(i))*rxyz(3,i)
end do
comx=comx/total_mass
comy=comy/total_mass
comz=comz/total_mass

do i=1,natoms
rxyz(1,i)=rxyz(1,i)-comx
rxyz(2,i)=rxyz(2,i)-comy
rxyz(3,i)=rxyz(3,i)-comz
end do

! get angular momentum and intertia tensor.  tensor code taken from NDDO code
angx=0.
angy=0.
angz=0.

tensor=0.
do i=1,natoms
zcharge=atom_mass(type_atom(i))
tensor(1,1)=tensor(1,1)+zcharge* ( rxyz(2,i)**2 + rxyz(3,i)**2)
tensor(2,2)=tensor(2,2)+zcharge* ( rxyz(1,i)**2 + rxyz(3,i)**2)
tensor(3,3)=tensor(3,3)+zcharge* ( rxyz(1,i)**2 + rxyz(2,i)**2)
tensor(2,1)=tensor(2,1)-zcharge*rxyz(1,i)*rxyz(2,i)
tensor(3,1)=tensor(3,1)-zcharge*rxyz(1,i)*rxyz(3,i)
tensor(3,2)=tensor(3,2)-zcharge*rxyz(3,i)*rxyz(2,i)
angx=angx+atom_mass(type_atom(i)) * ( rxyz(2,i) * velocity(3,i) - rxyz(3,i) * velocity(2,i))
angy=angy-atom_mass(type_atom(i)) * ( rxyz(1,i) * velocity(3,i) - rxyz(3,i) * velocity(1,i))
angz=angz+atom_mass(type_atom(i)) * ( rxyz(1,i) * velocity(2,i) - rxyz(2,i) * velocity(1,i))
end do

tensor(1,2)=tensor(2,1)
tensor(1,3)=tensor(3,1)
tensor(2,3)=tensor(3,2)


call migs(tensor,3,teninv,indx)

!compute angular speed around center of mass

wx=teninv(1,1)*angx + teninv(1,2)*angy + teninv(1,3)*angz
wy=teninv(2,1)*angx + teninv(2,2)*angy + teninv(2,3)*angz
wz=teninv(3,1)*angx + teninv(3,2)*angy + teninv(3,3)*angz

! compute linear speed v = w x r and subtract from velocity

do i=1,natoms
xlin = wy * rxyz(3,i) - wz * rxyz(2,i)
ylin = wx * rxyz(3,i) - wz * rxyz(1,i)
ylin=-ylin
zlin = wx * rxyz(2,i) - wy * rxyz(1,i)
velocity(1,i)=velocity(1,i)-xlin
velocity(2,i)=velocity(2,i)-ylin
velocity(3,i)=velocity(3,i)-zlin
end do


! move atoms back
do i=1,natoms
rxyz(1,i)=rxyz(1,i)+comx
rxyz(2,i)=rxyz(2,i)+comy
rxyz(3,i)=rxyz(3,i)+comz
end do



return
end subroutine com_momentum
