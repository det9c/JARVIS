subroutine rotate_cell_stress(cell,fxyz,nat,cell_rot,virial)
use control_module
use constants
use rotations
implicit double precision (a-h,o-z)

interface
 subroutine ryrz(x,y,z,phi,theta,xout,yout,zout,natoms)
 double precision,intent(in)::phi,theta
 double precision,dimension(:),intent(inout)::xout,yout,zout
 double precision,intent(in),dimension(:)::x,y,z
 integer,intent(in)::natoms
 end subroutine ryrz

 subroutine rz(x,y,z,phi,xout,yout,zout,natoms)
 double precision,intent(in)::phi
 double precision,dimension(:),intent(inout)::xout,yout,zout
 double precision,intent(in),dimension(:)::x,y,z
 integer,intent(in)::natoms
 end subroutine rz

 subroutine ry(x,y,z,phi,xout,yout,zout,natoms,clock)
 logical::clock
 double precision,intent(in)::phi
 double precision,dimension(:),intent(inout)::xout,yout,zout
 double precision,intent(in),dimension(:)::x,y,z
 integer,intent(in)::natoms
 end subroutine ry

 subroutine angles(x1,y1,z1,theta,phi)
double precision,intent(in)::x1,y1,z1
double precision,intent(inout)::theta,phi
end subroutine angles

subroutine quadrant(x,y,z,phi,sinth)
double precision, intent(inout)::phi,x,y,z,sinth
end subroutine quadrant

end interface


double precision,dimension(3,nat)::fxyz,fxyz_rot
double precision,dimension(nat+3)::fx,fy,fz,fxout,fyout,fzout
integer::nat
double precision,dimension(3,3)::cell,cell_rot,virial,virial_tmp
logical::clock
character(4)::ch

ndim=nat+3
factor=180./pi

do i=1,3
fx(i)=cell_rot(1,i)
fy(i)=cell_rot(2,i)
fz(i)=cell_rot(3,i)
end do
do i=4,ndim
fx(i)=fxyz(1,i-3)
fy(i)=fxyz(2,i-3)
fz(i)=fxyz(3,i-3)
end do

theta=2*pi-270./factor
phi=0.
call ryrz(fx,fy,fz,phi,theta,fxout,fyout,fzout,ndim)
fx=fxout
fy=fyout
fz=fzout
virial_tmp=matmul(rotmat,virial)
virial=matmul(virial_tmp,transpose(rotmat))







!undo second part

call angles(cell(1,2),cell(2,2),cell(3,2),theta,phi)
theta=0.
phi=2.0*pi-rot4
call ryrz(fx,fy,fz,phi,theta,fxout,fyout,fzout,ndim)
fx=fxout
fy=fyout
fz=fzout
virial_tmp=matmul(rotmat,virial)
virial=matmul(virial_tmp,transpose(rotmat))






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute theta and phi for a vector and rotate so that a vector is on z axis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



theta=2.0*pi-rot1
phi=0.  !2.0*pi-rot2
call ryrz(fx,fy,fz,phi,theta,fxout,fyout,fzout,ndim)
fx=fxout
fy=fyout
fz=fzout
virial_tmp=matmul(rotmat,virial)
virial=matmul(virial_tmp,transpose(rotmat))



theta=0.
phi=2.0*pi-rot2
call ryrz(fx,fy,fz,phi,theta,fxout,fyout,fzout,ndim)
fx=fxout
fy=fyout
fz=fzout
virial_tmp=matmul(rotmat,virial)
virial=matmul(virial_tmp,transpose(rotmat))



!do i=1,3
!print*, fx(i)*au_to_ang,fy(i)*au_to_ang,fz(i)*au_to_ang
!end do




do i=4,ndim
fxyz(1,i-3)=fx(i)
fxyz(2,i-3)=fy(i)
fxyz(3,i-3)=fz(i)
end do












end subroutine rotate_cell_stress

