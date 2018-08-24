subroutine rotate_cell(cell,rxyz0,nat,cell_rot,rxyz_rot)
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


double precision,dimension(3,nat)::rxyz0,rxyz_rot
double precision,dimension(nat+3)::x,y,z,xout,yout,zout
integer::nat
double precision,dimension(3,3)::cell,cell_rot
character(4)::ch

ndim=nat+3
factor=180./pi

do i=1,3
x(i)=cell(1,i)
y(i)=cell(2,i)
z(i)=cell(3,i)
end do


do i=4,ndim
x(i)=rxyz0(1,i-3)
y(i)=rxyz0(2,i-3)
z(i)=rxyz0(3,i-3)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute theta and phi for a vector and rotate so that a vector is on z axis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call angles(x(1),y(1),z(1),theta,phi)
call ryrz(x,y,z,phi,theta,xout,yout,zout,ndim)
x=xout
y=yout
z=zout
rot1=theta
rot2=phi







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute phi for b vector and rotate to +y side of x-axis
! rotating ccw about z so add 270 to complete orientation 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call angles(x(2),y(2),z(2),theta,phi)
theta=0.
phi=phi+270./factor
call ryrz(x,y,z,phi,theta,xout,yout,zout,ndim)
x=xout
y=yout
z=zout
rot3=theta
rot4=phi










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! rotating cw about y so do 270 degree to put a vector on x axis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

theta=270./factor
phi=0.
call ryrz(x,y,z,phi,theta,xout,yout,zout,ndim)
x=xout
y=yout
z=zout
rot5=theta
rot6=phi









do i=1,3
cell_rot(1,i)=x(i)
cell_rot(2,i)=y(i)
cell_rot(3,i)=z(i)
end do




do i=4,ndim
rxyz_rot(1,i-3)=x(i)
rxyz_rot(2,i-3)=y(i)
rxyz_rot(3,i-3)=z(i)
end do




end subroutine rotate_cell







subroutine angles(x1,y1,z1,theta,phi)
implicit double precision (a-h,o-z)
double precision,intent(in)::x1,y1,z1
double precision,intent(inout)::theta,phi

pi=3.14159d0
pi=dacos(-1.0d0)
one=1.0d0
zero=0.0d0
two=2.0d0
three=3.0d0

rsq=x1*x1 + y1*y1 + z1*z1
r=dsqrt(rsq)
costh=z1/r
if(abs(costh).gt.one)costh=dsign(one,costh)
 theta=dacos(costh)
 sinth=dsin(theta)
 if(abs(sinth)<1D-15)sinth=zero
 if(abs(sinth)==zero)then  ! there is a problem here with machine precision. 1d-20 effectively turns
! this check off
     phi=zero
     if(z1.gt.zero)then
     theta=zero
     else
     theta=pi
     end if
! if not on z axis then...
else
         cosphi=x1/(r*sinth)
         if(abs(cosphi).gt.one)cosphi=dsign(one,cosphi)
         phi=dacos(cosphi)
         cosphi=dcos(phi)
!         print*,'so phi is',phi*180.0d0/pi

! if on y axis or in plane of y axis then
!        rdiff=abs(phi)-(pi/two)
!    if(abs(phi).eq.(pi/two))then
!        if(abs(rdiff).lt.1.0d-6)then !BBBBBBBBBBBBBb
         if(abs(cosphi)<1D-15)cosphi=zero
         if(cosphi==zero)then
              if(y1.gt.zero)then
              phi=pi/two
              else
              phi=three*pi/two
              end if
!  elseif(abs(cosphi).lt.0.99999999999999999999d0)then
         elseif(y1==zero.and.z1==zero)then
               if(x1<zero)then
                  phi=pi
                  else
                  phi=zero
                  end if
         else
!            print*,'call quad',x1,y1,z1,phi*180.0d0/pi
         call quadrant(x1,y1,z1,phi,sinth)
!         print*,'after',phi*180.0d0/pi
          end if

end if

end subroutine angles

subroutine quadrant(x,y,z,phi,sinth)
implicit double precision (a-h,o-z)
double precision, intent(inout)::phi,x,y,z,sinth
one=1.0d0
zero=0.0d0
two=2.0d0
pi=3.14159d0
pi=dacos(-1.0d0)
!print*,'in quadrant',phi*180.0d0/pi
!check to see if in x-y plane
!if(abs(sinth).gt.0.99999999999999999999999999d0)then

if(abs(sinth)==one)then
!   print*,'in xy plane'
        if((x.gt.zero).and.(y.gt.zero))then
        return
        elseif((x.lt.zero).and.(y.gt.zero))then
        return
        else
        phi=two*pi-phi
        end if

else

   if((x.gt.zero).and.(y.gt.zero).and.(z.gt.zero))then
   return
   elseif((x.lt.zero).and.(y.gt.zero).and.(z.gt.zero))then
   return
   elseif((x.lt.zero).and.(y.lt.zero).and.(z.gt.zero))then
   iquad=3
   elseif((x.gt.zero).and.(y.lt.zero).and.(z.gt.zero))then
   iquad=4
   elseif((x.gt.zero).and.(y.gt.zero).and.(z.lt.zero))then
   iquad=5
   elseif((x.lt.zero).and.(y.gt.zero).and.(z.lt.zero))then
   iquad=6
   elseif((x.lt.zero).and.(y.lt.zero).and.(z.lt.zero))then
   iquad=7
   elseif((x.gt.zero).and.(y.lt.zero).and.(z.lt.zero))then
   iquad=8
   end if

   if(iquad==3)then
   phi=two*pi-phi
   elseif(iquad==4)then
   phi=two*pi-phi
   elseif(iquad==7)then
   phi=two*pi-phi
   elseif(iquad==8)then
   phi=two*pi-phi
   end if

end if

end subroutine quadrant

 subroutine ryrz(x,y,z,phi,theta,xout,yout,zout,natoms)
 use rotations
 implicit double precision (a-h,o-z)
 double precision,intent(in)::phi,theta
 double precision,dimension(:),intent(inout)::xout,yout,zout
! double precision,dimension(3,3)::rotmat
 double precision,intent(in),dimension(:)::x,y,z
 integer,intent(in)::natoms
zero=0.0d0
xout=zero
yout=zero
zout=zero
! compute rotation matrix
cosphi=dcos(phi)

sinphi=dsin(phi)
costh=dcos(theta)
sinth=dsin(theta)
rotmat(1,1)=cosphi*costh
rotmat(1,2)=-sinphi
rotmat(1,3)=cosphi*sinth
rotmat(2,1)=sinphi*costh
rotmat(2,2)=cosphi
rotmat(2,3)=sinphi*sinth
rotmat(3,1)=-sinth
rotmat(3,2)=zero
rotmat(3,3)=costh
!i borrowed this code from my nddo code which was rotating integrals
!from local framework to molecular framework.  here is am doing
!the inverse so to use same code, i must take transpose
rotmat=transpose(rotmat)

do i=1,natoms
xout(i)=rotmat(1,1)*x(i) + rotmat(1,2)*y(i)+rotmat(1,3)*z(i)
yout(i)=rotmat(2,1)*x(i) + rotmat(2,2)*y(i)+rotmat(2,3)*z(i)
zout(i)=rotmat(3,1)*x(i) + rotmat(3,2)*y(i)+rotmat(3,3)*z(i)
end do

end subroutine ryrz

 subroutine rz(x,y,z,phi,xout,yout,zout,natoms)
 implicit double precision (a-h,o-z)
 double precision,intent(in)::phi
 double precision,dimension(:),intent(inout)::xout,yout,zout
 double precision,dimension(3,3)::rotmat
 double precision,intent(in),dimension(:)::x,y,z
 integer,intent(in)::natoms
zero=0.0d0
xout=zero
yout=zero
zout=zero
! compute rotation matrix
cosphi=dcos(phi)
sinphi=dsin(phi)
rotmat(1,1)=cosphi
rotmat(1,2)=sinphi
rotmat(1,3)=zero
rotmat(2,1)=-sinphi
rotmat(2,2)=cosphi
rotmat(2,3)=zero
rotmat(3,1)=zero
rotmat(3,2)=zero
rotmat(3,3)=1.0d0

do i=1,natoms
xout(i)=rotmat(1,1)*x(i) + rotmat(1,2)*y(i)+rotmat(1,3)*z(i)
yout(i)=rotmat(2,1)*x(i) + rotmat(2,2)*y(i)+rotmat(2,3)*z(i)
zout(i)=rotmat(3,1)*x(i) + rotmat(3,2)*y(i)+rotmat(3,3)*z(i)
end do

end subroutine rz

 subroutine ry(x,y,z,phi,xout,yout,zout,natoms,clock)
 implicit double precision (a-h,o-z)
 logical::clock
 double precision,intent(in)::phi
 double precision,dimension(:),intent(inout)::xout,yout,zout
 double precision,dimension(3,3)::rotmat
 double precision,intent(in),dimension(:)::x,y,z
 integer,intent(in)::natoms
zero=0.0d0
xout=zero
yout=zero
zout=zero
! compute rotation matrix for ccw rotation
cosphi=dcos(phi)
sinphi=dsin(phi)
rotmat(1,1)=cosphi
rotmat(1,2)=zero
rotmat(1,3)=-sinphi
rotmat(2,1)=zero
rotmat(2,2)=1.0d0
rotmat(2,3)=zero
rotmat(3,1)=sinphi
rotmat(3,2)=zero
rotmat(3,3)=cosphi

if(clock)rotmat=transpose(rotmat)

do i=1,natoms
xout(i)=rotmat(1,1)*x(i) + rotmat(1,2)*y(i)+rotmat(1,3)*z(i)
yout(i)=rotmat(2,1)*x(i) + rotmat(2,2)*y(i)+rotmat(2,3)*z(i)
zout(i)=rotmat(3,1)*x(i) + rotmat(3,2)*y(i)+rotmat(3,3)*z(i)
end do

end subroutine ry





