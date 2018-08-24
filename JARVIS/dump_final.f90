subroutine dump_final(cell_in,volume,stress_copy)
use cell_module
use constants
implicit double precision (a-h,o-z)
double precision,dimension(3,3),intent(in)::cell_in,stress_copy
 double precision,intent(inout)::volume

!volume computed as A dot (B X C)



! do B X C
bxci=cell_in(2,2)*cell_in(3,3)-cell_in(3,2)*cell_in(2,3)  !xj*yk- xk*yj
bxcj=cell_in(1,2)*cell_in(3,3)-cell_in(3,2)*cell_in(1,3)  !xi*yk- xk*yi
bxcj=-1.0d0*bxcj
bxck=cell_in(1,2)*cell_in(2,3)-cell_in(2,2)*cell_in(1,3)  !xi*yj- xj*yi
! dot components into A lattice vector
volume=cell_in(1,1)*bxci + cell_in(2,1)*bxcj + cell_in(3,1)*bxck



do i=1,3
tot=0.
do j=1,3
tot=tot+cell_in(j,i)**2
end do
cell_lengths(i)=dsqrt(tot)
end do
tot=0.
!angles
kk=0
do i=1,2
do j=i+1,3
kk=kk+1
tot=cell_in(1,i)*cell_in(1,j)+cell_in(2,i)*cell_in(2,j)+cell_in(3,i)*cell_in(3,j)
cell_angles(kk)=tot/ ( cell_lengths(i) * cell_lengths(j) )
end do
end do

cell_angles=acos(cell_angles)

z=au_to_ang

 write(761,*)'                  Final Optimized Lattice Vectors               '
 write(761,*)' Vector            X-Comp            Y-Comp         Z-Comp          Length'
 write(761,900)'    A->          ',cell_in(1,1)*z,'        ' ,cell_in(2,1)*z,'     ',cell_in(3,1)*z,'      ',cell_lengths(1)*z
 write(761,900)'    B->          ',cell_in(1,2)*z,'        ' ,cell_in(2,2)*z,'     ',cell_in(3,2)*z,'      ',cell_lengths(2)*z
 write(761,900)'    C->          ',cell_in(1,3)*z,'        ' ,cell_in(2,3)*z,'     ',cell_in(3,3)*z,'      ',cell_lengths(3)*z
 write(761,901)'Vector Angles    ',cell_angles(3)*180./pi,'        ',cell_angles(2)*180./pi,'     ',cell_angles(1)*180./pi
 write(761,*)'Volume Ang**(3)',volume*au_to_ang**3
900 format(A,f10.5,A,f10.5,A,f10.5,A,f10.5)
901 format(A,f10.5,A,f10.5,A,f10.5)


open(unit=156,file='cell')
write(156,*)1
write(156,*)cell_lengths(1)*z
write(156,*)cell_lengths(2)*z
write(156,*)cell_lengths(3)*z
write(156,*)cell_angles(1)*180./pi
write(156,*)cell_angles(2)*180./pi
write(156,*)cell_angles(3)*180./pi
do i=1,3
do j=1,3
write(156,*)stress_copy(j,i)
end do 
end do




if(cell_stats)then
av_cell_lengths=av_cell_lengths+cell_lengths
av_cell_angles=av_cell_angles+cell_angles
av_cell_volume=av_cell_volume+volume
end if





end subroutine dump_final
