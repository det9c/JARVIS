subroutine neighbor_list(natoms,rxyz,ax,ay,az)
use cell_module
use control_module
use neighbors_module

implicit double precision (a-h,o-z) 
double precision,dimension(3,natoms)::rxyz
integer,intent(in)::natoms
integer,allocatable,dimension(:)::old_list

! allocate num_neighbors,ifirst_neighbor in input.f90


if(allocated(ipair_list))deallocate(ipair_list)
allocate(ipair_list(1))
allocate(old_list(1))
irow=0

do i=1,natoms
  near=0
  ifirst_neighbor(i)=irow+1
   do j=1,natoms 


             rxij=rxyz(1,i)-rxyz(1,j)
             ryij=rxyz(2,i)-rxyz(2,j)
             rzij=rxyz(3,i)-rxyz(3,j)
             if(pbc)then
               if(cubic)then
                 rxij=rxij-ax*anint(rxij/ax)
                 ryij=ryij-ay*anint(ryij/ay)
                 rzij=rzij-az*anint(rzij/az)
              else
         
                 sxij=fracs(1,i)-fracs(1,j)
                 syij=fracs(2,i)-fracs(2,j)
                 szij=fracs(3,i)-fracs(3,j)
                 sxij=sxij-anint(sxij)
                 syij=syij-anint(syij)
                 szij=szij-anint(szij)
                 rxij=cell(1,1)*sxij+cell(1,2)*syij+cell(1,3)*szij
                 ryij=cell(2,1)*sxij+cell(2,2)*syij+cell(2,3)*szij
                 rzij=cell(3,1)*sxij+cell(3,2)*syij+cell(3,3)*szij
              end if
            end if
         rijsq=rxij*rxij + ryij*ryij + rzij*rzij


         if(rijsq.lt.cutoff_sq)then
            near=near+1 
            irow=irow+1
            deallocate(ipair_list)
            allocate(ipair_list(irow))
            ipair_list(1:irow-1)=old_list(1:irow-1)
            ipair_list(irow)=j
            deallocate(old_list)
            allocate(old_list(irow))
            old_list=ipair_list
         end if
   end do
num_neighbors(i)=near
end do

deallocate(old_list)


!do j=1,natoms
!print*,'atom j',j
!print*,'-----------------'
!do i=1,num_neighbors(j)
!print*,ipair_list(ifirst_neighbor(j)+i-1)
!end do
!end do



end subroutine neighbor_list
























