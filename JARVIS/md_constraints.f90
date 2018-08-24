subroutine md_constraints
use atom_types_module
use control_module
use constants
use constraints_module
implicit double precision (a-h,o-z)
interface
   subroutine ke_and_t(velocity,kinetic,temperature,scale,scaled_t,pbcflag)
     double precision,dimension(:,:)::velocity
     double precision,intent(inout)::kinetic,temperature
     double precision,intent(in)::scaled_t
     logical,intent(inout)::scale,pbcflag
   end subroutine ke_and_t

   subroutine com_momentum(velocity,total_mass,rxyz)
     double precision,dimension(:,:),intent(inout)::velocity,rxyz
     double precision,intent(in)::total_mass
   end subroutine com_momentum
end interface
double precision,dimension(3,natoms)::velocity,rtdt,rt,acc,veltdt,forces
double precision::kinetic,k,kx,ky,kz
!!! test stuff
double precision,dimension(3,natoms)::q,gforce
double precision,dimension(100)::evec
!!! end test stuff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   set initial conditions

idum=-5
call rninit(idum)
! assign random initial velocities and compute total mass
   total_mass=zero
   do i=1,natoms
   do j=1,3
   velocity(j,i)=2.0d0*urand()-1.0d0 ! units are chosen to be Bohr/ps 
   end do
   total_mass=total_mass+atom_mass(type_atom(i))
   end do

xyz=xyz/au_to_ang !convert input coordinates to bohr
call com_momentum(velocity,total_mass,xyz)

scale=.true.
call ke_and_t(velocity,kinetic,temperature,scale,temp_input,pbc)

!                done with initial conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



open(unit=15,file='md_history.xyz') ! file to hold trajectory 



call space(5)

print*,'***************************************************************'
print*,'***************************************************************'
print*,'        VELOCITY-VERLET  MOLECULAR DYNAMICS ALGORITHM        '
print*,'                    Velocity Scaling'
print*,'               with intermolecular constraints               '
call space(2)
write(*,100)'              -> Number of timesteps:     ',isteps
write(*,110)'              -> Time step (ps):          ',tstep
write(*,110)'              -> Initial temperature (K): ',temperature
write(*,110)'              -> Cutoff (Angstroms):      ',cutoff       
print*,'             -> PBC flag                 ',pbc
if(pbc)write(*,120)'              -> Constraining center of mass translation'
if(.not. pbc)write(*,120)'              -> Constraining center of mass translation and rotation'
print*,'***************************************************************'
print*,'***************************************************************'
100 format(A,I10)
110 format(A,F10.4)
120 format(A)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up for MD loop. 
rt=xyz ! copy input coordinates into MD working array
rtdt=zero !initialize some MD arrays 
acc=zero
veltdt=zero 
call get_forces(rt, forces, ener, virial) !get initial forces and energy (hartrees/bohr and hartree)
forces=forces/hartree_per_eint !convert to internal units
! get initial accelerations from initial forces
do i=1,natoms
do j=1,3
acc(j,i)=forces(j,i)/atom_mass(type_atom(i))


end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




constraint_values=constraint_values/au_to_ang


do i=1,isteps ! MD loop
    rtdt=rt+velocity*tstep+.5d0*acc*tstep*tstep ! compute new coordinates
    gforce=rtdt


! now do constraint correction
do jj=1,1000000

do ii=1,num_constraints

 !evaluate error
 qx=0.
 qy=0.
 qz=0.
 rx=rt(1,iconstraint_pairs(ii,1))-rt(1,iconstraint_pairs(ii,2))
 ry=rt(2,iconstraint_pairs(ii,1))-rt(2,iconstraint_pairs(ii,2))
 rz=rt(3,iconstraint_pairs(ii,1))-rt(3,iconstraint_pairs(ii,2))
 xitemp=rtdt(1,iconstraint_pairs(ii,1))
 yitemp=rtdt(2,iconstraint_pairs(ii,1))
 zitemp=rtdt(3,iconstraint_pairs(ii,1))
 xjtemp=rtdt(1,iconstraint_pairs(ii,2))
 yjtemp=rtdt(2,iconstraint_pairs(ii,2))
 zjtemp=rtdt(3,iconstraint_pairs(ii,2))


 do icon=1,100
 xitemp=xitemp-tstep*qx/atom_mass(type_atom(iconstraint_pairs(ii,1)))
 yitemp=yitemp-tstep*qy/atom_mass(type_atom(iconstraint_pairs(ii,1)))
 zitemp=zitemp-tstep*qz/atom_mass(type_atom(iconstraint_pairs(ii,1)))
 xjtemp=xjtemp+tstep*qx/atom_mass(type_atom(iconstraint_pairs(ii,2)))
 yjtemp=yjtemp+tstep*qy/atom_mass(type_atom(iconstraint_pairs(ii,2)))
 zjtemp=zjtemp+tstep*qz/atom_mass(type_atom(iconstraint_pairs(ii,2)))
 sx=xitemp-xjtemp
 sy=yitemp-yjtemp
 sz=zitemp-zjtemp
 snorm=sx**2+sy**2+sz**2
 error=snorm-constraint_values(ii)**2

 if(dabs(error).lt.1d-8)then
 rtdt(1,iconstraint_pairs(ii,1))=xitemp
 rtdt(2,iconstraint_pairs(ii,1))=yitemp
 rtdt(3,iconstraint_pairs(ii,1))=zitemp
 rtdt(1,iconstraint_pairs(ii,2))=xjtemp
 rtdt(2,iconstraint_pairs(ii,2))=yjtemp
 rtdt(3,iconstraint_pairs(ii,2))=zjtemp
 exit
 end if
 sdotr=sx*rx+sy*ry+sz*rz
 g=(snorm-constraint_values(ii)**2)/(2.0*tstep*sdotr)
 g=g/(  (1./atom_mass(type_atom(iconstraint_pairs(ii,1)))) + (1./atom_mass(type_atom(iconstraint_pairs(ii,2))))   ) 
 qx=g*rx
 qy=g*ry
 qz=g*rz
 end do


end do

!evaluate constraint vector
rnorm=0.d0
do ik=1,num_constraints
rx=rtdt(1,iconstraint_pairs(ik,1))-rtdt(1,iconstraint_pairs(ik,2))
ry=rtdt(2,iconstraint_pairs(ik,1))-rtdt(2,iconstraint_pairs(ik,2))
rz=rtdt(3,iconstraint_pairs(ik,1))-rtdt(3,iconstraint_pairs(ik,2))
evec(ik)=rx**2 + ry**2 + rz**2
evec(ik)=evec(ik)-constraint_values(ik)**2
rnorm=rnorm+evec(ik)**2
end do


if(dsqrt(rnorm) .lt. 1d-3)then
call space(1)
do ik=1,num_constraints
rx=rtdt(1,iconstraint_pairs(ik,1))-rtdt(1,iconstraint_pairs(ik,2))
ry=rtdt(2,iconstraint_pairs(ik,1))-rtdt(2,iconstraint_pairs(ik,2))
rz=rtdt(3,iconstraint_pairs(ik,1))-rtdt(3,iconstraint_pairs(ik,2))
rnum=rx**2 + ry**2 + rz**2
rnum=dsqrt(rnum)
!print*,ik,rnum*au_to_ang,constraint_values(ik)*au_to_ang
end do
!print*,'error vec norm',dsqrt(rnorm)
rnorm2=rnorm
exit
end if

end do
!stop

! get the g matrix for virial and velocity update
gforce=rtdt-gforce ! this is .5*tstep**2*(G/m)
gforce=gforce/(.5d0*tstep*tstep) ! G/m. will have to scale for virial  


!    veltdt=velocity+0.5d0*tstep*acc             ! velocity at half timestep 
    veltdt=velocity+0.5d0*tstep*(acc+gforce) 
    call get_forces(rtdt, forces, ener, virial)         !hartrees/bohr and energy in hartree
    forces=forces/hartree_per_eint              !convert forces to Eo/Bohr
                                                ! convert forces to accelerations in this loop 
    do ii=1,natoms
    do j=1,3
      acc(j,ii)=forces(j,ii)/atom_mass(type_atom(ii))
   end do
   end do
                                                ! accelerations are done 
   veltdt=veltdt+.5d0 * tstep * acc             ! update half time step velocities to current velocities 

! so veltdt are the unconstrained velocities
do jj=1,1000000

do ii=1,num_constraints

 !evaluate error
 kx=0.
 ky=0.
 kz=0.
 rx=rtdt(1,iconstraint_pairs(ii,1))-rtdt(1,iconstraint_pairs(ii,2))
 ry=rtdt(2,iconstraint_pairs(ii,1))-rtdt(2,iconstraint_pairs(ii,2))
 rz=rtdt(3,iconstraint_pairs(ii,1))-rtdt(3,iconstraint_pairs(ii,2))
 vxitemp=veltdt(1,iconstraint_pairs(ii,1))
 vyitemp=veltdt(2,iconstraint_pairs(ii,1))
 vzitemp=veltdt(3,iconstraint_pairs(ii,1))
 vxjtemp=veltdt(1,iconstraint_pairs(ii,2))
 vyjtemp=veltdt(2,iconstraint_pairs(ii,2))
 vzjtemp=veltdt(3,iconstraint_pairs(ii,2))


 do icon=1,100
 vxitemp=vxitemp-kx/atom_mass(type_atom(iconstraint_pairs(ii,1)))
 vyitemp=vyitemp-ky/atom_mass(type_atom(iconstraint_pairs(ii,1)))
 vzitemp=vzitemp-kz/atom_mass(type_atom(iconstraint_pairs(ii,1)))
 vxjtemp=vxjtemp+kx/atom_mass(type_atom(iconstraint_pairs(ii,2)))
 vyjtemp=vyjtemp+ky/atom_mass(type_atom(iconstraint_pairs(ii,2)))
 vzjtemp=vzjtemp+kz/atom_mass(type_atom(iconstraint_pairs(ii,2)))
 sx=vxitemp-vxjtemp
 sy=vyitemp-vyjtemp
 sz=vzitemp-vzjtemp
 snorm=sx*rx + sy*ry+ sz*rz
 error=snorm
 if(dabs(error).lt.1d-4)then
 veltdt(1,iconstraint_pairs(ii,1))=vxitemp
 veltdt(2,iconstraint_pairs(ii,1))=vyitemp
 veltdt(3,iconstraint_pairs(ii,1))=vzitemp
 veltdt(1,iconstraint_pairs(ii,2))=vxjtemp
 veltdt(2,iconstraint_pairs(ii,2))=vyjtemp
 veltdt(3,iconstraint_pairs(ii,2))=vzjtemp
 exit
 end if


 k=snorm/ constraint_values(ii)**2
 k=k/ ( 1./atom_mass(type_atom(iconstraint_pairs(ii,1))) + 1./atom_mass(type_atom(iconstraint_pairs(ii,2))))  
 kx=k*rx
 ky=k*ry
 kz=k*rz
 end do


end do

!evaluate constraint vector
rnorm=0.d0
do ik=1,num_constraints
rx=rtdt(1,iconstraint_pairs(ik,1))-rtdt(1,iconstraint_pairs(ik,2))
ry=rtdt(2,iconstraint_pairs(ik,1))-rtdt(2,iconstraint_pairs(ik,2))
rz=rtdt(3,iconstraint_pairs(ik,1))-rtdt(3,iconstraint_pairs(ik,2))
vrx=veltdt(1,iconstraint_pairs(ik,1))-veltdt(1,iconstraint_pairs(ik,2))
vry=veltdt(2,iconstraint_pairs(ik,1))-veltdt(2,iconstraint_pairs(ik,2))
vrz=veltdt(3,iconstraint_pairs(ik,1))-veltdt(3,iconstraint_pairs(ik,2))
evec(ik)=rx*vrx + ry*vry + rz*vrz
rnorm=rnorm+evec(ik)**2
end do


if(dsqrt(rnorm) .lt. 1d-2)then
rnorm3=rnorm
exit
end if

end do






   call  com_momentum(veltdt,total_mass,rtdt)   ! remove translation/rotation

   if(mod(i,iscale_counter).eq.0)scale=.true.   !scale to T (subroutine sets scale to false if scaling is done)
   call ke_and_t(veltdt,kinetic,temperature,scale,temp_input,pbc)

   call space(2)                                ! print information
   print*,'                           Timestep # ',i
   print*,'           **************************************************'
   print*,'            Simulation Time   (fs):',float(i)*tstep*1000.0d0
   print*,'            Potential Energy (a.u):',ener
   print*,'            Kinetic Energy   (a.u):' ,kinetic*hartree_per_eint
   print*,'            Temperature        (K):',temperature
   print*,'            Constraint Vector Norm(Rij):',dsqrt(rnorm2)
   print*,'            Constraint Vector Norm(Vij):',dsqrt(rnorm3)
   velocity=veltdt                              ! save new velocities as old velocities
   rt=rtdt                                      ! save new coordinates as old coordinates

  if(i.eq.1 .or. mod(i,iwrite).eq.0)then
   write(15,*)natoms
   write(15,*)natoms
   do kk=1,natoms
   write(15,50)input_names(kk),rt(1,kk)*au_to_ang,rt(2,kk)*au_to_ang,rt(3,kk)*au_to_ang
   end do
   end if
   50 format(A,f10.5,f10.5,f10.5)


end do ! end loop over time steps



call space(2)
print*,'               ******************************************'
print*,'               ******************************************'
print*,'               ******************************************'
print*,'               **********MD run has completed************'
print*,'               ******************************************'
print*,'               ******************************************'
print*,'               ******************************************'
call space(1)
print*,'      Summary of final bond lengths and constraints'
print*,'     Constraint    Actual Value           Constraint Value'

do ik=1,num_constraints
rx=rtdt(1,iconstraint_pairs(ik,1))-rtdt(1,iconstraint_pairs(ik,2))
ry=rtdt(2,iconstraint_pairs(ik,1))-rtdt(2,iconstraint_pairs(ik,2))
rz=rtdt(3,iconstraint_pairs(ik,1))-rtdt(3,iconstraint_pairs(ik,2))
rnum=rx**2 + ry**2 + rz**2
rnum=dsqrt(rnum)
print*,ik,'   ',rnum*au_to_ang,constraint_values(ik)*au_to_ang
end do








end subroutine md_constraints



!   if(i.eq.1 .or. mod(i,100).eq.0)then
!   write(15,*)natoms
!   write(15,*)natoms
!   do kk=1,natoms
!   write(15,50)input_names(kk),rt(1,kk)*au_to_ang,rt(2,kk)*au_to_ang,rt(3,kk)*au_to_ang
!   end do
!   end if
!   50 format(A,f10.5,f10.5,f10.5)
