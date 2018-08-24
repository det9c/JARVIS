subroutine md_constraints_nose
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
double precision,dimension(3,natoms)::velocity,rtdt,rt,acc,veltdt,forces,vtemp,h,b,c,delv
double precision::kinetic

!!! test stuff
!double precision,dimension(3,natoms)::q
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
print*,'                   NOSE-HOOVER THERMOSTAT                    '    
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





!  do the dynamics

!notes for later. take tstep*tstep out of loop i.e. precompute it and 
! only compute acc-zeta*velocity once

a=-1.0d0
gfactor=3.0d0*natoms*boltzman*temp_input
zeta=0.
maxiter=6
Q=1.0d-1
constraint_values=constraint_values/au_to_ang

do i=1,isteps ! MD loop
    rtdt=rt+velocity*tstep+.5d0*(acc-zeta*velocity)*tstep*tstep ! compute new coordinates

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

if(dsqrt(rnorm) .lt. 1d-4)then
call space(1)
do ik=1,num_constraints
rx=rtdt(1,iconstraint_pairs(ik,1))-rtdt(1,iconstraint_pairs(ik,2))
ry=rtdt(2,iconstraint_pairs(ik,1))-rtdt(2,iconstraint_pairs(ik,2))
rz=rtdt(3,iconstraint_pairs(ik,1))-rtdt(3,iconstraint_pairs(ik,2))
rnum=rx**2 + ry**2 + rz**2
rnum=dsqrt(rnum)
!print*,ik,rnum*au_to_ang,constraint_values(ik)*au_to_ang
end do
rnorm2=rnorm
!print*,'error vec norm',dsqrt(rnorm)

exit
end if

end do
!stop


    veltdt=velocity+0.5d0*tstep*(acc-zeta*velocity)             ! velocity at half timestep 
    call get_forces(rtdt, forces, ener, virial)         !hartrees/bohr and energy in hartree
    forces=forces/hartree_per_eint              !convert forces to Eo/Bohr
                                                ! convert forces to accelerations in this loop 
    do ii=1,natoms
    do j=1,3
      acc(j,ii)=forces(j,ii)/atom_mass(type_atom(ii))
   end do
   end do
                                                ! accelerations are done 

! so veltdt is v'
! compute zeta'

   vx=0.
   vy=0.
   vz=0.
   total=0.
   do j=1,natoms ! can probably get this from kinetic energy from previous step
      vx=vx+atom_mass(type_atom(j))*velocity(1,j)**2
      vy=vy+atom_mass(type_atom(j))*velocity(2,j)**2
      vz=vz+atom_mass(type_atom(j))*velocity(3,j)**2
   end do
   total=vx+vy+vz
   zetap=zeta+0.5d0*tstep*(total-gfactor)/Q


! with v' and zeta', solve the equations

do k=1,maxiter

   h=veltdt+0.5d0*tstep*(acc-zeta*velocity)-velocity
   hn1=zetap+(total-gfactor)*0.5d0*tstep/Q - zeta


   do j=1,natoms
      do jj=1,3
         b(jj,j)=atom_mass(type_atom(j))*velocity(jj,j)
      end do
   end do
   b=b*tstep/Q
   c=-velocity*0.5d0*tstep
   d=-zeta*0.5d0*tstep-1.0d0
   term1=ddot(3*natoms,h,1,b,1)
   term2=ddot(3*natoms,b,1,c,1)
   delzeta=hn1*d-term1
   delzeta=delzeta/(-a*d +term2)
   delv=(-h-c*delzeta)
   delv=delv/d
   velocity=velocity+delv
   zeta=zeta+delzeta



if(k.lt.maxiter)then
   vx=0.
   vy=0.
   vz=0.
   total=0.
   do j=1,natoms ! can probably get this from kinetic energy from previous step
      vx=vx+atom_mass(type_atom(j))*velocity(1,j)**2
      vy=vy+atom_mass(type_atom(j))*velocity(2,j)**2
      vz=vz+atom_mass(type_atom(j))*velocity(3,j)**2
   end do
   total=vx+vy+vz

else
rnorm=ddot(3*natoms,h,1,h,1)
if(rnorm.gt.1d-6 .or. dabs(hn1).gt.1d-6)then
call space(1)
          print*,'        ********WARNING********'
          print*,'NOSE HOOVER ITERATES HAVE LOOSE CONVERGENCE'
          print*,'Current errors are',rnorm,' and',dabs(hn1)
          print*,'Try increasing parameter maxiter in md_nose.f90'
end if 
end if
end do ! end loop over newton raphson




   call  com_momentum(velocity,total_mass,rtdt)   ! remove translation/rotation

!   if(mod(i,iscale_counter).eq.0)scale=.true.   !scale to T (subroutine sets scale to false if scaling is done)
   scale=.false.
   call ke_and_t(velocity,kinetic,temperature,scale,temp_input,pbc)

   call space(2)                                ! print information
   print*,'                           Timestep # ',i
   print*,'           **************************************************'
   print*,'            Simulation Time   (fs):',float(i)*tstep*1000.0d0
   print*,'            Potential Energy (a.u):',ener
   print*,'            Kinetic Energy   (a.u):' ,kinetic*hartree_per_eint
   print*,'            Temperature        (K):',temperature
   print*,'            Constraint Vector Norm:',dsqrt(rnorm2)
!   velocity=veltdt                              ! save new velocities as old velocities
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





end subroutine md_constraints_nose



!   if(i.eq.1 .or. mod(i,100).eq.0)then
!   write(15,*)natoms
!   write(15,*)natoms
!   do kk=1,natoms
!   write(15,50)input_names(kk),rt(1,kk)*au_to_ang,rt(2,kk)*au_to_ang,rt(3,kk)*au_to_ang
!   end do
!   end if
!   50 format(A,f10.5,f10.5,f10.5)
