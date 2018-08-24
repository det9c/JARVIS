subroutine md_berendsen
use atom_types_module
use control_module
use cell_module
use constants
implicit double precision (a-h,o-z)
interface
   subroutine ke_and_t(velocity,kinetic,temperature,scale,scaled_t,pbcflag,kinetic_virial)
     double precision,dimension(:,:)::velocity
     double precision,intent(inout)::kinetic,temperature
     double precision,intent(in)::scaled_t
     double precision,intent(inout),dimension(3,3)::kinetic_virial
     logical,intent(inout)::scale,pbcflag
   end subroutine ke_and_t

   subroutine com_momentum(velocity,total_mass,rxyz)
     double precision,dimension(:,:),intent(inout)::velocity,rxyz
     double precision,intent(in)::total_mass
   end subroutine com_momentum

   subroutine cell_volume(cell_in,volume)
   double precision,dimension(3,3),intent(in)::cell_in
   double precision,intent(inout)::volume
   end subroutine cell_volume





end interface
double precision,dimension(3,natoms)::velocity,rtdt,rt,acc,veltdt,forces,rtdt2
double precision::kinetic
double precision,dimension(3,3)::stress_virial,kinetic_virial,pressure_virial


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
call ke_and_t(velocity,kinetic,temperature,scale,temp_input,pbc,kinetic_virial)

!                done with initial conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



open(unit=15,file='md_history.xyz') ! file to hold trajectory 



call space(5)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up for MD loop. 

rtdt=xyz !initialize some MD arrays 
acc=zero
veltdt=zero 
rtdt2=rtdt
    call get_forces(rtdt2, forces, ener, stress_virial) !forces hartrees/bohr. energy in hartree. virial hartree (sum r*f)
    forces=forces/hartree_per_eint              !convert forces to Eo/Bohr
    do ii=1,natoms
    do j=1,3
      acc(j,ii)=forces(j,ii)/atom_mass(type_atom(ii))
   end do
   end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

cell_flex=cell
call cell_volume(cell_flex,volume)
pressure_virial=(2.*kinetic_virial*hartree_per_eint+stress_virial)/3.0d0/volume
press=pressure_virial(1,1)+pressure_virial(2,2)+pressure_virial(3,3)





kinetic=kinetic*hartree_per_eint







print*,'***************************************************************'
print*,'***************************************************************'
print*,'          LEAPFROG  MOLECULAR DYNAMICS ALGORITHM               '
print*,'        BERENDSEN BAROSTAT/THERMOSTAT (UNCOUPLED)              '
call space(2)
write(*,100)'              -> Number of timesteps:     ',isteps
write(*,110)'              -> Time step (ps):          ',tstep
write(*,110)'              -> Initial temperature (K): ',temperature
write(*,*)'              -> Initial pressure (bar):  ',press*auforce_to_bar,press*auforce_to_bar/1.01325
write(*,110)'              -> Target pressure (bar):  ',pressure_input
write(*,110)'              -> Cutoff (Angstroms):      ',cutoff
print*,'             -> PBC flag                 ',pbc
if(pbc)write(*,120)'              -> Constraining center of mass translation'
if(.not. pbc)write(*,120)'              -> Constraining center of mass translation and rotation'
print*,'***************************************************************'
print*,'***************************************************************'
100 format(A,I10)
110 format(A,F15.4)
120 format(A)






!  do the dynamics


pressure_input=pressure_input/auforce_to_bar
paverage=0.
paverage2=0.
beta=4.6d-5*auforce_to_bar
av_cell_lengths=0.
av_cell_angles=0.



!print*,'tau is',tau,beta,pressure_input,press


do i=1,isteps ! MD loop

    rtdt2=rtdt
    call get_forces(rtdt2, forces, ener, stress_virial) !forces hartrees/bohr. energy in hartree. virial hartree (sum r*f)
    forces=forces/hartree_per_eint              !convert forces to Eo/Bohr
    do ii=1,natoms
    do j=1,3
      acc(j,ii)=forces(j,ii)/atom_mass(type_atom(ii))
   end do
   end do


! pressure=(2.0d0*kinetic+virial)/3.0d0/volume ! old kinetic energy

pressure_virial=(2.*kinetic_virial*hartree_per_eint+stress_virial)/3.0d0/volume
pressure=pressure_virial(1,1)+pressure_virial(2,2)+pressure_virial(3,3)


factor=1.0 - (beta*tstep/tau)*(pressure_input - pressure)
if(factor.lt.0.0)then
write(*,*)'scale factor should be > 0. Change barostat time constant Tau'
stop
end if

 dmu=(factor )**(1./3.) 

                                                ! accelerations are done 
   veltdt=velocity+  tstep * acc             ! update half time step velocities to current velocities 
   call  com_momentum(veltdt,total_mass,rtdt)   ! remove translation/rotation
   if(mod(i,iscale_counter).eq.0)scale=.true.   !scale to T (subroutine sets scale to false if scaling is done)
   call ke_and_t(veltdt,kinetic,temperature,scale,temp_input,pbc,kinetic_virial)
   kinetic=kinetic*hartree_per_eint !convert ke to hartree

   rt=rtdt
   velocity=veltdt 
    rtdt=rt+veltdt*tstep

! scale coordinates


  rtdt=rtdt*dmu
  cell=cell*dmu
! call matprt(cell,3,3,3,3)
! call cell_volume(cell,volume)

paverage=pressure+paverage



   call space(2)                                ! print information
   print*,'                           Timestep # ',i
   print*,'           **************************************************'
   print*,'            Simulation Time   (fs):',float(i)*tstep*1000.0d0
   print*,'            Potential Energy (a.u):',ener
   print*,'            Kinetic Energy   (a.u):' ,kinetic
   print*,'            Temperature        (K):',temperature
   print*,'            Pressure         (bar):',pressure*auforce_to_bar
   print*,'            Pressure-avg     (bar):',(paverage/float(i))*auforce_to_bar  

if(i .ge. istat_start)then
cell_stats=.true.
paverage2=pressure+paverage2
end if

call cell_volume(cell,volume)









   if(i.eq.1 .or. mod(i,iwrite).eq.0)then
   write(15,*)natoms
   write(15,*)'timestep',i 
   do kk=1,natoms
   write(15,50)input_names(kk),rt(1,kk)*au_to_ang,rt(2,kk)*au_to_ang,rt(3,kk)*au_to_ang
   end do
   end if
   50 format(A,f10.5,f10.5,f10.5)


end do ! end loop over time steps

mm=isteps-istat_start+1
write(*,*)'FINAL AVERAGE PRESSURE IS',auforce_to_bar*paverage2/float(mm),' OVER',mm,' STEPS'
av_cell_lengths=av_cell_lengths*0.5291772108d0/float(mm)
av_cell_angles=av_cell_angles*  180./pi/float(mm)
av_cell_volume=av_cell_volume/float(mm)
write(*,*)'Average Vector Lengths    ',av_cell_lengths(1),'        ',av_cell_lengths(2),'     ',av_cell_lengths(3)
write(*,*)'Vector Angles    ',av_cell_angles(1),'        ',av_cell_angles(2),'     ',av_cell_angles(3)
write(*,*)'Average Volume    ',av_cell_volume*au_to_ang**3








call space(2)
print*,'               ******************************************'
print*,'               ******************************************'
print*,'               ******************************************'
print*,'               **********MD run has completed************'
print*,'               ******************************************'
print*,'               ******************************************'
print*,'               ******************************************'




return

end subroutine md_berendsen



