subroutine md
use atom_types_module
use control_module
use cell_module
use constants
use hugoniot_module
implicit double precision (a-h,o-z)
interface
    subroutine ke_and_t(velocity,kinetic,temperature,scale,scaled_t,pbcflag, kinetic_virial,temp_tol)
     double precision,dimension(:,:)::velocity
     double precision,intent(inout)::kinetic,temperature
     double precision,intent(in)::scaled_t,temp_tol
     logical,intent(inout)::scale,pbcflag
     double precision,intent(inout),dimension(3,3)::kinetic_virial
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
double precision,dimension(3,natoms)::velocity,rtdt,rt,acc,veltdt,forces,rtdt2,force_average,force_average2
double precision::kinetic,temperature_av
double precision,dimension(3,3)::stress_virial,kinetic_virial,pressure_virial,press_average_mat,press_average_mat2,pwrite


open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'                      #     # ######'
write(761,*)'                      ##   ## #     #'
write(761,*)'                      # # # # #     #'
write(761,*)'                      #  #  # #     #'
write(761,*)'                      #     # #     #'
write(761,*)'                      #     # #     #'
write(761,*)'                      #     # ######'
write(761,*)''
write(761,*)' #####  ####### #     #  #####  #######    #    #     # #######'
write(761,*)'#     # #     # ##    # #     #    #      # #   ##    #    #'
write(761,*)'#       #     # # #   # #          #     #   #  # #   #    #'
write(761,*)'#       #     # #  #  #  #####     #    #     # #  #  #    #'
write(761,*)'#       #     # #   # #       #    #    ####### #   # #    #'
write(761,*)'#     # #     # #    ## #     #    #    #     # #    ##    #'
write(761,*)' #####  ####### #     #  #####     #    #     # #     #    #'
write(761,*)''
write(761,*)'             ####### ####### #     # ######'
write(761,*)'                #    #       ##   ## #     #'
write(761,*)'                #    #       # # # # #     #'
write(761,*)'                #    #####   #  #  # ######'
write(761,*)'                #    #       #     # #         ###'
write(761,*)'                #    #       #     # #         ###'
write(761,*)'                #    ####### #     # #         ###'
write(761,*)'#     #'
write(761,*)'#     #  ######  #        ####    ####      #     #####   #   #'
write(761,*)'#     #  #       #       #    #  #    #     #       #      # #'
write(761,*)'#     #  #####   #       #    #  #          #       #       #'
write(761,*)' #   #   #       #       #    #  #          #       #       #'
write(761,*)'  # #    #       #       #    #  #    #     #       #       #'
write(761,*)'   #     ######  ######   ####    ####      #       #       #'
write(761,*)''
write(761,*)' #####'
write(761,*)'#     #   ####     ##    #          #    #    #   ####'
write(761,*)'#        #    #   #  #   #          #    ##   #  #    #'
write(761,*)' #####   #       #    #  #          #    # #  #  #'
write(761,*)'      #  #       ######  #          #    #  # #  #  ###'
write(761,*)'#     #  #    #  #    #  #          #    #   ##  #    #'
write(761,*)' #####    ####   #    #  ######     #    #    #   ####'

close(761)



















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
call ke_and_t(velocity,kinetic,temperature,scale,temp_input,pbc,kinetic_virial,temp_tol)
!                done with initial conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call space(5)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up for MD loop. 
rt=xyz ! copy input coordinates into MD working array
rtdt=zero !initialize some MD arrays 
acc=zero
veltdt=zero 
rtdt2=rt
rnorm=0.

if(restart)then
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'***********READING RESTART FILE************'
open(unit=15,file='RESTARTOLD')
read(15,*)junk
write(761,*)'Restarting from timestep ',junk
read(15,*)cell
read(15,*)rt
rtdt2=rt
xyz=rt
read(15,*)velocity
close(15)
close(761)
call ke_and_t(velocity,kinetic,temperature,scale,temp_input,pbc,kinetic_virial,temp_tol)
end if







call get_forces(rtdt2, forces, ener, stress_virial) !get initial forces and energy (hartrees/bohr and hartree)
forces=forces/hartree_per_eint !convert to internal units
! get initial accelerations from initial forces
do i=1,natoms
do j=1,3
acc(j,i)=forces(j,i)/atom_mass(type_atom(i))
rnorm=rnorm+forces(j,i)*forces(j,i)*hartree_per_eint**2
end do
end do
rnorm=dsqrt(rnorm / float (3*natoms))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!call cell_volume(cell,volume)
!pressure_virial=(2.*kinetic_virial*hartree_per_eint+stress_virial)/3.0d0/volume
!press=pressure_virial(1,1)+pressure_virial(2,2)+pressure_virial(3,3)

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'                 |----------------------------------|'
write(761,*)'                 |     Initial Cell Parameters      |'
write(761,*)'                 |----------------------------------|'
call cell_volume(cell,volume)
call space(1)
pressure_virial=(2.*kinetic_virial*hartree_per_eint+stress_virial)/volume
!press=pressure_virial(1,1)+pressure_virial(2,2)+pressure_virial(3,3)
call space(2)
write(761,*)'                 |----------------------------------|'
write(761,*)'                 |     Initial Stress Tensor (GPa)  |'
write(761,*)'                 |----------------------------------|'
ab=29421.0107637093d0
write(761,17)pressure_virial(1,1)*ab,pressure_virial(1,2)*ab,pressure_virial(1,3)*ab
write(761,17)pressure_virial(2,1)*ab,pressure_virial(2,2)*ab,pressure_virial(2,3)*ab
write(761,17)pressure_virial(3,1)*ab,pressure_virial(3,2)*ab,pressure_virial(3,3)*ab
call space(2)

17 format(f20.6,f20.6,f20.6)


if(strain_is_on)then
call strain_lattice
rtdt2=xyz ! strained coordinates from strain_lattice call
rt=xyz
call get_forces(rtdt2, forces, ener, stress_virial) !get initial forces and energy (hartrees/bohr and hartree)
pressure_virial=(2.*kinetic_virial*hartree_per_eint+stress_virial)/volume_strain
forces=forces/hartree_per_eint !convert to internal units
! get initial accelerations from initial forces
rnorm=0.
do i=1,natoms
do j=1,3
acc(j,i)=forces(j,i)/atom_mass(type_atom(i))
rnorm=rnorm+forces(j,i)*forces(j,i)*hartree_per_eint**2
end do
end do
rnorm=dsqrt(rnorm / float (3*natoms))
write(761,*)'                 |-------------------------------------|'
write(761,*)'                 |  Stress Tensor (GPa)- After Strain  |'
write(761,*)'                 |-------------------------------------|'
ab=29421.0107637093d0
write(761,17)pressure_virial(1,1)*ab,pressure_virial(1,2)*ab,pressure_virial(1,3)*ab
write(761,17)pressure_virial(2,1)*ab,pressure_virial(2,2)*ab,pressure_virial(2,3)*ab
write(761,17)pressure_virial(3,1)*ab,pressure_virial(3,2)*ab,pressure_virial(3,3)*ab
call space(2)
end if










!press=(2.*kinetic*hartree_per_eint+virial)/3.0d0/volume
!print*,'initial press is',press*auforce_to_bar/1.01325
write(761,*)'***************************************************************'
write(761,*)'***************************************************************'
write(761,*)'        VELOCITY-VERLET  MOLECULAR DYNAMICS ALGORITHM        '
write(761,*)'                    Velocity Scaling'
call space(2)
write(761,100)'              -> Number of timesteps:     ',isteps
write(761,110)'              -> Time step (ps):          ',tstep
write(761,110)'              -> Initial temperature (K): ',temperature
!write(761,110)'              -> Initial pressure (bar):  ',press*auforce_to_bar,press*auforce_to_bar/1.01325
write(761,110)'              -> Cutoff (Angstroms):      ',cutoff
write(761,110)'              -> Initial Gradient Norm      ',rnorm
write(761,*)'             -> PBC flag                 ',pbc
if(pbc)write(761,120)'              -> Constraining center of mass translation'
if(.not. pbc)write(761,120)'              -> Constraining center of mass translation and rotation'
write(761,*)'***************************************************************'
write(761,*)'***************************************************************'
100 format(A,I10)
110 format(A,F15.8)
120 format(A)
close(761)





!  do the dynamics
temperature_av=0.
press_average_mat=0.
paverage=0.
etotal_av=0.
conav=0.0d0
force_average=0.0d0


open(unit=23,file='Instantaneous',access='append')
write(23,*)'Time   Temperature Press Stress Tensor '
close(23)
open(unit=24,file='Averages',access='append')
write(24,*)'Time   Temperature Press Stress Tensor '
close(24)




do i=1,isteps ! MD loop

    rtdt=rt+velocity*tstep+.5d0*acc*tstep*tstep ! compute new coordinates
    veltdt=velocity+0.5d0*tstep*acc             ! velocity at half timestep 
    rtdt2=rtdt
    call get_forces(rtdt2, forces, ener, stress_virial) !forces hartrees/bohr. energy in hartree. virial hartree (sum r*f)
    force_average=force_average+forces

    forces=forces/hartree_per_eint              !convert forces to Eo/Bohr
                                                ! convert forces to accelerations in this loop 
    rnorm=0.
    do ii=1,natoms
    do j=1,3
      acc(j,ii)=forces(j,ii)/atom_mass(type_atom(ii))
      rnorm=rnorm+forces(j,ii)*forces(j,ii)*hartree_per_eint**2
   end do
   end do
   rnorm=dsqrt(rnorm / float(3*natoms))



                                                ! accelerations are done 
   veltdt=veltdt+.5d0 * tstep * acc             ! update half time step velocities to current velocities 
   call  com_momentum(veltdt,total_mass,rtdt)   ! remove translation/rotation

   if(mod(i,iscale_counter).eq.0)scale=.true.   !scale to T (subroutine sets scale to false if scaling is done)
   call ke_and_t(veltdt,kinetic,temperature,scale,temp_input,pbc,kinetic_virial, temp_tol)
   kinetic=kinetic*hartree_per_eint !convert ke to hartree
 etotal_av=etotal_av+kinetic+ener
conav=conav+ener
   pressure_virial=(2.*kinetic_virial*hartree_per_eint+stress_virial)/volume
   pressure=(pressure_virial(1,1)+pressure_virial(2,2)+pressure_virial(3,3))/3.0d0

   temperature_av=temperature_av+temperature
   press_average_mat=pressure_virial+press_average_mat
   paverage=pressure+paverage

  force_average2=force_average/float(i)
  rnormpr=ddot(3*natoms,force_average2,1,force_average2,1)
  rnormpr=dsqrt(rnormpr / float (3*natoms))

   call space(2)                                ! print information
   open(unit=761,file='OUTPUT.DTPOLY',access='append')
   write(761,*)'                           Timestep # ',i
   write(761,*)'           **************************************************'
   write(761,*)'            Simulation Time   (fs):',float(i)*tstep*1000.0d0
   write(761,*)'            Potential Energy (a.u):',ener
   write(761,*)'            Kinetic Energy   (a.u):' ,kinetic
   write(761,*)'            Avg. Etotal      (a.u):',etotal_av/float(i)
   write(761,*)'            Temperature(unscaled)(K):',temp_unscaled
   write(761,*)'            Temperature(scaled)  (K):',temperature
   write(761,*)'            Avg. Temp.         (K):',temperature_av/float(i)
   write(761,*)'            Pressure         (GPa):',pressure*auforce_to_bar*1d-4
   write(761,*)'            Pressure-avg     (GPa):',(paverage/float(i))*auforce_to_bar*1d-4
   write(761,*)'            RMS Norm         :',rnorm
   write(761,*)'            Avg Force RMS Norm    :',rnormpr
hugoniot_pressure= (press_average_mat(1,1)+press_average_mat(2,2)+press_average_mat(3,3))*29421.0107637093d0/(3.0d0*float(i))
hugoniot_energy=etotal_av/float(i)
specific_energy=627.51d0*hugoniot_energy/(total_masskg*avogadro) ! kcal/kg
density=1.0D30 * total_masskg/(au_to_ang**3 * hugoniot_volume) ! in kg/m^3
specific_volume=1.0d0/density ! m^3/kg
hg=(specific_energy-specific_energy_ref) + 0.5d0 * (hugoniot_pressure_ref + hugoniot_pressure) * (specific_volume - specific_volume_ref) * rfactor
    write(761,*)'            Hugoniot Function     :',hg


hg_kelvin=avogadro*hg*total_masskg/(627.51d0*hartree_per_eint * boltzman* (3*natoms-3.0d0) )
   
  write(761,*)'            Hugoniot Function(K)     :',hg_kelvin

rfactor4=(1.0d0/29421.0107637093)*(627.51/avogadro)*(1.0d0/(0.52917D-10))**3*4184.0d0
shock_vel=dsqrt ( rfactor4* specific_volume_ref * (hugoniot_pressure - hugoniot_pressure_ref) / ( 1.0d0 - hugoniot_volume/hugoniot_volume_ref) )

write(761,*)'Shock Velocity (km/s)',shock_vel/1000.0d0
part_vel=dsqrt ( rfactor4* specific_volume_ref * (hugoniot_pressure - hugoniot_pressure_ref) * ( 1.0d0 - hugoniot_volume/hugoniot_volume_ref) )
write(761,*)'Particle Velocity (km/s)',part_vel/1000.0d0










   call space(2)
   write(761,*)'                 |----------------------------------|'
   write(761,*)'                 |         Stress Tensor (GPa)      |'
   write(761,*)'                 |----------------------------------|'
   write(761,17)pressure_virial(1,1)*29421.0107637093d0,pressure_virial(1,2)*29421.0107637093d0,pressure_virial(1,3)*29421.0107637093d0
   write(761,17)pressure_virial(2,1)*29421.0107637093d0,pressure_virial(2,2)*29421.0107637093d0,pressure_virial(2,3)*29421.0107637093d0
   write(761,17)pressure_virial(3,1)*29421.0107637093d0,pressure_virial(3,2)*29421.0107637093d0,pressure_virial(3,3)*29421.0107637093d0
   call space(2)

   press_average_mat2=press_average_mat*29421.0107637093d0/float(i)
   call space(2)
   write(761,*)'                 |----------------------------------|'
   write(761,*)'                 |     Average Stress Tensor (GPa)  |'
   write(761,*)'                 |----------------------------------|'
   write(761,17)press_average_mat2(1,1),press_average_mat2(1,2),press_average_mat2(1,3)
   write(761,17)press_average_mat2(2,1),press_average_mat2(2,2),press_average_mat2(2,3)
   write(761,17)press_average_mat2(3,1),press_average_mat2(3,2),press_average_mat2(3,3)
   call space(2)
   close(761)


pwrite=pressure_virial*29421.0107637093d0
open(unit=23,file='Instantaneous',access='append')
write(23,24)i,temperature,pressure*auforce_to_bar*1.0d-4,pwrite(1,1),pwrite(2,2),pwrite(3,3),pwrite(2,3),pwrite(1,3),pwrite(1,2),ener,kinetic+ener
close(23)
pwrite=press_average_mat2
open(unit=23,file='Averages',access='append')
write(23,24)i,temperature_av/float(i),(paverage/float(i))*auforce_to_bar*1.0d-4,pwrite(1,1),pwrite(2,2),pwrite(3,3),pwrite(2,3),pwrite(1,3),pwrite(1,2),conav/float(i),etotal_av/float(i)
close(23)
24 format(i8,f12.3,f12.3,f12.3,f12.3,f12.3,f12.3,f12.3,f12.3,f18.5,f18.5)

open(unit=23,file='FORCES.AV',access='append')
write(23,34)i,rnorm,rnormpr
close(23)
34 format(i8,f15.8,f15.8)





if(i.eq.500 .and. abs(hg).gt.5.0d0)return




   velocity=veltdt                              ! save new velocities as old velocities
   rt=rtdt                                      ! save new coordinates as old coordinates

   if(i.eq.1 .or. mod(i,iwrite).eq.0)then
 open(unit=15,file='OUTPUT.TRAJECTORY',access='append')
   write(15,*)natoms
   write(15,*)'timestep',i 
   do kk=1,natoms
   write(15,50)input_names(kk),rt(1,kk)*au_to_ang,rt(2,kk)*au_to_ang,rt(3,kk)*au_to_ang
   end do
 CLOSE(15)

   open(unit=15,file='RESTART')
   write(15,*)i
   write(15,*)cell
   write(15,*)rt
   write(15,*)velocity
   close(15)


open(unit=315,file='OUTPUT.TRAJECTORY.CELL',access='append')
  write(315,*)'timestep',iii
  write(315,*)cell*au_to_ang
close(315)








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




return

end subroutine md



