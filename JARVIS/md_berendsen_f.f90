subroutine md_berendsen_f
use atom_types_module
use control_module
use cell_module
use constants
use virial_mod
implicit double precision (a-h,o-z)
interface
   subroutine ke_and_t(velocity,kinetic,temperature,scale,scaled_t,pbcflag,kinetic_virial,temp_tol)
     double precision,dimension(:,:)::velocity
     double precision,intent(inout)::kinetic,temperature
     double precision,intent(in)::scaled_t,temp_tol
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
double precision,dimension(3,3)::stress_virial,kinetic_virial,pressure_virial,factor &
,pressure_in_mat,cell2,virial_av,press_average_mat,press_average_mat2,scaled_pressure_virial &
,pwrite

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'                  #     # ######'
write(761,*)'                  ##   ## #     #'
write(761,*)'                  # # # # #     #'
write(761,*)'                  #  #  # #     #'
write(761,*)'                  #     # #     #'
write(761,*)'                  #     # #     #'
write(761,*)'                  #     # ######'
write(761,*)''
write(761,*)' #####  ####### #     #  #####  #######    #    #     # #######'
write(761,*)'#     # #     # ##    # #     #    #      # #   ##    #    #'
write(761,*)'#       #     # # #   # #          #     #   #  # #   #    #'
write(761,*)'#       #     # #  #  #  #####     #    #     # #  #  #    #'
write(761,*)'#       #     # #   # #       #    #    ####### #   # #    #'
write(761,*)'#     # #     # #    ## #     #    #    #     # #    ##    #'
write(761,*)' #####  ####### #     #  #####     #    #     # #     #    #'
write(761,*)''
write(761,*)'        #####  ####### ######  #######  #####   #####'
write(761,*)'        #     #    #    #     # #       #     # #     #'
write(761,*)'        #          #    #     # #       #       #'
write(761,*)'         #####     #    ######  #####    #####   #####'
write(761,*)'              #    #    #   #   #             #       #'
write(761,*)'              #    #    #    #  #  #     #    # #     #'
write(761,*)'         #####     #    #     # #######  #####   #####'







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



!open(unit=15,file='md_history.xyz') ! file to hold trajectory 



call space(5)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up for MD loop. 

rtdt=xyz !initialize some MD arrays 
acc=zero
veltdt=zero 
rtdt2=rtdt
rnorm=0.

if(restart)then
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'***********READING RESTART FILE************'
open(unit=15,file='RESTARTOLD')
read(15,*)junk
write(761,*)'Restarting from timestep ',junk
read(15,*)cell
read(15,*)rtdt
rtdt2=rtdt
xyz=rtdt
read(15,*)velocity
close(761)
scale=.false.
call ke_and_t(velocity,kinetic,temperature,scale,temp_input,pbc,kinetic_virial,temp_tol)
end if

    call get_forces(rtdt2, forces, ener, stress_virial) !forces hartrees/bohr. energy in hartree. virial hartree (sum r*f)
    forces=forces/hartree_per_eint              !convert forces to Eo/Bohr
    do ii=1,natoms
    do j=1,3
      acc(j,ii)=forces(j,ii)/atom_mass(type_atom(ii))
rnorm=rnorm+forces(j,ii)*forces(j,ii)*hartree_per_eint**2 
   end do
   end do
rnorm=dsqrt(rnorm)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pressure_in_mat=stress_in
cell_flex=cell
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'                 |----------------------------------|'
write(761,*)'                 |     Initial Cell Parameters      |'
write(761,*)'                 |----------------------------------|'
call cell_volume(cell_flex,volume)
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
write(761,*)'                 |----------------------------------|'
write(761,*)'                 |     Target  Stress Tensor (GPa)  |'
write(761,*)'                 |----------------------------------|'
write(761,17)pressure_in_mat(1,1),pressure_in_mat(1,2),pressure_in_mat(1,3)
write(761,17)pressure_in_mat(2,1),pressure_in_mat(2,2),pressure_in_mat(2,3)
write(761,17)pressure_in_mat(3,1),pressure_in_mat(3,2),pressure_in_mat(3,3)
call space(2)
close(761)
17 format(f20.6,f20.6,f20.6)




if(strain_is_on)then
open(unit=761,file='OUTPUT.DTPOLY',access='append')
call strain_lattice
rtdt2=xyz ! strained coordinates from strain_lattice call
rtdt=xyz
call get_forces(rtdt2, forces, ener, stress_virial) !get initial forces and energy (hartrees/bohr and hartree)
pressure_virial=(2.*kinetic_virial*hartree_per_eint+stress_virial)/volume_strain
forces=forces/hartree_per_eint !convert to internal units
! get initial accelerations from initial forces
do i=1,natoms
do j=1,3
acc(j,i)=forces(j,i)/atom_mass(type_atom(i))
rnorm=rnorm+forces(j,i)*forces(j,i)*hartree_per_eint**2
end do
end do
rnorm=dsqrt(rnorm)
write(761,*)'                 |-------------------------------------|'
write(761,*)'                 |  Stress Tensor (GPa)- After Strain  |'
write(761,*)'                 |-------------------------------------|'
ab=29421.0107637093d0
write(761,17)pressure_virial(1,1)*ab,pressure_virial(1,2)*ab,pressure_virial(1,3)*ab
write(761,17)pressure_virial(2,1)*ab,pressure_virial(2,2)*ab,pressure_virial(2,3)*ab
write(761,17)pressure_virial(3,1)*ab,pressure_virial(3,2)*ab,pressure_virial(3,3)*ab
call space(2)
close(761)
end if



kinetic=kinetic*hartree_per_eint






open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'***************************************************************'
write(761,*)'***************************************************************'
write(761,*)'          LEAPFROG  MOLECULAR DYNAMICS ALGORITHM               '
write(761,*)'        BERENDSEN BAROSTAT/THERMOSTAT (UNCOUPLED)              '
write(761,*)'                 Fully Flexible Cell                           '
call space(2)
write(761,100)'              -> Number of timesteps:     ',isteps
write(761,110)'              -> Time step (ps):          ',tstep
write(761,110)'              -> Initial temperature (K): ',temperature
!write(761,110)'              -> Initial pressure (bar):  ',press*auforce_to_bar,press*auforce_to_bar/1.01325
!write(761,110)'              -> Target pressure (bar):  ',pressure_input
write(761,110)'              -> Cutoff (Angstroms):      ',cutoff
write(761,110)'              -> Initial Gradient Norm      ',rnorm
write(761,*)'             -> PBC flag                 ',pbc
if(pbc)write(761,120)'              -> Constraining center of mass translation'
if(.not. pbc)write(761,120)'              -> Constraining center of mass translation and rotation'
write(761,*)'***************************************************************'
write(761,*)'***************************************************************'
100 format(A,I10)
110 format(A,F15.4)
120 format(A)
close(761)





!  do the dynamics


!pressure_input=pressure_input/auforce_to_bar
!pressure_in_mat=unit_matrix*pressure_input
pressure_in_mat=stress_in/29421.0107637093d0


temperature_av=0.
paverage=0.
paverage2=0.
press_average_mat=0.
beta=4.6d-5*auforce_to_bar
av_cell_lengths=0.
av_cell_angles=0.
av_cell_volume=0.
virial_av=0.


!print*,'tau is',tau,beta,pressure_input,press
call cpusec(start_time)

open(unit=23,file='Instantaneous',access='append')
write(23,*)'Time   Temperature Press Stress Tensor '
close(23)
open(unit=24,file='Averages',access='append')
write(24,*)'Time   Temperature Press Stress Tensor '
close(24)

do i=1,isteps ! MD loop




    rtdt2=rtdt
    call get_forces(rtdt2, forces, ener, stress_virial) !forces hartrees/bohr. energy in hartree. virial hartree (sum r*f)
    forces=forces/hartree_per_eint              !convert forces to Eo/Bohr
    rnorm=0.
    do ii=1,natoms
    do j=1,3
      acc(j,ii)=forces(j,ii)/atom_mass(type_atom(ii))
      rnorm=rnorm+forces(j,ii)*forces(j,ii)*hartree_per_eint**2
   end do
   end do
   rnorm=dsqrt(rnorm)


pressure_virial=(2.*kinetic_virial*hartree_per_eint+stress_virial)/volume
pressure=( pressure_virial(1,1)+pressure_virial(2,2)+pressure_virial(3,3) )/3.0d0


scaled_pressure_virial=pressure_virial*scale_factor
!factor=unit_matrix - (beta*tstep/tau)*(pressure_in_mat - pressure_virial)
factor=unit_matrix - (beta*tstep/tau)*(pressure_in_mat - scaled_pressure_virial)
                                ! accelerations are done 
   veltdt=velocity+  tstep * acc             ! update half time step velocities to current velocities 
   call  com_momentum(veltdt,total_mass,rtdt)   ! remove translation/rotation
   if(mod(i,iscale_counter).eq.0)scale=.true.   !scale to T (subroutine sets scale to false if scaling is done)
   call ke_and_t(veltdt,kinetic,temperature,scale,temp_input,pbc,kinetic_virial,temp_tol)
   kinetic=kinetic*hartree_per_eint !convert ke to hartree


   velocity=veltdt
    rt=rtdt+veltdt*tstep



! scale coordinates - old spot
rtdt=matmul(factor,rt)
cell2=matmul(factor,cell)
cell=cell2




!  rtdt=rtdt*dmu
!  cell=cell*dmu



! call matprt(cell,3,3,3,3)
! call cell_volume(cell,volume)
temperature_av=temperature_av+temperature
paverage=pressure+paverage
press_average_mat=pressure_virial+press_average_mat

open(unit=761,file='OUTPUT.DTPOLY',access='append')
   call space(2)                                ! print information
   write(761,*)'                           Timestep # ',i
   write(761,*)'                         ****************'
   write(761,*)'            Simulation Time   (fs):',float(i)*tstep*1000.0d0
   write(761,*)'            Potential Energy (a.u):',ener
   write(761,*)'            Kinetic Energy   (a.u):' ,kinetic
   write(761,*)'            Temperature        (K):',temperature
   write(761,*)'            Avg. Temp.         (K):',temperature_av/float(i)
   write(761,*)'            Pressure         (GPa):',pressure*auforce_to_bar*1d-4
   write(761,*)'            Pressure-avg     (GPa):',(paverage/float(i))*auforce_to_bar*1d-4  
   write(761,*)'            Gradient Norm         :',rnorm


call space(2)
write(761,*)'                |----------------------------------|'
write(761,*)'                |         Stress Tensor (GPa)      |'
write(761,*)'                |----------------------------------|'
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

pwrite=pressure_virial*29421.0107637093d0
open(unit=23,file='Instantaneous',access='append')
write(23,24)i,temperature,pressure*auforce_to_bar*1.0d-4,pwrite(1,1),pwrite(2,2),pwrite(3,3),pwrite(2,3),pwrite(1,3),pwrite(1,2)
close(23)
pwrite=press_average_mat2
open(unit=23,file='Averages',access='append')
write(23,24)i,temperature_av/float(i),paverage/float(i)*auforce_to_bar*1.0d-4,pwrite(1,1),pwrite(2,2),pwrite(3,3),pwrite(2,3),pwrite(1,3),pwrite(1,2)
close(23)




24 format(i8,f12.3,f12.3,f12.3,f12.3,f12.3,f12.3,f12.3,f12.3)

if(i .ge. istat_start)then
cell_stats=.true.
paverage2=pressure+paverage2
virial_av=virial_av+pressure_virial
end if

call cell_volume(cell,volume)
call space(2)
write(761,46)'Average Vector Lengths ',av_cell_lengths(1)*0.5291772108d0/float(i),av_cell_lengths(2)*0.5291772108d0/float(i),av_cell_lengths(3)*0.5291772108d0/float(i)
write(761,46)'Average Vector Angles',av_cell_angles(3)*180./pi/float(i),av_cell_angles(2)*180./pi/float(i),av_cell_angles(1)*180./pi/float(i)
write(761,47)'Average Volume    ',av_cell_volume*au_to_ang**3/float(i)
write(761,*)'*******************************************************************************'
46 format(A30,f20.10,f20.10,f20.10)
47 format(A30,f20.10)
close(761)







if(i.eq.1 .or. mod(i,iwrite).eq.0)then
!   open(unit=215,file='OUTPUT.TRAJECTORY',access='append')
open(unit=215,file='OUTPUT.TRAJECTORY',access='append')
open(unit=315,file='OUTPUT.TRAJECTORY.CELL',access='append')
   write(215,*)natoms
   write(215,*)'timestep',i 
   do kk=1,natoms
   write(215,50)input_names(kk),rtdt(1,kk)*au_to_ang,rtdt(2,kk)*au_to_ang,rtdt(3,kk)*au_to_ang
   end do
close(215)
  write(315,*)'timestep',i
  write(315,*)cell*au_to_ang
close(315)



! dump restart file
   open(unit=15,file='RESTART')
   write(15,*)i
   write(15,*)cell
   write(15,*)rtdt
   write(15,*)velocity
   close(15)


end if
   50 format(A,f20.5,f20.5,f20.5)




end do ! end loop over time steps

mm=isteps-istat_start+1
write(*,*)'FINAL AVERAGE PRESSURE IS',auforce_to_bar*paverage2/float(mm),' OVER',mm,' STEPS'
av_cell_lengths=av_cell_lengths*0.5291772108d0/float(mm)
av_cell_angles=av_cell_angles*  180./pi/float(mm)
av_cell_volume=av_cell_volume/float(mm)
write(*,*)'Average Vector Lengths    ',av_cell_lengths(1),'        ',av_cell_lengths(2),'     ',av_cell_lengths(3)
write(*,*)'Vector Angles    ',av_cell_angles(1),'        ',av_cell_angles(2),'     ',av_cell_angles(3)
write(*,*)'Average Volume    ',av_cell_volume*au_to_ang**3

write(*,*)'Average Theta    ',theta_av/float(mm)


write(*,*)'average stress tensor'
virial_av=auforce_to_bar*virial_av/float(mm)
call matprt(virial_av,3,3,3,3)



open(unit=12,file='FINALXYZ')
   do kk=1,natoms
   write(12,50)input_names(kk),rtdt(1,kk)*au_to_ang,rtdt(2,kk)*au_to_ang,rtdt(3,kk)*au_to_ang
   end do
close(12)

open(unit=12,file='cell')
write(12,*)av_cell_lengths(1)
write(12,*)av_cell_angles(1)
write(12,*)av_cell_volume*au_to_ang**3
write(12,*)theta_av/float(mm)
write(12,*)av_cell_angles(1)
close(12)










call space(2)
print*,'               ******************************************'
print*,'               ******************************************'
print*,'               ******************************************'
print*,'               **********MD run has completed************'
print*,'               ******************************************'
print*,'               ******************************************'
print*,'               ******************************************'




return

end subroutine md_berendsen_f



