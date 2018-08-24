subroutine md_nose_f
use atom_types_module
use control_module
use constants
use cell_module
use virial_mod
use hugoniot_module
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
double precision,dimension(3,natoms)::velocity,rtdt,rt,acc,veltdt,forces,vtemp,h,b,c,delv,rtdt2,force_average,force_average2 
double precision::kinetic
double precision,dimension(3,3)::stress_virial,kinetic_virial,pressure_virial,factor &
,pressure_in_mat,cell2,virial_av,press_average_mat,press_average_mat2,scaled_pressure_virial &
,pwrite,press_average_mat_unsymm,pressure_notsymm


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
write(761,*)'             #     # #######  #####  #######'
write(761,*)'             ##    # #     # #     # #'
write(761,*)'             # #   # #     # #       #'
write(761,*)'             #  #  # #     #  #####  #####'
write(761,*)'             #   # # #     #       # #'
write(761,*)'             #    ## #     # #     # #'
write(761,*)'             #     # #######  #####  #######'
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




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up for MD loop. 
rt=xyz ! copy input coordinates into MD working array
rtdt=zero !initialize some MD arrays 
acc=zero
veltdt=zero 


if(restart)then
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'***********READING RESTART FILE************'
open(unit=15,file='RESTARTOLD')
read(15,*)junk
write(761,*)'Restarting from timestep ',junk
read(15,*)cell
read(15,*)rt
xyz=rt
!rtdt2=rtdt
read(15,*)velocity
close(15)
close(761)
scale=.false.
call ke_and_t(velocity,kinetic,temperature,scale,temp_input,pbc,kinetic_virial,temp_tol)
end if


if(restart_cell_opt)then
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'***********READING BEST_POINT FILE************'
open(unit=15,file='best_point')
read(15,*)svalue,svalue2
write(761,*)'gradient norm of restart point',svalue
read(15,*)cell
read(15,*)rt
xyz=rt
close(15)
close(761)
end if






rtdt2=rt
call get_forces(rtdt2, forces, ener, stress_virial) !get initial forces and energy (hartrees/bohr and hartree)
forces=forces/hartree_per_eint !convert to internal units
! get initial accelerations from initial forces
do i=1,natoms
do j=1,3
acc(j,i)=forces(j,i)/atom_mass(type_atom(i))
end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



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
rt=xyz ! strained coordinates from strain_lattice call
rtdt=zero
rtdt2=rt
write(761,*)'                 |-------------------------------------|'
write(761,*)'                 |  Stress Tensor Constraint Matrix    |'
write(761,*)'                 |-------------------------------------|'
write(761,17)scale_factor(1,1) ,scale_factor(1,2) ,scale_factor(1,3)
write(761,17)scale_factor(2,1) ,scale_factor(2,2) ,scale_factor(2,3)
write(761,17)scale_factor(3,1) ,scale_factor(3,2) ,scale_factor(3,3)
call space(2)


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





close(761)
end if

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'***************************************************************'
write(761,*)'***************************************************************'
write(761,*)'          VELOCITY-VERLET  MOLECULAR DYNAMICS ALGORITHM        '
write(761,*)'        BERENDSEN BAROSTAT/HOOVER THERMOSTAT (UNCOUPLED)       '
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

!notes for later. take tstep*tstep out of loop i.e. precompute it and 
! only compute acc-zeta*velocity once

a=-1.0d0
gfactor=3.0d0*natoms*boltzman*temp_input
zeta=0.
maxiter=6
Q=1.0d1
etotal_av=0.
temperature_av=0.
press_average_mat=0.
paverage=0.
beta=4.6d-5*auforce_to_bar
!explanation of beta for umpteenth time b/c i keep having to "refigure" it out
!the equation from the dlpoly manual is [beta*tstep/tau] *(Pext-P)
!tau, in my input file has units of time (ps). beta is set arbitrarily to the isothermal  
!compressibility of water which is 4.6d-10 INVERSE PASCALS. This converts to 4.6d-5 inverse bar.
! In the code below, where beta is used, I have
! (beta*tstep/tau)*(pressure_in_mat - scaled_pressure_virial)
! The tstep and tau units cancel and beta is in inverse bar.  So the auforce_to_bar that I include in the
! expression for beta above is actually there to convert the (pressure_in_mat - scaled_pressure_virial)
! from atomic units to bar and everything is correct as all units are now bar.  


pressure_in_mat=stress_in/29421.0107637093d0
av_cell_lengths=0.
av_cell_angles=0.
av_cell_volume=0.
virial_av=0.
hg_av=0.0d0
conav=0.
press_average_mat_unsymm=0.
force_average=0.0d0


!open(unit=23,file='Instantaneous',access='append')
!write(23,*)'Time   Temperature Press Stress Tensor '
!close(23)
!open(unit=24,file='Averages',access='append')
!write(24,*)'Time   Temperature Press Stress Tensor '
!close(24)






do iii=1,isteps ! MD loop

    rtdt=rt+velocity*tstep+.5d0*(acc-zeta*velocity)*tstep*tstep ! compute new coordinates
    veltdt=velocity+0.5d0*tstep*(acc-zeta*velocity)             ! velocity at half timestep 
    rtdt2=rtdt
    call get_forces(rtdt2, forces, ener, stress_virial)         !hartrees/bohr and energy in hartree
   force_average=force_average+forces
    rnormpr=ddot(3*natoms,forces,1,forces,1)
    rnormpr=dsqrt(rnormpr / float (3*natoms) )

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
open(unit=761,file='OUTPUT.DTPOLY',access='append')
          write(761,*)'        ********WARNING********'
          write(761,*)'NOSE HOOVER ITERATES HAVE LOOSE CONVERGENCE'
          write(761,*)'Current errors are',rnorm,' and',dabs(hn1)
          write(761,*)'Try increasing parameter maxiter in md_nose_f.f90'
close(761)
end if 
end if
end do ! end loop over newton raphson




   call  com_momentum(velocity,total_mass,rtdt)   ! remove translation/rotation

!   if(mod(i,iscale_counter).eq.0)scale=.true.   !scale to T (subroutine sets scale to false if scaling is done)
!   scale=.false.
!   call ke_and_t(velocity,kinetic,temperature,scale,temp_input,pbc)
   call ke_and_t(velocity,kinetic,temperature,scale,temp_input,pbc,kinetic_virial,temp_tol)

  kinetic=kinetic*hartree_per_eint !convert ke to hartree
  pressure_virial=(2.*kinetic_virial*hartree_per_eint+stress_virial)/volume
  pressure=( pressure_virial(1,1)+pressure_virial(2,2)+pressure_virial(3,3) )/3.0d0
if(dft)pressure_notsymm=(2.*kinetic_virial*hartree_per_eint+virial_unsymm)/volume

  
 temperature_av=temperature_av+temperature
 etotal_av=etotal_av+kinetic+ener
 conav=conav+ener

 paverage=pressure+paverage
press_average_mat=pressure_virial+press_average_mat
if(dft)press_average_mat_unsymm=press_average_mat_unsymm+pressure_notsymm

scaled_pressure_virial=pressure_virial*scale_factor
!factor=unit_matrix - (beta*tstep/tau)*(pressure_in_mat - pressure_virial)
factor=unit_matrix - (beta*tstep/tau)*(pressure_in_mat - scaled_pressure_virial)






  force_average2=force_average/float(iii)
  rnormpr2=ddot(3*natoms,force_average2,1,force_average2,1)
  rnormpr2=dsqrt(rnormpr2 / float (3*natoms))



open(unit=761,file='OUTPUT.DTPOLY',access='append')
   call space(2)                                ! print information
   write(761,*)'                           Timestep # ',iii
   write(761,*)'                         ****************'
   write(761,*)'            Simulation Time   (fs):',float(iii)*tstep*1000.0d0
   write(761,*)'            Potential Energy (a.u):',ener
   write(761,*)'            Kinetic Energy   (a.u):' ,kinetic
   write(761,*)'            Avg. Etotal      (a.u):',etotal_av/float(iii)
   write(761,*)'            Temperature(unscaled)(K):',temp_unscaled
   write(761,*)'            Temperature(scaled)  (K):',temperature
   write(761,*)'            Avg. Temp.         (K):',temperature_av/float(iii)
   write(761,*)'            Pressure         (GPa):',pressure*auforce_to_bar*1d-4
   write(761,*)'            Pressure-avg     (GPa):',(paverage/float(iii))*auforce_to_bar*1d-4
   write(761,*)'            RMS      Norm         :',rnormpr
   write(761,*)'            Avg Force RMS Norm    :',rnormpr2
   write(761,*)'            Hugoniot Function     :',hg,hg_av/float(iii)
   write(761,*)'            Hugoniot Function(K)  :',hg_kelvin
   write(761,*)'            Shock Vel (km/s)      :',shock_vel/1000.0d0
   write(761,*)'            Particle Vel (km/s)   :',part_vel/1000.0d0



call space(2)
write(761,*)'                |----------------------------------|'
write(761,*)'                |         Stress Tensor (GPa)      |'
write(761,*)'                |----------------------------------|'
write(761,17)pressure_virial(1,1)*29421.0107637093d0,pressure_virial(1,2)*29421.0107637093d0,pressure_virial(1,3)*29421.0107637093d0
write(761,17)pressure_virial(2,1)*29421.0107637093d0,pressure_virial(2,2)*29421.0107637093d0,pressure_virial(2,3)*29421.0107637093d0
write(761,17)pressure_virial(3,1)*29421.0107637093d0,pressure_virial(3,2)*29421.0107637093d0,pressure_virial(3,3)*29421.0107637093d0
call space(2)



press_average_mat2=press_average_mat*29421.0107637093d0/float(iii)
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
write(23,24)iii,temperature,pressure*auforce_to_bar*1.0d-4,pwrite(1,1),pwrite(2,2),pwrite(3,3),pwrite(2,3),pwrite(1,3),pwrite(1,2),ener,kinetic+ener
close(23)
pwrite=press_average_mat2
open(unit=23,file='Averages',access='append')
write(23,24)iii,temperature_av/float(iii),(paverage/float(iii))*auforce_to_bar*1.0d-4,pwrite(1,1),pwrite(2,2),pwrite(3,3),pwrite(2,3),pwrite(1,3),pwrite(1,2),conav/float(iii),etotal_av/float(iii)
close(23)
24 format(i8,f12.3,f12.3,f12.3,f12.3,f12.3,f12.3,f12.3,f12.3,f18.5,f18.5)

open(unit=23,file='FORCES.AV',access='append')
write(23,34)iii,rnormpr,rnormpr2
close(23)
34 format(i8,f15.8,f15.8)




if(dft)then
pwrite=press_average_mat_unsymm*29421.0107637093d0/float(iii)
open(unit=33,file='Averages-unsymm')
rewind 33
write(33,*)pwrite(1,1),pwrite(1,2),pwrite(1,3)
write(33,*)pwrite(2,1),pwrite(2,2),pwrite(2,3)
write(33,*)pwrite(3,1),pwrite(3,2),pwrite(3,3)
close(33)
end if




if(iii .ge. istat_start)then
cell_stats=.true.
end if


rt=matmul(factor,rtdt)
cell2=matmul(factor,cell)
cell=cell2
call cell_volume(cell,volume)




call space(2)
write(761,46)'Average Vector Lengths ',av_cell_lengths(1)*0.5291772108d0/float(iii),av_cell_lengths(2)*0.5291772108d0/float(iii),av_cell_lengths(3)*0.5291772108d0/float(iii)
write(761,46)'Average Vector Angles',av_cell_angles(3)*180./pi/float(iii),av_cell_angles(2)*180./pi/float(iii),av_cell_angles(1)*180./pi/float(iii)
write(761,47)'Average Volume    ',av_cell_volume*au_to_ang**3/float(iii)
hugoniot_volume=av_cell_volume/float(iii)

open(unit=23,file='OUTPUT.CELL',access='append')
write(23,28)iii,av_cell_lengths(1)*0.5291772108d0/float(iii),av_cell_lengths(2)*0.5291772108d0/float(iii),av_cell_lengths(3)*0.5291772108d0/float(iii),&
av_cell_angles(3)*180./pi/float(iii),av_cell_angles(2)*180./pi/float(iii),av_cell_angles(1)*180./pi/float(iii),volume*au_to_ang**3,&
av_cell_volume*au_to_ang**3/float(iii)
close(23)

open(unit=23,file='OUTPUT.CELL.INSTANT',access='append')
write(23,48)iii,cell_lengths_write(1),cell_lengths_write(2),cell_lengths_write(3),cell_angles_write(3),cell_angles_write(2),cell_angles_write(1),volume*au_to_ang**3
close(23)
48 format(i8,f12.5,f12.5,f12.5,f12.5,f12.5,f12.5,f18.5)


28 format(i8,f12.5,f12.5,f12.5,f12.5,f12.5,f12.5,f18.5,f18.5)

write(761,*)'*******************************************************************************'
46 format(A30,f20.10,f20.10,f20.10)
47 format(A30,f20.10)
close(761)


if(hydrostatic)then
hugoniot_pressure= (press_average_mat(1,1)+press_average_mat(2,2)+press_average_mat(3,3))*29421.0107637093d0/(3.0d0*float(iii))
else
hugoniot_pressure= (press_average_mat(iuniaxial,iuniaxial))*29421.0107637093d0/(float(iii))
end if
hugoniot_energy=etotal_av/float(iii)
specific_energy=627.51d0*hugoniot_energy/(total_masskg*avogadro) ! kcal/kg
density=1.0D30 * total_masskg/(au_to_ang**3 * hugoniot_volume) ! in kg/m^3
specific_volume=1.0d0/density ! m^3/kg
hg=(specific_energy-specific_energy_ref) + 0.5d0 * (hugoniot_pressure_ref + hugoniot_pressure) * (specific_volume - specific_volume_ref) * rfactor



if(hydrostatic)then
 dhugoniot_pressure=( pressure_virial(1,1)+pressure_virial(2,2)+pressure_virial(3,3) )/3.0d0
else
 dhugoniot_pressure=( pressure_virial(iuniaxial,iuniaxial) )/1.0d0
end if
dspecific_energy=627.51d0*dhugoniot_energy/(total_masskg*avogadro) ! kcal/kg
ddensity=1.0D30 * total_masskg/(au_to_ang**3 * volume) ! in kg/m^3
dspecific_volume=1.0d0/ddensity ! m^3/kg
dhg=(dspecific_energy-specific_energy_ref) + 0.5d0 * (hugoniot_pressure_ref + dhugoniot_pressure) * (dspecific_volume - specific_volume_ref) * rfactor
hg_av=hg_av+dhg
!!    write(761,*)'            Hugoniot Function     :',hg

hg_kelvin=avogadro*hg*total_masskg/(627.51d0*hartree_per_eint * boltzman* (3*natoms-3.0d0) )


rfactor4=(1.0d0/29421.0107637093)*(627.51/avogadro)*(1.0d0/(0.52917D-10))**3*4184.0d0

shock_vel=dsqrt ( rfactor4* specific_volume_ref * (hugoniot_pressure - hugoniot_pressure_ref) / ( 1.0d0 - hugoniot_volume/hugoniot_volume_ref) )

!write(761,*)'Shock Velocity (km/s)',shock_vel/1000.0d0
part_vel=dsqrt ( rfactor4* specific_volume_ref * (hugoniot_pressure - hugoniot_pressure_ref) * ( 1.0d0 - hugoniot_volume/hugoniot_volume_ref) )
!write(761,*)'Particle Velocity (km/s)',part_vel/1000.0d0

if(hugoniot_run)then
if(iii.eq.ihugsteps .and. abs(hg_kelvin).gt.3.0d0)return
end if













if(iii.eq.1 .or. mod(iii,iwrite).eq.0)then
!   open(unit=215,file='OUTPUT.TRAJECTORY',access='append')
open(unit=215,file='OUTPUT.TRAJECTORY',access='append')
open(unit=315,file='OUTPUT.TRAJECTORY.CELL',access='append')
   write(215,*)natoms
   write(215,*)'timestep',iii
   do kk=1,natoms
   write(215,50)input_names(kk),rt(1,kk)*au_to_ang,rt(2,kk)*au_to_ang,rt(3,kk)*au_to_ang
   end do
close(215)
  write(315,*)'timestep',iii
  write(315,*)cell*au_to_ang
close(315)



! dump restart file
   open(unit=15,file='RESTART')
   write(15,*)iii
   write(15,*)cell
   write(15,*)rt
   write(15,*)velocity
   write(15,*)etotal_av/float(iii)
   write(15,*)(paverage/float(iii))*auforce_to_bar*1d-4
   write(15,*)av_cell_volume/float(iii)
   close(15)





end if
   50 format(A,f20.5,f20.5,f20.5)



!rt=matmul(factor,rtdt)
!cell2=matmul(factor,cell)
!cell=cell2
!call cell_volume(cell,volume)


end do ! end loop over time steps


call space(2)
print*,'               ******************************************'
print*,'               ******************************************'
print*,'               ******************************************'
print*,'               **********MD run has completed************'
print*,'               ******************************************'
print*,'               ******************************************'
print*,'               ******************************************'







end subroutine md_nose_f



!   if(i.eq.1 .or. mod(i,100).eq.0)then
!   write(15,*)natoms
!   write(15,*)natoms
!   do kk=1,natoms
!   write(15,50)input_names(kk),rt(1,kk)*au_to_ang,rt(2,kk)*au_to_ang,rt(3,kk)*au_to_ang
!   end do
!   end if
!   50 format(A,f10.5,f10.5,f10.5)
