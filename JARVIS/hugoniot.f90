subroutine hugoniot
use atom_types_module
use control_module
use constants
use cell_module
use virial_mod
use hugoniot_module
implicit double precision (a-h,o-z)
interface
   subroutine cell_volume(cell_in,volume)
   double precision,dimension(3,3),intent(in)::cell_in
   double precision,intent(inout)::volume
   end subroutine cell_volume




end interface
double precision,dimension(3,natoms)::xyz_original,velocity_original
double precision,dimension(3,3)::cell_original

ihugsteps=1000
tracker=1d10
avogadro=6.022E23
rfactor=(1.0d0/29421.0107637093)*(627.51/avogadro)*(1.0d0/(0.52917D-10))**3

open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,50)'    /__/\        /__/\        /  /\        /  /\        /__/\       ___         /  /\         ___    '
write(761,50)'    \  \:\       \  \:\      /  /:/_      /  /::\       \  \:\     /  /\       /  /::\       /  /\   '
write(761,50)'     \__\:\       \  \:\    /  /:/ /\    /  /:/\:\       \  \:\   /  /:/      /  /:/\:\     /  /:/   '
write(761,50)' ___ /  /::\  ___  \  \:\  /  /:/_/::\  /  /:/  \:\  _____\__\:\ /__/::\     /  /:/  \:\   /  /:/    '
write(761,50)'/__/\  /:/\:\/__/\  \__\:\/__/:/__\/\:\/__/:/ \__\:\/__/::::::::\\__\/\:\__ /__/:/ \__\:\ /  /::\    '
write(761,50)'\  \:\/:/__\/\  \:\ /  /:/\  \:\ /~~/:/\  \:\ /  /:/\  \:\~~\~~\/   \  \:\/\\  \:\ /  /://__/:/\:\   '
write(761,50)' \  \::/      \  \:\  /:/  \  \:\  /:/  \  \:\  /:/  \  \:\  ~~~     \__\::/ \  \:\  /:/ \__\/  \:\  '
write(761,50)'  \  \:\       \  \:\/:/    \  \:\/:/    \  \:\/:/    \  \:\         /__/:/   \  \:\/:/       \  \:\ '
write(761,50)'   \  \:\       \  \::/      \  \::/      \  \::/      \  \:\        \__\/     \  \::/         \__\/ '
write(761,50)'    \__\/        \__\/        \__\/        \__\/        \__\/                   \__\/                '
close(761)
50 format(A103)
if(restart)then
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'***********READING RESTART FILE************'
open(unit=15,file='RESTARTOLD')
read(15,*)junk
write(761,*)'Restarting from timestep ',junk
read(15,*)cell_original
read(15,*)xyz_original
!rtdt2=rtdt
read(15,*)velocity_original
read(15,*)hugoniot_energy_ref
read(15,*)hugoniot_pressure_ref
read(15,*)hugoniot_volume_ref
close(15)
close(761)
end if

! get values for the reference state
!open(unit=15,file='RESTART')
!write(15,*)junk
!write(15,*)cell_original
!write(15,*)xyz_original
!write(15,*)velocity_original
!close(15)




!compute reference values
   total_mass=zero
   do i=1,natoms
   total_mass=total_mass+atom_mass(type_atom(i))
   end do
total_mass=total_mass/(1000D0*avogadro) ! in kg
total_masskg=total_mass
density=1.0D30 * total_mass/(au_to_ang**3 * hugoniot_volume_ref) ! in kg/m^3
specific_energy_ref=627.51d0*hugoniot_energy_ref/(total_mass*avogadro) ! kcal/kg
specific_volume_ref=1.0d0/density ! m^3/kg
!print*,total_mass,density,specific_energy,specific_volume

open(unit=762,file='OUTPUT.DTPOLY.HUGONIOT',access='append')
call spaces(3)
write(762,*)'             Hugoniot Reference Point'
write(762,*)'            -------------------------            '
write(762,*)'          Eo (a.u.)=',hugoniot_energy_ref
write(762,*)'          Vo (Ang**3)=',hugoniot_volume_ref*au_to_ang**3
write(762,*)'          Po (GPa)=',hugoniot_pressure_ref
write(762,*)'          Specific Energy (kcal/kg)=',specific_energy_ref
write(762,*)'          Specific Volume (m^3/kg)=',specific_volume_ref
write(762,*)'          Density (kg/m^3)=',density
if(hydrostatic)then
write(762,*)'             Hydrostatic algorithm     '
else
write(762,*)'             Uniaxial algorithm for diagonal tensor component    ',iuniaxial
end if
write(762,*)'***************************************************************'

call spaces(3)
close(762)


if(hug_restart)goto 400


! get guess for low T

100 continue
temp_input=tlow
open(unit=15,file='RESTARTOLD')
write(15,*)junk
write(15,*)cell_original
write(15,*)xyz_original
write(15,*)velocity_original
write(15,*)hugoniot_energy_ref
write(15,*)hugoniot_pressure_ref
write(15,*)hugoniot_volume_ref
close(15)
call md_nose_f
specific_energy=627.51d0*hugoniot_energy/(total_mass*avogadro) ! kcal/kg
density=1.0D30 * total_mass/(au_to_ang**3 * hugoniot_volume) ! in kg/m^3
specific_volume=1.0d0/density ! m^3/kg
tlow_hg=(specific_energy-specific_energy_ref) + 0.5d0 * (hugoniot_pressure_ref + hugoniot_pressure) * (specific_volume - specific_volume_ref) * rfactor
open(unit=762,file='OUTPUT.DTPOLY.HUGONIOT',access='append')
call spaces(1)
write(762,*)'Hg using guess lower temperature',tlow_hg
close(762)
if(tlow_hg .gt. 0.0d0)then
tlow=tlow-50.0d0
open(unit=762,file='OUTPUT.DTPOLY.HUGONIOT',access='append')
write(762,*)'Lower bound is positive. Reducing temperature to',tlow
close(762)
goto 100
end if



! get guess for high T

200 continue
temp_input=thigh
open(unit=15,file='RESTARTOLD')
write(15,*)junk
write(15,*)cell_original
write(15,*)xyz_original
write(15,*)velocity_original
write(15,*)hugoniot_energy_ref
write(15,*)hugoniot_pressure_ref
write(15,*)hugoniot_volume_ref
close(15)
call md_nose_f
specific_energy=627.51d0*hugoniot_energy/(total_mass*avogadro) ! kcal/kg
density=1.0D30 * total_mass/(au_to_ang**3 * hugoniot_volume) ! in kg/m^3
specific_volume=1.0d0/density ! m^3/kg
thigh_hg=(specific_energy-specific_energy_ref) + 0.5d0 * (hugoniot_pressure_ref + hugoniot_pressure) * (specific_volume - specific_volume_ref) * rfactor
open(unit=762,file='OUTPUT.DTPOLY.HUGONIOT',access='append')
call spaces(1)
write(762,*)'Hg using guess higher temperature',thigh_hg
close(762)
if(thigh_hg .lt. 0.0d0)then
thigh=thigh+50.0d0
open(unit=762,file='OUTPUT.DTPOLY.HUGONIOT',access='append')
write(762,*)'Upper bound is negative. Increasing temperature to',thigh
close(762)
goto 200
end if

open(unit=762,file='OUTPUT.DTPOLY.HUGONIOT',access='append')
call spaces(2)
write(762,*)'A suitable interval has been located. Commencing search for root.'
write(762,*)'Search Interval:',tlow, 'to',thigh
call spaces(2)
close(762)


400 continue
if(hug_restart)then
open(unit=15,file='RESTART.HUGONIOT')
read(15,*)tlow,thigh
close(15)
open(unit=762,file='OUTPUT.DTPOLY.HUGONIOT',access='append')
call spaces(2)
write(762,*)'Restarting search.' 
write(762,*)'Search Interval read from restart file:',tlow, 'to',thigh
call spaces(2)
close(762)
end if






iteration=0.


300 continue
iteration=iteration+1
tmid=(tlow+thigh)/2.0d0
temp_input=tmid
open(unit=15,file='RESTARTOLD')
write(15,*)junk
write(15,*)cell_original
write(15,*)xyz_original
write(15,*)velocity_original
write(15,*)hugoniot_energy_ref
write(15,*)hugoniot_pressure_ref
write(15,*)hugoniot_volume_ref
close(15)

open(unit=15,file='RESTART.HUGONIOT')
write(15,*)tlow,thigh
close(15)

open(unit=762,file='OUTPUT.DTPOLY.HUGONIOT',access='append')
write(762,*)'Iteration number:',iteration
write(762,*)'Current MD temperature',tmid
write(762,*)'Current lower bound of root',tlow
write(762,*)'Current upper bound of root',thigh
write(762,*)'Width of search interval',abs(thigh-tlow)
close(762)
call md_nose_f
specific_energy=627.51d0*hugoniot_energy/(total_mass*avogadro) ! kcal/kg
density=1.0D30 * total_mass/(au_to_ang**3 * hugoniot_volume) ! in kg/m^3
specific_volume=1.0d0/density ! m^3/kg
tmid_hg=(specific_energy-specific_energy_ref) + 0.5d0 * (hugoniot_pressure_ref + hugoniot_pressure) * (specific_volume - specific_volume_ref) * rfactor
open(unit=762,file='OUTPUT.DTPOLY.HUGONIOT',access='append')
 write(762,*)'Current Hugoniot function value',tmid_hg
 write(762,*)'Current Hugoniot function value(K)',hg_kelvin
 call spaces(2)
 write(762,*)'E (a.u.)=',hugoniot_energy
 write(762,*)'V (Ang**3)=',hugoniot_volume*au_to_ang**3
 write(762,*)'P (GPa)=',hugoniot_pressure
 write(762,*)'Specific Energy (kcal/kg)=',specific_energy
 write(762,*)'Specific Volume (m^3/kg)=',specific_volume
 write(762,*)'Density (kg/m^3)=',density
rfactor4=(1.0d0/29421.0107637093)*(627.51/avogadro)*(1.0d0/(0.52917D-10))**3*4184.0d0
shock_vel=dsqrt ( rfactor4* specific_volume_ref * (hugoniot_pressure - hugoniot_pressure_ref) / ( 1.0d0 - hugoniot_volume/hugoniot_volume_ref) )

write(762,*)'Shock Velocity (km/s)',shock_vel/1000.0d0
part_vel=dsqrt ( rfactor4* specific_volume_ref * (hugoniot_pressure - hugoniot_pressure_ref) * ( 1.0d0 - hugoniot_volume/hugoniot_volume_ref) )
write(762,*)'Particle Velocity (km/s)',part_vel/1000.0d0

if(abs(tmid_hg).lt.tracker)then
tracker=abs(tmid_hg)
open(unit=234,file='lowest.hugoniot')
write(234,235)tmid_hg,hg_kelvin,hugoniot_pressure,hugoniot_volume*au_to_ang**3,tmid,part_vel/1000.0d0,shock_vel/1000.0d0
close(234)
235 format(f15.5,f15.5,f15.5,f15.5,f15.5,f15.5,f15.5)
end if


call spaces(2)
write(762,*)'------------------------------------------------------'
close(762)
40 format(A4,I3,A,f20.10,A,f20.10,A,f20.10,A,f20.10)
if(abs(hg_kelvin) .lt. 1.0d0)then
call spaces(2)
open(unit=762,file='OUTPUT.DTPOLY.HUGONIOT',access='append')
 write(762,*)'Hugoniot has converged:'
 call spaces(2)
 close(762)
 stop

 else

if(tmid_hg .lt. 0.0)tlow=tmid
if(tmid_hg .gt. 0.0)thigh=tmid
goto 300


end if














stop


























end subroutine hugoniot
subroutine spaces(n)
implicit integer (a-z)
do i=1,n
write(762,*)''
end do
end subroutine spaces
