subroutine get_stress(cell_bfgs,xyz_bfgs, stress_bfgs, energy_bfgs,vol)
use atom_types_module
use control_module
use cell_module
use constants
use virial_mod
use elastic_module

implicit double precision (a-h,o-z)
double precision,dimension(3,natoms)::xyz_bfgs
double precision,dimension(3,3)::cell_bfgs,stress_bfgs
double precision,intent(inout)::energy_bfgs
double precision,intent(in)::vol


vol_global=vol

cell=cell_bfgs
cell_elastic=cell_bfgs


if(elas_is_on)xyz_bfgs=xyz_elastic

if(opt_bfgs)then
call bfgs(xyz_bfgs,natoms,afrc,arms,input_names,maxcyc,energy_bfgs,vol,stress_bfgs)
else
call cg_xyz(xyz_bfgs,natoms,frcmax,rmsmax,input_names,maxcyc,energy_bfgs,vol,stress_bfgs)
end if

xyz_elastic=xyz_bfgs


!stress_bfgs=virial_global
stress_elastic=stress_bfgs  !*29421.0107637093d0
stress_bfgs=stress_bfgs*29421.0107637093d0/vol !GPA




return
end subroutine get_stress
