subroutine get_stress_forces(cell_bfgs,xyz_bfgs, stress_bfgs, energy_bfgs,vol,fboth)
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
double precision,dimension(3,natoms),intent(inout)::fboth
double precision,intent(in)::vol


vol_global=vol

cell=cell_bfgs
cell_elastic=cell_bfgs


if(elas_is_on)xyz_bfgs=xyz_elastic

call get_forces(xyz_bfgs, fboth, energy_bfgs, stress_bfgs)

xyz_elastic=xyz_bfgs


!stress_bfgs=virial_global
stress_elastic=stress_bfgs  !*29421.0107637093d0
stress_bfgs=stress_bfgs*29421.0107637093d0/vol !GPA




return
end subroutine get_stress_forces
