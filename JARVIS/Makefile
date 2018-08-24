objects=header.o input.o forces.o bks.o control.o atom_types.o parse_comma.o md.o random.o ke_and_t.o matprt.o com_momentum.o arl_charm.o invert.o bfgs.o routines.o dblas1.o dblas2.o dblas3.o cpusec.o md_nose.o md_constraints.o constraints.o md_constraints_nose.o tersoff_key.o tersoff.o diff_theta.o lj_key.o lj_pot.o exp6_key.o exp6_pot.o morse_key.o morse_pot.o cell_key.o cell_volume.o md_berendsen.o neighbor_list.o md_berendsen_f.o angular.o bfgs_cell.o get_stress.o stress_key.o cg.o fcn.o cgsearch.o cgfam.o dump_final.o elastic.o eig.o strain_key.o strain_lattice.o md_nose_f.o lammps_pot.o rotate_cell.o rotate_cell_stress.o test_stress.o elastic_hessian.o hessian.o tred3.o tql3.o parallel_key.o cg_xyz.o hugoniot_key.o hugoniot.o hugoniotstat.o scale_key.o bfgs_all.o get_stress_forces.o

opt=    

compile90=ifort
compile77=ifort
md.exe: $(objects)
	$(compile90)  -o md.exe $(opt)   $(objects) 

header.o:header.f90
	$(compile90) -c  $(opt)  header.f90

input.o:input.f90 atom_types.f90 control.f90
	$(compile90) -c  $(opt)  input.f90

forces.o:forces.f90 input.f90 atom_types.f90
	$(compile90) -c  $(opt)  forces.f90

bks.o:bks.f90
	$(compile90) -c  $(opt)  bks.f90

control.o:control.f90
	$(compile90) -c  $(opt)  control.f90

atom_types.o:atom_types.f90
	$(compile90) -c  $(opt)  atom_types.f90

parse_comma.o:parse_comma.f90
	$(compile90) -c  $(opt)  parse_comma.f90

md.o:md.f90 atom_types.f90 input.f90 
	$(compile90) -c  $(opt)  md.f90

ke_and_t.o:ke_and_t.f90 atom_types.f90 input.f90
	$(compile90) -c  $(opt)  ke_and_t.f90

random.o:random.f
	$(compile77) -c  $(opt)  random.f

matprt.o:matprt.f
	$(compile77) -c  $(opt)  matprt.f

com_momentum.o:com_momentum.f90 atom_types.f90 input.f90 control.f90
	$(compile90) -c  $(opt)  com_momentum.f90

arl_charm.o:arl_charm.f
	$(compile77) -c  $(opt)  arl_charm.f

invert.o:invert.f
	$(compile77) -c  $(opt)  invert.f

bfgs.o: bfgs.f
	$(compile77) -c $(opt) bfgs.f

routines.o: routines.f
	$(compile77) -c $(opt)  routines.f

dblas1.o: dblas1.f
	$(compile77) -c  $(opt) dblas1.f

dblas2.o: dblas2.f
	$(compile77) -c  $(opt) dblas2.f

dblas3.o: dblas3.f
	$(compile77) -c  $(opt) dblas3.f

cpusec.o:cpusec.f
	$(compile77) -c  $(opt) cpusec.f

md_nose.o:md_nose.f90 input.f90 
	$(compile90) -c  $(opt)  md_nose.f90

md_constraints.o:md_constraints.f90 atom_types.f90 input.f90 control.f90 constraints.f90
	$(compile90) -c  $(opt)  md_constraints.f90

constraints.o:constraints.f90
	$(compile90) -c  $(opt)  constraints.f90

md_constraints_nose.o:md_constraints_nose.f90 atom_types.f90 input.f90 control.f90 constraints.f90
	$(compile90) -c  $(opt)  md_constraints_nose.f90

tersoff_key.o:tersoff_key.f90 input.f90
	$(compile90) -c  $(opt)  tersoff_key.f90

tersoff.o:tersoff.f90 input.f90 
	$(compile90) -c  $(opt)  tersoff.f90

diff_theta.o:diff_theta.f90
	$(compile90) -c  $(opt)  diff_theta.f90

lj_key.o:lj_key.f90 input.f90
	$(compile90) -c  $(opt)  lj_key.f90

lj_pot.o:lj_pot.f90 input.f90
	$(compile90) -c  $(opt)  lj_pot.f90

exp6_key.o:exp6_key.f90 input.f90
	$(compile90) -c  $(opt)  exp6_key.f90

exp6_pot.o:exp6_pot.f90 input.f90
	$(compile90) -c  $(opt)  exp6_pot.f90

morse_key.o:morse_key.f90 input.f90
	$(compile90) -c  $(opt)  morse_key.f90

morse_pot.o:morse_pot.f90 input.f90
	$(compile90) -c  $(opt)  morse_pot.f90

cell_key.o:cell_key.f90 input.f90
	$(compile90) -c  $(opt)  cell_key.f90

cell_volume.o::cell_volume.f90 input.f90
	$(compile90) -c  $(opt)  cell_volume.f90

md_berendsen.o:md_berendsen.f90 input.f90
	$(compile90) -c  $(opt)  md_berendsen.f90

neighbor_list.o:neighbor_list.f90 input.f90
	$(compile90) -c  $(opt)  neighbor_list.f90

md_berendsen_f.o:md_berendsen_f.f90 input.f90 control.f90
	$(compile90) -c  $(opt)  md_berendsen_f.f90

angular.o:angular.f90 input.f90
	$(compile90) -c  $(opt)  angular.f90

bfgs_cell.o: bfgs_cell.f
	$(compile77) -c $(opt) bfgs_cell.f

get_stress.o:get_stress.f90 input.f90 atom_types.f90
	$(compile90) -c  $(opt)  get_stress.f90

stress_key.o:stress_key.f90 input.f90
	$(compile90) -c  $(opt)  stress_key.f90

cg.o: cg.f
	$(compile77) -c $(opt) cg.f

fcn.o: fcn.f
	$(compile77) -c $(opt) fcn.f

cgsearch.o: cgsearch.f
	$(compile77) -c $(opt) cgsearch.f

cgfam.o: cgfam.f
	$(compile77) -c $(opt) cgfam.f

dump_final.o::dump_final.f90 input.f90
	$(compile90) -c  $(opt)  dump_final.f90

elastic.o::elastic.f90 input.f90
	 $(compile90) -c  $(opt) elastic.f90

eig.o:eig.f
	$(compile77) -c $(opt) eig.f

strain_key.o:strain_key.f90 input.f90
	$(compile90) -c  $(opt)  strain_key.f90

strain_lattice.o::strain_lattice.f90 input.f90
	$(compile90) -c  $(opt) strain_lattice.f90

md_nose_f.o:md_nose_f.f90 input.f90
	$(compile90) -c  $(opt)  md_nose_f.f90

lammps_pot.o:lammps_pot.f90 input.f90
	$(compile90) -c  $(opt)  lammps_pot.f90

rotate_cell.o:rotate_cell.f90 input.f90
	$(compile90) -c  $(opt)  rotate_cell.f90

rotate_cell_stress.o:rotate_cell_stress.f90 input.f90 rotate_cell.f90
	$(compile90) -c  $(opt)  rotate_cell_stress.f90

test_stress.o::test_stress.f90 input.f90
	$(compile90) -c  $(opt) test_stress.f90

elastic_hessian.o::elastic_hessian.f90 input.f90 
	$(compile90) -c  $(opt) elastic_hessian.f90

hessian.o::hessian.f90 input.f90
	$(compile90) -c  $(opt) hessian.f90

tql3.o: tql3.f
	$(compile77) -c    tql3.f

tred3.o: tred3.f
	$(compile77) -c   tred3.f

parallel_key.o:parallel_key.f90
	$(compile90) -c  $(opt)  parallel_key.f90

cg_xyz.o: cg_xyz.f
	$(compile77) -c $(opt) cg_xyz.f

hugoniot_key.o::hugoniot_key.f90 input.f90
	$(compile90) -c  $(opt) hugoniot_key.f90

hugoniot.o::hugoniot.f90 input.f90
	$(compile90) -c  $(opt) hugoniot.f90

hugoniotstat.o:hugoniotstat.f90 input.f90
	$(compile90) -c  $(opt)  hugoniotstat.f90

scale_key.o:scale_key.f90 input.f90
	$(compile90) -c  $(opt)  scale_key.f90

bfgs_all.o: bfgs_all.f
	$(compile77) -c $(opt) bfgs_all.f

get_stress_forces.o:get_stress_forces.f90 input.f90 atom_types.f90
	$(compile90) -c  $(opt)  get_stress_forces.f90
