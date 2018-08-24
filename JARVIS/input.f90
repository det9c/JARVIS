!////////////////////////////////////////////////////
!          modules for ARL stuff         //////
!/////////////////////////////////////////////////////

module hugoniot_module
   logical::hugoniot_run,hug_restart,hugoniot_stat,hydrostatic
   double precision::hugoniot_volume,hugoniot_energy,hugoniot_pressure,tlow,thigh  &
,specific_energy_ref,hugoniot_pressure_ref,specific_volume_ref,rfactor,total_masskg,avogadro,tau_stat,hg_kelvin, &
hugoniot_volume_ref,temp_tol_hug,eref_input,pref_input,vref_input,dampfac
 integer::iuniaxial,ihugsteps
end module hugoniot_module



module rotations
  double precision::rot1,rot2,rot3,rot4,rot5,rot6
  double precision,dimension(3,3)::rotmat
end module rotations



module parameters_bks
  double precision:: asio,bsio,csio,aoo,boo,coo,asisi,bsisi, &
                 csisi,qsi,qo 
end module parameters_bks

module cell_module
double precision,dimension(3,3)::cell,cell_copy,cell_flex,scale_factor
double precision,dimension(3)::cell_lengths,cell_angles,av_cell_lengths,av_cell_angles,&
cell_lengths_write,cell_angles_write

double precision,dimension(:,:),allocatable::fracs
logical::cubic,cell_stats,just_volume
double precision::cell_vol,av_cell_volume,vol_global,volume_strain
end module cell_module


module potentials
  logical:: bks,sherwin, eddq,charm
end module potentials

module control_module
  logical:: PBC,scale,predictor,verlet,bfgs_opt,dynamics,nose,berendsen,berendsen_f,angular_potential,cell_opt,cgopt, &
            elastic,refcell,dft,restart,strain_is_on,nose_f,lammps,metal,lreals,test_tens,linear,usenrm,restart_cell_md, &
            restart_cell_opt,oldcp2k,hessian,parallel_elastic,opt_bfgs,magic,fixed_frame,espresso,cp2k,mopac,skipelas, &
            aprun,lsda

  double precision::A_VEC,B_VEC,C_VEC,temp_input,tstep,cutoff,pressure_input,tau,cutoff_sq,strength,temp_tol, &
rmsmax,frcmax,delta_hessian,cellrms,afrc,cfrc,arms,crms
  integer::iscale_counter,isteps,iwrite,maxcyc,istat_start,iupdate_list,ncell,maxcyc2,optalgo
end module control_module

module atom_types_module
character(4),allocatable,dimension(:)::type_name,input_names
double precision,allocatable,dimension(:)::atom_mass,atom_charge
integer:: ntypes,natoms
double precision,allocatable,dimension(:,:)::xyz
integer,allocatable,dimension(:)::type_atom,num_type
end module atom_types_module

module constants
 ! module for internal constants
 ! mass -> amu
 !velocity bohr/ps
 ! time -> picosecond
 ! length ->bohr
 ! energy -> quasi-hartree amu*bohr^2/ps^2
double precision,parameter:: au_to_ang = 0.5291772108d0 !convert bohr to angstrom
double precision,parameter:: hartree_per_eint = 1.066572345d-6    ! convert hartree to internal unit of energy
double precision,parameter:: boltzman = 2.969151809d0           ! boltzman constant in internal energy units (Eo/K)  
double precision,parameter:: intvel_to_mpers = 1.889726125d-2 ! convert bohr/ps to m/s
double precision,parameter::zero=0.0d0  !boltzman=3.166d-6 atomic units Hartree/Kelvin
double precision,parameter::twothirds=.6666666666666666d0,half=.5d0
double precision::pi
double precision::auforce_to_bar=294210107.637093d0,temp_unscaled
end module constants

module constraints_module
integer::num_constraints
integer,allocatable,dimension(:,:)::iconstraint_pairs
double precision,allocatable,dimension(:)::constraint_values
logical::constrained_md
end module constraints_module


module tersoff_module
integer::num_tersoff
double precision,allocatable,dimension(:,:)::tersparams,ters_pair,terspairval

character(5),allocatable,dimension(:)::tersname
character(5),allocatable,dimension(:,:)::terspair2
logical::tersoff_potential
integer,allocatable,dimension(:)::terstype,offset
end module tersoff_module

module lj_module
integer::num_lj
character(5),allocatable,dimension(:,:)::ljpair
double precision,allocatable,dimension(:,:)::ljpairval,lj_paramval
logical::lj_potential
integer,allocatable,dimension(:)::offset_lj
end module lj_module

module exp6_module
integer::num_exp6
character(5),allocatable,dimension(:,:)::exp6pair
double precision,allocatable,dimension(:,:)::exp6pairval,exp6_paramval
logical::exp6_potential
integer,allocatable,dimension(:)::offset_exp6
end module exp6_module

module morse_module
integer::num_morse
character(5),allocatable,dimension(:,:)::morsepair
double precision,allocatable,dimension(:,:)::morsepairval,morse_paramval
logical::morse_potential
integer,allocatable,dimension(:)::offset_morse
end module morse_module

module neighbors_module
integer,allocatable,dimension(:)::num_neighbors,ifirst_neighbor,ipair_list
end module neighbors_module

module virial_mod
double precision,dimension(3,3)::stress_tens,press_tens,unit_matrix,virial_global,stress_in,strain_in,virial_unsymm
double precision,allocatable,dimension(:,:)::fxyz_global
end module virial_mod

module elastic_module
double precision,dimension(3,3)::cell_elastic,stress_elastic
double precision,allocatable,dimension(:,:)::xyz_elastic,hessmat,hessold
logical::elas_is_on,e1_only,e2_only,e3_only,e4_only,e5_only,e6_only
double precision::eps_elas
integer::ifirst_hess_atom,ilast_hess_atom
end module elastic_module


!/////////////////////////////////////////////////////
!          end of modules for ARL stuff         //////
!/////////////////////////////////////////////////////







subroutine arl_input
use potentials
use atom_types_module
use control_module
use constraints_module 
use tersoff_module
use lj_module
use exp6_module
use morse_module
use constants
use cell_module
use neighbors_module
use virial_mod
use elastic_module
use hugoniot_module
!/////////////////////////////////////////////////////
!  This subroutine loops over the supplementary
!  input file and reads the user defined keywords
!  for each potential. 
!/////////////////////////////////////////////////////
 implicit double precision (a-h,o-z)

      interface

         subroutine bks_key(keyword,ikey)
         character(200),intent(in):: keyword(:)
         integer,intent(in)::ikey
         end subroutine bks_key


         subroutine control(keyword,ikey)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      end subroutine control

      subroutine atom_types(keyword,ikey)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      end subroutine atom_types

      subroutine constraints(keyword,ikey)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      end subroutine constraints

         subroutine tersoff_key(keyword,ikey)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      end subroutine tersoff_key


         subroutine lj_key(keyword,ikey)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      end subroutine lj_key

         subroutine exp6_key(keyword,ikey)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      end subroutine exp6_key

      subroutine morse_key(keyword,ikey)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      end subroutine morse_key

      subroutine cell_key(keyword,ikey)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      end subroutine cell_key

     subroutine stress_key(keyword,ikey)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      end subroutine stress_key


     subroutine strain_key(keyword,ikey)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      end subroutine strain_key

     subroutine parallel_key(keyword,ikey)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      end subroutine parallel_key


     subroutine hugoniot_key(keyword,ikey)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      end subroutine hugoniot_key


     subroutine scale_key(keyword,ikey)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      end subroutine scale_key




      end interface
!************************************************************
 character(1000)line
 character(200),dimension(:),allocatable::keyword
 character(200)ch
 integer num_keywords(20)
 double precision,dimension(3,3)::stress_input,dummy
 double precision,dimension(34)::genetic

!*************************************************************

call space(2)

!set some defaults and set up some constants
scale_factor=1.0d0
strain_is_on=.false.
bks=.false.
constrained_md=.false.
tersoff_potential=.false.
lj_potential=.false.
exp6_potential=.false.
cubic=.false.
elas_is_on=.false.
rmsmax=1d-4
frcmax=1d-4
cellrms=1d-4
just_volume=.false.
elastic=.false.
parallel_elastic=.false.
e1_only=.false.
e2_only=.false.
e3_only=.false.
e4_only=.false.
e5_only=.false.
e6_only=.false.
hugoniot_run=.false.





pi=dacos(-1.0d0)
unit_matrix=zero
do i=1,3
unit_matrix(i,i)=1.0d0
end do


open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'                 PROGRAM CONTROLS SET BY USER:'
write(761,*)'          Character String        Internal Value'
write(761,*)'          ----------------        --------------'

open(unit=765,file='INPUT')
 isection=0 
 line(1:1000)=' '

do ! start of loop over input file sections
   read(765,*,iostat=io)line
   if(io<0)exit
   k=len_trim(line)
   ikey=1
   do i=1,k 
      ch=line(i:i)
      if(ch==' ')then
         ikey=ikey+1
      end if
   end do

  allocate(keyword(ikey))

  ibegin=1
  do i=1,ikey
     icount=0
     ch=' '
     do j=ibegin,k
        if(line(j:j)==' ')then
          ibegin=j+1
          exit
          else
          icount=icount+1
          ch(icount:icount)=line(j:j)
        end if
     end do
     keyword(i)=ch
   end do

!process input section by type

   k=len_trim(keyword(1))
   iexist=0
   if(keyword(1)(1:k).eq.'*POTENTIAL*')then
      iexist=1
      k=len_trim(keyword(2))
      select case(keyword(2)(1:k))
          case('BKS')
               call bks_key(keyword,ikey)
      end select
   end if


   if(keyword(1)(1:k).eq.'*CONTROL*')then
      iexist=1
      call control(keyword,ikey)
   end if

   if(keyword(1)(1:k).eq.'*ATOM_TYPES*')then
      iexist=1
      call atom_types(keyword,ikey)
   end if

   if(keyword(1)(1:k).eq.'*CONSTRAINTS*')then
      iexist=1
      call constraints(keyword,ikey)
   end if

   if(keyword(1)(1:k).eq.'*TERSOFF*')then
      iexist=1
      call tersoff_key(keyword,ikey)
   end if

   if(keyword(1)(1:k).eq.'*LJ*')then
      iexist=1
      call lj_key(keyword,ikey)
   end if

   if(keyword(1)(1:k).eq.'*EXP6*')then
      iexist=1
      call exp6_key(keyword,ikey)
   end if

   if(keyword(1)(1:k).eq.'*MORSE*')then
      iexist=1
      call morse_key(keyword,ikey)
   end if

   if(keyword(1)(1:k).eq.'*CELL*')then
      iexist=1
      call cell_key(keyword,ikey)
   end if

  if(keyword(1)(1:k).eq.'*STRESS*')then
      iexist=1
      call stress_key(keyword,ikey)
   end if

  if(keyword(1)(1:k).eq.'*STRAIN*')then
      iexist=1
      call strain_key(keyword,ikey)
   end if

  if(keyword(1)(1:k).eq.'*PARALLEL*')then
      iexist=1
      call parallel_key(keyword,ikey)
   end if
  if(keyword(1)(1:k).eq.'*HUGONIOT*')then
      iexist=1
      call hugoniot_key(keyword,ikey)
   end if

  if(keyword(1)(1:k).eq.'*SCALE*')then
      iexist=1
      call scale_key(keyword,ikey)
   end if




   if(iexist.eq.0)then
    write(*,*)'Error in INPUT file. Section ',keyword(1)(1:k),' not recognized'
    stop
    end if




   if(allocated(keyword))deallocate(keyword)

end do ! end outer loop over input file sections. 





!  print atom types defined by user. info provided here may 
! grow as code develops
call space(3)
write(761,*)'                     SUMMARY OF INPUT ATOM TYPES           '
write(761,*)'-----------------------------------------------------------------'
write(761,*)'Type      Atom Name               Mass                   Charge'
do i=1,ntypes 
write(761,20)i,'            ',type_name(i),'       ',atom_mass(i),'          ',atom_charge(i)
end do
20 format(I3,A,A,A,f14.3,A,f14.3)






! read coordinate file
natoms=0
open(unit=20,file='XYZ')
do 
read(20,*,iostat=io)ch,rjunk,rjunk,rjunk
if(io<0)exit
natoms=natoms+1
end do





allocate( xyz(3,natoms))
if(elastic)allocate(xyz_elastic(3,natoms))
if(hessian)then
allocate(fxyz_global(3,natoms))
allocate(hessmat(3*natoms,3*natoms))
allocate(hessold(3*natoms,3*natoms))
end if
allocate( input_names(natoms))
allocate(type_atom(natoms))
rewind 20
do i=1,natoms
read(20,*,iostat=io)input_names(i),xyz(1,i),xyz(2,i),xyz(3,i)
end do

allocate(num_type(ntypes))
num_type=0
do i=1,natoms
ch=input_names(i)
do j=1,ntypes
if(ch==type_name(j))exit
end do
type_atom(i)=j
num_type(j)=num_type(j)+1
end do


call space(3)
write(761,*)'                                 INPUT GEOMETRY                         '
write(761,*)'       -----------------------------------------------------------------'
write(761,*)'Atom Number   Atom Name       Type               X            Y            Z '
do i=1,natoms
write(761,30)' ',i,'            ',input_names(i),'        ',type_atom(i),'        ',xyz(1,i),xyz(2,i),xyz(3,i)
end do
30 format(A,I6,A,A,A,I3,A,f13.5,f13.5,f13.5)

call space(1)

isum=0
chg=0.
do i=1,ntypes
write(761,*)'             There are ',num_type(i),' atoms of kind ',type_name(i)
isum=isum+num_type(i)
chg=chg+num_type(i)*atom_charge(i)
end do
write(761,*)'             There are ',isum,' atoms total'
write(761,*)'             System Charge is ',chg


!compute reference values
   total_mass=zero
   do i=1,natoms
   total_mass=total_mass+atom_mass(type_atom(i))
   end do
write(761,*)'             Mass is ',total_mass,'g/mol'




write(761,*)'-------------------------------------------------------------------------------'
if(constrained_md)then
call space(2)
write(761,*)'               Summary of user defined bond constraints                  '
write(761,*)'       Constraint  Atom 1       Atom 2      Value'
do i=1,num_constraints
write(*,*)i,iconstraint_pairs(i,1),iconstraint_pairs(i,2),'  ',constraint_values(i)
end do
end if

call space(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!lattice vectors!!!!!!!!!!!!!!!!!!!!1 
do i=1,3
tot=0.
do j=1,3
tot=tot+cell(j,i)**2
end do
cell_lengths(i)=dsqrt(tot)
end do
tot=0.
!angles
kk=0
do i=1,2
do j=i+1,3
kk=kk+1
tot=cell(1,i)*cell(1,j)+cell(2,i)*cell(2,j)+cell(3,i)*cell(3,j)
cell_angles(kk)=tot/ ( cell_lengths(i) * cell_lengths(j) )
end do
end do
 write(761,*)'                       Lattice Vectors               '
 write(761,*)' Vector            X-Comp            Y-Comp         Z-Comp          Length' 
 write(761,900)'    A->          ',cell(1,1),'        ' ,cell(2,1),'     ',cell(3,1),'      ',cell_lengths(1)
 write(761,900)'    B->          ',cell(1,2),'        ' ,cell(2,2),'     ',cell(3,2),'      ',cell_lengths(2)
 write(761,900)'    C->          ',cell(1,3),'        ' ,cell(2,3),'     ',cell(3,3),'      ',cell_lengths(3)
cell=cell/au_to_ang
cell_lengths=cell_lengths/au_to_ang
cell_angles=acos(cell_angles)
 write(761,901)'Vector Angles    ',cell_angles(3)*180./pi,'        ',cell_angles(2)*180./pi,'     ',cell_angles(1)*180./pi
!if(cell_angles(1)*180./pi.eq.90.0 .and. cell_angles(2)*180./pi.eq.90.0 .and. cell_angles(3)*180./pi.eq.90.0 .and. PBC)then
!cubic=.true.
!call space(1)
!write(761,*)'***************Cubic Periodic Boundary Conditions Will Be Used****************'
!else
allocate(fracs(3,natoms))
!end if


900 format(A,f10.5,A,f10.5,A,f10.5,A,f10.5)
901 format(A,f10.5,A,f10.5,A,f10.5)











call space(3)





if(tersoff_potential)then

write(761,*)'reading parameters from genetic algorithm'
!open(unit=23,file='params')
!do i=1,34
!read(23,*)genetic(i)
!end do
!close(23)





write(761,*)'                Summary of Tersoff Potential Parameters               '
write(761,*)'              (See Phys. Rev. B Vol. 39, Number 8. (1989)                 '

 allocate(terstype(natoms))
 do i=1,natoms
 ch=input_names(i)                  
 do j=1,num_tersoff
 if(ch==tersname(j))exit         
 end do                           
 terstype(i)=j                    
 end do            
 call space(2)
 do i=1,num_tersoff
 write(761,*)'                        Atom Name:',tersname(i)
 write(761,*)'  --------------------------------------------------------------'
 write(761,*)'         A                      B                     lambda  '
 write(761,*)tersparams(1,i),tersparams(2,i),tersparams(3,i)
 write(761,*)'         mu                    beta                      n'
 write(761,*)tersparams(4,i),tersparams(5,i),tersparams(6,i)
 write(761,*)'         c                      d                        h '
 write(761,*)tersparams(7,i),tersparams(8,i),tersparams(9,i)
 write(761,*)'         R                      S'
 write(761,*)tersparams(10,i),tersparams(11,i)
 call space(2)
 tersparams(10,i)=tersparams(10,i)/au_to_ang
 tersparams(11,i)=tersparams(11,i)/au_to_ang
write(761,78)tersname(i),tersparams(1,i)*27.21,tersparams(3,i)/.52917,tersparams(2,i)*27.21,tersparams(4,i)/.52917,tersparams(10,i)*.52917
write(761,79)tersparams(11,i)*.52917,tersparams(5,i),tersparams(6,i),tersparams(7,i),tersparams(8,i),tersparams(9,i)

print*,tersparams(1,i)*27.21138386d0
print*,tersparams(3,i)/au_to_ang
print*,tersparams(2,i)*27.21138386d0
print*,tersparams(4,i)/au_to_ang
print*,tersparams(10,i)*au_to_ang
print*,tersparams(11,i)*au_to_ang
print*,tersparams(5,i)
print*,tersparams(6,i)
print*,tersparams(7,i)
print*,tersparams(8,i)
print*,tersparams(9,i)


 end do
78 format(A,f13.5,f13.5,f13.5,f13.5,f13.5)
79 format(f13.5,f13.5,f13.5,f13.5,f13.5,f13.5)


!compute pair parameters and report
 iterms=num_tersoff*(num_tersoff+1)/2
 allocate(ters_pair(8,iterms))

 allocate(offset(num_tersoff))
 do i=1,num_tersoff
 offset(i)=i*(i-1)/2
 end do


 do i=1,num_tersoff
 do j=1,i
 ij=j+offset(i)
 ters_pair(1,ij)=(tersparams(3,i)+tersparams(3,j))/2.0d0 ! lambda for pair
 ters_pair(2,ij)=(tersparams(4,i)+tersparams(4,j))/2.0d0 ! mu for pair 
 ters_pair(3,ij)=dsqrt(tersparams(1,i)*tersparams(1,j))  ! A for pair
 ters_pair(4,ij)=dsqrt(tersparams(2,i)*tersparams(2,j))  ! B for pair
 ters_pair(5,ij)=dsqrt(tersparams(10,i)*tersparams(10,j))! R for pair 
 ters_pair(6,ij)=dsqrt(tersparams(11,i)*tersparams(11,j))! S for pair
 end do 
 end do

! load pair parameters chi and omega into working arrays
do i=1,iterms
ch=terspair2(i,1)
do j=1,num_tersoff
if(ch==tersname(j))exit
end do
indi=j
ch=terspair2(i,2)
do j=1,num_tersoff
if(ch==tersname(j))exit
end do
indj=j
ii=max(indi,indj)
jj=min(indi,indj)
ij=jj+offset(ii)
ters_pair(7,ij)=terspairval(i,1)
ters_pair(8,ij)=terspairval(i,2)
print*,ters_pair(7,ij)
print*,ters_pair(8,ij)
end do

!ters_pair(3,1)=genetic(1)
!ters_pair(4,1)=genetic(2)
!ters_pair(1,1)=genetic(3)
!ters_pair(2,1)=genetic(4)
!ters_pair(3,2)=genetic(5)
!ters_pair(4,2)=genetic(6)
!ters_pair(1,2)=genetic(7)
!ters_pair(2,2)=genetic(8)
!ters_pair(3,3)=genetic(9)
!ters_pair(4,3)=genetic(10)
!ters_pair(1,3)=genetic(11)
!ters_pair(2,3)=genetic(12)




!ters_pair(7,1)=genetic(13)
!ters_pair(8,1)=genetic(14)
!ters_pair(7,2)=genetic(15)
!ters_pair(8,2)=genetic(16)
!ters_pair(7,3)=genetic(17)
!ters_pair(8,3)=genetic(18)


!tersparams(5,1)=genetic(19)
!tersparams(6,1)=genetic(20)
!tersparams(7,1)=genetic(21)
!tersparams(8,1)=genetic(22)
!tersparams(9,1)=genetic(23)
!tersparams(5,2)=genetic(24)
!tersparams(6,2)=genetic(25)
!tersparams(7,2)=genetic(26)
!tersparams(8,2)=genetic(27)
!tersparams(9,2)=genetic(28)


!ters_pair(5,1)=genetic(29)/au_to_ang
!ters_pair(6,1)=genetic(30)/au_to_ang
!ters_pair(5,2)=genetic(31)/au_to_ang
!ters_pair(6,2)=genetic(32)/au_to_ang
!ters_pair(5,3)=genetic(33)/au_to_ang
!ters_pair(6,3)=genetic(34)/au_to_ang




!ters_pair(3,1)=genetic(1)
!ters_pair(4,1)=genetic(2)
!ters_pair(1,1)=genetic(3)
!ters_pair(2,1)=genetic(4)
!ters_pair(3,2)=genetic(5)
!ters_pair(4,2)=genetic(6)
!ters_pair(1,2)=genetic(7)
!ters_pair(2,2)=genetic(8)
!ters_pair(3,3)=genetic(9)
!ters_pair(4,3)=genetic(10)
!ters_pair(1,3)=genetic(11)
!ters_pair(2,3)=genetic(12)




!ters_pair(7,1)=genetic(13)
!ters_pair(8,1)=genetic(14)
!ters_pair(7,2)=genetic(15)
!ters_pair(8,2)=genetic(16)
!ters_pair(7,3)=genetic(17)
!ters_pair(8,3)=genetic(18)

!write lammps parameters 
 do i=1,num_tersoff
 do j=1,num_tersoff
 do k=1,num_tersoff
         indi=i
         indj=j
         ii=max(indi,indj)
         jj=min(indi,indj)
         ij=jj+offset(ii)

         indk=k
         ii=max(indi,indk)
         jj=min(indi,indk)
         ik=jj+offset(ii)
         Rmink=ters_pair(5,ik)
         Smink=ters_pair(6,ik)
         womegak=ters_pair(8,ik)




          beta=tersparams(5,i)
          d_n=tersparams(6,i)
          c=tersparams(7,i)
          d=tersparams(8,i)
          h=tersparams(9,i)
         A=ters_pair(3,ij)
         B=ters_pair(4,ij)
         dlam=ters_pair(1,ij)
         dmu=ters_pair(2,ij)
         Rmin=ters_pair(5,ij)
         Smin=ters_pair(6,ij)
         chi=ters_pair(7,ij)
         womega=ters_pair(8,ij)

write(761,*)tersname(i),tersname(j),tersname(k), 1.0,womegak,0.0,c,d,h,d_n,beta,dmu/au_to_ang,27.21138386d0*B*chi,au_to_ang*(Rmink+Smink)/2.0d0,au_to_ang*(Smink-Rmink)/2.0d0,dlam/au_to_ang,A*27.21138386d0
end do
end do
end do







 write(761,*)'Tersoff Pair Parameters:'
 call space(1)
 do i=1,num_tersoff
 do j=1,i
 ij=j+offset(i)
 write(761,*)'                        Atom Pair:',tersname(j),'<->    ',tersname(i)
 write(761,*)'  --------------------------------------------------------------'
 write(761,*)'         A                      B                     lambda  '
 write(761,*)ters_pair(3,ij),ters_pair(4,ij),ters_pair(1,ij)
 write(761,*)'         mu                     R                        S'
 write(761,*)ters_pair(2,ij),ters_pair(5,ij),ters_pair(6,ij)
 write(761,*)'         chi                     omega                 '
 write(761,*)ters_pair(7,ij),ters_pair(8,ij)
 call space(2)
write(761,28)tersname(i),tersname(j),ters_pair(3,ij)*27.21,ters_pair(4,ij)*27.21,ters_pair(1,ij)/.52917,ters_pair(2,ij)/.52917,&
ters_pair(5,ij)*.52917,ters_pair(6,ij)*.52917
 end do
 end do
end if
28 format(A,A,f20.10,f20.10,f20.10,f20.10,f20.10,f20.10)




!!!!!!!!!!!!!!!!!!
if(lj_potential)then
 print*,'              Summary of Lennard-Jones Potential Parameters               '
 call space(2)
!compute pair parameters and report
 iterms=ntypes*(ntypes+1)/2
 allocate(lj_paramval(2,iterms))
 lj_paramval=0.
 allocate(offset_lj(ntypes))
 do i=1,ntypes
 offset_lj(i)=i*(i-1)/2
 end do
! load pair parameters epsilon and sigma into working arrays
do i=1,num_lj
ch=ljpair(i,1)
do j=1,ntypes
if(ch==type_name(j))exit
end do
indi=j
ch=ljpair(i,2)
do j=1,ntypes
if(ch==type_name(j))exit
end do
indj=j
ii=max(indi,indj)
jj=min(indi,indj)
ij=jj+offset_lj(ii)
lj_paramval(1,ij)=ljpairval(i,1)
lj_paramval(2,ij)=ljpairval(i,2)
end do
 print*,'Lennard-Jones Pair Parameters:'
 call space(1)
 do i=1,ntypes
 do j=1,i
 ij=j+offset_lj(i)
 print*,'                        Atom Pair:',type_name(j),'<->    ',type_name(i)
 print*,'  --------------------------------------------------------------'
 print*,'         eps                    sigma '
 write(*,*)lj_paramval(1,ij),lj_paramval(2,ij)
 call space(2)
 end do
 end do
end if



!!!!!!!!!!!!!!!!!!
if(exp6_potential)then
 print*,'                        Summary of Exp-6 Potential Parameters               '
 call space(2)
!compute pair parameters and report
 iterms=ntypes*(ntypes+1)/2
 allocate(exp6_paramval(3,iterms))
 exp6_paramval=0.
 allocate(offset_exp6(ntypes))
 do i=1,ntypes
 offset_exp6(i)=i*(i-1)/2
 end do
! load pair parameters a,b,c into working arrays
do i=1,num_exp6
ch=exp6pair(i,1)
do j=1,ntypes
if(ch==type_name(j))exit
end do
indi=j
ch=exp6pair(i,2)
do j=1,ntypes
if(ch==type_name(j))exit
end do
indj=j
ii=max(indi,indj)
jj=min(indi,indj)
ij=jj+offset_exp6(ii)
exp6_paramval(1,ij)=exp6pairval(i,1)
exp6_paramval(2,ij)=exp6pairval(i,2)
exp6_paramval(3,ij)=exp6pairval(i,3)
end do
 print*,'Exp-6 Pair Parameters:'
 call space(1)
 do i=1,ntypes
 do j=1,i
 ij=j+offset_exp6(i)
 print*,'                        Atom Pair:',type_name(j),'<->    ',type_name(i)
 print*,'  --------------------------------------------------------------'
 print*,'         A                      B                        C  '
 write(*,*)exp6_paramval(1,ij),exp6_paramval(2,ij),exp6_paramval(3,ij)
 call space(2)
 end do
 end do
end if


!!!!!!!!!!!!!!!!!!
if(morse_potential)then
 print*,'                        Summary of Morse Potential Parameters               '
 call space(2)
!compute pair parameters and report
 iterms=ntypes*(ntypes+1)/2
 allocate(morse_paramval(3,iterms))
 morse_paramval=0.
 allocate(offset_morse(ntypes))
 do i=1,ntypes
 offset_morse(i)=i*(i-1)/2
 end do
! load pair parameters a,b,c into working arrays
do i=1,num_morse
ch=morsepair(i,1)
do j=1,ntypes
if(ch==type_name(j))exit
end do
indi=j
ch=morsepair(i,2)
do j=1,ntypes
if(ch==type_name(j))exit
end do
indj=j
ii=max(indi,indj)
jj=min(indi,indj)
ij=jj+offset_morse(ii)
morse_paramval(1,ij)=morsepairval(i,1)
morse_paramval(2,ij)=morsepairval(i,2)
morse_paramval(3,ij)=morsepairval(i,3)
end do
 print*,'Morse Pair Parameters:'
 call space(1)
 do i=1,ntypes
 do j=1,i
 ij=j+offset_morse(i)
 print*,'                        Atom Pair:',type_name(j),'<->    ',type_name(i)
 print*,'  --------------------------------------------------------------'
 print*,'         Eo                     k                      R(eq)  '
 write(*,*)morse_paramval(1,ij),morse_paramval(2,ij),morse_paramval(3,ij)
 call space(2)
 end do
 end do
end if








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  move the atoms by bfgs or MD or both
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(.not. dynamics .and. .not. bfgs_opt)then
write(761,*)'Error: User has not selected BFGS or MD algorithm'
stop
end if 


! !
!                     chop users CP2K.INPUT file into 3 pieces
if(cp2k)then
   call system('wc -l INPUT.CP2K>num_lines')
   open(unit=1,file='num_lines')
   read(1,*)num_lines
   close(1)
   call system('rm -f num_lines')
   open(unit=1,file='INPUT.CP2K')
   do i=1,num_lines
      read(1,*,iostat=io)line
      if(line.eq.'*CELL')i1=i
      if(line.eq.'*XYZ')i2=i
   end do
   close(1)

   ihead=i1-1
   open(unit=1,file='chop')
   write(1,*)'head -n ',ihead,' INPUT.CP2K > part1'
   close(1)
   call system('chmod 755 chop ')
   call system(' ./chop ')

   ibetween=i2-i1-1
   itail=num_lines-i2
   open(unit=1,file='chop')
   write(1,*)'tail -n ',itail,' INPUT.CP2K > part3'
   close(1)
   call system('chmod 755 chop ')
   call system(' ./chop ')

   ihead=i2-1
   open(unit=1,file='chop')
   write(1,*)'head -n ',ihead,' INPUT.CP2K > part4'
   close(1)
   call system('chmod 755 chop ')
   call system(' ./chop ')

   open(unit=1,file='chop')
   write(1,*)'tail -n ',ibetween,' part4 > part2'
   close(1)
   call system('chmod 755 chop ')
   call system(' ./chop ')
   call system('rm part4')
end if


if(espresso)then
   call system('wc -l INPUT.QE>num_lines')
   open(unit=1,file='num_lines')
   read(1,*)num_lines
   close(1)
   call system('rm -f num_lines')
   open(unit=1,file='INPUT.QE')
   do i=1,num_lines
      line=''
      read(1,*,iostat=io)line
!!      print*,i,line,i1,i2,i3
      if(line.eq.'*GUESS')i1=i
      if(line.eq.'*CELL')i2=i
      if(line.eq.'*XYZ')i3=i
!!   print*,i,line(1:10)
   end do
   close(1)


   ihead=i1-1
   open(unit=1,file='chop')
   write(1,*)'head -n ',ihead,' INPUT.QE > part1'
   close(1)
   call system('chmod 755 chop ')
   call system(' ./chop ')


   ibetween=i2-i1-1
   ihead=i2-1
   open(unit=1,file='chop')
   write(1,*)'head -n ',ihead,' INPUT.QE > part4'
   close(1)
   call system('chmod 755 chop ')
   call system(' ./chop ')

   open(unit=1,file='chop')
   write(1,*)'tail -n ',ibetween,' part4 > part2'
   close(1)
   call system('chmod 755 chop ')
   call system(' ./chop ')
   call system('rm part4')



   ibetween=i3-i2-1
   ihead=i3-1
   open(unit=1,file='chop')
   write(1,*)'head -n ',ihead,' INPUT.QE > part4'
   close(1)
   call system('chmod 755 chop ')
   call system(' ./chop ')

   open(unit=1,file='chop')
   write(1,*)'tail -n ',ibetween,' part4 > part3'
   close(1)
   call system('chmod 755 chop ')
   call system(' ./chop ')
   call system('rm part4')


   itail=num_lines-i3
   open(unit=1,file='chop')
   write(1,*)'tail -n ',itail,' INPUT.QE > part4'
   close(1)
   call system('chmod 755 chop ')
   call system(' ./chop ')

end if




























close(761)


















!  set up pair lists
allocate(ifirst_neighbor(natoms))
allocate(num_neighbors(natoms))






if(test_tens)then
call test_stress
stop
end if


!if(hessian)then
!xyz=xyz/0.5291772108d0
!call get_hessian
!stop
!end if



if(bfgs_opt )then
xyz=xyz/0.5291772108d0


if(restart_cell_md)then
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'***********READING RESTART FILE************'
open(unit=15,file='RESTARTOLD')
read(15,*)junk
write(761,*)'Restarting from timestep ',junk
read(15,*)cell
read(15,*)xyz
close(15)
close(761)
end if

if(restart_cell_opt)then
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'***********READING BEST_POINT FILE************'
open(unit=15,file='best_point')
read(15,*)svalue,svalue2
write(761,*)'gradient norm of restart point',svalue,svalue2
read(15,*)cell
read(15,*)xyz
close(15)
close(761)
end if












if(.not. cell_opt)then
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'IF DOING A PLAIN OPTIMIZATION, HAVE TO PUT VOLUME IN CALL '
write(761,*)'TO SUBROUTINE IN INPUT.F90 '
close(761)
STOP
call bfgs(xyz,natoms,frcmax,rmsmax,input_names,maxcyc,dplace,1.0,dummy)
else
!call matprt(cell,3,3,3,3)

!stress_input=unit_matrix*pressure_input*1D-4 ! pressure in GPa

!stress_input(1,1)=0.
!stress_input(2,2)=4.0



do i=1,3
do j=1,3
if(stress_in(j,i).ne.0.0d0)usenrm=.true.
end do
end do




if(.not. cgopt)then
    if(optalgo.eq.3)then
      call bfgs_all(xyz,natoms,input_names,stress_in,usenrm,maxcyc2,transpose(cell),afrc,cfrc,arms,crms,fixed_frame,scale_factor)
    else
      call  bfgs_cell(xyz,natoms,input_names,maxcyc2,transpose(cell),stress_in,refcell,usenrm,optalgo,fixed_frame,cfrc,crms &
                      ,scale_factor)
    end if
else
call conjugate(xyz,natoms,frcmax,input_names,maxcyc2,transpose(cell),stress_in,refcell,usenrm)
end if


if(elastic)then
 if(hessian)then
 call elascon_hess
 else
 call elascon
 end if
end if

end if

end if

if(dynamics)then


!make sure cutoff does not exceed 1/2 box length for PBC runs before beginning
! vectors are still angstrom at this point. fix later
if(pbc)then
halfa=cell_lengths(1)/2.0d0
halfb=cell_lengths(2)/2.0d0
halfc=cell_lengths(3)/2.0d0
istop=0
if(cutoff.gt.halfa)then
write(761,*)'Error: Potential cutoff exceeds 1/2 length of A lattice vector'
istop=1
end if
if(cutoff.gt.halfb)then
write(761,*)'Error: Potential cutoff exceeds 1/2 length of B lattice vector'
istop=1
end if
if(cutoff.gt.halfc)then
write(761,*)'Error: Potential cutoff exceeds 1/2 length of C lattice vector'
istop=1
end if
if(istop.eq.1)stop
end if



close(761)
!  set up pair lists
!allocate(ifirst_neighbor(natoms))
!allocate(num_neighbors(natoms))

if(restart_cell_opt)then
open(unit=761,file='OUTPUT.DTPOLY',access='append')
write(761,*)'***********READING BEST_POINT FILE TO START MD************'
open(unit=15,file='best_point')
read(15,*)svalue,svalue2
write(761,*)'gradient norm of restart point',svalue,svalue2
read(15,*)cell
read(15,*)xyz
xyz=xyz*au_to_ang
close(15)
close(761)
end if



if(hugoniot_stat)call hugoniotstat
if(hugoniot_run)call hugoniot
if(.not. nose .and. .not. constrained_md .and. .not. berendsen .and. .not. berendsen_f .and. .not. nose_f .and. .not. hugoniot_run)call MD
if(.not. nose .and. constrained_md) call md_constraints
if(nose .and. .not. constrained_md .and. .not. nosef)call MD_NOSE
if(nose .and.  constrained_md .and. .not. nosef)call md_constraints_nose
if(berendsen)call md_berendsen
if(berendsen_f)call md_berendsen_f
if(nose_f)call md_nose_f



end if



200 continue
return






end subroutine arl_input




subroutine bks_key(keyword,ikey)
      !/////////////////////////////////////////////////////
      !Processes BKS keywords
      !/////////////////////////////////////////////////////
      use parameters_bks
      use potentials
      implicit double precision (a-h,o-z)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      character(200)word


! set defaults



! loop over possible keywords and process them




     do i=1,ikey
      
      k=len_trim(keyword(i))
!      write(*,*)keyword(i)(1:k)
      word=keyword(i)(1:k)
      value=0.0d0

!
!
! probably should have just called char_to_float
! from outside if/then block then it would not
! have appeared on every line below and would be easier
! to read. would have been a problem for T/F keywords though.  
!
!

        if(word(1:4)=='ASIO')then
            call char_to_float(word,value,k)
            asio=value

        elseif(word(1:4)=='BSIO')then
            call char_to_float(word,value,k)
            bsio=value

        elseif(word(1:4)=='CSIO')then
            call char_to_float(word,value,k)
            csio=value

        elseif(word(1:3)=='AOO')then
            call char_to_float(word,value,k)
            aoo=value

        elseif(word(1:3)=='BOO')then
            call char_to_float(word,value,k)
            boo=value

        elseif(word(1:3)=='COO')then
            call char_to_float(word,value,k)
            coo=value

        elseif(word(1:5)=='ASISI')then
            call char_to_float(word,value,k)
            asisi=value

        elseif(word(1:5)=='BSISI')then
            call char_to_float(word,value,k)
            bsisi=value

        elseif(word(1:5)=='CSISI')then
            call char_to_float(word,value,k)
            csisi=value

        elseif(word(1:3)=='QSI')then
            call char_to_float(word,value,k)
            qsi=value

        elseif(word(1:2)=='QO')then
            call char_to_float(word,value,k)
            qo=value

        elseif(word(1:3)=='BKS')then
            bks=.true.
 


                
      end if

      if(keyword(i)(1:k).ne.' ')then
      write(*,10)keyword(i)(1:k),value
      end if
     end do

10 format(A15,F25.10)
end subroutine bks_key




subroutine char_to_float(word,value,length)
!/////////////////////////////////////////////////////
!convert a character string to a floating point number
!/////////////////////////////////////////////////////
      implicit double precision (a-h,o-z)
      character(50),intent(in)::word
      integer,intent(in)::length
      double precision,intent(inout)::value


      nstart=scan(word,'=')
      nstart=nstart+1
      sign=1.0d0
      negative=scan(word,'-')
      if(negative/=0)then
      sign=-1.0d0   
      nstart=nstart+1
      end if
      ndecimal=scan(word,'.')

      if(ndecimal.eq.0)then
       write(*,*)'Floating point input required'
       write(*,*)'Error:',word
       stop
      end if 

! compute floating point total of keyword
      value=0.0d0
      nbefore=ndecimal-nstart
      iend=ndecimal-1
      do i=iend,nstart,-1
      j=ndecimal-i-1
      power=float(j)
      call atof(word(i:i),s)
      value=value+s*10**power
      end do


      istart=ndecimal+1
      do i=istart,length
      j=i-ndecimal
      power=float(j)*(-1.0d0)
      call atof(word(i:i),s)
      value=value+s*10**power
      end do


      value=value*sign
     
return
end subroutine char_to_float

subroutine atof(input,s)
character,intent(in)::input
double precision,intent(inout)::s
if(input=='0')then
s=0.
elseif(input=='1')then
s=1.
elseif(input=='2')then
s=2.
elseif(input=='3')then
s=3.
elseif(input=='4')then
s=4.
elseif(input=='5')then
s=5.
elseif(input=='6')then
s=6.
elseif(input=='7')then
s=7.
elseif(input=='8')then
s=8.
elseif(input=='9')then
s=9.
end if
end subroutine atof




