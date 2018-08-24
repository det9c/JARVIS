subroutine control(keyword,ikey)
use control_module
use elastic_module
      !/////////////////////////////////////////////////////
      !Processes control keywords
      !/////////////////////////////////////////////////////
      implicit double precision (a-h,o-z)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      character(200)word

! set defaults
predictor=.false.
verlet=.false.
bfgs_opt=.false.
cell_opt=.false.
cgopt=.false.
dynamics=.false.
nose=.false.
maxcyc=100
berendsen=.false.
berendsen_f=.false.
pressure_input=0.
iupdate_list=1000.
angular_potential=.false.
ncell=0
strength=0.
refcell=.false.
dft=.false.
espresso=.false.
cp2k=.false.
mopac=.false.
restart=.false.
cell_stats=.false.
nose_f=.false.

lammps=.false.
metal=.false.
lreals=.false.
test_tens=.false.
linear=.false.
usenrm=.false.
restart_cell_md=.false.
restart_cell_opt=.false.
oldcp2k=.false.
maxcyc2=100
hessian=.false.
delta_hessian=1.0d-3
opt_bfgs=.true.
optalgo=0 !0=Fixed Fracs |  1=Updated fracs | 2=XYZ algorithm | 3=Cell and cart in one shot
fixed_frame=.false.
skipelas=.false.
aprun=.false.
lsda = .false.
! loop over possible keywords and process them


     do i=1,ikey

      k=len_trim(keyword(i))
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

        if(word(1:3)=='PBC')then
            pbc=.true.

!        elseif(word(1:5)=='A_VEC')then
!            call char_to_float(word,value,k)
!            a_vec=value
!        elseif(word(1:5)=='B_VEC')then
!            call char_to_float(word,value,k)
!            b_vec=value
!        elseif(word(1:5)=='C_VEC')then
!            call char_to_float(word,value,k)
!            c_vec=value
        elseif(word(1:4)=='TEMP')then
            call char_to_float(word,value,k)
            temp_input=value
        elseif(word(1:11)=='SCALE_EVERY')then
            call char_to_float(word,value,k)
            scale=.true.
            iscale_counter=int(value)
        elseif(word(1:5)=='STEPS')then
            call char_to_float(word,value,k)
            isteps=int(value)
         elseif(word(1:9)=='PREDICTOR')then
            predictor=.true.
         elseif(word(1:9)=='TIME_STEP')then
           call char_to_float(word,value,k)
           tstep=value
         elseif(word(1:6)=='VERLET')then
            verlet=.true.
         elseif(word(1:6)=='CUTOFF')then
            call char_to_float(word,value,k)
            cutoff=value
            cutoff=cutoff/0.5291772108d0
            cutoff_sq=cutoff*cutoff
         elseif(word(1:4)=='BFGS')then
            bfgs_opt=.true.
         elseif(word(1:2)=='MD')then
            dynamics=.true.
         elseif(word(1:4)=='NOSE')then
            nose=.true.
         elseif(word(1:10)=='TRAJECTORY')then
            call char_to_float(word,value,k)
            iwrite=int(value)
         elseif(word(1:6)=='MAXCYC')then
            call char_to_float(word,value,k)
            maxcyc=int(value)
         elseif(word(1:11)=='F_BERENDSEN')then
            berendsen_f=.true.

         elseif(word(1:9)=='BERENDSEN')then
            berendsen=.true.
         elseif(word(1:8)=='PRESSURE')then
            call char_to_float(word,value,k)
            pressure_input=value
         elseif(word(1:3)=='TAU')then
            call char_to_float(word,value,k)
            tau=value
         elseif(word(1:9)=='START_AVG')then
            call char_to_float(word,value,k)
            istat_start=int(value)
            
         elseif(word(1:6)=='UPDATE')then
            call char_to_float(word,value,k)
            iupdate_list=int(value)
       
         elseif(word(1:7)=='ANGULAR')then
            angular_potential=.true.

         elseif(word(1:5)=='NCELL')then
            call char_to_float(word,value,k)
            ncell=int(value)
         elseif(word(1:8)=='STRENGTH')then
            call char_to_float(word,value,k)
            strength=value


         elseif(word(1:9)=='CONJUGATE')then
            cgopt=.true.
         elseif(word(1:7)=='ELASTIC')then
            elastic=.true.
         elseif(word(1:3)=='EPS')then
            call char_to_float(word,value,k)
            eps_elas=value
         elseif(word(1:7)=='REFCELL')then
            refcell=.true.

         elseif(word(1:7)=='DFT')then
            dft=.true.

        elseif(word(1:8)=='TOL_TEMP')then
            call char_to_float(word,value,k)
            temp_tol=value
        elseif(word(1:7)=='RESTART')then
            restart=.true.
         elseif(word(1:6)=='F_NOSE')then
            nose_f=.true.

            elseif(word(1:6)=='LAMMPS')then
            lammps=.true.
         elseif(word(1:6)=='METALS')then
            metal=.true.
            lreals=.false.
            elseif(word(1:5)=='REALS')then
             metal=.false.
            lreals=.true.
         elseif(word(1:9)=='TESTTENS')then
            test_tens=.true.
         elseif(word(1:6)=='LINEAR')then
            linear=.true.
         elseif(word(1:8)=='USE_NORM')then
         usenrm=.true.
         elseif(word(1:15)=='CELL_MD_RESTART')then
         restart_cell_md=.true.
         elseif(word(1:16)=='CELL_OPT_RESTART')then
         restart_cell_opt=.true.
         elseif(word(1:8)=='CELL_OPT')then
            cell_opt=.true.
         elseif(word(1:6)=='RMSMAX')then
            call char_to_float(word,value,k)
            rmsmax=value
         elseif(word(1:6)=='FRCMAX')then
            call char_to_float(word,value,k)
            frcmax=value
        elseif(word(1:7)=='OLDCP2K')then
           oldcp2k=.true.
        elseif(word(1:10)=='CELL_STEPS')then
            call char_to_float(word,value,k)
            maxcyc2=int(value)
        elseif(word(1:7)=='HESSIAN')then
           hessian=.true.
         elseif(word(1:2)=='DX')then
            call char_to_float(word,value,k)
            delta_hessian=value
          elseif(word(1:6)=='USE_CG')then
            opt_bfgs=.false.
        elseif(word(1:7)=='OPTALGO')then
            call char_to_float(word,value,k)
            optalgo=int(value)
        elseif(word(1:11)=='FIXED_FRAME')then
           fixed_frame=.true.
         elseif(word(1:7)=='CELLRMS')then
            call char_to_float(word,value,k)
            cellrms=value
         elseif(word(1:4)=='AFRC')then
            call char_to_float(word,value,k)
            afrc=value
         elseif(word(1:4)=='CFRC')then
            call char_to_float(word,value,k)
            cfrc=value
         elseif(word(1:4)=='ARMS')then
            call char_to_float(word,value,k)
            arms=value
         elseif(word(1:4)=='CRMS')then
            call char_to_float(word,value,k)
            crms=value
          elseif(word(1:8)=='ESPRESSO')then
            espresso=.true.
          elseif(word(1:4)=='CP2K')then
            cp2k=.true.
          elseif(word(1:5)=='MOPAC')then
            mopac=.true.
          elseif(word(1:8)=='SKIPELAS')then
            skipelas=.true.
          elseif(word(1:5)=='APRUN')then
            aprun=.true.
          elseif(word(1:5)=='LSDA')then
            lsda=.true.



      end if

      if(keyword(i)(1:k).ne.' ')then
      write(761,10)keyword(i)(1:k),value
      end if
     end do


     if(cgopt)opt_bfgs=.true.



10 format(A20,F25.10)
end subroutine control
