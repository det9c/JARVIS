subroutine hugoniot_key(keyword,ikey)
use control_module
use hugoniot_module
      !/////////////////////////////////////////////////////
      !Processes parellel options for elastic constants
      !/////////////////////////////////////////////////////
      implicit double precision (a-h,o-z)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      character(200)word



! loop over possible keywords and process them

hugoniot_run=.true.
hugoniot_stat=.false.
hug_restart=.false.
hydrostatic=.true.
dampfac=1.0d0
     do i=1,ikey

      k=len_trim(keyword(i))
      word=keyword(i)(1:k)
      value=0.0d0
           if(word(1:4)=='TLOW')then
            call char_to_float(word,value,k)
            tlow=value
            elseif(word(1:5)=='THIGH')then
               call char_to_float(word,value,k)
            thigh=value
            elseif(word(1:7)=='RESTART')then
            hug_restart=.true.
            elseif(word(1:7)=='STAT')then
            hugoniot_stat=.true.
            hugoniot_run=.false.
            elseif(word(1:3)=='TAU')then
            call char_to_float(word,value,k)
            tau_stat=value
            elseif(word(1:8)=='TOL_TEMP')then
             call char_to_float(word,value,k)
             temp_tol_hug=value
            elseif(word(1:4)=='EREF')then
             call char_to_float(word,value,k)
             eref_input=value
            elseif(word(1:4)=='PREF')then
             call char_to_float(word,value,k)
             pref_input=value
            elseif(word(1:4)=='VREF')then
             call char_to_float(word,value,k)
             vref_input=value
            elseif(word(1:5)=='HYDRO')then
               hydrostatic=.true.
            elseif(word(1:3)=='UNI')then
               call char_to_float(word,value,k)
            iuniaxial=int(value)
            hydrostatic=.false.
           elseif(word(1:4)=='DAMP')then
            call char_to_float(word,value,k)
            dampfac=value

           end if


      if(keyword(i)(1:k).ne.' ')then
      write(761,10)keyword(i)(1:k),value
      end if

 
     end do

10 format(A20,F25.10)

end subroutine hugoniot_key
