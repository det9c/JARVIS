subroutine parallel_key(keyword,ikey)
use control_module
use elastic_module
      !/////////////////////////////////////////////////////
      !Processes parellel options for elastic constants
      !/////////////////////////////////////////////////////
      implicit double precision (a-h,o-z)
      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      character(200)word



! loop over possible keywords and process them

parallel_elastic=.true.

     do i=1,ikey

      k=len_trim(keyword(i))
      word=keyword(i)(1:k)
      value=0.0d0

           if(word(1:2)=='E1')then
            e1_only=.true.
       elseif(word(1:2)=='E2')then
            e2_only=.true.
       elseif(word(1:2)=='E3')then
            e3_only=.true.
       elseif(word(1:2)=='E4')then
            e4_only=.true.
       elseif(word(1:2)=='E5')then
            e5_only=.true.
       elseif(word(1:2)=='E6')then
            e6_only=.true.
       elseif(word(1:5)=='FIRST')then
            call char_to_float(word,value,k)
            ifirst_hess_atom=int(value)
       elseif(word(1:4)=='LAST')then
            call char_to_float(word,value,k)
            ilast_hess_atom=int(value)


      end if

 
     end do



end subroutine parallel_key
