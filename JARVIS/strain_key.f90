subroutine strain_key(keyword,ikey)
use control_module
use virial_mod
implicit double precision (a-h,o-z)

      !/////////////////////////////////////////////////////
      !Process STRAIN section of input file
      !/////////////////////////////////////////////////////
interface
subroutine parse_comma(word,value,length)
      character(200),intent(in)::word
      integer,intent(in)::length
      double precision,intent(inout)::value
end subroutine parse_comma
end interface



      character(200),intent(in):: keyword(:)
      integer,intent(in)::ikey
      character(200)word,ch1,ch
       integer,dimension(8)::icomma
! set defaults

strain_is_on=.true.
do i=2,ikey  
    k=len_trim(keyword(i))
    word=keyword(i)(1:k)

! find the 8 commas and store positions in icomma array
      istart=1 
      do j=1,8
         do m=istart,k
          if(word(m:m).eq.',')then
             icomma(j)=m
             istart=m+1
             exit
          end if   
         end do
      end do    

! parse string based on comma location


ch=word(1:icomma(1)-1)
call parse_comma(ch,value,icomma(1)-1)
strain_in(1,1)=value


ch=word(icomma(1)+1:icomma(2)-1)
call parse_comma(ch,value,icomma(2)-icomma(1)-1)
strain_in(2,1)=value

ch=word(icomma(2)+1:icomma(3)-1)
call parse_comma(ch,value,icomma(3)-icomma(2)-1)
strain_in(3,1)=value


ch=word(icomma(3)+1:icomma(4)-1)
call parse_comma(ch,value,icomma(4)-icomma(3)-1)
strain_in(1,2)=value

ch=word(icomma(4)+1:icomma(5)-1)
call parse_comma(ch,value,icomma(5)-icomma(4)-1)
strain_in(2,2)=value

ch=word(icomma(5)+1:icomma(6)-1)
call parse_comma(ch,value,icomma(6)-icomma(5)-1)
strain_in(3,2)=value

ch=word(icomma(6)+1:icomma(7)-1)
call parse_comma(ch,value,icomma(7)-icomma(6)-1)
strain_in(1,3)=value

ch=word(icomma(7)+1:icomma(8)-1)
call parse_comma(ch,value,icomma(8)-icomma(7)-1)
strain_in(2,3)=value

ch=word(icomma(8)+1:k)
call parse_comma(ch,value,k-icomma(8)-1)
strain_in(3,3)=value





end do








end subroutine strain_key
