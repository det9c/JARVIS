subroutine scale_key(keyword,ikey)
use cell_module
implicit double precision (a-h,o-z)

      !/////////////////////////////////////////////////////
      !Process STRESS section of input file
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

! loop over possible keywords and process them


! need to parse string of format xx,yx,zx,xy,yy,zy,xz,yz,zz e.g load by columns
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
scale_factor(1,1)=value


ch=word(icomma(1)+1:icomma(2)-1)
call parse_comma(ch,value,icomma(2)-icomma(1)-1)
scale_factor(1,2)=value

ch=word(icomma(2)+1:icomma(3)-1)
call parse_comma(ch,value,icomma(3)-icomma(2)-1)
scale_factor(1,3)=value


ch=word(icomma(3)+1:icomma(4)-1)
call parse_comma(ch,value,icomma(4)-icomma(3)-1)
scale_factor(2,1)=value

ch=word(icomma(4)+1:icomma(5)-1)
call parse_comma(ch,value,icomma(5)-icomma(4)-1)
scale_factor(2,2)=value

ch=word(icomma(5)+1:icomma(6)-1)
call parse_comma(ch,value,icomma(6)-icomma(5)-1)
scale_factor(2,3)=value

ch=word(icomma(6)+1:icomma(7)-1)
call parse_comma(ch,value,icomma(7)-icomma(6)-1)
scale_factor(3,1)=value

ch=word(icomma(7)+1:icomma(8)-1)
call parse_comma(ch,value,icomma(8)-icomma(7)-1)
scale_factor(3,2)=value



ch=word(icomma(8)+1:k)
call parse_comma(ch,value,k-icomma(8)-1)
scale_factor(3,3)=value





end do







end subroutine scale_key
