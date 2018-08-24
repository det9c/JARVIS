subroutine exp6_key(keyword,ikey)
use constraints_module
use exp6_module

implicit double precision (a-h,o-z)

      !/////////////////////////////////////////////////////
      !Process exp6 section of input file
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
      character(200)word,ch
       integer,dimension(11)::icomma
! set defaults
! count number of pair parameter sections 

! loop over possible keywords and process them

num_exp6=ikey-1
allocate(exp6pair(num_exp6,2))
allocate(exp6pairval(num_exp6,3))
exp6_potential=.true.


! need to parse string of format Name,Name,a,b,c

jcounter=0
do i=2,ikey     ! loop over all atomic parameters
    k=len_trim(keyword(i))
    word=keyword(i)(1:k)
! find the commas and store positions in icomma array
jcounter=jcounter+1
      istart=1
      do j=1,4
         do m=istart,k
          if(word(m:m).eq.',')then
             icomma(j)=m
             istart=m+1
             exit
          end if
         end do
      end do
! parse string based on comma location
! atom name first

ch=word(1:icomma(1)-1)
exp6pair(jcounter,1)=ch
ch=word(icomma(1)+1:icomma(2)-1)
exp6pair(jcounter,2)=ch

ch=word(icomma(2)+1:icomma(3)-1)
print*,ch
call parse_comma(ch,value,icomma(3)-icomma(2)-1)
exp6pairval(jcounter,1)=value

ch=word(icomma(3)+1:icomma(4)-1)
print*,ch
call parse_comma(ch,value,icomma(4)-icomma(3)-1)
exp6pairval(jcounter,2)=value

ch=word(icomma(4)+1:k)
print*,ch
call parse_comma(ch,value,k-icomma(4)-1)
exp6pairval(jcounter,3)=value





end do


end subroutine exp6_key
