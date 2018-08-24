subroutine lj_key(keyword,ikey)
use constraints_module
use lj_module

implicit double precision (a-h,o-z)

      !/////////////////////////////////////////////////////
      !Process lj section of input file
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
print*,'in lj'

num_lj=ikey-1
allocate(ljpair(num_lj,2))
allocate(ljpairval(num_lj,2))
lj_potential=.true.


! need to parse string of format Name,Name,epsilon,sigma

jcounter=0
do i=2,ikey     ! loop over all atomic parameters
    k=len_trim(keyword(i))
    word=keyword(i)(1:k)
! find the commas and store positions in icomma array
jcounter=jcounter+1
      istart=1
      do j=1,3
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
ljpair(jcounter,1)=ch
ch=word(icomma(1)+1:icomma(2)-1)
ljpair(jcounter,2)=ch

ch=word(icomma(2)+1:icomma(3)-1)
call parse_comma(ch,value,icomma(3)-icomma(2)-1)
ljpairval(jcounter,1)=value

ch=word(icomma(3)+1:k)
call parse_comma(ch,value,k-icomma(3)-1)
ljpairval(jcounter,2)=value




end do

print*,'out'

end subroutine lj_key
