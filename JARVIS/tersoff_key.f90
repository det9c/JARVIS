subroutine tersoff_key(keyword,ikey)
use constraints_module
use tersoff_module
implicit double precision (a-h,o-z)

      !/////////////////////////////////////////////////////
      !Process tersoff section of input file
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
kk=0
do i=2,ikey
k=len_trim(keyword(i))
word=keyword(i)(1:k)
icom=0
do j=1,k
if(word(j:j).eq.',')icom=icom+1
end do
if(icom.eq.3)kk=kk+1
end do

! loop over possible keywords and process them

num_tersoff=ikey-1
num_tersoff=num_tersoff-kk
allocate(tersname(num_tersoff))
allocate(tersparams(11,num_tersoff))

allocate(terspair2(kk,2))
allocate(terspairval(kk,2))




tersoff_potential=.true.


! need to parse string of format Name,A,B,lambda,mu,beta,n,c,d,h,R,S
icounter=0
jcounter=0
do i=2,ikey     ! loop over all atomic parameters
    k=len_trim(keyword(i))
    word=keyword(i)(1:k)

! see if this is parameter string or chi/omega string
icom=0
do j=1,k
if(word(j:j).eq.',')icom=icom+1
end do


if(icom.eq.11)then
icounter=icounter+1
! find the commas and store positions in icomma array
      istart=1 
      do j=1,11
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
tersname(icounter)=word(1:icomma(1)-1)


do j=1,10
ch=word(icomma(j)+1:icomma(j+1)-1)
call parse_comma(ch,value,icomma(j+1)-icomma(j)-1)
tersparams(j,icounter)=value
end do

ch=word(icomma(11)+1:k)
call parse_comma(ch,value,k-icomma(11)-1)
tersparams(11,icounter)=value


else !process pair parameter string

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
terspair2(jcounter,1)=ch
ch=word(icomma(1)+1:icomma(2)-1)
terspair2(jcounter,2)=ch

ch=word(icomma(2)+1:icomma(3)-1)
call parse_comma(ch,value,icomma(3)-icomma(2)-1)
terspairval(jcounter,1)=value

ch=word(icomma(3)+1:k)
call parse_comma(ch,value,k-icomma(3)-1)
terspairval(jcounter,2)=value






end if







end do



end subroutine tersoff_key
