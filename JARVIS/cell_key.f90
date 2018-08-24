subroutine cell_key(keyword,ikey)
use cell_module
implicit double precision (a-h,o-z)

      !/////////////////////////////////////////////////////
      !Process CELL section of input file
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
       integer,dimension(3)::icomma
! set defaults

! loop over possible keywords and process them


! need to parse string of format VEC,ax,yx,zx
do i=2,ikey  
    k=len_trim(keyword(i))
    word=keyword(i)(1:k)

! find the 3 commas and store positions in icomma array
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
ch1=word(1:icomma(1)-1)

ch=word(icomma(1)+1:icomma(2)-1)
call parse_comma(ch,value,icomma(2)-icomma(1)-1)
xcomp=value


ch=word(icomma(2)+1:icomma(3)-1)
call parse_comma(ch,value,icomma(3)-icomma(2)-1)
ycomp=value

ch=word(icomma(3)+1:k)
call parse_comma(ch,value,k-icomma(3)-1)
zcomp=value


if(ch1.eq.'A')then
cell(1,1)=xcomp
cell(2,1)=ycomp
cell(3,1)=zcomp
elseif(ch1.eq.'B')then
cell(1,2)=xcomp
cell(2,2)=ycomp
cell(3,2)=zcomp
elseif(ch1.eq.'C')then
cell(1,3)=xcomp
cell(2,3)=ycomp
cell(3,3)=zcomp
end if




end do


open(unit=10,file='cellin')
!read(10,*)cell(1,1)
!read(10,*)cell(2,1)
!read(10,*)cell(3,1)
!read(10,*)cell(1,2)
!read(10,*)cell(2,2)
!read(10,*)cell(3,2)
!read(10,*)cell(1,3)
!read(10,*)cell(2,3)
!read(10,*)cell(3,3)
close(10)







end subroutine cell_key
