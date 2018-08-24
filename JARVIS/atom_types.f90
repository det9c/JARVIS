subroutine atom_types(keyword,ikey)
use atom_types_module
implicit double precision (a-h,o-z)

      !/////////////////////////////////////////////////////
      !Process ATOM_TYPES section of input file
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
       integer,dimension(2)::icomma
! set defaults

! loop over possible keywords and process them

ntypes=ikey-1
allocate(type_name(ntypes))
allocate(atom_mass(ntypes))
allocate(atom_charge(ntypes))



! need to parse string of format NAME,MASS,CHARGE.
do i=2,ikey     ! loop over all atom types
    k=len_trim(keyword(i))
    word=keyword(i)(1:k)

! find the 2 commas and store positions in icomma array
      istart=1 
      do j=1,2
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
type_name(i-1)=word(1:icomma(1)-1)

!atom mass next
ch=word(icomma(1)+1:icomma(2)-1)
call parse_comma(ch,value,icomma(2)-icomma(1)-1)
atom_mass(i-1)=value

!atom charge last
ch=word(icomma(2)+1:k)
call parse_comma(ch,value,k-icomma(2))
atom_charge(i-1)=value



end do

end subroutine atom_types
