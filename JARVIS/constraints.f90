subroutine constraints(keyword,ikey)
use constraints_module
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

num_constraints=ikey-1
if(num_constraints .ne.0)constrained_md=.true.
allocate(iconstraint_pairs(num_constraints,2))
allocate(constraint_values(num_constraints))




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
ch=word(1:icomma(1)-1)
call parse_comma(ch,value,icomma(1)-1)
iconstraint_pairs(i-1,1)=int(value)

!atom mass next
ch=word(icomma(1)+1:icomma(2)-1)
call parse_comma(ch,value,icomma(2)-icomma(1)-1)
iconstraint_pairs(i-1,2)=int(value)

!atom charge last
ch=word(icomma(2)+1:k)
call parse_comma(ch,value,k-icomma(2))
constraint_values(i-1)=value



end do

!call space(2)
!print*,'               Summary of user defined bond constraints                  '
!print*,'       Constraint  Atom 1       Atom 2      Value'                        
!do i=1,num_constraints
!write(*,*)i,iconstraint_pairs(i,1),iconstraint_pairs(i,2),'  ',constraint_values(i)
!end do


end subroutine constraints
