subroutine parse_comma(word,value,length)
!/////////////////////////////////////////////////////
!convert a character string to a floating point number
!/////////////////////////////////////////////////////
      implicit double precision (a-h,o-z)
      character(200),intent(in)::word
      integer,intent(in)::length
      double precision,intent(inout)::value

      nstart=1
      negative=scan(word,'-')
      sign=1.0d0
      if(negative/=0)then
      sign=-1.0d0
      nstart=nstart+1
      end if
      ndecimal=scan(word,'.')

      if(ndecimal.eq.0)then
       write(*,*)'Floating point input required'
       write(*,*)'Error:',word
       stop
      end if

! compute floating point total of keyword
      value=0.0d0
      nbefore=ndecimal-nstart
      iend=ndecimal-1
      do i=iend,nstart,-1
      j=ndecimal-i-1
      power=float(j)
      call atof(word(i:i),s)
      value=value+s*10**power
      end do

      istart=ndecimal+1
      do i=istart,length
      j=i-ndecimal
      power=float(j)*(-1.0d0)
      call atof(word(i:i),s)
      value=value+s*10**power
      end do

      value=value*sign
return
end subroutine parse_comma
