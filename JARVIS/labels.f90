subroutine arl_labels(ax,ay,az,atmsym,nat,rxyz0, fxyz, ener)
  use parameters_sherwin !module with parameters read from supp. input file 
  implicit double precision (a-h,o-z)
  integer::nlabel


  open(unit=100,file='labels')
  read(100,*)nlabel
  allocate(labels(nlabel))
  do i=1,nlabel
  read(100,*)labels(i)
  end do 


  num_si=0
  num_o=0
  num_h=0
  num_o1=0

! count silicon atoms
  do i=1,nlabel
  if(labels(i).eq.'SI')num_si=num_si+1
  if(labels(i).eq.'O ' .or. labels(i).eq.'O1')num_o=num_o+1
  if(labels(i).eq.'H ')num_h=num_h+1
  if(labels(i).eq.'O1')num_o1=num_o1+1
  end do

  print*,'Number of silicon atoms:',num_si
  print*,'Number of general oxygen atoms:',num_o
  print*,'Number of silanol oxygen atoms:',num_o1
  print*,'Number of hydrogen atoms:',num_h 

  allocate(si_list(num_si))
  allocate(o_list(num_o))
  allocate(h_list(num_h))
  allocate(o1_list(num_o1))

  irow=0
  jrow=0
  krow=0
  lrow=0


   do i=1,nlabel
     if(labels(i).eq.'SI')then
        irow=irow+1
        si_list(irow)=i

      elseif(labels(i).eq.'O ')then
         jrow=jrow+1
         o_list(jrow)=i

      elseif(labels(i).eq.'H ')then
         krow=krow+1
         h_list(krow)=i


      elseif(labels(i).eq.'O1')then
         lrow=lrow+1
         o1_list(lrow)=i
         jrow=jrow+1
         o_list(jrow)=i

      end if   
   end do


print*,'si',si_list

print*,'o list',o_list

print*,'h list',h_list


print*,'o1 list',o1_list








end subroutine arl_labels

