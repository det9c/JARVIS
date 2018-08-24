      subroutine fcn( n, x, f, g ,natoms,xyz,input_names,force,
     $stress_in,itimes,fracs,refcell,usenrm)
      implicit double precision (a-h,o-z)
      dimension x(n),g(n),xyz(3,natoms),stress_in(3,3),cell(3,3)
     $,cell_copy(3,3),stress_copy(3,3)
     $,fracs(3,natoms),stress(3,3),cinv(3,3)
     $ ,smat(3,3),indx(3),cell_deriv(3,3)     
     $,sinv(3,3)
      character*4 input_names(natoms)
      integer itimes
      logical refcell,usenrm
      save flowest

      factor=.5291772108d0
      delstr=0.001

        irow=0
        do i=1,3
        do j=1,3
           irow=irow+1
        cell(i,j)=x(irow)
        end do
        end do




        open(unit=761,file='OUTPUT.DTPOLY',access='append')
         write(761,*)''
         write(761,*)' -------------------------------'
         write(761,*)'|  Iteration Number',itimes,'   |'
         write(761,*)' -------------------------------'
         f=0.0D0


         if(itimes.eq.1)then
            flowest=1d10
                     do i=1,3
         do j=1,3
c         cell(j,i)=cell(j,i)*factor
         cell_copy(j,i)=cell(j,i)
         end do
         end do


         call migs(cell_copy,3,sinv,indx)

                      do i=1,natoms
          ssx=sinv(1,1)*xyz(1,i)+sinv(2,1)*xyz(2,i)+sinv(3,1)*xyz(3,i)
          ssy=sinv(1,2)*xyz(1,i)+sinv(2,2)*xyz(2,i)+sinv(3,2)*xyz(3,i)
          ssz=sinv(1,3)*xyz(1,i)+sinv(2,3)*xyz(2,i)+sinv(3,3)*xyz(3,i)
          xss=ssx
          yss=ssy
          zss=ssz
          fracs(1,i)=xss
          fracs(2,i)=yss
          fracs(3,i)=zss
          end do
          end if





         open(unit=133,file='HISTORY.CELL')
         write(133,*)'---------------------------------'
          write(133,*)'optimization step ',itimes
          write(133,*)'cell is:'
         write(133,33)'A',x(1)*factor,x(2)*factor,x(3)*factor
         write(133,33)'B',x(4)*factor,x(5)*factor,x(6)*factor
         write(133,33)'C',x(7)*factor,x(8)*factor,x(9)*factor
         close(133)





         do i=1,3
         do j=1,3
         cell_copy(j,i)=cell(i,j)
         end do
         end do

       call cell_volume(cell_copy,volume)

 33               format(A2,f20.15,f20.15,f20.15)
         open(unit=133,file='HISTORY.CELL')
          do i=1,natoms
          xss=fracs(1,i)
          yss=fracs(2,i)
          zss=fracs(3,i)
          xyz(1,i)=x(1)*xss+x(4)*yss+x(7)*zss
          xyz(2,i)=x(2)*xss+x(5)*yss+x(8)*zss
          xyz(3,i)=x(3)*xss+x(6)*yss+x(9)*zss
          write(133,33)input_names(i),xyz(1,i)*factor,xyz(2,i)*factor,
     $         xyz(3,i)*factor
          end do
          close(133)

          close(761)
          eold=enew
       call get_stress(cell_copy,xyz, stress, enew,volume)
       open(unit=761,file='OUTPUT.DTPOLY',access='append')
          write(761,*)''
          write(761,*)'energy is',enew
          write(761,*)'delta E is',enew-eold
          write(761,*)''
          write(761,*)'Current Stress Tensor (GPa)'
          write(761,23)stress(1,1),stress(1,2),stress(1,3)
          write(761,23)stress(2,1),stress(2,2),stress(2,3)
          write(761,23)stress(3,1),stress(3,2),stress(3,3)
 23                     format(f25.15,f25.15,f25.15)

          press=(stress(1,1)+stress(2,2)+stress(3,3))/3.0d0
          write(761,*)'Current internal pressure',press



          do i=1,3
          do j=1,3
          stress(j,i)=-stress(j,i)+stress_in(j,i)
          stress_copy(j,i)=stress(j,i)
          end do
          end do


c compute strain energy
          estrold=estrain
          estrain=0.0
          do i=1,3
          do j=1,i
          estrain=estrain+delstr*stress(j,i)*volume/29421.0107637093d0
          end do
          end do
          write(761,*)''
          write(761,*)'Strain energy',estrain
          write(761,*)'delta Estrain is',estrain-estrold




c convert stress tensor to cell derivatives
        do i=1,3
        do j=1,3
          stress(i,j)=1.0d0*stress(i,j)/29421.0107637093d0
        end do
        end do

        call migs(cell,3,cinv,indx)

        call mmult(stress,cinv,smat,3)
        do i=1,3
        do j=1,3
        cell_deriv(j,i)=smat(j,i)*volume
        end do
        end do


        write(761,*)''
          write(761,*)'Current Cell Derivatives'
          write(761,23)cell_deriv(1,1),cell_deriv(1,2),cell_deriv(1,3)
          write(761,23)cell_deriv(2,1),cell_deriv(2,2),cell_deriv(2,3)
          write(761,23)cell_deriv(3,1),cell_deriv(3,2),cell_deriv(3,3)


c copy gradient to g vector
          rnorm=0.
         irow=0
        do i=1,3
        do j=1,3
           irow=irow+1
        g(irow)=cell_deriv(i,j)
        rnorm=rnorm+g(irow)**2
        end do
        end do


        write(761,*)'norm',dsqrt(rnorm)

        if(dsqrt(rnorm).lt.flowest)then
           flowest=dsqrt(rnorm)
           open(unit=34,file='best_point')
           write(34,*)flowest
           write(34,*)cell_copy
           write(34,*)xyz
           close(34)
        end if




        f=estrain
        f=enew
        if(usenrm)f=dsqrt(rnorm)
        write(761,*)'f returned =',f
        close(761)


        if(dsqrt(rnorm).lt.force)then

           open(unit=761,file='OUTPUT.DTPOLY',access='append')
           write(761,*)''
           write(761,*)''
           write(761,*)''
           write(761,*)''
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
           write(761,*)'    CELL OPTIMIZATION HAS CONVERGED  '
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'
           write(761,*)'*****************************************'

        call dump_final(cell_copy,volume,stress_copy)
           write(761,*)''
           write(761,*)''
           write(761,*)''
           write(761,*)''

        close(761)
        stop
        end if

c compute new fractionals for new cell above





         
        if(.not. refcell)then
           open(unit=761,file='OUTPUT.DTPOLY',access='append')
           write(761,*)'updating fractionals aka refcell is off'
           close(761)
          do i=1,natoms
          ssx=cinv(1,1)*xyz(1,i)+cinv(2,1)*xyz(2,i)+cinv(3,1)*xyz(3,i)
          ssy=cinv(1,2)*xyz(1,i)+cinv(2,2)*xyz(2,i)+cinv(3,2)*xyz(3,i)
          ssz=cinv(1,3)*xyz(1,i)+cinv(2,3)*xyz(2,i)+cinv(3,3)*xyz(3,i)
          xss=ssx
          yss=ssy
          zss=ssz
          fracs(1,i)=xss
          fracs(2,i)=yss
          fracs(3,i)=zss
          end do

          end if



      
      

      return
      end

