subroutine angular(ax,ay,az,atmsym,nat,rxyz0, fxyz, ener, virial)
  use parameters_bks !module with parameters read from supp. input file 
  use control_module
  use tersoff_module 
  use constants 
  use cell_module
  use neighbors_module
  use atom_types_module
  implicit double precision (a-h,o-z)

interface
subroutine diff_theta(y,rij,rik,rjk,rijsq,riksq,rjksq, &
                      rxij,ryij,rzij,rxik,ryik,rzik,   &
                      rxjk,ryjk,rzjk,dtdxi,dtdyi,dtdzi, &
                      dtdxj,dtdyj,dtdzj,dtdxk,dtdyk,dtdzk)
double  precision,intent(in)::y,rij,rik,rjk,rijsq,riksq,rjksq, &
                      rxij,ryij,rzij,rxik,ryik,rzik,   &
                      rxjk,ryjk,rzjk
double  precision,intent(inout)::dtdxi,dtdyi,dtdzi, &
                      dtdxj,dtdyj,dtdzj,dtdxk,dtdyk,dtdzk
end subroutine diff_theta
end interface

  double precision,dimension(3,nat)::rxyz0,fxyz,zetapr
  integer::nat
  character(4),dimension(nat)::atmsym
  double precision::ener,ax,ay,az
   integer,dimension(nat)::ilabel
   double precision,dimension(3,3)::virial,virterm
  integer,save:: icounter=10
  integer,dimension(8)::icarbon

rcut=cutoff

       fxyz=0.0d0
       ener=0.0d0
       virial=0.d0




icarbon(1)=1
icarbon(2)=100
icarbon(3)=103
icarbon(4)=106
icarbon(5)=109
icarbon(6)=112
icarbon(7)=115
icarbon(8)=118


          
       do icell=1,ncell
          

        ibelow=(icell-1)*120  
          
      do jgroup=1,8
         i=ibelow+icarbon(jgroup)
         j=i+1
         k=i+2


             rxij=rxyz0(1,i)-rxyz0(1,j)
             ryij=rxyz0(2,i)-rxyz0(2,j)
             rzij=rxyz0(3,i)-rxyz0(3,j)
             rxik=rxyz0(1,i)-rxyz0(1,k)
             ryik=rxyz0(2,i)-rxyz0(2,k)
             rzik=rxyz0(3,i)-rxyz0(3,k)



         if(pbc)then
            if(cubic)then
             rxij=rxij-ax*anint(rxij/ax)
             ryij=ryij-ay*anint(ryij/ay)
             rzij=rzij-az*anint(rzij/az)
             rxik=rxik-ax*anint(rxik/ax)
             ryik=ryik-ay*anint(ryik/ay)
             rzik=rzik-az*anint(rzik/az)
            else
             sxij=fracs(1,i)-fracs(1,j)
             syij=fracs(2,i)-fracs(2,j)
             szij=fracs(3,i)-fracs(3,j)
             sxij=sxij-anint(sxij)
             syij=syij-anint(syij)
             szij=szij-anint(szij)
             rxij=cell(1,1)*sxij+cell(1,2)*syij+cell(1,3)*szij
             ryij=cell(2,1)*sxij+cell(2,2)*syij+cell(2,3)*szij
             rzij=cell(3,1)*sxij+cell(3,2)*syij+cell(3,3)*szij

             sxik=fracs(1,i)-fracs(1,k)
             syik=fracs(2,i)-fracs(2,k)
             szik=fracs(3,i)-fracs(3,k)
             sxik=sxik-anint(sxik)
             syik=syik-anint(syik)
             szik=szik-anint(szik)
             rxik=cell(1,1)*sxik+cell(1,2)*syik+cell(1,3)*szik
             ryik=cell(2,1)*sxik+cell(2,2)*syik+cell(2,3)*szik
             rzik=cell(3,1)*sxik+cell(3,2)*syik+cell(3,3)*szik

            end if 
          end if



         rijsq=rxij*rxij + ryij*ryij + rzij*rzij
         rij=dsqrt(rijsq)

         riksq=rxik*rxik + ryik*ryik + rzik*rzik
         rik=dsqrt(riksq)

         rxjk=rxik-rxij
         ryjk=ryik-ryij
         rzjk=rzik-rzij
         rjksq=rxjk*rxjk + ryjk*ryjk + rzjk*rzjk
         rjk=dsqrt(rjksq)

        y     =  (-1.0D0*rjksq + rijsq + riksq)/(2.0D0*rij*rik)
        
        if(y.ge.1.0d0)y=.9999999999d0
        if(y.le.-1.0d0)y=-.9999999999d0
        theta =  acos(y)

!        angmax=3.054326
!        if(theta.lt.angmax)then
!        print*,'bad angle set',strength*(y+1.0d0)**2,theta*180/pi
!        print*,'rij,rik,rjk',rij,rik,rjk
!        print*,'atoms i,j,k',i,j,k
!        call matprt(cell,3,3,3,3)
!        print*,'i frac',fracs(1,i),fracs(2,i),fracs(3,i)
!        print*,'j frac',fracs(1,j),fracs(2,j),fracs(3,j)
!        print*,'k frac',fracs(1,k),fracs(2,k),fracs(3,k)
!        print*,'sxij',sxij,syij,szij
!        print*,'sxik',sxik,syik,szik
!        end if

        ener=ener+strength*(y+1.0d0)**2


         call diff_theta(y,rij,rik,rjk,rijsq,riksq,rjksq, &
         rxij,ryij,rzij,rxik,ryik,rzik,   &
         rxjk,ryjk,rzjk,dtdxi,dtdyi,dtdzi, &
         dtdxj,dtdyj,dtdzj,dtdxk,dtdyk,dtdzk)

         dgdt=-2.0d0*strength*(y+1.0d0)*sin(theta)

         fxyz(1,i)=fxyz(1,i)+dgdt*dtdxi
         fxyz(2,i)=fxyz(2,i)+dgdt*dtdyi
         fxyz(3,i)=fxyz(3,i)+dgdt*dtdzi
         fxyz(1,j)=fxyz(1,j)+dgdt*dtdxj
         fxyz(2,j)=fxyz(2,j)+dgdt*dtdyj
         fxyz(3,j)=fxyz(3,j)+dgdt*dtdzj
         fxyz(1,k)=fxyz(1,k)+dgdt*dtdxk
         fxyz(2,k)=fxyz(2,k)+dgdt*dtdyk
         fxyz(3,k)=fxyz(3,k)+dgdt*dtdzk



         virial(1,1)=virial(1,1)+(dgdt*dtdxj*rxij + dgdt*dtdxk*rxik)
         virial(2,1)=virial(2,1)+(dgdt*dtdxj*ryij + dgdt*dtdxk*ryik)
         virial(3,1)=virial(3,1)+(dgdt*dtdxj*rzij + dgdt*dtdxk*rzik)
         virial(2,2)=virial(2,2)+(dgdt*dtdyj*ryij + dgdt*dtdyk*ryik)
         virial(3,2)=virial(3,2)+(dgdt*dtdyj*rzij + dgdt*dtdyk*rzik)
         virial(3,3)=virial(3,3)+(dgdt*dtdzj*rzij + dgdt*dtdzk*rzik)

         end do
         end do


         fxyz=-fxyz




end subroutine angular


