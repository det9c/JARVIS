subroutine j_tersoff(ax,ay,az,atmsym,nat,rxyz0, fxyz, ener, virial)
  use parameters_bks !module with parameters read from supp. input file 
  use control_module
  use tersoff_module 
  use constants 
  use cell_module
  use neighbors_module
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


rcut=cutoff

       fxyz=0.0d0
       ener=0.0d0
       virial=0.d0



          



       do i=1,nat

          beta=tersparams(5,terstype(i))
          d_n=tersparams(6,terstype(i))
          c=tersparams(7,terstype(i))
          d=tersparams(8,terstype(i))
          h=tersparams(9,terstype(i))
	  dxi=0.0d0
	  dyi=0.0d0
	  dzi=0.0d0
	  dedr=0.0d0

       do j1=1,num_neighbors(i)
          j=ipair_list(ifirst_neighbor(i)+j1-1)



         if(i.eq.j)goto 123
             rxij=rxyz0(1,i)-rxyz0(1,j)
             ryij=rxyz0(2,i)-rxyz0(2,j)
             rzij=rxyz0(3,i)-rxyz0(3,j)
         if(pbc)then
            if(cubic)then
             rxij=rxij-ax*anint(rxij/ax)
             ryij=ryij-ay*anint(ryij/ay)
             rzij=rzij-az*anint(rzij/az)
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
            end if 
          end if


!        print*,'i frac',fracs(1,1),fracs(2,1),fracs(3,1)
!        print*,'j frac',fracs(1,23),fracs(2,23),fracs(3,23)
!        print*,sxij,syij,szij



         rijsq=rxij*rxij + ryij*ryij + rzij*rzij
         rij=dsqrt(rijsq)
!         if(rij.gt.rcut)goto 123
! load parameters for easier coding
         indi=terstype(i)
         indj=terstype(j)
         ii=max(indi,indj)
         jj=min(indi,indj)
         ij=jj+offset(ii)
         A=ters_pair(3,ij)
         B=ters_pair(4,ij)
         dlam=ters_pair(1,ij)
         dmu=ters_pair(2,ij)
         Rmin=ters_pair(5,ij)
         Smin=ters_pair(6,ij)
         chi=ters_pair(7,ij)



! compute fc

         if(rij.lt.Rmin)then
          term_fc=1.0d0
          term1=0. ! derivative term used later
         elseif(Rmin .lt. rij .and. rij.lt.Smin)then
          t1= pi * (rij - Rmin) / (Smin - Rmin)
          term_fc=0.5d0 + 0.5d0 * cos(t1)

          term1=-0.5d0 * sin(t1)*pi/(Smin - Rmin)

         else
          goto 123
         end if 

         F_r=A*exp(-dlam*rij)
         F_a=-B*exp(-dmu*rij)

! compute zeta looping over k
         zeta=0.0d0
         zetapr=0.0d0
         virterm=0.

       do k1=1,num_neighbors(i)
          k=ipair_list(ifirst_neighbor(i)+k1-1)

         if(k.eq.i .or. k.eq.j)goto 130
         rxik=rxyz0(1,i)-rxyz0(1,k)
         ryik=rxyz0(2,i)-rxyz0(2,k)
         rzik=rxyz0(3,i)-rxyz0(3,k)

         if(pbc)then
            if(cubic)then
             rxik=rxik-ax*anint(rxik/ax)
             ryik=ryik-ay*anint(ryik/ay)
             rzik=rzik-az*anint(rzik/az)
            else
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

         riksq=rxik*rxik + ryik*ryik + rzik*rzik
         rik=dsqrt(riksq)
!	 if(rik.gt.rcut)goto 130
         rxjk=rxik-rxij
         ryjk=ryik-ryij
         rzjk=rzik-rzij
         rjksq=rxjk*rxjk + ryjk*ryjk + rzjk*rzjk
         rjk=dsqrt(rjksq)

         indk=terstype(k)
         ii=max(indi,indk)
         jj=min(indi,indk)
         ij=jj+offset(ii)
         Rmink=ters_pair(5,ij)
         Smink=ters_pair(6,ij)
         womega=ters_pair(8,ij)



          if(rik.lt.Rmink)then
          term_fck=womega
          term1k=0.
         elseif(Rmink .lt. rik .and. rik.lt.Smink)then
          t1= pi * (rik - Rmink) / (Smink- Rmink)
          term_fck=womega * (0.5d0 + 0.5d0 * cos(t1))
          term1k=-0.5d0*womega * sin(t1)*pi/(Smink - Rmink)
         else
          goto 130
         end if

        y     =  (-1.0D0*rjksq + rijsq + riksq)/(2.0D0*rij*rik)
        
!        if(y.ge.1.0d0)y=1.0d0
!        if(y.le.-1.0d0)y=-1.0d0
         if(y.ge.1.0d0) y= .9999999999999d0
         if(y.le.-1.0d0)y=-.9999999999999d0
        theta =  acos(y)

         gtheta=1.0d0 + c*c / (d*d) - c*c /(d*d+ (h-cos(theta))**2)
         zeta=zeta+term_fck*gtheta

         dxi=term1k*gtheta*rxik/rik
         dyi=term1k*gtheta*ryik/rik
         dzi=term1k*gtheta*rzik/rik


!         virterm=virterm-dxi*rxik-dyi*ryik-dzi*rzik

         virterm(1,1)=virterm(1,1)-dxi*rxik
         virterm(2,1)=virterm(2,1)-dxi*ryik
         virterm(3,1)=virterm(3,1)-dxi*rzik
         virterm(2,2)=virterm(2,2)-dyi*ryik
         virterm(3,2)=virterm(3,2)-dyi*rzik
         virterm(3,3)=virterm(3,3)-dzi*rzik







         zetapr(1,i)=zetapr(1,i)+dxi
         zetapr(2,i)=zetapr(2,i)+dyi
         zetapr(3,i)=zetapr(3,i)+dzi
         zetapr(1,k)=zetapr(1,k)-dxi
         zetapr(2,k)=zetapr(2,k)-dyi
         zetapr(3,k)=zetapr(3,k)-dzi


    !     if(y.eq.1.0d0 .or. y.eq.-1.0d0)goto 130
         call diff_theta(y,rij,rik,rjk,rijsq,riksq,rjksq, &
         rxij,ryij,rzij,rxik,ryik,rzik,   &
         rxjk,ryjk,rzjk,dtdxi,dtdyi,dtdzi, &
         dtdxj,dtdyj,dtdzj,dtdxk,dtdyk,dtdzk)
         

         dgdt=(c*c)*2.0d0*(h-cos(theta))*sin(theta)/( d*d+ (h-cos(theta))**2)**2
         dgdt=dgdt*term_fck
         zetapr(1,i)=zetapr(1,i)+dgdt*dtdxi
         zetapr(2,i)=zetapr(2,i)+dgdt*dtdyi
         zetapr(3,i)=zetapr(3,i)+dgdt*dtdzi
         zetapr(1,j)=zetapr(1,j)+dgdt*dtdxj
         zetapr(2,j)=zetapr(2,j)+dgdt*dtdyj
         zetapr(3,j)=zetapr(3,j)+dgdt*dtdzj
         zetapr(1,k)=zetapr(1,k)+dgdt*dtdxk
         zetapr(2,k)=zetapr(2,k)+dgdt*dtdyk
         zetapr(3,k)=zetapr(3,k)+dgdt*dtdzk



         virterm(1,1)=virterm(1,1)+(dgdt*dtdxj*rxij + dgdt*dtdxk*rxik)
         virterm(2,1)=virterm(2,1)+(dgdt*dtdxj*ryij + dgdt*dtdxk*ryik)
         virterm(3,1)=virterm(3,1)+(dgdt*dtdxj*rzij + dgdt*dtdxk*rzik)
         virterm(2,2)=virterm(2,2)+(dgdt*dtdyj*ryij + dgdt*dtdyk*ryik)
         virterm(3,2)=virterm(3,2)+(dgdt*dtdyj*rzij + dgdt*dtdyk*rzik)
         virterm(3,3)=virterm(3,3)+(dgdt*dtdzj*rzij + dgdt*dtdzk*rzik)










         
!         virterm=virterm+(dgdt*dtdxj*rxij+dgdt*dtdyj*ryij+dgdt*dtdzj*rzij&
!+dgdt*dtdxk*rxik+dgdt*dtdyk*ryik+dgdt*dtdzk*rzik)
!         virterm=virterm-dgdt*dtdxi*rxij-dgdt*dtdyi*ryij-dgdt*dtdzi*rzij
!virterm=virterm-dgdt*dtdxi*rxik-dgdt*dtdyi*ryik-dgdt*dtdzi*rzik





130      continue
         end do

         b=chi*(1.0d0 + (beta * zeta)**d_n)**(-0.5d0/d_n)



         ener=ener+term_fc*(F_r + b*F_a)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         get the derivative                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         term2=-dlam*F_r
         term4=-dmu*F_a

         dedr=term_fc*term2 + F_r*term1
         dedr=dedr+term1*b*F_a+term_fc*b*term4

         dxi=dedr*rxij/rij
         dyi=dedr*ryij/rij
         dzi=dedr*rzij/rij

!         virial=virial-dxi*rxij-dyi*ryij-dzi*rzij
         virial(1,1)=virial(1,1)-dxi*rxij
         virial(2,1)=virial(2,1)-dxi*ryij
         virial(3,1)=virial(3,1)-dxi*rzij
         virial(2,2)=virial(2,2)-dyi*ryij
         virial(3,2)=virial(3,2)-dyi*rzij
         virial(3,3)=virial(3,3)-dzi*rzij






         fxyz(1,i)=dxi+fxyz(1,i)
         fxyz(2,i)=dyi+fxyz(2,i)
         fxyz(3,i)=dzi+fxyz(3,i)
         fxyz(1,j)=fxyz(1,j)-dxi
         fxyz(2,j)=fxyz(2,j)-dyi
         fxyz(3,j)=fxyz(3,j)-dzi


! three body contribution
         if(zeta.ne.0)then
         dv=F_a*term_fc*(-0.5d0/d_n)*(1.0d0 + (beta * zeta)**d_n)**(-0.5d0/d_n - 1.0d0)
         dv=dv*chi*d_n * (beta**d_n) * zeta**(d_n-1.0)
         else
         dv=0.
         end if

         zetapr=zetapr*dv
         

         virterm=virterm*dv

         virial=virial+virterm 
         fxyz=fxyz+zetapr





 123   continue


       end do
       end do
       fxyz=-fxyz*0.5d0
       ener=ener*0.5d0
       virial=virial*.5d0
!      do i=1,nat
!       print*,'ters',fxyz(1,i),fxyz(2,i),fxyz(3,i)
!       end do
!       print*,'total',ener,ener*27.211


      open(unit=1,file='energy.out')
       write(1,*)ener
        close(1)




end subroutine j_tersoff


