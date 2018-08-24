subroutine diff_theta(y,rij,rik,rjk,rijsq,riksq,rjksq, &
                      rxij,ryij,rzij,rxik,ryik,rzik,   &
                      rxjk,ryjk,rzjk,dtdxi,dtdyi,dtdzi, &
                      dtdxj,dtdyj,dtdzj,dtdxk,dtdyk,dtdzk)
implicit double precision (a-h,o-z) 
double  precision,intent(in)::y,rij,rik,rjk,rijsq,riksq,rjksq, &
                      rxij,ryij,rzij,rxik,ryik,rzik,   &
                      rxjk,ryjk,rzjk
double  precision,intent(inout)::dtdxi,dtdyi,dtdzi, &
                      dtdxj,dtdyj,dtdzj,dtdxk,dtdyk,dtdzk 



      dtdxi=-1.0D0/sqrt(1.0D0-y**2)
      term4=(1.0D0/(rik*rij))*((rxij)+(rxik))
      term5=-0.5D0*(rjksq-riksq-rijsq)
      term6=-1.0D0*(rxij)/(rik*(rij**3))
      term7=-1.0D0*(rxik)/(rij*(rik**3))
      term8=term5*(term6+term7)
      term4=term4+term8
      dtdxi=dtdxi*term4

      dtdyi=-1.0D0/sqrt(1.0D0-y**2)
      term4=(1.0D0/(rik*rij))*((ryij)+(ryik))
      term5=-0.5D0*(rjksq-riksq-rijsq)
      term6=-1.0D0*(ryij)/(rik*(rij**3))
      term7=-1.0D0*(ryik)/(rij*(rik**3))
      term8=term5*(term6+term7)
      term4=term4+term8
      dtdyi=dtdyi*term4

      dtdzi=-1.0D0/sqrt(1.0D0-y**2)
      term4=(1.0D0/(rik*rij))*((rzij)+(rzik))
      term5=-0.5D0*(rjksq-riksq-rijsq)
      term6=-1.0D0*(rzij)/(rik*(rij**3))
      term7=-1.0D0*(rzik)/(rij*(rik**3))
      term8=term5*(term6+term7)
      term4=term4+term8
      dtdzi=dtdzi*term4

! atom j

      dtdxj=-1.0D0/sqrt(1.0D0-y**2)
      term4=(1.0D0/(rik*rij))*(-1.0d0*(rxjk)+(-1.*rxij))
      term5=-0.5D0*(rjksq-riksq-rijsq)
      term6=-1.0D0*(-1.*rxij)/(rik*(rij**3))
      term8=term5*(term6)
      term4=term4+term8
      dtdxj=dtdxj*term4

      dtdyj=-1.0D0/sqrt(1.0D0-y**2)
      term4=(1.0D0/(rik*rij))*(-1.0d0*(ryjk)+(-1.*ryij))
      term5=-0.5D0*(rjksq-riksq-rijsq)
      term6=-1.0D0*(-1.*ryij)/(rik*(rij**3))
      term8=term5*(term6)
      term4=term4+term8
      dtdyj=dtdyj*term4

      dtdzj=-1.0D0/sqrt(1.0D0-y**2)
      term4=(1.0D0/(rik*rij))*(-1.0d0*(rzjk)+(-1.*rzij))
      term5=-0.5D0*(rjksq-riksq-rijsq)
      term6=-1.0D0*(-1.*rzij)/(rik*(rij**3))
      term8=term5*(term6)
      term4=term4+term8
      dtdzj=dtdzj*term4

! atom k

      dtdxk=-1.0D0/sqrt(1.0D0-y**2)
      term4=(1.0D0/(rik*rij))*(-1.0d0*(-1.*rxjk)+(-1.*rxik))
      term5=-0.5D0*(rjksq-riksq-rijsq)
      term6=-1.0D0*(-1.*rxik)/(rij*(rik**3))
      term8=term5*(term6)
      term4=term4+term8
      dtdxk=dtdxk*term4

      dtdyk=-1.0D0/sqrt(1.0D0-y**2)                           
      term4=(1.0D0/(rik*rij))*(-1.0d0*(-1.*ryjk)+(-1.*ryik))  
      term5=-0.5D0*(rjksq-riksq-rijsq)                        
      term6=-1.0D0*(-1.*ryik)/(rij*(rik**3))                  
      term8=term5*(term6)                                     
      term4=term4+term8                                       
      dtdyk=dtdyk*term4        


      dtdzk=-1.0D0/sqrt(1.0D0-y**2)                           
      term4=(1.0D0/(rik*rij))*(-1.0d0*(-1.*rzjk)+(-1.*rzik))  
      term5=-0.5D0*(rjksq-riksq-rijsq)                        
      term6=-1.0D0*(-1.*rzik)/(rij*(rik**3))                  
      term8=term5*(term6)                                     
      term4=term4+term8                                       
      dtdzk=dtdzk*term4        



return
end subroutine diff_theta 

