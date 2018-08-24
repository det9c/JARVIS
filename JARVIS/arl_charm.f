       subroutine arl_charm(ax,ay,az,atmsym,nat,rxyz0,fxyz,ener,pbc)
c----------------------------------------------------------------------
c     ax,ay,az: length of lattice vectors in simulation cell
c     Atmsym: char*2 vector of length nat with atom labels
c     Nat: number of atoms
c     Rxyz0: (3,nat) array of atomic coordinates
c     Forces: (3,nat) array of forces in eV/angstrom
c     Ener: potential energy in eV
c     Pbc: logical variable. True if pbc is on, false if no pbc
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      double precision l_j
      parameter(nsurf=9731, nbnd=12643, nang=24968, ndi=9947)
      character*2 atmsym(*)
      logical pbc
      common /mgc/ pbl(2,3),pba(2,5),pdh(3,4),
     .             eV_per_au,ang_per_au,deg_per_rad,ck_per_au,
     .             ibl(3,nbnd),iba(4,nang),idh(5,ndi)

      dimension rxyz0(*),fxyz(3,*),rxyz(3,nsurf)
      print*,'char'


      open(unit=1,file='bonds')
      read(1,*)ibl
      close(1)
      open(unit=1,file='angles')
      read(1,*)iba
      close(1)
      open(unit=1,file='dihedrals')
      read(1,*)idh
      close(1)





      do i=1, nsurf
         do j=1,3
           fxyz(j,i) = 0.0d0
           rxyz(j,i) = rxyz0((i-1)*3 + j) / ang_per_au  ! Starting Point Bohr
         enddo
      enddo

      ener = 0.0d0

c
c     Inter-molecular forces
c
c     Loop over the Individual Probe and Surface Bond Lengths
c

      
c carlos old code      do i= 1, nbnd
        cvrt = ang_per_au**2 / ck_per_au  ! Converts kcal/mol-A**2  to au/bohr**2
        do i= 1, nbnd
        fbl = pbl(1,ibl(3,i)) * cvrt          !  kb in Hartree/bohr**2

c carlos old line        fbl = pbl(1,ibl(3,i)) / ck_per_au        !  kb in Hartree/bohr**2
        r0  = pbl(2,ibl(3,i)) / ang_per_au       !  r0 in Bohr
        ener = ener + hookbl(ibl(1,i),ibl(2,i),rxyz,fxyz,fbl,r0)
      enddo

c
c     Loop over the Individual Probe and Surface Bond Angles
c

      do i= 1, nang
        fa = pba(1,iba(4,i)) / ck_per_au      !  ktheta in Hartree/rad**2
        a0 = pba(2,iba(4,i)) / deg_per_rad    !  theta0 in Radians
        ener= ener + hookang(iba(1,i),iba(2,i),iba(3,i),rxyz,fxyz,fa,a0) 
      enddo

c
c     Loop over the Individual Probe and Surface Dihedral Angles
c

      do i=1, ndi
        fd    = pdh(1,idh(5,i)) / ck_per_au    !  k-dihed in Hartree
        phi0  = pdh(2,idh(5,i)) / deg_per_rad  !  Phase shift in Radians
        cmult = pdh(3,idh(5,i))                !  Multiplicity of the Dihedral
        ener = ener + dihed(idh(1,i),idh(2,i),idh(3,i),idh(4,i),fd,phi0,
     .                                                  cmult,rxyz,fxyz)
      enddo

c
c     Array fxyz returns the Gradient of the Energy wrt the Nuclear Coords.
c
c     Convert back from au to eV and Ang.
c

      ener = ener * eV_per_au   ! Energy in eV 

      do i=1, nsurf
         do j=1,3
           fxyz(j,i) = fxyz(j,i) * eV_per_au / ang_per_au  ! Forces in eV/Ang.
         enddo
      enddo


      return
      end



      double precision function hookang(i,j,k,rxyz0,fxyz,ka,a0)
c
c     Hook's Law --- Bond Angle 
c
      implicit double precision (a-h,o-z)
      double precision ka
      dimension rxyz0(3,*),fxyz(3,*)
      hookang = 0.0d0

c     ka in units of au/rad**2
c     a0 in units of radians

      rij = drij(i,j,rxyz0,rxij,ryij,rzij)
      rik = drij(i,k,rxyz0,rxik,ryik,rzik)
      rjk = drij(j,k,rxyz0,rxjk,ryjk,rzjk)

      y = (rij**2 + rjk**2 - rik**2)/(2.0d0*rij*rjk) ! Cosine(aijk)  (Radians)
      aijk = dacos(y)                                ! Angle i--j--k (Radians)
      hookang = ka * (aijk - a0)**2                  ! Bend Energy   (au)

c
c     Three-body Forces
c
c     dV(hookang)/dq , q = x, y, z
c
      t1 = -ka * 2.0d0 * (aijk - a0)                       ! d[V(aijk)]/d[aijk]
      t2 = -1.0d0 / dsin(aijk)                             ! d[aijk]/dy
      t3 = 0.5d0*(1.0d0/rjk - rjk/rij**2 + rik**2/(rjk*rij**2)) ! dy/d[rij]
      t4 = -rik /(rij * rjk)                                    ! dy/d[rik]
      t5 = 0.5d0*(1.0d0/rij - rij/rjk**2 + rik**2/(rij*rjk**2)) ! dy/d[rjk]
      t1t2 = t1 * t2

      fxyz(1,i) = fxyz(1,i) -  t1t2 * (t3 * rxij + t4 * rxik)  ! -d[V(aijk)]/d[x(i)]
      fxyz(2,i) = fxyz(2,i) -  t1t2 * (t3 * ryij + t4 * ryik)  ! -d[V(aijk)]/d[y(i)]
      fxyz(3,i) = fxyz(3,i) -  t1t2 * (t3 * rzij + t4 * rzik)  ! -d[V(aijk)]/d[z(i)]

      fxyz(1,j) = fxyz(1,j) -  t1t2 * (t5 * rxjk - t3 * rxij)  ! -d[V(aijk)]/d[x(j)]
      fxyz(2,j) = fxyz(2,j) -  t1t2 * (t5 * ryjk - t3 * ryij)  ! -d[V(aijk)]/d[y(j)]
      fxyz(3,j) = fxyz(3,j) -  t1t2 * (t5 * rzjk - t3 * rzij)  ! -d[V(aijk)]/d[z(j)]

      fxyz(1,k) = fxyz(1,k) +  t1t2 * (t4 * rxik + t5 * rxjk)  ! -d[V(aijk)]/d[x(k)]
      fxyz(2,k) = fxyz(2,k) +  t1t2 * (t4 * ryik + t5 * ryjk)  ! -d[V(aijk)]/d[y(k)]
      fxyz(3,k) = fxyz(3,k) +  t1t2 * (t4 * rzik + t5 * rzjk)  ! -d[V(aijk)]/d[z(k)]

      return
      end



      double precision function hookbl(i,j,rxyz0,fxyz,kb,r0)
c
c     Hook's Law --- Bond Length
c
      implicit double precision (a-h,o-z)
      double precision kb
      dimension rxyz0(3,*),fxyz(3,*)

      hookbl = 0.0d0

c
c     E(hookbl) = kb(r[Ni]-r[Hj])**2
c

c     kb in units of au/Bohr**2
c     r0 in units of Bohr
      rij = drij(i,j,rxyz0,rxij,ryij,rzij)

      hookbl = kb * (rij - r0)**2

c
c     Two-body Forces 
c
c     dV(hookbl)/dq , q = x, y, z
c

      y = kb * 2.0d0 * (rij - r0)    ! dV[(hookbl)]/d[rij]

      prod = rxij * y
      fxyz(1,j) = fxyz(1,j) - prod
      fxyz(1,i) = fxyz(1,i) + prod

      prod = ryij * y
      fxyz(2,j) = fxyz(2,j) - prod
      fxyz(2,i) = fxyz(2,i) + prod

      prod = rzij * y
      fxyz(3,j) = fxyz(3,j) - prod
      fxyz(3,i) = fxyz(3,i) + prod

      return
      end



      double precision function dihed(i,j,k,l,kt,phi,cmult,rxyz0,fxyz)
c
c     CHARMM Dihedral Potential: kt(1+cos(mult*theta - phi))
c     kt = force constant --- input in kcal
c     mult = multiplicity --- input as an integer
c     phi = phase factor --- input in degrees
c
      implicit double precision (a-h,o-z)
      double precision kt
      dimension rxyz0(3,*),fxyz(3,*),vrijxkj(3),vrkjxkl(3)

      dihed = 0.0d0

c     kt  in units of au
c     phi in radians

      rkj = drij(k,j,rxyz0,rx23,ry23,rz23)
      call crx(i,j,k,j,rxyz0,vrijxkj,rijkj)
      call crx(k,j,k,l,rxyz0,vrkjxkl,rkjkl)
      call cdot(vrijxkj,vrkjxkl,dp)

      y = dp /(rijkj * rkjkl) ! = Cosine theta
       theta = sgn(vrijxkj,vrkjxkl,rxyz0,j,k) * dacos(y)


      dihed = kt * (1.0d0 + dcos(cmult*theta - phi)) ! Dihedral (torsion) energy

c
c     four-body Forces
c
c     dV(dihed)/d(q) , q = x, y, z
c

      dpijkj = fdot(i,j,k,j,rxyz0)
      dpklkj = fdot(k,l,k,j,rxyz0)
      dpijkj = fdot(i,j,k,j,rxyz0)
      dpklkj = fdot(k,l,k,j,rxyz0)

      t1   = -cmult * kt * dsin(cmult*theta - phi) ! dV(di-hed)/d(theta)
      t2   =  rkj/rijkj**2
      t3   = -rkj/rkjkl**2
      t4   = dpijkj/rkj**2
      t4mo = t4 - 1.0d0
      t5   = dpklkj/rkj**2
      t5mo = t5 - 1.0d0

      dydxi = t2 * vrijxkj(1)             ! d[arccos(y)]/d[xi]
      dydyi = t2 * vrijxkj(2)             ! d[arccos(y)]/d[yi]
      dydzi = t2 * vrijxkj(3)             ! d[arccos(y)]/d[zi]

      dydxl = t3 * vrkjxkl(1)             ! d[arccos(y)]/d[xl]
      dydyl = t3 * vrkjxkl(2)             ! d[arccos(y)]/d[yl]
      dydzl = t3 * vrkjxkl(3)             ! d[arccos(y)]/d[zl]

      dydxj = t4mo * dydxi - t5 * dydxl   ! d[arccos(y)]/d[xj]
      dydyj = t4mo * dydyi - t5 * dydyl   ! d[arccos(y)]/d[yj]
      dydzj = t4mo * dydzi - t5 * dydzl   ! d[arccos(y)]/d[zj]

      dydxk = t5mo * dydxl - t4 * dydxi   ! d[arccos(y)]/d[xk]
      dydyk = t5mo * dydyl - t4 * dydyi   ! d[arccos(y)]/d[xk]
      dydzk = t5mo * dydzl - t4 * dydzi   ! d[arccos(y)]/d[xk]


      fxyz(1,i) = fxyz(1,i) + t1 * dydxi ! dV(theta)/d(ri[x])
      fxyz(2,i) = fxyz(2,i) + t1 * dydyi ! dV(theta)/d(ri[y])
      fxyz(3,i) = fxyz(3,i) + t1 * dydzi ! dV(theta)/d(ri[z])

      fxyz(1,j) = fxyz(1,j) + t1 * dydxj ! dV(theta)/d(rj[x])
      fxyz(2,j) = fxyz(2,j) + t1 * dydyj ! dV(theta)/d(rj[y])
      fxyz(3,j) = fxyz(3,j) + t1 * dydzj ! dV(theta)/d(rj[z])

      fxyz(1,k) = fxyz(1,k) + t1 * dydxk ! dV(theta)/d(rk[x])
      fxyz(2,k) = fxyz(2,k) + t1 * dydyk ! dV(theta)/d(rk[y])
      fxyz(3,k) = fxyz(3,k) + t1 * dydzk ! dV(theta)/d(rk[z])

      fxyz(1,l) = fxyz(1,l) + t1 * dydxl ! dV(theta)/d(rl[x])
      fxyz(2,l) = fxyz(2,l) + t1 * dydyl ! dV(theta)/d(rl[y])
      fxyz(3,l) = fxyz(3,l) + t1 * dydzl ! dV(theta)/d(rl[z])

      return
      end



      double precision function drij(i,j,rxyz0,rx,ry,rz)
      implicit double precision (a-h,o-z)
      dimension rxyz0(3,*)

      rx = rxyz0(1,i) - rxyz0(1,j)
      ry = rxyz0(2,i) - rxyz0(2,j)
      rz = rxyz0(3,i) - rxyz0(3,j)

      drij = dsqrt(rx**2 + ry**2 + rz**2)   ! rij

      rx = rx / drij        ! d[rij]/dx(i) = -d[rij]/dx(j)
      ry = ry / drij        ! d[rij]/dy(i) = -d[rij]/dy(j)
      rz = rz / drij        ! d[rij]/dz(i) = -d[rij]/dz(j)

      return
      end


      subroutine crx(i,j,k,l,r0,xp,d)
      implicit double precision (a-h,o-z)
      dimension r0(3,*),xp(3),rij(3),rkl(3)


      rij(1) = r0(1,i) - r0(1,j)
      rij(2) = r0(2,i) - r0(2,j)
      rij(3) = r0(3,i) - r0(3,j)

      rkl(1) = r0(1,k) - r0(1,l)
      rkl(2) = r0(2,k) - r0(2,l)
      rkl(3) = r0(3,k) - r0(3,l)

      xp(1) = rij(2) * rkl(3) - rij(3) * rkl(2)
      xp(2) = rij(3) * rkl(1) - rij(1) * rkl(3)
      xp(3) = rij(1) * rkl(2) - rij(2) * rkl(1)

      d = dsqrt(xp(1)**2 + xp(2)**2 + xp(3)**2)

      return
      end


      subroutine cdot(v1,v2,dp)
      implicit double precision (a-h,o-z)
      dimension v1(*),v2(*)

      dp = v1(1)* v2(1) + v1(2)* v2(2) + v1(3)* v2(3)

      return
      end


      double precision function fdot(i,j,k,l,r0)
      implicit double precision (a-h,o-z)
      dimension r0(3,*),v1(3),v2(3)

      v1(1) = r0(1,i) - r0(1,j)
      v1(2) = r0(2,i) - r0(2,j)
      v1(3) = r0(3,i) - r0(3,j)

      v2(1) = r0(1,k) -r0(1,l)
      v2(2) = r0(2,k) -r0(2,l)
      v2(3) = r0(3,k) -r0(3,l)

      fdot = v1(1)* v2(1) + v1(2)* v2(2) + v1(3)* v2(3)

      return
      end


      double precision function sgn(v1,v2,r0,j,i)
      implicit double precision (a-h,o-z)
      dimension v1(*),v2(*),r0(3,*),vp(3),rij(3)

      rij(1) = r0(1,i) - r0(1,j)
      rij(2) = r0(2,i) - r0(2,j)
      rij(3) = r0(3,i) - r0(3,j)

      vp(1) = v1(2) * v2(3) - v1(3) * v2(2)
      vp(2) = v1(3) * v2(1) - v1(1) * v2(3)
      vp(3) = v1(1) * v2(2) - v1(2) * v2(1)

      call cdot(rij,vp,rp)

      sgn = dsign(1.0d0,rp)

      return
      end


      block data
      implicit double precision(a-h,o-z)
      parameter(nsurf=9731, nbnd=12643, nang=24968, ndi=9947)
      common /mgc/ pbl(2,3),pba(2,5),pdh(3,4),
     .             eV_per_au,ang_per_au,deg_per_rad,ck_per_au,
     .             ibl(3,nbnd),iba(4,nang),idh(5,ndi)

      data ck_per_au   /627.5095d0/        ! (kcal/mol)/Hartree
      data eV_per_au   /27.2116d0/         ! eV/Hartree
      data ang_per_au  /0.529177249d0/     ! Angstrom/Bohr
      data deg_per_rad /57.29577951d0/     ! Degrees per Radian


      data pbl /    ! paired by type: Kbl(kcal/mol-A**2), b0(A), ...
     .       566.0d0, 0.9750d0,  325.0d0, 1.680d0,  302.0d0, 1.698d0 /

      data pba /    ! paired by type:  Ktheta(kcal/rad**2), theta0(deg(, ...
     .       34.0d0, 150.0d0, 32.0d0, 126.0d0, 30.0d0, 121.5d0, 34.0d0,
     .      122.5d0, 30.0d0, 117.0d0 /

      data pdh / !ensemble by type: Kphi(kcal/deg),phi0(phase-deg),mult ...
     .      0.35d0, 0.0d0, 3.0d0,  0.35d0, 0.0d0, 3.0d0,  0.18d0, 0.0d0,
     .      5.0d0,  0.12d0, 0.0d0, 5.0d0 /


      end

