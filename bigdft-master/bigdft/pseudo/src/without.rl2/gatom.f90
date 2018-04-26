!> @file
!! Generate atomic electronic configuration (without rl2)
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by matthias krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/
!!
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Updated version with spin polarization and libxc support
!! some notes about input variables:
!! nspol are the spin channels for the xc function.
!! nspin are spin components of orbitals with l > 0.
!! nspin may differ from nspol in the relativistic case.
!! the xc treatment and gaussian basis are initialized in the calling routine.
!! energ requests a total energy calculation.
!! verbose requests detailed output after the wfn is converged.
!! above two logical variables have been moved from a common block.
!! @note 
!!   A note about the r_r variable: For now, it does nothing, but it will
!!   be an experimental feature to test Gaussian type projectors with two
!!   different length scales r_l and r_r. The corresponding elements are
!!   named pp and qq here, ppr and qqr for the residues.
!!   There is no convention yet to read, pack and fit r_r(l).

      subroutine gatom(nspol,energ,verbose,  &
     &     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,  &
     &     occup,aeval,chrg,dhrg,ehrg,res,wght,wfnode,psir0,  &
     &     rcov,rprb,rcore,gcore,znuc,zion,rloc,gpot,r_l,hsep,  &
     &     vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,igrad,  &
     &     rr,rw,rd,ntime,itertot,etotal)


!          updated version with spin polarization and libXC support
!          ________________________________________________________


!          Some notes about input variables:

!          nspol are the spin channels for the XC function.
!          nspin are spin components of orbitals with l > 0.
!          nspin may differ from nspol in the relativistic case.
!          the XC treatment and gaussian basis are initialized in the calling routine.
!          energ requests a total energy calculation.
!          verbose requests detailed output after the wfn is converged.
!          Above two logical variables have been moved from a common block.

      
      implicit real*8 (a-h,o-z)
      dimension occup(noccmx,lmx,nsmx),  &
     &     aeval(noccmx,lmx,nsmx),chrg(noccmx,lmx,nsmx),  &
     &     dhrg(noccmx,lmx,nsmx),ehrg(noccmx,lmx,nsmx),  &
     &     res(noccmx,lmx,nsmx),wght(noccmx,lmx,nsmx,8),  &
     &     wfnode(noccmx,lmx,nsmx,3),  &
     &     gpot(4),r_l(lmx),hsep(6,lpmx,nsmx),  &
     &     gcore(4),  &
     &     vh(((ng+1)*(ng+2))/2,lcx+1,((ng+1)*(ng+2))/2,lmax+1),  &
     &     xp(0:ng), rmt(nint,((ng+1)*(ng+2))/2,lmax+1),  &
     &     rmtg(nint,((ng+1)*(ng+2))/2,lmax+1),  &
     &     ud(nint,((ng+1)*(ng+2))/2,lcx+1),  &
     &     psi(0:ngmx,noccmx,lmx,nsmx)
!     local arrays
!     notice nspin and nspol differ
      dimension  hh(0:ng,0:ng,lmax+1,nspin),ss(0:ng,0:ng,lmax+1),  &
     &     hht(0:ng,0:ng),sst(0:ng,0:ng),hhsc(((ng+1)*(ng+2))/2,lmax+1),  &
!          new work array for XC, used only in the polarized case
     &     hhxc(((ng+1)*(ng+2))/2,lmax+1,nspin),  &
     &     eval(0:ng),evec(0:ng,0:ng),pp1(0:ng,lpx+1),  &
     &     pp2(0:ng,lpx+1),pp3(0:ng,lpx+1),potgrd(nint),  &
     &     rho(((ng+1)*(ng+2))/2,lmax+1,nspol),  &
     &     rhoold(((ng+1)*(ng+2))/2,lmax+1,nspol),excgrd(nint),  &
     &     vxcgrd(nint,nspol),rr(nint),rw(nint),rd(nint),pexgrd(nint),  &
     &     ppr1(nint,lmax+1),ppr2(nint,lmax+1),ppr3(nint,lmax+1),  &
     &     aux1(nint),aux2(nint,0:ng,lmax+1),  &
     &     expxpr(0:ng,nint), tts(nspol)

      dimension rhogrd(nint,nspol),drhogrd(nint,nspol),  &
     &                             rhocore(nint,nspol)


!          two lines for experimental feature: Separable term with two r_l
      real(8):: r_r(lmx),qq1(0:ng,lpx+1),qq2(0:ng,lpx+1),qq3(0:ng,lpx+1), &
     &         qqr1(nint,lmax+1),qqr2(nint,lmax+1),qqr3(nint,lmax+1)  

      dimension y1(nint),y2(nint),y3(nint),  &
     &     rlist(0:nint),drlist(0:nint),ddrlist(0:nint)

       character*10 is(2)
      real*8 gamma
      logical energ,igrad, verbose, pol

      external gamma,wave,dwave,ddwave
      save nscf,nscfo,delta,odelta
      fourpi = 16.d0*atan(1.d0)
!      print*,'entered gatom, nspin and nspol are',nspin,nspol
!      is(1)= 'so=+0.5'
!      is(2)= 'so=-0.5'
!      print*,'rcov,rprb,zion,rloc,gpot'
!      print*,rcov,rprb,zion,rloc
!      print*,gpot(1)
!      print*,gpot(2)
!      print*,gpot(3)
!      print*,gpot(4)
!      print*,'lpx,lmax,lcx,noccmax,nspin'
!      print*,lpx,lmax,lcx,noccmax,nspin
!      print*,'ng,nint:',ng,nint
!            do l=0,lpx
!               write(6,*) 'l=',l
!               write(6,'(f7.3,t8,6e11.3,t76,a)') r_l(l+1),
!     :              (hsep(i,l+1,1),i=1,6),'r_l(),hsep(), '//is(1)
!               if (l.gt.0 .and. nspin.eq.2)
!     :              write(6,'(t8,6e11.3,t76,a)')
!     :              (hsep(i,l+1,2),i=1,6),'       hsep(), '//is(2)
!            enddo
!       print*,'xp',xp
!c----------------------------------------------end of debug output




      if (ntime .eq. 0) then
!     c. hartwig    modified scf mixing
!     initialize variables for first iteration rmix
         delta= 0.25d0
         odelta = 0.0d0
         nscf = 0
         nscfo=2**10
!     initial guess for wavefunction
         do l=0,lmax
            do iocc=1,noccmax
               do ispin=1,nspin
                  do i=0,ng
                     psi(i,iocc,l+1,ispin)=1.d-8
                  enddo
               enddo
            enddo
         enddo
      endif
      ntime=ntime+1


!***********************************************************************

!      just to be sure
       rhocore=0d0

       if(rcore>0d0)then
!      ***************************
!      *Nonlinear Core Correction*
!      ***************************

!      set up a simplest model charge distrubution
!      to take into account a frozen core charge for
!      the evaluation of exchange correlation terms. 

!      To be consistent with the convention in BigDFT
!      let us use a factor of four pi for the NLCC here.
       fourpi = 16.d0*atan(1.d0)

!      Careful: In BigDFT, the polynomial part is
!               NOT scaled by rcore. This results
!               easily in large coeffs and is not
!               optimal for fitting. Therefore,
!               we have slightly different conventions.
!               the corresponding transformation is
!               done once in pseudo for i/o of nlcc.
!
!               here we have powers of r/rcore, not of r.

          do k= 1,nint
              r2=(rr(k)/rcore)**2
              rhocore(k,1) =  &
     &        exp(-.5d0*r2) /fourpi  *(  &
     &        gcore(1)  &
     &      + gcore(2)*r2   &
     &      + gcore(3)*r2**2   &
     &      + gcore(4)*r2**3 )
!           write(17,*)rr(k),rhocore(k,1)
          end do

          if(nspol==2)then
!           split the charge equally among the two channels
!           even though the core charge should not be polarized,
!           it is stored in two spin channels for now.
            do k= 1,nint
               rhocore(k,1)=rhocore(k,1)*.5d0
               rhocore(k,2)=rhocore(k,1)
            end do
          end if
        end if

!***********************************************************************


! set up all quantities that are not density dependent
!
      do l=lcx+1,lmax
!        no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
         do ispin=1,max(min(2*l+1,nspin),nspol)
            if (occup(1,l+1,ispin).gt.1.d-10) stop 'lcx too small'
         enddo
      enddo

! projectors
      do l=0,lpx
!     gaussians
         gml1=sqrt( gamma(l+1.5d0) / (2.d0*r_l(l+1)**(2*l+3)) )  
         gml2=sqrt( gamma(l+3.5d0) / (2.d0*r_l(l+1)**(2*l+7)) )  &
     &        /(l+2.5d0)
         gml3=sqrt( gamma(l+5.5d0) / (2.d0*r_l(l+1)**(2*l+11)) )  &
     &        /((l+3.5d0)*(l+4.5d0))
         tt=1.d0/(2.d0*r_l(l+1)**2)
         do i=0,ng
            ttt=1.d0/(xp(i)+tt)
            pp1(i,l+1)=gml1*(sqrt(ttt)**(2*l+3))
            pp2(i,l+1)=gml2*ttt*(sqrt(ttt)**(2*l+3))
            pp3(i,l+1)=gml3*ttt**2*(sqrt(ttt)**(2*l+3))
         enddo
!     radial grid
        rnrm1=1.d0/sqrt(.5d0*gamma(l+1.5d0)*r_l(l+1)**(2*l+3))
        rnrm2=1.d0/sqrt(.5d0*gamma(l+3.5d0)*r_l(l+1)**(2*l+7))
        rnrm3=1.d0/sqrt(.5d0*gamma(l+5.5d0)*r_l(l+1)**(2*l+11))
        do k=1,nint
           r=rr(k)
           ppr1(k,l+1)=rnrm1*r**l*exp(-.5d0*(r/r_l(l+1))**2)
           ppr2(k,l+1)=rnrm2*r**(l+2)*exp(-.5d0*(r/r_l(l+1))**2)
           ppr3(k,l+1)=rnrm3*r**(l+4)*exp(-.5d0*(r/r_l(l+1))**2)
        enddo
      enddo

!     experimental feature: Generalize the seprable part to Gaussian
!     type projectors with a second length scale r_r 
!     -- not functional for now.

!     r_r is disabled:
      r_r = r_l

      do l=0,lpx
!        do not recompute qq when r_l is equal to r_r
!        its a waste of memory, though
         if( (r_r(l+1)-r_l(l+1))**2 <1d-6) then
            do i=0,ng
               qq1(i,l+1)=pp1(i,l+1)
               qq2(i,l+1)=pp2(i,l+1)
               qq3(i,l+1)=pp3(i,l+1)
            end do
            do k=1,nint
               qqr1(k,l+1)=ppr1(k,l+1)
               qqr2(k,l+1)=ppr2(k,l+1)
               qqr3(k,l+1)=ppr3(k,l+1)
            end do
            cycle ! next l
         end if
         
!     Gaussians
         gml1=sqrt( gamma(l+1.5d0) / (2.d0*r_r(l+1)**(2*l+3)) )
         gml2=sqrt( gamma(l+3.5d0) / (2.d0*r_r(l+1)**(2*l+7)) )  &
     &        /(l+2.5d0)
         gml3=sqrt( gamma(l+5.5d0) / (2.d0*r_r(l+1)**(2*l+11)) )  &
     &        /((l+3.5d0)*(l+4.5d0))
         tt=1.d0/(2.d0*r_r(l+1)**2)
         do i=0,ng
            ttt=1.d0/(xp(i)+tt)
            qq1(i,l+1)=gml1*(sqrt(ttt)**(2*l+3))
            qq2(i,l+1)=gml2*ttt*(sqrt(ttt)**(2*l+3))
            qq3(i,l+1)=gml3*ttt**2*(sqrt(ttt)**(2*l+3))
         enddo
!     radial grid
        rnrm1=1.d0/sqrt(.5d0*gamma(l+1.5d0)*r_l(l+1)**(2*l+3))
        rnrm2=1.d0/sqrt(.5d0*gamma(l+3.5d0)*r_l(l+1)**(2*l+7))
        rnrm3=1.d0/sqrt(.5d0*gamma(l+5.5d0)*r_l(l+1)**(2*l+11))
        do k=1,nint
           r=rr(k)
           qqr1(k,l+1)=rnrm1*r**l*exp(-.5d0*(r/r_r(l+1))**2)
           qqr2(k,l+1)=rnrm2*r**(l+2)*exp(-.5d0*(r/r_r(l+1))**2)
           qqr3(k,l+1)=rnrm3*r**(l+4)*exp(-.5d0*(r/r_r(l+1))**2)
        enddo
      enddo
!     end of experimental feature: r_r --> qq, qqr 



!   external potential on grid
      do k=1,nint
         r=rr(k)
         pexgrd(k)=.5d0*(r/rprb**2)**2-zion*Derf(r/(sqrt(2.d0)*rloc))/r  &
     &        + exp(-.5d0*(r/rloc)**2)*  &
     &        ( gpot(1) + gpot(2)*(r/rloc)**2 + gpot(3)*(r/rloc)**4 +  &
     &        gpot(4)*(r/rloc)**6 )
      enddo

!     store exp(-xp(i)*r**2) in expxpr()
      do k=1,nint
         r=rr(k)
         do i=0,ng
            expxpr(i,k)= exp(-xp(i)*r**2)
         enddo
      enddo

!     auxillary grids for resid:
      do k=1,nint
         r=rr(k)
         aux1(k)=fourpi/rw(k)
         do ll=0,lmax
            do i=0,ng
               aux2(k,i,ll+1)=(xp(i)*(3.d0+2.d0  &
     &              *ll-2.d0*xp(i)*r**2)*expxpr(i,k))
            enddo
         enddo
      enddo
!
! set up charge independent part of hamiltonian
!
      do l=0,lmax
!        no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
         do ispin=1,max(min(2*l+1,nspin),nspol)
            gml=0.5d0*gamma(l+0.5d0)
!

!     lower triangles only
!
            do j=0,ng
               do i=j,ng
                  hhij=0.0d0
                  d=xp(i)+xp(j)
                  sxp=1.d0/d
                  const=gml*sqrt(sxp)**(2*l+1)
!     overlap
                  ss(i,j,l+1)=const*sxp*(l+.5d0)
!     kinetic energy
                  hhij=.5d0*const*sxp**2*(3.d0*xp(i)*xp(j)+  &
     &                 l*(6.d0*xp(i)*xp(j)-xp(i)**2-xp(j)**2) -  &
     &                 l**2*(xp(i)-xp(j))**2  )+ .5d0*l*(l+1.d0)*const
!     potential energy from parabolic potential
                  hhij=hhij+  &
     &                 .5d0*const*sxp**2*(l+.5d0)*(l+1.5d0)/rprb**4
!     hartree potential from ionic core charge
                  tt=sqrt(1.d0+2.d0*rloc**2*d)
                  if (l.eq.0) then
                     hhij=hhij-zion/(2.d0*d*tt)
                  else if (l.eq.1) then
                     hhij=hhij-zion*  &
     &                    (1.d0 + 3.d0*rloc**2*d)/(2.d0*d**2*tt**3)
                  else if (l.eq.2) then
                     hhij=hhij-zion*  &
     &                    (2.d0+10.d0*rloc**2*d+15.d0*rloc**4*d**2)  &
     &                    /(2.d0*d**3*tt**5)
                  else if (l.eq.3) then
                     hhij=hhij-zion*3.d0*  &
     &                    (2.d0+14.d0*rloc**2*d+35.d0*rloc**4*d**2  &
     &                    +35.d0*rloc**6*d**3)/(2.d0*d**4*tt**7)
                  else
                     stop 'l too big'
                  endif
!     potential from repulsive gauss potential
                  tt=rloc**2/(.5d0+d*rloc**2)

                  pw1=1.5d0+dble(l)
                  pw2=2.5d0+dble(l)
                  pw3=3.5d0+dble(l)
                  pw4=4.5d0+dble(l)
                  hhij=hhij+gpot(1)*.5d0*gamma(pw1)*tt**pw1  &
     &                 + (gpot(2)/rloc**2)*.5d0*gamma(pw2)*tt**pw2  &
     &                 + (gpot(3)/rloc**4)*.5d0*gamma(pw3)*tt**pw3  &
     &                 + (gpot(4)/rloc**6)*.5d0*gamma(pw4)*tt**pw4
!     separabel terms
!     Experimental feature: Use two length scales r_l and r_r.
!     So we have two factors pp and qq.
!     For now, the two are identical, so this does not have an effect.
                  if (l.le.lpx) then
                     hhij = hhij  &
     &                    + pp1(i,l+1)*hsep(1,l+1,ispin)*qq1(j,l+1)  &
     &                    + pp1(i,l+1)*hsep(2,l+1,ispin)*qq2(j,l+1)  &
     &                    + pp2(i,l+1)*hsep(2,l+1,ispin)*qq1(j,l+1)  &
     &                    + pp2(i,l+1)*hsep(3,l+1,ispin)*qq2(j,l+1)  &
     &                    + pp1(i,l+1)*hsep(4,l+1,ispin)*qq3(j,l+1)  &
     &                    + pp3(i,l+1)*hsep(4,l+1,ispin)*qq1(j,l+1)  &
     &                    + pp2(i,l+1)*hsep(5,l+1,ispin)*qq3(j,l+1)  &
     &                    + pp3(i,l+1)*hsep(5,l+1,ispin)*qq2(j,l+1)  &
     &                    + pp3(i,l+1)*hsep(6,l+1,ispin)*qq3(j,l+1)
                  endif
                  hh(i,j,l+1,ispin)=hhij
               enddo
            enddo
         enddo
      enddo
!     hhxc is kept constant at zero in the unpolarized case
      hhxc=0d0

!
! finished setup of hh()
!


! initial charge and ev
      do l=0,lmax
         ij=0
         do j=0,ng
            do i=j,ng
               ij=ij+1
               rho(ij,l+1,:)=0.d0
            enddo
         enddo
      enddo
      evsum=1.d30



!ccccccccccccccccccccccccccccccc
!     begin the SCF cycles     c
!ccccccccccccccccccccccccccccccc


      do it=1,200
         evsumold=evsum
         evsum=0.d0
!
!     coefficients of charge density
!
!     the index for spin polarization shall be
         isp=1
!     in the unpolarized case, then
!     rho(i,l,1)  holds the total charge,
!     while in the polarized case nspol=2
!     rho(i,l,isp) are the two spin channels

         do l=0,lmax
            ij=0
            do j=0,ng
               do i=j,ng
                  ij=ij+1
                  rhoold(ij,l+1,:)=rho(ij,l+1,:)
                  rho(ij,l+1,:)=0.d0
               enddo
            enddo
         enddo
         do l=0,lmax
!           no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
            do ispin=1,max(min(2*l+1,nspin),nspol)
!              isp is always one in the unpolarized case
!              it determines which spin index of rho is addressed
               isp=min(ispin,nspol)
               do iocc=1,noccmax
                  if (occup(iocc,l+1,ispin).gt.1.d-10) then
                     ij=0
                     do j=0,ng
                        i=j
                        ij=ij+1
                        rho(ij,l+1,isp)=rho(ij,l+1,isp) +  &
     &                       psi(i,iocc,l+1,ispin)  &
     &                       *psi(j,iocc,l+1,ispin)  &
     &                       *occup(iocc,l+1,ispin)
                        do i=j+1,ng
                           ij=ij+1
                           rho(ij,l+1,isp)=rho(ij,l+1,isp) +  &
     &                          psi(i,iocc,l+1,ispin)  &
     &                          *psi(j,iocc,l+1,ispin)  &
     &                          *(2.d0*occup(iocc,l+1,ispin))
                        enddo
                     enddo
                  endif
               enddo
            enddo
         enddo

! determine optimal value for rmix
!     for minimum number of scf iterations
         if ( mod(ntime,20).eq.0 .and. it.eq.2) then
            tttt = delta
            if (nscf .lt. nscfo) then
               if (delta .gt. odelta) then
                  delta = delta + 0.05d0
               else
                  delta = delta - 0.05d0
               endif
            else
               if (delta .gt. odelta) then
                  delta = delta - 0.05d0
               else
                  delta = delta + 0.05d0
               endif
            endif
            delta = max(0.1d0,delta)
            delta = min(0.9d0,delta)
            odelta = tttt
            nscfo = nscf
            nscf = 0
         endif
!     F90 intrinsic
         call random_number(rmix)
         rmix = delta + (.5d0-delta/2.d0)*rmix
!     Intel (ifc)
!        rmix = delta + (.5d0-delta/2.d0)*dble(rand(0.0d0))
!     IBM/DEC/PGI
!        rmix = delta + (.5d0-delta/2.d0)*dble(rand())
!     CRAY
!        rmix = delta + (.5d0-delta/2.d0)*ranf()
!     rmix = delta
         if (it.eq.1) rmix=1.d0
         do l=0,lmax
            ij=0
            do j=0,ng
               do i=j,ng
                  ij=ij+1
!                 the : is over nspol components, 1 or 2
                  tts=rmix*rho(ij,l+1,:) + (1.d0-rmix)*rhoold(ij,l+1,:)
                  rho(ij,l+1,:)=tts
               enddo
            enddo
         enddo
!
!     calc. gradient only if xc-func. with gradient-corrections
!     rho on grid ij=1,nint:
!     rhogrd(k) =+ rmt(k,i,j,l+1)*rho(i,j,l+1)/(4*pi)
!     corresponding gradient:
!     drhogrd(k)=+ rmtg(k,i,j,l+1)*rho(i,j,l+1)/(4*pi)





         tt=1.d0/(16.d0*atan(1.d0))
         call DGEMV('N',nint,((ng+1)*(ng+2))/2*(lcx+1),  &
     &        tt,rmt,nint,rho(:,:,1),1,0.d0,rhogrd(:,1),1)
!        try yo keep it simple. Same procedure for spin down charge
         if(nspol==2)then
         call DGEMV('N',nint,((ng+1)*(ng+2))/2*(lcx+1),  &
     &        tt,rmt,nint,rho(:,:,2),1,0.d0,rhogrd(:,2),1)
         end if
!     for ggaenergy15, we don't need the gradient, as that driver
!     provides the derivative by finite differences on the radial grid.
!     Therefore, calculation of drhogrid is commented out and not
!     generalized to the spin polarized case.

!        if(igrad) call DGEMV('N',nint,((ng+1)*(ng+2))/2*(lcx+1),
!    &           tt,rmtg,nint,rho,1,0.d0,drhogrd,1)



!      The Nonlinear Core Correction (NLCC)
!      is added to the charge density
!      prior to calling the XC drivers
       if(rcore>0d0)then
           do k=1,nint
               rhogrd(k,1) = rhogrd(k,1) + rhocore(k,1)
           end do
       end if

       if(rcore>0d0.and.nspol==2)then
           do k=1,nint
               rhogrd(k,2) = rhogrd(k,2) + rhocore(k,2)
           end do
       end if


!     hutter
!      call evxc(nint,rr,rhogrd,drhogrd,vxcgrd,excgrd)
!     goedecker
!     libXC wrapper
!     call ggaenergy_15(nspol,nint,rw,rd,rhogrd,enexc,vxcgrd,excgrd)
!     call ggaenergy_15(nspol,nint,rr,rw,rd,rhogrd,enexc,vxcgrd,excgrd)
      call driveXC(nspol,nint,rr,rw,rd,rhogrd,enexc,vxcgrd,excgrd)
!        multiply with dr*r^2 to speed up calculation of matrix elements
!        open(11,file='rhogrd')
         do k=1,nint
            vxcgrd(k,:)=vxcgrd(k,:)*rw(k)/fourpi
         enddo
!        close(11)

!      IMPORTANT: since the real space representation
!      rhogrd will be used again only after recomputing
!      it from the gaussian representation rho,
!      there is no need to subtract the core charge here.

!      EXCEPTION: rhogrd will be passed to etot, 
!      which needs the core charge to calculate 
!      the EXC energy termi and then subtracts 
!      rhocore from rhogrd to compute the other
!      energy functionals (including VXC).

!
!     charge dependent part of hamiltonian
!     hartree potential from valence charge distribution

!     if nspol=2, the up and down potential differs due to Vxc 

!     do 4982,lp=0,lcx
!     do 4982,jp=0,ng
!     do 4982,ip=0,ng
!     hhsc(i,j,l+1) =+ vh(ip,jp,lp+1,i,j,l+1)*rho(ip,jp,lp+1)

         call DGEMV('T',(lcx+1)*((ng+1)*(ng+2))/2,  &
     &        (lmax+1)*((ng+1)*(ng+2))/2,1.d0,  &
     &        vh,(lcx+1)*((ng+1)*(ng+2))/2,rho(:,:,1),1,0.d0,hhsc,1)

!     if nspol=2, add the hartree potential from the 2nd charge channel
!     Note: It seems easier to add the charges first and then do one dgemv.
         if(nspol==2)   &
     &   call DGEMV('T',(lcx+1)*((ng+1)*(ng+2))/2,  &
     &        (lmax+1)*((ng+1)*(ng+2))/2,1.d0,  &
     &        vh,(lcx+1)*((ng+1)*(ng+2))/2,rho(:,:,2),1,1.d0,hhsc,1)
!                                                  ^    ^ 

!     potential from XC libraries 

!     do 8049,k=1,nint
!     8049 hhsc(i,j,l+1) =+ vxcgrd(k)*rmt(k,i,j,l+1)

!     MODIFICATION: if spin polarized, add this term to hhxc, not hhsc.
!                   hxc is a spin polarized matrix only for that purpose.
         hhxc=0d0
         if(nspol==1)then
              call DGEMV('T',nint,(lmax+1)*((ng+1)*(ng+2))/2,1.0d0,  &
     &        rmt,nint,vxcgrd(:,1),1,1.d0,hhsc,1)
         else
              call DGEMV('T',nint,(lmax+1)*((ng+1)*(ng+2))/2,1.0d0,  &
     &        rmt,nint,vxcgrd(:,1),1,0.d0,hhxc(:,:,1),1)
              call DGEMV('T',nint,(lmax+1)*((ng+1)*(ng+2))/2,1.0d0,  &
     &        rmt,nint,vxcgrd(:,2),1,1.d0,hhxc(:,:,2),1)

!        spin polarized XC term end
         end if


!     DIAGONALIZE
         do l=0,lmax
            do ispin=1,max(min(2*l+1,nspin), nspol)
!     LAPACK
               ij=0
               do j=0,ng
                  do i=j,ng
                     ij=ij+1
                     hht(i,j)=hh(i,j,l+1,ispin)+hhsc(ij,l+1)  &
     &                       +hhxc(ij,l+1,ispin) 
                     sst(i,j)=ss(i,j,l+1)
                  enddo
               enddo
!     IBM/DEC
               call DSYGV(1,'V','L',ng+1,hht,ng+1,sst,ng+1,  &
     &              eval,evec,(ng+1)**2,info)
!     the routine DSYGV is also included in sub_lapack.f
!     CRAY:
!     call SSYGV(1,'V','L',ng+1,hht,ng+1,sst,ng+1,
!     1       eval,evec,(ng+1)**2,info)
               if (info.ne.0) write(6,*) 'LAPACK',info
               do iocc=0,noccmax-1
                  do i=0,ng
                     evec(i,iocc)=hht(i,iocc)
                  enddo
               enddo
!     end LAPACK
               do iocc=1,noccmax
!                 write(6,*)'DEBUG: E(iocc,ispin,l,it)',
!    :                       eval(iocc-1),iocc-1,ispin,l,it
                  evsum=evsum+eval(iocc-1)
                  aeval(iocc,l+1,ispin)=eval(iocc-1)
                  do i=0,ng
                     psi(i,iocc,l+1,ispin)=evec(i,iocc-1)
                  enddo
               enddo
!     write(6,*) 'eval',l
!     55         format(5(e14.7))
!     write(6,55) eval
!     write(6,*) 'evec',l
!     do i=0,ng
!     33            format(10(e9.2))
!     write(6,33) (evec(i,iocc),iocc=0,noccmax-1)
!     enddo
           enddo
         enddo
         tt=abs(evsum-evsumold)
!        write(6,*)'DEBUG: residue=',tt,it
         if (tt.lt.1.d-8) goto 3000
      enddo
      write(6,*) 'WARNING: NO SC CONVERGENCE',tt
 3000 continue

!ccccccccccccccccccccccccccccccc
!    end of the SCF cycles     c
!ccccccccccccccccccccccccccccccc

!     write(*,*)'DEBUG: KS eigenvalues',aeval
      itertot=itertot+it
      nscf = nscf +it
      call resid(nspol,  &
     &     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,  &
     &     aeval,res,  &
     &     hsep,  &
     &     ud,nint,ng,ngmx,psi,rho,pp1,pp2,pp3,  &
     &     potgrd,pexgrd,vxcgrd,rr,rw,ppr1,ppr2,ppr3,aux1,aux2,  &
     &     qqr1,qqr2,qqr3,&
     &     expxpr)
!     etot evaluates Ehartree using rhogrd,
      if (energ) call etot(verbose,nspol,  &
     &     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,  &
     &     aeval,  &
     &     rprb,zion,rloc,gpot,r_l,hsep,  &
     &     xp,ud,nint,ng,ngmx,psi,rho,pp1,pp2,pp3,  &
     &     vxcgrd,excgrd,rhogrd,rhocore,occup,rr,rw,  &
     &     expxpr,etotal)
!
!     charge up to radius rcov or infinity
!

!     MODIFICATION: Can one integrate up to a different rcov for semicore states?
!     problem: loop over l is implicit here, one can not tell easily where nl<occmax(l)

      if (lmax.gt.3) stop 'cannot calculate chrg'
      do l=0,lmax
!        no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
         do ispin=1,max(min(2*l+1,nspin), nspol)
            do iocc=1,noccmax
               chrg(iocc,l+1,ispin)=0.d0
               dhrg(iocc,l+1,ispin)=0.d0
               ehrg(iocc,l+1,ispin)=0.d0
            enddo
         enddo
      enddo
      do ispin=1,max(min(2*l+1,nspin),nspol)
!     here, l=lmax+1, so do ispin=1,2 if lmax>0 and nspin=2
         do 3711,iocc=1,noccmax
!        loop over all nl(l)
            do 3762,j=0,ng
               do 3762,i=0,ng
                  d=xp(i)+xp(j)
                  sd=sqrt(d)
                  terf=Derf(sd*rcov)
                  texp=exp(-d*rcov**2)
                  tt0=0.4431134627263791d0*terf/sd**3-0.5d0*rcov*texp/d
                  tt1=0.6646701940895686d0*terf/sd**5 +  &  
     &                 (-0.75d0*rcov*texp - 0.5d0*d*rcov**3*texp)/d**2
                  chrg(iocc,1,ispin)=chrg(iocc,1,ispin) +  &
     &                 psi(i,iocc,1,ispin)*psi(j,iocc,1,ispin)*tt0
!     integrate up to rcov
!               dhrg(iocc,1,ispin)=dhrg(iocc,1,ispin) +
!     1              psi(i,iocc,1,ispin)*psi(j,iocc,1,ispin)*tt1
!     integrate up to inf
                  dhrg(iocc,1,ispin)=dhrg(iocc,1,ispin) +  &
     &                 psi(i,iocc,1,ispin)*psi(j,iocc,1,ispin)  &
     &                 *0.6646701940895686d0/sd**5
                  ehrg(iocc,1,ispin)=ehrg(iocc,1,ispin) +  &
     &                 psi(i,iocc,1,ispin)*psi(j,iocc,1,ispin)  &
     &                 *1.66167548522392d0/sd**7
                  if (lmax.eq.0) goto 3762

                  tt2=1.661675485223921d0*terf/sd**7 +  &
     &                 (-1.875d0*rcov*texp-1.25d0*d*rcov**3*texp-  &
     &                 0.5d0*d**2*rcov**5*texp)/d**3
                  chrg(iocc,2,ispin)=chrg(iocc,2,ispin) +  &
     &                 psi(i,iocc,2,ispin)*psi(j,iocc,2,ispin)*tt1
!     integrate up to rcov
!               dhrg(iocc,2,ispin)=dhrg(iocc,2,ispin) +
!     1              psi(i,iocc,2,ispin)*psi(j,iocc,2,ispin)*tt2
!     integrate up to inf
                  dhrg(iocc,2,ispin)=dhrg(iocc,2,ispin) +  &
     &                 psi(i,iocc,2,ispin)*psi(j,iocc,2,ispin)  &
     &                 *1.661675485223921d0/sd**7
                  ehrg(iocc,2,ispin)=ehrg(iocc,2,ispin) +  &
     &                 psi(i,iocc,2,ispin)*psi(j,iocc,2,ispin)  &
     &                 *5.815864198283725d0/sd**9
                  if (lmax.eq.1) goto 3762

                  tt3=5.815864198283725d0*terf/sd**9 +  &
     &                 (-6.5625d0*rcov*texp-4.375d0*d*rcov**3*texp-  &
     &                 1.75d0*d**2*rcov**5*texp -  &
     &                 0.5d0*d**3*rcov**7*texp)/d**4
                  chrg(iocc,3,ispin)=chrg(iocc,3,ispin) +  &
     &                 psi(i,iocc,3,ispin)*psi(j,iocc,3,ispin)*tt2
!     integrate up to rcov
!               dhrg(iocc,3,ispin)=dhrg(iocc,3,ispin) +
!     1              psi(i,iocc,3,ispin)*psi(j,iocc,3,ispin)*tt3
!     integrate up to inf
                  dhrg(iocc,3,ispin)=dhrg(iocc,3,ispin) +  &
     &                 psi(i,iocc,3,ispin)*psi(j,iocc,3,ispin)  &
     &                 *5.815864198283725d0/sd**9
                  ehrg(iocc,3,ispin)=ehrg(iocc,3,ispin) +  &
     &                 psi(i,iocc,3,ispin)*psi(j,iocc,3,ispin)  &
     &                 *26.17138889227676d0/sd**11

                  if (lmax.eq.2) goto 3762

                  chrg(iocc,4,ispin)=chrg(iocc,4,ispin) +  &
     &                 psi(i,iocc,4,ispin)*psi(j,iocc,4,ispin)*tt3
!     integrate up to rcov
!                  tt4=26.17138889227676d0*terf/sd**11+(-29.53125d0*
!     :                 rcov*texp-19.6875d0*d*rcov**3*texp-7.875d0*d**2
!     :                 *rcov**5*texp-2.25d0*d**3*rcov**7*texp-  &
!     &                 0.5d0*d**4*rcov**9*texp)/d**5
!               dhrg(iocc,4,ispin)=dhrg(iocc,4,ispin) +
!     1              psi(i,iocc,4,ispin)*psi(j,iocc,4,ispin)*tt4
!     integrate up to inf
                  dhrg(iocc,4,ispin)=dhrg(iocc,4,ispin) +  &
     &                 psi(i,iocc,4,ispin)*psi(j,iocc,4,ispin)  &
     &                 *26.17138889227676d0/sd**11
                  ehrg(iocc,4,ispin)=dhrg(iocc,4,ispin) +  &
     &                 psi(i,iocc,4,ispin)*psi(j,iocc,4,ispin)  &
     &                 *143.9426389075222d0/sd**13

 3762       continue
 3711    continue
      enddo

!
!     value at origin
!
      psir0=0.d0
      do i=0,ng
         psir0=psir0+psi(i,1,1,1)
      enddo
      psir0=psir0**2

!     node locations of psi*r
!     n-1 nodes allowed!
!     search only for nodes if the corresponding weights are <> zero
!     to avoid bumpy wavefunctions: no node of the first derivative of
!     the pseudowavefunction*r  and only one node
!     of the second derivative  up to the rmax of the lowest valence state
!
      tol =1.0d-12
! initialize ALL elements of wfnode to zero, also unused ones. 
! this is only to test whether this helps to get rid of a bug
      wfnode=0d0
      do l=0,lmax
!        no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
         do ispin=1, min(min(2*l+1,nspin),nspol)
            do nocc=1,noccmax
               if ( (wght(nocc,l+1,ispin,6).ne.0.d0)  &
     &              .or.  (wght(nocc,l+1,ispin,7).ne.0.d0)  &
     &              .or.  (wght(nocc,l+1,ispin,8).ne.0.d0) ) then

!       print*,'node search, l, nocc:',l,' ',nocc
                  wfnode(nocc,l+1,ispin,1)=0.d0
                  wfnode(nocc,l+1,ispin,2)=0.d0
                  wfnode(nocc,l+1,ispin,3)=0.d0
                  nnode=0
                  ndnode=0
                  nddnode=0
                  rnode=0.d0
                  rrnode=0.0d0
                  rrdnode=0.d0
                  dnode=0.d0
                  ddnode=0.0d0
!     find outer max of psi, search from ~10 bohr down
                  call detnp(nint,rr,10.0d0,kout)
                  ttrmax=rr(kout)
                  ra=ttrmax
                  ttmax= dabs(wave2(ng,l,psi(0,nocc,l+1,ispin),  &
     &                 expxpr,ra,kout,nint))
!      print*,'ttmax=',ttmax
                  do k=kout,1, -1
                     ra= rr(k)
                     ttpsi= dabs(wave2(ng,l,psi(0,nocc,l+1,ispin),  &
     &                 expxpr,ra,k,nint))
                     if ( ttpsi .gt. ttmax  &
     &                    .and. ttpsi .gt. 1.0d-4 ) then
                        ttmax=ttpsi
                        ttrmax=ra
                     endif
                  if (ttpsi.lt.ttmax .and. ttpsi.gt.1.0d-4) goto 3456
                  enddo
 3456             continue
!     search up to 90% of rmax
                  ttrmax=max(0.90d0*ttrmax,rr(1))
                  call detnp(nint,rr,ttrmax,kout)
                  ttrmax=rr(kout)
!       print*,'search up to ',ttrmax,ttmax
!     calc wavefunction and it's first two derivatives on the grid
!
                  do k=1,kout
                     call wave3(ng,l,xp,psi(0,nocc,l+1,ispin),  &
     &                    expxpr,rr(k),k,nint,y1(k),y2(k),y3(k))
                  enddo

                     do k = 2,kout
!     nodes of wavefunction
                     if (y1(k)*y1(k-1).lt.0.d0) then
                        nnode = nnode +1
                        x1=rr(k-1)
                        x2=rr(k)
                        rrnode = zbrent(wave,ng,ngmx,l,lmx,xp,psi,  &
     &                       nocc,noccmx,ispin,nsmx,  &
     &                       X1,X2,TOL)
                        if (nnode .ge.nocc) then
                           rnode=rnode+rrnode
!                          print*,'found rnode at:',rrnode
                        endif
                        rlist(nnode)=rrnode
                     endif
!     nodes of first derivative
                     if (y2(k)*y2(k-1).lt.0.d0) then
                        ndnode = ndnode +1
                        x1=rr(k-1)
                        x2=rr(k)
                        rrnode = zbrent(dwave,ng,ngmx,l,lmx,xp,psi,  &
     &                       nocc,noccmx,ispin,nsmx,  &
     &                       X1,X2,TOL)
                        if (ndnode .ge.nocc) then
                           dnode=dnode+rrnode
!                        print*,'found dnode at:',rrnode
                        endif
                        drlist(ndnode)=rrnode
                     endif
!     second derivative test:
                     if (y3(k)*y3(k-1).lt.0.d0) then
                        nddnode = nddnode + 1
                        x1=rr(k-1)
                        x2=rr(k)
                        rrnode = zbrent(ddwave,ng,ngmx,l,lmx,xp,psi,  &
     &                       nocc,noccmx,ispin,nsmx,  &
     &                       X1,X2,TOL)
!     only add the lowest node! (this one shoud dissapear)
                        if (nddnode .ge. nocc +1 ) then
                           ddnode = ddnode + rrdnode
!                          print*,'found ddnode at:',rrnode
                        else
                           rrdnode=rrnode
                        endif
                        ddrlist(nddnode)=rrnode
                     endif
                  enddo

!     print*,'rnode,dnode,ddnode',rnode,dnode,ddnode,nnode

!     new version: use integral of the relevant functions between the nodes
!     not the node-locations!
!     calc. necessary integrals:
                  sum1=0.0d0
                  sum2=0.0d0
                  sum3=0.0d0
!     rnodes:
                  do i=nnode+1-nocc,1,-2
                     aa=Wwav(ng,l,xp,psi(0,nocc,l+1,ispin),rlist(i))  &
     &                 -Wwav(ng,l,xp,psi(0,nocc,l+1,ispin),rlist(i-1))
                     sum1 = sum1+aa
                  enddo
!     dnodes
                  do i=ndnode+1-nocc,1,-2
                     aa=wave(ng,l,xp,psi(0,nocc,l+1,ispin),drlist(i))  &
     &                 -wave(ng,l,xp,psi(0,nocc,l+1,ispin),drlist(i-1))
                     sum2 = sum2+aa
                  enddo
!     ddnodes
                  do i=nddnode+1-nocc,1,-2
                     aa=dwave(ng,l,xp,psi(0,nocc,l+1,ispin),ddrlist(i))  &
     &                -dwave(ng,l,xp,psi(0,nocc,l+1,ispin),ddrlist(i-1))
!                    this test line is quite slow, for debuging purposes
                     sum3 = sum3+aa
                  enddo
!     old version for nodes as used in the paper:
!                  wfnode(nocc,l+1,ispin,1)=rnode
!                  wfnode(nocc,l+1,ispin,2)=dnode
!                  wfnode(nocc,l+1,ispin,3)=ddnode
!     new version, using the integrals of the function between the nodes
                  wfnode(nocc,l+1,ispin,1)=sum1
                  wfnode(nocc,l+1,ispin,2)=sum2
                  wfnode(nocc,l+1,ispin,3)=sum3

!                 Some lines for bugfixing of wfnode 
                  if(.not. sum1*sum1+sum2*sum2+sum3*sum3 >= 0d0)then
!                   let us use this condition as an isNaN function
!                   that should work with all fortran compilers.
!                   indeed NaN was observed sometimes for penalty terms from nodes
                    write(6,*)'Ouch! Some node integral is NaN for'
                    write(6,*)'l=',l,' s=',ispin,' n-ncore(l)=',nocc
                    write(6,*)'(NaN?)   node=',sum1,'   rnode=',rnode
                    write(6,*)'(NaN?)  dnode=',sum2,'  drnode=',drnode
                    write(6,*)'(NaN?) ddnode=',sum3,' ddrnode=',ddrnode
                    if(.not. sum1*sum1 >=0d0 )wfnode(nocc,l+1,ispin,1)=0d0
                    if(.not. sum2*sum2 >=0d0 )wfnode(nocc,l+1,ispin,2)=0d0
                    if(.not. sum3*sum3 >=0d0 )wfnode(nocc,l+1,ispin,3)=0d0
                  end if
               endif
            enddo
         enddo
      enddo

!     print*,'leave gatom'

      end

