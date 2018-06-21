!> @file
!!  Routines for the ionic energy contribution
!! @author
!!    Copyright (C) 2007-2015 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Calculate the ionic contribution to the energy and the forces
subroutine IonicEnergyandForces(iproc,nproc,dpbox,at,elecfield,&
     & rxyz,eion,fion,dispersion,edisp,fdisp,ewaldstr,&
     & pot_ion,pkernel,psoffset)
  use dynamic_memory
  use dictionaries
  use wrapper_mpi
  use f_utils
  use module_defs
  use module_base, only: bigdft_mpi
  use numerics
  use f_enums
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use multipole_preserving
  use module_atoms
  use module_dpbox
  use abi_interfaces_geometry, only: abi_metric
  use abi_interfaces_common, only: abi_ewald, abi_ewald2
  use m_paw_numeric, only: paw_splint
  use abi_interfaces_numeric, only: abi_derf_ab
  use vdwcorrection
  use yaml_output
  use public_enums, only: PSPCODE_PAW, PSPCODE_GTH, PSPCODE_HGH, PSPCODE_HGH_K, PSPCODE_HGH_K_NLCC
  use bounds, only: ext_buffers
  use sparsematrix_init, only: distribute_on_tasks
  use gaussians
  use box
  use bounds
  implicit none
  !Arguments
  type(denspot_distribution), intent(inout) :: dpbox
  type(atoms_data), intent(in) :: at
  integer, intent(in) :: iproc,nproc,dispersion
  real(gp), dimension(3), intent(in) :: elecfield
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  type(coulomb_operator), intent(inout) :: pkernel
  real(gp), intent(out) :: eion,edisp,psoffset
  real(dp), dimension(6),intent(out) :: ewaldstr
  real(gp), dimension(:,:), pointer :: fion,fdisp
  real(dp), dimension(*), intent(out) :: pot_ion
  !Local variables
  real(gp), parameter :: mp_tiny = 1.e-30_gp
  logical :: slowion=.false.!,use_iterator=.false.
  logical :: perx,pery,perz
  integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,n3i,n3pi,i3s
  integer :: n1i,n2i,i,iat,ii,ityp,jat,jtyp,natp,isat,iiat
  integer :: mpnx,mpny,mpnz
  real(gp) :: ucvol,rloc
  real(gp) :: twopitothreehalf,atint,shortlength,charge,eself,rx,ry,rz
  real(gp) :: fxion,fyion,fzion,dist,fxerf,fyerf,fzerf
  real(gp) :: hxh,hyh,hzh
  real(gp) :: chgprod
  real(gp) :: xp,Vel,prefactor,ehart,de
  !real(gp) :: Mz,cmassy
  real(gp), dimension(3,3) :: gmet,rmet,rprimd,gprimd
  !other arrays for the ewald treatment
  real(gp), dimension(:,:), allocatable :: fewald,xred
  real(gp), dimension(:), allocatable  :: mpx,mpy,mpz
  real(gp), dimension(3) :: cc
  type(atoms_iterator) :: atit
  type(gaussian_real_space) :: g

  call f_routine(id='IonicEnergyandForces')
  call timing(iproc,'ionic_energy','ON')
  fion = f_malloc_ptr((/ 3, at%astruct%nat /),id='fion')
  fdisp = f_malloc_ptr((/ 3, at%astruct%nat /),id='fdisp')

  !Initialize the work arrays needed to integrate with isf
  !Determine the number of points depending on the min rloc
  if (at%multipole_preserving) &
     call initialize_real_space_conversion(isf_m=at%mp_isf,rlocs=at%psppar(0,0,:))

  ! Aliasing
  hxh = dpbox%mesh%hgrids(1)
  hyh = dpbox%mesh%hgrids(2)
  hzh = dpbox%mesh%hgrids(3)
  n1i = dpbox%mesh%ndims(1)
  n2i = dpbox%mesh%ndims(2)
  n3i = dpbox%mesh%ndims(3)
  i3s = dpbox%i3s + dpbox%i3xcsh
  n3pi = dpbox%n3pi

  psoffset=0.0_gp
  ewaldstr=0.0_gp
  if (at%astruct%geocode == 'P') then
     !here we insert the calculation of the ewald forces
     fewald = f_malloc((/ 3, at%astruct%nat /),id='fewald')
     xred = f_malloc((/ 3, at%astruct%nat /),id='xred')

     !calculate rprimd
     rprimd(:,:)=0.0_gp

     rprimd(1,1)=at%astruct%cell_dim(1)
     rprimd(2,2)=at%astruct%cell_dim(2)
     rprimd(3,3)=at%astruct%cell_dim(3)

     !calculate the metrics and the volume
     call abi_metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

     !calculate reduced coordinates
     !$omp parallel if (at%astruct%nat>1000) &
     !$omp default(none) shared(at, xred, gprimd, rxyz) private(iat, ii)
     !$omp do schedule(static)
     do iat=1,at%astruct%nat
        do ii=1,3
           xred(ii,iat)= gprimd(1,ii)*rxyz(1,iat)+gprimd(2,ii)*rxyz(2,iat)+&
                gprimd(3,ii)*rxyz(3,iat)
        end do
     end do
     !$omp end do
     !$omp end parallel

     !calculate ewald energy and forces + stress
     call abi_ewald(iproc,nproc,bigdft_mpi%mpi_comm,&
          eion,gmet,fewald,at%astruct%nat,at%astruct%ntypes,rmet,at%astruct%iatype,ucvol,&
          xred,real(at%nelpsp,gp))
     ewaldstr=0.0_dp
     call abi_ewald2(iproc,nproc,bigdft_mpi%mpi_comm,&
          gmet,at%astruct%nat,at%astruct%ntypes,rmet,rprimd,ewaldstr,at%astruct%iatype,&
          ucvol,xred,real(at%nelpsp,gp))

     ! our sequence of strten elements : 11 22 33 12 13 23
     ! abinit output                   : 11 22 33 23 13 12

     !make forces dimensional
     !$omp parallel if (at%astruct%nat>1000) &
     !$omp default(none) shared(at, fion, gprimd, fewald) private(iat, ii)
     !$omp do schedule(static)
     do iat=1,at%astruct%nat
        do ii=1,3
           fion(ii,iat)= - (gprimd(ii,1)*fewald(1,iat)+&
                gprimd(ii,2)*fewald(2,iat)+&
                gprimd(ii,3)*fewald(3,iat))
        end do
        !if (nproc==1 .and. slowion) print *,'iat,fion',iat,(fion(j1,iat),j1=1,3)
     end do
     !$omp end do
     !$omp end parallel

     call f_free(xred)
     call f_free(fewald)

     !now calculate the integral of the local psp
     !this is the offset to be applied in the Poisson Solver to have a neutralizing background
     psoffset=0.0_gp
     shortlength=0.0_gp
     charge=0.0_gp
     twopitothreehalf=2.0_gp*pi*sqrt(2.0_gp*pi)
!!$     !$omp parallel if (at%astruct%nat>1000) &
!!$     !$omp default(none) shared(at, psoffset, shortlength, charge) &
!!$     !$omp private(iat, ityp, rloc, atint)
!!$     !$omp do reduction(+: psoffset, shortlength, charge)
     do iat=1,at%astruct%nat
        ityp=at%astruct%iatype(iat)
        if (at%npspcode(ityp) == PSPCODE_PAW) then
           shortlength = shortlength + at%epsatm(ityp)
        else if (at%npspcode(ityp) == PSPCODE_GTH .or. &
             & at%npspcode(ityp) == PSPCODE_HGH .or. &
             & at%npspcode(ityp) == PSPCODE_HGH_K .or. &
             & at%npspcode(ityp) == PSPCODE_HGH_K_NLCC) then
           rloc=at%psppar(0,0,ityp)
           atint=at%psppar(0,1,ityp)+3.0_gp*at%psppar(0,2,ityp)+&
                15.0_gp*at%psppar(0,3,ityp)+105.0_gp*at%psppar(0,4,ityp)
           psoffset=psoffset+twopitothreehalf*rloc**3*atint
           shortlength=shortlength+real(at%nelpsp(ityp),gp)*rloc**2*2.d0*pi
        end if
        charge=charge+real(at%nelpsp(ityp),gp)
     end do
!!$     !$omp end do
!!$     !$omp end parallel

     !print *,'psoffset',psoffset,'pspcore', &
     !    (psoffset+shortlength)*charge/(at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3))
     !if (iproc ==0) print *,'eion',eion,charge/ucvol*(psoffset+shortlength)
     !correct ionic energy taking into account the PSP core correction
     !write(*,*) "EION", eion
     !write(*,*) "PSP CORE", charge/ucvol*(psoffset+shortlength), psoffset, shortlength
     eion=eion+charge/ucvol*(psoffset+shortlength)

     !symmetrization of ewald stress (probably not needed)
     if (at%astruct%sym%symObj >= 0) call symm_stress(ewaldstr,at%astruct%sym%symObj)
     !PSP core correction of the stress tensor (diag.)
     ewaldstr(1:3)=ewaldstr(1:3)-charge*(psoffset+shortlength)/ucvol/ucvol

!!$     if (iproc == 0) then
!!$        write(*,*) 'STRESS TENSOR: EWALD + PSP-CORE'
!!$        write(*,*) ewaldstr(1:3)
!!$        write(*,*) ewaldstr(6),ewaldstr(5),ewaldstr(4)
!!$     end if

     !!-     !in the surfaces case, correct the energy term following (J.Chem.Phys. 111(7)-3155, 1999)
     !!-     if (at%astruct%geocode == 'S') then
     !!-        !calculate the Mz dipole component (which in our case corresponds to y direction)
     !!-        !first calculate the center of mass
     !!-        cmassy=0.0_gp
     !!-        do iat=1,at%astruct%nat
     !!-           cmassy=cmassy+rxyz(2,iat)
     !!-        end do
     !!-
     !!-        Mz=0.0_gp
     !!-        do iat=1,at%astruct%nat
     !!-           ityp=at%astruct%iatype(iat)
     !!-           Mz=Mz+real(at%nelpsp(ityp),gp)*(rxyz(2,iat)-cmassy)
     !!-        end do
     !!-
     !!-        !correct energy and forces in the y direction
     !!-        eion=eion+0.5_gp/ucvol*Mz**2
     !!-        do iat=1,at%astruct%nat
     !!-           ityp=at%astruct%iatype(iat)
     !!-           fion(2,iat)=fion(2,iat)-real(at%nelpsp(ityp),gp)/ucvol*Mz
     !!-           if (nproc==1 .and. slowion) print *,'iat,fion',iat,(fion(j1,iat),j1=1,3)
     !!-        end do
     !!-
     !!-     end if

  else if (at%astruct%geocode == 'F') then

     eion=0.0_gp
     eself=0.0_gp

     !LR: commented hessian as not currently using it

     ! Parallelization over the atoms.
!!! First distribute evenly...
     !!natp = at%astruct%nat/nproc
     !!isat = iproc*natp
!!! ... and now distribute the remaining atoms.
     !!ii = at%astruct%nat-nproc*natp
     !!if (iproc<ii) natp = natp + 1
     !!isat = isat + min(iproc,ii)
     call distribute_on_tasks(at%astruct%nat, iproc, nproc, natp, isat)


!!!$omp parallel default(none) &
!!!$omp private(iat,ityp,rx,ry,rz,fxion,fyion,fzion,jtyp,chgprod,dist,jat) &
!!!$omp shared(at,rxyz,fion,eself,eion)
!!!$omp do reduction(+:eself,eion)
     call f_zero(fion)
     !do iat=1,at%astruct%nat
     do iat=1,natp
        iiat=iat+isat
        ityp=at%astruct%iatype(iiat)
        rx=rxyz(1,iiat)
        ry=rxyz(2,iiat)
        rz=rxyz(3,iiat)
        !initialization of the forces
        fxion=0.0_gp
        fyion=0.0_gp
        fzion=0.0_gp
        !initialisation of the hessian
        !hxx=0.0_gp
        !hxy=0.0_gp
        !hxz=0.0_gp
        !hyy=0.0_gp
        !hyz=0.0_gp
        !hzz=0.0_gp

        ! Ion-ion interaction
        !$omp parallel if (iiat>1000) &
        !$omp default(none) &
        !$omp shared(iiat,rx,ry,rz,rxyz,at,ityp,eion,fxion,fyion,fzion) &
        !$omp private(jat,dist,jtyp,chgprod)
        !$omp do schedule(static) reduction(+:eion,fxion,fyion,fzion)
        do jat=1,iiat-1
           dist=sqrt((rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2)
           jtyp=at%astruct%iatype(jat)
           chgprod=real(at%nelpsp(jtyp),gp)*real(at%nelpsp(ityp),gp)
           eion=eion+chgprod/dist
           !forces
           fxion=fxion+chgprod/(dist**3)*(rx-rxyz(1,jat))
           fyion=fyion+chgprod/(dist**3)*(ry-rxyz(2,jat))
           fzion=fzion+chgprod/(dist**3)*(rz-rxyz(3,jat))
           !hessian matrix
           !hxx=hxx+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))**2-chgprod/(dist**3)
           !hxy=hxy+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(ry-rxyz(2,jat))
           !hxz=hxz+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(rz-rxyz(3,jat))
           !hyy=hyy+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))**2-chgprod/(dist**3)
           !hyz=hyz+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))*(rz-rxyz(3,jat))
           !hzz=hzz+3.0_gp*chgprod/(dist**5)*(rz-rxyz(3,jat))**2-chgprod/(dist**3)
        end do
        !$omp end do
        !$omp end parallel
        !$omp parallel if (at%astruct%nat-iiat>1000) &
        !$omp default(none) &
        !$omp shared(iiat,rx,ry,rz,rxyz,at,ityp,eion,fxion,fyion,fzion) &
        !$omp private(jat,dist,jtyp,chgprod)
        !$omp do schedule(static) reduction(+:fxion,fyion,fzion)
        do jat=iiat+1,at%astruct%nat
           dist=sqrt((rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2)
           jtyp=at%astruct%iatype(jat)
           chgprod=real(at%nelpsp(jtyp),gp)*real(at%nelpsp(ityp),gp)

           !forces
           fxion=fxion+chgprod/(dist**3)*(rx-rxyz(1,jat))
           fyion=fyion+chgprod/(dist**3)*(ry-rxyz(2,jat))
           fzion=fzion+chgprod/(dist**3)*(rz-rxyz(3,jat))
           !hessian matrix
           !hxx=hxx+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))**2-chgprod/(dist**3)
           !hxy=hxy+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(ry-rxyz(2,jat))
           !hxz=hxz+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(rz-rxyz(3,jat))
           !hyy=hyy+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))**2-chgprod/(dist**3)
           !hyz=hyz+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))*(rz-rxyz(3,jat))
           !hzz=hzz+3.0_gp*chgprod/(dist**5)*(rz-rxyz(3,jat))**2-chgprod/(dist**3)
        end do
        !$omp end do
        !$omp end parallel

        fion(1,iiat)=fxion
        fion(2,iiat)=fyion
        fion(3,iiat)=fzion

        !if (nproc==1 .and. slowion) print *,'iat,fion',iat,(fion(j1,iat),j1=1,3)
        !energy which comes from the self-interaction of the spread charge
        eself=eself+real(at%nelpsp(ityp)**2,gp)*0.5_gp*sqrt(1.d0/pi)/at%psppar(0,0,ityp)
     end do
!!!$omp end do
!!!$omp end parallel
     call fmpi_allreduce(fion, FMPI_SUM, comm=bigdft_mpi%mpi_comm)
     call fmpi_allreduce(eion, 1, FMPI_SUM, comm=bigdft_mpi%mpi_comm)
     call fmpi_allreduce(eself, 1, FMPI_SUM, comm=bigdft_mpi%mpi_comm)

     !if (nproc==1 .and. slowion) print *,'eself',eself

  end if

  !for the surfaces BC,
  !activate for the moment only the slow calculation of the ionic energy and forces
  !the slowion command has also to be activated for the cavity calculation
  !if (at%astruct%geocode == 'S' .or. at%astruct%geocode == 'P') slowion=.true.
  if (at%astruct%geocode == 'S' .or. pkernel%method /= 'VAC') slowion=.true.

  slowion_if: if (slowion) then

     !case of slow ionic calculation
     !conditions for periodicity in the three directions
     perx=(at%astruct%geocode /= 'F')
     pery=(at%astruct%geocode == 'P')
     perz=(at%astruct%geocode /= 'F')

     call ext_buffers(perx,nbl1,nbr1)
     call ext_buffers(pery,nbl2,nbr2)
     call ext_buffers(perz,nbl3,nbr3)
     !the ions corresponds to gaussian charges disposed in the same way as the pseudopotentials

     !first calculate the self-energy and the forces
     !(the latter are zero for a symmetric grid distribution)

     !self energy initialisation
     eself=0.0_gp
     do iat=1,at%astruct%nat

        fion(1,iat)=0.0_gp
        fion(2,iat)=0.0_gp
        fion(3,iat)=0.0_gp

        ityp=at%astruct%iatype(iat)
        rloc=at%psppar(0,0,ityp)
        charge=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
        prefactor=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**5)

        !calculate the self energy of the isolated bc
        eself=eself+real(at%nelpsp(ityp),gp)**2/rloc

     end do

     eself=0.5_gp/sqrt(pi)*eself

     !if (nproc==1)
     !print *,'iproc,eself',iproc,eself
     call f_zero(dpbox%mesh%ndims(1)*dpbox%mesh%ndims(2)*dpbox%n3pi,pot_ion(1))

     !if (dpbox%n3pi >0 ) then
     !this section is important only if we are not doing a cavity calculation
     if (pkernel%method == 'VAC' .and. dpbox%n3pi >0 ) then

        !then calculate the hartree energy and forces of the charge distributions
        !(and save the values for the ionic potential)

        !Determine the maximal bounds for mpx, mpy, mpz (1D-integral)
        call mp_range(at%multipole_preserving,at%mp_isf,1,&
             hxh,hyh,hzh,maxval(at%psppar(0,0,:)),mpnx,mpny,mpnz)
        !Separable function: do 1-D integrals before and store it.
        mpx = f_malloc( (/ 0 .to. mpnx /),id='mpx')
        mpy = f_malloc( (/ 0 .to. mpny /),id='mpy')
        mpz = f_malloc( (/ 0 .to. mpnz /),id='mpz')

        atit = atoms_iter(at%astruct)
        do while(atoms_iter_next(atit))

           !!-        do iat=1,at%astruct%nat
           !!-           ityp=at%astruct%iatype(iat)
           !!-          rx=rxyz(1,iat)
           !!-          ry=rxyz(2,iat)
           !!-          rz=rxyz(3,iat)


           rx=rxyz(1,atit%iat)
           ry=rxyz(2,atit%iat)
           rz=rxyz(3,atit%iat)

           !!-           rloc=at%psppar(0,0,ityp)
           !!-           charge=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
!!$           rloc=at%psppar(0,0,atit%ityp)
!!$           rlocinv2sq=0.5_gp/rloc**2
!!$           charge=real(at%nelpsp(atit%ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
!!$           call gaussian_real_space_set(g,rlocinv2sq,1,[charge],[0,0,0])
           call atomic_charge_density(g,at,atit)
           call three_dimensional_density(dpbox%bitp,g,-1.0_dp,rxyz(1,atit%iat),pot_ion)

!!$           !cutoff of the range
!!$           cutoff=10.0_gp*rloc
!!$           if (at%multipole_preserving) then
!!$              !We want to have a good accuracy of the last point rloc*10
!!$              !cutoff=cutoff+max(hxh,hyh,hzh)*real(16,kind=gp)
!!$              cutoff=cutoff+max(hxh,hyh,hzh)*real(at%mp_isf,kind=gp)
!!$           end if
!!$
!!$           isx=floor((rx-cutoff)/hxh)
!!$           isy=floor((ry-cutoff)/hyh)
!!$           isz=floor((rz-cutoff)/hzh)
!!$
!!$           iex=ceiling((rx+cutoff)/hxh)
!!$           iey=ceiling((ry+cutoff)/hyh)
!!$           iez=ceiling((rz+cutoff)/hzh)
!!$
!!$           !Separable function: do 1-D integrals before and store it.
!!$           !call mp_calculate(rx,ry,rz,hxh,hyh,hzh,cutoff,rlocinv2sq,at%multipole_preserving,mpx,mpy,mpz)
!!$           !!mpx = f_malloc( (/ isx.to.iex /),id='mpx')
!!$           !!mpy = f_malloc( (/ isy.to.iey /),id='mpy')
!!$           !!mpz = f_malloc( (/ isz.to.iez /),id='mpz')
!!$           do i1=isx,iex
!!$              mpx(i1-isx) = mp_exp(hxh,rx,rlocinv2sq,i1,0,at%multipole_preserving)
!!$           end do
!!$           do i2=isy,iey
!!$              mpy(i2-isy) = mp_exp(hyh,ry,rlocinv2sq,i2,0,at%multipole_preserving)
!!$           end do
!!$           do i3=isz,iez
!!$              mpz(i3-isz) = mp_exp(hzh,rz,rlocinv2sq,i3,0,at%multipole_preserving)
!!$           end do
!!$
!!$           !these nested loops will be used also for the actual ionic forces, to be recalculated
!!$           do i3=isz,iez
!!$              zp = mpz(i3-isz)
!!$              if (abs(zp) < mp_tiny) cycle
!!$              z=real(i3,gp)*hzh-rz
!!$              !call ind_positions(perz,i3,n3,j3,goz)
!!$              call ind_positions_new(perz,i3,n3i,j3,goz)
!!$              j3=j3+nbl3+1
!!$              do i2=isy,iey
!!$                 yp = zp*mpy(i2-isy)
!!$                 if (abs(yp) < mp_tiny) cycle
!!$                 y=real(i2,gp)*hyh-ry
!!$                 !call ind_positions(pery,i2,n2,j2,goy)
!!$                 call ind_positions_new(pery,i2,n2i,j2,goy)
!!$                 do i1=isx,iex
!!$                    xp = yp*mpx(i1-isx)
!!$                    if (abs(xp) < mp_tiny) cycle
!!$                    x=real(i1,gp)*hxh-rx
!!$                    !call ind_positions(perx,i1,n1,j1,gox)
!!$                    call ind_positions_new(perx,i1,n1i,j1,gox)
!!$                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
!!$                       ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s)*n1i*n2i
!!$                       pot_ion(ind)=pot_ion(ind)-xp*charge
!!$                    endif
!!$                 enddo
!!$              enddo
!!$           enddo
!!$
!!$           !De-allocate the 1D temporary arrays for separability
!!$           !call f_free(mpx,mpy,mpz)

        enddo

        !De-allocate the 1D temporary arrays for separability
        call f_free(mpx,mpy,mpz)

     end if

  end if slowion_if

  !in the case of cavity the ionic energy is only considered as the self energy
  nocavity_if: if (pkernel%method /= 'VAC') then
     eion=-eself
  else if (slowion) then
     !now call the Poisson Solver for the global energy forces
     call H_potential('D',pkernel,pot_ion,pot_ion,ehart,-2.0_gp*psoffset,.false.)

     eion=ehart-eself

     !print *,'ehart,eself',iproc,ehart,eself

     !if (nproc==1)
     !print *,'iproc,eion',iproc,eion

     !calculate the forces near the atom due to the error function part of the potential
     !calculate forces for all atoms only in the distributed part of the simulation box
     if (dpbox%n3pi >0 ) then

        !Determine the maximal bounds for mpx, mpy, mpz (1D-integral)
        call mp_range(at%multipole_preserving,at%mp_isf,1,&
            hxh,hyh,hzh,maxval(at%psppar(0,0,:)),mpnx,mpny,mpnz)
        !Separable function: do 1-D integrals before and store it.
        mpx = f_malloc( (/ 0 .to. mpnx /),id='mpx')
        mpy = f_malloc( (/ 0 .to. mpny /),id='mpy')
        mpz = f_malloc( (/ 0 .to. mpnz /),id='mpz')

        atit = atoms_iter(at%astruct)
        do while(atoms_iter_next(atit))

           !!-        do iat=1,at%astruct%nat
           !!-           ityp=at%astruct%iatype(iat)
           !!-          rx=rxyz(1,iat)
           !!-          ry=rxyz(2,iat)
           !!-          rz=rxyz(3,iat)

           rx=rxyz(1,atit%iat)
           ry=rxyz(2,atit%iat)
           rz=rxyz(3,atit%iat)

           !inizialization of the forces
           fxerf=0.0_gp
           fyerf=0.0_gp
           fzerf=0.0_gp

           call atomic_charge_density(g,at,atit)
           call set_box_around_gaussian(dpbox%bitp,g,rxyz(1,atit%iat))
           !call box_iter_set_nbox(dpbox%bitp,nbox=gaussian_nbox(rxyz(1,atit%iat),dpbox%bitp%mesh,g))
           do while(box_next_point(dpbox%bitp))
              !r = distance(dpbox%bitp%mesh,dpbox%bitp%rxyz,rxyz(1,atit%iat))
              xp=gaussian_radial_value(g,rxyz(1,atit%iat),dpbox%bitp)/g%sigma**2
              !error fucntion per volume element for the integration
              dpbox%bitp%tmp=closest_r(dpbox%bitp%mesh,dpbox%bitp%rxyz,rxyz(1,atit%iat))
              Vel=pot_ion(dpbox%bitp%ind)*dpbox%bitp%mesh%volume_element
              fxerf=fxerf+xp*Vel*dpbox%bitp%tmp(1) !warning, this should be absolute x
              fyerf=fyerf+xp*Vel*dpbox%bitp%tmp(2)
              fzerf=fzerf+xp*Vel*dpbox%bitp%tmp(3)
           end do
           call box_iter_expand_nbox(dpbox%bitp)
           fion(1,atit%iat) = fion(1,atit%iat) + fxerf
           fion(2,atit%iat) = fion(2,atit%iat) + fyerf
           fion(3,atit%iat) = fion(3,atit%iat) + fzerf

           !print *,'fxerf',fxerf,fyerf,fzerf


!!$           !local part
!!$           !!-           rloc=at%psppar(0,0,ityp)
!!$           !!-           prefactor=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**5)
!!$           rloc=at%psppar(0,0,atit%ityp)
!!$           prefactor=real(at%nelpsp(atit%ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**5)
!!$           rlocinv2sq=0.5_gp/rloc**2
!!$           !maximum extension of the gaussian
!!$           cutoff=10.0_gp*rloc
!!$           if (at%multipole_preserving) then
!!$              !We want to have a good accuracy of the last point rloc*10
!!$              !cutoff=cutoff+max(hxh,hyh,hzh)*real(16,kind=gp)
!!$              cutoff=cutoff+max(hxh,hyh,hzh)*real(at%mp_isf,kind=gp)
!!$           end if
!!$
!!$           !final result of the forces
!!$
!!$           isx=floor((rx-cutoff)/hxh)
!!$           isy=floor((ry-cutoff)/hyh)
!!$           isz=floor((rz-cutoff)/hzh)
!!$
!!$           iex=ceiling((rx+cutoff)/hxh)
!!$           iey=ceiling((ry+cutoff)/hyh)
!!$           iez=ceiling((rz+cutoff)/hzh)
!!$
!!$           !Separable function: do 1-D integrals before and store it.
!!$           !call mp_calculate(rx,ry,rz,hxh,hyh,hzh,cutoff,rlocinv2sq,at%multipole_preserving,mpx,mpy,mpz)
!!$           !!mpx = f_malloc( (/ isx.to.iex /),id='mpx')
!!$           !!mpy = f_malloc( (/ isy.to.iey /),id='mpy')
!!$           !!mpz = f_malloc( (/ isz.to.iez /),id='mpz')
!!$           do i1=isx,iex
!!$              mpx(i1-isx) = mp_exp(hxh,rx,rlocinv2sq,i1,0,at%multipole_preserving)
!!$           end do
!!$           do i2=isy,iey
!!$              mpy(i2-isy) = mp_exp(hyh,ry,rlocinv2sq,i2,0,at%multipole_preserving)
!!$           end do
!!$           do i3=isz,iez
!!$              mpz(i3-isz) = mp_exp(hzh,rz,rlocinv2sq,i3,0,at%multipole_preserving)
!!$           end do
!!$
!!$           do i3=isz,iez
!!$              z=real(i3,gp)*hzh-rz
!!$              zp = mpz(i3-isz)
!!$              if (abs(zp) < mp_tiny) cycle
!!$              !call ind_positions(perz,i3,n3,j3,goz)
!!$              call ind_positions_new(perz,i3,n3i,j3,goz)
!!$              j3=j3+nbl3+1
!!$              do i2=isy,iey
!!$                 yp = zp*mpy(i2-isy)
!!$                 if (abs(yp) < mp_tiny) cycle
!!$                 y=real(i2,gp)*hyh-ry
!!$                 !call ind_positions(pery,i2,n2,j2,goy)
!!$                 call ind_positions_new(pery,i2,n2i,j2,goy)
!!$                 do i1=isx,iex
!!$                    xp = yp*mpx(i1-isx)
!!$                    if (abs(xp) < mp_tiny) cycle
!!$                    x=real(i1,gp)*hxh-rx
!!$                    !call ind_positions(perx,i1,n1,j1,gox)
!!$                    call ind_positions_new(perx,i1,n1i,j1,gox)
!!$                    r2=x**2+y**2+z**2
!!$                    arg=r2/rloc**2
!!$                    xp=exp(-.5_gp*arg)
!!$                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
!!$                       ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s)*n1i*n2i
!!$                       !error function part
!!$                       Vel=pot_ion(ind)
!!$                       fxerf=fxerf+xp*Vel*x
!!$                       fyerf=fyerf+xp*Vel*y
!!$                       fzerf=fzerf+xp*Vel*z
!!$                    endif
!!$                 end do
!!$              end do
!!$           end do
!!$           !!-           fion(1,iat)=fion(1,iat)+(hxh*hyh*hzh*prefactor)*fxerf
!!$           !!-           fion(2,iat)=fion(2,iat)+(hxh*hyh*hzh*prefactor)*fyerf
!!$           !!-           fion(3,iat)=fion(3,iat)+(hxh*hyh*hzh*prefactor)*fzerf
!!$
!!$           fion(1,atit%iat) = fion(1,atit%iat) + (hxh*hyh*hzh*prefactor)*fxerf
!!$           fion(2,atit%iat) = fion(2,atit%iat) + (hxh*hyh*hzh*prefactor)*fyerf
!!$           fion(3,atit%iat) = fion(3,atit%iat) + (hxh*hyh*hzh*prefactor)*fzerf
!!$
!!$           !if (nproc==1) print *,'iat,fion',iat,(fion(j1,iat),j1=1,3)
!!$
!!$           !!-        write(10+iat,'(1x,f8.3,i5,(1x,3(1x,1pe12.5)))',advance='no') &
!!$           !!-             hxh,iat,(fion(j1,iat),j1=1,3)
!!$

        end do !do iat

        !De-allocate the 1D temporary arrays for separability
        call f_free(mpx,mpy,mpz)

     end if !if dpbox%n3pi

     if (pkernel%mpi_env%nproc > 1) then
        call fmpi_allreduce(fion,FMPI_SUM,comm=pkernel%mpi_env%mpi_comm)
     end if

     !if (iproc ==0) print *,'eion',eion,psoffset,shortlength

  end if nocavity_if


  ! Add contribution from constant electric field to the forces
  if(any(elecfield(1:3)/=0._gp)) then
     call center_of_charge(at,rxyz,cc)
     do iat=1,at%astruct%nat
        ityp=at%astruct%iatype(iat)
        charge=real(at%nelpsp(ityp),gp)
        fion(1:3,iat)=fion(1:3,iat)+(charge*elecfield(1:3))
        !ry=rxyz(2,iat)
        !eion=eion-(charge*elecfield)*ry
        de=0.0_gp
        do i=1,3
           de=de+elecfield(i)*(rxyz(i,iat)-cc(i))
        end do
        !eion=eion-charge*sum(elecfield(1:3)*rxyz(1:3,iat))
        eion=eion-charge*de
     enddo
  end if

  if (iproc == 0) then
     if(all(elecfield(1:3)==0._gp)) then
        !      write(*,'(1x,a,1pe22.14)') 'ion-ion interaction energy',eion
        call yaml_map('Ion-Ion interaction energy',eion,fmt='(1pe22.14)')
     else
        call yaml_map('Ion-electric field interaction energy',eion,fmt='(1pe22.14)')
     endif
  end if

  ! Add empiric correction for Van der Waals forces and energy.
  call vdwcorrection_calculate_energy(edisp,rxyz,at,dispersion)
  if (iproc == 0 .and. edisp /= 0.0_gp) then
     call yaml_map('Dispersion Correction Energy (Ha)',edisp,fmt='(1pe22.14)', advance = "no")
     call yaml_comment(trim(vdw_correction_names(dispersion + 1)))
  end if

  call vdwcorrection_calculate_forces(fdisp,rxyz,at,dispersion)
  call vdwcorrection_freeparams()

  if (at%multipole_preserving) call finalize_real_space_conversion()
  call f_release_routine()
  call timing(iproc,'ionic_energy','OF')

END SUBROUTINE IonicEnergyandForces


!> Create the effective ionic potential (main ionic + counter ions)
subroutine createEffectiveIonicPotential(iproc, verb, input, atoms, rxyz, shift, &
     & dpbox, pkernel, pot_ion, rho_ion, elecfield, psoffset)

  use module_base
  use module_dpbox, only: denspot_distribution
  use module_types

  implicit none

  !Arguments
  integer, intent(in) :: iproc
  logical, intent(in) :: verb
  !!-  real(gp), intent(in) :: hxh,hyh,hzh
  real(gp), intent(in) :: psoffset
  type(atoms_data), intent(in) :: atoms
  !!-  type(locreg_descriptors), intent(in) :: Glr
  type(input_variables), intent(in) :: input
  type(denspot_distribution), intent(inout) :: dpbox
  real(gp), dimension(3), intent(in) :: elecfield
  real(gp), dimension(3), intent(in) :: shift
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  type(coulomb_operator), intent(inout) :: pkernel
  real(wp), dimension(*), intent(inout) :: pot_ion
  real(wp), dimension(*), intent(inout) :: rho_ion

  !Local variables
  logical :: counterions
  real(dp), dimension(:), allocatable :: counter_ions
  integer :: ncounter_ions

  call f_routine(id='createEffectiveIonicPotential')

  ! Compute the main ionic potential.
  call createIonicPotential(iproc, verb, atoms, rxyz, &
       & elecfield, dpbox, pkernel, pot_ion, rho_ion, psoffset)

  !inquire for the counter_ion potential calculation (for the moment only xyz format)
  inquire(file='posinp_ci.xyz',exist=counterions)
  if (counterions) then
     if (dpbox%n3pi > 0) then
        !counter_ions = f_malloc(dpbox%ndims(1)*dpbox%ndims(2)*dpbox%n3pi,id='counter_ions')
        ncounter_ions = dpbox%mesh%ndims(1)*dpbox%mesh%ndims(2)*dpbox%n3pi
     else
        !counter_ions = f_malloc(1,id='counter_ions')
        ncounter_ions = 1
     end if
     counter_ions = f_malloc(ncounter_ions,id='counter_ions')

     call CounterIonPotential(iproc,input,shift,dpbox,pkernel,ncounter_ions,counter_ions)

     !sum that to the ionic potential
     call axpy(dpbox%mesh%ndims(1)*dpbox%mesh%ndims(2)*dpbox%n3pi,1.0_dp,counter_ions(1),1,&
          &   pot_ion(1),1)

     call f_free(counter_ions)
  end if

  call f_release_routine()

END SUBROUTINE createEffectiveIonicPotential


!> Create the ionic potential
subroutine createIonicPotential(iproc,verb,at,rxyz,&
     elecfield,dpbox,pkernel,pot_ion,rho_ion,psoffset)

  use module_base
  use module_types
  use yaml_output
  use m_paw_numeric, only: paw_splint
  use multipole_preserving
  use module_atoms
  use ao_inguess, only: atomic_info
  use module_dpbox
  !  use module_interfaces, only: mp_calculate
  !  use module_interfaces, except_this_one => createIonicPotential
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use abi_interfaces_numeric, only: abi_derf_ab
  use public_enums, only: PSPCODE_PAW, PSPCODE_PSPIO, PSPCODE_GTH, PSPCODE_HGH, &
       & PSPCODE_HGH_K, PSPCODE_HGH_K_NLCC
  use pspiof_m, only: pspiof_pspdata_get_vlocal, pspiof_potential_t, &
       & pspiof_potential_eval_deriv, pspiof_potential_eval_deriv2, &
       & pspiof_potential_eval, pspiof_pspdata_get_zvalence
  use bounds, only: ext_buffers
  use box
  use gaussians
  implicit none

  !Arguments
  !!-  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: iproc
  !!-  integer, intent(in) :: n1,n2,n3
  logical, intent(in) :: verb
  real(gp), intent(in) :: psoffset
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3), intent(in) :: elecfield
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  type(denspot_distribution), intent(inout) :: dpbox
  type(coulomb_operator), intent(inout) :: pkernel
  real(wp), dimension(*), intent(inout) :: pot_ion
  real(dp), dimension(*), intent(out) :: rho_ion
  !Local variables
  character(len = 3) :: quiet
  character(len=20) :: rprbstr
  !  logical, parameter :: efwrite=.false.
  real(gp), parameter :: mp_tiny = 1.e-30_gp !<put zero when the value is lower than this
  logical :: perx,pery,perz,gox,goy,goz,yespar
  logical :: htoobig=.false.,check_potion=.false.!use_iterator=.false.
  integer :: i1,i2,i3,ierr !n(c) nspin
  integer :: nloc,iloc
  integer  :: i3s,n3pi,nbl1,nbr1,nbl2,nbl3,nbr2,nbr3
  real(kind=8) :: raux1(1),rr1(1)
  integer :: iex,iey,iez,ind,indj3,indj23,isx,isy,isz,j1,j2,j3
  integer :: n1i,n2i,n3i,mpnx,mpny,mpnz
  real(gp) :: hxh,hyh,hzh
  real(gp) :: rloc,cutoff,r2,arg,xp,tt,rx,ry,rz
  real(gp) :: tt_tot,potxyz,aval
  real(gp) :: rlocinvsq,rlocinv2sq
  real(gp) :: x,y,z,yp,zp,zsq,yzsq,rholeaked,rholeaked_tot,psoff,psoff_tot
  !  real(gp), dimension(1) :: raux,rr
  real(wp) :: maxdiff
  real(gp) :: ehart
  logical, dimension(3) :: peri
  real(dp), dimension(3) :: charges_mpi
  real(dp), dimension(:), allocatable :: potion_corr, pot_add
  logical, parameter :: pawErfCorrection = .true.
  real(gp), dimension(:), allocatable  :: mpx,mpy,mpz
  !real(dp), dimension(:), allocatable :: den_aux
  type(atoms_iterator) :: atit
  type(dpbox_iterator) :: boxit
  type(box_iterator) :: bitp
  type(gaussian_real_space) :: g
  type(pspiof_potential_t) :: pot
  real(gp), dimension(3) :: center_of_charge_ions

  call f_routine(id='createIonicPotential')
  call timing(iproc,'CrtLocPot     ','ON')

  !Initialize the work arrays needed to integrate with isf
  !Determine the number of points depending on the min rloc
  if (at%multipole_preserving) &
     call initialize_real_space_conversion(isf_m=at%mp_isf,&
          rlocs=at%psppar(0,0,:),verbose=(iproc==0))

  ! Aliasing
  hxh = dpbox%mesh%hgrids(1)
  hyh = dpbox%mesh%hgrids(2)
  hzh = dpbox%mesh%hgrids(3)
  n1i = dpbox%mesh%ndims(1)
  n2i = dpbox%mesh%ndims(2)
  n3i = dpbox%mesh%ndims(3)
  i3s = dpbox%i3s+dpbox%i3xcsh
  n3pi = dpbox%n3pi

  ! Ionic charge (must be calculated for the PS active processes)
  rholeaked=0.d0
  ! Ionic energy (can be calculated for all the processors)

  ! Recalculated PSoffset on grid for real space potentials.
  psoff = 0._gp

  !Creates charge density arising from the ionic PSP cores
  call f_zero(n1i*n2i*dpbox%n3pi,pot_ion(1))
  if (any(at%npspcode == PSPCODE_PAW) .or. any(at%npspcode == PSPCODE_PSPIO)) then
     pot_add = f_malloc0(n1i*n2i*dpbox%n3pi, id = "pot_add")
  end if

  !conditions for periodicity in the three directions
  peri=cell_periodic_dims(dpbox%mesh)
  perx=peri(1)!(dpbox%geocode /= 'F')
  pery=peri(2)!(dpbox%geocode == 'P')
  perz=peri(3)!(dpbox%geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  if (dpbox%n3pi >0 .and. .not. htoobig) then

     !Determine the maximal bounds for mpx, mpy, mpz (1D-integral)
     call mp_range(at%multipole_preserving,at%mp_isf,at%astruct%nat,&
         hxh,hyh,hzh,maxval(at%psppar(0,0,:)),mpnx,mpny,mpnz)
     !Separable function: do 1-D integrals before and store it.
     mpx = f_malloc( (/ 0 .to. mpnx /),id='mpx')
     mpy = f_malloc( (/ 0 .to. mpny /),id='mpy')
     mpz = f_malloc( (/ 0 .to. mpnz /),id='mpz')

     atit = atoms_iter(at%astruct)
     do while(atoms_iter_next(atit))

!!$        !!-     do iat=1,at%astruct%nat
!!$        !!-        ityp=at%astruct%iatype(iat)
!!$        !!-       rx=rxyz(1,iat)
!!$        !!-       ry=rxyz(2,iat)
!!$        !!-       rz=rxyz(3,iat)
!!$        rx=rxyz(1,atit%iat)
!!$        ry=rxyz(2,atit%iat)
!!$        rz=rxyz(3,atit%iat)
!!$
!!$        !!-        rloc=at%psppar(0,0,ityp)
!!$        !!-        charge=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
!!$        rloc=at%psppar(0,0,atit%ityp)
!!$        rlocinv2sq=0.5_gp/rloc**2
!!$        charge=real(at%nelpsp(atit%ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
!!$
!!$        !cutoff of the range
!!$        cutoff=10.0_gp*rloc
!!$        if (at%multipole_preserving) then
!!$           !We want to have a good accuracy of the last point rloc*10
!!$           cutoff=cutoff+max(hxh,hyh,hzh)*real(at%mp_isf,kind=gp)
!!$        end if
!!$
!!$        isx=floor((rx-cutoff)/hxh)
!!$        isy=floor((ry-cutoff)/hyh)
!!$        isz=floor((rz-cutoff)/hzh)
!!$
!!$        iex=ceiling((rx+cutoff)/hxh)
!!$        iey=ceiling((ry+cutoff)/hyh)
!!$        iez=ceiling((rz+cutoff)/hzh)
!!$
!!$        !Separable function: do 1-D integrals before and store it.
!!$        !call mp_calculate(rx,ry,rz,hxh,hyh,hzh,cutoff,rlocinv2sq,at%multipole_preserving,mpx,mpy,mpz)
!!$        !!mpx = f_malloc( (/ isx.to.iex /),id='mpx')
!!$        !!mpy = f_malloc( (/ isy.to.iey /),id='mpy')
!!$        !!mpz = f_malloc( (/ isz.to.iez /),id='mpz')
!!$        do i1=isx,iex
!!$           mpx(i1-isx) = mp_exp(hxh,rx,rlocinv2sq,i1,0,at%multipole_preserving)
!!$        end do
!!$        do i2=isy,iey
!!$           mpy(i2-isy) = mp_exp(hyh,ry,rlocinv2sq,i2,0,at%multipole_preserving)
!!$        end do
!!$        do i3=isz,iez
!!$           mpz(i3-isz) = mp_exp(hzh,rz,rlocinv2sq,i3,0,at%multipole_preserving)
!!$        end do

!!$           !Calculate Ionic Density using HGH parameters.
!!$           !Eq. 1.104, T. Deutsch and L. Genovese, JDN. 12, 2011
!!$           do i3=isz,iez
!!$              zp = mpz(i3-isz)
!!$              !call ind_positions(perz,i3,n3,j3,goz)
!!$              call ind_positions_new(perz,i3,n3i,j3,goz)
!!$              j3=j3+nbl3+1
!!$              if ( goz .and. (j3<i3s.or.j3>i3s+n3pi-1) ) cycle
!!$              indj3=(j3-i3s)*n1i*n2i
!!$              do i2=isy,iey
!!$                 yp = zp*mpy(i2-isy)
!!$                 !call ind_positions(pery,i2,n2,j2,goy)
!!$                 call ind_positions_new(pery,i2,n2i,j2,goy)
!!$                 if (goz.and.(.not.goy)) cycle
!!$                 indj23=1+nbl1+(j2+nbl2)*n1i+indj3
!!$                 do i1=isx,iex
!!$                    xp = yp*mpx(i1-isx)
!!$                    !call ind_positions(perx,i1,n1,j1,gox)
!!$                    call ind_positions_new(perx,i1,n1i,j1,gox)
!!$                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1 .and. goy .and. gox) then
!!$                       ind=j1+indj23
!!$                       pot_ion(ind)=pot_ion(ind)-xp*charge
!!$                       !write(*,'(4(i0,1x),2(1pe20.10))') i1,i2,i3,ind,xp,pot_ion(ind)
!!$                    else if (.not. goz ) then
!!$                       rholeaked=rholeaked+xp*charge
!!$                    endif
!!$                 enddo
!!$              enddo
!!$           enddo

           if (at%npspcode(atit%ityp) == PSPCODE_PSPIO) &
                & pot = pspiof_pspdata_get_vlocal(at%pspio(atit%ityp))

           call atomic_charge_density(g,at,atit)
           call set_box_around_gaussian(dpbox%bitp,g,rxyz(1,atit%iat))
           do while(box_next_point(dpbox%bitp))
              !Take the HGH form for rho_L (long range)
              raux1(1)=-gaussian_radial_value(g,rxyz(1,atit%iat),dpbox%bitp)
              pot_ion(dpbox%bitp%ind)=pot_ion(dpbox%bitp%ind)+raux1(1)

              !Create the additional part, should be done also for HGH in fact.
              if (at%npspcode(atit%ityp) == PSPCODE_PAW .or. &
                   & at%npspcode(atit%ityp) == PSPCODE_PSPIO) then
                 rr1(1)=distance(dpbox%bitp%mesh,dpbox%bitp%rxyz,rxyz(1,atit%iat))
                 if (at%npspcode(atit%ityp) == PSPCODE_PAW) then
                    call paw_splint(at%pawtab(atit%ityp)%wvl%rholoc%msz, &
                         & at%pawtab(atit%ityp)%wvl%rholoc%rad, &
                         & at%pawtab(atit%ityp)%wvl%rholoc%d(:,3), &
                         & at%pawtab(atit%ityp)%wvl%rholoc%d(:,4), &
                         & 1,rr1,raux1,ierr)
                 else
                    raux1(1) = pspiof_potential_eval(pot, rr1(1))
                 end if
                 if (rr1(1) < 1e-6_dp) then
                    tt = -at%nelpsp(atit%ityp) / sqrt(2._dp * pi) / at%psppar(0,0,atit%ityp) * 2._dp
                 else
                    !call abi_derf_ab(tt, rr1(1) / sqrt(2._dp) / at%psppar(0,0,atit%ityp))
                    tt = safe_erf(rr1(1) / sqrt(2._dp) / at%psppar(0,0,atit%ityp))
                    tt = -at%nelpsp(atit%ityp) * tt / rr1(1)
                 end if
                 pot_add(dpbox%bitp%ind) = pot_add(dpbox%bitp%ind) + raux1(1) - tt
                 if (cell_geocode(dpbox%mesh) == 'P') psoff = psoff + raux1(1) - tt
              end if
           end do
           call box_iter_expand_nbox(dpbox%bitp)
           rholeaked=0.0_dp

     enddo

     !De-allocate for multipole preserving
     call f_free(mpx,mpy,mpz)

  end if

  ! Check
  tt=0.d0
  do j3=1,n3pi
     indj3=(j3-1)*n1i*n2i
     do i2= -nbl2,n2i-nbl2-1!2*n2+1+nbr2
        indj23=1+nbl1+(i2+nbl2)*n1i+indj3
        do i1= -nbl1,n1i-nbl1-1!2*n1+1+nbr1
           ind=i1+indj23
           tt=tt+pot_ion(ind)
        enddo
     enddo
  enddo

  tt=tt*dpbox%mesh%volume_element!hxh*hyh*hzh
  rholeaked=rholeaked*dpbox%mesh%volume_element!hxh*hyh*hzh
  psoff=psoff*dpbox%mesh%volume_element!hxh*hyh*hzh

  !print *,'test case input_rho_ion',iproc,i3start,i3end,n3pi,2*n3+16,tt
  !if rho_ion is needed for the SCF cycle copy in the array
  if (pkernel%method /= 'VAC') then
     call f_memcpy(n=n1i*n2i*dpbox%n3pi,src=pot_ion(1),dest=rho_ion(1))
  end if

  if (pkernel%mpi_env%nproc > 1) then
     charges_mpi(1)=tt
     charges_mpi(2)=rholeaked
     charges_mpi(3)=psoff

     !Reduce from all mpi proc
     call fmpi_allreduce(charges_mpi,FMPI_SUM,comm=pkernel%mpi_env%mpi_comm)

     tt_tot=charges_mpi(1)
     rholeaked_tot=charges_mpi(2)
     psoff_tot=charges_mpi(3)
  else
     tt_tot=tt
     rholeaked_tot=rholeaked
     psoff_tot=psoff
  end if

  if (verb) then
     call yaml_comment('Ionic Potential Creation',hfill='-')
     call yaml_map('Total ionic charge',tt_tot,fmt='(f26.12)')
     if (rholeaked_tot /= 0.0_gp) call yaml_map('Leaked charge',rholeaked_tot,fmt='(1pe10.3)')
     quiet = "no "
  else
     quiet = "yes"
  end if

  if (.not. htoobig) then
     call timing(iproc,'CrtLocPot     ','OF')
     !here the value of the datacode must be kept fixed
     !n(c) nspin=1

     !if (nproc > 1) call MPI_BARRIER(bigdft_mpi%mpi_env%mpi_comm,ierr)

     !in the case of vacuum the pot_ion treatment is as usual
     !otherwise the pot_ion array is set to zero and can be filled with external potentials
     !like the gaussian part of the PSP ad/or external electric fields
     if (pkernel%method /= 'VAC') then
        call f_zero(n1i*n2i*dpbox%n3pi,pot_ion(1))
     else
        if (allocated(pot_add)) then
           call H_potential('D',pkernel,pot_ion,pot_add,ehart,-psoff_tot,.true.,quiet=quiet)
        else
           call H_potential('D',pkernel,pot_ion,pot_ion,ehart,-psoffset,.false.,quiet=quiet)
        end if
     end if
     
     call timing(iproc,'CrtLocPot     ','ON')

     if (check_potion) then
        !if (iproc == 0) write(*,'(1x,a)',advance='no') 'Check the ionic potential...'

        potion_corr = f_malloc0(n1i*n2i*dpbox%n3pi,id='potion_corr')

        !call to_zero(n1i*n2i*n3pi,potion_corr)
        do while(box_next_point(dpbox%bitp))
           x=dpbox%bitp%rxyz(1)!we should use absolute coordinate
           y=dpbox%bitp%rxyz(2)
           z=dpbox%bitp%rxyz(3)
           call sum_erfcr(at%astruct%nat,at%astruct%ntypes,x,y,z,at%astruct%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
           potion_corr(ind)=potion_corr(ind)+potxyz
        end do
        maxdiff=0.0_dp
        !call f_diff(n1i*n2i*n3pi,potion_corr(1),pot_ion(1),maxdiff)

!!$        !calculate pot_ion with an explicit error function to correct in the case of big grid spacings
!!$        !for the moment works only in the isolated BC case
!!$        do i3=1,n3pi
!!$           z=real(i3+i3s-1-nbl3-1,gp)*hzh
!!$           do i2=1,n2i
!!$              y=real(i2-nbl2-1,gp)*hyh
!!$              do i1=1,n1i
!!$                 x=real(i1-nbl1-1,gp)*hxh
!!$                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
!!$                 !if (i1==49 .and. i2==46 .and. i3==44) then
!!$                 call sum_erfcr(at%astruct%nat,at%astruct%ntypes,x,y,z,at%astruct%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
!!$                 !   stop
!!$                 !end if
!!$                 potion_corr(ind)=potion_corr(ind)+potxyz
!!$                 !write(18,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind)
!!$              end do
!!$           end do
!!$        end do

!!$        !then calculate the maximum difference in the sup norm
!!$        maxdiff=0.0_wp
!!$        do i3=1,n3pi
!!$           do i2=1,n2i
!!$              do i1=1,n1i
!!$                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
!!$                 maxdiff=max(maxdiff,abs(potion_corr(ind)-pot_ion(ind)))
!!$                 !write(17,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind),maxdiff
!!$              end do
!!$           end do
!!$        end do

        if (pkernel%mpi_env%nproc > 1) then
           call fmpi_allreduce(maxdiff,1,FMPI_MAX,comm=pkernel%mpi_env%mpi_comm)
        end if

        if (iproc == 0) call yaml_map('Check the ionic potential',maxdiff,fmt='(1pe24.17)')

        stop

        call f_free(potion_corr)
     end if

  end if

  if (allocated(pot_add)) call f_free(pot_add)

  !!-  !calculate the value of the offset to be put
  !!-  tt_tot=0.d0
  !!-  do ind=1,n1i*n2i*n3i
  !!-     tt_tot=tt_tot+pot_ion(ind)
  !!-  end do
  !!-  print *,'previous offset',tt_tot*hxh*hyh*hzh

  if (dpbox%n3pi > 0) then

     !Determine the maximal bounds for mpx, mpy, mpz (1D-integral)
     call mp_range(at%multipole_preserving,at%mp_isf,at%astruct%nat,&
         hxh,hyh,hzh,maxval(at%psppar(0,0,:)),mpnx,mpny,mpnz)
     !Separable function: do 1-D integrals before and store it.
     mpx = f_malloc( (/ 0 .to. mpnx /),id='mpx')
     mpy = f_malloc( (/ 0 .to. mpny /),id='mpy')
     mpz = f_malloc( (/ 0 .to. mpnz /),id='mpz')

     ! Only for HGH pseudos
     atit = atoms_iter(at%astruct)
     do while(atoms_iter_next(atit))

        !!-     do iat=1,at%astruct%nat
        !!-        ityp=at%astruct%iatype(iat)
        !!-        rx=rxyz(1,iat)
        !!-        ry=rxyz(2,iat)
        !!-        rz=rxyz(3,iat)

        rx=rxyz(1,atit%iat)
        ry=rxyz(2,atit%iat)
        rz=rxyz(3,atit%iat)

        rloc=at%psppar(0,0,atit%ityp)
        rlocinvsq=1.0_gp/rloc**2
        rlocinv2sq=0.5_gp/rloc**2
        cutoff=10.0_gp*rloc

        if (at%multipole_preserving) then
           !We want to have a good accuracy of the last point rloc*10
           !cutoff=cutoff+max(hxh,hyh,hzh)*real(16,kind=gp)
           cutoff=cutoff+max(hxh,hyh,hzh)*real(at%mp_isf,kind=gp)
        end if

        !Separable function: do 1-D integrals before and store it.
        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)
        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)

        do i1=isx,iex
          mpx(i1-isx) = mp_exp(hxh,rx,rlocinv2sq,i1,0,at%multipole_preserving)
        end do
        do i2=isy,iey
          mpy(i2-isy) = mp_exp(hyh,ry,rlocinv2sq,i2,0,at%multipole_preserving)
        end do
        do i3=isz,iez
          mpz(i3-isz) = mp_exp(hzh,rz,rlocinv2sq,i3,0,at%multipole_preserving)
        end do

        if (at%npspcode(atit%ityp) == PSPCODE_GTH .or. &
             & at%npspcode(atit%ityp) == PSPCODE_HGH .or. &
             & at%npspcode(atit%ityp) == PSPCODE_HGH_K .or. &
             & at%npspcode(atit%ityp) == PSPCODE_HGH_K_NLCC) then

           ! Add the remaining local terms of Eq. (9) in JCP 129, 014109(2008)

           ! Determine the number of local terms
           nloc=0
           do iloc=1,4
              !!-              if (at%psppar(0,iloc,ityp) /= 0.d0) nloc=iloc
              if (at%psppar(0,iloc,atit%ityp) /= 0.d0) nloc=iloc
           enddo

           !do not add the local part for the vacancy
           if (nloc /= 0) then

              do i3=isz,iez
                 !call ind_positions(perz,i3,n3,j3,goz)
                 call ind_positions_new(perz,i3,n3i,j3,goz)
                 j3=j3+nbl3+1
                 indj3=(j3-i3s)*n1i*n2i
                 if (goz .and. j3 >= i3s .and. j3 <=  i3s+n3pi-1) then
                    zp = mpz(i3-isz)
                    if (abs(zp) < mp_tiny) cycle
                    z=real(i3,gp)*hzh-rz
                    zsq=z**2
                    do i2=isy,iey
                       !call ind_positions(pery,i2,n2,j2,goy)
                       call ind_positions_new(pery,i2,n2i,j2,goy)
                       indj23=1+nbl1+(j2+nbl2)*n1i+indj3
                       if (goy) then
                          yp = zp*mpy(i2-isy)
                          if (abs(yp) < mp_tiny) cycle
                          y=real(i2,gp)*hyh-ry
                          yzsq=y**2+zsq
                          do i1=isx,iex
                             !call ind_positions(perx,i1,n1,j1,gox)
                             call ind_positions_new(perx,i1,n1i,j1,gox)
                             if (gox) then
                                xp = yp*mpx(i1-isx)
                                if (abs(xp) < mp_tiny) cycle
                                x=real(i1,gp)*hxh-rx
                                r2=x**2+yzsq
                                arg=r2*rlocinvsq
                                tt=at%psppar(0,nloc,atit%ityp)
                                do iloc=nloc-1,1,-1
                                   tt=arg*tt+at%psppar(0,iloc,atit%ityp)
                                enddo
                                ind=j1+indj23
                                pot_ion(ind)=pot_ion(ind)+xp*tt
                             end if
                          enddo
                       end if
                    enddo
                 end if
              end do
           end if !nloc
        end if ! at%npspcode(iat) /= PSPCODE_PAW
     end do !iat

     !De-allocate the 1D temporary arrays for separability
     call f_free(mpx,mpy,mpz)

     !debug exit

     if (htoobig) then
        !add to pot_ion an explicit error function to correct in the case of big grid spacing
        !for the moment works only in the isolated BC case
        !!-        do i3=1,n3pi
        !!-           z=real(i3+i3s-1-nbl3-1,gp)*hzh
        !!-           do i2=1,n2i
        !!-              y=real(i2-nbl2-1,gp)*hyh
        !!-              do i1=1,n1i
        !!-                 x=real(i1-nbl1-1,gp)*hxh
        !!-                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
        !!-                 call sum_erfcr(at%astruct%nat,at%astruct%ntypes,x,y,z,at%astruct%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
        !!-                 pot_ion(ind)=pot_ion(ind)+potxyz
        !!-              end do
        !!-           end do
        !!-        end do

        boxit = dpbox_iter(dpbox,DPB_POT_ION)
        do while(dpbox_iter_next(boxit))
           call sum_erfcr(at%astruct%nat,at%astruct%ntypes, &
                &         boxit%x,boxit%y,boxit%z, &
                &         at%astruct%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
           pot_ion(boxit%ind) = pot_ion(boxit%ind) + potxyz
        end do

     end if

  end if

  !!-  !calculate the value of the offset to be put
  !!-  tt_tot=0.d0
  !!-  do ind=1,n1i*n2i*n3i
  !!-     tt_tot=tt_tot+pot_ion(ind)
  !!-  end do
  !!-  print *,'actual offset',tt_tot*hxh*hyh*hzh

  !use rhopotential to calculate the potential from a constant electric field along y direction
  if (any(elecfield(1:3) /= 0.0_gp)) then
     !constant electric field allowed only for surface and free BC
     if (cell_geocode(dpbox%mesh) == 'P') then
        !if (iproc == 0)
        call f_err_throw('The constant electric field is not allowed for Fully Periodic BC.', &
             err_name='BIGDFT_RUNTIME_ERROR')
        !constant electric field allowed for surface BC only normal to the surface
     elseif (cell_geocode(dpbox%mesh) == 'S' .and. (elecfield(1) /= 0.0_gp .or. elecfield(3) /= 0.0_gp) ) then
        !if (iproc == 0)
        call f_err_throw('Only normal constant electric field (Ex=Ez=0) is allowed for Surface BC.', &
             err_name='BIGDFT_RUNTIME_ERROR')
     end if
     if (verb) call yaml_map('Constant electric field (Ha/Bohr)',elecfield(1:3),fmt='(es10.2)')
     !if (verb) write(*,'(1x,a,"(",es10.2,", ",es10.2,", ",es10.2,") ", a)') &
     !     'Constant electric field ',elecfield(1:3),' Ha/Bohr'
     !or         'Parabolic confining potential: rprb=',elecfield,&
     !           ';  v_conf(r)= 1/(2*rprb**4) * r**2'

     !write or not electric field in a separate file

     call center_of_charge(at,rxyz,center_of_charge_ions)
     do while(box_next_point(dpbox%bitp))
        dpbox%bitp%tmp=dpbox%bitp%rxyz-center_of_charge_ions
        pot_ion(dpbox%bitp%ind)=pot_ion(dpbox%bitp%ind)+&
             elecfield(1)*dpbox%bitp%tmp(1)+elecfield(2)*dpbox%bitp%tmp(2)+elecfield(3)*dpbox%bitp%tmp(3)
     end do

!!$        if (dpbox%n3pi > 0) then
!!$           do i3=1,n3pi
!!$              z=real(i3+i3s-1-nbl3-1,gp)*hzh
!!$              do i2=1,n2i
!!$                 y=real(i2-nbl2-1,gp)*hyh
!!$                 do i1=1,n1i
!!$                    x=real(i1-nbl1-1,gp)*hxh
!!$                    ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
!!$                    pot_ion(ind)=pot_ion(ind)+elecfield(1)*x+elecfield(2)*y+elecfield(3)*z
!!$!                    parabola: these two lines replace the above line comment out the if case and calculate x, z
!!$!                    r2=(x-rx)**2+(y-ry)**2+(z-rz)**2
!!$!                    pot_ion(ind)=pot_ion(ind)+0.5_gp/(elecfield**4)*r2
!!$                 end do
!!$              end do
!!$           end do
!!$        end if

     !        if (efwrite .and. iproc == 0) then
     !           open(unit=17,file='elecpotential_x',status='unknown')
     !           write(17,*) "# x , external electric potential(x,y=0,z=0)"
     !           do i1=nbl1+1,n1i-nbr1-1
     !              x=real(i1-nbl1-1,gp)*hxh
     !              write(17,*)x,-elecfield(1)*x
     !           end do
     !           close(17)
     !           open(unit=17,file='elecpotential_y',status='unknown')
     !           write(17,*) "# y , external electric potential(x=0,y,z=0)"
     !           do i2=nbl2+1,n2i-nbr2-1
     !              y=real(i2-nbl2-1,gp)*hyh
     !              write(17,*)y,-elecfield(2)*y
     !           end do
     !           close(17)
     !           open(unit=17,file='elecpotential_z',status='unknown')
     !           write(17,*) "# z , external electric potential(x=0,y=0,z)"
     !           do i3=1,n3pi
     !              z=real(i3+i3s-1-nbl3-1,gp)*hzh
     !              write(17,*)z,-elecfield(3)*z
     !           end do
     !           close(17)
     !        end if
     !
     !!-     end if
  end if

  !also includes extra treatment for atom-centered parabolic potentials

  !watch if any of the atoms has the "x2" attribute
  yespar=.false.
  atit = atoms_iter(at%astruct)
  do while(atoms_iter_next(atit))
     if ('x2' .notin. atit%attrs) cycle
     !get the value of the parabola and add it to the external potential
     rprbstr=atit%attrs//'x2'
     if (trim(rprbstr) .eqv. 'rprb') then
        call atomic_info(at%nzatom(atit%ityp),at%nelpsp(atit%ityp),rprb=aval)
        aval=0.5_gp/aval**4
     else
        aval=atit%attrs//'x2'
     end if
     if (verb) then
        if (.not. yespar) then
           call yaml_sequence_open('Parabolic confining potentials')
           yespar=.true.
        end if
        call yaml_sequence(advance='no')
        call yaml_mapping_open(flow=.true.)
        call yaml_map('Atom id',atit%iat)
        call yaml_map('Name',atit%name)
        call yaml_map('Prefactor',aval,fmt='(es10.2)')
        call yaml_mapping_close()
     end if
     bitp=dpbox%bitp !shallow copy to avoid intent(out) in dpbox
     do while(box_next_point(bitp))
        r2=square_gd(dpbox%mesh,bitp%rxyz-rxyz(:,atit%iat))
        pot_ion(bitp%ind)=pot_ion(bitp%ind)+aval*r2
     end do
  end do
  if (yespar) call yaml_sequence_close()

  if (at%multipole_preserving) call finalize_real_space_conversion()

  call timing(iproc,'CrtLocPot     ','OF')
  call f_release_routine()
END SUBROUTINE createIonicPotential


!> Create the ionic potential plus the terms describing the environmental information
subroutine external_potential(iproc,verb,at,rxyz,&
     elecfield,dpbox,pkernel,pot_ion,rho_ion,psoffset)
  use module_base
  use module_types
  use yaml_output
  use m_paw_numeric, only: paw_splint
  use multipole_preserving
  use module_atoms
  use module_dpbox
  !  use module_interfaces, only: mp_calculate
  !  use module_interfaces, except_this_one => createIonicPotential
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use abi_interfaces_numeric, only: abi_derf_ab
  use public_enums, only: PSPCODE_PAW
  use bounds, only: ext_buffers
  use box
  implicit none

  !Arguments
  !!-  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: iproc
  !!-  integer, intent(in) :: n1,n2,n3
  logical, intent(in) :: verb
  real(gp), intent(in) :: psoffset
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3), intent(in) :: elecfield
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  type(denspot_distribution), intent(inout) :: dpbox
  type(coulomb_operator), intent(inout) :: pkernel
  real(wp), dimension(*), intent(inout) :: pot_ion
  real(dp), dimension(*), intent(out) :: rho_ion
  !Local variables
  real(gp), parameter :: mp_tiny = 1.e-30_gp
  character(len = 3) :: quiet
  logical :: perx,pery,perz,gox,goy,goz
  logical :: htoobig=.false.,check_potion=.false.
  integer :: i1,i2,i3,ierr !n(c) nspin
  integer :: nloc,iloc
  integer  :: i3s,n3pi,nbl1,nbr1,nbl2,nbl3,nbr2,nbr3
  real(kind=8) :: raux1(1),rr1(1)
  integer :: iex,iey,iez,ind,indj3,indj23,isx,isy,isz,j1,j2,j3
  integer :: n1i,n2i,n3i,mpnx,mpny,mpnz
  real(gp) :: hxh,hyh,hzh
  real(gp) :: rloc,charge,cutoff,r2,arg,xp,tt,rx,ry,rz
  real(gp) :: tt_tot,potxyz
  real(gp) :: raux2,rlocinvsq,rlocinv2sq
  real(gp) :: x,y,z,yp,zp,zsq,yzsq,rholeaked,rholeaked_tot
  !  real(gp), dimension(1) :: raux,rr
  real(wp) :: maxdiff
  real(gp) :: ehart
  logical, dimension(3) :: peri
  real(dp), dimension(2) :: charges_mpi
  real(dp), dimension(:), allocatable :: potion_corr
  logical, parameter :: pawErfCorrection = .true.
  real(gp), dimension(:), allocatable  :: mpx,mpy,mpz
  !real(dp), dimension(:), allocatable :: den_aux
  type(atoms_iterator) :: atit
  type(dpbox_iterator) :: boxit
  real(gp), dimension(3) :: center_of_charge_ions

  call f_routine(id='createIonicPotential')
  call timing(iproc,'CrtLocPot     ','ON')

  !Initialize the work arrays needed to integrate with isf
  !Determine the number of points depending on the min rloc
  if (at%multipole_preserving) &
     call initialize_real_space_conversion(isf_m=at%mp_isf,rlocs=at%psppar(0,0,:))

  ! Aliasing
  hxh = dpbox%mesh%hgrids(1)
  hyh = dpbox%mesh%hgrids(2)
  hzh = dpbox%mesh%hgrids(3)
  n1i = dpbox%mesh%ndims(1)
  n2i = dpbox%mesh%ndims(2)
  n3i = dpbox%mesh%ndims(3)
  i3s = dpbox%i3s+dpbox%i3xcsh
  n3pi = dpbox%n3pi

  ! Ionic charge (must be calculated for the PS active processes)
  rholeaked=0.d0
  ! Ionic energy (can be calculated for all the processors)

  !Creates charge density arising from the ionic PSP cores
  call f_zero(n1i*n2i*dpbox%n3pi,pot_ion(1))

  !conditions for periodicity in the three directions
  peri=cell_periodic_dims(dpbox%mesh)
  perx=peri(1)!(dpbox%geocode /= 'F')
  pery=peri(2)!(dpbox%geocode == 'P')
  perz=peri(3)!(dpbox%geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  if (dpbox%n3pi >0 .and. .not. htoobig) then

    !Determine the maximal bounds for mpx, mpy, mpz (1D-integral)
    call mp_range(at%multipole_preserving,at%mp_isf,at%astruct%nat,&
         hxh,hyh,hzh,maxval(at%psppar(0,0,:)),mpnx,mpny,mpnz)
     !Separable function: do 1-D integrals before and store it.
     mpx = f_malloc( (/ 0 .to. mpnx /),id='mpx')
     mpy = f_malloc( (/ 0 .to. mpny /),id='mpy')
     mpz = f_malloc( (/ 0 .to. mpnz /),id='mpz')

     atit = atoms_iter(at%astruct)
     do while(atoms_iter_next(atit))

        rx=rxyz(1,atit%iat)
        ry=rxyz(2,atit%iat)
        rz=rxyz(3,atit%iat)

        rloc=at%psppar(0,0,atit%ityp)
        rlocinv2sq=0.5_gp/rloc**2
        charge=real(at%nelpsp(atit%ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)

        !cutoff of the range
        cutoff=10.0_gp*rloc
        if (at%multipole_preserving) then
           !We want to have a good accuracy of the last point rloc*10
           cutoff=cutoff+max(hxh,hyh,hzh)*real(at%mp_isf,kind=gp)
        end if

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)

        do i1=isx,iex
           mpx(i1-isx) = mp_exp(hxh,rx,rlocinv2sq,i1,0,at%multipole_preserving)
        end do
        do i2=isy,iey
           mpy(i2-isy) = mp_exp(hyh,ry,rlocinv2sq,i2,0,at%multipole_preserving)
        end do
        do i3=isz,iez
           mpz(i3-isz) = mp_exp(hzh,rz,rlocinv2sq,i3,0,at%multipole_preserving)
        end do

        if ( .not. any(at%npspcode == PSPCODE_PAW) ) then
           !Calculate Ionic Density using HGH parameters.
           !Eq. 1.104, T. Deutsch and L. Genovese, JDN. 12, 2011
           do i3=isz,iez
              zp = mpz(i3-isz)
              !call ind_positions(perz,i3,n3,j3,goz)
              call ind_positions_new(perz,i3,n3i,j3,goz)
              j3=j3+nbl3+1
              if ( goz .and. (j3<i3s.or.j3>i3s+n3pi-1) ) cycle
              indj3=(j3-i3s)*n1i*n2i
              do i2=isy,iey
                 yp = zp*mpy(i2-isy)
                 !call ind_positions(pery,i2,n2,j2,goy)
                 call ind_positions_new(pery,i2,n2i,j2,goy)
                 if (goz.and.(.not.goy)) cycle
                 indj23=1+nbl1+(j2+nbl2)*n1i+indj3
                 do i1=isx,iex
                    xp = yp*mpx(i1-isx)
                    !call ind_positions(perx,i1,n1,j1,gox)
                    call ind_positions_new(perx,i1,n1i,j1,gox)
                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1 .and. goy .and. gox) then
                       ind=j1+indj23
                       pot_ion(ind)=pot_ion(ind)-xp*charge
                       !write(*,'(4(i0,1x),2(1pe20.10))') i1,i2,i3,ind,xp,pot_ion(ind)
                    else if (.not. goz ) then
                       rholeaked=rholeaked+xp*charge
                    endif
                 enddo
              enddo
           enddo
        else
           !Calculate Ionic Density using splines, PAW case
           !r2paw=at%pawtab(atit%ityp)%rpaw**2
           do i3=isz,iez
              zp = mpz(i3-isz)
              if (abs(zp) < mp_tiny) cycle
              !call ind_positions(perz,i3,n3,j3,goz)
              call ind_positions_new(perz,i3,n3i,j3,goz)
              j3=j3+nbl3+1
              indj3=(j3-i3s)*n1i*n2i
              z=real(i3,gp)*hzh-rz
              zsq=z**2
              do i2=isy,iey
                 yp = zp*mpy(i2-isy)
                 if (abs(yp) < mp_tiny) cycle
                 !call ind_positions(pery,i2,n2,j2,goy)
                 call ind_positions_new(pery,i2,n2i,j2,goy)
                 indj23=1+nbl1+(j2+nbl2)*n1i+indj3
                 y=real(i2,gp)*hyh-ry
                 yzsq=y**2+zsq
                 do i1=isx,iex
                    xp = yp*mpx(i1-isx)
                    if (abs(xp) < mp_tiny) cycle
                    !call ind_positions(perx,i1,n1,j1,gox)
                    call ind_positions_new(perx,i1,n1i,j1,gox)
                    x=real(i1,gp)*hxh-rx
                    r2=x**2+yzsq
                    !if(r2>r2paw) cycle
                    if(.not. pawErfCorrection) then
                       !This converges very slowly
                       rr1(1)=sqrt(r2)
                       call paw_splint(at%pawtab(atit%ityp)%wvl%rholoc%msz, &
                            & at%pawtab(atit%ityp)%wvl%rholoc%rad, &
                            & at%pawtab(atit%ityp)%wvl%rholoc%d(:,1), &
                            & at%pawtab(atit%ityp)%wvl%rholoc%d(:,2), &
                            & 1,rr1,raux1,ierr)
                    else
                       !Take the HGH form for rho_L (long range)
                       raux1(1)=-xp*charge
                    end if
                    !raux=-4.d0**(3.0d0/2.0d0)*exp(-4.d0*pi*r2)

                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                       ind=j1+indj23
                       pot_ion(ind)=pot_ion(ind)+raux1(1)
                    else if (.not. goz) then
                       rholeaked=rholeaked-raux1(1)
                    endif
                 enddo
              enddo
           enddo
        end if
     enddo
     !De-allocate for multipole preserving
     call f_free(mpx,mpy,mpz)
  end if

  ! Check
  tt=0.d0
  do j3=1,n3pi
     indj3=(j3-1)*n1i*n2i
     do i2= -nbl2,n2i-nbl2-1!2*n2+1+nbr2
        indj23=1+nbl1+(i2+nbl2)*n1i+indj3
        do i1= -nbl1,n1i-nbl1-1!2*n1+1+nbr1
           ind=i1+indj23
           tt=tt+pot_ion(ind)
        enddo
     enddo
  enddo

  tt=tt*hxh*hyh*hzh
  rholeaked=rholeaked*hxh*hyh*hzh

  if (pkernel%method /= 'VAC') &
       call f_memcpy(n=n1i*n2i*dpbox%n3pi,src=pot_ion(1),dest=rho_ion(1))

  if (pkernel%mpi_env%nproc > 1) then
     charges_mpi(1)=tt
     charges_mpi(2)=rholeaked

     !Reduce from all mpi proc
     call fmpi_allreduce(charges_mpi,FMPI_SUM,comm=pkernel%mpi_env%mpi_comm)

     tt_tot=charges_mpi(1)
     rholeaked_tot=charges_mpi(2)
  else
     tt_tot=tt
     rholeaked_tot=rholeaked
  end if

  if (verb) then
     call yaml_comment('Ionic Potential Creation',hfill='-')
     call yaml_map('Total ionic charge',tt_tot,fmt='(f26.12)')
     if (rholeaked_tot /= 0.0_gp) call yaml_map('Leaked charge',rholeaked_tot,fmt='(1pe10.3)')
     quiet = "no "
  else
     quiet = "yes"
  end if

  if (.not. htoobig) then
     call timing(iproc,'CrtLocPot     ','OF')
     if (pkernel%method /= 'VAC') then
        call f_zero(n1i*n2i*dpbox%n3pi,pot_ion(1))
     else
        call H_potential('D',pkernel,pot_ion,pot_ion,ehart,-psoffset,.false.,quiet=quiet)
     end if
     call timing(iproc,'CrtLocPot     ','ON')

     if (check_potion) then
        potion_corr = f_malloc0(n1i*n2i*dpbox%n3pi,id='potion_corr')
        !calculate pot_ion with an explicit error function to correct in the case of big grid spacings
        !for the moment works only in the isolated BC case
        do i3=1,n3pi
           z=real(i3+i3s-1-nbl3-1,gp)*hzh
           do i2=1,n2i
              y=real(i2-nbl2-1,gp)*hyh
              do i1=1,n1i
                 x=real(i1-nbl1-1,gp)*hxh
                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                 call sum_erfcr(at%astruct%nat,at%astruct%ntypes,x,y,z,at%astruct%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
                 potion_corr(ind)=potion_corr(ind)+potxyz
              end do
           end do
        end do
        !then calculate the maximum difference in the sup norm
        maxdiff=0.0_wp
        do i3=1,n3pi
           do i2=1,n2i
              do i1=1,n1i
                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                 maxdiff=max(maxdiff,abs(potion_corr(ind)-pot_ion(ind)))
                 !write(17,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind),maxdiff
              end do
           end do
        end do

        if (pkernel%mpi_env%nproc > 1) then
           call fmpi_allreduce(maxdiff,1,FMPI_MAX,comm=pkernel%mpi_env%mpi_comm)
        end if

        if (iproc == 0) call yaml_map('Check the ionic potential',maxdiff,fmt='(1pe24.17)')
        stop
        call f_free(potion_corr)
     end if
  end if

  if (dpbox%n3pi > 0) then

     !Determine the maximal bounds for mpx, mpy, mpz (1D-integral)
     call mp_range(at%multipole_preserving,at%mp_isf,at%astruct%nat,&
          hxh,hyh,hzh,maxval(at%psppar(0,0,:)),mpnx,mpny,mpnz)
     !Separable function: do 1-D integrals before and store it.
     mpx = f_malloc( (/ 0 .to. mpnx /),id='mpx')
     mpy = f_malloc( (/ 0 .to. mpny /),id='mpy')
     mpz = f_malloc( (/ 0 .to. mpnz /),id='mpz')

     ! Only for HGH pseudos
     atit = atoms_iter(at%astruct)
     do while(atoms_iter_next(atit))
        rx=rxyz(1,atit%iat)
        ry=rxyz(2,atit%iat)
        rz=rxyz(3,atit%iat)

        rloc=at%psppar(0,0,atit%ityp)
        rlocinvsq=1.0_gp/rloc**2
        rlocinv2sq=0.5_gp/rloc**2
        cutoff=10.0_gp*rloc

        if (at%multipole_preserving) then
           !We want to have a good accuracy of the last point rloc*10
           !cutoff=cutoff+max(hxh,hyh,hzh)*real(16,kind=gp)
           cutoff=cutoff+max(hxh,hyh,hzh)*real(at%mp_isf,kind=gp)
        end if
        !Separable function: do 1-D integrals before and store it.
        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)
        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)

        do i1=isx,iex
           mpx(i1-isx) = mp_exp(hxh,rx,rlocinv2sq,i1,0,at%multipole_preserving)
        end do
        do i2=isy,iey
           mpy(i2-isy) = mp_exp(hyh,ry,rlocinv2sq,i2,0,at%multipole_preserving)
        end do
        do i3=isz,iez
           mpz(i3-isz) = mp_exp(hzh,rz,rlocinv2sq,i3,0,at%multipole_preserving)
        end do

        if( at%npspcode(atit%ityp) /= PSPCODE_PAW) then

           ! Add the remaining local terms of Eq. (9) in JCP 129, 014109(2008)

           ! Determine the number of local terms
           nloc=0
           do iloc=1,4
              if (at%psppar(0,iloc,atit%ityp) /= 0.d0) nloc=iloc
           enddo

           !do not add the local part for the vacancy
           if (nloc /= 0) then
              do i3=isz,iez
                 call ind_positions_new(perz,i3,n3i,j3,goz)
                 j3=j3+nbl3+1
                 indj3=(j3-i3s)*n1i*n2i
                 if (goz .and. j3 >= i3s .and. j3 <=  i3s+n3pi-1) then
                    zp = mpz(i3-isz)
                    if (abs(zp) < mp_tiny) cycle
                    z=real(i3,gp)*hzh-rz
                    zsq=z**2
                    do i2=isy,iey
                       !call ind_positions(pery,i2,n2,j2,goy)
                       call ind_positions_new(pery,i2,n2i,j2,goy)
                       indj23=1+nbl1+(j2+nbl2)*n1i+indj3
                       if (goy) then
                          yp = zp*mpy(i2-isy)
                          if (abs(yp) < mp_tiny) cycle
                          y=real(i2,gp)*hyh-ry
                          yzsq=y**2+zsq
                          do i1=isx,iex
                             !call ind_positions(perx,i1,n1,j1,gox)
                             call ind_positions_new(perx,i1,n1i,j1,gox)
                             if (gox) then
                                xp = yp*mpx(i1-isx)
                                if (abs(xp) < mp_tiny) cycle
                                x=real(i1,gp)*hxh-rx
                                r2=x**2+yzsq
                                arg=r2*rlocinvsq
                                tt=at%psppar(0,nloc,atit%ityp)
                                do iloc=nloc-1,1,-1
                                   tt=arg*tt+at%psppar(0,iloc,atit%ityp)
                                enddo
                                ind=j1+indj23
                                pot_ion(ind)=pot_ion(ind)+xp*tt
                             end if
                          enddo
                       end if
                    enddo
                 end if
              end do
           end if
        else if (pawErfCorrection) then
           ! For PAW, add V^PAW-V_L^HGH
           charge=real(at%nelpsp(atit%ityp),gp)
           !!-           charge=real(at%nelpsp(ityp),gp)
           do i3=isz,iez
              z=real(i3,gp)*hzh-rz
              !call ind_positions(perz,i3,n3,j3,goz)
              call ind_positions_new(perz,i3,n3i,j3,goz)
              j3=j3+nbl3+1
              indj3=(j3-i3s)*n1i*n2i
              zsq=z**2
              if (goz .and. j3 >= i3s .and. j3 <=  i3s+n3pi-1) then
                 do i2=isy,iey
                    y=real(i2,gp)*hyh-ry
                    !call ind_positions(pery,i2,n2,j2,goy)
                    call ind_positions_new(pery,i2,n2i,j2,goy)
                    indj23=1+nbl1+(j2+nbl2)*n1i+indj3
                    yzsq=y**2+zsq
                    if (goy) then
                       do i1=isx,iex
                          x=real(i1,gp)*hxh-rx
                          !call ind_positions(perx,i1,n1,j1,gox)
                          call ind_positions_new(perx,i1,n1i,j1,gox)
                          if (gox) then
                             r2=x**2+yzsq
                             rr1(1)=sqrt(r2)
                             !1) V_L^HGH
                             if(rr1(1)>0.01d0) then
                                arg=rr1(1)/(sqrt(2.0)*rloc)
                                call abi_derf_ab(tt,arg)
                                raux2=-charge/rr1(1)*tt
                             else
                                !In this case we deduce the values
                                !from a quadratic interpolation (due to 1/rr factor)
                                call interpol_vloc(rr1(1),rloc,charge,raux2)
                             end if
                             !2) V^PAW from splines
                             call paw_splint(at%pawtab(atit%ityp)%wvl%rholoc%msz, &
                                  & at%pawtab(atit%ityp)%wvl%rholoc%rad, &
                                  & at%pawtab(atit%ityp)%wvl%rholoc%d(:,3), &
                                  & at%pawtab(atit%ityp)%wvl%rholoc%d(:,4), &
                                  & 1,rr1,raux1,ierr)
                             ind=j1+indj23
                             pot_ion(ind)=pot_ion(ind)+raux1(1)-raux2
                          end if
                       enddo
                    end if
                 enddo
              end if
           end do
        end if ! at%npspcode(iat) /= PSPCODE_PAW
     end do !iat

     !De-allocate the 1D temporary arrays for separability
     call f_free(mpx,mpy,mpz)

     !debug exit

     if (htoobig) then
        !add to pot_ion an explicit error function to correct in the case of big grid spacing
        !for the moment works only in the isolated BC case
        !!-        do i3=1,n3pi
        !!-           z=real(i3+i3s-1-nbl3-1,gp)*hzh
        !!-           do i2=1,n2i
        !!-              y=real(i2-nbl2-1,gp)*hyh
        !!-              do i1=1,n1i
        !!-                 x=real(i1-nbl1-1,gp)*hxh
        !!-                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
        !!-                 call sum_erfcr(at%astruct%nat,at%astruct%ntypes,x,y,z,at%astruct%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
        !!-                 pot_ion(ind)=pot_ion(ind)+potxyz
        !!-              end do
        !!-           end do
        !!-        end do

        boxit = dpbox_iter(dpbox,DPB_POT_ION)
        do while(dpbox_iter_next(boxit))
           call sum_erfcr(at%astruct%nat,at%astruct%ntypes, &
                &         boxit%x,boxit%y,boxit%z, &
                &         at%astruct%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
           pot_ion(boxit%ind) = pot_ion(boxit%ind) + potxyz
        end do

     end if

  end if

  !use rhopotential to calculate the potential from a constant electric field along y direction
  if (any(elecfield(1:3) /= 0.0_gp)) then
     !constant electric field allowed only for surface and free BC
     if (cell_geocode(dpbox%mesh) == 'P') then
        call f_err_throw('The constant electric field is not allowed for Fully Periodic BC.', &
             err_name='BIGDFT_RUNTIME_ERROR')
        !constant electric field allowed for surface BC only normal to the surface
     elseif (cell_geocode(dpbox%mesh) == 'S' .and. (elecfield(1) /= 0.0_gp .or. elecfield(3) /= 0.0_gp) ) then
        call f_err_throw('Only normal constant electric field (Ex=Ez=0) is allowed for Surface BC.', &
             err_name='BIGDFT_RUNTIME_ERROR')
     end if
     if (verb) call yaml_map('Constant electric field (Ha/Bohr)',elecfield(1:3),fmt='(es10.2)')
     !if (verb) write(*,'(1x,a,"(",es10.2,", ",es10.2,", ",es10.2,") ", a)') &
     !     'Constant electric field ',elecfield(1:3),' Ha/Bohr'
     !or         'Parabolic confining potential: rprb=',elecfield,&
     !           ';  v_conf(r)= 1/(2*rprb**4) * r**2'

     !the iterator here is on the potential distribution
     call center_of_charge(at,rxyz,center_of_charge_ions)
     do while(box_next_point(dpbox%bitp))
        dpbox%bitp%tmp=dpbox%bitp%rxyz-center_of_charge_ions
        pot_ion(dpbox%bitp%ind)=pot_ion(dpbox%bitp%ind)+&
             elecfield(1)*dpbox%bitp%tmp(1)+elecfield(2)*dpbox%bitp%tmp(2)+elecfield(3)*dpbox%bitp%tmp(3)

     end do
!!$     !write or not electric field in a separate file
!!$     if (dpbox%n3pi > 0) then
!!$        do i3=1,n3pi
!!$           z=real(i3+i3s-1-nbl3-1,gp)*hzh
!!$           do i2=1,n2i
!!$              y=real(i2-nbl2-1,gp)*hyh
!!$              do i1=1,n1i
!!$                 x=real(i1-nbl1-1,gp)*hxh
!!$                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
!!$                 pot_ion(ind)=pot_ion(ind)+elecfield(1)*x+elecfield(2)*y+elecfield(3)*z
!!$                 !                    parabola: these two lines replace the above line comment out the if case and calculate x, z
!!$                 !                    r2=(x-rx)**2+(y-ry)**2+(z-rz)**2
!!$                 !                    pot_ion(ind)=pot_ion(ind)+0.5_gp/(elecfield**4)*r2
!!$              end do
!!$           end do
!!$        end do
!!$     end if
  end if

  if (at%multipole_preserving) call finalize_real_space_conversion()

  call timing(iproc,'CrtLocPot     ','OF')
  call f_release_routine()

contains

  ! We use a quadratic interpolation to get vloc(x)
  ! useful for small values of x
  SUBROUTINE interpol_vloc(xx,rloc,charge,yy)
    implicit none
    real(dp),intent(in)  :: xx,rloc,charge
    real(dp),intent(out) :: yy
    !   local variables
    real(dp)::l0,l1,l2,x0,x1,x2,y0,y1,y2

    !   Find 3 points (x0,y0), (x1,y1), (x2,y2).
    x0=0.01d0; x1=0.02d0; x2=0.03d0
    call calcVloc(y0,x0,rloc,charge)
    call calcVloc(y1,x1,rloc,charge)
    call calcVloc(y2,x2,rloc,charge)

    !   Find a polynomial of the form:
    !   P(x)=y0L0(x) + y1L1(x) + y2L2(x)

    !   L0(x) = (x-x1)(x-x2)/((x0-x1)(x0-x2))
    l0=(xx-x1)*(xx-x2)/((x0-x1)*(x0-x2))
    !   L1(x) = (x-x0)(x-x2)/((x1-x0)(x1-x2))
    l1=(xx-x0)*(xx-x2)/((x1-x0)*(x1-x2))
    !   L2(x) = (x-x0)(x-x1)/((x2-x0)(x2-x1))
    l2=(xx-x0)*(xx-x1)/((x2-x0)*(x2-x1))

    yy=y0*l0+y1*l1+y2*l2

  END SUBROUTINE interpol_vloc

  subroutine calcVloc(yy,xx,rloc,Z)
    use abi_interfaces_numeric, only: abi_derf_ab
    implicit none
    !Arguments
    real(dp),intent(in)  :: xx,rloc,Z
    real(dp),intent(out) :: yy
    !Local variables
    !integer, parameter   :: dp = kind(1.0d0) !< double precision
    real(dp):: arg,tt

    arg=xx/(sqrt(2.0_dp)*rloc)
    call abi_derf_ab(tt,arg)
    yy=-Z/xx*tt

  end subroutine calcVloc

END SUBROUTINE external_potential



!> Determine the index in which the potential must be inserted, following the BC
!! Determine also whether the index is inside or outside the box for free BC
!!- subroutine ind_positions(periodic,i,n,j,go)
!!-   implicit none
!!-   logical, intent(in) :: periodic
!!-   integer, intent(in) :: i,n
!!-   logical, intent(out) :: go
!!-   integer, intent(out) :: j
!!-
!!-   if (periodic) then
!!-      go=.true.
!!-      j=modulo(i,2*n+2)
!!-   else
!!-      j=i
!!-      if (i >= -14 .and. i <= 2*n+16) then
!!-         go=.true.
!!-      else
!!-         go=.false.
!!-      end if
!!-   end if
!!-
!!- END SUBROUTINE ind_positions


!> Determine the index in which the potential must be inserted, following the BC
!! Determine also whether the index is inside or outside the box for free BC
subroutine ind_positions_new(periodic,i,ni,j,go)
  implicit none
  logical, intent(in) :: periodic
  integer, intent(in) :: i,ni
  logical, intent(out) :: go
  integer, intent(out) :: j

  if (periodic) then
     go=.true.
     j=modulo(i,ni)
  else
     j=i
     if (i >= -14 .and. i <= ni-15) then
        go=.true.
     else
        go=.false.
     end if
  end if

END SUBROUTINE ind_positions_new


!> Calculate the 1D separable integral for the exponential function used in the local potential
!!- subroutine mp_calculate(rx,ry,rz,hxh,hyh,hzh,cutoff,rlocinv2sq,mp,mpx,mpy,mpz)
!!-   use module_base
!!-   use gaussians, only: mp_exp
!!-   !Arguments
!!-   real(gp), intent(in) :: rx,ry,rz,hxh,hyh,hzh,cutoff,rlocinv2sq
!!-   logical, intent(in) :: mp
!!-   real(gp), dimension(:), allocatable, intent(out) :: mpx,mpy,mpz
!!-
!!-   !Local variable
!!-   integer :: isx,iex,isy,iey,isz,iez
!!-   logical :: lc
!!-
!!-   isx=floor((rx-cutoff)/hxh)
!!-   isy=floor((ry-cutoff)/hyh)
!!-   isz=floor((rz-cutoff)/hzh)
!!-
!!-   iex=ceiling((rx+cutoff)/hxh)
!!-   iey=ceiling((ry+cutoff)/hyh)
!!-   iez=ceiling((rz+cutoff)/hzh)
!!-
!!-   mpx = f_malloc( (/ isx.to.iex /),id='mpx')
!!-   mpy = f_malloc( (/ isy.to.iey /),id='mpy')
!!-   mpz = f_malloc( (/ isz.to.iez /),id='mpz')
!!-
!!- !  lc = .true.
!!-   do i1=isx,iex
!!-      mpx(i1) = mp_exp(hxh,rx,rlocinv2sq,i1,0,mp)
!!- !     if (abs(mpx(i1)) < mp_tiny) then
!!- !        if (lc) istart = i1+1
!!- !     else
!!- !        lc = .false.
!!- !        iend = i1
!!- !     end if
!!-   end do
!!- !  isx = istart
!!- !  iex = iend
!!- !  lc = .true.
!!-   do i2=isy,iey
!!-      mpy(i2) = mp_exp(hyh,ry,rlocinv2sq,i2,0,mp)
!!- !     if (abs(mpy(i2)) < mp_tiny) then
!!- !        if (lc) istart = i2+1
!!- !     else
!!- !        lc = .false.
!!- !        iend = i2
!!- !     end if
!!-   end do
!!- !  isy = istart
!!- !  iey = iend
!!- !  lc = .true.
!!-   do i3=isz,iez
!!-      mpz(i3) = mp_exp(hzh,rz,rlocinv2sq,i3,0,mp)
!!- !     if (abs(mpz(i3)) < mp_tiny) then
!!- !        if (lc) istart = i3+1
!!- !     else
!!- !        lc = .false.
!!- !        iend = i3
!!- !     end if
!!-   end do
!!- !  isz = istart
!!- !  iez = iend
!!-
!!- end subroutine mp_calculate


subroutine sum_erfcr(nat,ntypes,x,y,z,iatype,nelpsp,psppar,rxyz,potxyz)
  use module_base
  use abi_interfaces_numeric, only: abi_derf_ab
  implicit none
  integer, intent(in) :: nat,ntypes
  real(gp) :: x,y,z
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(gp), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(wp), intent(out) :: potxyz
  !local variables
  integer :: iat,ityp
  real(wp) :: charge
  real(gp) :: r,sq2rl,rx,ry,rz,derf_val

  potxyz =0.0_wp

  do iat=1,nat

     ityp=iatype(iat)
     sq2rl=sqrt(2.0_gp)*psppar(0,0,ityp)
     charge=real(nelpsp(ityp),wp)

     rx=rxyz(1,iat)-x
     ry=rxyz(2,iat)-y
     rz=rxyz(3,iat)-z

     r=sqrt(rx**2+ry**2+rz**2)

     if (r == 0.0_gp) then
        potxyz = potxyz - charge*2.0_wp/(sqrt(pi)*real(sq2rl,wp))
     else
        call abi_derf_ab(derf_val,r/sq2rl)
        potxyz = potxyz - charge*real(derf_val/r,wp)
     end if

     !write(*,'(a,1x,i0,3(1pe24.17))')'iat,r,derf_val,derf_val/r',iat,r/sq2rl,derf_val,derf_val/r

  end do

END SUBROUTINE sum_erfcr


!> Read and initialize counter-ions potentials (read psp files)
subroutine CounterIonPotential(iproc,in,shift,dpbox,pkernel,npot_ion,pot_ion)
  use module_base
  use module_types
  !use module_interfaces, except_this_one => CounterIonPotential
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use module_input_dicts
  use public_keys, only: IG_OCCUPATION
  use dictionaries
  use yaml_output
  use multipole_preserving
  use module_atoms
  use module_dpbox
!!$  use multipole, only: gaussian_density
  use bounds, only: ext_buffers
  use pseudopotentials, only: psp_dict_fill_all
  use box, only: cell_periodic_dims
  use gaussians
  implicit none
  !Arguments
  integer, intent(in) :: iproc, npot_ion
  !!-  real(gp), intent(in) :: hxh,hyh,hzh
  real(gp), dimension(3), intent(in) :: shift
  type(input_variables), intent(in) :: in
  !!-  type(grid_dimensions), intent(in) :: grid
  type(denspot_distribution), intent(in) :: dpbox
  type(coulomb_operator), intent(inout) :: pkernel
  real(wp), dimension(npot_ion), intent(inout) :: pot_ion
  !Local variables
  real(gp), parameter :: mp_tiny = 1.e-30_gp
  logical, parameter :: htoobig=.false.,check_potion=.false.,use_iterator=.false.
  logical :: perx,pery,perz
  integer :: j3,nmpx,nmpy,nmpz
  integer :: i1,i2,i3,ityp,nspin,indj3,indj23,n1i,n2i,n3i
  integer :: ind,nbl1,nbr1,nbl2,nbr2,n3pi,nbl3,nbr3,i3s
  real(kind=8) :: cutoff,tt,rx,ry,rz
  real(kind=8) :: x,y,z,rholeaked,rholeaked_tot
  real(kind=8) :: hxh,hyh,hzh,tt_tot,potxyz
  real(wp) :: maxdiff
  real(gp) :: ehart
  type(atoms_data) :: at
  type(dictionary), pointer :: dict
  logical, dimension(3) :: peri
  real(dp), dimension(2) :: charges_mpi
  real(dp), dimension(:), allocatable :: potion_corr
  real(dp), dimension(:), allocatable  :: mpx,mpy,mpz
  !  real(gp), dimension(:,:), allocatable :: radii_cf
  type(atoms_iterator) :: atit
  type(dpbox_iterator) :: boxit
  type(gaussian_real_space) :: g


  call timing(iproc,'CrtLocPot     ','ON')

  n3pi = dpbox%n3pi
  i3s = dpbox%i3s + dpbox%i3xcsh
  hxh = dpbox%mesh%hgrids(1)
  hyh = dpbox%mesh%hgrids(2)
  hzh = dpbox%mesh%hgrids(3)
  n1i=dpbox%mesh%ndims(1)
  n2i=dpbox%mesh%ndims(2)
  n3i=dpbox%mesh%ndims(3)


  if (iproc.eq.0) then
     !write(*,'(1x,a)')&
     !     '--------------------------------------------------- Counter Ionic Potential Creation'
     call yaml_comment('Counter Ionic Potential Creation',hfill='-')
  end if

  at = atoms_data_null()
  !read the positions of the counter ions from file
  call dict_init(dict)
  call astruct_file_merge_to_dict(dict, "posinp", 'posinp_ci')
  call astruct_set_from_dict(dict // "posinp", at%astruct)

  call atoms_file_merge_to_dict(dict)
  do ityp = 1, at%astruct%ntypes, 1
     call psp_dict_fill_all(dict, at%astruct%atomnames(ityp), in%ixc, in%projrad, in%crmult, in%frmult)
  end do
  call psp_dict_analyse(dict, at)
  ! Read associated pseudo files.
  call atomic_data_set_from_dict(dict,IG_OCCUPATION, at, in%nspin)
  call dict_free(dict)

  !read the specifications of the counter ions from pseudopotentials
  !  radii_cf = f_malloc((/ at%astruct%ntypes, 3 /),id='radii_cf')
  !  radii_cf = at%radii_cf
  if (iproc == 0) call print_atomic_variables(at, max(in%hx,in%hy,in%hz), in%ixc)


  !Initialize the work arrays needed to integrate with isf
  !Determine the number of points depending on the min rloc
  if (at%multipole_preserving) &
     call initialize_real_space_conversion(isf_m=at%mp_isf,rlocs=at%psppar(0,0,:))

  ! Ionic charge (must be calculated for the PS active processes)
  rholeaked=0.d0
  ! Ionic energy (can be calculated for all the processors)

  !Creates charge density arising from the ionic PSP cores
  call f_zero(dpbox%mesh%ndims(1)*dpbox%mesh%ndims(2)*dpbox%n3pi,pot_ion(1))


  if (.not. use_iterator) then
     peri=cell_periodic_dims(dpbox%mesh)
     !conditions for periodicity in the three directions
     perx=peri(1)!(dpbox%geocode /= 'F')
     pery=peri(2)!(dpbox%geocode == 'P')
     perz=peri(3)!(dpbox%geocode /= 'F')

     !Calculate external buffers for each direction
     call ext_buffers(perx,nbl1,nbr1)
     call ext_buffers(pery,nbl2,nbr2)
     call ext_buffers(perz,nbl3,nbr3)
  end if

  if (dpbox%n3pi >0 .and. .not. htoobig) then

     !Determine the maximal bounds for mpx, mpy, mpz (1D-integral)
     cutoff=10.0_gp*maxval(at%psppar(0,0,:))
     if (at%multipole_preserving) then
        !We want to have a good accuracy of the last point rloc*10
        cutoff=cutoff+max(hxh,hyh,hzh)*real(at%mp_isf,kind=gp)
     end if
     !Separable function: do 1-D integrals before and store it.
     nmpx = (ceiling(cutoff/hxh) - floor(-cutoff/hxh)) + 1
     nmpy = (ceiling(cutoff/hyh) - floor(-cutoff/hyh)) + 1
     nmpz = (ceiling(cutoff/hzh) - floor(-cutoff/hzh)) + 1
     mpx = f_malloc( (/ 0 .to. nmpx /),id='mpx')
     mpy = f_malloc( (/ 0 .to. nmpy /),id='mpy')
     mpz = f_malloc( (/ 0 .to. nmpz /),id='mpz')

     atit = atoms_iter(at%astruct)
     if (iproc==0) then
        call yaml_sequence_open('Counter ions')
     end if
     do while(atoms_iter_next(atit))

        !!-     do iat=1,at%astruct%nat
        !!-        ityp=at%astruct%iatype(iat)

        !shift the positions of the counter_ion wrt the box
        !!-        rx=at%astruct%rxyz(1,iat)-shift(1)
        !!-        ry=at%astruct%rxyz(2,iat)-shift(2)
        !!-        rz=at%astruct%rxyz(3,iat)-shift(3)
        rx=at%astruct%rxyz(1,atit%iat)-shift(1)
        ry=at%astruct%rxyz(2,atit%iat)-shift(2)
        rz=at%astruct%rxyz(3,atit%iat)-shift(3)

        if (iproc == 0) then
           !write(*,'(1x,a,i6,3(1pe14.7))')'counter ion No. ',atit%iat,rx,ry,rz
           call yaml_sequence(advance='no')
           call yaml_mapping_open(flow=.true.)
           call yaml_map(trim(at%astruct%atomnames(at%astruct%iatype(atit%iat))),(/rx,ry,rz/),fmt='(1es20.12)')
           call yaml_mapping_close(advance='no')
           call yaml_comment(trim(yaml_toa(atit%iat,fmt='(i4.4)')))
        end if

        ! SM: COPY FROM HERE... #############################################################

!!$        call gaussian_density(perx, pery, perz, n1i, n2i, n3i, nbl1, nbl2, nbl3, i3s, n3pi, hxh, hyh, hzh, rx, ry, rz, &
!!$             at%psppar(0,0,atit%ityp), at%nelpsp(atit%ityp), at%multipole_preserving, use_iterator, at%mp_isf, &
!!$             dpbox, nmpx, nmpy, nmpz, mpx, mpy, mpz, npot_ion, pot_ion, rholeaked)

        call atomic_charge_density(g,at,atit)
        call three_dimensional_density(dpbox%bitp,g,-1.0_dp,[rx,ry,rz],pot_ion)

!!!!!-        rloc=at%psppar(0,0,ityp)
!!!!!-        charge=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
!!!        rloc=at%psppar(0,0,atit%ityp)
!!!        rlocinv2sq=0.5_gp/rloc**2
!!!        charge=real(at%nelpsp(atit%ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
!!!
!!!        write(*,*) 'at%nelpsp(atit%ityp),rloc,charge',at%nelpsp(atit%ityp),rloc,charge
!!!
!!!        !cutoff of the range
!!!        cutoff=10.0_gp*rloc
!!!        if (at%multipole_preserving) then
!!!           !We want to have a good accuracy of the last point rloc*10
!!!           !cutoff=cutoff+max(hxh,hyh,hzh)*real(16,kind=gp)
!!!           cutoff=cutoff+max(hxh,hyh,hzh)*real(at%mp_isf,kind=gp)
!!!        end if
!!!        write(*,*) 'cutoff',cutoff
!!!
!!!        if (use_iterator) then
!!!           nbox(1,1)=floor((rx-cutoff)/hxh)
!!!           nbox(1,2)=floor((ry-cutoff)/hyh)
!!!           nbox(1,3)=floor((rz-cutoff)/hzh)
!!!           nbox(2,1)=ceiling((rx+cutoff)/hxh)
!!!           nbox(2,2)=ceiling((ry+cutoff)/hyh)
!!!           nbox(2,3)=ceiling((rz+cutoff)/hzh)
!!!
!!!           !Separable function: do 1-D integrals before and store it.
!!!           !mpx = f_malloc( (/ nbox(1,1).to.nbox(2,1) /),id='mpx')
!!!           !mpy = f_malloc( (/ nbox(1,2).to.nbox(2,2) /),id='mpy')
!!!           !mpz = f_malloc( (/ nbox(1,3).to.nbox(2,3) /),id='mpz')
!!!           do i1=nbox(1,1),nbox(2,1)
!!!              mpx(i1-nbox(1,1)) = mp_exp(hxh,rx,rlocinv2sq,i1,0,at%multipole_preserving)
!!!           end do
!!!           do i2=nbox(1,2),nbox(2,2)
!!!              mpy(i2-nbox(1,2)) = mp_exp(hyh,ry,rlocinv2sq,i2,0,at%multipole_preserving)
!!!           end do
!!!           do i3=nbox(1,3),nbox(2,3)
!!!              mpz(i3-nbox(1,3)) = mp_exp(hzh,rz,rlocinv2sq,i3,0,at%multipole_preserving)
!!!           end do
!!!           boxit = dpbox_iter(dpbox,DPB_POT_ION,nbox)
!!!
!!!
!!!           do while(dpbox_iter_next(boxit))
!!!              xp = mpx(boxit%ibox(1)-nbox(1,1)) * mpy(boxit%ibox(2)-nbox(1,2)) * mpz(boxit%ibox(3)-nbox(1,3))
!!!              pot_ion(boxit%ind) = pot_ion(boxit%ind) - xp*charge
!!!           end do
!!!
!!!        else
!!!           isx=floor((rx-cutoff)/hxh)
!!!           isy=floor((ry-cutoff)/hyh)
!!!           isz=floor((rz-cutoff)/hzh)
!!!
!!!           iex=ceiling((rx+cutoff)/hxh)
!!!           iey=ceiling((ry+cutoff)/hyh)
!!!           iez=ceiling((rz+cutoff)/hzh)
!!!
!!!           !Separable function: do 1-D integrals before and store it.
!!!           !call mp_calculate(rx,ry,rz,hxh,hyh,hzh,cutoff,rlocinv2sq,at%multipole_preserving,mpx,mpy,mpz)
!!!           !mpx = f_malloc( (/ isx.to.iex /),id='mpx')
!!!           !mpy = f_malloc( (/ isy.to.iey /),id='mpy')
!!!           !mpz = f_malloc( (/ isz.to.iez /),id='mpz')
!!!           write(*,*) 'before exps'
!!!           do i1=isx,iex
!!!              mpx(i1-isx) = mp_exp(hxh,rx,rlocinv2sq,i1,0,at%multipole_preserving)
!!!           end do
!!!           do i2=isy,iey
!!!              mpy(i2-isy) = mp_exp(hyh,ry,rlocinv2sq,i2,0,at%multipole_preserving)
!!!           end do
!!!           do i3=isz,iez
!!!              mpz(i3-isz) = mp_exp(hzh,rz,rlocinv2sq,i3,0,at%multipole_preserving)
!!!           end do
!!!           do i3=isz,iez
!!!           write(*,*) 'i3',i3
!!!              zp = mpz(i3-isz)
!!!              if (abs(zp) < mp_tiny) cycle
!!!              !call ind_positions(perz,i3,grid%n3,j3,goz)
!!!              call ind_positions_new(perz,i3,n3i,j3,goz)
!!!              j3=j3+nbl3+1
!!!              do i2=isy,iey
!!!                 yp = zp*mpy(i2-isy)
!!!                 if (abs(yp) < mp_tiny) cycle
!!!                 !call ind_positions(pery,i2,grid%n2,j2,goy)
!!!                 call ind_positions_new(pery,i2,n2i,j2,goy)
!!!                 do i1=isx,iex
!!!                    xp = yp*mpx(i1-isx)
!!!                    if (abs(xp) < mp_tiny) cycle
!!!                    !call ind_positions(perx,i1,grid%n1,j1,gox)
!!!                    call ind_positions_new(perx,i1,n1i,j1,gox)
!!!                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
!!!                       ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s)*n1i*n2i
!!!                       pot_ion(ind)=pot_ion(ind)-xp*charge
!!!                    else if (.not. goz ) then
!!!                       rholeaked=rholeaked+xp*charge
!!!                    endif
!!!                 enddo
!!!              enddo
!!!           enddo
!!!
!!!        end if
!!!
!!!        ! SM: TO HERE... #############################################################

        !!De-allocate for multipole preserving
        !call f_free(mpx,mpy,mpz)

     end do

     if (iproc==0) then
        call yaml_sequence_close()
     end if
     call f_free(mpx,mpy,mpz)

  end if

  ! Check
  tt=0.d0
  if (use_iterator) then
     boxit = dpbox_iter(dpbox,DPB_POT_ION)
     do while(dpbox_iter_next(boxit))
        tt = tt + pot_ion(boxit%ind)
     end do
  else
     do j3=1,n3pi
        indj3=(j3-1)*n1i*n2i
        do i2= -nbl2,n2i-nbl2-1!2*n2+1+nbr2
           indj23=1+nbl1+(i2+nbl2)*n1i+indj3
           do i1= -nbl1,n1i-nbl1-1!2*n1+1+nbr1
              ind=i1+indj23
              tt=tt+pot_ion(ind)
           enddo
        enddo
     enddo
  end if

  tt=tt*hxh*hyh*hzh
  rholeaked=rholeaked*hxh*hyh*hzh

  if (pkernel%mpi_env%nproc > 1) then
     charges_mpi(1)=tt
     charges_mpi(2)=rholeaked

     call fmpi_allreduce(charges_mpi,FMPI_SUM,comm=pkernel%mpi_env%mpi_comm)

     tt_tot=charges_mpi(1)
     rholeaked_tot=charges_mpi(2)
  else
     tt_tot=tt
     rholeaked_tot=rholeaked
  end if

  !!-  if (iproc == 0) write(*,'(1x,a,f26.12,2x,1pe10.3)') &
  !!-       'total ionic charge, leaked charge ',tt_tot,rholeaked_tot
  if (iproc == 0) call yaml_map('total ionic charge',tt_tot)

  if (.not. htoobig) then
     call timing(iproc,'CrtLocPot     ','OF')
     !here the value of the datacode must be kept fixed
     nspin=1

     call H_potential('D',pkernel,pot_ion,pot_ion,ehart,0.0_gp,.false.)

     call timing(iproc,'CrtLocPot     ','ON')

     if (check_potion) then
        !if (iproc == 0) write(*,'(1x,a)',advance='no') 'Check the ionic potential...'
        potion_corr = f_malloc0(dpbox%mesh%ndims(1)*dpbox%mesh%ndims(2)*dpbox%n3pi,id='potion_corr')

        !call to_zero(grid%n1i*grid%n2i*n3pi,potion_corr)

        if (use_iterator) then
           boxit = dpbox_iter(dpbox,DPB_POT_ION)
           do while(dpbox_iter_next(boxit))
              call sum_erfcr(at%astruct%nat,at%astruct%ntypes, &
                   &         boxit%x,boxit%y,boxit%z, &
                   &         at%astruct%iatype,at%nelpsp,at%psppar,at%astruct%rxyz,potxyz)
              potion_corr(boxit%ind) = potion_corr(boxit%ind )+ potxyz
           end do
           !then calculate the maximum difference in the sup norm
           maxdiff=0.0_wp
           boxit = dpbox_iter(dpbox,DPB_POT_ION)
           do while(dpbox_iter_next(boxit))
              maxdiff=max(maxdiff,abs(potion_corr(boxit%ind)-pot_ion(boxit%ind)))
           end do
        else
           !calculate pot_ion with an explicit error function to correct in the case of big grid spacings
           !for the moment works only in the isolated BC case
           do i3=1,n3pi
              z=real(i3+i3s-1-nbl3-1,gp)*hzh
              do i2=1,n2i
                 y=real(i2-nbl2-1,gp)*hyh
                 do i1=1,n1i
                    x=real(i1-nbl1-1,gp)*hxh
                    ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                    !if (i1==49 .and. i2==46 .and. i3==44) then
                    call sum_erfcr(at%astruct%nat,at%astruct%ntypes,x,y,z,at%astruct%iatype,at%nelpsp,at%psppar,&
                         at%astruct%rxyz,potxyz)
                    !   stop
                    !end if
                    potion_corr(ind)=potion_corr(ind)+potxyz
                    !write(18,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind)
                 end do
              end do
           end do
           maxdiff=0.0_wp
           do i3=1,n3pi
              do i2=1,n2i
                 do i1=1,n1i
                    ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                    maxdiff=max(maxdiff,abs(potion_corr(ind)-pot_ion(ind)))
                    !write(17,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind),maxdiff
                 end do
              end do
           end do
        end if


        if (pkernel%mpi_env%nproc > 1) then
           call fmpi_allreduce(maxdiff,1,FMPI_MAX,comm=pkernel%mpi_env%mpi_comm)
        end if

        if (iproc == 0) call yaml_map('Check the ionic potential',maxdiff,fmt='(1pe24.17)')
        !if (iproc == 0) write(*,'(1x,a,1pe24.17)')'...done. MaxDiff=',maxdiff

        stop

        call f_free(potion_corr)

     end if

  end if

  if (dpbox%n3pi > 0 .and. htoobig) then
     if (use_iterator) then
        boxit = dpbox_iter(dpbox,DPB_POT_ION)
        do while(dpbox_iter_next(boxit))
           call sum_erfcr(at%astruct%nat,at%astruct%ntypes, &
                &         boxit%x,boxit%y,boxit%z, &
                &         at%astruct%iatype,at%nelpsp,at%psppar,at%astruct%rxyz,potxyz)
           pot_ion(boxit%ind) = pot_ion(boxit%ind) + potxyz
        end do
     else
        !add to pot_ion an explicit error function to correct in the case of big grid spacing
        !for the moment works only in the isolated BC case
        do i3=1,n3pi
           z=real(i3+i3s-1-nbl3-1,gp)*hzh
           do i2=1,n2i
              y=real(i2-nbl2-1,gp)*hyh
              do i1=1,n1i
                 x=real(i1-nbl1-1,gp)*hxh
                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                 call sum_erfcr(at%astruct%nat,at%astruct%ntypes,x,y,z,at%astruct%iatype,at%nelpsp,at%psppar,at%astruct%rxyz,potxyz)
                 pot_ion(ind)=pot_ion(ind)+potxyz
              end do
           end do
        end do
     end if
  end if

  !deallocations
  call deallocate_atoms_data(at)

  !  call f_free(radii_cf)

  call f_free_ptr(at%astruct%rxyz)


  if (at%multipole_preserving) call finalize_real_space_conversion()

  call timing(iproc,'CrtLocPot     ','OF')

END SUBROUTINE CounterIonPotential
