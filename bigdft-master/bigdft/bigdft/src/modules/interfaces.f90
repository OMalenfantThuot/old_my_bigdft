!> @file
!! Define the module module_interfaces containing all interfaces
!!
!! @author
!!    Copyright (C) 2007-2011 BigDFT group (LG,DC)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!>  Modules which contains all interfaces
module module_interfaces

   implicit none

      interface
        subroutine copy_old_wavefunctions(nproc,orbs,psi,&
            &   wfd_old,psi_old)
        use module_defs, only: wp
        use module_types
        use compression
         implicit none
         integer, intent(in) :: nproc
         type(orbitals_data), intent(in) :: orbs
         type(wavefunctions_descriptors), intent(in) :: wfd_old
         real(wp), dimension(:), pointer :: psi,psi_old
        END SUBROUTINE copy_old_wavefunctions
      end interface


      interface
        subroutine orbitals_descriptors(iproc,nproc,norb,norbu,norbd,nspin,nspinor, &
                 nkpt,kpt,wkpt,orbs,linear_partition,basedist,basedistu,basedistd)
         use module_defs, only: gp
         use module_types
         implicit none
         integer, intent(in) :: linear_partition !< repartition mode for the linear scaling version
         integer, intent(in) :: iproc,nproc,norb,norbu,norbd,nkpt,nspin
         integer, intent(in) :: nspinor
         type(orbitals_data), intent(inout) :: orbs
         real(gp), dimension(nkpt), intent(in) :: wkpt
         real(gp), dimension(3,nkpt), intent(in) :: kpt
         integer, dimension(0:nproc-1,nkpt), intent(in), optional :: basedist
         integer, dimension(0:nproc-1,nkpt), intent(in), optional :: basedistu
         integer, dimension(0:nproc-1,nkpt), intent(in), optional :: basedistd
        END SUBROUTINE orbitals_descriptors
      end interface


      interface
        subroutine createWavefunctionsDescriptors(iproc,hx,hy,hz,atoms,rxyz,&
           &   crmult,frmult,calculate_bounds,Glr,output_denspot)
       use module_defs, only: gp
       use module_types
       use locregs
        implicit none
        !Arguments
        type(atoms_data), intent(in) :: atoms
        integer, intent(in) :: iproc
        real(gp), intent(in) :: hx,hy,hz,crmult,frmult
        real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
        !real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
        logical,intent(in) :: calculate_bounds
        type(locreg_descriptors), intent(inout) :: Glr
        logical, intent(in), optional :: output_denspot
        END SUBROUTINE createWavefunctionsDescriptors
      end interface

!!$      interface
!!$        subroutine createProjectorsArrays(lr,rxyz,at,ob,&
!!$          cpmult,fpmult,hx,hy,hz,dry_run,nlpsp,&
!!$          init_projectors_completely)
!!$       !n(c) use module_base
!!$       use module_types
!!$       use orbitalbasis
!!$       implicit none
!!$       type(atoms_data), intent(in) :: at
!!$       !type(orbitals_data), intent(in) :: orbs
!!$       type(orbital_basis), intent(in) :: ob
!!$       real(kind=8), intent(in) :: cpmult,fpmult,hx,hy,hz
!!$       type(locreg_descriptors),intent(in) :: lr
!!$       real(kind=8), dimension(3,at%astruct%nat), intent(in) :: rxyz
!!$       !real(kind=8), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
!!$       logical, intent(in) :: dry_run
!!$       type(DFT_PSP_projectors), intent(out) :: nlpsp
!!$       logical,intent(in),optional :: init_projectors_completely
!!$        END SUBROUTINE createProjectorsArrays
!!$      end interface


       interface
         subroutine IonicEnergyandForces(iproc,nproc,dpbox,at,elecfield,&
          & rxyz,eion,fion,dispersion,edisp,fdisp,ewaldstr,&
          & pot_ion,pkernel,psoffset)
       use module_defs, only: gp,dp
       use module_dpbox
       use module_types
       implicit none
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
         END SUBROUTINE IonicEnergyandForces
       end interface

!!$      interface
!     subroutine mp_calculate(rx,ry,rz,hxh,hyh,hzh,cutoff,rlocinv2sq,mp,mpx,mpy,mpz)
!       use module_base
!       use gaussians, only: mp_exp
!       !Arguments
!       real(gp), intent(in) :: rx,ry,rz,hxh,hyh,hzh,cutoff,rlocinv2sq
!       logical, intent(in) :: mp
!       real(gp), dimension(:), allocatable, intent(out) :: mpx,mpy,mpz
!     end subroutine mp_calculate
!!$     subroutine density_and_hpot(dpbox,symObj,orbs,Lzd,pkernel,rhodsc,GPU,xc,psi,rho,vh,rho_ion,hstrten)
!!$      use module_defs, only: gp,wp,dp
!!$      use module_dpbox
!!$      use module_types
!!$      use module_atoms, only: symmetry_data
!!$      use module_xc
!!$      implicit none
!!$      type(denspot_distribution), intent(in) :: dpbox
!!$      type(rho_descriptors),intent(inout) :: rhodsc
!!$      type(orbitals_data), intent(in) :: orbs
!!$      type(local_zone_descriptors), intent(in) :: Lzd
!!$      type(symmetry_data), intent(in) :: symObj
!!$        type(coulomb_operator), intent(inout) :: pkernel
!!$      type(xc_info), intent(in) :: xc
!!$      real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
!!$      type(GPU_pointers), intent(inout) :: GPU
!!$      real(gp), dimension(6), intent(out) :: hstrten
!!$      real(dp), dimension(:), pointer :: rho,vh
!!$        real(dp), dimension(:,:,:,:), pointer :: rho_ion
!!$        END SUBROUTINE density_and_hpot
!!$      end interface

      interface
        subroutine sumrho(dpbox,orbs,Lzd,GPU,symObj,rhodsc,xc,psi,rho_p,mapping)
      use module_defs, only: wp,gp,dp
      use module_atoms, only: symmetry_data
      use module_dpbox
      use module_types
      use module_xc
      implicit none
      !Arguments
      type(denspot_distribution), intent(in) :: dpbox
      type(rho_descriptors),intent(in) :: rhodsc
      type(orbitals_data), intent(in) :: orbs
      type(local_zone_descriptors), intent(in) :: Lzd
      type(symmetry_data), intent(in) :: symObj
      type(xc_info), intent(in) :: xc
      real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
      real(dp), dimension(:,:), pointer :: rho_p
      type(GPU_pointers), intent(inout) :: GPU
      integer,dimension(orbs%norb),intent(in),optional:: mapping
        END SUBROUTINE sumrho
      end interface

     !starting point for the communication routine of the density
      interface
        subroutine communicate_density(dpbox,nspin,rhodsc,rho_p,rho,keep_rhop)
        use module_defs, only: gp,dp,wp
      use module_dpbox
      use module_types
      implicit none
      logical, intent(in) :: keep_rhop !< preserves the total density in the rho_p array
      integer, intent(in) :: nspin
      type(rho_descriptors),intent(in) :: rhodsc
      type(denspot_distribution), intent(in) :: dpbox
      real(dp), dimension(:,:), pointer :: rho_p !< partial density in orbital distribution scheme
      real(dp), dimension(max(dpbox%mesh%ndims(1)*dpbox%mesh%ndims(2)*dpbox%n3d,1),nspin), intent(out) :: rho
        END SUBROUTINE communicate_density
      end interface


       interface
         subroutine LocalHamiltonianApplication(iproc,nproc,at,npsidim_orbs,orbs,&
          Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
          energs,SIC,GPU,PotOrKin,xc,pkernel,orbsocc,psirocc,dpbox,potential,comgp,hpsi_noconf,econf)
         use module_defs, only: gp,dp,wp
       use module_dpbox
       use module_types
       use module_xc
         use locreg_operations, only: confpot_data
       implicit none
       integer, intent(in) :: PotOrKin !< if true, only the potential operator is applied
       integer, intent(in) :: iproc,nproc,npsidim_orbs
       type(atoms_data), intent(in) :: at
       type(orbitals_data), intent(in) :: orbs
       type(local_zone_descriptors), intent(in) :: Lzd
       type(SIC_data), intent(in) :: SIC
       type(xc_info), intent(in) :: xc
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
       real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
       type(confpot_data), dimension(orbs%norbp) :: confdatarr
       real(wp), dimension(:), pointer :: pot
       !real(wp), dimension(*) :: pot
       type(energy_terms), intent(inout) :: energs
       real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(inout) :: hpsi
       type(GPU_pointers), intent(inout) :: GPU
       type(coulomb_operator), intent(inout), optional :: pkernel
       type(orbitals_data), intent(in), optional :: orbsocc
       real(wp), dimension(:), pointer, optional :: psirocc
       type(denspot_distribution),intent(in),optional :: dpbox
       real(wp), dimension(*), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
       type(p2pComms),intent(inout), optional:: comgp
       real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(inout),optional :: hpsi_noconf
       real(gp),intent(out),optional :: econf
         END SUBROUTINE LocalHamiltonianApplication
       end interface


      interface
        subroutine SynchronizeHamiltonianApplication(nproc,npsidim_orbs,orbs,Lzd,GPU,xc,hpsi,&
           energs,energs_work)
        use module_defs, only: gp,wp
        use module_types
        use module_xc
        implicit none
        integer, intent(in) :: nproc,npsidim_orbs
        type(orbitals_data),  intent(in) :: orbs
        type(local_zone_descriptors), intent(in) :: Lzd
        type(GPU_pointers), intent(inout) :: GPU
        type(xc_info), intent(in) :: xc
        type(energy_terms), intent(inout) :: energs
        !real(gp), intent(inout) :: ekin_sum,epot_sum,eproj_sum,eSIC_DC,eexctX
        real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: hpsi
        type(work_mpiaccumulate),optional,intent(inout) :: energs_work
        END SUBROUTINE SynchronizeHamiltonianApplication
      end interface

      interface
        subroutine hpsitopsi(iproc,nproc,iter,idsx,wfn,&
           at,nlpsp,eproj_sum)
        use module_defs, only: gp
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc,idsx,iter
         type(DFT_wavefunction), intent(inout) :: wfn
         type(atoms_data), intent(in) :: at
         type(DFT_PSP_projectors), intent(inout) :: nlpsp
         real(gp),optional, intent(out) :: eproj_sum
        END SUBROUTINE hpsitopsi
      end interface

      interface
        subroutine last_orthon(iproc,nproc,iter,wfn,evsum,opt_keeppsit)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        integer, intent(in) :: iproc,nproc,iter
        real(wp), intent(out) :: evsum
        type(DFT_wavefunction), intent(inout) :: wfn
        logical, optional :: opt_keeppsit
        END SUBROUTINE last_orthon
      end interface

      interface
        subroutine preconditionall2(iproc,nproc,orbs,Lzd,hx,hy,hz,ncong,npsidim,hpsi,confdatarr,gnrm,gnrm_zero, &
                 linear_precond_convol_workarrays, linear_precond_workarrays)
        use module_defs, only: gp,dp,wp
        use module_types
        use locreg_operations, only: workarrays_quartic_convolutions,workarr_precond,confpot_data
        implicit none
        integer, intent(in) :: iproc,nproc,ncong,npsidim
        real(gp), intent(in) :: hx,hy,hz
        type(local_zone_descriptors), intent(in) :: Lzd
        type(orbitals_data), intent(in) :: orbs
        real(dp), intent(out) :: gnrm,gnrm_zero
        real(wp), dimension(npsidim), intent(inout) :: hpsi
        type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
        type(workarrays_quartic_convolutions),dimension(orbs%norbp),intent(inout),optional :: linear_precond_convol_workarrays !< convolution workarrays for the linear case
        type(workarr_precond),dimension(orbs%norbp),intent(inout),optional :: linear_precond_workarrays !< workarrays for the linear case
        END SUBROUTINE preconditionall2
      end interface

      interface
        subroutine partial_density_free(rsflag,nproc,n1i,n2i,n3i,npsir,nspinn,nrhotot,&
            &   hfac,nscatterarr,spinsgn,psir,rho_p,ibyyzz_r) !ex-optional argument
         use module_defs, only: gp,dp,wp
         implicit none
         logical, intent(in) :: rsflag
         integer, intent(in) :: nproc,n1i,n2i,n3i,nrhotot,nspinn,npsir
         real(gp), intent(in) :: hfac,spinsgn
         integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
         real(wp), dimension(n1i,n2i,n3i,nspinn), intent(in) :: psir
         real(dp), dimension(n1i,n2i,nrhotot,nspinn), intent(inout) :: rho_p
         integer, dimension(:,:,:), pointer :: ibyyzz_r
        END SUBROUTINE partial_density_free
      end interface

      interface
        subroutine parse_cp2k_files(iproc,basisfile,orbitalfile,nat,ntypes,orbs,iatype,rxyz,&
            &   CP2K,wfn_cp2k)
         use module_defs, only: gp, wp
         use module_types, only: orbitals_data, gaussian_basis
         implicit none
         character(len=*), intent(in) :: basisfile,orbitalfile
         integer, intent(in) :: iproc,nat,ntypes
         type(orbitals_data), intent(in) :: orbs
         integer, dimension(nat), intent(in) :: iatype
         real(gp), dimension(3,nat), target, intent(in) :: rxyz
         type(gaussian_basis), intent(out) :: CP2K
         real(wp), dimension(:,:), pointer :: wfn_cp2k
        END SUBROUTINE parse_cp2k_files
      end interface

      interface
        subroutine read_gaussian_information(orbs,G,coeffs,filename, opt_fillrxyz)
          use module_defs, only: wp
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         type(orbitals_data), intent(inout) :: orbs
         type(gaussian_basis), intent(out) :: G
         real(wp), dimension(:,:), pointer :: coeffs
         logical, optional :: opt_fillrxyz
        END SUBROUTINE read_gaussian_information
      end interface

      interface
        subroutine restart_from_gaussians(iproc,nproc,orbs,Lzd,hx,hy,hz,psi,G,coeffs)
        use module_defs, only: gp,wp
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc
         real(gp), intent(in) :: hx,hy,hz
         type(orbitals_data), intent(in) :: orbs
         type(local_zone_descriptors), intent(in) :: Lzd
         type(gaussian_basis), intent(inout) :: G
         real(wp), dimension(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%norbp), intent(out) :: psi
         real(wp), dimension(:,:), pointer :: coeffs
        END SUBROUTINE restart_from_gaussians
      end interface

      interface
        subroutine inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,nvirt,nspin,&
            orbs,orbse,norbsc_arr,locrad,G,psigau,eks,iversion,mapping,quartic_prefactor)
        use module_defs, only: gp,wp
         use module_types
         implicit none
         integer, intent(in) :: iproc,nproc,nspin
         integer, intent(inout) :: nvirt
         type(atoms_data), intent(in) :: at
         type(orbitals_data), intent(in) :: orbs
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         real(gp), intent(out) :: eks
         integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
         real(gp), dimension(at%astruct%nat), intent(out) :: locrad
         type(orbitals_data), intent(out) :: orbse
         type(gaussian_basis), intent(out) :: G
         real(wp), dimension(:,:,:), pointer :: psigau
         integer,intent(in) :: iversion !< 1:cubic, 2:linear
         integer,dimension(orbs%norb),intent(in),optional:: mapping
         real(gp),dimension(at%astruct%ntypes),intent(in),optional:: quartic_prefactor
        END SUBROUTINE inputguess_gaussian_orbitals
      end interface


     interface
       subroutine inputguess_gaussian_orbitals_forLinear(iproc,nproc,norb,at,rxyz,nvirt,nspin,&
          nlr, norbsPerAt, mapping, &
          orbs,orbse,norbsc_arr,locrad,G,psigau,eks,quartic_prefactor)
       use module_defs, only: gp,dp,wp
       use module_types
       implicit none
       integer, intent(in) :: iproc,nproc,nspin,nlr,norb
       integer, intent(inout) :: nvirt
       type(atoms_data), intent(in) :: at
       type(orbitals_data), intent(in) :: orbs
       real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
       integer,dimension(norb),intent(in):: mapping
       integer,dimension(at%astruct%nat),intent(in):: norbsPerAt
       real(gp), intent(out) :: eks
       integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
       real(gp), dimension(at%astruct%nat), intent(out) :: locrad
       type(orbitals_data), intent(inout) :: orbse
       type(gaussian_basis), intent(out) :: G
       real(wp), dimension(:,:,:), pointer :: psigau
       real(gp),dimension(at%astruct%ntypes),intent(in),optional:: quartic_prefactor
       END SUBROUTINE inputguess_gaussian_orbitals_forLinear
     end interface

     interface
       subroutine AtomicOrbitals(iproc,at,rxyz,norbe,orbse,norbsc,&
          nspin,eks,G,gaucoeff,iorbtolr,mapping,quartic_prefactor)
       use module_defs, only: gp,dp,wp
       use module_types
       implicit none
       integer, intent(in) :: norbe,iproc
       integer, intent(in) :: norbsc,nspin
       type(atoms_data), intent(in) :: at
       !logical, dimension(4,2,at%natsc), intent(in) :: scorb
       real(gp), dimension(3,at%astruct%nat), intent(in), target :: rxyz
       type(orbitals_data), intent(inout) :: orbse
       type(gaussian_basis), intent(out) :: G
       real(gp), intent(out) :: eks
       integer, dimension(orbse%norbp), intent(out) :: iorbtolr !assign the localisation region
       real(wp), intent(out) :: gaucoeff !norbe=G%ncoeff !fake interface for passing address
       integer,dimension(orbse%norb), optional, intent(in):: mapping
       real(gp),dimension(at%astruct%ntypes),intent(in),optional:: quartic_prefactor
        END SUBROUTINE AtomicOrbitals
      end interface

      interface
        subroutine apply_potential(n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot,psir,pot,epot,&
            &   ibyyzz_r) !optional
         use module_defs, only: gp,dp,wp
         implicit none
         integer, intent(in) :: n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot
         real(wp), dimension(-nl1:2*n1+2+nl1,-nl2:2*n2+2+nl2,-nl3:2*n3+2+nl3,nspinor), intent(inout) :: psir
         real(wp), dimension(-nl1:2*n1+2+nl1-4*nbuf,-nl2:2*n2+2+nl2-4*nbuf,-nl3:2*n3+2+nl3-4*nbuf,npot), intent(in) :: pot
         integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
         real(gp), intent(out) :: epot
        END SUBROUTINE apply_potential
      end interface


      interface
        subroutine read_pw_waves(filename, iproc, nproc, at, rxyz, Glr, orbs, psig, rhoij)
        use module_defs, only: gp, wp
        use module_atoms
        use locregs
        use module_types
        implicit none
        character(len = *), intent(in) :: filename
        integer, intent(in) :: iproc, nproc
        type(atoms_data), intent(in) :: at
        real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
        type(locreg_descriptors), intent(in) :: Glr
        type(orbitals_data), intent(in) :: orbs
        real(wp), dimension(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, orbs%norbp), intent(out) :: psig
        real(wp), dimension(:,:,:), pointer, optional :: rhoij
        END SUBROUTINE read_pw_waves
      end interface
      interface
        subroutine read_potfile4b2B(filename,n1i,n2i,n3i, rho, alat1, alat2, alat3)
         use module_defs, only: gp,dp,wp
         implicit none
         character(len=*), intent(in) :: filename
         integer, intent(out) :: n1i,n2i,n3i
         real(gp) alat1, alat2, alat3, dum, dum1
         ! real(dp), dimension(n1i*n2i*n3d), intent(out) :: rho
         real(gp), pointer :: rho(:)
        END SUBROUTINE read_potfile4b2B
      end interface

      interface
        subroutine gaussian_pswf_basis(ng,enlargerprb,iproc,nspin,at,rxyz,G,Gocc, gaenes, &
            &   iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg )
        use module_defs, only: gp,wp
         use module_types
         implicit none
         logical, intent(in) :: enlargerprb
         integer, intent(in) :: iproc,nspin,ng
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%astruct%nat), target, intent(in) :: rxyz
         type(gaussian_basis), intent(out) :: G
         real(wp), dimension(:), pointer :: Gocc
         real(gp), pointer, optional :: gaenes(:)
         integer, pointer, optional :: iorbtolr(:)
         integer, pointer, optional :: iorbto_l(:)
         integer, pointer, optional :: iorbto_m(:)
         integer, pointer, optional :: iorbto_ishell(:)
         integer, pointer, optional :: iorbto_iexpobeg(:)
        END SUBROUTINE gaussian_pswf_basis
      end interface

      interface
        subroutine gaussian_pswf_basis_for_paw(at,rxyz,G,  &
            &   iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg, iorbto_paw_nchannels,&
         iorbto_imatrixbeg )
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%astruct%nat), target, intent(in) :: rxyz
         type(gaussian_basis_c), intent(inout) :: G

         integer, pointer :: iorbtolr(:)
         integer, pointer :: iorbto_l(:)
         integer, pointer :: iorbto_paw_nchannels(:)
         integer, pointer :: iorbto_m(:)
         integer, pointer :: iorbto_ishell(:)
         integer, pointer :: iorbto_iexpobeg(:)
         integer, pointer :: iorbto_imatrixbeg(:)

         !local variables
        END SUBROUTINE gaussian_pswf_basis_for_paw
      end interface


      interface
        subroutine local_analysis(iproc,nproc,hx,hy,hz,at,rxyz,lr,orbs,orbsv,psi,psivirt)
        use module_defs, only: gp,wp
        use module_types
        use locregs
         implicit none
         integer, intent(in) :: iproc,nproc
         real(gp), intent(in) :: hx,hy,hz
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(in) :: orbs,orbsv
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         real(wp), dimension(:), pointer :: psi,psivirt
        END SUBROUTINE local_analysis
      end interface

      interface
        subroutine plot_gatom_basis(filename,iat,ngx,G,Gocc,rhocoeff,rhoexpo)
        use module_defs, only: wp
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         integer, intent(in) :: iat,ngx
         type(gaussian_basis), intent(in) :: G
         real(wp), dimension(:), pointer :: Gocc
         real(wp), dimension((ngx*(ngx+1))/2), intent(out) :: rhoexpo
         real(wp), dimension((ngx*(ngx+1))/2,4), intent(out) :: rhocoeff
        END SUBROUTINE plot_gatom_basis
      end interface

      interface
        subroutine calculate_rhocore(at,rxyz,dpbox,rhocore)
        use module_defs, only: gp,dp,wp
        use module_atoms, only: atoms_data
        use module_dpbox, only: denspot_distribution
        implicit none
        type(atoms_data), intent(in) :: at
        type(denspot_distribution), intent(in) :: dpbox
        real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
        real(wp), dimension(:,:,:,:), pointer :: rhocore
        END SUBROUTINE calculate_rhocore
      end interface

      interface
        subroutine XC_potential(geocode,datacode,iproc,nproc,mpi_comm,n01,n02,n03,xc,hgrids,&
           rho,exc,vxc,nspin,rhocore,rhohat,potxc,xcstr,dvxcdrho)
        use module_defs, only: gp,dp,wp
        use module_xc
        implicit none
        character(len=1), intent(in) :: geocode  !< @copydoc poisson_solver::doc::geocode
        character(len=1), intent(in) :: datacode !< @copydoc poisson_solver::doc::datacode
        integer, intent(in) :: iproc,nproc,n01,n02,n03,nspin,mpi_comm
        type(xc_info), intent(in) :: xc
        real(gp), dimension(3), intent(in) :: hgrids
        real(gp), intent(out) :: exc,vxc
        real(dp), dimension(*), intent(inout) :: rho
        real(wp), dimension(:,:,:,:), pointer :: rhocore !associated if useful
        real(wp), dimension(:,:,:,:), pointer :: rhohat
        real(wp), dimension(*), intent(out) :: potxc
        real(dp), dimension(6), intent(out) :: xcstr
        real(dp), dimension(:,:,:,:), target, intent(out), optional :: dvxcdrho
        END SUBROUTINE XC_potential
      end interface

      interface
        subroutine xc_energy(geocode,m1,m3,md1,md2,md3,nxc,nwb,nxt,nwbl,nwbr,&
           nxcl,nxcr,xc,hx,hy,hz,rhopot,pot_ion,sumpion,zf,zfionxc,exc,vxc,nproc,nspden)
        use module_defs, only: gp,dp,wp
        use module_xc
        use abi_interfaces_xc_lowlevel, only: abi_drivexc,abi_size_dvxc
        implicit none
        character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
        logical, intent(in) :: sumpion
        integer, intent(in) :: m1,m3,nxc,nwb,nxcl,nxcr,nxt,md1,md2,md3,nproc,nspden
        integer, intent(in) :: nwbl,nwbr
        real(gp), intent(in) :: hx,hy,hz
        type(xc_info), intent(in) :: xc
        real(dp), dimension(m1,m3,nxt,nspden), intent(inout) :: rhopot
        real(wp), dimension(*), intent(in) :: pot_ion
        real(dp), dimension(md1,md3,md2/nproc), intent(out) :: zf
        real(wp), dimension(md1,md3,md2/nproc,nspden), intent(out) :: zfionxc
        real(dp), intent(out) :: exc,vxc
        END SUBROUTINE xc_energy
      end interface

      interface
        subroutine write_eigenvalues_data(etol,orbs,mom_vec)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        real(gp), intent(in) :: etol
        type(orbitals_data), intent(in) :: orbs
        real(gp), dimension(:,:,:), intent(in), pointer :: mom_vec
        END SUBROUTINE write_eigenvalues_data
      end interface

      interface
        subroutine write_eigen_objects(iproc,occorbs,nspin,nvirt,nplot,hx,hy,hz,at,rxyz,lr,orbs,orbsv,psi,psivirt)
        use module_defs, only: gp,wp
        use module_types
        use locregs
         implicit none
         logical, intent(in) :: occorbs
         integer, intent(in) :: iproc,nspin,nvirt,nplot
         real(gp), intent(in) :: hx,hy,hz
         type(atoms_data), intent(in) :: at
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(in) :: orbs,orbsv
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         real(wp), dimension(:), pointer :: psi,psivirt
         END SUBROUTINE write_eigen_objects
       end interface

      interface
        subroutine free_full_potential(nproc,flag,xc,pot)
         use module_defs, only: gp,dp,wp
         use module_xc
         implicit none
         integer, intent(in) :: nproc,flag
         type(xc_info), intent(in) :: xc
         real(wp), dimension(:), pointer :: pot
        END SUBROUTINE free_full_potential
      end interface

      interface
        subroutine orthoconstraint(iproc,nproc,orbs,comms,symm,tr_min,&
            psi,hpsi,scprsum,spsi) !n(c) wfd (arg:5)
        use module_defs, only: gp,dp,wp
        use module_types
        use communications_base, only: comms_cubic
        implicit none
        logical, intent(in) :: symm !< symmetrize the lagrange multiplier after calculation
        logical, intent(in) :: tr_min !< optimize the trace of the hamiltonian instead of the energy
        integer, intent(in) :: iproc,nproc
        type(orbitals_data), intent(in) :: orbs
        type(comms_cubic), intent(in) :: comms
        !n(c) type(wavefunctions_descriptors), intent(in) :: wfd
        real(wp), dimension(orbs%npsidim_comp), intent(in) :: psi
        real(wp), dimension(orbs%npsidim_comp), intent(inout) :: hpsi
        real(dp), intent(out) :: scprsum
        real(wp), dimension(orbs%npsidim_comp), optional, intent(in) :: spsi
        END SUBROUTINE orthoconstraint
      end interface


      interface
        subroutine constrained_davidson(iproc,nproc,in,at,&
           orbs,orbsv,nvirt,Lzd,comms,commsv,&
           hx,hy,hz,rxyz,rhopot,psi,v,dpbox,xc,GPU)
        use module_defs, only: gp,dp,wp
        use module_dpbox
        use module_types
        use communications_base, only: comms_cubic
        use module_xc
        implicit none
        integer, intent(in) :: iproc,nproc
        integer, intent(in) :: nvirt
        type(input_variables), intent(in) :: in
        type(atoms_data), intent(in) :: at
        type(local_zone_descriptors), intent(in) :: Lzd
        type(orbitals_data), intent(in) :: orbs
        type(comms_cubic), intent(in) :: comms, commsv
        type(denspot_distribution), intent(in) :: dpbox
        type(xc_info), intent(in) :: xc
        real(gp), intent(in) :: hx,hy,hz
        real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
        real(dp), dimension(*), intent(in) :: rhopot
        type(orbitals_data), intent(inout) :: orbsv
        type(GPU_pointers), intent(inout) :: GPU
        real(wp), dimension(:), pointer :: psi,v!=psivirt(nvctrp,nvirtep*nproc)
        !v, that is psivirt, is transposed on input and direct on output
        END SUBROUTINE constrained_davidson
      end interface

      interface
        subroutine NK_SIC_potential(lr,orbs,xc,fref,hgrids,pkernel,psi,poti,eSIC_DC,potandrho,wxdsave)
        use module_defs, only: gp,wp,dp
        use module_types
        use module_xc
        use locregs
         implicit none
         real(gp), intent(in) :: fref
         type(locreg_descriptors), intent(in) :: lr
         type(orbitals_data), intent(in) :: orbs
         type(coulomb_operator), intent(inout) :: pkernel
         type(xc_info), intent(in) :: xc
         real(gp), dimension(3), intent(in) :: hgrids
         real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
         !real(wp), dimension((lr%d%n1i*lr%d%n2i*lr%d%n3i*((orbs%nspinor/3)*3+1)),max(orbs%norbp,orbs%nspin)), intent(inout) :: poti
         real(wp), intent(inout) :: poti
         real(gp), intent(out) :: eSIC_DC
         real(dp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,2*orbs%nspin), intent(in), optional :: potandrho
         real(dp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspin), intent(out), optional :: wxdsave
        END SUBROUTINE NK_SIC_potential
      end interface

      interface
        subroutine readmywaves(iproc,filename,iformat,orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  &
         wfd,psi,orblist,pawrhoij)
         use module_defs, only: gp,dp,wp
         use module_types
         use compression
         use m_pawrhoij
         implicit none
         integer, intent(in) :: iproc,n1,n2,n3, iformat
         real(gp), intent(in) :: hx,hy,hz
         type(wavefunctions_descriptors), intent(in) :: wfd
         type(orbitals_data), intent(inout) :: orbs
         type(atoms_data), intent(in) :: at
         real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
         integer, dimension(orbs%norb), optional :: orblist
         real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
         real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(out) :: psi
         character(len=*), intent(in) :: filename
         type(pawrhoij_type), dimension(at%astruct%nat), intent(inout), optional :: pawrhoij
        END SUBROUTINE readmywaves
      end interface

      interface
        subroutine open_filename_of_iorb(unitfile,lbin,filename,orbs,iorb,ispinor,iorb_out,iorb_shift,iiorb)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         logical, intent(in) :: lbin
         integer, intent(in) :: iorb,ispinor
         integer, intent(inout) :: unitfile
         type(orbitals_data), intent(in) :: orbs
         integer, intent(out) :: iorb_out
         integer,intent(in),optional :: iorb_shift,iiorb
        END SUBROUTINE open_filename_of_iorb
      end interface

      interface
        subroutine filename_of_iorb(lbin,filename,orbs,iorb,ispinor,filename_out,iorb_out,iorb_shift,iiorb)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         character(len=*), intent(in) :: filename
         logical, intent(in) :: lbin
         integer, intent(in) :: iorb,ispinor
         type(orbitals_data), intent(in) :: orbs
         character(len=*), intent(out) :: filename_out
         integer, intent(out) :: iorb_out
         integer,intent(in),optional :: iorb_shift,iiorb
        END SUBROUTINE filename_of_iorb
      end interface

      interface
        subroutine verify_file_presence(filerad,orbs,iformat,nproc,nforb)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        integer, intent(in) :: nproc
        character(len=*), intent(in) :: filerad
        type(orbitals_data), intent(in) :: orbs
        integer, intent(out) :: iformat
        integer, optional, intent(in) :: nforb
        END SUBROUTINE verify_file_presence
      end interface

      interface
        subroutine readwavetoisf(lstat, filename, formatted, hx, hy, hz, &
           & n1, n2, n3, nspinor, psiscf)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        character(len = *), intent(in) :: filename
        logical, intent(in) :: formatted
        integer, intent(out) :: n1, n2, n3, nspinor
        real(gp), intent(out) :: hx, hy, hz
        real(wp), dimension(:,:,:,:), pointer :: psiscf
        logical, intent(out) :: lstat
        END SUBROUTINE readwavetoisf
      end interface
      interface
        subroutine readwavetoisf_etsf(lstat, filename, iorbp, hx, hy, hz, &
           & n1, n2, n3, nspinor, psiscf)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        character(len = *), intent(in) :: filename
        integer, intent(in) :: iorbp
        integer, intent(out) :: n1, n2, n3, nspinor
        real(gp), intent(out) :: hx, hy, hz
        real(wp), dimension(:,:,:,:), pointer :: psiscf
        logical, intent(out) :: lstat
        END SUBROUTINE readwavetoisf_etsf
      end interface

      interface
        subroutine read_wave_to_isf(lstat, filename, ln, iorbp, hx, hy, hz, &
           & n1, n2, n3, nspinor, psiscf)
        use module_defs, only: gp,dp,wp
        use module_types
        implicit none
        integer, intent(in) :: ln
        character(len = ln), intent(in) :: filename
        integer, intent(in) :: iorbp
        integer, intent(out) :: n1, n2, n3, nspinor
        real(gp), intent(out) :: hx, hy, hz
        real(wp), dimension(:,:,:,:), pointer :: psiscf
        logical, intent(out) :: lstat
        END SUBROUTINE read_wave_to_isf
      end interface
      interface
        subroutine free_wave_to_isf(psiscf)
        use module_defs, only: gp,dp,wp
        implicit none
        real(wp), dimension(:,:,:,:), pointer :: psiscf
        END SUBROUTINE free_wave_to_isf
      end interface

    interface
      subroutine inputguessConfinement(iproc, nproc, at, input, hx, hy, hz, &
         rxyz, nlpsp, GPU, orbs, kswfn, tmb, denspot, rhopotold, energs,&
         locregcenters)
      ! Input wavefunctions are found by a diagonalization in a minimal basis set
      ! Each processors write its initial wavefunctions into the wavefunction file
      ! The files are then read by readwave
      use module_defs, only: gp,dp,wp
      use module_types
      implicit none
      !Arguments
      integer, intent(in) :: iproc,nproc
      real(gp), intent(in) :: hx, hy, hz
      type(atoms_data), intent(inout) :: at
      type(DFT_PSP_projectors), intent(inout) :: nlpsp
      type(GPU_pointers), intent(inout) :: GPU
      type(input_variables),intent(in) :: input
      real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
      type(orbitals_data),intent(inout) :: orbs
      type(DFT_wavefunction),intent(inout) :: kswfn, tmb
      type(DFT_local_fields), intent(inout) :: denspot
      real(dp), dimension(max(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin), intent(inout) ::  rhopotold
      type(energy_terms),intent(inout) :: energs
      real(kind=8),dimension(3,at%astruct%nat),intent(in),optional :: locregcenters
      END SUBROUTINE inputguessConfinement
    end interface

     interface
       subroutine LDiagHam(iproc,nproc,natsc,nspin,orbs,Lzd,Lzde,comms,&
          psi,hpsi,psit,orthpar,passmat,mixing,Tel,occopt,& !mandatory
          orbse,commse,etol,norbsc_arr) !optional
       use module_defs, only: gp,dp,wp
       use module_types
       use communications_base, only: comms_cubic
       implicit none
       logical, intent(in) :: mixing
       integer, intent(in) :: iproc,nproc,natsc,nspin,occopt
       real(gp), intent(in) :: Tel
       type(local_zone_descriptors) :: Lzd        !< Information about the locregs after LIG
       type(local_zone_descriptors) :: Lzde       !< Information about the locregs for LIG
       type(comms_cubic), intent(in) :: comms
       type(orbitals_data), intent(inout) :: orbs
       type(orthon_data), intent(in):: orthpar
       real(wp), dimension(*), intent(out) :: passmat !< passage matrix for building the eigenvectors (the size depends of the optional arguments)
       real(wp), dimension(:), pointer :: psi,hpsi,psit
       real(gp), intent(in) :: etol
       type(orbitals_data), intent(inout) :: orbse
       type(comms_cubic), intent(in) :: commse
       integer, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
       END SUBROUTINE LDiagHam
     end interface

     interface
       subroutine FullHamiltonianApplication(iproc,nproc,at,orbs,&
          Lzd,nlpsp,confdatarr,ngatherarr,Lpot,psi,hpsi,paw,&
          energs,SIC,GPU,xc,pkernel,orbsocc,psirocc)
       use module_defs, only: gp,dp,wp
       use module_types
       use module_xc
       use gaussians, only: gaussian_basis
       use locreg_operations, only: confpot_data
       implicit none
       integer, intent(in) :: iproc,nproc!,nspin
       type(atoms_data), intent(in) :: at
       type(orbitals_data), intent(in) :: orbs
       type(local_zone_descriptors),intent(in) :: Lzd
       type(DFT_PSP_projectors), intent(inout) :: nlpsp
       type(SIC_data), intent(in) :: SIC
       type(xc_info), intent(in) :: xc
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
       real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
       type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
       !real(wp), dimension(lzd%ndimpotisf) :: Lpot
       real(wp), dimension(:),pointer :: Lpot
       type(energy_terms), intent(inout) :: energs
       real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(out) :: hpsi
       type(GPU_pointers), intent(inout) :: GPU
       type(coulomb_operator), intent(in), optional :: pkernel
       type(orbitals_data), intent(in), optional :: orbsocc
       real(wp), dimension(:), pointer, optional :: psirocc
       !PAW variables:
       type(paw_objects),intent(inout)::paw
       END SUBROUTINE FullHamiltonianApplication
     end interface

       interface
         subroutine export_grids(fname, atoms, rxyz, hx, hy, hz, n1, n2, n3, logrid_c, logrid_f)
         use module_defs, only: gp
         use module_types
         implicit none
         character(len = *), intent(in) :: fname
         type(atoms_data), intent(in) :: atoms
         real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
         real(gp), intent(in) :: hx, hy, hz
         integer, intent(in) :: n1, n2, n3
         logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid_c
         logical, dimension(0:n1,0:n2,0:n3), intent(in), optional :: logrid_f
         END SUBROUTINE export_grids
       end interface

       interface
         subroutine system_initialization(iproc,nproc,dump,inputpsi,input_wf_format,&
            & dry_run,in,atoms,rxyz,OCLconv,&
            orbs,lnpsidim_orbs,lnpsidim_comp,lorbs,Lzd,Lzd_lin,nlpsp,comms,&
            ref_frags, denspot, locregcenters, inwhichlocreg_old, onwhichatom_old, &
            norb_par_ref, norbu_par_ref, norbd_par_ref,output_grid)
         use module_defs, only: gp,dp,wp
         use f_enums, only: f_enumerator
         use module_types
         use module_fragments
         use communications_base, only: comms_cubic
         implicit none
         integer, intent(in) :: iproc,nproc
         integer, intent(out) :: input_wf_format,lnpsidim_orbs,lnpsidim_comp
         type(f_enumerator), intent(inout) :: inputpsi
         type(input_variables), intent(inout) :: in
         type(atoms_data), intent(inout) :: atoms
         real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
         logical, intent(in) :: OCLconv
         type(orbitals_data), intent(inout) :: orbs,lorbs
         type(local_zone_descriptors), intent(inout) :: Lzd, Lzd_lin
         type(DFT_local_fields), intent(out), optional :: denspot
         type(DFT_PSP_projectors), intent(out) :: nlpsp
         type(comms_cubic), intent(out) :: comms
         !real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
         type(system_fragment), dimension(:), pointer :: ref_frags
         real(kind=8),dimension(3,atoms%astruct%nat),intent(inout),optional :: locregcenters
         integer,dimension(:),pointer,optional:: inwhichlocreg_old, onwhichatom_old
         integer,dimension(0:nproc-1),optional:: norb_par_ref, norbu_par_ref, norbd_par_ref !< support function distribution to be used as a reference
         logical, intent(in) :: dry_run, dump
         logical, intent(in), optional :: output_grid
         END SUBROUTINE system_initialization
       end interface

       interface
         subroutine input_check_psi_id(inputpsi, input_wf_format, dir_output, orbs, lorbs, iproc, nproc, nfrag, frag_dir, ref_frags)
         use module_types
         use f_enums
         use module_fragments
         implicit none
         integer, intent(out) :: input_wf_format         !< (out) Format of WF
         type(f_enumerator), intent(inout) :: inputpsi   !< (in) indicate how check input psi, (out) give how to build psi
         integer, intent(in) :: iproc                    !< (in)  id proc
         integer, intent(in) :: nproc                    !< (in)  #proc
         integer, intent(in) :: nfrag                    !< number of fragment directories which need checking
         type(system_fragment), dimension(:), pointer :: ref_frags  !< number of orbitals for each fragment
         character(len=100), dimension(nfrag), intent(in) :: frag_dir !< label for fragment subdirectories (blank if not a fragment calculation)
         character(len = *), intent(in) :: dir_output
         type(orbitals_data), intent(in) :: orbs, lorbs
         END SUBROUTINE input_check_psi_id
       end interface

       interface
         subroutine calc_gradient(geocode,n1,n2,n3,n3grad,deltaleft,deltaright,rhoinp,nspden,hx,hy,hz,&
            gradient,rhocore)
         use module_defs, only: gp,dp,wp
         implicit none
         !Arguments
         character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
         integer, intent(in) :: n1,n2,n3,n3grad,deltaleft,deltaright,nspden
         real(dp), intent(in) :: hx,hy,hz
         real(dp), dimension(n1,n2,n3,nspden), intent(inout) :: rhoinp
         real(dp), dimension(n1,n2,n3grad,2*nspden-1,0:3), intent(out) :: gradient
         real(dp), dimension(:,:,:,:), pointer :: rhocore
         END SUBROUTINE calc_gradient
       end interface


       interface
         subroutine update_locreg(iproc, nproc, nlr, locrad, locrad_kernel, locrad_mult, locregCenter, glr_tmp, &
                  useDerivativeBasisFunctions, nscatterarr, hx, hy, hz, astruct, input, &
                  orbs_KS, orbs, lzd, npsidim_orbs, npsidim_comp, lbcomgp, lbcollcom, lfoe, lice, lbcollcom_sr)
         use module_defs, only: gp,dp,wp
         use module_types
         use foe_base, only: foe_data
         use communications_base, only: p2pComms
         use locregs
         implicit none
         integer,intent(in):: iproc, nproc, nlr
         integer,intent(out) :: npsidim_orbs, npsidim_comp
         logical,intent(in):: useDerivativeBasisFunctions
         integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
         real(8),intent(in):: hx, hy, hz
         type(atomic_structure),intent(in) :: astruct
         type(input_variables),intent(in) :: input
         real(8),dimension(nlr),intent(in):: locrad, locrad_kernel, locrad_mult
         type(orbitals_data),intent(in):: orbs_KS, orbs
         real(8),dimension(3,nlr),intent(in):: locregCenter
         type(locreg_descriptors),intent(in):: glr_tmp
         type(local_zone_descriptors),intent(inout):: lzd
         type(p2pComms),intent(inout):: lbcomgp
         type(foe_data),intent(inout),optional :: lfoe
         type(comms_linear),intent(inout):: lbcollcom
         type(comms_linear),intent(inout),optional :: lbcollcom_sr
         type(foe_data),intent(inout),optional :: lice
         END SUBROUTINE update_locreg
       end interface

       interface
         subroutine init_orbitals_data_for_linear(iproc, nproc, nspinor, input, astruct, rxyz, lorbs, &
           norb_par_ref, norbu_par_ref, norbd_par_ref)
         use module_defs, only: gp,dp,wp
         use module_types
         implicit none
         integer,intent(in):: iproc, nproc, nspinor
         type(input_variables),intent(in):: input
         type(atomic_structure),intent(in):: astruct
         real(8),dimension(3,astruct%nat),intent(in):: rxyz
         type(orbitals_data),intent(out):: lorbs
         integer,dimension(0:nproc-1),intent(in),optional :: norb_par_ref, norbu_par_ref, norbd_par_ref
         END SUBROUTINE init_orbitals_data_for_linear
       end interface

       interface
         subroutine psi_to_vlocpsi(iproc,npsidim_orbs,orbs,Lzd,&
            ipotmethod,confdatarr,pot,psi,vpsi,pkernel,xc,alphaSIC,epot_sum,evSIC,comgp,vpsi_noconf,econf_sum)
         use module_defs, only: gp,dp,wp
         use module_types
         use module_xc
         use locreg_operations, only: confpot_data
         implicit none
         integer, intent(in) :: iproc,ipotmethod,npsidim_orbs
         real(gp), intent(in) :: alphaSIC
         type(xc_info), intent(in) :: xc
         type(orbitals_data), intent(in) :: orbs
         type(local_zone_descriptors), intent(in) :: Lzd
         type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
         real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi !this dimension will be modified
         real(wp), dimension(*) :: pot !< the potential, with the dimension compatible with the ipotmethod flag
         real(gp), intent(out) :: epot_sum,evSIC
         real(wp), dimension(orbs%npsidim_orbs), intent(inout) :: vpsi
         type(coulomb_operator), intent(inout) ::  pkernel !< the PSolver kernel which should be associated for the SIC schemes
         type(p2pcomms),intent(in) :: comgp
         real(wp), dimension(orbs%npsidim_orbs), intent(inout),optional :: vpsi_noconf
         real(gp),intent(out),optional :: econf_sum
         END SUBROUTINE psi_to_vlocpsi
       end interface

       interface
         subroutine initialize_linear_from_file(iproc,nproc,input_frag,astruct,rxyz,orbs,Lzd,&
              iformat,dir_output,filename,ref_frags,orblist)
         use module_defs, only: gp,dp,wp
         use module_types
         use module_fragments
         implicit none
         integer, intent(in) :: iproc, nproc, iformat
         type(orbitals_data), intent(inout) :: orbs  !< orbs related to the basis functions, inwhichlocreg generated in this routine
         type(atomic_structure), intent(in) :: astruct
         real(gp), dimension(3,astruct%nat), intent(in) :: rxyz
         character(len=*), intent(in) :: filename, dir_output
         type(local_zone_descriptors), intent(inout) :: Lzd !< must already contain Glr and hgrids
         type(fragmentInputParameters), intent(in) :: input_frag
         type(system_fragment), dimension(input_frag%nfrag_ref), intent(inout) :: ref_frags
         integer, dimension(orbs%norb), optional :: orblist
         END SUBROUTINE initialize_linear_from_file
       end interface

        interface
          subroutine readmywaves_linear_new(iproc,nproc,dir_output,filename,iformat,at,tmb,rxyz,&
               ref_frags,input_frag,frag_calc,kernel_restart,max_nbasis_env,frag_env_mapping,orblist)
          use module_defs, only: gp,dp,wp
          use module_types
          use module_fragments
          implicit none
          integer, intent(in) :: iproc, nproc
          integer, intent(in) :: iformat
          type(atoms_data), intent(in) :: at
          type(DFT_wavefunction), intent(inout) :: tmb
          real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
          !real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
          character(len=*), intent(in) :: dir_output, filename
          type(fragmentInputParameters), intent(in) :: input_frag
          type(system_fragment), dimension(:), pointer :: ref_frags
          logical, intent(in) :: frag_calc, kernel_restart
          integer, intent(in) :: max_nbasis_env
          integer, dimension(input_frag%nfrag,max_nbasis_env,3), intent(inout) :: frag_env_mapping
          integer, dimension(tmb%orbs%norb), intent(in), optional :: orblist
          END SUBROUTINE readmywaves_linear_new
        end interface

        interface
          subroutine allocate_auxiliary_basis_function(npsidim, subname, lphi, lhphi)
          use module_defs, only: gp,dp,wp
          implicit none
          integer,intent(in):: npsidim
          real(8),dimension(:),pointer,intent(out):: lphi, lhphi
          character(len=*),intent(in):: subname
          END SUBROUTINE allocate_auxiliary_basis_function
        end interface

        interface
          subroutine deallocate_auxiliary_basis_function(subname, lphi, lhphi)
          use module_defs, only: gp,dp,wp
          implicit none
          real(8),dimension(:),pointer:: lphi, lhphi
          character(len=*),intent(in):: subname
          END SUBROUTINE deallocate_auxiliary_basis_function
        end interface

        interface
          subroutine denspot_set_history(denspot, iscf, &
               npulayit)
            use module_types, only: DFT_local_fields
          use f_enums, only: f_enumerator
          implicit none
          type(DFT_local_fields), intent(inout) :: denspot
          type(f_enumerator), intent(in) :: iscf
          integer,intent(in),optional :: npulayit
          END SUBROUTINE denspot_set_history
        end interface

        interface
          subroutine cholesky(iproc, nspin,norbIn, psi, &
          orbs, comms, ndim_ovrlp, ovrlp, norbTot, block1, &
          ispinIn, paw)
          !use module_defs, only: gp,dp,wp
          use module_types
          use communications_base, only: comms_cubic
          implicit none

          integer:: iproc,nvctrp,norbIn, nspin, block1, ispinIn
          type(orbitals_data), intent(in) :: orbs
          type(comms_cubic):: comms
          real(kind=8),dimension(orbs%npsidim_comp),intent(in out):: psi
          integer,dimension(nspin,0:orbs%nkpts):: ndim_ovrlp
          real(kind=8),dimension(ndim_ovrlp(nspin,orbs%nkpts),1):: ovrlp
          integer,dimension(nspin):: norbTot
          type(paw_objects),optional,intent(inout)::paw
          END SUBROUTINE cholesky
        end interface

        interface
          subroutine gsChol(iproc, nproc, psi, orthpar, nspinor,&
          orbs, nspin,ndim_ovrlp,norbArr,comms,paw)
          use module_defs, only: wp
          use module_types
          use communications_base, only: comms_cubic
          implicit none
          integer, intent(in) :: iproc, nproc,nspin
          integer, intent(inout) ::  nspinor
          type(orthon_data), intent(in):: orthpar
          type(orbitals_data):: orbs
          type(comms_cubic), intent(in) :: comms
          integer, dimension(nspin), intent(in) :: norbArr
          integer, dimension(nspin,0:orbs%nkpts), intent(inout) :: ndim_ovrlp
          real(wp),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(inout):: psi
          type(paw_objects),optional,intent(inout)::paw
          END SUBROUTINE gsCHol
        end interface

        interface
          subroutine loewdin(iproc, norbIn, block1, ispinIn,&
          orbs, comms, nspin, psit, ovrlp, ndim_ovrlp, norbTot, paw)
          !use module_defs, only: gp,dp,wp
          use module_types
          use communications_base, only: comms_cubic
          implicit none
          integer,intent(in):: iproc,norbIn, nspin, block1, ispinIn
          type(orbitals_data),intent(in):: orbs
          type(comms_cubic),intent(in):: comms
          real(kind=8),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(in out):: psit
          integer,dimension(nspin,0:orbs%nkpts):: ndim_ovrlp
          real(kind=8),dimension(ndim_ovrlp(nspin,orbs%nkpts)):: ovrlp
          integer,dimension(nspin):: norbTot
          type(paw_objects),optional,intent(inout)::paw
          END SUBROUTINE loewdin
        end interface

        interface
          subroutine gramschmidt(iproc, norbIn, psit, ndim_ovrlp, ovrlp, orbs, nspin,&
          nspinor, comms, norbTot, block1, block2, ispinIn,paw)
          use module_defs, only: wp !module_base
          use module_types
          use communications_base, only: comms_cubic
          implicit none
          integer,intent(in):: iproc, norbIn, nspin, block1, block2, ispinIn
          integer, intent(out) :: nspinor
          type(orbitals_data):: orbs
          type(comms_cubic), intent(in) :: comms
          type(paw_objects),optional,intent(inout)::paw
          real(wp),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(inout):: psit
          integer,dimension(nspin,0:orbs%nkpts):: ndim_ovrlp
          real(wp),dimension(ndim_ovrlp(nspin,orbs%nkpts)):: ovrlp
          integer,dimension(nspin):: norbTot
          END SUBROUTINE gramschmidt
        end interface

        interface
          subroutine orthogonalize(iproc,nproc,orbs,comms,psi,orthpar,paw)
          use module_defs, only: wp
          use module_types
          use communications_base, only: comms_cubic
          implicit none
          integer, intent(in) :: iproc,nproc
          type(orbitals_data), intent(in) :: orbs
          type(comms_cubic), intent(in) :: comms
          type(orthon_data), intent(in) :: orthpar
          real(wp), dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb), intent(inout) :: psi
          type(paw_objects),optional,intent(inout) :: paw
          END SUBROUTINE orthogonalize
        end interface

  interface
     subroutine copy_old_supportfunctions(iproc,orbs,lzd,phi,lzd_old,phi_old)
       use module_defs, only: gp,dp,wp
       use module_types
       implicit none
       integer,intent(in) :: iproc
       type(orbitals_data), intent(in) :: orbs
       type(local_zone_descriptors), intent(in) :: lzd
       type(local_zone_descriptors), intent(inout) :: lzd_old
       real(wp), dimension(:), pointer :: phi,phi_old
     END SUBROUTINE copy_old_supportfunctions
  end interface

  interface
     subroutine copy_old_coefficients(norb_tmb, nfvctr, coeff, coeff_old)
       use module_defs, only: gp,dp,wp
       implicit none
       integer,intent(in):: norb_tmb, nfvctr
       real(8),dimension(:,:),pointer:: coeff, coeff_old
     END SUBROUTINE copy_old_coefficients
  end interface

  interface
     subroutine copy_old_inwhichlocreg(norb_tmb, inwhichlocreg, inwhichlocreg_old, onwhichatom, onwhichatom_old)
       use module_defs, only: gp,dp,wp
       implicit none
       integer,intent(in):: norb_tmb
       integer,dimension(:),pointer:: inwhichlocreg, inwhichlocreg_old, onwhichatom, onwhichatom_old
     END SUBROUTINE copy_old_inwhichlocreg
  end interface

  interface
     subroutine reformat_supportfunctions(iproc,nproc,at,rxyz_old,rxyz,add_derivatives,tmb,ndim_old,lzd_old,&
          frag_trans,psi_old,input_dir,input_frag,ref_frags,max_shift,phi_array_old)
       use module_defs, only: gp,dp,wp
       use module_types
       use module_fragments
       use rototranslations, only: rototranslation
       implicit none
       integer, intent(in) :: iproc,nproc
       integer, intent(in) :: ndim_old
       type(atoms_data), intent(in) :: at
       real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz,rxyz_old
       type(DFT_wavefunction), intent(inout) :: tmb
       type(local_zone_descriptors), intent(in) :: lzd_old
       type(rototranslation), dimension(tmb%orbs%norbp), intent(in) :: frag_trans
       real(wp), dimension(:), pointer :: psi_old
       type(phi_array), dimension(tmb%orbs%norbp), optional, intent(in) :: phi_array_old
       logical, intent(in) :: add_derivatives
       character(len=*), intent(in) :: input_dir
       type(fragmentInputParameters), intent(in) :: input_frag
       type(system_fragment), dimension(:), pointer :: ref_frags
       real(gp),intent(out) :: max_shift
     END SUBROUTINE reformat_supportfunctions
  end interface

  interface
     subroutine integral_equation(iproc,nproc,atoms,wfn,ngatherarr,local_potential,GPU,xc,nlpsp,rxyz,paw)
       use module_defs, only: gp,dp,wp
       use module_types
       use module_xc
       implicit none
       integer, intent(in) :: iproc,nproc
       type(atoms_data), intent(in) :: atoms
       type(DFT_wavefunction), intent(in) :: wfn
       type(GPU_pointers), intent(inout) :: GPU
       type(DFT_PSP_projectors), intent(inout) :: nlpsp
       type(xc_info), intent(in) :: xc
       type(paw_objects), intent(inout) :: paw
       integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
       real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
       real(dp), dimension(:), pointer :: local_potential
     END SUBROUTINE integral_equation
  end interface

  interface
     subroutine atoms_new(atoms)
       use module_types
       implicit none
       type(atoms_data), pointer :: atoms
     END SUBROUTINE atoms_new
  end interface

  interface
     subroutine inputs_new(in)
       use module_input_keys, only: input_variables
       implicit none
       type(input_variables), pointer :: in
     END SUBROUTINE inputs_new
  end interface

  interface
     subroutine calculate_residue_ks(iproc, nproc, num_extra, ksorbs, tmb, hpsit_c, hpsit_f)
       use module_types
       implicit none

       ! Calling arguments
       integer, intent(in) :: iproc, nproc, num_extra
       type(dft_wavefunction), intent(inout) :: tmb
       type(orbitals_data), intent(in) :: ksorbs
       real(kind=8),dimension(:),pointer :: hpsit_c, hpsit_f
     END SUBROUTINE calculate_residue_ks
  end interface


  interface
     subroutine plot_wf(units_provided,orbname,nexpo,at,factor,lr,hgrids,rxyz,psi, &
          unit0_, unitx_, unity_, unitz_)
       use module_defs, only: gp,dp,wp
       use locregs, only: locreg_descriptors
       use module_types, only: atoms_data
       implicit none
       logical,intent(in) :: units_provided
       character(len=*) :: orbname
       integer, intent(in) :: nexpo
       real(dp), intent(in) :: factor
       real(gp), dimension(3), intent(in) :: hgrids
       type(atoms_data), intent(in) :: at
       real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
       type(locreg_descriptors), intent(in) :: lr
       real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
       integer,intent(in),optional :: unit0_, unitx_, unity_, unitz_
     END SUBROUTINE plot_wf
  end interface

  interface
     subroutine write_orbital_density(iproc, transform_to_global, iformat, &
          filename, npsidim, psi, orbs, lzd_g, at, rxyz, dens, iorb_shift, lzd_l, in_which_locreg)
       use module_defs, only: gp,dp,wp
       use module_types
       implicit none
       logical,intent(in) :: transform_to_global
       character(len=*),intent(in) :: filename
       integer,intent(in) :: iproc, npsidim, iformat
       real(kind=8),dimension(npsidim),intent(in),target :: psi
       type(orbitals_data),intent(in) :: orbs !< orbitals descriptors
       type(local_zone_descriptors),intent(inout) :: lzd_g !< global descriptors
       type(atoms_data),intent(in) :: at
       real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
       integer,intent(in),optional :: iorb_shift
       type(local_zone_descriptors),intent(in),optional :: lzd_l !< local descriptors
       logical,intent(in) :: dens !< density of wavefunctions or just wavefunctions
       integer,dimension(orbs%norb),intent(in),optional :: in_which_locreg
     END SUBROUTINE write_orbital_density
  end interface

END MODULE module_interfaces
