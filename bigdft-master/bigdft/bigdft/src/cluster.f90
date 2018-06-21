!> @file
!!   Routines to use BigDFT as a blackbox
!! @author
!!   Copyright (C) 2005-2015 BigDFT group
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS


!>  Main routine which does self-consistent loop.
!!  Does not parse input file and no geometry optimization.
!!  Does an electronic structure calculation.
!!  Output is the total energy and the forces
!!   @warning psi, keyg, keyv and eval should be freed after use outside of the routine.
subroutine cluster(nproc,iproc,atoms,rxyz,energy,energs,fxyz,strten,fnoise,pressure,&
     KSwfn,tmb,rxyz_old,in,GPU,infocode)
  use module_base
  use locregs, only: deallocate_locreg_descriptors
  use module_dpbox, only: dpbox_free
  use module_types
  use module_interfaces, only: IonicEnergyandForces, &
       & XC_potential, communicate_density, copy_old_wavefunctions, &
       denspot_set_history, &
       & gaussian_pswf_basis, local_analysis, &
       & orbitals_descriptors, sumrho, system_initialization,readmywaves,FullHamiltonianApplication
  use gaussians, only: deallocate_gwf
  use module_fragments
  use constrained_dft
  use f_ternary
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use module_xc
  use m_libpaw_libxc, only: libxc_functionals_init, libxc_functionals_end
  use communications_init, only: orbitals_communicators, write_memory_requirements_collcom
  use communications_base, only: deallocate_comms
  !  use vdwcorrection
  use yaml_output
  use psp_projectors_base, only: free_dft_psp_projectors
  use rhopotential, only: full_local_potential
  use sparsematrix_base, only: sparse_matrix_null, matrices_null, allocate_matrices, &
                               SPARSE_TASKGROUP, sparsematrix_malloc_ptr, assignment(=), &
                               DENSE_PARALLEL, DENSE_MATMUL, SPARSE_FULL, sparse_matrix_metadata_null
  use sparsematrix_init, only: write_sparsematrix_info
  use bigdft_matrices, only: check_local_matrix_extents, init_matrixindex_in_compressed_fortransposed, &
                             init_bigdft_matrices
  use sparsematrix_wrappers, only: init_sparse_matrix_wrapper, init_sparse_matrix_for_KSorbs, check_kernel_cutoff
  use sparsematrix, only: check_matrix_compression
  use communications_base, only: comms_linear_null
  use unitary_tests, only: check_communication_potential, check_communication_sumrho, &
                           check_communications_locreg
  use multipole, only: potential_from_charge_multipoles, interaction_multipoles_ions, &
                       ionic_energy_of_external_charges
  use public_enums
  use module_input_keys, only: SIC_data_null,print_dft_parameters,inputpsiid_set_policy,set_inputpsiid
  use orbitalbasis
  use io, only: plot_density,io_files_exists, writemywaves
  use PSbox, only: PS_gather
  use foe_common, only: init_foe
  use test_mpi_wrappers, only: test_mpi_alltoallv
  implicit none
  !Arguments
  integer, intent(in) :: nproc,iproc
  type(input_variables), intent(inout) :: in !SM: I had to change this to inout due to the new CheSS parameters within input_variables, maybe to be changed again later...
  type(atoms_data), intent(inout) :: atoms
  type(GPU_pointers), intent(inout) :: GPU
  type(DFT_wavefunction), intent(inout) :: KSwfn, tmb
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz_old
  real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
  type(energy_terms), intent(out) :: energs
  real(gp), intent(out) :: energy,fnoise,pressure
  real(gp), dimension(6), intent(out) :: strten
  real(gp), dimension(3,atoms%astruct%nat), intent(out) :: fxyz
  !> Encloses some information about the status of the run
  !!   - 0 run successfully succeded
  !!   - 1 the run ended after the allowed number of minimization steps. gnrm_cv not reached
  !!       forces may be meaningless
  !!   - 2 (present only for inputPsiId=INPUT_PSI_MEMORY_WVL) gnrm of the first iteration > 1 AND growing in
  !!       the second iteration OR grnm 1st >2.
  !!       Input wavefunctions need to be recalculated. Routine exits.
  !!   - 3 (present only for inputPsiId=INPUT_PSI_LCAO) gnrm > 4. SCF error. Routine exits.
  integer, intent(out) :: infocode
  !Local variables
  character(len=*), parameter :: subname='cluster'
  character(len=5) :: gridformat, wfformat
  logical :: refill_proj, calculate_dipole, calculate_quadrupole !,potential_from_disk=.false.
  logical :: DoDavidson,DoLastRunThings
  integer :: nvirt,norbv
  integer :: i, input_wf_format
  integer :: n1,n2,n3
  integer :: ncount0,ncount1,ncount_rate,ncount_max,n1i,n2i,n3i
  integer :: igroup,ikpt,ifrag
  !integer :: nproctiming
  real :: tcpu0,tcpu1
  real(kind=8) :: tel
  type(local_zone_descriptors) :: lzd_old
  type(DFT_PSP_projectors) :: nlpsp
  type(DFT_wavefunction) :: VTwfn !< Virtual wavefunction
  type(DFT_wavefunction) :: tmb_old
  !!type(DFT_wavefunction) :: tmb
  type(system_fragment), dimension(:), pointer :: ref_frags
  type(cdft_data) :: cdft
  real(dp), dimension(6) :: ewaldstr,xcstr
  real(gp), dimension(:,:), allocatable :: thetaphi,band_structure_eval,rxyz_tmp
  real(gp), dimension(:,:), pointer :: fdisp,fion,fpulay
  ! Charge density/potential,ionic potential, pkernel
  type(DFT_local_fields) :: denspot
  type(DFT_optimization_loop) :: optLoop
  real(gp), dimension(:), allocatable:: denspot0
  !wavefunction gradients, hamiltonian on vavefunction
  !transposed  wavefunction
  ! Pointers and variables to store the last psi
  ! before reformatting if useFormattedInput is .true.
  real(wp), dimension(:), pointer :: psi_old
  type(memory_estimation) :: mem
  !real(gp) :: energy_constrained
  ! PSP projectors
  real(kind=8), dimension(:), pointer :: gbd_occ!,rhocore
  ! Variables for the virtual orbitals and band diagram.
  integer :: nkptv, nvirtu, nvirtd
  real(gp), dimension(:), allocatable :: wkptv,psi_perturbed,hpsi_perturbed
  real(gp), dimension(:), allocatable :: h2psi_perturbed!,hpsi_tmp
  type(f_enumerator) :: inputpsi,output_denspot
  type(dictionary), pointer :: dict_timing_info
  type(orbital_basis) :: ob,ob_occ,ob_virt,ob_prime
  type(energy_terms) :: energs_fake
  real(kind=8),dimension(:,:),allocatable :: locreg_centers
  real(kind=8),dimension(:),allocatable :: charge_fake
  !Variable for TDHF
  real(gp) :: exc_fac !Factor in front of the term to subtract when doing TDHF or TDDFT with hybrid functionals

  ! testing
  real(kind=8),dimension(:,:),pointer :: locregcenters
  integer :: ilr, nlr, ioffset
  type(f_enumerator) :: linear_iscf
  integer,dimension(2) :: irow, icol, iirow, iicol
  character(len=20) :: comment
  !integer :: i1,i2,i3p,n1m,n2m,n3m !MM

  integer :: ishift, extra_states, i1, i2, i3, ii
  integer :: ind_min_s, ind_mas_s
  integer :: ind_min_m, ind_mas_m
  integer :: ind_min_l, ind_mas_l

  !debug
  !real(kind=8) :: ddot
  !proprietary interfaces
  interface
     subroutine input_wf(iproc,nproc,in,GPU,atoms,rxyz,&
          denspot,denspot0,nlpsp,KSwfn,tmb,energs,inputpsi,input_wf_format,norbv,&
          lzd_old,psi_old,rxyz_old,tmb_old,ref_frags,cdft,&
          locregcenters)
       use module_defs, only: gp,wp
       use f_enums, only: f_enumerator
       use module_types
       use module_fragments
       use constrained_dft
       implicit none
       integer, intent(in) :: iproc, nproc, input_wf_format
       type(f_enumerator), intent(in) :: inputpsi
       type(input_variables), intent(in) :: in
       type(GPU_pointers), intent(inout) :: GPU
       type(atoms_data), intent(inout) :: atoms
       real(gp), dimension(3, atoms%astruct%nat), target, intent(in) :: rxyz
       type(DFT_local_fields), intent(inout) :: denspot
       type(DFT_wavefunction), intent(inout) :: KSwfn,tmb,tmb_old !<input wavefunction
       type(energy_terms), intent(inout) :: energs !<energies of the system
       real(gp), dimension(*), intent(out) :: denspot0 !< Initial density / potential, if needed
       real(wp), dimension(:), pointer :: psi_old
       integer, intent(out) :: norbv
       type(DFT_PSP_projectors), intent(inout) :: nlpsp
       real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz_old
       type(local_zone_descriptors),intent(in):: lzd_old
       type(system_fragment), dimension(:), pointer :: ref_frags
       type(cdft_data), intent(out) :: cdft
       real(kind=8),dimension(3,atoms%astruct%nat),intent(in),optional :: locregcenters
     END SUBROUTINE input_wf
  end interface

  interface
     subroutine davidson(iproc,nproc,in,at,&
          orbs,orbsv,nvirt,Lzd,comms,commsv,&
          rxyz,rhopot,nlpsp,pkernel,psi,v,dpbox,xc,GPU)
       use module_defs, only: gp,dp,wp
       use module_types
       use communications_base, only: comms_cubic
       use module_xc
       use module_dpbox
       use locregs
       implicit none
       integer, intent(in) :: iproc,nproc
       integer, intent(in) :: nvirt
       type(input_variables), intent(in) :: in
       type(atoms_data), intent(in) :: at
       type(DFT_PSP_projectors), intent(inout) :: nlpsp
       type(local_zone_descriptors), intent(inout) :: Lzd
       type(orbitals_data), intent(inout) :: orbs
       type(comms_cubic), intent(in) :: comms, commsv
       type(denspot_distribution), intent(in) :: dpbox
       real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
       type(coulomb_operator), intent(inout) :: pkernel
       real(dp), dimension(*), intent(in) :: rhopot
       type(orbitals_data), intent(inout) :: orbsv
       type(GPU_pointers), intent(inout) :: GPU
       type(xc_info), intent(in) :: xc
       real(wp), dimension(:), pointer :: psi,v!=psivirt(nvctrp,nvirtep*nproc)
     END SUBROUTINE davidson
  end interface
  interface
     subroutine direct_minimization(iproc,nproc,in,at,nvirt,rxyz,&
          rhopot,nlpsp,pkernel,dpbox,xc,GPU,KSwfn,VTwfn)
       use module_defs, only: gp,dp,wp
       use module_types
       use module_xc
       use module_dpbox
       implicit none
       integer, intent(in) :: iproc,nproc,nvirt
       type(input_variables), intent(in) :: in
       type(atoms_data), intent(in) :: at
       type(DFT_PSP_projectors), intent(inout) :: nlpsp
       type(denspot_distribution), intent(in) :: dpbox
       type(DFT_wavefunction), intent(inout) :: KSwfn,VTwfn
       real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
       type(coulomb_operator), intent(inout) :: pkernel
       real(dp), dimension(*), intent(in), target :: rhopot
       type(GPU_pointers), intent(inout) :: GPU
       type(xc_info), intent(in) :: xc
     END SUBROUTINE direct_minimization
  end interface
  interface
     subroutine linearScaling(iproc,nproc,KSwfn,tmb,at,input,rxyz,denspot,rhopotold,nlpsp,GPU,&
          energs,energy,fpulay,infocode,ref_frags,cdft, &
          fdisp, fion)
       use module_defs, only: gp,dp,wp
       use module_types
       use module_fragments
       use constrained_dft
       implicit none
       integer,intent(in):: iproc, nproc
       type(atoms_data),intent(inout):: at
       type(input_variables),intent(in):: input
       real(8),dimension(3,at%astruct%nat),intent(inout):: rxyz
       real(8),dimension(3,at%astruct%nat),intent(out):: fpulay
       type(DFT_local_fields), intent(inout) :: denspot
       real(gp), dimension(*), intent(inout) :: rhopotold
       type(DFT_PSP_projectors),intent(inout):: nlpsp
       type(GPU_pointers),intent(in out):: GPU
       type(energy_terms),intent(inout) :: energs
       real(gp), dimension(:), pointer :: rho,pot
       real(8),intent(out):: energy
       type(DFT_wavefunction),intent(inout),target:: tmb
       type(DFT_wavefunction),intent(inout),target:: KSwfn
       integer,intent(out):: infocode
       type(system_fragment), dimension(:), pointer :: ref_frags
       type(cdft_data), intent(inout) :: cdft
       real(kind=8),dimension(3,at%astruct%nat),intent(in) :: fdisp, fion
     END SUBROUTINE linearScaling
  end interface

  call f_routine(id=subname)

  energs = energy_terms_null()
  energs_fake= energy_terms_null()

  DoLastRunThings=.false. !to avoid the implicit save attribute

  !copying the input variables for readability
  !this section is of course not needed
  !note that this procedure is convenient ONLY in the case of scalar variables
  !an array would have been copied, thus occupying more memory space
  !Hence WARNING: these variables are copied, in case of an update the new value should be
  !reassigned inside the structure

  if (iproc == 0) then
     !start a new document in the beginning of the output, if the document is closed before
     call yaml_new_document()
     call print_dft_parameters(in,atoms)
  end if

  !Time initialization
  call f_timing_reset(filename=trim(in%dir_output)//'time'+in%naming_id+'.yaml',master=iproc==0,&
       verbose_mode=get_verbose_level()>2 .and. nproc>1)
  call cpu_time(tcpu0)
  call system_clock(ncount0,ncount_rate,ncount_max)

  !!! Test MPI wrappers
  !!call test_mpi_alltoallv(iproc, nproc, bigdft_mpi%mpi_comm, &
  !!     maxsize_local=2000000000, ntest=5)

  !Nullify for new input guess
  call nullify_local_zone_descriptors(lzd_old)
  ! We save the variables that defined the previous psi if the restart is active
  inputpsi = in%inputPsiId
  if (in%inputPsiId == 'INPUT_PSI_MEMORY_WVL') then
     if (associated(KSwfn%psi)) then
        call copy_local_zone_descriptors(KSwfn%Lzd, lzd_old, subname)
        !if the history is bigger than two, create the workspace to store the wavefunction
        if (in%wfn_history > 0) then
!           if(iproc==0)print *, "NNdbg: KSwfn%Psi copied to Psi(",in%wfn_history+1,")"
           call old_wavefunction_set(KSwfn%oldpsis(in%wfn_history+1),&
                atoms%astruct%nat,KSwfn%orbs%norbp*KSwfn%orbs%nspinor,&
                KSwfn%Lzd,rxyz_old,KSwfn%psi)
           call f_free_ptr(KSwfn%psi)
        else
           call copy_old_wavefunctions(nproc,KSwfn%orbs,&
                KSwfn%psi,lzd_old%Glr%wfd,psi_old)
        end if
        !to maintain the same treatment destroy wfd afterwards (to be unified soon)
        !deallocation
        call deallocate_locreg_descriptors(KSwfn%Lzd%Glr)
     else
        !inputpsi = INPUT_PSI_LCAO
        call inputpsiid_set_policy(ENUM_SCRATCH,inputpsi)
     end if
  else if(in%inputPsiId == 'INPUT_PSI_MEMORY_GAUSS') then
     if (associated(KSwfn%psi)) then
        !deallocate wavefunction and descriptors for placing the gaussians
        call deallocate_locreg_descriptors(KSwfn%Lzd%Glr)
        call f_free_ptr(KSwfn%psi)
     else
        !inputpsi = INPUT_PSI_LCAO
        call set_inputpsiid(INPUT_PSI_LCAO,inputpsi) !special treatment here
     end if
  else if ( in%inputPsiId == 'INPUT_PSI_MEMORY_LINEAR' .and. associated(KSwfn%psi)) then !in%inputPsiId == INPUT_PSI_MEMORY_LINEAR_
     if (associated(tmb%psi)) then
        !this is a tmb_null constructor
        tmb_old%lzd = local_zone_descriptors_null()
        !!tmb_old%linmat%smmd = sparse_matrix_metadata_null()
        !!tmb_old%linmat%smat(1) = sparse_matrix_null()
        !!tmb_old%linmat%smat(2) = sparse_matrix_null()
        !!tmb_old%linmat%smat(3) = sparse_matrix_null()
        !!!!tmb_old%linmat%ks = sparse_matrix_null()
        !!!!tmb_old%linmat%ks_e = sparse_matrix_null()
        !!nullify(tmb_old%linmat%ks)
        !!nullify(tmb_old%linmat%ks_e)
        !!tmb_old%linmat%ovrlp_ = matrices_null()
        !!tmb_old%linmat%ham_ = matrices_null()
        !!tmb_old%linmat%kernel_ = matrices_null()
        !!do i=1,size(tmb_old%linmat%ovrlppowers_)
        !!    tmb_old%linmat%ovrlppowers_(i) = matrices_null()
        !!end do
        tmb_old%linmat = linear_matrices_null()
        tmb_old%collcom = comms_linear_null()
        call copy_tmbs(iproc, tmb, tmb_old, subname)
        call destroy_DFT_wavefunction(tmb)
        call f_free_ptr(KSwfn%psi)
        call deallocate_locreg_descriptors(KSwfn%Lzd%Glr)
        !call deallocate_wfd(KSwfn%Lzd%Glr%wfd)
     else
        !inputpsi = INPUT_PSI_LINEAR_AO
        call inputpsiid_set_policy(ENUM_SCRATCH,inputpsi)
     end if
  end if

  ! Setup all descriptors and allocate what should be.
  if (inputpsi .hasattr. 'LINEAR') then
     locregcenters=f_malloc_ptr((/3,atoms%astruct%nat/),id=' locregcenters')
      if (in%explicit_locregcenters) then
          open(unit=123, file='locregcenters.xyz')
          read(123,*) nlr
          if (nlr/=atoms%astruct%nat) stop 'ERROR: wrong nlr'
          read(123,*) comment
          do ilr=1,nlr
              read(123,*) comment, locregcenters(1,ilr), locregcenters(2,ilr), locregcenters(3,ilr)
          end do
      else
          locregcenters = rxyz
      end if
  end if

  if(f_int(inputpsi) == INPUT_PSI_MEMORY_LINEAR) then
    call system_initialization(iproc,nproc,.true.,inputpsi,input_wf_format,.false.,in,atoms,rxyz,GPU%OCLconv,&
         KSwfn%orbs,tmb%npsidim_orbs,tmb%npsidim_comp,tmb%orbs,KSwfn%Lzd,tmb%Lzd,nlpsp,&
         KSwfn%comms,ref_frags,denspot=denspot,locregcenters=locregcenters,&
         inwhichlocreg_old=tmb_old%orbs%inwhichlocreg,onwhichatom_old=tmb_old%orbs%onwhichatom,&
         norb_par_ref=tmb_old%orbs%norb_par, norbu_par_ref=tmb_old%orbs%norbu_par, norbd_par_ref=tmb_old%orbs%norbd_par)
  else if(inputpsi == 'INPUT_PSI_LINEAR_AO' .or. inputpsi == 'INPUT_PSI_DISK_LINEAR') then
    call system_initialization(iproc,nproc,.true.,inputpsi,input_wf_format,.false.,in,atoms,rxyz,GPU%OCLconv,&
         KSwfn%orbs,tmb%npsidim_orbs,tmb%npsidim_comp,tmb%orbs,KSwfn%Lzd,tmb%Lzd,nlpsp,&
         KSwfn%comms,ref_frags,denspot=denspot,locregcenters=locregcenters)
  else
    call system_initialization(iproc,nproc,.true.,inputpsi,input_wf_format,&
         & .false.,in,atoms,rxyz,GPU%OCLconv,&
         KSwfn%orbs,tmb%npsidim_orbs,tmb%npsidim_comp,tmb%orbs,KSwfn%Lzd,tmb%Lzd,nlpsp,&
         KSwfn%comms,ref_frags,denspot=denspot)
  end if

  ! Ugly here to be moved elsewhere.
  if (associated(atoms%pawtab)) then
     call libxc_functionals_init(denspot%xc%ixc, denspot%dpbox%nrhodim)
  end if

  !memory estimation, to be rebuilt in a more modular way
  call MemoryEstimator(nproc,in%idsx,KSwfn%Lzd%Glr,&
       KSwfn%orbs%norb,KSwfn%orbs%nspinor,KSwfn%orbs%nkpts,&
       nlpsp%nprojel,in%nspin,in%itrpmax,f_int(in%scf),mem)
  if (.not.(inputpsi .hasattr. 'LINEAR') .and. iproc==0 .and. get_verbose_level() > 0) then
      call print_memory_estimation(mem)
  end if

  if (in%lin%fragment_calculation .and. inputpsi == 'INPUT_PSI_DISK_LINEAR') then
     call output_fragment_rotations(iproc,atoms%astruct%nat,rxyz,1,trim(in%dir_output),in%frag,ref_frags)
  end if

  ! temporary, really want to just initialize it here rather than copy
  ! but still need to move all cubic references to KSwfn%orbs%npsidim to just KSwfn%npsidim
  KSwfn%npsidim_orbs = KSwfn%orbs%npsidim_orbs
  KSwfn%npsidim_comp = KSwfn%orbs%npsidim_comp

  ! We complete here the definition of DFT_wavefunction structures.
  if (inputpsi .hasattr. 'LINEAR') then

     !!call init_p2p_tags(nproc)
     !!tag=0

     call kswfn_init_comm(tmb, denspot%dpbox, iproc, nproc, in%nspin, in%imethod_overlap)

     !!locreg_centers = f_malloc((/3,tmb%lzd%nlr/),id='locreg_centers')
     !!do ilr=1,tmb%lzd%nlr
     !!    locreg_centers(1:3,ilr)=tmb%lzd%llr(ilr)%locregcenter(1:3)
     !!end do
     call init_foe_wrapper(iproc, nproc, in, KSwfn%orbs, 0.d0, tmb%foe_obj)
     ! Do the same for the object which handles the calculation of the inverse.
     charge_fake = f_malloc0(in%nspin,id='charge_fake')
     call init_foe(iproc, nproc, in%nspin, charge_fake, tmb%ice_obj, &
          tmprtr=0.d0, evbounds_nsatur=in%cp%foe%evbounds_nsatur, &
          evboundsshrink_nsatur=in%cp%foe%evboundsshrink_nsatur, &
          evlow=0.5d0, evhigh=1.5d0, fscale=in%cp%foe%fscale, &
          ef_interpol_det=in%cp%foe%ef_interpol_det, &
          ef_interpol_chargediff=in%cp%foe%ef_interpol_chargediff, &
          fscale_lowerbound=in%cp%foe%fscale_lowerbound, &
          fscale_upperbound=in%cp%foe%fscale_upperbound, &
          eval_multiplicator=1.d0, &
          accuracy_function=in%cp%foe%accuracy_ice, accuracy_penalty=in%cp%foe%accuracy_penalty, &
          betax=in%cp%foe%betax_ice, occupation_function=in%cp%foe%occupation_function, &
          adjust_fscale=in%cp%foe%adjust_fscale, &
          fscale_ediff_low=in%cp%foe%fscale_ediff_low, &
          fscale_ediff_up=in%cp%foe%fscale_ediff_up)
     call f_free(charge_fake)

     !!call f_free(locreg_centers)
     !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, atoms%astruct, in, KSwfn%orbs, tmb%orbs, tmb%foe_obj, .true.)

     call create_large_tmbs(iproc, nproc, KSwfn, tmb, denspot,nlpsp,in, atoms, rxyz, .false.)

     if (iproc==0) then
         call yaml_mapping_open('Workarray memory requirements for transposed communication')
     end if
     call write_memory_requirements_collcom(iproc, nproc, bigdft_mpi%mpi_comm, 'Normal locregs', tmb%collcom)
     call write_memory_requirements_collcom(iproc, nproc, bigdft_mpi%mpi_comm, 'Normal locregs sumrho', tmb%collcom_sr)
     call write_memory_requirements_collcom(iproc, nproc, bigdft_mpi%mpi_comm, 'Large locregs', tmb%ham_descr%collcom)
     if (iproc==0) then
         call yaml_mapping_close()
     end if

     call init_bigdft_matrices(iproc, nproc, atoms, tmb%orbs, &
          tmb%collcom, tmb%collcom_sr, tmb%ham_descr%collcom, &
          tmb%lzd, tmb%ham_descr%lzd, in, tmb%linmat)


     if (iproc==0) then
         call yaml_mapping_open('Matrices')
         call write_sparsematrix_info(tmb%linmat%smat(1), 'Overlap matrix')
         call write_sparsematrix_info(tmb%linmat%smat(2), 'Hamiltonian matrix')
         call write_sparsematrix_info(tmb%linmat%smat(3), 'Density kernel matrix')
         call yaml_mapping_close()
     end if

     tmb%linmat%kernel_ = matrices_null()
     tmb%linmat%ham_ = matrices_null()
     tmb%linmat%ovrlp_ = matrices_null()
     !!call allocate_matrices(tmb%linmat%smat(3), allocate_full=.false., &
     !!     matname='tmb%linmat%kernel_', mat=tmb%linmat%kernel_)
     !!call allocate_matrices(tmb%linmat%smat(2), allocate_full=.false., &
     !!     matname='tmb%linmat%ham_', mat=tmb%linmat%ham_)
     !!call allocate_matrices(tmb%linmat%smat(1), allocate_full=.false., &
     !!     matname='tmb%linmat%ovrlp_', mat=tmb%linmat%ovrlp_)
     tmb%linmat%kernel_%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%smat(3), &
         iaction=SPARSE_TASKGROUP,id='tmb%linmat%kernel_%matrix_compr')
     ! SM: Probably not used anymore
     !!if (tmb%linmat%smat(3)%smmm%nfvctrp>tmb%linmat%smat(3)%nfvctrp) then
     !!    tmb%linmat%kernel_%matrixp = sparsematrix_malloc_ptr(tmb%linmat%smat(3), &
     !!        iaction=DENSE_MATMUL,id='tmb%linmat%kernel_%matrixp')
     !!else
     !!    tmb%linmat%kernel_%matrixp = sparsematrix_malloc_ptr(tmb%linmat%smat(3), &
     !!        iaction=DENSE_PARALLEL,id='tmb%linmat%kernel_%matrixp')
     !!end if
     tmb%linmat%ham_%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%smat(2), &
         iaction=SPARSE_TASKGROUP,id='tmb%linmat%ham_%matrix_compr')
     tmb%linmat%ovrlp_%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%smat(1), &
         iaction=SPARSE_TASKGROUP,id='tmb%linmat%ovrlp_%matrix_compr')

     if (iproc==0) then
         call yaml_comment('Unitary tests',hfill='-')
         call yaml_mapping_open('Results of unitary tests')
     end if

     if (in%check_matrix_compression) then
         if (iproc==0) call yaml_mapping_open('Checking Compression/Uncompression of small sparse matrices')
         !call check_matrix_compression(iproc,tmb%linmat%ham)
         !call check_matrix_compression(iproc,tmb%linmat%ovrlp)
         call check_matrix_compression(iproc, nproc, tmb%linmat%smat(2), tmb%linmat%ham_)
         call check_matrix_compression(iproc, nproc, tmb%linmat%smat(1), tmb%linmat%ovrlp_)
         if (iproc ==0) call yaml_mapping_close()
     end if

     do i=1,size(tmb%linmat%ovrlppowers_)
         tmb%linmat%ovrlppowers_(i) = matrices_null()
         !call allocate_matrices(tmb%linmat%smat(3), allocate_full=.false., &
         !     matname='tmb%linmat%ovrlppowers_(i)', mat=tmb%linmat%ovrlppowers_(i))
         tmb%linmat%ovrlppowers_(i)%matrix_compr = &
             sparsematrix_malloc_ptr(tmb%linmat%smat(3), iaction=SPARSE_TASKGROUP, id='tmb%linmat%ovrlppowers_(i)%matrix_compr')
     end do

     !call nullify_sparse_matrix(tmb%linmat%inv_ovrlp_large)
     !tmb%linmat%inv_ovrlp_large=sparse_matrix_null()
     !call sparse_copy_pattern(tmb%linmat%smat(3), tmb%linmat%inv_ovrlp_large, iproc, subname)

     ! Initializes a sparse matrix type compatible with the ditribution of the
     ! KS orbitals. This is required for the re-orthonromalization of the
     ! KS espansion coefficients, so it is not necessary for FOE.
     nullify(tmb%linmat%ks)
     nullify(tmb%linmat%ks_e)
     if (in%lin%scf_mode/=LINEAR_FOE .or. &
         (mod(in%lin%plotBasisFunctions,10) /= WF_FORMAT_NONE) .or. in%lin%diag_end .or. in%output_wf /= ENUM_EMPTY  &
         .or. f_int(inputpsi) == INPUT_PSI_DISK_LINEAR) then

         if (in%lin%fragment_calculation) then
            if (in%lin%extra_states/=0) then
               !if (iproc==0)
               call f_err_throw('ERROR: extra_states must be zero for fragment calculation',&
                    err_name='BIGDFT_INPUT_VARIABLES_ERROR')
               !call MPI_ABORT(bigdft_mpi%mpi_comm,10,ierr)
            end if
            extra_states = tmb%orbs%norb - KSwfn%orbs%norb
         else
            extra_states = in%lin%extra_states
         end if
         call init_sparse_matrix_for_KSorbs(iproc, nproc, &
              KSwfn%orbs, in, atoms%astruct%geocode, &
              atoms%astruct%cell_dim, extra_states, tmb%linmat%ks, tmb%linmat%ks_e)
     end if


     if (in%check_matrix_compression) then
         if (iproc==0) then
             call yaml_mapping_open('Checking Compression/Uncompression of large sparse matrices')
         end if
         call check_matrix_compression(iproc, nproc, tmb%linmat%smat(3), tmb%linmat%kernel_)
         if (iproc ==0) then
             call yaml_mapping_close()
         end if
     end if

     if (in%check_sumrho>0) then
         call check_communication_potential(iproc,nproc,denspot,tmb)
         call check_communication_sumrho(iproc, nproc, tmb%orbs, tmb%lzd, tmb%collcom_sr, &
              denspot, tmb%linmat%smat(3), tmb%linmat%auxl, tmb%linmat%kernel_, in%check_sumrho)
     end if

     if (iproc==0) then
         call yaml_mapping_open('Checking Communications of Minimal Basis')
     end if
     call check_communications_locreg(iproc,nproc,tmb%orbs,in%nspin,tmb%lzd, &
          tmb%collcom,tmb%linmat%smat(1),tmb%linmat%auxs,tmb%linmat%ovrlp_, &
          tmb%npsidim_orbs,tmb%npsidim_comp,in%check_overlap)
     if (iproc==0) then
         call yaml_mapping_close()
     end if

     if (iproc==0) then
         call yaml_mapping_open('Checking Communications of Enlarged Minimal Basis')
     end if
     call check_communications_locreg(iproc,nproc,tmb%orbs,in%nspin,tmb%ham_descr%lzd, &
          tmb%ham_descr%collcom,tmb%linmat%smat(2),tmb%linmat%auxm,tmb%linmat%ham_, &
          tmb%ham_descr%npsidim_orbs,tmb%ham_descr%npsidim_comp,in%check_overlap)
     if (iproc ==0) then
         call yaml_mapping_close()
         call yaml_mapping_close()
     end if

     if (in%lin%scf_mode/=LINEAR_FOE .or. &
         (mod(in%lin%plotBasisFunctions,10) /= WF_FORMAT_NONE) .or. in%lin%diag_end .or. in%output_wf /= ENUM_EMPTY &
          .or. f_int(inputpsi) == INPUT_PSI_DISK_LINEAR .or. mod(in%lin%output_coeff_format,10) /= WF_FORMAT_NONE) then
        tmb%coeff = f_malloc_ptr((/ tmb%linmat%smat(2)%nfvctr , tmb%orbs%norb /),id='tmb%coeff')
     else
        nullify(tmb%coeff)
     end if

     denspot0 = f_malloc(max(denspot%dpbox%ndimrhopot, denspot%dpbox%nrhodim),id='denspot0')
  else
     denspot0 = f_malloc(1,id='denspot0')
  end if

  !the lookup tables for the application of the nonlocal potential can be created from now on
  optLoop%scf = in%scf
  optLoop%itrpmax = in%itrpmax
  optLoop%nrepmax = in%nrepmax
  optLoop%itermax = in%itermax
  optLoop%itermin = in%itermin!Bastian
  optLoop%gnrm_cv = in%gnrm_cv
  optLoop%rpnrm_cv = in%rpnrm_cv
  optLoop%gnrm_startmix = in%gnrm_startmix
  optLoop%itrp = 0
  optLoop%itrep = 0
  optLoop%iter = 0
  optLoop%infocode = 0

  call system_signaling(iproc,in%signaling,in%gmainloop,&
       KSwfn,tmb,energs,denspot,optloop,&
       atoms%astruct%ntypes,atoms%radii_cf,in%crmult,in%frmult)

  !variables substitution for the PSolver part
  n1=KSwfn%Lzd%Glr%d%n1
  n2=KSwfn%Lzd%Glr%d%n2
  n3=KSwfn%Lzd%Glr%d%n3

  !calculate the rhocore contribution to the energy value
  if (associated(denspot%rho_C)) then
     !calculate the XC energy of rhocore, use the rhov array as a temporary variable
     !use Vxc and other quantities as local variables
     call xc_init_rho(denspot%xc, denspot%dpbox%nrhodim,denspot%rhov,1)
     denspot%rhov=1.d-16
     call XC_potential(atoms%astruct%geocode,'D',denspot%pkernel%mpi_env%iproc,denspot%pkernel%mpi_env%nproc,&
          denspot%pkernel%mpi_env%mpi_comm,&
          denspot%dpbox%mesh%ndims(1),denspot%dpbox%mesh%ndims(2),denspot%dpbox%mesh%ndims(3),denspot%xc,&
          denspot%dpbox%mesh%hgrids,&
          denspot%rhov,energs%excrhoc,tel,KSwfn%orbs%nspin,denspot%rho_C,&
          denspot%rhohat,denspot%V_XC,xcstr)
     if (iproc==0) call yaml_map('Value for Exc[rhoc]',energs%excrhoc)
     !if (iproc==0) write(*,*)'value for Exc[rhoc]',energs%excrhoc
  end if

  !Calculate the ionic energy and forces accordingly
  call IonicEnergyandForces(iproc,nproc,denspot%dpbox,atoms,in%elecfield,rxyz,&
       energs%eion,fion,in%dispersion,energs%edisp,fdisp,ewaldstr,&
       denspot%V_ext,denspot%pkernel,denspot%psoffset)
  !Calculate effective ionic potential, including counter ions if any.
  call createEffectiveIonicPotential(iproc,(iproc == 0),in,atoms,rxyz,atoms%astruct%shift,&
       denspot%dpbox,denspot%pkernel,denspot%V_ext,denspot%rho_ion,in%elecfield,denspot%psoffset)
  call potential_from_charge_multipoles(iproc, nproc, atoms, denspot, in%ep, 1, denspot%dpbox%mesh%ndims(1), 1,&
       denspot%dpbox%mesh%ndims(2), &
       denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+1, &
       denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,2), &
       denspot%dpbox%mesh%hgrids(1),denspot%dpbox%mesh%hgrids(2),denspot%dpbox%mesh%hgrids(3), atoms%astruct%shift, verbosity=1, &
       ixc=in%ixc, lzd=tmb%lzd, pot=denspot%V_ext, &
       rxyz=rxyz, ixyz0=in%plot_mppot_axes, write_directory=trim(in%dir_output))
  call interaction_multipoles_ions(bigdft_mpi%iproc, KSwfn%Lzd%Glr%mesh, in%ep, atoms, energs%eion, fion)
!!$  call interaction_multipoles_ions(bigdft_mpi%iproc, in%ep, atoms, energs%eion, fion)
  !write(*,*) 'eion before', energs%eion
  call ionic_energy_of_external_charges(bigdft_mpi%iproc, KSwfn%Lzd%Glr%mesh, in%ep, atoms, energs%eion)
!!$  call ionic_energy_of_external_charges(bigdft_mpi%iproc, in%ep, atoms, energs%eion)
  !write(*,*) 'eion after', energs%eion
  !call yaml_map('rxyz',rxyz)
  !call yaml_map('atoms%astruct%rxyz',atoms%astruct%rxyz)
  !call yaml_map('hgrids',(/denspot%dpbox%mesh%hgrids(1),denspot%dpbox%mesh%hgrids(2),denspot%dpbox%mesh%hgrids(3)/))
  !!call potential_from_multipoles(in%ep, 1, denspot%dpbox%mesh%ndims(1), 1, denspot%dpbox%mesh%ndims(2), &
  !!     denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+1, &
  !!     denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,3)+denspot%dpbox%nscatterarr(denspot%dpbox%mpi_env%iproc,2), &
  !!     1, denspot%dpbox%mesh%ndims(3), &
  !!     denspot%dpbox%mesh%hgrids(1),denspot%dpbox%mesh%hgrids(2),denspot%dpbox%mesh%hgrids(3), shift, denspot%V_ext)
  if (denspot%c_obj /= 0) then
     call denspot_emit_v_ext(denspot, iproc, nproc)
  end if

  !ii = 0
  !do i3=1,denspot%dpbox%mesh%ndims(3)
  !    do i2=1,denspot%dpbox%mesh%ndims(2)
  !        do i1=1,denspot%dpbox%mesh%ndims(1)
  !            ii = ii + 1
  !            write(300,*) 'vals', i1, i2, i3, denspot%V_ext(i1,i2,i3,1)
  !        end do
  !    end do
  !end do

  norbv=abs(in%norbv)
  if (in%inputPsiId .hasattr. 'LINEAR') then
     ! Setup the mixing, if necessary -- NEW
     if (in%lin%mixHist_lowaccuracy /= in%lin%mixHist_highaccuracy) then
         ! This must be fixed later
        call f_err_throw('in%lin%mixHist_lowaccuracy /= in%lin%mixHist_highaccuracy',&
             err_name='BIGDFT_RUNTIME_ERROR')
     end if

     if (in%lin%mixHist_lowaccuracy==0) then
        ! simple mixing
        linear_iscf=f_enumerator('METHOD',AB7_MIXING_SIMPLE,null())
     else
        ! Pulay mixing
        linear_iscf=f_enumerator('METHOD',AB7_MIXING_PULAY,null())
     end if

     select case (in%lin%scf_mode)
     case (LINEAR_MIXPOT_SIMPLE)
        call f_enum_attr(linear_iscf,POT_MIX_ENUM)
     case (LINEAR_FOE,LINEAR_PEXSI,LINEAR_DIRECT_MINIMIZATION,&
          LINEAR_MIXDENS_SIMPLE)
        call f_enum_attr(linear_iscf,DEN_MIX_ENUM)
     case default
         stop 'ERROR: wrong in%lin%scf_mode'
     end select
     call f_enum_attr(linear_iscf,RSPACE_MIX_ENUM)
     call f_enum_attr(linear_iscf,MIX_ENUM)
     call denspot_set_history(denspot,linear_iscf, &
          & npulayit=in%lin%mixHist_lowaccuracy)
     tmb%damping_factor_confinement = 1.d0 !will be modified in case of a restart
     call input_wf(iproc,nproc,in,GPU,atoms,rxyz,denspot,denspot0,nlpsp,KSwfn,tmb,energs,&
          inputpsi,input_wf_format,norbv,lzd_old,psi_old,rxyz_old,tmb_old,ref_frags,cdft,&
          locregcenters)
      call f_free_ptr(locregcenters)
  else
      call input_wf(iproc,nproc,in,GPU,atoms,rxyz,denspot,denspot0,nlpsp,KSwfn,tmb,energs,&
           inputpsi,input_wf_format,norbv,lzd_old,psi_old,rxyz_old,tmb_old,ref_frags,cdft)
   end if

  nvirt=max(in%nvirt,norbv)

  ! modified by SM
  call deallocate_local_zone_descriptors(lzd_old)

  !end of the initialization part
  call timing(bigdft_mpi%mpi_comm,'INIT','PR')

  !start the optimization
  energs%eexctX=0.0_gp
  ! Skip the following part in the linear scaling case.
  skip_if_linear: if (inputpsi .hasattr. 'CUBIC') then
     call kswfn_optimization_loop(iproc, nproc, optLoop, &
     & in%alphamix, in%idsx, inputpsi, KSwfn, denspot, nlpsp, energs, atoms, GPU, xcstr, &
     & in)
     infocode = optLoop%infocode

     !if we are in the last_run case, validate the last_run only for the last cycle
     !do the last_run things regardless of infocode
     !nrepmax=0 is needed for the Band Structure calculations
     DoLastRunThings=(in%last_run == 1 .and. optLoop%nrepmax == 0) .or. &
          & (in%last_run == 1 .and. optLoop%itrep >= optLoop%nrepmax)
              !print the energies only if they are meaningful
     energy = energs%energy
     !Davidson is set to false first because used in deallocate_before_exiting
     DoDavidson= .false.

     ! Treat the info code from the optimization routine.
     if (infocode == 2 .or. infocode == 3) then
        !call deallocate_convolutions_bounds(KSwfn%lzd%glr%bounds)
        call deallocate_before_exiting
        return
     end if
  else

     ! Allocation of array for Pulay forces (only needed for linear version)
     fpulay = f_malloc_ptr((/ 3, atoms%astruct%nat /),id='fpulay')

     call linearScaling(iproc,nproc,KSwfn,tmb,atoms,in,&
          rxyz,denspot,denspot0,nlpsp,GPU,energs,energy,fpulay,infocode,ref_frags,cdft,&
          fdisp,fion)

     ! Clean denspot parts only needed in the SCF loop -- NEW
     call denspot_free_history(denspot)

     ! maybe not the best place to keep it - think about it!
     if (in%lin%calc_transfer_integrals) then
        if (in%lin%constrained_dft) then
           ! switch excess charge to other fragment, recalculate kernel and density and reset lagrange multiplier
           if (iproc==0) write(*,*) '--------------------------------------------------------------------------------------'
           if (iproc==0) write(*,*) 'Warning: site-energy/transfer integral calculation not yet working for constrained DFT'
           if (iproc==0) write(*,*) '--------------------------------------------------------------------------------------'

           !in_frag_charge=f_malloc_ptr(in%frag%nfrag,id='in_frag_charge')
           !call vcopy(in%frag%nfrag,in%frag%charge(1),1,in_frag_charge(1),1)
           !! assume all other fragments neutral, use total system charge to get correct charge for the other fragment
           !in_frag_charge(cdft%ifrag_charged(1))=in%ncharge - in_frag_charge(cdft%ifrag_charged(2))
           !overlap_calculated=.true.
           !call fragment_coeffs_to_kernel(iproc,in%frag,in_frag_charge,ref_frags,tmb,KSwfn%orbs,overlap_calculated)
           !call f_free_ptr(in_frag_charge)
           !cdft%charge=-cdft%charge

           !call reconstruct_kernel(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orthpar%blocksize_pdsyev, tmb%orthpar%blocksize_pdgemm, &
           !     KSwfn%orbs, tmb, overlap_calculated)
           !tmb%can_use_transposed=.false. ! - do we really need to deallocate here?
           !i_all = -product(shape(tmb%psit_c))*kind(tmb%psit_c)
           !deallocate(tmb%psit_c,stat=i_stat)
           !call memocc(i_stat,i_all,'tmb%psit_c',subname)
           !i_all = -product(shape(tmb%psit_f))*kind(tmb%psit_f)
           !deallocate(tmb%psit_f,stat=i_stat)
           !call memocc(i_stat,i_all,'tmb%psit_f',subname)

           !! Now need to calculate the charge density and the potential related to this inputguess
           !call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, max(tmb%npsidim_orbs,tmb%npsidim_comp), &
           !     tmb%orbs, tmb%psi, tmb%collcom_sr)

           !call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%mesh%hgrids(1), KSwfn%Lzd%mesh%hgrids(2), KSwfn%Lzd%mesh%hgrids(3), &
           !     tmb%collcom_sr, tmb%linmat%denskern, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)

           !! Must initialize rhopotold (FOR NOW... use the trivial one)
           !call vcopy(max(denspot%dpbox%mesh%ndims(1)*denspot%dpbox%mesh%ndims(2)*denspot%dpbox%n3p,1)*in%nspin, &
           !     denspot%rhov(1), 1, denspot0(1), 1)
           !!!call deallocateCommunicationbufferSumrho(tmb%comsr, subname)
           !call updatePotential(in%ixc,in%nspin,denspot,energs%eh,energs%exc,energs%evxc)
           !call local_potential_dimensions(tmb%lzd,tmb%orbs,denspot%dpbox%ngatherarr(0,1))

           !! keep a copy of previous wavefunctions and energies...
           !allocate(psi_constrained(tmb%npsidim_orbs), stat=i_stat)
           !call memocc(i_stat, psi_constrained, 'psi_constrained', subname)
           !call vcopy(tmb%npsidim_orbs,tmb%psi(1),1,psi_constrained(1),1)
           !energy_constrained=energy

           !call linearScaling(iproc,nproc,KSwfn,tmb,atoms,in,&
           !     rxyz,denspot,denspot0,nlpsp,GPU,energs,energy,fpulay,infocode,ref_frags,cdft)

           !! calculate matrix elements here...

           !i_all=-product(shape(psi_constrained))*kind(psi_constrained)
           !deallocate(psi_constrained, stat=i_stat)
           !call memocc(i_stat, i_all, 'psi_constrained', subname)

        else
           if (.not. in%lin%fragment_calculation) stop 'Error, fragment calculation needed for transfer integral calculation'
           !if (input%frag%nfrag==2) call calc_transfer_integrals_old(iproc,nproc,input%frag,ref_frags,tmb%orbs,&
           !     tmb%linmat%ham,tmb%linmat%ovrlp)
           call calc_site_energies_transfer_integrals(iproc,nproc,in%lin%order_taylor,&
                in%frag,ref_frags,tmb%orbs,tmb%linmat%smat(2),tmb%linmat%ham_,tmb%linmat%smat(1),tmb%linmat%ovrlp_,&
                tmb%linmat%ks)
        end if
     end if

     ! deallocate fragments
      if (in%lin%fragment_calculation) then ! we really need to deallocate
         do ifrag=1,in%frag%nfrag_ref
            call fragment_free(ref_frags(ifrag))
         end do
        deallocate(ref_frags)
     else if (f_int(inputpsi) == INPUT_PSI_DISK_LINEAR) then! we haven't actually allocated anything, so can just nullify - should make this more robust/general
         do ifrag=1,in%frag%nfrag_ref
            call fragment_free(ref_frags(ifrag))
            ref_frags(ifrag)%astruct_frg%nat=-1
            ref_frags(ifrag)%fbasis%forbs=minimal_orbitals_data_null()
            !ref_frags(ifrag)=fragment_null()
         end do
        deallocate(ref_frags)
     end if

     denspot%rho_work = f_malloc_ptr(max(denspot%dpbox%ndimrhopot,1),id='denspot%rho_work')
     !then calculate the contribution coming from the Hartree part
     ! keep only the essential part of the density, without the GGA bufffers
     ioffset=kswfn%lzd%glr%d%n1i*kswfn%lzd%glr%d%n2i*denspot%dpbox%i3xcsh
     !write(*,*) 'ioffset', ioffset
     if (denspot%dpbox%ndimrhopot>0) then
        call f_memcpy(n=denspot%dpbox%ndimpot,src=denspot%rhov(ioffset+1),&
             dest=denspot%rho_work(1))
        ! add the spin down part if present
        if (denspot%dpbox%nrhodim==2) then
           ishift=denspot%dpbox%ndimrhopot/denspot%dpbox%nrhodim !start of the spin down part
           call axpy(denspot%dpbox%ndimpot, 1.d0, &
                denspot%rhov(ioffset+ishift+1), &
                1, denspot%rho_work(1),1)
        end if
     end if

     if (infocode==2) then
        !!! Allocate this array since it will be deallocated in deallocate_before_exiting
        call f_free_ptr(fpulay)
        call destroy_DFT_wavefunction(tmb)
        call f_free_ptr(KSwfn%psi)
        !call deallocate_wfd(KSwfn%Lzd%Glr%wfd)
        call deallocate_locreg_descriptors(KSwfn%Lzd%Glr)
        call f_free_ptr(denspot%rho_work)
        call f_free_ptr(KSwfn%orbs%eval)
        call deallocate_before_exiting()
        return
     end if

     !infocode = 0
  end if skip_if_linear

  !last run things has to be done:
  !if it is the last run and the infocode is zero
  !if infocode is not zero but the last run has been done for nrepmax times

  DoLastRunThings= (in%last_run == 1 .and. infocode == 0) .or. DoLastRunThings

  !analyse the possibility to calculate Davidson treatment
  !(nvirt > 0 .and. in%inputPsiId == 0)
  DoDavidson= (abs(in%norbv) > 0 .or. in%nkptv > 0) .and. DoLastRunThings

  !project the wavefunctions on a gaussian basis and keep in memory
  !if (in%gaussian_help) then
  if (in%inputPsiId .hasattr. 'GAUSSIAN') then
     call timing(iproc,'gauss_proj','ON')
     if (iproc == 0.and.get_verbose_level() >1) then
        call yaml_comment('Gaussian Basis Projection',hfill='-')
     end if

     !extract the gaussian basis from the pseudowavefunctions
     call gaussian_pswf_basis(21,.false.,iproc,in%nspin,atoms,rxyz,KSwfn%gbd,gbd_occ)
     call f_free_ptr(gbd_occ)

     if (.not. associated(KSwfn%gaucoeffs)) then
        KSwfn%gaucoeffs = f_malloc_ptr((/ KSwfn%gbd%ncoeff, KSwfn%orbs%nspinor*KSwfn%orbs%norbp /),id='KSwfn%gaucoeffs')
     end if

     thetaphi = f_malloc0((/ 2, KSwfn%gbd%nat /),id='thetaphi')
     !thetaphi=0.0_gp
!!$     call wavelets_to_gaussians(atoms%astruct%geocode,KSwfn%orbs%norbp,KSwfn%orbs%nspinor,&
     call wavelets_to_gaussians(KSwfn%Lzd%Glr,KSwfn%orbs%norbp,KSwfn%orbs%nspinor,&
          n1,n2,n3,KSwfn%gbd,thetaphi,&
          KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
          KSwfn%Lzd%Glr%wfd,KSwfn%psi,KSwfn%gaucoeffs)

     call f_free(thetaphi)
     call timing(iproc,'gauss_proj','OF')
  end if

  !  write all the wavefunctions into files
  !if (in%output_wf_format /= WF_FORMAT_NONE .and. DoLastRunThings) then
  if (in%output_wf /= ENUM_EMPTY .and. DoLastRunThings) then
     !add flag for writing waves in the gaussian basis form
     if ((in%inputPsiId .hasattr. 'GAUSSIAN') .and. (in%inputPsiId .hasattr. 'CUBIC')) then

        !write the coefficients and the basis on a file
        if (iproc ==0) call yaml_map('Writing wavefunctions in file','wavefunction.gau')
        call write_gaussian_information(iproc,nproc,KSwfn%orbs,KSwfn%gbd,KSwfn%gaucoeffs,trim(in%dir_output) // 'wavefunctions.gau')

        !build dual coefficients
        call dual_gaussian_coefficients(KSwfn%orbs%norbp*KSwfn%orbs%nspinor,KSwfn%gbd,KSwfn%gaucoeffs)

        !control the accuracy of the expansion
        call check_gaussian_expansion(iproc,nproc,KSwfn%orbs,KSwfn%Lzd,KSwfn%psi,KSwfn%gbd,KSwfn%gaucoeffs)

        call deallocate_gwf(KSwfn%gbd)
        call f_free_ptr(KSwfn%gaucoeffs)
        nullify(KSwfn%gbd%rxyz)

     else if (inputpsi .hasattr. 'CUBIC') then
        call writemywaves(iproc,trim(in%dir_output) // trim(in%outputpsiid), f_int(in%output_wf), &
             KSwfn%orbs,n1,n2,n3,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
             atoms,rxyz,KSwfn%Lzd%Glr%wfd,KSwfn%psi,paw = KSwfn%paw)
     end if
  end if

  call f_zero(gridformat)
  if (in%output_denspot .hasattr. 'ETSF') then
     call f_strcpy(src='.etsf',dest=gridformat)
  else if (in%output_denspot .hasattr. 'CUBE') then
     call f_strcpy(src='.cube',dest=gridformat)
  end if

  !plot the ionic potential, if required by output_denspot
  if (in%output_denspot == 'DENSPOT' .and. DoLastRunThings) then
     if (all(in%plot_pot_axes>=0)) then
        if (iproc == 0) call yaml_map('Writing external potential in file', 'external_potential'//gridformat)
        call plot_density(iproc,nproc,trim(in%dir_output)//'external_potential' // gridformat,&
             atoms,rxyz,denspot%pkernel,1,denspot%V_ext,ixyz0=in%plot_pot_axes)
        if (iproc == 0) call yaml_map('Writing local potential in file','local_potential'//gridformat)
        call plot_density(iproc,nproc,trim(in%dir_output)//'local_potential' // gridformat,&
             atoms,rxyz,denspot%pkernel,in%nspin,denspot%rhov,ixyz0=in%plot_pot_axes)
     else if (any(in%plot_pot_axes>=0)) then
        call f_err_throw('The coordinates of the point through which '//&
             'the potential shall be plotted must all be non-zero')
     else
        if (iproc == 0) call yaml_map('Writing external potential in file', 'external_potential'//gridformat)
        call plot_density(iproc,nproc,trim(in%dir_output)//'external_potential' // gridformat,&
             atoms,rxyz,denspot%pkernel,1,denspot%V_ext)
        if (iproc == 0) call yaml_map('Writing local potential in file','local_potential'//gridformat)
        call plot_density(iproc,nproc,trim(in%dir_output)//'local_potential' // gridformat,&
             atoms,rxyz,denspot%pkernel,in%nspin,denspot%rhov)
     end if
  end if

  call f_free_ptr(denspot%V_ext)

  !variables substitution for the PSolver part
  n1i=KSwfn%Lzd%Glr%d%n1i
  n2i=KSwfn%Lzd%Glr%d%n2i
  n3i=KSwfn%Lzd%Glr%d%n3i

  if (f_int(inputpsi) /= INPUT_PSI_EMPTY) then
     !------------------------------------------------------------------------
     ! here we start the calculation of the forces
     if (iproc == 0) then
        call yaml_comment('Forces Calculation',hfill='-')
     end if

     !refill projectors for tails, davidson
     refill_proj=((in%rbuf > 0.0_gp) .or. DoDavidson .or. in%sdos) .and. DoLastRunThings
     if (inputpsi .hasattr. 'CUBIC') then
        fpulay = f_malloc0_ptr((/ 3, atoms%astruct%nat /),id='fpulay')
     end if

     if (DoLastRunThings) then
        if (inputpsi .hasattr. 'CUBIC') then
           calculate_dipole=.true.
           calculate_quadrupole=.true.
        else
           calculate_dipole = in%lin%calc_dipole
           calculate_quadrupole=in%lin%calc_quadrupole
        end if
        output_denspot = in%output_denspot
     else
        output_denspot = ENUM_EMPTY
        calculate_dipole = .false.
        calculate_quadrupole = .false.
     end if

     call kswfn_post_treatments(iproc, nproc, KSwfn, tmb, &
          inputpsi .hasattr. 'LINEAR',&
          fxyz, fnoise, fion, fdisp, fpulay, &
          strten, pressure, ewaldstr, xcstr, GPU, denspot, atoms, rxyz, nlpsp, &
          output_denspot, in%dir_output, gridformat, refill_proj, calculate_dipole, calculate_quadrupole, &
          in%calculate_strten,in%nspin, in%plot_pot_axes)

     call f_free_ptr(fpulay)
  end if

  call f_free_ptr(fion)
  call f_free_ptr(fdisp)

  !if (nvirt > 0 .and. in%inputPsiId == 0) then
  if (DoDavidson .and. (inputpsi .hasattr. 'CUBIC')) then

     !for a band structure calculation allocate the array in which to put the eigenvalues
     if (associated(in%kptv) .and. in%nkptv > 0) then
        band_structure_eval = f_malloc((/ KSwfn%orbs%norbu+KSwfn%orbs%norbd+in%nspin*norbv, in%nkptv /),id='band_structure_eval')
     end if

     !calculate Davidson procedure for all the groups of k-points which are chosen
     ikpt=1
     do igroup=1,in%ngroups_kptv

        ! Set-up number of states and shifting values.
        nvirtu = norbv
        nvirtd = 0
        if (in%nspin==2) nvirtd=nvirtu
        ! Create the orbitals.
        if (associated(in%kptv) .and. in%nkptv > 0) then
           nvirtu = nvirtu + KSwfn%orbs%norbu
           nvirtd = nvirtd + KSwfn%orbs%norbd
           nvirt  = nvirtu+nvirtd

           !number of k-points for this group
           nkptv = in%nkptsv_group(igroup) !size(in%kptv, 2)

           wkptv = f_malloc(nkptv,id='wkptv')
           wkptv(:) = real(1.0, gp) / real(nkptv, gp)
           call orbitals_descriptors(iproc,nproc,nvirtu+nvirtd,nvirtu,nvirtd, &
                KSwfn%orbs%nspin,KSwfn%orbs%nspinor,nkptv, &
                in%kptv(:,sum(in%nkptsv_group(1:igroup - 1)) + 1:sum(in%nkptsv_group(1:igroup))), &
                wkptv,VTwfn%orbs,LINEAR_PARTITION_NONE)
           !allocate communications arrays for virtual orbitals
           call orbitals_communicators(iproc,nproc,KSwfn%Lzd%Glr,VTwfn%orbs,VTwfn%comms)

           call f_free(wkptv)

           !free projectors
           call free_DFT_PSP_projectors(nlpsp)

           ! Calculate all projectors, or allocate array for on-the-fly calculation
           call timing(iproc,'CrtProjectors ','ON')
           call orbital_basis_associate(ob,orbs=VTwfn%orbs,Lzd=KSwfn%Lzd)
           call createProjectorsArrays(iproc,nproc,KSwfn%Lzd%Glr,rxyz,atoms,ob,&
                in%frmult,in%frmult,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
                in%projection,.false.,nlpsp,.true.)
           call orbital_basis_release(ob)
           call timing(iproc,'CrtProjectors ','OF')
           if (iproc == 0) call print_nlpsp(nlpsp)

           nvirt = KSwfn%orbs%norbu + in%nvirt
        else
           !the virtual orbitals should be in agreement with the traditional k-points
           call orbitals_descriptors(iproc,nproc,nvirtu+nvirtd,nvirtu,nvirtd, &
                KSwfn%orbs%nspin,KSwfn%orbs%nspinor,KSwfn%orbs%nkpts,&
                KSwfn%orbs%kpts,KSwfn%orbs%kwgts,VTwfn%orbs,LINEAR_PARTITION_NONE,&
                basedist=KSwfn%orbs%norb_par(0:,1:),basedistu=KSwfn%orbs%norbu_par(0:,1:),&
                basedistd=KSwfn%orbs%norbd_par(0:,1:))
           !allocate communications arrays for virtual orbitals
           call orbitals_communicators(iproc,nproc,KSwfn%Lzd%Glr,VTwfn%orbs,VTwfn%comms,&
                basedist=KSwfn%comms%nvctr_par(0:,1:))

           nvirt = .if. (in%nvirt > 0) .then. in%nvirt .else. norbv
        end if

        !allocate psivirt pointer (note the orbs dimension)
        VTwfn%psi = f_malloc_ptr(max(VTwfn%orbs%npsidim_comp, VTwfn%orbs%npsidim_orbs),id='VTwfn%psi')
        !to avoid problems with the bindings
        VTwfn%c_obj=0
        !no paw for the Virtual Wavefunction (and this for a while)
        call nullify_paw_objects(VTwfn%paw)
        !define Local zone descriptors, should pass from initializers
        VTwfn%Lzd = KSwfn%Lzd
        VTwfn%orthpar=KSwfn%orthpar
        VTwfn%SIC=SIC_data_null() !then fill it if needed
        allocate(VTwfn%confdatarr(VTwfn%orbs%norbp))
        call default_confinement_data(VTwfn%confdatarr,VTwfn%orbs%norbp)

        if (in%norbv < 0) then
           call direct_minimization(iproc,nproc,in,atoms,&
                nvirt,rxyz,denspot%rhov,nlpsp, &
                denspot%pkernelseq,denspot%dpbox,denspot%xc,GPU,KSwfn,VTwfn)

           if(abs(in%nplot)>KSwfn%orbs%norb+nvirt) then
              if(iproc==0) call yaml_warning('More plots requested than orbitals calculated')
           end if
        else if (in%norbv > 0 .or. in%nkptv > 0) then
           call davidson(iproc,nproc,in,atoms,&
                KSwfn%orbs,VTwfn%orbs,nvirt,VTwfn%Lzd,&
                KSwfn%comms,VTwfn%comms,&
                rxyz,denspot%rhov,nlpsp, &
                denspot%pkernelseq,KSwfn%psi,VTwfn%psi,denspot%dpbox,denspot%xc,GPU)
!!$           call constrained_davidson(iproc,nproc,in,atoms,&
!!$                orbs,orbsv,in%nvirt,Lzd%Glr,comms,VTwfn%comms,&
!!$                hx,hy,hz,rxyz,denspot%rhov,nlpsp, &
!!$                psi,VTwfn%psi,nscatterarr,ngatherarr,GPU)
           if(abs(in%nplot)>KSwfn%orbs%norb+in%nvirt) then
              if(iproc==0) call yaml_warning('More plots requested than orbitals calculated')
           end if
        end if
        if(in%output_wf /= ENUM_EMPTY .and. abs(in%norbv)>0 ) then
           call dump_eigenfunctions(trim(in%dir_output),in%nplot,atoms,VTwfn%Lzd%hgrids,VTwfn%Lzd%Glr,&
                KSwfn%orbs,VTwfn%orbs,rxyz,KSwfn%psi,VTwfn%psi)
        end if

        deallocate(VTwfn%confdatarr)

        ! Write virtual wavefunctions in ETSF format
        if (in%output_wf /= ENUM_EMPTY  .and. abs(in%norbv) > 0) then
           call  writemywaves(iproc,trim(in%dir_output) // "virtuals", f_int(in%output_wf), &
                VTwfn%orbs,n1,n2,n3,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
                atoms,rxyz,KSwfn%Lzd%Glr%wfd,VTwfn%psi)
        end if

        if (trim(in%dir_perturbation)/='none/') then !input variable for the Completeness relation expression
           if (io_files_exists(in%dir_perturbation,'wavefunction',&
                Kswfn%orbs)) then
              rxyz_tmp=f_malloc([3,atoms%astruct%nat],id='rxyz_old')
              psi_perturbed=f_malloc(max(KSwfn%orbs%npsidim_orbs,KSwfn%orbs%npsidim_comp),id='psi_perturbed')
              hpsi_perturbed=f_malloc(max(KSwfn%orbs%npsidim_orbs,KSwfn%orbs%npsidim_comp),id='hpsi_perturbed')
              h2psi_perturbed=f_malloc(max(KSwfn%orbs%npsidim_orbs,KSwfn%orbs%npsidim_comp),id='hpsi_perturbed')
              call readmywaves(iproc,trim(in%dir_perturbation)// "wavefunction",input_wf_format,KSwfn%orbs,&
                   n1,n2,n3,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),atoms,rxyz_tmp,rxyz,  &
                   KSwfn%Lzd%Glr%wfd,psi_perturbed)

              !apply the hamiltonian to the perturbed wavefunctions
              !allocate the potential in the full box
              call full_local_potential(iproc,nproc,KSwfn%orbs,KSwfn%Lzd,0,&
                   denspot%dpbox,denspot%xc,denspot%rhov,denspot%pot_work)
              call FullHamiltonianApplication(iproc,nproc,atoms,KSwfn%orbs,&
                   KSwfn%Lzd,nlpsp,KSwfn%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,&
                   psi_perturbed,hpsi_perturbed,&
                   KSwfn%PAW,energs_fake,in%SIC,GPU,denspot%xc,&
                   denspot%pkernelseq)
              call FullHamiltonianApplication(iproc,nproc,atoms,KSwfn%orbs,&
                   KSwfn%Lzd,nlpsp,KSwfn%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,&
                   hpsi_perturbed,h2psi_perturbed,&
                   KSwfn%PAW,energs_fake,in%SIC,GPU,denspot%xc,&
                   denspot%pkernelseq)

              call orbital_basis_associate(ob_occ,orbs=KSwfn%orbs,&
                   Lzd=KSwfn%Lzd,phis_wvl=KSwfn%psi,comms=KSwfn%comms)
              call orbital_basis_associate(ob_virt,orbs=VTwfn%orbs,&
                   Lzd=VTwfn%Lzd,phis_wvl=VTwfn%psi,comms=VTwfn%comms)
              call orbital_basis_associate(ob_prime,orbs=KSwfn%orbs,&
                   Lzd=KSwfn%Lzd,phis_wvl=psi_perturbed,comms=KSwfn%comms)

              call evaluate_completeness_relation(ob_occ,ob_virt,ob_prime,hpsi_perturbed,h2psi_perturbed)

              call f_free(h2psi_perturbed)
              call f_free(hpsi_perturbed)
              call f_free(psi_perturbed)
              call f_free(rxyz_tmp)
              if (nproc > 1) then
                 call f_free_ptr(denspot%pot_work)
              else
                 nullify(denspot%pot_work)
              end if
              call orbital_basis_release(ob_prime)
              call orbital_basis_release(ob_virt)
              call orbital_basis_release(ob_occ)
           else
              if (bigdft_mpi%iproc ==0) &
                   call yaml_warning('Files for perturbation expansion "'&
                   +in%dir_perturbation+'" do not exist')
           end if
        end if

        !start the Casida's treatment
        !if (in%tddft_approach=='TDA') then
        
        if (in%tddft_approach .ne. 'none') then

           !does it make sense to use GPU only for a one-shot sumrho?
           if (GPU%OCLconv) then
              call allocate_data_OCL(KSwfn%Lzd%Glr%d%n1,KSwfn%Lzd%Glr%d%n2,KSwfn%Lzd%Glr%d%n3,&
                   KSwfn%Lzd%Glr%mesh_coarse,&
                   in%nspin,KSwfn%Lzd%Glr%wfd,KSwfn%orbs,GPU)
           end if

           !this could have been calculated before
           ! Potential from electronic charge density
           !WARNING: this is good just because the TDDFT is done with LDA
           call sumrho(denspot%dpbox,KSwfn%orbs,KSwfn%Lzd,&
                GPU,atoms%astruct%sym,denspot%rhod,denspot%xc,KSwfn%psi,denspot%rho_psi)
           call communicate_density(denspot%dpbox,KSwfn%orbs%nspin,&
                denspot%rhod,denspot%rho_psi,denspot%rhov,.false.)
           call denspot_set_rhov_status(denspot, ELECTRONIC_DENSITY, -1,iproc,nproc)

           if (GPU%OCLconv) then
              call free_gpu_OCL(GPU,KSwfn%orbs,in%nspin)
           end if

           !Allocate second Exc derivative
           if (denspot%dpbox%n3p >0) then
              !initialize fxc to zero in the case of HF calculation
              denspot%f_XC = f_malloc0_ptr((/ n1i , n2i , denspot%dpbox%n3p , in%nspin+1 /),id='denspot%f_XC')
           else
              denspot%f_XC = f_malloc_ptr((/ 1 , 1 , 1 , in%nspin+1 /),id='denspot%f_XC')
           end if

!!$$!MM test
!!$$ if (iproc==1) open(unit=201,file='density.dat')
!!$$ n1m=KSwfn%Lzd%Glr%d%n1i
!!$$ n2m=KSwfn%Lzd%Glr%d%n2i
!!$$ n3m=denspot%dpbox%n3p
!!$$ i1=n1m/2
!!$$ i2=n2m/2
!!$$ if (in%nspin==1) then
!!$$    do i3p=1,n3m
!!$$       if (iproc==1) write(201,'(E16.9E2,2x)') denspot%rhov(i1+((i2-1)+(i3p-1)*n2m)*n1m)
!!$$    end do
!!$$ else
!!$$    do i3p=1,n3m
!!$$       if (iproc==1) write(201,'(2(E16.9E2,2x))') denspot%rhov(i1+((i2-1)+(i3p-1)*n2m)*n1m),&
!!$$                             denspot%rhov(i1+((i2-1)+(i3p-1)*n2m)*n1m+n1m*n2m*n3m)
!!$$    end do
!!$$ end if
!!$$ if (iproc==1) close(201)
!!$$!MM test END


           call XC_potential(atoms%astruct%geocode,'D',iproc,nproc,bigdft_mpi%mpi_comm,&
                KSwfn%Lzd%Glr%d%n1i,KSwfn%Lzd%Glr%d%n2i,KSwfn%Lzd%Glr%d%n3i,denspot%xc,&
                denspot%dpbox%mesh%hgrids,&
                denspot%rhov,energs%exc,energs%evxc,in%nspin,denspot%rho_C,&
                denspot%rhohat,denspot%V_XC,xcstr,denspot%f_XC)

!!$$!MM test 
!!$$           if (iproc==0) open(unit=201,file='toplot0.dat')
!!$$           if (iproc==2) open(unit=203,file='toplot2.dat')
!!$$           n1m=KSwfn%Lzd%Glr%d%n1i
!!$$           n2m=KSwfn%Lzd%Glr%d%n2i
!!$$           n3m=denspot%dpbox%n3p
!!$$           i1=n1m/2
!!$$           i2=n2m/2
!!$$
!!$$           if (in%nspin==1) then
!!$$
!!$$              do i3p=1,n3m
!!$$                 if (iproc==0) write(201,'(i4,4(2x,E16.9E2))') i3p,&
!!$$                                  denspot%rhov(i1+((i2-1)+(i3p-1)*n2m)*n1m),&
!!$$                                  denspot%V_XC(i1,i2,i3p,1), &
!!$$                                  denspot%f_XC(i1,i2,i3p,1), denspot%f_XC(i1,i2,i3p,2)
!!$$                 if (iproc==2) write(203,'(i4,4(2x,E16.9E2))') i3p,&
!!$$                                  denspot%rhov(i1+((i2-1)+(i3p-1)*n2m)*n1m),&
!!$$                                  denspot%V_XC(i1,i2,i3p,1), &
!!$$                                  denspot%f_XC(i1,i2,i3p,1), denspot%f_XC(i1,i2,i3p,2)
!!$$              end do
!!$$
!!$$           else if (in%nspin==2) then
!!$$
!!$$              do i3p=1,n3m
!!$$                 if (iproc==0) write(201,'(i4,7(2x,E16.9E2))') i3p,&
!!$$                                  denspot%rhov(i1+((i2-1)+(i3p-1)*n2m)*n1m),& 
!!$$                                  denspot%rhov(i1+((i2-1)+(i3p-1)*n2m)*n1m+n1m*n2m*n3m),&
!!$$                                  denspot%V_XC(i1,i2,i3p,1), denspot%V_XC(i1,i2,i3p,2), &
!!$$                                  denspot%f_XC(i1,i2,i3p,1), denspot%f_XC(i1,i2,i3p,2), denspot%f_XC(i1,i2,i3p,3)
!!$$                 if (iproc==2) write(203,'(i4,7(2x,E16.9E2))') i3p,&
!!$$                                  denspot%rhov(i1+((i2-1)+(i3p-1)*n2m)*n1m),& 
!!$$                                  denspot%rhov(i1+((i2-1)+(i3p-1)*n2m)*n1m+n1m*n2m*n3m),&
!!$$                                  denspot%V_XC(i1,i2,i3p,1), denspot%V_XC(i1,i2,i3p,2), &
!!$$                                  denspot%f_XC(i1,i2,i3p,1), denspot%f_XC(i1,i2,i3p,2), denspot%f_XC(i1,i2,i3p,3)
!!$$              end do
!!$$
!!$$           end if
!!$$
!!$$           if (iproc==0) close(unit=201)
!!$$           if (iproc==2) close(unit=203)
!!$$!MM end test 


           call denspot_set_rhov_status(denspot, CHARGE_DENSITY, -1,iproc,nproc)

           !select the active space if needed
           exc_fac = 0.0
           if (in%ixc==100) then 
              !call yaml_comment('HF detected',hfill='-')
              exc_fac=1.0 !if HF, exc_fac=1.0 
           end if

           call tddft_casida(iproc,nproc,trim(in%dir_output),atoms,rxyz,&
                denspot%dpbox%n3p,denspot%dpbox%ngatherarr(0,1),&
                KSwfn%Lzd%Glr,in%tddft_approach,KSwfn%orbs,VTwfn%orbs,&
                denspot%f_XC,denspot%pkernelseq,KSwfn%psi,VTwfn%psi,&
                exc_fac,denspot%dpbox%bitp)

           call f_free_ptr(denspot%f_XC)

        end if

        call deallocate_comms(VTwfn%comms)
        call deallocate_orbs(VTwfn%orbs)

        !in the case of band structure calculation, copy the values of the eigenvectors
        !into a new array to write them afterwards
        if (associated(in%kptv) .and. in%nkptv > 0) then
           call vcopy(VTwfn%orbs%norb*nkptv,VTwfn%orbs%eval(1),1,band_structure_eval(1,ikpt),1)
           !increment the value of ikpt
           ikpt=ikpt+in%nkptsv_group(igroup)
        end if

        call f_free_ptr(VTwfn%orbs%eval)

        !if the local analysis has to be performed the deallocation should not be done
        call f_free_ptr(VTwfn%psi)

     end do

     if (associated(in%kptv) .and. in%nkptv > 0) then
        !dump the band structure eigenvalue on a file and deallocate it
        if (iproc == 0) then
           open(unit=11,file='band_structure.dat',status='unknown')
           do ikpt=1,in%nkptv
              write(11,'(i5,3(f12.6),10000(1pe12.4))')ikpt,&
                   (in%kptv(i,ikpt),i=1,3),(band_structure_eval(i,ikpt),i=1,VTwfn%orbs%norb)
           end do
           !tentative gnuplot string for the band structure file
           write(11,'(a,9999(a,i6,a))')&
                "#plot 'band_structure.dat' u 1:5 w l t ''",&
                (",'' u 1:",5+i-1," w l t ''" ,i=2,VTwfn%orbs%norb)
           close(unit=11)
        end if
        call f_free(band_structure_eval)
     end if

  else if (in%nplot > 0 .and. DoLastRunThings) then
     !plot the occupied wavefunctions for visualization purposes
     call dump_eigenfunctions(trim(in%dir_output),-in%nplot,atoms,KSwfn%Lzd%hgrids,KSwfn%Lzd%Glr,&
          KSwfn%orbs,KSwfn%orbs,rxyz,KSwfn%psi,KSwfn%psi)
  end if

  if (DoLastRunThings) then
     !perform here the mulliken charge and density of states
     !localise them on the basis of gatom of a number of atoms
     if ((inputpsi .hasattr. 'GAUSSIAN') .and. (inputpsi .hasattr. 'CUBIC')) then
        !here one must check if psivirt should have been kept allocated
        if (.not. DoDavidson) then
           VTwfn%orbs%norb=0
           VTwfn%orbs%norbp=0
        end if
        call local_analysis(iproc,nproc,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
             atoms,rxyz,KSwfn%Lzd%Glr,KSwfn%orbs,VTwfn%orbs,KSwfn%psi,VTwfn%psi)
     else if (optLoop%itrpmax /= 1 .and. get_verbose_level() >= 2) then
        ! Do a full DOS calculation.
        if (iproc == 0) call global_analysis(KSwfn%orbs, in%Tel,in%occopt,trim(in%dir_output) // "dos.gnuplot")
     end if
     if (in%sdos) then !spatially resolved DOS
        if (iproc==0) call yaml_comment('Calculating Spatially Resolved DOS',hfill='-')
!!$        !apply the hamiltonian to the perturbed wavefunctions
!!$        hpsi_tmp=f_malloc(max(KSwfn%orbs%npsidim_orbs,KSwfn%orbs%npsidim_comp),id='hpsi_tmp')
!!$        !allocate the potential in the full box
!!$        call full_local_potential(iproc,nproc,KSwfn%orbs,KSwfn%Lzd,0,&
!!$             denspot%dpbox,denspot%xc,denspot%rhov,denspot%pot_work)
!!$        call FullHamiltonianApplication(iproc,nproc,atoms,KSwfn%orbs,&
!!$             KSwfn%Lzd,nlpsp,KSwfn%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,&
!!$             KSwfn%psi,hpsi_tmp,&
!!$             KSwfn%PAW,energs_fake,in%SIC,GPU,denspot%xc,&
!!$             denspot%pkernelseq)
        call orbital_basis_associate(ob_occ,orbs=KSwfn%orbs,&
             Lzd=KSwfn%Lzd,phis_wvl=KSwfn%psi,comms=KSwfn%comms)
        call spatially_resolved_dos(ob_occ,trim(in%dir_output))
        call orbital_basis_release(ob_occ)
        if (nproc > 1) then
           call f_free_ptr(denspot%pot_work)
        else
           nullify(denspot%pot_work)
        end if       
!!$        call f_free(hpsi_tmp)
     end if
  end if
  if (((in%exctxpar == 'OP2P' .and. xc_exctXfac(denspot%xc) /= 0.0_gp) &
       .or. in%SIC%alpha /= 0.0_gp) .and. nproc >1) then
     if (.not. associated(denspot%pkernelseq%kernel,target=denspot%pkernel%kernel) .and. &
          associated(denspot%pkernelseq%kernel)) then
        call pkernel_free(denspot%pkernelseq)
     end if
  else if (nproc == 1 .and. (in%exctxpar == 'OP2P' .or. in%SIC%alpha /= 0.0_gp)) then
     nullify(denspot%pkernelseq%kernel)
  end if

  call deallocate_paw_objects(KSwfn%paw)

  !------------------------------------------------------------------------
  if ((in%rbuf > 0.0_gp) .and. atoms%astruct%geocode == 'F' .and. DoLastRunThings ) then
     if (in%SIC%alpha /= 0.0_gp) call f_err_throw('Tail correction not admitted with SIC corrections for the moment',&
             err_name='BIGDFT_RUNTIME_ERROR')

     call timing(iproc,'Tail          ','ON')
     !    Calculate energy correction due to finite size effects
     !    ---reformat potential
     !!denspot%pot_work = f_malloc_ptr(n1i*n2i*n3i*in%nspin,id='denspot%pot_work')
     call full_local_potential(iproc,nproc,KSwfn%orbs,KSwfn%Lzd,0,&
          & denspot%dpbox,denspot%xc,denspot%rhov,denspot%pot_work)
     call PS_gather(denspot%rhov,denspot%pkernel,dest=denspot%pot_work,&
          nsrc=in%nspin)

     call pkernel_free(denspot%pkernel)

     call dpbox_free(denspot%dpbox)
     call f_free_ptr(denspot%V_XC)
     if (nproc > 1) call f_free_ptr(denspot%rhov)

     call CalculateTailCorrection(iproc,nproc,atoms,in%rbuf,KSwfn%orbs,&
          KSwfn%Lzd%Glr,nlpsp,in%ncongt,denspot%pot_work,KSwfn%Lzd%hgrids,&
          rxyz,in%crmult,in%frmult,in%nspin,&
          KSwfn%psi,(in%output_denspot /= ENUM_EMPTY),energs%ekin,energs%epot,energs%eproj,KSwfn%paw)

     call f_free_ptr(denspot%pot_work)

     energs%ebs=energs%ekin+energs%epot+energs%eproj
     energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%evsic+energs%eion+energs%edisp-energs%eTS+energs%ePV

     if (iproc == 0) then
        call yaml_mapping_open('Corrected Energies', flow=.true.)
        call yaml_map('Ekin', energs%ekin, fmt='(1pe18.11)')
        call yaml_map('Epot', energs%epot, fmt='(1pe18.11)')
        call yaml_map('Eproj',energs%eproj,fmt='(1pe18.11)')
        call yaml_mapping_close()
        call yaml_map('Total energy with tail correction',energy,fmt='(1pe24.17)')
        call yaml_mapping_close()
     endif

     call timing(iproc,'Tail          ','OF')
  else
     call pkernel_free(denspot%pkernel)
     !    No tail calculation
     if (nproc > 1) call fmpi_barrier(comm=bigdft_mpi%mpi_comm)!,ierr)
     call f_free_ptr(denspot%rhov)
     call f_free_ptr(denspot%V_XC)
     call dpbox_free(denspot%dpbox)
  endif
  ! --- End if of tail calculation

  !?!   !Finally, we add the entropic contribution to the energy from non-integer occnums
  !?!   if(orbs%eTS>0_gp) then
  !?!      energy=energy - orbs%eTS
  !?!
  !?!      if (iproc == 0) then
  !?!         write( *,'(1x,a,1(1x,1pe18.11))')&
  !?!              '  Entropic correction due to electronic temperature',orbs%eTS
  !?!         write( *,'(1x,a,1x,1pe24.17)')&
  !?!              'Free energy (= total energy - T*S)  ',energy
  !?!      endif
  !?!    endif

  call deallocate_before_exiting

!START debug code added by bastian
!write(*,*)"BIGDFTbastian debug exit sub. cluster, iproc",iproc
!call f_utils_flush(6)
!END debug code added by bastian

contains

  !> Routine which deallocate the pointers and the arrays before exiting
  subroutine deallocate_before_exiting
    use communications_base, only: deallocate_comms
    use module_cfd, only: cfd_free
    implicit none
    external :: gather_timings
  !when this condition is verified we are in the middle of the SCF cycle
    if (infocode /=0 .and. infocode /=1 .and. inputpsi /= 'INPUT_PSI_EMPTY') then
       call f_free_ptr(denspot%V_ext)

       if (((in%exctxpar == 'OP2P' .and. xc_exctXfac(denspot%xc) /= 0.0_gp) &
            .or. in%SIC%alpha /= 0.0_gp) .and. nproc >1) then
          if (.not. associated(denspot%pkernelseq%kernel,target=denspot%pkernel%kernel) .and. &
               associated(denspot%pkernelseq%kernel)) then
             call pkernel_free(denspot%pkernelseq)
          end if
       else if (nproc == 1 .and. (in%exctxpar == 'OP2P' .or. in%SIC%alpha /= 0.0_gp)) then
          nullify(denspot%pkernelseq%kernel)
       end if
       call pkernel_free(denspot%pkernel)

       ! calc_tail false
       call f_free_ptr(denspot%rhov)
       call f_free_ptr(denspot%V_XC)

       call dpbox_free(denspot%dpbox)
       call f_free_ptr(fion)
       call f_free_ptr(fdisp)
    end if
    call xc_end(denspot%xc)
    !@todo fix a method for freeing the denspot structure
    call f_free_ptr(denspot%rho_ion)
    !@todo Ugly here to be moved elsewhere.
    if (associated(atoms%pawtab)) then
       call libxc_functionals_end()
    end if

    !free GPU if it is the case
    if (GPU%OCLconv .and. .not.(DoDavidson)) then
       call free_gpu_OCL(GPU,KSwfn%orbs,in%nspin)
    end if

    ! Free all remaining parts of denspot
    call deallocate_rho_descriptors(denspot%rhod)
    if(associated(denspot%rho_C)) then
       call f_free_ptr(denspot%rho_C)
    end if
    call f_free_ptr(denspot%rhohat)
    call cfd_free(denspot%cfd)
       
    call f_free(denspot0)

    ! Free all remaining parts of KSwfn
    call deallocate_Lzd_except_Glr(KSwfn%Lzd)

    call deallocate_comms(KSwfn%comms)
    call deallocate_orbs(KSwfn%orbs)
    if (inputpsi .hasattr. 'CUBIC') deallocate(KSwfn%confdatarr)

    ! Free projectors.
    call free_DFT_PSP_projectors(nlpsp)

    ! Stop signals
    if (in%signaling .and. iproc == 0) then
       call bigdft_signals_rm_denspot(in%gmainloop)
       call bigdft_signals_rm_energs(in%gmainloop)
       call bigdft_signals_rm_wf(in%gmainloop)
       call bigdft_signals_rm_optloop(in%gmainloop)
       call localfields_free_wrapper(denspot%c_obj)
       !call energs_free_wrapper(energs%c_obj)
       call optloop_free_wrapper(optLoop%c_obj)
       call wf_free_wrapper(KSwfn%c_obj)
       call wf_free_wrapper(tmb%c_obj)
    end if

     if (iproc == 0 .and. (in%inputPsiId .hasattr. 'CUBIC') .and. infocode==1) then
        call yaml_warning('Self-consistent cycle did not meet convergence criteria')
     end if
    !release the yaml document
    call yaml_release_document()

    call f_release_routine()

    !end of wavefunction minimisation
    call timing(bigdft_mpi%mpi_comm,'LAST','PR')
    call dict_init(dict_timing_info)
    if (DoLastRunThings) call f_malloc_dump_status(dict_summary=dict_timing_info)
    call mpi_environment_dict(bigdft_mpi,dict_timing_info)
    !call build_dict_info(dict_timing_info)
    call f_timing_stop(mpi_comm=bigdft_mpi%mpi_comm,nproc=bigdft_mpi%nproc,&
         gather_routine=gather_timings,dict_info=dict_timing_info)
    call dict_free(dict_timing_info)
    call cpu_time(tcpu1)
    call system_clock(ncount1,ncount_rate,ncount_max)
    tel=dble(ncount1-ncount0)/dble(ncount_rate)
    if (iproc == 0) then
       call yaml_comment('Timing for root process',hfill='-')
       call yaml_mapping_open('Timings for root process')
       call yaml_map('CPU time (s)',tcpu1-tcpu0,fmt='(f12.2)')
       call yaml_map('Elapsed time (s)',tel,fmt='(f12.2)')
       call yaml_mapping_close()
       call yaml_flush_document()
    end if

  END SUBROUTINE deallocate_before_exiting

END SUBROUTINE cluster


!> Kohn-Sham wavefunction optimization loop
subroutine kswfn_optimization_loop(iproc, nproc, opt, &
     & alphamix, idsx, inputpsi, KSwfn, denspot, nlpsp, energs, atoms, GPU, xcstr, &
     & in)
  use module_base
  use module_types
  use module_interfaces, only: denspot_set_history, hpsitopsi, last_orthon
  use io, only: write_energies
  use module_xc, only: XC_NO_HARTREE
  use yaml_output
  use public_enums
  implicit none
  real(dp), dimension(6), intent(out) :: xcstr
  integer, intent(in) :: iproc, nproc, idsx
  type(f_enumerator), intent(in) :: inputpsi
  real(gp), intent(in) :: alphamix
  type(DFT_optimization_loop), intent(inout) :: opt
  type(DFT_wavefunction), intent(inout) :: KSwfn
  type(DFT_local_fields), intent(inout) :: denspot
  type(energy_terms), intent(inout) :: energs
  type(atoms_data), intent(in) :: atoms
  type(GPU_pointers), intent(inout) :: GPU
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(input_variables), intent(in) :: in !< @todo Remove me

  character(len = *), parameter :: subname = "kswfn_optimization_loop"
  logical :: endloop, scpot, endlooprp, lcs
  integer :: ndiis_sd_sw, idsx_actual_before, linflag, ierr,iter_for_diis
  integer :: ikpt_homo,ikpt_lumo,ispin_homo,ispin_lumo
  real(gp) :: gnrm_zero,homo,lumo,occup_lumo,minres_gpe
  character(len=5) :: final_out
  !temporary variables for PAPI computation
  ! real(kind=4) :: rtime, ptime,  mflops
  ! integer(kind=8) ::flpops
  call f_routine(id=subname)
!  !start PAPI counting
!  if (iproc==0) call PAPIF_flops(rtime, ptime, flpops, mflops,ierr)

  ! Setup the mixing, if necessary
  call denspot_set_history(denspot,opt%scf)

  ! allocate arrays necessary for DIIS convergence acceleration
  call allocate_diis_objects(idsx,in%alphadiis,sum(KSwfn%comms%ncntt(0:nproc-1)),&
       KSwfn%orbs%nkptsp,KSwfn%orbs%nspinor,KSwfn%diis)

  !number of switching betweed DIIS and SD during self-consistent loop
  ndiis_sd_sw=0
  !previous value of idsx_actual to control if switching has appeared
  idsx_actual_before=KSwfn%diis%idsx

  gnrm_zero=0.0d0
  opt%gnrm=1.d10
  opt%rpnrm=1.d10
  endlooprp=.false.
  energs%e_prev=0.0_gp

  !normal opt%infocode, if everything go through smoothly we should keep this
  opt%infocode=0
  !yaml output
  if (iproc==0 .and. opt%itrpmax > 0 ) then
     call yaml_comment('Self-Consistent Cycle',hfill='-')
     call yaml_sequence_open('Ground State Optimization')
  end if
  opt%itrp=1
  rhopot_loop: do
     KSwfn%diis%energy_old=1.d100
     if (opt%itrp > opt%itrpmax) exit
     !yaml output
     if (iproc==0) then
        call yaml_sequence(advance='no')
        call yaml_sequence_open("Hamiltonian Optimization",label=&
             'itrp'//trim(adjustl(yaml_toa(opt%itrp,fmt='(i3.3)'))))

     end if
     !set the opt%infocode to the value it would have in the case of no convergence
     opt%infocode=1
     opt%itrep=1
     iter_for_diis=0 !initialize it here for keeping the history also after a subspace diagonalization
     subd_loop: do
        if (opt%itrep > opt%nrepmax) exit subd_loop
        !yaml output
        if (iproc==0) then
           call yaml_sequence(advance='no')
           call yaml_mapping_open("Subspace Optimization",label=&
                'itrep'//trim(adjustl(yaml_toa(opt%itrp,fmt='(i3.3)')))//'-'//&
                trim(adjustl(yaml_toa(opt%itrep,fmt='(i2.2)'))))
        end if

        !yaml output
        if (iproc==0) then
           call yaml_sequence_open("Wavefunctions Iterations")
        end if
        opt%iter=1
        iter_for_diis=0
        wfn_loop: do
           if (opt%iter > opt%itermax) exit wfn_loop

           !control whether the minimisation iterations should end after the hamiltionian application
           endloop= (opt%gnrm <= opt%gnrm_cv .or. opt%iter == opt%itermax) .and. opt%iter >= opt%itermin !Bastian

           if (iproc == 0) then
              !yaml output
              if (endloop .and. (opt%itrpmax==1 .or. opt%itrpmax >1 .and. endlooprp)) then
                 call yaml_sequence(label='FINAL'//trim(adjustl(yaml_toa(opt%itrep,fmt='(i3.3)'))),advance='no')
              else if (endloop .and. opt%itrep == opt%nrepmax) then
                 call yaml_sequence(label='final'//trim(adjustl(yaml_toa(opt%itrp,fmt='(i4.4)'))),&
                      advance='no')
              else
                 call yaml_sequence(advance='no')
              end if
              call yaml_mapping_open(flow=.true.)
              if (get_verbose_level() > 0) &
                   call yaml_comment('iter:'//yaml_toa(opt%iter,fmt='(i6)'),hfill='-')
           endif

           !control how many times the DIIS has switched into SD
           if (KSwfn%diis%idsx /= idsx_actual_before) ndiis_sd_sw=ndiis_sd_sw+1

           !let SD runs if the DIIS did not work the second time
           if (ndiis_sd_sw > 1) then
              KSwfn%diis%switchSD=.false.
           end if

           !stop the partial timing counter if necessary
           if (endloop .and. opt%itrpmax==1 .and. (opt%itrep == opt%nrepmax .or. opt%gnrm < opt%gnrm_cv)) &
                call timing(bigdft_mpi%mpi_comm,'WFN_OPT','PR')
           !logical flag for the self-consistent potential
           scpot=(( (opt%scf .hasattr. 'MIXING') .and. opt%iter==1 .and. opt%itrep==1) .or. & !mixing to be done
                (.not. (opt%scf .hasattr. 'MIXING'))) .and. & !direct minimisation
                .not. (denspot%xc%ixc == XC_NO_HARTREE) ! Need to calculate the scp pot (i.e. Hartree + XC)
           !allocate the potential in the full box
           !temporary, should change the use of flag in full_local_potential2
           linflag = 1
           if(in%linear == INPUT_IG_OFF) linflag = 0
           if(in%linear == INPUT_IG_TMO) linflag = 2

           !if (opt%iter == 1) minres_gpe=denspot%pkernel%minres
           !denspot%pkernel%minres=max(min(1.e-6_gp,opt%gnrm**2) ,minres_gpe)!!opt%gnrm_cv**2)

           nlpsp%apply_gamma_target=((opt%scf .hasattr. 'MIXING') .and. opt%itrp <= in%occupancy_control_itermax &
                .and. opt%itrep <= in%occupancy_control_nrepmax) .or. &
                (.not. (opt%scf .hasattr. 'MIXING') .and. opt%iter <= in%occupancy_control_itermax &
                .and. opt%itrep <= in%occupancy_control_nrepmax)

           !Calculates the application of the Hamiltonian on the wavefunction
           call psitohpsi(iproc,nproc,atoms,scpot,denspot,opt%itrp,opt%iter,opt%scf,alphamix,&
                nlpsp,linflag,in%unblock_comms,GPU,KSwfn,energs,opt%rpnrm,xcstr)

           endlooprp= (opt%itrp > 1 .and. opt%rpnrm <= opt%rpnrm_cv) .or. opt%itrp == opt%itrpmax

           call total_energies(energs, opt%iter, iproc)

           !check for convergence or whether max. numb. of iterations exceeded
           if (endloop) then
              if (opt%gnrm < opt%gnrm_cv) opt%infocode=0
              exit wfn_loop
           endif

           !evaluate the functional of the wavefunctions and put it into the diis structure
           !the energy values is printed out in this routine
           call calculate_energy_and_gradient(opt%iter,iproc,nproc,GPU,in%ncong,opt%scf,&
                energs,KSwfn,opt%gnrm,gnrm_zero)

           !control the previous value of idsx_actual
           idsx_actual_before=KSwfn%diis%idsx
           iter_for_diis=iter_for_diis+1
           call hpsitopsi(iproc,nproc,iter_for_diis,idsx,KSwfn,atoms,nlpsp,energs%eproj)

           if (inputpsi == 'INPUT_PSI_LCAO') then !might insert the policy here
              if ((opt%gnrm > 4.d0 .and. KSwfn%orbs%norbu /= KSwfn%orbs%norbd) .or. &
                   &   (KSwfn%orbs%norbu == KSwfn%orbs%norbd .and. opt%gnrm > 10.d0)) then
                 opt%infocode=3
              end if
           else if (inputpsi == 'INPUT_PSI_MEMORY_WVL') then
              if (opt%gnrm > 1.d0) then
                 opt%infocode=2
              end if
           end if
           !flush all writings on standard output
           if (iproc==0) then
              !yaml output
              call yaml_mapping_close() !iteration
              call yaml_flush_document()
           end if
           ! Emergency exit case
           if (opt%infocode == 2 .or. opt%infocode == 3) then
              if (nproc > 1) call fmpi_barrier(comm=bigdft_mpi%mpi_comm)
              if (iproc==0) then
                 !call yaml_mapping_close()
                 call yaml_sequence_close() !wfn iterations
                 call yaml_mapping_close()
                 call yaml_sequence_close() !itrep
                 if (opt%infocode==2) then
                    call yaml_warning('The norm of the residue is too large, need to recalculate input wavefunctions')
                 else if (opt%infocode ==3) then
                    call yaml_warning('The norm of the residue is too large also with input wavefunctions.')
                 end if
              end if
              exit rhopot_loop
           end if

           if (opt%c_obj /= 0) then
              call optloop_emit_iter(opt, OPTLOOP_WAVEFUNCTIONS, energs, iproc, nproc)
           end if

           opt%iter = opt%iter + 1
        !if (opt%iter == 2) stop
        end do wfn_loop


        if (opt%c_obj /= 0) then
           call optloop_emit_done(opt, OPTLOOP_WAVEFUNCTIONS, energs, iproc, nproc)
        end if

        if (iproc == 0) then
           !if (verbose > 1) write( *,'(1x,a,i0,a)')'done. ',opt%iter,' minimization iterations required'
           !write( *,'(1x,a)') &
           !     &   '--------------------------------------------------- End of Wavefunction Optimisation'
           if ((opt%itrpmax >1 .and. endlooprp) .or. opt%itrpmax == 1) then
              write(final_out, "(A5)") "FINAL"
           else
              write(final_out, "(A5)") "final"
           end if
           call write_energies(opt%iter,energs,opt%gnrm,gnrm_zero,final_out)
           call yaml_mapping_close()
           call yaml_flush_document()
           if (opt%itrpmax >1) then
              if ( KSwfn%diis%energy - KSwfn%diis%energy_min > 1.e-10_gp) &
                   call yaml_warning('Found an energy value lower than the ' // final_out // ' energy (delta=' // &
                   trim(yaml_toa(KSwfn%diis%energy-KSwfn%diis%energy_min,fmt='(1pe9.2)')) // ')')
           else
              !write this warning only if the system is closed shell
              call check_closed_shell(KSwfn%orbs,lcs)
              if (lcs) then
                 if ( energs%eKS - KSwfn%diis%energy_min > 1.e-10_gp) &
                      call yaml_warning('Found an energy value lower than the FINAL energy (delta=' // &
                      trim(yaml_toa(energs%eKS-KSwfn%diis%energy_min,fmt='(1pe9.2)')) // ')')
              end if
           end if
        end if

        if (iproc==0) then
           call yaml_sequence_close() !wfn iterations
           if (opt%iter == opt%itermax .and. opt%infocode/=0) &
                call yaml_warning('No convergence within the allowed number of minimization steps')
        end if

        call last_orthon(iproc,nproc,opt%iter,KSwfn,energs%evsum,.true.) !never deallocate psit and hpsi

!!$        !EXPERIMENTAL
!!$        !check if after convergence the integral equation associated with Helmholtz' Green function is satisfied
!!$        !note: valid only for negative-energy eigenstates
!!$        call integral_equation(iproc,nproc,atoms,KSwfn,denspot%dpbox%ngatherarr,denspot%rhov,GPU,proj,nlpspd,rxyz,KSwfn%paw)

        !exit if the opt%infocode is correct
        if (opt%infocode /= 0) then
           if(iproc==0) then
              if (opt%itrp == opt%itrpmax .and. opt%gnrm_cv > 0.0_gp) &
                   call yaml_warning('Wavefunctions not converged after cycle '// trim(yaml_toa(opt%itrep,fmt='(i0)')))
              if (opt%itrpmax > 1 .and. opt%itrp == opt%itrpmax .and. opt%gnrm > sqrt(opt%rpnrm)) &
                   call yaml_warning('Wavefunction residue is not consistent with density convergence (T_el too small?)')
              if (opt%itrep < opt%nrepmax) call yaml_comment('restart after diagonalisation')
              ! write(*,*) ' WARNING: Wavefunctions not converged after cycle',opt%itrep
              ! if (opt%itrep < opt%nrepmax) write(*,*)' restart after diagonalisation'
           end if
           opt%gnrm=1.d10

           if (opt%itrpmax == 1 .and. in%norbsempty > 0) then
              !recalculate orbitals occupation numbers
              call evaltoocc(iproc,nproc,.false.,in%Tel,KSwfn%orbs,in%occopt)

              !opt%gnrm =1.d10
              KSwfn%diis%energy_min=1.d10
              !KSwfn%diis%alpha=2.d0
              KSwfn%diis%alpha=in%alphadiis
           end if
        end if

        if (opt%itrpmax ==1) then
           call eigensystem_info(iproc,nproc,opt%gnrm,&
             KSwfn%Lzd%Glr%wfd%nvctr_c+7*KSwfn%Lzd%Glr%wfd%nvctr_f,&
             KSwfn%orbs,KSwfn%psi)
           if (opt%infocode /=0) then
              opt%gnrm =1.d10
           end if
        end if

        if (iproc==0) then
           call yaml_mapping_close()
           call yaml_flush_document()
        end if

        if (opt%infocode ==0) exit subd_loop

        if (opt%c_obj /= 0) then
           call optloop_emit_iter(opt, OPTLOOP_SUBSPACE, energs, iproc, nproc)
        end if

        opt%itrep = opt%itrep + 1
     end do subd_loop
     if (opt%c_obj /= 0) then
        call optloop_emit_done(opt, OPTLOOP_SUBSPACE, energs, iproc, nproc)
     end if

     if (iproc==0) then
        call yaml_sequence_close() !itrep
     end if


     if (opt%itrpmax > 1) then

        !recalculate orbitals occupation numbers
        call evaltoocc(iproc,nproc,.false.,in%Tel,KSwfn%orbs,in%occopt)

        call eigensystem_info(iproc,nproc,opt%gnrm,&
             KSwfn%Lzd%Glr%wfd%nvctr_c+7*KSwfn%Lzd%Glr%wfd%nvctr_f,&
             KSwfn%orbs,KSwfn%psi)

!!$        !this section has to be improved to determine whether the system has a semiconductor-like occupation
!!$        !or not
!!$        call orbs_get_gap(KSwfn%orbs,ikpt_homo,ikpt_lumo,ispin_homo,ispin_lumo,homo,lumo,occup_lumo)
!!$
!!$        if (iproc==0) then
!!$           call yaml_mapping_open('Highest Occupied Orbital',flow=.true.)
!!$           call yaml_map('e',homo,fmt='(1pe15.7)')
!!$           call yaml_map('kpt',ikpt_homo)
!!$           call yaml_map('s',ispin_homo)
!!$           call yaml_mapping_close()
!!$           call yaml_mapping_open('Lowest Unoccupied Orbital',flow=.true.)
!!$           call yaml_map('e',lumo,fmt='(1pe15.7)')
!!$           call yaml_map('kpt',ikpt_lumo)
!!$           call yaml_map('s',ispin_lumo)
!!$           call yaml_mapping_close()
!!$           call yaml_map('Gap (Ha, eV)',(lumo-homo)*[1.0_gp,Ha_eV],fmt='(1pe15.7)')
!!$        end if
!!$
!!$        !if the gap is too high for the electronic temperature inform that Direct minimization is available
!!$        if (occup_lumo < 1.d-6) then !lumo-homo > 5.0_gp*in%Tel) then
!!$           if (iproc==0) call yaml_warning('The gap for this system is very high, consider to perform direct minimization')
!!$        end if
        !stop the partial timing counter if necessary
        if (endlooprp) then
           call timing(bigdft_mpi%mpi_comm,'WFN_OPT','PR')
           exit rhopot_loop
        end if

        opt%gnrm =1.d10
        KSwfn%diis%energy_min=1.d10
        ! this line can be commented
        !KSwfn%diis%alpha=in%alphadiis
     end if

     if (iproc == 0) then
        !yaml output
        !summarize the key elements in the opt%itrp element
        if (opt%itrp >1) then
           call yaml_map('RhoPot Delta','*rpnrm'//trim(adjustl(yaml_toa(opt%itrp,fmt='(i4.4)'))))
           call yaml_map('Energies','*final'//trim(adjustl(yaml_toa(opt%itrp,fmt='(i4.4)'))))
!!$           call yaml_comment('End RhoPot Iterations, itrp:'//&
!!$                yaml_toa(opt%itrp,fmt='(i6)'))
        end if
     end if
     if (opt%c_obj /= 0) then
        call optloop_emit_iter(opt, OPTLOOP_HAMILTONIAN, energs, iproc, nproc)
     end if

     opt%itrp = opt%itrp + 1
  end do rhopot_loop

!!$  if (iproc ==0) then
!!$     call PAPIF_flops(rtime, ptime, flpops, mflops,ierr)
!!$
!!$     write (*,90) rtime, ptime, flpops, mflops
!!$
!!$90   format('           Real time (secs) :', f15.3, &
!!$          /'            CPU time (secs) :', f15.3,&
!!$          /'Floating point instructions :', i15,&
!!$          /'                     MFLOPS :', f15.3)
!!$
!!$
!!$  end if


  if (opt%c_obj /= 0) then
     call optloop_emit_done(opt, OPTLOOP_HAMILTONIAN, energs, iproc, nproc)
  end if
  if (iproc==0) call yaml_sequence_close() !opt%itrp
  !recuperate the information coming from the last iteration (useful for post-processing of the document)
  !only if everything got OK
  if (iproc==0 .and. opt%infocode == BIGDFT_SUCCESS .and. opt%itrpmax > 0) &
       call yaml_map('Last Iteration','*FINAL'//trim(adjustl(yaml_toa(opt%itrep,fmt='(i3.3)'))))

  !!do i_all=1,size(rhopot)
  !!    write(10000+iproc,*) rhopot(i_all)
  !!end do
  !!do i_all=1,size(psi)
  !!    write(11000+iproc,*) psi(i_all)
  !!end do
  !!do i_all=1,size(psi)
  !!    write(12000+iproc,*) psi(i_all)
  !!end do

  !this warning can be deplaced in write_energies
  if (inputpsi /= 'INPUT_PSI_EMPTY') then
     energs%ebs=energs%ekin+energs%epot+energs%eproj !the potential energy contains also exctX
     !write this warning only if the system is closed shell
     call check_closed_shell(KSwfn%orbs,lcs)
     if (abs(energs%evsum-energs%ebs) > 1.d-8 .and. iproc==0 .and. lcs) then
        call yaml_newline()
        call yaml_mapping_open('Energy inconsistencies')
        call yaml_map('Band Structure Energy',energs%ebs,fmt='(1pe22.14)')
        call yaml_map('Sum of Eigenvalues',energs%evsum,fmt='(1pe22.14)')
        if (energs%evsum /= 0.0_gp) call yaml_map('Relative inconsistency',(energs%ebs-energs%evsum)/energs%evsum,fmt='(1pe9.2)')
        call yaml_mapping_close()
        !write( *,'(1x,a,2(1x,1pe20.13))')&
        !  &   'Difference:evsum,energybs',energs%evsum,energs%ebs
     end if
  end if
  ! Clean KSwfn parts only needed in the SCF loop.
  call kswfn_free_scf_data(KSwfn, (nproc > 1))
  ! Clean denspot parts only needed in the SCF loop.
  call denspot_free_history(denspot)

  call f_release_routine()
END SUBROUTINE kswfn_optimization_loop


subroutine kswfn_post_treatments(iproc, nproc, KSwfn, tmb, linear, &
     & fxyz, fnoise, fion, fdisp, fpulay, &
     & strten, pressure, ewaldstr, xcstr, &
     & GPU, denspot, atoms, rxyz, nlpsp, &
     & output_denspot, dir_output, gridformat, refill_proj, &
     & calculate_dipole, calculate_quadrupole, calculate_strten,nspin, &
     & plot_pot_axes)
  use module_base
  use module_types
  use module_interfaces, only: XC_potential,sumrho,communicate_density
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use yaml_output
  use communications_base, only: deallocate_comms_linear, deallocate_p2pComms
  use communications, only: synchronize_onesided_communication
  use sparsematrix_base, only: deallocate_matrices, deallocate_sparse_matrix
  use multipole, only: calculate_dipole_moment
  use public_enums
  use orbitalbasis
  use io, only: plot_density
  use module_xc, only: XC_NO_HARTREE
  use PSbox
  use box
  implicit none
  !Arguments
  type(DFT_wavefunction), intent(in) :: KSwfn
  type(DFT_wavefunction), intent(inout) :: tmb
  type(GPU_pointers), intent(inout) :: GPU
  type(DFT_local_fields), intent(inout) :: denspot
  type(atoms_data), intent(in) :: atoms
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  logical, intent(in) :: linear, refill_proj, calculate_dipole, calculate_quadrupole, calculate_strten
  integer, intent(in) :: iproc, nproc, nspin
  type(f_enumerator), intent(in) :: output_denspot
  character(len = *), intent(in) :: dir_output
  character(len = *), intent(in) :: gridformat
  real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(3, atoms%astruct%nat), intent(in) :: fdisp,fpulay
  real(gp), dimension(3, atoms%astruct%nat), intent(inout) :: fion
  real(dp), dimension(6), intent(in) :: ewaldstr
  real(dp), dimension(6), intent(inout) :: xcstr
  real(gp), intent(out) :: fnoise, pressure
  real(gp), dimension(6), intent(out) :: strten
  real(gp), dimension(3, atoms%astruct%nat), intent(out) :: fxyz
  integer,dimension(3),intent(in) :: plot_pot_axes

  !Local variables
  character(len = *), parameter :: subname = "kswfn_post_treatments"
  integer ::  jproc, nsize_psi, imode, i3xcsh_old,ioffset,ishift
  !integer :: i,ispin
  real(dp), dimension(6) :: hstrten
  real(dp), dimension(:,:), pointer :: rho_p
  real(dp), dimension(:), allocatable :: fpcm
  real(gp) :: ehart_fake, exc_fake, evxc_fake
  type(orbital_basis) :: ob
  type(PSolver_energies) :: PSenergies

  call f_routine(id=subname)
  !manipulate scatter array for avoiding the GGA shift
!!$     call dpbox_repartition(denspot%dpbox%iproc,denspot%dpbox%nproc,atoms%astruct%geocode,'D',1,denspot%dpbox)
  !n3d=n3p
  denspot%dpbox%n3d=denspot%dpbox%n3p
  !i3xcsh=0
  denspot%dpbox%i3s=denspot%dpbox%i3s+denspot%dpbox%i3xcsh
  !save the value for future reference for the core density
  i3xcsh_old=denspot%dpbox%i3xcsh
  denspot%dpbox%i3xcsh=0
  do jproc=0,denspot%dpbox%mpi_env%nproc-1
     !n3d=n3p
     denspot%dpbox%nscatterarr(jproc,1)=denspot%dpbox%nscatterarr(jproc,2)
     !i3xcsh=0
     denspot%dpbox%nscatterarr(jproc,4)=0
     !the same for the density
     denspot%dpbox%ngatherarr(:,3)=denspot%dpbox%ngatherarr(:,1)
  end do
  !change communication scheme to LDA case
  !only in the case of no PSolver tasks
  if (denspot%dpbox%mpi_env%nproc < nproc) then
     denspot%rhod%icomm=0
     denspot%rhod%nrhotot=denspot%dpbox%mesh%ndims(3)
  else
     denspot%rhod%icomm=1
     denspot%rhod%nrhotot=sum(denspot%dpbox%nscatterarr(:,1))
  end if

  if (linear) then
!!$     !then calculate the contribution coming from the Hartree part
!!$     ! keep only the essential part of the density, without the GGA bufffers
!!$     ioffset=kswfn%lzd%glr%d%n1i*kswfn%lzd%glr%d%n2i*denspot%dpbox%i3xcsh
!!$     if (denspot%dpbox%ndimrhopot>0) then
!!$        call f_memcpy(n=denspot%dpbox%ndimpot,src=denspot%rhov(ioffset+1),&
!!$             dest=denspot%rho_work(1))
!!$        ! add the spin down part if present
!!$        if (denspot%dpbox%nrhodim==2) then
!!$           ishift=denspot%dpbox%ndimrhopot/denspot%dpbox%nrhodim !start of the spin down part
!!$           call axpy(denspot%dpbox%ndimpot, 1.d0, &
!!$                denspot%rhov(ioffset+ishift+1), &
!!$                1, denspot%rho_work(1),1)
!!$        end if
!!$     end if

     ! This seems to be necessary to get the correct value of xcstr.
     call XC_potential(atoms%astruct%geocode,'D',iproc,nproc,bigdft_mpi%mpi_comm,&
          KSwfn%Lzd%Glr%d%n1i,KSwfn%Lzd%Glr%d%n2i,KSwfn%Lzd%Glr%d%n3i,denspot%xc,&
          denspot%dpbox%mesh%hgrids,&
          denspot%rhov,exc_fake,evxc_fake,nspin,denspot%rho_C,denspot%rhohat,&
          denspot%V_XC,xcstr)
  else
     denspot%rho_work = f_malloc_ptr(max(denspot%dpbox%ndimrhopot,1),id='denspot%rho_work')
     nullify(rho_p)
     call sumrho(denspot%dpbox,KSwfn%orbs,KSwfn%Lzd,GPU,atoms%astruct%sym,denspot%rhod,denspot%xc,&
          KSwfn%psi,rho_p)
     call communicate_density(denspot%dpbox,KSwfn%orbs%nspin,denspot%rhod,rho_p,denspot%rho_work,.false.)
     !Calculate the total density in the case of nspin==2
     if (KSwfn%orbs%nspin==2) then
        call axpy(denspot%dpbox%ndimpot,1.0_dp,denspot%rho_work(1+denspot%dpbox%ndimpot),1,denspot%rho_work(1),1)
     end if
  end if
  !xc stress, diagonal for the moment
  if (atoms%astruct%geocode=='P' .and. atoms%astruct%sym%symObj >= 0) &
       call symm_stress(xcstr,atoms%astruct%sym%symObj)

  denspot%pot_work = f_malloc_ptr(max(denspot%dpbox%ndimpot,1),id='denspot%pot_work')
  if (denspot%xc%id(1) /= XC_NO_HARTREE) then
     !Calculate electrostatic potential
     call f_memcpy(n=denspot%dpbox%ndimpot,src=denspot%rho_work(1),dest=denspot%pot_work(1))
     !the ionic denisty is given in the case of the embedded solver
     call PS_set_options(denspot%pkernel,calculate_strten=.true.,&
          final_call=denspot%pkernel%method .hasattr. 'rigid',&
          cavity_info=output_denspot /= ENUM_EMPTY)
     if (denspot%pkernel%method /= 'VAC') then
        call Electrostatic_Solver(denspot%pkernel,denspot%pot_work,PSenergies,rho_ion=denspot%rho_ion)
        !here the ionic forces can be corrected with the value coming from the cavity
        if (denspot%pkernel%method .hasattr. 'rigid') then
           fpcm=f_malloc0(3*atoms%astruct%nat,id='fpcm')
           call ps_soft_PCM_forces(denspot%pkernel,fpcm)
           call PS_reduce(fpcm,denspot%pkernel)
           !call f_axpy(a=f_1,x=fpcm,y=fion)
           !fion=fion .plus. fpcm
           call axpy(3*atoms%astruct%nat,1.0_dp,fpcm(1),1,fion(1,1),1)
           call f_free(fpcm)
        end if
     else
        call Electrostatic_Solver(denspot%pkernel,denspot%pot_work,PSenergies)
     end if
  else
     call f_zero(denspot%pot_work)
  end if
  hstrten=PSenergies%strten

  !In principle symmetrization of the stress tensor is not needed since the density has been 
  !already symmetrized
  if (atoms%astruct%sym%symObj >= 0 .and. cell_geocode(denspot%pkernel%mesh)=='P') &
       call symm_stress(hstrten,atoms%astruct%sym%symObj)

  !SM: for a spin polarized calculation, rho_work already contains the full
  !density in the first half of the array. Therefore I think that calc_dipole should be
  !called with nspin=1 and not nspin=2 as it used to be.
  if (calculate_dipole) then
     ! calculate dipole moment associated to the charge density
     !call calc_dipole(denspot%dpbox,denspot%dpbox%nrhodim,atoms,rxyz,denspot%rho_work,.false.)
     call calculate_dipole_moment(denspot%dpbox,1,atoms,rxyz,denspot%rho_work,&
          calculate_quadrupole)
  end if
  !plot the density on the cube file
  !to be done either for post-processing or if a restart is to be done with mixing enabled
  if (output_denspot /= ENUM_EMPTY) then
     if (iproc == 0) call yaml_map('Writing electronic density in file','electronic_density'//gridformat)

     if (all(plot_pot_axes>=0)) then
        call plot_density(iproc,nproc,trim(dir_output)//'electronic_density' // gridformat,&
             atoms,rxyz,denspot%pkernel,denspot%dpbox%nrhodim,denspot%rho_work, &
             ixyz0=plot_pot_axes)
     else if (any(plot_pot_axes>=0)) then
        call f_err_throw('The coordinates of the point through which '//&
             'the density shall be plotted must all be non-zero',&
             err_name='BIGDFT_RUNTIME_ERROR')
     else
         call plot_density(iproc,nproc,trim(dir_output)//'electronic_density' // gridformat,&
              atoms,rxyz,denspot%pkernel,denspot%dpbox%nrhodim,denspot%rho_work)
      end if

     call PS_dump_coulomb_operator(denspot%pkernel,trim(dir_output))
!---------------------------------------------------

     ! SM: the check whether denspot%dpbox%n3d>0 might lead to deadlocks (plot_density contains
     ! an MPI gather. To avoid out of bounds errors, rho_C is now allocted min max(n3d,1).
     !if (associated(denspot%rho_C) .and. denspot%dpbox%n3d>0) then
     if (associated(denspot%rho_C)) then
        if (iproc == 0) call yaml_map('Writing core density in file','core_density'//gridformat)
        call plot_density(iproc,nproc,trim(dir_output)//'core_density' // gridformat,&
             atoms,rxyz,denspot%pkernel,1,denspot%rho_C(1:,1:,i3xcsh_old+1:,1:))
     end if
     if (associated(denspot%rhohat)) then
        if (iproc == 0) call yaml_map('Writing compensation density in file', 'hat_density'//gridformat)
        call plot_density(iproc,nproc,trim(dir_output)//'hat_density' // gridformat,&
             atoms,rxyz,denspot%pkernel,1,denspot%rhohat)
     end if
  end if
  !plot also the electrostatic potential
  if (output_denspot == 'DENSPOT') then
     if (iproc == 0) call yaml_map('Writing Hartree potential in file','hartree_potential'//gridformat)
     if (all(plot_pot_axes>=0)) then
        call plot_density(iproc,nproc,trim(dir_output)//'hartree_potential' // gridformat, &
             atoms,rxyz,denspot%pkernel,1,denspot%pot_work, &
             ixyz0=plot_pot_axes)
     else if (any(plot_pot_axes>=0)) then
        call f_err_throw('The coordinates of the point through which '//&
             'the potential shall be plotted must all be non-zero',&
             err_name='BIGDFT_RUNTIME_ERROR')
     else 
        call plot_density(iproc,nproc,trim(dir_output)//'hartree_potential' // gridformat, &
             atoms,rxyz,denspot%pkernel,1,denspot%pot_work)
     end if
  end if

  call timing(iproc,'Forces        ','ON')

  ! Calculate the forces. Pass the Pulay forces in the linear scaling case.
  if (linear) then
     imode = 1
     nsize_psi=1
     ! This is just to save memory, since calculate_forces will require quite a lot
!     call deallocate_comms_linear(tmb%collcom)
!     call deallocate_comms_linear(tmb%ham_descr%collcom)
!     call deallocate_comms_linear(tmb%collcom_sr)
!     call deallocate_p2pcomms(tmb%comgp)
!     call deallocate_p2pcomms(tmb%ham_descr%comgp)
!     do i=1,size(tmb%linmat%ovrlppowers_)
!         call deallocate_matrices(tmb%linmat%ovrlppowers_(i))
!     end do
!     call deallocate_matrices(tmb%linmat%ham_)
!     call deallocate_matrices(tmb%linmat%ovrlp_)
!     call deallocate_sparse_matrix(tmb%linmat%smat(1))
!     call deallocate_sparse_matrix(tmb%linmat%smat(2))
!     if (associated(tmb%linmat%ks)) then
!         do ispin=1,tmb%linmat%smat(3)%nspin
!             call deallocate_sparse_matrix(tmb%linmat%ks(ispin))
!         end do
!         deallocate(tmb%linmat%ks)
!     end if
!     if (associated(tmb%linmat%ks_e)) then
!         do ispin=1,tmb%linmat%smat(3)%nspin
!             call deallocate_sparse_matrix(tmb%linmat%ks_e(ispin))
!         end do
!         deallocate(tmb%linmat%ks_e)
!     end if

  else
     imode = 0
     nsize_psi = (KSwfn%Lzd%Glr%wfd%nvctr_c+7*KSwfn%Lzd%Glr%wfd%nvctr_f)*KSwfn%orbs%nspinor*KSwfn%orbs%norbp
  end if
  if (.not. KSwfn%paw%usepaw) then
     !to be checked, as the description might change for the linear case
     if (linear) then
        call orbital_basis_associate(ob,orbs=KSwfn%orbs,Lzd=KSwfn%Lzd)
     else
        call orbital_basis_associate(ob,orbs=KSwfn%orbs,phis_wvl=KSwfn%psi,Lzd=KSwfn%Lzd)        
     end if
  call calculate_forces(iproc,nproc,denspot%pkernel%mpi_env%nproc,KSwfn%Lzd%Glr,atoms,ob,nlpsp,rxyz,&
          KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
       denspot%dpbox, &
          denspot%dpbox%i3s+denspot%dpbox%i3xcsh,denspot%dpbox%n3p,&
          denspot%dpbox%nrhodim,refill_proj,denspot%dpbox%ngatherarr,denspot%rho_work,&
          denspot%pot_work,denspot%V_XC,nsize_psi,KSwfn%psi,fion,fdisp,fxyz,&
       calculate_strten,ewaldstr,hstrten,xcstr,strten,pressure,denspot%psoffset,imode,tmb,KSwfn%paw,fpulay)
  call orbital_basis_release(ob)
  end if

  !here we might printout the atomic density matrix

  call f_free_ptr(denspot%rho_work)
  call f_free_ptr(denspot%pot_work)
  if (linear) then
     ! to eventually be better sorted
!     call synchronize_onesided_communication(iproc, nproc, tmb%ham_descr%comgp)
!     call deallocate_p2pComms(tmb%ham_descr%comgp)
!     call deallocate_local_zone_descriptors(tmb%ham_descr%lzd)
!     call deallocate_comms_linear(tmb%ham_descr%collcom)
!     call deallocate_auxiliary_basis_function(subname, tmb%ham_descr%psi, tmb%hpsi)

!!!! TEST ##################
     !!fxyz=0.d0
     !!tmb%psi(1:KSwfn%orbs%npsidim_orbs)=KSwfn%psi(1:KSwfn%orbs%npsidim_orbs)
     !!tmb%wfnmd%density_kernel=0.d0
     !!do i_stat=1,KSwfn%orbs%norb
     !!    tmb%wfnmd%density_kernel(i_stat,i_stat)=1.d0
     !!end do
     !!call  nonlocal_forces(tmb%lzd%glr,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
     !! atoms,rxyz,&
     !! KSwfn%orbs,nlpsp,proj,tmb%lzd%glr%wfd,KSwfn%psi,fxyz,refill_proj,strten)
     !!call nonlocal_forces_linear(iproc,nproc,tmb%lzd%glr,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),&
     !!     KSwfn%Lzd%hgrids(3),atoms,rxyz,&
     !!     tmb%orbs,nlpsp,proj,tmb%lzd,tmb%psi,tmb%wfnmd%density_kernel,fxyz,refill_proj,strten)
     !!call nonlocal_forces_linear(iproc,nproc,tmb%lzd%glr,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),&
     !!     KSwfn%Lzd%hgrids(3),atoms,rxyz,&
     !!     tmb%orbs,nlpsp,proj,tmb%ham_descr%lzd,tmb%ham_descr%psi,tmb%wfnmd%density_kernel,fxyz,refill_proj,strten)
     !!if (nproc > 1) then
     !!   call mpiallred(fxyz(1,1),3*atoms%astruct%nat,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
     !!end if
     !!if (iproc==0) then
     !!     do iat=1,atoms%astruct%nat
     !!         write(*,'(a,3es18.8)') 'new forces',fxyz(1,iat), fxyz(2,iat), fxyz(3,iat)
     !!     end do
     !!end if
!!!! #######################
  end if

  !!stop
  call timing(iproc,'Forces        ','OF')
  call f_release_routine()
END SUBROUTINE kswfn_post_treatments
