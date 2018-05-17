!> @file
!!  Routines related to system properties
!! @author
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Initialize the objects needed for the computation: basis sets, allocate required space
subroutine system_initialization(iproc,nproc,dump,inputpsi,input_wf_format,dry_run,&
     & in,atoms,rxyz,OCLconv,&
     orbs,lnpsidim_orbs,lnpsidim_comp,lorbs,Lzd,Lzd_lin,nlpsp,comms,&
     ref_frags, denspot, locregcenters, inwhichlocreg_old, onwhichatom_old, &
     norb_par_ref, norbu_par_ref, norbd_par_ref,output_grid)
  use module_base
  use module_types
  use module_interfaces, only: createWavefunctionsDescriptors, &
       & init_orbitals_data_for_linear, orbitals_descriptors
  use module_xc
  use module_fragments
  use vdwcorrection
  use yaml_output
  use module_atoms, only: set_symmetry_data
  use communications_base, only: comms_cubic
  use communications_init, only: orbitals_communicators
  use Poisson_Solver, only: pkernel_allocate_cavity
  use public_enums
  use f_enums
  use locreg_operations
  use locregs_init, only: initLocregs,lr_set
  use orbitalbasis
  use chess_base, only: chess_init
  use module_dpbox, only: dpbox_set
  use rhopotential, only: set_cfd_data
  implicit none
  integer, intent(in) :: iproc,nproc 
  logical, intent(in) :: dry_run, dump
  integer, intent(out) :: input_wf_format, lnpsidim_orbs, lnpsidim_comp
  type(f_enumerator), intent(inout) :: inputpsi
  type(input_variables), intent(inout) :: in !SM: I had to change this to inout due to the new CheSS parameters within input_variables, maybe to be changed again later...
  type(atoms_data), intent(inout) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
  logical, intent(in) :: OCLconv
  type(orbitals_data), intent(inout) :: orbs, lorbs
  type(local_zone_descriptors), intent(inout) :: Lzd, Lzd_lin
  type(DFT_PSP_projectors), intent(out) :: nlpsp
  type(comms_cubic), intent(out) :: comms
  !real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
  type(system_fragment), dimension(:), pointer :: ref_frags
  real(kind=8),dimension(3,atoms%astruct%nat),intent(inout),optional :: locregcenters
  integer,dimension(:),pointer,optional:: inwhichlocreg_old, onwhichatom_old
  integer,dimension(0:nproc-1),optional:: norb_par_ref, norbu_par_ref, norbd_par_ref !< support function distribution to be used as a reference
  type(DFT_local_fields), intent(out), optional :: denspot
  logical, intent(in), optional :: output_grid
  !local variables
  character(len = *), parameter :: subname = "system_initialization"
  integer :: nB,nKB,nMB,ii,iat,iorb,nspin_ig,norbe,norbsc,ifrag,nspinor
  real(gp), dimension(3) :: h_input
  logical:: present_inwhichlocreg_old, present_onwhichatom_old, output_grid_, frag_allocated, calculate_bounds
  integer, dimension(:,:), allocatable :: norbsc_arr
  integer, dimension(:), allocatable :: norb_par, norbu_par, norbd_par
  real(kind=8), dimension(:), allocatable :: locrad, times_convol
  integer :: ilr, iilr
  real(kind=8),dimension(:),allocatable :: totaltimes, locrads
  real(kind=8),dimension(2) :: time_max, time_average
  !real(kind=8) :: ratio_before, ratio_after
  logical :: init_projectors_completely
  type(orbital_basis) :: ob
  call f_routine(id=subname)


  output_grid_ = .false.
  if (present(output_grid)) output_grid_ = output_grid

  if (iproc == 0 .and. dump) &
       & call print_atomic_variables(atoms, max(in%hx,in%hy,in%hz), in%ixc)

  !grid spacings of the zone descriptors (not correct, the set is done by system size)
  Lzd=default_lzd()
  !h_input=(/ in%hx, in%hy, in%hz /)
  !call lzd_set_hgrids(Lzd,h_input) 
  Lzd%hgrids=(/ in%hx, in%hy, in%hz /) !to be adjusted with the constraints of the box
  ! Create wavefunctions descriptors and allocate them inside the global locreg desc.
  calculate_bounds = .not. (inputpsi .hasattr. 'LINEAR')
  call lr_set(lzd%Glr,iproc,OCLconv,dump,in%crmult,in%frmult,lzd%hgrids,rxyz,atoms,&
       calculate_bounds,output_grid_)

  if (present(locregcenters)) then
      do iat=1,atoms%astruct%nat
          locregcenters(1:3,iat)=locregcenters(1:3,iat)-atoms%astruct%shift(1:3)
          if (any(locregcenters(:,iat) < 0.0_gp .or. &
               locregcenters(:,iat) > lzd%glr%mesh_coarse%ndims*lzd%glr%mesh_coarse%hgrids)) then
!!$          if (locregcenters(1,iat)<dble(0)*lzd%hgrids(1) .or. locregcenters(1,iat)>dble(lzd%glr%d%n1+1)*lzd%hgrids(1) .or. &
!!$              locregcenters(2,iat)<dble(0)*lzd%hgrids(2) .or. locregcenters(2,iat)>dble(lzd%glr%d%n2+1)*lzd%hgrids(2) .or. &
!!$              locregcenters(3,iat)<dble(0)*lzd%hgrids(3) .or. locregcenters(3,iat)>dble(lzd%glr%d%n3+1)*lzd%hgrids(3)) then
             call f_err_throw('locregcenter outside of global box!', err_name='BIGDFT_RUNTIME_ERROR')
          end if
       end do
  end if

  ! Initialize the object holding the CheSS parameters
  call chess_init(in%chess_dict, in%cp)

  if (present(denspot)) then
     call initialize_DFT_local_fields(denspot, in%ixc, in%nspin, in%alpha_hartree_fock)

     !here the initialization of dpbox can be set up
     call dpbox_set(denspot%dpbox,Lzd%Glr%mesh,denspot%xc,iproc,nproc,bigdft_mpi%mpi_comm, &
          in%SIC%approach, in%nspin)

     ! Create the Poisson solver kernels.
     call system_initKernels(.true.,iproc,nproc,atoms%astruct%geocode,in,denspot)
     call system_createKernels(denspot, (get_verbose_level() > 1))
     if (denspot%pkernel%method .hasattr. 'rigid') then
        call epsilon_cavity(atoms,rxyz,denspot%pkernel)
        !allocate cavity, in the case of nonvacuum treatment
     else if (denspot%pkernel%method /= 'VAC') then 
          call pkernel_allocate_cavity(denspot%pkernel,&
          vacuum=.not. (denspot%pkernel%method .hasattr. 'sccs'))

          call epsinnersccs_cavity(atoms,rxyz,denspot%pkernel)

     end if
  end if

  ! Create global orbs data structure.
  if(in%nspin==4) then
     nspinor=4
  else
     nspinor=1
  end if
  call orbitals_descriptors(iproc, nproc,in%gen_norb,in%gen_norbu,in%gen_norbd,in%nspin,nspinor,&
       in%gen_nkpt,in%gen_kpt,in%gen_wkpt,orbs,LINEAR_PARTITION_NONE)
  !!write(*,*) 'orbs%norbu', orbs%norbu
  !!write(*,*) 'orbs%norbd', orbs%norbd
  !!write(*,*) 'orbs%norb', orbs%norb
  !!write(*,*) 'orbs%norbup', orbs%norbup
  !!write(*,*) 'orbs%norbdp', orbs%norbdp
  !!write(*,*) 'orbs%norbp', orbs%norbp
  orbs%occup(1:orbs%norb*orbs%nkpts) = in%gen_occup
  if (dump .and. iproc==0) call print_orbitals(orbs, atoms%astruct%geocode)

  ! See if linear scaling should be activated and build the correct Lzd 
  call check_linear_and_create_Lzd(iproc,nproc,in%linear,Lzd,atoms,orbs,in%nspin,rxyz)
  lzd_lin=default_lzd()
  call nullify_local_zone_descriptors(lzd_lin)
  lzd_lin%nlr = 0

  ! Create linear orbs data structure.
  !if (inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR &
  !    .or. inputpsi == INPUT_PSI_MEMORY_LINEAR) then
  if (inputpsi .hasattr. 'LINEAR') then
     if (.not. (inputpsi .hasattr. 'MEMORY')) then ! == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR) then
        ! First do a simple redistribution
        call init_linear_orbs(LINEAR_PARTITION_SIMPLE)
     else
        ! Directly used the reference distribution
        norb_par = f_malloc(0.to.nproc-1,id='norb_par')
        norbu_par = f_malloc(0.to.nproc-1,id='norbu_par')
        norbd_par = f_malloc(0.to.nproc-1,id='norbd_par')
        if (.not.present(norb_par_ref)) then
           call f_err_throw('norb_par_ref not present', err_name='BIGDFT_RUNTIME_ERROR')
        end if
        call vcopy(nproc, norb_par_ref(0), 1, norb_par(0), 1)
        if (.not.present(norbu_par_ref)) then
           call f_err_throw('norbu_par_ref not present', err_name='BIGDFT_RUNTIME_ERROR')
        end if
        call vcopy(nproc, norbu_par_ref(0), 1, norbu_par(0), 1)
        if (.not.present(norbd_par_ref)) then
           call f_err_throw('norbd_par_ref not present', err_name='BIGDFT_RUNTIME_ERROR')
        end if
        call vcopy(nproc, norbd_par_ref(0), 1, norbd_par(0), 1)
        call init_linear_orbs(LINEAR_PARTITION_OPTIMAL)
        call f_free(norb_par)
        call f_free(norbu_par)
        call f_free(norbd_par)
     end if
     call fragment_stuff()
     call init_lzd_linear()
     ! For restart calculations, the suport function distribution must not be modified
     !if (inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR .or. in%lin%fragment_calculation) then
     !SM: added the ".or. fin%lin%fragment_calculation", as this came from a merge with Laura...
     if (.not. (inputpsi .hasattr. 'MEMORY') .or. in%lin%fragment_calculation) then
         times_convol = f_malloc(lorbs%norb,id='times_convol')
         call test_preconditioning()
         time_max(1) = sum(times_convol(lorbs%isorb+1:lorbs%isorb+lorbs%norbp))
         time_average(1) = time_max(1)/real(nproc,kind=8)
         norb_par = f_malloc(0.to.nproc-1,id='norb_par')
         norbu_par = f_malloc(0.to.nproc-1,id='norbu_par')
         norbd_par = f_malloc(0.to.nproc-1,id='norbd_par')
         call optimize_loadbalancing2()
         call deallocate_orbitals_data(lorbs)
         call deallocate_local_zone_descriptors(lzd_lin)
         lzd_lin=default_lzd()
         call nullify_local_zone_descriptors(lzd_lin)
         lzd_lin%nlr = 0
         ! Deallocate here fragment stuff
         !if (.not.(frag_allocated .and. (.not. in%lin%fragment_calculation) .and. inputpsi /= INPUT_PSI_DISK_LINEAR)) then
         !if (frag_allocated) then
         !if (inputpsi == INPUT_PSI_DISK_LINEAR .or. in%lin%fragment_calculation) then
         if (frag_allocated) then
             do ifrag=1,in%frag%nfrag_ref
                call fragment_free(ref_frags(ifrag))
                ref_frags(ifrag)%astruct_frg%nat=-1
                ref_frags(ifrag)%fbasis%forbs=minimal_orbitals_data_null()
                !ref_frags(ifrag)=fragment_null()
                !!call f_free_ptr(ref_frags(ifrag)%astruct_frg%iatype)
                !!call f_free_ptr(ref_frags(ifrag)%astruct_frg%ifrztyp)
                !!call f_free_ptr(ref_frags(ifrag)%astruct_frg%input_polarization)
                !!call f_free_ptr(ref_frags(ifrag)%astruct_frg%rxyz)
                !!call f_free_ptr(ref_frags(ifrag)%astruct_frg%rxyz_int)
                !!call f_free_ptr(ref_frags(ifrag)%astruct_frg%ixyz_int)
             end do
            deallocate(ref_frags)
         end if
         call init_linear_orbs(LINEAR_PARTITION_OPTIMAL)
         totaltimes = f_malloc0(nproc,id='totaltimes')
         call fragment_stuff()
         call init_lzd_linear()
         call test_preconditioning()
         time_max(2) = sum(times_convol(lorbs%isorb+1:lorbs%isorb+lorbs%norbp))
         time_average(2) = time_max(2)/real(nproc,kind=8)
         totaltimes(iproc+1) = time_max(2)
         if (nproc>1) then
             call fmpi_allreduce(time_max, FMPI_MAX, comm=bigdft_mpi%mpi_comm)
             call fmpi_allreduce(time_average, FMPI_SUM, comm=bigdft_mpi%mpi_comm)
             call fmpi_allreduce(totaltimes, FMPI_SUM, comm=bigdft_mpi%mpi_comm)
         end if
         !ratio_before = real(time_max(1),kind=8)/real(max(1.d0,time_min(1)),kind=8) !max to prevent divide by zero
         !ratio_after = real(time_max(2),kind=8)/real(max(1.d0,time_min(2)),kind=8) !max to prevent divide by zero
         !if (iproc==0) call yaml_map('preconditioning load balancing min/max before',(/time_min(1),time_max(1)/),fmt='(es9.2)')
         !if (iproc==0) call yaml_map('preconditioning load balancing min/max after',(/time_min(2),time_max(2)/),fmt='(es9.2)')
         if (iproc==0) call yaml_map('preconditioning load balancing before',time_max(1)/time_average(1),fmt='(es9.2)')
         if (iproc==0) call yaml_map('preconditioning load balancing after',time_max(2)/time_average(2),fmt='(es9.2)')
         if (iproc==0) call yaml_map('task with max load',maxloc(totaltimes)-1)
         call f_free(norb_par)
         call f_free(norbu_par)
         call f_free(norbd_par)
         call f_free(times_convol)
         call f_free(totaltimes)
     end if
  end if

  !In the case in which the number of orbitals is not "trivial" check whether they are too many
  if (inputpsi /= 'INPUT_PSI_RANDOM') then

     ! Allocations for readAtomicOrbitals (check inguess.dat and psppar files)
     norbsc_arr = f_malloc((/ atoms%natsc+1, in%nspin /),id='norbsc_arr')
     locrad = f_malloc(atoms%astruct%nat,id='locrad')

     !calculate the inputguess orbitals
     !spin for inputguess orbitals
     if (in%nspin==4) then
        nspin_ig=1
     else
        nspin_ig=in%nspin
     end if

     ! Read the inguess.dat file or generate the input guess via the inguess_generator
     call readAtomicOrbitals(atoms,norbe,norbsc,nspin_ig,orbs%nspinor,&
          norbsc_arr,locrad)

     if (in%nspin==4) then
        !in that case the number of orbitals doubles
        norbe=2*norbe
     end if

     ! De-allocations
     call f_free(locrad)
     call f_free(norbsc_arr)

     !Check if orbitals and electrons are present
     if (orbs%norb*orbs%nkpts == 0) &
        & call f_err_throw('No electrons in the system! Check your input variables or atomic positions.', &
        & err_id=BIGDFT_INPUT_VARIABLES_ERROR)
     ! Check the maximum number of orbitals
     if (in%nspin==1 .or. in%nspin==4) then
        if (orbs%norb>norbe) then
           call f_err_throw('The number of orbitals ('+yaml_toa(orbs%norb)// &
                &   ') must not be greater than the number of orbitals ('+yaml_toa(norbe)// &
                &   ') generated from the input guess.',err_id=BIGDFT_INPUT_VARIABLES_ERROR)
        end if
     else if (in%nspin == 2) then
        if (orbs%norbu > norbe) then
           call f_err_throw('The number of orbitals up ('+yaml_toa(orbs%norbu)// &
                &   ') must not be greater than the number of orbitals ('+yaml_toa(norbe)// &
                &   ') generated from the input guess.',err_id=BIGDFT_INPUT_VARIABLES_ERROR)
        end if
        if (orbs%norbd > norbe) then
           call f_err_throw('The number of orbitals down ('+yaml_toa(orbs%norbd) //&
                &   ') must not be greater than the number of orbitals ('+yaml_toa(norbe) //&
                &   ') generated from the input guess.',err_id=BIGDFT_INPUT_VARIABLES_ERROR)
        end if
     end if
  end if

  !allocate communications arrays (allocate it before Projectors because of the definition
  !of iskpts and nkptsp)
  call orbitals_communicators(iproc,nproc,Lzd%Glr,orbs,comms)  
  !if (inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR &
  !.or. inputpsi == INPUT_PSI_MEMORY_LINEAR) then
  if ((inputpsi .hasattr. 'LINEAR') .and. .not. dry_run) then
     if(iproc==0 .and. dump) call print_orbital_distribution(iproc, nproc, lorbs)
  end if

  if (.not.(inputpsi .hasattr. 'LINEAR') .and. iproc == 0 .and. dump) then
     nB=max(orbs%npsidim_orbs,orbs%npsidim_comp)*8
     nMB=nB/1024/1024
     nKB=(nB-nMB*1024*1024)/1024
     nB=modulo(nB,1024)
     call yaml_map('Wavefunctions memory occupation for root MPI process',&
          trim(yaml_toa(nMB,fmt='(i5)'))//' MB'//trim(yaml_toa(nKB,fmt='(i5)'))//&
          ' KB'//trim(yaml_toa(nB,fmt='(i5)'))//' B')
 !!$     write(*,'(1x,a,3(i5,a))') &
 !!$       'Wavefunctions memory occupation for root MPI process: ',&
 !!$       nMB,' MB ',nKB,' KB ',nB,' B'
  end if
  ! Done orbs

  !!! fragment initializations - if not a fragment calculation, set to appropriate dummy values
  !if (inputpsi /= INPUT_PSI_LINEAR_AO .and. inputpsi /= INPUT_PSI_DISK_LINEAR &
  !    .and. inputpsi /= INPUT_PSI_MEMORY_LINEAR) then
  if (.not. (inputpsi .hasattr. 'LINEAR')) then
      call fragment_stuff()
  end if
  !!frag_allocated=.false.
  !!if (inputpsi == INPUT_PSI_DISK_LINEAR .or. in%lin%fragment_calculation) then
  !!   allocate(ref_frags(in%frag%nfrag_ref))
  !!   do ifrag=1,in%frag%nfrag_ref
  !!      ref_frags(ifrag)=fragment_null()
  !!   end do
  !!   call init_fragments(in,lorbs,atoms%astruct,ref_frags)
  !!  frag_allocated=.true.
  !!else
  !!   nullify(ref_frags)
  !!end if

  !!call input_check_psi_id(inputpsi, input_wf_format, in%dir_output, &
  !!     orbs, lorbs, iproc, nproc, in%frag%nfrag_ref, in%frag%dirname, ref_frags)

  !!! we need to deallocate the fragment arrays we just allocated as not a restart calculation so this is no longer needed
  !!if (frag_allocated .and. (.not. in%lin%fragment_calculation) .and. inputpsi /= INPUT_PSI_DISK_LINEAR) then
  !!    do ifrag=1,in%frag%nfrag_ref
  !!       ref_frags(ifrag)%astruct_frg%nat=-1
  !!       ref_frags(ifrag)%fbasis%forbs=minimal_orbitals_data_null()
  !!       call fragment_free(ref_frags(ifrag))
  !!       !ref_frags(ifrag)=fragment_null()
  !!    end do
  !!   deallocate(ref_frags)
  !!end if

  ! Calculate all projectors, or allocate array for on-the-fly calculation
  ! SM: For a linear scaling calculation, some parts can be done later.
  ! SM: The following flag is false for linear scaling and true otherwise.
  init_projectors_completely =   .not.(inputpsi .hasattr. 'LINEAR')
  !(inputpsi /= INPUT_PSI_LINEAR_AO .and. &
  !                              inputpsi /= INPUT_PSI_DISK_LINEAR .and. &
  !                              inputpsi /= INPUT_PSI_MEMORY_LINEAR)
  call orbital_basis_associate(ob,orbs=orbs,Lzd=Lzd,id='system_initialization')
  call createProjectorsArrays(iproc,nproc,Lzd%Glr,rxyz,atoms,ob,&
       in%frmult,in%frmult,Lzd%hgrids(1),Lzd%hgrids(2),&
       Lzd%hgrids(3),dry_run,nlpsp,init_projectors_completely)
  call orbital_basis_release(ob)
  if (iproc == 0 .and. dump) call print_nlpsp(nlpsp)
  if (iproc == 0 .and. .not. nlpsp%on_the_fly .and. .false.) then
     call writemyproj("proj",WF_FORMAT_BINARY,orbs,Lzd%hgrids(1),Lzd%hgrids(2),&
       Lzd%hgrids(3),atoms,rxyz,nlpsp)
  end if
  !the complicated part of the descriptors has not been filled
  if (dry_run) then
     call f_release_routine()
     return
  end if
  !calculate the partitioning of the orbitals between the different processors
!  print *,'here the localization regions should have been filled already'
!  stop

  if (present(denspot)) then
     !here dpbox can be put as input
     call density_descriptors(iproc,nproc,denspot%xc,in%nspin,in%crmult,in%frmult,atoms,&
          denspot%dpbox,in%rho_commun,rxyz,denspot%rhod)
     !allocate the arrays.
     call allocateRhoPot(Lzd%Glr,in%nspin,atoms,rxyz,denspot)
     !here insert the conditional for the constrained field dynamics
     if (in%calculate_magnetic_torque) then
        call set_cfd_data(denspot%cfd,Lzd%Glr%mesh,atoms%astruct,rxyz)
     end if
     if (in%do_spin_dynamics) then
        if (in%calculate_magnetic_torque) then
           print *,'spin dynamics is go for launch'
           !call asd_allocate(asd)
        else
           print *,'spin dynamics needs constraining fields..'
        end if
     end if
  end if

  !calculate the irreductible zone for this region, if necessary.
  call set_symmetry_data(atoms%astruct%sym,atoms%astruct%geocode, &
       & Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i, in%nspin)

  ! A message about dispersion forces.
  call vdwcorrection_initializeparams(in%ixc, in%dispersion)

  !check the communication distribution
  !if(inputpsi /= INPUT_PSI_LINEAR_AO .and. inputpsi /= INPUT_PSI_DISK_LINEAR &
  !   .and. inputpsi /= INPUT_PSI_MEMORY_LINEAR) then
  if (.not.(inputpsi .hasattr. 'LINEAR')) then
     call check_communications(iproc,nproc,orbs,Lzd,comms)
  else
      ! Do not call check_communication, since the value of orbs%npsidim_orbs is wrong
      !if(iproc==0) call yaml_warning('Do not call check_communications in the linear scaling version!')
      !if(iproc==0) write(*,*) 'WARNING: do not call check_communications in the linear scaling version!'
  end if

  call f_release_routine()
  !---end of system definition routine


  contains


    subroutine init_linear_orbs(linear_partition)
      use module_interfaces, only: init_orbitals_data_for_linear
     implicit none
     integer,intent(in) :: linear_partition

     if (linear_partition==LINEAR_PARTITION_SIMPLE) then
         if (present(locregcenters)) then
             call init_orbitals_data_for_linear(iproc, nproc, orbs%nspinor, in, atoms%astruct, locregcenters, lorbs)
         else
             call init_orbitals_data_for_linear(iproc, nproc, orbs%nspinor, in, atoms%astruct, atoms%astruct%rxyz, lorbs)
         end if
     else if (linear_partition==LINEAR_PARTITION_OPTIMAL) then
         if (present(locregcenters)) then
             call init_orbitals_data_for_linear(iproc, nproc, orbs%nspinor, in, atoms%astruct, locregcenters, lorbs, &
                  norb_par, norbu_par, norbd_par)
         else
             call init_orbitals_data_for_linear(iproc, nproc, orbs%nspinor, in, atoms%astruct, atoms%astruct%rxyz, lorbs, &
                  norb_par, norbu_par, norbd_par)
         end if
     else
         !stop 'init_linear_orbs: wrong value of linear_partition'
         call f_err_throw('init_linear_orbs: wrong value of linear_partition',err_name='BIGDFT_RUNTIME_ERROR')
     end if

       ! There are needed for the restart (at least if the atoms have moved...)
       present_inwhichlocreg_old = present(inwhichlocreg_old)
       present_onwhichatom_old = present(onwhichatom_old)
       if (present_inwhichlocreg_old .and. .not.present_onwhichatom_old &
           .or. present_onwhichatom_old .and. .not.present_inwhichlocreg_old) then
           call f_err_throw('inwhichlocreg_old and onwhichatom_old should be present at the same time', &
           & err_name='BIGDFT_INPUT_VARIABLES_ERROR')
           !stop 
       end if
       if (present_inwhichlocreg_old .and. present_onwhichatom_old) then
           call vcopy(lorbs%norb, onwhichatom_old(1), 1, lorbs%onwhichatom(1), 1)
           !call vcopy(lorbs%norb, inwhichlocreg_old(1), 1, lorbs%inwhichlocreg(1), 1)
           !use onwhichatom to build the new inwhichlocreg (because the old inwhichlocreg can be ordered differently)
           ii = 0
           do iat=1, atoms%astruct%nat
              !only want to do this for spin=1!
              do iorb=1,lorbs%norbu
                 if(iat ==  lorbs%onwhichatom(iorb)) then
                    ii = ii + 1
                    lorbs%inwhichlocreg(iorb)= ii
                 end if
              end do 
           end do

           ! LR: not sure if this is the best way to do this but it seems to work...
           ! Correction for spin polarized systems. For non polarized systems, norbu=norb and the loop does nothing.
           do iorb=lorbs%norbu+1,lorbs%norb
               lorbs%inwhichlocreg(iorb)=lorbs%inwhichlocreg(iorb-lorbs%norbu)+lorbs%norbu
           end do

           !i_all=-product(shape(inwhichlocreg_old))*kind(inwhichlocreg_old)
           !deallocate(inwhichlocreg_old,stat=i_stat)
           !call memocc(i_stat,i_all,'inwhichlocreg_old',subname)
           !i_all=-product(shape(onwhichatom_old))*kind(onwhichatom_old)
           !deallocate(onwhichatom_old,stat=i_stat)
           !call memocc(i_stat,i_all,'onwhichatom_old',subname)
       end if
     end subroutine init_linear_orbs

     subroutine init_lzd_linear()
       use module_interfaces, only: initialize_linear_from_file
       use locregs_init, only: initLocregs
       use locregs, only: copy_locreg_descriptors
       use sparsematrix_wrappers, only: check_kernel_cutoff
       implicit none
       call copy_locreg_descriptors(Lzd%Glr, lzd_lin%glr)
       !call lzd_set_hgrids(lzd_lin, Lzd%hgrids)
       lzd_lin%hgrids=Lzd%hgrids
       if (inputpsi == 'INPUT_PSI_LINEAR_AO' .or. inputpsi == 'INPUT_PSI_MEMORY_LINEAR') then
           !!write(*,*) 'rxyz',rxyz
           !!write(*,*) 'locregcenters',locregcenters
           if (present(locregcenters)) then
              call lzd_init_llr(iproc, nproc, in, atoms%astruct, locregcenters, lorbs, lzd_lin)
          else
              call lzd_init_llr(iproc, nproc, in, atoms%astruct, atoms%astruct%rxyz, lorbs, lzd_lin)
          end if

       else
          call initialize_linear_from_file(iproc,nproc,in%frag,atoms%astruct,rxyz,lorbs,lzd_lin,&
               input_wf_format,in%dir_output,'minBasis',ref_frags)
          !what to do with derivatives?
          ! These values are not read from file, not very nice this way
          do ilr=1,lzd_lin%nlr
              iilr=mod(ilr-1,lorbs%norbu)+1 !correct value for a spin polarized system
              lzd_lin%llr(ilr)%locrad_kernel=in%lin%locrad_kernel(iilr)
              lzd_lin%llr(ilr)%locrad_mult=in%lin%locrad_mult(iilr)
          end do
       end if

       ! Make sure that the cutoff for the multiplications
       ! is larger than the kernel cutoff
       call check_kernel_cutoff(iproc, lorbs, atoms, in%hamapp_radius_incr, lzd_lin)
       do ilr=1,lzd_lin%nlr
          if (lzd_lin%llr(ilr)%locrad_mult < lzd_lin%llr(ilr)%locrad_kernel) then
             call f_err_throw('rloc_kernel_foe ('//trim(yaml_toa(lzd_lin%llr(ilr)%locrad_mult,fmt='(f5.2)'))//&
                  &') too small, must be at least as big as rloc_kernel('&
                  &//trim(yaml_toa(lzd_lin%llr(ilr)%locrad_kernel,fmt='(f5.2)'))//')', err_id=BIGDFT_RUNTIME_ERROR)
          end if
       end do

       if (.not. dry_run) then
          locrads = f_malloc(lzd_lin%nlr)
          locrads = lzd_lin%llr(:)%locrad
          call initLocregs(iproc, nproc, lzd_lin, Lzd_lin%hgrids(1), Lzd_lin%hgrids(2),Lzd_lin%hgrids(3), &
               atoms%astruct%rxyz,locrads, lorbs, Lzd_lin%Glr, 's')
          call update_wavefunctions_size(lzd_lin,lnpsidim_orbs,lnpsidim_comp,lorbs,iproc,nproc)
          call f_free(locrads)
       else
          call yaml_warning("Locregs not initialized for linear.")
       end if
     end subroutine init_lzd_linear


     !this routine should go in the locreg_operation module
     subroutine test_preconditioning()
       implicit none

       !Local variables
       integer :: iorb, iiorb, ilr, ncplx, ist, i, ierr, ii, jj
       logical :: with_confpot
       real(gp) :: kx, ky, kz
       type(workarrays_quartic_convolutions),dimension(:),allocatable :: precond_convol_workarrays
       type(workarr_precond),dimension(:),allocatable :: precond_workarrays
       real(kind=8),dimension(:),allocatable :: phi
       real(kind=8) :: t1, t2, time, tt, maxtime
       integer,parameter :: nit=5
       real(kind=8),dimension(2*nit+1) :: times

      call f_zero(times_convol)
      do iorb=1,lorbs%norbp
          iiorb=lorbs%isorb+iorb
          ilr=lorbs%inwhichlocreg(iiorb)
          ii = (lzd_lin%llr(ilr)%d%n1+1)*(lzd_lin%llr(ilr)%d%n2+1)*(lzd_lin%llr(ilr)%d%n3+1)
          jj = 7*(lzd_lin%llr(ilr)%d%nfu1-lzd_lin%llr(ilr)%d%nfl1+1)*&
                 (lzd_lin%llr(ilr)%d%nfu2-lzd_lin%llr(ilr)%d%nfl2+1)*&
                 (lzd_lin%llr(ilr)%d%nfu3-lzd_lin%llr(ilr)%d%nfl3+1)
          times_convol(iiorb) = real(ii+jj,kind=8)
      end do
      if (nproc>1) then
          call fmpi_allreduce(times_convol, FMPI_SUM, comm=bigdft_mpi%mpi_comm)
      end if

      return !###############################################3

       phi = f_malloc(lnpsidim_orbs,id='phi')
       phi=1.d-5

      allocate(precond_convol_workarrays(lorbs%norbp))
      allocate(precond_workarrays(lorbs%norbp))
      do iorb=1,lorbs%norbp
          iiorb=lorbs%isorb+iorb
          ilr=lorbs%inwhichlocreg(iiorb)
          with_confpot = .true.
          call init_local_work_arrays(lzd_lin%llr(ilr)%d%n1, lzd_lin%llr(ilr)%d%n2, lzd_lin%llr(ilr)%d%n3, &
               lzd_lin%llr(ilr)%d%nfl1, lzd_lin%llr(ilr)%d%nfu1, &
               lzd_lin%llr(ilr)%d%nfl2, lzd_lin%llr(ilr)%d%nfu2, &
               lzd_lin%llr(ilr)%d%nfl3, lzd_lin%llr(ilr)%d%nfu3, &
               with_confpot, precond_convol_workarrays(iorb))
          kx=lorbs%kpts(1,lorbs%iokpt(iorb))
          ky=lorbs%kpts(2,lorbs%iokpt(iorb))
          kz=lorbs%kpts(3,lorbs%iokpt(iorb))
          if (kx**2+ky**2+kz**2 > 0.0_gp .or. lorbs%nspinor==2 ) then
             ncplx=2
          else
             ncplx=1
          end if
          call allocate_work_arrays(lzd_lin%llr(ilr)%geocode, lzd_lin%llr(ilr)%hybrid_on, &
               ncplx, lzd_lin%llr(ilr)%d, precond_workarrays(iorb))
      end do


      call f_zero(times_convol)

      if (nproc>1) then
          call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
      end if
       ist=0
       tt = 0.d0
       do iorb=1,lorbs%norbp
           iiorb = lorbs%isorb + iorb
           ilr = lorbs%inwhichlocreg(iiorb)
           kx=lorbs%kpts(1,lorbs%iokpt(iorb))
           ky=lorbs%kpts(2,lorbs%iokpt(iorb))
           kz=lorbs%kpts(3,lorbs%iokpt(iorb))
           if (kx**2+ky**2+kz**2 > 0.0_gp .or. lorbs%nspinor==2 ) then
              ncplx=2
           else
              ncplx=1
           end if
           do i=1,2*nit+1
               t1 = mpi_wtime()
               call solvePrecondEquation(iproc, nproc, lzd_lin%llr(ilr), ncplx, 6, -0.5d0, &
                    lzd_lin%hgrids(1), lzd_lin%hgrids(2), lzd_lin%hgrids(3), &
                    lorbs%kpts(1,lorbs%iokpt(iorb)), lorbs%kpts(1,lorbs%iokpt(iorb)), lorbs%kpts(1,lorbs%iokpt(iorb)), &
                    phi(1+ist), lzd_lin%llr(ilr)%locregCenter,&
                    1.d-3, 4, precond_convol_workarrays(iorb), precond_workarrays(iorb))
               t2 = mpi_wtime()
               times(i) = t2-t1
           end do
           ! Take the median
           write(20000+iproc,'(a,i7,2es13.2,5x,11es12.2)') 'iiorb, time, tottime, alltimes', iiorb, time, tt, times
           maxtime = maxval(times)
           do i=1,nit
               times(minloc(times,1)) = maxtime
           end do
           do i=1,2*nit
               times(maxloc(times,1)) = 0.d0
           end do
           time = maxval(times)
           tt = tt + time
           write(20000+iproc,'(a,i7,2es13.2,5x,11es12.2)') 'iiorb, time, tottime, alltimes', iiorb, time, tt, times
           times_convol(iiorb) = time
           ist = ist + (lzd_lin%llr(ilr)%wfd%nvctr_c+7*lzd_lin%llr(ilr)%wfd%nvctr_f)*ncplx
       end do
       write(20000+iproc,'(a)') '==========================='

       do iorb=1,lorbs%norbp
           iiorb=lorbs%isorb+iorb
           ilr=lorbs%inwhichlocreg(iiorb)
           call deallocate_workarrays_quartic_convolutions(precond_convol_workarrays(iorb))
           kx=lorbs%kpts(1,lorbs%iokpt(iorb))
           ky=lorbs%kpts(2,lorbs%iokpt(iorb))
           kz=lorbs%kpts(3,lorbs%iokpt(iorb))
           if (kx**2+ky**2+kz**2 > 0.0_gp .or. lorbs%nspinor==2 ) then
              ncplx=2
           else
              ncplx=1
           end if
           call deallocate_work_arrays(lzd_lin%llr(ilr)%geocode, lzd_lin%llr(ilr)%hybrid_on, &
                ncplx, precond_workarrays(iorb))
       end do
       deallocate(precond_convol_workarrays)
       deallocate(precond_workarrays)

       call f_free(phi)

       if (nproc>1) then
           call fmpi_allreduce(times_convol, FMPI_SUM, comm=bigdft_mpi%mpi_comm)
       end if

     end subroutine test_preconditioning

     !!subroutine optimize_loadbalancing()
     !!  implicit none

     !!  ! Local variables
     !!  integer(kind=8) :: isize, isize_ideal

     !!  ! Sum up the total size of all support functions
     !!  isize = int(lnpsidim_orbs,kind=8)
     !!  if (nproc>1) then
     !!      call fmpi_allreduce(isize, 1, FMPI_SUM, bigdft_mpi%mpi_comm)
     !!  end if

     !!  ! Ideal size per task (integer division)
     !!  isize_ideal = isize/int(nproc,kind=8)

     !!  ! Redistribute the support functions such that the load balancing is optimal
     !!  call redistribute(lorbs%norb, isize_ideal, norb_par)

     !!  ! The same for the up and down orbitals
     !!  isize_ideal = isize_ideal/int(in%nspin,kind=8)
     !!  call redistribute(lorbs%norbu, isize_ideal, norbu_par)
     !!  call redistribute(lorbs%norbd, isize_ideal, norbd_par)

     !!end subroutine optimize_loadbalancing

     !!subroutine redistribute(norb, isize_ideal, norb_par)
     !!  implicit none
     !!  integer,intent(in) :: norb
     !!  integer(kind=8),intent(in) :: isize_ideal
     !!  integer,dimension(0:nproc-1),intent(out) :: norb_par
     !!  integer :: jjorbtot, jjorb, jproc, jlr, jorb
     !!  integer(kind=8) :: ncount

     !!  call f_zero(norb_par)
     !!  ncount = int(0,kind=8)
     !!  jproc = 0
     !!  jjorb = 0
     !!  jjorbtot = 0
     !!  do jorb=1,norb
     !!      jjorb = jjorb + 1
     !!      jlr = lorbs%inwhichlocreg(jorb)
     !!      ncount = ncount + int(lzd_lin%llr(jlr)%wfd%nvctr_c+7*lzd_lin%llr(jlr)%wfd%nvctr_f,kind=8)
     !!      if (ncount>=isize_ideal*int(jproc+1,kind=8)) then
     !!          norb_par(jproc) = jjorb
     !!          jjorbtot = jjorbtot + jjorb
     !!          jjorb = 0
     !!          jproc = jproc + 1
     !!      end if
     !!      if (jproc==nproc-1) exit
     !!  end do
     !!  norb_par(nproc-1) = (norb - jjorbtot) !+jjorb -> why was this needed?! !take the rest
     !!  do jproc=0,nproc-1
     !!      !if (iproc==0) write(*,*) 'jproc, norb_par(jproc)', jproc, norb_par(jproc)
     !!  end do
     !!end subroutine redistribute


     subroutine optimize_loadbalancing2()
       use sparsematrix_init, only: redistribute
       implicit none

       ! Local variables
       real(kind=8) :: time, time_ideal

       ! Sum up the total size of all support functions
       time = sum(times_convol)

       ! Ideal size per task (integer division)
       time_ideal = time/int(nproc,kind=8)

       ! Redistribute the support functions such that the load balancing is optimal
       !!call redistribute2(lorbs%norb, time_ideal, norb_par)
       call redistribute(nproc, lorbs%norb, times_convol, time_ideal, norb_par)

       ! The same for the up and down orbitals
       time_ideal = time_ideal/int(in%nspin,kind=8)
       !!call redistribute2(lorbs%norbu, time_ideal, norbu_par)
       !!call redistribute2(lorbs%norbd, time_ideal, norbd_par)
       call redistribute(nproc, lorbs%norbu, times_convol, time_ideal, norbu_par)
       call redistribute(nproc, lorbs%norbd, times_convol, time_ideal, norbd_par)

     end subroutine optimize_loadbalancing2

     !!subroutine redistribute2(norb, time_ideal, norb_par)
     !!  implicit none
     !!  integer,intent(in) :: norb
     !!  real(kind=8),intent(in) :: time_ideal
     !!  integer,dimension(0:nproc-1),intent(out) :: norb_par
     !!  integer :: jjorbtot, jjorb, jproc, jlr, jorb
     !!  real(kind=8) :: tcount
       !integer :: jlr

     !!  call f_zero(norb_par)
     !!  if (norb>=nproc) then
     !!      tcount = 0.d0
     !!      jproc = 0
     !!      jjorb = 0
     !!      jjorbtot = 0
     !!      do jorb=1,norb
     !!          if (jproc==nproc-1) exit
     !!          jjorb = jjorb + 1
     !!          if(jorb==norb) exit !just to besure that no out of bound happens
     !!          tcount = tcount + times_convol(jorb)
     !!          !if (iproc==0) write(*,'(a,2i8,2es14.5)') 'jorb, jproc, tcount, diff to target', jorb, jproc, tcount, abs(tcount-time_ideal*real(jproc+1,kind=8))
     !!          if (abs(tcount-time_ideal*real(jproc+1,kind=8)) <= &
     !!                  abs(tcount+times_convol(jorb+1)-time_ideal*real(jproc+1,kind=8))) then
     !!              norb_par(jproc) = jjorb
     !!              jjorbtot = jjorbtot + jjorb
     !!              jjorb = 0
     !!              jproc = jproc + 1
     !!          end if
     !!      end do
     !!      norb_par(nproc-1) = jjorb + (norb - jjorbtot) !take the rest
     !!      !do jproc=0,nproc-1
     !!      !    if (iproc==0) write(*,*) 'jproc, norb_par(jproc)', jproc, norb_par(jproc)
     !!      !end do
     !!  else
     !!      ! Equal distribution
     !!      norb_par(0:norb-1) = 1
     !!  end if
     !!end subroutine redistribute2


     subroutine fragment_stuff()
       use module_interfaces, only: input_check_psi_id
       implicit none
       integer :: iorbn
       integer, allocatable, dimension(:) :: inwhichlocreg_tmp, onwhichatom_tmp

       frag_allocated=.false.
       if (inputpsi == 'INPUT_PSI_DISK_LINEAR' .or. in%lin%fragment_calculation) then
          allocate(ref_frags(in%frag%nfrag_ref))
          do ifrag=1,in%frag%nfrag_ref
             ref_frags(ifrag)=fragment_null()
          end do
          call init_fragments(in,lorbs,atoms%astruct,ref_frags)
         frag_allocated=.true.
       else
          nullify(ref_frags)
         frag_allocated=.false.
       end if

       call input_check_psi_id(inputpsi, input_wf_format, in%dir_output, &
            orbs, lorbs, iproc, nproc, in%frag%nfrag_ref, in%frag%dirname, ref_frags)

       ! we need to deallocate the fragment arrays we just allocated as not a restart calculation so this is no longer needed
       if (frag_allocated .and. (.not. in%lin%fragment_calculation) .and. inputpsi /= 'INPUT_PSI_DISK_LINEAR') then
           do ifrag=1,in%frag%nfrag_ref
              call fragment_free(ref_frags(ifrag))
              ref_frags(ifrag)%astruct_frg%nat=-1
              ref_frags(ifrag)%fbasis%forbs=minimal_orbitals_data_null()
              !ref_frags(ifrag)=fragment_null()
           end do
          deallocate(ref_frags)
         frag_allocated=.false.
       end if

       ! not sure if this will mess up load balancing but we want the tmbs to stay in input file order
       ! for reading for disk this is taken care of automatically
       if (in%lin%fragment_calculation .and. inputpsi /= 'INPUT_PSI_DISK_LINEAR') then
          inwhichlocreg_tmp = f_malloc(lorbs%norb,id='inwhichlocreg_tmp')
          onwhichatom_tmp = f_malloc(lorbs%norb,id='onwhichatom_tmp') 

          call vcopy(lorbs%norb,lorbs%inwhichlocreg(1),1,inwhichlocreg_tmp(1),1)
          call vcopy(lorbs%norb,lorbs%onwhichatom(1),1,onwhichatom_tmp(1),1)

          !!print*,lorbs%norb
          !!print*,lorbs%inwhichlocreg
          !!print*,lorbs%onwhichatom

          ! might be a better way of doing this...
          ! assume that norbu==norbd etc, check if this is always true...
          iorbn=0
          do iat=1,atoms%astruct%nat
             do iorb=1,lorbs%norbu
                if (onwhichatom_tmp(iorb)/=iat) cycle
                iorbn = iorbn+1
                lorbs%onwhichatom(iorbn) = iat
                lorbs%inwhichlocreg(iorbn) = inwhichlocreg_tmp(iorb)
                ! perform the same change for spin down
                if (lorbs%nspin==2) then
                   lorbs%onwhichatom(iorbn+lorbs%norbu) = iat
                   lorbs%inwhichlocreg(iorbn+lorbs%norbu) = inwhichlocreg_tmp(iorb+lorbs%norbu)
                end if
                ! double check - might not always be true in future
                if (iorbn /= inwhichlocreg_tmp(iorb)) then
                   write(*,*) 'Error reordering tmbs for fragment calculation (1)',iorbn,inwhichlocreg_tmp(iorb)
                   stop
                end if
             end do
          end do
          if (iorbn /= lorbs%norbu) then
             write(*,*) 'Error reordering tmbs for fragment calculation (2)',iorbn,lorbs%norb
             stop
          end if


          call f_free(inwhichlocreg_tmp)
          call f_free(onwhichatom_tmp)
       end if

     end subroutine fragment_stuff

END SUBROUTINE system_initialization


subroutine system_initKernels(verb, iproc, nproc, geocode, in, denspot)
  use module_types
  use module_xc
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use module_base
  implicit none
  logical, intent(in) :: verb
  integer, intent(in) :: iproc, nproc
  character, intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  type(input_variables), intent(in) :: in
  type(DFT_local_fields), intent(inout) :: denspot

  integer, parameter :: ndegree_ip = 16

  denspot%pkernel=pkernel_init(iproc,nproc,in%PS_dict,&
       geocode,denspot%dpbox%mesh%ndims,denspot%dpbox%mesh%hgrids,&
       mpi_env=denspot%dpbox%mpi_env)

  !create the sequential kernel if the exctX parallelisation scheme requires it
  if ((xc_exctXfac(denspot%xc) /= 0.0_gp .and. in%exctxpar=='OP2P' .or. in%SIC%alpha /= 0.0_gp)&
       .and. denspot%dpbox%mpi_env%nproc > 1) then
     !the communicator of this kernel is bigdft_mpi%mpi_comm
     !this might pose problems when using SIC or exact exchange with taskgroups
     denspot%pkernelseq=pkernel_init(0,1,in%PS_dict_seq,&
          geocode,denspot%dpbox%mesh%ndims,denspot%dpbox%mesh%hgrids)
  else 
     denspot%pkernelseq = denspot%pkernel
  end if

END SUBROUTINE system_initKernels

subroutine system_createKernels(denspot, verb)
  use module_base
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp
  implicit none
  logical, intent(in) :: verb
  type(DFT_local_fields), intent(inout) :: denspot
  call pkernel_set(denspot%pkernel,verbose=verb)
    !create the sequential kernel if pkernelseq is not pkernel
  if (denspot%pkernelseq%mpi_env%nproc == 1 .and. denspot%pkernel%mpi_env%nproc /= 1) then
     call pkernel_set(denspot%pkernelseq,verbose=.false.)
  else
     denspot%pkernelseq = denspot%pkernel
  end if

END SUBROUTINE system_createKernels

!> calculate the dielectric function for the cavitBy
subroutine epsilon_cavity(atoms,rxyz,pkernel)
  use dynamic_memory
  use Poisson_Solver
  use module_atoms
  use ao_inguess, only: atomic_info
  use yaml_output
  use numerics, only : Bohr_Ang
  use module_base, only: bigdft_mpi,f_zero
  use module_defs, only: UNINITIALIZED
  use f_enums, f_str => str
  use yaml_output
  use dictionaries, only: f_err_throw
  use box
  use bounds, only: locreg_mesh_origin
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  type(coulomb_operator), intent(inout) :: pkernel
  !local variables
  real(gp), parameter :: fact_Pau=1.36d0 ! Multiplying factor to enlarge the rigid cavity, Pauling radii.
  real(gp), parameter :: fact_Bondi=1.32d0 ! Multiplying factor to enlarge the rigid cavity, Bondi radii.
  real(gp), parameter :: fact_UFF=1.12d0 ! Multiplying factor to enlarge the rigid cavity, UFF radii.
  integer :: i,iat
  integer :: i1,i2,i3,unt,i3s,i23
  real(gp) :: IntSur,IntVol,noeleene,Cavene,Repene,Disene,IntSurt,IntVolt,diffSur,diffVol
  real(gp) :: radii_Pau,radii_Bondi,radii_UFF,fact,rcav
  type(atoms_iterator) :: it
  real(gp), dimension(:), allocatable :: radii
  real(gp), dimension(:,:), allocatable :: rxyz_shifted
  real(gp), dimension(:,:,:), allocatable :: eps,oneoeps,oneosqrteps,corr
  real(gp), dimension(:,:,:,:), allocatable :: dlogeps
!!$  real(gp), dimension(:,:,:), allocatable :: epst,oneoepst,oneosqrtepst,corrt
  real(gp), dimension(:,:,:,:), allocatable :: dlogepst
  !type(cell) :: mesh
  real(dp), dimension(3) :: origin
  integer, parameter :: radii_cav = 3 ! 1 for Pauling, 2 for Bondi, 3 for UFF.
  character(2) :: atname

  !set the vdW radii for the cavity definition
  !iterate above atoms

  radii=f_malloc(atoms%astruct%nat,id='radii')
  !radii_nofact=f_malloc(atoms%astruct%nat,id='radii_nofact')
  eps=f_malloc(pkernel%mesh%ndims,id='eps')
  dlogeps=f_malloc([3,pkernel%mesh%ndims(1),pkernel%mesh%ndims(2),pkernel%mesh%ndims(3)],id='dlogeps')
  oneoeps=f_malloc(pkernel%mesh%ndims,id='oneoeps')
  oneosqrteps=f_malloc(pkernel%mesh%ndims,id='oneosqrteps')
  corr=f_malloc(pkernel%mesh%ndims,id='corr')

  it=atoms_iter(atoms%astruct)
  !python metod
  if (bigdft_mpi%iproc==0) call yaml_mapping_open('Covalent radii',flow=.true.)
  do while(atoms_iter_next(it))
     !only amu is extracted here
     call atomic_info(atoms%nzatom(it%ityp),atoms%nelpsp(it%ityp),&
          rcov=radii(it%iat))
     call astruct_at_from_dict(it%attrs,cavity_radius=rcav)
     if (rcav==UNINITIALIZED(rcav)) then
        radii(it%iat)=pkernel_get_radius(pkernel,it%name)
     else
        radii(it%iat)=rcav
     end if
     if (bigdft_mpi%iproc==0) call yaml_map(it%name,radii(it%iat))
  end do
  if (bigdft_mpi%iproc==0) call yaml_mapping_close()

!  if (bigdft_mpi%iproc==0) call yaml_map('Covalent radii',radii)

  fact=pkernel%cavity%fact_rigid
  it=atoms_iter(atoms%astruct)
  do while(atoms_iter_next(it))
     radii(it%iat)=fact*radii(it%iat)/Bohr_Ang
  end do

  !here the pkernel_set_epsilon routine should been modified to accept
  !already the radii and the atoms

  !mesh=cell_new(atoms%astruct%geocode,pkernel%mesh%ndims,pkernel%mesh%hgrids)
  origin=locreg_mesh_origin(pkernel%mesh)
  rxyz_shifted=f_malloc([3,atoms%astruct%nat],id='rxyz_shifted')
  do iat=1,atoms%astruct%nat
     rxyz_shifted(:,iat)=rxyz(:,iat)+origin
  end do
  call pkernel_set_epsilon(pkernel,nat=atoms%astruct%nat,rxyz=rxyz_shifted,radii=radii)
  call f_free(rxyz_shifted)
  call f_free(radii)
  call f_free(eps)
  call f_free(dlogeps)
  call f_free(oneoeps)
  call f_free(oneosqrteps)
  call f_free(corr)
end subroutine epsilon_cavity

!> Check the difference of two 3 dimensional vectors
subroutine check_accuracy_3d(n01,n02,n03,i,r1,r2)
  use yaml_output
  use dynamic_memory
  use f_utils
  use module_base, only: bigdft_mpi
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: i
  real(kind=8), dimension(n01,n02,n03), intent(in) :: r1,r2
  !automatic array, to be check is stack poses problem
  real(kind=8), dimension(:,:,:), allocatable :: re
  integer :: i1,i2,i3,j,i1_max,i2_max,i3_max,jj,unt
  real(kind=8) :: max_val,fact
  character(len=20) :: str

  re=f_malloc([n01,n02,n03],id='re')

      max_val = 0.d0
      i1_max = 1
      i2_max = 1
      i3_max = 1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               re(i1,i2,i3) = r1(i1,i2,i3) - r2(i1,i2,i3)
               fact=abs(re(i1,i2,i3))
               if (max_val < fact) then
                  max_val = fact
                  i1_max = i1
                  i2_max = i2
                  i3_max = i3
               end if
            end do
         end do
      end do
  if (bigdft_mpi%iproc==0) then
      if (max_val == 0.d0) then
         call yaml_map('Inf. Norm difference with reference',0.d0)
      else
         call yaml_mapping_open('Inf. Norm difference with reference')
         call yaml_map('Value',max_val,fmt='(1pe22.15)')
         call yaml_map('Point',[i1_max,i2_max,i3_max],fmt='(i4)')
         call yaml_map('Some values',[re(n01/2,n02/2,n03/2),re(2,n02/2,n03/2),re(10,n02/2,n03/2)],&
              fmt='(1pe22.15)')
         call yaml_mapping_close()
      end if
  end if
  call f_free(re)

end subroutine check_accuracy_3d

!> Check the difference of two 4 dimensional vectors
subroutine check_accuracy_4d(n01,n02,n03,i,r1,r2)
  use yaml_output
  use dynamic_memory
  use f_utils
  use module_base, only: bigdft_mpi
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: i
  real(kind=8), dimension(3,n01,n02,n03), intent(in) :: r1,r2
  !automatic array, to be check is stack poses problem
  real(kind=8), dimension(:,:,:,:), allocatable :: re
  integer :: i1,i2,i3,j,i1_max,i2_max,i3_max,jj,unt
  real(kind=8) :: max_val,fact
  character(len=20) :: str
  
  re=f_malloc([3,n01,n02,n03],id='re')

      max_val = 0.d0
      i1_max = 1
      i2_max = 1
      i3_max = 1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               do j=1,3
               re(j,i1,i2,i3) = r1(j,i1,i2,i3) - r2(j,i1,i2,i3)
               fact=abs(re(j,i1,i2,i3))
               if (max_val < fact) then
                  max_val = fact
                  i1_max = i1
                  i2_max = i2
                  i3_max = i3
               end if
               end do
            end do
         end do
      end do
  if (bigdft_mpi%iproc==0) then
      if (max_val == 0.d0) then
         call yaml_map('Inf. Norm difference with reference',0.d0)
      else
         call yaml_mapping_open('Inf. Norm difference with reference')
         call yaml_map('Value',max_val,fmt='(1pe22.15)')
         call yaml_map('Point',[i1_max,i2_max,i3_max],fmt='(i4)')
         call yaml_map('Some values',[re(1,n01/2,n02/2,n03/2),re(1,2,n02/2,n03/2),re(1,10,n02/2,n03/2)],&
              fmt='(1pe22.15)')
         call yaml_mapping_close()
      end if
  end if
  call f_free(re)

end subroutine check_accuracy_4d

!> Calculate the inner cavity for a sccs run to avoit discontinuity in epsilon
!! due to near-zero edens near atoms
subroutine epsinnersccs_cavity(atoms,rxyz,pkernel)
  use dynamic_memory
  use Poisson_Solver
  use module_atoms
  use ao_inguess, only: atomic_info
  !use yaml_output
  use numerics, only : Bohr_Ang
  use module_base, only: bigdft_mpi
  use f_enums, f_str => str
  use yaml_output
  use dictionaries, only: f_err_throw
  use PStypes, only: epsilon_inner_cavity
  use box
  use bounds, only: locreg_mesh_origin
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  type(coulomb_operator), intent(inout) :: pkernel

  !local variables
  integer :: i,n1,n23,i3s
  real(gp) :: delta
  !type(atoms_iterator) :: it
  real(gp), dimension(:), allocatable :: radii
  real(gp), dimension(:,:,:), allocatable :: eps

  radii=f_malloc(atoms%astruct%nat,id='radii')

  delta=2.0*maxval(pkernel%mesh%hgrids)
  do i=1,atoms%astruct%nat
   radii(i) = 0.5d0/Bohr_Ang
  end do
  call epsilon_inner_cavity(pkernel,atoms%astruct%nat,rxyz,radii,delta,locreg_mesh_origin(pkernel%mesh))

  call f_free(radii)
end subroutine epsinnersccs_cavity

!> Calculate the important objects related to the physical properties of the system
subroutine system_properties(iproc,nproc,in,atoms,orbs)!,radii_cf)
  use module_base
  use module_types
  use module_interfaces, only: orbitals_descriptors
  use public_enums
  implicit none
  integer, intent(in) :: iproc,nproc
  type(input_variables), intent(in) :: in
  type(atoms_data), intent(in) :: atoms
  type(orbitals_data), intent(inout) :: orbs
!  real(gp), dimension(atoms%astruct%ntypes,3), intent(out) :: radii_cf
  !local variables
  !n(c) character(len=*), parameter :: subname='system_properties'
  integer :: nspinor

!  radii_cf = atoms%radii_cf
 !!$  call read_radii_variables(atoms, radii_cf, in%crmult, in%frmult, in%projrad)
 !!$  call read_atomic_variables(atoms, trim(in%file_igpop),in%nspin)
  if (iproc == 0) call print_atomic_variables(atoms, max(in%hx,in%hy,in%hz), in%ixc, in%dispersion)
  if(in%nspin==4) then
     nspinor=4
  else
     nspinor=1
  end if
  call orbitals_descriptors(iproc, nproc,in%gen_norb,in%gen_norbu,in%gen_norbd,in%nspin,nspinor,&
       in%gen_nkpt,in%gen_kpt,in%gen_wkpt,orbs,LINEAR_PARTITION_NONE)
  orbs%occup(1:orbs%norb*orbs%nkpts) = in%gen_occup
  if (iproc==0) call print_orbitals(orbs, atoms%astruct%geocode)

  !Check if orbitals and electrons
  if (orbs%norb*orbs%nkpts == 0) &
     & call f_err_throw('No electrons in the system. Check your input variables or atomic positions', &
     & err_id=BIGDFT_INPUT_VARIABLES_ERROR)

END SUBROUTINE system_properties


!> Check for the need of a core density and fill the rhocore array which
!! should be passed at the rhocore pointer
subroutine calculate_rhocore(at,rxyz,dpbox,rhocore)
  use module_base
  use module_types
  use module_dpbox, only: denspot_distribution
  use public_enums, only: PSPCODE_PAW
  use m_pawrad,  only : pawrad_type, pawrad_init, pawrad_free
  use yaml_output
  implicit none
  type(atoms_data), intent(in) :: at
  type(denspot_distribution), intent(in) :: dpbox
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(wp), dimension(:,:,:,:), pointer :: rhocore
  !local variables
  character(len=*), parameter :: subname='calculate_rhocore'
  integer :: ityp,iat,j3,i1,i2,ncmax,i !,ierr,ind
  real(wp) :: tt
  real(gp) :: rx,ry,rz,rloc,cutoff,chgat
  real(gp), dimension(:), allocatable :: chg_at
  logical :: donlcc
  type(pawrad_type)::core_mesh
  !allocatable arrays
  integer,allocatable :: ifftsph(:)
  real(gp),allocatable:: rr(:), raux(:), rcart(:,:)

  !check for the need of a nonlinear core correction
  donlcc=.false.
  chk_nlcc: do ityp=1,at%astruct%ntypes
     if (at%nlcc_ngv(ityp) /= UNINITIALIZED(1) .or. &
          & at%nlcc_ngc(ityp) /= UNINITIALIZED(1) .or. &
          & at%npspcode(ityp) == PSPCODE_PAW) then
        donlcc = .true.
        exit
     end if
  end do chk_nlcc
  
  if (.not. donlcc) then
     !No NLCC needed, nullify the pointer 
     nullify(rhocore)
     return
  end if

  !allocate pointer rhocore
     chg_at=f_malloc0(at%astruct%ntypes,id='chg_at')
  rhocore = f_malloc0_ptr((/ dpbox%mesh%ndims(1) , dpbox%mesh%ndims(2) ,max(dpbox%n3d,1) , 10 /), id = 'rhocore')
  !perform the loop on any of the atoms which have this feature
  do ityp = 1, at%astruct%ntypes
     if (at%npspcode(ityp) == PSPCODE_PAW) then
        !  Create mesh_core object
        !  since core_mesh_size can be bigger than pawrad%mesh_size, 
        call pawrad_init(core_mesh, mesh_size = at%pawtab(ityp)%core_mesh_size, &
             & mesh_type = at%pawrad(ityp)%mesh_type, &
             & rstep = at%pawrad(ityp)%rstep, &
             & lstep = at%pawrad(ityp)%lstep)

        !  Set radius size:
        do i = 1, size(at%pawtab(ityp)%tcoredens, 1) - 1
           if (core_mesh%rad(i) > at%pawtab(ityp)%rpaw .and. &
                & at%pawtab(ityp)%tcoredens(i, 1) < 1e-10) exit
        end do
        rloc = core_mesh%rad(i)
        cutoff = rloc * 1.1d0

        !  allocate arrays
        if (dpbox%n3p > 0) then
           ! ncmax=1+int(1.1_dp*nfft*four_pi/(three*ucvol)*rshp**3)
           ! ncmax=1+int(1.1d0*((rshp/dpbox%hgrids(1))*(rshp/dpbox%hgrids(2))*pi_param))
           ! Calculation is done per z plane.
           ncmax = 1 + int(4._gp * pi_param * rloc ** 2 / dpbox%mesh%hgrids(1) / dpbox%mesh%hgrids(2))
        else
           ncmax = 1
        end if
        ifftsph = f_malloc(ncmax, id = "ifftsph_tmp")
        rr = f_malloc(ncmax, id = "rr")
        raux = f_malloc(ncmax, id = "raux")
        rcart = f_malloc((/ 3, ncmax /), id = "rcart")
     else
        rloc = at%psppar(0,0,ityp)
        cutoff = 10.d0 * rloc
     end if

     do iat = 1, at%astruct%nat
        if (at%astruct%iatype(iat) /= ityp) cycle

        rx=rxyz(1,iat) 
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)
        
        if (at%npspcode(ityp) == PSPCODE_PAW) then
           call mkcore_paw_iat(bigdft_mpi%iproc,at,ityp,rx,ry,rz,cutoff,&
                & dpbox%mesh%hgrids(1),dpbox%mesh%hgrids(2),dpbox%mesh%hgrids(3), &
                & dpbox%mesh%ndims(1), dpbox%mesh%ndims(2),dpbox%mesh%ndims(3), &
                & dpbox%i3s,dpbox%n3d,core_mesh, rhocore, ncmax, ifftsph, &
                & rr, rcart, raux)
        else
           call calc_rhocore_iat(bigdft_mpi%iproc,at,ityp,rx,ry,rz,cutoff,&
                & dpbox%mesh%hgrids(1),dpbox%mesh%hgrids(2),dpbox%mesh%hgrids(3), &
                & dpbox%mesh%ndims(1), dpbox%mesh%ndims(2),dpbox%mesh%ndims(3), &
                & dpbox%i3s,dpbox%n3d,chgat,rhocore)
           chg_at(ityp)=chg_at(ityp)+chgat
        end if
     end do
     !  Deallocate
     if (at%npspcode(ityp) == PSPCODE_PAW) then
        call pawrad_free(core_mesh)
        call f_free(ifftsph)
        call f_free(rr)
        call f_free(raux)
        call f_free(rcart)
     end if
  end do

  !calculate total core charge in the grid
  !In general this should be really bad

!!$     do j3=1,n3d
!!$        tt=0.0_wp
!!$        do i2=1,d%n2i
!!$           do i1=1,d%n1i
!!$              !ind=i1+(i2-1)*d%n1i+(j3+i3xcsh-1)*d%n1i*d%n2i
!!$              tt=tt+rhocore(i1,i2,j3,1)
!!$           enddo
!!$        enddo
!!$        write(17+iproc,*)j3+i3s-1,tt
!!$     enddo
!!$call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
!!$stop
  tt=0.0_wp
  do j3=1,dpbox%n3p
     do i2=1,dpbox%mesh%ndims(2)
        do i1=1,dpbox%mesh%ndims(1)
           !ind=i1+(i2-1)*d%n1i+(j3+i3xcsh-1)*d%n1i*d%n2i
           tt=tt+rhocore(i1,i2,j3+dpbox%i3xcsh,1)
        enddo
     enddo
  enddo

  if (bigdft_mpi%nproc > 1) call fmpi_allreduce(tt,1,FMPI_SUM,comm=bigdft_mpi%mpi_comm)
  tt=tt*dpbox%mesh%volume_element
  if (bigdft_mpi%iproc == 0) then
     call yaml_mapping_open('Analytic core charges for atom species')
     do ityp=1,at%astruct%ntypes
        if (chg_at(ityp) /= 0.0_gp) &
             call yaml_map(trim(at%astruct%atomnames(ityp)),chg_at(ityp),fmt='(f15.7)')
     end do
     call yaml_mapping_close()
     call yaml_map('Total core charge',sum(chg_at),fmt='(f15.7)')
     call yaml_map('Total core charge on the grid', tt,fmt='(f15.7)', advance = "no")
     call yaml_comment('To be compared with analytic one')
  end if
  call f_free(chg_at)

END SUBROUTINE calculate_rhocore




!> Calculate the number of electrons and check the polarisation (mpol)
subroutine read_n_orbitals(iproc, qelec_up, qelec_down, norbe, &
     & atoms, qcharge, nspin, mpol, norbsempty)
  use module_atoms, only: atoms_data
  use ao_inguess, only: charge_and_spol
  use module_defs, only: gp
  use dictionaries, only: f_err_throw
  use yaml_strings
  use yaml_output, only: yaml_warning, yaml_comment
  use dynamic_memory
  !use ao_inguess, only : count_atomic_shells
  implicit none
  !Arguments
  type(atoms_data), intent(in) :: atoms
  integer, intent(out) :: norbe
  real(gp), intent(out) :: qelec_up, qelec_down
  real(gp), intent(in) :: qcharge
  integer, intent(in) :: nspin, mpol, norbsempty, iproc
  !Local variables
  logical :: int_charge
  integer :: nel, nel_up,nel_dwn,nchg,iat, ityp, ispinsum, ichgsum, ichg, ispol,iabspol!, nspinor
  real(gp) :: qelec

  call f_routine(id='read_n_orbitals')

  !calculate number of electrons and orbitals
  ! Number of electrons and number of semicore atoms
  qelec=0
  do iat=1,atoms%astruct%nat
     ityp=atoms%astruct%iatype(iat)
     qelec=qelec+real(atoms%nelpsp(ityp),gp)
  enddo
  qelec=qelec-qcharge
  !roundoff of the charge
  nel=nint(qelec)
  if (qelec - real(nel,gp) > 1.e-12_gp) nel=nel+1
  nchg=nint(-qcharge)
  if (-qcharge - real(nchg,gp) > 1.e-12_gp) nchg=nchg+1
  nchg=-nchg

  int_charge = real(nint(qelec),gp) == qelec

  if(qelec < 0.0_gp ) then
    call f_err_throw('Number of electrons is negative:' // trim(yaml_toa(qelec)) // &
      & '. FIX: decrease value of qcharge.', err_name='BIGDFT_RUNTIME_ERROR')
  end if

  ! Number of orbitals
  if (nspin==1) then
     qelec_up=qelec
     qelec_down=0.0_gp
  else if(nspin==4) then
     qelec_up=qelec
     qelec_down=0.0_gp
  else 
     if (mod(nel+mpol,2) /=0 .and. int_charge) then
          call f_err_throw('The mpol polarization should have the same parity of the (rounded) number of electrons. ' // &
            & '(mpol='+mpol+' and qelec='+qelec+')', &
            & err_name='BIGDFT_INPUT_VARIABLES_ERROR')

     end if
     !put the charge according to the polarization.
     !non-integer part always goes to the upper spin shell
     !nelec_up=min((nelec+mpol)/2,nelec) !this is the integer part (rounded)
     nel_up=min((nel+mpol)/2,nel)
     nel_dwn=nel-nel_up
     qelec_down=real(nel_dwn,gp)
     !then the elec_up part is redefined with the actual charge
     qelec_up=qelec-qelec_down

     !test if the spin is compatible with the input guess polarisations
     ispinsum=0
     ichgsum=0
     iabspol=0
     do iat=1,atoms%astruct%nat
        call charge_and_spol(atoms%astruct%input_polarization(iat),ichg,ispol)
        ispinsum=ispinsum+ispol
        ichgsum=ichgsum+ichg
        iabspol=iabspol+abs(ispol)
     end do

     if (ispinsum /= nel_up-nel_dwn .and. int_charge) then
        call f_err_throw('Total polarisation for the input guess (found ' // &
             trim(yaml_toa(ispinsum)) // &
             ') must be equal to rounded nel_up-nel_dwn ' // &
             '(nelec=' // trim(yaml_toa(qelec)) // ', mpol=' // trim(yaml_toa(mpol)) // &
             ', nel_up-nel_dwn=' // trim((yaml_toa(nel_up-nel_dwn))) // &
             ', nel_up=' // trim((yaml_toa(nel_up))) // &
             ', nel_dwn=' // trim((yaml_toa(nel_dwn))) // &
             '). Use the keyword "IGSpin" or add a spin component for the input guess per atom.', &
             err_name='BIGDFT_INPUT_VARIABLES_ERROR')
     end if

     if (ichgsum /= nchg .and. ichgsum /= 0) then
        call f_err_throw('Total input charge (found ' // trim(yaml_toa(ichgsum)) // &
             & ') cannot be different than rounded charge. With charge =' // trim(yaml_toa(qcharge)) // &
             & ' and input charge=' // trim(yaml_toa(ichgsum)), &
             & err_name='BIGDFT_INPUT_VARIABLES_ERROR')
     end if

     !now warn if there is no input guess spin polarisation
!!$     ispinsum=0
!!$     do iat=1,atoms%astruct%nat
!!$        call charge_and_spol(atoms%astruct%input_polarization(iat),ichg,ispol)
!!$        ispinsum=ispinsum+abs(ispol)
!!$     end do
!!$     if (ispinsum == 0) then
     if (iabspol == 0 .and. iproc==0 .and. norbsempty == 0) &
          call yaml_warning('Found no input polarisation, add it for a correct input guess')
  end if

  norbe = 0
  do iat=1,atoms%astruct%nat
     norbe=norbe+atoms%aoig(iat)%nao
  end do
  if (norbe == 0) norbe = nel !elec_up + nelec_down ! electron gas case

  call f_release_routine()

end subroutine read_n_orbitals


!> Find the correct position of the nlcc parameters
subroutine nlcc_start_position(ityp,atoms,ngv,ngc,islcc)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ityp
  type(atoms_data), intent(in) :: atoms
  integer, intent(out) :: ngv,ngc,islcc
  !local variables
  integer :: ilcc,jtyp

  ilcc=0
  do jtyp=1,ityp-1
     ngv=atoms%nlcc_ngv(jtyp)
     if (ngv /= UNINITIALIZED(ngv)) ilcc=ilcc+(ngv*(ngv+1)/2)
     ngc=atoms%nlcc_ngc(jtyp)
     if (ngc /= UNINITIALIZED(ngc)) ilcc=ilcc+(ngc*(ngc+1))/2
  end do
  islcc=ilcc

  ngv=atoms%nlcc_ngv(ityp)
  if (ngv==UNINITIALIZED(1)) ngv=0
  ngc=atoms%nlcc_ngc(ityp)
  if (ngc==UNINITIALIZED(1)) ngc=0
END SUBROUTINE nlcc_start_position


!!!!> Define the descriptors of the orbitals from a given norb
!!!!! It uses the cubic strategy for partitioning the orbitals
!!!subroutine orbitals_descriptors_forLinear(iproc,nproc,norb,norbu,norbd,nspin,nspinor,nkpt,kpt,wkpt,orbs)
!!!  use module_base
!!!  use module_types
!!!  implicit none
!!!  integer, intent(in) :: iproc,nproc,norb,norbu,norbd,nkpt,nspin
!!!  integer, intent(in) :: nspinor
!!!  type(orbitals_data), intent(out) :: orbs
!!!  real(gp), dimension(nkpt), intent(in) :: wkpt
!!!  real(gp), dimension(3,nkpt), intent(in) :: kpt
!!!  !local variables
!!!  character(len=*), parameter :: subname='orbitals_descriptors'
!!!  integer :: iorb,jproc,norb_tot,ikpt,i_stat,jorb,ierr,i_all,iiorb
!!!  integer :: mpiflag
!!!  logical, dimension(:), allocatable :: GPU_for_orbs
!!!  integer, dimension(:,:), allocatable :: norb_par !(with k-pts)
!!!
!!!
!!!  allocate(orbs%norb_par(0:nproc-1,0:nkpt),stat=i_stat)
!!!  call memocc(i_stat,orbs%norb_par,'orbs%norb_par',subname)
!!!
!!!  !assign the value of the k-points
!!!  orbs%nkpts=nkpt
!!!  !allocate vectors related to k-points
!!!  allocate(orbs%kpts(3,orbs%nkpts),stat=i_stat)
!!!  call memocc(i_stat,orbs%kpts,'orbs%kpts',subname)
!!!  allocate(orbs%kwgts(orbs%nkpts),stat=i_stat)
!!!  call memocc(i_stat,orbs%kwgts,'orbs%kwgts',subname)
!!!  orbs%kpts(:,1:nkpt) = kpt(:,:)
!!!  orbs%kwgts(1:nkpt) = wkpt(:)
!!!
!!!  ! Change the wavefunctions to complex if k-points are used (except gamma).
!!!  orbs%nspinor=nspinor
!!!  if (nspinor == 1) then
!!!     if (maxval(abs(orbs%kpts)) > 0._gp) orbs%nspinor=2
!!!     !nspinor=2 !fake, used for testing with gamma
!!!  end if
!!!  orbs%nspin = nspin
!!!
!!!  !initialise the array
!!!  do jproc=0,nproc-1
!!!     orbs%norb_par(jproc,0)=0 !size 0 nproc-1
!!!  end do
!!!
!!!  !create an array which indicate which processor has a GPU associated 
!!!  !from the viewpoint of the BLAS routines (deprecated, not used anymore)
!!!  if (.not. GPUshare) then
!!!     allocate(GPU_for_orbs(0:nproc-1),stat=i_stat)
!!!     call memocc(i_stat,GPU_for_orbs,'GPU_for_orbs',subname)
!!!     
!!!     if (nproc > 1) then
!!!        call MPI_ALLGATHER(GPUconv,1,MPI_LOGICAL,GPU_for_orbs(0),1,MPI_LOGICAL,&
!!!             bigdft_mpi%mpi_comm,ierr)
!!!     else
!!!        GPU_for_orbs(0)=GPUconv
!!!     end if
!!!     
!!!     i_all=-product(shape(GPU_for_orbs))*kind(GPU_for_orbs)
!!!     deallocate(GPU_for_orbs,stat=i_stat)
!!!     call memocc(i_stat,i_all,'GPU_for_orbs',subname)
!!!  end if
!!!
!!!  allocate(norb_par(0:nproc-1,orbs%nkpts),stat=i_stat)
!!!  call memocc(i_stat,norb_par,'norb_par',subname)
!!!
!!!  !old system for calculating k-point repartition
!!!!!$  call parallel_repartition_with_kpoints(nproc,orbs%nkpts,norb,orbs%norb_par)
!!!!!$
!!!!!$  !check the distribution
!!!!!$  norb_tot=0
!!!!!$  do jproc=0,iproc-1
!!!!!$     norb_tot=norb_tot+orbs%norb_par(jproc)
!!!!!$  end do
!!!!!$  !reference orbital for process
!!!!!$  orbs%isorb=norb_tot
!!!!!$  do jproc=iproc,nproc-1
!!!!!$     norb_tot=norb_tot+orbs%norb_par(jproc)
!!!!!$  end do
!!!!!$
!!!!!$  if(norb_tot /= norb*orbs%nkpts) then
!!!!!$     write(*,*)'ERROR: partition of orbitals incorrect, report bug.'
!!!!!$     write(*,*)orbs%norb_par(:),norb*orbs%nkpts
!!!!!$     stop
!!!!!$  end if
!!!!!$
!!!!!$  !calculate the k-points related quantities
!!!!!$  allocate(mykpts(orbs%nkpts),stat=i_stat)
!!!!!$  call memocc(i_stat,mykpts,'mykpts',subname)
!!!!!$
!!!!!$  call parallel_repartition_per_kpoints(iproc,nproc,orbs%nkpts,norb,orbs%norb_par,&
!!!!!$       orbs%nkptsp,mykpts,norb_par)
!!!!!$  if (orbs%norb_par(iproc) >0) then
!!!!!$     orbs%iskpts=mykpts(1)-1
!!!!!$  else
!!!!!$     orbs%iskpts=0
!!!!!$  end if
!!!!!$  i_all=-product(shape(mykpts))*kind(mykpts)
!!!!!$  deallocate(mykpts,stat=i_stat)
!!!!!$  call memocc(i_stat,i_all,'mykpts',subname)
!!!
!!!  !new system for k-point repartition
!!!  call kpts_to_procs_via_obj(nproc,orbs%nkpts,norb,norb_par)
!!!  !assign the values for norb_par and check the distribution
!!!  norb_tot=0
!!!  do jproc=0,nproc-1
!!!     if (jproc==iproc) orbs%isorb=norb_tot
!!!     do ikpt=1,orbs%nkpts
!!!        orbs%norb_par(jproc,0)=orbs%norb_par(jproc,0)+norb_par(jproc,ikpt)
!!!     end do
!!!     norb_tot=norb_tot+orbs%norb_par(jproc,0)
!!!  end do
!!!
!!!  if(norb_tot /= norb*orbs%nkpts) then
!!!     write(*,*)'ERROR: partition of orbitals incorrect, report bug.'
!!!     write(*,*)orbs%norb_par(:,0),norb*orbs%nkpts
!!!     stop
!!!  end if
!!!
!!!
!!!  !allocate(orbs%ikptsp(orbs%nkptsp+ndebug),stat=i_stat)
!!!  !call memocc(i_stat,orbs%ikptsp,'orbs%ikptsp',subname)
!!!  !orbs%ikptsp(1:orbs%nkptsp)=mykpts(1:orbs%nkptsp)
!!!
!!!  !this array will be reconstructed in the orbitals_communicators routine
!!!  i_all=-product(shape(norb_par))*kind(norb_par)
!!!  deallocate(norb_par,stat=i_stat)
!!!  call memocc(i_stat,i_all,'norb_par',subname)
!!!
!!!  !assign the values of the orbitals data
!!!  orbs%norb=norb
!!!  orbs%norbp=orbs%norb_par(iproc,0)
!!!  orbs%norbu=norbu
!!!  orbs%norbd=norbd
!!!
!!!  ! Modify these values
!!!  call repartitionOrbitals2(iproc, nproc, orbs%norb, orbs%norb_par, orbs%norbp, orbs%isorb)
!!!
!!!
!!!  allocate(orbs%iokpt(orbs%norbp+ndebug),stat=i_stat)
!!!  call memocc(i_stat,orbs%iokpt,'orbs%iokpt',subname)
!!!
!!!  !assign the k-point to the given orbital, counting one orbital after each other
!!!  jorb=0
!!!  do ikpt=1,orbs%nkpts
!!!     do iorb=1,orbs%norb
!!!        jorb=jorb+1 !this runs over norb*nkpts values
!!!        if (jorb > orbs%isorb .and. jorb <= orbs%isorb+orbs%norbp) then
!!!           orbs%iokpt(jorb-orbs%isorb)=ikpt
!!!        end if
!!!     end do
!!!  end do
!!!
!!!  !allocate occupation number and spinsign
!!!  !fill them in normal way
!!!  allocate(orbs%occup(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
!!!  call memocc(i_stat,orbs%occup,'orbs%occup',subname)
!!!  allocate(orbs%spinsgn(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
!!!  call memocc(i_stat,orbs%spinsgn,'orbs%spinsgn',subname)
!!!  orbs%occup(1:orbs%norb*orbs%nkpts)=1.0_gp 
!!!  do ikpt=1,orbs%nkpts
!!!     do iorb=1,orbs%norbu
!!!        orbs%spinsgn(iorb+(ikpt-1)*orbs%norb)=1.0_gp
!!!     end do
!!!     do iorb=1,orbs%norbd
!!!        orbs%spinsgn(iorb+orbs%norbu+(ikpt-1)*orbs%norb)=-1.0_gp
!!!     end do
!!!  end do
!!!
!!!  !put a default value for the fermi energy
!!!  orbs%efermi = UNINITIALIZED(orbs%efermi)
!!!  !and also for the gap
!!!  orbs%HLgap = UNINITIALIZED(orbs%HLgap)
!!!
!!!  ! allocate inwhichlocreg
!!!
!!!  allocate(orbs%inwhichlocreg(orbs%norb*orbs%nkpts),stat=i_stat)
!!!  call memocc(i_stat,orbs%inwhichlocreg,'orbs%inwhichlocreg',subname)
!!!  ! default for inwhichlocreg
!!!  orbs%inwhichlocreg = 1
!!!
!!!  !nullify(orbs%inwhichlocregP)
!!!
!!!  !allocate the array which assign the k-point to processor in transposed version
!!!  allocate(orbs%ikptproc(orbs%nkpts+ndebug),stat=i_stat)
!!!  call memocc(i_stat,orbs%ikptproc,'orbs%ikptproc',subname)
!!!
!!!  !initialize the starting point of the potential for each orbital (to be removed?)
!!!  allocate(orbs%ispot(orbs%norbp),stat=i_stat)
!!!  call memocc(i_stat,orbs%ispot,'orbs%ispot',subname)
!!!
!!!
!!!  ! Define two new arrays:
!!!  ! - orbs%isorb_par is the same as orbs%isorb, but every process also knows
!!!  !   the reference orbital of each other process.
!!!  ! - orbs%onWhichMPI indicates on which MPI process a given orbital
!!!  !   is located.
!!!  allocate(orbs%isorb_par(0:nproc-1), stat=i_stat)
!!!  call memocc(i_stat, orbs%isorb_par, 'orbs%isorb_par', subname)
!!!  allocate(orbs%onWhichMPI(sum(orbs%norb_par(:,0))), stat=i_stat)
!!!  call memocc(i_stat, orbs%onWhichMPI, 'orbs%onWhichMPI', subname)
!!!  iiorb=0
!!!  orbs%isorb_par=0
!!!  do jproc=0,nproc-1
!!!      do iorb=1,orbs%norb_par(jproc,0)
!!!          iiorb=iiorb+1
!!!          orbs%onWhichMPI(iiorb)=jproc
!!!      end do
!!!      if(iproc==jproc) then
!!!          orbs%isorb_par(jproc)=orbs%isorb
!!!      end if
!!!  end do
!!!  call MPI_Initialized(mpiflag,ierr)
!!!  if(mpiflag /= 0 .and. nproc > 1) call fmpi_allreduce(orbs%isorb_par(0), nproc, FMPI_SUM, bigdft_mpi%mpi_comm, ierr)
!!!
!!!END SUBROUTINE orbitals_descriptors_forLinear


!> Routine which assigns to each processor the repartition of nobj*nkpts objects
subroutine kpts_to_procs_via_obj(nproc,nkpts,nobj,nobj_par)
  use module_base
  implicit none
  integer, intent(in) :: nproc !< No. of proc
  integer, intent(in) :: nkpts !< No. K points
  integer, intent(in) :: nobj  !< Object number (i.e. nvctr)
  integer, dimension(0:nproc-1,nkpts), intent(out) :: nobj_par !< iresult of the partition
  !local varaibles
  logical :: intrep
  integer :: jproc,ikpt,iobj,nobjp_max_kpt,nprocs_with_floor,jobj,nobjp
  integer :: jkpt,nproc_per_kpt,nproc_left,kproc,nkpt_per_proc,nkpts_left
  real(gp) :: robjp,rounding_ratio

  call f_routine(id='kpts_to_procs_via_obj')

  !decide the naive number of objects which should go to each processor.
  robjp=real(nobj,gp)*real(nkpts,gp)/real(nproc,gp)
  !print *,'hereweare',robjp,nobj   
  !maximum number of objects which has to go to each processor per k-point
  nobjp_max_kpt=ceiling(modulo(robjp-epsilon(1.0_gp),real(nobj,gp)))


!see the conditions for the integer repartition of k-points
  if (nobjp_max_kpt == nobj .or. (nobjp_max_kpt==1 .and. robjp < 1.0_gp)) then
     intrep=.true.
     rounding_ratio=0.0_gp
     nprocs_with_floor=0
  else
     intrep=.false.
     !the repartition is not obvious, some processors take nobj_max_kpt objects, others take the previous integer.
     !to understand how many, we round the percentage of processors which is given by
     rounding_ratio=(robjp-real(floor(robjp), gp))
     !then this is the number of processors which will take the floor
     nprocs_with_floor=ceiling((1.0_gp-rounding_ratio)*real(nproc,gp))!nproc-(nobj*nkpts-floor(robjp)*nproc)
     !print *,'rounding_ratio,nprocs_with_floor',rounding_ratio,nprocs_with_floor
     if (nprocs_with_floor > nproc) stop 'ERROR: should not happen'
     !if (nprocs_with_floor == nproc) nprocs_with_floor=nproc-1
  end if

  !start separating the objects for the repartition which is suggested by rounding_ratio and nprocs_with_floor
  nobj_par(0:nproc-1,1:nkpts)=0
  !integer repartition
  if (intrep) then 
     !strategy for the repartition
     if (nproc >= nkpts) then
        !decide in how many processors a single k-point can be partitioned
        nproc_per_kpt=max((nproc-1),1)/nkpts !this is the minimum
        !count how many processors are left that way
        !distribute the k-point among these
        nproc_left=nproc-nproc_per_kpt*nkpts
        ikpt=0
        jproc=0
        !print *,'here',nproc_left,nproc_per_kpt
        do kproc=0,nproc_left-1
           ikpt=ikpt+1
           if (ikpt > nkpts) stop 'ERROR: also this should not happen3'
           do iobj=0,nobj-1
              nobj_par(jproc+modulo(iobj,nproc_per_kpt+1),ikpt)=nobj_par(jproc+modulo(iobj,nproc_per_kpt+1),ikpt)+1
           end do
           jproc=jproc+nproc_per_kpt+1
        end do
        !print *,'debug'
        if ((nproc_per_kpt+1)*nproc_left < nproc) then
           do jproc=(nproc_per_kpt+1)*nproc_left,nproc-1,nproc_per_kpt
              ikpt=ikpt+1
              !print *,'passed through here',modulo(nproc,nkpts),nkpts,ikpt,nproc_per_kpt,nproc,jproc,nproc_left
              if (ikpt > nkpts .or. jproc > nproc-1) stop 'ERROR: also this should not happen3b'
              do iobj=0,nobj-1
                 nobj_par(jproc+modulo(iobj,nproc_per_kpt),ikpt)=nobj_par(jproc+modulo(iobj,nproc_per_kpt),ikpt)+1
              end do
           end do
        end if
        !print *,'passed through here',modulo(nproc,nkpts),nkpts,ikpt,nproc_per_kpt,nproc
     else
        !decide in how many kpoints single processor can be partitioned
        nkpt_per_proc=max((nkpts-1),1)/nproc !this is the minimum
        !count how many k-points are left that way
        !distribute the processors among these
        nkpts_left=nkpts-nkpt_per_proc*nproc
        ikpt=1
        jproc=-1
        !print *,'hello',nkpts_left,nkpts_per_proc
        do jkpt=1,nkpts_left
           jproc=jproc+1
           if (jproc > nproc-1) stop 'ERROR: also this should not happen4'
           do iobj=0,(nobj)*(nkpt_per_proc+1)-1
              nobj_par(jproc,ikpt+modulo(iobj,nkpt_per_proc+1))=nobj_par(jproc,ikpt+modulo(iobj,nkpt_per_proc+1))+1
           end do
           ikpt=ikpt+nkpt_per_proc+1
        end do
        !print *,'ciao'
        if ((nkpt_per_proc+1)*nkpts_left < nkpts) then
           do ikpt=(nkpt_per_proc+1)*nkpts_left+1,nkpts,nkpt_per_proc
              jproc=jproc+1
              !print *,'passed through here',modulo(nproc,nkpts),nkpts,ikpt,nproc_per_kpt,nproc,jproc,nproc_left
              if (ikpt > nkpts .or. jproc > nproc-1) stop 'ERROR: also this should not happen4b'
              do iobj=0,(nobj)*(nkpt_per_proc)-1
                 nobj_par(jproc,ikpt+modulo(iobj,nkpt_per_proc))=nobj_par(jproc,ikpt+modulo(iobj,nkpt_per_proc))+1
              end do
           end do
        end if
           !print *,'passed through here',modulo(nproc,nkpts),nkpts,ikpt,nproc_per_kpt,nproc
     end if
  else
     !non-integer repartition
     iobj=0
     ikpt=0
     do jproc=0,nproc-2 !leave the last processor at the end
        nobjp=floor(robjp)
        !respect the rounding ratio
        if (nproc-jproc > nprocs_with_floor) nobjp=nobjp+1
        !print *,'jproc,nobjp',jproc,nobjp,nkpts,nobj,nkpts*nobj,iobj,nprocs_with_floor
        do jobj=1,nobjp
           if (modulo(iobj,nobj) ==0) ikpt=ikpt+1
           iobj=iobj+1
           if (iobj > nobj*nkpts) stop 'ERROR: also this should not happen'
           nobj_par(jproc,ikpt)=nobj_par(jproc,ikpt)+1
        end do
     end do
     !in the last processor we put the objects which are lacking
     nobjp=nobj*nkpts-iobj
     do jobj=1,nobjp
        if (modulo(iobj,nobj) ==0) ikpt=ikpt+1
        iobj=iobj+1
        !print *,'finished',jobj,nobjp,iobj,nobj*nkpts,jproc,ikpt
        if (iobj > nobj*nkpts) stop 'ERROR: also this should not happen2'
        nobj_par(nproc-1,ikpt)=nobj_par(nproc-1,ikpt)+1
     end do
  end if

  call f_release_routine()

END SUBROUTINE kpts_to_procs_via_obj


subroutine components_kpt_distribution(nproc,nkpts,norb,nvctr,norb_par,nvctr_par)
  use module_base, only: gp, f_err_throw, f_zero,BIGDFT_RUNTIME_ERROR,&
       UNINITIALIZED
  implicit none
  !Arguments
  integer, intent(in) :: nproc,nkpts,nvctr,norb
  integer, dimension(0:nproc-1,nkpts), intent(in) :: norb_par
  integer, dimension(0:nproc-1,nkpts), intent(out) :: nvctr_par
  !local variables
  integer :: ikpt,jsproc,jeproc,kproc,icount,ivctr,jproc,numproc
  real(gp) :: strprc,endprc

  !for any of the k-points find the processors which have such k-point associated
  call f_zero(nvctr_par)

  !Loop over each k point
  do ikpt=1,nkpts
     jsproc=UNINITIALIZED(1)
     jeproc=UNINITIALIZED(1)
     find_start: do jproc=0,nproc-1
        if(norb_par(jproc,ikpt) > 0) then 
           jsproc=jproc
           exit find_start
        end if
     end do find_start
     if (jsproc == UNINITIALIZED(1)) call f_err_throw('ERROR in kpt assignments',err_id=BIGDFT_RUNTIME_ERROR)
     if(norb_par(jsproc,ikpt) /= norb) then
        strprc=real(norb_par(jsproc,ikpt),gp)/real(norb,gp)     
     else
        strprc=1.0_gp
     end if
     if (ikpt < nkpts) then
        find_end: do jproc=jsproc,nproc-1
           if(norb_par(jproc,ikpt+1) > 0) then
              if (norb_par(jproc,ikpt)==0) then
                 jeproc=jproc-1
              else
                 jeproc=jproc
              end if
              exit find_end
           end if
        end do find_end
        if (jeproc == UNINITIALIZED(1)) call f_err_throw('ERROR in kpt assignments',err_id=BIGDFT_RUNTIME_ERROR)
     else
        jeproc=nproc-1
     end if
     if (jeproc /= jsproc) then
        endprc=real(norb_par(jeproc,ikpt),gp)/real(norb,gp)     
     else
        endprc=0.0_gp
     end if
     !if the number of processors is bigger than the number of orbitals this means 
     !that strprc and endprc are not correctly evaluated
     !evaluate the percentage on the number of components
     if (jeproc-jsproc+1 > norb) then
        strprc=1.0_gp/real(jeproc-jsproc+1,gp)
        endprc=strprc
     end if

     !assign the number of components which corresponds to the same orbital distribution
     numproc=jeproc-jsproc+1
     icount=0

     !print *,'kpoint',ikpt,jsproc,jeproc,strprc,endprc,ceiling(strprc*real(nvctr,gp)),nvctr
     !start filling the first processor
     ivctr=min(ceiling(strprc*real(nvctr,gp)-epsilon(1.0_gp)),nvctr)
     nvctr_par(jsproc,ikpt)=ivctr!min(ceiling(strprc*real(nvctr,gp)),nvctr)
     fill_array: do 
        if (ivctr==nvctr) exit fill_array
        icount=icount+1
        kproc=jsproc+modulo(icount,numproc)
        !put the floor of the components to the first processor
        if (strprc /= 1.0_gp .and. kproc==jsproc .and. &
             nvctr_par(kproc,ikpt)==ceiling(strprc*real(nvctr,gp)-epsilon(1.0_gp))) then
           !do nothing, skip away
        else
           nvctr_par(kproc,ikpt) = nvctr_par(kproc,ikpt)+1
           ivctr=ivctr+1
        end if
     end do fill_array
     !print '(a,i3,i3,i6,2(1pe25.17),i7,20i5)','here',ikpt,jsproc,jeproc,strprc,endprc,sum(nvctr_par(:,ikpt)),nvctr_par(:,ikpt)
     !print '(a,i3,i3,i6,2(1pe25.17),i7,20i5)','there',ikpt,jsproc,jeproc,strprc,endprc,sum(nvctr_par(:,ikpt)),norb_par(:,ikpt)
  end do

END SUBROUTINE components_kpt_distribution


!> Check the distribution of k points over the processors
subroutine check_kpt_distributions(nproc,nkpts,norb,ncomp,norb_par,ncomp_par,info,lub_orbs,lub_comps)
  use module_base
  implicit none
  integer, intent(in) :: nproc,nkpts,norb,ncomp
  integer, dimension(0:nproc-1,nkpts), intent(in) :: norb_par
  integer, dimension(0:nproc-1,nkpts), intent(in) :: ncomp_par
  integer, intent(inout) :: info
  integer, intent(out) :: lub_orbs,lub_comps
  !local variables
  character(len=*), parameter :: subname='check_kpt_distributions'
  logical :: notcompatible,couldbe
  integer :: ikpt,jproc,norbs,ncomps,kproc,ieproc,isproc,jkpt
  integer, dimension(:,:), allocatable :: load_unbalancing
  !before printing the distribution schemes, check that the two distributions contain
  !the same k-points
  if (info == 0) call print_distribution_schemes(nproc,nkpts,norb_par,ncomp_par)

  do ikpt=1,nkpts
     isproc=UNINITIALIZED(1)
     find_isproc : do kproc=0,nproc-1
        if (ncomp_par(kproc,ikpt) > 0) then
           isproc=kproc
           exit find_isproc
        end if
     end do find_isproc
     !if (isproc == UNINITIALIZED(1)) stop 'ERROR(check_kpt_distributions): isproc cannot be found'
     if (isproc == UNINITIALIZED(1)) call f_err_throw( &
        & 'isproc cannot be found',err_name='BIGDFT_RUNTIME_ERROR')
     ieproc=UNINITIALIZED(1)
     find_ieproc : do kproc=nproc-1,0,-1
        if (ncomp_par(kproc,ikpt) > 0) then
           ieproc=kproc
           exit find_ieproc
        end if
     end do find_ieproc
     !if (ieproc == UNINITIALIZED(1)) stop 'ERROR(check_kpt_distributions): ieproc cannot be found'
     if (ieproc == UNINITIALIZED(1)) call f_err_throw( &
        & 'ieproc cannot be found', err_name='BIGDFT_RUNTIME_ERROR')

     norbs=0
     ncomps=0
     do jproc=0,nproc-1
        !count the total number of components
        norbs=norbs+norb_par(jproc,ikpt)
        ncomps=ncomps+ncomp_par(jproc,ikpt)
        notcompatible=(ncomp_par(jproc,ikpt) == 0 .neqv. norb_par(jproc,ikpt) == 0) 
        !check whether there are only 0 orbitals
        if (notcompatible .and. norb_par(jproc,ikpt)==0) then
           !if the processor is the last one then there should not be other k-points on this processors
           couldbe=.false.
           if (jproc == ieproc) then
              couldbe=.true.
              do jkpt=ikpt+1,nkpts
                 couldbe=couldbe .and. (norb_par(jproc,jkpt) ==0 .and. ncomp_par(jproc,jkpt)==0)
              end do
           end if
           if ((isproc < jproc .and. jproc < ieproc) .or. couldbe) notcompatible=.false.
        end if
        if (notcompatible) then     
           if (info == 0) write(*,*)' ERROR: processor ', jproc,' kpt,',ikpt,&
                'have components and orbital distributions not compatible'
           info=1
           return
           !call MPI_ABORT(bigdft_mpi%mpi_comm, ierr)
        end if
     end do
     if (norb/=norbs .or. ncomps /= ncomp) then
        if (info == 0) write(*,*)' ERROR: kpt,',ikpt,&
             'has components or orbital distributions not correct'
        info=2
        return
        !call MPI_ABORT(bigdft_mpi%mpi_comm, ierr)
     end if
  end do

  load_unbalancing = f_malloc((/ 0.to.nproc-1, 1.to.2 /),id='load_unbalancing')

  do jproc=0,nproc-1
     load_unbalancing(jproc,:)=0
     do ikpt=1,nkpts
        load_unbalancing(jproc,1)=load_unbalancing(jproc,1)+norb_par(jproc,ikpt)
        load_unbalancing(jproc,2)=load_unbalancing(jproc,2)+ncomp_par(jproc,ikpt)
     end do
  end do

  !calculate the maximum load_unbalancing
  lub_orbs=0
  lub_comps=0
  do jproc=0,nproc-1
     do kproc=0,nproc-1
        lub_orbs=max(lub_orbs,load_unbalancing(jproc,1)-load_unbalancing(kproc,1))
        lub_comps=max(lub_comps,load_unbalancing(jproc,2)-load_unbalancing(kproc,2))
     end do
  end do

  if (info==0) write(*,*)' Kpoints Distribuitions are compatible, load unbalancings, orbs,comps:',lub_orbs,&
       '/',max(minval(load_unbalancing(:,1)),1),lub_comps,'/',minval(load_unbalancing(:,2))
  info=0
  call f_free(load_unbalancing)


END SUBROUTINE check_kpt_distributions

!> Routine which associates to any of the processor a given number of objects
!! depending of the number of processors and k-points
subroutine parallel_repartition_with_kpoints(nproc,nkpts,nobj,nobj_par)
  use module_base
  implicit none
  integer, intent(in) :: nkpts,nobj,nproc
  integer, dimension(0:nproc-1), intent(out) :: nobj_par
  !local variables
  integer :: n_i,n_ip,rs_i,N_a,N_b,N_c,ikpt,jproc,i,ntmp
 !!$  real(gp) :: rtmp

  ! Strategy to divide between k points.
  ! There is an nproc length to divide into orbs%nkpts segments.
  ! Segment (ikpt - 1) expand in 0 <= r_i < r_ip <= nproc.
  ! where r_i and r_ip are real values. There are two possibilities:
  !  - We can write r_i <= n_i <= n_ip <= r_ip with n_i and n_ip integers ;
  !  - or r_i <= n_i and n_ip <= r_ip and n_i = n_ip + 1.
  ! For both cases, we can divide nobj into the partition (real values):
  !  - N_a = (n_i - r_i)*nobj*nkpts/nproc (the initial part);
  !  - N_b = max((n_ip - n_i)*nobj*nkpts / nproc, 0) (the naive part, the only one if nkpts is a multiple of nproc);
  !  - N_c = (r_ip - n_ip) * nobj * orbs%nkpts / nproc (the final part);
  ! Before going to integer values, we have r_i = (ikpt - 1) * nproc / orbs%nkpts (the naive division)
  ! and r_ip = (ikpt) * nproc / orbs%nkpts (the segment endpoint)
  ! So N_a and N_b can be simplified and written instead:
  !  - N_a = int(nobj * (n_i * orbs%nkpts - (ikpt - 1) * nproc) / nproc);
  !  - N_c = int(nobj * ((ikpt) * nproc - n_ip * orbs%nkpts) / nproc)
  !  - N_b = nobj - N_a - N_c 
  ! After, if N_a > 0, we put this quantity to proc n_i - 1, if N_c > 0
  ! we put its quantity to proc n_ip ; and finally N_b is distributed
  ! among [n_i;n_ip[ procs.

  nobj_par(:)=0
  do ikpt=1,nkpts
     ! Calculation of n_i and n_ip, rs_i = r_i * orbs%nkpts to avoid rounding.
     rs_i=(ikpt-1)*nproc !integer variable for rounding purposes

     if (mod(rs_i,nkpts) == 0) then
        n_i=rs_i/nkpts 
     else
        n_i=rs_i/nkpts+1
     end if

     rs_i=ikpt*nproc
     n_ip=rs_i/nkpts
 !!$     print *,'ikpt,ni,nip',ikpt,n_i,n_ip
     ! Calculation of N_a, N_b and N_c from given n_i and n_ip.
     if (n_ip >= n_i) then
        ntmp = (n_i*nkpts-(ikpt-1)*nproc) * nobj
        if (modulo(ntmp, nproc) == 0) then
           N_a = ntmp / nproc
        else
           N_a = (ntmp - modulo(ntmp, nproc) + nproc) / nproc
        end if
 !!$        ntmp=n_i*nkpts-(ikpt-1)*nproc
 !!$        rtmp=real(nobj,gp)/real(nproc,gp)
 !!$        rtmp=rtmp*real(ntmp,gp)
 !!$        N_a=floor(rtmp)
 !!$        if (iproc == 0) print *,'ikpts,rtmp',ikpt,rtmp
        ntmp = (ikpt*nproc-n_ip*nkpts) * nobj
        if (modulo(ntmp, nproc) == 0) then
           N_c = ntmp / nproc
        else
           N_c = (ntmp - modulo(ntmp, nproc) + nproc) / nproc
        end if

 !!$        ntmp=ikpt*nproc-n_ip*nkpts
 !!$        rtmp=real(nobj,gp)/real(nproc,gp)
 !!$        rtmp=rtmp*real(ntmp,gp)
 !!$        N_c=ceiling(rtmp)
 !!$        if (iproc == 0) print *,'ikpts,rtmp2',ikpt,rtmp,N_a,N_c
        !the corrections above are to avoid the 32 bit integer overflow
        !N_a=nint(real(nobj*(n_i*nkpts-(ikpt-1)*nproc),gp)/real(nproc,gp))
        !N_c=nint(real(nobj*(ikpt*nproc-n_ip*nkpts),gp)/real(nproc,gp))
     else
        N_c=nobj/2
        N_a=nobj-N_c
     end if
     N_b=nobj-N_a-N_c
     if (N_b == -1) then
        N_c = N_c - 1
        N_b = 0
     end if
 !!$     write(*,*) ikpt, N_a, N_b, N_c
     if (nkpts > 1 .and. N_b < n_ip - n_i) stop 'ERROR:parallel_repartion_with_kpoints'
     !assign to procs the objects.
     if (N_a>0) nobj_par(n_i-1)=nobj_par(n_i-1)+N_a
     if (N_b>0) then
        do i=0,N_b-1
           jproc=n_i+mod(i,n_ip-n_i)
           nobj_par(jproc)=nobj_par(jproc)+1
        end do
     end if
     if (N_c>0) nobj_par(n_ip)=nobj_par(n_ip)+N_c
  end do
END SUBROUTINE parallel_repartition_with_kpoints


subroutine parallel_repartition_per_kpoints(iproc,nproc,nkpts,nobj,nobj_par,&
     nkptsp,mykpts,nobj_pkpt)
  implicit none
  integer, intent(in) :: iproc,nproc,nkpts,nobj
  integer, dimension(0:nproc-1), intent(in) :: nobj_par
  integer, intent(out) :: nkptsp
  integer, dimension(nkpts), intent(out) :: mykpts
  integer, dimension(0:nproc-1,nkpts), intent(out) :: nobj_pkpt
  !local variables
  integer :: ikpts,jproc,jobj,norb_tot,iorbp

  !initialise the array
  do ikpts=1,nkpts
     do jproc=0,nproc-1
        nobj_pkpt(jproc,ikpts)=0 
     end do
  end do

  !assign the k-point, counting one object after each other
  jobj=1
  ikpts=1
  !print *,'here',nobj_par(:)
  do jproc=0,nproc-1
     do iorbp=1,nobj_par(jproc)
        nobj_pkpt(jproc,ikpts)=nobj_pkpt(jproc,ikpts)+1
        if (mod(jobj,nobj)==0) then
           ikpts=ikpts+1
        end if
        jobj=jobj+1
     end do
  end do
  !some checks
  if (nobj /= 0) then
     !check the distribution
     do ikpts=1,nkpts
        !print *,'partition',ikpts,orbs%nkpts,'ikpts',nobj_pkpt(:,ikpts)
        norb_tot=0
        do jproc=0,nproc-1
           norb_tot=norb_tot+nobj_pkpt(jproc,ikpts)
        end do
        if(norb_tot /= nobj) then
           write(*,*)'ERROR: partition of objects incorrect, kpoint:',ikpts
           stop
        end if
     end do
  end if

  !calculate the number of k-points treated by each processor in both
  ! the component distribution and the orbital distribution.
  nkptsp=0
  do ikpts=1,nkpts
     if (nobj_pkpt(iproc,ikpts) /= 0) then
        nkptsp=nkptsp+1
        mykpts(nkptsp) = ikpts
     end if
  end do

END SUBROUTINE parallel_repartition_per_kpoints

subroutine pawpatch_from_file( filename, atoms,ityp, paw_tot_l, &
     paw_tot_q, paw_tot_coefficients, paw_tot_matrices, &
     storeit)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  type(atoms_data), intent(inout) :: atoms
  integer , intent(in):: ityp 
  integer , intent(inout):: paw_tot_l, paw_tot_q, paw_tot_coefficients, paw_tot_matrices
  logical, intent(in) :: storeit

!! local variables  
  character(len=*), parameter :: subname='pawpatch_from_file'
  integer :: npawl, ipawl, paw_l
  integer :: paw_nofgaussians, paw_nofchannels, il, ierror, ig
  real(gp) :: paw_greal, paw_gimag, paw_ccoeff, paw_scoeff, dumpaw
  character(len=100) :: string

  !parameters for abscalc-paw

  if(.not. storeit) then
     !if(ityp == 1) then !this implies that the PSP are all present
     if (.not. associated(atoms%paw_NofL)) then
        atoms%paw_NofL = f_malloc_ptr(atoms%astruct%ntypes,id='atoms%paw_NofL')
     end if
     ! if (iproc.eq.0) write(*,*) 'opening PSP file ',filename
     open(unit=11,file=trim(filename),status='old',iostat=ierror)
     !Check the open statement
     if (ierror /= 0) then
        write(*,*) ': Failed to open the PAWpatch file "',&
             trim(filename),'"'
        stop
     end if
     
     !! search for paw_patch informations
     
     atoms%paw_NofL(ityp)=0
     do while(.true.)
        read(11,'(a)',iostat=ierror, END=110)  string
        if ( trim(string).eq. 'PAWPATCH') then
           exit
        endif
     end do
     !! explain_paw_psp_terms_in_atom_data()
     
     read(11,*) npawl
     
     atoms%paw_NofL(ityp) = npawl
     
     paw_tot_l = paw_tot_l + npawl
     do ipawl=1,npawl
        read(11,*) paw_l
        read(11,*) paw_greal
        read(11,*) paw_nofgaussians
        read(11,*) paw_nofchannels
        read(11,*)  string  !!  follow 300 PAW_Gimag factors
        paw_tot_q = paw_tot_q+paw_nofgaussians
        paw_tot_coefficients = paw_tot_coefficients + paw_nofchannels*paw_nofgaussians*2
        paw_tot_matrices=paw_tot_matrices+paw_nofchannels**2
        
        
        do ig=1, paw_nofgaussians
           read(11,*)  paw_gimag
        enddo
        read(11,*)  string  !!  !!  follow for each of the 7 channels 300 (cos_factor, sin_factor)  pairs
        do il=1, paw_nofchannels
           do ig=1, paw_nofgaussians
              read(11,*)  paw_ccoeff, paw_scoeff
           enddo
        enddo
        read(11,*)  string  !! pawpatch matrix
        do il=1, paw_nofchannels
           do ig=1, paw_nofchannels
              read(11,*)  dumpaw
           end do
        enddo
        read(11,*)  string  !! S  matrix
        do il=1, paw_nofchannels
           do ig=1, paw_nofchannels
              read(11,*)  dumpaw
           end do
        enddo
        
        read(11,*)  string  !! Sm1  matrix
        do il=1, paw_nofchannels
           do ig=1, paw_nofchannels
              read(11,*)  dumpaw
           end do
        enddo
     enddo
110  close(11)

  else
     if(ityp.eq.1) then
        atoms%paw_l   = f_malloc_ptr(paw_tot_l,id='atoms%paw_l  ')
        atoms%paw_nofchannels   = f_malloc_ptr(paw_tot_l,id='atoms%paw_nofchannels  ')
        atoms%paw_nofgaussians   = f_malloc_ptr(paw_tot_l,id='atoms%paw_nofgaussians  ')
        atoms%paw_Greal   = f_malloc_ptr(paw_tot_l,id='atoms%paw_Greal  ')
        atoms%paw_Gimag  = f_malloc_ptr(paw_tot_q   ,id='atoms%paw_Gimag ')
        atoms%paw_Gcoeffs  = f_malloc_ptr(paw_tot_coefficients  ,id='atoms%paw_Gcoeffs ')
        atoms%paw_H_matrices = f_malloc_ptr(paw_tot_matrices,id='atoms%paw_H_matrices')
        atoms%paw_S_matrices  = f_malloc_ptr(paw_tot_matrices  ,id='atoms%paw_S_matrices ')
        atoms%paw_Sm1_matrices  = f_malloc_ptr(paw_tot_matrices  ,id='atoms%paw_Sm1_matrices ')
        
        
        paw_tot_l=0
        paw_tot_q = 0
        paw_tot_coefficients = 0
        paw_tot_matrices= 0
     endif
     
     if( atoms%paw_NofL(ityp).gt.0) then

        open(unit=11,file=trim(filename),status='old')

        do while(.true.)
           read(11,'(a)', END=220)  string
           if ( trim(string).eq. 'PAWPATCH') then
              exit
           endif
        end do
220     continue
        if(trim(string) .ne. 'PAWPATCH') then
           print *, "paw section not found re-reading  file ", filename
           close(11)
           stop
        end if
        
        read(11,*) npawl
        atoms%paw_NofL(ityp) = npawl
        
        
        do ipawl=1,npawl
           !! explain_paw_psp_terms_in_atom_data()
           read(11,*) atoms%paw_l(paw_tot_l+ipawl  )
           read(11,*) atoms%paw_greal(paw_tot_l+ipawl  )
           read(11,*) atoms%paw_nofgaussians(paw_tot_l+ipawl  )
           read(11,*) atoms%paw_nofchannels(paw_tot_l+ipawl  )
           paw_nofchannels = atoms%paw_nofchannels(paw_tot_l+ipawl  )
           paw_nofgaussians = atoms%paw_nofgaussians(paw_tot_l+ipawl  )
           read(11,'(a)')  string  !!  follow  paw_nofchannels PAW_Gimag factors
           
           do ig=1, paw_nofgaussians
              read(11,*)  atoms%paw_Gimag(paw_tot_q + ig )
           enddo
           read(11,'(a)')  string  !!  !!  follow for each of the Npaw channels  nofgaussians (cos_factor, sin_factor)  pairs
           !!print *, string, " reading  " , filename
           
           do il=1, paw_nofchannels
              do ig=1, paw_nofgaussians
                 read(11,*)  atoms%paw_Gcoeffs( paw_tot_coefficients + 2*(il-1)*paw_nofgaussians+2*ig-1) , &
                      atoms%paw_Gcoeffs( paw_tot_coefficients + 2*(il-1)*paw_nofgaussians+2*ig)
              enddo
           enddo
           read(11,'(a)')  string  !! pawpatch matrix
           !!print *, string, " reading  " , filename
           do il=1, paw_nofchannels
              do ig=1, paw_nofchannels
                 read(11,*)  dumpaw 
                 atoms%paw_H_matrices(paw_tot_matrices +(il-1)*paw_nofchannels+ig   )=dumpaw
              end do
           enddo
           read(11,'(a)')  string  !! S  matrix
           !!print *, string, " reading >>>>>  " , filename
           do il=1, paw_nofchannels
              do ig=1, paw_nofchannels
                 read(11,*)  dumpaw
                 atoms%paw_S_matrices(paw_tot_matrices +(il-1)*paw_nofchannels+ig   )=dumpaw
              end do
           enddo
           read(11,'(a)')  string  !! Sm1  matrix
           !!print *, string, " reading  " , filename
           do il=1, paw_nofchannels
              do ig=1, paw_nofchannels
                 read(11,*)  dumpaw
                 atoms%paw_Sm1_matrices(paw_tot_matrices +(il-1)*paw_nofchannels+ig   )=dumpaw
                 !!print *, dumpaw
              end do
           enddo
           paw_tot_q = paw_tot_q+paw_nofgaussians
           paw_tot_coefficients = paw_tot_coefficients + paw_nofchannels*paw_nofgaussians*2
           paw_tot_matrices=paw_tot_matrices+paw_nofchannels**2
        end do
        paw_tot_l = paw_tot_l + npawl
     end if
     close(11)
  endif
end subroutine pawpatch_from_file

subroutine paw_init(iproc, paw, at, rxyz, d, dpbox, nspinor, npsidim, norb, nkpts)
  !use module_base
  use module_defs, only: gp
  use module_types, only: atoms_data, paw_objects, nullify_paw_objects
  use module_dpbox, only: denspot_distribution
  use locregs, only: grid_dimensions
  use dynamic_memory
  use m_paw_an, only: paw_an_init
  use m_paw_ij, only: paw_ij_init
  use m_pawcprj, only: pawcprj_alloc, pawcprj_getdim
  use m_pawfgrtab, only: pawfgrtab_init
  use abi_interfaces_add_libpaw, only: abi_initrhoij, abi_wvl_nhatgrid
  use abi_interfaces_geometry, only: abi_metric
  implicit none
  type(paw_objects), intent(out) :: paw
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3, at%astruct%nat), intent(in) :: rxyz
  type(grid_dimensions), intent(in) :: d
  type(denspot_distribution), intent(in) :: dpbox
  integer, intent(in) :: iproc, nspinor, npsidim, norb, nkpts

  integer, parameter :: pawxcdev = 1, pawspnorb = 0, nspinor_ = 1, cplex = 1
  integer :: i
  integer, dimension(:), allocatable :: nlmn, nattyp, l_size_atm, lexexch, lpawu, atindx1
  real(gp), dimension(:,:), allocatable :: spinat
  integer, parameter :: optcut = 0, optgr0 = 1, optgr1 = 0, optgr2 = 0, optrad = 0

  call nullify_paw_objects(paw)
  if (.not. associated(at%pawtab)) return

  paw%usepaw = .true.
  paw%ntypes = at%astruct%ntypes
  paw%natom  = at%astruct%nat
  paw%lmnmax = maxval(at%pawtab(:)%lmn_size)

  allocate(paw%paw_an(at%astruct%nat))
  call paw_an_init(paw%paw_an, at%astruct%nat, at%astruct%ntypes, 0, dpbox%nrhodim, 1, pawxcdev, &
       & at%astruct%iatype, at%pawang, at%pawtab, &
       & has_vxc=1, has_vxc_ex=1, has_vhartree=1) !, &
  !&   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

  allocate(paw%paw_ij(at%astruct%nat))
  call paw_ij_init(paw%paw_ij,cplex,nspinor_,dpbox%nrhodim,dpbox%nrhodim,&
       &   0,at%astruct%nat,at%astruct%ntypes,at%astruct%iatype,at%pawtab,&
       &   has_dij=1,has_dijhartree=1,has_dijso=0,has_dijhat=0,&
       &   has_pawu_occ=1,has_exexch_pot=1) !,&
  !&   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

  allocate(paw%cprj(at%astruct%nat, norb * nkpts * nspinor_))
  nlmn = f_malloc(at%astruct%nat, id = "nlmn")
  nattyp = f_malloc0(at%astruct%ntypes, id = "nattyp")
  do i = 1, at%astruct%nat
     nattyp(at%astruct%iatype(i)) = nattyp(at%astruct%iatype(i)) + 1
  end do
  call pawcprj_getdim(nlmn, at%astruct%nat, nattyp, at%astruct%ntypes, &
       & at%astruct%iatype, at%pawtab, 'O')
  call pawcprj_alloc(paw%cprj, 0, nlmn)
  call f_free(nlmn)

  allocate(paw%fgrtab(at%astruct%nat))
  l_size_atm = f_malloc(at%astruct%nat, id = "l_size_atm")
  do i = 1, at%astruct%nat
     l_size_atm(i) = at%pawtab(at%astruct%iatype(i))%lcut_size
  end do
  call pawfgrtab_init(paw%fgrtab, cplex, l_size_atm, dpbox%nrhodim, at%astruct%iatype) !,&
  !&     mpi_atmtab=mpi_enreg%my_atmtab,mpi_comm_atom=mpi_enreg%comm_atom)
  call f_free(l_size_atm)
  atindx1 = f_malloc(at%astruct%nat, id = "atindx1")
  do i = 1, at%astruct%nat
     atindx1(i) = i
  end do
  !ucvol = product(denspot%dpbox%mesh%ndims) * product(denspot%dpbox%mesh%hgrids)
  call abi_wvl_nhatgrid(atindx1, at%astruct%geocode, dpbox%mesh%hgrids, dpbox%i3s + dpbox%i3xcsh, &
       & size(at%pawtab), at%astruct%nat, nattyp, at%astruct%ntypes, &
       & d%n1, d%n1i, d%n2, d%n2i, d%n3, dpbox%n3pi, &
       & optcut, optgr0, optgr1, optgr2, optrad, &
       & paw%fgrtab, at%pawtab, 0, rxyz) !,&
  !&         mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
  !&         mpi_comm_fft=spaceComm_fft,distribfft=mpi_enreg%distribfft)
  call f_free(atindx1)
  call f_free(nattyp)

  allocate(paw%pawrhoij(at%astruct%nat))
  spinat = f_malloc0((/3, at%astruct%nat/), id = "spinat")

  lexexch = f_malloc(at%astruct%ntypes, id = "lexexch")
  lexexch = -1
  lpawu = f_malloc(at%astruct%ntypes, id = "lpawu")
  lpawu = -1
  call abi_initrhoij(cplex, lexexch, lpawu, &
       & at%astruct%nat, at%astruct%nat, dpbox%nrhodim, nspinor_, dpbox%nrhodim, &
       & at%astruct%ntypes, paw%pawrhoij, pawspnorb, at%pawtab, spinat, at%astruct%iatype) !,&
  !&   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
  call f_free(lpawu)
  call f_free(lexexch)
  call f_free(spinat)

  paw%spsi = f_malloc_ptr(npsidim,id='paw%spsi')
end subroutine paw_init

subroutine system_signaling(iproc, signaling, gmainloop, KSwfn, tmb, energs, denspot, optloop, &
       & ntypes, radii_cf, crmult, frmult)
  use module_defs, only: gp, UNINITIALIZED
  use module_types
  implicit none
  integer, intent(in) :: iproc, ntypes
  logical, intent(in) :: signaling
  double precision, intent(in) :: gmainloop
  type(DFT_wavefunction), intent(inout) :: KSwfn, tmb
  type(DFT_local_fields), intent(inout) :: denspot
  type(DFT_optimization_loop), intent(inout) :: optloop
  type(energy_terms), intent(inout) :: energs
  real(gp), dimension(ntypes,3), intent(in) :: radii_cf
  real(gp), intent(in) :: crmult, frmult

  if (signaling) then
     ! Only iproc 0 has the C wrappers.
     if (iproc == 0) then
        call wf_new_wrapper(KSwfn%c_obj, KSwfn, 0)
        call wf_copy_from_fortran(KSwfn%c_obj, radii_cf, crmult, frmult)
        call wf_new_wrapper(tmb%c_obj, tmb, 1)
        call wf_copy_from_fortran(tmb%c_obj, radii_cf, crmult, frmult)
        call bigdft_signals_add_wf(gmainloop, KSwfn%c_obj, tmb%c_obj)
        !call energs_new_wrapper(energs%c_obj, energs)
        call bigdft_signals_add_energs(gmainloop, energs%c_obj)
        call localfields_new_wrapper(denspot%c_obj, denspot)
        call bigdft_signals_add_denspot(gmainloop, denspot%c_obj)
        call optloop_new_wrapper(optLoop%c_obj, optLoop)
        call bigdft_signals_add_optloop(gmainloop, optLoop%c_obj)
     else
        KSwfn%c_obj   = UNINITIALIZED(KSwfn%c_obj)
        tmb%c_obj     = UNINITIALIZED(tmb%c_obj)
        denspot%c_obj = UNINITIALIZED(denspot%c_obj)
        optloop%c_obj = UNINITIALIZED(optloop%c_obj)
     end if
  else
     KSwfn%c_obj  = 0
     tmb%c_obj    = 0
  end if
END SUBROUTINE system_signaling
