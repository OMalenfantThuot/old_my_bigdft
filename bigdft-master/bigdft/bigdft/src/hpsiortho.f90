!> @file
!!  Application of the Hamiltonian + orthonormalize constraints
!! @author
!!    Copyright (C) 2007-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Calculates the application of the Hamiltonian on the wavefunction. The hamiltonian can be self-consistent or not.
!! In the latter case, the potential should be given in the rhov array of denspot structure.
!! Otherwise, rhov array is filled by the self-consistent density
subroutine psitohpsi(iproc,nproc,atoms,scf,denspot,itrp,itwfn,scf_mode,alphamix,&
     nlpsp,linflag,unblock_comms,GPU,wfn,&
     energs,rpnrm,xcstr)
  use module_base
  use module_types
  use module_interfaces, only: LocalHamiltonianApplication, SynchronizeHamiltonianApplication, &
       & XC_potential, communicate_density, free_full_potential, sumrho
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use module_mixing
  use yaml_output
  use rhopotential, only: full_local_potential
  use public_enums
  use rhopotential, only: updatePotential,exchange_and_correlation
  use module_cfd, only: cfd_dump_info, cfd_field, cfd_is_converged
  use f_enums, only: toi
  use module_asd, only: asd_data, asd_wrapper, asd
  !
  implicit none
  !Arguments
  logical, intent(in) :: scf  !< If .false. do not calculate the self-consistent potential
  integer, intent(in) :: iproc,nproc,itrp,linflag,itwfn
  type(f_enumerator), intent(in) :: scf_mode
  character(len=3), intent(in) :: unblock_comms
  real(gp), intent(in) :: alphamix
  type(atoms_data), intent(in) :: atoms
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(DFT_local_fields), intent(inout) :: denspot
  type(energy_terms), intent(inout) :: energs
  type(DFT_wavefunction), intent(inout) :: wfn
  type(GPU_pointers), intent(inout) :: GPU
  !real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
  real(gp), intent(inout) :: rpnrm
  real(gp), dimension(6), intent(out) :: xcstr
  !real(wp), dimension(orbs%npsidim_orbs), intent(out) :: hpsi
  !local variables
  character(len=*), parameter :: subname='psitohpsi'
  logical :: unblock_comms_den,unblock_comms_pot,whilepot,savefields
  integer :: nthread_max,ithread,nthread,irhotot_add,irho_add,ispin,correcth,ishift
  real(gp) :: ehart_ps
  !integer :: ii,jj
  !$ integer :: omp_get_max_threads,omp_get_thread_num,omp_get_num_threads
  real(gp) :: compch_sph
  !$ real(wp), dimension(:), allocatable :: temp,m_norm,temp2 !to be removed
  !!! Temporary asd data type. To be moved..
  !!type(asd_data) :: asd


  call f_routine(id=subname)
  !in the default case, non local hamiltonian is done after potential creation
  whilepot=.true.
  !if (wfn%paw%usepaw) whilepot = .false.
  !flag for saving the local fields (rho,vxc,vh)
  savefields= (toi(scf_mode)==SCF_KIND_GENERALIZED_DIRMIN)
  correcth=1
  !do not do that if rho_work is already associated
  if (savefields .and. associated(denspot%rho_work)) then
     !flag for correcting the hamiltonian (either is false or toggles savefields)
     correcth=2
     savefields=.false.
  end if

  nthread_max=1
  ithread=0
  nthread=1
  !control if we have the possibility of using OMP to de-synchronize communication
  !$ nthread_max=omp_get_max_threads()


  !decide the communication strategy
  unblock_comms_den=mpi_thread_funneled_is_supported .and. unblock_comms=='DEN' .and. &
      nthread_max > 1 .and. scf .and. .not. GPU%OCLconv! density is done only if scf is present
  if (unblock_comms_den) whilepot=.false. !anticipate the NlHamiltonian if density should be overlapped

  unblock_comms_pot=mpi_thread_funneled_is_supported .and. unblock_comms=='POT' .and. &
      nthread_max > 1 .and. whilepot .and. .not. GPU%OCLconv

  if ((unblock_comms_den .or. unblock_comms_pot) .and. iproc==0) then
     if (unblock_comms_den) call yaml_map('Overlapping communication of','Density')
     if (unblock_comms_pot) call yaml_map('Overlapping communication of','Potential')
     call yaml_map('No. of OMP Threads',nthread_max)
     call yaml_newline()
  end if

  !zero the output wavefunction
  if (wfn%orbs%npsidim_orbs > 0) call f_zero(wfn%orbs%npsidim_orbs,wfn%hpsi(1))

  !calculate the self-consistent potential
  if (scf) then
     !update the entropic energy
     energs%eTS=wfn%orbs%eTS
     !safe the previous value of the energy
     energs%e_prev=energs%energy
     !print *,'here',savefields,correcth,energs%ekin,energs%epot,dot(wfn%orbs%npsidim_orbs,wfn%psi(1),1,wfn%psi(1),1)
     ! Potential from electronic charge density
     call sumrho(denspot%dpbox,wfn%orbs,wfn%Lzd,GPU,atoms%astruct%sym,denspot%rhod,denspot%xc,wfn%psi,denspot%rho_psi)
     !print *,'here',wfn%orbs%occup(:),'there',savefields,correcth,energs%ekin,energs%epot,&
     !     dot(wfn%Lzd%Glr%d%n1i*wfn%Lzd%Glr%d%n2i*denspot%dpbox%n3p*wfn%orbs%nspin,&
     !     denspot%rho_psi(1,1),1,denspot%rho_psi(1,1),1)

     !initialize nested approach
     !this has always to be done for using OMP parallelization in the
     !projector case
     !if nesting is not supported, a bigdft_nesting routine should not be called
     !$ if (unblock_comms_den) then
     !$ call timing(iproc,'UnBlockDen    ','ON')
     !$ call bigdft_open_nesting(2)
     !$ end if
     !print *,'how many threads ?' ,nthread_max
     !$OMP PARALLEL IF(unblock_comms_den) DEFAULT(shared), PRIVATE(ithread,nthread)
     !$ ithread=omp_get_thread_num()
     !$ nthread=omp_get_num_threads() !this should be 2 if active
     !print *,'hello, I am thread no.',ithread,' out of',nthread,'of iproc', iproc
     ! thread 0 does mpi communication
     if (ithread == 0) then
        !$ if (unblock_comms_den) call OMP_SET_NUM_THREADS(1)
        !communicate density
        !the rho_p pointer, allocated aoutside form the nested region, is by default
        !freed by communicate_density routine in the nested region
        call communicate_density(denspot%dpbox,wfn%orbs%nspin,denspot%rhod,&
             denspot%rho_psi,denspot%rhov,unblock_comms_den)
        !write(*,*) 'node:', iproc, ', thread:', ithread, 'mpi communication finished!!'
     end if

     if (wfn%paw%usepaw) then
        !write(*,*) "DENSITY", sum(denspot%rhov) * product(denspot%dpbox%hgrids)
        call paw_compute_rhoij(wfn%paw, wfn%orbs, atoms)
        call paw_update_rho(wfn%paw, denspot, atoms)
        !write(*,*) "DENSITY", sum(denspot%rhov) * product(denspot%dpbox%hgrids)
     end if

     !in case of GPU do not overlap density communication and projectors
     if ((ithread > 0 .or. nthread==1) .and. .not. whilepot .and. .not. GPU%OCLconv) then
        ! Only the remaining threads do computations (if active)
        !$ if (unblock_comms_den) call OMP_SET_NUM_THREADS(nthread_max-1)

        !nonlocal hamiltonian
        !$ if (get_verbose_level() > 2 .and. iproc==0 .and. unblock_comms_den)&
        !$ & print *,'NonLocalHamiltonian with nthread:, out to:' ,omp_get_max_threads(),nthread_max
        call NonLocalHamiltonianApplication(iproc,atoms,wfn%orbs%npsidim_orbs,wfn%orbs,&
             wfn%Lzd,nlpsp,wfn%psi,wfn%hpsi,energs%eproj,wfn%paw)
     end if
     !$OMP END PARALLEL !if unblock_comms_den
     !$ if (unblock_comms_den) then
     !$ call bigdft_close_nesting(nthread_max)
     !$ call f_free_ptr(denspot%rho_psi) !now the pointer can be freed
     !$ call timing(iproc,'UnBlockDen    ','OF')
     !$ end if

     ithread=0
     nthread=1

     !here the density can be mixed
     if (scf_mode .hasattr. 'MIXING') then
        if (denspot%mix%kind == AB7_MIXING_DENSITY) then
           call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,alphamix,denspot%mix,&
                denspot%rhov,itrp,denspot%dpbox%mesh%ndims(1),denspot%dpbox%mesh%ndims(2),denspot%dpbox%mesh%ndims(3),&
                atoms%astruct%cell_dim(1)*atoms%astruct%cell_dim(2)*atoms%astruct%cell_dim(3),&
                rpnrm,denspot%dpbox%nscatterarr)

           if (iproc == 0 .and. itrp > 1) then
              call yaml_newline()
              call yaml_map('itrp',itrp,fmt='(i4)')
              call yaml_map('Mixing on','Density')
              call yaml_map('RhoPot delta per volume unit',rpnrm,fmt='(1pe9.2)',&
                   label='rpnrm'//trim(adjustl(yaml_toa(itrp,fmt='(i4.4)'))))
              call yaml_newline()
              !write( *,'(1x,a,i6,2x,(1x,1pe9.2))') &
              !     &   'DENSITY iteration,Delta : (Norm 2/Volume)',itrp,rpnrm
           end if
           ! xc_init_rho should be put in the mixing routines
           denspot%rhov = abs(denspot%rhov) + 1.0d-20
        end if
     end if

     call denspot_set_rhov_status(denspot, ELECTRONIC_DENSITY, itwfn, iproc, nproc)

     !before creating the potential, save the density in the second part
     !in the case of NK SIC, so that the potential can be created afterwards
     !copy the density contiguously since the GGA is calculated inside the NK routines
     !with the savefield scheme, this can be avoided in the future
     if (wfn%SIC%approach=='NK') then !here the density should be copied somewhere else
        irhotot_add=denspot%dpbox%mesh%ndims(1)*denspot%dpbox%mesh%ndims(2)*denspot%dpbox%i3xcsh+1
        irho_add=denspot%dpbox%mesh%ndims(1)*denspot%dpbox%mesh%ndims(2)*denspot%dpbox%n3d*wfn%orbs%nspin+1
        do ispin=1,wfn%orbs%nspin
           call vcopy(denspot%dpbox%ndimpot,&
                denspot%rhov(irhotot_add),1,denspot%rhov(irho_add),1)
           irhotot_add=irhotot_add+denspot%dpbox%mesh%ndims(1)*denspot%dpbox%mesh%ndims(2)*denspot%dpbox%n3d
           irho_add=irho_add+denspot%dpbox%mesh%ndims(1)*denspot%dpbox%mesh%ndims(2)*denspot%dpbox%n3p
        end do
     end if

!!$     if(wfn%orbs%nspinor==4) then
!!$        !this wrapper can be inserted inside the XC_potential routine
!!$        call PSolverNC(atoms%astruct%geocode,'D',denspot%pkernel%mpi_env%iproc,denspot%pkernel%mpi_env%nproc,&
!!$             denspot%dpbox%mesh%ndims(1),denspot%dpbox%mesh%ndims(2),denspot%dpbox%mesh%ndims(3),&
!!$             denspot%dpbox%n3d,denspot%xc,&
!!$             denspot%dpbox%mesh%hgrids,&
!!$             denspot%rhov,denspot%pkernel%kernel,denspot%V_ext,&
!!$             energs%eh,energs%exc,energs%evxc,0.d0,.true.,4)
!!$     else

        !!           denspot%rhov denspot%rhov+2.e-7   STEFAN Goedecker

!!$     !to reproduce the previous behaviour save the density in the temporary array
!!$     if(wfn%orbs%nspinor==4) then
!!$        temp=f_malloc(src=denspot%rhov,id='temp')
!!$        m_norm=f_malloc0(size(denspot%rhov),id='m_norm')
!!$        temp2=f_malloc0(size(denspot%rhov),id='temp2')
!!$        call get_local_magnetization(denspot%dpbox%ndimrho,temp,m_norm,temp2)
!!$        call f_free(temp2)
!!$     end if

     !here we need to calculate the atomic magnetic moments of the provided density
     !the conditional should be added
     if (denspot%cfd%nat >0) then
        call atomic_magnetic_moments(denspot%dpbox%bitp,denspot%dpbox%ndimpot,atoms%astruct%nat,denspot%cfd%rxyz,denspot%cfd%radii,&
             denspot%rhov(1+denspot%dpbox%mesh%ndims(1)*denspot%dpbox%mesh%ndims(2)*denspot%dpbox%i3xcsh),&
             denspot%cfd%rho_at,denspot%cfd%m_at)
        if (nproc > 1) then
           call fmpi_allreduce(denspot%cfd%rho_at,op=FMPI_SUM,comm=bigdft_mpi%mpi_comm)
           call fmpi_allreduce(denspot%cfd%m_at,op=FMPI_SUM,comm=bigdft_mpi%mpi_comm)
        end if
        if(iproc==0) call cfd_dump_info(denspot%cfd)
!!!>
!!!>     !here there should be a call to the CFD routines
!!!>     call cfd_magnetic_field(cfd)
!!!         call constrain(denspot%cfd%m_at,denspot%cfd%m_at,denspot%cfd%B_at,denspot%cfd%nat)
!!!>
     end if

     call exchange_and_correlation(denspot%xc,denspot%dpbox,&
          denspot%rhov,energs%exc,energs%evxc,wfn%orbs%nspin,denspot%rho_C,&
          denspot%rhohat,denspot%V_XC,xcstr)

     if (denspot%cfd%nat >0) then
!!$        !here the constraining magnetic field is added on top of the local xc potential
           ! First calculate the new constraining field
           if(iproc==0) call cfd_field(denspot%cfd,iproc)
           call mpibcast(denspot%cfd%B_at,root=0,comm=bigdft_mpi%mpi_comm)!,check=.true.)
           call mpibcast(denspot%cfd%constrained_mom_err,root=0,comm=bigdft_mpi%mpi_comm)!,check=.true.)
!!$        call f_zero(denspot%cfd%B_at)
!!$        denspot%cfd%B_at(1,1)=0.5_gp
!!$        denspot%cfd%B_at(1,2)=0.5_gp
           ! Then apply the field
        call atomic_magnetic_field(denspot%dpbox%bitp,denspot%dpbox%ndimpot,atoms%astruct%nat,&
             denspot%cfd%rxyz,denspot%cfd%radii,denspot%cfd%B_at,denspot%V_XC)
!          call mpiallred(denspot%cfd%B_at,op=MPI_SUM,comm=bigdft_mpi%mpi_comm)
          if (cfd_is_converged(denspot%cfd)) then
             if (iproc==0) then
                call yaml_mapping_open('Calling ASD routines')
                call yaml_map('Moments in',denspot%cfd%m_at_ref,fmt='(f12.6)')
!!$             if(iproc==0) print *, 'Calling ASD! moments_in:'
!!$             if(iproc==0) print '(3f12.6)', denspot%cfd%m_at_ref
                !if(iproc==0) 
                call asd_wrapper(asd,denspot%cfd%m_at_ref,denspot%cfd%B_at,denspot%cfd%nat)
             end if
             call mpibcast(denspot%cfd%m_at_ref,root=0,comm=bigdft_mpi%mpi_comm)!,check=.true.)
             if (iproc==0) then
                call yaml_map('Moments out',denspot%cfd%m_at_ref,fmt='(f12.6)')
                call yaml_mapping_close()
             end if
!!$             if(iproc==0) print *, 'ASD called... moments_out' 
!!$             if(iproc==0) print '(3f12.6)', denspot%cfd%m_at_ref
          end if
     end if

!!$        call XC_potential(atoms%astruct%geocode,'D',denspot%pkernel%mpi_env%iproc,denspot%pkernel%mpi_env%nproc,&
!!$             denspot%pkernel%mpi_env%mpi_comm,&
!!$             denspot%dpbox%mesh%ndims(1),denspot%dpbox%mesh%ndims(2),denspot%dpbox%mesh%ndims(3),denspot%xc,&
!!$             denspot%dpbox%mesh%hgrids,&
!!$             denspot%rhov,energs%exc,energs%evxc,wfn%orbs%nspin,denspot%rho_C,&
!!$             denspot%rhohat,denspot%V_XC,xcstr)
        call denspot_set_rhov_status(denspot, CHARGE_DENSITY, itwfn, iproc, nproc)

        call H_potential('D',denspot%pkernel,&
             denspot%rhov,denspot%V_ext,ehart_ps,0.0_dp,.true.,&
             quiet=denspot%PSquiet,rho_ion=denspot%rho_ion) !optional argument

        if (denspot%pkernel%method /= 'VAC') then
           energs%eelec=ehart_ps
           energs%eh=0.0_gp
        else
           energs%eelec=0.0_gp
           energs%eh=ehart_ps
        end if

        !this is not true, there is also Vext
        call denspot_set_rhov_status(denspot, HARTREE_POTENTIAL, itwfn, iproc, nproc)

        !sum the two potentials in rhopot array
        !fill the other part, for spin, polarised
        ishift=denspot%dpbox%ndimpot
        if (wfn%orbs%nspin ==4) then
           call f_zero(2*denspot%dpbox%ndimpot,denspot%rhov(1+ishift))
           call f_memcpy(n=denspot%dpbox%ndimpot,src=denspot%rhov(1),&
                dest=denspot%rhov(1+3*ishift))
        else
           do ispin=2,wfn%orbs%nspin
              call f_memcpy(n=denspot%dpbox%ndimpot,src=denspot%rhov(1),&
                   dest=denspot%rhov(1+ishift))
              ishift=ishift+denspot%dpbox%ndimpot
           end do
        end if
!!$        if (wfn%orbs%nspin == 2) then
!!$           call vcopy(denspot%dpbox%ndimpot,denspot%rhov(1),1,&
!!$                denspot%rhov(1+denspot%dpbox%ndimpot),1)
!!$        end if

        !AB ALERT
        !spin up and down together with the XC part
        call axpy(denspot%dpbox%ndimpot*wfn%orbs%nspin,&
              1.0_dp,denspot%V_XC(1,1,1,1),1,&
             denspot%rhov(1),1)


        !put this term for the implementation of the (presumably incorrect) previous version
!!$        if (wfn%orbs%nspinor==4) then
!!$           !here temp is (still) the original density whereas rhov contains the full V_HXC
!!$          call get_spinorial_potential(denspot%dpbox%ndimpot,denspot%rhov,m_norm,temp)
!!$          call f_memcpy(src=temp,dest=denspot%rhov)
!!$          call f_free(temp,m_norm)
!!$       end if


!!$        !here a external potential with spinorial indices can be added
!!$     end if

     !this part has to be replaced by the updatepotential routine
!     call updatePotential(wfn%orbs%nspinor,denspot,energs)

     !here the potential can be mixed
     if (scf_mode .hasattr. 'MIXING') then
        if (denspot%mix%kind == AB7_MIXING_POTENTIAL) then
           call mix_rhopot(iproc,nproc,denspot%mix%nfft*denspot%mix%nspden,alphamix,denspot%mix,&
                denspot%rhov,itrp,denspot%dpbox%mesh%ndims(1),denspot%dpbox%mesh%ndims(2),denspot%dpbox%mesh%ndims(3),&
                atoms%astruct%cell_dim(1)*atoms%astruct%cell_dim(2)*atoms%astruct%cell_dim(3),&!volume should be used
                rpnrm,denspot%dpbox%nscatterarr)
           if (iproc == 0 .and. itrp > 1) then
              call yaml_newline()
              call yaml_map('itrp',itrp,fmt='(i4)')
              call yaml_map('Mixing on','Potential')
              call yaml_map('RhoPot delta per volume unit',rpnrm,fmt='(1pe9.2)',&
                   label='rpnrm'//trim(adjustl(yaml_toa(itrp,fmt='(i4.4)'))))
              call yaml_newline()
              !write( *,'(1x,a,i6,2x,(1x,1pe9.2))') &
              !     &   'POTENTIAL iteration,Delta P (Norm 2/Volume)',itrp,rpnrm
           end if
        end if
     end if
     call denspot_set_rhov_status(denspot, KS_POTENTIAL, itwfn, iproc, nproc)

     if (savefields) then
        if (associated(denspot%rho_work)) then
           call f_err_throw('The reference potential should be empty to correct the hamiltonian!',&
                err_name='BIGDFT_RUNTIME_ERROR')
           return
        end if

        denspot%rho_work = f_malloc_ptr(denspot%dpbox%ndimpot*denspot%dpbox%nrhodim,id='denspot%rho_work')
        call vcopy(denspot%dpbox%ndimpot*denspot%dpbox%nrhodim,denspot%rhov(1),1,&
             denspot%rho_work(1),1)
     end if

     !if the hamiltonian should have to be just updated, perform the difference between potentials
     if (correcth==2) then
        if (.not. associated(denspot%rho_work)) then
           write(*,*)'ERROR: need a reference potential to correct the hamiltonian!'
           stop
        end if
        call yaml_newline()
        call yaml_comment('Calculating potential delta')
        !subtract the previous potential from the new one
        call axpy(denspot%dpbox%ndimpot*denspot%dpbox%nrhodim,&
             -1.0_dp,denspot%rho_work(1),1,denspot%rhov(1),1)

        !deallocation should be deplaced
        call f_free_ptr(denspot%rho_work)
        nullify(denspot%rho_work)
     end if

     if (wfn%paw%usepaw) then
        call paw_compute_dij(wfn%paw, atoms, denspot, denspot%V_XC, &
             & energs%epaw, energs%epawdc, compch_sph)
     end if
  end if !if (scf)

  !debug
  !call MPI_BARRIER(MPI_COMM_WORLD,i_stat)
  !write(*,*)'hpsiortho l325,erase me'
  !ifile=iproc+140
  !jj=size(wfn%psi)
  !ii=size(wfn%hpsi)
  !write(ifile,*)'# ',jj,ii
  !do ii=1,jj
  ! write(ifile,'(2f20.8)')wfn%psi(ii),wfn%hpsi(ii)
  !end do
  !ifile=iproc+150
  !jj=size(denspot%rhov)
  !write(ifile,*)'# ',jj
  !do ii=1,jj
  ! write(ifile,'(2f20.8)')denspot%rhov(ii)
  !end do
  !call MPI_BARRIER(MPI_COMM_WORLD,i_stat)
  !end debug

  ! self-consistent case or not: rhov should be the total potential
  if (denspot%rhov_is /= KS_POTENTIAL) &
       call f_err_throw('psitohpsi: KS_potential not available, control the operations on rhov',&
       err_name='BIGDFT_RUNTIME_ERROR')
  

  !temporary, to be corrected with comms structure
  if (wfn%exctxpar == 'OP2P') energs%eexctX = UNINITIALIZED(1.0_gp)

  !initialize nested approach
  !this has always to be done for using OMP parallelization in the
  !projector case
  !if nesting is not supported, bigdft_nesting routine should not be called
  !$ if (unblock_comms_pot) then
  !$ call timing(iproc,'UnBlockPot    ','ON')
  !$ call bigdft_open_nesting(2)
  !$ end if
  !print *,'how many threads ?' ,nthread_max
  !$OMP PARALLEL IF (unblock_comms_pot) DEFAULT(shared), PRIVATE(ithread,nthread)
  !$ ithread=omp_get_thread_num()
  !$ nthread=omp_get_num_threads() !this should be 2 if active
  !print *,'hello, I am thread no.',ithread,' out of',nthread,'of iproc', iproc
  ! thread 0 does mpi communication
  if (ithread == 0) then
     !$ if (unblock_comms_pot) call OMP_SET_NUM_THREADS(1)
     call full_local_potential(iproc,nproc,wfn%orbs,wfn%Lzd,linflag,&
          denspot%dpbox,denspot%xc,denspot%rhov,denspot%pot_work)
     !write(*,*) 'node:', iproc, ', thread:', ithread, 'mpi communication finished!!'
  end if

  if ((ithread > 0 .or. nthread==1) .and. whilepot .and. .not. GPU%OCLconv) then
     ! Only the remaining threads do computations (if active)
     !$ if (unblock_comms_pot) call OMP_SET_NUM_THREADS(nthread_max-1)

     !nonlocal hamiltonian
     !$ if (get_verbose_level() > 2 .and. iproc==0 .and. unblock_comms_pot)&
     !$ & call yaml_map('NonLocalHamiltonian with nthread out to',[omp_get_max_threads(),nthread_max])
     call NonLocalHamiltonianApplication(iproc,atoms,wfn%orbs%npsidim_orbs,wfn%orbs,&
          wfn%Lzd,nlpsp,wfn%psi,wfn%hpsi,energs%eproj,wfn%paw)
  end if
  !$OMP END PARALLEL !if unblock_comms_pot
  !$ if (unblock_comms_pot) then
  !$ call bigdft_close_nesting(nthread_max)
  !$ call timing(iproc,'UnBlockPot    ','OF')
  !$ end if

  ithread=0
  nthread=1

  !here exctxpar might be passed
  !choose to just add the potential if needed
  call LocalHamiltonianApplication(iproc,nproc,atoms,wfn%orbs%npsidim_orbs,wfn%orbs,&
       wfn%Lzd,wfn%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,wfn%psi,wfn%hpsi,&
       energs,wfn%SIC,GPU,correcth,denspot%xc,pkernel=denspot%pkernelseq)

  !in the case of OCL GPU the nonlocal hamiltonian can run after the local hamiltonian to overlap GPU-CPU computation
  if (GPU%OCLconv) then
     call NonLocalHamiltonianApplication(iproc,atoms,wfn%orbs%npsidim_orbs,wfn%orbs,&
          wfn%Lzd,nlpsp,wfn%psi,wfn%hpsi,energs%eproj,wfn%paw)
  end if

  call SynchronizeHamiltonianApplication(nproc,wfn%orbs%npsidim_orbs,wfn%orbs,wfn%Lzd,&
       & GPU,denspot%xc,wfn%hpsi,&
       energs)!%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)

  !here we can reduce and output the density matrix if required
  if (associated(nlpsp%gamma_mmp) .and. nproc > 1) &
       call fmpi_allreduce(nlpsp%gamma_mmp,op=FMPI_SUM,comm=bigdft_mpi%mpi_comm)

  if (iproc==0 .and. get_verbose_level() > 1) call write_atomic_density_matrix(wfn%orbs%nspin,atoms%astruct,nlpsp)

  !here we might rework the value of gamma in case we would like to apply some extra
  !term in the following iteration
  !in case there is any target occupancy then calculate the delta
  if (nlpsp%apply_gamma_target) then
     call atomic_density_matrix_delta(iproc==0,wfn%orbs%nspin,atoms%astruct,nlpsp,&
          atoms%gamma_targets)
  end if


!!$  if (iproc ==0) then
!!$     !compute the proper hartree energy in the case of a cavity calculation
!!$     !energs%eh=energs%epot-energs%eh-energs%evxc
!!$     call yaml_map('Potential',energs%epot)
!!$     call yaml_map('Ionic',energs%eion)
!!$     call yaml_map('Evxc',energs%evxc)
!!$     call yaml_map('Old EH',energs%eh)
!!$     call yaml_map('rho Vext+static',energs%epot-2*energs%eh-energs%evxc)
!!$     call yaml_map('EHtot',energs%epot-energs%eh-energs%evxc)
!!$  end if

  ! Emit that hpsi are ready.
  if (wfn%c_obj /= 0) then
     call kswfn_emit_psi(wfn, itwfn, 1, iproc, nproc)
  end if


  !deallocate potential
  call free_full_potential(denspot%dpbox%mpi_env%nproc,linflag,denspot%xc,denspot%pot_work)
  !----
  if (iproc==0 .and. get_verbose_level() > 0) then
     if (correcth==2) then
        call yaml_map('Hamiltonian Applied',.false.)
        call yaml_comment('Only local potential')
     else
        call yaml_map('Hamiltonian Applied',.true.)
     end if
  end if

  call f_release_routine()
end subroutine psitohpsi


!> Application of the Full Hamiltonian
subroutine FullHamiltonianApplication(iproc,nproc,at,orbs,&
     Lzd,nlpsp,confdatarr,ngatherarr,pot,psi,hpsi,paw,&
     energs,SIC,GPU,xc,pkernel,orbsocc,psirocc)
  use module_base
  use module_types
  use module_interfaces, only: LocalHamiltonianApplication, SynchronizeHamiltonianApplication
  use module_xc
  use public_enums, only: PSPCODE_PAW
  use yaml_output
  use locreg_operations, only: confpot_data
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: Lzd
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(SIC_data), intent(in) :: SIC
  type(xc_info), intent(in) :: xc
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
  type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
  !Local variables
  !real(wp), dimension(lzd%ndimpotisf) :: pot
  real(wp), dimension(:),pointer :: pot
  type(energy_terms), intent(inout) :: energs
  !real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum,evsic
  real(wp), target, dimension(max(1,orbs%npsidim_orbs)), intent(out) :: hpsi
  type(GPU_pointers), intent(inout) :: GPU
  type(coulomb_operator), intent(inout), optional :: pkernel
  type(orbitals_data), intent(in), optional :: orbsocc
  real(wp), dimension(:), pointer, optional :: psirocc
  !PAW variables:
  type(paw_objects),intent(inout)::paw

  call f_routine(id='FullHamiltonianApplication')

  !put to zero hpsi array (now important since any of the pieces of the hamiltonian is accumulating)
  if (orbs%npsidim_orbs > 0) call f_zero(orbs%npsidim_orbs,hpsi(1))

  !write(*,*) 'lzd%ndimpotisf', lzd%ndimpotisf
  !do i=1,lzd%ndimpotisf
  !    write(210,*) pot(i)
  !end do

  if (.not. present(pkernel)) then
     call LocalHamiltonianApplication(iproc,nproc,at,orbs%npsidim_orbs,orbs,&
          Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
          energs,SIC,GPU,1,xc)
  else if (present(pkernel) .and. .not. present(orbsocc)) then
     call LocalHamiltonianApplication(iproc,nproc,at,orbs%npsidim_orbs,orbs,&
          Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
          energs,SIC,GPU,1,xc,pkernel=pkernel)
  else if (present(pkernel) .and. present(orbsocc) .and. present(psirocc)) then
     call LocalHamiltonianApplication(iproc,nproc,at,orbs%npsidim_orbs,orbs,&
          Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
          energs,SIC,GPU,1,xc,pkernel,orbsocc,psirocc)
  else
     call f_err_throw('HamiltonianApplication, argument error')
  end if

  !these two sections have to be inverted to profit of overlapping in GPU accelerated case
  call NonLocalHamiltonianApplication(iproc,at,orbs%npsidim_orbs,orbs,&
       Lzd,nlpsp,psi,hpsi,energs%eproj,paw)

  call SynchronizeHamiltonianApplication(nproc,orbs%npsidim_orbs,orbs,Lzd,GPU,xc,hpsi,&
       energs)!%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)

  !to be adjusted
!!$  if (trim(denspot%pkernel%method) /= 'VAC') then
  !     energs%eh=energs%epot-energs%eh-energs%evxc
!!$  end if

  call f_release_routine()
END SUBROUTINE FullHamiltonianApplication


!> Application of the Local Hamiltonian
subroutine LocalHamiltonianApplication(iproc,nproc,at,npsidim_orbs,orbs,&
     Lzd,confdatarr,ngatherarr,pot,psi,hpsi,&
     energs,SIC,GPU,PotOrKin,xc,pkernel,orbsocc,psirocc,dpbox,potential,comgp,hpsi_noconf,econf)
  use module_base
  use module_dpbox, only: denspot_distribution
  use module_types
  use module_xc
  use module_interfaces, only: NK_SIC_potential, psi_to_vlocpsi
  use orbitalbasis
  use yaml_output
  use communications_base, only: p2pComms
  use locreg_operations
  use overlap_point_to_point
  implicit none
  !logical, intent(in) :: onlypot !< if true, only the potential operator is applied
  integer, intent(in) :: PotOrKin
  integer, intent(in) :: iproc,nproc,npsidim_orbs
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  type(SIC_data), intent(in) :: SIC
  type(xc_info), intent(in) :: xc
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(wp), dimension(npsidim_orbs), intent(in) :: psi
  type(confpot_data), dimension(orbs%norbp) :: confdatarr
  real(wp), dimension(:), pointer :: pot
  !real(wp), dimension(*) :: pot
  type(energy_terms), intent(inout) :: energs
  real(wp), target, dimension(max(1,npsidim_orbs)), intent(inout) :: hpsi
  type(GPU_pointers), intent(inout) :: GPU
  type(coulomb_operator), intent(inout), optional :: pkernel
  type(orbitals_data), intent(in), optional :: orbsocc
  real(wp), dimension(:), pointer, optional :: psirocc
  type(denspot_distribution),intent(in),optional :: dpbox
  !!real(wp), dimension(max(dpbox%ndimrhopot,orbs%nspin)), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
  real(wp), dimension(*), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
  type(p2pComms),intent(inout), optional:: comgp
  real(wp), target, dimension(max(1,npsidim_orbs)), intent(inout),optional :: hpsi_noconf
  real(gp),intent(out),optional :: econf
  !local variables
  character(len=*), parameter :: subname='LocalHamiltonianApplication'
  logical :: exctX,op2p_flag, symmetric
  integer :: n3p,ispot,ipotmethod,ngroup,prc,isorb,jproc,ndim,norbp
  integer :: igpu,gpudirect,i_stat,nsize,nspinor
  real(gp) :: evsic_tmp, ekin, epot,sfac
  real(f_double) :: tel,trm
  type(coulomb_operator) :: pkernelSIC
  type(ket) :: psi_it
  type(orbital_basis) :: psi_ob
  type(workarr_locham) :: wrk_lh
  type(workarr_sumrho) :: w
  type(OP2P_data) :: OP2P
  type(OP2P_iterator) :: iter
  integer, dimension(:,:), allocatable :: nobj_par
  real(wp), dimension(:,:), allocatable :: vsicpsir
  real(wp), dimension(:,:), allocatable :: psir
  !real(wp), dimension(:,:), allocatable :: vpsi_tmp
  real(wp), dimension(:), allocatable :: rp_ij
  real(wp), dimension(:), pointer :: hpsi_ptr
  real(gp) :: eSIC_DCi,fi

  call f_routine(id=subname)

  ! local potential and kinetic energy for all orbitals belonging to iproc
  !if (iproc==0 .and. get_verbose_level() > 1) then
  !call yaml_comment('Hamiltonian application, ',advance='no')
  !write(*,'(1x,a)',advance='no')&
  !     'Hamiltonian application...'
  !end if

  !initialise exact exchange energy
  op2p_flag=(energs%eexctX == UNINITIALIZED(1.0_gp))
  energs%eexctX=0.0_gp
  energs%evsic=0.0_gp
  evsic_tmp=0.0_gp

  exctX = xc_exctXfac(xc) /= 0.0_gp


  ispot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%nspin+1

  !potential method
  !traditional case
  ipotmethod=0
  if (exctX) ipotmethod=1


  !the PZ-SIC correction does not makes sense for virtual orbitals procedure
  !if alphaSIC is zero no SIC correction
  if (SIC%approach == 'PZ' .and. .not. present(orbsocc) .and. SIC%alpha /= 0.0_gp ) ipotmethod=2
  if (SIC%approach == 'NK' .and. SIC%alpha /= 0.0_gp) ipotmethod=3

  !the poisson kernel should be present and associated in the case of SIC
  if ((ipotmethod /= 0) .and. present(pkernel)) then
     if (.not. associated(pkernel%kernel)) then
        if (iproc ==0) write(*,*)&
             &   'ERROR(LocalHamiltonianApplication): Poisson Kernel must be associated in SIC case'
        stop
     end if
  end if

  !associate the poisson kernel pointer in case of SIC
  if (ipotmethod == 2 .or. ipotmethod == 3) then
     pkernelSIC = pkernel
  else
     nullify(pkernelSIC%kernel)
  end if


  !fill the rest of the potential with the exact-exchange terms
  if (ipotmethod==1) then
     n3p=ngatherarr(iproc,1)/(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i)
     !exact exchange for virtual orbitals (needs psirocc)

     !here we have to add the round part
     if (present(psirocc) .and. present(orbsocc)) then
        call exact_exchange_potential_virt(iproc,nproc,at%astruct%geocode,orbs%nspin,&
             Lzd%Glr,orbsocc,orbs,ngatherarr(0,1),n3p,&
             0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
             pkernel,psirocc,psi,pot(ispot))
        energs%eexctX = 0._gp
     else
        !here the condition for the scheme should be chosen
        if (.not. op2p_flag) then
           call exact_exchange_potential(iproc,nproc,at%astruct%geocode,xc,orbs%nspin,&
                Lzd%Glr,orbs,ngatherarr(0,1),n3p,&
                0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
                pkernel,psi,pot(ispot),energs%eexctX)
        else

!!$            call exact_exchange_potential_round(iproc,nproc,xc,orbs%nspin,Lzd%Glr,orbs,&
!!$                0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
!!$                pkernel,psi,pot(ispot),energs%eexctX)

           !the psi should be transformed in real space, do it within the orbital basis iterators
           psir = f_malloc0([Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i, orbs%norbp],id='psir')
           !initialize the orbital basis object, for psi and hpsi
           call orbital_basis_associate(psi_ob,orbs=orbs,phis_wvl=psi,Lzd=Lzd,id='LocalHamiltonianApplication')
           !iterate over the orbital_basis
           psi_it=orbital_basis_iterator(psi_ob)
           do while(ket_next_locreg(psi_it))
              call initialize_work_arrays_sumrho(psi_it%lr,.true.,w)
              do while(ket_next(psi_it,ilr=psi_it%ilr))
                 call daub_to_isf(psi_it%lr,w,psi_it%phi_wvl,psir(1,psi_it%iorbp))
              end do
              !deallocations of work arrays
              call deallocate_work_arrays_sumrho(w)
           end do
           call orbital_basis_release(psi_ob)

           !starting new approach for the exact exchange
           !vpsi_tmp = f_malloc0([Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i, orbs%norbp],id='vpsi_tmp')
           if (orbs%nspin==2) then
              sfac=1.0_gp
              ngroup=2
           else
              sfac=0.5_gp
              ngroup=1
           end if
           !construct the OP2P scheme and test it
           !use temporaryly tyhe nvrct_par array
           nobj_par = f_malloc0((/ 0.to.nproc-1, 1.to.ngroup /),id='nobj_par')

           call fill_nobj_par_for_OP2P(nproc,ngroup,orbs,nobj_par)

           ndim=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i

           symmetric=.true.
           call f_zero(ndim*orbs%norbp,pot(ispot))
           !if (iproc==0) call yaml_map('Orbital repartition',nobj_par)
           !Comment because the test is not fully operational.
           !!$ call OP2P_unitary_test(bigdft_mpi%mpi_comm,iproc,nproc,ngroup,ndim,nobj_par,symmetric)
           if(pkernel%igpu==1 .and. pkernel%initCufftPlan==0) then
              igpu=0
              if (iproc==0) call yaml_warning("not enough memory to allocate cuFFT plans on the GPU - no GPUDirect either")
           else
              igpu=pkernel%igpu
           end if

           gpudirect=0
           !check that we did not deactivate gpudirect manually
           if(igpu==1 .and. pkernel%use_gpu_direct) gpudirect=1

           call initialize_OP2P_data(OP2P,bigdft_mpi%mpi_comm,iproc,nproc,ngroup,ndim,nobj_par,gpudirect,symmetric)

           !initialization deactivates gpudirect if not enough memory (these parts should go via accessors in the PS_set_options somehow)
           if(gpudirect==1 .and. OP2P%gpudirect==1) pkernel%stay_on_gpu=1

           !allocate work array for the internal exctx calculation
           rp_ij = f_malloc(ndim,id='rp_ij')
           energs%eexctX=0.0_gp
           !initialize the OP2P descriptor for the communication
           call set_OP2P_iterator(iproc,OP2P,iter,orbs%norbp,psir,pot(ispot:ispot+OP2P%ndim*orbs%norbp-1))!vpsi_tmp)
           !main loop
           if (iproc == 0) call yaml_newline()
           OP2P_exctx_loop: do
              call OP2P_communication_step(iproc,OP2P,iter)
              if(igpu==1) call synchronize()
              if (iter%event == OP2P_EXIT) exit OP2P_exctx_loop
              call internal_calculation_exctx(iter%istep,sfac,pkernel,orbs%norb,orbs%occup,orbs%spinsgn,&
                   iter%remote_result,iter%nloc_i,iter%nloc_j,iter%isloc_i,iter%isloc_j,&
                   iter%phi_i,iter%phi_j,energs%eexctX,rp_ij)
              if (iproc == 0) then
                 call OP2P_info(iter,OP2P,prc,tel,trm)
                 call yaml_comment('Exact exchange calculation: '+prc**'(i3)'+&
                      '%; Time (s): Elapsed='+tel**'(1pg12.2)'&
                      +', Remaining='+trm**'(1pg12.2)')
              end if
           end do OP2P_exctx_loop

           !we need to get the result back from the card (and synchronize with it, finally)
           if(pkernel%igpu==1 .and. pkernel%stay_on_gpu==1) then
              call get_gpu_data(1, energs%eexctX, pkernel%w%eexctX_GPU )
              call cudamemset(pkernel%w%eexctX_GPU,0,1,i_stat)
              if (i_stat /= 0) call f_err_throw('error cudamalloc eexctX_GPU (GPU out of memory ?) ')
              pkernel%stay_on_gpu=0
           end if
           call free_OP2P_data(OP2P)
           if (nproc>1) call fmpi_allreduce(energs%eexctX,1,FMPI_SUM,comm=bigdft_mpi%mpi_comm)
           !the exact exchange energy is half the Hartree energy (which already has another half)
           energs%eexctX=-xc_exctXfac(xc)*energs%eexctX
           if (iproc == 0) call yaml_map('Exact Exchange Energy',energs%eexctX,fmt='(1pe18.11)')
           !call f_memcpy(n=ndim*orbs%norbp,src=vpsi_tmp(1,1),dest=pot(ispot))
           call f_free(nobj_par)
           call f_free(rp_ij)
           !call f_free(vpsi_tmp)


!!$            call exact_exchange_potential_round_clean(iproc,nproc,xc,orbs%nspin,&
!!$                 Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,orbs,&
!!$                 !0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
!!$                 pkernel,psir,pot(ispot),energs%eexctX)

           call f_free(psir)


           !call exact_exchange_potential_op2p(iproc,nproc,xc,Lzd%Glr,orbs,pkernel,psi,pot(ispot),energs%eexctX)

        end if
     end if
     !print *,'iproc,eexctX',iproc,eexctX
  else if (ipotmethod==3) then
     !put fref=1/2 for the moment
     if (present(orbsocc) .and. present(psirocc)) then
        call NK_SIC_potential(Lzd%Glr,orbs,xc,SIC%fref,&
             0.5_gp*Lzd%hgrids,&
             pkernelSIC,psi,pot(ispot),evsic_tmp,&
             potandrho=psirocc)
     else
        call NK_SIC_potential(Lzd%Glr,orbs,xc,SIC%fref,&
             0.5_gp*Lzd%hgrids,&
             pkernelSIC,psi,pot(ispot),evsic_tmp)
     end if
  end if

  !GPU are supported only for ipotmethod=0
  if (GPU%OCLconv .and. ipotmethod /=0) then
     call f_err_throw('ERROR(HamiltonianApplication): '//&
          'Accelerated hamiltonian are possible only with ipotmethod==0)',&
          err_name='BIGDFT_RUNTIME_ERROR')
  end if

!!!

  !apply the local hamiltonian for each of the orbitals
  !given to each processor
  !pot=0.d0
  !psi=1.d0
  !switch between GPU/CPU treatment
  !  do i=1,(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp
  !       call random_number(psi(i))
  !  end do
  if(GPU%OCLconv) then! needed also in the non_ASYNC since now NlPSP is before .and. ASYNCconv)) then
     GPU%hpsi_ASYNC = f_malloc_ptr(max(1, (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp),id='GPU%hpsi_ASYNC')
     !      call f_zero((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp,GPU%hpsi_ASYNC(1))!hpsi(1))
     !else if (GPU%OCLconv) then
     !   GPU%hpsi_ASYNC => hpsi
  end if
  if (GPU%OCLconv) then
     !pin potential
     !call timing(iproc,'ApplyLocPotKin','ON')
     call local_hamiltonian_OCL(orbs,Lzd%Glr,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
          orbs%nspin,pot,psi,GPU%hpsi_ASYNC,energs%ekin,energs%epot,GPU)
     !call timing(iproc,'ApplyLocPotKin','OF')
  else

     !local hamiltonian application for different methods
     !print *,'here',ipotmethod,associated(pkernelSIC)
     if (PotOrKin==1) then ! both
        call timing(iproc,'ApplyLocPotKin','ON')
!!!here we can branch into the new ket-based application of the hamiltonian
        !initialize the orbital basis object, for psi and hpsi
        call orbital_basis_associate(psi_ob,orbs=orbs,confdatarr=confdatarr,&
             phis_wvl=psi,Lzd=Lzd,id='LochamPotKin')

        energs%ekin=0.0_gp
        energs%epot=0.0_gp


!!$         call test_iterator(psi_ob)
!!$!         stop

        !iterate over the orbital_basis
        psi_it=orbital_basis_iterator(psi_ob)
        ! Get the maximal size of the work array psir
        nsize = 0
        nspinor = 0
        loop_lr_size: do while(ket_next_locreg(psi_it))
           nsize = max(nsize,psi_it%lr%d%n1i*psi_it%lr%d%n2i*psi_it%lr%d%n3i)
           nspinor = max(nspinor,psi_it%nspinor)
        end do loop_lr_size
        psir = f_malloc([nsize,nspinor],id='psir')
        psi_it=orbital_basis_iterator(psi_ob)
        ! print *,'orbs0',psi_it%iorb,psi_it%ilr
        loop_lr: do while(ket_next_locreg(psi_it))
           ! print *,'orbs1',psi_it%iorb,psi_it%ilr,psi_it%nspinor,associated(psi_it%lr)
           !psir = f_malloc0([psi_it%lr%d%n1i*psi_it%lr%d%n2i*psi_it%lr%d%n3i,psi_it%nspinor],id='psir')
           call f_zero(psir)
           call initialize_work_arrays_locham(psi_it%lr,psi_it%nspinor,.true.,wrk_lh)
           ! wavefunction after application of the self-interaction potential
           if (ipotmethod == 2 .or. ipotmethod == 3) then
              vsicpsir = f_malloc([psi_it%lr%d%n1i*psi_it%lr%d%n2i*psi_it%lr%d%n3i,psi_it%nspinor],id='vsicpsir')
           end if
           !print *,'orbs',psi_it%iorb,psi_it%ilr

           loop_psi_lr: do while(ket_next(psi_it,ilr=psi_it%ilr))
              fi=psi_it%kwgt*psi_it%occup
              hpsi_ptr => ob_ket_map(hpsi,psi_it)
              call local_hamiltonian_ket(psi_it,Lzd%hgrids,ipotmethod,xc,&
                   pkernelSIC,wrk_lh,psir,vsicpsir,hpsi_ptr,pot,eSIC_DCi,SIC%alpha,epot,ekin)
              energs%ekin=energs%ekin+fi*ekin
!!if(iproc==0)then
!!write(*,'(a,3(x,es14.7))')'Ekin per orbital',fi,ekin,fi*ekin
!!endif
              energs%epot=energs%epot+fi*epot
              energs%evsic=energs%evsic+SIC%alpha*eSIC_DCi
              !print *,'test',psi_it%iorb,sum(hpsi_ptr)
              !  print *,'orbs2',psi_it%iorbp,psi_it%iorb,psi_it%ikpt,psi_it%kwgt,psi_it%occup,epot,ekin,psi_it%ispsi,psi_it%nspinor
           end do loop_psi_lr
           !deallocations of work arrays
           !call f_free(psir)
           if (ipotmethod == 2 .or. ipotmethod ==3) then
              call f_free(vsicpsir)
           end if
           call deallocate_work_arrays_locham(wrk_lh)
        end do loop_lr
        call f_free(psir)
        !call mpibarrier()
        !stop
        call orbital_basis_release(psi_ob)

!!$         if(present(dpbox) .and. present(potential) .and. present(comgp)) then
!!$            call local_hamiltonian(iproc,nproc,npsidim_orbs,orbs,Lzd,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
!!$                 ipotmethod,confdatarr,pot,psi,hpsi,pkernelSIC,&
!!$                 xc,SIC%alpha,energs%ekin,energs%epot,energs%evsic,&
!!$                 dpbox,potential,comgp)
!!$         else
!!$            call local_hamiltonian(iproc,nproc,npsidim_orbs,orbs,Lzd,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
!!$                 ipotmethod,confdatarr,pot,psi,hpsi,pkernelSIC,&
!!$                 xc,SIC%alpha,energs%ekin,energs%epot,energs%evsic)
!!$         end if
        call timing(iproc,'ApplyLocPotKin','OF')

     else if (PotOrKin==2) then !only pot
        call timing(iproc,'ApplyLocPot','ON')
        if (.not.present(comgp)) then
            call f_err_throw('comgp not present')
        end if
        if (present(hpsi_noconf)) then
           if (.not.present(econf)) then
              call f_err_throw('econf must be present when hpsi_noconf is present')
           end if
           call psi_to_vlocpsi(iproc,npsidim_orbs,orbs,Lzd,&
                ipotmethod,confdatarr,pot,psi,hpsi,pkernelSIC,&
                xc,SIC%alpha,energs%epot,energs%evsic,comgp,hpsi_noconf,econf)
        else
           call psi_to_vlocpsi(iproc,npsidim_orbs,orbs,Lzd,&
                ipotmethod,confdatarr,pot,psi,hpsi,pkernelSIC,&
                xc,SIC%alpha,energs%epot,energs%evsic,comgp)
        end if
        call timing(iproc,'ApplyLocPot','OF')
     else if (PotOrKin==3) then !only kin
        call timing(iproc,'ApplyLocKin','ON')
        call psi_to_kinpsi(iproc,npsidim_orbs,orbs,lzd,psi,hpsi,energs%ekin)
        call timing(iproc,'ApplyLocKin','OF')
     end if

     !sum the external and the BS double counting terms
     energs%evsic=energs%evsic-SIC%alpha*evsic_tmp
  end if

  if (ipotmethod == 2 .or. ipotmethod==3) then
     nullify(pkernelSIC%kernel)
  end if

  call f_release_routine()

END SUBROUTINE LocalHamiltonianApplication

subroutine fill_nobj_par_for_OP2P(nproc,ngroup,orbs,nobj_par)
  use module_types, only: orbitals_data
  implicit none
  integer, intent(in) :: nproc,ngroup
  type(orbitals_data), intent(in) :: orbs
  integer, dimension(0:nproc-1,ngroup), intent(out) :: nobj_par
  !local variables
  integer :: isorb,jproc,norbp

  isorb=0
  do jproc=0,nproc-1
     norbp=orbs%norb_par(jproc,0)
     !transition region
     if (isorb+norbp > orbs%norbu .and. isorb < orbs%norbu) then
        nobj_par(jproc,1)=orbs%norbu-isorb
        if (ngroup==2) nobj_par(jproc,2)=isorb+norbp-orbs%norbu
     else if (isorb >= orbs%norbu .and. ngroup==2) then
        nobj_par(jproc,2)=norbp
     else
        nobj_par(jproc,1)=norbp
     end if
     isorb=isorb+norbp
  end do
end subroutine fill_nobj_par_for_OP2P

subroutine NonLocalHamiltonianApplication(iproc,at,npsidim_orbs,orbs,&
     Lzd,nl,psi,hpsi,eproj_sum,paw)
  use module_base
  use module_types
  use yaml_output
  !  use module_interfaces
  use psp_projectors, only: projector_has_overlap,&
       get_proj_locreg,hgh_psp_application
  use public_enums, only: PSPCODE_PAW
  use module_atoms
  use orbitalbasis
  use ao_inguess, only: lmax_ao
  use pseudopotentials, only: atomic_proj_coeff,nullify_atomic_proj_coeff,f_free_prj_ptr
  implicit none
  integer, intent(in) :: iproc, npsidim_orbs
  type(atoms_data), intent(in) :: at
  type(orbitals_data),  intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  type(DFT_PSP_projectors), intent(inout) :: nl
  real(wp), dimension(npsidim_orbs), intent(in) :: psi
  real(wp), dimension(npsidim_orbs), intent(inout) :: hpsi
  type(paw_objects),intent(inout)::paw
  real(gp), intent(out) :: eproj_sum
  !local variables
  character(len=*), parameter :: subname='NonLocalHamiltonianApplication'
  logical :: overlap
  integer :: istart_ck,nwarnings
  integer :: iproj,istart_c,mproj,iilr
  type(ket) :: psi_it
  type(orbital_basis) :: psi_ob
  type(atoms_iterator) :: atit
  real(wp), dimension(:), pointer :: hpsi_ptr,spsi_ptr

  !integer :: ierr
  !real(kind=4) :: tr0, tr1, t0, t1
  !real(kind=8) :: time0, time1, time2, time3, time4, time5, ttime
  !real(kind=8), dimension(0:4) :: times

  call f_routine(id=subname)

  eproj_sum=0.0_gp

  ! apply all PSP projectors for all orbitals belonging to iproc
  call timing(iproc,'ApplyProj     ','ON')

  !initialize the orbital basis object, for psi and hpsi
  call orbital_basis_associate(psi_ob,orbs=orbs,phis_wvl=psi,Lzd=Lzd,id='nonlocalham')
  !should we calculate the density matrix we have to zero it
  if (associated(nl%iagamma)) call f_zero(nl%gamma_mmp)
  !here we might rework the value of gamma in case we would like to apply some extra

  nwarnings=0
  if(paw%usepaw) call f_zero(orbs%npsidim_orbs, paw%spsi(1))

  !here the localisation region should be changed, temporary only for cubic approach

  !apply the projectors following the strategy (On-the-fly calculation or not)

  !apply the projectors  k-point of the processor
  !starting k-point
  istart_ck=1
  !iterate over the orbital_basis
  psi_it=orbital_basis_iterator(psi_ob)

  loop_kpt: do while(ket_next_kpt(psi_it))
     loop_lr: do while(ket_next_locreg(psi_it,ikpt=psi_it%ikpt))
        if (nl%on_the_fly) then
           iproj=0
           atit = atoms_iter(at%astruct)
           loop_atoms: do while(atoms_iter_next(atit))
              ! Check whether the projectors of this atom have an overlap with locreg ilr
              overlap = projector_has_overlap(atit%iat, psi_it%ilr,psi_it%lr,lzd%glr, nl)
              if(.not. overlap) cycle loop_atoms
              istart_c=1
              mproj=nl%pspd(atit%iat)%mproj
              call atom_projector(nl,atit%ityp,atit%iat, atit%name, &
                   at%astruct%geocode, 0, Lzd%Glr,&
                   Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3), &
                   psi_it%kpoint(1),psi_it%kpoint(2),psi_it%kpoint(3),&
                   istart_c, iproj, nwarnings)

              iilr=get_proj_locreg(nl,atit%iat,psi_it%ilr)
              loop_psi_kpt: do while(ket_next(psi_it,ikpt=psi_it%ikpt,ilr=psi_it%ilr))
                 istart_c=1
                 call nl_psp_application()
              end do loop_psi_kpt
           end do loop_atoms
        else
           loop_psi_kpt2: do while(ket_next(psi_it,ikpt=psi_it%ikpt,ilr=psi_it%ilr))
              atit = atoms_iter(at%astruct)
              istart_c=istart_ck !TO BE CHANGED IN ONCE-AND-FOR-ALL
              loop_atoms2: do while(atoms_iter_next(atit))
                 overlap = projector_has_overlap(atit%iat,psi_it%ilr,psi_it%lr, lzd%glr, nl)
                 if(.not. overlap) cycle loop_atoms2
                 mproj=nl%pspd(atit%iat)%mproj

                 iilr=get_proj_locreg(nl,atit%iat,psi_it%ilr)
                 call nl_psp_application()
              end do loop_atoms2
           end do loop_psi_kpt2
        end if
     end do loop_lr
     istart_ck=istart_c !only needed for the once-and-for-all
  end do loop_kpt

  call orbital_basis_release(psi_ob)

  if (.not. nl%on_the_fly .and. Lzd%nlr==1) then !TO BE REMOVED WITH NEW PROJECTOR APPLICATION
     if (istart_ck-1 /= nl%nprojel .and. orbs%norbp>0) then
        call f_err_throw('Incorrect once-and-for-all psp application, istart_ck, nprojel= '+&
             yaml_toa([istart_ck,nl%nprojel]),err_name='BIGDFT_RUNTIME_ERROR')
     end if
  end if

  !used on the on-the-fly projector creation
  if (nwarnings /= 0 .and. iproc == 0) then
     call yaml_map('Calculating wavelets expansion of projectors, found warnings',nwarnings,fmt='(i0)')
     call yaml_newline()
     call yaml_warning('Projectors too rough: Consider modifying hgrid and/or the localisation radii.')
  end if

  call timing(iproc,'ApplyProj     ','OF')

  call f_release_routine()

contains

  !>code factorization useful for routine restructuring
  subroutine nl_psp_application()
    implicit none
    !local variables
    integer :: ncplx_p,nvctr_p
    !real(gp), dimension(3,3,4) :: hij
    type(atomic_proj_coeff), dimension(:,:,:), pointer :: prj
    real(gp) :: eproj
    real(wp), dimension(:), pointer :: proj_ptr

    hpsi_ptr => ob_ket_map(hpsi,psi_it)

    if (paw%usepaw) then
       spsi_ptr => ob_ket_map(paw%spsi,psi_it)
       call apply_atproj_iorb_paw(atit%iat,psi_it%iorbp,istart_c,&
            at,psi_it%ob%orbs,psi_it%lr%wfd,nl,&
            psi_it%phi_wvl,hpsi_ptr,spsi_ptr,eproj_sum,&
            paw)
    else
       if (all(psi_it%kpoint == 0.0_gp)) then
          ncplx_p=1
       else
          ncplx_p=2
       end if
!!$       if (psi_it%nspinor > 1) then !which means 2 or 4
!!$          ncplx_w=2
!!$          n_w=psi_it%nspinor/2
!!$       else
!!$          ncplx_w=1
!!$          n_w=1
!!$       end if

       !extract hij parameters
!!$       call hgh_hij_matrix(at%npspcode(atit%ityp),at%psppar(0,0,atit%ityp),hij)

       call allocate_prj_ptr(atit%iat,atit%ityp,psi_it%ispin,at,nl,prj)

       nvctr_p=nl%pspd(atit%iat)%plr%wfd%nvctr_c+7*nl%pspd(atit%iat)%plr%wfd%nvctr_f

!!$       call NL_HGH_application(hij,&
!!$            ncplx_p,mproj,nl%pspd(atit%iat)%plr%wfd,nl%proj(istart_c),&
!!$            psi_it%ncplx,psi_it%n_ket,psi_it%lr%wfd,nl%pspd(atit%iat)%tolr(iilr),&
!!$            nl%wpack,nl%scpr,nl%cproj,nl%hcproj,&
!!$            psi_it%phi_wvl,hpsi_ptr,eproj)

       proj_ptr => nl%proj(istart_c:istart_c-1+nvctr_p*mproj*ncplx_p)

       call hgh_psp_application(prj,ncplx_p,mproj,nl%pspd(atit%iat)%plr%wfd,&
            proj_ptr,&
            psi_it%ncplx,psi_it%n_ket,psi_it%lr%wfd,nl%pspd(atit%iat)%tolr(iilr),&
            nl%wpack,nl%scpr,nl%cproj,nl%hcproj,&
            psi_it%phi_wvl,hpsi_ptr,eproj)

       !here the cproj can be extracted to update the density matrix for the atom iat
       if (associated(nl%iagamma)) then
          call cproj_to_gamma(atit%iat,nl%proj_G,psi_it%n_ket,mproj,lmax_ao,&
               max(psi_it%ncplx,ncplx_p),nl%cproj,psi_it%kwgt*psi_it%occup,&
               nl%iagamma(0,atit%iat),&
               nl%gamma_mmp(1,1,1,1,psi_it%ispin))
       end if

       call f_free_prj_ptr(prj)

       istart_c=istart_c+nvctr_p*ncplx_p*mproj

       eproj_sum=eproj_sum+psi_it%kwgt*psi_it%occup*eproj
    end if

  end subroutine nl_psp_application


  subroutine allocate_prj_ptr(iat,ityp,ispin,at,nl,prj)
    integer, intent(in) :: iat,ispin,ityp
    type(atoms_data), intent(in) :: at
    type(DFT_PSP_projectors), intent(inout) :: nl
    type(atomic_proj_coeff), dimension(:,:,:), pointer :: prj
    !local variables
    integer, parameter :: LMAX=3,IMAX=3
    logical :: occ_ctrl
    integer :: i,l,j,igamma,m
    real(gp), dimension(3,3,4) :: hij

    igamma=0
    occ_ctrl=.false.

    allocate(prj(3,3,4))

    !first extract the hij matrix as usual
    call hgh_hij_matrix(at%npspcode(ityp),at%psppar(0,0,ityp),hij)

    !then allocate the structure as needed
    do l=1,LMAX+1
       !the matrix is used only for i=1 at present
       if (associated(nl%iagamma)) igamma=nl%iagamma(l-1,iat)
       !if the matrix is available search for the target
       if (associated(at%gamma_targets)) &
            !if the given point need a target then associate the actual potential
            occ_ctrl= associated(at%gamma_targets(l-1,ispin,iat)%ptr) .and. nl%apply_gamma_target
       do i=1,IMAX
          call nullify_atomic_proj_coeff(prj(i,i,l))
          prj(i,i,l)%hij=hij(i,i,l)
          if (i==1 .and. occ_ctrl) then
             !it has to be discussed if the coefficient should change in to one
             !and we have to add the h11 term in the diagonal
             prj(i,i,l)%mat=&
                  f_malloc_ptr(src_ptr=at%gamma_targets(l-1,ispin,iat)%ptr,id='prjmat')
             do m=1,2*l-1
                prj(i,i,l)%mat(m,m)=prj(i,i,l)%mat(m,m)+prj(i,i,l)%hij
             end do
             prj(i,i,l)%hij=1.0_gp
          end if
          do j=i+1,IMAX
             call nullify_atomic_proj_coeff(prj(i,j,l))
             call nullify_atomic_proj_coeff(prj(j,i,l))
             !allocate upper triangular and associate lower triangular
             prj(i,j,l)%hij=hij(i,j,l)

             prj(j,i,l)%mat => prj(i,j,l)%mat
             prj(j,i,l)%hij=hij(j,i,l)
          end do
       end do
    end do

  end subroutine allocate_prj_ptr


end subroutine NonLocalHamiltonianApplication


!!$subroutine NonLocalHamiltonianApplication_new(iproc,at,npsidim_orbs,orbs,&
!!$     Lzd,nl,psi,hpsi,eproj_sum,paw)
!!$  use module_base
!!$  use module_types
!!$  use yaml_output
!!$  !  use module_interfaces
!!$  use psp_projectors_base, only: PSP_APPLY_SKIP
!!$  use psp_projectors, only: projector_has_overlap,get_proj_locreg
!!$  use public_enums, only: PSPCODE_PAW
!!$  use module_atoms
!!$  use orbitalbasis
!!$  implicit none
!!$  integer, intent(in) :: iproc, npsidim_orbs
!!$  type(atoms_data), intent(in) :: at
!!$  type(orbitals_data),  intent(in) :: orbs
!!$  type(local_zone_descriptors), intent(in) :: Lzd
!!$  type(DFT_PSP_projectors), intent(inout) :: nl
!!$  real(wp), dimension(npsidim_orbs), intent(in) :: psi
!!$  real(wp), dimension(npsidim_orbs), intent(inout) :: hpsi
!!$  type(paw_objects),intent(inout)::paw
!!$  real(gp), intent(out) :: eproj_sum
!!$  !local variables
!!$  character(len=*), parameter :: subname='NonLocalHamiltonianApplication'
!!$  logical :: overlap
!!$  integer :: istart_ck,nwarnings
!!$  integer :: iproj,istart_c,mproj,iilr
!!$  type(ket) :: psi_it
!!$  type(orbital_basis) :: psi_ob
!!$  type(atoms_iterator) :: atit
!!$  real(wp), dimension(:), pointer :: hpsi_ptr,spsi_ptr
!!$
!!$  !integer :: ierr
!!$  !real(kind=4) :: tr0, tr1, t0, t1
!!$  !real(kind=8) :: time0, time1, time2, time3, time4, time5, ttime
!!$  !real(kind=8), dimension(0:4) :: times
!!$
!!$  call f_routine(id='NonLocalHamiltonianApplication')
!!$
!!$  eproj_sum=0.0_gp
!!$
!!$  ! apply all PSP projectors for all orbitals belonging to iproc
!!$  call timing(iproc,'ApplyProj     ','ON')
!!$
!!$  !initialize the orbital basis object, for psi and hpsi
!!$  call orbital_basis_associate(psi_ob,orbs=orbs,phis_wvl=psi,Lzd=Lzd)
!!$
!!$  nwarnings=0
!!$  if(paw%usepaw) call f_zero(orbs%npsidim_orbs, paw%spsi(1))
!!$
!!$  !here the localisation region should be changed, temporary only for cubic approach
!!$
!!$  !apply the projectors following the strategy (On-the-fly calculation or not)
!!$
!!$  !apply the projectors  k-point of the processor
!!$  !iterate over the orbital_basis
!!$  psi_it=orbital_basis_iterator(psi_ob)
!!$
!!$  !iterate over the projectors
!!$  prit=projector_iterator(nl)
!!$
!!$  loop_kpt: do while(ket_next_kpt(psi_it))
!!$     loop_lr: do while(ket_next_locreg(psi_it,ikpt=psi_it%ikpt))
!!$        if (nl%on_the_fly) then
!!$           loop_proj: do while(projector_next(prit,ilr=psi_it%ilr,lr=psi_it%lr))
!!$              loop_psi_kpt: do while(ket_next(psi_it,ikpt=psi_it%ikpt,ilr=psi_it%ilr))
!!$                 hpsi_ptr => ob_ket_map(hpsi,psi_it)
!!$                 if (paw%usepaw) then
!!$                    spsi_ptr => ob_ket_map(paw%spsi,psi_it)
!!$                    call apply_atproj_iorb_paw(prit%iat,psi_it%iorbp,prit%istart_c,&
!!$                         at,psi_it%ob%orbs,psi_it%lr%wfd,nl,&
!!$                         psi_it%phi_wvl,hpsi_ptr,spsi_ptr,eproj_sum,&
!!$                         paw)
!!$                 else
!!$                    if (psi_it%nspinor > 1) then !which means 2 or 4
!!$                       ncplx_w=2
!!$                       n_w=psi_it%nspinor/2
!!$                    else
!!$                       ncplx_w=1
!!$                       n_w=1
!!$                    end if
!!$                    call NL_HGH_application(prit%hij,&
!!$                         prit%ncplx,prit%mproj,prit%lr%wfd,prit%proj,&
!!$                         ncplx_w,n_w,psi_it%lr%wfd,prit%tolr,prit%nl%wpack,prit%nl%scpr,prit%nl%cproj,prit%nl%hcproj,&
!!$                         psi_it%phi_wvl,hpsi_ptr,eproj)
!!$                    eproj_sum=eproj_sum+psi_it%kwgt*psi_it%occup*eproj
!!$                 end if
!!$              end do loop_psi_kpt
!!$           end do loop_proj
!!$        else
!!$           loop_psi_kpt2: do while(ket_next(psi_it,ikpt=psi_it%ikpt,ilr=psi_it%ilr))
!!$              loop_proj_2: do while(projector_next(prit,ilr=psi_it%ilr,lr=psi_it%lr))
!!$                 hpsi_ptr => ob_ket_map(hpsi,psi_it)
!!$                 if (paw%usepaw) then
!!$                    spsi_ptr => ob_ket_map(paw%spsi,psi_it)
!!$                    call apply_atproj_iorb_paw(prit%iat,psi_it%iorbp,prit%istart_c,&
!!$                         at,psi_it%ob%orbs,psi_it%lr%wfd,nl,&
!!$                         psi_it%phi_wvl,hpsi_ptr,spsi_ptr,eproj_sum,&
!!$                         paw)
!!$                 else
!!$                    if (psi_it%nspinor > 1) then !which means 2 or 4
!!$                       ncplx_w=2
!!$                       n_w=psi_it%nspinor/2
!!$                    else
!!$                       ncplx_w=1
!!$                       n_w=1
!!$                    end if
!!$                    call NL_HGH_application(prit%hij,&
!!$                         prit%ncplx,prit%mproj,prit%lr%wfd,prit%proj,&
!!$                         ncplx_w,n_w,psi_it%lr%wfd,prit%tolr,prit%nl%wpack,prit%nl%scpr,prit%nl%cproj,prit%nl%hcproj,&
!!$                         psi_it%phi_wvl,hpsi_ptr,eproj)
!!$                    eproj_sum=eproj_sum+psi_it%kwgt*psi_it%occup*eproj
!!$                 end if
!!$              end do loop_proj_2
!!$           end do loop_psi_kpt2
!!$        end if
!!$     end do loop_lr
!!$  end do loop_kpt
!!$
!!$  call orbital_basis_release(psi_ob)
!!$
!!$  if (.not. nl%on_the_fly .and. Lzd%nlr==1) then !TO BE REMOVED WITH NEW PROJECTOR APPLICATION
!!$     if (prit%istart_ck-1 /= nl%nprojel .and. orbs%norbp>0) then
!!$        call f_err_throw('Incorrect once-and-for-all psp application, istart_ck, nprojel= '+&
!!$             yaml_toa([prit%istart_ck,nl%nprojel]),err_name='BIGDFT_RUNTIME_ERROR')
!!$     end if
!!$  end if
!!$
!!$  !used on the on-the-fly projector creation
!!$  if (prit%nwarnings /= 0 .and. iproc == 0) then
!!$     call yaml_map('Calculating wavelets expansion of projectors, found warnings',nwarnings,fmt='(i0)')
!!$     call yaml_newline()
!!$     call yaml_warning('Projectors too rough: Consider modifying hgrid and/or the localisation radii.')
!!$  end if
!!$
!!$  call timing(iproc,'ApplyProj     ','OF')
!!$
!!$  call f_release_routine()
!!$
!!$end subroutine NonLocalHamiltonianApplication_new


!> Routine which calculates the application of nonlocal projectors on the wavefunctions
!! Reduce the wavefunction in case it is needed
subroutine NonLocalHamiltonianApplication_old(iproc,at,npsidim_orbs,orbs,&
     Lzd,nl,psi,hpsi,eproj_sum,paw)
  use module_base
  use module_types
  use yaml_output
!  use module_interfaces
  use psp_projectors, only: projector_has_overlap
  use f_utils
  use public_enums, only: PSPCODE_PAW
  implicit none
  integer, intent(in) :: iproc, npsidim_orbs
  type(atoms_data), intent(in) :: at
  type(orbitals_data),  intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  type(DFT_PSP_projectors), intent(inout) :: nl
  real(wp), dimension(npsidim_orbs), intent(in) :: psi
  real(wp), dimension(npsidim_orbs), intent(inout) :: hpsi
  type(paw_objects),intent(inout)::paw
  real(gp), intent(out) :: eproj_sum
  !local variables
  logical :: newmethod
  character(len=*), parameter :: subname='NonLocalHamiltonianApplication'
  logical :: dosome, overlap
  !logical :: goon
  integer :: ikpt,istart_ck,ispsi_k,isorb,ieorb,nspinor,iorb,iat,nwarnings
  integer :: iproj,ispsi,istart_c,ilr,ilr_skip,mproj,iatype,ispinor,iilr,jlr
  real(wp) :: hp,eproj
  real(wp), dimension(:), allocatable :: scpr
  !integer :: ierr
  !real(kind=4) :: tr0, tr1, t0, t1
  !real(kind=8) :: time0, time1, time2, time3, time4, time5, ttime
  !real(kind=8), dimension(0:4) :: times

  call f_routine(id='NonLocalHamiltonianApplication')

  newmethod=.true.

  eproj_sum=0.0_gp
  !quick return if no orbitals on this task
  if (orbs%norbp == 0) then
     if (paw%usepaw) then
        call gather_cprj(orbs, paw)
     end if
     call f_release_routine()
     return
  end if

  !FOR VESTA
  !call mpi_barrier(bigdft_mpi%mpi_comm,ierr)

  ! apply all PSP projectors for all orbitals belonging to iproc
  call timing(iproc,'ApplyProj     ','ON')

  !time1=0.0d0
  !time2=0.0d0
  !time3=0.0d0
  !time4=0.0d0
  !times=0.0d0
  !call cpu_time(t0)

  !array of the scalar products
  scpr=f_malloc(orbs%norbp*orbs%nspinor,id='scpr')


  nwarnings=0

  if(paw%usepaw) then
     call f_zero(orbs%npsidim_orbs, paw%spsi(1))
     !call f_memcpy(n = npsidim_orbs, dest = paw%spsi(1), src = psi(1))
     newmethod=.false.
  end if
  !here the localisation region should be changed, temporary only for cubic approach

  !apply the projectors following the strategy (On-the-fly calculation or not)

  !apply the projectors  k-point of the processor
  !starting k-point
  ikpt=orbs%iokpt(1)
  istart_ck=1
  ispsi_k=1
  ispsi=1 !to initialize the value in case of no projectors
  loop_kpt: do

     call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)
     !localisation regions loop
     loop_lr: do ilr=1,Lzd%nlr
        !do something only if at least one of the orbitals lives in the ilr
        dosome=.false.
        do iorb=isorb,ieorb
           !dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr)
           !dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr .and. lzd%doHamAppl(ilr))
           dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr)
           if (dosome) exit
        end do
        if (.not. dosome) cycle loop_lr

        if (nl%on_the_fly) then
           !first create a projector ,then apply it for everyone
           iproj=0
           loop_atoms_1: do iat=1,at%astruct%nat

              ! Check whether the projectors of this atom have an overlap with locreg ilr
              !!!goon=.false.
              !!!do jlr=1,nl%pspd(iat)%noverlap
              !!!    if (nl%pspd(iat)%lut_tolr(jlr)==ilr) then
              !!!        goon=.true.
              !!!        iilr=jlr
              !!!        exit
              !!!    end if
              !!!end do
              !!!if (.not.goon) cycle loop_atoms_1

              !!!iatype=at%astruct%iatype(iat)

              !!!mproj=nl%pspd(iat)%mproj
              !!!!no projector on this atom
              !!!if(mproj == 0) cycle
              !!!!projector not overlapping with the locreg
              !!!!!iilr=nl%pspd(iat)%lut_tolr(ilr)
              !!!!!if (iilr==PSP_APPLY_SKIP) cycle
              !!!if(nl%pspd(iat)%tolr(iilr)%strategy == PSP_APPLY_SKIP) cycle
              !!!!check if the atom projector intersect with the given localisation region
              !!!!this part can be moved at the place of the analysis between psp and lrs
              !!!!call cpu_time(tr0)

              !!!call check_overlap(Lzd%Llr(ilr), nl%pspd(iat)%plr, Lzd%Glr, overlap)

              ! Check whether the projectors of this atom have an overlap with locreg ilr
              overlap = projector_has_overlap(iat, ilr, lzd%llr(ilr), lzd%glr, nl)
              if(.not. overlap) cycle
              !call cpu_time(tr1)
              !time1=time1+real(tr1-tr0,kind=8)
              ! Now create the projector
              istart_c=1
              iatype=at%astruct%iatype(iat)
              do jlr=1,nl%pspd(iat)%noverlap
                  if (nl%pspd(iat)%lut_tolr(jlr)==ilr) then
                      iilr=jlr
                      exit
                  end if
              end do
              mproj=nl%pspd(iat)%mproj
              call atom_projector(nl, iatype, iat, at%astruct%atomnames(iatype), &
                   & at%astruct%geocode, 0, Lzd%Glr, Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3), &
                   & orbs%kpts(1,ikpt), orbs%kpts(2,ikpt), orbs%kpts(3,ikpt), &
                   & istart_c, iproj, nwarnings)
              !call cpu_time(tr0)
              !time2=time2+real(tr0-tr1,kind=8)
              !apply the projector to all the orbitals belonging to the processor
              !this part can be factorized somewhere else
              if (mproj ==1 .and. all(orbs%kpts(:,ikpt) == 0.0_gp) .and. &
                   newmethod .and. Lzd%nlr == 1) then
                 hp=at%psppar(1,1,iatype) !it is supposed that the only projector is the i=1 case
                 ispsi=ispsi_k
                 istart_c=1

                 call apply_oneproj_operator(nl%pspd(iat)%plr%wfd,nl%proj(istart_c),hp,&
                      (ieorb-isorb+1)*nspinor,Lzd%Llr(ilr)%wfd,psi(ispsi),hpsi(ispsi),scpr)

              !call cpu_time(tr1)
              !time3=time3+real(tr1-tr0,kind=8)

                 istart_c=istart_c+nl%pspd(iat)%plr%wfd%nvctr_c+7*nl%pspd(iat)%plr%wfd%nvctr_f
!                 call f_malloc_dump_status()
                 ispsi=ispsi+&
                      (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*nspinor*(ieorb-isorb+1)
                 do iorb=isorb,ieorb
                    do ispinor=1,nspinor
                       eproj=hp*(scpr((iorb-isorb)*nspinor+ispinor)**2)
                       eproj_sum=eproj_sum+orbs%kwgts(orbs%iokpt(iorb))*&
                            orbs%occup(iorb+orbs%isorb)*eproj
                    end do
                 end do
              else
                 ispsi=ispsi_k
                 do iorb=isorb,ieorb
                    if (orbs%inwhichlocreg(iorb+orbs%isorb) /= ilr) then
                       !increase ispsi to meet orbital index
                       ilr_skip=orbs%inwhichlocreg(iorb+orbs%isorb)
                       ispsi=ispsi+(Lzd%Llr(ilr_skip)%wfd%nvctr_c+7*Lzd%Llr(ilr_skip)%wfd%nvctr_f)*nspinor
                       cycle
                    end if
                    istart_c=1

                    call nl_psp_application_old()

                    !                print *,'iorb,iat,eproj',iorb+orbs%isorb,ispsi,iat,eproj_sum
                    ispsi=ispsi+&
                         (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*nspinor
                 end do
              end if
           end do loop_atoms_1
           !call cpu_time(tr0)
           !time4=time4+real(tr0-tr1,kind=8)
           !for the moment, localization region method is not tested with
           !once-and-for-all calculation
        else if (Lzd%nlr == 1) then

           !loop over the interesting my orbitals, and apply all the projectors over all orbitals
           ispsi=ispsi_k
           do iorb=isorb,ieorb
              if (orbs%inwhichlocreg(iorb+orbs%isorb) /= ilr) then
                 !increase ispsi to meet orbital index
                 ilr_skip=orbs%inwhichlocreg(iorb+orbs%isorb)
                 ispsi=ispsi+(Lzd%Llr(ilr_skip)%wfd%nvctr_c+7*Lzd%Llr(ilr_skip)%wfd%nvctr_f)*nspinor
                 cycle
              end if

              istart_c=istart_ck !TO BE CHANGED IN ONCE-AND-FOR-ALL
              loop_atoms_2: do iat=1,at%astruct%nat


                 !! ! Check whether the projectors of this atom have an overlap with locreg ilr
                 !! goon=.false.
                 !! do jlr=1,nl%pspd(iat)%noverlap
                 !!     if (nl%pspd(iat)%lut_tolr(jlr)==ilr) then
                 !!         goon=.true.
                 !!         iilr=jlr
                 !!         exit
                 !!     end if
                 !! end do
                 !! if (.not.goon) cycle loop_atoms_2

                 !!iatype=at%astruct%iatype(iat)
                 !!! Check if atom has projectors, if not cycle
                 !!mproj=nl%pspd(iat)%mproj
                 !!if(mproj == 0) cycle
                 !!!projector not overlapping with the locreg
                 !!!!iilr=nl%pspd(iat)%lut_tolr(ilr)
                 !!!!if (iilr==PSP_APPLY_SKIP) cycle
                 !!if(nl%pspd(iat)%tolr(iilr)%strategy == PSP_APPLY_SKIP) cycle

                 !!!check if the atom intersect with the given localisation region
                 !!call check_overlap(Lzd%Llr(ilr), nl%pspd(iat)%plr, Lzd%Glr, overlap)
                 overlap = projector_has_overlap(iat, ilr, lzd%llr(ilr), lzd%glr, nl)
                 if(.not. overlap) cycle loop_atoms_2 !stop 'ERROR all atoms should be in global'

                 iatype=at%astruct%iatype(iat)
                 do jlr=1,nl%pspd(iat)%noverlap
                     if (nl%pspd(iat)%lut_tolr(jlr)==ilr) then
                         iilr=jlr
                         exit
                     end if
                 end do
                 mproj=nl%pspd(iat)%mproj
                 call nl_psp_application_old()

                 !print *,'iorb,iat,eproj',iorb+orbs%isorb,iat,eproj_sum
              end do loop_atoms_2
              !print *,'TOTALPSI',iorb+orbs%isorb,sum(psi(ispsi:&
              !    ispsi+(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor-1)),&
              !     dot((Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor,&
              !     hpsi(ispsi),1,hpsi(ispsi),1)
              ispsi=ispsi+(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*nspinor
           end do
           istart_ck=istart_c !TO BE CHANGED IN THIS ONCE-AND-FOR-ALL
        else
           ! COULD CHANGE THIS NOW !!
           call yaml_warning('Localization Regions not allowed in once-and-for-all')
           stop
        end if

     end do loop_lr

     !last k-point has been treated
     if (ieorb == orbs%norbp) exit loop_kpt

     ikpt=ikpt+1
     ispsi_k=ispsi

  end do loop_kpt

  if(paw%usepaw) then
     call gather_cprj(orbs, paw)
  end if

  if (.not. nl%on_the_fly .and. Lzd%nlr==1) then !TO BE REMOVED WITH NEW PROJECTOR APPLICATION
     if (istart_ck-1 /= nl%nprojel) then
        call yaml_warning('Incorrect once-and-for-all psp application')
        stop
     end if
  end if
  !for the moment it has to be removed. A number of components in orbital distribution should be defined
  !if (ispsi-1 /= (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp) stop 'incorrect V_nl psi application'

  !used on the on-the-fly projector creation
  if (nwarnings /= 0 .and. iproc == 0) then
     call yaml_map('Calculating wavelets expansion of projectors, found warnings',nwarnings,fmt='(i0)')
      call yaml_newline()
      call yaml_warning('Projectors too rough: Consider modifying hgrid and/or the localisation radii.')
  end if

  call f_free(scpr)
  !call cpu_time(t1)
  !time0=real(t1-t0,kind=8)

  !print*,'iproc,times,sum,ttime',iproc,time1,time2,time3,time4,time1+time2+time3+time4,time0
  !call mpi_barrier(bigdft_mpi%mpi_comm,ierr)

  call timing(iproc,'ApplyProj     ','OF')

  call f_release_routine()

contains

!!$  !> count the number of projectors given a set of psppar
!!$  pure function nproj(psppar)
!!$    use module_base, only: gp
!!$    implicit none
!!$    real(gp), dimension(0:4,0:6), intent(in) :: psppar
!!$    integer :: nproj
!!$    !local variables
!!$    integer :: i,l
!!$    nproj = 0
!!$    !count over all the channels
!!$    do l=1,4
!!$       !loop over all the projectors of the channel
!!$       do i=1,3
!!$          !loop over all the components of the projector
!!$          if (psppar(l,i) /= 0.0_gp) nproj=nproj+2*l-1
!!$       end do
!!$    end do
!!$  end function nproj

  !>code factorization useful for routine restructuring
  subroutine nl_psp_application_old()
    implicit none
    !local variables
    integer :: ncplx_p,ncplx_w,n_w,nvctr_p
    real(gp), dimension(3,3,4) :: hij
    real(gp) :: eproj

    if (newmethod) then

       if (all(orbs%kpts(:,ikpt) == 0.0_gp)) then
          ncplx_p=1
       else
          ncplx_p=2
       end if
       if (nspinor > 1) then !which means 2 or 4
          ncplx_w=2
          n_w=nspinor/2
       else
          ncplx_w=1
          n_w=1
       end if

       !extract hij parameters
       call hgh_hij_matrix(at%npspcode(iatype),at%psppar(0,0,iatype),hij)

       call NL_HGH_application(hij,&
            ncplx_p,mproj,nl%pspd(iat)%plr%wfd,nl%proj(istart_c),&
            ncplx_w,n_w,Lzd%Llr(ilr)%wfd,nl%pspd(iat)%tolr(iilr),nl%wpack,nl%scpr,nl%cproj,nl%hcproj,&
            psi(ispsi),hpsi(ispsi),eproj)

       nvctr_p=nl%pspd(iat)%plr%wfd%nvctr_c+7*nl%pspd(iat)%plr%wfd%nvctr_f
       istart_c=istart_c+nvctr_p*ncplx_p*mproj

       eproj_sum=eproj_sum+&
            orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*eproj
    else
       if(paw%usepaw) then
          call apply_atproj_iorb_paw(iat,iorb,istart_c,&
               at,orbs,Lzd%Llr(ilr)%wfd,nl,&
               psi(ispsi),hpsi(ispsi),paw%spsi(ispsi),eproj_sum,&
               paw)
       else
          call apply_atproj_iorb_new(iat,iorb,istart_c,&
               nl%nprojel,&
               at,orbs,Lzd%Llr(ilr)%wfd,nl%pspd(iat)%plr,&
               nl%proj,psi(ispsi),hpsi(ispsi),eproj_sum)
       end if
    end if
  end subroutine nl_psp_application_old

END SUBROUTINE NonLocalHamiltonianApplication_old

subroutine inspect_hpsi(orbs,Lzd,hpsi)
  use module_types
  use module_defs, only: wp
  implicit none
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  real(wp), dimension(*), intent(in) :: hpsi
  !local variables
  integer :: ispsi,iorb,nvctr,ilr
  ispsi=1
  do iorb=1,orbs%norbp
     ilr=orbs%inWhichLocreg(orbs%isorb+iorb)
     nvctr=(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
     print *,'sum for orb',iorb,sum(hpsi(ispsi:ispsi+nvctr-1))
     ispsi=ispsi+nvctr
  end do
end subroutine inspect_hpsi


!> routine which puts a barrier to ensure that both local and nonlocal hamiltonians have been applied
!! in the GPU case puts a barrier to end the overlapped Local and nonlocal applications
subroutine SynchronizeHamiltonianApplication(nproc,npsidim_orbs,orbs,Lzd,GPU,xc,&
     hpsi,energs,energs_work)
   use module_base
   use module_types
   use module_xc
   implicit none
   integer, intent(in) :: nproc,npsidim_orbs
   type(orbitals_data), intent(in) :: orbs
   type(local_zone_descriptors), intent(in) :: Lzd
   type(GPU_pointers), intent(inout) :: GPU
   type(xc_info), intent(in) :: xc
   type(energy_terms), intent(inout) :: energs
   !real(gp), intent(inout) :: ekin_sum,epot_sum,eproj_sum,evsic,eexctX
   real(wp), dimension(npsidim_orbs), intent(inout) :: hpsi
   type(work_mpiaccumulate),optional,intent(inout) :: energs_work
   !local variables
   character(len=*), parameter :: subname='SynchronizeHamiltonianApplication'
   logical :: exctX
   integer :: iorb,ispsi,ilr
   real(gp), dimension(4) :: wrkallred

   call f_routine(id='SynchronizeHamiltonianApplication')

   if(GPU%OCLconv) then! needed also in the non_ASYNC since now NlPSP is before .and. ASYNCconv)) then
      call finish_hamiltonian_OCL(orbs,energs%ekin,energs%epot,GPU)
      ispsi=1
      do iorb=1,orbs%norbp
         ilr=orbs%inWhichLocreg(orbs%isorb+iorb)
         call axpy((Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor,&
              1.0_wp,GPU%hpsi_ASYNC(ispsi),1,hpsi(ispsi),1)
         ispsi=ispsi+&
             (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
      end do
      call f_free_ptr(GPU%hpsi_ASYNC)
   endif



   exctX = xc_exctXfac(xc) /= 0.0_gp

   !up to this point, the value of the potential energy is
   !only taking into account the local potential part
   !whereas it should consider also the value coming from the
   !exact exchange operator (twice the exact exchange energy)
   !this operation should be done only here since the exctX energy is already reduced
   !SM: Divide by nproc due to the reduction later on
   if (exctX) energs%epot=energs%epot+2.0_gp*energs%eexctX/real(nproc,kind=8)
   !energies reduction
   if (nproc > 1) then
      if (present(energs_work)) then
         energs_work%sendbuf(1) = energs%ekin
         energs_work%sendbuf(2) = energs%epot
         energs_work%sendbuf(3) = energs%eproj
         energs_work%sendbuf(4) = energs%evsic
         energs_work%receivebuf(:) = 0.d0
         !LG: why the window is opened with only one element whereas we communicate 4?
         !! I correct this point as it seems a bug to me
         !energs_work%window = mpiwindow(1, energs_work%receivebuf(1), bigdft_mpi%mpi_comm)
         call fmpi_win_create(energs_work%window,energs_work%receivebuf(1),4,comm=bigdft_mpi%mpi_comm)
         call fmpi_win_fence(energs_work%window,FMPI_WIN_OPEN)
         call fmpi_accumulate(energs_work%sendbuf(1),target_rank=0,count=4,op=FMPI_SUM,&
              target_disp=int(0,fmpi_address),win=energs_work%window)
         !call mpiaccumulate(energs_work%sendbuf(1), 4, 0, &
         !     int(0,kind=mpi_address_kind), FMPI_SUM, energs_work%window)
      else
         wrkallred(1)=energs%ekin
         wrkallred(2)=energs%epot
         wrkallred(3)=energs%eproj
         wrkallred(4)=energs%evsic

         call fmpi_allreduce(wrkallred,FMPI_SUM,comm=bigdft_mpi%mpi_comm)

         energs%ekin=wrkallred(1)
         energs%epot=wrkallred(2)
         energs%eproj=wrkallred(3)
         energs%evsic=wrkallred(4)
      end if
   else if (present(energs_work)) then
       ! Do a "fake communication"
       energs_work%receivebuf(1) = energs%ekin
       energs_work%receivebuf(2) = energs%epot
       energs_work%receivebuf(3) = energs%eproj
       energs_work%receivebuf(4) = energs%evsic
   end if

   !define the correct the hartree energy in the case of electrostatic contribution
   if (energs%eelec /= 0.0_gp) then
      energs%eh=energs%epot-energs%eelec-energs%evxc
   end if
   call f_release_routine()

END SUBROUTINE SynchronizeHamiltonianApplication


!> Nullify potential (pot)
subroutine free_full_potential(nproc,flag,xc,pot)
   use module_base
   use module_xc
   implicit none
   integer, intent(in) :: nproc, flag
   type(xc_info), intent(in) :: xc
   real(wp), dimension(:), pointer :: pot
   !local variables
   logical :: odp
   call f_routine(id='free_full_potential')
   odp = xc_exctXfac(xc) /= 0.0_gp
   if (nproc > 1 .or. odp .or. flag > 0 ) then
      !call f_free_ptr(pot)
      call f_free_ptr(pot)
      nullify(pot)
   else
      nullify(pot)
   end if

   call f_release_routine()
END SUBROUTINE free_full_potential


!> Calculate total energies from the energy terms
subroutine total_energies(energs, iter, iproc)
  use module_base
  use module_types
  implicit none
  type(energy_terms), intent(inout) :: energs
  integer, intent(in) :: iter, iproc
  call f_routine(id='total_energies')
  !band structure energy calculated with occupation numbers
  energs%ebs=energs%ekin+energs%epot+energs%eproj !the potential energy contains also exctX
  !this is the Kohn-Sham energy
  energs%eKS=energs%ebs-energs%eh+energs%exc-energs%evxc-&
       energs%eexctX-energs%evsic+energs%eion+energs%edisp+&
       energs%epawdc!-energs%excrhoc

  ! Gibbs Free Energy
  energs%energy=energs%eKS-energs%eTS+energs%ePV

  if (energs%c_obj /= 0) then
     call timing(iproc,'energs_signals','ON')
     call energs_emit(energs%c_obj, iter, 0) ! 0 is for BIGDFT_E_KS in C.
     call timing(iproc,'energs_signals','OF')
  end if
  call f_release_routine()
end subroutine total_energies


!> Extract the energy (the quantity which has to be minimised by the wavefunction)
!! and calculate the corresponding gradient.
!! The energy can be the actual Kohn-Sham energy or the trace of the hamiltonian,
!! depending of the functional we want to calculate. The gradient wrt the wavefunction
!! is put in hpsi accordingly to the functional
subroutine calculate_energy_and_gradient(iter,iproc,nproc,GPU,ncong,scf_mode,&
     energs,wfn,gnrm,gnrm_zero)
  use module_base
  use module_types
  use module_interfaces, only: orthoconstraint, preconditionall2
  use io, only: write_energies
  use yaml_output
  use communications, only: transpose_v, untranspose_v
  use communications, only: toglobal_and_transpose
  implicit none
  integer, intent(in) :: iproc,nproc,ncong,iter
  type(f_enumerator), intent(in) :: scf_mode
  type(energy_terms), intent(inout) :: energs
  type(GPU_pointers), intent(in) :: GPU
  type(DFT_wavefunction), intent(inout) :: wfn
  real(gp), intent(out) :: gnrm,gnrm_zero
  !local variables
  character(len=*), parameter :: subname='calculate_energy_and_gradient'
  logical :: lcs,tr_min
  integer :: ikpt,iorb,k
  real(gp) :: rzeroorbs,tt,garray(2)
  real(wp), dimension(:,:,:), pointer :: mom_vec

  call f_routine(id=subname)
  !calculate orbital polarisation directions
  if(wfn%orbs%nspinor==4) then
     mom_vec = f_malloc_ptr((/ 4, wfn%orbs%norb, min(nproc, 2) /),id='mom_vec')

     call calc_moments(iproc,nproc,wfn%orbs%norb,wfn%orbs%norb_par,&
          wfn%Lzd%Glr%wfd%nvctr_c+7*wfn%Lzd%Glr%wfd%nvctr_f,wfn%orbs%nspinor,wfn%psi,mom_vec)
  else
     nullify(mom_vec)
  end if


  !transpose the hpsi wavefunction
   call toglobal_and_transpose(iproc,nproc,wfn%orbs,wfn%Lzd,wfn%comms,wfn%hpsi,wfn%psi)
  !transpose the spsi wavefunction
  if(wfn%paw%usepaw) then
     call toglobal_and_transpose(iproc,nproc,wfn%orbs,wfn%Lzd,wfn%comms,wfn%paw%spsi,wfn%psi)
  end if

  if (nproc == 1) then
     !associate psit pointer for orthoconstraint and transpose it (for the non-collinear case)
     wfn%psit => wfn%psi
     ! work array not used for nproc==1, so pass the same address
     call transpose_v(iproc,nproc,wfn%orbs,wfn%lzd%glr%wfd,wfn%comms,wfn%psit,wfn%psit)
  end if

!!$  if (iproc==0 .and. get_verbose_level() > 0) then
!!$     call yaml_map('Symmetrize Lagr. Multiplier',wfn%SIC%alpha/=0.0_gp)
!!$  end if

  ! Apply  orthogonality constraints to all orbitals belonging to iproc
  !takes also into account parallel k-points distribution
  !here the orthogonality with respect to other occupied functions should be
  !passed as an optional argument
  energs%trH_prev=energs%trH
  !disable correction for exact exhange computations, as the density matrix is the object to minimize.
  tr_min=(scf_mode .hasattr. 'MIXING') .or. energs%eexctX /=0.0_gp
  if(wfn%paw%usepaw) then
    !PAW: spsi is used.
    call orthoconstraint(iproc,nproc,wfn%orbs,wfn%comms,wfn%SIC%alpha/=0.0_gp,tr_min,&
         wfn%psit,wfn%hpsi,energs%trH,wfn%paw%spsi)
  else
    !NC:
    call orthoconstraint(iproc,nproc,wfn%orbs,wfn%comms,wfn%SIC%alpha/=0.0_gp,tr_min,&
         wfn%psit,wfn%hpsi,energs%trH)
  end if

  !retranspose the hpsi wavefunction
  call untranspose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,wfn%hpsi,wfn%psi)
  if(wfn%paw%usepaw) then
   !retranspose the spsi wavefunction
   call untranspose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,wfn%paw%spsi,wfn%psi)
  end if


  !after having calcutated the trace of the hamiltonian, the functional have to be defined
  !new value without the trace, to be added in hpsitopsi
  if (scf_mode .hasattr. 'MIXING') then
     wfn%diis%energy=energs%trH
  else
     wfn%diis%energy=energs%eKS!trH-eh+exc-evxc-eexctX+eion+edisp(not correct for non-integer occnums)
  end if

  if (iproc==0 .and. get_verbose_level() > 0) then
     call yaml_map('Orthoconstraint',.true.)
  end if


  !check that the trace of the hamiltonian is compatible with the
  !band structure energy
  !this can be done only if the occupation numbers are all equal
  tt=(energs%ebs-energs%trH)/energs%trH
!print *,'tt,energybs,trH',tt,energybs,trH
  if (abs(tt) > 1.d-10 .and. iproc==0) then
     !write this warning only if the system is closed shell
     call check_closed_shell(wfn%orbs,lcs)
     if (lcs) then
        call yaml_newline()
        call yaml_mapping_open('Energy inconsistencies')
           call yaml_map('Band Structure Energy',energs%ebs,fmt='(1pe22.14)')
           call yaml_map('Trace of the Hamiltonian',energs%trH,fmt='(1pe22.14)')
           call yaml_map('Relative inconsistency',tt,fmt='(1pe9.2)')
        call yaml_mapping_close()
     end if
  endif

  call timing(iproc,'Precondition  ','ON')

  !Preconditions all orbitals belonging to iproc
  !and calculate the partial norm of the residue
  !switch between CPU and GPU treatment
  if (GPU%OCLconv) then
     call preconditionall_OCL(wfn%orbs,wfn%Lzd%Glr,&
          wfn%Lzd%hgrids(1),wfn%Lzd%hgrids(2),wfn%Lzd%hgrids(3),ncong,&
          wfn%hpsi,gnrm,gnrm_zero,GPU)
  else
     !this is the final routine, the confining potential has to be passed to
     !switch between the global and delocalized preconditioner
     call preconditionall2(iproc,nproc,wfn%orbs,wfn%Lzd,&
          wfn%Lzd%hgrids(1),wfn%Lzd%hgrids(2),wfn%Lzd%hgrids(3),&
          ncong,wfn%orbs%npsidim_orbs,wfn%hpsi,wfn%confdatarr,gnrm,gnrm_zero)
  end if

  call timing(iproc,'Precondition  ','OF')

  !sum over all the partial residues
  if (nproc > 1) then
      garray(1)=gnrm
      garray(2)=gnrm_zero
     call fmpi_allreduce(garray,FMPI_SUM,comm=bigdft_mpi%mpi_comm)
      gnrm     =garray(1)
      gnrm_zero=garray(2)
  endif

  !count the number of orbitals which have zero occupation number
  !weight this with the corresponding k point weight
  rzeroorbs=0.0_gp
  do ikpt=1,wfn%orbs%nkpts
     do iorb=1,wfn%orbs%norb
        if (wfn%orbs%occup(iorb+(ikpt-1)*wfn%orbs%norb) == 0.0_gp) then
           rzeroorbs=rzeroorbs+wfn%orbs%kwgts(ikpt)
        end if
     end do
  end do
  !commented out, the kwgts sum already to one
  !if (orbs%nkpts > 1) nzeroorbs=nint(real(nzeroorbs,gp)/real(orbs%nkpts,gp))

  gnrm=sqrt(gnrm/(real(wfn%orbs%norb,gp)-rzeroorbs))

  if (rzeroorbs /= 0.0_gp) then
     gnrm_zero=sqrt(gnrm_zero/rzeroorbs)
  else
     gnrm_zero=0.0_gp
  end if

  !if (iproc==0 .and. get_verbose_level() > 1) write(*,'(1x,a)') 'done.'

  if (iproc==0  .and. get_verbose_level() > 0) then
     call yaml_map('Preconditioning',.true.)
     call yaml_newline()
  end if

  if (wfn%orbs%nspinor == 4) then
     !only the root process has the correct array
     ! AB supress this printout
!!!      if(iproc==0 .and. get_verbose_level() > 0) then
!!!         call yaml_sequence_open('Magnetic polarization per orbital')
!!!         call yaml_newline()
!!!         !write(*,'(1x,a)')&
!!!         !     &   'Magnetic polarization per orbital'
!!! !!$        write(*,'(1x,a)')&
!!! !!$             &   '  iorb    m_x       m_y       m_z'
!!!         do iorb=1,wfn%orbs%norb
!!!            call yaml_sequence(advance='no')
!!!            call yaml_mapping_open(flow=.true.)
!!!            call yaml_map('iorb',iorb,fmt='(i5)')
!!!            call yaml_map('M',(/(mom_vec(k,iorb,1)/mom_vec(1,iorb,1),k=2,4)/),fmt='(3f10.5)')
!!!            call yaml_mapping_close()
!!!            if (iorb < wfn%orbs%norb)call yaml_newline()
!!!            !write(*,'(1x,i5,3f10.5)') &
!!!            !     &   iorb,(mom_vec(k,iorb,1)/mom_vec(1,iorb,1),k=2,4)
!!!         end do
!!!         call yaml_sequence_close()
!!!      end if
     call f_free_ptr(mom_vec)
  end if

  !write the energy information
  if (iproc == 0) call write_energies(iter,energs,gnrm,gnrm_zero,' ',scf_mode)
  call f_release_routine()
END SUBROUTINE calculate_energy_and_gradient


!> Operations after h|psi>
!! (transposition, orthonormalisation, inverse transposition)
subroutine hpsitopsi(iproc,nproc,iter,idsx,wfn,&
   at,nlpsp,eproj_sum)
   use module_base
   use module_types
   use module_interfaces, only: orthogonalize
   use yaml_output
   use communications, only: transpose_v, untranspose_v
   use public_enums, only: PSPCODE_PAW
   implicit none
   !Arguments
   integer, intent(in) :: iproc,nproc,idsx,iter
   type(DFT_wavefunction), intent(inout) :: wfn
   type(atoms_data), intent(in) :: at
   type(DFT_PSP_projectors), intent(inout) :: nlpsp
   real(gp),optional, intent(out) :: eproj_sum
   !local variables
   !character(len=*), parameter :: subname='hpsitopsi'
   !debug
!!$   integer :: jorb,iat
   !end debug
   call f_routine(id='hpsitopsi')

   if(wfn%paw%usepaw) then
     if ( .not. present(eproj_sum)) then
        call f_err_throw('hpsitopsi for PAW needs the following optional variables: eproj',err_name='BIGDFT_RUNTIME_ERROR')
        !write(*,*)'ERROR: hpsitopsi for PAW needs the following optional variables::'
        !write(*,*)'       eproj'
        !stop
     end if
   end if

   !adjust the save variables for DIIS/SD switch
   if (iter == 1) then
      wfn%diis%ids=0
      wfn%diis%mids=1
      wfn%diis%idiistol=0
   end if
   !update variables at each iteration step
   if (idsx > 0) then
      wfn%diis%mids=mod(wfn%diis%ids,idsx)+1
      wfn%diis%ids=wfn%diis%ids+1
   end if

   wfn%diis%energy_min=min(wfn%diis%energy_min,wfn%diis%energy)

   !transpose the hpsi wavefunction
   call transpose_v(iproc,nproc,wfn%orbs,wfn%lzd%glr%wfd,wfn%comms,&
        wfn%hpsi,wfn%psi)

   !!experimental, orthogonalize the preconditioned gradient wrt wavefunction
   !call orthon_virt_occup(iproc,nproc,orbs,orbs,comms,comms,psit,hpsi,(get_verbose_level() > 2))

   !apply the minimization method (DIIS or steepest descent)
   if (iter > 0) then
      call timing(iproc,'Diis          ','ON')
      call psimix(iproc,nproc,sum(wfn%comms%ncntt(0:nproc-1)),wfn%orbs,wfn%comms,wfn%diis,wfn%hpsi,wfn%psit)
      call timing(iproc,'Diis          ','OF')
   end if

   !Update spsi, since psi has changed
   !Pending: make this with the transposed wavefunctions:
   if(wfn%paw%usepaw) then
     !retranspose psit
     call untranspose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,&
        &   wfn%psit,wfn%hpsi,out_add=wfn%psi)

     !Calculate  hpsi,spsi and cprj with new psi
     if (wfn%orbs%npsidim_orbs >0) then
       call f_zero(wfn%orbs%npsidim_orbs,wfn%hpsi(1))
       !call f_zero(wfn%orbs%npsidim_orbs,wfn%paw%spsi(1))
     end if
     call NonLocalHamiltonianApplication(iproc,at,wfn%orbs%npsidim_orbs,wfn%orbs,&
          wfn%Lzd,nlpsp,wfn%psi,wfn%hpsi,eproj_sum,wfn%paw)

!    Transpose spsi:
     call transpose_v(iproc,nproc,wfn%orbs,wfn%lzd%glr%wfd,wfn%comms,wfn%paw%spsi,wfn%hpsi)
     if (nproc == 1) &
          & call transpose_v(iproc,nproc,wfn%orbs,wfn%lzd%glr%wfd,wfn%comms,wfn%psi,wfn%hpsi)
   end if

   if (iproc == 0 .and. get_verbose_level() > 1) then
      !write(*,'(1x,a)',advance='no')&
      !&   'Orthogonalization...'
      call yaml_map('Orthogonalization Method',wfn%orthpar%methortho,fmt='(i3)')
   end if

   if(wfn%paw%usepaw) then
     call orthogonalize(iproc,nproc,wfn%orbs,wfn%comms,wfn%psit,wfn%orthpar,wfn%paw)
   else
     call orthogonalize(iproc,nproc,wfn%orbs,wfn%comms,wfn%psit,wfn%orthpar)
   end if

   !call checkortho_p(iproc,nproc,norb,nvctrp,psit)
!!$   if(wfn%paw%usepaw) then
!!$      !debug:
!!$      call checkortho_paw(iproc,wfn%orbs%norb*wfn%orbs%nspinor,&
!!$           wfn%comms%nvctr_par(iproc,0),wfn%psit,wfn%paw%spsi)
!!$      if (iproc == 0 .and. get_verbose_level() > 1) then
!!$         call yaml_sequence_open('cprj(:,1) (5 first orbitals)')
!!$         do iat=1,wfn%paw%natom
!!$            call yaml_mapping_open("atom" // trim(yaml_toa(iat, fmt = "(I0)")))
!!$            do jorb=1,min(5, wfn%orbs%norbu)
!!$               call yaml_map("iorb" // trim(yaml_toa(jorb, fmt = "(I0)")), yaml_toa(wfn%paw%cprj(iat,jorb)%cp(:,1)))
!!$               !write(*,'(a,2(i4,1x),1000f20.12)')' cprj(iat,jorb)%cp(:,:)=',iat,jorb,wfn%paw%cprj(iat,jorb)%cp(:,:)
!!$               !call yaml_comment("iorb #" // yaml_toa(jorb, fmt = "(I0)"))
!!$            end do
!!$            call yaml_mapping_close()
!!$         end do
!!$         call yaml_sequence_close()
!!$      end if
!!$   end if

   call untranspose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,&
        wfn%psit,wfn%hpsi,out_add=wfn%psi)

   if (nproc == 1) then
      nullify(wfn%psit)
   end if

   !if (iproc == 0 .and. get_verbose_level() > 1) write(*,'(1x,a)') 'done.'

   ! Emit that new wavefunctions are ready.
   if (wfn%c_obj /= 0) then
      call kswfn_emit_psi(wfn, iter, 0, iproc, nproc)
   end if

   if (iter > 0) then
      call diis_or_sd(iproc,idsx,wfn%diis)
   end if

   !previous value already filled
   wfn%diis%energy_old=wfn%diis%energy

   !
   !DEBUG hpsi
   !if(paw%usepaw==1 .and. nproc==1) then
   !   call debug_hpsi(wfn,at)

   !   !
   !   !Recalculate  hpsi,spsi and cprj with new psi
   !   !PENDING: Here we only should recalculate cprj, hpsi and spsi are already updated.
   !   !         make a routine calculate_cprj based on NonLocalHam..?
   !   write(*,*)'Recalculate hpsi,spsi and cprj with new psi'
   !   if (wfn%orbs%npsidim_orbs >0) call f_zero(wfn%orbs%npsidim_orbs,wfn%hpsi(1))
   !   if (wfn%orbs%npsidim_orbs >0.and. paw%usepaw==1) call f_zero(wfn%orbs%npsidim_orbs,paw%spsi(1))
   !   call NonLocalHamiltonianApplication(iproc,at,wfn%orbs,rxyz,&
   !          proj,wfn%Lzd,nlpspd,wfn%psi,wfn%hpsi,eproj_sum,proj_G,paw)

   !   !check ortho again, to see if psi and spsi are correct
   !   if(paw%usepaw==1) call checkortho_paw(iproc,wfn%orbs%norb*wfn%orbs%nspinor,wfn%comms%nvctr_par(iproc,0),wfn%psi,paw%spsi)

   !end if
   !END DEBUG
  call f_release_routine()
END SUBROUTINE hpsitopsi



!>   First orthonormalisation
subroutine first_orthon(iproc,nproc,orbs,lzd,comms,psi,hpsi,psit,orthpar,paw)
   use module_base
   use module_types
   use module_interfaces, only: orthogonalize
   use communications_base, only: comms_cubic
   use communications, only: transpose_v, untranspose_v
   implicit none
   integer, intent(in) :: iproc,nproc
   type(orbitals_data), intent(in) :: orbs
   type(local_zone_descriptors),intent(in) :: lzd
   type(comms_cubic), intent(in) :: comms
   type(orthon_data):: orthpar
   real(wp), dimension(:) , pointer :: psi,hpsi,psit
   type(paw_objects),optional,intent(inout)::paw
   !local variables
   character(len=*), parameter :: subname='first_orthon'
   logical :: usepaw

   call f_routine(id=subname)
   usepaw = .false.
   if(present(paw)) usepaw = paw%usepaw
   !!!  if(nspin==4) then
   !!!     nspinor=4
   !!!  else
   !!!     nspinor=1
   !!!  end if

   if (nproc > 1) then
      !allocate hpsi array (used also as transposed)
      !allocated in the transposed way such as
      !it can also be used as the transposed hpsi
      hpsi =f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='hpsi')
      !allocate transposed principal wavefunction
      psit = f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='psit')
   else
      psit => psi
   end if

   !to be substituted, must pass the wavefunction descriptors to the routine
   if (nproc>1) then
       call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi,&
          &   hpsi,recvbuf=psit)
       if (usepaw) call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,&
            & paw%spsi, hpsi)
   else
       ! work array not nedded for nproc==1, so pass the same address
       call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi,&
          &   psi,recvbuf=psit)
       if (usepaw) call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,&
            & paw%spsi, paw%spsi)
   end if

   if(usepaw) then
     call orthogonalize(iproc,nproc,orbs,comms,psit,orthpar,paw)
   else
     call orthogonalize(iproc,nproc,orbs,comms,psit,orthpar)
   end if

   !call checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)

   if (nproc>1) then
       call untranspose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psit,&
          &   hpsi,out_add=psi)
   else
       ! work array not nedded for nproc==1, so pass the same address
       call untranspose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psit,&
          &   psit,out_add=psi)
   end if

   if (nproc == 1) then
      nullify(psit)
      !allocate hpsi array
      hpsi = f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='hpsi')
   end if

   call f_release_routine()
END SUBROUTINE first_orthon


!>   Transform to KS orbitals and deallocate hpsi wavefunction (and also psit in parallel)
subroutine last_orthon(iproc,nproc,iter,wfn,evsum,opt_keeppsit)
   use module_base
   use module_types
   use communications, only: transpose_v, untranspose_v
   implicit none
   integer, intent(in) :: iproc,nproc,iter
   real(wp), intent(out) :: evsum
   type(DFT_wavefunction), intent(inout) :: wfn
   logical, optional :: opt_keeppsit
   !local variables
   logical :: keeppsit
   character(len=*), parameter :: subname='last_orthon'

   call f_routine(id=subname)
   if (present(opt_keeppsit)) then
      keeppsit=opt_keeppsit
   else
      keeppsit=.false.
   end if

   call transpose_v(iproc,nproc,wfn%orbs,wfn%lzd%glr%wfd,wfn%comms,wfn%hpsi,wfn%psi)
   if (nproc==1) then
      wfn%psit => wfn%psi
      ! workarray not used for nporc==1, so pass the sa,e address
      call transpose_v(iproc,nproc,wfn%orbs,wfn%lzd%glr%wfd,wfn%comms,wfn%psit,wfn%psit)
   end if

   call subspace_diagonalisation(iproc,nproc,wfn%orbs,wfn%comms,wfn%psit,wfn%hpsi,evsum)

   !here we should preserve hpsi and transpose it if we are in ensemble mimimization scheme

   call untranspose_v(iproc,nproc,wfn%orbs,wfn%Lzd%Glr%wfd,wfn%comms,&
        wfn%psit,wfn%hpsi,out_add=wfn%psi)
   ! Emit that new wavefunctions are ready.
   if (wfn%c_obj /= 0) then
      call kswfn_emit_psi(wfn, iter, 0, iproc, nproc)
   end if

   if(.not.  keeppsit) then
      if (nproc > 1  ) then
         call f_free_ptr(wfn%psit)
      else
         nullify(wfn%psit)
      end if

      call f_free_ptr(wfn%hpsi)

   endif

   !call eigensystem_info(iproc,nproc,wfn%Lzd%Glr%wfd%nvctr_c+7*wfn%Lzd%Glr%wfd%nvctr_f,&
   !     wfn%orbs,wfn%psi)
   call f_release_routine()
END SUBROUTINE last_orthon

subroutine eigensystem_info(iproc,nproc,tolerance,nvctr,orbs,psi)
  use module_base
  use module_types
  use module_interfaces, only: write_eigenvalues_data
  implicit none
  integer, intent(in) :: iproc,nproc,nvctr
  real(gp), intent(in) :: tolerance !< threshold to classify degenerate eigenvalues
  type(orbitals_data), intent(inout) :: orbs
  real(wp), dimension(nvctr,orbs%nspinor,orbs%norbp), intent(in) :: psi
  !local variables
  character(len=*), parameter :: subname='eigensystem_info'
  real(wp), dimension(:,:,:), pointer :: mom_vec

  call f_routine(id=subname)

  !for a non-collinear treatment,
  !we add the calculation of the moments for printing their value
  !close to the corresponding eigenvector
  if(orbs%nspinor==4) then
     mom_vec = f_malloc_ptr((/ 4, orbs%norb, min(nproc, 2) /),id='mom_vec')

     call calc_moments(iproc,nproc,orbs%norb,orbs%norb_par,&
          nvctr,&
          orbs%nspinor,psi,mom_vec)
  else
     nullify(mom_vec)
  end if

  ! Send all eigenvalues to all procs.
  !! needed only if nkpts > 1
  call broadcast_kpt_objects(nproc,orbs%nkpts,orbs%norb, &
       orbs%eval,orbs%ikptproc)

  !here the new occupation numbers should be recalculated for future needs

  !print the found eigenvalues
  if (iproc == 0) then
     call write_eigenvalues_data(tolerance,orbs,mom_vec)
  end if

  if (orbs%nspinor ==4) then
     call f_free_ptr(mom_vec)
  end if
  call f_release_routine()
end subroutine eigensystem_info

!> find the gap once the fermi level has been found
 subroutine orbs_get_gap(orbs,ikpt_homo,ikpt_lumo,ispin_homo,ispin_lumo,homo,lumo,occup_lumo)
   use module_base
   use module_defs, only: gp
   use module_types, only: orbitals_data
   use dictionaries, only: f_err_throw
   implicit none
   integer, intent(out) :: ikpt_lumo,ikpt_homo,ispin_homo,ispin_lumo
   real(gp), intent(out) :: homo,lumo,occup_lumo
   type(orbitals_data), intent(inout) :: orbs
   !local variables
   integer :: ikpt,iorb,jorb

   call f_routine(id='orbs_get_gap')
   homo=-1.e100_gp
   lumo=1.e100_gp
   occup_lumo=2.d0
   do ikpt=1,orbs%nkpts
      find_up: do iorb=1,orbs%norbu
         jorb=iorb+(ikpt-1)*orbs%norb
         if (orbs%occup(jorb) < 0.5_gp) then
            if (iorb ==1) call f_err_throw(&
                 'Fermi level badly calculated, first orbital is already above',&
                 err_name='BIGDFT_RUNTIME_ERROR')
            if (orbs%eval(jorb-1) > homo) then
               homo=orbs%eval(jorb-1)
               ikpt_homo=ikpt
               ispin_homo=1
            end if
            if (orbs%eval(jorb) < lumo) then
               lumo=orbs%eval(jorb)
               occup_lumo=min(occup_lumo,orbs%occup(jorb))
               ikpt_lumo=ikpt
               ispin_lumo=1
            end if
            exit find_up
         end if
      end do find_up
      find_down: do iorb=1,orbs%norbd
         jorb=orbs%norbu+iorb+(ikpt-1)*orbs%norb
         if (orbs%occup(jorb) < 0.5_gp) then
            if (iorb ==1) call f_err_throw(&
                 'Fermi level badly calculated (down spin), first orbital is already above',&
                 err_name='BIGDFT_RUNTIME_ERROR')
            if (orbs%eval(jorb-1) > homo) then
               homo=orbs%eval(jorb-1)
               ikpt_homo=ikpt
               ispin_homo=-1
            end if
            if (orbs%eval(jorb) < lumo) then
               lumo=orbs%eval(jorb)
               occup_lumo=min(occup_lumo,orbs%occup(jorb))
               ikpt_lumo=ikpt
               ispin_lumo=-1
            end if
            exit find_down
         end if
      end do find_down
   end do
   !now verify that the gap has been found
   !if (lumo < homo) call f_err_throw('Error in determining homo-lumo gap',&
   !     err_name='BIGDFT_RUNTIME_ERROR')
   call f_release_routine()
 end subroutine orbs_get_gap


!!subroutine eFermi_nosmearing(iproc,orbs)
!!   use module_base
!!   use module_types
!!   use yaml_output
!!   implicit none
!!   integer, intent(in) :: iproc
!!   type(orbitals_data), intent(inout) :: orbs
!!   !local variables
!!   integer :: iu,id,n,nzeroorbs,ikpt,iorb
!!   real(gp) :: charge
!!   real(wp) :: eF
!!
!!   call f_routine(id='eFermi_nosmearing')
!!   !SM: I think iu and id should be initialized to these values, in case the
!!   ! large if will not be executed.
!!   iu=orbs%norbu
!!   id=orbs%norbd
!!   eF = 0._wp
!!   do ikpt=1,orbs%nkpts
!!      !number of zero orbitals for the given k-point
!!      nzeroorbs=0
!!      !overall charge of the system
!!      charge=0.0_gp
!!      do iorb=1,orbs%norb
!!         if (orbs%occup(iorb+(ikpt-1)*orbs%norb) == 0.0_gp) then
!!            nzeroorbs=nzeroorbs+1
!!         else
!!            charge=charge+orbs%occup(iorb+(ikpt-1)*orbs%norb)
!!         end if
!!      end do
!!      if (nzeroorbs /= 0 .and. orbs%norbd .gt.0) then
!!         do iorb=1,orbs%norbu-1
!!            if (orbs%eval((ikpt-1)*orbs%norb+iorb) > orbs%eval((ikpt-1)*orbs%norb+iorb+1)) &
!!               &   write(*,*) 'wrong ordering of up EVs',iorb,iorb+1
!!         end do
!!         do iorb=1,orbs%norbd-1
!!            if (orbs%eval((ikpt-1)*orbs%norb+iorb+orbs%norbu) > orbs%eval((ikpt-1)*orbs%norb+iorb+1+orbs%norbu))&
!!               &   write(*,*) 'wrong ordering of dw EVs',iorb+orbs%norbu,iorb+1+orbs%norbu
!!         enddo
!!
!!         iu=0
!!         id=0
!!         n=0
!!         do while (real(n,gp) < charge)
!!            if (orbs%eval((ikpt-1)*orbs%norb+iu+1) <= orbs%eval((ikpt-1)*orbs%norb+id+1+orbs%norbu)) then
!!               iu=iu+1
!!               eF=orbs%eval((ikpt-1)*orbs%norb+iu+1)
!!            else
!!               id=id+1
!!               eF=orbs%eval((ikpt-1)*orbs%norb+id+1+orbs%norbu)
!!            endif
!!            n=n+1
!!         enddo
!!         if (iproc==0) then
!!            !write(*,'(1x,a,1pe21.14,a,i4)') 'Suggested Homo energy level',eF,', Spin polarization',iu-id
!!            call yaml_map('Suggested Fermi Level',ef,fmt='(1pe21.14)')
!!            call yaml_map('Suggested Spin pol.',iu-id,fmt='(i4)')
!!         end if
!!         !write(*,*) 'up,down, up-down',iu,id,iu-id
!!      end if
!!   end do
!!   orbs%efermi=eF
!!   !assign the values for the occupation numbers
!!   do iorb=1,iu
!!      orbs%occup(iorb)=1.0_gp
!!   end do
!!   do iorb=iu+1,orbs%norbu
!!      orbs%occup(iorb)=0.0_gp
!!   end do
!!   do iorb=1,id
!!      orbs%occup(iorb+orbs%norbu)=1.0_gp
!!   end do
!!   do iorb=id+1,orbs%norbd
!!      orbs%occup(iorb+orbs%norbu)=0.0_gp
!!   end do
!!
!!   call f_release_routine()
!!
!!END SUBROUTINE eFermi_nosmearing



!>   Calculate magnetic moments
subroutine calc_moments(iproc,nproc,norb,norb_par,nvctr,nspinor,psi,mom_vec)
   use module_base
   use module_types
   implicit none
   integer, intent(in) :: iproc,nproc,norb,nvctr,nspinor
   integer, dimension(0:nproc-1), intent(in) :: norb_par
   real(wp), dimension(nvctr,norb*nspinor), intent(in) :: psi
   real(wp), dimension(4,norb,min(nproc,2)), intent(out) :: mom_vec
   !local variables
   character(len=*), parameter :: subname='calc_moments'
   integer :: ierr,iorb,jproc
   integer :: ndim,oidx
   integer, dimension(:), allocatable :: norb_displ
   real(wp) :: m00,m11,m13,m24,m14,m23
   !real(wp), dimension(:,:,:), allocatable :: mom_vec

   call f_routine(id=subname)
   ndim=2
   if (nproc==1) ndim=1

   if(nspinor==4) then

      call f_zero(mom_vec)

      do iorb=1,norb_par(iproc)
         oidx=(iorb-1)*nspinor+1
         m00=dot(2*nvctr,psi(1,oidx),1,psi(1,oidx),1)
         m11=dot(2*nvctr,psi(1,oidx+2),1,psi(1,oidx+2),1)
         m13=dot(nvctr,psi(1,oidx),1,psi(1,oidx+2),1)
         m24=dot(nvctr,psi(1,oidx+1),1,psi(1,oidx+3),1)
         !        m12=dot(nvctr,psi(1,oidx),1,psi(1,oidx+1),1)
         !        m34=dot(nvctr,psi(1,oidx+2),1,psi(1,oidx+3),1)
         m14=dot(nvctr,psi(1,oidx),1,psi(1,oidx+3),1)
         m23=dot(nvctr,psi(1,oidx+1),1,psi(1,oidx+2),1)

         mom_vec(1,iorb,ndim)=(m00+m11) !rho
         mom_vec(2,iorb,ndim)=2.0d0*(m13+m24)       !m_x
         !        mom_vec(3,iorb,ndim)=2.0d0*(m12-m34)       !m_y
         mom_vec(3,iorb,ndim)=2.0d0*(m14-m23)       !m_y
         mom_vec(4,iorb,ndim)=(m00-m11)             !m_z
      end do

      if(nproc>1) then
         norb_displ = f_malloc(0.to.nproc-1,id='norb_displ')

         norb_displ(0)=0
         do jproc=1,nproc-1
            norb_displ(jproc)=norb_displ(jproc-1)+norb_par(jproc-1)
         end do

         call MPI_GATHERV(mom_vec(1,1,2),4*norb_par(iproc),mpidtypw,&
            &   mom_vec(1,1,1),4*norb_par,4*norb_displ,mpidtypw,&
         0,bigdft_mpi%mpi_comm,ierr)

         call f_free(norb_displ)
      end if

   end if
   call f_release_routine()
END SUBROUTINE calc_moments


!> Check communications (transpose and untranspose orbitals)
subroutine check_communications(iproc,nproc,orbs,lzd,comms)
   use module_base
   use module_types
   use yaml_output
   use communications_base, only: comms_cubic
   use communications, only: transpose_v, untranspose_v
   implicit none
   integer, intent(in) :: iproc,nproc
   type(orbitals_data), intent(in) :: orbs
   type(local_zone_descriptors), intent(in) :: lzd
   type(comms_cubic), intent(in) :: comms
   !local variables
   character(len=*), parameter :: subname='check_communications'
   integer :: i,ispinor,iorb,indspin,indorb,jproc,iscomp,idsx,index,ikptsp
   integer :: ikpt,ispsi,nspinor,nvctrp,ierr
   real(wp) :: psival,maxdiff
   real(wp), dimension(:), allocatable :: psi
   real(wp), dimension(:), pointer :: pwork
   real(wp) :: epsilon
   character(len = 25) :: filename
   logical :: abort

   call f_routine(id=subname)
   if(bigdft_mpi%iproc==0) call yaml_mapping_open('Communication checks')

   !allocate the "wavefunction" and fill it, and also the workspace
   psi = f_malloc(max(orbs%npsidim_orbs, orbs%npsidim_comp),id='psi')
   pwork = f_malloc_ptr(max(orbs%npsidim_orbs, orbs%npsidim_comp),id='pwork')

   do iorb=1,orbs%norbp
      ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
      !valkpt=real(512*ikpt,wp)
      !valorb=real(orbs%isorb+iorb-(ikpt-1)*orbs%norb,wp)+valkpt
      indorb=(iorb-1)*(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f)*orbs%nspinor
      do ispinor=1,orbs%nspinor
         indspin=(ispinor-1)*(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f)
         do i=1,lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f
            !vali=real(i,wp)/512.0_wp  ! *1.d-5
            call test_value(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
            psi(i+indspin+indorb)=psival!(valorb+vali)*(-1)**(ispinor-1)
         end do
      end do
   end do

   !transpose the hpsi wavefunction
   call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi,pwork)

   !check the results of the transposed wavefunction
   maxdiff=0.0_wp
   ispsi=0
   do ikptsp=1,orbs%nkptsp
      ikpt=orbs%iskpts+ikptsp!orbs%ikptsp(ikptsp)
      !valkpt=real(512*ikpt,wp)
      !calculate the starting point for the component distribution
      iscomp=0
      do jproc=0,iproc-1
         iscomp=iscomp+comms%nvctr_par(jproc,ikpt)
      end do
      nvctrp=comms%nvctr_par(iproc,ikpt)
      nspinor=orbs%nspinor

      do iorb=1,orbs%norb
         !valorb=real(iorb,wp)+valkpt
         indorb=(iorb-1)*nvctrp*nspinor
         do idsx=1,(nspinor-1)/2+1
            do i=1,nvctrp
               !vali=real(i+iscomp,wp)/512.d0  ! *1.d-5
               do ispinor=1,((2+nspinor)/4+1)
                  !psival=(-1)**(ispinor-1)*(valorb+vali)
                  call test_value(ikpt,iorb,ispinor,i+iscomp,psival)
                  !this is just to force the IEEE representation of psival
                  !              if (psival .lt. 0.d0) then
                  !              write(321,*) psival,psival**2
                  !              endif
                  index=ispinor+(i-1)*((2+nspinor)/4+1)+&
                     &   (idsx-1)*((2+nspinor)/4+1)*nvctrp+indorb+ispsi
                  maxdiff=max(abs(psi(index)-psival),maxdiff)
               end do
            end do
         end do
      end do
      ispsi=ispsi+nvctrp*orbs%norb*nspinor
   end do

   abort = .false.
   if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp) .and. maxdiff < 1.e-6_wp) then
      call yaml_warning('Transposition error of'//&
           trim(yaml_toa(maxdiff,fmt='(1pe15.7)'))//', check whether results are meaningful!')
   else if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp)) then
      write(*,*)'ERROR: process',iproc,'does not transpose wavefunctions correctly!'
      write(*,*)'       found an error of',maxdiff,'cannot continue.'
      write(*,*)'       data are written in the file transerror.log, exiting...'

      write(filename, "(A,I0,A)") 'transerror', iproc, '.log'
      open(unit=22,file=trim(filename),status='unknown')
      ispsi=0
      do ikptsp=1,orbs%nkptsp
         ikpt=orbs%iskpts+ikptsp!orbs%ikptsp(ikptsp)
         !valkpt=real(512*ikpt,wp)
         !calculate the starting point for the component distribution
         iscomp=0
         do jproc=0,iproc-1
            iscomp=iscomp+comms%nvctr_par(jproc,ikpt)
         end do
         nvctrp=comms%nvctr_par(iproc,ikpt)
         nspinor=orbs%nspinor

         do iorb=1,orbs%norb
            !valorb=real(iorb,wp)+valkpt
            indorb=(iorb-1)*nvctrp*nspinor
            do idsx=1,(nspinor-1)/2+1
               do i=1,nvctrp
                  !vali=real(i+iscomp,wp)/512.d0  !*1.d-5
                  do ispinor=1,((2+nspinor)/4+1)
                     !psival=(-1)**(ispinor-1)*(valorb+vali)
                     call test_value(ikpt,iorb,ispinor,i+iscomp,psival)
                     index=ispinor+(i-1)*((2+nspinor)/4+1)+&
                        &   (idsx-1)*((2+nspinor)/4+1)*nvctrp+indorb+ispsi
                     maxdiff=abs(psi(index)-psival)
                     if (maxdiff > 0.d0) then
                        write(22,'(i3,i6,2i4,3(1x,1pe13.6))')ispinor,i+iscomp,iorb,ikpt,psival,&
                           &   psi(index),maxdiff
                     end if
                  end do
               end do
            end do
         end do
         ispsi=ispsi+nvctrp*orbs%norb*nspinor
      end do
      close(unit=22)
      abort = .true.
      write(filename, "(A,I0,A)") 'distscheme', iproc, '.log'
      !open(unit=22,file=trim(filename),status='unknown')
      call print_distribution_schemes(nproc,orbs%nkpts,orbs%norb_par(0,1),comms%nvctr_par(0,1))
      !close(unit=22)
   end if

   if(bigdft_mpi%iproc==0) call yaml_map('Transpositions', .not. abort)
   call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
   if (abort) then
      if (iproc == 0) call print_distribution_schemes(nproc,orbs%nkpts,orbs%norb_par(0,1),comms%nvctr_par(0,1))
      call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)
   end if

   !retranspose the hpsi wavefunction
   call untranspose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,&
      &   psi,pwork)

   maxdiff=0.0_wp
   do iorb=1,orbs%norbp
      ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
      !valkpt=real(512*ikpt,wp)
      !valorb=real(orbs%isorb+iorb-(ikpt-1)*orbs%norb,wp)+valkpt
      indorb=(iorb-1)*(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f)*orbs%nspinor
      do ispinor=1,orbs%nspinor
         indspin=(ispinor-1)*(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f)
         do i=1,lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f
            !vali=real(i,wp)/512.d0  !*1.d-5
            !psival=(valorb+vali)*(-1)**(ispinor-1)
            call test_value(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
            maxdiff=max(abs(psi(i+indspin+indorb)-psival),maxdiff)
         end do
      end do
   end do

   abort = .false.
   if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp) .and. maxdiff < 1.e-6_wp) then
      call yaml_warning('Inverse transposition error of'//&
           trim(yaml_toa(maxdiff,fmt='(1pe15.7)'))//', check whether results are meaningful!')
   else if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp)) then
      write(*,*)'ERROR: process',iproc,'does not untranspose wavefunctions correctly!'
      write(*,*)'       found an error of',maxdiff,'cannot continue.'
      write(*,*)'       data are written in the file transerror.log, exiting...'

      write(filename, "(A,I0,A)") 'transerror', iproc, '.log'
      open(unit=22,file=trim(filename),status='unknown')
      maxdiff=0.0_wp
      do iorb=1,orbs%norbp
         ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
         !valkpt=real(512*ikpt,wp)
         !valorb=real(orbs%isorb+iorb-(ikpt-1)*orbs%norb,wp)+valkpt
         indorb=(iorb-1)*(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f)*orbs%nspinor
         do ispinor=1,orbs%nspinor
            indspin=(ispinor-1)*(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f)
            do i=1,lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f
               !vali=real(i,wp)/512.d0  !*1.d-5
               !psival=(valorb+vali)*(-1)**(ispinor-1)
               call test_value(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
               maxdiff=abs(psi(i+indspin+indorb)-psival)
               if (maxdiff > 0.d0) then
                  write(22,'(i3,i6,2i4,3(1x,1pe13.6))')ispinor,i,iorb,orbs%isorb,psival,&
                     &   psi(ispinor+(i-1)*orbs%nspinor+indorb),maxdiff
               end if
            end do
         end do
      end do
      close(unit=22)
      abort = .true.
   end if

   if(bigdft_mpi%iproc==0) call yaml_map('Reverse transpositions', .not. abort)
   call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
   if (abort) call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)

   call f_free(psi)
   call f_free_ptr(pwork)

   if(bigdft_mpi%iproc==0) call yaml_mapping_close()

   call f_release_routine()
END SUBROUTINE check_communications


!> define a value for the wavefunction which is dependent of the indices
subroutine test_value(ikpt,iorb,ispinor,icomp,val)
   use module_base
   implicit none
   integer, intent(in) :: ikpt,icomp,iorb,ispinor
   real(wp), intent(out) :: val
   !local variables
   integer, parameter :: ilog=6
   real(wp) :: valkpt,valorb,vali

   ! recognizable pattern, for debugging
    valkpt=real(10**ilog*(ikpt-1),wp)!real(512*ikpt,wp)
    valorb=real(iorb,wp)+valkpt
    vali=real(icomp,wp)*10.0_wp**(-ilog)  !real(icomp,wp)/512.0_wp  ! *1.d-5
   !
   ! val=(valorb+vali)*(-1)**(ispinor-1)

   !valkpt=real(512*ikpt,wp)
   !valorb=real(iorb,wp)+valkpt
   !vali=real(icomp,wp)/512.0_wp  ! *1.d-5

   val=(valorb+vali)*(-1)**(ispinor-1)

END SUBROUTINE test_value

!> Determine the components which were not communicated correctly
!! works only with the recognizable pattern of test function
subroutine wrong_components(psival,ikpt,iorb,icomp)
   use module_base
  implicit none
  real(wp), intent(in) :: psival
  integer, intent(out) :: ikpt,iorb,icomp
  integer, parameter :: ilog=6

  icomp=nint((psival-real(floor(psival),wp))*10.0_wp**ilog)
  ikpt=floor(psival)/(10**ilog)
  iorb=floor(psival)-(ikpt-1)*(10**ilog)

end subroutine wrong_components


subroutine debug_hpsi(wfn,at)!,proj_G,paw,nlpspd)
   use module_base
   use module_types
   implicit none
   !Arguments
   !type(nonlocal_psp_descriptors), intent(in) :: nlpspd
   type(DFT_wavefunction), intent(in) :: wfn
   !type(paw_objects),intent(in) :: paw
   type(atoms_data), intent(in) :: at
   !type(gaussian_basis),dimension(at%astruct%ntypes),intent(in)::proj_G !projectors in gaussian basis (for PAW)
   !Local variables
   integer :: ispsi_a,ispsi_b,ia,ib,iatype,ispinor
   integer :: nvctr_c,nvctr_f,nvctr_tot
   integer :: ncplx,nspinor
   !integer :: i,istart_j,jlmn,j_shell,j_l,j_m
   !integer :: mbseg_c,mbseg_f,mbvctr_c,mbvctr_f,istart_ck
   integer :: ikpt,ispsi_k,iat,ilr,ilr_skip,iorb
   integer :: isorb,ieorb,ispsi
   logical :: dosome
   real(gp),dimension(2) :: scpr
   real(gp) :: scalprod(2,2)
   real(dp) :: ddot

   write(*,*)'Calculate <PSI|H|PSI> for debugging'

   nvctr_c=wfn%Lzd%glr%wfd%nvctr_c
   nvctr_f=wfn%Lzd%glr%wfd%nvctr_f
   nvctr_tot=nvctr_c+7*nvctr_f

   ikpt=wfn%orbs%iokpt(1)

   ispsi_k=1
   loop_kpt: do

      call orbs_in_kpt(ikpt,wfn%orbs,isorb,ieorb,nspinor)

      loop_lr: do ilr=1,wfn%Lzd%nlr
         !do something only if at least one of the orbitals lives in the ilr
         dosome=.false.
         do iorb=isorb,ieorb
            dosome = (wfn%orbs%inwhichlocreg(iorb+wfn%orbs%isorb) == ilr)
            if (dosome) exit
         end do
         if (.not. dosome) cycle loop_lr

      do iat=1,at%astruct%nat
           iatype=at%astruct%iatype(iat)
              ispsi=ispsi_k
            do iorb=isorb,ieorb
               if (wfn%orbs%inwhichlocreg(iorb+wfn%orbs%isorb) /= ilr) then
                  !increase ispsi to meet orbital index
                  ilr_skip=wfn%orbs%inwhichlocreg(iorb+wfn%orbs%isorb)
                  ispsi=ispsi+(wfn%Lzd%Llr(ilr_skip)%wfd%nvctr_c+7*wfn%Lzd%Llr(ilr_skip)%wfd%nvctr_f)*nspinor
                  cycle
               end if

               !call plr_segs_and_vctrs(nlpspd%plr(iat),mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
               call ncplx_kpt(wfn%orbs%iokpt(iorb),wfn%orbs,ncplx)
               !loop over all the components of the wavefunction
               do ispinor=1,nspinor,ncplx
                 ispsi_a=ispsi
                 ispsi_b=ispsi
                 do ia=1,ncplx
                    do ib=1,ncplx
                       scalprod(ia,ib)=ddot(nvctr_tot,wfn%psi(ispsi_a),1,wfn%hpsi(ispsi_b),1)
                       ispsi_b=ispsi_b+nvctr_tot
                    end do
                    ispsi_a=ispsi_a+nvctr_tot
                 end do
                 !then define the result
                 if (ncplx == 1) then
                    scpr(1)=scalprod(1,1)
                 else if (ncplx == 2) then
                    scpr(1)=scalprod(1,1)+scalprod(2,2)
                    scpr(2)=scalprod(1,2)-scalprod(2,1)
                 end if
              write(*,'("<psi|H|psi> for kpt=,",i5," iat=",i5,"=>",2(f15.5,1x))')ikpt,iat,scpr
                 write(*,*)'ispsi_a','ispsi_b','nvctr_tot',ispsi_a,ispsi_b,nvctr_tot
               end do !ispinor
            end do !iorb
         end do !iat
      end do loop_lr

      !last k-point has been treated
      if (ieorb == wfn%orbs%norbp) exit loop_kpt

      ikpt=ikpt+1
      ispsi_k=ispsi

   end do  loop_kpt !ikpt

   !if (.not. DistProjApply) then !TO BE REMOVED WITH NEW PROJECTOR APPLICATION
   !   if (istart_ck-1 /= nlpspd%nprojel) &
   !      &   stop 'incorrect once-and-for-all psp application'
   !end if

END SUBROUTINE debug_hpsi


subroutine broadcast_kpt_objects(nproc, nkpts, ndata, data, ikptproc)
   use module_base
   use module_types
   implicit none
   integer, intent(in) :: nproc, nkpts, ndata
   integer, dimension(nkpts), intent(in) :: ikptproc
   real(gp), dimension(ndata,nkpts), intent(inout) :: data

   integer :: ikpt


   if (nproc > 1) then
      call mpibarrier(comm=bigdft_mpi%mpi_comm)
      do ikpt = 1, nkpts
         !print *,'data(:),ikpt',data(:,ikpt),ikpt,bigdft_mpi%iproc,ikptproc(ikpt)
         call mpibcast(data(:,ikpt),root=ikptproc(ikpt),&
              comm=bigdft_mpi%mpi_comm)!,check=.true.)
         !call MPI_BCAST(data(1,ikpt), ndata,mpidtypg, &
         !   &   ikptproc(ikpt), bigdft_mpi%mpi_comm, ierr)
         !redundant barrier
         call mpibarrier(comm=bigdft_mpi%mpi_comm)
         !call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
      end do
   end if
END SUBROUTINE broadcast_kpt_objects


!!subroutine minimize_by_orthogonal_transformation(iproc, nproc, orbs, wfd, comms, orthpar, E0, El, stepsize, hpsi, psi, derivative)
!!use module_base
!!use module_types
!!use module_interfaces
!!implicit none
!!
!!  integer, intent(in) :: iproc,nproc
!!  type(orbitals_data), intent(in) :: orbs
!!  type(wavefunctions_descriptors),intent(in):: wfd
!!  type(comms_cubic), intent(in) :: comms
!!  type(orthon_data), intent(in) :: orthpar
!!  !n(c) type(wavefunctions_descriptors), intent(in) :: wfd
!!  !real(wp), dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb), intent(in) :: hpsi
!!  !real(wp), dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb), intent(inout) :: psi
!!  real(8),intent(in):: E0, El
!!  real(8),intent(inout):: stepsize
!!  real(wp), dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)), intent(inout) :: hpsi
!!  real(wp), dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)), intent(inout) :: psi
!!  real(8),intent(out):: derivative
!!  !local variables
!!  character(len=*), parameter :: subname='orthoconstraint'
!!  integer :: i_stat,i_all,ierr,iorb !n(c) ise
!!  integer :: ispin,nspin,ikpt,norb,norbs,ncomp,nvctrp,ispsi,ikptp,nspinor,iiorb,jjorb,jorb,lwork,info,k
!!  real(dp) :: occ !n(c) tt
!!  integer, dimension(:,:), allocatable :: ndimovrlp
!!  real(wp), dimension(:), allocatable :: alag, work, rwork, eval
!!  real(wp),dimension(:,:),allocatable:: gmat
!!  complex(8),dimension(:,:),allocatable:: gmatc, tempmatc, tempmat2c, omatc
!!  complex(8),dimension(:),allocatable:: expDc
!!  real(8):: lstep, dfactorial
!!  complex(8):: ttc
!!  real(8),dimension(:),pointer:: psiwork
!!  real(8),dimension(:,:),allocatable:: omat
!!
!!  write(*,*) 'iproc, orbs%npsidim',iproc,max(orbs%npsidim_orbs,orbs%npsidim_comp)
!!
!!  !separate the orthogonalisation procedure for up and down orbitals
!!  !and for different k-points
!!  call timing(iproc,'LagrM_comput  ','ON')
!!
!!  !number of components of the overlap matrix for parallel case
!!  !calculate the dimension of the overlap matrix for each k-point
!!  if (orbs%norbd > 0) then
!!     nspin=2
!!  else
!!     nspin=1
!!  end if
!!
!!  !number of components for the overlap matrix in wp-kind real numbers
!!
!!  allocate(ndimovrlp(nspin,0:orbs%nkpts),stat=i_stat)
!!  call memocc(i_stat,ndimovrlp,'ndimovrlp',subname)
!!
!!  call dimension_ovrlp(nspin,orbs,ndimovrlp)
!!
!!  allocate(alag(ndimovrlp(nspin,orbs%nkpts) ),stat=i_stat)
!!  call memocc(i_stat,alag,'alag',subname)
!!
!!  !put to zero all the k-points which are not needed
!!  call f_zero(ndimovrlp(nspin,orbs%nkpts),alag)
!!
!!  ! Transpose orbitals
!!  allocate(psiwork(size(psi)), stat=i_stat)
!!  call transpose_v(iproc, nproc, orbs, wfd, comms, psi, work=psiwork)
!!  call transpose_v(iproc, nproc, orbs, wfd, comms, hpsi, work=psiwork)
!!
!!  !do it for each of the k-points and separate also between up and down orbitals in the non-collinear case
!!  ispsi=1
!!  do ikptp=1,orbs%nkptsp
!!     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
!!
!!     do ispin=1,nspin
!!
!!        call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
!!             nvctrp,norb,norbs,ncomp,nspinor)
!!        if (nvctrp == 0) cycle
!!
!!        if(nspinor==1) then
!!           call gemm('T','N',norb,norb,nvctrp,1.0_wp,psi(ispsi),&
!!                max(1,nvctrp),hpsi(ispsi),max(1,nvctrp),0.0_wp,&
!!                alag(ndimovrlp(ispin,ikpt-1)+1),norb)
!!        else
!!           !this part should be recheck in the case of nspinor == 2
!!           call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(ispsi),&
!!                max(1,ncomp*nvctrp), &
!!                hpsi(ispsi),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
!!                alag(ndimovrlp(ispin,ikpt-1)+1),norb)
!!        end if
!!        ispsi=ispsi+nvctrp*norb*nspinor
!!     end do
!!  end do
!!
!!  if (nproc > 1) then
!!     call timing(iproc,'LagrM_comput  ','OF')
!!     call timing(iproc,'LagrM_commun  ','ON')
!!     call fmpi_allreduce(alag(1),ndimovrlp(nspin,orbs%nkpts),FMPI_SUM,bigdft_mpi%mpi_comm,ierr)
!!     call timing(iproc,'LagrM_commun  ','OF')
!!     call timing(iproc,'LagrM_comput  ','ON')
!!  end if
!!
!!  !now each processors knows all the overlap matrices for each k-point
!!  !even if it does not handle it.
!!  !this is somehow redundant but it is one way of reducing the number of communications
!!  !without defining group of processors
!!
!!  i_stat=0
!!  do iorb=1,orbs%norb
!!      do jorb=1,orbs%norb
!!          i_stat=i_stat+1
!!          if(iproc==0) write(1750,'(a,3i8,es20.8)') 'iorb, jorb, i_stat, alag(i_stat)', iorb, jorb, i_stat, alag(i_stat)
!!      end do
!!  end do
!!
!!  ! Build the antisymmetric matrix "alag(iorb,jorb)-alag(jorb,iorb)"
!!  allocate(gmat(orbs%norb,orbs%norb), stat=i_stat)
!!  call memocc(i_stat, gmat, 'gmat', subname)
!!  do iorb=1,orbs%norb
!!      do jorb=1,orbs%norb
!!          iiorb=(iorb-1)*orbs%norb+jorb
!!          jjorb=(jorb-1)*orbs%norb+iorb
!!          gmat(jorb,iorb)=alag(iiorb)-alag(jjorb)
!!      end do
!!  end do
!!
!!
!!
!!
!!  !Build the complex matrix -igmat
!!  allocate(gmatc(orbs%norb,orbs%norb), stat=i_stat)
!!  !call memocc(i_stat, gmatc, 'gmatc', subname) !memocc not working for complex
!!  do iorb=1,orbs%norb
!!      do jorb=1,orbs%norb
!!          gmatc(jorb,iorb)=cmplx(0.d0,-gmat(jorb,iorb),kind=8)
!!      end do
!!  end do
!!
!!
!!
!!  ! Check whether hermitian
!!  do iorb=1,orbs%norb
!!      do jorb=1,orbs%norb
!!          if(real(gmatc(jorb,iorb))/=real(gmatc(iorb,jorb)) .or. aimag(gmatc(jorb,iorb))/=-aimag(gmatc(iorb,jorb))) then
!!              write(*,'(a,4es16.7)') 'ERROR: (Gmatc(1,jorb,iorb)), (Gmatc(1,iorb,jorb)), (Gmatc(2,jorb,iorb)), (Gmatc(2,iorb,jorb))', Gmatc(jorb,iorb), Gmatc(iorb,jorb)
!!          end if
!!          if(iproc==0) write(1710,'(a,2i8,2es20.12)') 'iorb, jorb, gmatc(jorb,iorb)', iorb, jorb, gmatc(jorb,iorb)
!!      end do
!!  end do
!!
!!
!!
!!  ! Diagonalize gmatc
!!  allocate(eval(orbs%norb), stat=i_stat)
!!  call memocc(i_stat, eval, 'eval', subname)
!!  lwork=10*orbs%norb
!!  allocate(work(2*lwork), stat=i_stat) ! factor of 2 since it is assumed to be complex
!!  allocate(rwork(lwork), stat=i_stat)
!!  call zheev('v', 'l', orbs%norb, gmatc(1,1), orbs%norb, eval(1), work, lwork, rwork, info)
!!  if(info/=0) stop 'ERROR in zheev'
!!  deallocate(work)
!!  deallocate(rwork)
!!
!!  do iorb=1,orbs%norb
!!      if(iproc==0) write(1720,'(a,i8,es20.8)') 'iorb, eval(iorb)', iorb, eval(iorb)
!!  end do
!!
!!
!!  if(stepsize>0.d0) then
!!      if(iproc==0) write(*,*) 'in first branch'
!!      ! Calculate the derivative
!!      ! Calculate optimal step size
!!      lstep=derivative*stepsize**2/(2.d0*(derivative*stepsize+E0-El))
!!      stepsize=-1.d0
!!      derivative=0.d0
!!  else
!!      if(iproc==0) write(*,*) 'in second branch'
!!      lstep=.1d0/(maxval(abs(eval)))
!!      stepsize=lstep
!!      derivative=0.d0
!!      do iorb=1,orbs%norb
!!          do jorb=1,orbs%norb
!!              derivative=derivative+gmat(jorb,iorb)**2
!!          end do
!!      end do
!!      derivative=-.5d0*derivative
!!  end if
!!  if(iproc==0) write(*,'(a,2es12.4)') '>> LSTEP:, derivative', lstep, derivative
!!
!!
!!  ! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
!!  ! This is also a diagonal matrix, so only calculate the diagonal part.
!!  allocate(expDc(orbs%norb), stat=i_stat)
!!  !call memocc(i_stat, expDc, 'expDc', subname)
!!  do iorb=1,orbs%norb
!!     ttc=cmplx(0.d0,-lstep*eval(iorb),kind=8)
!!     expDc(iorb)=cmplx(0.d0,0.d0,kind=8)
!!      do k=0,100
!!          expDc(iorb)=expDc(iorb)+ttc**k/dfactorial(k)
!!      end do
!!  end do
!!  do iorb=1,orbs%norb
!!     if(iproc==0) write(1740,'(a,i8,2es20.8)') 'iorb, expDc(iorb)', iorb, expDc(iorb)
!!  end do
!!
!!
!!
!!  ! Calculate the matrix O
!!  allocate(tempmatc(orbs%norb,orbs%norb), stat=i_stat)
!!  allocate(tempmat2c(orbs%norb,orbs%norb), stat=i_stat)
!!  allocate(omatc(orbs%norb,orbs%norb), stat=i_stat)
!!  do iorb=1,orbs%norb
!!     do jorb=1,orbs%norb
!!         if(iorb==jorb) then
!!             tempmat2c(jorb,iorb)=expDc(iorb)
!!         else
!!             tempmat2c(jorb,iorb)=cmplx(0.d0,0.d0,kind=8)
!!         end if
!!     end do
!!  end do
!!  call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmat2c(1,1), orbs%norb, &
!!      gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1), orbs%norb)
!!  do iorb=1,orbs%norb
!!     do jorb=1,orbs%norb
!!         if(iproc==0) write(1730,'(a,2i8,2es20.8)') 'iorb, jorb, tempmatc(jorb,iorb)', iorb, jorb, tempmatc(jorb,iorb)
!!     end do
!!  end do
!!
!!  call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
!!      tempmatc(1,1), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
!!
!!  if(iproc==0) then
!!     do iorb=1,orbs%norb
!!         do jorb=1,orbs%norb
!!             write(1700,'(2i8,2es20.8)') iorb,jorb,omatc(jorb,iorb)
!!         end do
!!     end do
!!  end if
!!
!!
!!  ! Build new linear combinations
!!  call memocc(i_stat, psiwork, 'psiwork', subname)
!!
!!  allocate(omat(orbs%norb,orbs%norb), stat=i_stat)
!!  call memocc(i_stat, omat, 'omat', subname)
!!
!!  do iorb=1,orbs%norb
!!      do jorb=1,orbs%norb
!!          omat(jorb,iorb)=real(omatc(jorb,iorb))
!!          !!if(iorb==jorb) then
!!          !!    omat(jorb,iorb)=1.d0
!!          !!else
!!          !!    omat(jorb,iorb)=0.d0
!!          !!end if
!!      end do
!!  end do
!!
!!  nvctrp=comms%nvctr_par(iproc,0)
!!  call gemm('n', 'n', nvctrp, orbs%norb, orbs%norb, 1.0_wp, psi(1), max(1,nvctrp), &
!!            omat(1,1), orbs%norb, 0.0_wp, psiwork(1), max(1,nvctrp))
!!  call vcopy(size(psi), psiwork(1), 1, psi(1), 1)
!!
!!  ! I think this is not required..
!!  call orthogonalize(iproc,nproc,orbs,comms,psi,orthpar)
!!
!!  call untranspose_v(iproc, nproc, orbs, wfd, comms, psi, work=psiwork)
!!  call untranspose_v(iproc, nproc, orbs, wfd, comms, hpsi, work=psiwork)
!!
!!
!!  i_all=-product(shape(omat))*kind(omat)
!!  deallocate(omat,stat=i_stat)
!!  call memocc(i_stat,i_all,'omat',subname)
!!
!!  i_all=-product(shape(psiwork))*kind(psiwork)
!!  deallocate(psiwork,stat=i_stat)
!!  call memocc(i_stat,i_all,'psiwork',subname)
!!
!!  i_all=-product(shape(alag))*kind(alag)
!!  deallocate(alag,stat=i_stat)
!!  call memocc(i_stat,i_all,'alag',subname)
!!
!!  i_all=-product(shape(ndimovrlp))*kind(ndimovrlp)
!!  deallocate(ndimovrlp,stat=i_stat)
!!  call memocc(i_stat,i_all,'ndimovrlp',subname)
!!
!!  i_all=-product(shape(gmat))*kind(gmat)
!!  deallocate(gmat,stat=i_stat)
!!  call memocc(i_stat,i_all,'gmat',subname)
!!
!!  i_all=-product(shape(eval))*kind(eval)
!!  deallocate(eval,stat=i_stat)
!!  call memocc(i_stat,i_all,'eval',subname)
!!
!!  !!memoc not working for complex
!!  !i_all=-product(shape(gmatc))*kind(gmatc)
!!  deallocate(gmatc,stat=i_stat)
!!  !call memocc(i_stat,i_all,'gmatc',subname)
!!
!!  !i_all=-product(shape(tempmatc))*kind(tempmatc)
!!  deallocate(tempmatc,stat=i_stat)
!!  !call memocc(i_stat,i_all,'tempmatc',subname)
!!
!!  !i_all=-product(shape(tempmat2c))*kind(tempmat2c)
!!  deallocate(tempmat2c,stat=i_stat)
!!  !call memocc(i_stat,i_all,'tempmat2c',subname)
!!
!!  !i_all=-product(shape(omatc))*kind(omatc)
!!  deallocate(omatc,stat=i_stat)
!!  !call memocc(i_stat,i_all,'omatc',subname)
!!
!!  !i_all=-product(shape(expDc))*kind(expDc)
!!  deallocate(expDc,stat=i_stat)
!!  !call memocc(i_stat,i_all,'expDc',subname)
!!
!!  call timing(iproc,'LagrM_comput  ','OF')
!!
!!
!!
!!end subroutine minimize_by_orthogonal_transformation


subroutine integral_equation(iproc,nproc,atoms,wfn,ngatherarr,local_potential,GPU,xc,nlpsp,rxyz,paw)
  use module_base
  use module_types
  use module_xc
  use module_interfaces, only: LocalHamiltonianApplication, plot_wf
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use yaml_output
  use yaml_parse, only: yaml_load
  use locreg_operations
  implicit none
  integer, intent(in) :: iproc,nproc
  type(atoms_data), intent(in) :: atoms
  type(DFT_wavefunction), intent(in) :: wfn
  type(GPU_pointers), intent(inout) :: GPU
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(xc_info), intent(in) :: xc
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  real(dp), dimension(:), pointer :: local_potential
  type(paw_objects), intent(inout) :: paw
  !local variables
  integer :: iorb,nbox,ilr,ist
  real(gp) :: eh_fake,eks
  type(energy_terms) :: energs_tmp
  type(coulomb_operator) :: G_Helmholtz
  type(dictionary), pointer :: dict
  type(workarr_sumrho) :: w
  real(wp), dimension(:), allocatable :: vpsi,vpsir

  call f_routine(id='helmholtz_equation')

  !first, apply the potential operator to the wavefunction
  vpsi=f_malloc0(wfn%orbs%npsidim_orbs,id='vpsi')

  call LocalHamiltonianApplication(iproc,nproc,atoms,wfn%orbs%npsidim_orbs,wfn%orbs,&
       wfn%Lzd,wfn%confdatarr,ngatherarr,local_potential,wfn%psi,vpsi,&
       energs_tmp,wfn%SIC,GPU,2,xc) !potential only

  call NonLocalHamiltonianApplication(iproc,atoms,wfn%orbs%npsidim_orbs,wfn%orbs,&
       wfn%Lzd,nlpsp,wfn%psi,vpsi,energs_tmp%eproj,paw)

  !now vpsi is a wavefunction array in orbitals parallelization scheme which is associated to Vpsi
  !rescale it to match with the Green's function treatment
  call vscal(wfn%orbs%npsidim_orbs,-0.5_gp/pi,vpsi(1),1)

  !helmholtz-based preconditioning
  ilr=1 !for the moment only cubic version
  call initialize_work_arrays_sumrho(wfn%Lzd%Llr(ilr),.true.,w)
  !box elements size
  nbox=wfn%Lzd%Llr(ilr)%d%n1i*wfn%Lzd%Llr(ilr)%d%n2i*wfn%Lzd%Llr(ilr)%d%n3i

  ! Wavefunction in real space
  vpsir=f_malloc0(nbox,id='vpsir')

  !case for helmholtz-based preconditioning
  ist=0
  do iorb=1,wfn%orbs%norbp*wfn%orbs%nspinor
     ilr = wfn%orbs%inwhichlocreg(iorb+wfn%orbs%isorb)
     !entering
     eks=wfn%orbs%eval(iorb+wfn%orbs%isorb)
     print *,'iorb,initial',iorb,nrm2(wfn%Lzd%Llr(ilr)%wfd%nvctr_c+7*wfn%Lzd%Llr(ilr)%wfd%nvctr_f,wfn%psi(1+ist),1)

!!$     call plot_wf('psi'//trim(adjustl(yaml_toa(iorb))),1,atoms,1.0_gp,wfn%Lzd%llr(ilr),&
!!$          wfn%Lzd%hgrids(1),wfn%Lzd%hgrids(2),wfn%Lzd%hgrids(3),rxyz,wfn%psi(1+ist:))

!     call axpy(wfn%Lzd%Llr(ilr)%wfd%nvctr_c+7*wfn%Lzd%Llr(ilr)%wfd%nvctr_f,-eks,wfn%psi(1+ist),1,vpsi(1+ist),1)

     call plot_wf(.false.,'Vpsi'//trim(adjustl(yaml_toa(iorb))),1,atoms,1.0_gp,wfn%Lzd%llr(ilr),&
          wfn%Lzd%hgrids,rxyz,vpsi(1+ist))


     call daub_to_isf(wfn%Lzd%Llr(ilr),w,vpsi(1+ist),vpsir(1))

     !sequential kernel

     dict => yaml_load('{kernel: {screening:'//sqrt(2.0_gp*abs(eks))//'},'//&
          'setup : { verbose: No}}')
     G_Helmholtz=pkernel_init(0,1,dict,wfn%Lzd%Llr(ilr)%geocode,&
          (/wfn%Lzd%Llr(ilr)%d%n1i,wfn%Lzd%Llr(ilr)%d%n2i,wfn%Lzd%Llr(ilr)%d%n3i/),&
          0.5_gp*wfn%Lzd%hgrids)
     call dict_free(dict)

     call pkernel_set(G_Helmholtz,verbose=.true.)

     !apply it to the gradient to smooth it
     call H_potential('D',G_Helmholtz,vpsir(1),vpsir(1),eh_fake,0.d0,.false.)

     !convert the gradient back to the locreg
     call isf_to_daub(wfn%Lzd%Llr(ilr),w,vpsir(1),vpsi(1+ist))
!     call vscal(ncplx*(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f),1.0_gp/(16.0_gp*atan(1.0_gp)),hpsi(1+ist),1)
     print *,'iorb,gradient',iorb,nrm2(wfn%Lzd%Llr(ilr)%wfd%nvctr_c+7*wfn%Lzd%Llr(ilr)%wfd%nvctr_f,vpsi(1+ist),1),eh_fake
!!$     call plot_wf('GVpsi'//trim(adjustl(yaml_toa(iorb))),1,atoms,1.0_gp,wfn%Lzd%llr(ilr),&
!!$          wfn%Lzd%hgrids(1),wfn%Lzd%hgrids(2),wfn%Lzd%hgrids(3),rxyz,vpsi(1+ist:))


     ist=ist+(wfn%Lzd%Llr(ilr)%wfd%nvctr_c+7*wfn%Lzd%Llr(ilr)%wfd%nvctr_f)

     call pkernel_free(G_Helmholtz)

  end do
  call deallocate_work_arrays_sumrho(w)

  call f_free(vpsi)
  call f_free(vpsir)

  call f_release_routine()

end subroutine integral_equation


!> Compute the Dij coefficients from the current KS potential.
subroutine paw_compute_dij(paw, at, denspot, vxc, e_paw, e_pawdc, compch_sph)
  use module_base, only: bigdft_mpi, f_err_throw
  use module_types, only: paw_objects, atoms_data, DFT_local_fields, &
       & TCAT_PAW_DIJ, TCAT_LIBPAW
  use public_enums, only: KS_POTENTIAL
  use module_defs, only: gp
  use numerics, only: Ha_eV
  use m_paw_an, only: paw_an_reset_flags
  use m_paw_ij, only: paw_ij_reset_flags
  use m_pawdij, only: pawdij
  use abi_defs_basis, only: tol14
  use abi_interfaces_add_libpaw, only: abi_pawdenpot
  use dynamic_memory
  use f_utils
  use time_profiling
  use yaml_output
  use module_xc, only: xc_isgga
  implicit none
  type(paw_objects), intent(inout) :: paw
  type(atoms_data), intent(in) :: at
  type(DFT_local_fields), intent(inout) :: denspot
  real(gp), dimension(denspot%dpbox%ndimpot, denspot%dpbox%nrhodim), intent(in) :: vxc
  real(gp), intent(out) :: e_paw, e_pawdc, compch_sph

  integer, parameter :: cplex = 1, pawprtvol = 0, pawspnorb = 0, pawxcdev = 1, &
       & enunit = 0, ipert = 0, nzlmopt = 0, option = 0
  integer :: nfft, nfftot, xclevel
  real(gp), parameter :: spnorbscl = 1._gp, charge = 0._gp, xc_denpos = tol14
  real(gp) :: ucvol
  real(gp), dimension(3), parameter :: qphon = (/ 0._gp, 0._gp, 0._gp /)
  real(gp), dimension(3,3) :: gprimd ! Used only for phonons.
  real(gp), dimension(:,:), allocatable :: xred ! Used only for phonons.
  real(gp), dimension(2) :: vpotzero

  call f_timing(TCAT_PAW_DIJ, "ON")

  call paw_an_reset_flags(paw%paw_an) ! Force the recomputation of on-site potentials
  call paw_ij_reset_flags(paw%paw_ij,self_consistent=.true.) ! Force the recomputation of Dij

  call f_zero(gprimd)
  xred = f_malloc((/ 3, at%astruct%nat /), id = "xred")

  nfft = denspot%dpbox%ndimpot
  nfftot = product(denspot%dpbox%mesh%ndims)
  ucvol = product(denspot%dpbox%mesh%ndims) *denspot%dpbox%mesh%volume_element
  xclevel = 1
  if (xc_isgga(denspot%xc)) xclevel = 2

  call f_timing(TCAT_LIBPAW, "ON")
  call abi_pawdenpot(compch_sph, e_paw, e_pawdc, ipert, denspot%xc%ixc, &
       & size(paw%pawrhoij), at%astruct%nat,&
       & denspot%dpbox%nrhodim, at%astruct%ntypes, nzlmopt, option, paw%paw_an, &
       & paw%paw_an, paw%paw_ij, at%pawang, pawprtvol, at%pawrad, &
       & paw%pawrhoij, pawspnorb, at%pawtab, pawxcdev, spnorbscl, xclevel, &
       & xc_denpos, ucvol, real(at%nzatom, gp), vpotzero = vpotzero) !, &
  !&     mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
  !&     vpotzero=vpotzero)
  call f_timing(TCAT_LIBPAW, "OF")

  !call yaml_map('average electrostatic smooth potential [Ha] , [eV]', &
  !     & (/ SUM(vpotzero(:)),SUM(vpotzero(:))*Ha_eV /))

  denspot%rhov(:) = denspot%rhov(:) + SUM(vpotzero(:))
  if(option/=1)then
     !      Fix the direct total energy (non-zero only for charged systems)
     e_paw = e_paw - SUM(vpotzero) * charge
     !      Fix the double counting total energy accordingly (for both charged AND
     !      neutral systems)
     e_pawdc = e_pawdc - SUM(vpotzero) * sum(at%nelpsp) + vpotzero(2) * charge
  end if
  if (bigdft_mpi%iproc == 0) then
     call yaml_newline()
     call yaml_map('Compensation charge on libPAW mesh', compch_sph)
     call yaml_newline()
  end if

  call f_timing(TCAT_LIBPAW, "ON")
  if (denspot%rhov_is /= KS_POTENTIAL) &
     & call f_err_throw("rhov must be KS potential here",err_name="BIGDFT_RUNTIME_ERROR")
     !stop "rhov must be KS pot here."
  call pawdij(cplex, enunit, gprimd, ipert, size(paw%pawrhoij), at%astruct%nat, nfft, nfftot, &
       & denspot%dpbox%nrhodim, at%astruct%ntypes, paw%paw_an, paw%paw_ij, at%pawang, &
       & paw%fgrtab, pawprtvol, at%pawrad, paw%pawrhoij, pawspnorb, at%pawtab, pawxcdev, &
       & qphon, spnorbscl, ucvol, charge, denspot%rhov, vxc, xred, &
       & mpi_comm_grid = bigdft_mpi%mpi_comm) !, &
  !&     natvshift=dtset%natvshift,atvshift=dtset%atvshift,fatvshift=fatvshift) !,&
  !&     mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
  call f_timing(TCAT_LIBPAW, "OF")

  call f_free(xred)

  call f_timing(TCAT_PAW_DIJ, "OF")
end subroutine paw_compute_dij


!> Compute the PAW quantities rhoij (augmentation occupancies)
!  Remember:for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
subroutine paw_compute_rhoij(paw, orbs, atoms)
  use module_base, only: bigdft_mpi, gp
  use module_types, only: orbitals_data, paw_objects, atoms_data, TCAT_LIBPAW
  use m_pawrhoij, only: pawrhoij_type, pawrhoij_init_unpacked, pawrhoij_mpisum_unpacked, &
       & pawrhoij_free, symrhoij, pawrhoij_free_unpacked
  use abi_defs_basis, only: AB7_NO_ERROR
  use m_ab6_symmetry
  use abi_interfaces_numeric, only: abi_mati3inv
  use dynamic_memory
  use time_profiling
  use abi_interfaces_add_libpaw, only: abi_pawaccrhoij
  implicit none
  type(paw_objects), intent(inout) :: paw
  type(orbitals_data), intent(in) :: orbs
  type(atoms_data), intent(in) :: atoms

  integer, parameter :: nspinor = 1, ipert = 0, pawprtvol = 0
  integer :: iatom, iorbp, isppol, ncplx, iorb
  integer, dimension(atoms%astruct%nat) :: atindx
  type(pawrhoij_type), pointer :: pawrhoij_all(:)
  integer :: errno, isym, nsym
  type(symmetry_type), pointer :: symObj
  integer, allocatable :: symrec(:,:,:)
  integer, pointer  :: sym(:,:,:), indsym(:,:,:)
  integer, pointer  :: symAfm(:)
  real(gp), pointer :: transNon(:,:)

  call f_timing(TCAT_LIBPAW, "ON")
  !Build and initialize unpacked rhoij (to be computed here)
  call pawrhoij_init_unpacked(paw%pawrhoij)

  !If pawrhoij is MPI-distributed over atomic sites, gather it
  if (size(paw%pawrhoij) < atoms%astruct%nat) then
     allocate(pawrhoij_all(atoms%astruct%nat))
  else
     pawrhoij_all => paw%pawrhoij
  end if

  do iatom = 1, atoms%astruct%nat
     atindx(iatom) = iatom
  end do

  ! Loop on band for this proc.
  do iorbp = 1, orbs%norbp
     iorb = orbs%isorb + iorbp
     isppol = 1 + modulo(iorb - 1, orbs%norb) / orbs%norbu
     !call ncplx_kpt(orbs%iokpt(iorbp),orbs,ncplx)
     ncplx = 2
     call abi_pawaccrhoij(atindx, ncplx, &
          & paw%cprj(:, iorb:iorb), paw%cprj(:, iorb:iorb),&
          & 0, isppol, size(paw%pawrhoij), atoms%astruct%nat, nspinor, &
          & orbs%occup(iorb), 1, pawrhoij_all, &
          & .false., orbs%kwgts(orbs%iokpt(iorbp)))
  end do

  !MPI: need to exchange rhoij_ between procs
  call pawrhoij_mpisum_unpacked(pawrhoij_all, bigdft_mpi%mpi_comm)

  if (atoms%astruct%sym%nsym > 0) then
     call symmetry_get_matrices_p(atoms%astruct%sym%symObj, nsym, sym, transNon, symAfm, indSym = indSym, errno = errno)
     if (errno /= AB7_NO_ERROR) stop
     call symmetry_get_from_id(symObj, atoms%astruct%sym%symObj, errno)
     if (errno /= AB7_NO_ERROR) stop
  else
     nsym = 1
     sym = f_malloc0_ptr((/ 3, 3, 1 /), id = "sym")
     sym(1,1,1) = 1
     sym(2,2,1) = 1
     sym(3,3,1) = 1
     indsym = f_malloc_ptr((/ 4, nsym, atoms%astruct%nat /), id = "indsym")
     symafm = f_malloc_ptr((/ nsym /), id = "symafm")
     allocate(symObj)
     symObj%xred = f_malloc0_ptr((/ 3, atoms%astruct%nat /), id = "xred")
  end if

  !Get the symmetry matrices in terms of reciprocal basis
  symrec = f_malloc((/ 3, 3, nsym /), id = "symrec")
  do isym = 1, nsym, 1
     call abi_mati3inv(sym(:,:,isym), symrec(:,:,isym))
  end do

!!$  if (present(pawang_sym)) then
!!$     if (present(comm_atom) .and. present(mpi_atmtab)) then
!!$        call symrhoij(pawrhoij,pawrhoij_all,choice,gprimd,indsym,ipert,natom,nsym,ntypat,&
!!$             &   option,pawang_sym,pawprtvol,pawtab,rprimd,symafm,symrec,typat,&
!!$             &   comm_atom=comm_atom,mpi_atmtab=mpi_atmtab,&
!!$             &   qphon=qphon)
!!$     else
!!$        call symrhoij(pawrhoij,pawrhoij_all,choice,gprimd,indsym,ipert,natom,nsym,ntypat,&
!!$             &   option,pawang_sym,pawprtvol,pawtab,rprimd,symafm,symrec,typat,&
!!$             &   qphon=qphon)
!!$     end if
!!$  else
!!$     if (size(paw%pawrhoij) < atoms%astruct%nat) then
!!$        call symrhoij(pawrhoij,pawrhoij_all,choice,gprimd,indsym,ipert,natom,nsym,ntypat,&
!!$             &   option,pawang,pawprtvol,pawtab,rprimd,symafm,symrec,typat,&
!!$             &   comm_atom=comm_atom,mpi_atmtab=mpi_atmtab,&
!!$             &   qphon=qphon)
!!$     else
        call symrhoij(pawrhoij_all, pawrhoij_all, 1, symObj%gprimd, indsym, ipert, &
             & atoms%astruct%nat, nsym, atoms%astruct%ntypes, 1, atoms%pawang, &
             & pawprtvol, atoms%pawtab, symObj%rprimd, symafm, symrec, &
             & atoms%astruct%iatype)
!!$     end if
!!$  end if

  call f_free(symrec)
  if (atoms%astruct%sym%nsym == 0) then
     call f_free_ptr(symObj%xred)
     deallocate(symObj)
     call f_free_ptr(indsym)
     call f_free_ptr(symafm)
     call f_free_ptr(sym)
  end if

  call pawrhoij_free_unpacked(pawrhoij_all)

  !In case of distribution over atomic sites, dispatch rhoij
  if (size(paw%pawrhoij) < atoms%astruct%nat) then
!!$     do iatom=1,size(paw%mpi_atmtab)
!!$        paw%pawrhoij(iatom)%rhoij_(:,:)=pawrhoij_all(paw%mpi_atmtab(iatom))%rhoij_(:,:)
!!$     end do
     call pawrhoij_free(pawrhoij_all)
     deallocate(pawrhoij_all)
  end if

  call f_timing(TCAT_LIBPAW, "OF")
END SUBROUTINE paw_compute_rhoij

subroutine paw_update_rho(paw, denspot, atoms)
  use module_defs, only: gp, dp
  use module_base, only: bigdft_mpi
  use module_types, only: paw_objects, DFT_local_fields, TCAT_LIBPAW, TCAT_PAW_RHOIJ
  use public_enums, only: ELECTRONIC_DENSITY
  use module_atoms
  use yaml_output
  use abi_interfaces_add_libpaw, only: abi_pawmknhat
  use dynamic_memory
  use time_profiling
  use wrapper_MPI
  use wrapper_linalg
  implicit none
  type(paw_objects), intent(inout) :: paw
  type(DFT_local_fields), intent(inout) :: denspot
  type(atoms_data), intent(in) :: atoms

  integer :: offset
  integer :: nfft, ngfft(18)
  real(gp) :: ucvol
  real(gp) :: compch_fft
  real(dp) :: rhodum(0,0,0)

  integer, parameter :: cplex = 1, pawprtvol = 0, usewvl = 1, pawxcdev = 1, ipert = 0, idir = 1
  real(gp), dimension(3), parameter :: qphon = (/ 0._gp, 0._gp, 0._gp /)
  real(gp), dimension(3,3) :: gprimd, rprimd ! Used only for phonons
  real(gp), dimension(3, atoms%astruct%nat) :: xred ! Idem

  call f_timing(TCAT_PAW_RHOIJ, "ON")

  nfft = denspot%dpbox%mesh%ndims(1) * denspot%dpbox%mesh%ndims(2) * denspot%dpbox%n3p
  ngfft(1:3) = denspot%dpbox%mesh%ndims
  ucvol = product(denspot%dpbox%mesh%ndims) *denspot%dpbox%mesh%volume_element

  !if (denspot%rhov_is /= ELECTRONIC_DENSITY) stop "rhov must be density here."
  if (bigdft_mpi%iproc == 0) then
     call yaml_newline()
  end if
  if (.not. associated(denspot%rhohat)) then
     denspot%rhohat = f_malloc_ptr((/ denspot%dpbox%mesh%ndims(1), denspot%dpbox%mesh%ndims(2), &
          & denspot%dpbox%n3p , denspot%dpbox%nrhodim /), id='denspot%rhohat')
  end if

  call f_timing(TCAT_LIBPAW, "ON")
  !Compute compensation charge density
  call abi_pawmknhat(compch_fft, cplex, 0, idir, ipert, 0, gprimd, &
       & size(paw%pawrhoij), atoms%astruct%nat, nfft, ngfft, 0, denspot%dpbox%nrhodim, &
       & atoms%astruct%ntypes, atoms%pawang, paw%fgrtab, rhodum, denspot%rhohat, &
       & paw%pawrhoij, paw%pawrhoij, atoms%pawtab, qphon, rprimd, ucvol, usewvl, &
       & xred, mpi_comm_wvl = bigdft_mpi%mpi_comm)
!!$  call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,"nhat0.cube",atoms,&
!!$       & symObj%xred,denspot%dpbox,denspot%dpbox%nrhodim,&
!!$       & denspot%rhohat(1, 1, denspot%dpbox%i3xcsh + 1, 1))
  call f_timing(TCAT_LIBPAW, "OF")
  if (bigdft_mpi%iproc == 0) then
     call yaml_map('Compensation charge on wavelet grid', compch_fft)
     call yaml_newline()
  end if

  if (denspot%dpbox%ndimpot > 0) then
     offset = denspot%dpbox%mesh%ndims(1) * denspot%dpbox%mesh%ndims(2) * denspot%dpbox%i3xcsh
     call axpy(denspot%dpbox%ndimpot, 1._dp, denspot%rhohat(1,1,1,1), 1, &
          & denspot%rhov(offset + 1), 1)
     if (denspot%dpbox%nrhodim == 2) then
        call axpy(denspot%dpbox%ndimpot, 1._dp, denspot%rhohat(1,1,1,2), 1, &
             & denspot%rhov(denspot%dpbox%ndimrhopot / denspot%dpbox%nrhodim + offset + 1), 1)
     end if
  end if

  call f_timing(TCAT_PAW_RHOIJ, "OF")
END SUBROUTINE paw_update_rho
