!> @file
!! H|psi> and orthonormalization (linear scaling version)
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine calculate_energy_and_gradient_linear(iproc, nproc, it, &
           ldiis, fnrmOldArr, alpha, trH, trHold, fnrm, fnrmMax, alpha_mean, alpha_max, &
           energy_increased, tmb, lhphiold, overlap_calculated, &
           energs, hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint, &
           hpsi_small, experimental_mode, correction_co_contra, hpsi_noprecond, &
           norder_taylor, method_updatekernel, precond_convol_workarrays, precond_workarrays, &
           kernel_old)
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, except_this_one => calculate_energy_and_gradient_linear
  use communications, only: transpose_localized
  use sparsematrix_base, only: matrices, matrices_null, deallocate_matrices, &
                               sparsematrix_malloc_ptr, assignment(=), SPARSE_FULL
  use sparsematrix_init, only: matrixindex_in_compressed
  use sparsematrix, only: transform_sparse_matrix
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, nproc, it, norder_taylor, method_updatekernel
  type(DFT_wavefunction), target, intent(inout):: tmb
  type(localizedDIISParameters), intent(inout) :: ldiis
  real(kind=8), dimension(tmb%orbs%norb), intent(inout) :: fnrmOldArr
  real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: alpha
  real(kind=8), intent(out):: trH, fnrm, fnrmMax, alpha_mean, alpha_max
  real(kind=8), intent(inout):: trHold
  logical,intent(out) :: energy_increased
  real(kind=8), dimension(tmb%npsidim_orbs), intent(inout):: lhphiold
  logical, intent(inout):: overlap_calculated
  type(energy_terms), intent(in) :: energs
  real(kind=8), dimension(:), pointer:: hpsit_c, hpsit_f
  integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint
  logical, intent(in) :: experimental_mode, correction_co_contra
  real(kind=8), dimension(tmb%npsidim_orbs), intent(out) :: hpsi_small
  real(kind=8), dimension(tmb%npsidim_orbs), optional,intent(out) :: hpsi_noprecond
  type(workarrays_quartic_convolutions),dimension(tmb%orbs%norbp),intent(inout) :: precond_convol_workarrays
  type(workarr_precond),dimension(tmb%orbs%norbp),intent(inout) :: precond_workarrays
  real(kind=8),dimension(tmb%linmat%smat(3)%nvctr),intent(inout) :: kernel_old

  ! Local variables
  integer :: iorb, iiorb, ilr, ncount, ierr, ist, ncnt, istat, iall, ii, jjorb, i
  integer :: lwork, info
  real(kind=8) :: ddot, tt, fnrmOvrlp_tot, fnrm_tot, fnrmold_tot
  character(len=*), parameter :: subname='calculate_energy_and_gradient_linear'
  real(kind=8), dimension(:), pointer :: hpsittmp_c, hpsittmp_f
  real(kind=8), dimension(:), allocatable :: fnrmOvrlpArr, fnrmArr
  real(kind=8), dimension(:), allocatable :: hpsi_conf, hpsi_tmp
  real(kind=8), dimension(:), pointer :: kernel_compr_tmp
  real(kind=8), dimension(:), allocatable :: prefac, hpsit_nococontra_c, hpsit_nococontra_f
  real(kind=8),dimension(3) :: reducearr
  real(wp), dimension(2) :: garray
  real(dp) :: gnrm,gnrm_zero,gnrmMax,gnrm_old ! for preconditional2, replace with fnrm eventually, but keep separate for now
  type(matrices) :: matrixm
  real(kind=8),dimension(:),allocatable :: kernel_old_small
  real(kind=8),dimension(:),pointer :: hpsi_old, hpsit_c_old, hpsit_f_old
  real(kind=8),dimension(:),allocatable :: kernel_tmp


 !!!!@ EXPERIMENTAL #################################
 !!!call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
 !!!      tmb%ham_descr%psi, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%lzd)
 !!!call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
 !!!     tmb%ham_descr%psit_c, hpsit_c, tmb%ham_descr%psit_f, hpsit_f, tmb%linmat%smat(2), tmb%linmat%ham_)
 !!!kernel_old_small = f_malloc(tmb%linmat%smat(2)%nvctr,id='kernel_old_small')
 !!!call transform_sparse_matrix(tmb%linmat%smat(2), tmb%linmat%smat(3), &
 !!!     kernel_old_small, kernel_old, 'large_to_small')
 !!!tt=ddot(tmb%linmat%smat(2)%nvctr, kernel_old_small(1), 1, tmb%linmat%ham_%matrix_compr(1), 1)
 !!!if (iproc==0) write(*,*) 'tt',tt
 !!!call f_free(kernel_old_small)
 !!!!@ END EXPERIMENTAL ##############################

  call f_routine(id='calculate_energy_and_gradient_linear')

  if (target_function==TARGET_FUNCTION_IS_HYBRID) then
      hpsi_conf = f_malloc(tmb%npsidim_orbs,id='hpsi_conf')
      call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
           tmb%orbs, tmb%hpsi, hpsi_conf)
      call timing(iproc,'eglincomms','ON')
      ist=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          tt=ddot(ncount, hpsi_conf(ist), 1, tmb%psi(ist), 1)
          call daxpy(ncount, -tt, tmb%psi(ist), 1, hpsi_conf(ist), 1)
          ist=ist+ncount
      end do
      call timing(iproc,'eglincomms','OF')
  end if

  ! by default no quick exit
  energy_increased=.false.


  hpsittmp_c = f_malloc_ptr(sum(tmb%ham_descr%collcom%nrecvcounts_c),id='hpsittmp_c')
  hpsittmp_f = f_malloc_ptr(7*sum(tmb%ham_descr%collcom%nrecvcounts_f),id='hpsittmp_f')

  hpsi_old= f_malloc_ptr(tmb%ham_descr%npsidim_orbs,id='hpsi_old')
  hpsit_c_old = f_malloc_ptr(tmb%ham_descr%collcom%ndimind_c,id='hpsit_c_old')
  hpsit_f_old = f_malloc_ptr(7*tmb%ham_descr%collcom%ndimind_f,id='hpsit_f_old')
  kernel_tmp = f_malloc(tmb%linmat%smat(3)%nvctr,id='kernel_tmp')

  if(target_function==TARGET_FUNCTION_IS_ENERGY .or. &
     target_function==TARGET_FUNCTION_IS_HYBRID) then

      if(sum(tmb%ham_descr%collcom%nrecvcounts_c)>0) &
          call vcopy(sum(tmb%ham_descr%collcom%nrecvcounts_c), hpsit_c(1), 1, hpsittmp_c(1), 1)
      if(sum(tmb%ham_descr%collcom%nrecvcounts_f)>0) &
          call vcopy(7*sum(tmb%ham_descr%collcom%nrecvcounts_f), hpsit_f(1), 1, hpsittmp_f(1), 1)

      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          ! Gradient with old kernel
          call vcopy(tmb%ham_descr%npsidim_orbs, tmb%hpsi(1), 1, hpsi_old(1), 1)
          call vcopy(tmb%ham_descr%collcom%ndimind_c, hpsit_c(1), 1, hpsit_c_old(1), 1)
          call vcopy(7*tmb%ham_descr%collcom%ndimind_f, hpsit_f(1), 1, hpsit_f_old(1), 1)
          call vcopy(tmb%linmat%smat(3)%nvctr, tmb%linmat%kernel_%matrix_compr(1), 1, kernel_tmp(1), 1)
          call get_hybrid_gradient(kernel_old, hpsi_old, hpsit_c_old, hpsit_f_old)
          call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
               tmb%ham_descr%psit_c, hpsit_c_old, tmb%ham_descr%psit_f, hpsit_f_old, tmb%linmat%smat(2), tmb%linmat%ham_)
          trH=0.d0
          call timing(iproc,'eglincomms','ON')
          do iorb=1,tmb%orbs%norb
             ii=matrixindex_in_compressed(tmb%linmat%smat(2),iorb,iorb)
             trH = trH + tmb%linmat%ham_%matrix_compr(ii)
          end do
          if (iproc==0) write(*,*) 'NEW: oldkernel trH',trH

          call vcopy(tmb%linmat%smat(3)%nvctr, kernel_tmp(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)


          ! Correct gradient
          !!call vcopy(tmb%ham_descr%collcom%ndimind_c, hpsittmp_c(1), 1, hpsit_c(1), 1)
          !!call vcopy(tmb%ham_descr%collcom%ndimind_f, hpsittmp_c(1), 1, hpsit_f(1), 1)
          call get_hybrid_gradient(tmb%linmat%kernel_%matrix_compr, tmb%hpsi, hpsit_c, hpsit_f)

          call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
               tmb%ham_descr%psit_c, hpsit_c, tmb%ham_descr%psit_f, hpsit_f, tmb%linmat%smat(2), tmb%linmat%ham_)
          trH=0.d0
          call timing(iproc,'eglincomms','ON')
          do iorb=1,tmb%orbs%norb
             ii=matrixindex_in_compressed(tmb%linmat%smat(2),iorb,iorb)
             trH = trH + tmb%linmat%ham_%matrix_compr(ii)
          end do
          if (iproc==0) write(*,*) 'NEW: newkernel trH',trH


          !!call timing(iproc,'eglincomms','ON')
          !!kernel_compr_tmp = f_malloc_ptr(tmb%linmat%smat(3)%nvctr,id='kernel_compr_tmp')
          !!call vcopy(tmb%linmat%smat(3)%nvctr, tmb%linmat%kernel_%matrix_compr(1), 1, kernel_compr_tmp(1), 1)
          !!    do ii=1,tmb%linmat%smat(3)%nvctr
          !!            iiorb = tmb%linmat%smat(3)%orb_from_index(1,ii)
          !!            jjorb = tmb%linmat%smat(3)%orb_from_index(2,ii)
          !!        if(iiorb==jjorb) then
          !!            tmb%linmat%kernel_%matrix_compr(ii)=0.d0
          !!        else
          !!            tmb%linmat%kernel_%matrix_compr(ii)=kernel_compr_tmp(ii)
          !!        end if
          !!    end do
          !!!end do

          !!ist=1
          !!do iorb=tmb%orbs%isorb+1,tmb%orbs%isorb+tmb%orbs%norbp
          !!    ilr=tmb%orbs%inwhichlocreg(iorb)
          !!    do ii=1,tmb%linmat%smat(3)%nvctr
          !!            iiorb = tmb%linmat%smat(3)%orb_from_index(1,ii)
          !!            jjorb = tmb%linmat%smat(3)%orb_from_index(2,ii)
          !!            if(iiorb==jjorb .and. iiorb==iorb) then
          !!                ncount=tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_f
          !!                call dscal(ncount, kernel_compr_tmp(ii), tmb%hpsi(ist), 1)
          !!                ist=ist+ncount
          !!            end if
          !!    end do
          !!end do
          !!call timing(iproc,'eglincomms','OF')
          !!call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
          !!     tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
          !!call build_linear_combination_transposed(tmb%ham_descr%collcom, &
          !!     tmb%linmat%smat(3), tmb%linmat%kernel_, hpsittmp_c, hpsittmp_f, .false., hpsit_c, hpsit_f, iproc)
          !!! copy correct kernel back
          !!call vcopy(tmb%linmat%smat(3)%nvctr, kernel_compr_tmp(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
          !!call f_free_ptr(kernel_compr_tmp)
      else
          call build_linear_combination_transposed(tmb%ham_descr%collcom, &
               tmb%linmat%smat(3), tmb%linmat%kernel_, hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f, iproc)
      end if
  end if

  hpsit_nococontra_c = f_malloc(tmb%ham_descr%collcom%ndimind_c,id='hpsit_nococontra_c')
  hpsit_nococontra_f = f_malloc(7*tmb%ham_descr%collcom%ndimind_f,id='hpsit_nococontra_f')

  call vcopy(tmb%ham_descr%collcom%ndimind_c, hpsit_c(1), 1, hpsit_nococontra_c(1), 1)
  call vcopy(7*tmb%ham_descr%collcom%ndimind_f, hpsit_f(1), 1, hpsit_nococontra_f(1), 1)

  if (correction_co_contra) then
      !@NEW correction for contra / covariant gradient

      if ((method_updatekernel/=UPDATE_BY_FOE .and. method_updatekernel/=UPDATE_BY_RENORMALIZATION) &
           .or. target_function/=TARGET_FUNCTION_IS_HYBRID) then
          call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
               tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
          tmb%can_use_transposed=.true.

          call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
               tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%smat(1), tmb%linmat%ovrlp_)
      end if
      call vcopy(tmb%ham_descr%collcom%ndimind_c, hpsit_c(1), 1, hpsittmp_c(1), 1)
      call vcopy(7*tmb%ham_descr%collcom%ndimind_f, hpsit_f(1), 1, hpsittmp_f(1), 1)
      matrixm = matrices_null()
      matrixm%matrix_compr = sparsematrix_malloc_ptr(tmb%linmat%smat(2), iaction=SPARSE_FULL, id='matrixm%matrix_compr')
      call transform_sparse_matrix(tmb%linmat%smat(1), tmb%linmat%smat(2), &
           tmb%linmat%ovrlp_%matrix_compr, matrixm%matrix_compr, 'small_to_large')
      call build_linear_combination_transposed(tmb%ham_descr%collcom, &
           tmb%linmat%smat(2), matrixm, hpsittmp_c, hpsittmp_f, .true., hpsit_c, hpsit_f, iproc)

      ! with old gradient ##############@######################
      call vcopy(tmb%ham_descr%collcom%ndimind_c, hpsit_c_old(1), 1, hpsittmp_c(1), 1)
      call vcopy(7*tmb%ham_descr%collcom%ndimind_f, hpsit_f_old(1), 1, hpsittmp_f(1), 1)
      call build_linear_combination_transposed(tmb%ham_descr%collcom, &
           tmb%linmat%smat(2), matrixm, hpsittmp_c, hpsittmp_f, .true., hpsit_c_old, hpsit_f_old, iproc)
      ! END with old gradient #################################

      call deallocate_matrices(matrixm)

      !@END NEW correction for contra / covariant gradient
  end if



  call f_free_ptr(hpsittmp_c)
  call f_free_ptr(hpsittmp_f)

  ! Calculate the overlap matrix if necessary
  if (.not.correction_co_contra) then
      call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
           tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
      tmb%can_use_transposed=.true.
      call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
           tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%smat(1), tmb%linmat%ovrlp_)
  end if




  if (.not.present(hpsi_noprecond)) stop 'hpsi_noprecond not present'
  call orthoconstraintNonorthogonal(iproc, nproc, tmb%ham_descr%lzd, &
       tmb%ham_descr%npsidim_orbs, tmb%ham_descr%npsidim_comp, &
       tmb%orbs, tmb%ham_descr%collcom, tmb%orthpar, correction_orthoconstraint, &
       tmb%linmat, tmb%ham_descr%psi, tmb%hpsi, &
       tmb%linmat%smat(2), tmb%linmat%ham_, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, &
       hpsit_c, hpsit_f, hpsit_c_old, hpsit_f_old, &
       hpsit_nococontra_c, hpsit_nococontra_f, tmb%ham_descr%can_use_transposed, &
       overlap_calculated, experimental_mode, norder_taylor, &
       tmb%npsidim_orbs, tmb%lzd, hpsi_noprecond)

  call f_free(hpsit_nococontra_c)
  call f_free(hpsit_nococontra_f)
  call f_free(kernel_tmp)
  call f_free_ptr(hpsi_old)
  call f_free_ptr(hpsit_c_old)
  call f_free_ptr(hpsit_f_old)


  call large_to_small_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
       tmb%orbs, tmb%hpsi, hpsi_small)


  ! Calculate trace (or band structure energy, resp.)
  trH=0.d0
  call timing(iproc,'eglincomms','ON')
  do iorb=1,tmb%orbs%norb
     ii=matrixindex_in_compressed(tmb%linmat%smat(2),iorb,iorb)
     trH = trH + tmb%linmat%ham_%matrix_compr(ii)
  end do
  if (iproc==0) write(*,*) 'trH',trH
  call timing(iproc,'eglincomms','OF')

  ! trH is now the total energy (name is misleading, correct this)
  ! Multiply by 2 because when minimizing trace we don't have kernel
  if(tmb%orbs%nspin==1 .and. target_function==TARGET_FUNCTION_IS_TRACE) trH=2.d0*trH
  trH=trH-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp


  ! Cycle if the trace increased (steepest descent only)
  if(.not. ldiis%switchSD .and. ldiis%isx==0) then
      if(trH > ldiis%trmin+1.d-12*abs(ldiis%trmin)) then !1.d-12 is here to tolerate some noise...
          !!if(iproc==0) write(*,'(1x,a,es18.10,a,es18.10)') &
          !!    'WARNING: the target function is larger than its minimal value reached so far:',trH,' > ', ldiis%trmin
          if (iproc==0) then
              call yaml_newline()
              call yaml_warning('target function larger than its minimal value reached so far, &
                  &D='//trim(yaml_toa(trH-ldiis%trmin,fmt='(1es10.3)')))!//'. &
                  !&Decrease step size and restart with previous TMBs')
          end if
          !if(iproc==0) write(*,'(1x,a)') 'Decrease step size and restart with previous TMBs'
          energy_increased=.true.
      end if
  end if

  fnrmOvrlpArr = f_malloc(tmb%orbs%norb,id='fnrmOvrlpArr')
  fnrmArr = f_malloc(tmb%orbs%norb,id='fnrmArr')

  ! Calculate the norm of the gradient (fnrmArr) and determine the angle between the current gradient and that
  ! of the previous iteration (fnrmOvrlpArr).
  call timing(iproc,'eglincomms','ON')
  ist=1
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      if(it>1) then
         fnrmOvrlpArr(iorb)=ddot(ncount, hpsi_small(ist), 1, lhphiold(ist), 1)
         !fnrmOldArr(iorb)=ddot(ncount, lhphiold(ist), 1, lhphiold(ist), 1)
      end if
      fnrmArr(iorb)=ddot(ncount, hpsi_small(ist), 1, hpsi_small(ist), 1)
      ist=ist+ncount
  end do



  ! Determine the gradient norm and its maximal component. In addition, adapt the
  ! step size for the steepest descent minimization (depending on the angle 
  ! between the current gradient and the one from the previous iteration).
  ! This is of course only necessary if we are using steepest descent and not DIIS.
  ! if newgradient is true, the angle criterion cannot be used and the choice whether to
  ! decrease or increase the step size is only based on the fact whether the trace decreased or increased.

  ! TEMPORARY #################################
  ! This is just for tests, can be improved
  fnrmOvrlp_tot=0.d0
  fnrm_tot=0.d0
  fnrmOld_tot=0.d0
  do iorb=1,tmb%orbs%norbp
      if (it>1) fnrmOvrlp_tot=fnrmOvrlp_tot+fnrmOvrlpArr(iorb)
      fnrm_tot=fnrm_tot+fnrmArr(iorb)
      if (it>1) fnrmOld_tot=fnrmOld_tot+fnrmOldArr(iorb)
  end do

  if (nproc > 1) then
     reducearr(1)=fnrm_tot
     reducearr(2)=fnrmOld_tot
     if (it>1) then
         reducearr(3)=fnrmOvrlp_tot
         call mpiallred(reducearr(1), 3, mpi_sum, bigdft_mpi%mpi_comm)
         fnrmOvrlp_tot=reducearr(3)
     else
         call mpiallred(reducearr(1), 2, mpi_sum, bigdft_mpi%mpi_comm)
     end if
     fnrm_tot=reducearr(1)
     fnrmOld_tot=reducearr(2)
  end if

  ! ###########################################

  fnrm=0.d0
  fnrmMax=0.d0
  do iorb=1,tmb%orbs%norbp
      fnrm=fnrm+fnrmArr(iorb)
      if(fnrmArr(iorb)>fnrmMax) fnrmMax=fnrmArr(iorb)
      if(it>1 .and. ldiis%isx==0 .and. .not.ldiis%switchSD) then
      ! Adapt step size for the steepest descent minimization.
          if (experimental_mode) then
              if (iproc==0 .and. iorb==1) then
                  !write(*,*) 'WARNING: USING SAME STEP SIZE'
                  call yaml_warning('Using same step size for all TMBs')
              end if
              tt=fnrmOvrlp_tot/sqrt(fnrm_tot*fnrmOld_tot)
          else
              tt=fnrmOvrlpArr(iorb)/sqrt(fnrmArr(iorb)*fnrmOldArr(iorb))
          end if
          ! apply thresholds so that alpha never goes below around 1.d-2 and above around 2
          if(tt>.6d0 .and. trH<trHold .and. alpha(iorb)<1.8d0) then
              alpha(iorb)=alpha(iorb)*1.1d0
          else if (alpha(iorb)>1.7d-3) then
              alpha(iorb)=alpha(iorb)*.6d0
          end if
          !!alpha(iorb)=min(alpha(iorb),1.5d0)
      end if
  end do
  
  if (nproc > 1) then
     call mpiallred(fnrm, 1, mpi_sum, bigdft_mpi%mpi_comm)
     call mpiallred(fnrmMax, 1, mpi_max, bigdft_mpi%mpi_comm)
  end if

  fnrm=sqrt(fnrm/dble(tmb%orbs%norb))
  fnrmMax=sqrt(fnrmMax)

  call vcopy(tmb%orbs%norb, fnrmArr(1), 1, fnrmOldArr(1), 1)

  call f_free(fnrmOvrlpArr)
  call f_free(fnrmArr)

  ! Determine the mean step size for steepest descent iterations.
  tt=sum(alpha)
  alpha_max=maxval(alpha)
  if (nproc > 1) then
     call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm)
     call mpiallred(alpha_max, 1, mpi_max, bigdft_mpi%mpi_comm)
  end if
  alpha_mean=tt/dble(tmb%orbs%norb)

  ! Copy the gradient (will be used in the next iteration to adapt the step size).
  call vcopy(tmb%npsidim_orbs, hpsi_small(1), 1, lhphiold(1), 1)
  call timing(iproc,'eglincomms','OF')

  ! if energy has increased or we only wanted to calculate the energy, not gradient, we can return here
  ! rather than calculating the preconditioning for nothing
  if ((energy_increased) .and. target_function/=TARGET_FUNCTION_IS_HYBRID) return



  if(target_function==TARGET_FUNCTION_IS_HYBRID) then
     hpsi_tmp = f_malloc(tmb%npsidim_orbs,id='hpsi_tmp')
     ist=1
     do iorb=1,tmb%orbs%norbp
        iiorb=tmb%orbs%isorb+iorb
        ilr = tmb%orbs%inWhichLocreg(iiorb)
        ncnt=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f

        tt=ddot(ncnt, hpsi_conf(ist), 1, hpsi_small(ist), 1)
        tt=tt/ddot(ncnt, hpsi_conf(ist), 1, hpsi_conf(ist), 1)
        do i=ist,ist+ncnt-1
           hpsi_tmp(i)=tt*hpsi_conf(i)
        end do
        call daxpy(ncnt, -tt, hpsi_conf(ist), 1, hpsi_small(ist), 1)

        ist=ist+ncnt
     end do

     call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
          tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
          nit_precond,tmb%npsidim_orbs,hpsi_tmp,tmb%confdatarr,gnrm,gnrm_zero, &
          precond_convol_workarrays, precond_workarrays)

     ! temporarily turn confining potential off...
     prefac = f_malloc(tmb%orbs%norbp,id='prefac')
     prefac(:)=tmb%confdatarr(:)%prefac
     tmb%confdatarr(:)%prefac=0.0d0
     call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
          tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
          nit_precond,tmb%npsidim_orbs,hpsi_small,tmb%confdatarr,gnrm,gnrm_zero, & ! prefac should be zero
          precond_convol_workarrays, precond_workarrays)
     call daxpy(tmb%npsidim_orbs, 1.d0, hpsi_tmp(1), 1, hpsi_small(1), 1)
     ! ...revert back to correct value
     tmb%confdatarr(:)%prefac=prefac

     call f_free(prefac)
     call f_free(hpsi_conf)
     call f_free(hpsi_tmp)
  else
     call preconditionall2(iproc,nproc,tmb%orbs,tmb%Lzd,&
          tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3),&
          nit_precond,tmb%npsidim_orbs,hpsi_small,tmb%confdatarr,gnrm,gnrm_zero,&
          precond_convol_workarrays, precond_workarrays)
  end if

  if (iproc==0) then
      call yaml_map('Preconditioning',.true.)
  end if

  call f_release_routine()


  contains


    subroutine get_hybrid_gradient(kernel, hpsi, hpsit_c, hpsit_f)
      real(kind=8),dimension(tmb%linmat%smat(3)%nvctr),intent(inout) :: kernel
      real(kind=8),dimension(tmb%ham_descr%npsidim_orbs),intent(inout) :: hpsi
      real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_c),intent(inout) :: hpsit_c
      real(kind=8),dimension(tmb%ham_descr%collcom%ndimind_f),intent(inout) :: hpsit_f

      if(sum(tmb%ham_descr%collcom%nrecvcounts_c)>0) &
          call vcopy(sum(tmb%ham_descr%collcom%nrecvcounts_c), hpsit_c(1), 1, hpsittmp_c(1), 1)
      if(sum(tmb%ham_descr%collcom%nrecvcounts_f)>0) &
          call vcopy(7*sum(tmb%ham_descr%collcom%nrecvcounts_f), hpsit_f(1), 1, hpsittmp_f(1), 1)

      call timing(iproc,'eglincomms','ON')
      kernel_compr_tmp = f_malloc_ptr(tmb%linmat%smat(3)%nvctr,id='kernel_compr_tmp')
      call vcopy(tmb%linmat%smat(3)%nvctr, kernel(1), 1, kernel_compr_tmp(1), 1)
          do ii=1,tmb%linmat%smat(3)%nvctr
                  iiorb = tmb%linmat%smat(3)%orb_from_index(1,ii)
                  jjorb = tmb%linmat%smat(3)%orb_from_index(2,ii)
              if(iiorb==jjorb) then
                  tmb%linmat%kernel_%matrix_compr(ii)=0.d0
              else
                  tmb%linmat%kernel_%matrix_compr(ii)=kernel_compr_tmp(ii)
              end if
          end do
      !end do

      ist=1
      do iorb=tmb%orbs%isorb+1,tmb%orbs%isorb+tmb%orbs%norbp
          ilr=tmb%orbs%inwhichlocreg(iorb)
          do ii=1,tmb%linmat%smat(3)%nvctr
                  iiorb = tmb%linmat%smat(3)%orb_from_index(1,ii)
                  jjorb = tmb%linmat%smat(3)%orb_from_index(2,ii)
                  if(iiorb==jjorb .and. iiorb==iorb) then
                      ncount=tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%ham_descr%lzd%llr(ilr)%wfd%nvctr_f
                      call dscal(ncount, kernel_compr_tmp(ii), hpsi(ist), 1)
                      ist=ist+ncount
                  end if
          end do
      end do
      call timing(iproc,'eglincomms','OF')
      call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
           hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
      call build_linear_combination_transposed(tmb%ham_descr%collcom, &
           tmb%linmat%smat(3), tmb%linmat%kernel_, hpsittmp_c, hpsittmp_f, .false., hpsit_c, hpsit_f, iproc)
      ! copy correct kernel back
      call vcopy(tmb%linmat%smat(3)%nvctr, kernel_compr_tmp(1), 1, kernel(1), 1)
      call f_free_ptr(kernel_compr_tmp)

    end subroutine get_hybrid_gradient

end subroutine calculate_energy_and_gradient_linear


subroutine calculate_residue_ks(iproc, nproc, num_extra, ksorbs, tmb, hpsit_c, hpsit_f)
  use module_base
  use module_types
  use module_interfaces, except_this_one => calculate_residue_ks
  use sparsematrix_base, only: sparse_matrix, sparse_matrix_null, deallocate_sparse_matrix, &
                               matrices_null, allocate_matrices, deallocate_matrices
  use sparsematrix, only: uncompress_matrix
  implicit none

  ! Calling arguments
  integer, intent(in) :: iproc, nproc, num_extra
  type(dft_wavefunction), intent(inout) :: tmb
  type(orbitals_data), intent(in) :: ksorbs
  real(kind=8),dimension(:),pointer :: hpsit_c, hpsit_f

  integer :: iorb, istat, iall, ierr
  real(kind=8) :: ksres_sum
  real(kind=8), dimension(:), allocatable :: ksres
  real(kind=8), dimension(:,:), allocatable :: coeff_tmp, grad_coeff
  type(sparse_matrix) :: grad_ovrlp
  type(matrices) :: grad_ovrlp_
  character(len=256) :: subname='calculate_residue_ks'


  call f_routine(id='calculate_residue_ks')

  ! want to calculate the residue of the KS states here, not just the tmbs
  ! for now just occupied, eventually would want occupied+num_extra
  ! probably better to calculate c_i^a|g_a> first but copying and pasting for now (INEFFICIENT but just testing)
  !!if(associated(hpsit_c)) then
  !!    iall=-product(shape(hpsit_c))*kind(hpsit_c)
  !!    deallocate(hpsit_c, stat=istat)
  !!    call memocc(istat, iall, 'hpsit_c', subname)
  !!end if
  !!if(associated(hpsit_f)) then
  !!    iall=-product(shape(hpsit_f))*kind(hpsit_f)
  !!    deallocate(hpsit_f, stat=istat)
  !!    call memocc(istat, iall, 'hpsit_f', subname)
  !!end if
  !!allocate(hpsit_c(sum(collcom%nrecvcounts_c)), stat=istat)
  !!call memocc(istat, hpsit_c, 'hpsit_c', subname)
  !!allocate(hpsit_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
  !!call memocc(istat, hpsit_f, 'hpsit_f', subname)

  ! should already be done in orthoconstraintnonorthogonal so can use directly
  !call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
  !     tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
  !!can_use_transposed=.true.

  !call nullify_sparse_matrix(grad_ovrlp)
  grad_ovrlp=sparse_matrix_null()
  call sparse_copy_pattern(tmb%linmat%smat(2), grad_ovrlp, iproc, subname)
  grad_ovrlp%matrix_compr=f_malloc_ptr(grad_ovrlp%nvctr,id='grad_ovrlp%matrix_compr')
  grad_ovrlp_ = matrices_null()
  call allocate_matrices(tmb%linmat%smat(2), allocate_full=.false., &
       matname='grad_ovrlp_', mat=grad_ovrlp_)

  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, hpsit_c, hpsit_c, &
       hpsit_f, hpsit_f, tmb%linmat%smat(2), grad_ovrlp_)
  ! This can then be deleted if the transition to the new type has been completed.
  grad_ovrlp%matrix_compr=grad_ovrlp_%matrix_compr

  call deallocate_matrices(grad_ovrlp_)

  grad_coeff = f_malloc((/ tmb%orbs%norb, tmb%orbs%norb /),id='grad_coeff')
  coeff_tmp = f_malloc((/ tmb%orbs%norbp, max(tmb%orbs%norb, 1) /),id='coeff_tmp')

  grad_ovrlp%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='grad_ovrlp%matrix')
  call uncompress_matrix(iproc,grad_ovrlp)

  ! can change this so only go up to ksorbs%norb...
  if (tmb%orbs%norbp>0) then
     call dgemm('n', 'n', tmb%orbs%norbp, tmb%orbs%norb, tmb%orbs%norb, 1.d0, grad_ovrlp%matrix(tmb%orbs%isorb+1,1), &
          tmb%orbs%norb, tmb%coeff(1,1), tmb%orbs%norb, 0.d0, coeff_tmp, tmb%orbs%norbp)
     call dgemm('t', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norbp, 1.d0, tmb%coeff(tmb%orbs%isorb+1,1), &
          tmb%orbs%norb, coeff_tmp, tmb%orbs%norbp, 0.d0, grad_coeff, tmb%orbs%norb)
  else
     call to_zero(tmb%orbs%norb**2, grad_coeff(1,1))
  end if

  call f_free(coeff_tmp)

  if (nproc>1) then
      call mpiallred(grad_coeff(1,1), tmb%orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  ! now calculate sqrt(<g_i|g_i>) and mean value
  ksres = f_malloc(ksorbs%norb+num_extra,id='ksres')
  
  ksres_sum=0.0d0
  do iorb=1,ksorbs%norb+num_extra
    ksres(iorb)=dsqrt(grad_coeff(iorb,iorb))
    !ksres_sum=ksres_sum+ksres(iorb)
    ksres_sum=ksres_sum+grad_coeff(iorb,iorb)
    if (iproc==0) write(*,*) 'KS residue',iorb,ksres(iorb)!,tmb%orbs%occup(iorb)
  end do
  if (iproc==0) write(*,*) 'Average KS residue',sqrt(ksres_sum/real(ksorbs%norb+num_extra,gp))


  !call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
  !     tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%ham)
  ! calculate Tr[Kg]  (not recalculating kernel as don't have the correct occs here)
  !call calculate_kernel_and_energy(iproc,nproc,denskern,grad_coeff,ksres_sum,tmb%coeff,orbs,tmb%orbs,.true.)
  grad_ovrlp_ = matrices_null()
  call allocate_matrices(tmb%linmat%smat(2), allocate_full=.false., matname='grad_ovrlp_', mat=grad_ovrlp_)
  grad_ovrlp_%matrix_compr=grad_ovrlp%matrix_compr
  call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%smat(3),grad_ovrlp,&
       tmb%linmat%kernel_, grad_ovrlp_, &
       ksres_sum,tmb%coeff,tmb%orbs,tmb%orbs,.false.)
  call deallocate_matrices(grad_ovrlp_)
  if (iproc==0) write(*,*) 'KS residue from trace',dsqrt(ksres_sum)/real(tmb%orbs%norb,gp) ! should update normalization as would only be occ here not extra?

  call deallocate_sparse_matrix(grad_ovrlp, subname)

  call f_free(grad_coeff)
  call f_free(ksres)

  call f_release_routine()

end subroutine calculate_residue_ks


subroutine hpsitopsi_linear(iproc, nproc, it, ldiis, tmb,  &
           lphiold, alpha, trH, alpha_mean, alpha_max, alphaDIIS, hpsi_small, ortho, psidiff, &
           experimental_mode, order_taylor,trH_ref, kernel_best, complete_reset)
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, except_this_one => hpsitopsi_linear
  use communications, only: transpose_localized, untranspose_localized
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, it, order_taylor
  type(localizedDIISParameters), intent(inout) :: ldiis
  type(DFT_wavefunction), target,intent(inout) :: tmb
  real(kind=8), dimension(tmb%npsidim_orbs), intent(inout) :: lphiold
  real(kind=8), intent(in) :: trH, alpha_mean, alpha_max
  real(kind=8), dimension(tmb%orbs%norbp), intent(inout) :: alpha, alphaDIIS
  real(kind=8), dimension(tmb%npsidim_orbs), intent(inout) :: hpsi_small
  real(kind=8), dimension(tmb%npsidim_orbs), optional,intent(out) :: psidiff
  logical, intent(in) :: ortho, experimental_mode
  real(kind=8),intent(out) :: trH_ref
  real(kind=8),dimension(tmb%linmat%smat(3)%nvctr),intent(out) :: kernel_best
  logical,intent(out) :: complete_reset

  ! Local variables
  integer :: i, iorb, ilr, ist, iiorb, ncount
  character(len=*), parameter :: subname='hpsitopsi_linear'
  real(kind=8), dimension(:), allocatable :: norm
  real(kind=8) :: ddot, dnrm2, tt

  call f_routine(id='hpsitopsi_linear')

  call DIISorSD(iproc, it, trH, tmb, ldiis, alpha, alphaDIIS, lphiold, trH_ref, kernel_best, complete_reset)
  if(iproc==0) then
      call yaml_newline()
      call yaml_mapping_open('Optimization',flow=.true.)
      if(ldiis%isx>0) then
          call yaml_map('algorithm','DIIS')
          call yaml_map('history length',ldiis%isx)
          call yaml_map('consecutive failures',ldiis%icountDIISFailureCons)
          call yaml_map('total failures',ldiis%icountDIISFailureTot)
      else
          call yaml_map('algorithm','SD')
          call yaml_map('mean alpha',alpha_mean,fmt='(es9.3)')
          call yaml_map('max alpha',alpha_max,fmt='(es9.3)')
          call yaml_map('consecutive successes',ldiis%icountSDSatur)
      end if
      call yaml_mapping_close()
      call yaml_newline()
  end if

  ! Improve the orbitals, depending on the choice made above.
  if (present(psidiff)) call vcopy(tmb%npsidim_orbs, tmb%psi(1), 1, psidiff(1), 1)
  if(.not.ldiis%switchSD) then
      call improveOrbitals(iproc, nproc, tmb, ldiis, alpha, hpsi_small, experimental_mode)
  else
      if (iproc==0) then
          call yaml_warning('no improvement of the orbitals, recalculate gradient')
          call yaml_newline()
      end if
  end if

  ! The transposed quantities can now not be used any more...
  if(tmb%can_use_transposed) then
      tmb%can_use_transposed=.false.
  end if


  if (.not.ortho .and. iproc==0) then
      call yaml_map('Orthogonalization',.false.)
  end if

  if(.not.ldiis%switchSD.and.ortho) then
      if (present(psidiff)) then
          do i=1,tmb%npsidim_orbs
              psidiff(i)=tmb%psi(i)-psidiff(i)
          end do 
      end if

      call orthonormalizeLocalized(iproc, nproc, order_taylor, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, &
           tmb%linmat%smat(1), tmb%linmat%smat(3), tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, &
           tmb%can_use_transposed, tmb%foe_obj)
      if (iproc == 0) then
          call yaml_map('Orthogonalization',.true.)
      end if
  else if (experimental_mode .or. .not.ldiis%switchSD) then
      ist=1
      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          ilr=tmb%orbs%inwhichlocreg(iiorb)
          ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          tt=dnrm2(ncount, tmb%psi(ist), 1)
          tt=1.d0/tt
          call dscal(ncount, tt, tmb%psi(ist), 1)
          ist=ist+ncount
      end do
      if (iproc == 0) then
          call yaml_map('Normalization',.true.)
          !call yaml_map('Normalization',.false.)
      end if
      if (present(psidiff)) then
          do i=1,tmb%npsidim_orbs
              psidiff(i)=tmb%psi(i)-psidiff(i)
          end do 
      end if
  end if

  ! Emit that new wavefunctions are ready.
  if (tmb%c_obj /= 0) then
     call kswfn_emit_psi(tmb, it, 0, iproc, nproc)
  end if

  call f_release_routine

end subroutine hpsitopsi_linear
