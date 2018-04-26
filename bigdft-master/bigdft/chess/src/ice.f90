!> @file
!!   File containing the main matrix power expansion routine
!! @author
!!   Copyright (C) 2016 CheSS developers
!!
!!   This file is part of CheSS.
!!   
!!   CheSS is free software: you can redistribute it and/or modify
!!   it under the terms of the GNU Lesser General Public License as published by
!!   the Free Software Foundation, either version 3 of the License, or
!!   (at your option) any later version.
!!   
!!   CheSS is distributed in the hope that it will be useful,
!!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU Lesser General Public License for more details.
!!   
!!   You should have received a copy of the GNU Lesser General Public License
!!   along with CheSS.  If not, see <http://www.gnu.org/licenses/>.


!> Module Inverse Chebyshev Expansion
module ice
  use sparsematrix_base
  use dynamic_memory
  use yaml_output
  use dictionaries, only: f_err_throw
  use time_profiling
  implicit none

  private

  !> Public routines
  !public :: inverse_chebyshev_expansion
  public :: inverse_chebyshev_expansion_new
  public :: calculate_fermi_function_entropy

  contains

!> Chebyshev expansion of the inverse overlap (Inverse Chebyshev Expansion)
!!    subroutine inverse_chebyshev_expansion(iproc, nproc, norder_polynomial, &
!!               ovrlp_smat, inv_ovrlp_smat, ncalc, ex, ovrlp_mat, inv_ovrlp, &
!!               verbosity, npl_auto)
!!      use module_base
!!      use yaml_output
!!      use sparsematrix_base, only: sparsematrix_malloc_ptr, sparsematrix_malloc, &
!!                                   sparsematrix_malloc0_ptr, assignment(=), &
!!                                   SPARSE_TASKGROUP, SPARSE_MATMUL_SMALL, &
!!                                   matrices, sparse_matrix, matrices_null, deallocate_matrices
!!      use sparsematrix_init, only: matrixindex_in_compressed
!!      use sparsematrix, only: compress_matrix, uncompress_matrix, compress_matrix_distributed_wrapper, &
!!                              transform_sparsity_pattern
!!      use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
!!                          foe_data_set_logical, foe_data_get_logical
!!      use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
!!                             fermilevel_get_real, fermilevel_get_logical
!!      use chebyshev, only: chebyshev_clean, chebyshev_fast
!!      use foe_common, only: scale_and_shift_matrix, &
!!                            evnoise, check_eigenvalue_spectrum_new, get_chebyshev_expansion_coefficients, &
!!                            get_chebyshev_polynomials
!!      use module_func
!!      implicit none
!!
!!      ! Calling arguments
!!      integer, intent(in) :: iproc, nproc, norder_polynomial, ncalc
!!      type(sparse_matrix), intent(in) :: ovrlp_smat, inv_ovrlp_smat
!!      real(kind=mp),dimension(ncalc),intent(in) :: ex
!!      type(matrices), intent(in) :: ovrlp_mat
!!      type(matrices), dimension(ncalc), intent(inout) :: inv_ovrlp
!!      integer, intent(in),optional :: verbosity
!!      logical, intent(in),optional :: npl_auto
!!
!!      ! Local variables
!!      integer :: npl, jorb, it, ii, iseg, verbosity_
!!      integer :: isegstart, isegend, iismall, nsize_polynomial
!!      integer :: iismall_ovrlp, iismall_ham, npl_boundaries, i, ipl
!!      integer, parameter :: nplx=50000
!!      real(kind=mp),dimension(:,:),pointer :: chebyshev_polynomials
!!      real(kind=mp),dimension(:,:,:),pointer :: inv_ovrlp_matrixp
!!      real(kind=mp),dimension(:,:,:),allocatable :: penalty_ev
!!      real(kind=mp),dimension(:,:,:),pointer :: cc
!!      real(kind=mp) :: anoise, scale_factor, shift_value
!!      real(kind=mp) :: evlow_old, evhigh_old, tt
!!      real(kind=mp) :: x_max_error_fake, max_error_fake, mean_error_fake
!!      real(kind=mp) :: tt_ovrlp, tt_ham, eval_multiplicator, eval_multiplicator_total
!!      logical :: restart, calculate_SHS
!!      logical, dimension(2) :: emergency_stop
!!      real(kind=mp),dimension(2) :: allredarr
!!      real(kind=mp),dimension(:),allocatable :: hamscal_compr
!!      logical, dimension(2) :: eval_bounds_ok
!!      integer, dimension(2) :: irowcol
!!      integer :: irow, icol, iflag, ispin, isshift, ilshift, ilshift2
!!      logical :: overlap_calculated, evbounds_shrinked, degree_sufficient, reached_limit, npl_auto_
!!      integer, parameter :: NPL_MIN=5
!!      real(kind=mp),parameter :: DEGREE_MULTIPLICATOR_MAX=20.d0
!!      real(kind=mp) :: degree_multiplicator
!!      integer, parameter :: SPARSE=1
!!      integer, parameter :: DENSE=2
!!      integer, parameter :: imode=SPARSE
!!      type(foe_data) :: foe_obj
!!      real(kind=mp),dimension(:),allocatable :: eval, work, x_max_error, max_error, mean_error
!!      real(kind=mp),dimension(:,:),allocatable :: tempmat
!!      integer :: lwork, info, j, icalc, iline, icolumn
!!      real(kind=mp),dimension(:,:),allocatable :: inv_ovrlp_matrixp_new
!!      real(kind=mp),dimension(:,:),allocatable :: penalty_ev_new
!!      real(kind=mp),dimension(:,:),allocatable :: inv_ovrlp_matrixp_small_new
!!      type(matrices) :: ovrlp_scaled
!!      character(len=3), parameter :: old='old'
!!      character(len=3), parameter :: new='new'
!!      character(len=3) :: mode=old
!!
!!      !!real(kind=mp),dimension(ovrlp_smat%nfvctr,ovrlp_smat%nfvctr) :: overlap
!!      !!real(kind=mp),dimension(ovrlp_smat%nfvctr) :: eval
!!      !!integer, parameter :: lwork=100000
!!      !!real(kind=mp),dimension(lwork) :: work
!!      !!integer :: info
!!
!!      call f_routine(id='inverse_chebyshev_expansion')
!!
!!      if (present(npl_auto)) then
!!          npl_auto_ = npl_auto
!!      else
!!          npl_auto_ = .false.
!!      end if
!!
!!      if (present(verbosity)) then
!!          verbosity_ = verbosity
!!      else
!!          verbosity_ = 1
!!      end if
!!
!!
!!      penalty_ev_new = f_malloc((/inv_ovrlp_smat%smmm%nvctrp,2/),id='penalty_ev_new')
!!      inv_ovrlp_matrixp_new = f_malloc((/max(inv_ovrlp_smat%smmm%nvctrp,1),ncalc/),id='inv_ovrlp_matrixp_new')
!!      inv_ovrlp_matrixp_small_new = f_malloc((/max(inv_ovrlp_smat%smmm%nvctrp_mm,1),ncalc/),id='inv_ovrlp_matrixp_small_new')
!!
!!
!!    !@ JUST FOR THE MOMENT.... ########################
!!         foe_obj%ef = f_malloc0_ptr(ovrlp_smat%nspin,id='(foe_obj%ef)')
!!         foe_obj%evlow = f_malloc0_ptr(ovrlp_smat%nspin,id='foe_obj%evlow')
!!         foe_obj%evhigh = f_malloc0_ptr(ovrlp_smat%nspin,id='foe_obj%evhigh')
!!         foe_obj%bisection_shift = f_malloc0_ptr(ovrlp_smat%nspin,id='foe_obj%bisection_shift')
!!         foe_obj%charge = f_malloc0_ptr(ovrlp_smat%nspin,id='foe_obj%charge')
!!         do ispin=1,ovrlp_smat%nspin
!!             call foe_data_set_real(foe_obj,"ef",0.d0,ispin)
!!             call foe_data_set_real(foe_obj,"evlow",0.3d0,ispin)
!!             call foe_data_set_real(foe_obj,"evhigh",2.2d0,ispin)
!!             call foe_data_set_real(foe_obj,"bisection_shift",1.d-1,ispin)
!!             call foe_data_set_real(foe_obj,"charge",0.d0,ispin)
!!         end do
!!
!!         call foe_data_set_real(foe_obj,"fscale",1.d-1)
!!         call foe_data_set_real(foe_obj,"ef_interpol_det",0.d0)
!!         call foe_data_set_real(foe_obj,"ef_interpol_chargediff",0.d0)
!!         call foe_data_set_int(foe_obj,"evbounds_isatur",0)
!!         call foe_data_set_int(foe_obj,"evboundsshrink_isatur",0)
!!         call foe_data_set_int(foe_obj,"evbounds_nsatur",10)
!!         call foe_data_set_int(foe_obj,"evboundsshrink_nsatur",10)
!!         call foe_data_set_real(foe_obj,"fscale_lowerbound",1.d-2)
!!         call foe_data_set_real(foe_obj,"fscale_upperbound",0.d0)
!!    !@ ################################################
!!
!!
!!      evbounds_shrinked = .false.
!!
!!      !@ TEMPORARY: eigenvalues of  the overlap matrix ###################
!!      !call get_minmax_eigenvalues(iproc, ovrlp_smat, ovrlp_mat)
!!
!!      !!tempmat = f_malloc0((/ovrlp_smat%nfvctr,ovrlp_smat%nfvctr/),id='tempmat')
!!      !!do iseg=1,ovrlp_smat%nseg
!!      !!    ii=ovrlp_smat%keyv(iseg)
!!      !!    do i=ovrlp_smat%keyg(1,1,iseg),ovrlp_smat%keyg(2,1,iseg)
!!      !!        tempmat(i,ovrlp_smat%keyg(1,2,iseg)) = ovrlp_mat%matrix_compr(ii)
!!      !!        ii = ii + 1
!!      !!    end do
!!      !!end do
!!      !!!!if (iproc==0) then
!!      !!!!    do i=1,ovrlp_smat%nfvctr
!!      !!!!        do j=1,ovrlp_smat%nfvctr
!!      !!!!            write(*,'(a,2i6,es17.8)') 'i,j,val',i,j,tempmat(j,i)
!!      !!!!        end do
!!      !!!!    end do
!!      !!!!end if
!!      !!eval = f_malloc(ovrlp_smat%nfvctr,id='eval')
!!      !!lwork=100*ovrlp_smat%nfvctr
!!      !!work = f_malloc(lwork,id='work')
!!      !!call dsyev('n','l', ovrlp_smat%nfvctr, tempmat, ovrlp_smat%nfvctr, eval, work, lwork, info)
!!      !!!if (iproc==0) write(*,*) 'eval',eval
!!      !!if (iproc==0) call yaml_map('eval max/min',(/eval(1),eval(ovrlp_smat%nfvctr)/),fmt='(es16.6)')
!!
!!      !!call f_free(tempmat)
!!      !!call f_free(eval)
!!      !!call f_free(work)
!!      !@ END TEMPORARY: eigenvalues of  the overlap matrix ###############
!!
!!
!!      call timing(iproc, 'FOE_auxiliary ', 'ON')
!!
!!
!!
!!      !!penalty_ev = f_malloc((/inv_ovrlp_smat%nfvctr,inv_ovrlp_smat%smmm%nfvctrp,2/),id='penalty_ev')
!!
!!
!!      if (npl_auto_) then
!!          ovrlp_scaled = matrices_null()
!!          ovrlp_scaled%matrix_compr = sparsematrix_malloc_ptr(ovrlp_smat, &
!!              iaction=SPARSE_TASKGROUP, id='ovrlp_scaled%matrix_compr')
!!          call f_memcpy(src=ovrlp_mat%matrix_compr,dest=ovrlp_scaled%matrix_compr)
!!      end if
!!      hamscal_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_TASKGROUP, id='hamscal_compr')
!!
!!
!!      ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
!!      nsize_polynomial = inv_ovrlp_smat%smmm%nvctrp_mm
!!
!!
!!      ! Fake allocation, will be modified later
!!      chebyshev_polynomials = f_malloc_ptr((/nsize_polynomial,1/),id='chebyshev_polynomials')
!!
!!
!!      !inv_ovrlp_matrixp = sparsematrix_malloc0_ptr(inv_ovrlp_smat, &
!!      !                         iaction=DENSE_MATMUL, id='inv_ovrlp_matrixp')
!!      !!inv_ovrlp_matrixp = f_malloc_ptr((/inv_ovrlp_smat%nfvctr,inv_ovrlp_smat%smmm%nfvctrp,ncalc/),&
!!      !!                                  id='inv_ovrlp_matrixp')
!!
!!
!!          spin_loop: do ispin=1,ovrlp_smat%nspin
!!
!!              degree_multiplicator = real(norder_polynomial,kind=mp)/ &
!!                                     (foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin))
!!              degree_multiplicator = min(degree_multiplicator,DEGREE_MULTIPLICATOR_MAX)
!!
!!              isshift=(ispin-1)*ovrlp_smat%nvctr
!!              ilshift=(ispin-1)*inv_ovrlp_smat%nvctr
!!              ilshift2=(ispin-1)*inv_ovrlp_smat%nvctr
!!
!!              evlow_old=1.d100
!!              evhigh_old=-1.d100
!!
!!              eval_multiplicator = 1.d0
!!              eval_multiplicator_total = 1.d0
!!
!!
!!                  !!calculate_SHS=.true.
!!
!!              !if (inv_ovrlp_smat%smmm%nfvctrp>0) then !LG: this conditional seems decorrelated
!!              !call f_zero(inv_ovrlp_smat%nfvctr*inv_ovrlp_smat%smmm%nfvctrp*ncalc, inv_ovrlp_matrixp(1,1,1))
!!              !end if
!!              !!    call f_zero(inv_ovrlp_matrixp)
!!
!!
!!                  it=0
!!                  eval_bounds_ok=.false.
!!                  !!bisection_bounds_ok=.false.
!!                  main_loop: do
!!
!!                      it=it+1
!!
!!
!!                      ! Scale the Hamiltonian such that all eigenvalues are in the intervall [0:1]
!!                      if (foe_data_get_real(foe_obj,"evlow",ispin)/=evlow_old .or. &
!!                          foe_data_get_real(foe_obj,"evhigh",ispin)/=evhigh_old) then
!!                          !!call scale_and_shift_matrix()
!!                          if (npl_auto_) then
!!                              call dscal(size(ovrlp_scaled%matrix_compr), eval_multiplicator, ovrlp_scaled%matrix_compr(1), 1)
!!                              eval_multiplicator_total = eval_multiplicator_total*eval_multiplicator
!!                              call scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, inv_ovrlp_smat, &
!!                                   ovrlp_smat, ovrlp_scaled, isshift, &
!!                                   matscal_compr=hamscal_compr, scale_factor=scale_factor, shift_value=shift_value)
!!                              if (iproc==0) then
!!                                  write(*,*) 'eval_multiplicator, eval_multiplicator_total', &
!!                                              eval_multiplicator, eval_multiplicator_total
!!                              end if
!!                          else
!!                              call scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, inv_ovrlp_smat, &
!!                                   ovrlp_smat, ovrlp_mat, isshift, &
!!                                   matscal_compr=hamscal_compr, scale_factor=scale_factor, shift_value=shift_value)
!!                          end if
!!                          calculate_SHS=.true.
!!                      else
!!                          calculate_SHS=.false.
!!                      end if
!!                      !!do i=1,size(ovrlp_mat%matrix_compr)
!!                      !!    write(900+iproc,*) i, ovrlp_mat%matrix_compr(i)
!!                      !!end do
!!                      !!do i=1,size(hamscal_compr)
!!                      !!    write(950+iproc,*) i, hamscal_compr(i)
!!                      !!end do
!!                      evlow_old=foe_data_get_real(foe_obj,"evlow",ispin)
!!                      evhigh_old=foe_data_get_real(foe_obj,"evhigh",ispin)
!!
!!
!!                      !call uncompress_matrix(iproc,ovrlp_smat,ovrlp_mat%matrix_compr,overlap)
!!                      !call dsyev('v', 'l', ovrlp_smat%nfvctr, overlap, ovrlp_smat%nfvctr, eval, work, lwork, info)
!!                      !if (iproc==0) write(*,*) 'ovrlp_mat%matrix_compr: eval low / high',eval(1), eval(ovrlp_smat%nfvctr)
!!                      !call uncompress_matrix(iproc,inv_ovrlp_smat,hamscal_compr,overlap)
!!                      !call dsyev('v', 'l', ovrlp_smat%nfvctr, overlap, ovrlp_smat%nfvctr, eval, work, lwork, info)
!!                      !if (iproc==0) write(*,*) 'hamscal_compr: eval low / high',eval(1), eval(ovrlp_smat%nfvctr)
!!
!!
!!                      ! Determine the degree of the polynomial
!!                      if (.not. npl_auto_) then
!!                          npl=nint(degree_multiplicator* &
!!                               (foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin)))
!!                          npl=max(npl,NPL_MIN)
!!                          npl_boundaries = nint(degree_multiplicator* &
!!                              (foe_data_get_real(foe_obj,"evhigh",ispin)-foe_data_get_real(foe_obj,"evlow",ispin)) &
!!                                  /foe_data_get_real(foe_obj,"fscale_lowerbound")) ! max polynomial degree for given eigenvalue boundaries
!!                          if (npl>npl_boundaries) then
!!                              npl=npl_boundaries
!!                              if (iproc==0) call yaml_warning('very sharp decay of error function, polynomial degree reached limit')
!!                              if (iproc==0) write(*,*) 'STOP SINCE THIS WILL CREATE PROBLEMS WITH NPL_CHECK'
!!                              stop
!!                          end if
!!                          if (npl>nplx) stop 'npl>nplx'
!!                      else
!!                          stop 'use the new wrapper for this'
!!                          !!call get_polynomial_degree(iproc, nproc, ispin, ncalc, ex, foe_obj, 5, 100, 1.d-10, &
!!                          !!     npl, cc, anoise)
!!                      end if
!!
!!                      ! Array that holds the Chebyshev polynomials. Needs to be recalculated
!!                      ! every time the Hamiltonian has been modified.
!!                      if (iproc==0 .and. verbosity_>0) then
!!                          call yaml_newline()
!!                          call yaml_mapping_open('ICE')
!!                          call yaml_map('eval bounds',&
!!                               (/foe_data_get_real(foe_obj,"evlow",ispin),foe_data_get_real(foe_obj,"evhigh",ispin)/),fmt='(f5.2)')
!!                          call yaml_map('mult.',degree_multiplicator,fmt='(f5.2)')
!!                          call yaml_map('pol. deg.',npl)
!!                      end if
!!                      if (calculate_SHS) then
!!                          call f_free_ptr(chebyshev_polynomials)
!!                          chebyshev_polynomials = f_malloc_ptr((/nsize_polynomial,npl/),id='chebyshev_polynomials')
!!                      end if
!!
!!
!!                      if (.not. npl_auto_) then
!!                          cc = f_malloc_ptr((/npl,3,ncalc/),id='cc')
!!
!!                          !!if (foe_data_get_real(foe_obj,"evlow")>=0.d0) then
!!                          !!    stop 'ERROR: lowest eigenvalue must be negative'
!!                          !!end if
!!                          if (foe_data_get_real(foe_obj,"evhigh",ispin)<=0.d0) then
!!                              stop 'ERROR: highest eigenvalue must be positive'
!!                          end if
!!
!!                          call timing(iproc, 'FOE_auxiliary ', 'OF')
!!                          call timing(iproc, 'chebyshev_coef', 'ON')
!!
!!                          max_error = f_malloc(ncalc,id='max_error')
!!                          x_max_error = f_malloc(ncalc,id='x_max_error')
!!                          mean_error = f_malloc(ncalc,id='mean_error')
!!                          do icalc=1,ncalc
!!                              call func_set(FUNCTION_POLYNOMIAL, powerx=ex(icalc))
!!                              call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
!!                                   foe_data_get_real(foe_obj,"evhigh",ispin), npl, func, cc(1:,1:,icalc), &
!!                                   x_max_error(icalc), max_error(icalc), mean_error(icalc))
!!                              !!##call cheb_exp(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
!!                              !!##     foe_data_get_real(foe_obj,"evhigh",ispin), npl, cc(1:,1:,icalc:), ex(icalc), &
!!                              !!##     x_max_error(icalc), max_error(icalc), mean_error(icalc))
!!                              !call chder(foe_data_get_real(foe_obj,"evlow",ispin), &
!!                              !     foe_data_get_real(foe_obj,"evhigh",ispin), cc(1:,1:,icalc:), cc(1:,2:,icalc:), npl)
!!                              call func_set(FUNCTION_EXPONENTIAL, betax=-40.d0, &
!!                                   muax=foe_data_get_real(foe_obj,"evlow",ispin), mubx=foe_data_get_real(foe_obj,"evhigh",ispin))
!!                              call get_chebyshev_expansion_coefficients(iproc, nproc, foe_data_get_real(foe_obj,"evlow",ispin), &
!!                                   foe_data_get_real(foe_obj,"evhigh",ispin), npl, func, cc(1:,2:,icalc), &
!!                                   x_max_error_fake, max_error_fake, mean_error_fake)
!!                              do ipl=1,npl
!!                                 cc(ipl,3,1) = -cc(ipl,2,1)
!!                              end do
!!                              !!##call chebyshev_coefficients_penalyfunction(foe_data_get_real(foe_obj,"evlow",ispin), &
!!                              !!##     foe_data_get_real(foe_obj,"evhigh",ispin), npl, cc(1:,2:,icalc:), max_error_fake)
!!                              call evnoise(npl, cc(1:,2:,icalc:), foe_data_get_real(foe_obj,"evlow",ispin), &
!!                                   foe_data_get_real(foe_obj,"evhigh",ispin), anoise)
!!                          end do
!!                          if (iproc==0 .and. verbosity_>0) then
!!                              call yaml_mapping_open('accuracy (x, max err, mean err)')
!!                              do icalc=1,ncalc
!!                                  call yaml_map('Operation '//trim(yaml_toa(icalc)), &
!!                                      (/x_max_error(icalc),max_error(icalc),mean_error(icalc)/),fmt='(es9.2)')
!!                              end do
!!                              call yaml_mapping_close()
!!                          end if
!!                          call f_free(mean_error)
!!                          call f_free(max_error)
!!                          call f_free(x_max_error)
!!
!!                          call timing(iproc, 'chebyshev_coef', 'OF')
!!                          call timing(iproc, 'FOE_auxiliary ', 'ON')
!!                      end if
!!                      if (iproc==0 .and. verbosity_>0) then
!!                          call yaml_mapping_close()
!!                      end if
!!
!!                      !!do j=1,npl
!!                      !!    write(*,*) 'in main: j, cc(j,1,1), cc(j,2,1)', j, cc(j,1,1), cc(j,2,1)
!!                      !!end do
!!
!!
!!                      call timing(iproc, 'FOE_auxiliary ', 'OF')
!!
!!                      emergency_stop=.false.
!!                      if (calculate_SHS) then
!!                          ! Passing inv_ovrlp(1)%matrix_compr as it will not be
!!                          ! used, to be improved...
!!                          call chebyshev_clean(iproc, nproc, npl, cc, &
!!                               inv_ovrlp_smat, hamscal_compr, &
!!                               .false., &
!!                               nsize_polynomial, ncalc, inv_ovrlp_matrixp_new, penalty_ev_new, chebyshev_polynomials, &
!!                               emergency_stop)
!!                           !write(*,*) 'sum(hamscal_compr)', sum(hamscal_compr)
!!                           !write(*,*) 'npl, sum(cc(:,1,1)), sum(chebyshev_polynomials)', &
!!                           !    npl, sum(cc(:,1,1)), sum(chebyshev_polynomials)
!!                           !write(*,*) 'sum(inv_ovrlp_matrixp_new)',sum(inv_ovrlp_matrixp_new)
!!                           !do i=1,size(inv_ovrlp_matrixp_new)
!!                           !    write(300,*) 'i, inv_ovrlp_matrixp_new(i)', i, inv_ovrlp_matrixp_new(i,1)
!!                           !end do
!!                           !!!@NEW#####################################################################################
!!                           !!if (mode==new) then
!!                           !!    call f_free_ptr(chebyshev_polynomials)
!!                           !!    call get_chebyshev_polynomials(iproc, nproc, 1, 2, npl, ovrlp_smat, inv_ovrlp_smat, &
!!                           !!                              ovrlp_scaled, foe_obj, chebyshev_polynomials, eval_bounds_ok)
!!                           !!    call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
!!                           !!         inv_ovrlp_smat%nfvctr, inv_ovrlp_smat%smmm%nfvctrp, &
!!                           !!         inv_ovrlp_smat, chebyshev_polynomials, ncalc, cc, inv_ovrlp_matrixp_new)
!!                           !!    penalty_ev_new = 0.d0
!!                           !!end if
!!                           !@END NEW#################################################################################
!!                           !!do i=1,size(inv_ovrlp_matrixp_new,1)
!!                           !!    write(400+iproc,*) i, inv_ovrlp_matrixp_new(i,1)
!!                           !!end do
!!                          if (inv_ovrlp_smat%smmm%nvctrp>0) then
!!                              do icalc=1,ncalc
!!                                  call transform_sparsity_pattern(inv_ovrlp_smat%nfvctr, &
!!                                       inv_ovrlp_smat%smmm%nvctrp_mm, inv_ovrlp_smat%smmm%isvctr_mm, &
!!                                       inv_ovrlp_smat%nseg, inv_ovrlp_smat%keyv, inv_ovrlp_smat%keyg, &
!!                                       inv_ovrlp_smat%smmm%line_and_column_mm, &
!!                                       inv_ovrlp_smat%smmm%nvctrp, inv_ovrlp_smat%smmm%isvctr, &
!!                                       inv_ovrlp_smat%smmm%nseg, inv_ovrlp_smat%smmm%keyv, inv_ovrlp_smat%smmm%keyg, &
!!                                       inv_ovrlp_smat%smmm%istsegline, 'large_to_small', &
!!                                       inv_ovrlp_matrixp_small_new(1,icalc), inv_ovrlp_matrixp_new(1,icalc))
!!                                 !!do i=1,size(inv_ovrlp_matrixp_small_new,1)
!!                                 !!    write(410+iproc,*) i, inv_ovrlp_matrixp_small_new(i,icalc)
!!                                 !!end do
!!                              end do
!!                          end if
!!                          !write(*,*) 'size(inv_ovrlp_matrixp_small_new), sum(inv_ovrlp_matrixp_small_new), ncalc', &
!!                          !     size(inv_ovrlp_matrixp_small_new), sum(inv_ovrlp_matrixp_small_new), ncalc
!!
!!                           !write(*,'(a,i5,2es24.8)') 'iproc, sum(inv_ovrlp_matrixp(:,:,1:2)', (sum(inv_ovrlp_matrixp(:,:,icalc)),icalc=1,ncalc)
!!                          !!do i=1,inv_ovrlp_smat%smmm%nvctrp
!!                          !!    ii = inv_ovrlp_smat%smmm%isvctr + i
!!                          !!    call get_line_and_column(ii, inv_ovrlp_smat%smmm%nseg, inv_ovrlp_smat%smmm%keyv, inv_ovrlp_smat%smmm%keyg, iline, icolumn)
!!                          !!    do icalc=1,ncalc
!!                          !!        inv_ovrlp_matrixp(icolumn,iline-inv_ovrlp_smat%smmm%isfvctr,icalc) = inv_ovrlp_matrixp_new(i,icalc)
!!                          !!    end do
!!                          !!    !!penalty_ev(icolumn,iline-inv_ovrlp_smat%smmm%isfvctr,1) = penalty_ev_new(i,1)
!!                          !!    !!penalty_ev(icolumn,iline-inv_ovrlp_smat%smmm%isfvctr,2) = penalty_ev_new(i,2)
!!                          !!end do
!!                      else
!!                          ! The Chebyshev polynomials are already available
!!                          !if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','from memory')
!!                          call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
!!                               inv_ovrlp_smat%nfvctr, inv_ovrlp_smat%smmm%nfvctrp, &
!!                               inv_ovrlp_smat, chebyshev_polynomials, ncalc, cc, inv_ovrlp_matrixp_new)
!!                          !do icalc=1,ncalc
!!                          !    write(*,*) 'sum(inv_ovrlp_matrixp_new(:,icalc))',sum(inv_ovrlp_matrixp_new(:,icalc))
!!                          !end do
!!                          !!do icalc=1,ncalc
!!                          !!    call uncompress_polynomial_vector(iproc, nproc, nsize_polynomial, &
!!                          !!         inv_ovrlp_smat, inv_ovrlp_matrixp_new, inv_ovrlp_matrixp(:,:,icalc))
!!                          !!end do
!!                      end if
!!
!!
!!
!!                     !!! Check for an emergency stop, which happens if the kernel explodes, presumably due
!!                     !!! to the eigenvalue bounds being too small.
!!                     !!call check_emergency_stop(nproc,emergency_stop)
!!                     !!if (emergency_stop) then
!!                     !!     eval_bounds_ok(1)=.false.
!!                     !!     call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)/1.2d0,ispin)
!!                     !!     eval_bounds_ok(2)=.false.
!!                     !!     call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*1.2d0,ispin)
!!                     !!     call f_free(cc)
!!                     !!     cycle main_loop
!!                     !!end if
!!
!!
!!                      call timing(iproc, 'FOE_auxiliary ', 'ON')
!!
!!
!!                      restart=.false.
!!
!!                      ! Check the eigenvalue bounds. Only necessary if calculate_SHS is true
!!                      ! (otherwise this has already been checked in the previous iteration).
!!                      if (calculate_SHS) then
!!                          !call check_eigenvalue_spectrum()
!!                          !!call check_eigenvalue_spectrum(nproc, inv_ovrlp_smat, ovrlp_smat, ovrlp_mat, 1, &
!!                          !!     0, 1.2d0, 1.d0/1.2d0, penalty_ev, anoise, .false., emergency_stop, &
!!                          !!     foe_obj, restart, eval_bounds_ok)
!!                          if (mode==old) then
!!                              call check_eigenvalue_spectrum_new(nproc, inv_ovrlp_smat, ispin, &
!!                                   0, 1.2d0, 1.d0/1.2d0, penalty_ev_new, anoise, .false., emergency_stop, &
!!                                   foe_obj, restart, eval_bounds_ok, &
!!                                   verbosity=verbosity_, eval_multiplicator=eval_multiplicator)
!!                          else if (mode==new) then
!!                              if (.not.eval_bounds_ok(1)) then
!!                                  ! lower bound too large
!!                                  call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)/1.2d0,ispin)
!!                                  restart=.true.
!!                                  eval_multiplicator = 2.0d0
!!                              else if (.not.eval_bounds_ok(2)) then
!!                                  ! upper bound too small
!!                                  call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*1.2d0,ispin)
!!                                  restart=.true.
!!                                  eval_multiplicator = 1.d0/2.0d0
!!                              end if
!!                          end if
!!                      end if
!!
!!                      call f_free_ptr(cc)
!!
!!                      if (restart) then
!!                          if(evbounds_shrinked) then
!!                              ! this shrink was not good, increase the saturation counter
!!                              call foe_data_set_int(foe_obj,"evboundsshrink_isatur", &
!!                                   foe_data_get_int(foe_obj,"evboundsshrink_isatur")+1)
!!                          end if
!!                          call foe_data_set_int(foe_obj,"evbounds_isatur",0)
!!                          cycle
!!                      end if
!!
!!                      ! eigenvalue bounds ok
!!                      if (calculate_SHS) then
!!                          call foe_data_set_int(foe_obj,"evbounds_isatur",foe_data_get_int(foe_obj,"evbounds_isatur")+1)
!!                      end if
!!
!!
!!                      exit
!!
!!
!!                  end do main_loop
!!
!!
!!
!!              !if (inv_ovrlp_smat%smmm%nvctrp>0) then
!!                  do icalc=1,ncalc
!!                      !!call compress_matrix_distributed(iproc, nproc, inv_ovrlp_smat, DENSE_MATMUL, inv_ovrlp_matrixp(1:,1:,icalc), &
!!                      !!     inv_ovrlp(icalc)%matrix_compr(ilshift2+1:))
!!                      call compress_matrix_distributed_wrapper(iproc, nproc, inv_ovrlp_smat, &
!!                           SPARSE_MATMUL_SMALL, inv_ovrlp_matrixp_small_new(:,icalc), &
!!                           inv_ovrlp(icalc)%matrix_compr(ilshift2+1:))
!!                      !write(*,*) 'sum(inv_ovrlp(icalc)%matrix_compr)', sum(inv_ovrlp(icalc)%matrix_compr)
!!                  end do
!!              !end if
!!
!!              if (npl_auto_) then
!!                  do icalc=1,ncalc
!!                      call dscal(inv_ovrlp_smat%nvctrp_tg, 1.d0/eval_multiplicator**ex(icalc), &
!!                           inv_ovrlp(icalc)%matrix_compr(ilshift2+1), 1)
!!                      !write(*,*) 'sum(inv_ovrlp(icalc)%matrix_compr)', sum(inv_ovrlp(icalc)%matrix_compr)
!!                  end do
!!              end if
!!
!!
!!          end do spin_loop
!!
!!      !call f_free_ptr(inv_ovrlp_matrixp)
!!      call f_free(inv_ovrlp_matrixp_small_new)
!!      call f_free(inv_ovrlp_matrixp_new)
!!      call f_free_ptr(chebyshev_polynomials)
!!      !!call f_free(penalty_ev)
!!      call f_free(penalty_ev_new)
!!      call f_free(hamscal_compr)
!!      if (npl_auto_) then
!!          call deallocate_matrices(ovrlp_scaled)
!!      end if
!!
!!      call f_free_ptr(foe_obj%ef)
!!      call f_free_ptr(foe_obj%evlow)
!!      call f_free_ptr(foe_obj%evhigh)
!!      call f_free_ptr(foe_obj%bisection_shift)
!!      call f_free_ptr(foe_obj%charge)
!!
!!      call timing(iproc, 'FOE_auxiliary ', 'OF')
!!
!!      call f_release_routine()
!!
!!
!!    end subroutine inverse_chebyshev_expansion




!!    subroutine get_minmax_eigenvalues(iproc, ovrlp_smat, ovrlp_mat)
!!      implicit none
!!
!!      ! Calling arguments
!!      integer, intent(in) :: iproc
!!      type(sparse_matrix), intent(in) :: ovrlp_smat
!!      type(matrices), intent(in) :: ovrlp_mat
!!
!!      ! Local variables
!!      integer :: iseg, ii, i, lwork, info
!!      real(kind=mp),dimension(:,:),allocatable :: tempmat
!!      real(kind=mp),dimension(:),allocatable :: eval, work
!!
!!      call f_routine(id='get_minmax_eigenvalues')
!!
!!      tempmat = f_malloc0((/ovrlp_smat%nfvctr,ovrlp_smat%nfvctr/),id='tempmat')
!!      do iseg=1,ovrlp_smat%nseg
!!          ii=ovrlp_smat%keyv(iseg)
!!          do i=ovrlp_smat%keyg(1,1,iseg),ovrlp_smat%keyg(2,1,iseg)
!!              tempmat(i,ovrlp_smat%keyg(1,2,iseg)) = ovrlp_mat%matrix_compr(ii)
!!              ii = ii + 1
!!          end do
!!      end do
!!      !!if (iproc==0) then
!!      !!    do i=1,ovrlp_smat%nfvctr
!!      !!        do j=1,ovrlp_smat%nfvctr
!!      !!            write(*,'(a,2i6,es17.8)') 'i,j,val',i,j,tempmat(j,i)
!!      !!        end do
!!      !!    end do
!!      !!end if
!!      eval = f_malloc(ovrlp_smat%nfvctr,id='eval')
!!      lwork=100*ovrlp_smat%nfvctr
!!      work = f_malloc(lwork,id='work')
!!      call dsyev('n','l', ovrlp_smat%nfvctr, tempmat, ovrlp_smat%nfvctr, eval, work, lwork, info)
!!      !if (iproc==0) write(*,*) 'eval',eval
!!      if (iproc==0) call yaml_map('eval max/min',(/eval(1),eval(ovrlp_smat%nfvctr)/),fmt='(es16.6)')
!!
!!      call f_free(tempmat)
!!      call f_free(eval)
!!      call f_free(work)
!!
!!      call f_release_routine()
!!
!!    end subroutine get_minmax_eigenvalues


    subroutine inverse_chebyshev_expansion_new(iproc, nproc, comm, &
               ovrlp_smat, inv_ovrlp_smat, ncalc, ex, ovrlp_mat, inv_ovrlp, &
               verbosity, npl_auto, ice_objx)
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: compress_matrix, uncompress_matrix, compress_matrix_distributed_wrapper, &
                              transform_sparsity_pattern, get_minmax_eigenvalues
      use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                          foe_data_set_logical, foe_data_get_logical, &
                          foe_data_deallocate
      use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
                             fermilevel_get_real, fermilevel_get_logical
      use chebyshev, only: chebyshev_clean, chebyshev_fast
      use foe_common, only: evnoise, get_chebyshev_expansion_coefficients, &
                            get_bounds_and_polynomials, &
                            init_foe
      use module_func
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, ncalc
      type(sparse_matrix), intent(in) :: ovrlp_smat, inv_ovrlp_smat
      real(kind=mp),dimension(ncalc),intent(in) :: ex
      type(matrices), intent(in) :: ovrlp_mat
      type(matrices),dimension(ncalc),intent(out) :: inv_ovrlp
      integer, intent(in),optional :: verbosity
      logical, intent(in),optional :: npl_auto
      type(foe_data),intent(inout),target,optional :: ice_objx

      ! Local variables
      integer :: npl
      integer :: nsize_polynomial
      !integer :: npl_boundaries
      integer, parameter :: nplx=50000
      real(kind=mp),dimension(:,:,:),pointer :: chebyshev_polynomials
      !real(kind=mp),dimension(:,:,:),pointer :: inv_ovrlp_matrixp
      !real(kind=mp),dimension(:,:,:),allocatable :: penalty_ev
      real(kind=mp),dimension(:,:,:),pointer :: cc
      real(kind=mp) :: scale_factor, shift_value, betax
      real(kind=mp) :: tt
      !real(kind=mp),dimension(ovrlp_smat%nspin) :: eval_min, eval_max
      !real(kind=mp) :: x_max_error_fake 
      real(kind=mp) :: eval_multiplicator, eval_multiplicator_total, bounds_limit
      !logical :: restart
      !logical, dimension(2) :: emergency_stop
      !real(kind=mp),dimension(2) :: allredarr
      real(kind=mp),dimension(:),allocatable :: hamscal_compr
      !logical, dimension(2) :: eval_bounds_ok
      !integer, dimension(2) :: irowcol
      integer :: ispin, isshift, ilshift, ilshift2, verbosity_, npl_min_fake
      logical :: evbounds_shrinked
      real(kind=mp),parameter :: DEGREE_MULTIPLICATOR_MAX=20.d0
      real(kind=mp) :: accuracy_function, accuracy_penalty
      integer, parameter :: SPARSE=1
      integer, parameter :: DENSE=2
      integer, parameter :: imode=SPARSE
      type(foe_data),pointer :: ice_obj
      type(foe_data),target :: ice_obj_
      real(kind=mp),dimension(:),allocatable :: x_max_error, max_error, mean_error, charge_fake
      !real(kind=mp),dimension(:,:),allocatable :: tempmat
      integer :: icalc
      real(kind=mp),dimension(:,:),allocatable :: inv_ovrlp_matrixp_new
      !real(kind=mp),dimension(:,:),allocatable :: penalty_ev_new
      real(kind=mp),dimension(:,:),allocatable :: inv_ovrlp_matrixp_small_new
      !type(matrices) :: ovrlp_scaled
      character(len=3), parameter :: old='old'
      character(len=3), parameter :: new='new'
      !character(len=3) :: mode=old
      !!integer,parameter :: NPL_MIN = 5
      !!integer,parameter :: NPL_MAX = 5000
      !!integer,parameter :: NPL_STRIDE = 5
      integer :: npl_min, npl_max, npl_stride

      call f_routine(id='inverse_chebyshev_expansion_new')
      call f_timing(TCAT_CME_AUXILIARY,'ON')


      if (.not.inv_ovrlp_smat%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      if (present(verbosity)) then
          verbosity_ = verbosity
      else
          verbosity_ = 100
      end if

      inv_ovrlp_matrixp_new = f_malloc((/inv_ovrlp_smat%smmm%nvctrp,ncalc/),id='inv_ovrlp_matrixp_new')
      inv_ovrlp_matrixp_small_new = f_malloc((/inv_ovrlp_smat%smmm%nvctrp_mm,ncalc/),id='inv_ovrlp_matrixp_small_new')
      hamscal_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_TASKGROUP, id='hamscal_compr')

      if (present(ice_objx)) then
          ice_obj => ice_objx
      else
          charge_fake = f_malloc0(ovrlp_smat%nspin,id='charge_fake')
          call init_foe(iproc, nproc, ovrlp_smat%nspin, charge=charge_fake, &
               evlow=0.5_mp, evhigh=1.5_mp, foe_obj=ice_obj_, &
               accuracy_function=1.e-8_mp, accuracy_penalty=1.e-5_mp)
          call f_free(charge_fake)
          ice_obj => ice_obj_
      end if

      npl_min = foe_data_get_int(ice_obj,"npl_min")
      npl_max = foe_data_get_int(ice_obj,"npl_max")
      npl_stride = foe_data_get_int(ice_obj,"npl_stride")
      betax = foe_data_get_real(ice_obj,"betax")

      !!@ TEMPORARY: eigenvalues of  the overlap matrix ###################
      !!call get_minmax_eigenvalues(iproc, nproc, comm, 'standard', -8, &
      !!     ovrlp_smat, ovrlp_mat, eval_min, eval_max)


      ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
      nsize_polynomial = inv_ovrlp_smat%smmm%nvctrp_mm

      max_error = f_malloc(ncalc,id='max_error')
      x_max_error = f_malloc(ncalc,id='x_max_error')
      mean_error = f_malloc(ncalc,id='mean_error')

      evbounds_shrinked=.false.
      ! try to decrease the eigenvalue spectrum a bit
      if (foe_data_get_int(ice_obj,"evbounds_isatur")>foe_data_get_int(ice_obj,"evbounds_nsatur") .and. &
          foe_data_get_int(ice_obj,"evboundsshrink_isatur")<=foe_data_get_int(ice_obj,"evboundsshrink_nsatur")) then
          do ispin=1,inv_ovrlp_smat%nspin
              ! Make sure that the lower bound is not larger than the upper bound
              bounds_limit = 0.5_mp*(foe_data_get_real(ice_obj,"evlow",ispin)+foe_data_get_real(ice_obj,"evhigh",ispin))
              tt = min(1.d0/0.9d0*foe_data_get_real(ice_obj,"evlow",ispin),bounds_limit-1.e-2_mp)
              call foe_data_set_real(ice_obj,"evlow",tt,ispin)
              tt = max(0.9d0*foe_data_get_real(ice_obj,"evhigh",ispin),bounds_limit+1.e-2_mp)
              call foe_data_set_real(ice_obj,"evhigh",tt,ispin)
          end do
          evbounds_shrinked=.true.
      else
          evbounds_shrinked=.false.
      end if



      eval_multiplicator = foe_data_get_real(ice_obj,"eval_multiplicator",1)
      eval_multiplicator_total = 1.d0


      npl_min_fake = NPL_MIN !since intent(inout)
      accuracy_function = foe_data_get_real(ice_obj,"accuracy_function")
      accuracy_penalty = foe_data_get_real(ice_obj,"accuracy_penalty")
      call get_bounds_and_polynomials(iproc, nproc, comm, 1, 1, NPL_MAX, NPL_STRIDE, &
           ncalc, FUNCTION_POLYNOMIAL, accuracy_function, accuracy_penalty, .true., 1.0_mp/1.2_mp, 1.2_mp, verbosity_, &
           ovrlp_smat, inv_ovrlp_smat, ovrlp_mat, ice_obj, npl_min_fake, &
           inv_ovrlp(1)%matrix_compr, chebyshev_polynomials, &
           npl, scale_factor, shift_value, hamscal_compr, &
           ex=ex, scaling_factor_low=1.2_mp, scaling_factor_up=1.0_mp/1.2_mp, &
           eval_multiplicator=eval_multiplicator, eval_multiplicator_total=eval_multiplicator_total, &
           cc=cc, max_errorx=max_error)



      spin_loop: do ispin=1,ovrlp_smat%nspin
          isshift=(ispin-1)*ovrlp_smat%nvctr
          ilshift=(ispin-1)*inv_ovrlp_smat%nvctr
          ilshift2=(ispin-1)*inv_ovrlp_smat%nvctrp_tg
          call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
               inv_ovrlp_smat%nfvctr, inv_ovrlp_smat%smmm%nfvctrp, &
               inv_ovrlp_smat, chebyshev_polynomials(:,:,ispin), ncalc, cc, inv_ovrlp_matrixp_small_new)
          do icalc=1,ncalc
              call compress_matrix_distributed_wrapper(iproc, nproc, inv_ovrlp_smat, &
                   SPARSE_MATMUL_SMALL, inv_ovrlp_matrixp_small_new(:,icalc), &
                   ONESIDED_FULL, inv_ovrlp(icalc)%matrix_compr(ilshift2+1:))
              call dscal(inv_ovrlp_smat%nvctrp_tg, 1.d0/eval_multiplicator_total**ex(icalc), &
                   inv_ovrlp(icalc)%matrix_compr(ilshift2+1), 1)
          end do
      end do spin_loop

      call f_free_ptr(cc)
      call f_free_ptr(chebyshev_polynomials)

      call foe_data_set_real(ice_obj,"eval_multiplicator",eval_multiplicator_total,1)


      call f_free(inv_ovrlp_matrixp_small_new)
      call f_free(inv_ovrlp_matrixp_new)
      call f_free(hamscal_compr)
      call f_free(max_error)
      call f_free(x_max_error)
      call f_free(mean_error)

      if (.not.present(ice_objx)) then
          call foe_data_deallocate(ice_obj_)
      end if

      call f_timing(TCAT_CME_AUXILIARY,'OF')
      call f_release_routine()

    end subroutine inverse_chebyshev_expansion_new


    !!subroutine adjust_eval_bounds
    !!  implicit none

    !!  ! Calling arguments
    !!  logical, dimension(2), intent(in) :: eval_bounds_ok
    !!  type(foe_data),intent(inout) :: ice_obj

    !!  if (.not. eval_bounds_ok(1)) then
    !!      ! Lower bounds too large

    !!end subroutine adjust_eval_bounds




    subroutine calculate_fermi_function_entropy(iproc, nproc, comm, &
               ovrlp_smat, kernel_smat, entropykernel_smat, &
               ovrlp_mat, kernel_mat, ovrlp_minus_one_half_mat, &
               accuracy_entropy, entropykernel_mat, eTS, eTS_check, &
               verbosity, ice_objx)
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: compress_matrix, uncompress_matrix, compress_matrix_distributed_wrapper, &
                              transform_sparsity_pattern, get_minmax_eigenvalues, &
                              trace_sparse_matrix
      use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                          foe_data_set_logical, foe_data_get_logical, &
                          foe_data_deallocate
      use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
                             fermilevel_get_real, fermilevel_get_logical
      use chebyshev, only: chebyshev_clean, chebyshev_fast
      use foe_common, only: evnoise, get_chebyshev_expansion_coefficients, &
                            get_bounds_and_polynomials, &
                            init_foe
      use module_func
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      type(sparse_matrix), intent(in) :: kernel_smat, entropykernel_smat, ovrlp_smat
      type(matrices), intent(in) :: kernel_mat, ovrlp_mat, ovrlp_minus_one_half_mat
      real(mp),intent(in) :: accuracy_entropy
      type(matrices),intent(out) :: entropykernel_mat
      real(mp),intent(out) :: eTS, eTS_check
      integer, intent(in),optional :: verbosity
      type(foe_data),intent(inout),target,optional :: ice_objx

      ! Local variables
      integer :: npl, jorb, it, ii, iseg
      integer :: isegstart, isegend, iismall, nsize_polynomial
      integer :: iismall_ovrlp, iismall_ham, npl_boundaries, i, ipl
      integer, parameter :: nplx=50000
      real(kind=mp),dimension(:,:,:),pointer :: chebyshev_polynomials
      real(kind=mp),dimension(:,:,:),pointer :: inv_ovrlp_matrixp
      real(kind=mp),dimension(:,:,:),allocatable :: penalty_ev
      real(kind=mp),dimension(:,:,:),pointer :: cc
      real(kind=mp) :: anoise, scale_factor, shift_value, betax
      real(kind=mp) :: evlow_old, evhigh_old, tt
      real(kind=mp),dimension(kernel_smat%nspin) :: eval_min, eval_max
      real(kind=mp) :: x_max_error_fake, max_error_fake, mean_error_fake
      real(kind=mp) :: tt_ovrlp, tt_ham, eval_multiplicator, eval_multiplicator_total, bounds_limit
      logical :: restart, calculate_SHS
      logical, dimension(2) :: emergency_stop
      real(kind=mp),dimension(2) :: allredarr
      real(kind=mp),dimension(:),allocatable :: hamscal_compr
      logical, dimension(2) :: eval_bounds_ok
      integer, dimension(2) :: irowcol
      integer :: irow, icol, iflag, ispin, isshift, ilshift, ilshift2, verbosity_, npl_min_fake
      logical :: overlap_calculated, evbounds_shrinked, degree_sufficient, reached_limit, npl_auto_
      real(kind=mp),parameter :: DEGREE_MULTIPLICATOR_MAX=20.d0
      real(kind=mp) :: degree_multiplicator, accuracy_function, accuracy_penalty, sumn, val
      integer, parameter :: SPARSE=1
      integer, parameter :: DENSE=2
      integer, parameter :: imode=SPARSE
      type(foe_data),pointer :: ice_obj
      type(foe_data),target :: ice_obj_
      real(kind=mp),dimension(:),allocatable :: eval, work, x_max_error, max_error, mean_error, charge_fake
      real(kind=mp),dimension(:,:),allocatable :: tempmat
      integer :: lwork, info, j, icalc, iline, icolumn
      real(kind=mp),dimension(:),allocatable :: inv_ovrlp_matrixp_new
      real(kind=mp),dimension(:,:),allocatable :: penalty_ev_new
      real(kind=mp),dimension(:),allocatable :: inv_ovrlp_matrixp_small_new
      type(matrices) :: ovrlp_scaled
      character(len=3), parameter :: old='old'
      character(len=3), parameter :: new='new'
      character(len=3) :: mode=old
      !!integer,parameter :: NPL_MIN = 5
      !!integer,parameter :: NPL_MAX = 5000
      !!integer,parameter :: NPL_STRIDE = 5
      integer :: npl_min, npl_max, npl_stride
      integer,parameter :: ncalc = 1
      real(mp),dimension(:),allocatable :: evals

      call f_routine(id='inverse_chebyshev_expansion_new')
      call f_timing(TCAT_CME_AUXILIARY,'ON')


      if (.not.entropykernel_smat%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      if (present(verbosity)) then
          verbosity_ = verbosity
      else
          verbosity_ = 100
      end if

      inv_ovrlp_matrixp_new = f_malloc((/entropykernel_smat%smmm%nvctrp/),id='inv_ovrlp_matrixp_new')
      inv_ovrlp_matrixp_small_new = f_malloc((/entropykernel_smat%smmm%nvctrp_mm/),id='inv_ovrlp_matrixp_small_new')
      hamscal_compr = sparsematrix_malloc(entropykernel_smat, iaction=SPARSE_TASKGROUP, id='hamscal_compr')

      if (present(ice_objx)) then
          ice_obj => ice_objx
      else
          charge_fake = f_malloc0(kernel_smat%nspin,id='charge_fake')
          call init_foe(iproc, nproc, kernel_smat%nspin, charge=charge_fake, &
               evlow=-0.1_mp, evhigh=1.1_mp, foe_obj=ice_obj_, &
               accuracy_function=accuracy_entropy, accuracy_penalty=1.e-5_mp)
          call f_free(charge_fake)
          ice_obj => ice_obj_
      end if

      npl_min = foe_data_get_int(ice_obj,"npl_min")
      npl_max = foe_data_get_int(ice_obj,"npl_max")
      npl_stride = foe_data_get_int(ice_obj,"npl_stride")
      betax = foe_data_get_real(ice_obj,"betax")

      !!@ TEMPORARY: eigenvalues of  the overlap matrix ###################
      !!evals = f_malloc(kernel_smat%nfvctr,id='evals')
      !!call get_minmax_eigenvalues(iproc, nproc, comm, 'generalized', -8, &
      !!     kernel_smat, kernel_mat, eval_min, eval_max, &
      !!     smat2=ovrlp_smat, mat2=ovrlp_mat, &
      !!     evals=evals)
      !!eTS_check = 0.d0
      !!sumn = 0.d0
      !!do i=1,kernel_smat%nfvctr
      !!    val = evals(i)
      !!    sumn = sumn + val
      !!    if (val<1.e-30_mp .or. (val-1._mp)>-1.e-30_mp) then
      !!        tt = 0._mp
      !!    else
      !!        tt = -(val*log(val) + (1-val)*log(1-val))
      !!    end if
      !!    write(*,*) 'i, val, tt', i, val, tt
      !!    eTS_check = eTS_check + tt
      !!end do
      !!write(*,*) 'sumn',sumn
      !!call f_free(evals)
      !! ##################################################################


      ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
      nsize_polynomial = entropykernel_smat%smmm%nvctrp_mm

      max_error = f_malloc(ncalc,id='max_error')
      x_max_error = f_malloc(ncalc,id='x_max_error')
      mean_error = f_malloc(ncalc,id='mean_error')

      evbounds_shrinked=.false.
      ! try to decrease the eigenvalue spectrum a bit
      if (foe_data_get_int(ice_obj,"evbounds_isatur")>foe_data_get_int(ice_obj,"evbounds_nsatur") .and. &
          foe_data_get_int(ice_obj,"evboundsshrink_isatur")<=foe_data_get_int(ice_obj,"evboundsshrink_nsatur")) then
          do ispin=1,entropykernel_smat%nspin
              ! Make sure that the lower bound is not larger than the upper bound
              bounds_limit = 0.5_mp*(foe_data_get_real(ice_obj,"evlow",ispin)+foe_data_get_real(ice_obj,"evhigh",ispin))
              tt = min(1.d0/0.9d0*foe_data_get_real(ice_obj,"evlow",ispin),bounds_limit-1.e-2_mp)
              call foe_data_set_real(ice_obj,"evlow",tt,ispin)
              tt = max(0.9d0*foe_data_get_real(ice_obj,"evhigh",ispin),bounds_limit+1.e-2_mp)
              call foe_data_set_real(ice_obj,"evhigh",tt,ispin)
          end do
          evbounds_shrinked=.true.
      else
          evbounds_shrinked=.false.
      end if



      !!eval_multiplicator = foe_data_get_real(ice_obj,"eval_multiplicator",1)
      !!eval_multiplicator_total = 1.d0

      eval_multiplicator = 1.00_mp
      eval_multiplicator_total = 1.d0

      npl_min_fake = NPL_MIN !since intent(inout)
      !accuracy_function = 1.d-4 !foe_data_get_real(ice_obj,"accuracy_function")
      accuracy_function = foe_data_get_real(ice_obj,"accuracy_function")
      accuracy_penalty = foe_data_get_real(ice_obj,"accuracy_penalty")
      !!write(*,*) 'verbosity_', verbosity_
      call get_bounds_and_polynomials(iproc, nproc, comm, 2, 1, NPL_MAX, NPL_STRIDE, &
           ncalc, FUNCTION_FERMIFUNCTION_ENTROPY, accuracy_function, accuracy_penalty, .false., 1.2_mp, 1.2_mp, verbosity_, &
           kernel_smat, entropykernel_smat, kernel_mat, ice_obj, npl_min_fake, &
           entropykernel_mat%matrix_compr, chebyshev_polynomials, &
           npl, scale_factor, shift_value, hamscal_compr, &
           scaling_factor_low=1.2_mp, scaling_factor_up=1.0_mp/1.2_mp, &
           eval_multiplicator=eval_multiplicator, eval_multiplicator_total=eval_multiplicator_total, &
           cc=cc, max_errorx=max_error, &
           smats=ovrlp_smat, ovrlp_=ovrlp_mat, ovrlp_minus_one_half_=ovrlp_minus_one_half_mat)



      spin_loop: do ispin=1,kernel_smat%nspin
          isshift=(ispin-1)*kernel_smat%nvctr
          ilshift=(ispin-1)*entropykernel_smat%nvctr
          ilshift2=(ispin-1)*entropykernel_smat%nvctrp_tg
          call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
               entropykernel_smat%nfvctr, entropykernel_smat%smmm%nfvctrp, &
               entropykernel_smat, chebyshev_polynomials(:,:,ispin), ncalc, cc, inv_ovrlp_matrixp_small_new)
          !do icalc=1,ncalc
              call compress_matrix_distributed_wrapper(iproc, nproc, entropykernel_smat, &
                   SPARSE_MATMUL_SMALL, inv_ovrlp_matrixp_small_new, &
                   ONESIDED_FULL, entropykernel_mat%matrix_compr(ilshift2+1:))
              !!call dscal(entropykernel_smat%nvctrp_tg, 1.d0/eval_multiplicator_total**ex(icalc), &
              !!     entropykernel_mat%matrix_compr(ilshift2+1), 1)
          !end do
      end do spin_loop

      call f_free_ptr(cc)
      call f_free_ptr(chebyshev_polynomials)

      !!call foe_data_set_real(ice_obj,"eval_multiplicator",eval_multiplicator_total,1)


      call f_free(inv_ovrlp_matrixp_small_new)
      call f_free(inv_ovrlp_matrixp_new)
      call f_free(hamscal_compr)
      call f_free(max_error)
      call f_free(x_max_error)
      call f_free(mean_error)

      sumn = trace_sparse_matrix(iproc, nproc, comm, kernel_smat, kernel_mat%matrix_compr)
      eTS = trace_sparse_matrix(iproc, nproc, comm, kernel_smat, entropykernel_mat%matrix_compr)
      !eTS = eTS/scale_factor+shift_value*sumn

      if (.not.present(ice_objx)) then
          call foe_data_deallocate(ice_obj_)
      end if

      call f_timing(TCAT_CME_AUXILIARY,'OF')
      call f_release_routine()

    end subroutine calculate_fermi_function_entropy


end module ice
