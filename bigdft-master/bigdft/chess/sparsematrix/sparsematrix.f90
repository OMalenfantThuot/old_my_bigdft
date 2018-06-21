!> @file
!!   Basic routines dealing with the sparse matrices
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


!> Module to deal with the sparse matrices
module sparsematrix
  use dictionaries, only: f_err_throw
  use yaml_strings
  use sparsematrix_base
  use time_profiling
  use wrapper_mpi
  use wrapper_linalg
  use f_utils
  implicit none

  private


  !> Public routines
  public :: compress_matrix, compress_matrix2
  public :: uncompress_matrix, uncompress_matrix2
  public :: check_matrix_compression
  public :: transform_sparse_matrix!, transform_sparse_matrix_local
  public :: compress_matrix_distributed_wrapper
  public :: uncompress_matrix_distributed2
  public :: sequential_acces_matrix_fast, sequential_acces_matrix_fast2
  public :: sparsemm_new
  public :: gather_matrix_from_taskgroups, gather_matrix_from_taskgroups_inplace
  public :: extract_taskgroup_inplace, extract_taskgroup
  public :: write_matrix_compressed
  public :: write_sparsematrix
  public :: write_sparsematrix_CCS
  public :: transform_sparsity_pattern
  public :: matrix_matrix_mult_wrapper
  public :: trace_sparse_matrix
  public :: trace_sparse_matrix_product
  public :: delete_coupling_terms
  public :: synchronize_matrix_taskgroups
  public :: max_asymmetry_of_matrix
  public :: symmetrize_matrix
  public :: check_deviation_from_unity_sparse
  !!public :: operation_using_dense_lapack
  public :: matrix_power_dense_lapack
  public :: diagonalizeHamiltonian2
  public :: get_minmax_eigenvalues
  public :: resize_matrix_to_taskgroup


  interface compress_matrix_distributed_wrapper
      module procedure compress_matrix_distributed_wrapper_1
      module procedure compress_matrix_distributed_wrapper_2
  end interface compress_matrix_distributed_wrapper

  contains

    !> subroutine to compress the matrix to sparse form
    subroutine compress_matrix(iproc,nproc,sparsemat,inmat,outmat)
      use dynamic_memory
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: iproc, nproc
      type(sparse_matrix),intent(in) :: sparsemat
      real(kind=mp),dimension(sparsemat%nfvctr,sparsemat%nfvctr,sparsemat%nspin),target,intent(in) :: inmat
      real(kind=mp),dimension(sparsemat%nvctr*sparsemat%nspin),target,intent(out) :: outmat
    
      ! Local variables
      integer :: iseg, j, jj, ishift, ispin
      !integer,dimension(2) :: irowcol
      !integer :: ierr, irow, jcol, jjj
      real(kind=mp),dimension(:,:,:),pointer :: inm
      real(kind=mp),dimension(:),pointer :: outm
      
      call f_routine(id='compress_matrix')
      !call timing(iproc,'compressd_mcpy','ON')
      call f_timing(TCAT_SMAT_COMPRESSION,'ON')

      !if (present(outmat)) then
      !    if (sparsemat%parallel_compression/=0 .and. bigdft_mpi%nproc>1) then
      !        stop 'outmat not allowed for the given options'
      !    end if
          outm => outmat
      !else
      !    outm => sparsemat%matrix_compr
      !end if

      !if (present(inmat)) then
      !    if (sparsemat%parallel_compression/=0 .and. bigdft_mpi%nproc>1) then
      !        stop 'in not allowed for the given options'
      !    end if
          inm => inmat
      !else
      !    inm => sparsemat%matrix
      !end if
    
    
      if (sparsemat%parallel_compression==0.or.nproc==1) then
         do ispin=1,sparsemat%nspin
             ishift=(ispin-1)*sparsemat%nvctr
             !OpenMP broken on Vesta
             !!! !$omp parallel default(none) private(irowcol) &
             !$omp parallel default(none) private(iseg,j,jj) &
             !$omp shared(sparsemat,inm,outm,ishift,ispin)
             !$omp do
             do iseg=1,sparsemat%nseg
                 jj=sparsemat%keyv(iseg)
                 ! A segment is always on one line, therefore no double loop
                 do j=sparsemat%keyg(1,1,iseg),sparsemat%keyg(2,1,iseg)
                    !irow = sparsemat%orb_from_index(1,jj)
                    !jcol = sparsemat%orb_from_index(2,jj)
                    !irowcol = orb_from_index(sparsemat, j)
                    !write(*,*) 'iseg, j, jj', iseg, j, jj
                    outm(jj+ishift)=inm(j,sparsemat%keyg(1,2,iseg),ispin)
                    jj=jj+1
                 end do
             end do
             !$omp end do
             !$omp end parallel
         end do
      else if (sparsemat%parallel_compression==1) then
         stop 'this needs to be fixed'
         !!#call to_zero(sparsemat%nvctr, sparsemat%matrix_compr(1))
         !!#!$omp parallel do default(private) shared(sparsemat)
         !!#do jj=1,sparsemat%nvctrp
         !!#   jjj=jj+sparsemat%isvctr
         !!#   irow = sparsemat%orb_from_index(1,jjj)
         !!#   jcol = sparsemat%orb_from_index(2,jjj)
         !!#   sparsemat%matrix_compr(jjj)=sparsemat%matrix(irow,jcol)
         !!#end do
         !!#!$omp end parallel do
         !!#if (bigdft_mpi%nproc > 1) then
         !!#   call fmpi_allreduce(sparsemat%matrix_compr(1), sparsemat%nvctr, FMPI_SUM, bigdft_mpi%mpi_comm)
         !!#end if
      else
         stop 'this needs to be fixed'
         !!#sparsemat%matrix_comprp=f_malloc_ptr((sparsemat%nvctrp),id='sparsemat%matrix_comprp')
         !!#!$omp parallel do default(private) shared(sparsemat)
         !!#do jj=1,sparsemat%nvctrp
         !!#   jjj=jj+sparsemat%isvctr
         !!#   irow = sparsemat%orb_from_index(1,jjj)
         !!#   jcol = sparsemat%orb_from_index(2,jjj)
         !!#   sparsemat%matrix_comprp(jj)=sparsemat%matrix(irow,jcol)
         !!#end do
         !!#!$omp end parallel do
         !!#if (bigdft_mpi%nproc > 1) &
         !!#   & call mpi_allgatherv(sparsemat%matrix_comprp, sparsemat%nvctrp, &
         !!#   &    mpi_double_precision, sparsemat%matrix_compr, sparsemat%nvctr_par(:), &
         !!#   &    sparsemat%isvctr_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
         !!#call f_free_ptr(sparsemat%matrix_comprp)
      end if
    
      !call timing(iproc,'compressd_mcpy','OF')
      call f_timing(TCAT_SMAT_COMPRESSION,'OF')
      call f_release_routine()
    
    end subroutine compress_matrix


    subroutine compress_matrix2(iproc, nproc, sparsemat, inmat, outmat)
      use dynamic_memory
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: iproc, nproc
      type(sparse_matrix),intent(in) :: sparsemat
      real(kind=mp),dimension(sparsemat%nfvctr,sparsemat%nfvctr,sparsemat%nspin),intent(in) :: inmat
      real(kind=mp),dimension(sparsemat%nvctrp_tg*sparsemat%nspin),intent(out) :: outmat

      ! Local variables
      real(kind=mp),dimension(:),allocatable :: tmparr

      tmparr = sparsematrix_malloc(sparsemat,iaction=SPARSE_FULL,id='tmparr')
      call compress_matrix(iproc,nproc,sparsemat,inmat,tmparr)
      call extract_taskgroup(sparsemat, tmparr, outmat)
      call f_free(tmparr)
    end subroutine compress_matrix2


    !> subroutine to uncompress the matrix from sparse form
    subroutine uncompress_matrix(iproc,nproc,sparsemat,inmat,outmat)
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: iproc, nproc
      type(sparse_matrix), intent(in) :: sparsemat
      real(kind=mp),dimension(sparsemat%nvctr*sparsemat%nspin),target,intent(in) :: inmat
      real(kind=mp),dimension(sparsemat%nfvctr,sparsemat%nfvctr,sparsemat%nspin),target,intent(inout) :: outmat
      
      ! Local variables
      integer :: iseg, i, ii, ishift, ispin
      !integer, dimension(2) :: irowcol
      !integer ::  jcol, irow, iii, ierr
      real(kind=mp),dimension(:),pointer :: inm
      real(kind=mp),dimension(:,:,:),pointer :: outm

      !!if (present(outmat)) then
      !!    if (sparsemat%parallel_compression/=0 .and. bigdft_mpi%nproc>1) then
      !!        stop 'outmat not allowed for the given options'
      !!    end if
          outm => outmat
      !!else
      !!    outm => sparsemat%matrix
      !!end if

      !!if (present(inmat)) then
      !!    if (sparsemat%parallel_compression/=0 .and. bigdft_mpi%nproc>1) then
      !!        stop 'inmat not allowed for the given options'
      !!    end if
          inm => inmat
      !!else
      !!    inm => sparsemat%matrix_compr
      !!end if
    
      !call timing(iproc,'compressd_mcpy','ON')
      call f_timing(TCAT_SMAT_COMPRESSION,'ON')
    
      if (sparsemat%parallel_compression==0.or.nproc==1) then
         !call to_zero(sparsemat%nfvctr**2*sparsemat%nspin, outm(1,1,1))
         call f_zero(outm)
         do ispin=1,sparsemat%nspin
             ishift=(ispin-1)*sparsemat%nvctr
             !openmp broken on vesta
             !!! !$omp parallel default(none) private(irowcol)
             !$omp parallel default(none) private(iseg,i,ii) shared(sparsemat,inm,outm,ispin,ishift)
             !$omp do
             do iseg=1,sparsemat%nseg
                 ii=sparsemat%keyv(iseg)
                 ! a segment is always on one line, therefore no double loop
                 do i=sparsemat%keyg(1,1,iseg),sparsemat%keyg(2,1,iseg)
                    !irow = sparsemat%orb_from_index(1,ii)
                    !jcol = sparsemat%orb_from_index(2,ii)
                    !irowcol = orb_from_index(sparsemat, i)
                    outm(i,sparsemat%keyg(1,2,iseg),ispin)=inm(ii+ishift)
                    ii=ii+1
                end do
             end do
             !$omp end do
             !$omp end parallel
         end do
      else if (sparsemat%parallel_compression==1) then
         stop 'needs to be fixed'
         !!call to_zero(sparsemat%nfvctr**2, sparsemat%matrix(1,1))
         !!!$omp parallel do default(private) shared(sparsemat)
         !!do ii=1,sparsemat%nvctrp
         !!   iii=ii+sparsemat%isvctr
         !!   irow = sparsemat%orb_from_index(1,iii)
         !!   jcol = sparsemat%orb_from_index(2,iii)
         !!   sparsemat%matrix(irow,jcol)=sparsemat%matrix_compr(iii)
         !!end do
         !!!$omp end parallel do
         !!if (bigdft_mpi%nproc > 1) then
         !!   call fmpi_allreduce(sparsemat%matrix(1,1), sparsemat%nfvctr**2,FMPI_SUM,bigdft_mpi%mpi_comm)
         !!end if
      else
         stop 'needs to be fixed'
         !!sparsemat%matrixp=f_malloc_ptr((/sparsemat%nfvctr,sparsemat%nfvctrp/),id='sparsemat%matrixp')
         !!call to_zero(sparsemat%nfvctr*sparsemat%nfvctrp, sparsemat%matrixp(1,1))
         !!!$omp parallel do default(private) shared(sparsemat)
         !!do ii=1,sparsemat%nvctrp
         !!   iii=ii+sparsemat%isvctr
         !!   irow = sparsemat%orb_from_index(1,iii)
         !!   jcol = sparsemat%orb_from_index(2,iii) - sparsemat%isfvctr
         !!   sparsemat%matrixp(irow,jcol)=sparsemat%matrix_compr(iii)
         !!end do
         !!!$omp end parallel do
         !!if (bigdft_mpi%nproc > 1) &
         !!   & call mpi_allgatherv(sparsemat%matrixp, sparsemat%nfvctr*sparsemat%nfvctrp, mpi_double_precision, &
         !!   &   sparsemat%matrix, sparsemat%nfvctr*sparsemat%nfvctr_par(:), sparsemat%nfvctr*sparsemat%isfvctr_par, &
         !!   &   mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
         !!call f_free_ptr(sparsemat%matrixp)
      end if
    
      !call timing(iproc,'compressd_mcpy','OF')
      call f_timing(TCAT_SMAT_COMPRESSION,'OF')
    
    end subroutine uncompress_matrix


    subroutine uncompress_matrix2(iproc, nproc, comm, smat, matrix_compr, matrix)
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      type(sparse_matrix),intent(in) :: smat
      real(kind=mp),dimension(smat%nvctrp_tg*smat%nspin),intent(in) :: matrix_compr
      real(kind=mp),dimension(smat%nfvctr,smat%nfvctr,smat%nspin),intent(out) :: matrix

      ! Local variables
      real(kind=mp),dimension(:),allocatable :: tmparr

      tmparr = sparsematrix_malloc(smat,iaction=SPARSE_FULL,id='tmparr')
      call gather_matrix_from_taskgroups(iproc, nproc, comm, smat, matrix_compr, tmparr)
      call uncompress_matrix(iproc, nproc, smat, inmat=tmparr, outmat=matrix)
      call f_free(tmparr)
    end subroutine uncompress_matrix2


    subroutine check_matrix_compression(iproc, nproc, sparsemat, mat)
      use yaml_output
      use dynamic_memory
      implicit none
      integer,intent(in) :: iproc, nproc
      type(sparse_matrix),intent(in) :: sparsemat
      type(matrices),intent(inout) :: mat
      !Local variables
      character(len=*), parameter :: subname='check_matrix_compression'
      real(kind=mp), parameter :: tol=1.e-10
      integer :: jorb, irow, icol, iseg, ii
      real(kind=mp) :: maxdiff
    
      call f_routine('check_matrix_compression')

      call f_free_ptr(mat%matrix_compr)
      mat%matrix_compr = sparsematrix_malloc_ptr(sparsemat,iaction=SPARSE_FULL,id='mat%matrix_compr')
    
      mat%matrix = sparsematrix_malloc_ptr(sparsemat, iaction=DENSE_FULL, id='mat%matrix')
    
      !call to_zero(sparsemat%nfvctr**2*sparsemat%nspin,mat%matrix(1,1,1))
      call f_zero(mat%matrix)
      do iseg = 1, sparsemat%nseg
         ! A segment is always on one line, therefore no double loop
         do jorb = sparsemat%keyg(1,1,iseg), sparsemat%keyg(2,1,iseg)
            !call get_indices(jorb,irow,icol)
            irow = jorb
            icol = sparsemat%keyg(1,2,iseg)
            !print *,'jorb, irow,icol',jorb, irow, icol,test_value_matrix(sparsemat%nfvctr, irow, icol)
            !SM: need to fix spin 
            mat%matrix(irow,icol,1) = test_value_matrix(sparsemat%nfvctr, irow, icol)
         end do
      end do
      
      call compress_matrix(iproc, nproc, sparsemat, inmat=mat%matrix, outmat=mat%matrix_compr)
      !write(*,*) 'mat%matrix',mat%matrix
      !write(*,*) 'mat%matrix_compr',mat%matrix_compr
    
      maxdiff = 0.d0
      do iseg = 1, sparsemat%nseg
         ii=0
         ! A segment is always on one line, therefore no double loop
         do jorb = sparsemat%keyg(1,1,iseg), sparsemat%keyg(2,1,iseg)
            !call get_indices(jorb,irow,icol)
            irow = jorb
            icol = sparsemat%keyg(1,2,iseg)
            !write(*,'(a,4i8,2es13.3)') 'jorb, irow, icol, sparsemat%keyv(iseg)+ii, val, ref', jorb, irow, icol, sparsemat%keyv(iseg)+ii, mat%matrix_compr(sparsemat%keyv(iseg)+ii), test_value_matrix(sparsemat%nfvctr, irow, icol)
            maxdiff = max(abs(mat%matrix_compr(sparsemat%keyv(iseg)+ii)&
                 -test_value_matrix(sparsemat%nfvctr, irow, icol)),maxdiff)
            ii=ii+1
         end do
      end do
    
      if (iproc==0) call yaml_map('Tolerances for this check',tol,fmt='(1pe25.17)')
    
      if(iproc==0) then
        if (maxdiff > tol) then
           call yaml_warning('COMPRESSION ERROR : difference of '//trim(yaml_toa(maxdiff,fmt='(1pe12.5)')))
        else
           call yaml_map('Maxdiff for compress', maxdiff,fmt='(1pe25.17)')
        end if
      end if
    
      call uncompress_matrix(iproc, nproc, sparsemat, inmat=mat%matrix_compr, outmat=mat%matrix)
    
      maxdiff = 0.d0
      do iseg = 1, sparsemat%nseg
         ! A segment is always on one line, therefore no double loop
         do jorb = sparsemat%keyg(1,1,iseg), sparsemat%keyg(2,1,iseg)
            !call get_indices(jorb,irow,icol)
            irow = jorb
            icol = sparsemat%keyg(1,2,iseg)
            maxdiff = max(abs(mat%matrix(irow,icol,1)-test_value_matrix(sparsemat%nfvctr, irow, icol)),maxdiff) 
         end do
      end do
    
      if(iproc==0) then
        if (maxdiff > tol) then
           call yaml_warning('UNCOMPRESSION ERROR : difference of '//trim(yaml_toa(maxdiff,fmt='(1pe12.5)')))
        else
           call yaml_map('Maxdiff for uncompress', maxdiff,fmt='(1pe25.17)')
        end if
      end if
    
      call f_free_ptr(mat%matrix)
      !!call f_free_ptr(sparsemat%matrix_compr)

      call f_free_ptr(mat%matrix_compr)
      mat%matrix_compr = sparsematrix_malloc_ptr(sparsemat,iaction=SPARSE_TASKGROUP,id='mat%matrix_compr')

      call f_release_routine()
    
    contains
       !> define a value for the wavefunction which is dependent of the indices
       function test_value_matrix(norb,iorb,jorb)
          implicit none
          integer, intent(in) :: norb,iorb,jorb
          real(kind=mp) :: test_value_matrix
    
          test_value_matrix = norb*(iorb-1)+jorb
          !print *,iorb,jorb,test_value_matrix
       END FUNCTION test_value_matrix
    
       subroutine get_indices(ind,irow,icol)
         implicit none
         integer, intent(in) :: ind
         integer, intent(out) :: irow, icol
    
         icol = (ind - 1) / sparsemat%nfvctr + 1
         irow = ind - (icol-1)*sparsemat%nfvctr
         !print *,'irow,icol',irow,icol
       END SUBROUTINE get_indices 
    end subroutine check_matrix_compression


    subroutine transform_sparse_matrix(iproc, smat, lmat, imode, direction, &
         smat_in, lmat_in, smat_out, lmat_out)
      use dynamic_memory
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, imode
      type(sparse_matrix),intent(in) :: smat, lmat
      character(len=14),intent(in) :: direction
      real(kind=mp),dimension(:),intent(in),optional :: smat_in
      real(kind=mp),dimension(:),intent(in),optional :: lmat_in
      real(kind=mp),dimension(:),intent(out),optional :: smat_out
      real(kind=mp),dimension(:),intent(out),optional :: lmat_out
    
      ! Local variables
      integer(kind=mp) :: isstart, isend, ilstart, ilend, iostart, ioend
      integer :: idir, icheck, isseg, ilseg, isoffset_tg, iloffset_tg, issegstartx, issegendx
      integer :: ilength, iscostart, ilcostart, nssize, nlsize, i, ilsegstartx, ilsegendx
      integer :: ilsegstart, ispin, isshift, ilshift, isoffset, iloffset
      integer,parameter :: SMALL_TO_LARGE=1
      integer,parameter :: LARGE_TO_SMALL=2
    
      call f_routine(id='transform_sparse_matrix')

      ! determine whether the routine should handle a matrix taskgroup or the full matrices
      if (imode==SPARSE_FULL) then
          ilsegstartx = 1
          issegstartx = 1
          ilsegendx = lmat%nseg
          issegendx = smat%nseg
          isoffset_tg = 0
          iloffset_tg = 0
          nssize = smat%nvctr*smat%nspin
          nlsize = lmat%nvctr*lmat%nspin
      else if (imode==SPARSE_TASKGROUP) then
          ilsegstartx = lmat%iseseg_tg(1)
          issegstartx = smat%iseseg_tg(1)
          ilsegendx = lmat%iseseg_tg(2)
          issegendx = smat%iseseg_tg(2)
          isoffset_tg = smat%isvctrp_tg
          iloffset_tg = lmat%isvctrp_tg
          nssize = smat%nvctrp_tg*smat%nspin
          nlsize = lmat%nvctrp_tg*lmat%nspin
      else if (imode==SPARSE_MATMUL_SMALL) then
          if (direction/='small_to_large' .and. direction/='SMALL_TO_LARGE') then
              call f_err_throw('SPARSE_MATMUL_SMALL only posible for small_to_large')
          end if
          if (.not.lmat%smatmul_initialized) then
              call f_err_throw('lmat%smatmul_initialized mus be true')
          end if
          ilsegstartx = lmat%smmm%istartendseg_mm(1)
          issegstartx = smat%iseseg_tg(1)
          ilsegendx = lmat%smmm%istartendseg_mm(2)
          issegendx = smat%iseseg_tg(2)
          isoffset_tg = smat%isvctrp_tg
          iloffset_tg = lmat%isvctrp_tg
          nssize = smat%nvctrp_tg*smat%nspin
          nlsize = lmat%nvctrp_tg*lmat%nspin
      else
          call f_err_throw('wrong imode')
      end if
    
      ! determine the case:
      ! SMALL_TO_LARGE -> transform from large sparsity pattern to small one
      ! LARGE_TO_SMALL -> transform from small sparsity pattern to large one
      if (direction=='small_to_large' .or. direction=='SMALL_TO_LARGE') then
          idir=SMALL_TO_LARGE
          if (.not.present(smat_in)) call f_err_throw('smat_in not present')
          if (.not.present(lmat_out)) call f_err_throw('lmat_out not present')
          if (size(smat_in)/=nssize) then
              call f_err_throw('the size of smat_in ('//trim(yaml_toa(size(smat_in)))//&
                   &') is not the correct one ('//trim(yaml_toa(nssize))//')')
          end if
          if (size(lmat_out)/=nlsize) then
              call f_err_throw('the size of lmat_out ('//trim(yaml_toa(size(lmat_out)))//&
                   &') is not the correct one ('//trim(yaml_toa(nlsize))//')')
          end if
      else if (direction=='large_to_small' .or. direction=='LARGE_TO_SMALL') then
          idir=LARGE_TO_SMALL
          if (.not.present(lmat_in)) call f_err_throw('lmat_in not present')
          if (.not.present(smat_out)) call f_err_throw('smat_out not present')
          if (size(lmat_in)/=nlsize) then
              call f_err_throw('the size of lmat_in ('//trim(yaml_toa(size(lmat_in)))//&
                   &') is not the correct one ('//trim(yaml_toa(nlsize))//')')
          end if
          if (size(smat_out)/=nssize) then
              call f_err_throw('the size of smat_out ('//trim(yaml_toa(size(smat_out)))//&
                   &') is not the correct one ('//trim(yaml_toa(nssize))//')')
          end if
      else
          call f_err_throw('wrong direction')
      end if
    
      select case (idir)
      case (SMALL_TO_LARGE)
         call f_zero(lmat_out)
      case (LARGE_TO_SMALL)
         call f_zero(smat_out)
      case default
          call f_err_throw('wrong idir')
      end select
    
      !call timing(iproc,'transform_matr','IR')
      call f_timing(TCAT_SMAT_TRANSFORMATION,'IR')


      icheck=0
      do ispin=1,smat%nspin

          isshift=(ispin-1)*smat%nvctr
          ilshift=(ispin-1)*lmat%nvctr
    
          ilsegstart = ilsegstartx
          !$omp parallel default(none) &
          !$omp shared(smat, lmat, issegstartx, issegendx, idir, icheck, isshift, ilshift, isoffset_tg, iloffset_tg) &
          !$omp shared(smat_in, lmat_in, smat_out, lmat_out, ilsegendx) &
          !$omp private(isseg, isstart, isend, ilstart, ilend, iostart, ioend) &
          !$omp private(isoffset, iloffset, iscostart, ilcostart, ilength,ilseg,i) &
          !$omp firstprivate(ilsegstart)
          !$omp do reduction(+:icheck)
          sloop: do isseg=issegstartx,issegendx
              isstart = int((smat%keyg(1,2,isseg)-1),kind=mp)*int(smat%nfvctr,kind=mp) + int(smat%keyg(1,1,isseg),kind=mp)
              isend = int((smat%keyg(2,2,isseg)-1),kind=mp)*int(smat%nfvctr,kind=mp) + int(smat%keyg(2,1,isseg),kind=mp)
              ! A segment is always on one line, therefore no double loop
              lloop: do ilseg=ilsegstart,ilsegendx
                  ilstart = int((lmat%keyg(1,2,ilseg)-1),kind=mp)*int(lmat%nfvctr,kind=mp) + int(lmat%keyg(1,1,ilseg),kind=mp)
                  ilend = int((lmat%keyg(2,2,ilseg)-1),kind=mp)*int(lmat%nfvctr,kind=mp) + int(lmat%keyg(2,1,ilseg),kind=mp)
    
                  ! check whether there is an overlap:
                  ! if not, increase loop counters
                  if (ilstart>isend) then
                      !ilsegstart=ilseg
                      exit lloop
                  end if
                  if (isstart>ilend) then
                      ilsegstart=ilseg
                      cycle lloop
                  end if
                  ! if yes, determine start end end of overlapping segment (in uncompressed form)
                  iostart=max(isstart,ilstart)
                  ioend=min(isend,ilend)
                  ilength=int(ioend-iostart+1,kind=4)
    
                  ! offset with respect to the starting point of the segment
                  isoffset = int(iostart - &
                             (int((smat%keyg(1,2,isseg)-1),kind=mp)*int(smat%nfvctr,kind=mp) &
                               + int(smat%keyg(1,1,isseg),kind=mp)),kind=4)
                  iloffset = int(iostart - &
                             (int((lmat%keyg(1,2,ilseg)-1),kind=mp)*int(lmat%nfvctr,kind=mp) &
                               + int(lmat%keyg(1,1,ilseg),kind=mp)),kind=4)
    
                  ! determine start end and of the overlapping segment in compressed form
                  iscostart=smat%keyv(isseg)+isoffset
                  ilcostart=lmat%keyv(ilseg)+iloffset
    
                  ! copy the elements
                  select case (idir)
                  case (SMALL_TO_LARGE) 
                      do i=0,ilength-1
                          lmat_out(ilcostart+i+ilshift-iloffset_tg)=smat_in(iscostart+i+isshift-isoffset_tg)
                      end do
                  case (LARGE_TO_SMALL) 
                      do i=0,ilength-1
                          smat_out(iscostart+i+isshift-isoffset_tg)=lmat_in(ilcostart+i+ilshift-iloffset_tg)
                      end do
                  case default
                      stop 'wrong idir'
                  end select
                  icheck=icheck+ilength
              end do lloop
          end do sloop
          !$omp end do 
          !$omp end parallel

      end do
    
      if (imode==SPARSE_MATMUL_SMALL) then
          !!if (icheck/=lmat%smmm%nvctrp_mm) then
          !!    call f_err_throw(trim(yaml_toa(icheck))//'=icheck /= lmat%smmm%nvctrp_mm='//trim(yaml_toa(lmat%smmm%nvctrp_mm)))
          !!end if
      else
          ! all elements of the small matrix must have been processed, no matter in
          ! which direction the transformation has been executed
          !if (icheck/=smat%nvctr*smat%nspin) then
          if (icheck/=nssize) then
              call f_err_throw(trim(yaml_toa(icheck))//'=icheck /= nssize='//trim(yaml_toa(nssize)))
          end if
      end if

      !call timing(iproc,'transform_matr','RS')
      call f_timing(TCAT_SMAT_TRANSFORMATION,'RS')
      call f_release_routine()
    
    end subroutine transform_sparse_matrix



!!    subroutine transform_sparse_matrix_local(iproc, smat, lmat, cmode, &
!!               smatrix_compr_in, lmatrix_compr_in, smatrix_compr_out, lmatrix_compr_out)
!!      implicit none
!!    
!!      ! Calling arguments
!!      integer,intent(in) :: iproc
!!      type(sparse_matrix),intent(in) :: smat, lmat
!!      character(len=14),intent(in) :: cmode
!!      real(kind=8),dimension(smat%nspin*smat%nvctrp_tg),intent(in),optional :: smatrix_compr_in
!!      real(kind=8),dimension(lmat%nspin*lmat%nvctrp_tg),intent(in),optional :: lmatrix_compr_in
!!      real(kind=8),dimension(smat%nspin*smat%nvctrp_tg),intent(out),optional :: smatrix_compr_out
!!      real(kind=8),dimension(lmat%nspin*lmat%nvctrp_tg),intent(out),optional :: lmatrix_compr_out
!!    
!!      ! Local variables
!!      real(kind=8),dimension(:),allocatable :: tmparrs, tmparrl
!!      integer :: ishift_src, ishift_dst, imode, ispin, i
!!      integer,parameter :: SMALL_TO_LARGE=1
!!      integer,parameter :: LARGE_TO_SMALL=2
!!    
!!      call f_routine(id='transform_sparse_matrix_local')
!!
!!
!!      ! determine the case:
!!      ! SMALL_TO_LARGE -> transform from large sparsity pattern to small one
!!      ! LARGE_TO_SMALL -> transform from small sparsity pattern to large one
!!      if (cmode=='small_to_large' .or. cmode=='SMALL_TO_LARGE') then
!!          imode=SMALL_TO_LARGE
!!      else if (cmode=='large_to_small' .or. cmode=='LARGE_TO_SMALL') then
!!          imode=LARGE_TO_SMALL
!!      else
!!          stop 'wrong cmode'
!!      end if
!!
!!    
!!      select case (imode)
!!      case (SMALL_TO_LARGE)
!!          ! check the aruments
!!          if (.not.present(smatrix_compr_in)) call f_err_throw("'smatrix_compr_in' not present")
!!          if (.not.present(lmatrix_compr_out)) call f_err_throw("'lmatrix_compr_out' not present")
!!          tmparrs = sparsematrix_malloc0(smat,iaction=SPARSE_FULL,id='tmparrs')
!!          tmparrl = sparsematrix_malloc(lmat,iaction=SPARSE_FULL,id='tmparrl')
!!          call transform_sparse_matrix_test(iproc, smat, lmat, cmode, &
!!               smat_in=smatrix_compr_in, lmat_out=lmatrix_compr_out)
!!          do i=1,lmat%nspin*lmat%nvctrp_tg
!!              write(1000+iproc,*) 'i, val', i, lmatrix_compr_out(i)
!!          end do
!!          do ispin=1,smat%nspin
!!              ishift_src = (ispin-1)*smat%nvctrp_tg
!!              ishift_dst = (ispin-1)*smat%nvctr
!!              call vcopy(smat%nvctrp_tg, smatrix_compr_in(ishift_src+1), 1, &
!!                   tmparrs(ishift_dst+smat%isvctrp_tg+1), 1)
!!              call transform_sparse_matrix(iproc, smat, lmat, cmode, smat_in=tmparrs, lmat_out=tmparrl)
!!              call extract_taskgroup(lmat, tmparrl, lmatrix_compr_out)
!!          end do
!!          do i=1,lmat%nspin*lmat%nvctrp_tg
!!              write(1100+iproc,*) 'i, val', i, lmatrix_compr_out(i)
!!          end do
!!      case (LARGE_TO_SMALL)
!!          if (.not.present(lmatrix_compr_in)) call f_err_throw("'lmatrix_compr_in' not present")
!!          if (.not.present(smatrix_compr_out)) call f_err_throw("'smatrix_compr_out' not present")
!!          tmparrs = sparsematrix_malloc(smat,iaction=SPARSE_FULL,id='tmparrs')
!!          tmparrl = sparsematrix_malloc0(lmat,iaction=SPARSE_FULL,id='tmparrl')
!!          call transform_sparse_matrix_test(iproc, smat, lmat, cmode, &
!!               lmat_in=lmatrix_compr_in, smat_out=smatrix_compr_out)
!!          do i=1,smat%nspin*smat%nvctrp_tg
!!              write(2000+iproc,*) 'i, val', i, smatrix_compr_out(i)
!!          end do
!!          do ispin=1,smat%nspin
!!              ishift_src = (ispin-1)*lmat%nvctrp_tg
!!              ishift_dst = (ispin-1)*lmat%nvctr
!!              call vcopy(lmat%nvctrp_tg, lmatrix_compr_in(ishift_src+1), 1, &
!!                   tmparrl(ishift_dst+lmat%isvctrp_tg+1), 1)
!!              call transform_sparse_matrix(iproc, smat, lmat, cmode, lmat_in=tmparrl, smat_out=tmparrs)
!!              call extract_taskgroup(smat, tmparrs, smatrix_compr_out)
!!          end do
!!          do i=1,smat%nspin*smat%nvctrp_tg
!!              write(2100+iproc,*) 'i, val', i, smatrix_compr_out(i)
!!          end do
!!      case default
!!          stop 'wrong imode'
!!      end select
!!
!!      call f_free(tmparrs)
!!      call f_free(tmparrl)
!!
!!      call f_release_routine()
!!    
!!  end subroutine transform_sparse_matrix_local




   subroutine compress_matrix_distributed_wrapper_1(iproc, nproc, smat, layout, matrixp, &
              onesided_action, matrix_compr, matrix_localx, windowsx)
     use dynamic_memory
     implicit none

     ! Calling arguments
     integer,intent(in) :: iproc, nproc, layout, onesided_action
     type(sparse_matrix),intent(in) :: smat
     real(kind=mp),dimension(:),target,intent(in) :: matrixp
     real(kind=mp),dimension(smat%nvctrp_tg),intent(out) :: matrix_compr
     real(kind=mp),dimension(:),intent(inout),target,optional :: matrix_localx
     type(fmpi_win),dimension(:),intent(inout),optional :: windowsx

     ! Local variables
     !integer :: nvctrp
     !integer :: ncount
     !integer :: window, sizeof
     !integer,dimension(:),pointer :: nvctr_par
     !integer,dimension(:),allocatable :: request, windows
     real(kind=mp),dimension(:),pointer :: matrix_local
     !real(kind=mp),dimension(:),allocatable :: recvbuf

     call f_routine(id='compress_matrix_distributed_wrapper_1')

     !call timing(iproc,'compressd_mcpy','ON')


     if (layout==SPARSE_MATMUL_SMALL) then
         if (.not.smat%smatmul_initialized) then
             call f_err_throw('sparse matrix multiplication not initialized', &
                  err_name='SPARSEMATRIX_RUNTIME_ERROR')
         end if
         if (size(matrixp)/=smat%smmm%nvctrp_mm) then
             call f_err_throw('Array matrixp has size '//trim(yaml_toa(size(matrixp),fmt='(i0)'))//&
                  &' instead of '//trim(yaml_toa(smat%smmm%nvctrp_mm,fmt='(i0)')), &
                  err_name='SPARSEMATRIX_MANIPULATION_ERROR')
         end if
         matrix_local => matrixp
     else if (layout==SPARSE_MATMUL_LARGE) then
         if (.not.smat%smatmul_initialized) then
             call f_err_throw('sparse matrix multiplication not initialized', &
                  err_name='SPARSEMATRIX_RUNTIME_ERROR')
         end if
         if (size(matrixp)/=smat%smmm%nvctrp) then
             call f_err_throw('Array matrixp has size '//trim(yaml_toa(size(matrixp),fmt='(i0)'))//&
                  &' instead of '//trim(yaml_toa(smat%smmm%nvctrp,fmt='(i0)')), &
                  err_name='SPARSEMATRIX_MANIPULATION_ERROR')
         end if
         select case (onesided_action)
         case (ONESIDED_POST,ONESIDED_GATHER)
             if (.not.present(matrix_localx)) then
                 call f_err_throw('matrix_localx not present')
             end if
             if (size(matrix_localx)/=smat%smmm%nvctrp_mm) then
                 call f_err_throw('Array matrix_localx has size '//trim(yaml_toa(size(matrix_localx),fmt='(i0)'))//&
                      &' instead of '//trim(yaml_toa(smat%smmm%nvctrp_mm,fmt='(i0)')), &
                      err_name='SPARSEMATRIX_MANIPULATION_ERROR')
             end if
             matrix_local => matrix_localx
         case (ONESIDED_FULL)
             matrix_local = f_malloc_ptr(smat%smmm%nvctrp_mm,id='matrix_local')
         case default
             call f_err_throw('wrong value for onesided_action')
         end select
         if (onesided_action==ONESIDED_POST .or. onesided_action==ONESIDED_FULL) then
             call transform_sparsity_pattern(iproc, smat%nfvctr, smat%smmm%nvctrp_mm, smat%smmm%isvctr_mm, &
                  smat%nseg, smat%keyv, smat%keyg, smat%smmm%line_and_column_mm, &
                  smat%smmm%nvctrp, smat%smmm%isvctr, &
                  smat%smmm%nseg, smat%smmm%keyv, smat%smmm%keyg, smat%smmm%istsegline, &
                  'large_to_small', matrix_s_out=matrix_local, matrix_l_in=matrixp)
         end if
         !call f_free_ptr(matrix_local)
     else
             call f_err_throw('layout has the value '//trim(yaml_toa(layout,fmt='(i0)'))//&
                  &'; allowed are '//trim(yaml_toa(SPARSE_MATMUL_SMALL,fmt='(i0)'))//&
                  &' and '//trim(yaml_toa(SPARSE_MATMUL_LARGE,fmt='(i0)')), &
                  err_name='SPARSEMATRIX_MANIPULATION_ERROR')
     end if

     !not needed as the interfaces are explicit
!!$     if (present(windowsx)) then
!!$         call compress_matrix_distributed_core(iproc, nproc, smat, SPARSE_MATMUL_SMALL, &
!!$              matrix_local, onesided_action, matrix_compr, windowsx)
!!$     else
!!$         call compress_matrix_distributed_core(iproc, nproc, smat, SPARSE_MATMUL_SMALL, &
!!$              matrix_local, onesided_action, matrix_compr)
!!$     end if
     call compress_matrix_distributed_core(iproc, nproc, smat, SPARSE_MATMUL_SMALL, &
          matrix_local, onesided_action, matrix_compr, windowsx=windowsx)

     if (layout==SPARSE_MATMUL_LARGE) then
         if (onesided_action==ONESIDED_FULL) then
             call f_free_ptr(matrix_local)
         end if
     end if

     !!call timing(iproc,'compressd_comm_new','OF')

     call f_release_routine()

  end subroutine compress_matrix_distributed_wrapper_1


   subroutine compress_matrix_distributed_wrapper_2(iproc, nproc, smat, layout, matrixp, &
              onesided_action, matrix_compr, matrix_localx, windowsx)
     !!use yaml_output
     use dynamic_memory
     implicit none

     ! Calling arguments
     integer,intent(in) :: iproc, nproc, layout, onesided_action
     type(sparse_matrix),intent(in) :: smat
     real(kind=mp),dimension(:,:),intent(in) :: matrixp
     real(kind=mp),dimension(smat%nvctrp_tg),intent(out) :: matrix_compr
     real(kind=mp),dimension(:),intent(inout),target,optional :: matrix_localx
     type(fmpi_win),dimension(:),intent(inout),optional :: windowsx

     ! Local variables
     integer :: isegstart, isegend, iseg, ii, jorb, nfvctrp, isfvctr, nvctrp, isvctr
     !integer :: ncount
     integer :: iorb
     !integer :: window
     integer,dimension(:),pointer :: isvctr_par, nvctr_par
     real(kind=mp),dimension(:),pointer :: matrix_local
     !real(kind=mp),dimension(:),allocatable :: recvbuf

     call f_routine(id='compress_matrix_distributed_wrapper_2')

     !call timing(iproc,'compressd_mcpy','ON')
     call f_timing(TCAT_SMAT_COMPRESSION,'ON')

     ! Check the dimensions of the input array and assign some values
     !if (size(matrixp,1)/=smat%nfvctr) stop 'size(matrixp,1)/=smat%nfvctr'
     if (size(matrixp,1)/=smat%nfvctr) then
         call f_err_throw('Array matrixp has size '//trim(yaml_toa(size(matrixp,1),fmt='(i0)'))//&
              &' instead of '//trim(yaml_toa(smat%nfvctr,fmt='(i0)')), &
              err_name='SPARSEMATRIX_MANIPULATION_ERROR')
     end if
     if (layout==DENSE_PARALLEL) then
         if (size(matrixp,2)/=smat%nfvctrp) stop '(ubound(matrixp,2)/=smat%nfvctrp'
         nfvctrp = smat%nfvctrp
         isfvctr = smat%isfvctr
         nvctrp = smat%nvctrp
         isvctr = smat%isvctr
         isvctr_par => smat%isvctr_par
         nvctr_par => smat%nvctr_par
     else if (layout==DENSE_MATMUL) then
         if (.not.smat%smatmul_initialized) then
             call f_err_throw('sparse matrix multiplication not initialized', &
                  err_name='SPARSEMATRIX_RUNTIME_ERROR')
         end if
         if (size(matrixp,2)/=smat%smmm%nfvctrp) stop '(ubound(matrixp,2)/=smat%smmm%nfvctrp'
         nfvctrp = smat%smmm%nfvctrp
         isfvctr = smat%smmm%isfvctr
         nvctrp = smat%smmm%nvctrp_mm
         isvctr = smat%smmm%isvctr_mm
         isvctr_par => smat%smmm%isvctr_mm_par
         nvctr_par => smat%smmm%nvctr_mm_par
     else
         call f_err_throw('layout has the value '//trim(yaml_toa(layout,fmt='(i0)'))//&
              &'; allowed are '//trim(yaml_toa(DENSE_PARALLEL,fmt='(i0)'))//&
              &' and '//trim(yaml_toa(DENSE_MATMUL,fmt='(i0)')), &
              err_name='SPARSEMATRIX_MANIPULATION_ERROR')
     end if



     !@ NEW #####################
     matrix_local = f_malloc_ptr(nvctrp,id='matrix_local')
     if (layout==DENSE_PARALLEL) then
         ii = 0
         if (nfvctrp>0) then
             isegstart=smat%istsegline(isfvctr+1)
             isegend=smat%istsegline(isfvctr+nfvctrp)+smat%nsegline(isfvctr+nfvctrp)-1
             do iseg=isegstart,isegend
                 ! A segment is always on one line, therefore no double loop
                 do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                     iorb = smat%keyg(1,2,iseg)
                     ii = ii + 1
                     matrix_local(ii) = matrixp(jorb,iorb-isfvctr)
                 end do
             end do
         end if
         if (ii/=nvctrp) stop 'compress_matrix_distributed: ii/=nvctrp'
     else if (layout==DENSE_MATMUL) then
         ii = 0
         if (nvctrp>0) then
             isegstart=smat%istsegline(isfvctr+1)
             isegend=smat%istsegline(isfvctr+nfvctrp)+smat%nsegline(isfvctr+nfvctrp)-1
             do iseg=isegstart,isegend
                 ! A segment is always on one line, therefore no double loop
                 do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                     iorb = smat%keyg(1,2,iseg)
                     ii = ii + 1
                     matrix_local(ii) = matrixp(jorb,iorb-isfvctr)
                 end do
             end do
         end if
         if (ii/=nvctrp) stop 'compress_matrix_distributed: ii/=nvctrp'
     else
         stop 'compress_matrix_distributed: wrong data_strategy'
     end if

     !call timing(iproc,'compressd_mcpy','OF')
     call f_timing(TCAT_SMAT_COMPRESSION,'OF')

     if (onesided_action/=ONESIDED_FULL) then
         ! For the other cases I have to pass a proper workarray matrix_local
         call f_err_throw('compress_matrix_distributed_wrapper_2 not yet functional for onesided_action/=ONESIDED_FULL')
     end if

     !if (present(windowsx)) then
     !no need to put the if for fortran norm
     call compress_matrix_distributed_core(iproc, nproc, smat, SPARSE_PARALLEL, &
              matrix_local, onesided_action, matrix_compr, windowsx)
     !else
     !    call compress_matrix_distributed_core(iproc, nproc, smat, SPARSE_PARALLEL, &
     !         matrix_local, onesided_action, matrix_compr)
     !end if
     call f_free_ptr(matrix_local)
     !@ END NEW #################

     call f_release_routine()

  end subroutine compress_matrix_distributed_wrapper_2



   !> Gathers together the matrix parts calculated by other tasks.
   !! Attention: Even if the output array has size smat%nvctrp_tg, only the
   !! relevant part (indicated by smat%istartend_local) is calculated
  subroutine compress_matrix_distributed_core(iproc, nproc, smat, layout, matrixp, &
             onesided_action, matrix_compr, windowsx)
    use dynamic_memory
     implicit none

     ! Calling arguments
     integer,intent(in) :: iproc, nproc, layout, onesided_action
     type(sparse_matrix),intent(in) :: smat
     !real(kind=mp),dimension(smat%smmm%nvctrp_mm),intent(in) :: matrixp
     real(kind=mp),dimension(:),intent(in) :: matrixp
     real(kind=mp),dimension(smat%nvctrp_tg),target,intent(out) :: matrix_compr
     type(fmpi_win),dimension(:),intent(inout),target,optional :: windowsx

     ! Local variables
     integer :: ii, nvctrp
     integer :: ncount, itg, iitg, ist_send, ist_recv
     integer :: jproc_send, jproc, nccomm
     real(kind=mp) :: window_fake
     !integer,dimension(:),pointer :: nvctr_par
     type(fmpi_win),dimension(:),pointer :: windows
     !real(kind=mp),dimension(:),pointer :: matrix_local
     !real(kind=mp),dimension(:),allocatable :: recvbuf
     integer,parameter :: ALLGATHERV=51, GET=52, GLOBAL_MATRIX=101, SUBMATRIX=102
     integer,parameter :: comm_strategy=GET
     integer,parameter :: data_strategy=SUBMATRIX!GLOBAL_MATRIX
     integer,dimension(:,:),pointer :: luccomm

     call f_routine(id='compress_matrix_distributed_core')

     !call timing(iproc,'compressd_mcpy','ON')
     call f_timing(TCAT_SMAT_COMPRESSION,'ON')

     ! Check the arguments
     select case (onesided_action)
     case (ONESIDED_POST,ONESIDED_GATHER)
         if (nproc>1) then
            if (.not.present(windowsx)) call f_err_throw('windowsx not present')
             if (size(windowsx)/=smat%ntaskgroup) then
                 call f_err_throw('size(windowsx)='//trim(yaml_toa(size(windowsx))) //&
                      &' /= smat%ntaskgroup='//trim(yaml_toa(smat%ntaskgroup)))
             end if
             !!!!write(*,*) 'windwos is present, iproc', iproc
             windows => windowsx
         end if
     case (ONESIDED_FULL)
         if (nproc>1) then
             ! Create a window for all taskgroups to which iproc belongs (max 2)
             windows = f_malloc_ptr(smat%ntaskgroup,id='windows')
         end if
     case default
         call f_err_throw('wrong value for onesided_action')
     end select


     if (data_strategy==GLOBAL_MATRIX) then
         stop 'compress_matrix_distributed: option GLOBAL_MATRIX is deprecated'
     else if (data_strategy==SUBMATRIX) then
         !if (layout==DENSE_PARALLEL) then
         if (layout==SPARSE_PARALLEL) then
             !stop 'layout==DENSE_PARALLEL not yet implemented'
             luccomm => smat%luccomm
             nvctrp = smat%nvctrp
             nccomm = smat%nccomm
         !else if (layout==DENSE_MATMUL) then
         else if (layout==SPARSE_MATMUL_SMALL) then
             if (.not.smat%smatmul_initialized) then
                 call f_err_throw('sparse matrix multiplication not initialized', &
                      err_name='SPARSEMATRIX_RUNTIME_ERROR')
             end if
             luccomm => smat%smmm%luccomm_smmm
             nvctrp = smat%smmm%nvctrp_mm
             nccomm = smat%smmm%nccomm_smmm
         end if
         if (size(matrixp)/=nvctrp) then
             call f_err_throw('Array matrixp has size '//trim(yaml_toa(size(matrixp),fmt='(i0)'))//&
                  &' instead of '//trim(yaml_toa(nvctrp,fmt='(i0)')), &
                  err_name='SPARSEMATRIX_MANIPULATION_ERROR')
         end if

             !call timing(iproc,'compressd_mcpy','OF')
             call f_timing(TCAT_SMAT_COMPRESSION,'OF')
             !call timing(iproc,'compressd_comm','ON')
             call f_timing(TCAT_SMAT_COMPRESSION_COMMUNICATION,'ON')

             if (onesided_action==ONESIDED_POST .or. onesided_action==ONESIDED_FULL) then
                 if (nproc>1) then
                    call f_zero(matrix_compr)

                     !!! Create a window for all taskgroups to which iproc belongs (max 2)
                     !!windows = f_malloc(smat%ntaskgroup,id='windows')
                     do itg=1,smat%ntaskgroupp
                         iitg = smat%taskgroupid(itg)
                         ! Use a fake window if nvctrp is zero
                         if (nvctrp>0) then
                            !windows(iitg) = mpiwindow(nvctrp, matrixp(1), smat%mpi_groups(iitg)%mpi_comm)
                            call fmpi_win_create(windows(iitg),matrixp(1),&
                                 int(nvctrp,f_long),&
                                 comm=smat%mpi_groups(iitg)%mpi_comm)
                            call fmpi_win_fence(windows(iitg),FMPI_WIN_OPEN)
                         else
                            !windows(iitg) = mpiwindow(1, window_fake, smat%mpi_groups(iitg)%mpi_comm)
                            call fmpi_win_create(windows(iitg),window_fake,&
                                 int(1,f_long),comm=smat%mpi_groups(iitg)%mpi_comm)
                            call fmpi_win_fence(windows(iitg),FMPI_WIN_OPEN)
                            !no fence here
                         end if
                     end do
                     do jproc=1,nccomm
                         jproc_send = luccomm(1,jproc)
                         ist_send = luccomm(2,jproc)
                         ist_recv = luccomm(3,jproc)
                         ncount = luccomm(4,jproc)
                         !write(*,'(5(a,i0))') 'task ',iproc,' gets ',ncount,' elements at position ',ist_recv,' from position ',ist_send,' on task ',jproc_send
                         iitg = get_taskgroup_id(iproc,jproc_send)
                         ! Now get the task ID on the taskgroup (subtract the ID of the first task)
                         !jproc_send = jproc_send - smat%isrank(iitg)
                         ii = jproc_send
                         jproc_send = get_rank_on_taskgroup(ii,iitg)
                         !call mpiget(matrix_compr(ist_recv), ncount, jproc_send, int(ist_send-1,kind=mpi_address_kind), window)
                         !!write(*,'(6(a,i0))') 'task ',iproc,' gets ',ncount,' elements at position ',ist_recv, &
                         !!                     ' from task ',jproc_send,' with offset ',ist_send-1,' on window ',iitg
                         call fmpi_get(matrix_compr(ist_recv),jproc_send,&
                              windows(iitg),ncount,int(ist_send-1,fmpi_address))
!!$                         call mpiget(matrix_compr(ist_recv), ncount, jproc_send, &
!!$                              int(ist_send-1,kind=mpi_address_kind), windows(iitg))
                     end do
                 else
                     ist_send = luccomm(2,1)
                     ist_recv = luccomm(3,1)
                     ncount = luccomm(4,1)
                     call vcopy(ncount, matrixp(ist_send), 1, matrix_compr(ist_recv), 1)
                 end if
             end if

             if (onesided_action==ONESIDED_GATHER .or. onesided_action==ONESIDED_FULL) then
                 if (nproc>1) then
                     ! Synchronize the communication
                     do itg=1,smat%ntaskgroupp
                         iitg = smat%taskgroupid(itg)
                         !!write(*,'(a,i0,a)') 'task ',iproc,' calls fence'
                         !call mpi_fenceandfree(windows(iitg))
                         call fmpi_win_fence(windows(iitg),FMPI_WIN_CLOSE)
                         call fmpi_win_free(windows(iitg))
                         !LG: why closing the windows (that is what fenceanfree does) if we have windowsx?
                     end do
                     if (onesided_action==ONESIDED_FULL) then
                         call free_fmpi_win_ptr(windows)
                     end if
                     !call mpi_fenceandfree(window)
                 end if
             end if

             !!call f_free_ptr(matrix_local)

         !!end if

         !call timing(iproc,'compressd_comm','OF')
         call f_timing(TCAT_SMAT_COMPRESSION_COMMUNICATION,'OF')
         call f_timing(TCAT_SMAT_COMPRESSION,'ON')
     else
         stop 'compress_matrix_distributed: wrong data_strategy'
     end if

     !!!call timing(iproc,'compressd_comm','OF')
     call f_timing(TCAT_SMAT_COMPRESSION,'OF')

     call f_release_routine()


     contains



       !> Get the taskgroup which should be used for the communication, i.e. the
       !! one to which both iproc and jproc belong
       integer function get_taskgroup_id(iproc,jproc)
         implicit none
         integer,intent(in) :: iproc, jproc

         ! Local variables
         integer :: itg, iitg, jtg, jjtg
         logical :: found

         ! A task never belongs to more than 2 taskgroups 
         found = .false.
         iloop: do itg=1,2
             iitg = smat%inwhichtaskgroup(itg,iproc)
             jloop: do jtg=1,2
                 jjtg = smat%inwhichtaskgroup(jtg,jproc)
                 if (iitg==jjtg) then
                     get_taskgroup_id = iitg
                     found = .true.
                     exit iloop
                 end if
             end do jloop
         end do iloop
         if (.not.found) then
             call f_err_throw('get_taskgroup_id did not suceed')
         end if
       end function get_taskgroup_id


       ! Get the ID of task iiproc on taskgroup iitg
       integer function get_rank_on_taskgroup(iiproc,iitg)
         implicit none
         ! Calling arguments
         integer,intent(in) :: iiproc, iitg
         ! Local variables
         integer :: jproc
         logical :: found

         found = .false.
         do jproc=0,smat%nranks(iitg)-1
             if (smat%tgranks(jproc,iitg) == iiproc) then
                 get_rank_on_taskgroup = jproc
                 found = .true.
                 exit
             end if
         end do
         if (.not.found) stop 'get_rank_on_taskgroup did not suceed'
       end function get_rank_on_taskgroup

  end subroutine compress_matrix_distributed_core



  subroutine uncompress_matrix_distributed2(iproc, smat, layout, matrix_compr, matrixp)
    use dynamic_memory
    implicit none

    ! Calling arguments
    integer,intent(in) :: iproc, layout
    type(sparse_matrix),intent(in) :: smat
    real(kind=mp),dimension(smat%nvctrp_tg),intent(in) :: matrix_compr
    real(kind=mp),dimension(:,:),intent(out) :: matrixp

    ! Local variables
    integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb, nfvctrp, isfvctr

     !call timing(iproc,'compressd_mcpy','ON')
     call f_timing(TCAT_SMAT_COMPRESSION,'ON')

     ! Check the dimensions of the output array and assign some values
     if (size(matrixp,1)/=smat%nfvctr) stop 'size(matrixp,1)/=smat%nfvctr'
     if (layout==DENSE_PARALLEL) then
         if (size(matrixp,2)/=smat%nfvctrp) stop '(ubound(matrixp,2)/=smat%nfvctrp'
         nfvctrp=smat%nfvctrp
         isfvctr=smat%isfvctr
     else if (layout==DENSE_MATMUL) then
         if (.not.smat%smatmul_initialized) then
             call f_err_throw('sparse matrix multiplication not initialized', &
                  err_name='SPARSEMATRIX_RUNTIME_ERROR')
         end if
         if (size(matrixp,2)/=smat%smmm%nfvctrp) stop '(ubound(matrixp,2)/=smat%smmm%nfvctrp'
         nfvctrp=smat%smmm%nfvctrp
         isfvctr=smat%smmm%isfvctr
     end if

       if (nfvctrp>0) then

          !call to_zero(smat%nfvctr*nfvctrp,matrixp(1,1))
          call f_zero(matrixp)

           isegstart=smat%istsegline(isfvctr+1)
           isegend=smat%istsegline(isfvctr+nfvctrp)+smat%nsegline(isfvctr+nfvctrp)-1
           !!isegstart=smat%istsegline(smat%isfvctr_par(iproc)+1)
           !!if (smat%isfvctr_par(iproc)+smat%nfvctrp<smat%nfvctr) then
           !!    isegend=smat%istsegline(smat%isfvctr_par(iproc+1)+1)-1
           !!else
           !!    isegend=smat%nseg
           !!end if
           !$omp parallel do default(private) &
           !$omp shared(isegstart, isegend, smat, matrixp, matrix_compr, isfvctr)
           do iseg=isegstart,isegend
               ii=smat%keyv(iseg)-1
               ! A segment is always on one line, therefore no double loop
               do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                   ii=ii+1
                   iiorb = smat%keyg(1,2,iseg)
                   jjorb = jorb
                   matrixp(jjorb,iiorb-isfvctr) = matrix_compr(ii-smat%isvctrp_tg)
               end do
           end do
           !$omp end parallel do
       end if

     !call timing(iproc,'compressd_mcpy','OF')
     call f_timing(TCAT_SMAT_COMPRESSION,'OF')

   end subroutine uncompress_matrix_distributed2



   subroutine sequential_acces_matrix_fast(smat, a, a_seq)
     use dynamic_memory
     implicit none
   
     ! Calling arguments
     type(sparse_matrix),intent(in) :: smat
     real(kind=mp),dimension(smat%nvctr),intent(in) :: a
     real(kind=mp),dimension(smat%smmm%nseq),intent(out) :: a_seq
   
     ! Local variables
     integer :: ii
     integer(kind=8) :: iseq

     call f_routine(id='sequential_acces_matrix_fast')

     if (.not.smat%smatmul_initialized) then
         call f_err_throw('sparse matrix multiplication not initialized', &
              err_name='SPARSEMATRIX_RUNTIME_ERROR')
     end if
   
     !$omp parallel do default(none) private(iseq, ii) &
     !$omp shared(smat, a_seq, a)
     do iseq=1,smat%smmm%nseq
         ii=smat%smmm%indices_extract_sequential(iseq)
         a_seq(iseq)=a(ii)
     end do
     !$omp end parallel do

     call f_release_routine()
   
   end subroutine sequential_acces_matrix_fast

   subroutine sequential_acces_matrix_fast2(smat, a, a_seq)
     use dynamic_memory
     implicit none
   
     ! Calling arguments
     type(sparse_matrix),intent(in) :: smat
     real(kind=mp),dimension(smat%nvctrp_tg),intent(in) :: a
     real(kind=mp),dimension(smat%smmm%nseq),intent(out) :: a_seq
   
     ! Local variables
     integer :: ii
     integer(kind=8) :: iseq

     call f_routine(id='sequential_acces_matrix_fast2')

     if (.not.smat%smatmul_initialized) then
         call f_err_throw('sparse matrix multiplication not initialized', &
              err_name='SPARSEMATRIX_RUNTIME_ERROR')
     end if
   
     !$omp parallel do schedule(guided) &
     !$omp default(none) private(iseq, ii) &
     !$omp shared(smat, a_seq, a)
     do iseq=1,smat%smmm%nseq
         ii=smat%smmm%indices_extract_sequential(iseq)
         a_seq(iseq)=a(ii-smat%isvctrp_tg)
     end do
     !$omp end parallel do

     call f_release_routine()
   
   end subroutine sequential_acces_matrix_fast2




   subroutine sparsemm_new(iproc, smat, a_seq, b, c)
     use yaml_output
     use dynamic_memory
     implicit none
   
     !Calling Arguments
     integer,intent(in) :: iproc
     type(sparse_matrix),intent(in) :: smat
     real(kind=mp), dimension(smat%smmm%nvctrp),intent(in) :: b
     real(kind=mp), dimension(smat%smmm%nseq),intent(in) :: a_seq
     real(kind=mp), dimension(smat%smmm%nvctrp), intent(out) :: c
   
     !Local variables
     !character(len=*), parameter :: subname='sparsemm'
     integer :: i,jorb,jjorb,iend,nblock, iblock, ncount
     integer :: ii, ilen, iout, iiblock, isblock, is,ie
     real(kind=mp) :: tt0
     integer :: n_dense
     real(kind=mp),dimension(:,:),allocatable :: a_dense, b_dense, c_dense
     !real(kind=mp),dimension(:),allocatable :: b_dense, c_dense
     !!integer,parameter :: MATMUL_NEW = 101
     !!integer,parameter :: MATMUL_OLD = 102
     !!integer,parameter :: matmul_version = MATMUL_NEW 
     logical,parameter :: count_flops = .false.
     real(kind=mp) :: ts, te, op, gflops
     real(kind=mp),parameter :: flop_per_op = 2.d0 !<number of FLOPS per operations
   
     call f_routine(id='sparsemm')

     if (.not.smat%smatmul_initialized) then
         call f_err_throw('sparse matrix multiplication not initialized', &
              err_name='SPARSEMATRIX_RUNTIME_ERROR')
     end if

     if (count_flops) then
         n_dense = nint(sqrt(real(smat%smmm%nseq,kind=mp)))
         !n_dense = smat%nfvctr
         a_dense = f_malloc((/n_dense,n_dense/),id='a_dense')
         b_dense = f_malloc((/n_dense,1/),id='b_dense')
         c_dense = f_malloc((/n_dense,1/),id='c_dense')
     end if
     !call timing(iproc, 'sparse_matmul ', 'IR')
     call f_timing(TCAT_SMAT_MULTIPLICATION,'IR')


     ! The choice for matmul_version can be made in sparsematrix_base
     if (matmul_version==MATMUL_NEW) then

         if (count_flops) then
             ! Start time
             ts = mpi_wtime()
         end if
         !$omp parallel default(private) shared(smat, a_seq, b, c)
         !$omp do schedule(guided)
         do iout=1,smat%smmm%nout
             i=smat%smmm%onedimindices_new(1,iout)
             nblock=smat%smmm%onedimindices_new(4,iout)
             isblock=smat%smmm%onedimindices_new(5,iout)
             tt0=0.d0

             is = isblock + 1
             ie = isblock + nblock
             !do iblock=1,nblock
             do iblock=is,ie
                 !iiblock = isblock + iblock
                 jorb = smat%smmm%consecutive_lookup(1,iblock)
                 jjorb = smat%smmm%consecutive_lookup(2,iblock)
                 ncount = smat%smmm%consecutive_lookup(3,iblock)
                 !tt0 = tt0 + ddot(ncount, b(jjorb), 1, a_seq(jorb), 1)
                 !avoid calling ddot from OpenMP region on BG/Q as too expensive
                 !tt0=tt0+my_dot(ncount,b(jjorb:jjorb+ncount-1),a_seq(jorb:jorb+ncount-1))
                 tt0=tt0+my_dot(ncount,b(jjorb),a_seq(jorb))
             end do

             c(i) = tt0
         end do 
         !$omp end do
         !$omp end parallel

         if (count_flops) then
             ! End time
             te = mpi_wtime()
             ! Count the operations
             op = 0.d0
             do iout=1,smat%smmm%nout
                 nblock=smat%smmm%onedimindices_new(4,iout)
                 isblock=smat%smmm%onedimindices_new(5,iout)
                 do iblock=1,nblock
                     iiblock = isblock + iblock
                     ncount = smat%smmm%consecutive_lookup(3,iiblock)
                     op = op + real(ncount,kind=mp)
                 end do
             end do
             gflops = 1.d-9*op/(te-ts)*flop_per_op
             call yaml_map('SPARSE: operations',op,fmt='(es9.3)')
             call yaml_map('SPARSE: time',te-ts,fmt='(es9.3)')
             call yaml_map('SPARSE: GFLOPS',gflops)

             ! Compare with dgemm of comparable size
             ts = mpi_wtime()
             call dgemv('n', n_dense, n_dense, 1.d0, a_dense, n_dense, b_dense, 1, 0.d0, c_dense, 1)
             !call dgemm('n', 'n', n_dense, n_dense, n_dense, 1.d0, a_dense, n_dense, &
             !     b_dense, n_dense, 0.d0, c_dense, n_dense)
             te = mpi_wtime()
             op = real(n_dense,kind=mp)*real(n_dense,kind=mp)
             gflops = 1.d-9*op/(te-ts)*flop_per_op
             call yaml_map('DGEMV: operations',op,fmt='(es9.3)')
             call yaml_map('DGEMV: time',te-ts,fmt='(es9.3)')
             call yaml_map('DGEMV: GFLOPS',gflops)
         end if


     else if (matmul_version==MATMUL_OLD) then

         !$omp parallel default(private) shared(smat, a_seq, b, c)
         !$omp do
         do iout=1,smat%smmm%nout
             i=smat%smmm%onedimindices_new(1,iout)
             ilen=smat%smmm%onedimindices_new(2,iout)
             ii=smat%smmm%onedimindices_new(3,iout)
             tt0=0.d0

             iend=ii+ilen-1

             do jorb=ii,iend
                jjorb=smat%smmm%ivectorindex_new(jorb)
                tt0 = tt0 + b(jjorb)*a_seq(jorb)
             end do
             c(i) = tt0
         end do 
         !$omp end do
         !$omp end parallel


     else

         stop 'wrong value of matmul_version'

     end if

   
     !call timing(iproc, 'sparse_matmul ', 'RS')
     call f_timing(TCAT_SMAT_MULTIPLICATION,'RS')
     call f_release_routine()

        contains

     pure function my_dot(n,x,y) result(tt)
       implicit none
       integer , intent(in) :: n
       double precision :: tt
       double precision, dimension(n), intent(in) :: x,y
       !local variables
       integer :: i

       tt=0.d0
       do i=1,n
          tt=tt+x(i)*y(i)
       end do
     end function
       
   end subroutine sparsemm_new


   !!function orb_from_index(smat, ival)
   !!  use sparsematrix_base, only: sparse_matrix
   !!  implicit none
   !!  ! Calling arguments
   !!  type(sparse_matrix),intent(in) :: smat
   !!  integer,intent(in) :: ival
   !!  integer,dimension(2) :: orb_from_index

   !!  orb_from_index(2) = (ival-1)/smat%nfvctr + 1
   !!  !orb_from_index(1) = ival - (orb_from_index_fn(2)-1)*smat%nfvctr
   !!  orb_from_index(1) = mod(ival-1,smat%nfvctr) + 1

   !!end function orb_from_index


   subroutine gather_matrix_spin(nproc,mat_tg,count,mat_global,ld_global,recvcounts,displs,comm)
     implicit none
     integer, intent(in) :: count,ld_global,comm,nproc
     real(mp), dimension(count), intent(in) :: mat_tg
     real(mp), dimension(ld_global), intent(out) :: mat_global
     integer, dimension(nproc), intent(in) :: recvcounts,displs
     
      call fmpi_allgather(sendbuf=mat_tg,sendcount=count,&
           recvbuf=mat_global,recvcounts=recvcounts, &
           displs=displs,comm=comm,algorithm=ONESIDED_ENUM)

    end subroutine gather_matrix_spin

    subroutine gather_matrix_from_taskgroups(iproc, nproc, comm, smat, mat_tg, mat_global)
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      type(sparse_matrix),intent(in) :: smat
      real(kind=mp),dimension(smat%nvctrp_tg*smat%nspin),intent(in) :: mat_tg !< matrix distributed over the taskgroups
      real(kind=mp),dimension(smat%nvctr*smat%nspin),intent(out) :: mat_global !< global matrix gathered together

      ! Local variables
      integer,dimension(:),allocatable :: recvcounts, recvdspls
      integer :: ncount, ist_send, jproc, ispin, ishift

      call f_routine(id='gather_matrix_from_taskgroups')

      if (nproc>1) then
         recvcounts = f_malloc0(0.to.nproc-1,id='recvcounts')
         recvdspls = f_malloc0(0.to.nproc-1,id='recvdspls')
         !call to_zero(nproc, recvcounts(0))
         !call to_zero(nproc, recvdspls(0))
         !ncount = smat%smmm%istartend_mm_dj(2) - smat%smmm%istartend_mm_dj(1) + 1
         ncount = smat%nvctrp
         recvcounts(iproc) = ncount
         call fmpi_allreduce(recvcounts(0), nproc, FMPI_SUM, comm=comm)
         recvdspls(0) = 0
         do jproc=1,nproc-1
            recvdspls(jproc) = recvdspls(jproc-1) + recvcounts(jproc-1)
         end do
         do ispin=1,smat%nspin
            ishift = (ispin-1)*smat%nvctr
            !ist_send = smat%smmm%istartend_mm_dj(1) - smat%isvctrp_tg + ishift
            ist_send = smat%isvctr + 1 - smat%isvctrp_tg + ishift
            ! The following condition is necessary for ncount=0, in order to avoid out of bound problems
            ist_send = min(ist_send,ispin*smat%nvctrp_tg)
            call gather_matrix_spin(nproc,mat_tg(ist_send),ncount,mat_global(ishift+1),&
                 sum(recvcounts),recvcounts,recvdspls,comm)
!!$            call fmpi_allgather(mat_tg(ist_send), ncount, &
!!$                 mat_global(ishift+1), recvcounts=recvcounts, &
!!$                 displs=recvdspls, comm=comm,&
!!$                 algorithm=ONESIDED_ENUM)
            !!call mpi_allgatherv(mat_tg(ist_send), ncount, mpi_double_precision, &
            !!                    mat_global(1), recvcounts, recvdspls, mpi_double_precision, &
            !!                    bigdft_mpi%mpi_comm, ierr)
         end do
         call f_free(recvcounts)
         call f_free(recvdspls)
      else
         call vcopy(smat%nvctrp*smat%nspin, mat_tg(1), 1, mat_global(1), 1)
      end if

      call f_release_routine()

    end subroutine gather_matrix_from_taskgroups


   subroutine gather_matrix_from_taskgroups_deprecated(iproc, nproc, comm, smat, mat_tg, mat_global)
     use dynamic_memory
     implicit none
   
     ! Calling arguments
     integer,intent(in) :: iproc, nproc, comm
     type(sparse_matrix),intent(in) :: smat
     real(kind=mp),dimension(smat%nvctrp_tg*smat%nspin),intent(in) :: mat_tg !< matrix distributed over the taskgroups
     real(kind=mp),dimension(smat%nvctr*smat%nspin),intent(out) :: mat_global !< global matrix gathered together
   
     ! Local variables
     integer,dimension(:),allocatable :: recvcounts, recvdspls
     integer :: ncount, ist_send, jproc, ispin, ishift

     call f_routine(id='gather_matrix_from_taskgroups')
   
     if (nproc>1) then
         recvcounts = f_malloc0(0.to.nproc-1,id='recvcounts')
         recvdspls = f_malloc0(0.to.nproc-1,id='recvdspls')
         !call to_zero(nproc, recvcounts(0))
         !call to_zero(nproc, recvdspls(0))
         !ncount = smat%smmm%istartend_mm_dj(2) - smat%smmm%istartend_mm_dj(1) + 1
         ncount = smat%nvctrp
         recvcounts(iproc) = ncount
         call fmpi_allreduce(recvcounts(0), nproc, FMPI_SUM, comm=comm)
         recvdspls(0) = 0
         do jproc=1,nproc-1
             recvdspls(jproc) = recvdspls(jproc-1) + recvcounts(jproc-1)
         end do
         do ispin=1,smat%nspin
             ishift = (ispin-1)*smat%nvctr
             !ist_send = smat%smmm%istartend_mm_dj(1) - smat%isvctrp_tg + ishift
             ist_send = smat%isvctr + 1 - smat%isvctrp_tg + ishift
             ! The following condition is necessary for ncount=0, in order to avoid out of bound problems
             ist_send = min(ist_send,ispin*smat%nvctrp_tg)
             call fmpi_allgather(mat_tg(ist_send), ncount, &
                  mat_global(ishift+1), recvcounts=recvcounts, &
                  displs=recvdspls, comm=comm,&
                  algorithm=ONESIDED_ENUM)
             !!call mpi_allgatherv(mat_tg(ist_send), ncount, mpi_double_precision, &
             !!                    mat_global(1), recvcounts, recvdspls, mpi_double_precision, &
             !!                    bigdft_mpi%mpi_comm, ierr)
         end do
         call f_free(recvcounts)
         call f_free(recvdspls)
     else
         call vcopy(smat%nvctrp*smat%nspin, mat_tg(1), 1, mat_global(1), 1)
     end if

     call f_release_routine()

   end subroutine gather_matrix_from_taskgroups_deprecated


   subroutine gather_matrix_from_taskgroups_inplace(iproc, nproc, comm, smat, mat)
     use dynamic_memory
     implicit none

   
     ! Calling arguments
     integer,intent(in) :: iproc, nproc, comm
     type(sparse_matrix),intent(in) :: smat
     type(matrices),intent(inout) :: mat
   
     ! Local variables
     integer,dimension(:),allocatable :: recvcounts, recvdspls
     integer :: ncount, ist_send, jproc, ispin, ishift
     real(kind=mp),dimension(:),allocatable :: mat_global

     call f_routine(id='gather_matrix_from_taskgroups_inplace')

     if (.not.smat%smatmul_initialized) then
         call f_err_throw('sparse matrix multiplication not initialized', &
              err_name='SPARSEMATRIX_RUNTIME_ERROR')
     end if
   
     mat_global = sparsematrix_malloc(smat,iaction=SPARSE_FULL,id='mat_global')
     if (nproc>1) then
         recvcounts = f_malloc0(0.to.nproc-1,id='recvcounts')
         recvdspls = f_malloc0(0.to.nproc-1,id='recvdspls')
         !call to_zero(nproc, recvcounts(0))
         !call to_zero(nproc, recvdspls(0))
         ncount = smat%smmm%istartend_mm_dj(2) - smat%smmm%istartend_mm_dj(1) + 1
         recvcounts(iproc) = ncount
         call fmpi_allreduce(recvcounts(0), nproc, FMPI_SUM, comm=comm)
         recvdspls(0) = 0
         do jproc=1,nproc-1
             recvdspls(jproc) = recvdspls(jproc-1) + recvcounts(jproc-1)
         end do
         do ispin=1,smat%nspin
             ishift = (ispin-1)*smat%nvctr
             ist_send = smat%smmm%istartend_mm_dj(1) - smat%isvctrp_tg + ishift
             call fmpi_allgather(mat%matrix_compr(ist_send), ncount, &
                  mat_global(ishift+1), recvcounts=recvcounts, &
                  displs=recvdspls, comm=comm,&
                  algorithm=ONESIDED_ENUM)
             !!call mpi_allgatherv(mat%matrix_compr(ist_send), ncount, mpi_double_precision, &
             !!                    mat_global(1), recvcounts, recvdspls, mpi_double_precision, &
             !!                    bigdft_mpi%mpi_comm, ierr)
         end do
         call f_free(recvcounts)
         call f_free(recvdspls)
     else
         call vcopy(smat%nvctrp_tg*smat%nspin, mat%matrix_compr(1), 1, mat_global(1), 1)
     end if

     call f_free_ptr(mat%matrix_compr)
     mat%matrix_compr = sparsematrix_malloc_ptr(smat,iaction=SPARSE_FULL,id='mat_global')
     call f_memcpy(src=mat_global, dest=mat%matrix_compr)
     call f_free(mat_global)

     call f_release_routine()

   end subroutine gather_matrix_from_taskgroups_inplace


   subroutine extract_taskgroup_inplace(smat, mat)
     implicit none
   
     ! Calling arguments
     type(sparse_matrix),intent(in) :: smat
     type(matrices),intent(inout) :: mat

     ! Local variables
     integer :: i, ispin, ishift_tg, ishift_glob

     do ispin=1,smat%nspin
         ishift_tg = (ispin-1)*smat%nvctrp_tg
         ishift_glob = (ispin-1)*smat%nvctr
         do i=1,smat%nvctrp_tg
             mat%matrix_compr(i+ishift_tg) = mat%matrix_compr(i+smat%isvctrp_tg+ishift_glob)
         end do
     end do

   end subroutine extract_taskgroup_inplace


   subroutine extract_taskgroup(smat, mat_glob, mat_tg)
     implicit none
   
     ! Calling arguments
     type(sparse_matrix),intent(in) :: smat
     real(kind=mp),dimension(smat%nvctr*smat%nspin),intent(in) :: mat_glob
     real(kind=mp),dimension(smat%nvctrp_tg*smat%nspin),intent(out) :: mat_tg

     ! Local variables
     integer :: i, ispin, ishift_tg, ishift_glob

     do ispin=1,smat%nspin
         ishift_tg = (ispin-1)*smat%nvctrp_tg
         ishift_glob = (ispin-1)*smat%nvctr
         do i=1,smat%nvctrp_tg
             mat_tg(i+ishift_tg) = mat_glob(i+smat%isvctrp_tg+ishift_glob)
         end do
     end do

   end subroutine extract_taskgroup


    subroutine write_matrix_compressed(iproc, nproc, comm, message, smat, mat)
      use yaml_output
      use dynamic_memory
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      character(len=*),intent(in) :: message
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
    
      ! Local variables
      !integer, dimension(2) :: irowcol
      integer :: iseg, i, ii
      real(kind=mp),dimension(:),allocatable :: matrix_compr
      !integer :: iorb, jorb
    
      !!call yaml_sequence_open(trim(message))
      !!do iseg=1,smat%nseg
      !!    call yaml_sequence(advance='no')
      !!    ilen=smat%keyg(2,iseg)-smat%keyg(1,iseg)+1
      !!    call yaml_mapping_open(flow=.true.)
      !!    call yaml_map('segment',iseg)
      !!    istart=smat%keyv(iseg)
      !!    iend=smat%keyv(iseg)+ilen
      !!    call yaml_map('values',smat%matrix_compr(istart:iend))
      !!    call yaml_mapping_close()
      !!    call yaml_newline()
      !!end do
      !!call yaml_sequence_close()

      matrix_compr = sparsematrix_malloc(smat,iaction=SPARSE_FULL,id='matrix_compr')
      call gather_matrix_from_taskgroups(iproc, nproc, comm, &
           smat, mat%matrix_compr, matrix_compr)
    
      if (iproc==0) then
          call yaml_sequence_open(trim(message))
          do iseg=1,smat%nseg
              ! A segment is always on one line, therefore no double loop
              call yaml_sequence(advance='no')
              !ilen=smat%keyg(2,iseg)-smat%keyg(1,iseg)+1
              call yaml_mapping_open(flow=.true.)
              call yaml_map('segment',iseg)
              call yaml_sequence_open('elements')
              !istart=smat%keyv(iseg)
              !iend=smat%keyv(iseg)+ilen-1
              !do i=istart,iend
              ii=smat%keyv(iseg)
              do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                  call yaml_newline()
                  call yaml_sequence(advance='no')
                  call yaml_mapping_open(flow=.true.)
                  !irowcol=orb_from_index(smat,i)
                  !iorb=orb_from_index(1,i)
                  !jorb=orb_from_index(2,i)
                  call yaml_map('coordinates',(/smat%keyg(1,2,iseg),i/))
                  call yaml_map('value',matrix_compr(ii))
                  call yaml_mapping_close()
                  ii=ii+1
              end do
              call yaml_sequence_close()
              !call yaml_map('values',smat%matrix_compr(istart:iend))
              call yaml_mapping_close()
              call yaml_newline()
          end do
          call yaml_sequence_close()
      end if

      call f_free(matrix_compr)
    
    end subroutine write_matrix_compressed


     !integer :: mp1, jjorb0, jjorb1, jjorb2, jjorb3, jjorb4, jjorb5, jjorb6


    !> Write a sparse matrix to a file
    subroutine write_sparsematrix(filename, smat, mat)
      use yaml_output
      use dynamic_memory
      implicit none
    
      ! Calling arguments
      character(len=*),intent(in) :: filename
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
    
      ! Local variables
      integer :: iseg, i, ii
      integer,parameter :: iunit=234

      call f_routine(id='write_sparsematrix')

      ! First check that no taskgroups are used. Otherwise this routine does not work
      if (smat%ntaskgroup>1) then
          call f_err_throw('write_sparsematrix has not yet been implememted for matrix taskgroups', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if
    
      open(unit=iunit,file=filename)

      write(iunit,*) smat%nfvctr, '# number of columns'
      write(iunit,*) smat%nseg, '# number of segments'
      write(iunit,*) smat%nvctr, '# number of non-zero elements'
      do iseg=1,smat%nseg
          if(iseg==1) then
              write(iunit,*) smat%keyv(iseg), '# values of keyv'
          else
              write(iunit,*) smat%keyv(iseg)
          end if
      end do
      do iseg=1,smat%nseg
          if(iseg==1) then
              write(iunit,*) smat%keyg(1:2,1:2,iseg), '# values of keyg'
          else
              write(iunit,*) smat%keyg(1:2,1:2,iseg)
          end if
      end do
    
      do iseg=1,smat%nseg
          ! A segment is always on one line, therefore no double loop
          ii=smat%keyv(iseg)
          do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
              if (i==1 .and. iseg==1) then
                  write(iunit,*) mat%matrix_compr(ii), '# values of matrix_compr'
              else
                  write(iunit,*) mat%matrix_compr(ii)
              end if
              ii=ii+1
          end do
      end do

      close(unit=iunit)

      call f_release_routine()
    
    end subroutine write_sparsematrix



    !> Write a sparse matrix to a file, using the CCS format
    subroutine write_sparsematrix_CCS(filename, smat, mat)
      use yaml_output
      use dynamic_memory
      implicit none
    
      ! Calling arguments
      character(len=*),intent(in) :: filename
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
    
      ! Local variables
      integer :: iseg, i, ii, icol, imat
      integer, dimension(:), allocatable :: col_ptr, row_ind
     ! logical, dimension(:,:), allocatable :: matg
      logical :: first_in_column_set
      !integer :: j
      !logical :: column_started
      real(kind=mp),dimension(:),allocatable :: val
      integer,parameter :: iunit=234, iunit2=235
      character(len=10) :: num
      character(len=100) :: frmt

      call f_routine(id='write_sparsematrix_CCS')

      col_ptr = f_malloc(smat%nfvctr,id='col_ptr')
      row_ind = f_malloc(smat%nvctr,id='row_ind')
      val = f_malloc(smat%nvctr,id='val')

      ii = 0
      do icol=1,smat%nfvctr
          imat = 0
          first_in_column_set = .false.
          do iseg=1,smat%nseg
              do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                  imat = imat + 1
                  if (i==icol) then
                      ! We are in column icol
                      ii = ii + 1
                      row_ind(ii) = smat%keyg(1,2,iseg) !row index
                      val(ii) = mat%matrix_compr(imat)
                      if (.not.first_in_column_set) then
                          col_ptr(icol) = ii
                          first_in_column_set = .true.
                      end if
                  end if
              end do
          end do
      end do
      if (ii/=smat%nvctr) stop 'ERROR in write_sparsematrix_CCS: ii/=smat%nvctr'

      !!matg = f_malloc((/smat%nfvctr,smat%nfvctr/),id='matg')
      !!matg = .false.
      !!do iseg=1,smat%nseg
      !!    ! A segment is always on one line, therefore no double loop
      !!    ii=smat%keyv(iseg)
      !!    do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
      !!        matg(smat%keyg(1,2,iseg),i) = .true.
      !!    end do
      !!end do

      !!ii = 0
      !!do i=1,smat%nfvctr
      !!    column_started = .false.
      !!    do j=1,smat%nfvctr
      !!       if(matg(j,i)) then
      !!           ii = ii + 1
      !!           row_ind(ii) = j
      !!           if (.not.column_started) then
      !!               col_ptr(i) = ii
      !!               column_started = .true.
      !!           end if
      !!       end if
      !!    end do
      !!end do
    
      open(unit=iunit,file=trim(filename))
      open(unit=iunit2,file=trim(filename)//'_2')

      write(iunit,*) smat%nfvctr, smat%nvctr, '# number of rows/columns, number of non-zero entries'
      write(iunit2,*) smat%nfvctr, smat%nfvctr, smat%nvctr
      do i=1,smat%nfvctr
          if (i==1) then
              write(iunit,*) col_ptr(i), '# col_ptr'
          else
              write(iunit,*) col_ptr(i)
          end if
      end do
      write(num,'(i0)') smat%nfvctr
      frmt='('//num//'(i0,1x))'
      write(iunit2,trim(frmt)) (col_ptr(i),i=1,smat%nfvctr)

      do i=1,smat%nvctr
          if (i==1) then
              write(iunit,*) row_ind(i), '# row_ind'
          else
              write(iunit,*) row_ind(i)
          end if
      end do
      write(num,'(i0)') smat%nvctr
      frmt='('//num//'(i0,1x))'
      write(iunit2,trim(frmt)) (row_ind(i),i=1,smat%nvctr)
      
      do i=1,smat%nvctr
          if(i==1) then
              !!write(iunit,*) mat%matrix_compr(i), '# values of matrix_compr' 
              !!write(iunit2,*) mat%matrix_compr(i)
              write(iunit,*) val(i), '# values of matrix_compr' 
              write(iunit2,*) val(i)
          else
              !write(iunit,*) mat%matrix_compr(i) 
              !write(iunit2,*) mat%matrix_compr(i) 
              write(iunit,*) val(i) 
              write(iunit2,*) val(i) 
          end if
      end do
      !write(num,'(i0)') smat%nvctr
      !frmt='('//num//'i9)'
      !write(iunit2,trim(frmt)) (mat%matrix_compr(i),i=1,smat%nvctr)

      close(unit=iunit)
      close(unit=iunit2)

      call f_free(col_ptr)
      call f_free(row_ind)
      call f_free(val)
      !!call f_free(matg)

      call f_release_routine()
    
    end subroutine write_sparsematrix_CCS


    !> Transform a matrix from a large parsity pattern *_l to a small sparsity pattern *_s or vice versa.
    !! The small pattern must be contained within the large one.
    subroutine transform_sparsity_pattern(iproc, nfvctr, nvctrp_s, isvctr_s, nseg_s, keyv_s, keyg_s, line_and_column_s, &
               nvctrp_l, isvctr_l, nseg_l, keyv_l, keyg_l, istsegline_l, direction, &
               matrix_s_in, matrix_l_in, matrix_s_out, matrix_l_out)
      use sparsematrix_init, only: matrixindex_in_compressed_lowlevel
      use dynamic_memory
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nfvctr, nvctrp_s, isvctr_s, nseg_s, nvctrp_l, isvctr_l, nseg_l
      integer,dimension(2,nvctrp_s),intent(in) :: line_and_column_s
      integer,dimension(nseg_s),intent(in) :: keyv_s
      integer,dimension(2,2,nseg_s),intent(in) :: keyg_s
      integer,dimension(nseg_l),intent(in) :: keyv_l
      integer,dimension(2,2,nseg_l),intent(in) :: keyg_l
      integer,dimension(nfvctr),intent(in) :: istsegline_l
      character(len=*),intent(in) :: direction
      real(kind=mp),dimension(nvctrp_l),intent(in),optional :: matrix_l_in
      real(kind=mp),dimension(nvctrp_s),intent(in),optional :: matrix_s_in
      real(kind=mp),dimension(nvctrp_l),intent(out),optional :: matrix_l_out
      real(kind=mp),dimension(nvctrp_s),intent(out),optional :: matrix_s_out
      ! Local variables
      integer :: i, ii, ind, iline, icolumn

      call f_routine(id='transform_sparsity_pattern')
      !call timing(iproc, 'transformspars', 'ON')
      call f_timing(TCAT_SMAT_TRANSFORMATION,'ON')

        if (direction=='large_to_small') then

            if (.not.present(matrix_l_in)) then
                call f_err_throw("'matrix_l_in' not present",err_name='BIGDFT_RUNTIME_ERROR')
            end if
            if (.not.present(matrix_s_out)) then
                call f_err_throw("'matrix_s_out' not present",err_name='BIGDFT_RUNTIME_ERROR')
            end if

            ! No need for f_zero since every value will be overwritten.
            !$omp parallel default(none) &
            !$omp shared(nvctrp_s, isvctr_s, isvctr_l, line_and_column_s) &
            !$omp shared(nfvctr, nseg_l, keyv_l, keyg_l, istsegline_l, matrix_s_out, matrix_l_in) &
            !$omp private(i, ii, iline, icolumn, ind)
            !$omp do schedule(guided)
            do i=1,nvctrp_s
                ii = isvctr_s + i
                !!call get_line_and_column(ii, nseg_s, keyv_s, keyg_s, iline, icolumn)
                iline = line_and_column_s(1,i)
                icolumn = line_and_column_s(2,i)
                ind = matrixindex_in_compressed_lowlevel(icolumn, iline, nfvctr, &
                      nseg_l, keyv_l, keyg_l, istsegline_l)
                ind = ind - isvctr_l
                matrix_s_out(i) = matrix_l_in(ind)
            end do
            !$omp end do
            !$omp end parallel

        else if (direction=='small_to_large') then
            if (.not.present(matrix_s_in)) then
                call f_err_throw("'matrix_s_in' not present",err_name='BIGDFT_RUNTIME_ERROR')
            end if
            if (.not.present(matrix_l_out)) then
                call f_err_throw("'matrix_l_out' not present",err_name='BIGDFT_RUNTIME_ERROR')
            end if
            call f_zero(matrix_l_out)
            !$omp parallel default(none) &
            !$omp shared(nvctrp_s, isvctr_s, isvctr_l, line_and_column_s) &
            !$omp shared(nfvctr, nseg_l, keyv_l, keyg_l, istsegline_l, matrix_s_in, matrix_l_out) &
            !$omp private(i, ii, iline, icolumn, ind)
            !$omp do schedule(guided)
            do i=1,nvctrp_s
                ii = isvctr_s + i
                !call get_line_and_column(ii, nseg_s, keyv_s, keyg_s, iline, icolumn)
                iline = line_and_column_s(1,i)
                icolumn = line_and_column_s(2,i)
                ind = matrixindex_in_compressed_lowlevel(icolumn, iline, nfvctr, &
                      nseg_l, keyv_l, keyg_l, istsegline_l)
                ind = ind - isvctr_l
                matrix_l_out(ind) = matrix_s_in(i)
            end do
            !$omp end do
            !$omp end parallel
        else
            stop 'wrong direction'
        end if

      !call timing(iproc, 'transformspars', 'OF')
      call f_timing(TCAT_SMAT_TRANSFORMATION,'OF')
      call f_release_routine()

    end subroutine transform_sparsity_pattern



    !> Calculates c = a*b for matrices a,b,c
    subroutine matrix_matrix_mult_wrapper(iproc, nproc, smat, a, b, c)
      use dynamic_memory
      implicit none



      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(sparse_matrix),intent(in) :: smat
      real(kind=mp),dimension(smat%nvctrp_tg),intent(in) :: a, b
      real(kind=mp),dimension(smat%nvctrp_tg),intent(inout) :: c

      ! Local variables
      real(kind=mp),dimension(:),allocatable :: b_exp, c_exp, a_seq

      call f_routine(id='matrix_matrix_mult_wrapper')

      if (.not.smat%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      b_exp = f_malloc(smat%smmm%nvctrp, id='b_exp')
      c_exp = f_malloc(smat%smmm%nvctrp, id='c_exp')
      a_seq = sparsematrix_malloc(smat, iaction=SPARSEMM_SEQ, id='a_seq')

      call sequential_acces_matrix_fast2(smat, a, a_seq)
      if (smat%smmm%nvctrp_mm>0) then !to avoid out of bounds error...
          call transform_sparsity_pattern(iproc, smat%nfvctr, smat%smmm%nvctrp_mm, smat%smmm%isvctr_mm, &
               smat%nseg, smat%keyv, smat%keyg, &
               smat%smmm%line_and_column_mm, &
               smat%smmm%nvctrp, smat%smmm%isvctr, &
               smat%smmm%nseg, smat%smmm%keyv, smat%smmm%keyg, smat%smmm%istsegline, &
               'small_to_large', matrix_s_in=b(smat%smmm%isvctr_mm-smat%isvctrp_tg+1), matrix_l_out=b_exp)
      end if
      call sparsemm_new(iproc, smat, a_seq, b_exp, c_exp)
      call compress_matrix_distributed_wrapper(iproc, nproc, smat, SPARSE_MATMUL_LARGE, &
           c_exp, ONESIDED_FULL, c)

      call f_free(b_exp)
      call f_free(c_exp)
      call f_free(a_seq)

      call f_release_routine()

    end subroutine matrix_matrix_mult_wrapper




    !< Calculates the trace of the sparse matrix mat
    function trace_sparse_matrix(iproc, nproc, comm, smat, mat) result(tr)
      use sparsematrix_init, only: matrixindex_in_compressed
      use dynamic_memory
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc,  nproc, comm
      type(sparse_matrix),intent(in) :: smat
      real(kind=mp),dimension(smat%nvctrp_tg),intent(in) :: mat
    
      ! Local variables
      integer :: irow, ind
      real(kind=mp) :: tr
    
      call f_routine(id='trace_sparse_matrix')

      tr = 0.0_mp 
      do irow=1,smat%nfvctr
          ind = matrixindex_in_compressed(smat, irow, irow)
          if (ind<smat%isvctr+1) cycle
          if (ind>smat%isvctr+smat%nvctrp) then
              exit
          end if
          tr = tr + mat(ind-smat%isvctrp_tg)
      end do
      if (nproc > 1) then
          call fmpi_allreduce(tr, 1, FMPI_SUM, comm=comm)
      end if
    
      call f_release_routine()
    
    end function trace_sparse_matrix





    !< Calculates the trace of the matrix product amat*bmat.
    !< WARNING: It is mandatory that the sparsity pattern of amat be contained
    !< within the sparsity pattern of bmat!
    function trace_sparse_matrix_product(iproc, nproc, comm, asmat, bsmat, amat, bmat) result(sumn)
      use sparsematrix_init, only: matrixindex_in_compressed
      use dynamic_memory
      implicit none

    
      ! Calling arguments
      integer,intent(in) :: iproc,  nproc, comm
      type(sparse_matrix),intent(in) :: asmat, bsmat
      real(kind=mp),dimension(asmat%nvctrp_tg),intent(in) :: amat
      real(kind=mp),dimension(bsmat%nvctrp_tg),intent(in) :: bmat
    
      ! Local variables
      integer :: iseg, ii, jorb, iiorb, jjorb, iilarge
      integer :: iashift, ibshift, iel
      real(kind=mp) :: sumn
    
    
      call f_routine(id='trace_sparse')
    
      iashift = 0!(ispin-1)*asmat%nvctr
      ibshift = 0!(ispin-1)*bsmat%nvctr
    
      !if (.not.asmat%smatmul_initialized) then
      !    call f_err_throw('The sparse matrix multiplications must &
      !         &be initialized to use the routine trace_sparse')
      !end if

      sumn=0.d0
      !$omp parallel default(none) &
      !$omp private(iseg, ii, jorb, iiorb, jjorb, iilarge, iel) &
      !$omp shared(bsmat, asmat, amat, bmat, iashift, ibshift, sumn)
      !$omp do reduction(+:sumn)
      do iseg=asmat%isseg,asmat%ieseg
          iel = asmat%keyv(iseg) - 1
          ii=iashift+asmat%keyv(iseg)-1
          ! A segment is always on one line, therefore no double loop
          do jorb=asmat%keyg(1,1,iseg),asmat%keyg(2,1,iseg)
              iel = iel + 1
              if (iel<asmat%isvctr+1) cycle
              if (iel>asmat%isvctr+asmat%nvctrp) then
                  exit
              end if
              ii=ii+1
              iiorb = asmat%keyg(1,2,iseg)
              jjorb = jorb
              iilarge = ibshift + matrixindex_in_compressed(bsmat, iiorb, jjorb)
              sumn = sumn + amat(ii-asmat%isvctrp_tg)*bmat(iilarge-bsmat%isvctrp_tg)
          end do  
      end do
      !$omp end do
      !$omp end parallel

      if (nproc > 1) then
          call fmpi_allreduce(sumn, 1, FMPI_SUM, comm=comm)
      end if
    
      call f_release_routine()
    
    end function trace_sparse_matrix_product


    !> Set to zero all term which couple different atoms
    subroutine delete_coupling_terms(iproc, nproc, comm, smmd, smat, mat_compr)
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      type(sparse_matrix_metadata),intent(in) :: smmd
      type(sparse_matrix),intent(in) :: smat
      real(kind=mp),dimension(smat%nvctrp_tg*smat%nspin),intent(inout) :: mat_compr

      ! Local variables
      integer :: ispin, ishift, iseg, ii, i, iiat, jjat
      real(kind=mp),dimension(:),allocatable :: fullmat_compr
      
      fullmat_compr = sparsematrix_malloc(smat,iaction=SPARSE_FULL,id='tmparr')
      call gather_matrix_from_taskgroups(iproc, nproc, comm, &
           smat, mat_compr, fullmat_compr)

      do ispin=1,smat%nspin
          ishift=(ispin-1)*smat%nvctr
          !!$omp parallel default(none) private(iseg,i,ii,irowcol) shared(sparsemat,inm,outm,ispin,ishift)
          !!$omp do
          do iseg=1,smat%nseg
              ii=smat%keyv(iseg)
              ! a segment is always on one line, therefore no double loop
              do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                 iiat = smmd%on_which_atom(i)
                 jjat = smmd%on_which_atom(smat%keyg(1,2,iseg))
                 if (iiat/=jjat) then
                     fullmat_compr(ii+ishift) = 0.d0
                 end if
                 ii=ii+1
             end do
          end do
          !!$omp end do
          !!$omp end parallel
      end do

      call extract_taskgroup(smat, fullmat_compr, mat_compr)

   end subroutine delete_coupling_terms


   subroutine synchronize_matrix_taskgroups(iproc, nproc, smat, mat)
     use dynamic_memory
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
    
      ! Local variables
      integer :: ncount, itg, iitg, ispin, ishift, ist_send, ist_recv
      integer,dimension(:),allocatable :: request
      real(kind=mp),dimension(:),allocatable :: recvbuf
    
      call f_routine(id='synchronize_matrix_taskgroups')
       
      if (nproc>1) then
          request = f_malloc(smat%ntaskgroupp,id='request')
          ncount = 0
          do itg=1,smat%ntaskgroupp
              iitg = smat%taskgroupid(itg)
              ncount = ncount + smat%taskgroup_startend(2,1,iitg)-smat%taskgroup_startend(1,1,iitg)+1
          end do
          recvbuf = f_malloc(ncount,id='recvbuf')
          do ispin=1,smat%nspin
              ishift = (ispin-1)*smat%nvctrp_tg
    
              ncount = 0
              do itg=1,smat%ntaskgroupp
                  iitg = smat%taskgroupid(itg)
                  ist_send = smat%taskgroup_startend(1,1,iitg) - smat%isvctrp_tg
                  ist_recv = ncount + 1
                  ncount = smat%taskgroup_startend(2,1,iitg)-smat%taskgroup_startend(1,1,iitg)+1
                  !!call mpi_iallreduce(mat%matrix_compr(ist_send), recvbuf(ist_recv), ncount, &
                  !!     mpi_double_precision, FMPI_SUM, smat%mpi_groups(iitg)%mpi_comm, request(itg), ierr)
                  if (smat%mpi_groups(iitg)%nproc>1) then
                     !call mpiiallred(mat%matrix_compr(ishift+ist_send), recvbuf(ist_recv), ncount, &
                     !      mpi_sum, smat%mpi_groups(iitg)%mpi_comm, request(itg))
                     call fmpi_allreduce(sendbuf=mat%matrix_compr(ishift+ist_send),&
                           recvbuf=recvbuf(ist_recv),count=ncount,op=FMPI_SUM,&
                           comm=smat%mpi_groups(iitg)%mpi_comm,request=request(itg))
                  else
                      call vcopy(ncount, mat%matrix_compr(ishift+ist_send), 1, recvbuf(ist_recv), 1)
                  end if
              end do
              if (smat%mpi_groups(iitg)%nproc > 1) then
                  call fmpi_waitall(smat%ntaskgroupp, request)
              end if
              ncount = 0
              do itg=1,smat%ntaskgroupp
                  iitg = smat%taskgroupid(itg)
                  ist_send = smat%taskgroup_startend(1,1,iitg) - smat%isvctrp_tg
                  ist_recv = ncount + 1
                  ncount = smat%taskgroup_startend(2,1,iitg)-smat%taskgroup_startend(1,1,iitg)+1
                  !call vcopy(ncount, recvbuf(ist_recv), 1, mat%matrix_compr(ishift+ist_send), 1)
                  call dcopy(ncount, recvbuf(ist_recv), 1, mat%matrix_compr(ishift+ist_send), 1)
              end do
          end do
          call f_free(request)
          call f_free(recvbuf)
      end if

      call f_release_routine()

    end subroutine synchronize_matrix_taskgroups



    subroutine max_asymmetry_of_matrix(iproc, nproc, comm, sparsemat, mat_tg, error_max)!, ispinx)
      use sparsematrix_init, only: matrixindex_in_compressed
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      type(sparse_matrix),intent(in) :: sparsemat
      real(kind=mp),dimension(sparsemat%nvctrp_tg),intent(in) :: mat_tg
      real(kind=mp),intent(out) :: error_max
      !integer,intent(in),optional :: ispinx

      ! Local variables
      !real(kind=mp),dimension(:),allocatable :: mat_full
      integer :: iseg, i, ind, ind_trans, iel
      real(kind=mp) :: val, val_trans, error

      call f_routine(id='max_asymmetry_of_matrix')

      !!! Gather together the matrix from the taskgroups
      !!mat_full = sparsematrix_malloc(sparsemat,iaction=SPARSE_FULL,id='mat_full')
      !!call gather_matrix_from_taskgroups(iproc, nproc, comm, sparsemat, mat_tg, mat_full)

      error_max = 0.0_mp
      !!do ispin=1,sparsemat%nspin
      !!    if (present(ispinx)) then
      !!        if (ispin/=ispinx) cycle
      !!    end if
      !!    ishift=(ispin-1)*sparsemat%nvctr
          ! SM: The function matrixindex_in_compressed is rather expensive, so I think OpenMP is always worth
          !$omp parallel default(none) &
          !$omp shared(sparsemat, mat_tg, error_max) &
          !$omp private(iseg, iel, i, ind, ind_trans, val, val_trans, error)
          !$omp do schedule(guided) reduction(max: error_max)
          do iseg=sparsemat%isseg,sparsemat%ieseg
              iel = sparsemat%keyv(iseg) - 1
              do i=sparsemat%keyg(1,1,iseg),sparsemat%keyg(2,1,iseg)
                  iel = iel + 1
                  if (iel<sparsemat%isvctr+1) cycle
                  if (iel>sparsemat%isvctr+sparsemat%nvctrp) then
                      exit
                  end if
                  ind = matrixindex_in_compressed(sparsemat, i, sparsemat%keyg(1,2,iseg)) - sparsemat%isvctrp_tg
                  ind_trans = matrixindex_in_compressed(sparsemat, sparsemat%keyg(1,2,iseg), i) - sparsemat%isvctrp_tg
                  !val = mat_full(ind)
                  !val_trans = mat_full(ind_trans)
                  val = mat_tg(ind)
                  val_trans = mat_tg(ind_trans)
                  error = abs(val-val_trans)
                  if (error>error_max) then
                      error_max = error
                  end if
              end do
          end do
          !$omp end do
          !$omp end parallel
      !!end do
      call fmpi_allreduce(error_max, 1, FMPI_MAX, comm=comm)
      !if (iproc==0) call yaml_map('max asymmetry',error_max)

      !!call f_free(mat_full)

      call f_release_routine()

    end subroutine max_asymmetry_of_matrix


    !!!!subroutine transform_sparse_matrix_test(iproc, smat, lmat, cmode, &
    !!!!           smat_in, lmat_in, smat_out, lmat_out)
    !!!!  implicit none
    !!!!
    !!!!  ! Calling arguments
    !!!!  integer,intent(in) :: iproc
    !!!!  type(sparse_matrix),intent(in) :: smat, lmat
    !!!!  character(len=14),intent(in) :: cmode
    !!!!  real(kind=8),dimension(smat%nspin*smat%nvctrp_tg),intent(in),optional :: smat_in
    !!!!  real(kind=8),dimension(lmat%nspin*lmat%nvctrp_tg),intent(in),optional :: lmat_in
    !!!!  real(kind=8),dimension(smat%nspin*smat%nvctrp_tg),intent(out),optional :: smat_out
    !!!!  real(kind=8),dimension(lmat%nspin*lmat%nvctrp_tg),intent(out),optional :: lmat_out
    !!!!
    !!!!  ! Local variables
    !!!!  integer(kind=8) :: isstart, isend, ilstart, ilend, iostart, ioend
    !!!!  integer :: imode, icheck, isseg, ilseg
    !!!!  integer :: ilength, iscostart, ilcostart, i
    !!!!  integer :: ilsegstart, ispin, isshift, ilshift, isoffset, iloffset
    !!!!  integer,parameter :: SMALL_TO_LARGE=1
    !!!!  integer,parameter :: LARGE_TO_SMALL=2
    !!!!
    !!!!  call f_routine(id='transform_sparse_matrix')
    !!!!
    !!!!  ! determine the case:
    !!!!  ! SMALL_TO_LARGE -> transform from large sparsity pattern to small one
    !!!!  ! LARGE_TO_SMALL -> transform from small sparsity pattern to large one
    !!!!  if (cmode=='small_to_large' .or. cmode=='SMALL_TO_LARGE') then
    !!!!      imode=SMALL_TO_LARGE
    !!!!      if (.not.present(smat_in)) call f_err_throw('smat_in not present')
    !!!!      if (.not.present(lmat_out)) call f_err_throw('lmat_out not present')
    !!!!  else if (cmode=='large_to_small' .or. cmode=='LARGE_TO_SMALL') then
    !!!!      imode=LARGE_TO_SMALL
    !!!!      if (.not.present(lmat_in)) call f_err_throw('lmat_in not present')
    !!!!      if (.not.present(smat_out)) call f_err_throw('smat_out not present')
    !!!!  else
    !!!!      call f_err_throw('wrong cmode')
    !!!!  end if
    !!!!
    !!!!  select case (imode)
    !!!!  case (SMALL_TO_LARGE)
    !!!!     !call to_zero(lmat%nvctr*lmat%nspin,lmatrix_compr(1))
    !!!!     call f_zero(lmat_out)
    !!!!  case (LARGE_TO_SMALL)
    !!!!     !call to_zero(smat%nvctr*lmat%nspin,smatrix_compr(1))
    !!!!     call f_zero(smat_out)
    !!!!  case default
    !!!!      call f_err_throw('wrong imode')
    !!!!  end select
    !!!!
    !!!!  call timing(iproc,'transform_matr','IR')


    !!!!  icheck=0
    !!!!  do ispin=1,smat%nspin

    !!!!      isshift=(ispin-1)*smat%nvctr
    !!!!      ilshift=(ispin-1)*lmat%nvctr
    !!!!
    !!!!      ilsegstart=lmat%iseseg_tg(1)
    !!!!      !$omp parallel default(private) &
    !!!!      !$omp shared(smat, lmat, imode, icheck, isshift, ilshift) &
    !!!!      !$omp shared(smat_in, lmat_in, smat_out, lmat_out) &
    !!!!      !$omp firstprivate(ilsegstart)
    !!!!      !$omp do reduction(+:icheck)
    !!!!      sloop: do isseg=smat%iseseg_tg(1),smat%iseseg_tg(2)
    !!!!          isstart = int((smat%keyg(1,2,isseg)-1),kind=8)*int(smat%nfvctr,kind=8) + int(smat%keyg(1,1,isseg),kind=8)
    !!!!          isend = int((smat%keyg(2,2,isseg)-1),kind=8)*int(smat%nfvctr,kind=8) + int(smat%keyg(2,1,isseg),kind=8)
    !!!!          ! A segment is always on one line, therefore no double loop
    !!!!          lloop: do ilseg=ilsegstart,lmat%iseseg_tg(2)
    !!!!              ilstart = int((lmat%keyg(1,2,ilseg)-1),kind=8)*int(lmat%nfvctr,kind=8) + int(lmat%keyg(1,1,ilseg),kind=8)
    !!!!              ilend = int((lmat%keyg(2,2,ilseg)-1),kind=8)*int(lmat%nfvctr,kind=8) + int(lmat%keyg(2,1,ilseg),kind=8)
    !!!!
    !!!!              ! check whether there is an overlap:
    !!!!              ! if not, increase loop counters
    !!!!              if (ilstart>isend) then
    !!!!                  !ilsegstart=ilseg
    !!!!                  exit lloop
    !!!!              end if
    !!!!              if (isstart>ilend) then
    !!!!                  ilsegstart=ilseg
    !!!!                  cycle lloop
    !!!!              end if
    !!!!              ! if yes, determine start end end of overlapping segment (in uncompressed form)
    !!!!              iostart=max(isstart,ilstart)
    !!!!              ioend=min(isend,ilend)
    !!!!              ilength=int(ioend-iostart+1,kind=4)
    !!!!
    !!!!              ! offset with respect to the starting point of the segment
    !!!!              isoffset = int(iostart - &
    !!!!                         (int((smat%keyg(1,2,isseg)-1),kind=8)*int(smat%nfvctr,kind=8) &
    !!!!                           + int(smat%keyg(1,1,isseg),kind=8)),kind=4)
    !!!!              iloffset = int(iostart - &
    !!!!                         (int((lmat%keyg(1,2,ilseg)-1),kind=8)*int(lmat%nfvctr,kind=8) &
    !!!!                           + int(lmat%keyg(1,1,ilseg),kind=8)),kind=4)
    !!!!
    !!!!              ! determine start end and of the overlapping segment in compressed form
    !!!!              iscostart=smat%keyv(isseg)+isoffset
    !!!!              ilcostart=lmat%keyv(ilseg)+iloffset
    !!!!
    !!!!              ! copy the elements
    !!!!              select case (imode)
    !!!!              case (SMALL_TO_LARGE) 
    !!!!                  do i=0,ilength-1
    !!!!                      lmat_out(ilcostart+i+ilshift-lmat%isvctrp_tg)=smat_in(iscostart+i+isshift-smat%isvctrp_tg)
    !!!!                  end do
    !!!!              case (LARGE_TO_SMALL) 
    !!!!                  do i=0,ilength-1
    !!!!                      smat_out(iscostart+i+isshift-smat%isvctrp_tg)=lmat_in(ilcostart+i+ilshift-lmat%isvctrp_tg)
    !!!!                  end do
    !!!!              case default
    !!!!                  stop 'wrong imode'
    !!!!              end select
    !!!!              icheck=icheck+ilength
    !!!!          end do lloop
    !!!!      end do sloop
    !!!!      !$omp end do 
    !!!!      !$omp end parallel

    !!!!  end do
    !!!!
    !!!!  ! all elements of the small matrix must have been processed, no matter in
    !!!!  ! which direction the transformation has been executed
    !!!!  if (icheck/=smat%nvctrp_tg*smat%nspin) then
    !!!!      write(*,'(a,2i8)') 'ERROR: icheck/=smat%nvctr*smat%nspin', icheck, smat%nvctr*smat%nspin
    !!!!      stop
    !!!!  end if

    !!!!  call timing(iproc,'transform_matr','RS')
    !!!!  call f_release_routine()
    !!!!
    !!!!end subroutine transform_sparse_matrix_test



    subroutine symmetrize_matrix(smat, csign, mat_in, mat_out, ispinx)
      use sparsematrix_init, only: matrixindex_in_compressed
      use dynamic_memory
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      character(len=*),intent(in) :: csign
      real(mp),dimension(smat%nvctrp_tg*smat%nspin),intent(in) :: mat_in
      real(mp),dimension(smat%nvctrp_tg*smat%nspin),intent(out) :: mat_out
      integer,intent(in),optional :: ispinx

      ! Local variables
      integer :: ispin, ishift, ishift_tg, ii, i, ii_trans
      logical :: minus
      real(mp) :: half
    
      call f_routine(id='symmetrize_matrix')

      if (trim(csign)=='plus') then
          minus = .false.
      else if (trim(csign)=='minus') then
          minus = .true.
      else
          call f_err_throw("wrong value of 'csign'", err_name='BIGDFT_RUNTIME_ERROR')
      end if

      half=0.5_mp
      if (minus) half=-half

      do ispin=1,smat%nspin
          if (present(ispinx)) then
              if (ispin/=ispinx) cycle
          end if
          ishift = (ispin-1)*smat%nvctrp_tg
          ishift_tg = ishift-smat%isvctrp_tg
          !!!!$omp parallel default(none) &
          !!!!$omp shared(smat,mat_in,mat_out,ishift_tg,half) &
          !!!!$omp private(iseg,ii,i,ii_trans)
          !!!!$omp do schedule(guided)
          !!!do iseg=smat%istartendseg_local(1),smat%istartendseg_local(2)
          !!!    ii = smat%keyv(iseg)
          !!!    ! A segment is always on one line, therefore no double loop
          !!!    do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg) !this is too much, but for the moment ok
          !!!        ii_trans = matrixindex_in_compressed(smat,smat%keyg(1,2,iseg),i)
          !!!            mat_out(ii+ishift_tg) = half*(&
          !!!                 mat_in(ii+ishift_tg)+&
          !!!                 mat_in(ii_trans+ishift_tg))
          !!!        ii=ii+1
          !!!    end do
          !!!end do
          !!!!$omp end do
          !!!!$omp end parallel
          !$omp parallel default(none) &
          !$omp shared(smat,mat_in,mat_out,ishift_tg,half) &
          !$omp private(ii,ii_trans)
          !$omp do schedule(guided)
          do ii=smat%istartend_local(1),smat%istartend_local(2)
              ii_trans = smat%transposed_lookup_local(ii)
              mat_out(ii+ishift_tg) = half*(&
                   mat_in(ii+ishift_tg)+&
                   mat_in(ii_trans+ishift_tg))
          end do
          !$omp end do
          !$omp end parallel
      end do

!!$      if (minus) then
!!$          ! There should be a scal wrapper...
!!$          call dscal(smat%nvctrp_tg*smat%nspin, -1.d0, mat_out(1), 1)
!!$      end if
    
      call f_release_routine()
    
    end subroutine symmetrize_matrix


    subroutine check_deviation_from_unity_sparse(iproc, smat, mat, max_error, mean_error)
      use sparsematrix_base, only: sparse_matrix,matrices
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
      real(kind=8),intent(out) :: mean_error, max_error

      ! Local variables
      integer :: iseg, ii, i, irow, icolumn
      real(kind=8) :: error

      call f_routine(id='check_deviation_from_unity_sparse')

      mean_error = 0.d0
      max_error = 0.d0
      do iseg=1,smat%nseg
          ii=smat%keyv(iseg)
          ! A segment is always on one line, therefore no double loop
          do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
             irow = i
             icolumn = smat%keyg(1,2,iseg)
             if (irow==icolumn) then
                 error = abs(mat%matrix_compr(ii)-1.d0)
             else
                 error = abs(mat%matrix_compr(ii))
             end if
             mean_error = mean_error + error
             max_error = max(max_error,error)
             ii=ii+1
         end do
      end do
      mean_error = mean_error/real(smat%nvctr,kind=8)

      !if (iproc==0) then
      !    call yaml_mapping_open('Check the deviation from unity of the operation S^x*S^-x')
      !    call yaml_map('max_error',max_error,fmt='(es10.3)')
      !    call yaml_map('mean_error',mean_error/real(smat%nvctr,kind=8),fmt='(es10.3)')
      !    call yaml_mapping_close()
      !end if

      call f_release_routine()

    end subroutine check_deviation_from_unity_sparse



    !!subroutine operation_using_dense_lapack(iproc, nproc, exp_power, smat_in, mat_in)
    !!  use parallel_linalg, only: dgemm_parallel
    !!  implicit none

    !!  ! Calling arguments
    !!  integer,intent(in) :: iproc, nproc
    !!  real(mp),intent(in) :: exp_power
    !!  type(sparse_matrix),intent(in) :: smat_in
    !!  type(matrices),intent(in) :: mat_in

    !!  ! Local variables
    !!  integer :: blocksize
    !!  real(kind=8),dimension(:,:),allocatable :: mat_in_dense, mat_out_dense
    !!  real(kind=8),dimension(:,:,:),allocatable :: mat_check_accur_dense
    !!  external :: gather_timings

    !!  call f_routine(id='operation_using_dense_lapack')

    !!  blocksize = -100
    !!  mat_in_dense = f_malloc((/smat_in%nfvctr,smat_in%nfvctr/),id='mat_in_dense')
    !!  mat_out_dense = f_malloc((/smat_in%nfvctr,smat_in%nfvctr/),id='mat_out_dense')
    !!  mat_check_accur_dense = f_malloc((/smat_in%nfvctr,smat_in%nfvctr,2/),id='mat_check_accur_dense')
    !!  call uncompress_matrix(iproc, nproc, &
    !!       smat_in, mat_in%matrix_compr, mat_in_dense)

    !!  call f_timing_checkpoint(ctr_name='INIT_CUBIC',mpi_comm=mpiworld(),nproc=nproc,&
    !!       gather_routine=gather_timings)

    !!  call matrix_power_dense(iproc, nproc, comm, blocksize, smat_in%nfvctr, &
    !!       mat_in_dense, exp_power, mat_out_dense)

    !!  call f_timing_checkpoint(ctr_name='CALC_CUBIC',mpi_comm=mpiworld(),nproc=nproc,&
    !!       gather_routine=gather_timings)

    !!  call matrix_power_dense(iproc, nproc, comm, blocksize, smat_in%nfvctr, &
    !!       mat_in_dense, -exp_power, mat_check_accur_dense)
    !!  call dgemm_parallel(iproc, nproc, blocksize, mpi_comm_world, 'n', 'n', &
    !!       smat_in%nfvctr, smat_in%nfvctr, smat_in%nfvctr, &
    !!       1.d0, mat_out_dense(1,1), smat_in%nfvctr, &
    !!       mat_check_accur_dense(1,1,1), smat_in%nfvctr, 0.d0, mat_check_accur_dense(1,1,2), smat_in%nfvctr)
    !!  call check_deviation_from_unity_dense(iproc, smat_in%nfvctr, mat_check_accur_dense(1,1,2))

    !!  call f_timing_checkpoint(ctr_name='CHECK_CUBIC',mpi_comm=mpiworld(),nproc=nproc,&
    !!       gather_routine=gather_timings)

    !!  call f_free(mat_check_accur_dense)
    !!  call f_free(mat_in_dense)
    !!  call f_free(mat_out_dense)

    !!  call f_release_routine()

    !!end subroutine operation_using_dense_lapack



    subroutine matrix_power_dense_lapack(iproc, nproc, comm, blocksize_diag, blocksize_matmul, keep_full_result, &
               exp_power, smat_in, smat_out, mat_in, mat_out, algorithm, overwrite)
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, blocksize_diag, blocksize_matmul
      logical,intent(in) :: keep_full_result
      real(mp),intent(in) :: exp_power
      type(sparse_matrix),intent(in) :: smat_in, smat_out
      type(matrices),intent(inout),target :: mat_in
      type(matrices),intent(out),target :: mat_out
      character(len=*),intent(in),optional :: algorithm
      logical,intent(in),optional :: overwrite

      ! Local variables
      real(kind=8),dimension(:,:),pointer :: mat_in_dense, mat_out_dense
      logical :: full_available, overwrite_

      call f_routine(id='operation_using_dense_lapack')

      ! Check the validity of the arguments
      full_available = .true.
      if (.not.associated(mat_in%matrix)) then
          if (keep_full_result) call f_err_throw('mat_in%matrix must be associated')
          full_available = .false.
      end if
      if (size(mat_in%matrix,1)/=smat_in%nfvctr) then
          if (keep_full_result) call f_err_throw('wrong first dimension of mat_in%matrix')
          full_available = .false.
      end if
      if (size(mat_in%matrix,2)/=smat_in%nfvctr) then
          if (keep_full_result) call f_err_throw('wrong second dimension of mat_in%matrix')
          full_available = .false.
      end if
      if (.not.associated(mat_out%matrix)) then
          if (keep_full_result) call f_err_throw('mat_out%matrix must be associated')
          full_available = .false.
      end if
      if (size(mat_out%matrix,1)/=smat_out%nfvctr) then
          if (keep_full_result) call f_err_throw('wrong first dimension of mat_out%matrix')
          full_available = .false.
      end if
      if (size(mat_out%matrix,2)/=smat_out%nfvctr) then
          if (keep_full_result) call f_err_throw('wrong second dimension of mat_out%matrix')
          full_available = .false.
      end if
      if (size(mat_in%matrix_compr)/=smat_in%nvctrp_tg) then
          call f_err_throw('wrong dimension of mat_in%matrix_compr')
      end if
      if (size(mat_out%matrix_compr)/=smat_out%nvctrp_tg) then
          call f_err_throw('wrong dimension of mat_out%matrix_compr')
      end if

      if (full_available) then
          mat_in_dense => mat_in%matrix(:,:,1)
          mat_out_dense => mat_out%matrix(:,:,1)
      else
          mat_in_dense = f_malloc_ptr((/smat_in%nfvctr,smat_in%nfvctr/),id='mat_in_dense')
          mat_out_dense = f_malloc_ptr((/smat_out%nfvctr,smat_out%nfvctr/),id='mat_out_dense')
      end if
      call uncompress_matrix2(iproc, nproc, comm, &
           smat_in, mat_in%matrix_compr, mat_in_dense)
      overwrite_ = .false.
      if (present(overwrite)) overwrite_ = overwrite
      if (present(algorithm)) then
          call matrix_power_dense(iproc, nproc, comm, blocksize_diag, blocksize_matmul, smat_in%nfvctr, &
               mat_in_dense, exp_power, mat_out_dense, algorithm=algorithm, overwrite=overwrite_)
      else
          call matrix_power_dense(iproc, nproc, comm, blocksize_diag, blocksize_matmul, smat_in%nfvctr, &
               mat_in_dense, exp_power, mat_out_dense, overwrite=overwrite_)
      end if
      call compress_matrix2(iproc, nproc, smat_out, mat_out_dense, mat_out%matrix_compr)
      !!if (keep_full_result) then
      !!    call f_memcpy(src=mat_in_dense, dest=mat_in%matrix)
      !!    call f_memcpy(src=mat_out_dense, dest=mat_out%matrix)
      !!end if
      if (.not.full_available) then
          call f_free_ptr(mat_in_dense)
          call f_free_ptr(mat_out_dense)
      end if

      call f_release_routine()

    end subroutine matrix_power_dense_lapack



    !> Calculate matrix**power, using the dense matrix and exact LAPACK operations
    subroutine matrix_power_dense(iproc, nproc, comm, blocksize_diag, blocksize_matmul, &
               n, mat_in, ex, mat_out, algorithm, overwrite)
      !use module_base
      use parallel_linalg, only: dgemm_parallel, dsyev_parallel
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, blocksize_diag, blocksize_matmul, n
      real(kind=8),dimension(n,n),intent(inout),target :: mat_in
      real(kind=8),intent(in) :: ex
      real(kind=8),dimension(n,n),intent(out) :: mat_out
      character(len=*),intent(in),optional :: algorithm
      logical,intent(in),optional :: overwrite

      ! Local variables
      integer :: i, j, info
      real(kind=8) :: tt
      real(kind=8),dimension(:,:),pointer :: mat_diag
      real(kind=8),dimension(:,:),pointer :: mat_tmp
      real(kind=8),dimension(:),allocatable :: eval
      logical :: positive_definite, overwrite_

      call f_routine(id='matrix_power_dense')

      overwrite_ = .false.
      if (present(overwrite)) overwrite_ = overwrite

      ! Diagonalize the matrix
      if (overwrite_) then
          mat_diag => mat_in
      else
          mat_diag = f_malloc_ptr((/n,n/),id='mat_diag')
          ! f_memcpy can cause segfault for large matrices (I assume integer overflows)
          !call f_memcpy(src=mat_in, dest=mat_tmp)
          do i=1,n
              do j=1,n
                  mat_diag(j,i) = mat_in(j,i)
              end do
          end do
      end if
      eval = f_malloc(n,id='eval')
      if (present(algorithm)) then
          call dsyev_parallel(iproc, nproc, blocksize_diag, comm, 'v', 'l', n, mat_diag, n, eval, info, algorithm=algorithm)
      else
          call dsyev_parallel(iproc, nproc, blocksize_diag, comm, 'v', 'l', n, mat_diag, n, eval, info)
      endif
      if (info /= 0) then
          if (iproc==0) then
              call f_err_throw('dsyev_parallel issued error code '//trim(yaml_toa(info)))
          end if
      end if

      ! Multiply a diagonal matrix containing the eigenvalues to the power ex with the diagonalized matrix.
      ! If all eigenvalues are positive, we can save some memory and take the square root of this value, 
      ! and later on multiply the resulting matrix with its own.
      if (minval(eval)>0.d0 .and. maxval(eval)>0.d0) then
          ! Matrix is positive definite
          positive_definite = .true.
          do i=1,n
              tt = sqrt(eval(i)**ex)
              call dscal(n, tt, mat_diag(1,i), 1)
          end do
      else
          ! Matrix is not positive definite
          positive_definite = .false.
          mat_tmp = f_malloc_ptr((/n,n/),id='mat_tmp')
          do i=1,n
              tt = eval(i)**ex
              do j=1,n
                  mat_tmp(j,i) = mat_diag(j,i)*tt
              end do
          end do
      end if

      call f_free(eval)

      ! Apply the diagonalized matrix to the matrix constructed above
      if (positive_definite) then
          call dgemm_parallel(iproc, nproc, blocksize_matmul, comm, 'n', 't', n, n, n, 1.d0, mat_diag, n, &
               mat_diag, n, 0.d0, mat_out, n)
      else
          call dgemm_parallel(iproc, nproc, blocksize_matmul, comm, 'n', 't', n, n, n, 1.d0, mat_diag, n, &
               mat_tmp, n, 0.d0, mat_out, n)
          call f_free_ptr(mat_tmp)
      end if
      !call dgemm_parallel(iproc, nproc, -1, comm, 'n', 't', n, n, n, 1.d0, mat_tmp(1:,1:,1), n, &
      !     mat_tmp(1:,1:,2), n, 0.d0, mat_out, n)

      if (.not.overwrite_) then
          call f_free_ptr(mat_diag)
      end if

      call f_release_routine()


    end subroutine matrix_power_dense


    subroutine check_deviation_from_unity_dense(iproc, n, mat)
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, n
      real(kind=8),dimension(n,n),intent(in) :: mat

      ! Local variables
      integer :: i, j
      real(kind=8) :: sum_error, max_error, error

      sum_error = 0.d0
      max_error = 0.d0
      do i=1,n
          do j=1,n
              if (j==i) then
                  error = abs(mat(j,i)-1.d0)
              else
                  error = abs(mat(j,i))
              end if
              sum_error = sum_error + error
              max_error = max(max_error,error)
          end do
      end do


      if (iproc==0) then
          call yaml_mapping_open('Check the deviation from unity of the operation S^x*S^-x')
          call yaml_map('max_error',max_error,fmt='(es10.3)')
          call yaml_map('sum_error',sum_error/real(n**2,kind=8),fmt='(es10.3)')
          call yaml_mapping_close()
      end if

    end subroutine check_deviation_from_unity_dense


    subroutine diagonalizeHamiltonian2(iproc, nproc, comm, itype, blocksize, norb, HamSmall, ovrlp, eval)
      use dynamic_memory
      !
      ! Purpose:
      ! ========
      !   Diagonalizes the Hamiltonian HamSmall and makes sure that all MPI processes give
      !   the same result. This is done by requiring that the first entry of each vector
      !   is positive.
      !
      ! Calling arguments:
      ! ==================
      !   Input arguments:
      !   ----------------
      !     iproc     process ID
      !     nproc     number of MPI processes
      !     comm      MPI communicator
      !     itype     type of eigenvalue problem to solve:
      !                 1: A*x = (lambda)*B*x
      !                 2: A*B*x = (lambda)*x
      !                 3: B*A*x = (lambda)*x
      !     blocksize ScaLAPACK block size (negative for sequential LAPACK)
      !     norb      size of the matrices
      !   Input / Putput arguments
      !     HamSmall  on input: the Hamiltonian
      !               on exit: the eigenvectors
      !     ovrlp     on input: overlap matrix
      !               on exit: overwritten
      !   Output arguments
      !     eval      the associated eigenvalues 
      !
      use yaml_output
      use parallel_linalg, only: dsygv_parallel
      implicit none
    
      ! Calling arguments
      integer, intent(in) :: iproc, nproc, comm, itype, blocksize, norb
      real(mp),dimension(norb, norb),intent(inout) :: HamSmall
      real(mp),dimension(norb, norb),intent(inout) :: ovrlp
      real(mp),dimension(norb),intent(out) :: eval
    
      ! Local variables
      integer :: lwork, info
      real(mp),dimension(:),allocatable :: work
      character(len=*),parameter :: subname='diagonalizeHamiltonian'


      call dsygv_parallel(iproc, nproc, comm, blocksize, nproc, &
           itype, 'v', 'l', norb, HamSmall, norb, ovrlp, norb, eval, info)
      if (info/=0) then
          call yaml_warning('dsygv_parallel returned non-zero error code, value = '//trim(yaml_toa(info)))
      end if

      !!!!real(8),dimension(:,:),pointer :: hamtmp, ovrlptmp, invovrlp, tmpmat, tmpmat2
      !!!!real(8) :: tt, tt2
      !!!!integer :: nproc
      !!!!real(8),dimension(norb,norb) :: kernel
    
      !!!!allocate(hamtmp(norb,norb))
      !!!!allocate(ovrlptmp(norb,norb))
      !!!!allocate(invovrlp(norb,norb))
      !!!!allocate(tmpmat(norb,norb))
      !!!!allocate(tmpmat2(norb,norb))
    
      !!!!call mpi_comm_size(mpi_comm_world,nproc,istat)
    
      !!!!hamtmp=HamSmall
      !!!!ovrlptmp=ovrlp
      !!!!call overlapPowerGeneral(iproc, nproc, 100, -2, -1, norb, ovrlptmp, invovrlp, tt)
    
      !!!!call dgemm('n', 'n', norb, norb, norb, 1.d0, invovrlp, norb, hamtmp, norb, 0.d0, tmpmat, norb)
      !!!!call dgemm('n', 'n', norb, norb, norb, 1.d0, tmpmat, norb, invovrlp, norb, 0.d0, tmpmat2, norb)
    
      !!!!lwork=10000
      !!!!allocate(work(lwork))
      !!!!call dsyev('v', 'l', norb, tmpmat2, norb, eval, work, lwork, info)
      !!!!deallocate(work)
    
      !!!!ovrlptmp=ovrlp
      !!!!tmpmat=tmpmat2
      !!!!call overlapPowerGeneral(iproc, nproc, 100, -2, -1, norb, ovrlptmp, invovrlp, tt)
      !!!!!call dgemm('n', 'n', norb, norb, norb, 1.d0, invovrlp, norb, tmpmat, norb, 0.d0, tmpmat2, norb)
      !!!!!if (iproc==0) then
      !!!!!    do istat=1,norb
      !!!!!        do iall=1,norb
      !!!!!            write(200,*) tmpmat2(iall,istat)
      !!!!!        end do
      !!!!!    end do
      !!!!!end if
    
      !!!!call dgemm('n', 't', norb, norb, 28, 1.d0, tmpmat2, norb, tmpmat2, norb, 0.d0, kernel, norb)
      !!!!if (iproc==0) then
      !!!!    tt=0.d0
      !!!!    tt2=0.d0
      !!!!    do istat=1,norb
      !!!!        do iall=1,norb
      !!!!            write(300,*) kernel(iall,istat)
      !!!!            if (istat==iall) tt=tt+kernel(iall,istat)
      !!!!            tt2=tt2+kernel(iall,istat)*ovrlp(iall,istat)
      !!!!        end do
      !!!!    end do
      !!!!    write(*,*) 'Before: trace(K)',tt
      !!!!    write(*,*) 'Before: trace(KS)',tt2
      !!!!end if
    
      !!!!call dgemm('n', 'n', norb, norb, norb, 1.d0, invovrlp, norb, kernel, norb, 0.d0, tmpmat, norb)
      !!!!call dgemm('n', 'n', norb, norb, norb, 1.d0, tmpmat, norb, invovrlp, norb, 0.d0, kernel, norb)
      !!!!if (iproc==0) then
      !!!!    tt=0.d0
      !!!!    tt2=0.d0
      !!!!    do istat=1,norb
      !!!!        do iall=1,norb
      !!!!            write(305,*) kernel(iall,istat)
      !!!!            if (istat==iall) tt=tt+kernel(iall,istat)
      !!!!            tt2=tt2+kernel(iall,istat)*ovrlp(iall,istat)
      !!!!        end do
      !!!!    end do
      !!!!    write(*,*) 'After: trace(K)',tt
      !!!!    write(*,*) 'After: trace(KS)',tt2
      !!!!end if
    
    
      !!!call timing(iproc,'diagonal_seq  ','ON')
      !!call f_timing(TCAT_SMAT_HL_DSYGV,'ON')
      !!call f_routine(id='diagonalizeHamiltonian2')
    
      !!! DEBUG: print hamiltonian and overlap matrices
      !!!if (iproc==0) then
      !!!   open(10)
      !!!   open(11)
      !!!   do iorb=1,orbs%norb
      !!!      do jorb=1,orbs%norb
      !!!         write(10,*) iorb,jorb,HamSmall(iorb,jorb)
      !!!         write(11,*) iorb,jorb,ovrlp(iorb,jorb)
      !!!      end do
      !!!      write(10,*) ''
      !!!      write(11,*) ''
      !!!   end do
      !!!   close(10)
      !!!   close(11)
      !!!end if
      !!! DEBUG: print hamiltonian and overlap matrices
    
      !!!call yaml_map('Hamiltonian before',HamSmall)
      !!! Get the optimal work array size
      !!lwork=-1 
      !!work = f_malloc(100,id='work')
      !!call dsygv(1, 'v', 'l', norb, HamSmall(1,1), norb, ovrlp(1,1), norb, eval(1), work(1), lwork, info) 
      !!lwork=int(work(1))
    
      !!! Deallocate the work array and reallocate it with the optimal size
      !!call f_free(work)
      !!work = f_malloc(lwork,id='work')
    
      !!! Diagonalize the Hamiltonian
      !!call dsygv(1, 'v', 'l', norb, HamSmall(1,1), norb, ovrlp(1,1), norb, eval(1), work(1), lwork, info) 
      !!if(info/=0)then
      !!  write(*,*) 'ERROR: dsygv in diagonalizeHamiltonian2, info=',info,'N=',norb
      !!end if
      !!!!if (iproc==0) then
      !!!!    do istat=1,norb
      !!!!        do iall=1,norb
      !!!!            write(201,*) hamsmall(iall,istat)
      !!!!        end do
      !!!!    end do
      !!!!end if
    
      !!call f_free(work)
    
      !!call f_release_routine()
      !!!call timing(iproc,'diagonal_seq  ','OF')
      !!call f_timing(TCAT_SMAT_HL_DSYGV,'OF')
    
    end subroutine diagonalizeHamiltonian2


    !> Get the minimal and maximal eigenvalue of a matrix
    subroutine get_minmax_eigenvalues(iproc, nproc, comm, mode, scalapack_blocksize, &
               smat, mat, eval_min, eval_max, &
               algorithm, quiet, smat2, mat2, evals)
      use parallel_linalg, only: dsyev_parallel, dsygv_parallel
      use dynamic_memory
      use yaml_output
      implicit none

      ! Calling arguments
      integer, intent(in) :: iproc, nproc, comm, scalapack_blocksize
      character(len=*),intent(in) :: mode
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
      real(mp),dimension(smat%nspin),intent(out) :: eval_min, eval_max
      character(len=*),intent(in),optional :: algorithm
      logical,intent(in),optional :: quiet
      type(sparse_matrix),intent(in),optional :: smat2
      type(matrices),intent(in),optional :: mat2
      real(mp),dimension(smat%nfvctr*smat%nspin),intent(out),optional :: evals

      ! Local variables
      integer :: info, ispin, ishift, imode
      real(kind=mp),dimension(:,:,:),allocatable :: tempmat, tempmat2
      real(kind=mp),dimension(:),allocatable :: eval
      logical :: quiet_

      call f_routine(id='get_minmax_eigenvalues')

      if (trim(mode)=='standard') then
          !Everything ok, standard eigenvalue problem
          imode = 1
      else if (trim(mode)=='generalized') then
          ! Generalized eigenvalue problem, check the presene of the required optional arguments
          if (.not.present(smat2)) call f_err_throw('smat2 must be present for a generalized eigenvalue problem')
          if (.not.present(mat2)) call f_err_throw('mat2 must be present for a generalized eigenvalue problem')
          imode = 2
      else
          call f_err_throw("wrong value for 'mode', must be 'standard' or 'generalized'")
      end if

      quiet_ = .false.
      if (present(quiet)) quiet_ = quiet

      eval = f_malloc(smat%nfvctr,id='eval')

      tempmat = f_malloc0((/smat%nfvctr,smat%nfvctr,smat%nspin/),id='tempmat')
      call uncompress_matrix2(iproc, nproc, comm, smat, mat%matrix_compr, tempmat)
      if (imode==2) then
          tempmat2 = f_malloc0((/smat%nfvctr,smat%nfvctr,smat%nspin/),id='tempmat2')
          call uncompress_matrix2(iproc, nproc, comm, smat2, mat2%matrix_compr, tempmat2)
      end if

      do ispin=1,smat%nspin


          ishift = (ispin-1)*smat%nvctr

          !!do iseg=1,smat%nseg
          !!    ii=smat%keyv(iseg)
          !!    do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
          !!        tempmat(i,smat%keyg(1,2,iseg)) = mat%matrix_compr(ishift+ii)
          !!        ii = ii + 1
          !!    end do
          !!end do

          !!if (imode==2) then
          !!    tempmat2 = f_malloc0((/smat%nfvctr,smat%nfvctr/),id='tempmat2')
          !!    do iseg=1,smat2%nseg
          !!        ii=smat2%keyv(iseg)
          !!        do i=smat2%keyg(1,1,iseg),smat2%keyg(2,1,iseg)
          !!            tempmat2(i,smat2%keyg(1,2,iseg)) = mat2%matrix_compr(ishift+ii)
          !!            ii = ii + 1
          !!        end do
          !!    end do
          !!end if

          if (imode==1) then
              if (present(algorithm)) then
                  call dsyev_parallel(iproc, nproc, scalapack_blocksize, comm, 'n', 'l', &
                       smat%nfvctr, tempmat(:,:,ispin), smat%nfvctr, eval, info, algorithm=algorithm)
              else
                  call dsyev_parallel(iproc, nproc, scalapack_blocksize, comm, 'n', 'l', &
                       smat%nfvctr, tempmat(:,:,ispin), smat%nfvctr, eval, info)
              end if
          else if (imode==2) then
              call dsygv_parallel(iproc, nproc, comm, scalapack_blocksize, nproc, 1, 'n', 'l', &
                   smat%nfvctr, tempmat(:,:,ispin), smat%nfvctr, tempmat2(:,:,ispin), smat%nfvctr, eval, info)
          end if
          if (info/=0) then
              if (iproc==0) then
                  if (imode==1) then
                      call f_err_throw('dsyev_parallel issued error code '//trim(yaml_toa(info)))
                  else if (imode==2) then
                      call f_err_throw('dsygv_parallel issued error code '//trim(yaml_toa(info)))
                  end if
              end if
          end if
          if (iproc==0 .and. .not.quiet_) then
              call yaml_map('eval max/min',(/eval(1),eval(smat%nfvctr)/),fmt='(es16.6)')
          end if
          eval_min(ispin) = eval(1)
          eval_max(ispin) = eval(smat%nfvctr)

          if (present(evals)) then
              call vcopy(smat%nfvctr, eval(1), 1, evals((ispin-1)*smat%nfvctr+1), 1)
          end if

      end do

      call f_free(tempmat)
      if (imode==2) then
          call f_free(tempmat2)
      end if
      call f_free(eval)

      call f_release_routine()

    end subroutine get_minmax_eigenvalues


    subroutine resize_matrix_to_taskgroup(smat, mat)
      use futile
      implicit none
      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(inout) :: mat
      ! Local variables
      real(kind=mp),dimension(:),allocatable :: mat_tg

      call f_routine(id='resize_matrix_to_taskgroup')

      mat_tg = sparsematrix_malloc(smat,iaction=SPARSE_TASKGROUP,id='mat_tg')
      call extract_taskgroup(smat, mat%matrix_compr, mat_tg)
      call f_free_ptr(mat%matrix_compr)
      mat%matrix_compr = sparsematrix_malloc_ptr(smat,iaction=SPARSE_TASKGROUP,id='mat%matrix_compr')
      call f_memcpy(src=mat_tg, dest=mat%matrix_compr)
      call f_free(mat_tg)

      call f_release_routine()

    end subroutine resize_matrix_to_taskgroup

end module sparsematrix
