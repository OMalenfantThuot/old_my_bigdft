!> @file
!!  
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


module sparsematrix_timing
  use time_profiling, only: TIMING_UNINITIALIZED
  implicit none

  private

  ! Timings categories
  integer,public,save :: TCAT_SMAT_COMPRESSION = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_COMPRESSION_COMMUNICATION = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_TRANSFORMATION = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_MULTIPLICATION = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_INITIALIZATION = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_CME_AUXILIARY = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_CME_POLYNOMIALS = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_CME_COEFFICIENTS = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_HL_MATRIX_OPERATIONS = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_HL_MATRIX_COMMUNICATIONS = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_HL_MATRIX_CHECKS = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_HL_DGEMM = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_HL_DSYEV = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_HL_DSYGV = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_HL_DGESV = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_HL_DPOTRF = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_HL_DPOTRI = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_PEXSI_INIT = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_PEXSI_DISTRIBUTE = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_PEXSI_KERNEL = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_PEXSI_FACTORIZATION = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_PEXSI_INVERSE = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_PEXSI_RETRIEVE = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_PEXSI_FINALIZE = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_PEXSI_COMMUNICATE = TIMING_UNINITIALIZED
  integer,public,save :: TCAT_SMAT_READ = TIMING_UNINITIALIZED

  !> Public routines
  public :: sparsematrix_initialize_timing_categories
  public :: build_dict_info

  contains

    !> Switch on the timing categories for the sparse matrices
    subroutine sparsematrix_initialize_timing_categories()
      use time_profiling, only: f_timing_category_group,f_timing_category
      implicit none
      character(len=*), parameter :: smat_manip = 'sparse matrix manipulation'
      character(len=*), parameter :: smat_init = 'sparse matrix initialization'
      character(len=*), parameter :: smat_comm = 'sparse matrix communications'
      character(len=*), parameter :: cme = 'chebyshev matrix expansion'
      character(len=*), parameter :: smat_blaslapack = 'BLAS / LAPACK'
      character(len=*), parameter :: smat_pexsi = 'PEXSI'
      character(len=*), parameter :: smat_io = 'sparse matrix I/O'
  
  
      ! Definition of the groups
      call f_timing_category_group(smat_manip, 'sparse matrix operations')
      call f_timing_category_group(smat_init, 'sparse matrix initialization')
      call f_timing_category_group(smat_comm, 'sparse matrix communications')
      call f_timing_category_group(cme, 'chebyshev matrix expansion')
      call f_timing_category_group(smat_blaslapack, 'BLAS/LAPACK')
      call f_timing_category_group(smat_pexsi, 'PEXSI')
      call f_timing_category_group(smat_io, 'sparse matrix I/O')
  
      ! Define the timing categories
  
      ! Initialization timing
      call f_timing_category('sparse matrix initialization', smat_init, &
           'sparse matrix initialization', TCAT_SMAT_INITIALIZATION)
  
      ! Low level timings
      call f_timing_category('Sparse matrix compression', smat_manip, &
           '(un)compression of sparse matrices', TCAT_SMAT_COMPRESSION)
      call f_timing_category('Sparse matrix compression communication', smat_comm, &
           '(un)compression communication of sparse matrices', TCAT_SMAT_COMPRESSION_COMMUNICATION)
      call f_timing_category('Sparse matrix transformation', smat_manip, &
           'sparsity pattern transformation of sparse matrices', TCAT_SMAT_TRANSFORMATION)
      call f_timing_category('Sparse matrix multiplication', smat_manip, &
           'sparse matrix matrix multiplication', TCAT_SMAT_MULTIPLICATION)

      ! I/O operations
      call f_timing_category('Sparse matrix read', smat_io, &
           'Reading of the sparse matrices', TCAT_SMAT_READ)
  
      ! Chebyshev Matrix Expansion timing
      call f_timing_category('CME auxiliary', cme, &
           'Chebyshev matrix expansion auxiliary', TCAT_CME_AUXILIARY)
      call f_timing_category('CME polynomials', cme, &
           'Chebyshev matrix expansion polynomials', TCAT_CME_POLYNOMIALS)
      call f_timing_category('CME coefficients', cme, &
           'Chebyshev matrix expansion coefficients', TCAT_CME_COEFFICIENTS)
  
      ! High level timings
      call f_timing_category('highlevel matrix operations', smat_manip, &
           'highlevel matrix operations', TCAT_HL_MATRIX_OPERATIONS)
      call f_timing_category('highlevel matrix communications', smat_comm, &
           'highlevel matrix communications', TCAT_HL_MATRIX_COMMUNICATIONS)
      call f_timing_category('highlevel matrix checks', smat_manip, &
           'highlevel matrix checks', TCAT_HL_MATRIX_CHECKS)
  
      ! BLAS / LAPACK timing
      call f_timing_category('DGEMM', smat_blaslapack, &
           '(Sca)LAPACK DGEMM', TCAT_HL_DGEMM)
      call f_timing_category('DSYEV', smat_blaslapack, &
           '(Sca)LAPACK DSYEV', TCAT_SMAT_HL_DSYEV)
      call f_timing_category('DSYGV', smat_blaslapack, &
           '(Sca)LAPACK DSYGV', TCAT_SMAT_HL_DSYGV)
      call f_timing_category('DGESV', smat_blaslapack, &
           '(Sca)LAPACK DGESV', TCAT_SMAT_HL_DGESV)
      call f_timing_category('DPOTRF', smat_blaslapack, &
           '(Sca)LAPACK DPOTRF', TCAT_SMAT_HL_DPOTRF)
      call f_timing_category('DPOTRI', smat_blaslapack, &
           '(Sca)LAPACK DPOTRI', TCAT_SMAT_HL_DPOTRI)

      ! PEXSI timing
      call f_timing_category('PEXSI Init', smat_pexsi, &
           'PEXSI initialization', TCAT_PEXSI_INIT)
      call f_timing_category('PEXSI Distribute', smat_pexsi, &
           'PEXSI matrix distribution', TCAT_PEXSI_DISTRIBUTE)
      call f_timing_category('PEXSI Kernel', smat_pexsi, &
           'PEXSI kernel calculation', TCAT_PEXSI_KERNEL)
      call f_timing_category('PEXSI Factorization', smat_pexsi, &
           'PEXSI factorization', TCAT_PEXSI_FACTORIZATION)
      call f_timing_category('PEXSI Inverse', smat_pexsi, &
           'PEXSI inverse calculation', TCAT_PEXSI_INVERSE)
      call f_timing_category('PEXSI Retrieve', smat_pexsi, &
           'PEXSI matrix retrieval', TCAT_PEXSI_RETRIEVE)
      call f_timing_category('PEXSI Finalize', smat_pexsi, &
           'PEXSI finalize', TCAT_PEXSI_FINALIZE)
      call f_timing_category('PEXSI Communicate', smat_pexsi, &
           'PEXSI communicate', TCAT_PEXSI_COMMUNICATE)
  
    end subroutine sparsematrix_initialize_timing_categories


    !> construct the dictionary needed for the timing information
    subroutine build_dict_info(dict_info)
      use wrapper_MPI
      use dynamic_memory
      use dictionaries
      implicit none
      !calling arguments
      !integer,intent(in) :: iproc, nproc
      type(dictionary), pointer :: dict_info
 
      call dict_init(dict_info)
      call f_malloc_dump_status(dict_summary=dict_info)
      call mpi_environment_dict(mpi_environment_comm(),dict_info)
 
    end subroutine build_dict_info

end module sparsematrix_timing
