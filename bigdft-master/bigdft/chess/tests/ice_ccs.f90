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


!> @file
!! Test of the sparsematrix library
!! @author
!!    Copyright (C) 2015-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


program driver_css
  use module_base
  use sparsematrix_base, only: sparse_matrix, matrices, &
                               deallocate_sparse_matrix, deallocate_matrices
  use bigdft_run, only: bigdft_init
  use sparsematrix_highlevel, only: sparse_matrix_and_matrices_init_from_file_ccs, &
                                    sparse_matrix_init_from_file_ccs, matrices_init, &
                                    matrices_set, sparse_matrix_init_from_data, &
                                    ccs_data_from_sparse_matrix, ccs_matrix_write, &
                                    matrix_matrix_multiplication, matrix_chebyshev_expansion
  use utilities, only: get_ccs_data_from_file
  use sparsematrix, only: write_matrix_compressed
  implicit none

  ! Variables
  integer :: i
  integer :: norder_polynomial
  type(sparse_matrix) :: smat1, smat2, smat3
  type(matrices) :: mat1
  type(matrices),dimension(2) :: mat2
  type(matrices),dimension(2) :: mat3
  integer :: nfvctr, nvctr
  integer,dimension(:),pointer :: row_ind, col_ptr

  ! Initialize flib
  call f_lib_initialize()

  ! General initialization, including MPI. Here we have:
  ! bigdft_mpi%iproc is the task ID
  ! bigdft_mpi%nproc is the total number of tasks
  ! PROBLEM: This is in src/modules...
  call bigdft_init()

  ! Read from matrix1.dat and create the type containing the sparse matrix descriptors (smat1) as well as
  ! the type which contains the matrix data (overlap). The matrix element are stored in mat1%matrix_compr.
  call sparse_matrix_and_matrices_init_from_file_ccs('matrix1.dat', bigdft_mpi%iproc, bigdft_mpi%nproc, &
       bigdft_mpi%mpi_comm, smat1, mat1)

  ! Read from matrix2.dat and creates the type containing the sparse matrix descriptors (smat2).
  call sparse_matrix_init_from_file_ccs('matrix2.dat', bigdft_mpi%iproc, bigdft_mpi%nproc, &
       bigdft_mpi%mpi_comm, smat2)

  ! Prepares the type containing the matrix data.
  call matrices_init(smat2, mat2(1))

  ! Calculate the square root of the matrix described by the pair smat1/mat1 and store the result in
  ! smat2/mat2. Attention: The sparsity pattern of smat2 must be contained within smat1.
  ! It is your responsabilty to assure this, the routine does only some minimal checks.
  ! The final result is contained in mat2(1)%matrix_compr.
  norder_polynomial = 30
  call matrix_chebyshev_expansion(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
       1, (/0.5d0/), &
       smat1, smat2, mat1, mat2(1))

  ! Write the result in YAML format to the standard output (required for non-regression tests).
  call write_matrix_compressed('Result of second matrix', smat2, mat2(1))

  ! Create another matrix type, this time directly with the CCS format descriptors.
  ! Get these descriptors from an auxiliary routine using again matrix2.dat
  call get_ccs_data_from_file('matrix2.dat', nfvctr, nvctr, row_ind, col_ptr)
  call sparse_matrix_init_from_data(bigdft_mpi%iproc, bigdft_mpi%nproc, nfvctr, nvctr, row_ind, col_ptr, smat3)

  ! Prepare a new matrix data type.
  call matrices_init(smat2, mat2(2))
  ! Set the contents of this new matrix type, using the above result.
  call matrices_set(smat2, mat2(1)%matrix_compr, mat2(2))

  ! Prepare two new matrix data types
  do i=1,2
      call matrices_init(smat3, mat3(i))
  end do

  ! Calculate at the same time the square root and the inverse square of the matrix described  by the pair smat2/mat2
  ! and store the result in smat3/mat3. The final results are thus in mat3(1)%matrix_compr and mat3(2)%matrix_compr.
  ! The same wraning as above applies.
  norder_polynomial = 40
  call matrix_chebyshev_expansion(bigdft_mpi%iproc, bigdft_mpi%nproc, bigdft_mpi%mpi_comm, &
       2, (/0.5d0,-0.5d0/), &
       smat2, smat3, mat2(2), mat3)

  ! Write the result in YAML format to the standard output (required for non-regression tests).
  call write_matrix_compressed('Result of second matrix', smat3, mat3(1))
  call write_matrix_compressed('Result of third matrix', smat3, mat3(2))

  ! Calculate the CCS descriptors from the sparse_matrix type.
  call ccs_data_from_sparse_matrix(smat3, row_ind, col_ptr)

  ! Write the two matrices to disk, using the CCS format
  call ccs_matrix_write('squareroot.dat', smat3, row_ind, col_ptr, mat3(1))
  call ccs_matrix_write('invsquareroot.dat', smat3, row_ind, col_ptr, mat3(2))

  ! Multiply the two matrices calculated above (i.e. the square root and the inverse square root) and
  ! store the result in mat2. The final result is thus contained in mat2%matrix_compr.
  call matrix_matrix_multiplication(bigdft_mpi%iproc, bigdft_mpi%nproc, smat3, &
       mat3(1), mat3(2), mat2(1))

  ! Write the result of the above multiplication to a file. Since we multiply the square root times the
  ! inverse square root, the result will be the unity matrix.
  call ccs_matrix_write('unity.dat', smat3, row_ind, col_ptr, mat2(1))

  ! Write the result also in YAML format to the standard output (required for non-regression tests).
  call write_matrix_compressed('Result of fourth matrix', smat3, mat2(1))

  ! Deallocate all the sparse matrix descriptrs types
  call deallocate_sparse_matrix(smat1)
  call deallocate_sparse_matrix(smat2)
  call deallocate_sparse_matrix(smat3)

  ! Deallocate all the matrix data types
  call deallocate_matrices(mat1)
  do i=1,2
      call deallocate_matrices(mat2(i))
      call deallocate_matrices(mat3(i))
  end do

  ! Deallocate the CCS format descriptors
  call f_free_ptr(row_ind)
  call f_free_ptr(col_ptr)

  ! Finalize flib
  call f_lib_finalize()


end program driver_css
