!> @file
!!   Test of the matrix power expansion
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


program driver_single
  ! The following module are part of the sparsematrix library
  use sparsematrix_base
  use sparsematrix_highlevel, only: sparse_matrix_and_matrices_init_from_file_bigdft, &
                                    sparse_matrix_init_from_file_bigdft, &
                                    matrices_init, &
                                    matrix_chebyshev_expansion, &
                                    matrix_matrix_multiplication
  use sparsematrix, only: check_deviation_from_unity_sparse, &
                          get_minmax_eigenvalues, matrix_power_dense_lapack
  use sparsematrix_init, only: write_sparsematrix_info
  use dynamic_memory
  use dictionaries
  use yaml_output
  use wrapper_MPI
  use time_profiling

  implicit none

  ! External routines
  external :: gather_timings

  ! Variables
  type(sparse_matrix) :: smat_in, smat_out
  type(matrices) :: mat_in
  type(matrices),dimension(1) :: mat_out
  type(matrices),dimension(3) :: mat_check_accur
  integer :: ierr, nthread, blocksize, iproc, nproc, i
  real(mp) :: exp_power, max_error, mean_error
  real(mp),dimension(1) :: eval_min, eval_max
  real(mp) :: max_error_rel, mean_error_rel, tt, tt_rel
  character(len=200) :: filename_in, filename_out, filename_out_matmul
  type(dictionary), pointer :: dict_timing_info
  !$ integer :: omp_get_max_threads

  ! Initialize flib
  call f_lib_initialize()

  ! MPI initialization; we have:
  ! iproc is the task ID
  ! nproc is the total number of tasks
  call mpiinit()
  iproc=mpirank()
  nproc=mpisize()

  call f_malloc_set_status(iproc=iproc)

  ! Initialize the sparsematrix error handling and timing.
  call sparsematrix_init_errors()
  call sparsematrix_initialize_timing_categories()

  if (iproc==0) then
      call yaml_new_document()
  end if
  call f_timing_reset(filename='time.yaml',master=iproc==0,verbose_mode=.true. .and. nproc>1)


  if (iproc==0) then
      call yaml_scalar('',hfill='~')
      call yaml_scalar('DRIVER FOR MATRIX CHEBYSHEV EXPANSION',hfill='~')
  end if

  if (iproc==0) then
      call yaml_mapping_open('Parallel environment')
      call yaml_map('MPI tasks',nproc)
      nthread = 1
      !$ nthread = omp_get_max_threads()
      call yaml_map('OpenMP threads',nthread)
      call yaml_mapping_close()
  end if



  ! Read in the parameters for the run:
  ! 1) name of the file with the input matrix
  ! 2) name of the file with the descriptors for the output matrix
  ! 3) exponent for the operation (mat**exponent)
  ! 4) the degree of the polynomial that shall be used
  if (iproc==0) then
      read(*,*) filename_in, filename_out, filename_out_matmul, exp_power

      call yaml_mapping_open('Input parameters')
      call yaml_map('File with input matrix',trim(filename_in))
      call yaml_map('File with output matrix descriptors',trim(filename_out))
      call yaml_map('File with output matrix multiplication descriptors',trim(filename_out_matmul))
      call yaml_map('Exponent',exp_power)
      call yaml_mapping_close()
  end if

  ! Send the input parameters to all MPI tasks
  call mpibcast(filename_in, root=0, comm=mpi_comm_world)
  call mpibcast(filename_out, root=0, comm=mpi_comm_world)
  call mpibcast(filename_out_matmul, root=0, comm=mpi_comm_world)
  call mpibcast(exp_power, root=0, comm=mpi_comm_world)

  ! Read the input matrix descriptors and the matrix itself, and create the correpsonding structures
  if (iproc==0) then
      call yaml_comment('Input matrix',hfill='-')
      !call yaml_mapping_open('Input matrix structure')
  end if
  call sparse_matrix_and_matrices_init_from_file_bigdft('serial', trim(filename_in), &
       iproc, nproc, mpi_comm_world, smat_in, mat_in)
  if (iproc==0) then
      call write_sparsematrix_info(smat_in, 'Input matrix')
      call yaml_mapping_close()
  end if
  !!mat_in%matrix_compr(:) = 100.d0*mat_in%matrix_compr(:)

  ! Read the output matrix descriptors, and create the correpsonding structures
  if (iproc==0) then
      call yaml_comment('Output matrix',hfill='-')
      !call yaml_mapping_open('Output matrix structure')
  end if
  call sparse_matrix_init_from_file_bigdft('serial_text', trim(filename_out), &
       iproc, nproc, mpi_comm_world, smat_out, &
       init_matmul=.true., filename_mult=trim(filename_out_matmul))
  if (iproc==0) then
      call write_sparsematrix_info(smat_out, 'Input matrix')
      call yaml_mapping_close()
  end if

  ! Prepare the output matrix structure
  call matrices_init(smat_out, mat_out(1))


  ! Calculate the minimal and maximal eigenvalue, to determine the condition number
  call get_minmax_eigenvalues(iproc, nproc, mpiworld(), 'standard', -1, &
       smat_in, mat_in, eval_min, eval_max, quiet=.true.)
  if (iproc==0) then
      call yaml_comment('Eigenvalue informations',hfill='-')
      call yaml_mapping_open('Eigenvalue properties')
      call yaml_map('Minimal',eval_min)
      call yaml_map('Maximal',eval_max)
      call yaml_map('Condition number',eval_max/eval_min)
      call yaml_mapping_close()
  end if

  if (iproc==0) then
      call yaml_comment('Performing Matrix Chebyshev Expansion',hfill='-')
  end if

  !call timing(mpi_comm_world,'INIT','PR')
  call f_timing_checkpoint(ctr_name='INIT',mpi_comm=mpiworld(),nproc=mpisize(),&
       gather_routine=gather_timings)

  ! Perform the operation mat_out = mat_in**exp_power
  call matrix_chebyshev_expansion(iproc, nproc, mpi_comm_world, &
       1, (/exp_power/), &
       smat_in, smat_out, mat_in, mat_out, npl_auto=.true.)

  !call timing(mpi_comm_world,'CALC_LINEAR','PR')
  call f_timing_checkpoint(ctr_name='CALC_LINEAR',mpi_comm=mpiworld(),nproc=mpisize(),&
       gather_routine=gather_timings)


  ! Now perform a check of the accuracy by calculating the matrix with the inverse power. The previously
  ! calculated matrix times this result should the give the identity.
  ! Prepare the structures needed for the check of the accuracy
  call matrices_init(smat_out, mat_check_accur(1))
  call matrices_init(smat_out, mat_check_accur(2))
  call matrices_init(smat_out, mat_check_accur(3))

  ! Perform the operation mat_out = mat_in**exp_power
  if (iproc==0) then
      call yaml_comment('Performing Matrix Chebyshev Expansion',hfill='-')
  end if
  call matrix_chebyshev_expansion(iproc, nproc, mpi_comm_world, &
       1, (/-exp_power/), &
       smat_in, smat_out, mat_in, mat_check_accur(1), npl_auto=.true.)
  ! Multiply the previously calculated one with this new none. The result should be the identity.
  call matrix_matrix_multiplication(iproc, nproc, smat_out, &
       mat_out(1), mat_check_accur(1), mat_check_accur(2))

  ! Check the deviation from unity
  if (iproc==0) then
      call yaml_comment('Check deviation from Unity')
  end if
  call check_deviation_from_unity_sparse(iproc, smat_out, mat_check_accur(2), &
       max_error, mean_error)
  if (iproc==0) then
      call yaml_mapping_open('Check the deviation from unity of the operation mat^x*mat^-x')
      call yaml_map('max_error',max_error,fmt='(es10.3)')
      call yaml_map('mean_error',mean_error,fmt='(es10.3)')
      call yaml_mapping_close()
  end if


  !call timing(mpi_comm_world,'CHECK_LINEAR','PR')
  call f_timing_checkpoint(ctr_name='CHECK_LINEAR',mpi_comm=mpiworld(),nproc=mpisize(),&
       gather_routine=gather_timings)


  ! Do the operation using exact LAPACK and the dense matrices
  if (iproc==0) then
      call yaml_comment('Do the same calculation using dense LAPACK',hfill='-')
  end if
  !call operation_using_dense_lapack(iproc, nproc, exp_power, smat_in, mat_in)
  call matrix_power_dense_lapack(iproc, nproc, mpiworld(), -1, -1, .false., &
       exp_power, smat_in, smat_out, mat_in, mat_check_accur(3))
  !!!call write_dense_matrix(iproc, nproc, mpi_comm_world, smat, mat_check_accur(1), 'resultchebyshev.dat', binary=.false.)
  !!!call write_dense_matrix(iproc, nproc, mpi_comm_world, smat, mat_check_accur(3), 'resultlapack.dat', binary=.false.)
  max_error = 0.0_mp
  mean_error = 0.0_mp
  max_error_rel = 0.0_mp
  mean_error_rel = 0.0_mp
  do i=1,smat_out%nvctr
      tt = abs(mat_out(1)%matrix_compr(i)-mat_check_accur(3)%matrix_compr(i))
      tt_rel = tt/abs(mat_check_accur(3)%matrix_compr(i))
      mean_error = mean_error + tt
      max_error = max(max_error,tt)
      mean_error_rel = mean_error_rel + tt_rel
      max_error_rel = max(max_error_rel,tt_rel)
  end do
  mean_error = mean_error/real(smat_out%nvctr,kind=8)
  mean_error_rel = mean_error_rel/real(smat_out%nvctr,kind=8)
  if (iproc==0) then
      call yaml_mapping_open('Check the deviation from the exact result using BLAS (only within the sparsity pattern)')
      call yaml_map('max error',max_error,fmt='(es10.3)')
      call yaml_map('mean error',mean_error,fmt='(es10.3)')
      call yaml_map('max error relative',max_error_rel,fmt='(es10.3)')
      call yaml_map('mean error relative',mean_error_rel,fmt='(es10.3)')
      call yaml_mapping_close()
  end if

  !call timing(mpi_comm_world,'CHECK_CUBIC','PR')
  call f_timing_checkpoint(ctr_name='CHECK_CUBIC',mpi_comm=mpiworld(),nproc=mpisize(),&
       gather_routine=gather_timings)

  ! Deallocate all structures
  call deallocate_sparse_matrix(smat_in)
  call deallocate_sparse_matrix(smat_out)
  call deallocate_matrices(mat_in)
  call deallocate_matrices(mat_out(1))
  call deallocate_matrices(mat_check_accur(1))
  call deallocate_matrices(mat_check_accur(2))
  call deallocate_matrices(mat_check_accur(3))

  !call timing(mpi_comm_world,'FINISH','PR')
  call f_timing_checkpoint(ctr_name='FINISH',mpi_comm=mpiworld(),nproc=mpisize(),&
       gather_routine=gather_timings)


  call build_dict_info(dict_timing_info)
  call f_timing_stop(mpi_comm=mpi_comm_world, nproc=nproc, &
       gather_routine=gather_timings, dict_info=dict_timing_info)
  call dict_free(dict_timing_info)

  ! Finalize MPI
  call mpifinalize()

  ! Finalize flib
  call f_lib_finalize()


  !!contains

   !!!!!> Construct the dictionary needed for the timing information
   !!!! subroutine build_dict_info(dict_info)
   !!!!   !use module_base
   !!!!   use dynamic_memory
   !!!!   use dictionaries
   !!!!   implicit none
   !!!!   type(dictionary), pointer :: dict_info
   !!!!   !local variables
   !!!!   integer :: ierr,namelen,nthreads
   !!!!   character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
   !!!!   character(len=MPI_MAX_PROCESSOR_NAME), dimension(:), allocatable :: nodename
   !!!!   type(dictionary), pointer :: dict_tmp
   !!!!   !$ integer :: omp_get_max_threads

   !!!!   call dict_init(dict_info)
!  !!!!bastian: comment out 4 followinf lines for debug purposes (7.12.2014)
   !!!!   !if (DoLastRunThings) then
   !!!!      call f_malloc_dump_status(dict_summary=dict_tmp)
   !!!!      call set(dict_info//'Routines timing and number of calls',dict_tmp)
   !!!!   !end if
   !!!!   nthreads = 0
   !!!!   !$  nthreads=omp_get_max_threads()
   !!!!   call set(dict_info//'CPU parallelism'//'MPI tasks',nproc)
   !!!!   if (nthreads /= 0) call set(dict_info//'CPU parallelism'//'OMP threads',&
   !!!!        nthreads)

   !!!!   nodename=f_malloc0_str(MPI_MAX_PROCESSOR_NAME,0.to.nproc-1,id='nodename')
   !!!!   if (nproc>1) then
   !!!!      call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)
   !!!!      !gather the result between all the process
   !!!!      call MPI_GATHER(nodename_local,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
   !!!!           nodename(0),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,0,&
   !!!!           mpi_comm_world,ierr)
   !!!!      if (iproc==0) call set(dict_info//'Hostnames',&
   !!!!              list_new(.item. nodename))
   !!!!   end if
   !!!!   call f_free_str(MPI_MAX_PROCESSOR_NAME,nodename)

   !!!! end subroutine build_dict_info





end program driver_single
