!> @file
!!   Test of the eigenvalue calculation (using FOE) using the CCS or the BigDFT format
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


program driver_eigenvalues
  ! The following module are part of the sparsematrix library
  use sparsematrix_base
  use foe_base, only: foe_data, foe_data_deallocate, foe_data_get_real
  use foe_common, only: init_foe
  use sparsematrix_highlevel, only: sparse_matrix_and_matrices_init_from_file_ccs, &
                                    sparse_matrix_init_from_file_ccs, matrices_init, &
                                    matrices_get_values, matrices_set_values, &
                                    sparse_matrix_init_from_data_ccs, &
                                    ccs_data_from_sparse_matrix, ccs_matrix_write, &
                                    matrix_matrix_multiplication, matrix_fermi_operator_expansion, &
                                    trace_A, trace_AB, sparse_matrix_metadata_init_from_file, &
                                    sparse_matrix_and_matrices_init_from_file_bigdft
  use sparsematrix, only: write_matrix_compressed, transform_sparse_matrix, get_minmax_eigenvalues
  use sparsematrix_init, only: matrixindex_in_compressed, write_sparsematrix_info, &
                               get_number_of_electrons, distribute_on_tasks
  ! The following module is an auxiliary module for this test
  use utilities, only: get_ccs_data_from_file, calculate_error
  use futile
  use wrapper_MPI
  use wrapper_linalg
  use pexsi, only: pexsi_wrapper
  use coeffs, only: get_coeffs_diagonalization, calculate_kernel_and_energy
  use fermi_level, only: eval_to_occ, SMEARING_DIST_ERF
  use highlevel_wrappers, only: calculate_eigenvalues, solve_eigensystem_lapack

  implicit none

  ! Variables
  type(sparse_matrix) :: smat_s, smat_h, smat_k
  type(matrices) :: mat_s, mat_h, mat_k, mat_ek
  type(matrices),dimension(1) :: mat_ovrlpminusonehalf
  type(sparse_matrix_metadata) :: smmd
  integer :: nfvctr, nvctr, ierr, iproc, nproc, nthread, ncharge, nfvctr_mult, nvctr_mult, scalapack_blocksize, icheck
  integer :: ispin, ihomo, imax, ntemp, npl_max, pexsi_npoles, norbu, norbd, ii, info, norbp, isorb, norb, iorb, pexsi_np_sym_fact
  integer :: iev_min, iev_max, iev, iev_minx, iev_maxx
  real(mp) :: pexsi_mumin, pexsi_mumax, pexsi_mu, pexsi_DeltaE, pexsi_temperature, pexsi_tol_charge, fscale
  integer,dimension(:),pointer :: row_ind, col_ptr, row_ind_mult, col_ptr_mult
  real(mp),dimension(:),pointer :: kernel, overlap, overlap_large, evals, evals_check
  real(mp),dimension(:),allocatable :: charge, eval_min, eval_max, eval_all, eval_occup, occup
  real(mp),dimension(:,:),allocatable :: coeff
  real(mp) :: energy, tr_KS, tr_KS_check, ef, energy_fake, efermi, eTS
  type(foe_data) :: foe_obj, ice_obj
  real(mp) :: tr, fscale_lowerbound, fscale_upperbound
  type(dictionary),pointer :: dict_timing_info, options
  type(yaml_cl_parse) :: parser !< command line parser
  character(len=1024) :: metadata_file, overlap_file, hamiltonian_file, kernel_file, kernel_matmul_file
  character(len=1024) :: sparsity_format, matrix_format, kernel_method
  logical :: check_spectrum, do_cubic_check
  integer,parameter :: nthreshold = 10 !< number of checks with threshold
  real(mp),dimension(nthreshold),parameter :: threshold = (/ 1.e-1_mp, &
                                                             1.e-2_mp, &
                                                             1.e-3_mp, &
                                                             1.e-4_mp, &
                                                             1.e-5_mp, &
                                                             1.e-6_mp, &
                                                             1.e-7_mp, &
                                                             1.e-8_mp,&
                                                             1.e-9_mp,&
                                                             1.e-10_mp /) !< threshold for the relative errror

  external :: gather_timings
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
      !call print_logo()
  end if

  !Time initialization
  call f_timing_reset(filename='time.yaml',master=(iproc==0),verbose_mode=.false.)

  if (iproc==0) then
      call yaml_scalar('',hfill='~')
      call yaml_scalar('CHESS EIGENVALUES TEST DRIVER',hfill='~')
  end if

  if (iproc==0) then
      call yaml_map('Timestamp of the run',yaml_date_and_time_toa())
      call yaml_mapping_open('Parallel environment')
      call yaml_map('MPI tasks',nproc)
      nthread = 1
      !$ nthread = omp_get_max_threads()
      call yaml_map('OpenMP threads',nthread)
      call yaml_mapping_close()
  end if

  ! Read in the parameters for the run and print them.
  if (iproc==0) then
      parser=yaml_cl_parse_null()
      call commandline_options(parser)
      call yaml_cl_parse_cmd_line(parser,args=options)
      call yaml_cl_parse_free(parser)

      metadata_file = options//'metadata_file'
      overlap_file = options//'overlap_file'
      hamiltonian_file = options//'hamiltonian_file'
      kernel_file = options//'kernel_file'
      kernel_matmul_file = options//'kernel_matmul_file'
      sparsity_format = options//'sparsity_format'
      matrix_format = options//'matrix_format'
      kernel_method = options//'kernel_method'
      scalapack_blocksize = options//'scalapack_blocksize'
      check_spectrum = options//'check_spectrum'
      fscale_lowerbound = options//'fscale_lowerbound'
      fscale_upperbound = options//'fscale_upperbound'
      ntemp = options//'ntemp'
      ef = options//'ef'
      npl_max = options//'npl_max'
      pexsi_npoles = options//'pexsi_npoles'
      pexsi_mumin = options//'pexsi_mumin'
      pexsi_mumax = options//'pexsi_mumax'
      pexsi_mu = options//'pexsi_mu'
      pexsi_temperature = options//'pexsi_temperature'
      pexsi_tol_charge = options//'pexsi_tol_charge'
      pexsi_np_sym_fact = options//'pexsi_np_sym_fact'
      pexsi_DeltaE = options//'pexsi_DeltaE'
      do_cubic_check = options//'do_cubic_check'
      fscale = options//'fscale'
      iev_min = options//'iev_min'
      iev_max = options//'iev_max'
     
      call dict_free(options)

      call yaml_mapping_open('Input parameters')
      call yaml_map('Sparsity format',trim(sparsity_format))
      call yaml_map('Matrix format',trim(matrix_format))
      call yaml_map('Metadata file',trim(metadata_file))
      call yaml_map('Overlap matrix file',trim(overlap_file))
      call yaml_map('Hamiltonian matrix file',trim(hamiltonian_file))
      call yaml_map('Density kernel matrix file',trim(kernel_file))
      call yaml_map('Density kernel matrix multiplication file',trim(kernel_matmul_file))
      call yaml_map('Kernel method',trim(kernel_method))
      call yaml_map('Blocksize for ScaLAPACK',scalapack_blocksize)
      call yaml_map('Check the Hamiltonian spectrum',check_spectrum)
      call yaml_map('Lower bound for decay length',fscale_lowerbound)
      call yaml_map('Upper bound for decay length',fscale_upperbound)
      call yaml_map('Iterations with varying temperatures',ntemp)
      call yaml_map('Guess for Fermi energy',ef)
      call yaml_map('Maximal polynomial degree',npl_max)
      call yaml_map('PEXSI number of poles',pexsi_npoles)
      call yaml_map('PEXSI mu min',pexsi_mumin)
      call yaml_map('PEXSI mu max',pexsi_mumax)
      call yaml_map('PEXSI mu',pexsi_mu)
      call yaml_map('PEXSI Delta E',pexsi_DeltaE)
      call yaml_map('PEXSI temperature',pexsi_temperature)
      call yaml_map('PEXSI charge tolerance',pexsi_tol_charge)
      call yaml_map('PEXSI number of procs for symbolic factorization',pexsi_np_sym_fact)
      call yaml_map('Do a check with cubic scaling (Sca)LAPACK',do_cubic_check)
      call yaml_map('Minimal eigenvalue to be calculated',iev_min)
      call yaml_map('Maximal eigenvalue to be calculated',iev_max)
      call yaml_map('Decay length of the error function used to find the eigenvalue',fscale)
      call yaml_mapping_close()
  end if

  ! Send the input parameters to all MPI tasks
  call mpibcast(sparsity_format, root=0, comm=mpi_comm_world)
  call mpibcast(matrix_format, root=0, comm=mpi_comm_world)
  call mpibcast(metadata_file, root=0, comm=mpi_comm_world)
  call mpibcast(overlap_file, root=0, comm=mpi_comm_world)
  call mpibcast(hamiltonian_file, root=0, comm=mpi_comm_world)
  call mpibcast(kernel_file, root=0, comm=mpi_comm_world)
  call mpibcast(kernel_matmul_file, root=0, comm=mpi_comm_world)
  call mpibcast(kernel_method, root=0, comm=mpi_comm_world)
  call mpibcast(scalapack_blocksize, root=0, comm=mpi_comm_world)
  call mpibcast(kernel_matmul_file, root=0, comm=mpi_comm_world)
  call mpibcast(fscale_lowerbound, root=0, comm=mpi_comm_world)
  call mpibcast(fscale_upperbound, root=0, comm=mpi_comm_world)
  call mpibcast(ntemp, root=0, comm=mpi_comm_world)
  call mpibcast(ef, root=0, comm=mpi_comm_world)
  call mpibcast(npl_max, root=0, comm=mpi_comm_world)
  call mpibcast(pexsi_npoles, root=0, comm=mpi_comm_world)
  call mpibcast(pexsi_mumin, root=0, comm=mpi_comm_world)
  call mpibcast(pexsi_mumax, root=0, comm=mpi_comm_world)
  call mpibcast(pexsi_mu, root=0, comm=mpi_comm_world)
  call mpibcast(pexsi_DeltaE, root=0, comm=mpi_comm_world)
  call mpibcast(pexsi_temperature, root=0, comm=mpi_comm_world)
  call mpibcast(pexsi_tol_charge, root=0, comm=mpi_comm_world)
  call mpibcast(pexsi_np_sym_fact, root=0, comm=mpi_comm_world)
  call mpibcast(iev_min, root=0, comm=mpi_comm_world)
  call mpibcast(iev_max, root=0, comm=mpi_comm_world)
  call mpibcast(fscale, root=0, comm=mpi_comm_world)
  ! Since there is no wrapper for logicals...
  if (iproc==0) then
      if (check_spectrum) then
          icheck = 1
      else
          icheck = 0
      end if
  end if
  call mpibcast(icheck, root=0, comm=mpi_comm_world)
  if (icheck==1) then
      check_spectrum = .true.
  else
      check_spectrum = .false.
  end if
  if (iproc==0) then
      if (do_cubic_check) then
          icheck = 1
      else
          icheck = 0
      end if
  end if
  call mpibcast(icheck, root=0, comm=mpi_comm_world)
  if (icheck==1) then
      do_cubic_check = .true.
  else
      do_cubic_check = .false.
  end if



  !!! Read in the overlap matrix and create the type containing the sparse matrix descriptors (smat_s) as well as
  !!! the type which contains the matrix data (overlap). The matrix element are stored in mat_s%matrix_compr.
  !!! Do the same also for the Hamiltonian matrix.
  !!if (trim(sparsity_format)=='ccs') then
  !!    if (trim(matrix_format)/='serial_text') then
  !!        call f_err_throw('Wrong matrix format')
  !!    end if
  !!    !if (iproc==0) then
  !!    !    call yaml_scalar('Initializing overlap matrix',hfill='-')
  !!    !    call yaml_map('Reading from file','overlap_ccs.txt')
  !!    !end if
  !!    call sparse_matrix_and_matrices_init_from_file_ccs(trim(overlap_file), &
  !!         iproc, nproc, mpi_comm_world, smat_s, mat_s)

  !!    !if (iproc==0) then
  !!    !    call yaml_scalar('Initializing Hamiltonian matrix',hfill='-')
  !!    !    call yaml_map('Reading from file','hamiltonian_ccs.txt')
  !!    !end if
  !!    call sparse_matrix_and_matrices_init_from_file_ccs(trim(hamiltonian_file), &
  !!         iproc, nproc, mpi_comm_world, smat_h, mat_h)

  !!    ! Create another matrix type, this time directly with the CCS format descriptors.
  !!    ! Get these descriptors from an auxiliary routine using matrix3.dat
  !!    !if (iproc==0) then
  !!    !    call yaml_scalar('Initializing density kernel matrix',hfill='-')
  !!    !    call yaml_map('Reading from file','density_kernel_ccs.txt')
  !!    !end if
  !!    call get_ccs_data_from_file(trim(kernel_file), nfvctr, nvctr, row_ind, col_ptr)
  !!    !call get_ccs_data_from_file('density_kernel_matmul_ccs.dat', nfvctr_mult, nvctr_mult, row_ind_mult, col_ptr_mult)
  !!    call get_ccs_data_from_file(trim(kernel_matmul_file), nfvctr_mult, nvctr_mult, row_ind_mult, col_ptr_mult)
  !!    if (nfvctr_mult/=nfvctr) then
  !!        call f_err_throw('nfvctr_mult/=nfvctr',err_name='SPARSEMATRIX_INITIALIZATION_ERROR')
  !!    end if
  !!    call sparse_matrix_init_from_data_ccs(iproc, nproc, mpi_comm_world, &
  !!         nfvctr, nvctr, row_ind, col_ptr, smat_k, &
  !!         init_matmul=.true., nvctr_mult=nvctr_mult, row_ind_mult=row_ind_mult, col_ptr_mult=col_ptr_mult)
  !!    call f_free_ptr(row_ind)
  !!    call f_free_ptr(col_ptr)
  !!    call f_free_ptr(row_ind_mult)
  !!    call f_free_ptr(col_ptr_mult)
  !!else if (trim(sparsity_format)=='bigdft') then
  !!    call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, overlap_file, &
  !!         iproc, nproc, mpi_comm_world, smat_s, mat_s, init_matmul=.false.)
  !!    call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, hamiltonian_file, &
  !!         iproc, nproc, mpi_comm_world, smat_h, mat_h, init_matmul=.false.)
  !!    call sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, kernel_file, &
  !!         iproc, nproc, mpi_comm_world, smat_k, mat_k, init_matmul=.true., filename_mult=trim(kernel_matmul_file))
  !!else
  !!    call f_err_throw('Wrong sparsity format')
  !!end if

  !!if (iproc==0) then
  !!    call yaml_mapping_open('Matrix properties')
  !!    call write_sparsematrix_info(smat_s, 'Overlap matrix')
  !!    call write_sparsematrix_info(smat_h, 'Hamiltonian matrix')
  !!    call write_sparsematrix_info(smat_k, 'Density kernel')
  !!    call yaml_mapping_close()
  !!end if

  !!call sparse_matrix_metadata_init_from_file(trim(metadata_file), smmd)
  !!call get_number_of_electrons(smmd, ncharge)
  !!if (iproc==0) then
  !!    call yaml_map('Number of electrons',ncharge)
  !!end if

  !!! Prepares the type containing the matrix data.
  !!if (trim(sparsity_format)=='ccs') then
  !!    ! Otherwise already done above
  !!    call matrices_init(smat_k, mat_k)
  !!end if
  !!call matrices_init(smat_k, mat_ek)
  !!call matrices_init(smat_k, mat_ovrlpminusonehalf(1))

  !!! Initialize the opaque object holding the parameters required for the Fermi Operator Expansion.
  !!! Only provide the mandatory values and take for the optional values the default ones.
  !!charge = f_malloc(smat_s%nspin,id='charge')
  !!charge(:) = real(ncharge,kind=mp)
  !!call init_foe(iproc, nproc, smat_s%nspin, charge, foe_obj, &
  !!     fscale_lowerbound=fscale_lowerbound, fscale_upperbound=fscale_upperbound, &
  !!     ntemp = ntemp, ef=ef, npl_max=npl_max)
  !!! Initialize the same object for the calculation of the inverse. Charge does not really make sense here...
  !!call init_foe(iproc, nproc, smat_s%nspin, charge, ice_obj, evlow=0.5_mp, evhigh=1.5_mp)

  !!call mpibarrier()
  !!call f_timing_checkpoint(ctr_name='INIT',mpi_comm=mpiworld(),nproc=mpisize(), &
  !!     gather_routine=gather_timings)

  !!if(check_spectrum) then
  !!    evals = f_malloc(smat_h%nfvctr,id='evals')
  !!    eval_min = f_malloc(smat_h%nspin,id='eval_min')
  !!    eval_max = f_malloc(smat_h%nspin,id='eval_max')
  !!    call get_minmax_eigenvalues(iproc, nproc, mpiworld(), 'generalized', scalapack_blocksize, &
  !!         smat_h, mat_h, eval_min, eval_max, quiet=.true., smat2=smat_s, mat2=mat_s, evals=evals)
  !!    if (iproc==0) then
  !!        call yaml_mapping_open('Hamiltonian spectrum')
  !!        do ispin=1,smat_h%nspin
  !!            if (smat_h%nspin>1) then
  !!                if (ispin==1) then
  !!                    call yaml_mapping_open('Spin up')
  !!                else if (ispin==2) then
  !!                    call yaml_mapping_open('Spin down')
  !!                end if
  !!                ihomo = ncharge
  !!            else
  !!                ihomo = ncharge/2
  !!            end if
  !!            call yaml_map('Lowest eigenvalue',evals(1))
  !!            call yaml_map('Highest eigenvalue',evals(smat_h%nfvctr))
  !!            call yaml_map('HOMO eigenvalue',evals(ihomo))
  !!            call yaml_map('LUMO eigenvalue',evals(ihomo+1))
  !!            call yaml_map('Spectral width',evals(smat_h%nfvctr)-evals(1))
  !!            call yaml_map('HOMO-lUMO gap',evals(ihomo+1)-evals(ihomo))
  !!            if (smat_h%nspin>1) then
  !!                call yaml_mapping_close()
  !!            end if
  !!            call yaml_mapping_close()
  !!        end do
  !!    end if
  !!    call f_free(evals)
  !!    call f_free(eval_min)
  !!    call f_free(eval_max)
  !!end if

  !!call mpibarrier()
  !!call f_timing_checkpoint(ctr_name='INFO',mpi_comm=mpiworld(),nproc=mpisize(), &
  !!     gather_routine=gather_timings)


  ! Calculate the eigenvalues
  call calculate_eigenvalues(iproc, nproc, matrix_format, metadata_file, &
       overlap_file, hamiltonian_file, kernel_file, kernel_matmul_file, &
       1, iev_min, iev_max, fscale, evals_out=evals)

  call mpibarrier()
  call f_timing_checkpoint(ctr_name='CALC',mpi_comm=mpiworld(),nproc=mpisize(), &
       gather_routine=gather_timings)

  if (do_cubic_check) then
      !!norbu = smat_h%nfvctr
      !!norbd = 0
      !!norb = norbu + norbd
      !!coeff = f_malloc((/smat_h%nfvctr,norb/),id='coeff')
      !!eval_all = f_malloc(smat_h%nfvctr,id='eval_all')
      !!eval_occup = f_malloc(norb,id='eval_occup')
      !!write(*,*) 'calling check...'
      !!write(*,*) 'smat_s%nvctrp_tg, smat_s%nspin', smat_s%nvctrp_tg, smat_s%nspin
      !!call get_coeffs_diagonalization(iproc, nproc, mpi_comm_world, &
      !!     smat_h%nfvctr, norbu, norbd, norb, scalapack_blocksize, &
      !!     smat_s, smat_h, mat_s, mat_h, coeff, &
      !!     eval_all, eval_occup, info)
      call solve_eigensystem_lapack(iproc, nproc, mpiworld(), 1, matrix_format, metadata_file, &
           overlap_file, hamiltonian_file, scalapack_blocksize, write_coeff=.false., write_eval=.false., &
           evals_out=evals_check)
      !!call yaml_map('evals',evals_check)

      call f_timing_checkpoint(ctr_name='CHECK',mpi_comm=mpiworld(),nproc=mpisize(), &
           gather_routine=gather_timings)

      if (iproc==0) then
          call yaml_sequence_open('Comparing the eigenvalues')
          iev_minx = lbound(evals,1)
          iev_maxx = ubound(evals,1)
          do iev=iev_minx,iev_maxx
              call yaml_sequence(advance='no')
              call yaml_mapping_open(flow=.true.)
              call yaml_map('ID',iev)
              call yaml_map('CheSS',evals(iev),fmt='(es14.7)')
              call yaml_map('LAPACK',evals_check(iev),fmt='(es14.7)')
              call yaml_map('difference',evals(iev)-evals_check(iev),fmt='(es15.8)')
              call yaml_mapping_close()
          end do
          call yaml_sequence_close
      end if
      call f_free_ptr(evals_check)
  end if

  call f_free_ptr(evals)

  !!! Calculate the density kernel for the system described by the pair smat_s/mat_s and smat_h/mat_h and 
  !!! store the result in smat_k/mat_k.
  !!! Attention: The sparsity pattern of smat_s must be contained within that of smat_h
  !!! and the one of smat_h within that of smat_k. It is your responsabilty to assure this, 
  !!! the routine does only some minimal checks.
  !!! The final result will be contained in mat_k%matrix_compr.
  !!if (trim(kernel_method)=='FOE') then
  !!    call matrix_fermi_operator_expansion(iproc, nproc, mpi_comm_world, &
  !!         foe_obj, ice_obj, smat_s, smat_h, smat_k, &
  !!         mat_s, mat_h, mat_ovrlpminusonehalf, mat_k, energy, &
  !!         calculate_minusonehalf=.true., foe_verbosity=1, symmetrize_kernel=.true., &
  !!         calculate_energy_density_kernel=.true., energy_kernel=mat_ek)
  !!else if (trim(kernel_method)=='PEXSI') then
  !!    call pexsi_wrapper(iproc, nproc, mpi_comm_world, smat_s, smat_h, smat_k, mat_s, mat_h, &
  !!         foe_data_get_real(foe_obj,"charge",1), pexsi_npoles, pexsi_mumin, pexsi_mumax, pexsi_mu, pexsi_DeltaE, &
  !!         pexsi_temperature, pexsi_tol_charge, pexsi_np_sym_fact, &
  !!         mat_k, energy, mat_ek)
  !!else if (trim(kernel_method)=='LAPACK') then
  !!    norbu = smat_h%nfvctr
  !!    norbd = 0
  !!    norb = norbu + norbd
  !!    coeff = f_malloc((/smat_h%nfvctr,norb/),id='coeff')
  !!    eval_all = f_malloc(smat_h%nfvctr,id='eval_all')
  !!    eval_occup = f_malloc(norb,id='eval_occup')
  !!    call get_coeffs_diagonalization(iproc, nproc, mpi_comm_world, &
  !!         smat_h%nfvctr, norbu, norbd, norb, scalapack_blocksize, &
  !!         smat_s, smat_h, mat_s, mat_h, coeff, &
  !!         eval_all, eval_occup, info)
  !!    occup = f_malloc0(norb,id='occup')
  !!    if (smat_h%nspin/=1) call f_err_throw('The kernel calculation with LAPACK is currently not possible for nspin/=1')
  !!    ii = nint(0.5_mp*foe_data_get_real(foe_obj,"charge",1))
  !!    occup(1:ii) = 2.0_mp
  !!    call eval_to_occ(iproc, nproc, norbu, norbd, norb, 1, (/1.0_mp/), &
  !!         eval_occup, occup, .false., .true., pexsi_temperature, SMEARING_DIST_ERF, efermi, eTS, &
  !!         norbu, norbd)
  !!    call distribute_on_tasks(norb, iproc, nproc, norbp, isorb)
  !!    ! Density kernel
  !!    call calculate_kernel_and_energy(iproc, nproc, mpi_comm_world, &
  !!         smat_k, smat_h, mat_k, mat_h, energy,&
  !!         coeff, norbp, isorb, norbu, norb, occup, .true.)
  !!    ! Energy density kernel
  !!    do iorb=1,norb
  !!        occup(iorb) = occup(iorb)*eval_occup(iorb)
  !!    end do
  !!    call calculate_kernel_and_energy(iproc, nproc, mpi_comm_world, &
  !!         smat_k, smat_h, mat_ek, mat_h, energy_fake,&
  !!         coeff, norbp, isorb, norbu, norb, occup, .true.)
  !!    call f_free(coeff)
  !!    call f_free(eval_all)
  !!    call f_free(eval_occup)
  !!    call f_free(occup)
  !!else
  !!    call f_err_throw("wrong value for 'kernel_method'; possible values are 'FOE', 'PEXSI' or 'LAPACK'")
  !!end if

  !!call mpibarrier()
  !!call f_timing_checkpoint(ctr_name='CALC',mpi_comm=mpiworld(),nproc=mpisize(), &
  !!     gather_routine=gather_timings)

  !!!tr = trace_A(iproc, nproc, mpi_comm_world, smat_k, mat_ek, 1)
  !!tr = trace_AB(iproc, nproc, mpi_comm_world, smat_s, smat_k, mat_s, mat_ek, 1)
  !!if (iproc==0) then
  !!    call yaml_map('Energy from FOE',energy)
  !!    call yaml_map('Trace of energy density kernel', tr)
  !!    call yaml_map('Difference',abs(energy-tr))
  !!end if

  !!!! Write the result in YAML format to the standard output (required for non-regression tests).
  !!!if (iproc==0) call write_matrix_compressed('Result of FOE', smat_k, mat_k)

  !!! Calculate trace(KS)
  !!!tr_KS = trace_sparse(iproc, nproc, smat_s, smat_k, mat_s%matrix_compr, mat_k%matrix_compr, 1)
  !!tr_KS = trace_AB(iproc, nproc, mpi_comm_world, smat_s, smat_k, mat_s, mat_k, 1)

  !!! Write the result
  !!if (iproc==0) call yaml_map('trace(KS)',tr_KS)

  !!! Extract the compressed kernel matrix from the data type.
  !!! The first routine allocates an array with the correct size, the second one extracts the result.
  !!kernel = sparsematrix_malloc_ptr(smat_k, iaction=SPARSE_FULL, id='kernel')
  !!call matrices_get_values(smat_k, mat_k, kernel)

  !!! Do the same also for the overlap matrix
  !!overlap = sparsematrix_malloc_ptr(smat_s, iaction=SPARSE_FULL, id='overlap')
  !!call matrices_get_values(smat_s, mat_s, overlap)

  !!! Transform the overlap matrix to the sparsity pattern of the kernel
  !!overlap_large = sparsematrix_malloc_ptr(smat_k, iaction=SPARSE_FULL, id='overlap_large')
  !!call transform_sparse_matrix(iproc, smat_s, smat_k, SPARSE_FULL, 'small_to_large', &
  !!     smat_in=overlap, lmat_out=overlap_large)

  !!! Again calculate trace(KS), this time directly with the array holding the data.
  !!! Since both matrices are symmetric and have now the same sparsity pattern, this is a simple ddot.
  !!tr_KS_check = dot(smat_k%nvctr, kernel(1), 1, overlap_large(1), 1)

  !!! Write the result
  !!if (iproc==0) call yaml_map('trace(KS) check',tr_KS_check)

  !!! Write the difference to the previous result
  !!if (iproc==0) call yaml_map('difference',tr_KS-tr_KS_check)

  !!if (do_cubic_check) then
  !!    norbu = smat_h%nfvctr
  !!    norbd = 0
  !!    norb = norbu + norbd
  !!    coeff = f_malloc((/smat_h%nfvctr,norb/),id='coeff')
  !!    eval_all = f_malloc(smat_h%nfvctr,id='eval_all')
  !!    eval_occup = f_malloc(norb,id='eval_occup')
  !!    call get_coeffs_diagonalization(iproc, nproc, mpi_comm_world, &
  !!         smat_h%nfvctr, norbu, norbd, norb, scalapack_blocksize, &
  !!         smat_s, smat_h, mat_s, mat_h, coeff, &
  !!         eval_all, eval_occup, info)
  !!    occup = f_malloc0(norb,id='occup')
  !!    if (smat_h%nspin/=1) call f_err_throw('The kernel calculation with LAPACK is currently not possible for nspin/=1')
  !!    ii = nint(0.5_mp*foe_data_get_real(foe_obj,"charge",1))
  !!    occup(1:ii) = 2.0_mp
  !!    call eval_to_occ(iproc, nproc, norbu, norbd, norb, 1, (/1.0_mp/), &
  !!         eval_occup, occup, .false., .true., pexsi_temperature, SMEARING_DIST_ERF, efermi, eTS, &
  !!         norbu, norbd)
  !!    call distribute_on_tasks(norb, iproc, nproc, norbp, isorb)
  !!    ! Density kernel... Use mat_ek to store this reference kernel.
  !!    call calculate_kernel_and_energy(iproc, nproc, mpi_comm_world, &
  !!         smat_k, smat_h, mat_ek, mat_h, energy,&
  !!         coeff, norbp, isorb, norbu, norb, occup, .true.)
  !!    call calculate_error(iproc, smat_k, mat_k, mat_ek, nthreshold, threshold, .false., &
  !!         'Check the deviation from the reference calculation using LAPACK')
  !!    call f_free(coeff)
  !!    call f_free(eval_all)
  !!    call f_free(eval_occup)
  !!    call f_free(occup)
  !!end if

  !!! Deallocate the object holding the FOE parameters
  !!call foe_data_deallocate(foe_obj)
  !!call foe_data_deallocate(ice_obj)

  ! Deallocate all the sparse matrix descriptors types
  !!call deallocate_sparse_matrix(smat_s)
  !!call deallocate_sparse_matrix(smat_h)
  !!call deallocate_sparse_matrix(smat_k)
  !!call deallocate_sparse_matrix_metadata(smmd)

  ! Deallocate all the matrix data types
  !!call deallocate_matrices(mat_s)
  !!call deallocate_matrices(mat_h)
  !!call deallocate_matrices(mat_k)
  !!call deallocate_matrices(mat_ek)
  !!call deallocate_matrices(mat_ovrlpminusonehalf(1))

  !!! Deallocate all the remaining arrays
  !!call f_free(charge)
  !!call f_free_ptr(kernel)
  !!call f_free_ptr(overlap)
  !!call f_free_ptr(overlap_large)

  call mpibarrier()
  call f_timing_checkpoint(ctr_name='LAST',mpi_comm=mpiworld(),nproc=mpisize(), &
       gather_routine=gather_timings)

  call build_dict_info(dict_timing_info)
  call f_timing_stop(mpi_comm=mpi_comm_world,nproc=nproc,&
       gather_routine=gather_timings,dict_info=dict_timing_info)
  call dict_free(dict_timing_info)

  if (iproc==0) then
      call yaml_release_document()
  end if

  ! Finalize MPI
  call mpifinalize()

  ! Finalize flib
  call f_lib_finalize()


  !!contains

  !!  !> construct the dictionary needed for the timing information
  !!  !! SM: This routine should go to a module
  !!  subroutine build_dict_info(dict_info)
  !!    use wrapper_MPI
  !!    use dynamic_memory
  !!    use dictionaries
  !!    implicit none

  !!    type(dictionary), pointer :: dict_info
  !!    !local variables
  !!    integer :: ierr,namelen,nthreads
  !!    character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
  !!    character(len=MPI_MAX_PROCESSOR_NAME), dimension(:), allocatable :: nodename
  !!    type(dictionary), pointer :: dict_tmp
  !!    !$ integer :: omp_get_max_threads

  !!    call dict_init(dict_info)
! !! bastian: comment out 4 followinf lines for debug purposes (7.12.2014)
  !!    !if (DoLastRunThings) then
  !!       call f_malloc_dump_status(dict_summary=dict_tmp)
  !!       call set(dict_info//'Routines timing and number of calls',dict_tmp)
  !!    !end if
  !!    nthreads = 0
  !!    !$  nthreads=omp_get_max_threads()
  !!    call set(dict_info//'CPU parallelism'//'MPI tasks',nproc)
  !!    if (nthreads /= 0) call set(dict_info//'CPU parallelism'//'OMP threads',&
  !!         nthreads)

  !!    nodename=f_malloc0_str(MPI_MAX_PROCESSOR_NAME,0.to.nproc-1,id='nodename')
  !!    if (nproc>1) then
  !!       call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)
  !!       !gather the result between all the process
  !!       call MPI_GATHER(nodename_local,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
  !!            nodename(0),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,0,&
  !!            mpi_comm_world,ierr)
  !!       if (iproc==0) call set(dict_info//'Hostnames',&
  !!               list_new(.item. nodename))
  !!    end if
  !!    call f_free_str(MPI_MAX_PROCESSOR_NAME,nodename)

  !!  end subroutine build_dict_info

end program driver_eigenvalues


subroutine commandline_options(parser)
  use yaml_parse
  use dictionaries, only: dict_new,operator(.is.)
  implicit none
  type(yaml_cl_parse),intent(inout) :: parser

  call yaml_cl_parse_option(parser,'metadata_file','sparsematrix_metadata.dat',&
       'input file with the matrix metadata',&
       help_dict=dict_new('Usage' .is. &
       'Input file with the matrix metadata',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'overlap_file','overlap_sparse.txt',&
       'input file with the overlap matrix',&
       help_dict=dict_new('Usage' .is. &
       'Input file with the overlap matrix',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'hamiltonian_file','hamiltonian_sparse.txt',&
       'input file with the Hamiltonian matrix',&
       help_dict=dict_new('Usage' .is. &
       'Input file with the Hamiltonian matrix',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'kernel_file','density_kernel_sparse.txt',&
       'input file with the density kernel matrix',&
       help_dict=dict_new('Usage' .is. &
       'Input file with the density kernel matrix',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'kernel_matmul_file','density_kernel_sparse_matmul.txt',&
       'input file with the density kernel matrix multiplication matrix',&
       help_dict=dict_new('Usage' .is. &
       'Input file with the density kernel matrix multiplication matrix',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'kernel_method','FOE',&
       'Indicate which kernel method should be used (FOE or PEXSI)',&
       help_dict=dict_new('Usage' .is. &
       'Indicate which kernel method should be used (FOE or PEXSI)',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'sparsity_format','bigdft',&
       'indicate the sparsity format',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the sparsity format of the input matrices',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'matrix_format','serial_text',&
       'indicate the matrix format',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the format of the input matrices',&
       'Allowed values' .is. &
       'String'))

  call yaml_cl_parse_option(parser,'scalapack_blocksize','-1',&
       'Indicate the blocksize for ScaLAPACK (negative for LAPACK)',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the blocksize for ScaLAPACK (negative for LAPACK)',&
       'Allowed values' .is. &
       'Integer'))

  call yaml_cl_parse_option(parser,'check_spectrum','.false.',&
       'Indicate whether the spectral properties of the Hamiltonian shall be calculated',&
       help_dict=dict_new('Usage' .is. &
       'Indicate whether the spectral properties of the Hamiltonian shall be calculated by a diagonalization',&
       'Allowed values' .is. &
       'Logical'))

  call yaml_cl_parse_option(parser,'fscale_lowerbound','5.e-3',&
       'Indicate the minimal value for the Fermi function decay length',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the minimal value for the Fermi function decay length that should be used to calculate the density kernel',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'fscale_upperbound','5.e-2',&
       'Indicate the maximal value for the Fermi function decay length',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the maximal value for the Fermi function decay length that should be used to calculate the density kernel',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'ntemp','4',&
       'Indicate the maximal number of FOE iterations with various temperatures',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the maximal number of FOE iterations with various temperatures, to determine dynamically the width of the gap',&
       'Allowed values' .is. &
       'Integer'))

  call yaml_cl_parse_option(parser,'ef','0.0',&
       'Indicate the initial guess for the Fermi energy',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the initial guess for the Fermi energy, will then be adjusted automatically',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'npl_max','5000',&
       'Indicate the maximal polynomial degree',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the maximal polynomial degree',&
       'Allowed values' .is. &
       'Integer'))

  call yaml_cl_parse_option(parser,'pexsi_npoles','40',&
       'Indicate the number of poles to be used by PEXSI',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the number of poles to be used by PEXSI',&
       'Allowed values' .is. &
       'Integer'))

  call yaml_cl_parse_option(parser,'pexsi_mumin','-1.0',&
       'Initial guess for the lower bound of the chemical potential used by PEXSI',&
       help_dict=dict_new('Usage' .is. &
       'Initial guesss for the lower bound of the chemical potential (in hartree?) used by PEXSI,&
       & will be adjusted automatically later',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'pexsi_mumax','1.0',&
       'Initial guess for the lower bound of the chemical potential used by PEXSI',&
       help_dict=dict_new('Usage' .is. &
       'Initial guesss for the upper bound of the chemical potential (in hartree?) used by PEXSI,&
       & will be adjusted automatically later',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'pexsi_mu','0.0',&
       'initial guess for the chemical potential used by PEXSI',&
       help_dict=dict_new('Usage' .is. &
       'Initial guesss for the chemical potential (in hartree?) used by PEXSI, will be adjusted automatically later',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'pexsi_DeltaE','10.0',&
       'upper bound for the spectral radius of S^-1H used by PEXSI',&
       help_dict=dict_new('Usage' .is. &
       'upper bound for the spectral radius of S^-1H (in hartree?) used by PEXSI',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'pexsi_temperature','1.e-3',&
       'Indicate the temperature used by PEXSI',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the temperature (in atomic units?) used by PEXSI',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'pexsi_tol_charge','1.e-3',&
       'Indicate the charge tolerance used PEXSI',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the tolerance on the number of electrons used by PEXSI',&
       'Allowed values' .is. &
       'Double'))

  call yaml_cl_parse_option(parser,'pexsi_np_sym_fact','16',&
       'Indicate the number of tasks used for the symbolic factorization within PEXSI',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the number of tasks used for the symbolic factorization within PEXSI',&
       'Allowed values' .is. &
       'Integer'))

   call yaml_cl_parse_option(parser,'do_cubic_check','.true.',&
       'perform a check using cubic scaling dense (Sca)LAPACK',&
       help_dict=dict_new('Usage' .is. &
       'Indicate whether a cubic scaling check using dense (Sca)LAPACK should be performed',&
       'Allowed values' .is. &
       'Logical'))

   call yaml_cl_parse_option(parser,'iev_min','1',&
       'Lowest eigenvalue to be calculated',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the lowest eigenvalue to be calculated',&
       'Allowed values' .is. &
       'Integer'))

   call yaml_cl_parse_option(parser,'iev_max','1000000',&
       'Higest eigenvalue to be calculated',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the highest eigenvalue to be calculated',&
       'Allowed values' .is. &
       'Integer'))

   call yaml_cl_parse_option(parser,'fscale','5.e-3',&
       'Decay length of the error function used to find the eigenvalue',&
       help_dict=dict_new('Usage' .is. &
       'Indicate the decay length of the error function used to find the eigenvalue',&
       'Allowed values' .is. &
       'Double'))


end subroutine commandline_options
