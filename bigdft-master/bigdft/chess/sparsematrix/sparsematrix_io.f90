
!> @file
!!   Sparse matrix I/O routines
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


module sparsematrix_io
  use sparsematrix_base
  use wrapper_mpi
  use dynamic_memory
  use f_utils
  use dictionaries, only: f_err_throw
  use yaml_output
  implicit none

  private

  public :: write_ccs_matrix
  public :: read_sparse_matrix
  public :: read_sparse_matrix_metadata
  public :: write_sparse_matrix
  public :: write_sparse_matrix_metadata
  public :: write_linear_coefficients
  public :: read_linear_coefficients
  public :: write_dense_matrix
  public :: read_dense_matrix
  public :: write_linear_eigenvalues

  contains

    subroutine write_ccs_matrix(filename, nfvctr, nvctr, row_ind, col_ptr, mat_compr)
      use dynamic_memory
      use f_utils
      implicit none
      !Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(in) :: nfvctr !number of rows/columns
      integer,intent(in) :: nvctr !number of non-zero elements
      integer,dimension(nvctr),intent(in) :: row_ind
      integer,dimension(nfvctr),intent(in) :: col_ptr
      real(kind=mp),dimension(nvctr),intent(in) :: mat_compr
      ! Local variables
      integer :: i, iunit

      call f_routine(id='write_ccs_matrix')

      iunit = 99
      call f_open_file(iunit, file=trim(filename), binary=.false.)

      write(iunit,'(4(i0,1x))') nfvctr, nfvctr, nvctr, 0
      write(iunit,'(10000(i0,1x))') (col_ptr(i),i=1,nfvctr),nvctr+1
      write(iunit,'(10000(i0,1x))') (row_ind(i),i=1,nvctr)
      do i=1,nvctr
          write(iunit,'(es24.15)') mat_compr(i)
      end do

      call f_close(iunit)

      call f_release_routine()

    end subroutine write_ccs_matrix


    subroutine read_sparse_matrix(mode, filename, iproc, nproc, comm, nspin, nfvctr, nseg, nvctr, keyv, keyg, mat_compr)
      use time_profiling
      use dynamic_memory
      use f_utils
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      character(len=*),intent(in) :: mode
      character(len=*),intent(in) :: filename
      integer,intent(out) :: nspin, nfvctr, nseg, nvctr
      integer,dimension(:),pointer,intent(out) :: keyv
      integer,dimension(:,:,:),pointer,intent(out) :: keyg
      real(kind=mp),dimension(:),pointer,intent(out) :: mat_compr

      ! Local variables
      integer :: iunit, iseg, icol, irow, jorb, ind, ispin, index_dot
      !real(kind=mp) :: dummy_double
      !character(len=20) :: dummy_char
      character(len=128) :: filename_extension
      logical :: file_present

      call f_routine(id='read_sparse_matrix')
      call f_timing(TCAT_SMAT_READ,'ON')

      if (iproc==0) call yaml_comment('Reading from file '//trim(filename),hfill='~')
      inquire(file=trim(filename),exist=file_present)
      if (.not.file_present) then
          call f_err_throw("File '"//trim(filename)//"' is not present", &
               err_name='SPARSEMATRIX_IO_ERROR')
      end if

      index_dot = index(filename,'.',back=.true.)
      filename_extension = filename(index_dot:)

      if (trim(mode)=='parallel_mpi-native') then
          if (trim(filename_extension)/='.mpi') then
              call f_err_throw("Wrong file extension; '.mpi' is required, but found "//trim(filename_extension)&
                   &//" ("//trim(filename)//")", &
                   err_name='SPARSEMATRIX_IO_ERROR')
          end if
          call read_sparse_matrix_parallel(filename, iproc, nproc, comm, nspin, nfvctr, nseg, nvctr, keyv, keyg, mat_compr)
      else if (trim(mode)=='serial_text') then
          if (trim(filename_extension)/='.txt') then
              call f_err_throw("Wrong file extension; '.txt' is required, but found "//trim(filename_extension)&
                   &//" ("//trim(filename)//")", &
                   err_name='SPARSEMATRIX_IO_ERROR')
          end if

          if (iproc==0) then

              ! Read in the matrix on task 0
              iunit = 99
              call f_open_file(iunit, file=trim(filename), binary=.false.)

              read(iunit,*) nspin, nfvctr, nseg, nvctr
              keyv = f_malloc_ptr(nseg,id='keyv')
              keyg = f_malloc_ptr((/2,2,nseg/),id='keyg')

              do iseg=1,nseg
                  read(iunit,*) keyv(iseg), keyg(1,1,iseg), keyg(2,1,iseg), keyg(1,2,iseg), keyg(2,2,iseg)
              end do

              mat_compr = f_malloc_ptr(nvctr*nspin,id='mat_compr')
              ind = 0
              do ispin=1,nspin
                  do iseg=1,nseg
                      icol = keyg(1,2,iseg)
                      do jorb=keyg(1,1,iseg),keyg(2,1,iseg)
                          irow = jorb
                          ind = ind + 1
                          read(iunit,*) mat_compr(ind)
                      end do
                  end do
              end do
          end if

          ! Communicate to the other tasks
          call mpibcast(nspin, count=1, root=0, comm=comm)
          call mpibcast(nfvctr, count=1, root=0, comm=comm)
          call mpibcast(nseg, count=1, root=0, comm=comm)
          call mpibcast(nvctr, count=1, root=0, comm=comm)
          if (iproc/=0) then
              keyv = f_malloc_ptr(nseg,id='keyv')
              keyg = f_malloc_ptr((/2,2,nseg/),id='keyg')
          end if
          call mpibcast(keyv, root=0, comm=comm)
          call mpibcast(keyg, root=0, comm=comm)
          if (iproc/=0) then
              mat_compr = f_malloc_ptr(nvctr*nspin,id='mat_compr')
          end if
          call mpibcast(mat_compr, root=0, comm=comm)

      else
          call f_err_throw("wrong value for 'mode'")
      end if

      call f_close(iunit)

      call f_timing(TCAT_SMAT_READ,'OF')
      call f_release_routine()

    end subroutine read_sparse_matrix


    subroutine read_sparse_matrix_parallel(filename, iproc, nproc, comm, nspin, nfvctr, nseg, nvctr, keyv, keyg, mat_compr)
      use sparsematrix_init, only: distribute_on_tasks
      implicit none
      
      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(in) :: iproc, nproc, comm
      integer,intent(out) :: nspin, nfvctr, nseg, nvctr
      integer,dimension(:),pointer,intent(out) :: keyv
      integer,dimension(:,:,:),pointer,intent(out) :: keyg
      real(kind=mp),dimension(:),pointer,intent(out) :: mat_compr

      ! Local variables
      integer :: i, ii, np, is, size_of_integer, size_of_double, ierr, thefile, ispin, ioffset
      !character(len=1024) :: filename_base, filename_extension, filename_matmul
      integer(kind=mpi_offset_kind) :: disp
      integer,dimension(4) :: workarr_header
      integer,dimension(:,:),allocatable :: workarr_keys
      integer(kind=f_long) :: is_long, nseg_long, four_long, five_long, size_of_integer_long, size_of_double_long

      call f_routine(id='read_sparse_matrix_parallel')

      call mpi_file_open(comm, trim(filename), & 
           mpi_mode_rdonly, & 
           mpi_info_null, thefile, ierr) 
      size_of_integer = mpitypesize(1)
      size_of_double = mpitypesize(1.0_mp)
      size_of_integer_long = int(size_of_integer,kind=f_long)
      size_of_double_long = int(size_of_double,kind=f_long)
      four_long = int(4,kind=f_long)
      five_long = int(5,kind=f_long)

      ! Read the header
      disp = int(0,kind=mpi_offset_kind)
      call mpi_file_set_view(thefile, disp, mpi_integer, mpi_integer, 'native', mpi_info_null, ierr) 
      !if (iproc==0) then
          call mpi_file_read(thefile, workarr_header, 4, mpi_integer, mpi_status_ignore, ierr)
          nspin = workarr_header(1)
          nfvctr = workarr_header(2)
          nseg =  workarr_header(3)
          nvctr =  workarr_header(4)
      !end if

      ! Read the matrix keys
      keyv = f_malloc0_ptr(nseg,id='keyv')
      keyg = f_malloc0_ptr((/2,2,nseg/),id='keyg')
      call distribute_on_tasks(nseg, iproc, nproc, np, is)
      workarr_keys = f_malloc((/5,np/),id='workarr_keys')
      disp = int((4+5*is)*size_of_integer,kind=mpi_offset_kind)
      call mpi_file_set_view(thefile, disp, mpi_integer, mpi_integer, 'native', mpi_info_null, ierr) 
      call mpi_file_read(thefile, workarr_keys, 5*np, mpi_integer, mpi_status_ignore, ierr)
      do i=1,np
          ii = is + i
          keyv(ii) = workarr_keys(1,i)
          keyg(1,1,ii) = workarr_keys(2,i)
          keyg(2,1,ii) = workarr_keys(3,i)
          keyg(1,2,ii) = workarr_keys(4,i)
          keyg(2,2,ii) = workarr_keys(5,i)
      end do
      call f_free(workarr_keys)
      call fmpi_allreduce(keyv, FMPI_SUM, comm=comm)
      call fmpi_allreduce(keyg, FMPI_SUM, comm=comm)

      ! Read the matrices
      mat_compr = f_malloc0_ptr(nvctr*nspin,id='mat_compr')
      call distribute_on_tasks(nvctr, iproc, nproc, np, is)
      do ispin=1,nspin
          ioffset = (ispin-1)*nvctr
          is = is + ioffset
          is_long = int(is,kind=f_long)
          nseg_long = int(nseg,kind=f_long)
          disp = int((four_long+five_long*nseg_long)*size_of_integer_long+is_long*size_of_double_long,kind=mpi_offset_kind)
          call mpi_file_set_view(thefile, disp, mpi_double_precision, mpi_double_precision, 'native', mpi_info_null, ierr) 
          if (np>1) then
              call mpi_file_read(thefile, mat_compr(is+1), np, mpi_double_precision, mpi_status_ignore, ierr)
          end if
      end do
      call fmpi_allreduce(mat_compr, FMPI_SUM, comm=comm)

      call mpi_file_close(thefile, ierr)      

      call f_release_routine()

    end subroutine read_sparse_matrix_parallel


    subroutine read_sparse_matrix_metadata(filename, nfvctr, nat, ntypes, units, geocode, cell_dim, shift, &
               nzatom, nelpsp, atomnames, iatype, rxyz, on_which_atom)
      use dynamic_memory
      use f_utils
      implicit none
      
      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(out) :: nfvctr
      character(len=20),intent(out) :: units
      character(len=1),intent(out) :: geocode
      real(kind=mp),dimension(3),intent(out) :: cell_dim, shift
      integer,intent(out) :: nat, ntypes
      integer,dimension(:),pointer,intent(inout) :: nzatom, nelpsp, iatype
      character(len=20),dimension(:),pointer,intent(inout) :: atomnames
      real(kind=mp),dimension(:,:),pointer,intent(inout) :: rxyz
      integer,dimension(:),pointer,intent(inout) :: on_which_atom

      ! Local variables
      integer :: iunit, iat, itype, i
      !real(kind=mp) :: dummy_double
      !character(len=20) :: dummy_char
      !logical :: read_rxyz, read_on_which_atom

      call f_routine(id='read_sparse_matrix_metadata')

      iunit = 99
      call f_open_file(iunit, file=trim(filename), binary=.false.)

      read(iunit,*) nfvctr, nat, ntypes
      read(iunit,*) units
      read(iunit,*) geocode, cell_dim
      read(iunit,*) shift
      nzatom = f_malloc_ptr(ntypes,id='nzatom')
      nelpsp = f_malloc_ptr(ntypes,id='nelpsp')
      atomnames = f_malloc0_str_ptr(len(atomnames),ntypes,id='atomnames')

      do itype=1,ntypes
          read(iunit,*) nzatom(itype), nelpsp(itype), atomnames(itype)
      end do
      rxyz = f_malloc_ptr((/3,nat/),id='rxyz')
      iatype = f_malloc_ptr(nat,id='iatype')
      do iat=1,nat
          read(iunit,*) iatype(iat), rxyz(1,iat), rxyz(2,iat), rxyz(3,iat)
      end do
      on_which_atom = f_malloc_ptr(nfvctr,id='on_which_atom')
      do i=1,nfvctr
          read(iunit,*) on_which_atom(i)
      end do

      call f_close(iunit)

      call f_release_routine()

    end subroutine read_sparse_matrix_metadata


    !> Write a sparse matrix to disk.
    !! ATTENTION: This routine must be called by all MPI tasks due to the fact that the matrix 
    !! is distributed among the matrix taksgroups
    subroutine write_sparse_matrix(mode, iproc, nproc, comm, smat, mat, filename)
      use dynamic_memory
      use f_utils
      use sparsematrix, only: gather_matrix_from_taskgroups
      implicit none
      
      ! Calling arguments
      character(len=*),intent(in) :: mode
      integer,intent(in) :: iproc, nproc, comm
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(in) :: mat
      character(len=*),intent(in) :: filename

      ! Local variables
      integer :: iunit, iseg, icol, irow, jorb, ind, ispin
      !integer :: size_of_integer, size_of_double
      real(kind=mp),dimension(:),allocatable :: matrix_compr
      character(len=1024) :: filename_matmul

      call f_routine(id='write_sparse_matrix')
      
      if (len(filename)>1024) then
          call f_err_throw('filename is too long')
      end if

      matrix_compr = sparsematrix_malloc(smat,iaction=SPARSE_FULL,id='matrix_compr')
      call gather_matrix_from_taskgroups(iproc, nproc, comm, &
           smat, mat%matrix_compr, matrix_compr)

      !!call write_sparse_matrix_parallel(iproc, nproc, comm, smat, matrix_compr, trim(filename)//'.mpi')
      if (trim(mode)=='parallel_mpi-native') then
          call write_sparse_matrix_parallel(iproc, nproc, comm, &
               smat%nspin, smat%nfvctr, smat%nseg, smat%nvctr, smat%keyv, smat%keyg, &
               matrix_compr, trim(filename)//'.mpi')
          call f_free(matrix_compr)
          if (smat%smatmul_initialized) then
              !index_dot = index(filename,'.',back=.true.)
              !filename_base = filename(1:index_dot-1)
              !filename_extension = filename(index_dot:)
              !filename_matmul = trim(filename_base)//'_matmul'//trim(filename_extension)
              filename_matmul = trim(filename)//'_matmul'
              ! Fill the array with zero, as the entries have no meaning
              matrix_compr = f_malloc0(sum(smat%smmm%nvctr_par),id='matrix_compr')
              call write_sparse_matrix_parallel(iproc, nproc, comm, &
                   smat%nspin, smat%nfvctr, smat%smmm%nseg, sum(smat%smmm%nvctr_par), &
                   smat%smmm%keyv, smat%smmm%keyg, &
                   matrix_compr, trim(filename_matmul)//'.mpi')
              call f_free(matrix_compr)
          end if
      else if (trim(mode)=='serial_text') then

          if (iproc==0) then

              iunit = 99
              call f_open_file(iunit, file=trim(filename)//'.txt', binary=.false.)

              !!write(iunit,'(i10,2i6,3x,a,3es24.16,a)') nat, ntypes, smat%nspin, smat%geocode, smat%cell_dim, &
              !!    '   # number of atoms, number of atom types, nspin, geocode, cell_dim'
              !!do itype=1,ntypes
              !!    write(iunit,'(2i8,3x,a,a)') nzatom(itype), nelpsp(itype), trim(atomnames(itype)), &
              !!        '   # nz, nelpsp, name'
              !!end do
              !!do iat=1,nat
              !!    write(iunit,'(i5, 3es24.16,a,i0)') iatype(iat), rxyz(1:3,iat), '   # atom no. ',iat
              !!end do
              write(iunit,'(4i12,a)') smat%nspin, smat%nfvctr, smat%nseg, smat%nvctr, '   # nspin, nfvctr, nseg, nvctr'
              do iseg=1,smat%nseg
                  write(iunit,'(5i12,a)') smat%keyv(iseg), smat%keyg(1,1,iseg), smat%keyg(2,1,iseg), &
                      smat%keyg(1,2,iseg), smat%keyg(2,2,iseg), '   # keyv, keyg(1,1), keyg(2,1), keyg(1,2), keyg(2,2)'
              end do
              ind = 0
              do ispin=1,smat%nspin
                  do iseg=1,smat%nseg
                      icol = smat%keyg(1,2,iseg)
                      !iat = smat%on_which_atom(icol)
                      do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                          irow = jorb
                          !jat = smat%on_which_atom(irow)
                          ind = ind + 1
                          !write(iunit,'(es24.16,2i12,a)') matrix_compr(ind), jat, iat, '   # matrix, jat, iat'
                          write(iunit,'(es24.16,a,i0,a)') matrix_compr(ind),'   # matrix_compr(',ind,')'
                      end do
                  end do
              end do

              call f_close(iunit)
              call f_free(matrix_compr)


              if (smat%smatmul_initialized) then
                  iunit = 99
                  !index_dot = index(filename,'.',back=.true.)
                  !filename_base = filename(1:index_dot-1)
                  !filename_extension = filename(index_dot:)
                  !filename_matmul = trim(filename_base)//'_matmul'//trim(filename_extension)
                  filename_matmul = trim(filename)//'_matmul'
                  !call f_open_file(iunit, file=trim(filename)//'_matmul', binary=.false.)
                  call f_open_file(iunit, file=trim(filename_matmul)//'.txt', binary=.false.)

                  write(iunit,'(4i12,a)') smat%nspin, smat%nfvctr, smat%smmm%nseg, sum(smat%smmm%nvctr_par), &
                      '   # nspin, nfvctr, nseg, nvctr'
                  do iseg=1,smat%smmm%nseg
                      write(iunit,'(5i12,a)') smat%smmm%keyv(iseg), smat%smmm%keyg(1,1,iseg), smat%smmm%keyg(2,1,iseg), &
                          smat%smmm%keyg(1,2,iseg), smat%smmm%keyg(2,2,iseg), &
                          '   # keyv, keyg(1,1), keyg(2,1), keyg(1,2), keyg(2,2)'
                  end do
                  ind = 0
                  do ispin=1,smat%nspin
                      do iseg=1,smat%smmm%nseg
                          icol = smat%smmm%keyg(1,2,iseg)
                          !iat = smat%on_which_atom(icol)
                          do jorb=smat%smmm%keyg(1,1,iseg),smat%smmm%keyg(2,1,iseg)
                              irow = jorb
                              !jat = smat%on_which_atom(irow)
                              ind = ind + 1
                              !write(iunit,'(es24.16,2i12,a)') matrix_compr(ind), jat, iat, '   # matrix, jat, iat'
                              write(iunit,'(es24.16,a,i0,a)') 0.0_mp,'   # matrix_compr(',ind,')'
                          end do
                      end do
                  end do

                  call f_close(iunit)
              end if

          end if

      else
          call f_err_throw("wrong value for 'mode'")
      end if


      call f_release_routine()

    end subroutine write_sparse_matrix



    subroutine write_sparse_matrix_parallel(iproc, nproc, comm, nspin, nfvctr, nseg, nvctr, &
               keyv, keyg, matrix_compr, filename)
      use sparsematrix_init, only: distribute_on_tasks
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      integer,intent(in) :: nspin, nfvctr, nseg, nvctr
      integer,dimension(nseg),intent(in) :: keyv
      integer,dimension(2,2,nseg),intent(in) :: keyg
      real(kind=mp),dimension(nvctr*nspin),intent(in) :: matrix_compr
      character(len=*),intent(in) :: filename

      ! Local variables
      integer :: i, ii, np, is, size_of_integer, size_of_double, ierr, thefile, ispin, ioffset
      !character(len=1024) :: filename_base, filename_extension, filename_matmul
      integer(kind=mpi_offset_kind) :: disp
      integer,dimension(4) :: workarr_header
      integer,dimension(:,:),allocatable :: workarr_keys
      integer(kind=f_long) :: is_long, nseg_long, four_long, five_long, size_of_integer_long, size_of_double_long

      call f_routine(id='write_sparse_matrix_parallel')

      call mpi_file_open(comm, trim(filename), & 
           mpi_mode_wronly + mpi_mode_create, & 
           mpi_info_null, thefile, ierr) 
      size_of_integer = mpitypesize(1)
      size_of_double = mpitypesize(1.0_mp)
      size_of_integer_long = int(size_of_integer,kind=f_long)
      size_of_double_long = int(size_of_double,kind=f_long)
      four_long = int(4,kind=f_long)
      five_long = int(5,kind=f_long)

      ! Write the header
      disp = int(0,kind=mpi_offset_kind)
      call mpi_file_set_view(thefile, disp, mpi_integer, mpi_integer, 'native', mpi_info_null, ierr) 
      if (iproc==0) then
          workarr_header(1) = nspin
          workarr_header(2) = nfvctr
          workarr_header(3) = nseg
          workarr_header(4) = nvctr
          call mpi_file_write(thefile, workarr_header, 4, mpi_integer, mpi_status_ignore, ierr)
      end if

      ! Write the matrix keys
      call distribute_on_tasks(nseg, iproc, nproc, np, is)
      workarr_keys = f_malloc((/5,np/),id='workarr_keys')
      do i=1,np
          ii = is + i
          workarr_keys(1,i) = keyv(ii)
          workarr_keys(2,i) = keyg(1,1,ii)
          workarr_keys(3,i) = keyg(2,1,ii)
          workarr_keys(4,i) = keyg(1,2,ii)
          workarr_keys(5,i) = keyg(2,2,ii)
      end do
      !disp = int((4+5*is)*size_of_integer,kind=mpi_offset_kind)
      is_long = int(is,kind=f_long)
      disp = int((four_long+five_long*is_long)*size_of_integer_long,kind=mpi_offset_kind)
      call mpi_file_set_view(thefile, disp, mpi_integer, mpi_integer, 'native', mpi_info_null, ierr) 
      call mpi_file_write(thefile, workarr_keys, 5*np, mpi_integer, mpi_status_ignore, ierr)
      call f_free(workarr_keys)

      ! Write the matrices
      call distribute_on_tasks(nvctr, iproc, nproc, np, is)
      do ispin=1,nspin
          ioffset = (ispin-1)*nvctr
          is = is + ioffset
          is_long = int(is,kind=f_long)
          nseg_long = int(nseg,kind=f_long)
          disp = int((four_long+five_long*nseg_long)*size_of_integer_long+is_long*size_of_double_long,kind=mpi_offset_kind)
          call mpi_file_set_view(thefile, disp, mpi_double_precision, mpi_double_precision, 'native', mpi_info_null, ierr)
          if (np>1) then
              call mpi_file_write(thefile, matrix_compr(is+1), np, mpi_double_precision, mpi_status_ignore, ierr)
          end if
      end do

      call mpi_file_close(thefile, ierr)      

      call f_release_routine()

    end subroutine write_sparse_matrix_parallel




    subroutine write_sparse_matrix_metadata(iproc, nfvctr, nat, ntypes, units, geocode, cell_dim, shift, iatype, &
               rxyz, nzatom, nelpsp, atomnames, on_which_atom, filename)
      use dynamic_memory
      use f_utils
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nfvctr, nat, ntypes
      character(len=*),intent(in) :: units
      character(len=1),intent(in) :: geocode
      real(kind=mp),dimension(3),intent(in) :: cell_dim, shift
      integer,dimension(nat),intent(in) :: iatype
      real(kind=mp),dimension(3,nat),intent(in) :: rxyz
      integer,dimension(ntypes),intent(in) :: nzatom, nelpsp
      character(len=*),dimension(ntypes),intent(in) :: atomnames
      integer,dimension(nfvctr),intent(in) :: on_which_atom
      character(len=*),intent(in) :: filename

      ! Local variables
      integer :: iunit, iat, itype, i
      !real(kind=mp),dimension(:),allocatable :: matrix_compr

      call f_routine(id='write_sparse_matrix_metadata')

      if (iproc==0) then

          iunit = 99
          call f_open_file(iunit, file=trim(filename), binary=.false.)

          write(iunit,'(i10,i8,i6,a)') nfvctr, nat, ntypes, &
              '   # matrix dimension, number of atoms, number of atom types'
          write(iunit,'(a,a)') trim(units), '  # units'
          write(iunit,'(a,3es24.16,a)') geocode, cell_dim, &
              '   # geocode, cell_dim'
          write(iunit,'(3es24.16,a)') -shift, &
              '   # atomic shift'
          do itype=1,ntypes
              write(iunit,'(2i8,3x,a,a)') nzatom(itype), nelpsp(itype), trim(atomnames(itype)), &
                  '   # nz, nelpsp, name'
          end do
          do iat=1,nat
              write(iunit,'(i5,3es24.16,a,i0)') iatype(iat), rxyz(1:3,iat), '   # atom no. ',iat
          end do
          do i=1,nfvctr
              write(iunit,'(i8,a,i0,a)') on_which_atom(i), '   # on_which_atom(',i,')'
          end do

          call f_close(iunit)

      end if

      call f_release_routine()

    end subroutine write_sparse_matrix_metadata


    subroutine write_linear_coefficients(mode, iproc, nproc, comm, iroot, filename, verbosity, &
               nfvctr, ntmb, nspin, coeff, eval)
      use yaml_output
      implicit none
      ! Calling arguments
      character(len=*),intent(in) :: mode
      character(len=*),intent(in) :: filename
      !type(atoms_data),intent(in) :: at
      integer,intent(in) :: iproc, nproc, comm, iroot, verbosity, nfvctr, ntmb, nspin
      real(mp), dimension(nfvctr,ntmb), intent(in) :: coeff
      real(mp), dimension(ntmb), intent(in) :: eval
      ! Local variables
      integer :: iunit, i, j
      logical :: scaled

      call f_routine(id='write_linear_coefficients')


      if (trim(mode)=='parallel_mpi-native') then
          call write_linear_coefficients_parallel(iproc, nproc, comm, trim(filename), &
               nspin, ntmb, nfvctr, eval, coeff)
      else if (trim(mode)=='serial_text') then
          if (iproc==iroot) then

              iunit = 99
              call f_open_file(iunit, file=trim(filename), binary=.false.)
    
              ! Write the Header
              !!!write(iunit,'(i10,2i6,a)') at%astruct%nat, at%astruct%ntypes, nspin, &
              !!write(iunit,'(i10,2i6,a)') nat, ntypes, nspin, &
              !!    '   # number of atoms, number of atom types, nspin'
              !!!do itype=1,at%astruct%ntypes
              !!do itype=1,ntypes
              !!    !write(iunit,'(2i8,3x,a,a)') at%nzatom(itype), at%nelpsp(itype), trim(at%astruct%atomnames(itype)), &
              !!    write(iunit,'(2i8,3x,a,a)') nzatom(itype), nelpsp(itype), trim(atomnames(itype)), &
              !!        '   # nz, nelpsp, name'
              !!end do
              !!!do iat=1,at%astruct%nat
              !!do iat=1,nat
              !!    !write(iunit,'(i5, 3es24.16,a,i0)') at%astruct%iatype(iat), rxyz(1:3,iat), '   # atom no. ',iat
              !!    write(iunit,'(i5, 3es24.16,a,i0)') iatype(iat), rxyz(1:3,iat), '   # atom no. ',iat
              !!end do
              write(iunit,'(3i12,a)') nspin, ntmb, nfvctr, '   # nspin, ntmb, nfvctr'
              do i=1,ntmb
                  write(iunit,'(es24.16,a,i0)') eval(i), '   # eval no. ', i
              enddo
    
              ! Now write the coefficients
              do i=1,ntmb
                 ! First element always positive, for consistency when using for transfer integrals;
                 ! unless 1st element below some threshold, in which case first significant element.
                 scaled = .false.
                 do j=1,nfvctr
                    if (abs(coeff(j,i))>1.0d-3) then
                       if (coeff(j,i)<0.0_mp) call dscal(ntmb,-1.0_mp,coeff(1,i),1)
                       scaled = .true.
                       exit
                    end if
                 end do
                 if (.not.scaled) then
                     call yaml_warning('Consistency between the written coefficients not guaranteed')
                 end if
    
                 do j = 1,nfvctr
                     write(iunit,'(es24.16,2i9,a)') coeff(j,i), j, i, '   # coeff, j, i'
                 end do
              end do  
              if (verbosity >= 2 .and. iproc==0) call yaml_map('Wavefunction coefficients written',.true.)

              call f_close(iunit)

          end if
      else
          call f_err_throw("wrong value for 'mode'")
      end if

      call f_release_routine()
    
    end subroutine write_linear_coefficients


    subroutine write_linear_coefficients_parallel(iproc, nproc, comm, filename, &
               nspin, ntmb, nfvctr, eval, coeff)
      use sparsematrix_init, only: distribute_on_tasks
      use wrapper_linalg, only: vcopy
      use wrapper_mpi, only: mpibcast
      use yaml_output
      implicit none

      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(in) :: iproc, nproc, comm, nspin, ntmb, nfvctr
      real(mp), dimension(ntmb), intent(in) :: eval
      real(mp), dimension(nfvctr,ntmb), intent(in) :: coeff

      ! Local variables
      integer :: thefile, size_of_integer, size_of_double, i, ii, j, is, np, ierr
      integer(kind=f_long) :: size_of_integer_long, size_of_double_long, disp, is_long
      integer(kind=f_long) :: one_long, three_long, five_long, nfvctr_long, ntmb_long
      logical :: scaled
      integer,dimension(3) :: workarr_header
      real(mp),dimension(:),allocatable :: workarr_eval
      real(mp),dimension(:,:),allocatable :: workarr_coeff

      call f_routine(id='write_linear_coefficients_parallel')

      ! Make sure that all processes have the same filename
      call mpibcast(filename, root=0, comm=comm)

      call mpi_file_open(comm, trim(filename), & 
           mpi_mode_wronly + mpi_mode_create, & 
           mpi_info_null, thefile, ierr) 
      size_of_integer = mpitypesize(1)
      size_of_double = mpitypesize(1.0_mp)
      size_of_integer_long = int(size_of_integer,kind=f_long)
      size_of_double_long = int(size_of_double,kind=f_long)
      one_long = int(1,kind=f_long)
      three_long = int(3,kind=f_long)
      five_long = int(5,kind=f_long)
      nfvctr_long = int(nfvctr,kind=f_long)
      ntmb_long = int(ntmb,kind=f_long)

      ! Write the header
      disp = int(0,kind=mpi_offset_kind)
      call mpi_file_set_view(thefile, disp, mpi_integer, mpi_integer, 'native', mpi_info_null, ierr) 
      if (iproc==0) then
          workarr_header(1) = nspin
          workarr_header(2) = ntmb
          workarr_header(3) = nfvctr
          call mpi_file_write(thefile, workarr_header, 3, mpi_integer, mpi_status_ignore, ierr)
      end if

      ! Write the eigenvalues
      call distribute_on_tasks(ntmb, iproc, nproc, np, is)
      workarr_eval = f_malloc(np,id='workarr_eval')
      do i=1,np
          ii = is + i
          workarr_eval(i) = eval(ii)
      end do
      is_long = int(is,kind=f_long)
      disp = int(three_long*size_of_integer_long+is_long*size_of_double_long,kind=mpi_offset_kind)
      call mpi_file_set_view(thefile, disp, mpi_double_precision, mpi_double_precision, 'native', mpi_info_null, ierr) 
      call mpi_file_write(thefile, workarr_eval, np, mpi_double_precision, mpi_status_ignore, ierr)
      call f_free(workarr_eval)

      workarr_coeff = f_malloc((/nfvctr,np/),id='workarr_coeff')
      do i=1,np
          ii = is + i
          ! First element always positive, for consistency when using for transfer integrals;
          ! unless 1st element below some threshold, in which case first significant element.
          scaled = .false.
          do j=1,nfvctr
             if (abs(coeff(j,ii))>1.e-3_mp) then
                if (coeff(j,ii)<0.0_mp) call dscal(ntmb,-1.0_mp,coeff(1,i),1)
                scaled = .true.
                exit
             end if
          end do
          if (.not.scaled) then
              call yaml_warning('Consistency between the written coefficients not guaranteed')
          end if
          !call f_memcpy(n=nfvctr, src=coeff(1:nfvctr,ii:ii), dest=workarr_coeff(1:nfvctr,i:i))
          call vcopy(nfvctr, coeff(1,ii), 1, workarr_coeff(1,i), 1)
      end do
      disp = int(three_long*size_of_integer_long+(ntmb_long+is_long*nfvctr_long)*size_of_double_long,kind=mpi_offset_kind)
      call mpi_file_set_view(thefile, disp, mpi_double_precision, mpi_double_precision, 'native', mpi_info_null, ierr) 
      call mpi_file_write(thefile, workarr_coeff, nfvctr*np, mpi_double_precision, mpi_status_ignore, ierr)
      call f_free(workarr_coeff)


      call f_release_routine()
    
    end subroutine write_linear_coefficients_parallel


    subroutine write_linear_eigenvalues(iproc, iroot, filename, nfvctr, ntmb, nspin, eval)
      use yaml_output
      implicit none
      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(in) :: iproc, iroot, nfvctr, ntmb, nspin
      real(mp), dimension(ntmb), intent(in) :: eval
      ! Local variables
      integer :: iunit, i

      call f_routine(id='write_linear_eigenvalues')

      if (iproc==iroot) then

          iunit = 99
          call f_open_file(iunit, file=trim(filename), binary=.false.)
    
          call yaml_map('nspin',nspin,unit=iunit)
          call yaml_map('ntmb',ntmb,unit=iunit)
          call yaml_map('nfvctr',nfvctr,unit=iunit)
          call yaml_sequence_open('Eigenvalues',unit=iunit)
          do i=1,ntmb
              call yaml_sequence(yaml_toa(eval(i)),advance='no',unit=iunit)
              !call yaml_sequence(advance='no',unit=iunit)
              !call yaml_map('e',eval(i),unit=iunit,advance='no')
              call yaml_comment(adjustl(trim(yaml_toa(i))),unit=iunit)
          end do
          call yaml_sequence_close()
          !!write(iunit,'(3i12,a)') nspin, ntmb, nfvctr, '   # nspin, ntmb, nfvctr'
          !!do i=1,ntmb
          !!    write(iunit,'(es24.16,a,i0)') eval(i), '   # eval no. ', i
          !!end do

          call f_close(iunit)

      end if
    
      call f_release_routine()
    
    end subroutine write_linear_eigenvalues


    subroutine read_linear_coefficients(mode, iproc, nproc, comm, filename, nspin, nfvctr, ntmb, coeff, eval)
      use yaml_output
      implicit none
      ! Calling arguments
      character(len=*),intent(in) :: mode
      character(len=*),intent(in) :: filename
      integer,intent(in) :: iproc, nproc, comm
      integer,intent(out) :: nspin, nfvctr, ntmb
      real(kind=8),dimension(:,:),pointer,intent(inout) :: coeff
      !integer,intent(out),optional :: nat, ntypes
      !integer,dimension(:),pointer,intent(inout),optional :: nzatom, nelpsp, iatype
      !character(len=20),dimension(:),pointer,intent(inout),optional :: atomnames
      !real(kind=8),dimension(:,:),pointer,intent(inout),optional :: rxyz
      real(kind=8),dimension(:),pointer,intent(inout),optional :: eval
      ! Local variables
      real(kind=8) :: dummy_double
      !character(len=20) :: dummy_char
      integer :: iunit, i, j, index_dot
      logical :: scaled, read_rxyz, read_eval, file_present
      character(len=1024) :: filename_extension

      call f_routine(id='read_linear_coefficients')

      !!if (present(nat) .and. present(ntypes) .and. present(nzatom) .and.  &
      !!    present(nelpsp) .and. present(atomnames) .and. present(iatype) .and. present(rxyz)) then
      !!    read_rxyz = .true.
      !!else if (present(nat) .or. present(ntypes) .or. present(nzatom) .or.  &
      !!    present(nelpsp) .or. present(atomnames) .or. present(iatype) .or. present(rxyz)) then
      !!    call f_err_throw("not all optional arguments were given", &
      !!         err_name='BIGDFT_RUNTIME_ERROR')
      !!else
          read_rxyz = .false.
      !!end if

      if (present(eval)) then
          read_eval = .true.
      else
          read_eval = .false.
      end if

      if (iproc==0) call yaml_comment('Reading from file '//trim(filename),hfill='~')
      inquire(file=trim(filename),exist=file_present)
      !write(*,*) 'file_present',file_present
      if (.not.file_present) then
          call f_err_throw("File '"//trim(filename)//"' is not present", &
               err_name='SPARSEMATRIX_IO_ERROR')
      end if

      index_dot = index(filename,'.',back=.true.)
      filename_extension = filename(index_dot:)

      if (trim(mode)=='parallel_mpi-native') then
          if (trim(filename_extension)/='.mpi') then
              call f_err_throw("Wrong file extension; '.mpi' is required, but found "//trim(filename_extension)&
                   &//" ("//trim(filename)//")", &
                   err_name='SPARSEMATRIX_IO_ERROR')
          end if
          if (read_eval) then
              call read_linear_coefficients_parallel(iproc, nproc, comm, trim(filename), &
                   nspin, ntmb, nfvctr, coeff, eval)
          else
              call read_linear_coefficients_parallel(iproc, nproc, comm, trim(filename), &
                   nspin, ntmb, nfvctr, coeff)
          end if
      else if (trim(mode)=='serial_text') then
          if (trim(filename_extension)/='.txt') then
              call f_err_throw("Wrong file extension; '.txt' is required, but found "//trim(filename_extension)&
                   &//" ("//trim(filename)//")", &
                   err_name='SPARSEMATRIX_IO_ERROR')
          end if

          iunit = 99
          call f_open_file(iunit, file=trim(filename), binary=.false.)
    
          ! Read the Header
          !!if (read_rxyz) then
          !!    read(iunit,*) nat, ntypes, nspin
          !!    nzatom = f_malloc_ptr(ntypes,id='nzatom')
          !!    nelpsp = f_malloc_ptr(ntypes,id='nelpsp')
          !!    atomnames = f_malloc0_str_ptr(len(atomnames),ntypes,id='atomnames')

          !!    do itype=1,ntypes
          !!        read(iunit,*) nzatom(itype), nelpsp(itype), atomnames(itype)
          !!    end do
          !!    rxyz = f_malloc_ptr((/3,nat/),id='rxyz')
          !!    iatype = f_malloc_ptr(nat,id='iatype')
          !!    do iat=1,nat
          !!        read(iunit,*) iatype(iat), rxyz(1,iat), rxyz(2,iat), rxyz(3,iat)
          !!    end do
          !!else
          !!    read(iunit,*) nat_, ntypes_, nspin
          !!    do itype=1,ntypes_
          !!        read(iunit,*) dummy_int, dummy_int, dummy_char
          !!    end do
          !!    do iat=1,nat_
          !!        read(iunit,*) dummy_int, dummy_double, dummy_double, dummy_double
          !!    end do
          !!end if

          read(iunit,*) nspin, ntmb, nfvctr

          if (read_eval) then
              eval = f_malloc_ptr(ntmb,id='eval')
              do i=1,ntmb
                  read(iunit,*) eval(i)
              end do
          else
              do i=1,ntmb
                  read(iunit,*) dummy_double
              end do
          end if
    
          ! Now read the coefficients
          coeff = f_malloc_ptr((/nfvctr,ntmb/),id='coeff')

          do i=1,ntmb

             do j = 1,nfvctr
                 read(iunit,*) coeff(j,i)
             end do

             ! First element always positive, for consistency when using for transfer integrals;
             ! unless 1st element below some threshold, in which case first significant element.
             scaled = .false.
             do j=1,nfvctr
                if (abs(coeff(j,i))>1.0d-3) then
                   if (coeff(j,i)<0.0_mp) call dscal(ntmb,-1.0_mp,coeff(1,i),1)
                   scaled = .true.
                   exit
                end if
             end do
             if (.not.scaled) then
                 call yaml_warning('Consistency between the written coefficients not guaranteed')
             end if
    
          end do  

          call f_close(iunit)

      else
          call f_err_throw("wrong value for 'mode'")
      end if


      call f_release_routine()
    
    end subroutine read_linear_coefficients


    subroutine read_linear_coefficients_parallel(iproc, nproc, comm, filename, &
               nspin, ntmb, nfvctr, coeff, eval)
      use sparsematrix_init, only: distribute_on_tasks
      use wrapper_linalg, only: vcopy
      use yaml_output
      implicit none

      ! Calling arguments
      character(len=*),intent(in) :: filename
      integer,intent(in) :: iproc, nproc, comm
      integer,intent(out) :: nspin, ntmb, nfvctr
      real(mp),dimension(:,:),pointer,intent(out) :: coeff
      real(mp),dimension(:),pointer,intent(out),optional :: eval

      ! Local variables
      integer :: thefile, size_of_integer, size_of_double, i, ii, is, np, ierr
      integer(kind=f_long) :: size_of_integer_long, size_of_double_long, disp, is_long
      integer(kind=f_long) :: one_long, three_long, five_long, nfvctr_long, ntmb_long
      !logical :: scaled
      integer,dimension(3) :: workarr_header
      real(mp),dimension(:),allocatable :: workarr_eval
      real(mp),dimension(:,:),allocatable :: workarr_coeff

      call f_routine(id='read_linear_coefficients_parallel')

      call mpi_file_open(comm, trim(filename), & 
           mpi_mode_rdonly, & 
           mpi_info_null, thefile, ierr) 
      size_of_integer = mpitypesize(1)
      size_of_double = mpitypesize(1.0_mp)
      size_of_integer_long = int(size_of_integer,kind=f_long)
      size_of_double_long = int(size_of_double,kind=f_long)
      one_long = int(1,kind=f_long)
      three_long = int(3,kind=f_long)
      five_long = int(5,kind=f_long)

      ! Read the header
      disp = int(0,kind=mpi_offset_kind)
      call mpi_file_set_view(thefile, disp, mpi_integer, mpi_integer, 'native', mpi_info_null, ierr) 
      call mpi_file_read(thefile, workarr_header, 3, mpi_integer, mpi_status_ignore, ierr)
      nspin = workarr_header(1)
      ntmb = workarr_header(2)
      nfvctr = workarr_header(3)

      nfvctr_long = int(nfvctr,kind=f_long)
      ntmb_long = int(ntmb,kind=f_long)

      ! Read the eigenvalues
      call distribute_on_tasks(ntmb, iproc, nproc, np, is)
      is_long = int(is,kind=f_long)
      if (present(eval)) then
          workarr_eval = f_malloc(np,id='workarr_eval')
          disp = int(three_long*size_of_integer_long+is_long*size_of_double_long,kind=mpi_offset_kind)
          call mpi_file_set_view(thefile, disp, mpi_double_precision, mpi_double_precision, 'native', mpi_info_null, ierr) 
          call mpi_file_read(thefile, workarr_eval, np, mpi_double_precision, mpi_status_ignore, ierr)
          eval = f_malloc0_ptr(ntmb,id='eval')
          do i=1,np
              ii = is + i
              eval(ii) = workarr_eval(i) 
          end do
          call f_free(workarr_eval)
          call fmpi_allreduce(eval, FMPI_SUM, comm=comm)
      end if

      workarr_coeff = f_malloc((/nfvctr,np/),id='workarr_coeff')
      disp = int(three_long*size_of_integer_long+(ntmb_long+is_long*nfvctr_long)*size_of_double_long,kind=mpi_offset_kind)
      call mpi_file_set_view(thefile, disp, mpi_double_precision, mpi_double_precision, 'native', mpi_info_null, ierr) 
      call mpi_file_read(thefile, workarr_coeff, nfvctr*np, mpi_double_precision, mpi_status_ignore, ierr)
      coeff = f_malloc0_ptr((/nfvctr,ntmb/),id='coeff')
      do i=1,np
          ii = is + i
          call vcopy(nfvctr, workarr_coeff(1,i), 1, coeff(1,ii), 1)
      end do
      call f_free(workarr_coeff)
      call fmpi_allreduce(coeff, FMPI_SUM, comm=comm)

      call f_release_routine()
    
    end subroutine read_linear_coefficients_parallel


    !> Write a sparse matrix to disk, but in dense format
    !! ATTENTION: This routine must be called by all MPI tasks due to the fact that the matrix 
    !! in distributed among the matrix taksgroups
    subroutine write_dense_matrix(iproc, nproc, comm, smat, mat, uncompress, filename, binary)
      use sparsematrix_base
      use sparsematrix, only: uncompress_matrix2
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      type(sparse_matrix),intent(in) :: smat
      type(matrices),intent(inout) :: mat
      logical,intent(in) :: uncompress
      character(len=*),intent(in) :: filename
      logical, intent(in) :: binary

      ! Local variables
      integer :: iunit, jorb, ispin, iorb
      !real(kind=8),dimension(:),allocatable :: matrix_compr

      call f_routine(id='write_dense_matrix')


      if (uncompress) then
          mat%matrix = sparsematrix_malloc_ptr(smat, iaction=DENSE_FULL, id='mat%matrix')
          call uncompress_matrix2(iproc, nproc, comm, &
               smat, mat%matrix_compr, mat%matrix)
      end if

      if (iproc==0) then

          iunit = 99
          call f_open_file(iunit, file=trim(filename), binary=binary)

          if (.not. binary) then
              write(iunit,'(2i10,a)') smat%nspin, smat%nfvctr, ' # nspin, nfvctr'
          else
              write(iunit) smat%nspin, smat%nfvctr, ' # nspin, nfvctr'
          end if
    
          do ispin=1,smat%nspin
             do iorb=1,smat%nfvctr
                do jorb=1,smat%nfvctr
                   if (.not. binary) then
                      write(iunit,'(2(i6,1x),es19.12,a)') iorb,jorb,mat%matrix(iorb,jorb,ispin), &
                          '   # i, j, mat(i,j)'
                   else
                      write(iunit) iorb,jorb,mat%matrix(iorb,jorb,ispin)
                   end if
                end do
             end do
          end do

          call f_close(iunit)

          if (uncompress) then
              call f_free_ptr(mat%matrix)
          end if
      end if

      call f_release_routine()

    end subroutine write_dense_matrix




    !> Read a sparse matrix from disk, but in dense format
    !! ATTENTION: This routine must be called by all MPI tasks
    subroutine read_dense_matrix(binary, filename, iproc, nproc, comm, nfvctr, nspin, matrix)
      implicit none
      
      ! Calling arguments
      logical, intent(in) :: binary
      character(len=*),intent(in) :: filename
      integer,intent(in) :: iproc, nproc, comm
      integer,intent(out) :: nfvctr, nspin
      real(kind=mp),dimension(:,:,:),pointer,intent(out) :: matrix

      ! Local variables
      integer :: iunit, ispin, jorb, iorb, idum, jdum

      call f_routine(id='read_dense_matrix')


      if (binary) then
          call f_err_throw('dense matrices cannot be read in binary mode')
      end if

      if (.not.binary) then

          ! Read in the matrix on task 0
          if (iproc==0) then

              iunit = 99
              call f_open_file(iunit, file=trim(filename), binary=binary)

              read(iunit,*) nspin, nfvctr

              matrix = f_malloc_ptr((/nfvctr,nfvctr,nspin/),id='matrix')
    
              do ispin=1,nspin
                  do iorb=1,nfvctr
                      do jorb=1,nfvctr
                          read(iunit,*) idum,jdum,matrix(iorb,jorb,ispin)
                      end do
                  end do
              end do

              call f_close(iunit)

          end if

          ! Communicate the matrix
          call mpibcast(nfvctr, root=0, comm=comm)
          call mpibcast(nspin, root=0, comm=comm)
          if (iproc/=0) then
              matrix = f_malloc_ptr((/nfvctr,nfvctr,nspin/),id='matrix')
          end if
          call mpibcast(matrix, root=0, comm=comm)

      end if

      call f_release_routine()

    end subroutine read_dense_matrix

end module sparsematrix_io
