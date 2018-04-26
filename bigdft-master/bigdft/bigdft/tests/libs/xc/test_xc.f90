!> @file
!! Test of the libxc library
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Program to test the libxc library used by BigDFT
program test_xc

  use module_base
  use module_xc
  use yaml_output

  implicit none

  integer, parameter :: n_funcs = 54
  integer, dimension(n_funcs), parameter :: funcs = (/ &
       & 1, -20, 2, -1009, 3, 4, -1002, 5, -1004, 6, -6006, &
       & 7, -1012, 8, -1, 9, -1003, -1005, -1007, -1008, -1010, &
       & -1011, -1013, -1014, &
       & 11, -101130, 12, -101, 13, -160012, 14, -102130, 15, -117130, &
       & 16, -161, 17, -162, 23, -118, 26, -163, 27, -164, &
       & -102, -103, -104, -105, -106, -107, -108, -109, -110, &
       & -406 /)
!!$  integer, parameter :: n_funcs = 1
!!$  integer, dimension(n_funcs), parameter :: funcs = (/ -101130 /)
  integer :: ifunc, ixc_prev, ierr, iproc, nproc
  real(dp) :: tt0
  real(dp) :: exc_(2, n_funcs), dt_(n_funcs),tt(2)
  real(dp) :: exc(2, n_funcs), dt(n_funcs)

  !Initialize f_lib
  call f_lib_initialize()
  
  !Initialize MPI environment
  !call MPI_INIT(ierr)
  !call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  !call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  call mpiinit()
  iproc=mpirank()
  nproc=mpisize()

  call f_malloc_set_status(iproc=iproc)
  exc_ = 0.d0
  dt_  = 0.d0
  do ifunc = 1, 1!,n_funcs, 1
     if (modulo(ifunc, nproc) == iproc) then
        !print *,'exc',ifunc,funcs(ifunc)
        call test(funcs(ifunc), tt, tt0)
        exc_(1,ifunc)=tt(1)
        exc_(2,ifunc)=tt(2)
        dt_(ifunc)=tt0
     end if
!!$     if (funcs(ifunc) < 0) then
!!$        call test(funcs(ifunc), exc, dt, option = XC_LIBXC)
!!$        write(*,"(1x,A,I7,3x,A,F17.8,1x,A,1x,A,F17.8,3x,A,F10.5,1x,A)") &
!!$             & "ixc = ", funcs(ifunc), "nosp = ", exc(1), "|", "scol = ", &
!!$             & exc(2), "time = ", dt, "s"
!!$     end if
  end do
  if (nproc == 1) then
     exc = exc_
     dt = dt_
  else
     call MPI_ALLREDUCE(exc_, exc, 2 * n_funcs, &
          & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(dt_, dt, n_funcs, &
          & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  end if

  if (iproc == 0) then
     call yaml_sequence_open('Test XC')
     ixc_prev = 1
     do ifunc = 1, 1!n_funcs, 1
        !if (ixc_prev * funcs(ifunc) > 0 .or. (ixc_prev < 0 .and. funcs(ifunc) > 0)) &
        !     & write(*,"(1x,A,A,A)") repeat("-", 41), "+", repeat("-", 44)
        !write(*,"(1x,A,I7,3x,A,F17.8,1x,A,1x,A,F17.8,3x,A,F10.5,1x,A)") &
        !     & "ixc = ", funcs(ifunc), "nosp = ", exc(1, ifunc), "|", "scol = ", &
        !     & exc(2, ifunc), "time = ", dt(ifunc), "s"
        call yaml_sequence(advance='no')
        call yaml_mapping_open(flow=.true.)
        call yaml_map('ixc',funcs(ifunc))
        call yaml_map('nosp',exc(1,ifunc),fmt='(f17.8)')
        call yaml_map('scol',exc(2,ifunc),fmt='(f17.8)')
        call yaml_map('time',dt(ifunc),fmt='(f10.5)')
        call yaml_mapping_close()
        ixc_prev = funcs(ifunc)
     end do
     call yaml_sequence_close()
     !write(*,"(1x,A,A,A)") repeat("-", 41), "+", repeat("-", 44)
  end if

  call mpifinalize()
  call f_lib_finalize()

contains

  subroutine test(ixc, excs, dt, option)
    
    implicit none

    integer, intent(in) :: ixc
    real(dp), intent(out) :: excs(2), dt
    integer, intent(in), optional :: option

    integer :: i, itmp, n, type,unt,j
    type(xc_info) :: xc
    integer, parameter :: n_rho = 100000, n_runs = 1
    real(dp), dimension(:,:), allocatable :: rho, vxc
    real(dp), dimension(:,:), allocatable :: rhogr, vxcgr,fxc,gradv,gradr
    real(dp), dimension(:), allocatable :: exc
    integer :: start, end, countPerSecond

    exc=f_malloc(n_rho,id='exc')
    rho=f_malloc([n_rho,2],id='rho')
    vxc=f_malloc0([n_rho,2],id='vxc')
    rhogr=f_malloc([n_rho,3],id='rhogr')
    vxcgr=f_malloc0([n_rho,3],id='vxcgr')
    gradv=f_malloc([n_rho,3],id='gradv')
    gradr=f_malloc([n_rho,2],id='gradr')
    fxc=f_malloc([n_rho,3],id='fxc')

    if (present(option)) then
       type = option
    else
       if (ixc < 0) then
          type = XC_MIXED
       else
          type = XC_ABINIT
       end if
    end if
    call system_clock(count_rate = countPerSecond)
    call system_clock(start)
    do n = 1, n_runs, 1
       do i = 1, 2, 1
          itmp=2
          call xc_init(xc, ixc, type, i)
!!$          if (i == 1 .and. n == 1) call xc_dump()

          call gauss(xc, rho, n_rho, i, type)
          if (xc_isgga(xc)) call gaussgr(rhogr, rho, n_rho)
          call xc_getvxc(xc, n_rho, exc, i, rho(1,1), vxc(1,1), rhogr, vxcgr,dvxci=fxc)
          !call gaussgr(rhogr, rho, n_rho)
          excs(i) = sum(exc)
          do j=1,itmp
             call fgr(gradv(1,j),vxc(1,j),n_rho)
             call fgr(gradr(1,j),rho(1,j),n_rho)
          end do
          call xc_end(xc)
          call f_open_file(unit=unt,file='test'+ixc+'-'+i+'.dat')
          !if (i==1) then
          !   do j=1,n_rho
          !      write(unt,'(i8,15(1pe25.17))')j,rho(j,1),exc(j),vxc(j,1),&
          !         fxc(j,1:2),gradv(j,1),(fxc(j,1)+fxc(j,2))*gradr(j,1)
          !         !fxc(j,1:2),gradv(j,1),2.d0*fxc(j,1)*gradr(j,1),2.d0*fxc(j,2)*gradr(j,1)
          !   end do
          !else if (i==2) then
             do j=1,n_rho
                write(unt,'(i8,15(1pe25.17))')j,rho(j,1:2),exc(j),vxc(j,1:2),&
                   fxc(j,1:3),gradv(j,1:2),gradr(j,1:2),fxc(j,1)*gradr(j,1)+fxc(j,2)*gradr(j,2),&
                   fxc(j,2)*gradr(j,1)+fxc(j,3)*gradr(j,2)
             end do
          !end if
          call f_close(unt)
       end do
    end do
    unt=11

    call system_clock(end)

    dt = real(end - start) / real(countPerSecond) / real(n_runs)

    call f_free(exc)
    call f_free(rho)
    call f_free(vxc)
    call f_free(rhogr)
    call f_free(vxcgr)
    call f_free(fxc)
    call f_free(gradv)
    call f_free(gradr)

  end subroutine test

  subroutine gauss(xc, rho, n_rho, nspin, type)
    use module_base
    use module_xc

    implicit none

    type(xc_info), intent(in) :: xc
    integer, intent(in) :: n_rho, nspin, type
    real(dp), intent(out) :: rho(n_rho * 2)

    integer :: j, delta
    real(dp) :: sigma

    call xc_init_rho(xc, n_rho * 2, rho, 1)
    delta = 0
    !if (nspin == 2 ) delta = int(real(n_rho) * 0.05)
    sigma = 1.d0 / (real(n_rho, dp) * 0.25d0)
    do j = 5, n_rho - 5
       if (type == XC_LIBXC .and. nspin == 2) then
          rho(2 * j - 1) = exp(-((sigma * (n_rho / 2 - j + delta)) ** 2)) / nspin / n_rho
          rho(2 * j    ) = exp(-((sigma * (n_rho / 2 - j - delta)) ** 2)) / nspin / n_rho
       else
          rho(j        ) = exp(-((sigma * (n_rho / 2 - j + delta)) ** 2)) / nspin / n_rho
          rho(j + n_rho) = exp(-((sigma * (n_rho / 2 - j - delta)) ** 2)) / nspin / n_rho
       end if
    end do
    call xc_clean_rho(xc, n_rho * 2, rho, 1)
    rho = rho / sum(rho(1:nspin * n_rho))
  end subroutine gauss

  subroutine gaussgr(rhogr, rho, n_rho)
    implicit none

    integer, intent(in) :: n_rho
    real(dp), intent(in) :: rho(n_rho, 2)
    real(dp), intent(out) :: rhogr(n_rho, 3)

    integer :: j

    rhogr(1, :) = 0.d0
    do j = 2, n_rho - 1
       rhogr(j, 1) = (rho(j + 1, 1) - rho(j - 1, 1)) / 2.d0 * n_rho
       rhogr(j, 2) = (rho(j + 1, 2) - rho(j - 1, 2)) / 2.d0 * n_rho
       rhogr(j, 3) = rhogr(j, 1) + rhogr(j, 2)
       rhogr(j, 1) = rhogr(j, 1) * rhogr(j, 1)
       rhogr(j, 2) = rhogr(j, 2) * rhogr(j, 2)
       rhogr(j, 3) = rhogr(j, 3) * rhogr(j, 3)
    end do
    rhogr(n_rho, :) = 0.d0
  end subroutine gaussgr

  subroutine fgr(rhogr, rho, n_rho)
    implicit none

    integer, intent(in) :: n_rho
    real(dp), intent(in) :: rho(n_rho)
    real(dp), intent(out) :: rhogr(n_rho)

    integer :: j

    rhogr(1) = 0.d0
    do j = 2, n_rho - 1
       rhogr(j) = (rho(j + 1) - rho(j - 1)) / 2.d0 * n_rho
    end do
    rhogr(n_rho) = 0.d0
  end subroutine fgr

end program test_xc
