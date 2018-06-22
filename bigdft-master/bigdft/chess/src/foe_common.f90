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


module module_func
  use sparsematrix_base
  use numerics
  use dynamic_memory
  use yaml_output
  use dictionaries, only: f_err_throw
  implicit none

  private

  ! Shared variables within the modules
  integer :: ifunc
  real(kind=mp) :: power, ef, fscale, beta, mua, mub

  ! Public routines
  public :: func_set
  public :: func

  ! Public parameters
  integer,parameter,public :: FUNCTION_POLYNOMIAL            = 101
  integer,parameter,public :: FUNCTION_ERRORFUNCTION         = 102
  integer,parameter,public :: FUNCTION_XTIMESERRORFUNCTION   = 103
  integer,parameter,public :: FUNCTION_EXPONENTIAL           = 104
  integer,parameter,public :: FUNCTION_ERRORFUNCTION_ENTROPY = 105
  integer,parameter,public :: FUNCTION_FERMIFUNCTION         = 106
  integer,parameter,public :: FUNCTION_FERMIFUNCTION_ENTROPY = 107

  contains

    subroutine func_set(ifuncx, powerx, efx, fscalex, betax, muax, mubx)
      implicit none
      integer,intent(in) :: ifuncx
      real(kind=mp),intent(in),optional :: powerx, efx, fscalex, betax, muax, mubx

      call f_routine(id='func_set')

      select case (ifuncx)
      case(FUNCTION_POLYNOMIAL)
          ifunc = FUNCTION_POLYNOMIAL
          if (.not.present(powerx)) call f_err_throw("'powerx' not present")
          power = powerx
      case(FUNCTION_ERRORFUNCTION)
          ifunc = FUNCTION_ERRORFUNCTION
          if (.not.present(efx)) call f_err_throw("'efx' not present")
          if (.not.present(fscalex)) call f_err_throw("'fscalex' not present")
          ef = efx
          fscale = fscalex
      case(FUNCTION_XTIMESERRORFUNCTION)
          ifunc = FUNCTION_XTIMESERRORFUNCTION
          if (.not.present(efx)) call f_err_throw("'efx' not present")
          if (.not.present(fscalex)) call f_err_throw("'fscalex' not present")
          ef = efx
          fscale = fscalex
      case(FUNCTION_EXPONENTIAL)
          ifunc = FUNCTION_EXPONENTIAL
          if (.not.present(betax)) call f_err_throw("'betax' not present")
          if (.not.present(muax)) call f_err_throw("'muax' not present")
          if (.not.present(mubx)) call f_err_throw("'mubx' not present")
          beta = betax
          mua = muax
          mub = mubx
      case(FUNCTION_ERRORFUNCTION_ENTROPY)
          ifunc = FUNCTION_ERRORFUNCTION_ENTROPY
          if (.not.present(efx)) call f_err_throw("'efx' not present")
          if (.not.present(fscalex)) call f_err_throw("'fscalex' not present")
          ef = efx
          fscale = fscalex
      case(FUNCTION_FERMIFUNCTION)
          ifunc = FUNCTION_FERMIFUNCTION
          if (.not.present(efx)) call f_err_throw("'efx' not present")
          if (.not.present(fscalex)) call f_err_throw("'fscalex' not present")
          ef = efx
          fscale = fscalex
      case(FUNCTION_FERMIFUNCTION_ENTROPY)
          ifunc = FUNCTION_FERMIFUNCTION_ENTROPY
      case default
          call f_err_throw("wrong value of 'ifuncx'")
      end select

      call f_release_routine()

    end subroutine func_set

    function func(x)
      use numerics, only: pi
      implicit none
      ! Calling arguments
      real(kind=mp),intent(in) :: x
      real(kind=mp) :: func
      ! Local parameters
      real(kind=mp) :: sqrt_pi = sqrt(pi)
      select case (ifunc)
      case(FUNCTION_POLYNOMIAL)
          func = x**power
      case(FUNCTION_ERRORFUNCTION)
          func = 0.5_mp*erfcc((x-ef)*(1._mp/fscale))
      case(FUNCTION_XTIMESERRORFUNCTION)
          func = x*0.5_mp*erfcc((x-ef)*(1._mp/fscale))
      case(FUNCTION_EXPONENTIAL)
          !func = safe_exp(beta*(x-mu))
          func = safe_exp(beta*(x-mua)) - safe_exp(-beta*(x-mub))
      case(FUNCTION_ERRORFUNCTION_ENTROPY)
          func = fscale/(2._mp*sqrt_pi)*safe_exp(-((x-ef)/fscale)**2)
      case(FUNCTION_FERMIFUNCTION)
          func = 1._mp/(1._mp+safe_exp((x-ef)/fscale))
      case(FUNCTION_FERMIFUNCTION_ENTROPY)
          ! We must be careful: This function is only properly defined in the interval (0:1),
          ! therefore set it manually to zero outside (including a small safety interval)
          if (x<1.e-30_mp .or. (x-1._mp)>-1.e-30_mp) then
              func = 0._mp
          else
              func = -(x*log(x) + (1.0_mp-x)*log(1._mp-x))
          end if
      case default
          call f_err_throw("wrong value of 'ifunc'")
      end select
    end function func



    !> Calculates the error function complement with an error of less than 1.2E-7
    function erfcc(x)
      implicit none

      ! Calling arguments
      real(8),intent(in) :: x
      real(8) :: erfcc

      ! Local variables
      real(8) :: z, t

      z=abs(x)
      t=1.d0/(1.+0.5d0*z)
      erfcc=t*safe_exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
            & t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
            & t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
      if (x.lt.0.) erfcc=2.D0-erfcc

    end function erfcc

end module module_func


module foe_common
  use foe_base
  use sparsematrix_base
  use dictionaries, only: f_err_throw
  use wrapper_mpi
  use yaml_output
  use numerics
  use f_utils
  use time_profiling
  implicit none

  private

  !> Public routines
  public :: get_chebyshev_expansion_coefficients
  public :: chder
  public :: evnoise
  public :: pltwght
  public :: pltexp
  public :: retransform_ext
  public :: init_foe
  public :: find_fermi_level
  !!public :: get_polynomial_degree
  public :: calculate_trace_distributed_new
  public :: get_bounds_and_polynomials


  contains


    subroutine get_chebyshev_expansion_coefficients(iproc, nproc, comm, A, B, N, func, cc, x_max_error,max_error,mean_error)
      use yaml_output
      use sparsematrix_init, only: distribute_on_tasks
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, n
      real(kind=mp),intent(in) :: A, B
      real(kind=mp),external :: func
      real(8),dimension(n),intent(out) :: cc
      real(kind=mp),intent(out) :: x_max_error, max_error, mean_error

      ! Local variables
      integer :: is, np

      call f_routine(id='get_chebyshev_expansion_coefficients')

      call f_routine(id='synchronize_mpi_tasks')
      call fmpi_barrier(comm)
      call f_release_routine()

      ! MPI parallelization... maybe only worth for large n?
      !call chebyshev_coefficients_init_parallelization(iproc, nproc, comm, n, np, is)
      ! Initialize the parallelization.
      call distribute_on_tasks(n, iproc, nproc, np, is)

      call chebyshev_coefficients_calculate(n, a, b, np, is, func, cc)

      call chebyshev_coefficients_communicate(comm, n, cc)

      call accuracy_of_chebyshev_expansion(iproc, nproc, comm, n, cc, A,B, &
           1.d-3, func, x_max_error, max_error, mean_error)

      call f_release_routine()

    end subroutine get_chebyshev_expansion_coefficients

    !!! Calculates chebychev expansion of fermi distribution.
    !!! Taken from numerical receipes: press et al
    !!subroutine chebft(iproc,nproc,A,B,N,cc,ef,fscale,tmprtr,x_max_error,max_error,mean_error)
    !!  use module_base
    !!  use module_func
    !!  use yaml_output
    !!  implicit none
    !!
    !!  ! Calling arguments
    !!  real(kind=mp),intent(in) :: A, B, ef, fscale, tmprtr
    !!  integer,intent(in) :: iproc, nproc, n
    !!  real(8),dimension(n),intent(out) :: cc
    !!  real(kind=mp),intent(out) :: x_max_error, max_error, mean_error
    !!
    !!  ! Local variables
    !!  integer :: k, j, is, np, ii, jj
    !!  real(kind=mp) :: bma, bpa, y, arg, fac, tt
    !!  real(kind=mp),dimension(50000) :: cf
    !!  !real(kind=mp),parameter :: pi=4.d0*atan(1.d0)
    !!
    !!  call f_routine(id='chebft')

    !!  if (tmprtr/=0.d0) call f_err_throw('tmprtr should be zero for the moment')

    !!  ! MPI parallelization... maybe only worth for large n?
    !!  ii = n/nproc
    !!  np = ii
    !!  is = iproc*ii
    !!  ii = n - nproc*ii
    !!  if (iproc<ii) then
    !!      np = np + 1
    !!  end if
    !!  is = is + min(iproc,ii)

    !!  !write(*,*) 'iproc, nproc, is, np, n', iproc, nproc, is, np, n
    !!  call f_zero(cc)
    !!
    !!  if (n>50000) stop 'chebft'
    !!  bma=0.5d0*(b-a)
    !!  bpa=0.5d0*(b+a)
    !!  fac=2.d0/n
    !!  !$omp parallel default(none) shared(bma,bpa,fac,n,tmprtr,cf,fscale,ef,cc,is,np,tt) &
    !!  !$omp private(k,y,arg,j,jj)
    !!  !$omp do
    !!  do k=1,n
    !!      y=cos(pi*(k-0.5d0)*(1.d0/n))
    !!      arg=y*bma+bpa
    !!      if (tmprtr.eq.0.d0) then
    !!          cf(k)=.5d0*erfcc((arg-ef)*(1.d0/fscale))
    !!      else
    !!          cf(k)=1.d0/(1.d0+safe_exp( (arg-ef)*(1.d0/tmprtr) ) )
    !!      end if
    !!  end do
    !!  !$omp end do
    !!  !$omp end parallel
    !!  do j=1,np
    !!      jj = j + is
    !!      tt=0.d0
    !!      !$omp parallel do default(none) shared(n,cf,jj) private(k) reduction(+:tt)
    !!      do  k=1,n
    !!          tt=tt+cf(k)*cos((pi*(jj-1))*((k-0.5d0)*(1.d0/n)))
    !!      end do
    !!      !$omp end parallel do
    !!      cc(jj)=fac*tt
    !!  end do

    !!  call fmpi_allreduce(cc, FMPI_SUM, comm=bigdft_mpi%mpi_comm)

    !!  call func_set(FUNCTION_ERRORFUNCTION, efx=ef, fscalex=fscale)
    !!  call accuracy_of_chebyshev_expansion(n, cc, (/A,B/), 1.d-3, func, x_max_error, max_error, mean_error)
    !!  !if (bigdft_mpi%iproc==0) call yaml_map('expected accuracy of Chebyshev expansion',max_error)
    !!
    !!  call f_release_routine()
    !!
    !!end subroutine chebft



    !!!! Calculates chebychev expansion of fermi distribution.
    !!!! Taken from numerical receipes: press et al
    !!!subroutine chebyshev_coefficients_penalyfunction(a,b,n,cc,max_error)
    !!!  use module_base
    !!!  use module_func
    !!!  implicit none
    !!!
    !!!  ! Calling arguments
    !!!  real(kind=mp),intent(in) :: a, b
    !!!  integer,intent(in) :: n
    !!!  real(kind=mp),dimension(n,2),intent(out) :: cc
    !!!  real(kind=mp),intent(out) :: max_error
    !!!
    !!!  ! Local variables
    !!!  integer :: k, j
    !!!  !real(kind=mp),parameter :: pi=4.d0*atan(1.d0)
    !!!  real(kind=mp) :: tt1, tt2, ttt, y, arg, fac, bma, bpa, x_max, max_err, mean_err
    !!!  real(kind=mp),dimension(50000) :: cf
    !!!
    !!!  call f_routine(id='chebyshev_coefficients_penalyfunction')
    !!!
    !!!  if (n>50000) stop 'chebyshev_coefficients_penalyfunction'
    !!!  bma=0.5d0*(b-a)
    !!!  bpa=0.5d0*(b+a)
    !!!  ! 3 gives broder safety zone than 4
    !!!  !ttt=3.0d0*n/(b-a)
    !!!  !ttt=4.d0*n/(b-a)
    !!!  ttt=40.d0
    !!!  fac=2.d0/n
    !!!  !$omp parallel default(none) shared(bma,bpa,ttt,fac,n,cf,a,b,cc) &
    !!!  !$omp private(k,y,arg,tt1,tt2,j)
    !!!  !$omp do
    !!!  do k=1,n
    !!!      y=cos(pi*(k-0.5d0)*(1.d0/n))
    !!!      arg=y*bma+bpa
    !!!      !write(*,*) 'arg, safe_exp(-(arg-a)*ttt)', arg, safe_exp(-(arg-a)*ttt)
    !!!      cf(k)= safe_exp(-(arg-a)*ttt)-safe_exp((arg-b)*ttt)
    !!!      !cf(k,2)=-safe_exp(-(arg-a)*ttt)+safe_exp((arg-b)*ttt)
    !!!      !cf(k,1)= safe_exp(-(arg-a)*ttt)
    !!!      !cf(k,2)= safe_exp((arg-b)*ttt)
    !!!  end do
    !!!  !$omp end do
    !!!  !$omp do
    !!!  do j=1,n
    !!!      tt1=0.d0
    !!!      tt2=0.d0
    !!!      do k=1,n
    !!!          tt1=tt1+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
    !!!          !tt2=tt2+cf(k,2)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
    !!!      end do
    !!!      cc(j,1)=fac*tt1
    !!!      !cc(j,2)=fac*tt2
    !!!  end do
    !!!  !$omp end do
    !!!  !$omp do
    !!!  do j=1,n
    !!!      cc(j,2) = -cc(j,1)
    !!!  end do
    !!!  !$omp end do
    !!!  !$omp end parallel
    !!!
    !!!  !!do j=1,n
    !!!  !!    write(*,*) 'j, cc(j,1), cc(j,2)', j, cc(j,1), cc(j,2)
    !!!  !!end do
    !!!  call func_set(FUNCTION_EXPONENTIAL, betax=-ttt, muax=a, mubx=b)
    !!!  call accuracy_of_chebyshev_expansion(n, cc(:,1), (/A,B/), 1.d-3, func, x_max, max_err, mean_err)
    !!!  max_error = max_err
    !!!  call func_set(FUNCTION_EXPONENTIAL, betax=ttt, muax=b, mubx=a)
    !!!  call accuracy_of_chebyshev_expansion(n, cc(:,2), (/A,B/), 1.d-3, func, x_max, max_err, mean_err)
    !!!  max_error = max(max_error,max_err)
    !!!  call f_release_routine()
    !!!
    !!!end subroutine chebyshev_coefficients_penalyfunction

    ! Calculates chebychev expansion of the derivative of Fermi distribution.
    subroutine chder(a,b,c,cder,n)
      use dynamic_memory
      implicit none

      ! Calling arguments
      real(kind=mp),intent(in) :: a, b
      integer,intent(in) :: n
      real(8),dimension(n),intent(in) :: c
      real(8),dimension(n),intent(out) :: cder

      ! Local variables
      integer :: j
      real(kind=mp) :: con

      call f_routine(id='chder')

      cder(n)=0.d0
      cder(n-1)=2*(n-1)*c(n)
      if(n>=3)then
          do j=n-2,1,-1
            cder(j)=cder(j+2)+2*j*c(j+1)
          end do
      end if
      con=2.d0/(b-a)
      do j=1,n
          cder(j)=cder(j)*con
      end do

      call f_release_routine()

    end subroutine chder


    !> Determine noise level
    subroutine evnoise(npl,cc,evlow,evhigh,anoise)
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: npl
      real(kind=mp),dimension(npl),intent(in) :: cc
      real(kind=mp),intent(in) :: evlow, evhigh
      real(kind=mp),intent(out) :: anoise

      ! Local variables
      integer :: i, n
      real(kind=mp) :: fact, dist, ddx, cent, tt, x

      call f_routine(id='evnoise')

      fact=1.d0
      dist=(fact*evhigh-fact*evlow)
      ddx=dist/(10*npl)
      cent=.5d0*(fact*evhigh+fact*evlow)
      !!tt=abs(chebev(evlow,evhigh,npl,cent,cc))
      !!do x=ddx,.25d0*dist,ddx
      !!    tt=max(tt,abs(chebev(evlow,evhigh,npl,cent+x,cc)), &
      !!       & abs(chebev(evlow,evhigh,npl,cent-x,cc)))
      !!end do
      ! Rewritten version of the above loop
      tt=abs(chebev(evlow,evhigh,npl,cent,cc))
      x=ddx
      n=ceiling((0.25d0*dist-ddx)/ddx)
      !$omp parallel default(none) shared(n,ddx,tt,evlow,evhigh,npl,cent,cc) private(i,x)
      !$omp do reduction(max:tt)
      do i=1,n
          x=real(i,kind=mp)*ddx
          tt=max(tt,abs(chebev(evlow,evhigh,npl,cent+x,cc)), &
             & abs(chebev(evlow,evhigh,npl,cent-x,cc)))
          !x=x+ddx
          !if (x>=.25d0*dist) exit
      end do
      !$omp end do
      !$omp end parallel
      !anoise=1.d0*tt
      anoise=20.d0*tt

      call f_release_routine()

    end subroutine evnoise





    !> Evaluates chebychev expansion
    function chebev(a,b,m,x,cc)
      implicit none

      ! Calling arguments
      real(kind=mp),intent(in) :: a, b, x
      integer,intent(in) :: m
      real(kind=mp),dimension(m),intent(in) :: cc
      real(kind=mp) :: chebev

      ! Local variables
      integer :: j
      real(kind=mp) :: d, dd, y, sv

      d=0.d0
      dd=0.d0
      y=2.d0*(2.d0*x-a-b)/(b-a)
      do j=m,2,-1
          sv=d
          d=y*d-dd+cc(j)
          dd=sv
      end do
      chebev= -dd + 0.5d0*(y*d+cc(1))

    end function chebev




    ! plots the approximate fermi distribution
            subroutine pltwght(npl,cc,cder,evlow,evhigh,ef,fscale,tmprtr)
              implicit none

              ! Calling arguments
              integer,intent(in) :: npl
              real(kind=mp),dimension(npl),intent(in) :: cc, cder
              real(kind=mp),intent(in) :: evlow, evhigh, ef, fscale, tmprtr

              ! Local variables
              integer :: ic
              real(kind=mp) :: ddx, x, tt, err

            open (unit=66,file='fermi',status='unknown')
    !     header for my favourite plotting program
            write(66,*) ' 3'
            write(66,*) ' #LINETYPE{132}'
    65        format(a,f5.2,a,i3,a)
            write(66,65) ' #TITLE{WEIGHT DISTR. for fscale=', fscale,' npl=',npl,'}'
            write(66,*) ' #XCAPT{ENERGY in eV}'
            write(66,*) ' #XYAPT{WEIGHT DISTR.}'
            write(66,*) ' #2YAXIS{2}'
            write(66,*) ' #YLOGSCALE2'
            write(66,*) ' #YCAPT2{ERROR}'
            write(66,*) ' $'
    !
    ! plot chebechev expansion of weight distribution function
    !
            ddx=(evhigh-evlow)/(10*npl)
    ! number of plot p[oints
            ic=0
            !!do x=evlow,evhigh,ddx
            !!    ic=ic+1
            !!end do
            x=evlow
            do
                ic=ic+1
                x=x+ddx
                if (x>=evhigh) exit
            end do
    ! weight distribution
            write(66,*) ic
            !!do x=evlow,evhigh,ddx
            !!    write(66,*) x,CHEBEV(evlow,evhigh,npl,x,cc)
            !!end do
            x=evlow
            do
                write(66,*) x,CHEBEV(evlow,evhigh,npl,x,cc)
                x=x+ddx
                if (x>=evhigh) exit
            end do
    ! derivative
            write(66,*) ic
            !!do x=evlow,evhigh,ddx
            !!    write(66,*) x,-CHEBEV(evlow,evhigh,npl,x,cder)
            !!end do
            x=evlow
            do
                write(66,*) x,-CHEBEV(evlow,evhigh,npl,x,cder)
                x=x+ddx
                if (x>=evhigh) exit
            end do
    ! error
            write(66,*) ic
            !!do x=evlow,evhigh,ddx
            !!    tt=tmprtr
            !!    if (tmprtr.eq.0.d0) tt=1.d-16
            !!    err=CHEBEV(evlow,evhigh,npl,x,cc) -1.d0/(1.d0+exp((x-ef)/tt))
            !!    write(66,*) x,err
            !!end do
            x=evlow
            do
                tt=tmprtr
                if (tmprtr.eq.0.d0) tt=1.d-16
                err=CHEBEV(evlow,evhigh,npl,x,cc) -1.d0/(1.d0+exp((x-ef)/tt))
                write(66,*) x,err
                x=x+ddx
                if (x>=evhigh) exit
            end do

            close(unit=66)
    end subroutine pltwght




    ! plots the approximate fermi distribution
    subroutine pltexp(anoise,npl,cc,evlow,evhigh)
            implicit none

            ! Calling arguments
            integer,intent(in) :: npl
            real(kind=mp),dimension(npl),intent(in) :: cc
            real(kind=mp),intent(in) :: anoise, evlow, evhigh

            ! Local variables
            integer :: ic
            real(kind=mp) :: fact, ddx, tt, x

            open (unit=66,file='exp',status='unknown')
    !     header for my favourite plotting program
            write(66,*) ' 2'
            write(66,*) ' #LINETYPE{12}'
            write(66,*) ' #TITLE{exponential}'
            write(66,*) ' #YLOGSCALE'
            write(66,*) ' #XCAPT{ENERGY in eV}'
            write(66,*) ' $'
    !
            fact=1.25d0
    ! plot chebechev expansion of weight distribution function
    !
            ddx=(fact*evhigh-fact*evlow)/(10*npl)
    ! number of plot p[oints
            ic=0
            !!do x=fact*evlow,fact*evhigh,ddx
            !!    ic=ic+1
            !!end do
            x=fact*evlow
            do
                ic=ic+1
                x=x+ddx
                if (x>=fact*evhigh) exit
            end do
    ! first curve
            write(66,*) ic
            !!do x=fact*evlow,fact*evhigh,ddx
            !!    tt=CHEBEV(evlow,evhigh,npl,x,cc)
            !!    if (abs(tt).lt.anoise) tt=anoise
            !!    write(66,*) x,tt
            !!end do
            x=fact*evlow
            do
                tt=CHEBEV(evlow,evhigh,npl,x,cc)
                if (abs(tt).lt.anoise) tt=anoise
                write(66,*) x,tt
                x=x+ddx
                if (x>=fact*evhigh) exit
            end do
    ! second curve
            write(66,*) ic
            !!do x=fact*evhigh,fact*evlow,-ddx
            !!    tt=CHEBEV(evlow,evhigh,npl,x,cc)
            !!    if (abs(tt).lt.anoise) tt=anoise
            !!    write(66,*) fact*evhigh-(x-fact*evlow),tt
            !!end do
            x=fact*evhigh
            do
                tt=CHEBEV(evlow,evhigh,npl,x,cc)
                if (abs(tt).lt.anoise) tt=anoise
                write(66,*) fact*evhigh-(x-fact*evlow),tt
                x=x-ddx
                if (x<=fact*evlow) exit
            end do

            close(unit=66)
    end subroutine pltexp



    subroutine check_eigenvalue_spectrum_new(iproc, nproc, comm, smat_l, ispin, isshift, &
               factor_high, factor_low, penalty_ev, anoise, trace_with_overlap, &
               emergency_stop, foe_obj, restart, eval_bounds_ok, &
               verbosity, eval_multiplicator, smat_s, mat)
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: transform_sparse_matrix
      use yaml_output
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, ispin, isshift
      type(sparse_matrix),intent(in) :: smat_l
      real(kind=mp),intent(in) :: factor_high, factor_low, anoise
      !real(kind=mp),dimension(smat_l%nfvctr,smat_l%smmm%nfvctrp,2),intent(in) :: penalty_ev
      real(kind=mp),dimension(smat_l%smmm%nvctrp),intent(in) :: penalty_ev
      logical,intent(in) :: trace_with_overlap
      logical,dimension(2) :: emergency_stop
      type(foe_data),intent(inout) :: foe_obj
      logical,intent(inout) :: restart
      logical,dimension(2),intent(out) :: eval_bounds_ok
      integer,intent(in),optional :: verbosity
      real(kind=mp),intent(inout),optional :: eval_multiplicator
      type(sparse_matrix),intent(in),optional :: smat_s
      type(matrices),intent(in),optional :: mat

      ! Local variables
      integer :: ii, iismall, i, iline, icolumn, verbosity_
      integer :: ishift
      real(mp),dimension(:),allocatable :: mat_large
      real(kind=mp) :: tt, noise, penalty
      !real(kind=mp),dimension(1) :: allredarr

      call f_routine(id='check_eigenvalue_spectrum_new')

      if (.not.smat_l%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      if (present(verbosity)) then
          verbosity_ = verbosity
      else
          verbosity_ = 1
      end if

      if (trace_with_overlap) then
          if (.not.present(smat_s) .or. .not.present(mat)) then
              call f_err_throw('not all required optional arguments are present')
          end if
          mat_large = sparsematrix_malloc(smat_l, iaction=sparse_taskgroup, id='mat_large')
          call transform_sparse_matrix(iproc, smat_s, smat_l, sparse_taskgroup, 'small_to_large', &
               smat_in=mat%matrix_compr, lmat_out=mat_large)
      end if

      penalty=0.d0
      ishift = (ispin-1)*smat_l%nvctrp_tg
      !$omp parallel if (smat_l%smmm%nvctrp>1000) &
      !$omp default(none) &
      !$omp shared(smat_l, penalty, trace_with_overlap, mat_large, ishift, penalty_ev) &
      !$omp private(i, ii, iline, icolumn, iismall, tt)
      !$omp do schedule(static) reduction(+:penalty)
      do i=1,smat_l%smmm%nvctrp
          ii = smat_l%smmm%isvctr + i
          iline = smat_l%smmm%line_and_column(1,i)
          icolumn = smat_l%smmm%line_and_column(2,i)
          if (trace_with_overlap) then
              ! Take the trace of the product matrix times overlap
              tt = mat_large(ishift+i)
          else
              ! Take the trace of the matrix alone, i.e. set the second matrix to the identity
              if (iline==icolumn) then
                  tt=1.d0
              else
                  tt=0.d0
              end if
          end if
          penalty = penalty + penalty_ev(i)*tt
      end do
      !$omp end do
      !$omp end parallel

      if (trace_with_overlap) then
          call f_free(mat_large)
      end if

      ! Divide the traces by the matrix dimension, to make them size independent
      penalty = penalty/real(smat_l%nfvctr,kind=mp)

      !!if (nproc > 1) then
      !!    call fmpi_allreduce(penalty, FMPI_SUM, comm=comm)
      !!end if
      call penalty_communicate(nproc, comm, penalty)


      !noise=10.d0*anoise
      noise = 1.d-5

      if (iproc==0 .and. verbosity_>0) then
          !call yaml_map('errors, noise',(/allredarr(1),allredarr(2),noise/),fmt='(es12.4)')
          !call yaml_map('pnlty',(/allredarr(1),allredarr(2)/),fmt='(es8.1)')
          call yaml_map('penalty',penalty,fmt='(es8.1)')
      end if

      eval_bounds_ok(1) = .true.
      eval_bounds_ok(2) = .true.
      if (any((/abs(penalty)>noise/))) then
          if (all((/abs(penalty)>noise/))) then
              if (penalty>0.d0) then
                  ! lower bound too large
                  eval_bounds_ok(1)=.false.
                  call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",1)*factor_low,1)
                  restart=.true.
                  if (present(eval_multiplicator)) then
                      !eval_multiplicator = eval_multiplicator*2.0d0
                      eval_multiplicator = 2.0d0
                  end if
              else if (penalty<0.d0) then
                  ! upper bound too small
                  eval_bounds_ok(2)=.false.
                  call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",1)*factor_high,1)
                  restart=.true.
                  if (present(eval_multiplicator)) then
                      !eval_multiplicator = eval_multiplicator/2.0d0
                      eval_multiplicator = 1.d0/2.0d0
                  end if
              else
                  call f_err_throw('The errors should have opposite signs')
              end if
          else
              call f_err_throw('The errors should have the same magnitude')
          end if
      end if

      call f_release_routine()

    end subroutine check_eigenvalue_spectrum_new


    subroutine scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, smatl, &
               smat1, mat1, i1shift, smat2, mat2, i2shift, &
               matscal_compr, scale_factor, shift_value)
      use sparsematrix_init, only: matrixindex_in_compressed
      use sparsematrix, only: transform_sparse_matrix
      use dynamic_memory
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, ispin, i1shift
      type(foe_data),intent(in) :: foe_obj
      type(sparse_matrix),intent(in) :: smatl, smat1
      type(matrices),intent(in) :: mat1
      type(sparse_matrix),intent(in),optional :: smat2
      type(matrices),intent(in),optional :: mat2
      integer,intent(in),optional :: i2shift
      real(kind=mp),dimension(smatl%nvctrp_tg),intent(out) :: matscal_compr
      real(kind=mp),intent(out) :: scale_factor, shift_value

      ! Local variables
      integer :: iseg, ii, i, jj
      integer :: j, ishift
      !integer,dimension(2) :: irowcol
      real(kind=mp) :: tt1, tt2
      logical :: with_overlap
      !real(kind=mp),dimension(:),pointer :: matscal_compr_local
      real(kind=mp),dimension(:),allocatable :: mat1_large, mat2_large
      integer,parameter :: ALLGATHERV=51, GET=52, GLOBAL_MATRIX=101, SUBMATRIX=102
      integer,parameter :: comm_strategy=GET
      integer,parameter :: data_strategy=SUBMATRIX!GLOBAL_MATRIX


      call f_routine(id='scale_and_shift_matrix')
      !call timing(iproc,'foe_aux_mcpy  ','ON')
      call f_timing(TCAT_CME_AUXILIARY,'ON')

      if (.not.smatl%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      call f_zero(matscal_compr)

      ! smat2 and mat2 must be present at the same time
      if (all((/present(smat2),present(mat2),present(i2shift)/))) then
          with_overlap = .true.
      else
          if (any((/present(smat2),present(mat2),present(i2shift)/))) then
              stop 'smat2, mat2 and i2shift must be present at the same time'
          end if
          with_overlap = .false.
      end if

      scale_factor=2.d0/(foe_data_get_real(foe_obj,"evhigh",1)-foe_data_get_real(foe_obj,"evlow",1))
      shift_value=.5d0*(foe_data_get_real(foe_obj,"evhigh",1)+foe_data_get_real(foe_obj,"evlow",1))

      if (data_strategy==GLOBAL_MATRIX) then
          stop 'scale_and_shift_matrix: data_strategy=GLOBAL_MATRIX is deprecated'
      else if (data_strategy==SUBMATRIX) then

          ishift = (ispin-1)*smatl%nvctrp_tg

          ! Transform all matrices to the large sparsity pattern.
          ! Takes some memory, but probably faster than the old way...
          mat1_large = sparsematrix_malloc(smatl, iaction=sparse_taskgroup, id='mat1_large')
          !!call transform_sparse_matrix(iproc, smat1, smatl, sparse_taskgroup, 'small_to_large', &
          !!     smat_in=mat1%matrix_compr, lmat_out=mat1_large)
          call transform_sparse_matrix(iproc, smat1, smatl, sparse_matmul_small, 'small_to_large', &
               smat_in=mat1%matrix_compr, lmat_out=mat1_large)
          if (with_overlap) then
              mat2_large = sparsematrix_malloc(smatl, iaction=sparse_taskgroup, id='mat2_large')
              !!call transform_sparse_matrix(iproc, smat2, smatl, sparse_taskgroup, 'small_to_large', &
              !!     smat_in=mat2%matrix_compr, lmat_out=mat2_large)
              call transform_sparse_matrix(iproc, smat2, smatl, sparse_matmul_small, 'small_to_large', &
                   smat_in=mat2%matrix_compr, lmat_out=mat2_large)
              !write(*,*) 'IN scale_and_shift_matrix: ispin, sum(H), sum(S)', &
              !    ispin, sum(mat1_large(ishift+1:ishift+smatl%nvctrp_tg)), &
              !    sum(mat2_large(ishift+1:ishift+smatl%nvctrp_tg))
          end if

          !$omp parallel if (smatl%smmm%istartendseg_mm(2)-smatl%smmm%istartendseg_mm(1)>100) &
          !$omp default(none) &
          !$omp shared(smatl, mat1_large, mat2_large, with_overlap, matscal_compr, scale_factor, shift_value) &
          !$omp shared(ishift) &
          !$omp private(iseg, j, ii, i, jj, tt1, tt2)
          !$omp do schedule(guided)
          do iseg=smatl%smmm%istartendseg_mm(1),smatl%smmm%istartendseg_mm(2)
              ! A segment is always on one line, therefore no double loop
              j = smatl%keyg(1,2,iseg)
              ii=smatl%keyv(iseg)-1
              do i=smatl%keyg(1,1,iseg),smatl%keyg(2,1,iseg)
                  ii = ii + 1
                  jj = ii-smatl%isvctrp_tg
                  tt1=mat1_large(ishift+jj)
                  if (with_overlap) then
                      tt2 = mat2_large(ishift+jj)
                  else
                      if (i==j) then
                          tt2 = 1.d0
                      else
                          tt2 = 0.d0
                      end if
                  end if
                  !!write(200+iproc,*) 'jj, tt1, tt2', jj, tt1, tt2
                  matscal_compr(jj)=scale_factor*(tt1-shift_value*tt2)
              end do
          end do
          !$omp end do
          !$omp end parallel

          call f_free(mat1_large)
          if (with_overlap) then
              call f_free(mat2_large)
          end if


          !!!write(*,*) 'smatl%smmm%istartendseg_mm',smatl%smmm%istartendseg_mm
          !!!$omp parallel default(none) private(ii,i,j,ii2,ii1,tt2,tt1,iseg) &
          !!!$omp shared(matscal_compr,scale_factor,shift_value,i2shift,i1shift,smatl,smat1,smat2,mat1,mat2,with_overlap)
          !!!$omp do
          !!do iseg=smatl%smmm%istartendseg_mm(1),smatl%smmm%istartendseg_mm(2)
          !!    !if (smatl%keyv(min(iseg+1,smatl%nseg))<smatl%smmm%istartend_mm(1)) cycle
          !!    !if (smatl%keyv(iseg)>smatl%smmm%istartend_mm(2)) exit
          !!    ! A segment is always on one line, therefore no double loop
          !!    j = smatl%keyg(1,2,iseg)
          !!    ii=smatl%keyv(iseg)-1
          !!    do i=smatl%keyg(1,1,iseg),smatl%keyg(2,1,iseg)
          !!        ii = ii + 1
          !!        ii1 = matrixindex_in_compressed(smat1, i, j)
          !!        if (ii1>0) then
          !!            tt1=mat1%matrix_compr(i1shift+ii1-smat1%isvctrp_tg)
          !!        else
          !!            tt1=0.d0
          !!        end if
          !!        if (with_overlap) then
          !!            ii2 = matrixindex_in_compressed(smat2, i, j)
          !!            if (ii2>0) then
          !!                tt2=mat2%matrix_compr(i2shift+ii2-smat2%isvctrp_tg)
          !!            else
          !!                tt2=0.d0
          !!            end if
          !!        else
          !!            if (i==j) then
          !!                tt2 = 1.d0
          !!            else
          !!                tt2 = 0.d0
          !!            end if
          !!        end if
          !!        !ii=matrixindex_in_compressed(smatl, i, j)
          !!        !write(*,*) 'i, ii, ii1, tt1, tt2', i, ii, ii1, tt1, tt2, i1shift, smat1%isvctrp_tg, i1shift+ii1-smat1%isvctrp_tg
          !!        write(300+iproc,*) 'jj, tt1, tt2', ii-smatl%isvctrp_tg, tt1, tt2
          !!        matscal_compr(ii-smatl%isvctrp_tg)=scale_factor*(tt1-shift_value*tt2)
          !!    end do
          !!end do
          !!!$omp end do
          !!!$omp end parallel
          !!!call timing(iproc,'foe_aux_mcpy  ','OF')
          call f_timing(TCAT_CME_AUXILIARY,'OF')
      else
          stop 'scale_and_shift_matrix: wrong data strategy'
      end if

      call f_release_routine()

    end subroutine scale_and_shift_matrix


    subroutine retransform_ext(iproc, nproc, smat, onesided_action, kernelpp_work, inv_ovrlp, kernel, &
               matrix_localx, windowsx)
        use sparsematrix, only: sequential_acces_matrix_fast, sequential_acces_matrix_fast2, &
                                compress_matrix_distributed_wrapper, &
                                sparsemm_new, transform_sparsity_pattern
        use dynamic_memory
        implicit none
        ! Calling arguments
        integer,intent(in) :: iproc, nproc, onesided_action
        type(sparse_matrix),intent(in) :: smat
        real(mp),dimension(smat%smmm%nvctrp),intent(inout) :: kernelpp_work
        real(kind=mp),dimension(smat%nvctrp_tg),intent(inout) :: inv_ovrlp
        real(kind=mp),dimension(smat%nvctrp_tg),intent(inout) :: kernel
        real(kind=mp),dimension(:),intent(inout),target,optional :: matrix_localx
        type(fmpi_win),dimension(:),target,intent(inout),optional :: windowsx

        ! Local variables
        real(kind=mp),dimension(:),pointer :: tempp_new, matrix_local
        real(kind=mp),dimension(:),allocatable :: mat_compr_seq
        type(fmpi_win),dimension(:),pointer :: windows

        call f_routine(id='retransform_ext')

        ! Check the arguments
        select case (onesided_action)
        case (ONESIDED_POST,ONESIDED_GATHER)
            if (nproc>1) then
               if (.not.present(windowsx)) call f_err_throw('windowsx not present')
                if (size(windowsx)/=smat%ntaskgroup) then
                    call f_err_throw('size(windowsx)='//trim(yaml_toa(size(windowsx))) //&
                         &' /= smat%ntaskgroup='//trim(yaml_toa(smat%ntaskgroup)))
                end if
                windows => windowsx
            end if
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
            !if (nproc>1) then
                ! Create a window for all taskgroups to which iproc belongs (max 2)
                windows = f_malloc_ptr(smat%ntaskgroup,id='windows')
            !end if
            matrix_local = f_malloc_ptr(smat%smmm%nvctrp_mm,id='matrix_local')
        case default
            call f_err_throw('wrong value for onesided_action')
        end select


    
        if (onesided_action==ONESIDED_POST .or. onesided_action==ONESIDED_FULL) then
            if (.not.smat%smatmul_initialized) then
                call f_err_throw('sparse matrix multiplication not initialized', &
                     err_name='SPARSEMATRIX_RUNTIME_ERROR')
            end if
            tempp_new = f_malloc_ptr(smat%smmm%nvctrp, id='tempp_new')
            mat_compr_seq = sparsematrix_malloc(smat, iaction=SPARSEMM_SEQ, id='mat_compr_seq')
            !!kernel_compr_seq = sparsematrix_malloc(smat, iaction=SPARSEMM_SEQ, id='kernel_compr_seq')
            if (smat%smmm%nvctrp_mm>0) then !to avoid an out of bounds error
                call transform_sparsity_pattern(iproc, smat%nfvctr, smat%smmm%nvctrp_mm, smat%smmm%isvctr_mm, &
                     smat%nseg, smat%keyv, smat%keyg, smat%smmm%line_and_column_mm, &
                     smat%smmm%nvctrp, smat%smmm%isvctr, &
                     smat%smmm%nseg, smat%smmm%keyv, smat%smmm%keyg, smat%smmm%istsegline, &
                     'small_to_large', &
                     matrix_s_in=inv_ovrlp(smat%smmm%isvctr_mm-smat%isvctrp_tg+1), matrix_l_out=kernelpp_work)
            end if
            call sequential_acces_matrix_fast2(smat, kernel, mat_compr_seq)
            call sparsemm_new(iproc, smat, mat_compr_seq, kernelpp_work, tempp_new)
            call sequential_acces_matrix_fast2(smat, inv_ovrlp, mat_compr_seq)
            call sparsemm_new(iproc, smat, mat_compr_seq, tempp_new, kernelpp_work)
            !kernelpp_work = 1.d0
            !write(*,*) 'after calc, iproc, sum(kernelpp_work)', iproc, sum(kernelpp_work)
            call f_free_ptr(tempp_new)
            call f_free(mat_compr_seq)
            call f_zero(kernel)
        end if
        !!if (onesided_action==ONESIDED_GATHER .or. onesided_action==ONESIDED_FULL) then
        !!if (onesided_action==ONESIDED_FULL) then
            !!write(*,*) 'before compress, iproc, sum(kernelpp_work)', iproc, sum(kernelpp_work)
            call compress_matrix_distributed_wrapper(iproc, nproc, smat, SPARSE_MATMUL_LARGE, &
                 kernelpp_work, onesided_action, kernel, &
                 matrix_localx=matrix_local, windowsx=windows)
            !!write(*,*) 'after compress, iproc, sum(kernel)', iproc, sum(kernel)
        !!end if

        if (onesided_action==ONESIDED_FULL) then
            call free_fmpi_win_ptr(windows)
            call f_free_ptr(matrix_local)
        end if

        call f_release_routine()

    end subroutine retransform_ext



    !!##! Calculates chebychev expansion of x**ex, where ex is any value (typically -1, -1/2, 1/2)
    !!##! Taken from numerical receipes: press et al
    !!##subroutine cheb_exp(iproc, nproc, A,B,N,cc,ex,x_max_error,max_error,mean_error)
    !!##  use module_base
    !!##  use module_func
    !!##  use yaml_output
    !!##  implicit none
    !!##
    !!##  ! Calling arguments
    !!##  integer,intent(in) :: iproc, nproc
    !!##  real(kind=mp),intent(in) :: A, B
    !!##  integer,intent(in) :: n
    !!##  real(kind=mp),intent(in) :: ex
    !!##  real(8),dimension(n),intent(out) :: cc
    !!##  real(kind=mp),intent(out) :: x_max_error, max_error,mean_error
    !!##
    !!##  ! Local variables
    !!##  integer :: k, j
    !!##  real(kind=mp) :: bma, bpa, y, arg, fac, tt
    !!##  real(kind=mp),dimension(50000) :: cf
    !!##  !real(kind=mp),parameter :: pi=4.d0*atan(1.d0)
    !!##
    !!##  call f_routine(id='chebft')

    !!##
    !!##  if (n>50000) stop 'chebft'
    !!##  bma=0.5d0*(b-a)
    !!##  bpa=0.5d0*(b+a)
    !!##  fac=2.d0/n
    !!##  !$omp parallel default(none) shared(bma,bpa,fac,n,cf,cc,ex) &
    !!##  !$omp private(k,y,arg,j,tt)
    !!##  !$omp do
    !!##  do k=1,n
    !!##      y=cos(pi*(k-0.5d0)*(1.d0/n))
    !!##      arg=y*bma+bpa
    !!##      cf(k)=arg**ex
    !!##  end do
    !!##  !$omp end do
    !!##  !$omp do
    !!##  do j=1,n
    !!##      tt=0.d0
    !!##      do  k=1,n
    !!##          tt=tt+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
    !!##      end do
    !!##      cc(j)=fac*tt
    !!##  end do
    !!##  !$omp end do
    !!##  !$omp end parallel

    !!##  call func_set(FUNCTION_POLYNOMIAL, powerx=ex)
    !!##  call accuracy_of_chebyshev_expansion(iproc, nproc, n, cc, (/A,B/), 1.d-3, func, x_max_error, max_error, mean_error)
    !!##  !if (bigdft_mpi%iproc==0) call yaml_map('expected accuracy of Chebyshev expansion',max_error)
    !!##
    !!##  call f_release_routine()
    !!##
    !!##end subroutine cheb_exp


    subroutine init_foe(iproc, nproc, nspin, charge, foe_obj, ef, tmprtr, evbounds_nsatur, evboundsshrink_nsatur, &
               evlow, evhigh, fscale, ef_interpol_det, ef_interpol_chargediff, &
               fscale_lowerbound, fscale_upperbound, eval_multiplicator, &
               npl_min, npl_max, npl_stride, betax, ntemp, accuracy_function, accuracy_penalty, occupation_function, &
               adjust_fscale, fscale_ediff_low, fscale_ediff_up)
      use foe_base, only: foe_data, foe_data_set_int, foe_data_set_real, foe_data_set_logical, foe_data_get_real, foe_data_null
      use dynamic_memory
      use chess_base, only: chess_params, chess_input_dict, chess_init
      use dictionaries
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nspin
      real(kind=mp),dimension(nspin),intent(in) :: charge
      type(foe_data),intent(out) :: foe_obj
      integer,intent(in),optional :: evbounds_nsatur
      integer,intent(in),optional :: evboundsshrink_nsatur
      real(kind=mp),intent(in),optional :: ef
      real(kind=mp),intent(in),optional :: evlow
      real(kind=mp),intent(in),optional :: evhigh
      real(kind=mp),intent(in),optional :: fscale
      real(kind=mp),intent(in),optional :: ef_interpol_det
      real(kind=mp),intent(in),optional :: ef_interpol_chargediff
      real(kind=mp),intent(in),optional :: fscale_lowerbound
      real(kind=mp),intent(in),optional :: fscale_upperbound
      real(kind=mp),intent(in),optional :: tmprtr
      real(kind=mp),intent(in),optional :: eval_multiplicator
      integer,intent(in),optional :: npl_min
      integer,intent(in),optional :: npl_max
      integer,intent(in),optional :: npl_stride
      real(kind=mp),intent(in),optional :: betax
      integer,intent(in),optional :: ntemp
      real(kind=mp),intent(in),optional :: accuracy_function
      real(kind=mp),intent(in),optional :: accuracy_penalty
      integer,intent(in),optional :: occupation_function
      logical,intent(in),optional :: adjust_fscale
      real(kind=mp),intent(in),optional :: fscale_ediff_low
      real(kind=mp),intent(in),optional :: fscale_ediff_up

      ! Local variables
      character(len=*), parameter :: subname='init_foe'
      integer :: ispin
      integer :: evbounds_nsatur_
      integer :: evboundsshrink_nsatur_
      real(kind=mp) :: ef_
      real(kind=mp) :: evlow_
      real(kind=mp) :: evhigh_
      real(kind=mp) :: fscale_
      real(kind=mp) :: ef_interpol_det_
      real(kind=mp) :: ef_interpol_chargediff_
      real(kind=mp) :: fscale_lowerbound_
      real(kind=mp) :: fscale_upperbound_
      real(kind=mp) :: tmprtr_
      real(kind=mp) :: eval_multiplicator_
      integer :: npl_min_
      integer :: npl_max_
      integer :: npl_stride_
      real(kind=mp) :: betax_
      integer :: ntemp_
      real(kind=mp) :: accuracy_function_
      real(kind=mp) :: accuracy_penalty_
      integer :: occupation_function_
      logical :: adjust_fscale_
      real(kind=mp) :: fscale_ediff_low_
      real(kind=mp) :: fscale_ediff_up_
      type(chess_params) :: cp
      type(dictionary),pointer :: dict

      !call timing(iproc,'init_matrCompr','ON')
      call f_timing(TCAT_CME_AUXILIARY,'ON')

      ! Define the default values... Is there a way to get them from input_variables_definition.yaml?
      nullify(dict)
      dict => dict_new()
      call chess_input_dict(dict)
      call chess_init(dict, cp)
      call dict_free(dict)


      ef_ = 0.0_mp
      evbounds_nsatur_ = cp%foe%evbounds_nsatur !3
      evboundsshrink_nsatur_ = cp%foe%evboundsshrink_nsatur !4
      evlow_ = cp%foe%eval_range_foe(1) !-0.5_mp
      evhigh_ = cp%foe%eval_range_foe(2) !0.5_mp
      fscale_ = cp%foe%fscale !2.d-2
      ef_interpol_det_ = cp%foe%ef_interpol_det !1.d-12
      ef_interpol_chargediff_ = cp%foe%ef_interpol_chargediff !1.0_mp
      fscale_lowerbound_ = cp%foe%fscale_lowerbound !5.d-3
      fscale_upperbound_ = cp%foe%fscale_upperbound !5.d-2
      occupation_function_ = cp%foe%occupation_function
      adjust_fscale_ = cp%foe%adjust_fscale
      tmprtr_ = 0.0_mp
      eval_multiplicator_ = 1.0_mp
      npl_min_ = 10
      npl_max_ = 5000
      npl_stride_ = 10
      betax_ = -1000.0_mp
      ntemp_ = 4
      accuracy_function_ = 1.e-5_mp
      accuracy_penalty_ = 1.e-5_mp
      fscale_ediff_low_ = 5.e-5
      fscale_ediff_up_ = 1.e-4

      if (present(ef)) ef_ = ef
      if (present(evbounds_nsatur)) evbounds_nsatur_ = evbounds_nsatur
      if (present(evboundsshrink_nsatur)) evboundsshrink_nsatur_ = evboundsshrink_nsatur
      if (present(evlow)) evlow_ = evlow
      if (present(evhigh)) evhigh_ = evhigh
      if (present(fscale)) fscale_ = fscale
      if (present(ef_interpol_det)) ef_interpol_det_ = ef_interpol_det
      if (present(ef_interpol_chargediff)) ef_interpol_chargediff_ = ef_interpol_chargediff
      if (present(fscale_lowerbound)) fscale_lowerbound_ = fscale_lowerbound
      if (present(fscale_upperbound)) fscale_upperbound_ = fscale_upperbound
      if (present(tmprtr)) tmprtr_ = tmprtr
      if (present(eval_multiplicator)) eval_multiplicator_ = eval_multiplicator
      if (present(npl_min)) npl_min_ = npl_min
      if (present(npl_max)) npl_max_ = npl_max
      if (present(npl_stride)) npl_stride_ = npl_stride
      if (present(betax)) betax_ = betax
      if (present(ntemp)) ntemp_ = ntemp
      if (present(accuracy_function)) accuracy_function_ = accuracy_function
      if (present(accuracy_penalty)) accuracy_penalty_ = accuracy_penalty
      if (present(occupation_function)) occupation_function_ = occupation_function
      if (present(adjust_fscale)) adjust_fscale_ = adjust_fscale
      if (present(fscale_ediff_low)) fscale_ediff_low_ = fscale_ediff_low
      if (present(fscale_ediff_low)) fscale_ediff_up_ = fscale_ediff_up
    
      foe_obj = foe_data_null()

      call foe_data_set_real(foe_obj,"ef",ef_)
      call foe_data_set_real(foe_obj,"fscale",fscale_)
      call foe_data_set_real(foe_obj,"ef_interpol_det",ef_interpol_det_)
      call foe_data_set_real(foe_obj,"ef_interpol_chargediff",ef_interpol_chargediff_)
      call foe_data_set_int(foe_obj,"evbounds_isatur",0)
      call foe_data_set_int(foe_obj,"evboundsshrink_isatur",0)
      call foe_data_set_int(foe_obj,"evbounds_nsatur",evbounds_nsatur_)
      call foe_data_set_int(foe_obj,"evboundsshrink_nsatur",evboundsshrink_nsatur_)
      call foe_data_set_real(foe_obj,"fscale_lowerbound",fscale_lowerbound_)
      call foe_data_set_real(foe_obj,"fscale_upperbound",fscale_upperbound_)
      call foe_data_set_real(foe_obj,"tmprtr",tmprtr_)
      call foe_data_set_int(foe_obj,"npl_min",npl_min_)
      call foe_data_set_int(foe_obj,"npl_max",npl_max_)
      call foe_data_set_int(foe_obj,"npl_stride",npl_stride_)
      call foe_data_set_real(foe_obj,"betax",betax_)
      call foe_data_set_int(foe_obj,"ntemp",ntemp_)
      call foe_data_set_real(foe_obj,"accuracy_function",accuracy_function_)
      call foe_data_set_real(foe_obj,"accuracy_penalty",accuracy_penalty_)
      call foe_data_set_int(foe_obj,"occupation_function",occupation_function_)
      call foe_data_set_logical(foe_obj,"adjust_fscale",adjust_fscale_)
      call foe_data_set_real(foe_obj,"fscale_ediff_low",fscale_ediff_low_)
      call foe_data_set_real(foe_obj,"fscale_ediff_up",fscale_ediff_up_)

      foe_obj%charge = f_malloc0_ptr(nspin,id='foe_obj%charge')
      foe_obj%evlow = f_malloc0_ptr(nspin,id='foe_obj%evlow')
      foe_obj%evhigh = f_malloc0_ptr(nspin,id='foe_obj%evhigh')
      foe_obj%bisection_shift = f_malloc0_ptr(nspin,id='foe_obj%bisection_shift')
      foe_obj%eval_multiplicator = f_malloc0_ptr(nspin,id='foe_obj%eval_multiplicator')
      do ispin=1,nspin
          call foe_data_set_real(foe_obj,"charge",charge(ispin),ispin)
          call foe_data_set_real(foe_obj,"evhigh",evhigh_,ispin)
          call foe_data_set_real(foe_obj,"evlow",evlow_,ispin)
          call foe_data_set_real(foe_obj,"bisection_shift",1.d-1,ispin)
          call foe_data_set_real(foe_obj,"eval_multiplicator",eval_multiplicator_,ispin)
      end do

      !call timing(iproc,'init_matrCompr','OF')
      call f_timing(TCAT_CME_AUXILIARY,'OF')


    end subroutine init_foe


    subroutine accuracy_of_chebyshev_expansion(iproc, nproc, comm, npl, coeff, bound_lower, bound_upper, &
               h, func, x_max_error, max_error, mean_error)
      use sparsematrix_init, only: distribute_on_tasks
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, npl
      real(kind=mp),dimension(npl),intent(in) :: coeff
      real(kind=mp),intent(in) :: bound_lower, bound_upper
      real(kind=mp),intent(in) :: h
      real(kind=mp),external :: func
      real(kind=mp),intent(out) :: x_max_error, max_error, mean_error

      ! Local variables
      integer :: isx, iex, i, ipl, n, iimin, iimax, ii, is, np, jproc, ithread, nthread
      real(kind=mp) :: x, xx, val_chebyshev, val_function, xxm1, xxm2, xxx, sigma, tau, error, val_chebyshev1
      real(kind=mp),dimension(:,:),allocatable :: max_errors
      real(kind=mp),dimension(:),allocatable :: max_error_arr, x_max_error_arr
      !$ integer :: omp_get_thread_num, omp_get_max_threads

      call f_routine(id='accuracy_of_chebyshev_expansion')

      sigma = 2.d0/(bound_upper-bound_lower)
      tau = (bound_lower+bound_upper)/2.d0

      isx = ceiling(bound_lower/h)
      iex = floor(bound_upper/h)
      n = iex - isx + 1

      ! MPI parallelization... maybe only worth for large n?
      !!ii = n/nproc
      !!np = ii
      !!is = iproc*ii
      !!ii = n - nproc*ii
      !!if (iproc<ii) then
      !!    np = np + 1
      !!end if
      !!is = is + min(iproc,ii)
      call distribute_on_tasks(n, iproc, nproc, np, is)

      is = is + isx - 1 !shift of the starting index
      !check
      ii = np
      call fmpi_allreduce(ii, 1, FMPI_SUM, comm=comm)
      if (ii/=n) then
          call f_err_throw('wrong partition: n='//trim(yaml_toa(n))//' /= '//trim(yaml_toa(ii))//'=ii &
               &(n='//trim(yaml_toa(n))//', np='//trim(yaml_toa(np))//')')
      end if
      iimin = 1 + is
      call fmpi_allreduce(iimin, 1, FMPI_MIN, comm=comm)
      if (iimin/=isx) then
          call f_err_throw('wrong starting index')
      end if
      iimax = np + is
      call fmpi_allreduce(iimax, 1, FMPI_MAX, comm=comm)
      if (iimax/=iex) then
          call f_err_throw('wrong ending index')
      end if

      nthread = 1
      !$ nthread = omp_get_max_threads()
      max_error_arr = f_malloc(0.to.nthread,id='max_error_arr')
      x_max_error_arr = f_malloc(0.to.nthread,id='x_max_error_arr')

      max_error_arr(:) = 0.d0
      mean_error = 0.d0
      val_chebyshev1 = 0.5d0*coeff(1)

      ithread = 0
      !$omp parallel if (np>1 .and. np*npl>1000) &
      !$omp default(shared) & ! (none) & !solve the problem of func function
      !!$omp shared(np, is, h, sigma, tau, val_chebyshev1, coeff, npl, mean_error, max_error_arr, x_max_error_arr) &
      !$omp private(i, ii, x, xx, val_chebyshev, xxm2, xxm1, ipl, xxx, val_function, error) &
      !$omp firstprivate(ithread)
      !$ ithread = omp_get_thread_num()
      !$omp do reduction(+: mean_error)
      do i=1,np
          ii = i + is
          x = real(ii,kind=mp)*h
          xx = sigma*(x-tau)
          val_chebyshev = val_chebyshev1 + coeff(2)*xx
          xxm2 = 1.d0
          xxm1 = xx
          xx=2.d0*xx
          do ipl=3,npl
              xxx = xx*xxm1 - xxm2
              val_chebyshev = val_chebyshev + coeff(ipl)*xxx
              xxm2 = xxm1
              xxm1 = xxx
          end do
          val_function = func(x)
          error = abs(val_chebyshev-val_function)
          if (error>max_error_arr(ithread)) then
              max_error_arr(ithread) = error
              x_max_error_arr(ithread) = x
          end if
          mean_error = mean_error + error
          !if (abs(bounds(1)-0.15d0)<1.d-1 .and. abs(bounds(2)-30.0d0)<1.d-1 .and. npl==100) then
          !    write(*,*) 'x, val_chebyshev, val_function, max_error', x, val_chebyshev, val_function, max_error
          !end if
          !write(*,*) 'x, val_chebyshev, exp(x-bounds(2))', x, val_chebyshev, exp(x-bounds(2))
      end do
      !$omp end do
      !$omp end parallel
      !write(*,*) 'max_error',max_error
      mean_error = mean_error/real(iex-isx+1,kind=mp)

      ! Get the maximum among the OpenMP threads
      max_error = 0.d0
      do ithread=0,nthread
          if (max_error_arr(ithread)>max_error) then
              max_error = max_error_arr(ithread)
              x_max_error = x_max_error_arr(ithread)
          end if
      end do
      call f_free(max_error_arr)
      call f_free(x_max_error_arr)

      ! Communicate the results... for the maximum in an array since also the position is required
      call fmpi_allreduce(mean_error, 1, FMPI_SUM, comm=comm)
      max_errors = f_malloc0((/1.to.2,0.to.nproc-1/),id='max_errors')
      max_errors(1,iproc) = max_error
      max_errors(2,iproc) = x_max_error
      call fmpi_allreduce(max_errors, FMPI_SUM, comm=comm)
      max_error = 0.d0
      do jproc=0,nproc-1
          if (max_errors(1,jproc)>max_error) then
              max_error = max_errors(1,jproc)
              x_max_error = max_errors(2,jproc)
          end if
      end do
      call f_free(max_errors)

      call f_release_routine()

    end subroutine accuracy_of_chebyshev_expansion


    !!pure function x_power(x, power)
    !!  implicit none
    !!  real(kind=mp),intent(in) :: x
    !!  real(kind=mp),intent(in) :: power
    !!  real(kind=mp) :: x_power
    !!  x_power = x**power
    !!end function x_power


    subroutine get_chebyshev_polynomials(iproc, nproc, comm, itype, foe_verbosity, npl, npl_penalty, smatm, smatl, &
               ham_, workarr_compr, foe_obj, chebyshev_polynomials, ispin, eval_bounds_ok, &
               hamscal_compr, scale_factor, shift_value, smats, ovrlp_, ovrlp_minus_one_half)
      use sparsematrix_init, only: analyze_unbalancing
      use sparsematrix, only: compress_matrix, uncompress_matrix, &
                              transform_sparsity_pattern, sequential_acces_matrix_fast2, sparsemm_new, &
                              compress_matrix_distributed_wrapper
      use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                          foe_data_get_logical
      use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
                             fermilevel_get_real, fermilevel_get_logical
      use chebyshev, only: chebyshev_clean, chebyshev_fast
      use module_func
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, itype, ispin
      integer,intent(in) :: foe_verbosity
      integer,intent(in) :: npl, npl_penalty
      type(sparse_matrix),intent(in) :: smatm, smatl
      type(matrices),intent(in) :: ham_
      real(kind=mp),dimension(smatl%nvctrp_tg),intent(inout) :: workarr_compr
      type(foe_data),intent(inout) :: foe_obj
      !real(kind=mp),dimension(:,:),pointer,intent(inout) :: chebyshev_polynomials
      real(kind=mp),dimension(smatl%smmm%nvctrp_mm,npl),intent(out) :: chebyshev_polynomials
      logical,dimension(2),intent(out) :: eval_bounds_ok
      real(kind=mp),dimension(smatl%nvctrp_tg),intent(out) :: hamscal_compr
      real(kind=mp),intent(out) :: scale_factor, shift_value
      type(sparse_matrix),intent(in),optional :: smats
      type(matrices),intent(in),optional :: ovrlp_
      real(kind=mp),dimension(smatl%nvctrp_tg),intent(in),optional :: ovrlp_minus_one_half

      ! Local variables
      integer :: ipl, it
      integer :: nsize_polynomial
      !integer :: it_shift
      integer,parameter :: nplx=50000
      real(kind=mp),dimension(:,:,:),allocatable :: cc
      !real(kind=mp),dimension(:,:),allocatable ::  fermip_check
      !real(kind=mp),dimension(:,:,:),allocatable :: penalty_ev
      real(kind=mp) :: anoise
      real(kind=mp) :: evlow_old, evhigh_old
      real(kind=mp) :: x_max_error_fake, max_error_fake, mean_error_fake
      !real(kind=mp) :: fscale, fscale_check, fscale_new
      logical :: restart, adjust_lower_bound, adjust_upper_bound, calculate_SHS, with_overlap
      logical,dimension(2) :: emergency_stop
      real(kind=mp),dimension(2) :: efarr
      !real(kind=mp) :: ebs_check, ef, ebsp
      integer :: isshift, imshift
      logical :: evbounds_shrinked
      integer,parameter :: NTEMP_ACCURATE=4
      integer,parameter :: NTEMP_FAST=1
      !real(kind=mp) :: x_max_error, x_max_error_check
      !real(kind=mp) :: mean_error, mean_error_check
      integer,parameter :: SPARSE=1
      integer,parameter :: DENSE=2
      integer,parameter :: imode=SPARSE
      type(fermi_aux) :: f
      !real(kind=mp),dimension(2) :: temparr
      real(kind=mp),dimension(:),allocatable :: penalty_ev_new, ham_eff, mat_seq, matmul_tmp, matrix_local
      real(kind=mp),dimension(:),allocatable :: fermi_new
      !integer :: icalc
      type(fmpi_win), dimension(:),allocatable :: windows
      logical :: measure_unbalance = .true.
      real(kind=mp) :: t1, t2, time
      real(kind=mp),dimension(:,:),allocatable :: vectors_new



      call f_routine(id='get_chebyshev_polynomials')

      !if (iproc==0) call yaml_comment('get Chebyshev polynomials',hfill='~')

      matrix_local = f_malloc(smatl%smmm%nvctrp_mm,id='matrix_local')
      windows = f_malloc(smatl%ntaskgroup,id='windows')


      !call timing(iproc, 'FOE_auxiliary ', 'ON')
      call f_timing(TCAT_CME_AUXILIARY,'ON')

      if (.not.smatl%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      ! Check the arguments
      select case (itype)
      case (1) !standard eigenvalue problem, i.e. the overlap matrix is the identity and is not required
          with_overlap = .false.
      case (2) !generalized eigenvalue problem, i.e. the overlap matrix must be provided
          if (.not.present(smats)) call f_err_throw('smats not present')
          if (.not.present(ovrlp_)) call f_err_throw('ovrlp_ not present')
          if (.not.present(ovrlp_minus_one_half)) call f_err_throw('ovrlp_minus_one_half not present')
          isshift = (ispin-1)*smats%nvctrp_tg
          with_overlap = .true.
      case default
          call f_err_throw('wrong value for itype')
      end select

      imshift = (ispin-1)*smatm%nvctrp_tg

      evbounds_shrinked=.false.



      penalty_ev_new = f_malloc((/smatl%smmm%nvctrp/),id='penalty_ev_new')
      fermi_new = f_malloc((/smatl%smmm%nvctrp/),id='fermi_new')

      ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
      nsize_polynomial = smatl%smmm%nvctrp_mm

      evlow_old=1.d100
      evhigh_old=-1.d100


      ! Don't let this value become too small.
      call foe_data_set_real(foe_obj, &
           "bisection_shift",max(foe_data_get_real(foe_obj,"bisection_shift",1),1.d-4), &
           1)

      efarr(1)=foe_data_get_real(foe_obj,"ef")-foe_data_get_real(foe_obj,"bisection_shift",1)
      efarr(2)=foe_data_get_real(foe_obj,"ef")+foe_data_get_real(foe_obj,"bisection_shift",1)

      call init_fermi_level(foe_data_get_real(foe_obj,"charge",1), foe_data_get_real(foe_obj,"ef"), f, &
           foe_data_get_real(foe_obj,"bisection_shift",1), foe_data_get_real(foe_obj,"ef_interpol_chargediff"), &
           foe_data_get_real(foe_obj,"ef_interpol_det"), foe_verbosity)
      !call foe_data_set_real(foe_obj,"ef",efarr(1),ispin)

      adjust_lower_bound=.true.
      adjust_upper_bound=.true.

      it=0
      eval_bounds_ok=.false.

      it=it+1


      ! Scale the Hamiltonian such that all eigenvalues are in the intervall [-1:1]
      if (with_overlap) then
          call scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, smatl, &
               smatm, ham_, i1shift=imshift, &
               smat2=smats, mat2=ovrlp_, i2shift=isshift, &
               matscal_compr=hamscal_compr, scale_factor=scale_factor, shift_value=shift_value)
      else
          call scale_and_shift_matrix(iproc, nproc, ispin, foe_obj, smatl, &
               smatm, ham_, i1shift=imshift, &
               matscal_compr=hamscal_compr, scale_factor=scale_factor, shift_value=shift_value)
      end if
      calculate_SHS=.true.
      evlow_old=foe_data_get_real(foe_obj,"evlow",1)
      evhigh_old=foe_data_get_real(foe_obj,"evhigh",1)

          if (with_overlap) then
              ham_eff = f_malloc0(smatl%smmm%nvctrp,id='ham_eff')
              if (smatl%smmm%nvctrp>0) then
                  mat_seq = sparsematrix_malloc(smatl, iaction=SPARSEMM_SEQ, id='mat_seq')
                  matmul_tmp = f_malloc0(smatl%smmm%nvctrp,id='matmul_tmp')
                  call prepare_matrix(smatl, ovrlp_minus_one_half, ham_eff)
                  call sequential_acces_matrix_fast2(smatl, hamscal_compr, mat_seq)
                  call sparsemm_new(iproc, smatl, mat_seq, ham_eff, matmul_tmp)
                  call f_zero(ham_eff)
                  call sequential_acces_matrix_fast2(smatl, ovrlp_minus_one_half, mat_seq)
                  call sparsemm_new(iproc, smatl, mat_seq, matmul_tmp, ham_eff)
                  call f_free(mat_seq)
                  call f_free(matmul_tmp)
              end if
              call compress_matrix_distributed_wrapper(iproc, nproc, smatl, SPARSE_MATMUL_LARGE, &
                   ham_eff, ONESIDED_POST, workarr_compr, matrix_localx=matrix_local, windowsx=windows)
              call f_free(ham_eff)
          else
              call f_memcpy(src=hamscal_compr, dest=workarr_compr)
          end if

      if (npl>nplx) then
          call f_err_throw('npl>nplx')
      end if


      cc = f_malloc((/npl,2,1/),id='cc')

      !call timing(iproc, 'FOE_auxiliary ', 'OF')
      call f_timing(TCAT_CME_AUXILIARY,'OF')
      !call timing(iproc, 'chebyshev_coef', 'ON')
      call f_timing(TCAT_CME_COEFFICIENTS,'ON')

      cc = 0.d0
      call func_set(FUNCTION_EXPONENTIAL, betax=foe_data_get_real(foe_obj,"betax"), &
           muax=foe_data_get_real(foe_obj,"evlow",1), mubx=foe_data_get_real(foe_obj,"evhigh",1))
      call get_chebyshev_expansion_coefficients(iproc, nproc, comm, foe_data_get_real(foe_obj,"evlow",1), &
           foe_data_get_real(foe_obj,"evhigh",1), npl, func, cc(1,2,1), &
           x_max_error_fake, max_error_fake, mean_error_fake)
      call evnoise(npl, cc(1,2,1), foe_data_get_real(foe_obj,"evlow",1), &
           foe_data_get_real(foe_obj,"evhigh",1), anoise)
      !call timing(iproc, 'chebyshev_coef', 'OF')
      call f_timing(TCAT_CME_COEFFICIENTS,'OF')
      !call timing(iproc, 'FOE_auxiliary ', 'ON')
      call f_timing(TCAT_CME_AUXILIARY,'ON')


      if (smatl%nspin==1) then
          do ipl=1,npl
              cc(ipl,1,1)=2.d0*cc(ipl,1,1)
              cc(ipl,2,1)=2.d0*cc(ipl,2,1)
              !!cc(ipl,3,1)=2.d0*cc(ipl,3,1)
          end do
      end if


      !call timing(iproc, 'FOE_auxiliary ', 'OF')
      call f_timing(TCAT_CME_AUXILIARY,'OF')

      emergency_stop=.false.
          if (with_overlap) then
              call compress_matrix_distributed_wrapper(iproc, nproc, smatl, SPARSE_MATMUL_LARGE, &
                   ham_eff, ONESIDED_GATHER, workarr_compr, matrix_localx=matrix_local, windowsx=windows)
          end if
          if (measure_unbalance) then
              call fmpi_barrier(comm=comm)
              t1 = mpi_wtime()
          end if
          mat_seq = sparsematrix_malloc(smatl, iaction=SPARSEMM_SEQ, id='mat_seq')
          if (smatl%smmm%nvctrp>0) then
              call sequential_acces_matrix_fast2(smatl, workarr_compr, mat_seq)
          end if
          vectors_new = f_malloc0((/smatl%smmm%nvctrp,4/),id='vectors_new')
          !!write(*,*) 'npl, npl_penalty', npl, npl_penalty
          call chebyshev_clean(iproc, nproc, comm, npl_penalty, cc(1:npl_penalty,1:2,1:1), &
               smatl, workarr_compr, &
               .false., &
               nsize_polynomial, 1, .false., mat_seq, vectors_new, fermi_new, penalty_ev_new, chebyshev_polynomials, &
               emergency_stop)
          !write(*,*) 'sum(penalty_ev_new)',sum(penalty_ev_new)
          !write(*,*) 'cc(:,2,1)', cc(:,2,1)
          if (measure_unbalance) then
              t2 = mpi_wtime()
              time = t2-t1
              !write(1000+iproc,*) time
              call analyze_unbalancing(iproc, nproc, comm, time)
              !!time_min = time
              !!time_max = time
              !!time_ideal = time/real(nproc,kind=8)
              !!call fmpi_allreduce(time_min, 1, FMPI_MIN, comm)
              !!call fmpi_allreduce(time_max, 1, FMPI_MAX, comm)
              !!call fmpi_allreduce(time_ideal, 1, FMPI_SUM, comm)
              !!if (iproc==0) then
              !!    call yaml_newline()
              !!    call yaml_mapping_open('Load unbalancing')
              !!    call yaml_map('Minimal time',time_min,fmt='(es9.2)')
              !!    call yaml_map('Maximal time',time_max,fmt='(es9.2)')
              !!    call yaml_map('Ideal time',time_ideal,fmt='(es9.2)')
              !!    call yaml_map('Unbalancing in %',(time_max-time_ideal)/time_ideal*100._mp,fmt='(f7.2)')
              !!    call yaml_mapping_close()
              !!end if
          end if
      !call timing(iproc, 'FOE_auxiliary ', 'ON')
      call f_timing(TCAT_CME_AUXILIARY,'ON')


      restart=.false.

      ! Check the eigenvalue bounds. Only necessary if calculate_SHS is true
      ! (otherwise this has already been checked in the previous iteration).
      call check_eigenvalue_spectrum_new(iproc, nproc, comm, smatl, ispin, &
            0, 1.0d0, 1.0d0, penalty_ev_new, anoise, .false., emergency_stop, &
            foe_obj, restart, eval_bounds_ok, foe_verbosity)


      !!write(*,*) 'restart',restart

      if (restart) then
          if(evbounds_shrinked) then
              ! this shrink was not good, increase the saturation counter
              call foe_data_set_int(foe_obj,"evboundsshrink_isatur", &
                   foe_data_get_int(foe_obj,"evboundsshrink_isatur")+1)
          end if
          call foe_data_set_int(foe_obj,"evbounds_isatur",0)
          if (iproc==0) then
              call yaml_map('npl calculated',npl_penalty)
          end if
      else
          call chebyshev_clean(iproc, nproc, comm, npl, cc, &
               smatl, workarr_compr, &
               .false., &
               nsize_polynomial, 1, .true., mat_seq, vectors_new, fermi_new, penalty_ev_new, chebyshev_polynomials, &
               emergency_stop, npl_resume=npl_penalty+1)
          if (iproc==0) then
              call yaml_map('npl calculated',npl)
          end if
      end if

      call f_free(cc)

      ! eigenvalue bounds ok
      if (calculate_SHS) then
          call foe_data_set_int(foe_obj,"evbounds_isatur",foe_data_get_int(foe_obj,"evbounds_isatur")+1)
      end if

      call f_free(mat_seq)
      call f_free(vectors_new)

      call f_free(penalty_ev_new)
      call f_free(fermi_new)
      call f_free(matrix_local)
      call free_fmpi_win_arr(windows)


      !call timing(iproc, 'FOE_auxiliary ', 'OF')
      call f_timing(TCAT_CME_AUXILIARY,'OF')

      call f_release_routine()

    end subroutine get_chebyshev_polynomials


    subroutine find_fermi_level(iproc, nproc, comm, npl, chebyshev_polynomials, &
               foe_verbosity, label, smatl, ispin, foe_obj, kernel_, calculate_spin_channels)
      use sparsematrix, only: compress_matrix, uncompress_matrix, &
                              transform_sparsity_pattern, compress_matrix_distributed_wrapper, &
                              max_asymmetry_of_matrix
      use foe_base, only: foe_data, foe_data_set_int, foe_data_get_int, foe_data_set_real, foe_data_get_real, &
                          foe_data_get_logical
      use fermi_level, only: fermi_aux, init_fermi_level, determine_fermi_level, &
                             fermilevel_get_real, fermilevel_get_logical
      use chebyshev, only: chebyshev_fast
      !!use foe_common, only: scale_and_shift_matrix, evnoise, &
      !!                      check_eigenvalue_spectrum_new, retransform_ext, get_chebyshev_expansion_coefficients
      use module_func
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, npl, ispin
      type(sparse_matrix),intent(in) :: smatl
      real(kind=mp),dimension(smatl%smmm%nvctrp_mm,npl,smatl%nspin),intent(in) :: chebyshev_polynomials
      integer,intent(in) :: foe_verbosity
      character(len=*),intent(in) :: label
      type(foe_data),intent(inout) :: foe_obj
      type(matrices),intent(inout) :: kernel_
      logical,dimension(smatl%nspin),intent(in) :: calculate_spin_channels

      ! Local variables
      integer :: ipl, it
      integer :: nsize_polynomial
      integer :: ilshift
      !integer,parameter :: nplx=50000
      real(kind=mp),dimension(:,:,:),allocatable :: cc
      !real(kind=mp),dimension(:,:),allocatable :: fermip_check
      !real(kind=mp),dimension(:,:,:),allocatable :: penalty_ev
      real(kind=mp) :: anoise, sumn, charge_diff
      real(kind=mp) :: evlow_old, evhigh_old, sumn_old, ef_old, tt
      real(kind=mp) :: x_max_error_fake, max_error_fake, mean_error_fake
      real(kind=mp) :: fscale, diff, fscale_new
      !logical :: calculate_SHS
      !logical,dimension(2) :: emergency_stop
      real(kind=mp),dimension(2) :: efarr, sumnarr
      !real(kind=mp),dimension(:),allocatable :: hamscal_compr
      !real(kind=mp),dimension(4,4) :: interpol_matrix
      !real(kind=mp),dimension(4) :: interpol_vector
      real(kind=mp),parameter :: charge_tolerance=1.d-6 ! exit criterion
      logical,dimension(2) :: eval_bounds_ok, bisection_bounds_ok
      real(kind=mp) :: temp_multiplicator, ef
      integer :: info
      logical :: evbounds_shrinked
      !real(kind=mp),parameter :: FSCALE_LOWER_LIMIT=5.d-3
      !real(kind=mp),parameter :: FSCALE_UPPER_LIMIT=5.d-2
      real(kind=mp),parameter :: DEGREE_MULTIPLICATOR_ACCURATE=3.d0
      real(kind=mp),parameter :: DEGREE_MULTIPLICATOR_FAST=2.d0
      real(kind=mp),parameter :: TEMP_MULTIPLICATOR_ACCURATE=1.d0
      real(kind=mp),parameter :: TEMP_MULTIPLICATOR_FAST=1.2d0 !2.d0 !1.2d0
      real(kind=mp),parameter :: CHECK_RATIO=1.25d0
      integer,parameter :: NPL_MIN=100
      !!type(matrices) :: inv_ovrlp
      integer,parameter :: NTEMP_ACCURATE=4
      integer,parameter :: NTEMP_FAST=1
      real(kind=mp) :: degree_multiplicator, x_max_error, max_error
      real(kind=mp) :: mean_error
      integer,parameter :: SPARSE=1
      integer,parameter :: DENSE=2
      integer,parameter :: imode=SPARSE
      type(fermi_aux) :: f
      !real(kind=mp),dimension(2) :: temparr
      real(kind=mp),dimension(:),allocatable :: occupations
      !real(kind=mp),dimension(:,:),allocatable :: penalty_ev_new
      real(kind=mp),dimension(:,:),allocatable :: fermi_small_new
      integer :: jspin
      type(fmpi_win) ,dimension(:), allocatable :: windowsx



      call f_routine(id='find_fermi_level')

      if (.not.smatl%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      !if (iproc==0) call yaml_comment('FOE calculation of kernel',hfill='~')


      !call timing(iproc, 'FOE_auxiliary ', 'ON')
      call f_timing(TCAT_CME_AUXILIARY,'ON')


      evbounds_shrinked=.false.

      windowsx = f_malloc(smatl%ntaskgroup,id='windowsx')

      fermi_small_new = f_malloc((/smatl%smmm%nvctrp_mm,smatl%nspin/),id='fermi_small_new')


      occupations = f_malloc0(smatl%nspin,id='occupations')


      !hamscal_compr = sparsematrix_malloc(smatl, iaction=SPARSE_TASKGROUP, id='hamscal_compr')


      ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
      nsize_polynomial = smatl%smmm%nvctrp_mm



      !ntemp = NTEMP_ACCURATE
      degree_multiplicator = DEGREE_MULTIPLICATOR_ACCURATE
      temp_multiplicator = TEMP_MULTIPLICATOR_ACCURATE

      fscale_new=1.d100


      !!spin_loop: do ispin=1,smatl%nspin

          !isshift=(ispin-1)*smats%nvctrp_tg
          !imshift=(ispin-1)*smatm%nvctrp_tg
          ilshift=(ispin-1)*smatl%nvctrp_tg
          !ilshift2=(ispin-1)*smatl%nvctrp_tg

          !call get_minmax_eigenvalues(iproc, smatm, ham_, imshift, smats, ovrlp_, isshift)


          fscale_new = temp_multiplicator*foe_data_get_real(foe_obj,"fscale")

          !temp_loop: do itemp=1,ntemp

              fscale = fscale_new
              !fscale = max(fscale,FSCALE_LOWER_LIMIT)
              !fscale = min(fscale,FSCALE_UPPER_LIMIT)

              evlow_old=1.d100
              evhigh_old=-1.d100


                  ! Don't let this value become too small.
                  call foe_data_set_real(foe_obj, &
                       "bisection_shift",max(foe_data_get_real(foe_obj,"bisection_shift",1),1.d-4), &
                       1)

                  efarr(1)=foe_data_get_real(foe_obj,"ef")-foe_data_get_real(foe_obj,"bisection_shift",1)
                  efarr(2)=foe_data_get_real(foe_obj,"ef")+foe_data_get_real(foe_obj,"bisection_shift",1)
                  !write(*,*) 'ef, efarr', foe_data_get_real(foe_obj,"ef",ispin), efarr

                  sumnarr(1)=0.d0
                  sumnarr(2)=1.d100
                  call init_fermi_level(foe_data_get_real(foe_obj,"charge",1), foe_data_get_real(foe_obj,"ef"), f, &
                       foe_data_get_real(foe_obj,"bisection_shift",1), foe_data_get_real(foe_obj,"ef_interpol_chargediff"), &
                       foe_data_get_real(foe_obj,"ef_interpol_det"), verbosity=foe_verbosity)
                  call foe_data_set_real(foe_obj,"ef",efarr(1),1)


                  if (iproc==0 .and. foe_verbosity>0) then
                      !if (foe_verbosity>=1) then
                      !    call yaml_sequence_open('FOE to determine density kernel',&
                      !         label='it_foe'//trim(label)//'-'//&
                      !         trim(adjustl(yaml_toa(itemp,fmt='(i2.2)')))//'-'//&
                      !         trim(adjustl(yaml_toa(ispin,fmt='(i2.2)'))))
                      !else
                      !    call yaml_sequence_open('FOE to determine density kernel')
                      !    if (iproc==0) call yaml_comment('FOE calculation of kernel',hfill='-')
                      !end if
                      call yaml_sequence_open('determine Fermi energy')
                  end if



                  it=0
                  eval_bounds_ok=.true.
                  bisection_bounds_ok=.false.
                  main_loop: do

                      it=it+1

                      if (iproc==0 .and. foe_verbosity>0) then
                          call yaml_newline()
                          call yaml_sequence(advance='no')
                          call yaml_mapping_open(flow=.true.)
                          if (foe_verbosity>=1) call yaml_comment('it FOE:'//yaml_toa(it,fmt='(i6)'),hfill='-')
                      end if

                      if (iproc==0) then
                          !if (foe_verbosity>=1) then
                          !    call yaml_map('bisec/eval bounds',&
                          !         (/fermilevel_get_real(f,"efarr(1)"),fermilevel_get_real(f,"efarr(2)"),&
                          !         foe_data_get_real(foe_obj,"evlow",ispin), &
                          !         foe_data_get_real(foe_obj,"evhigh",ispin)/),fmt='(f5.2)')
                          !else
                          !    call yaml_map('eval bounds',&
                          !         (/foe_data_get_real(foe_obj,"evlow",ispin), &
                          !         foe_data_get_real(foe_obj,"evhigh",ispin)/),fmt='(f5.2)')
                          !end if
                          !call yaml_map('pol deg',npl,fmt='(i0)')
                          !if (foe_verbosity>=1) call yaml_map('eF',foe_data_get_real(foe_obj,"ef",ispin),fmt='(es16.9)')
                      end if


                      cc = f_malloc((/npl,1,3/),id='cc')

                      !call timing(iproc, 'FOE_auxiliary ', 'OF')
                      call f_timing(TCAT_CME_AUXILIARY,'OF')
                      !call timing(iproc, 'chebyshev_coef', 'ON')
                      call f_timing(TCAT_CME_COEFFICIENTS,'ON')

                      call func_set(foe_data_get_int(foe_obj,"occupation_function"), &
                           efx=foe_data_get_real(foe_obj,"ef"), fscalex=fscale)
                      call get_chebyshev_expansion_coefficients(iproc, nproc, comm, foe_data_get_real(foe_obj,"evlow",1), &
                           foe_data_get_real(foe_obj,"evhigh",1), npl, func, cc(1,1,1), &
                           x_max_error, max_error, mean_error)
                      call func_set(FUNCTION_EXPONENTIAL, betax=foe_data_get_real(foe_obj,"betax"), &
                           muax=foe_data_get_real(foe_obj,"evlow",1), mubx=foe_data_get_real(foe_obj,"evhigh",1))
                      call get_chebyshev_expansion_coefficients(iproc, nproc, comm, foe_data_get_real(foe_obj,"evlow",1), &
                           foe_data_get_real(foe_obj,"evhigh",1), npl, func, cc(1,1,2), &
                           x_max_error_fake, max_error_fake, mean_error_fake)
                      do ipl=1,npl
                         cc(ipl,1,3) = -cc(ipl,1,2)
                      end do
                      call evnoise(npl, cc(1,1,2), foe_data_get_real(foe_obj,"evlow",1), &
                           foe_data_get_real(foe_obj,"evhigh",1), anoise)


                      !if (iproc==0 .and. foe_verbosity>=1) then
                      !    call yaml_newline()
                      !    call yaml_mapping_open('accuracy (x, max err, mean err)')
                      !    call yaml_map('main',(/x_max_error,max_error,mean_error/),fmt='(es9.2)')
                      !    !call yaml_map('check',(/x_max_error_check,max_error_check,max_error/),fmt='(es9.2)')
                      !    call yaml_mapping_close()
                      !    call yaml_newline()
                      !end if

                      !call timing(iproc, 'chebyshev_coef', 'OF')
                      call f_timing(TCAT_CME_COEFFICIENTS,'OF')
                      !call timing(iproc, 'FOE_auxiliary ', 'ON')
                      call f_timing(TCAT_CME_AUXILIARY,'ON')


                      if (smatl%nspin==1) then
                          !write(*,*) 'ef',foe_data_get_real(foe_obj,"ef",1)
                          !write(*,*) 'fscale',fscale
                          do ipl=1,npl
                              !write(*,*) 'cc, ipl, cc(ipl,1,1)', ipl, cc(ipl,1,1)
                              cc(ipl,1,1)=2.d0*cc(ipl,1,1)
                              cc(ipl,1,2)=2.d0*cc(ipl,1,2)
                              cc(ipl,1,3)=2.d0*cc(ipl,1,3)
                          end do
                      end if


                      sumn = 0.0_mp
                      do jspin=1,smatl%nspin
                      
                          if (.not. calculate_spin_channels(jspin)) cycle

                          !call timing(iproc, 'FOE_auxiliary ', 'OF')
                          call f_timing(TCAT_CME_AUXILIARY,'OF')

                          !if (foe_verbosity>=1 .and. iproc==0) call yaml_map('polynomials','from memory')

                          if (smatl%smmm%nvctrp_mm>0) then
                              call chebyshev_fast(iproc, nproc, nsize_polynomial, npl, &
                                   smatl%nfvctr, smatl%smmm%nfvctrp, &
                                   smatl, chebyshev_polynomials(:,:,jspin), 1, cc, fermi_small_new(:,jspin))
                          else
                              call f_zero(fermi_small_new(:,jspin))
                          end if

                          ilshift=(jspin-1)*smatl%nvctrp_tg
                          !!call compress_matrix_distributed_wrapper(iproc, nproc, smatl, SPARSE_MATMUL_SMALL, &
                          !!     fermi_small_new(:,jspin), ONESIDED_POST, &
                          !!     kernel_%matrix_compr(ilshift+1:), windowsx=windowsx)

                          !call timing(iproc, 'FOE_auxiliary ', 'ON')
                          call f_timing(TCAT_CME_AUXILIARY,'ON')




                          call calculate_trace_distributed_new(iproc, nproc, comm, smatl, fermi_small_new(:,jspin), tt)
                          occupations(jspin) = tt
                          !write(*,*) 'jspin, trace', jspin, tt
                          sumn = sumn + tt
                          !write(*,*) 'sumn',sumn

                      end do

                      call f_free(cc)


                      if (all(eval_bounds_ok) .and. all(bisection_bounds_ok)) then
                          ! Print these informations already now if all entries are true.
                          if (iproc==0) then
                              !!if (foe_verbosity>=1) call yaml_map('eval/bisection bounds ok',&
                              !!     (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                          end if
                      end if

                      if (iproc==0 .and. foe_verbosity>0) then
                          !call yaml_newline()
                          !call yaml_map('iter',it)
                          call yaml_map('eF',foe_data_get_real(foe_obj,"ef"),fmt='(es13.6)')
                          !call yaml_map('bisec bounds ok',&
                          !     (/bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                          if (smatl%nspin>1) call yaml_map('spin occupations',occupations,fmt='(es14.7)')
                          call yaml_map('Tr(K)',sumn,fmt='(es14.7)')
                          call yaml_map('D Tr(K)',sumn-foe_data_get_real(foe_obj,"charge",ispin),fmt='(es9.2)')
                      end if


                      call determine_fermi_level(iproc, f, sumn, ef, info)
                      bisection_bounds_ok(1) = fermilevel_get_logical(f,"bisection_bounds_ok(1)")
                      bisection_bounds_ok(2) = fermilevel_get_logical(f,"bisection_bounds_ok(2)")

                      charge_diff = sumn-foe_data_get_real(foe_obj,"charge",ispin)

                      ! If the charge difference is smaller than the threshold, there is no need to cycle even though we
                      ! are in principle still looking for the bisection bounds.
                      if (info<0 .and. abs(charge_diff)>=charge_tolerance) then
                          if (iproc==0 .and. foe_verbosity>0) then
                              !if (foe_verbosity>=1) call yaml_map('eval/bisection bounds ok',&
                              !     (/eval_bounds_ok(1),eval_bounds_ok(2),bisection_bounds_ok(1),bisection_bounds_ok(2)/))
                              call yaml_mapping_close()
                          end if
                          !!call f_free(cc_check)
                          ! Save the new fermi energy in the foe_obj structure
                          call foe_data_set_real(foe_obj,"ef",ef,ispin)
                          cycle
                      end if

                      !!! Save the new fermi energy and bisection_shift in the foe_obj structure
                      !!call foe_data_set_real(foe_obj,"ef",ef,ispin)
                      !!call foe_data_set_real(foe_obj,"bisection_shift",fermilevel_get_real(f,"bisection_shift"),1)

                      !!ef_old=foe_data_get_real(foe_obj,"ef",ispin)
                      !!sumn_old=sumn




                      if (iproc==0 .and. foe_verbosity>0) then
                          call yaml_mapping_close()
                      end if

                      if (abs(charge_diff)<charge_tolerance) then
                          if (iproc==0 .and. foe_verbosity>0) call yaml_sequence_close()
                          diff=0.d0

                          if (nproc > 1) then
                              call fmpi_allreduce(diff, 1, FMPI_SUM, comm=comm)
                          end if

                          diff=sqrt(diff)
                          !if (iproc==0) call yaml_map('diff from reference kernel',diff,fmt='(es10.3)')
                          exit
                      end if


                      ! Save the new fermi energy and bisection_shift in the foe_obj structure
                      call foe_data_set_real(foe_obj,"ef",ef,ispin)
                      call foe_data_set_real(foe_obj,"bisection_shift",fermilevel_get_real(f,"bisection_shift"),1)

                      ef_old=foe_data_get_real(foe_obj,"ef")
                      sumn_old=sumn

                  end do main_loop

                  if (iproc==0) then
                      call yaml_newline()
                      call yaml_mapping_open('summary',flow=.true.)
                      call yaml_map('nit',it)
                      call yaml_map('eF',foe_data_get_real(foe_obj,"ef"),fmt='(es13.6)')
                      call yaml_map('Tr(K)',sumn,fmt='(es14.7)')
                      call yaml_map('D Tr(K)',sumn-foe_data_get_real(foe_obj,"charge",ispin),fmt='(es9.2)')
                      call yaml_mapping_close()
                  end if
            
         do jspin=1,smatl%nspin
             ilshift=(jspin-1)*smatl%nvctrp_tg
             !!call compress_matrix_distributed_wrapper(iproc, nproc, smatl, SPARSE_MATMUL_SMALL, &
             !!     fermi_small_new(:,jspin), ONESIDED_GATHER, &
             !!     kernel_%matrix_compr(ilshift+1:), windowsx=windowsx)
             call compress_matrix_distributed_wrapper(iproc, nproc, smatl, SPARSE_MATMUL_SMALL, &
                  fermi_small_new(:,jspin), ONESIDED_FULL, &
                  kernel_%matrix_compr(ilshift+1:), windowsx=windowsx)
             !!tt = 0.d0
             !!do i=1,smatl%nfvctr
             !!    tt = tt + kernel_%matrix_compr(ilshift+(i-1)*smatl%nfvctr+1)
             !!end do
             !!write(*,*) 'jspin, tt kernel', jspin, tt
             !!write(*,*) 'ispin, sum(F), sum(K)', &
             !!    ispin, sum(fermi_small_new(:,jspin)), sum(kernel_%matrix_compr(ilshift+1:ilshift+smatl%nvctr))
         end do

      !end do spin_loop




      !if (iproc==0) call yaml_comment('FOE calculation of kernel finished',hfill='~')


      call f_free(occupations)
      call f_free(fermi_small_new)
      call free_fmpi_win_arr(windowsx)

      !call timing(iproc, 'FOE_auxiliary ', 'OF')
      call f_timing(TCAT_CME_AUXILIARY,'OF')

      call f_release_routine()




    end subroutine find_fermi_level


    subroutine calculate_trace_distributed_new(iproc, nproc, comm, smatl, matrixp, trace)
      use dynamic_memory
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      type(sparse_matrix),intent(in) :: smatl
      real(kind=mp),dimension(smatl%smmm%nvctrp_mm),intent(in) :: matrixp
      real(kind=mp),intent(out) :: trace
      integer :: i, iline, icolumn

      call f_routine(id='calculate_trace_distributed_new')

      if (.not.smatl%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      trace = 0.d0
      !$omp parallel default(none) &
      !$omp shared(trace, smatl, matrixp) &
      !$omp private(i, iline, icolumn)
      !$omp do reduction(+:trace)
      do i=1,smatl%smmm%nvctrp_mm
          iline = smatl%smmm%line_and_column_mm(1,i)
          icolumn = smatl%smmm%line_and_column_mm(2,i)
          if (iline==icolumn) then
              trace = trace + matrixp(i)
          end if
      end do
      !$omp end do
      !$omp end parallel

      if (nproc > 1) then
          call fmpi_allreduce(trace, 1, FMPI_SUM, comm=comm)
      end if

      call f_release_routine()
    end subroutine calculate_trace_distributed_new


    !> Determine the polynomial degree which yields the desired precision
    subroutine get_polynomial_degree(iproc, nproc, comm, ispin, ncalc, fun, foe_obj, &
               npl_min, npl_max, npl_stride, accuracy_function, accuracy_penalty, verbosity, npl, npl_penalty, cc, &
               max_error, x_max_error, mean_error, anoise, increase_degree_for_penaltyfunction, &
               ex, ef, fscale)
      use foe_base, only: foe_data, foe_data_get_real
      use yaml_output
      use module_func
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, ispin, ncalc, fun, verbosity
      integer,intent(in) :: npl_min, npl_max, npl_stride
      type(foe_data),intent(in) :: foe_obj
      real(kind=mp),intent(in) :: accuracy_function, accuracy_penalty
      integer,intent(out) :: npl, npl_penalty
      real(kind=mp),dimension(:,:,:),pointer,intent(inout) :: cc
      real(kind=mp),dimension(ncalc),intent(out) :: max_error, x_max_error, mean_error
      real(kind=mp),intent(out) :: anoise
      logical,intent(out) :: increase_degree_for_penaltyfunction
      real(kind=mp),dimension(ncalc),intent(in),optional :: ex, ef, fscale

      ! Local variables
      integer :: ipl, icalc, j, jpl
      logical :: error_ok, found_degree, found_penalty_degree
      real(kind=mp),dimension(:,:,:),allocatable :: cc_trial
      real(kind=mp) :: x_max_error_penaltyfunction, max_error_penaltyfunction, mean_error_penaltyfunction

      call f_routine(id='get_polynomial_degree')

      ! Check the arguments
      select case (fun)
      case (FUNCTION_POLYNOMIAL)
          if (.not. present(ex)) call f_err_throw("arguments 'ex' is not present")
      case (FUNCTION_ERRORFUNCTION, FUNCTION_FERMIFUNCTION)
          if (.not. present(ef)) call f_err_throw("arguments 'ef' is not present")
          if (.not. present(fscale)) call f_err_throw("arguments 'fscale' is not present")
          !write(*,*) 'iproc, ef, fscale, evlow, evhigh', &
          !    iproc, ef, fscale, foe_data_get_real(foe_obj,"evlow",ispin), foe_data_get_real(foe_obj,"evhigh",ispin)
      case (FUNCTION_FERMIFUNCTION_ENTROPY)
      case default
          call f_err_throw("wrong value of argument 'fun'")
      end select

      !max_error = f_malloc(ncalc,id='max_error')
      !x_max_error = f_malloc(ncalc,id='x_max_error')
      !mean_error = f_malloc(ncalc,id='mean_error')

      if (npl_min<3) then
          call f_err_throw('npl_min must be at least 3')
      end if
      if (npl_min>npl_max) then
          call f_err_throw('npl_min must be smaller or equal than npl_max')
      end if

      if (iproc==0 .and. verbosity>0) then
          call yaml_sequence_open('Determine polynomial degree')
      end if

      cc_trial = f_malloc0((/npl_max,ncalc,3/),id='cc_trial')

      found_degree = .false.
      increase_degree_for_penaltyfunction = .false.
      found_penalty_degree = .false.
      npl_penalty = huge(1)
      degree_loop: do ipl=npl_min,npl_max,npl_stride

          if (foe_data_get_real(foe_obj,"evhigh",ispin)<=0.d0) then
              stop 'ERROR: highest eigenvalue must be positive'
          end if

          !call timing(iproc, 'FOE_auxiliary ', 'OF')
          call f_timing(TCAT_CME_AUXILIARY,'OF')
          !call timing(iproc, 'chebyshev_coef', 'ON')
          call f_timing(TCAT_CME_COEFFICIENTS,'ON')

          do icalc=1,ncalc
              select case (fun)
              case (FUNCTION_POLYNOMIAL)
                  call func_set(FUNCTION_POLYNOMIAL, powerx=ex(icalc))
              case (FUNCTION_ERRORFUNCTION)
                  call func_set(FUNCTION_ERRORFUNCTION, efx=ef(icalc), fscalex=fscale(icalc))
              case (FUNCTION_FERMIFUNCTION)
                  call func_set(FUNCTION_FERMIFUNCTION, efx=ef(icalc), fscalex=fscale(icalc))
              case (FUNCTION_FERMIFUNCTION_ENTROPY)
                  call func_set(FUNCTION_FERMIFUNCTION_ENTROPY)
              case default
                  call f_err_throw('wrong value for fun')
              end select
              call get_chebyshev_expansion_coefficients(iproc, nproc, comm, foe_data_get_real(foe_obj,"evlow",ispin), &
                   foe_data_get_real(foe_obj,"evhigh",ispin), ipl, func, cc_trial(1:ipl,icalc,1), &
                   x_max_error(icalc), max_error(icalc), mean_error(icalc))
              !write(*,*) 'icalc, sum(cc_trial(:,icalc,1))', icalc, sum(cc_trial(:,icalc,1)), ex(icalc)
          end do

          !call timing(iproc, 'chebyshev_coef', 'OF')
          call f_timing(TCAT_CME_COEFFICIENTS,'OF')
          !call timing(iproc, 'FOE_auxiliary ', 'ON')
          call f_timing(TCAT_CME_AUXILIARY,'ON')

          if (iproc==0 .and. verbosity>0) then
              call yaml_mapping_open(flow=.true.)
              call yaml_map('ipl',ipl)
              do icalc=1,ncalc
                  !call yaml_map('Operation '//trim(yaml_toa(icalc)), &
                  !    (/x_max_error(icalc),max_error(icalc),mean_error(icalc),max_error_penaltyfunction/),fmt='(es9.2)')
                  !if (iproc==0) write(*,*) 'accuracy_function, accuracy_penalty', accuracy_function, accuracy_penalty
                  call yaml_map('Operation '//trim(yaml_toa(icalc)), &
                      (/x_max_error(icalc),max_error(icalc),mean_error(icalc)/),fmt='(es9.2)')
              end do
              call yaml_mapping_close()
          end if

          ! See whether for the penalty function we have already reached the required precision
          call func_set(FUNCTION_EXPONENTIAL, betax=foe_data_get_real(foe_obj,"betax"), &
               muax=foe_data_get_real(foe_obj,"evlow",ispin), mubx=foe_data_get_real(foe_obj,"evhigh",ispin))
          call get_chebyshev_expansion_coefficients(iproc, nproc, comm, foe_data_get_real(foe_obj,"evlow",ispin), &
               foe_data_get_real(foe_obj,"evhigh",ispin), ipl, func, cc_trial(1:ipl,icalc,2), &
               x_max_error_penaltyfunction, max_error_penaltyfunction, mean_error_penaltyfunction)
          if (.not. found_penalty_degree .and. max_error_penaltyfunction<=accuracy_penalty) then
              npl_penalty = ipl
              !!write(*,*) 'set npl_penalty to ', npl_penalty
              found_penalty_degree = .true.
          end if

          error_ok = .true.
          do icalc=1,ncalc
              if (max_error(icalc)>accuracy_function) then
                  error_ok = .false.
                  exit
              end if
          end do
          if (error_ok) then
              do icalc=1,ncalc
                  call func_set(FUNCTION_EXPONENTIAL, betax=foe_data_get_real(foe_obj,"betax"), &
                       muax=foe_data_get_real(foe_obj,"evlow",ispin), mubx=foe_data_get_real(foe_obj,"evhigh",ispin))
                  call get_chebyshev_expansion_coefficients(iproc, nproc, comm, foe_data_get_real(foe_obj,"evlow",ispin), &
                       foe_data_get_real(foe_obj,"evhigh",ispin), ipl, func, cc_trial(1:ipl,icalc,2), &
                       x_max_error_penaltyfunction, max_error_penaltyfunction, mean_error_penaltyfunction)
                  do jpl=1,ipl
                      cc_trial(jpl,icalc,3) = -cc_trial(jpl,icalc,2)
                  end do
                  if (max_error_penaltyfunction>accuracy_penalty) then
                      error_ok = .false.
                      increase_degree_for_penaltyfunction = .true.
                  end if
              end do
          end if
          if (error_ok) then
              npl = ipl
              found_degree = .true.
              exit degree_loop
          end if


      end do degree_loop

      if (.not.found_degree) then
          if (iproc==0) then
              call yaml_warning('Not possible to reach desired accuracy, using highest available polynomial degree')
          end if
          npl = npl_max
      end if

      if (iproc==0 .and. verbosity>0) then
          call yaml_sequence_close()
      end if

      cc = f_malloc_ptr((/npl,ncalc,3/),id='cc')
      !write(*,*) 'ef',ef
      !write(*,*) 'fscale',fscale
      do j=1,3
          do icalc=1,ncalc
              do ipl=1,npl
                  cc(ipl,icalc,j)=cc_trial(ipl,icalc,j)
                  !if (j==1) write(*,*) 'cc, ipl, cc(ipl,1,1)', ipl, cc(ipl,1,1)
                  !write(*,*) 'icalc, ipl, cc(ipl,icalc,1)', icalc, ipl, cc(ipl,icalc,1)
              end do
          end do
      end do

      ! npl_penalty must be at most npl
      npl_penalty = min(npl_penalty,npl)

      call f_free(cc_trial)
      !call f_free(mean_error)
      !call f_free(max_error)
      !call f_free(x_max_error)

      call f_release_routine

    end subroutine get_polynomial_degree


    subroutine chebyshev_coefficients_init_parallelization(iproc, nproc, comm, n, np, is)
      use dynamic_memory
      implicit none
      ! Caling arguments
      integer,intent(in) :: iproc, nproc, comm, n
      integer,intent(out) :: np, is

      ! Local variables
      integer :: ii

      call f_routine(id='chebyshev_coefficients_init_parallelization')

      ii = n/nproc
      np = ii
      is = iproc*ii
      ii = n - nproc*ii
      if (iproc<ii) then
          np = np + 1
      end if
      is = is + min(iproc,ii)
      !check
      ii = np
      call fmpi_allreduce(ii, 1, FMPI_SUM, comm=comm)
      if (ii/=n) then
          call f_err_throw('wrong partition of n')
      end if

      call f_release_routine()

    end subroutine chebyshev_coefficients_init_parallelization


    subroutine chebyshev_coefficients_calculate(n, a, b, np, is, func, cc)
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: n, np, is
      real(kind=mp),intent(in) :: a, b
      real(kind=mp),external :: func
      real(kind=mp),dimension(n),intent(out) :: cc

      ! Local variables
      integer :: k, j, jj
      real(kind=mp) :: bma, bpa, y, arg, fac, tt, one_over_n
      real(kind=mp),dimension(:),allocatable :: cf

      call f_routine(id='chebyshev_coefficients_calculate')


      call f_zero(cc)
      cf = f_malloc0(n,id='cf')

      bma=0.5d0*(b-a)
      bpa=0.5d0*(b+a)
      fac=2.d0/real(n,kind=mp)
      one_over_n = 1.d0/real(n,kind=mp)
      !$omp parallel default(shared) & !none) 
      !!$omp shared(bma,bpa,fac,n,cf,cc,is,np,tt,one_over_n) !&
      !$omp private(k,y,arg,j,jj)
      !$omp do
      do k=1,n
          y=cos(pi*(real(k,kind=mp)-0.5d0)*(one_over_n))
          arg=y*bma+bpa
          cf(k)=func(arg)
      end do
      !$omp end do
      !$omp end parallel

      do j=1,np
          jj = j + is
          tt=0.d0
          !$omp parallel do default(none) shared(n,cf,jj,one_over_n) private(k) reduction(+:tt)
          do  k=1,n
              tt=tt+cf(k)*cos((pi*real(jj-1,kind=mp))*((real(k,kind=mp)-0.5d0)*(one_over_n)))
          end do
          !$omp end parallel do
          cc(jj)=fac*tt
      end do
      call f_free(cf)

      call f_release_routine()

    end subroutine chebyshev_coefficients_calculate


    ! This routine is basically just here to get the profiling...
    subroutine chebyshev_coefficients_communicate(comm, n, cc)
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: comm, n
      real(kind=mp),dimension(n),intent(inout) :: cc

      call f_routine(id='chebyshev_coefficients_communicate')

      call fmpi_allreduce(cc, FMPI_SUM, comm=comm)

      call f_release_routine()

    end subroutine chebyshev_coefficients_communicate


    ! This routine is basically just here to get the profiling...
    subroutine penalty_communicate(nproc, comm, penalty)
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: nproc, comm
      real(mp),intent(inout) :: penalty

      call f_routine(id='penalty_communicate')

      if (nproc > 1) then
          call fmpi_allreduce(penalty, 1, FMPI_SUM, comm=comm)
      end if

      call f_release_routine()

    end subroutine penalty_communicate



    subroutine get_bounds_and_polynomials(iproc, nproc, comm, itype, ispin, &
               npl_max, npl_stride, ncalc, func_name, accuracy_function, accuracy_penalty, &
               do_scaling, bounds_factor_low, bounds_factor_up, foe_verbosity, &
               smatm, smatl, ham_, foe_obj, npl_min, workarr_compr, chebyshev_polynomials, &
               npl, scale_factor, shift_value, hamscal_compr, &
               smats, ovrlp_, ovrlp_minus_one_half_, efarr, fscale_arr, ex, &
               scaling_factor_low, scaling_factor_up, eval_multiplicator, eval_multiplicator_total, cc, max_errorx)
      use module_func
      use dynamic_memory
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, itype, ispin, npl_max, npl_stride, ncalc, func_name, foe_verbosity
      type(sparse_matrix),intent(in) :: smatm, smatl
      type(matrices),intent(in) :: ham_
      logical,intent(in) :: do_scaling
      real(mp),intent(in) :: accuracy_function, accuracy_penalty
      real(mp),intent(in),optional :: bounds_factor_low, bounds_factor_up
      type(foe_data),intent(inout) :: foe_obj
      integer,intent(inout) :: npl_min
      real(kind=mp),dimension(smatl%nspin*smatl%nvctrp_tg),intent(inout) :: workarr_compr
      real(mp),dimension(:,:,:),pointer,intent(inout) :: chebyshev_polynomials
      integer,intent(out) :: npl
      real(mp),intent(out) :: scale_factor, shift_value
      real(kind=mp),dimension(smatl%nvctrp_tg*smatl%nspin),intent(out) :: hamscal_compr
      type(sparse_matrix),intent(in),optional :: smats
      type(matrices),intent(in),optional :: ovrlp_, ovrlp_minus_one_half_
      real(kind=mp),dimension(ncalc),intent(in),optional :: efarr
      real(kind=mp),dimension(ncalc),intent(in),optional :: fscale_arr
      real(kind=mp),dimension(ncalc),intent(in),optional :: ex
      real(mp),intent(in),optional :: scaling_factor_low, scaling_factor_up
      real(mp),intent(inout),optional :: eval_multiplicator, eval_multiplicator_total
      real(kind=mp),dimension(:,:,:),pointer,intent(out),optional :: cc
      real(mp),dimension(ncalc),intent(out),optional :: max_errorx

      ! Local variables
      integer :: ilshift, i, jspin, imshift, npl_penalty
      real(mp),dimension(:),allocatable :: max_error, x_max_error, mean_error
      real(mp) :: anoise, tt
      real(kind=mp),dimension(:,:,:),pointer :: cc_
      logical,dimension(2) :: eval_bounds_ok, eval_bounds_ok_allspins, scale_matrix
      logical :: increase_degree_for_penaltyfunction
      type(matrices) :: ham_scaled

      call f_routine(id='get_bounds_and_polynomials')

      ! Check the arguments
      select case (itype)
      case (1) !standard eigenvalue problem, i.e. the overlap matrix is the identity and is not required
      case (2) !generalized eigenvalue problem, i.e. the overlap matrix must be provided
          if (.not.present(smats)) call f_err_throw('smats not present')
          if (.not.present(ovrlp_)) call f_err_throw('ovrlp_ not present')
          if (.not.present(ovrlp_minus_one_half_)) call f_err_throw('ovrlp_minus_one_half_ not present')
      case default
          call f_err_throw('wrong value for itype')
      end select

      select case (func_name)
      case (FUNCTION_ERRORFUNCTION, FUNCTION_FERMIFUNCTION) 
          if (.not.present(efarr)) call f_err_throw('efarr not present')
          if (.not.present(fscale_arr)) call f_err_throw('fscale_arr not present')
      case (FUNCTION_POLYNOMIAL) !generalized eigenvalue problem, i.e. the overlap matrix must be provided
          if (.not.present(ex)) call f_err_throw('ex not present')
      case (FUNCTION_FERMIFUNCTION_ENTROPY) 
      case default
          call f_err_throw('wrong value for func_name')
      end select

      if (do_scaling) then
          if (.not.present(scaling_factor_low)) call f_err_throw('scaling_factor_low not present')
          if (.not.present(scaling_factor_up)) call f_err_throw('scaling_factor_up not present')
          if (.not.present(eval_multiplicator)) call f_err_throw('eval_multiplicator not present')
          if (.not.present(eval_multiplicator_total)) call f_err_throw('eval_multiplicator_total not present')
      end if


      max_error = f_malloc(ncalc,id='max_error')
      x_max_error = f_malloc(ncalc,id='x_max_error')
      mean_error = f_malloc(ncalc,id='mean_error')
    
      ham_scaled = matrices_null()
      if (do_scaling) then
          ham_scaled%matrix_compr = sparsematrix_malloc_ptr(smatm, &
              iaction=SPARSE_TASKGROUP, id='ham_scaled%matrix_compr')
          call f_memcpy(src=ham_%matrix_compr,dest=ham_scaled%matrix_compr)
      else
          ham_scaled%matrix_compr => ham_%matrix_compr
      end if
    
    
      if (iproc==0 .and. foe_verbosity>0) then
          call yaml_map('beta for penaltyfunction',foe_data_get_real(foe_obj,"betax"),fmt='(f7.1)')
          call yaml_sequence_open('determine eigenvalue bounds')
      end if
      scale_matrix(1:2) = .false.
      bounds_loop: do
          eval_bounds_ok_allspins(:) = .true.
          if (do_scaling) then
              call dscal(size(ham_scaled%matrix_compr), eval_multiplicator, ham_scaled%matrix_compr(1), 1)
              eval_multiplicator_total = eval_multiplicator_total*eval_multiplicator
          end if
          spin_loop: do jspin=1,smatl%nspin 


              ilshift = (jspin-1)*smatl%nvctrp_tg
              imshift = (jspin-1)*smatm%nvctrp_tg
              !efarr(1) = foe_data_get_real(foe_obj,"ef",ispin)
              !fscale_arr(1) = foe_data_get_real(foe_obj,"fscale",ispin)

              if (jspin==1) then
    
                  if (func_name==FUNCTION_ERRORFUNCTION .or. func_name==FUNCTION_FERMIFUNCTION) then
                      call get_polynomial_degree(iproc, nproc, comm, 1, ncalc, &
                           foe_data_get_int(foe_obj,"occupation_function"), foe_obj, &
                           npl_min, npl_max, npl_stride, accuracy_function, accuracy_penalty, 0, npl, npl_penalty, cc_, &
                           max_error, x_max_error, mean_error, anoise, increase_degree_for_penaltyfunction, &
                           ef=efarr, fscale=fscale_arr)
                  else if (func_name==FUNCTION_POLYNOMIAL) then
                      call get_polynomial_degree(iproc, nproc, comm, 1, ncalc, FUNCTION_POLYNOMIAL, foe_obj, &
                           npl_min, npl_max, npl_stride, accuracy_function, accuracy_penalty, 0, npl, npl_penalty, cc_, &
                           max_error, x_max_error, mean_error, anoise, increase_degree_for_penaltyfunction, &
                           ex=ex)
                  else if (func_name==FUNCTION_FERMIFUNCTION_ENTROPY) then
                      !!write(*,*) 'npl_min, npl_max, npl_stride', npl_min, npl_max, npl_stride
                      call get_polynomial_degree(iproc, nproc, comm, 1, ncalc, func_name, foe_obj, &
                           npl_min, npl_max, npl_stride, accuracy_function, accuracy_penalty, 0, npl, npl_penalty, cc_, &
                           max_error, x_max_error, mean_error, anoise, increase_degree_for_penaltyfunction)
                  end if
                  !!npl_min = npl !to be used to speed up the search for npl in a following iteration
                  npl_min = min(npl_penalty,npl) !to be used to speed up the search for npl in a following iteration
                  if (iproc==0 .and. foe_verbosity>0) then
                      call yaml_newline()
                      call yaml_sequence(advance='no')
                      call yaml_mapping_open(flow=.true.)
                      call yaml_map('npl',npl)
                      call yaml_map('npl penalty',npl_penalty)
                      if (increase_degree_for_penaltyfunction) then
                          call yaml_map('npl determined by','penalty')
                      else
                          call yaml_map('npl determined by','function')
                      end if
                      if (do_scaling) call yaml_map('scale',eval_multiplicator_total,fmt='(es9.2)')
                      call yaml_newline()
                      call yaml_map('bounds', &
                           (/foe_data_get_real(foe_obj,"evlow",1),foe_data_get_real(foe_obj,"evhigh",1)/),fmt='(f7.3)')
                      call yaml_map('exp accur',max_error,fmt='(es8.2)')
                  end if
                  chebyshev_polynomials = f_malloc_ptr((/smatl%smmm%nvctrp_mm,npl,smatl%nspin/),id='chebyshev_polynomials')
              end if

              if (iproc==0 .and. smatl%nspin>1 .and. foe_verbosity>0) then
                  if (jspin==1) then
                      call yaml_mapping_open(trim(yaml_toa('up')))
                  else
                      call yaml_mapping_open(trim(yaml_toa('down')))
                  end if
              end if
    
              if (itype==2) then
                  call get_chebyshev_polynomials(iproc, nproc, comm, &
                       itype, foe_verbosity, npl, npl_penalty, smatm, smatl, &
                       ham_scaled, workarr_compr(ilshift+1:), foe_obj, &
                       chebyshev_polynomials(:,:,jspin), jspin, eval_bounds_ok, hamscal_compr(ilshift+1:), &
                       scale_factor, shift_value, &
                       smats=smats, ovrlp_=ovrlp_, &
                       ovrlp_minus_one_half=ovrlp_minus_one_half_%matrix_compr(ilshift+1:))
              else if (itype==1) then
                  call get_chebyshev_polynomials(iproc, nproc, comm, &
                       itype, foe_verbosity, npl, npl_penalty, smatm, smatl, &
                       ham_scaled, workarr_compr(ilshift+1:), foe_obj, &
                       chebyshev_polynomials(:,:,jspin), jspin, eval_bounds_ok, hamscal_compr(ilshift+1:), &
                       scale_factor, shift_value)
              end if
              !write(*,*) 'jspin, ilshift, sum(abs(hamscal_compr(ilshift+1:)))', &
              !    jspin, ilshift, sum(abs(hamscal_compr(ilshift+1:ilshift+smatl%nvctr)))
              !write(*,*) 'hamscal_compr',hamscal_compr(ilshift+1:smatl%nvctr)
              !write(*,*) 'jspin, sum(chebyshev_polynomials(:,:,jspin))', jspin, sum(chebyshev_polynomials(:,:,jspin))
              if (iproc==0 .and. foe_verbosity>0) then
                  call yaml_map('ok',eval_bounds_ok)
                  call yaml_mapping_close()
                  if (smatl%nspin>1 .and. jspin==smatl%nspin) then
                      call yaml_mapping_close()
                  end if
              end if

              do i=1,2
                  if (eval_bounds_ok_allspins(i)) then
                      eval_bounds_ok_allspins(i) = eval_bounds_ok(i)
                  end if 
              end do

          end do spin_loop

          !!write(*,*) 'eval_bounds_ok_allspins',eval_bounds_ok_allspins

          if (all(eval_bounds_ok_allspins)) then
              exit bounds_loop
          else
              if (.not.eval_bounds_ok_allspins(1)) then
                  ! lower bound not ok
                  !!call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",ispin)*1.2d0,ispin)
                  !!eval_multiplicator = 2.0d0
                  call foe_data_set_real(foe_obj,"evlow",foe_data_get_real(foe_obj,"evlow",1)*bounds_factor_low,1)
                  if (do_scaling) then
                      ! Scale by a smaller amount if there are oscillations
                      if (.not.scale_matrix(2)) then
                          tt = scaling_factor_low
                      else
                          tt = 0.5_mp*(1.0_mp+scaling_factor_low)
                      end if
                      eval_multiplicator = tt
                      scale_matrix(1) = .true.
                      scale_matrix(2) = .false.
                  end if
              else if (.not.eval_bounds_ok_allspins(2)) then
                  ! upper bound not ok
                  !!call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",ispin)*1.2d0,ispin)
                  !!eval_multiplicator = 1.d0/2.0d0
                  call foe_data_set_real(foe_obj,"evhigh",foe_data_get_real(foe_obj,"evhigh",1)*bounds_factor_up,1)
                  if (do_scaling) then
                      ! Scale by a smaller amount if there are oscillations
                      if (.not.scale_matrix(1)) then
                          tt = scaling_factor_up
                      else
                          tt = 0.5_mp*(1.0_mp+scaling_factor_up)
                      end if
                      eval_multiplicator = tt
                      scale_matrix(2) = .true.
                      scale_matrix(1) = .false.
                  end if
              end if
          end if
          call f_free_ptr(cc_)
          call f_free_ptr(chebyshev_polynomials)
      end do bounds_loop
      if (iproc==0 .and. foe_verbosity>0) then
          call yaml_sequence_close()
      end if

      if (iproc==0) then
          call yaml_mapping_open('summary',flow=.true.)
          call yaml_map('npl',npl)
          if (increase_degree_for_penaltyfunction) then
              call yaml_map('npl determined by','penalty')
          else
              call yaml_map('npl determined by','function')
          end if
          if (do_scaling) call yaml_map('scale',eval_multiplicator_total,fmt='(es9.2)')
          call yaml_newline()
          call yaml_map('bounds', &
               (/foe_data_get_real(foe_obj,"evlow",1),foe_data_get_real(foe_obj,"evhigh",1)/),fmt='(f7.3)')
          call yaml_map('exp accur',max_error,fmt='(es8.2)')
          call yaml_mapping_close()
      end if


      if (do_scaling) then
          call deallocate_matrices(ham_scaled)
      end if
      if (present(cc)) then
          !f_malloc((/npl,ncalc,3/),id='cc')
           cc = f_malloc_ptr((/npl,ncalc,3/), id='cc')
           call f_memcpy(src=cc_, dest=cc)
      end if
      call f_free_ptr(cc_)
      if (present(max_errorx)) then
          call f_memcpy(src=max_error, dest=max_errorx)
      end if
      call f_free(max_error)
      call f_free(x_max_error)
      call f_free(mean_error)

      call f_release_routine()

    end subroutine get_bounds_and_polynomials


    subroutine prepare_matrix(smat, invovrlp_compr, matrix)
      use sparsematrix_init, only: matrixindex_in_compressed
      use dynamic_memory
      implicit none

      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      real(kind=mp),dimension(smat%nvctrp_tg),intent(in) :: invovrlp_compr
      real(kind=mp),dimension(smat%smmm%nvctrp),intent(inout) :: matrix

      ! Local variables
      integer :: i, ii, iline, icolumn, jj

      call f_routine(id='prepare_matrix')

      if (.not.smat%smatmul_initialized) then
          call f_err_throw('sparse matrix multiplication not initialized', &
               err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if

      !$omp parallel &
      !$omp default(none) &
      !$omp shared(smat, matrix, invovrlp_compr) &
      !$omp private(i, ii, iline, icolumn, jj)
      !$omp do schedule(guided)
      do i=1,smat%smmm%nvctrp
          ii = smat%smmm%isvctr + i
          iline = smat%smmm%line_and_column(1,i)
          icolumn = smat%smmm%line_and_column(2,i)
          jj=matrixindex_in_compressed(smat, icolumn, iline)
          if (jj>0) then
              matrix(i) = invovrlp_compr(jj-smat%isvctrp_tg)
          else
              matrix(i) = 0.d0
          end if
      end do
      !$omp end do
      !$omp end parallel

      call f_release_routine()

    end subroutine prepare_matrix



end module foe_common
