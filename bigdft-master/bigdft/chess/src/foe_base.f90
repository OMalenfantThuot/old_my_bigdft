!> @file
!!   File containing basic routines for FOE
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


module foe_base
  use sparsematrix_base
  use dictionaries, only: f_err_throw
  implicit none

  private

  type,public :: foe_data
    real(kind=mp) :: ef     !< Fermi energy for FOE
    real(kind=mp),dimension(:),pointer :: evlow, evhigh !< Eigenvalue bounds for FOE (up/down spin)
    real(kind=mp),dimension(:),pointer :: bisection_shift !< Bisection shift to find Fermi energy (FOE) (up/down spin)
    real(kind=mp) :: fscale                      !< Length scale for complementary error function (FOE)
    real(kind=mp) :: ef_interpol_det             !< FOE: max determinant of cubic interpolation matrix
    real(kind=mp) :: ef_interpol_chargediff      !< FOE: max charge difference for interpolation
    real(kind=mp),dimension(:),pointer :: charge !< Total charge of the system (up/down spin)
    real(kind=mp) :: fscale_lowerbound           !< lower bound for the error function decay length
    real(kind=mp) :: fscale_upperbound           !< upper bound for the error function decay length
    real(kind=mp) :: tmprtr                      !< temperature
    integer :: evbounds_isatur, evboundsshrink_isatur, evbounds_nsatur, evboundsshrink_nsatur !< variables to check whether the eigenvalue bounds might be too big
    real(kind=mp) :: evlow_min, evhigh_max
    real(kind=mp),dimension(:),pointer :: eval_multiplicator !< multiplicative factor to scale the eigenvalue spectrum
    integer :: npl_min !< minimal polynomial degree
    integer :: npl_max !< maximal polynomial degree
    integer :: npl_stride !< stride to increase the polynomial order when searching for the degree yielding the desired precision
    real(kind=mp) :: betax !< exponent to be used for the exponential penalty function
    integer :: ntemp !< Number of FOE iterations with adapted temperature
    real(mp) :: accuracy_function !< Accuracy of the Chebyshev fit for the function to be approximated
    real(mp) :: accuracy_penalty !< Accuracy of the Chebyshev fit for the penalty function to estimate the eigenvalue bounds
    integer :: occupation_function !< Function to determine the occupation numbers
    logical :: adjust_fscale !< dynamically adjust fscale or not
    real(kind=mp) :: fscale_ediff_low !< lower bound for the optimal relative energy difference between the kernel and the control kernel
    real(kind=mp) :: fscale_ediff_up  !< upper bound for the optimal relative energy difference between the kernel and the control kernel
  end type foe_data


  public :: foe_data_null, nullify_foe_data
  public :: foe_data_deallocate
  public :: foe_data_set_int
  public :: foe_data_get_int
  public :: foe_data_set_logical
  public :: foe_data_set_real
  public :: foe_data_get_real
  public :: foe_data_get_logical
  public :: copy_foe_data


  contains

    pure subroutine nullify_foe_data(foe_obj)
      use f_utils, only: f_none,assignment(=)
      implicit none
      type(foe_data), intent(out) :: foe_obj
      nullify(foe_obj%evlow)
      nullify(foe_obj%evhigh)
      nullify(foe_obj%bisection_shift)
      nullify(foe_obj%charge)        
      nullify(foe_obj%eval_multiplicator)
      foe_obj%ef                     =f_none()
      foe_obj%fscale                 =f_none()! uninitialized(foe_obj%fscale)
      foe_obj%ef_interpol_det        =f_none()! uninitialized(foe_obj%ef_interpol_det)
      foe_obj%ef_interpol_chargediff =f_none()! uninitialized(foe_obj%ef_interpol_chargediff)
      foe_obj%fscale_lowerbound      =f_none()! uninitialized(foe_obj%fscale_lowerbound)
      foe_obj%fscale_upperbound      =f_none()! uninitialized(foe_obj%fscale_upperbound)
      foe_obj%tmprtr                 =f_none()! uninitialized(foe_obj%tmprtr)
      foe_obj%evbounds_isatur        =f_none()! uninitialized(foe_obj%evbounds_isatur)
      foe_obj%evboundsshrink_isatur  =f_none()! uninitialized(foe_obj%evboundsshrink_isatur)
      foe_obj%evbounds_nsatur        =f_none()! uninitialized(foe_obj%evbounds_nsatur)
      foe_obj%evboundsshrink_nsatur  =f_none()! uninitialized(foe_obj%evboundsshrink_nsatur)
      foe_obj%evlow_min              =f_none()! uninitialized(foe_obj%evlow_min)
      foe_obj%evhigh_max             =f_none()! uninitialized(foe_obj%evhigh_max)
      foe_obj%npl_min                =f_none()
      foe_obj%npl_max                =f_none()
      foe_obj%npl_stride             =f_none()
      foe_obj%betax                  =f_none()
      foe_obj%ntemp                  =f_none()
      foe_obj%accuracy_function      =f_none()
      foe_obj%accuracy_penalty       =f_none()
      foe_obj%occupation_function    =f_none()
      foe_obj%adjust_fscale          =f_none()
      foe_obj%fscale_ediff_low       =f_none()
      foe_obj%fscale_ediff_up        =f_none()
    end subroutine nullify_foe_data

    pure function foe_data_null() result(foe_obj)
      implicit none
      type(foe_data) :: foe_obj
      call nullify_foe_data(foe_obj)
    end function foe_data_null


    subroutine copy_foe_data(foe_obj_in, foe_obj_out)
      use dynamic_memory
      implicit none
      type(foe_data),intent(in) :: foe_obj_in
      type(foe_data),intent(out) :: foe_obj_out
      foe_obj_out%evlow              = &
          f_malloc_ptr(lbound(foe_obj_in%evlow).to.ubound(foe_obj_in%evlow))
      foe_obj_out%evhigh             = &
          f_malloc_ptr(lbound(foe_obj_in%evhigh).to.ubound(foe_obj_in%evhigh))
      foe_obj_out%bisection_shift    = &
          f_malloc_ptr(lbound(foe_obj_in%bisection_shift).to.ubound(foe_obj_in%bisection_shift))
      foe_obj_out%charge             = &
          f_malloc_ptr(lbound(foe_obj_in%charge).to.ubound(foe_obj_in%charge))
      foe_obj_out%eval_multiplicator = &
          f_malloc_ptr(lbound(foe_obj_in%eval_multiplicator).to.ubound(foe_obj_in%eval_multiplicator))
      call f_memcpy(src=foe_obj_in%evlow,              dest=foe_obj_out%evlow)
      call f_memcpy(src=foe_obj_in%evhigh,             dest=foe_obj_out%evhigh)
      call f_memcpy(src=foe_obj_in%bisection_shift,    dest=foe_obj_out%bisection_shift)
      call f_memcpy(src=foe_obj_in%charge,             dest=foe_obj_out%charge)
      call f_memcpy(src=foe_obj_in%eval_multiplicator, dest=foe_obj_out%eval_multiplicator)
      foe_obj_out%ef                     = foe_obj_in%ef
      foe_obj_out%fscale                 = foe_obj_in%fscale
      foe_obj_out%ef_interpol_det        = foe_obj_in%ef_interpol_det
      foe_obj_out%ef_interpol_chargediff = foe_obj_in%ef_interpol_chargediff
      foe_obj_out%fscale_lowerbound      = foe_obj_in%fscale_lowerbound
      foe_obj_out%fscale_upperbound      = foe_obj_in%fscale_upperbound
      foe_obj_out%tmprtr                 = foe_obj_in%tmprtr
      foe_obj_out%evbounds_isatur        = foe_obj_in%evbounds_isatur
      foe_obj_out%evboundsshrink_isatur  = foe_obj_in%evboundsshrink_isatur
      foe_obj_out%evbounds_nsatur        = foe_obj_in%evbounds_nsatur
      foe_obj_out%evboundsshrink_nsatur  = foe_obj_in%evboundsshrink_nsatur
      foe_obj_out%evlow_min              = foe_obj_in%evlow_min
      foe_obj_out%evhigh_max             = foe_obj_in%evhigh_max
      foe_obj_out%npl_min                = foe_obj_in%npl_min
      foe_obj_out%npl_max                = foe_obj_in%npl_max
      foe_obj_out%npl_stride             = foe_obj_in%npl_stride
      foe_obj_out%betax                  = foe_obj_in%betax
      foe_obj_out%ntemp                  = foe_obj_in%ntemp
      foe_obj_out%accuracy_function      = foe_obj_in%accuracy_function
      foe_obj_out%accuracy_penalty       = foe_obj_in%accuracy_penalty
      foe_obj_out%occupation_function    = foe_obj_in%occupation_function
      foe_obj_out%adjust_fscale          = foe_obj_in%adjust_fscale
      foe_obj_out%fscale_ediff_low       = foe_obj_in%fscale_ediff_low
      foe_obj_out%fscale_ediff_up        = foe_obj_in%fscale_ediff_up
    end subroutine copy_foe_data


    subroutine foe_data_deallocate(foe_obj)
      use dynamic_memory
      implicit none
      type(foe_data) :: foe_obj
      call f_free_ptr(foe_obj%evlow)
      call f_free_ptr(foe_obj%evhigh)
      call f_free_ptr(foe_obj%bisection_shift)
      call f_free_ptr(foe_obj%charge)
      call f_free_ptr(foe_obj%eval_multiplicator)
    end subroutine foe_data_deallocate


    subroutine foe_data_set_int(foe_obj, fieldname, val)
      type(foe_data) :: foe_obj
      character(len=*),intent(in) :: fieldname
      integer,intent(in) :: val

      select case (fieldname)
      case ("nseg")
          !!foe_obj%nseg = val
      case ("evbounds_isatur")
          foe_obj%evbounds_isatur = val
      case ("evboundsshrink_isatur")
          foe_obj%evboundsshrink_isatur = val
      case ("evbounds_nsatur")
          foe_obj%evbounds_nsatur = val
      case ("evboundsshrink_nsatur")
          foe_obj%evboundsshrink_nsatur = val
      case ("npl_min")
          foe_obj%npl_min = val
      case ("npl_max")
          foe_obj%npl_max = val
      case ("npl_stride")
          foe_obj%npl_stride = val
      case ("ntemp")
          foe_obj%ntemp = val
      case ("occupation_function")
          foe_obj%occupation_function = val
      case default
          call f_err_throw("wrong argument for "//trim(fieldname))
      end select

    end subroutine foe_data_set_int


    integer function foe_data_get_int(foe_obj, fieldname) result(val)
      type(foe_data) :: foe_obj
      character(len=*),intent(in) :: fieldname
      !integer,intent(in) :: val

      select case (fieldname)
      case ("nseg")
          !!val = foe_obj%nseg
      case ("evbounds_isatur")
          val = foe_obj%evbounds_isatur
      case ("evboundsshrink_isatur")
          val = foe_obj%evboundsshrink_isatur
      case ("evbounds_nsatur")
          val = foe_obj%evbounds_nsatur
      case ("evboundsshrink_nsatur")
          val = foe_obj%evboundsshrink_nsatur
      case ("npl_min")
          val = foe_obj%npl_min
      case ("npl_max")
          val = foe_obj%npl_max
      case ("npl_stride")
          val = foe_obj%npl_stride
      case ("ntemp")
          val = foe_obj%ntemp
      case ("occupation_function")
          val = foe_obj%occupation_function
      case default
          call f_err_throw("wrong argument for "//trim(fieldname))
      end select

    end function foe_data_get_int


    subroutine foe_data_set_real(foe_obj, fieldname, val, ind)
      type(foe_data) :: foe_obj
      character(len=*),intent(in) :: fieldname
      real(kind=mp),intent(in) :: val
      integer,intent(in),optional :: ind

      select case (fieldname)
      case ("ef")
          foe_obj%ef = val
      case ("evlow")
          if (.not.present(ind)) then
              call f_err_throw('ind not present')
          end if
          if (.not.associated(foe_obj%evlow)) then
              call f_err_throw('evlow not associated')
          end if
          foe_obj%evlow(ind) = val
      case ("evhigh")
          if (.not.present(ind)) then
              call f_err_throw('ind not present')
          end if
          if (.not.associated(foe_obj%evhigh)) then
              call f_err_throw('evhigh not associated')
          end if
          foe_obj%evhigh(ind) = val
      case ("bisection_shift")
          if (.not.present(ind)) then
              call f_err_throw('ind not present')
          end if
          if (.not.associated(foe_obj%bisection_shift)) then
              call f_err_throw('bisection_shift not associated')
          end if
          foe_obj%bisection_shift(ind) = val
      case ("fscale")
          foe_obj%fscale = val
      case ("ef_interpol_det")
          foe_obj%ef_interpol_det = val
      case ("ef_interpol_chargediff")
          foe_obj%ef_interpol_chargediff = val
      case ("charge")
          if (.not.present(ind)) then
              call f_err_throw('ind not present')
          end if
          if (.not.associated(foe_obj%charge)) then
              call f_err_throw('charge not associated')
          end if
          foe_obj%charge(ind) = val
      case ("fscale_lowerbound")
          foe_obj%fscale_lowerbound = val
      case ("fscale_upperbound")
          foe_obj%fscale_upperbound = val
      case ("tmprtr")
          foe_obj%tmprtr = val
      case ("evlow_min")
          foe_obj%evlow_min = val
      case ("evhigh_max")
          foe_obj%evhigh_max = val
      case ("eval_multiplicator")
          if (.not.present(ind)) then
              call f_err_throw('ind not present')
          end if
          if (.not.associated(foe_obj%eval_multiplicator)) then
              call f_err_throw('eval_multiplicator not associated')
          end if
          foe_obj%eval_multiplicator(ind) = val
      case ("betax")
          foe_obj%betax = val
      case ("accuracy_function")
          foe_obj%accuracy_function = val
      case ("accuracy_penalty")
          foe_obj%accuracy_penalty = val
      case ("fscale_ediff_low")
          foe_obj%fscale_ediff_low = val
      case ("fscale_ediff_up")
          foe_obj%fscale_ediff_up = val
      case default
          call f_err_throw("wrong argument for "//trim(fieldname))
      end select

    end subroutine foe_data_set_real


    real(kind=mp) function foe_data_get_real(foe_obj, fieldname, ind) result(val)
      type(foe_data) :: foe_obj
      character(len=*),intent(in) :: fieldname
      integer,intent(in),optional :: ind

      select case (fieldname)
      case ("ef")
          val = foe_obj%ef
      case ("evlow")
          if (.not.present(ind)) then
              call f_err_throw('ind not present')
          end if
          if (.not.associated(foe_obj%evlow)) then
              call f_err_throw('evlow not associated')
          end if
          val = foe_obj%evlow(ind)
      case ("evhigh")
          if (.not.present(ind)) then
              call f_err_throw('ind not present')
          end if
          if (.not.associated(foe_obj%evhigh)) then
              call f_err_throw('evhigh not associated')
          end if
          val = foe_obj%evhigh(ind)
      case ("bisection_shift")
          if (.not.present(ind)) then
              call f_err_throw('ind not present')
          end if
          if (.not.associated(foe_obj%bisection_shift)) then
              call f_err_throw('bisection_shift not associated')
          end if
          val = foe_obj%bisection_shift(ind)
      case ("fscale")
          val = foe_obj%fscale
      case ("ef_interpol_det")
          val = foe_obj%ef_interpol_det
      case ("ef_interpol_chargediff")
          val = foe_obj%ef_interpol_chargediff
      case ("charge")
          if (.not.present(ind)) then
              call f_err_throw('ind not present')
          end if
          if (.not.associated(foe_obj%charge)) then
              call f_err_throw('charge not associated')
          end if
          val = foe_obj%charge(ind)
      case ("fscale_lowerbound")
          val = foe_obj%fscale_lowerbound
      case ("fscale_upperbound")
          val = foe_obj%fscale_upperbound
      case ("tmprtr")
          val = foe_obj%tmprtr
      case ("evlow_min")
          val = foe_obj%evlow_min
      case ("evhigh_max")
          val = foe_obj%evhigh_max
      case ("eval_multiplicator")
          if (.not.present(ind)) then
              call f_err_throw('ind not present')
          end if
          if (.not.associated(foe_obj%eval_multiplicator)) then
              call f_err_throw('eval_multiplicator not associated')
          end if
          val = foe_obj%eval_multiplicator(ind)
      case ("betax")
          val = foe_obj%betax
      case ("accuracy_function")
          val = foe_obj%accuracy_function
      case ("accuracy_penalty")
          val = foe_obj%accuracy_penalty
      case ("fscale_ediff_low")
          val = foe_obj%fscale_ediff_low
      case ("fscale_ediff_up")
          val = foe_obj%fscale_ediff_up
      case default
          call f_err_throw("wrong argument for "//trim(fieldname))
      end select

    end function foe_data_get_real


    subroutine foe_data_set_logical(foe_obj, fieldname, val)
      type(foe_data) :: foe_obj
      character(len=*),intent(in) :: fieldname
      logical,intent(in) :: val

      select case (fieldname)
      case ("adjust_fscale")
          foe_obj%adjust_fscale = val
      case default
          call f_err_throw("wrong argument for "//trim(fieldname))
      end select

    end subroutine foe_data_set_logical


    logical function foe_data_get_logical(foe_obj, fieldname) result(val)
      type(foe_data) :: foe_obj
      character(len=*),intent(in) :: fieldname

      select case (fieldname)
      case ("adjust_fscale")
          val = foe_obj%adjust_fscale
      case default
          call f_err_throw("wrong argument for "//trim(fieldname))
      end select

    end function foe_data_get_logical

end module foe_base
