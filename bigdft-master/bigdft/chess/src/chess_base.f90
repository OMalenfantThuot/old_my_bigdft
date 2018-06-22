!> @file
!!   File containing the main FOE routine
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

module chess_base
  use f_precisions
  use dynamic_memory
  use f_input_file, only: ATTRS
  use sparsematrix_base
  implicit none

  private

  character(len=*), parameter :: DICT_COMPLETED          = '__dict_has_been_checked__'//ATTRS

  !> Public routines
  public :: chess_input_dict
  public :: chess_init

  !> Private types
  type, public :: foe_params
    real(mp) :: ef_interpol_det
    real(mp) :: ef_interpol_chargediff
    integer :: evbounds_nsatur
    integer :: evboundsshrink_nsatur
    integer :: occupation_function
    real(mp) :: fscale
    real(mp) :: fscale_lowerbound
    real(mp) :: fscale_upperbound
    real(mp),dimension(2) :: eval_range_foe
    real(mp) :: accuracy_foe, accuracy_ice, accuracy_penalty, accuracy_entropy
    real(mp) :: betax_foe, betax_ice
    logical :: adjust_fscale
    logical :: matmul_optimize_load_balancing
    real(mp) :: fscale_ediff_low
    real(mp) :: fscale_ediff_up
  end type foe_params

  type, public :: lapack_params
    integer :: blocksize_pdsyev
    integer :: blocksize_pdgemm
    integer :: maxproc_pdsyev
    integer :: maxproc_pdgemm
  end type lapack_params

  type, public :: pexsi_params
    integer :: pexsi_npoles
    integer :: pexsi_nproc_per_pole
    real(mp) :: pexsi_mumin
    real(mp) :: pexsi_mumax
    real(mp) :: pexsi_mu
    real(mp) :: pexsi_temperature
    real(mp) :: pexsi_tol_charge
    integer :: pexsi_np_sym_fact
    real(mp) :: pexsi_DeltaE
    logical :: pexsi_do_inertia_count
    integer :: pexsi_max_iter
    integer :: pexsi_verbosity
  end type pexsi_params

  !> Public types
  type,public :: chess_params
    type(foe_params) :: foe
    type(lapack_params) :: lapack
    type(pexsi_params) :: pexsi
  end type chess_params


  !> Parameters
  character(len=*),parameter :: FOE_PARAMETERS         = "foe"
  character(len=*),parameter :: LAPACK_PARAMETERS      = "lapack"
  character(len=*),parameter :: PEXSI_PARAMETERS       = "pexsi"

  character(len=*),parameter :: EF_INTERPOL_DET        = "ef_interpol_det"
  character(len=*),parameter :: EF_INTERPOL_CHARGEDIFF = "ef_interpol_chargediff"
  character(len=*),parameter :: EVBOUNDS_NSATUR        = "evbounds_nsatur"
  character(len=*),parameter :: EVBOUNDSSHRINK_NSATUR  = "evboundsshrink_nsatur"
  character(len=*),parameter :: FSCALE                 = "fscale"
  character(len=*),parameter :: FSCALE_LOWERBOUND      = "fscale_lowerbound"
  character(len=*),parameter :: FSCALE_UPPERBOUND      = "fscale_upperbound"
  character(len=*),parameter :: EVAL_RANGE_FOE         = "eval_range_foe"
  character(len=*),parameter :: BLOCKSIZE_PDSYEV       = "blocksize_pdsyev"
  character(len=*),parameter :: BLOCKSIZE_PDGEMM       = "blocksize_pdgemm"
  character(len=*),parameter :: MAXPROC_PDSYEV         = "maxproc_pdsyev"
  character(len=*),parameter :: MAXPROC_PDGEMM         = "maxproc_pdgemm"
  character(len=*),parameter :: PEXSI_NPOLES           = "pexsi_npoles"
  character(len=*),parameter :: PEXSI_NPROC_PER_POLE   = "pexsi_nproc_per_pole"
  character(len=*),parameter :: PEXSI_MUMIN            = "pexsi_mumin"
  character(len=*),parameter :: PEXSI_MUMAX            = "pexsi_mumax"
  character(len=*),parameter :: PEXSI_MU               = "pexsi_mu"
  character(len=*),parameter :: PEXSI_TEMPERATURE      = "pexsi_temperature"
  character(len=*),parameter :: PEXSI_TOL_CHARGE       = "pexsi_tol_charge"
  character(len=*),parameter :: PEXSI_NP_SYM_FACT      = "pexsi_np_sym_fact"
  character(len=*),parameter :: PEXSI_DELTAE           = "pexsi_DeltaE"
  character(len=*),parameter :: PEXSI_MAX_ITER         = "pexsi_max_iter"
  character(len=*),parameter :: PEXSI_VERBOSITY        = "pexsi_verbosity"
  character(len=*),parameter :: PEXSI_DO_INERTIA_COUNT = "pexsi_do_inertia_count"
  character(len=*),parameter :: ACCURACY_FOE           = "accuracy_foe"
  character(len=*),parameter :: ACCURACY_ICE           = "accuracy_ice"
  character(len=*),parameter :: ACCURACY_PENALTY       = "accuracy_penalty"
  character(len=*),parameter :: ACCURACY_ENTROPY       = "accuracy_entropy"
  character(len=*),parameter :: BETAX_FOE              = "betax_foe"
  character(len=*),parameter :: BETAX_ICE              = "betax_ice"
  character(len=*),parameter :: OCCUPATION_FUNCTION    = "occupation_function"
  character(len=*),parameter :: ADJUST_FSCALE          = "adjust_fscale"
  character(len=*),parameter :: MATMUL_OPTIMIZE_LOAD_BALANCING = "matmul_optimize_load_balancing"
  character(len=*),parameter :: FSCALE_EDIFF_LOW       = "fscale_ediff_low"
  character(len=*),parameter :: FSCALE_EDIFF_UP        = "fscale_ediff_up"



  contains

    pure function chess_params_null() result(cp)
      implicit none
      type(chess_params) :: cp
      cp%foe = foe_params_null()
      cp%lapack = lapack_params_null()
      cp%pexsi = pexsi_params_null()
    end function chess_params_null

    pure function foe_params_null() result(fp)
      implicit none
      type(foe_params) :: fp
      fp%ef_interpol_det = 0.0_mp
      fp%ef_interpol_chargediff = 0.0_mp
      fp%evbounds_nsatur = 0
      fp%evboundsshrink_nsatur = 0
      fp%fscale = 0.0_mp
      fp%fscale_lowerbound = 0.0_mp
      fp%fscale_upperbound = 0.0_mp
      fp%eval_range_foe(1:2) = 0.0_mp
      fp%accuracy_foe = 0.0_mp
      fp%accuracy_ice = 0.0_mp
      fp%accuracy_penalty = 0.0_mp
      fp%accuracy_entropy = 0.0_mp
      fp%betax_foe = 0.0_mp
      fp%betax_ice = 0.0_mp
      fp%occupation_function = 0
      fp%adjust_fscale = .false.
      fp%matmul_optimize_load_balancing = .false.
      fp%fscale_ediff_low = 0.0_mp
      fp%fscale_ediff_up = 0.0_mp
    end function foe_params_null

    pure function lapack_params_null() result(lp)
      implicit none
      type(lapack_params) :: lp
      lp%blocksize_pdsyev = 0
      lp%blocksize_pdgemm = 0
      lp%maxproc_pdsyev = 0
      lp%maxproc_pdgemm = 0
    end function lapack_params_null

    pure function pexsi_params_null() result(pp)
      implicit none
      type(pexsi_params) :: pp
      pp%pexsi_npoles = 0
      pp%pexsi_nproc_per_pole = 0
      pp%pexsi_mumin = 0.0_mp
      pp%pexsi_mumax = 0.0_mp
      pp%pexsi_mu = 0.0_mp
      pp%pexsi_temperature = 0.0_mp
      pp%pexsi_tol_charge = 0.0_mp
      pp%pexsi_np_sym_fact = 0
      pp%pexsi_DeltaE = 0.0_mp
      pp%pexsi_do_inertia_count = .false.
      pp%pexsi_max_iter = 0
      pp%pexsi_verbosity = 0
    end function pexsi_params_null



    subroutine chess_init(dict, cp)
      use dictionaries
      implicit none
      type(dictionary),pointer,intent(in) :: dict
      type(chess_params),intent(inout) :: cp
      cp = chess_params_null()
      call chess_fill_variables(dict, cp)
    end subroutine chess_init



    !>routine to fill the CheSS input variables
    subroutine chess_input_dict(dict,dict_minimal)
      use dictionaries
      use f_input_file
      use yaml_parse

      use yaml_output
      implicit none
      !>input dictionary, a copy of the user input, to be filled
      !!with all the variables on exit
      type(dictionary), pointer :: dict
      type(dictionary), pointer, optional :: dict_minimal
      !local variables
      integer(f_integer) :: params_size
      !integer(kind = 8) :: cbuf_add !< address of c buffer
      character, dimension(:), allocatable :: params
      type(dictionary), pointer :: parameters
      type(dictionary), pointer :: parsed_parameters
      type(dictionary), pointer :: profiles
      type(dictionary), pointer :: nested,asis
  
      call f_routine(id='chess_input_dict')
  
      nullify(parameters,parsed_parameters,profiles)
  
      !alternative filling of parameters from hard-coded source file
      !call getstaticinputdef(cbuf_add,params_size)
      call getchessinputdefsize(params_size)
      !allocate array
      params=f_malloc_str(1,params_size,id='params')
      !fill it and parse dictionary
      call getchessinputdef(params)
  
      call yaml_parse_from_char_array(parsed_parameters,params)
      !there is only one document in the input variables specifications
      parameters=>parsed_parameters//0
      profiles => parsed_parameters//1
      call f_free_str(1,params)
  
      call input_file_complete(parameters,dict)!,imports=profiles)
  
      if (present(dict_minimal)) then
         nullify(nested,asis)
         call input_file_minimal(parameters,dict,dict_minimal,nested,asis)
      end if
  
      if (associated(parsed_parameters)) then
         call dict_free(parsed_parameters)
         nullify(parameters)
         nullify(profiles)
      else
         call dict_free(parameters)
      end if
  
      !write in the dictionary that it has been completed
      call set(dict//DICT_COMPLETED,.true.)
  
      call f_release_routine()
  
    end subroutine chess_input_dict



    subroutine chess_fill_variables(dict, cp)
      use dictionaries
      implicit none
      type(dictionary),pointer,intent(in) :: dict
      type(chess_params),intent(inout) :: cp
      !local variables
      type(dictionary), pointer :: lvl,var

      call f_routine(id='chess_fill_variables')

      ! Transfer dict values into input_variables structure.
      lvl => dict_iter(dict)
      do while(associated(lvl))
         var => dict_iter(lvl)
         do while(associated(var))
            call chess_input_fill(var,dict_key(lvl),cp)
            var => dict_next(var)
         end do
         lvl => dict_next(lvl)
      end do

      call f_release_routine()

    end subroutine chess_fill_variables



    !> Set the dictionary from the input variables
    subroutine chess_input_fill(val, level, cp)
      use dictionaries
      use yaml_output
      implicit none
      type(dictionary),pointer,intent(in) :: val
      character(len = *),intent(in) :: level
      type(chess_params),intent(inout) :: cp

      if (index(dict_key(val), "_attributes") > 0) return

      select case(trim(level))
      case(FOE_PARAMETERS)
          select case (trim(dict_key(val)))
          case(EF_INTERPOL_DET)
              cp%foe%ef_interpol_det = val
          case(EF_INTERPOL_CHARGEDIFF)
              cp%foe%ef_interpol_chargediff = val
          case(EVBOUNDS_NSATUR)
              cp%foe%evbounds_nsatur = val
          case(EVBOUNDSSHRINK_NSATUR)
              cp%foe%evboundsshrink_nsatur = val
          case(FSCALE)
              cp%foe%fscale = val
          case(FSCALE_LOWERBOUND)
              cp%foe%fscale_lowerbound = val
          case(FSCALE_UPPERBOUND)
              cp%foe%fscale_upperbound = val
          case(EVAL_RANGE_FOE)
              cp%foe%eval_range_foe = val
          case(ACCURACY_FOE)
              cp%foe%accuracy_foe = val
          case(ACCURACY_ICE)
              cp%foe%accuracy_ice = val
          case(ACCURACY_PENALTY)
              cp%foe%accuracy_penalty = val
          case(ACCURACY_ENTROPY)
              cp%foe%accuracy_entropy = val
          case(BETAX_FOE)
              cp%foe%betax_foe = val
          case(BETAX_ICE)
              cp%foe%betax_ice = val
          case(OCCUPATION_FUNCTION)
              cp%foe%occupation_function = val
          case(ADJUST_FSCALE)
              cp%foe%adjust_fscale = val
          case(MATMUL_OPTIMIZE_LOAD_BALANCING)
              cp%foe%matmul_optimize_load_balancing = val
          case(FSCALE_EDIFF_LOW)
              cp%foe%fscale_ediff_low = val
          case(FSCALE_EDIFF_UP)
              cp%foe%fscale_ediff_up = val
          case default
              call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
          end select
      case(LAPACK_PARAMETERS)
          select case (trim(dict_key(val)))
          case(BLOCKSIZE_PDSYEV)
              cp%lapack%blocksize_pdsyev = val
          case(BLOCKSIZE_PDGEMM)
              cp%lapack%blocksize_pdgemm = val
          case(MAXPROC_PDSYEV)
              cp%lapack%maxproc_pdsyev = val
          case(MAXPROC_PDGEMM)
              cp%lapack%maxproc_pdgemm = val
          case default
              call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
          end select
      case(PEXSI_PARAMETERS)
          select case (trim(dict_key(val)))
          case(PEXSI_NPOLES)
              cp%pexsi%pexsi_npoles = val
          case(PEXSI_NPROC_PER_POLE)
              cp%pexsi%pexsi_nproc_per_pole = val
          case(PEXSI_MUMIN)
              cp%pexsi%pexsi_mumin = val
          case(PEXSI_MUMAX)
              cp%pexsi%pexsi_mumax = val
          case(PEXSI_MU)
              cp%pexsi%pexsi_mu = val
          case(PEXSI_TEMPERATURE)
              cp%pexsi%pexsi_temperature = val
          case(PEXSI_TOL_CHARGE)
              cp%pexsi%pexsi_tol_charge = val
          case(PEXSI_NP_SYM_FACT)
              cp%pexsi%pexsi_np_sym_fact = val
          case(PEXSI_DELTAE)
              cp%pexsi%pexsi_DeltaE = val
          case(PEXSI_DO_INERTIA_COUNT)
              cp%pexsi%pexsi_do_inertia_count = val
          case(PEXSI_MAX_ITER)
              cp%pexsi%pexsi_max_iter = val
          case(PEXSI_VERBOSITY)
              cp%pexsi%pexsi_verbosity = val
          case default
              call yaml_warning("unknown input key '" // trim(level) // "/" // trim(dict_key(val)) // "'")
          end select
      case default
          call yaml_warning("unknown input key '" // trim(level) // "'")
      end select

    end subroutine chess_input_fill


end module chess_base
