!> @file
!! Routines associated to the generation of data for atoms of the system
!! @author
!!    Copyright (C) 2007-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Handling of input guess creation from basis of atomic orbitals
!! and also pseudopotentials
module module_atoms
  use module_defs, only: dp,gp
  use ao_inguess, only: aoig_data
  use m_pawrad, only: pawrad_type
  use m_pawtab, only: pawtab_type
  use m_pawang, only: pawang_type
  use f_refcnts, only: f_reference_counter,nullify_f_ref,f_ref_new
  use public_keys, only : ASTRUCT_CELL,ASTRUCT_POSITIONS, &
       & ASTRUCT_PROPERTIES,ASTRUCT_UNITS, ASTRUCT_ATT_FROZEN, &
       & ASTRUCT_ATT_IGSPIN, ASTRUCT_ATT_IGCHRG, ASTRUCT_ATT_IXYZ_1, &
       & ASTRUCT_ATT_IXYZ_2, ASTRUCT_ATT_IXYZ_3, &
       & ASTRUCT_ATT_RXYZ_INT_1, ASTRUCT_ATT_RXYZ_INT_2, &
       & ASTRUCT_ATT_RXYZ_INT_3, ASTRUCT_ATT_MODE, ASTRUCT_ATT_CAVRAD
  use dictionaries, only: dictionary
  use f_trees, only: f_tree
  use f_arrays
  implicit none

  private

  !> Quantities used for the symmetry operators. To be used in atomic_structure derived type.
  type, public :: symmetry_data
     integer :: symObj                             !< The symmetry object from ABINIT
     integer :: nSym                               !< Number ofsymmetry (0 if disable)
     character(len=15) :: spaceGroup               !< Space group (disabled if not useful)
     integer, dimension(:,:,:), pointer :: irrzon
     real(dp), dimension(:,:,:), pointer :: phnons
  end type symmetry_data

  !> Stores a list of neighbours.
  type, public :: atomic_neighbours
     integer :: nat
     integer, dimension(:,:), pointer :: keynei
     integer, dimension(:), pointer :: nei

     ! Iterator part (to avoid to create yet another type.
     integer :: iat, ind
  end type atomic_neighbours

  !> Structure of the system. This derived type contains the information about the physical properties
  type, public :: atomic_structure
     character(len=1) :: geocode         !shoud become integer, dimension(3) :: bc
     character(len=5) :: inputfile_format  !< Can be xyz, ascii, int or yaml
     character(len=256) :: source          !< Name of the file or origin
     character(len=20) :: units            !< Can be angstroem or bohr
     character(len=20) :: angle            !< Can be radian or degree
     integer :: nat                        !< Number of atoms
     integer :: ntypes                     !< Number of atomic species in the structure
     real(gp), dimension(3) :: cell_dim    !< Dimensions of the simulation domain (each one periodic or free according to geocode)
     real(gp), dimension(3) :: shift       !< Rigid shift applied to the atomic positions
     !pointers
     real(gp), dimension(:,:), pointer :: rxyz             !< Atomic positions (always in AU, units variable is considered for I/O only)
     real(gp), dimension(:,:), pointer :: rxyz_int         !< Atomic positions in internal coordinates (Z matrix)
     integer, dimension(:,:), pointer :: ixyz_int          !< Neighbor list for internal coordinates (Z matrix)
     character(len=20), dimension(:), pointer :: atomnames !< Atomic species names
     integer, dimension(:), pointer :: iatype              !< Atomic species id
     integer, dimension(:), pointer :: ifrztyp             !< Freeze atoms while updating structure
     integer, dimension(:), pointer :: input_polarization  !< Used in AO generation for WFN input guess
     type(symmetry_data) :: sym                            !< The symmetry operators
     type(f_tree), dimension(:), pointer :: attributes     !< Static attributes per atom
     type(dictionary), pointer :: properties               !< Global properties
  end type atomic_structure

  !> Data containing the information about the atoms in the system
  type, public :: atoms_data
     !>reference counter
     type(f_reference_counter) :: refcnt
     type(atomic_structure) :: astruct                   !< Atomic structure (positions and so on)
     type(aoig_data), dimension(:), pointer :: aoig      !< Contains the information needed for generating the AO inputguess data for each atom
     integer :: natsc                                    !< Number of atoms with semicore occupations at the input guess
     integer, dimension(:), pointer :: nelpsp
     integer, dimension(:), pointer :: npspcode          !< PSP codes (see @link psp_projectors::pspcode_hgh @endlink)
     integer, dimension(:), pointer :: ixcpsp            !< PSP ixc code
     integer, dimension(:), pointer :: nzatom            !< Atomic number
     integer, dimension(:), pointer :: iradii_source     !< Source of the radii coefficients (Hard-Coded, PSP File, ...)
     real(gp), dimension(:,:), pointer :: radii_cf       !< User defined radii_cf, overridden in sysprop.f90
     real(gp), dimension(:), pointer :: amu              !< Amu(ntypes)  Atomic Mass Unit for each type of atoms
     !real(gp), dimension(:,:), pointer :: rloc          !< Localization regions for parameters of linear, to be moved somewhere else
     real(gp), dimension(:,:,:), pointer :: psppar       !< Pseudopotential parameters (HGH SR section)
     logical :: donlcc                                   !< Activate non-linear core correction treatment
     logical :: multipole_preserving                     !< Activate preservation of the multipole moment for the ionic charge
     integer :: mp_isf                                   !< Interpolating scaling function order for the multipole preserving
     integer, dimension(:), pointer :: nlcc_ngv,nlcc_ngc !< Number of valence and core gaussians describing NLCC
     real(gp), dimension(:,:), pointer :: nlccpar        !< Parameters for the non-linear core correction, if present
     logical, dimension(:,:), pointer :: dogamma         !< density matrix is calculated
     type(f_matrix), dimension(:,:,:), pointer :: gamma_targets !<array containing target gammas
     logical :: apply_gamma_target !< apply the gamma target to the NL PSP operator
     type(pawrad_type), dimension(:), pointer :: pawrad  !< PAW radial objects.
     type(pawtab_type), dimension(:), pointer :: pawtab  !< PAW objects for something.
     type(pawang_type) :: pawang                         !< PAW angular mesh definition.
     real(gp), dimension(:), pointer :: epsatm !< PAW pseudoatom energy for each type of atom

     !> for abscalc with pawpatch
     integer, dimension(:), pointer ::  paw_NofL, paw_l, paw_nofchannels
     integer, dimension(:), pointer ::  paw_nofgaussians
     real(gp), dimension(:), pointer :: paw_Greal, paw_Gimag, paw_Gcoeffs
     real(gp), dimension(:), pointer :: paw_H_matrices, paw_S_matrices, paw_Sm1_matrices
     integer :: iat_absorber
  end type atoms_data

  !> Iterator on atoms object
  type, public :: atoms_iterator
     integer :: iat                                 !< Atom number
     integer :: ityp                                !< Atom type
     character(len=20) :: name                      !< Atom name
     logical, dimension(3) :: frz                   !< Array of frozen coordinate of the atom
     real(gp), dimension(3) :: rxyz                 !< Atom positions
     type(atomic_structure), pointer :: astruct_ptr !< Private pointer to the atomic structure from which it derives
     type(dictionary), pointer :: attrs
  end type atoms_iterator


  !Public routines
  public :: atoms_data_null,nullify_atoms_data,deallocate_atoms_data,atoms_fill
  public :: atomic_structure_null,nullify_atomic_structure,deallocate_atomic_structure
  public :: astruct_merge_to_dict, astruct_at_from_dict
  public :: deallocate_symmetry_data,set_symmetry_data
  public :: set_astruct_from_file,astruct_dump_to_file,set_astruct_from_openbabel
  public :: allocate_atoms_data,move_this_coordinate,frozen_itof
  public :: rxyz_inside_box,check_atoms_positions
  public :: atomic_data_set_from_dict,atoms_iter,atoms_iter_next,atoms_iter_ensure_attr
  public :: nullify_atomic_neighbours, deallocate_atomic_neighbours
  public :: astruct_neighbours_iter, astruct_neighbours_next
  ! Dictionary inquire
  public :: astruct_dict_get_types
  ! Types from dictionaries
  public :: astruct_set_from_dict
  public :: astruct_file_merge_to_dict,atoms_file_merge_to_dict
  public :: psp_dict_analyse, nlcc_set_from_dict,atoms_gamma_from_dict
  public :: astruct_constraints,astruct_set,atomic_charge_density


contains


  !> Nullify an iterator over atoms
  pure subroutine nullify_atoms_iterator(it)
    implicit none
    type(atoms_iterator), intent(out) :: it
    it%iat=-1
    it%ityp=-1
    it%name=repeat(' ',len(it%name))
    it%frz=.false.
    it%rxyz=0.0_gp
    nullify(it%astruct_ptr)
    nullify(it%attrs)
  end subroutine nullify_atoms_iterator


  !> Function to nullify an iterator avor atoms
  pure function atoms_iterator_null() result(it)
    implicit none
    type(atoms_iterator) :: it

    call nullify_atoms_iterator(it)
  end function atoms_iterator_null


  !> Create an iterator over atoms (atomic_structure)
  function atoms_iter(astruct) result(it)
    implicit none
    type(atomic_structure), intent(in), target :: astruct
    type(atoms_iterator) :: it

    it=atoms_iterator_null()
    ! Start counting
    it%astruct_ptr => astruct
    it%iat=0
    it%ityp=0
  end function atoms_iter


  !> Increment a valid iterator
  !! the control for validity has to be done outside
  pure subroutine atoms_refresh_iterator(it)
    use yaml_strings, only: f_strcpy
    implicit none
    type(atoms_iterator), intent(inout) :: it
    it%ityp=it%astruct_ptr%iatype(it%iat)
    call f_strcpy(src=it%astruct_ptr%atomnames(it%ityp),dest=it%name)
    it%frz(1)=move_this_coordinate(it%astruct_ptr%ifrztyp(it%iat),1)
    it%frz(2)=move_this_coordinate(it%astruct_ptr%ifrztyp(it%iat),2)
    it%frz(3)=move_this_coordinate(it%astruct_ptr%ifrztyp(it%iat),3)
    it%rxyz=it%astruct_ptr%rxyz(:,it%iat)
    it%attrs=>it%astruct_ptr%attributes(it%iat)%d
  end subroutine atoms_refresh_iterator


  !> Increment, and nullify if ended
  !! if the iterator is nullified, it does nothing
  pure subroutine increment_atoms_iter(it)
    implicit none
    type(atoms_iterator), intent(inout) :: it

    if (associated(it%astruct_ptr)) then
       if (it%iat < it%astruct_ptr%nat) then
          it%iat=it%iat+1
          call atoms_refresh_iterator(it)
       else
          !end iteration, the iterator is destroyed
          call nullify_atoms_iterator(it)
       end if
    end if
  end subroutine increment_atoms_iter

  !> Logical function, returns .true. if the iterator is still valid
  pure function atoms_iter_is_valid(it)
    implicit none
    type(atoms_iterator), intent(in) :: it
    logical :: atoms_iter_is_valid

    atoms_iter_is_valid=associated(it%astruct_ptr)
  end function atoms_iter_is_valid

  !> Logical function for iterating above atoms
  function atoms_iter_next(it)
    implicit none
    type(atoms_iterator), intent(inout) :: it
    logical :: atoms_iter_next

    call increment_atoms_iter(it)
    atoms_iter_next=atoms_iter_is_valid(it)
  end function atoms_iter_next

  subroutine atoms_iter_ensure_attr(it)
    use dictionaries, only: dict_init
    implicit none
    type(atoms_iterator), intent(inout) :: it

    if (.not. associated(it%attrs)) then
       call dict_init(it%astruct_ptr%attributes(it%iat)%d)
       it%attrs => it%astruct_ptr%attributes(it%iat)%d
    end if
  end subroutine atoms_iter_ensure_attr

  !> Creators and destructors
  pure function symmetry_data_null() result(sym)
    implicit none
    type(symmetry_data) :: sym
    call nullify_symmetry_data(sym)
  end function symmetry_data_null

  pure subroutine nullify_symmetry_data(sym)
    use f_utils, only: f_zero
    implicit none
    type(symmetry_data), intent(out) :: sym
    sym%symObj=-1
    sym%nsym=-1
    call f_zero(sym%spaceGroup)
    nullify(sym%irrzon)
    nullify(sym%phnons)
  end subroutine nullify_symmetry_data

  pure subroutine nullify_atomic_neighbours(nei)
    type(atomic_neighbours), intent(out) :: nei
    nei%nat=-1
    nullify(nei%keynei)
    nullify(nei%nei)
  end subroutine nullify_atomic_neighbours

  !> Initialize the structure atomic_neighbours
  pure function atomic_structure_null() result(astruct)
    implicit none
    type(atomic_structure) :: astruct
    call nullify_atomic_structure(astruct)
  end function atomic_structure_null


  !> Initialize the structure atomic_structure
  pure subroutine nullify_atomic_structure(astruct)
    implicit none
    type(atomic_structure), intent(out) :: astruct
    astruct%geocode='X'
    astruct%inputfile_format=repeat(' ',len(astruct%inputfile_format))
    astruct%source=repeat(' ',len(astruct%source))
    astruct%units=repeat(' ',len(astruct%units))
    astruct%angle=repeat(' ',len(astruct%angle))
    astruct%nat=-1
    astruct%ntypes=-1
    astruct%cell_dim=0.0_gp
    astruct%shift=0.0_gp
    nullify(astruct%input_polarization)
    nullify(astruct%ifrztyp)
    nullify(astruct%atomnames)
    nullify(astruct%iatype)
    nullify(astruct%rxyz)
    nullify(astruct%rxyz_int)
    nullify(astruct%ixyz_int)
    call nullify_symmetry_data(astruct%sym)
    nullify(astruct%attributes)
    nullify(astruct%properties)
  end subroutine nullify_atomic_structure


  !> Nullify atoms_data structure (function)
  pure function atoms_data_null() result(at)
    type(atoms_data) :: at
    call nullify_atoms_data(at)
  end function atoms_data_null


  !> Nullify atoms_data structure (routine)
  pure subroutine nullify_atoms_data(at)
    !use m_pawang, only: pawang_nullify
    implicit none
    type(atoms_data), intent(out) :: at
    call nullify_f_ref(at%refcnt)
    call nullify_atomic_structure(at%astruct)
    nullify(at%aoig)
    at%natsc=0
    at%donlcc=.false.
    at%multipole_preserving=.false.
    at%mp_isf=0
    at%iat_absorber=-1
    at%apply_gamma_target=.false.
    !     nullify(at%iasctype)
    nullify(at%nelpsp)
    nullify(at%npspcode)
    nullify(at%ixcpsp)
    nullify(at%nzatom)
    nullify(at%radii_cf)
    nullify(at%iradii_source)
    nullify(at%amu)
    !nullify(at%rloc)
    nullify(at%psppar)
    nullify(at%dogamma)
    nullify(at%gamma_targets)
    nullify(at%nlcc_ngv)
    nullify(at%nlcc_ngc)
    nullify(at%nlccpar)
    nullify(at%paw_NofL)
    nullify(at%paw_l)
    nullify(at%paw_nofchannels)
    nullify(at%paw_nofgaussians)
    nullify(at%paw_Greal)
    nullify(at%paw_Gimag)
    nullify(at%paw_Gcoeffs)
    nullify(at%paw_H_matrices)
    nullify(at%paw_S_matrices)
    nullify(at%paw_Sm1_matrices)
    nullify(at%pawrad)
    nullify(at%pawtab)
    nullify(at%epsatm)
    !call pawang_nullify(at%pawang) !not needed in fact
  end subroutine nullify_atoms_data


  !> Destructor of symmetry data operations
  subroutine deallocate_symmetry_data(sym)
    use dynamic_memory, only: f_free_ptr
    use m_ab6_symmetry, only: symmetry_free
    implicit none
    type(symmetry_data), intent(inout) :: sym
    !character(len = *), intent(in) :: subname
    !      integer :: i_stat, i_all

    if (sym%symObj >= 0) then
       call symmetry_free(sym%symObj)
    end if
    call f_free_ptr(sym%irrzon)
    call f_free_ptr(sym%phnons)
  end subroutine deallocate_symmetry_data

  !> Deallocate the structure atomic_neighbours.
  subroutine deallocate_atomic_neighbours(nei)
    use dynamic_memory, only: f_free_ptr
    implicit none
    type(atomic_neighbours), intent(inout) :: nei
    if (nei%nat > 0) then
       call f_free_ptr(nei%keynei)
       call f_free_ptr(nei%nei)
    end if
  end subroutine deallocate_atomic_neighbours

  !> Deallocate the structure atoms_data.
  subroutine deallocate_atomic_structure(astruct)!,subname)
    use dynamic_memory, only: f_free_ptr,f_free_str_ptr
    use dictionaries, only: dict_free
    implicit none
    !character(len=*), intent(in) :: subname
    type(atomic_structure), intent(inout) :: astruct
    !local variables
    character(len=*), parameter :: subname='deallocate_atomic_structure' !remove
    integer :: iat
    ! Deallocations for the geometry part.
    if (astruct%nat >= 0) then
       call f_free_ptr(astruct%ifrztyp)
       call f_free_ptr(astruct%iatype)
       call f_free_ptr(astruct%input_polarization)
       call f_free_ptr(astruct%rxyz)
       call f_free_ptr(astruct%rxyz_int)
       call f_free_ptr(astruct%ixyz_int)
       if (associated(astruct%attributes)) then
          do iat = 1, astruct%nat
             call dict_free(astruct%attributes(iat)%d)
          end do
          deallocate(astruct%attributes)
          nullify(astruct%attributes)
       end if
    end if
    if (astruct%ntypes >= 0) then
       call f_free_str_ptr(len(astruct%atomnames),astruct%atomnames)
    end if
    ! Free additional stuff.
    call deallocate_symmetry_data(astruct%sym)
    if (associated(astruct%properties)) call dict_free(astruct%properties)
  END SUBROUTINE deallocate_atomic_structure


  !> Deallocate the structure atoms_data.
  subroutine deallocate_atoms_data(atoms)
    use module_base
    use m_pawtab, only: pawtab_free
    use m_pawrad, only: pawrad_free
    use m_pawang, only: pawang_free
    implicit none
    type(atoms_data), intent(inout) :: atoms
    !local variables
    character(len=*), parameter :: subname='dellocate_atoms_data' !remove
    integer :: ityp

    !check if the freeing is possible
    call f_ref_free(atoms%refcnt)

    ! Deallocate atomic structure
    call deallocate_atomic_structure(atoms%astruct)

    ! Deallocations related to pseudos.
    call f_free_ptr(atoms%nzatom)
    call f_free_ptr(atoms%psppar)
    call f_free_ptr(atoms%dogamma)
    call f_array_ptr_free(atoms%gamma_targets)
    atoms%apply_gamma_target=.false.
    call f_free_ptr(atoms%nelpsp)
    call f_free_ptr(atoms%ixcpsp)
    call f_free_ptr(atoms%npspcode)
    call f_free_ptr(atoms%nlcc_ngv)
    call f_free_ptr(atoms%nlcc_ngc)
    call f_free_ptr(atoms%radii_cf)
    call f_free_ptr(atoms%iradii_source)
    call f_free_ptr(atoms%amu)
    if (associated(atoms%aoig)) then
       !no flib calls for derived types for the moment
       deallocate(atoms%aoig)
       nullify(atoms%aoig)
    end if
    call f_free_ptr(atoms%nlccpar)

    !  Free data for pawpatch
    call f_free_ptr(atoms%paw_l)
    call f_free_ptr(atoms%paw_NofL)
    call f_free_ptr(atoms%paw_nofchannels)
    call f_free_ptr(atoms%paw_nofgaussians)
    call f_free_ptr(atoms%paw_Greal)
    call f_free_ptr(atoms%paw_Gimag)
    call f_free_ptr(atoms%paw_Gcoeffs)
    call f_free_ptr(atoms%paw_H_matrices)
    call f_free_ptr(atoms%paw_S_matrices)
    call f_free_ptr(atoms%paw_Sm1_matrices)

    ! Free PAW data.
    if (associated(atoms%pawrad)) then
       do ityp = 1, size(atoms%pawrad)
          call pawrad_free(atoms%pawrad(ityp))
       end do
       deallocate(atoms%pawrad)
    end if
    if (associated(atoms%pawtab)) then
       do ityp = 1, size(atoms%pawtab)
          call pawtab_free(atoms%pawtab(ityp))
       end do
       deallocate(atoms%pawtab)
    end if
    call pawang_free(atoms%pawang)
    if (associated(atoms%epsatm)) then
       call f_free_ptr(atoms%epsatm)
    end if
    END SUBROUTINE deallocate_atoms_data

    !> Start the iterator of an astruct_neighbours structure.
    !  Only one iterator is possible at a time.
    subroutine astruct_neighbours_iter(neighb, iat, n)
      implicit none
      type(atomic_neighbours), intent(inout) :: neighb
      integer, intent(in) :: iat
      integer, intent(out), optional :: n

      neighb%iat = iat
      neighb%ind = 0
      if (present(n)) n = neighb%keynei(1, neighb%iat)
    END SUBROUTINE astruct_neighbours_iter
    !> Return the next neighbour of a given atom, as initialised by
    !  astruct_neighbours_iter(). Return 0 if there is no next neighbours.
    function astruct_neighbours_next(neighb, inei)
      implicit none
      type(atomic_neighbours), intent(inout) :: neighb
      integer, intent(out) :: inei

      logical :: astruct_neighbours_next

      astruct_neighbours_next = .false.
      inei = 0

      neighb%ind = neighb%ind + 1
      if (neighb%ind > neighb%keynei(1, neighb%iat)) return

      inei = neighb%nei(neighb%keynei(2, neighb%iat) + neighb%ind - 1)
      astruct_neighbours_next = .true.
    END FUNCTION astruct_neighbours_next


    !> Determine the gaussian structure related to the given atom
    subroutine atomic_charge_density(g,at,atit)
      use gaussians
      use numerics, only: twopi
      implicit none
      type(atoms_data), intent(in) :: at
      type(atoms_iterator), intent(in) :: atit
      type(gaussian_real_space), intent(out) :: g
      !local variables
      integer :: mp_isf
      real(gp) :: rloc
      real(gp), dimension(1) :: charge
      integer, dimension(1) :: zero1
      integer, dimension(3) :: zeros

      rloc=at%psppar(0,0,atit%ityp)
      charge(1)=real(at%nelpsp(atit%ityp),gp)/(twopi*sqrt(twopi)*rloc**3)
      zeros=0
      zero1=0
      mp_isf=at%mp_isf
      if (.not. at%multipole_preserving) mp_isf=0
      call gaussian_real_space_set(g,rloc,1,charge,zeros,zero1,mp_isf)

    end subroutine atomic_charge_density

!!!    subroutine local_psp_terms(g,at,atit)
!!!      implicit none
!!!      type(gaussian_real_space), intent(out) :: g
!!!      !local variables
!!!      integer :: nloc,iloc
!!!      integer, dimension(4) :: pow
!!!      real(gp), dimension(4) :: factors
!!!
!!!      ! Determine the number of local terms
!!!      nloc=0
!!!      do iloc=1,4
!!!         if (at%psppar(0,iloc,atit%ityp) /= 0.d0) nloc=iloc
!!!      enddo
!!!
!!!      pow=0
!!!      factors=0.0_gp
!!!      sigma=at%psppar(0,0,atit%ityp)
!!!      rlocinvsq=onehalf/sigma**2
!!!      do iloc=1,nloc
!!!         pow(iloc)=2*(iloc-1)
!!!         factor(iloc)=at%psppar(0,iloc,atit%ityp)*rlocinvsq**(2*(iloc-1))
!!!      end do
!!!!!$      arg=r2*rlocinvsq
!!!!!$      tt=at%psppar(0,nloc,atit%ityp)
!!!!!$      do iloc=nloc-1,1,-1
!!!!!$         tt=arg*tt+at%psppar(0,iloc,atit%ityp)
!!!!!$      enddo
!!!      call gaussian_real_space_set(g,sigma,nloc,factor,pow,0,0.0_gp)
!!!
!!!    end subroutine local_psp_terms

    !> determine the atomic dictionary from the key
    subroutine get_atomic_dict(dict_tmp,dict,iat,name,fromname)
      use dictionaries
      use yaml_strings
      implicit none
      integer, intent(in) :: iat
      character(len=*), intent(in) :: name
      type(dictionary), pointer :: dict_tmp,dict
      logical, intent(out), optional :: fromname
      !local variables
      character(len = max_field_length) :: at

      nullify(dict_tmp)

      call f_strcpy(src="Atom "//trim(adjustl(yaml_toa(iat))),dest=at)
      dict_tmp = dict .get. at
      if (present(fromname)) fromname=.not. associated(dict_tmp)
      if (.not. associated(dict_tmp)) then
         dict_tmp = dict .get. name
         if (present(fromname)) fromname=associated(dict_tmp)
      end if

    end subroutine get_atomic_dict

    subroutine atoms_gamma_from_dict(dict,dict_targets,&
         lmax,astruct, dogamma,gamma_targets)
      use module_defs, only: gp
      use ao_inguess, only: shell_toi
      use dictionaries
      use yaml_strings, only: f_strcpy, yaml_toa
      use dynamic_memory
      use f_blas
      implicit none
      integer, intent(in) :: lmax
      type(dictionary), pointer :: dict,dict_targets
      type(atomic_structure), intent(in) :: astruct
      logical, dimension(0:lmax,astruct%nat), intent(out) :: dogamma
      type(f_matrix), dimension(:,:,:), pointer :: gamma_targets
      !local variables
      logical :: sometarget,fromname
      integer :: ln,i,ishell,ispin
      character(len=1) :: shell
      character(len=1), dimension(lmax+1) :: shells
      type(dictionary), pointer :: dict_tmp,iter
      type(atoms_iterator) :: it
      type(f_matrix), dimension(:,:,:), pointer :: gamma_ntypes

      call f_routine(id='atomic_gamma_from_dict')

      sometarget=.false.

      gamma_ntypes=f_malloc_ptr([0.to.lmax,1.to.2,1.to.astruct%ntypes],&
           id='gamma_ntypes')

      !iterate above atoms
      it=atoms_iter(astruct)
      do while(atoms_iter_next(it))
         dogamma(:,it%iat)=.false.
         ! Possible overwrite, if the dictionary has the item
         !get the particular species
         call get_atomic_dict(dict_tmp,dict,it%iat,it%name)
         !according to the dictionary
         if (associated(dict_tmp)) then
            !list case
            ln=dict_len(dict_tmp)
            if (ln > 0) then
               shells(1:ln)=dict_tmp
            else
               ln=1
               shells(1)=dict_tmp
            end if
            do i=1,ln
               dogamma(shell_toi(shells(i)),it%iat)=.true.
            end do
         end if
         !get the particular occupancy control
         call get_atomic_dict(dict_tmp,dict_targets,it%iat,it%name,fromname=fromname)
         if (.not. associated(dict_tmp)) cycle
         nullify(iter)
         do while(iterating(iter,on=dict_tmp))
            if (len_trim(dict_key(iter))/=1) cycle
            call f_strcpy(src=dict_key(iter),dest=shell)
            ishell=shell_toi(shell)
            !let us now see if the matrix has been already extracted
            do ispin=1,2 !HERE we should generalize to spinorial indices
               if (.not. fromname .or. .not. associated(gamma_ntypes(ishell,ispin,it%ityp)%ptr)) &
                    call get_gamma_target(iter,ishell,ispin,gamma_targets(ishell,ispin,it%iat))
               if (fromname) then
                  if (associated(gamma_ntypes(ishell,ispin,it%ityp)%ptr)) then
                     gamma_targets(ishell,ispin,it%iat)=&
                          gamma_ntypes(ishell,ispin,it%ityp)
                  else
                     !store it for future reference
                     gamma_ntypes(ishell,ispin,it%ityp)=&
                          gamma_targets(ishell,ispin,it%iat)
                  end if
               end if
            end do
            dogamma(ishell,it%iat)=.true.
         end do
         sometarget=.true.
      end do

      call f_array_ptr_free(gamma_ntypes)
      if (.not. sometarget) call f_array_ptr_free(gamma_targets)

   call f_release_routine()

  end subroutine atoms_gamma_from_dict

  subroutine get_gamma_target(dict,l,ispin,mat)
    use yaml_strings
    use dictionaries
    use f_blas
    use dynamic_memory
    implicit none
    integer, intent(in) :: ispin,l
    type(dictionary), pointer :: dict
    type(f_matrix), intent(inout) :: mat
    !local variables
    character(len=4) :: aspin
    type(dictionary), pointer :: iter

    select case(ispin)
    case(1)
       call f_strcpy(src='up',dest=aspin)
    case(2)
       call f_strcpy(src='down',dest=aspin)
    end select

    if (aspin .in. dict) then
       iter => dict // aspin
    else
       iter => dict
    end if
    !call f_matrix_allocate(mat,n=2*l+1)
    mat=f_malloc_ptr([2*l+1,2*l+1],id='mat')
    call get_matrix_from_dict(iter,l,mat%ptr)
  end subroutine get_gamma_target

  subroutine get_matrix_from_dict(dict,l,mat)
    use dictionaries
    use yaml_strings
    use f_utils, only: f_zero
    implicit none
    type(dictionary), pointer :: dict
    integer, intent(in) :: l
    real(gp), dimension(2*l+1,2*l+1), intent(out) :: mat
    !local variables
    integer :: i,ln
    real(gp) :: tr

    call f_zero(mat)
    ln=dict_len(dict)
    !first check if the matrix is given row by row
    if (ln > 0) then
       if (ln /= 2*l+1) then
          call f_err_throw(&
            'The number of lines of the matrix should be '//&
            trim(yaml_toa(2*l+1))//', found '//trim(yaml_toa(ln)),&
            err_name='BIGDFT_INPUT_FILE_ERROR')
          return
       end if
       do i=0,2*l
          mat(:,i+1)=dict//i
       end do
    else
       !the matrix is given as a trace
       tr=dict
       do i=1,2*l+1
          mat(i,i)=tr/real(2*l+1,gp)
       end do
    end if
  end subroutine get_matrix_from_dict


    subroutine atomic_data_set_from_dict(dict, key, atoms, nspin)
      use module_defs, only: gp
      use ao_inguess, only: ao_ig_charge,atomic_info,aoig_set_from_dict,&
           print_eleconf,aoig_set
      use dictionaries
      use yaml_output, only: yaml_warning, yaml_dict_dump
      use yaml_strings, only: f_strcpy, yaml_toa
      use module_base, only: bigdft_mpi
      use dynamic_memory
      implicit none
      type(dictionary), pointer :: dict
      type(atoms_data), intent(inout) :: atoms
      character(len = *), intent(in) :: key
      integer, intent(in) :: nspin

      !integer :: iat, ityp
      real(gp) :: elec
      type(dictionary), pointer :: dict_tmp
      type(atoms_iterator) :: it

      call f_routine(id='atomic_data_set_from_dict')

      !number of atoms with semicore channels
      atoms%natsc = 0
      !iterate above atoms
      it=atoms_iter(atoms%astruct)
      !python method
      do while(atoms_iter_next(it))
!!$      !fortran method
!!$      call increment_atoms_iter(it)
!!$      do while(atoms_iter_is_valid(it))
         !only amu is extracted here
         call atomic_info(atoms%nzatom(it%ityp),atoms%nelpsp(it%ityp),&
              amu=atoms%amu(it%ityp))
         !fill the atomic IG configuration from the input_polarization
         !standard form
         atoms%aoig(it%iat)=aoig_set(atoms%nzatom(it%ityp),atoms%nelpsp(it%ityp),&
              atoms%astruct%input_polarization(it%iat),nspin)

         ! Possible overwrite, if the dictionary has the item
         if (key .in. dict) then
            !get the particular species
            call get_atomic_dict(dict_tmp,dict//key,it%iat,it%name)
            !therefore the aiog structure has to be rebuilt
            !according to the dictionary
            if (associated(dict_tmp)) then
               atoms%aoig(it%iat)=aoig_set_from_dict(dict_tmp,nspin,atoms%aoig(it%iat))
               !check the total number of electrons
               elec=ao_ig_charge(nspin,atoms%aoig(it%iat)%aocc)
               if (nint(elec) /= atoms%nelpsp(it%ityp) .and. bigdft_mpi%iproc==0) then
                  call print_eleconf(nspin,atoms%aoig(it%iat)%aocc,atoms%aoig(it%iat)%nl_sc)
                  call yaml_warning('The total atomic charge '//trim(yaml_toa(elec))//&
                       ' is different from the PSP charge '//trim(yaml_toa(atoms%nelpsp(it%ityp))))
               end if
            end if
         end if
         if (atoms%aoig(it%iat)%nao_sc /= 0) atoms%natsc=atoms%natsc+1
!!$         !fortran method
!!$         call increment_atoms_iter(it)
      end do

      call f_release_routine()

    end subroutine atomic_data_set_from_dict


    !> Set irreductible Brillouin zone
    subroutine set_symmetry_data(sym, geocode, n1i, n2i, n3i, nspin)
      use module_base
      use m_ab6_kpoints
      use m_ab6_symmetry
      implicit none
      type(symmetry_data), intent(inout) :: sym
      integer, intent(in) :: n1i, n2i, n3i, nspin
      character, intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode

      character(len = *), parameter :: subname = "symmetry_set_irreductible_zone"
      integer :: i_stat, nsym, i_third
!!$      integer :: i_all
      integer, dimension(:,:,:), allocatable :: irrzon
      real(dp), dimension(:,:,:), allocatable :: phnons

      call f_free_ptr(sym%irrzon)
      call f_free_ptr(sym%phnons)


      if (sym%symObj >= 0) then
         call symmetry_get_n_sym(sym%symObj, nsym, i_stat)
         if (nsym > 1) then
            ! Current third dimension is set to 1 always
            ! since nspin == nsppol always in BigDFT
            i_third = 1
            if (geocode == "S") i_third = n2i
            sym%irrzon=f_malloc_ptr((/n1i*(n2i - i_third + 1)*n3i,2,i_third/),id='sym%irrzon')
            sym%phnons=f_malloc_ptr((/2,n1i*(n2i - i_third + 1)*n3i,i_third/),id='sym%phnons')
            if (geocode /= "S") then
               call kpoints_get_irreductible_zone(sym%irrzon, sym%phnons, &
                    &   n1i, n2i, n3i, nspin, nspin, sym%symObj, i_stat)
            else
               irrzon=f_malloc((/n1i*n3i,2,1/),id='irrzon')
               phnons=f_malloc((/2,n1i*n3i,1/),id='phnons')
               do i_third = 1, n2i, 1
                  call kpoints_get_irreductible_zone(irrzon, phnons, n1i, 1, n3i, &
                       & nspin, nspin, sym%symObj, i_stat)
                  sym%irrzon(:,:,i_third:i_third) = irrzon
                  !call vcopy(2*n1i*n3i, phnons(1,1,1), 1, sym%phnons(1,1,i_third), 1)
                  call f_memcpy(src=phnons,dest=sym%phnons(:,:,i_third))
               end do
               call f_free(irrzon)
               call f_free(phnons)
            end if
         end if
      end if

      if (.not. associated(sym%irrzon)) then
         ! Allocate anyway to small size otherwise the bounds check does not pass.
         sym%irrzon=f_malloc_ptr((/1,2,1/),id='sym%irrzon')
         sym%phnons=f_malloc_ptr((/2,1,1/),id='sym%phnons')
      end if
    END SUBROUTINE set_symmetry_data


    ! allocations, and setters
    !> fill the atomic structure datatype
    subroutine astruct_set(astruct,dict_posinp,randdis,disableSym,symTol,elecfield,nspin,simplify)
      implicit none
      type(atomic_structure), intent(out) :: astruct      !< Contains all info
      type(dictionary), pointer :: dict_posinp
      logical, intent(in) :: simplify
      real(gp), intent(in) :: randdis !< random displacement
      logical, intent(in) :: disableSym
      real(gp), intent(in) :: symtol
      real(gp), dimension(3), intent(in) :: elecfield
      integer, intent(in) :: nspin

      call nullify_atomic_structure(astruct)
      !astruct initialization
      call astruct_set_from_dict(dict_posinp, astruct)
      ! Shake atoms, if required.
      call astruct_set_displacement(astruct, randdis)
      call astruct_set_symmetries(astruct, disableSym,symTol,elecfield,nspin)
      call check_atoms_positions(astruct,simplify)
    end subroutine astruct_set


    !> Read atomic file
    subroutine set_astruct_from_file(file,iproc,astruct,comment,energy,fxyz,disableTrans,pos_format)
      use module_base
      use dictionaries, only: set, dictionary
      use yaml_strings
      use f_utils, only: f_zero
      use internal_coordinates, only: internal_to_cartesian
      use yaml_output, only: yaml_dict_dump
      use yaml_parse, only: yaml_parse_from_file,yaml_load
      implicit none
      !Arguments
      character(len=*), intent(in) :: file                  !< File name containing the atomic positions
      integer, intent(in) :: iproc
      type(atomic_structure), intent(inout) :: astruct      !< Contains all info
      real(gp), intent(out), optional :: energy
      real(gp), dimension(:,:), pointer, optional :: fxyz
      character(len=*), intent(in), optional :: pos_format  !< Explicit format of the file
      character(len = 1024), intent(out), optional :: comment
      logical, intent(in), optional :: disableTrans
      !Local variables
      integer :: iunit
      integer :: l, extract
      logical :: file_exists, archive
      character(len = 128) :: filename
      character(len = 15) :: arFile
      character(len = 6) :: ext
      real(gp) :: energy_
      type(dictionary), pointer :: yaml_file_dict,dict_extensions,iter
      real(gp), dimension(:,:), pointer :: fxyz_
      character(len = 1024) :: comment_, files
      external :: openNextCompress
      real(gp), parameter :: degree = 57.295779513d0

      iunit=99
      file_exists = .false.
      call f_zero(files)
      archive = .false.
      nullify(fxyz_)

      ! Extract from archive (posout_)
      if (index(file, "posout_") == 1 .or. index(file, "posmd_") == 1) then
         write(arFile, "(A)") "posout.tar.bz2"
         if (index(file, "posmd_") == 1) write(arFile, "(A)") "posmd.tar.bz2"
         call f_file_exists(arFile,file_exists)
         !arFile tested
         if (file_exists) then
!!$     call extractNextCompress(trim(arFile), len(trim(arFile)), &
!!$          & trim(file), len(trim(file)), extract, ext)
            call openNextCompress(trim(arFile), len(trim(arFile)), &
                 & trim(file), len(trim(file)), extract, ext)
            if (f_err_raise(extract == 0, &
                  & "Can't find '"//trim(file) //"' in archive.", &
                  & err_id=BIGDFT_INPUT_FILE_ERROR)) return
            archive = .true.
            write(filename, "(A)") file//'.'//trim(ext)
            write(astruct%inputfile_format, "(A)") trim(ext)
            write(astruct%source, "(A)") trim(filename)
         end if
      end if

!!$      ! Test file//'.xyz'
!!$      if (.not. file_exists) then
!!$         call f_file_exists(file//'.xyz',file_exists)
!!$         !files = trim(files) // "'" // trim(file)//".xyz'"
!!$         files = files+"'"+file+".xyz'"
!!$         if (file_exists) then
!!$            ext='xyz'
!!$            !write(filename, "(A)") file//'.xyz'!"posinp.xyz"
!!$            !write(astruct%inputfile_format, "(A)") "xyz"
!!$            !write(astruct%source, "(A)") trim(filename)
!!$            call f_open_file(unit=iunit,file=trim(filename),status='old')
!!$         end if
!!$      end if
!!$
!!$      ! Test file//'.ascii'
!!$      if (.not. file_exists) then
!!$         call f_file_exists(file//'.ascii',file_exists)
!!$         !files = trim(files) // ", '" //trim(file)//".ascii'"
!!$         files = files+"'"+file+".ascii'"
!!$         if (file_exists) then
!!$            ext='ascii'
!!$            !write(filename, "(A)") file//'.ascii'!"posinp.ascii"
!!$            !write(astruct%inputfile_format, "(A)") "ascii"
!!$            !write(astruct%source, "(A)") trim(filename)
!!$            call f_open_file(unit=iunit,file=trim(filename),status='old')
!!$         end if
!!$      end if
!!$      ! Test file//'.int'
!!$      if (.not. file_exists) then
!!$         call f_file_exists(file//'.int',file_exists)
!!$         files = files+"'"+file+".int'"
!!$         if (file_exists) then
!!$            ext='int'
!!$            !write(filename, "(A)") file//'.int'!"posinp.int
!!$            !write(astruct%inputfile_format, "(A)") "int"
!!$            !write(astruct%source, "(A)") trim(filename)
!!$            call f_open_file(unit=iunit,file=trim(filename),status='old')
!!$         end if
!!$      end if
!!$      ! Test file//'.yaml'
!!$      if (.not. file_exists) then
!!$         call f_file_exists(file//'.yaml',file_exists)
!!$         files = trim(files) // ", '" //trim(file)//".yaml'"
!!$         if (file_exists) then
!!$            write(filename, "(A)") file//'.yaml'!"posinp.yaml
!!$            write(astruct%inputfile_format, "(A)") "yaml"
!!$            write(astruct%source, "(A)") trim(filename)
!!$            ! Pb if toto.yaml because means that there is no key posinp!!
!!$         end if
!!$      end if

      if (.not. file_exists) then
         !general search for the file extension
         dict_extensions=>yaml_load('[xyz,ascii,int,yaml]')
         nullify(iter)
         find_file: do while(iterating(iter,on=dict_extensions))
            ext=trim(dict_value(iter))
            files = files+"'"+file+"."+ext+"',"
            call f_file_exists(file+'.'+ext,file_exists)
            if (file_exists) then
               write(filename, "(A)") file//'.'//trim(ext)
               write(astruct%inputfile_format, "(A)") trim(ext)
               write(astruct%source, "(A)") trim(filename)
               exit find_file
            end if
         end do find_file
         call dict_free(dict_extensions)
      end if


      ! Check if the format of the file detected corresponds to the specified format
      if (file_exists .and. present(pos_format)) then
         if (astruct%inputfile_format /= pos_format) &
              call f_err_throw("The detected filename '"//trim(filename)//"' has not the specified format '" // &
              trim(pos_format) // "'.", &
              err_name='BIGDFT_INPUT_VARIABLES_ERROR')
      end if

      ! Test the name directly
      if (.not. file_exists) then
         files = files+"'"+file+"',"
         call f_file_exists(file,file_exists)
         if (file_exists) then
            write(filename, "(A)") file
            write(astruct%source, "(A)") trim(filename)
            if (.not.present(pos_format)) then
               l = len(file)
               if (file(l-3:l) == ".xyz") then
                  write(astruct%inputfile_format, "(A)") "xyz"
               else if (file(l-5:l) == ".ascii") then
                  write(astruct%inputfile_format, "(A)") "ascii"
               else if (file(l-3:l) == ".int") then
                  write(astruct%inputfile_format, "(A)") "int"
               else if (file(l-4:l) == ".yaml") then
                  write(astruct%inputfile_format, "(A)") "yaml"
               else
                  !We assume that the format of the file is 'xyz'
                  write(astruct%inputfile_format, "(A)") "xyz"
               end if
            else
               ! The format is specified
               write(astruct%inputfile_format, "(A)") trim(pos_format)
            end if
!!$            if (trim(astruct%inputfile_format) /= "yaml") then
!!$               call f_open_file(unit=iunit,file=trim(filename),status='old')
!!$            end if
         end if
      end if

      if (f_err_raise(.not.file_exists, &
           "Atomic input file not found. Files looked for were: "//&
           trim(files) //"but none matched.", &
           err_id=BIGDFT_INPUT_FILE_ERROR)) return

      !We found a file
      select case (astruct%inputfile_format)
      case("xyz")
         call f_open_file(unit=iunit,file=trim(filename),status='old')
         !read atomic positions
         if (.not.archive) then
            call read_xyz_positions(iunit,filename,astruct,comment_,energy_,fxyz_,directGetLine,disableTrans)
         else
            call read_xyz_positions(iunit,filename,astruct,comment_,energy_,fxyz_,archiveGetLine,disableTrans)
         end if

      case("ascii")
         call f_open_file(unit=iunit,file=trim(filename),status='old')
         !read atomic positions
         if (.not.archive) then
            call read_ascii_positions(iunit,filename,astruct,comment_,energy_,fxyz_,directGetLine,disableTrans)
         else
            call read_ascii_positions(iunit,filename,astruct,comment_,energy_,fxyz_,archiveGetLine,disableTrans)
         end if

      case("int")
         call f_open_file(unit=iunit,file=trim(filename),status='old')
         !read atomic positions
         if (.not.archive) then
            call read_int_positions(iproc,iunit,astruct,comment_,energy_,fxyz_,directGetLine,disableTrans)
         else
            call read_int_positions(iproc,iunit,astruct,comment_,energy_,fxyz_,archiveGetLine,disableTrans)
         end if
         ! Fill the ordinary rxyz array
         !!! convert to rad
         !!astruct%rxyz_int(2:3,1:astruct%nat) = astruct%rxyz_int(2:3,1:astruct%nat) / degree
         ! The bond angle must be modified (take 180 degrees minus the angle)
         astruct%rxyz_int(2:2,1:astruct%nat) = pi_param - astruct%rxyz_int(2:2,1:astruct%nat)
         call internal_to_cartesian(astruct%nat, astruct%ixyz_int(1,:), astruct%ixyz_int(2,:), astruct%ixyz_int(3,:), &
              astruct%rxyz_int, astruct%rxyz)

       case("yaml")
          nullify(yaml_file_dict)
          call f_zero(energy_)
          call yaml_parse_from_file(yaml_file_dict,trim(filename))
          call astruct_set_from_dict(yaml_file_dict//0, astruct, comment_)
          call dict_free(yaml_file_dict)
      case default
         call f_err_throw(err_msg="The specified format '" // trim(astruct%inputfile_format) // "' is not recognised."// &
            & " The format should be 'yaml', 'int', 'ascii' or 'xyz'.",err_id=BIGDFT_INPUT_FILE_ERROR)

      end select
      !if an error has been produced return
      if (f_err_check()) return

      !Check the number of atoms (should be more a f_assert that a error raise)
      if (f_err_raise(astruct%nat < 0, &
           "In the file '"//trim(filename)//"' the number of atoms ("// &
           trim(yaml_toa(astruct%nat))//") should be >= 0.",err_id=BIGDFT_INPUT_VARIABLES_ERROR)) return

      ! We delay the calculation of the symmetries.
      !this should be already in the atoms_null routine
      astruct%sym=symmetry_data_null()

      ! close open file.
      if (.not.archive .and. trim(astruct%inputfile_format) /= "yaml") then
         close(iunit)
!!$  else
!!$     call unlinkExtract(trim(filename), len(trim(filename)))
      end if

      ! We transfer optionals.
      if (present(energy)) then
         energy = energy_
      end if
      if (present(comment)) then
         write(comment, "(A)") comment_
      end if
      if (present(fxyz)) then
         if (associated(fxyz_)) then
            !fxyz = f_malloc_ptr(src = fxyz_, id = "fxyz") not needed anymore
            fxyz => fxyz_
         else
            nullify(fxyz)
         end if
         !call f_free_ptr(fxyz_)
      else if (associated(fxyz_)) then
         call f_free_ptr(fxyz_)
      end if

    END SUBROUTINE set_astruct_from_file

    subroutine set_astruct_from_openbabel(astruct, obfile)
      use dictionaries
      implicit none
      type(atomic_structure), intent(out) :: astruct
      character(len = *), intent(in) :: obfile

      type(dictionary), pointer :: dict

      interface
         subroutine openbabel_load(d, f)
           use dictionaries
           implicit none
           type(dictionary), pointer :: d
           character(len = *), intent(in) :: f
         end subroutine openbabel_load
      end interface

      call dict_init(dict)
      call openbabel_load(dict, obfile)
      call astruct_set_from_dict(dict, astruct)
      call dict_free(dict)
    end subroutine set_astruct_from_openbabel

    !> Write an atomic file
    !! Yaml output included
    subroutine astruct_dump_to_file(astruct,filename,comment,energy,rxyz,forces,&
         fmt,unit,position)
      use module_base
      use yaml_output
      use yaml_strings, only: f_strcpy
      use internal_coordinates, only : xyzint
      implicit none
      character(len=*), intent(in) :: filename,comment
      type(atomic_structure), intent(in) :: astruct
      real(gp), intent(in), optional :: energy
      real(gp), dimension(3,astruct%nat), intent(in), target, optional :: rxyz
      real(gp), dimension(3,astruct%nat), intent(in), optional :: forces
      !> unit of file writing. Has to be already opened
      !! if filename is 'stdout' then this value is ignored
      !! otherwise, the filename is ignored
      integer, intent(in), optional :: unit
      !> force the format of the output
      !! the default is otherwise used as defined in inputfile_format
      character(len=*), intent(in), optional :: fmt
      !> position of the file at opening. Default is "rewind"
      character(len=*), intent(in), optional :: position

      !local variables
      character(len=6) :: pos
      character(len = 15) :: arFile
      integer :: iunit
      character(len = 1024) :: fname
      real(gp), dimension(3), parameter :: dummy = (/ 0._gp, 0._gp, 0._gp /)
      type(dictionary), pointer :: dict
      real(gp),dimension(:,:),allocatable :: rxyz_int
      real(kind=8),parameter :: degree=1.d0
      real(gp) :: energy_
      real(gp), dimension(:,:), pointer :: rxyz_
      character(len=len(astruct%inputfile_format)) :: formt

      call f_routine(id='astruct_dump_to_file')

      energy_=0.0_gp
      if (present(energy)) energy_ = energy
      rxyz_ => astruct%rxyz
      if (present(rxyz)) rxyz_ => rxyz
      formt=astruct%inputfile_format
      if (present(fmt)) call f_strcpy(src=fmt,dest=formt)
      iunit=f_get_free_unit(9)
      if (present(unit)) iunit=unit
      pos='rewind'
      if (present(position)) call f_strcpy(src=position,dest=pos)

      if (trim(filename) == "stdout") then
         iunit = 6
      else if (.not. present(unit)) then
         !also unit opening should be checked
         write(fname,"(A)") trim(filename)//'.'//trim(formt)
         if (formt == 'yaml') then
            call yaml_set_stream(unit = iunit, filename = trim(fname), &
                 & record_length = 92, tabbing = 0, setdefault = .false.,position=pos)
         else
            call yaml_set_stream(unit = iunit, filename = trim(fname), &
                 & record_length = 4096, tabbing = 0, setdefault = .false.,position=pos)
         end if
      end if

      select case(formt)
      case('xyz')
         call wtxyz(iunit,energy_,rxyz_,astruct,comment)
         if (present(forces)) call wtxyz_forces(iunit,forces,astruct)
      case('ascii')
         call wtascii(iunit,energy_,rxyz_,astruct,comment)
         if (present(forces)) call wtascii_forces(iunit,forces,astruct)
      case ('int')
         rxyz_int = f_malloc((/3,astruct%nat/),id='rxyz_int')
         !call wtint(iunit,energy_,rxyz,astruct,comment,ixyz(1,:),ixyz(2,:),ixyz(3,:))
         call xyzint(rxyz_, astruct%nat, &
              astruct%ixyz_int(1,:), astruct%ixyz_int(2,:),astruct%ixyz_int(3,:), &
              degree, rxyz_int)
         call wtint(iunit,energy_,rxyz_int,astruct,comment,astruct%ixyz_int(1,:),astruct%ixyz_int(2,:),astruct%ixyz_int(3,:))
         call f_free(rxyz_int)
      case('yaml')
         call yaml_new_document(unit = iunit)
         call dict_init(dict)
         call astruct_merge_to_dict(dict, astruct, rxyz_, comment)
         call yaml_dict_dump(dict, unit = iunit)
         call dict_free(dict)
      case default
         call f_err_throw('Writing the atomic file. Error, unknown file format ("'//&
              trim(astruct%inputfile_format)//'")', &
              err_name='BIGDFT_RUNTIME_ERROR')
      end select

      if (iunit /= 6 .and. .not. present(unit)) then
         call yaml_close_stream(unit = iunit)
         ! Add to archive
         if (index(filename, "posout_") == 1 .or. index(filename, "posmd_") == 1) then
            write(arFile, "(A)") "posout.tar.bz2"
            if (index(filename, "posmd_") == 1) write(arFile, "(A)") "posmd.tar.bz2"
            call addToCompress(trim(arFile), len(trim(arFile)), trim(fname), len(trim(fname)))
         end if
      end if

      call f_release_routine()

    END SUBROUTINE astruct_dump_to_file


    !> Convert astruct to dictionary for later dump.
    subroutine astruct_merge_to_dict(dict, astruct, rxyz, comment)
      use dynamic_memory
      use module_defs, only: gp, UNINITIALIZED
      use numerics, only: Bohr_Ang
      use dictionaries
      use yaml_strings
      use box, only: geocode_to_bc,bc_periodic_dims
      use ao_inguess, only: charge_and_spol
      use yaml_output !tmp
      
      implicit none
      type(dictionary), pointer :: dict
      type(atomic_structure), intent(in) :: astruct
      real(gp), dimension(3, astruct%nat), intent(in) :: rxyz
      character(len=*), intent(in), optional :: comment
      !local variables
      type(dictionary), pointer :: pos, at, last
      integer :: iat,ichg,ispol,i
      real(gp) :: factor(3)
      logical, dimension(3) :: peri
      logical :: reduced
      character(len = 4) :: frzstr

      call f_routine(id='astruct_merge_to_dict')

      !call dict_init(dict)

      reduced = .false.
      factor=1.0_gp
      Units: select case(trim(astruct%units))
      case('angstroem','angstroemd0')
         call set(dict // ASTRUCT_UNITS, 'angstroem')
         factor=Bohr_Ang
      case('reduced')
         call set(dict // ASTRUCT_UNITS, 'reduced')
         reduced = .true.
      case('atomic','atomicd0','bohr','bohrd0')
         ! Default, store nothing
      end select Units

      peri=bc_periodic_dims(geocode_to_bc(astruct%geocode))
      do i=1,3
         if (peri(i)) then         
            call set(dict // ASTRUCT_CELL // (i-1), yaml_toa(astruct%cell_dim(i)*factor(i)))
            if (reduced) factor(i) = 1._gp / astruct%cell_dim(i)
         else
            call set(dict // ASTRUCT_CELL // (i-1), '.inf')
         end if
      end do

      if (.not. any(peri)) then
         ! Default, store nothing and erase key if already exist.
         if (ASTRUCT_CELL .in. dict) then
            last=>dict .pop. ASTRUCT_CELL
            call dict_free(last)
            if (.not. associated(dict)) call dict_init(dict)
         end if
!!$         if (has_key(dict, ASTRUCT_CELL)) call dict_remove(dict, ASTRUCT_CELL,destroy=.false.)
      end if

      !cell information
!!$      BC :select case(astruct%geocode)
!!$      case('S')
!!$         call set(dict // ASTRUCT_CELL // 0, yaml_toa(astruct%cell_dim(1)*factor(1)))
!!$         call set(dict // ASTRUCT_CELL // 1, '.inf')
!!$         call set(dict // ASTRUCT_CELL // 2, yaml_toa(astruct%cell_dim(3)*factor(3)))
!!$         !angdeg to be added
!!$         if (reduced) then
!!$            factor(1) = 1._gp / astruct%cell_dim(1)
!!$            factor(3) = 1._gp / astruct%cell_dim(3)
!!$         end if
!!$      case('W')
!!$         call set(dict // ASTRUCT_CELL // 0, '.inf')
!!$         call set(dict // ASTRUCT_CELL // 1, '.inf')
!!$         call set(dict // ASTRUCT_CELL // 2, yaml_toa(astruct%cell_dim(3)*factor(3)))
!!$         if (reduced) then
!!$            factor(3) = 1._gp / astruct%cell_dim(3)
!!$         end if
!!$      case('P')
!!$         call set(dict // ASTRUCT_CELL // 0, yaml_toa(astruct%cell_dim(1)*factor(1)))
!!$         call set(dict // ASTRUCT_CELL // 1, yaml_toa(astruct%cell_dim(2)*factor(2)))
!!$         call set(dict // ASTRUCT_CELL // 2, yaml_toa(astruct%cell_dim(3)*factor(3)))
!!$         !angdeg to be added
!!$         if (reduced) then
!!$            factor(1) = 1._gp / astruct%cell_dim(1)
!!$            factor(2) = 1._gp / astruct%cell_dim(2)
!!$            factor(3) = 1._gp / astruct%cell_dim(3)
!!$         end if
!!$      case('F')
!!$         ! Default, store nothing and erase key if already exist.
!!$         if (has_key(dict, ASTRUCT_CELL)) call dict_remove(dict, ASTRUCT_CELL)
!!$      end select BC

      !if (has_key(dict, ASTRUCT_POSITIONS)) call dict_remove(dict, ASTRUCT_POSITIONS,destroy=.false.)
      if (ASTRUCT_POSITIONS .in. dict) then
         last=>dict .pop. ASTRUCT_POSITIONS
         call dict_free(last)
         if (.not. associated(dict)) call dict_init(dict)
      end if
      if (astruct%nat > 0) pos => dict // ASTRUCT_POSITIONS
      nullify(last)
      do iat=1,astruct%nat
         call dict_init(at)
         call add(at // astruct%atomnames(astruct%iatype(iat)), rxyz(1,iat) * factor(1))
         call add(at // astruct%atomnames(astruct%iatype(iat)), rxyz(2,iat) * factor(2))
         call add(at // astruct%atomnames(astruct%iatype(iat)), rxyz(3,iat) * factor(3))
         if (astruct%ifrztyp(iat) /= 0) then
            call frozen_itof(astruct%ifrztyp(iat), frzstr)
            call set(at // ASTRUCT_ATT_FROZEN, adjustl(frzstr))
         end if
         call charge_and_spol(astruct%input_polarization(iat),ichg,ispol)
         if (ichg /= 0) call set(at // ASTRUCT_ATT_IGCHRG, ichg)
         if (ispol /= 0) call set(at // ASTRUCT_ATT_IGSPIN, ispol)
         ! information for internal coordinates
         if (astruct%inputfile_format=='int') then
            call set(at // ASTRUCT_ATT_IXYZ_1, astruct%ixyz_int(1,iat))
            call set(at // ASTRUCT_ATT_IXYZ_2, astruct%ixyz_int(2,iat))
            call set(at // ASTRUCT_ATT_IXYZ_3, astruct%ixyz_int(3,iat))
            call set(at // ASTRUCT_ATT_RXYZ_INT_1, astruct%rxyz_int(1,iat))
            call set(at // ASTRUCT_ATT_RXYZ_INT_2, astruct%rxyz_int(2,iat))
            call set(at // ASTRUCT_ATT_RXYZ_INT_3, astruct%rxyz_int(3,iat))
         end if
         call dict_update(at, astruct%attributes(iat)%d)
         call add(pos, at, last)
      end do

      if (associated(astruct%properties)) then
         call dict_update(dict // ASTRUCT_PROPERTIES, astruct%properties)
      end if

      if (present(comment)) then
         if (len_trim(comment) > 0) then
             call add(dict // ASTRUCT_PROPERTIES // "info", comment)
         end if
      end if

      if (len_trim(astruct%inputfile_format) > 0) &
           & call set(dict // ASTRUCT_PROPERTIES // "format", astruct%inputfile_format)

      call f_release_routine()

    end subroutine astruct_merge_to_dict


    subroutine astruct_at_from_dict(dict, symbol, rxyz, rxyz_add, ifrztyp, igspin, igchrg, &
               ixyz, ixyz_add, rxyz_int, rxyz_int_add, cavity_radius,mode)
      use dictionaries
      use module_defs, only: UNINITIALIZED
      use dynamic_memory
      use f_enums, only: operator(==), f_int => toi
      implicit none
      type(dictionary), pointer :: dict
      character(len = max_field_length), intent(out), optional :: symbol !< Symbol
      integer, intent(out), optional :: ifrztyp !< Frozen id
      integer, intent(out), optional :: igspin  !< Spin for input guess
      integer, intent(out), optional :: igchrg  !< Charge for input guess
      integer, dimension(3), intent(out), optional :: ixyz !< Reference atom for internal coordinates
      integer, intent(out), optional :: ixyz_add !< Reference atom for internal coordinates address
      real(gp), dimension(3), intent(out), optional :: rxyz !< Coordinates.
      real(gp), intent(out), optional :: rxyz_add !< Coordinates address.
      real(gp), dimension(3), intent(out), optional :: rxyz_int !< Internal coordinates.
      real(gp), intent(out), optional :: rxyz_int_add !< Internal coordinates address.
      real(gp), intent(out), optional :: cavity_radius !< radius of the cavity
      character(len = max_field_length), intent(out), optional :: mode !< QM/MM treatment.

      type(dictionary), pointer :: atData
      character(len = max_field_length) :: str
      integer :: ierr
      integer, dimension(3) :: icoord
      real(gp), dimension(3) :: rcoord
      real(gp), dimension(3) :: rcoord_int

      ! Default values.
      if (present(symbol)) symbol = 'X'
      if (present(rxyz)) rxyz = UNINITIALIZED(rxyz(1))
      if (present(ixyz)) ixyz = UNINITIALIZED(ixyz(1))
      if (present(ifrztyp)) ifrztyp = 0
      if (present(igspin))  igspin = 0
      if (present(igchrg))  igchrg = 0
      if (present(mode)) write(mode, "(A)") ""
      if (present(cavity_radius)) cavity_radius=UNINITIALIZED(cavity_radius)

      atData => dict_iter(dict)
      do while(associated(atData))
         str = dict_key(atData)
         select case(trim(str))
         case(ASTRUCT_ATT_FROZEN)
            !if (trim(str) == ASTRUCT_ATT_FROZEN) then
            str = dict_value(atData)
            if (present(ifrztyp)) call frozen_ftoi(str(1:4), ifrztyp, ierr)
         case(ASTRUCT_ATT_IGSPIN)
            if (present(igspin)) igspin = atData
         case(ASTRUCT_ATT_IGCHRG)
            if (present(igchrg)) igchrg = atData
         case(ASTRUCT_ATT_IXYZ_1)
            if (present(ixyz)) ixyz(1) = atData
            if (present(ixyz_add)) then
               call f_memcpy(icoord(1), ixyz_add, 3)
               icoord(1) = atData
               call f_memcpy(ixyz_add, icoord(1), 3)
            end if
         case(ASTRUCT_ATT_IXYZ_2)
            if (present(ixyz)) ixyz(2) = atData
            if (present(ixyz_add)) then
               call f_memcpy(icoord(1), ixyz_add, 3)
               icoord(2) = atData
               call f_memcpy(ixyz_add, icoord(1), 3)
            end if
         case(ASTRUCT_ATT_IXYZ_3)
            if (present(ixyz)) ixyz(3) = atData
            if (present(ixyz_add)) then
               call f_memcpy(icoord(1), ixyz_add, 3)
               icoord(3) = atData
               call f_memcpy(ixyz_add, icoord(1), 3)
            end if
         case(ASTRUCT_ATT_MODE)
            if (present(mode)) then
               mode = dict_value(atData)
            end if
         case(ASTRUCT_ATT_RXYZ_INT_1)
            if (present(rxyz_int)) rxyz_int(1) = atData
            if (present(rxyz_int_add)) then
               call f_memcpy(rcoord_int(1), rxyz_int_add, 3)
               rcoord_int(1) = atData
               call f_memcpy(rxyz_int_add, rcoord_int(1), 3)
            end if
         case(ASTRUCT_ATT_RXYZ_INT_2)
            if (present(rxyz_int)) rxyz_int(2) = atData
            if (present(rxyz_int_add)) then
               call f_memcpy(rcoord_int(1), rxyz_int_add, 3)
               rcoord_int(2) = atData
               call f_memcpy(rxyz_int_add, rcoord_int(1), 3)
            end if
         case(ASTRUCT_ATT_RXYZ_INT_3)
            if (present(rxyz_int)) rxyz_int(3) = atData
            if (present(rxyz_int_add)) then
               call f_memcpy(rcoord_int(1), rxyz_int_add, 3)
               rcoord_int(3) = atData
               call f_memcpy(rxyz_int_add, rcoord_int(1), 3)
            end if
         case(ASTRUCT_ATT_CAVRAD)
            if (present(cavity_radius)) cavity_radius= atData
         case default
            ! Heuristic to find symbol and coordinates.
            if (dict_len(atData) == 3 .and. &
                 & (len_trim(str) < 3 .or. (index(str, "_") > 0 .and. &
                 & index(str, "_") < 4))) then
               if (present(symbol)) symbol = str
               if (present(rxyz)) rxyz = atData
               if (present(rxyz_add)) then
                  rcoord = atData
                  call f_memcpy(rxyz_add, rcoord(1), 3)
               end if
            end if
         end select
         atData => dict_next(atData)
      end do
    end subroutine astruct_at_from_dict


    !> Get types of the atoms in the dictionary
    subroutine astruct_dict_get_types(dict, types)
      use dictionaries
      implicit none
      !Arguments
      type(dictionary), pointer :: dict, types
      !Local variables
      type(dictionary), pointer :: atoms
      character(len = max_field_length) :: str
      integer :: ityp

      if (ASTRUCT_POSITIONS .notin. dict) then
         nullify(types)
         return
      end if
      call dict_init(types)
      ityp = 0
      atoms => dict_iter(dict // ASTRUCT_POSITIONS)
      do while(associated(atoms))
         call astruct_at_from_dict(atoms, symbol = str)
         !if (.not. has_key(types, str)) then
         if (str .notin. types) then
            ityp = ityp + 1
            call set(types // str, ityp)
         end if
         atoms => dict_next(atoms)
      end do
    end subroutine astruct_dict_get_types


    !> complete the atoms structure with the remaining information.
    !! must be called after the call to astruct_set
    subroutine atoms_fill(atoms,dict,nspin,multipole_preserving,mp_isf,ixc,alpha_hartree_fock)
      use ao_inguess, only: lmax_ao
      use public_keys, only: IG_OCCUPATION,OUTPUT_VARIABLES,ATOMIC_DENSITY_MATRIX,&
           DFT_VARIABLES,OCCUPANCY_CONTROL
      use public_enums, only: PSPCODE_PAW
      use abi_interfaces_add_libpaw, only : abi_pawinit
      use module_xc
      use numerics, only: pi_param => pi
      use m_ab6_symmetry
      use dynamic_memory, only: f_routine,f_release_routine
      use dictionaries
      implicit none
      integer, intent(in) :: nspin
      integer, intent(in) :: mp_isf               !< Interpolating scaling function order for multipole preserving
      integer, intent(in) :: ixc         !< XC functional Id
      logical, intent(in) :: multipole_preserving !< Preserve multipole for ionic charge (integrated isf)
      real(gp), intent(in) :: alpha_hartree_fock !< exact exchange contribution
      type(atoms_data), intent(inout) :: atoms
      type(dictionary), pointer :: dict
      !local variables
      integer, parameter :: pawlcutd = 10, pawlmix = 10, pawnphi = 13, pawntheta = 12, pawxcdev = 1
      integer, parameter :: usepotzero = 0
      integer :: iat,ierr,mpsang,nsym,xclevel
      real(gp) :: gsqcut_shp
      type(xc_info) :: xc
      type(dictionary), pointer :: dict_targets,dict_gamma

      call f_routine(id='atoms_fill')

      !fill the rest of the atoms structure
      call psp_dict_analyse(dict, atoms)
      call atomic_data_set_from_dict(dict,IG_OCCUPATION, atoms,nspin)
      !use get as the keys might not be there
      dict_gamma= dict .get. OUTPUT_VARIABLES
      dict_gamma= dict_gamma .get. ATOMIC_DENSITY_MATRIX
      dict_targets= dict .get. DFT_VARIABLES
      dict_targets= dict_targets .get. OCCUPANCY_CONTROL
      call atoms_gamma_from_dict(dict_gamma,dict_targets,&
           lmax_ao,atoms%astruct,atoms%dogamma,atoms%gamma_targets)
      ! Add multipole preserving information
      atoms%multipole_preserving = multipole_preserving
      atoms%mp_isf = mp_isf
      ! Complement PAW initialisation.
      if (any(atoms%npspcode == PSPCODE_PAW)) then
         call xc_init(xc, ixc, XC_MIXED, 1, alpha_hartree_fock)
         xclevel = 1 ! xclevel=XC functional level (1=LDA, 2=GGA)
         if (xc_isgga(xc)) xclevel = 2
         call xc_end(xc)
         !gsqcut_shp = two*abs(dtset%diecut)*dtset%dilatmx**2/pi**2
         gsqcut_shp = 2._gp * 2.2_gp / pi_param ** 2
         nsym = 0
         call symmetry_get_n_sym(atoms%astruct%sym%symObj, nsym, ierr)
         mpsang = -1
         do iat = 1, atoms%astruct%nat
            mpsang = max(mpsang, maxval(atoms%pawtab(iat)%orbitals))
         end do
         call abi_pawinit(1, gsqcut_shp, pawlcutd, pawlmix, mpsang + 1, &
              & pawnphi, nsym, pawntheta, atoms%pawang, atoms%pawrad, 0, &
              & atoms%pawtab, pawxcdev, xclevel, usepotzero)
      end if

      call f_release_routine()

    end subroutine atoms_fill


    !> Read old psppar file (check if not already in the dictionary) and merge to dict
    subroutine atoms_file_merge_to_dict(dict)
      use dictionaries
      use yaml_output, only: yaml_warning
      use public_keys, only: POSINP,SOURCE_KEY,PSPXC_KEY
      use pseudopotentials!, only: update_psp_dict
      implicit none
      type(dictionary), pointer :: dict
      !local variables
      type(dictionary), pointer :: types,iter
      !character(len = max_field_length) :: str
      !integer :: iat, stypes
      !character(len=max_field_length), dimension(:), allocatable :: keys
      !character(len=27) :: key
      !logical :: exists

      if (POSINP .notin. dict) return
      ! Loop on types for atomic data.
      call astruct_dict_get_types(dict // POSINP, types)
      nullify(iter)
      do while(iterating(iter,on=types))
         call update_psp_dict(dict,dict_key(iter))
      end do
!!$      if ( .not. associated(types)) return
!!$      allocate(keys(dict_size(types)))
!!$      keys = dict_keys(types)
!!$      stypes = dict_size(types)
!!$      do iat = 1, stypes, 1
!!$         key = 'psppar.' // trim(keys(iat))
!!$
!!$         exists = has_key(dict, key)
!!$
!!$         if (exists) then
!!$            if (has_key(dict // key, SOURCE_KEY)) then
!!$               str = dict_value(dict // key // SOURCE_KEY)
!!$            else
!!$               str = dict_value(dict // key)
!!$            end if
!!$            if (trim(str) /= "" .and. trim(str) /= TYPE_DICT) then
!!$               !Read the PSP file and merge to dict
!!$               if (trim(str) /= TYPE_LIST) then
!!$                  call psp_file_merge_to_dict(dict, key, filename = trim(str))
!!$               else
!!$                  call psp_file_merge_to_dict(dict, key, lstring = dict // key)
!!$               end if
!!$               if (.not. has_key(dict // key, PSPXC_KEY)) then
!!$                  call yaml_warning("Pseudopotential file '" // trim(str) // &
!!$                       & "' not found. Fallback to file '" // trim(key) // &
!!$                       & "' or hard-coded pseudopotential.")
!!$               end if
!!$            end if
!!$            exists = has_key(dict // key, PSPXC_KEY)
!!$         end if
!!$         if (.not. exists) call psp_file_merge_to_dict(dict, key, key)
!!$
!!$         exists = has_key(dict, key)
!!$         if (exists) exists = has_key(dict // key, 'Non Linear Core Correction term')
!!$         if (.not. exists) call nlcc_file_merge_to_dict(dict, key, 'nlcc.' // trim(keys(iat)))
!!$      end do
!!$      deallocate(keys)
      call dict_free(types)
    end subroutine atoms_file_merge_to_dict


    !> Read Atomic positions and merge into dict
    subroutine astruct_file_merge_to_dict(dict, key, filename, pos_format)
      use module_base, only: gp, UNINITIALIZED, bigdft_mpi,f_routine,f_release_routine, &
          & BIGDFT_INPUT_FILE_ERROR,f_free_ptr
      use public_keys, only: POSINP,GOUT_FORCES,GOUT_ENERGY,POSINP_SOURCE
      use yaml_strings
      use dictionaries
      use module_input_dicts, only: dict_get_run_properties
      implicit none
      !Arguments
      type(dictionary), pointer :: dict                      !< Contains (out) all the information
      character(len = *), intent(in) :: key                  !< Key of the dictionary where it should be have the information
      character(len = *), intent(in) :: filename             !< Name of the filename where the astruct should be read
      character(len = *), optional, intent(in) :: pos_format !< Explicit specified format
      !Local variables
      type(atomic_structure) :: astruct
      !type(DFT_global_output) :: outs
      character(len=max_field_length) :: msg,radical
      integer :: ierr,iat
      real(gp) :: energy
      real(gp), dimension(:,:), pointer :: fxyz
      type(dictionary), pointer :: dict_tmp,pos

      call f_routine(id='astruct_file_merge_to_dict')

      call nullify_atomic_structure(astruct)

      !Try to read the atomic coordinates from files (old way)
      call f_err_open_try()
      nullify(fxyz)
      ! if (present(pos_format)) then !not needed by fortran spec
      call set_astruct_from_file(filename, bigdft_mpi%iproc, astruct, &
           energy = energy, fxyz = fxyz, pos_format=pos_format)
      !else
      !call set_astruct_from_file(filename, bigdft_mpi%iproc, astruct, &
      !       energy = energy, fxyz = fxyz)
      !end if
      !print *,'test2',associated(fxyz)
      !Check if BIGDFT_INPUT_FILE_ERROR
      ierr = f_get_last_error(msg)
      call f_err_close_try()

      if (ierr == 0) then
        dict_tmp => dict // key
        !No errors: we have all information in astruct and put into dict
        call astruct_merge_to_dict(dict_tmp, astruct, astruct%rxyz)
        call set(dict_tmp // ASTRUCT_PROPERTIES // POSINP_SOURCE, astruct%source)

        if (GOUT_FORCES .in. dict_tmp) call dict_remove(dict_tmp, GOUT_FORCES)
        if (associated(fxyz)) then
          pos => dict_tmp // GOUT_FORCES
          do iat=1,astruct%nat
            call add(pos, dict_new(astruct%atomnames(astruct%iatype(iat)) .is. fxyz(:,iat)))
          end do
        end if

        if (GOUT_ENERGY .in. dict_tmp) call dict_remove(dict_tmp, GOUT_ENERGY)
        if (energy /= UNINITIALIZED(energy)) call set(dict_tmp // GOUT_ENERGY, energy)
        !call global_output_merge_to_dict(dict // key, outs, astruct)
        call deallocate_atomic_structure(astruct)

      else if (ierr == BIGDFT_INPUT_FILE_ERROR) then
        !Found no file, raise an error
        call f_strcpy(src='input',dest=radical)
        !modify the radical name if it exists
        call dict_get_run_properties(dict, input_id = radical)
        msg = "No section 'posinp' for the atomic positions in the file '"//&
              & trim(radical) // ".yaml'. " // trim(msg)
        call f_err_throw(err_msg=msg,err_id=ierr)
      else
        ! Raise an error
        call f_err_throw(err_msg=msg,err_id=ierr)
      end if

      call f_free_ptr(fxyz)
      call f_release_routine()

    end subroutine astruct_file_merge_to_dict


    !> Allocate the astruct variable from the dictionary of input data
    !! retrieve also other information like the energy and the forces if requested
    !! and presend in the dictionary
    subroutine astruct_set_from_dict(dict, astruct, comment)
      use module_base, only: bigdft_mpi
      use module_defs, only: gp,  UNINITIALIZED
      use numerics, only: Bohr_Ang
      use dynamic_memory
      use dictionaries
      use yaml_output, only: yaml_warning
      implicit none
      !Arguments
      type(dictionary), pointer :: dict !< dictionary of the input variables
      !! the keys have to be declared like input_dicts module
      type(atomic_structure), intent(out) :: astruct          !< Structure created from the file
      character(len = 1024), intent(out), optional :: comment !< Extra comment retrieved from the file if present
      !local variables
      character(len=*), parameter :: subname='astruct_set_from_dict'
      type(dictionary), pointer :: pos, at, types
      character(len = max_field_length) :: str
      integer :: iat, ityp, units, igspin, igchrg, ntyp

      call f_routine(id='astruct_set_from_dict')

      call nullify_atomic_structure(astruct)
      astruct%nat = -1
      if (present(comment)) write(comment, "(A)") " "

      ! The units
      units = 0
      write(astruct%units, "(A)") "bohr"
      if (has_key(dict, ASTRUCT_UNITS)) astruct%units = dict // ASTRUCT_UNITS
      select case(trim(astruct%units))
      case('atomic','atomicd0','bohr','bohrd0')
         units = 0
      case('angstroem','angstroemd0')
         units = 1
      case('reduced')
         units = 2
      end select
      ! The cell
      astruct%cell_dim = 0.0_gp
      if (.not. has_key(dict, ASTRUCT_CELL)) then
         astruct%geocode = 'F'
      else
         astruct%geocode = 'P'
         ! z
         astruct%cell_dim(3) = dict // ASTRUCT_CELL // 2
         ! y
         str = dict // ASTRUCT_CELL // 1
         if (trim(str) == ".inf") then
            astruct%geocode = 'S'
         else
            astruct%cell_dim(2) = dict // ASTRUCT_CELL // 1
         end if
         ! x
         str = dict // ASTRUCT_CELL // 0
         if (trim(str) == ".inf") then
            astruct%geocode = 'W'
         else
            astruct%cell_dim(1) = dict // ASTRUCT_CELL // 0
         end if
      end if
      if (units == 1) astruct%cell_dim = astruct%cell_dim / Bohr_Ang
      ! The types
      call astruct_dict_get_types(dict, types)
      ntyp = max(dict_size(types),0) !if types is nullified ntyp=-1
      call astruct_set_n_types(astruct, ntyp)
      ! astruct%atomnames = dict_keys(types)
      ityp = 1
      at => dict_iter(types)
      do while (associated(at))
         astruct%atomnames(ityp) = trim(dict_key(at))
         ityp = ityp + 1
         at => dict_next(at)
      end do
      ! The atoms
      if (ASTRUCT_POSITIONS .in. dict) then
         pos => dict // ASTRUCT_POSITIONS
         call astruct_set_n_atoms(astruct, dict_len(pos))
         at => dict_iter(pos)
         do while(associated(at))
            iat = dict_item(at) + 1

            call dict_copy(astruct%attributes(iat)%d, at)
            call astruct_at_from_dict(at, str, rxyz_add = astruct%rxyz(1, iat), &
                 & ifrztyp = astruct%ifrztyp(iat), igspin = igspin, igchrg = igchrg, &
                 & ixyz_add = astruct%ixyz_int(1,iat), rxyz_int_add = astruct%rxyz_int(1,iat))
            astruct%iatype(iat) = types // str
            astruct%input_polarization(iat) = 1000 * igchrg + sign(1, igchrg) * 100 + igspin

            if (units == 1) then
               astruct%rxyz(1,iat) = astruct%rxyz(1,iat) / Bohr_Ang
               astruct%rxyz(2,iat) = astruct%rxyz(2,iat) / Bohr_Ang
               astruct%rxyz(3,iat) = astruct%rxyz(3,iat) / Bohr_Ang
            endif
            if (units == 2) then !add treatment for reduced coordinates
               if (astruct%cell_dim(1) > 0.) astruct%rxyz(1,iat)=&
                    modulo(astruct%rxyz(1,iat),1.0_gp) * astruct%cell_dim(1)
               if (astruct%cell_dim(2) > 0.) astruct%rxyz(2,iat)=&
                    modulo(astruct%rxyz(2,iat),1.0_gp) * astruct%cell_dim(2)
               if (astruct%cell_dim(3) > 0.) astruct%rxyz(3,iat)=&
                    modulo(astruct%rxyz(3,iat),1.0_gp) * astruct%cell_dim(3)
            else if (astruct%geocode == 'P') then
               astruct%rxyz(1,iat)=modulo(astruct%rxyz(1,iat),astruct%cell_dim(1))
               astruct%rxyz(2,iat)=modulo(astruct%rxyz(2,iat),astruct%cell_dim(2))
               astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
            else if (astruct%geocode == 'S') then
               astruct%rxyz(1,iat)=modulo(astruct%rxyz(1,iat),astruct%cell_dim(1))
               astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
            else if (astruct%geocode == 'W') then
               astruct%rxyz(3,iat)=modulo(astruct%rxyz(3,iat),astruct%cell_dim(3))
            end if

            ! Copy attributes, except the coordinates.
            call dict_copy(astruct%attributes(iat)%d, at)
            call dict_remove(astruct%attributes(iat)%d, &
                 & trim(astruct%atomnames(astruct%iatype(iat))))

            at => dict_next(at)
         end do
      else
         call astruct_set_n_atoms(astruct,0)
      end if

      if (ASTRUCT_PROPERTIES .in. dict) then
         pos => dict // ASTRUCT_PROPERTIES
         if (("info" .in. pos) .and. present(comment)) comment = pos // "info"
         if ("format" .in. pos) then
            astruct%inputfile_format = pos // "format"
         else
            if (bigdft_mpi%iproc==0) &
                 call yaml_warning('Format not specified in the posinp dictionary, assuming yaml')
            astruct%inputfile_format='yaml'
         end if
         if ("source" .in. pos) astruct%source = pos // "source"
         call dict_copy(astruct%properties, pos)
      else
         if (bigdft_mpi%iproc==0) &
              call yaml_warning('Format not specified in the posinp dictionary, assuming yaml')
         astruct%inputfile_format='yaml'
      end if

      call dict_free(types)

      call f_release_routine()

    end subroutine astruct_set_from_dict

    include 'astruct-inc.f90'


    !> Terminate the allocation of the memory in the pointers of atoms
    subroutine allocate_atoms_data(atoms)
      implicit none
      type(atoms_data), intent(inout) :: atoms
      external :: allocate_atoms_nat,allocate_atoms_ntypes

      call allocate_atoms_nat(atoms)
      call allocate_atoms_ntypes(atoms)
    end subroutine allocate_atoms_data


    !> Fill up the atoms structure from dict
    subroutine psp_dict_analyse(dict, atoms)
      use module_defs, only: gp
      !use m_pawrad, only: pawrad_type !, pawrad_nullify
      !use m_pawtab, only: pawtab_type, pawtab_nullify
      use public_enums, only: PSPCODE_PAW
      !use public_keys, only: SOURCE_KEY
      use dynamic_memory
      use dictionaries
!      use m_libpaw_libxc, only: libxc_functionals_init, libxc_functionals_end
      use pseudopotentials, only: get_psp
      implicit none
      !Arguments
      type(dictionary), pointer :: dict        !< Input dictionary
      type(atoms_data), intent(inout) :: atoms !< Atoms structure to fill up
      !Local variables
      integer :: ityp
      character(len = 27) :: filename
      real(gp), dimension(3) :: radii_cf
      real(gp), dimension(0:4,0:6) :: psppar
      logical :: pawpatch
      integer :: paw_tot_l,  paw_tot_q, paw_tot_coefficients, paw_tot_matrices

      call f_routine(id='psp_dict_analyse')

      if (.not. associated(atoms%nzatom)) then
         call allocate_atoms_data(atoms)
      end if

      pawpatch = .true.
      do ityp=1,atoms%astruct%ntypes
         filename = 'psppar.'//atoms%astruct%atomnames(ityp)
         psppar(0:4,0:6) = atoms%psppar(0:4,0:6,ityp)
         radii_cf(:) = atoms%radii_cf(ityp,:) !To avoid a copy
         call get_psp(dict//filename,ityp,atoms%astruct%ntypes,&
              atoms%nzatom(ityp), atoms%nelpsp(ityp), atoms%npspcode(ityp), &
              atoms%ixcpsp(ityp), atoms%iradii_source(ityp),psppar,&
              radii_cf,pawpatch,&
              atoms%pawrad,atoms%pawtab,atoms%epsatm)
         atoms%radii_cf(ityp,1:3) = radii_cf(1:3)
         atoms%psppar(0:4,0:6,ityp) = psppar(0:4,0:6)
      end do
      call nlcc_set_from_dict(dict, atoms)

      !For PAW psp
      if (pawpatch.and. any(atoms%npspcode /= PSPCODE_PAW)) then
         paw_tot_l=0
         paw_tot_q=0
         paw_tot_coefficients=0
         paw_tot_matrices=0
         do ityp=1,atoms%astruct%ntypes
            filename = 'psppar.'//atoms%astruct%atomnames(ityp)
            call pawpatch_from_file( filename, atoms,ityp,&
                 paw_tot_l,  paw_tot_q, paw_tot_coefficients, paw_tot_matrices, .false.)
         end do
         do ityp=1,atoms%astruct%ntypes
            filename = 'psppar.'//atoms%astruct%atomnames(ityp)
            !! second time allocate and then store
            call pawpatch_from_file( filename, atoms,ityp,&
                 paw_tot_l, paw_tot_q, paw_tot_coefficients, paw_tot_matrices, .true.)
         end do
      else
         nullify(atoms%paw_l,atoms%paw_NofL,atoms%paw_nofchannels)
         nullify(atoms%paw_nofgaussians,atoms%paw_Greal,atoms%paw_Gimag)
         nullify(atoms%paw_Gcoeffs,atoms%paw_H_matrices,atoms%paw_S_matrices,atoms%paw_Sm1_matrices)
      end if

      call f_release_routine()

    end subroutine psp_dict_analyse


    subroutine nlcc_set_from_dict(dict, atoms)
      use module_defs, only: gp, UNINITIALIZED
      use dynamic_memory
      use dictionaries
      implicit none
      type(dictionary), pointer :: dict
      type(atoms_data), intent(inout) :: atoms
      !local variables
      type(dictionary), pointer :: nloc, coeffs
      integer :: ityp, nlcc_dim, n, i
      character(len=27) :: filename
      intrinsic :: int

      nlcc_dim = 0
      do ityp = 1, atoms%astruct%ntypes, 1
         atoms%nlcc_ngc(ityp)=UNINITIALIZED(atoms%nlcc_ngc(ityp))
         atoms%nlcc_ngv(ityp)=UNINITIALIZED(atoms%nlcc_ngv(ityp))
         filename = 'psppar.' // trim(atoms%astruct%atomnames(ityp))
         if (.not. has_key(dict, filename)) cycle
         if (.not. has_key(dict // filename, 'Non Linear Core Correction term')) cycle
         nloc => dict // filename // 'Non Linear Core Correction term'
         if (has_key(nloc, "Valence") .or. has_key(nloc, "Core")) then
            n = 0
            if (has_key(nloc, "Valence")) n = dict_len(nloc // "Valence")
            nlcc_dim = nlcc_dim + n
            atoms%nlcc_ngv(ityp) = int((sqrt(real(1 + 8 * n)) - 1) / 2)
            n = 0
            if (has_key(nloc, "Core")) n = dict_len(nloc // "Core")
            nlcc_dim = nlcc_dim + n
            atoms%nlcc_ngc(ityp) = int((sqrt(real(1 + 8 * n)) - 1) / 2)
         end if
         if (has_key(nloc, "Rcore") .and. has_key(nloc, "Core charge")) then
            nlcc_dim=nlcc_dim+1
            atoms%nlcc_ngc(ityp)=1
            atoms%nlcc_ngv(ityp)=0
         end if
      end do
      atoms%donlcc = (nlcc_dim > 0)
      atoms%nlccpar = f_malloc_ptr((/ 0 .to. 4, 1 .to. max(nlcc_dim,1) /), id = "nlccpar")
      !start again the file inspection to fill nlcc parameters
      if (atoms%donlcc) then
         nlcc_dim=0
         fill_nlcc: do ityp=1,atoms%astruct%ntypes
            filename = 'psppar.' // trim(atoms%astruct%atomnames(ityp))
            !ALEX: These are preferably read from psppar.Xy, as stored in the
            !local variables rcore and qcore
            nloc => dict // filename // 'Non Linear Core Correction term'
            if (has_key(nloc, "Valence") .or. has_key(nloc, "Core")) then
               n = 0
               if (has_key(nloc, "Valence")) n = dict_len(nloc // "Valence")
               do i = 1, n, 1
                  coeffs => nloc // "Valence" // (i - 1)
                  atoms%nlccpar(:, nlcc_dim + i) = coeffs
               end do
               nlcc_dim = nlcc_dim + n
               n = 0
               if (has_key(nloc, "Core")) n = dict_len(nloc // "Core")
               do i = 1, n, 1
                  coeffs => nloc // "Core" // (i - 1)
                  atoms%nlccpar(0, nlcc_dim + i) = coeffs // 0
                  atoms%nlccpar(1, nlcc_dim + i) = coeffs // 1
                  atoms%nlccpar(2, nlcc_dim + i) = coeffs // 2
                  atoms%nlccpar(3, nlcc_dim + i) = coeffs // 3
                  atoms%nlccpar(4, nlcc_dim + i) = coeffs // 4
               end do
               nlcc_dim = nlcc_dim + n
            end if
            if (has_key(nloc, "Rcore") .and. has_key(nloc, "Core charge")) then
               nlcc_dim=nlcc_dim+1
               atoms%nlcc_ngc(ityp)=1
               atoms%nlcc_ngv(ityp)=0
               atoms%nlccpar(0,nlcc_dim)=nloc // "Rcore"
               atoms%nlccpar(1,nlcc_dim)=nloc // "Core charge"
               atoms%nlccpar(2:4,nlcc_dim)=0.0_gp
            end if
         end do fill_nlcc
      end if
    end subroutine nlcc_set_from_dict

    !identify and inform about defined constraints on the system
    subroutine astruct_constraints(astruct)
      use dictionaries
      use yaml_output
      use numerics, only: Bohr_Ang
      use module_base, only: bigdft_mpi
      implicit none
      type(atomic_structure), intent(in) :: astruct
      !local variables
      !parameters to avoid typos in the keys while writing
      character(len=*), parameter :: CONSTRAINTS='constraints'
      character(len=*), parameter :: BONDS='bonds'
      character(len=*), parameter :: ANGLES='angles'
      character(len=*), parameter :: DIHEDRAL='dihedral'
      character(len=*), parameter :: TOLERANCE='tolerance'
      !fortran variables to store results
      logical :: msg
      integer :: iat,jat
      real(gp) :: tol,bnd,actual_bond,factor
      type(dictionary), pointer :: dict_cst,iter
      real(gp), dimension(3) :: dxyz

      !quick return if constraints are not defined according to API
      if (CONSTRAINTS .notin. astruct%properties) return

      msg=bigdft_mpi%iproc==0
      dict_cst=>astruct%properties//CONSTRAINTS
      select case (astruct%units)
      case('angstroem','angstroemd0')
         factor=Bohr_Ang
      case default
         factor=1.0_gp
      end select

      !show how to retrieve different values if needed
      tol=1.e-6_gp !default value
      tol=dict_cst .get. TOLERANCE !override if key 'tolerance' is specified
      !now inspect the different constraints
      if (BONDS .in. dict_cst) then
         if (msg) call yaml_sequence_open('Bonds for constraints')
         !iterate on the various bonds specified
         nullify(iter)
         do while(iterating(iter,on=dict_cst//BONDS))
            !for each loop cycle the iterator points
            !to the list [iat,jat,bond]
            !retrieve the values
            iat=iter//0
            jat=iter//1
            bnd=iter//2
            !and check if the values are in agreement with
            !the positions
            dxyz=astruct%rxyz(:,iat)-astruct%rxyz(:,jat)
            actual_bond=factor*sqrt(dxyz(1)**2+dxyz(2)**2+dxyz(3)**2)
            !inform the user about the bonds
            if (msg) then
            call yaml_sequence(advance='no')
             call yaml_mapping_open('Bond constraint')
              call yaml_map('Atoms',[iat,jat])
              call yaml_map('Constraint, actual bond',[bnd,actual_bond])
             call yaml_mapping_close()
            end if
         end do
         if (msg) call yaml_sequence_close()
      end if
      !similar information can be retrieved for other options and
      !constraints
      !...
    end subroutine astruct_constraints



END MODULE module_atoms

!> Allocation of the arrays inside the structure atoms_data, considering the part which is associated to astruct%nat
!! this routine is external to the module as it has to be called from C
subroutine astruct_set_n_atoms(astruct, nat)
  use module_base
  use module_atoms, only: atomic_structure
  implicit none
  type(atomic_structure), intent(inout) :: astruct
  integer, intent(in) :: nat
  !local variables
  character(len=*), parameter :: subname='astruct_set_n_atoms' !<remove
  integer :: iat

  astruct%nat = nat

  ! Allocate geometry related stuff.
  astruct%iatype = f_malloc_ptr(astruct%nat,id='astruct%iatype')
  astruct%ifrztyp = f_malloc_ptr(astruct%nat,id='astruct%ifrztyp')
  astruct%input_polarization = f_malloc_ptr(astruct%nat,id='astruct%input_polarization')
  astruct%rxyz = f_malloc0_ptr((/ 3,astruct%nat /),id='astruct%rxyz')
  astruct%rxyz_int = f_malloc_ptr((/ 3,astruct%nat /),id='astruct%rxyz_int')
  astruct%ixyz_int = f_malloc_ptr((/ 3,astruct%nat /),id='astruct%ixyz_int')

  allocate(astruct%attributes(astruct%nat))
  do iat = 1, astruct%nat
     nullify(astruct%attributes(iat)%d)
  end do

  !this array is useful for frozen atoms, no atom is frozen by default
  astruct%ifrztyp(:)=0
  !also the spin polarisation and the charge are is fixed to zero by default
  !this corresponds to the value of 100
  !RULE natpol=charge*1000 + 100 + spinpol
  astruct%input_polarization(:)=100

  !if (astruct%nat > 0) call to_zero(3 * astruct%nat, astruct%rxyz(1,1))
END SUBROUTINE astruct_set_n_atoms


!> Allocation of the memory space associated to the number of types astruct%ntypes
subroutine astruct_set_n_types(astruct, ntypes)
  use module_base
  use module_atoms, only: atomic_structure
  implicit none
  type(atomic_structure), intent(inout) :: astruct
  integer, intent(in) :: ntypes
  !character(len = *), intent(in) :: subname
  !local variables
  character(len=*), parameter :: subname='astruct_set_n_types' !<remove
  ! integer :: i_stat

  astruct%ntypes = ntypes

  ! Allocate geometry related stuff.
  astruct%atomnames=f_malloc0_str_ptr(len(astruct%atomnames),astruct%ntypes,&
       id='atomnames')
!!$  allocate(astruct%atomnames(astruct%ntypes),stat=i_stat)
!!$  call memocc(i_stat,astruct%atomnames,'astruct%atomnames',subname)

!!$  do i = 1, astruct%ntypes, 1
!!$     write(astruct%atomnames(i), "(A)") " "
!!$  end do
END SUBROUTINE astruct_set_n_types


!> Initialize the astruct variable from a file
!! Call the routine set_astruct_from_file with few arguments
subroutine astruct_set_from_file(lstat, astruct, filename)
  use module_base
  use module_atoms, only: atomic_structure,read_atomic_file=>set_astruct_from_file
  implicit none
  !Arguments
  logical, intent(out) :: lstat                     !< Error state
  type(atomic_structure), intent(inout) :: astruct  !< atomic positions structure
  character(len = *), intent(in) :: filename        !< Name of the file

  call f_err_open_try()
  call read_atomic_file(filename, 0, astruct)
  call f_err_close_try()
  lstat = (f_err_pop(BIGDFT_INPUT_VARIABLES_ERROR) == 0)

END SUBROUTINE astruct_set_from_file


!> Calculate the symmetries and update
subroutine astruct_set_symmetries(astruct, disableSym, tol, elecfield, nspin)
  use module_base
  use module_atoms, only: atomic_structure,deallocate_symmetry_data
!  use abi_defs_basis
  use m_ab6_symmetry
  implicit none
  type(atomic_structure), intent(inout) :: astruct
  logical, intent(in) :: disableSym
  real(gp), intent(in) :: tol
  real(gp), intent(in) :: elecfield(3)
  integer, intent(in) :: nspin
  !local variables
  character(len=*), parameter :: subname='astruct_set_symmetries'
  integer :: ierr
  real(gp), dimension(3,3) :: rprimd
  real(gp), dimension(:,:), allocatable :: xRed
  integer, dimension(3, 3, AB6_MAX_SYMMETRIES) :: sym
  integer, dimension(AB6_MAX_SYMMETRIES) :: symAfm
  real(gp), dimension(3, AB6_MAX_SYMMETRIES) :: transNon
  real(gp), dimension(3) :: genAfm
  integer :: spaceGroupId, pointGroupMagn

  ! Calculate the symmetries, if needed (for periodic systems only)
  if (astruct%geocode /= 'F') then
     if (astruct%sym%symObj < 0) then
        call symmetry_new(astruct%sym%symObj)
     end if
     ! Adjust tolerance
     if (tol > 0._gp) call symmetry_set_tolerance(astruct%sym%symObj, tol, ierr)
     ! New values
     rprimd(:,:) = 0.0_gp
     rprimd(1,1) = astruct%cell_dim(1)
     rprimd(2,2) = astruct%cell_dim(2)
     if (astruct%geocode == 'S') rprimd(2,2) = 1000._gp
     rprimd(3,3) = astruct%cell_dim(3)
     call symmetry_set_lattice(astruct%sym%symObj, rprimd, ierr)
     xRed = f_malloc((/ 3 , astruct%nat /),id='xRed')
     xRed(1,:) = modulo(astruct%rxyz(1, :) / rprimd(1,1), 1._gp)
     xRed(2,:) = modulo(astruct%rxyz(2, :) / rprimd(2,2), 1._gp)
     xRed(3,:) = modulo(astruct%rxyz(3, :) / rprimd(3,3), 1._gp)
     call symmetry_set_structure(astruct%sym%symObj, astruct%nat, astruct%iatype, xRed, ierr)
     call f_free(xRed)
     if (astruct%geocode == 'S') then
        call symmetry_set_periodicity(astruct%sym%symObj, &
             & (/ .true., .false., .true. /), ierr)
     else if (astruct%geocode == 'F') then
        call symmetry_set_periodicity(astruct%sym%symObj, &
             & (/ .false., .false., .false. /), ierr)
     end if
     !if (all(in%elecfield(:) /= 0)) then
     !     ! I'm not sure what this subroutine does!
     !   call symmetry_set_field(astruct%sym%symObj, (/ in%elecfield(1) , in%elecfield(2),in%elecfield(3) /), ierr)
     !elseif (in%elecfield(2) /= 0) then
     !   call symmetry_set_field(astruct%sym%symObj, (/ 0._gp, in%elecfield(2), 0._gp /), ierr)
     if (elecfield(2) /= 0) then
        call symmetry_set_field(astruct%sym%symObj, (/ 0._gp, elecfield(2), 0._gp /), ierr)
     end if
     if (nspin == 2) then
        call symmetry_set_collinear_spin(astruct%sym%symObj, astruct%nat, &
             & astruct%input_polarization, ierr)
!!$     else if (in%nspin == 4) then
!!$        call symmetry_set_spin(atoms%astruct%sym%symObj, atoms%astruct%nat, &
!!$             & atoms%astruct%input_polarization, ierror)
     end if
     if (disableSym) then
        call symmetry_set_n_sym(astruct%sym%symObj, 1, &
          & reshape((/ 1, 0, 0, 0, 1, 0, 0, 0, 1 /), (/ 3 ,3, 1 /)), &
          & reshape((/ 0.d0, 0.d0, 0.d0 /), (/ 3, 1/)), (/ 1 /), ierr)
     end if
  else
     call deallocate_symmetry_data(astruct%sym)
     astruct%sym%symObj = -1
  end if

  ! Generate symmetries for atoms
  if (.not. disableSym) then
     call symmetry_get_matrices(astruct%sym%symObj, astruct%sym%nSym, sym, transNon, symAfm, ierr)
     call symmetry_get_group(astruct%sym%symObj, astruct%sym%spaceGroup, &
          & spaceGroupId, pointGroupMagn, genAfm, ierr)
!     if (ierr == AB7_ERROR_SYM_NOT_PRIMITIVE) write(astruct%sym%spaceGroup, "(A)") "not prim."
  else
     astruct%sym%nSym = 0
     astruct%sym%spaceGroup = 'disabled'
  end if

END SUBROUTINE astruct_set_symmetries


!> Allocation of the arrays inside the structure atoms_data
subroutine allocate_atoms_nat(atoms)
  use module_base
  use module_atoms, only: atoms_data
  use ao_inguess, only : aoig_data_null,lmax_ao
  use f_arrays
  implicit none
  type(atoms_data), intent(inout) :: atoms
  integer :: iat

  atoms%refcnt=f_ref_new('atoms')

  allocate(atoms%aoig(atoms%astruct%nat))
  do iat=1,atoms%astruct%nat
     atoms%aoig(iat)=aoig_data_null()
  end do

  atoms%dogamma=f_malloc_ptr([0.to.lmax_ao,1.to.atoms%astruct%nat],&
       id='dogamma')

  !put also spin in the allocations
  !@todo here we have to decide if we need occupancy control in the
  !noncollinear case
  atoms%gamma_targets=f_malloc_ptr([0.to.lmax_ao,1.to.2,1.to.atoms%astruct%nat],id='gamma_targets')

END SUBROUTINE allocate_atoms_nat


!> Allocate the arrays in atoms_data structure.
subroutine allocate_atoms_ntypes(atoms)
  use module_base
  use module_atoms, only: atoms_data
  implicit none
  type(atoms_data), intent(inout) :: atoms
  !local variables
  character(len = *), parameter :: subname='allocate_atoms_ntypes'

  ! Allocate pseudo related stuff.
  ! store PSP parameters, modified to accept both GTH and HGHs pseudopotential types
  atoms%amu = f_malloc_ptr(atoms%astruct%ntypes,id='atoms%amu')
  atoms%psppar = f_malloc_ptr((/ 0.to.4 , 0.to.6 , 1.to.atoms%astruct%ntypes /),id='atoms%psppar')
  atoms%nelpsp = f_malloc_ptr(atoms%astruct%ntypes,id='atoms%nelpsp')
  atoms%npspcode = f_malloc_ptr(atoms%astruct%ntypes,id='atoms%npspcode')
  atoms%nzatom = f_malloc_ptr(atoms%astruct%ntypes,id='atoms%nzatom')
  atoms%ixcpsp = f_malloc_ptr(atoms%astruct%ntypes,id='atoms%ixcpsp')
  atoms%radii_cf = f_malloc_ptr((/ atoms%astruct%ntypes , 3 /),id='atoms%radii_cf')
  atoms%iradii_source = f_malloc_ptr(atoms%astruct%ntypes,id='atoms%iradii_source')
  ! parameters for NLCC
  atoms%nlcc_ngv = f_malloc_ptr(atoms%astruct%ntypes,id='atoms%nlcc_ngv')
  atoms%nlcc_ngc = f_malloc_ptr(atoms%astruct%ntypes,id='atoms%nlcc_ngc')
END SUBROUTINE allocate_atoms_ntypes


! Init and free routines


!> Allocate a new atoms_data type, for bindings.
subroutine atoms_new(atoms)
  use module_atoms, only: atoms_data,nullify_atoms_data
  use dynamic_memory
  implicit none
  type(atoms_data), pointer :: atoms
  type(atoms_data), pointer, save :: intern

  allocate(intern)
  call nullify_atoms_data(intern)! = atoms_data_null()
  atoms => intern
END SUBROUTINE atoms_new


!> Free an allocated atoms_data type.
subroutine atoms_free(atoms)
  use module_atoms, only: atoms_data,deallocate_atoms_data
  use f_refcnts, only: f_ref_count, f_ref_new
  implicit none
  type(atoms_data), pointer :: atoms

  if (f_ref_count(atoms%refcnt) < 0) then
     ! Trick here to be sure that the deallocate won't complain in case of not
     ! fully initialised atoms.
     atoms%refcnt=f_ref_new('atoms')
  end if
  call deallocate_atoms_data(atoms)
  deallocate(atoms)
END SUBROUTINE atoms_free


!> Add a displacement of atomic positions and put in the box
subroutine astruct_set_displacement(astruct, randdis)
  use module_defs, only: gp
  use module_atoms, only: rxyz_inside_box,atomic_structure,&
       move_this_coordinate
  implicit none
  type(atomic_structure), intent(inout) :: astruct
  real(gp), intent(in) :: randdis !< random displacement

  integer :: iat,i
  real :: tt

  !Shake atoms if required.
  if (randdis > 0.d0) then
     do iat=1,astruct%nat
        do i=1,3
           if (move_this_coordinate(astruct%ifrztyp(iat),i)) then
              call random_number(tt)
              astruct%rxyz(i,iat)=astruct%rxyz(i,iat)+randdis*real(tt,gp)
            end if
         end do
      end do
   end if

   call rxyz_inside_box(astruct)

 END SUBROUTINE astruct_set_displacement

 !> this routine should be merged with the cell definition in the box module
 subroutine astruct_distance(astruct, rxyz, dxyz, iat1, iat2)
   use module_defs, only: gp
   use module_atoms, only: atomic_structure
   implicit none
   type(atomic_structure), intent(in) :: astruct
   real(gp), dimension(3, astruct%nat), intent(in) :: rxyz
   real(gp), dimension(3), intent(out) :: dxyz
   integer, intent(in) :: iat1, iat2

   logical, dimension(3) :: per

   select case(astruct%geocode)
   case ("P")
      per = (/ .true., .true., .true. /)
   case ("S")
      per = (/ .true., .false., .true. /)
   case ("W")
      per = (/ .false., .true., .false. /)
   case default
      per = (/ .false., .false., .false. /)
   end select

   dxyz(:) = rxyz(:, iat2) - rxyz(:, iat1)
   where (per) dxyz = dxyz - astruct%cell_dim * nint(dxyz / astruct%cell_dim)
 END SUBROUTINE astruct_distance

!> Compute a list of neighbours for the given structure.
subroutine astruct_neighbours(astruct, rxyz, neighb)
  use module_defs, only: gp
  use module_atoms, only: atomic_structure, atomic_neighbours, nullify_atomic_neighbours
  use dynamic_memory
  use ao_inguess, only: atomic_info,atomic_z
  implicit none
  type(atomic_structure), intent(in) :: astruct
  real(gp), dimension(3, astruct%nat), intent(in) :: rxyz
  type(atomic_neighbours), intent(out) :: neighb

  integer :: maxnei, i, j, nnei
  integer, dimension(:,:), allocatable :: tmp_nei
  logical, dimension(3) :: per
  real(gp), dimension(3) :: dxyz
  real(gp), dimension(:), allocatable :: rcuts

  call nullify_atomic_neighbours(neighb)

  neighb%nat = astruct%nat
  neighb%keynei = f_malloc0_ptr((/ 2, neighb%nat /), id = "neighb%keynei")

  maxnei = min(astruct%nat, 50)
  tmp_nei = f_malloc((/ maxnei, astruct%nat /), id = "tmp_nei")

  select case(astruct%geocode)
  case ("P")
     per = (/ .true., .true., .true. /)
  case ("S")
     per = (/ .true., .false., .true. /)
  case ("W")
     per = (/ .false., .true., .false. /)
  case default
     per = (/ .false., .false., .false. /)
  end select

  rcuts = f_malloc(astruct%ntypes, id = "rcuts")
  do i = 1, astruct%ntypes, 1
     call atomic_info(atomic_z(trim(astruct%atomnames(i))), rcov = rcuts(i))
     rcuts(i) = rcuts(i) * 1.2_gp ! add 20% in case.
  end do

  nnei = 0
  do i = 1, astruct%nat
     do j = i + 1, astruct%nat

        dxyz(:) = rxyz(:, j) - rxyz(:, i)
        where (per) dxyz = dxyz - astruct%cell_dim * nint(dxyz / astruct%cell_dim)

        if (dxyz(1) * dxyz(1) + dxyz(2) * dxyz(2) + dxyz(3) * dxyz(3) < &
             & (rcuts(astruct%iatype(i)) + rcuts(astruct%iatype(j)))**2) then
           neighb%keynei(1, i) = neighb%keynei(1, i) + 1
           neighb%keynei(1, j) = neighb%keynei(1, j) + 1
           tmp_nei(neighb%keynei(1, i), i) = j
           tmp_nei(neighb%keynei(1, j), j) = i
        endif

     end do
     nnei = nnei + neighb%keynei(1, i)
  end do

  call f_free(rcuts)

  neighb%nei = f_malloc_ptr(nnei, id = "neighb%nei")

  nnei = 1
  do i = 1, neighb%nat
     neighb%keynei(2, i) = nnei
     if (neighb%keynei(1, i) > 0) &
          & neighb%nei(nnei: nnei + neighb%keynei(1, i) - 1) = &
          & tmp_nei(1:neighb%keynei(1, i), i)
     nnei = nnei + neighb%keynei(1, i)
  end do

  call f_free(tmp_nei)
END SUBROUTINE astruct_neighbours

subroutine astruct_from_subset(asub, astruct, rxyz, mask, passivate, bufNeighbours, frozen)
  use module_defs, only: gp
  use module_atoms, only: atomic_structure, nullify_atomic_structure, &
       & atomic_neighbours, astruct_neighbours_iter, astruct_neighbours_next, &
       & deallocate_atomic_neighbours
  use dynamic_memory
  use dictionaries
  use public_keys, only: ASTRUCT_ATT_ORIG_ID, ASTRUCT_ATT_FROZEN
  use ao_inguess, only: atomic_z, atomic_info
  implicit none
  type(atomic_structure), intent(out) :: asub
  type(atomic_structure), intent(in) :: astruct
  real(gp), dimension(3, astruct%nat), intent(in) :: rxyz
  logical, dimension(astruct%nat), intent(in) :: mask !< .true. for atoms in the subset.
  logical, intent(in) :: passivate
  integer, intent(in) :: bufNeighbours !< add a buffer layer made from the nth neighbours
  character(len = *), intent(in) :: frozen !< frozen attribute to apply to passivation and buffer layer.

  type(atomic_neighbours) :: nei
  integer :: i, iat, jat, nsub
  real(gp) :: rcutH, fact
  type(dictionary), pointer :: hlist, s, types
  logical, dimension(3) :: per
  real(gp), dimension(3) :: dxyz
  real(gp), dimension(:), allocatable :: rcuts
  logical, dimension(:), allocatable :: work
  logical, dimension(:), allocatable :: mask_

  call nullify_atomic_structure(asub)

  if (passivate .or. bufNeighbours > 0) then
     call astruct_neighbours(astruct, rxyz, nei)
  end if

  mask_ = f_malloc(astruct%nat, id = "mask_")
  mask_ = mask
  if (bufNeighbours > 0) then
     work = f_malloc(astruct%nat, id = "nmask")
     do i = 1, bufNeighbours
        call expand(astruct%nat, mask_, work)
     end do
     call f_free(work)
  end if

  call dict_init(hlist)
  if (passivate) then
     ! In case of passivation, every old neighbours that are cut, are replaced
     ! by an hydrogen.

     select case(astruct%geocode)
     case ("P")
        per = (/ .true., .true., .true. /)
     case ("S")
        per = (/ .true., .false., .true. /)
     case ("W")
        per = (/ .false., .true., .false. /)
     case default
        per = (/ .false., .false., .false. /)
     end select

     rcuts = f_malloc(astruct%ntypes, id = "rcuts")
     do i = 1, astruct%ntypes, 1
        call atomic_info(atomic_z(trim(astruct%atomnames(i))), rcov = rcuts(i))
     end do
     call atomic_info(1, rcov = rcutH)

     nullify(s)
     do iat = 1, astruct%nat
        if (mask_(iat)) then
           call astruct_neighbours_iter(nei, iat)
           do while(astruct_neighbours_next(nei, jat))
              if ( .not. mask_(jat)) then
                 dxyz(:) = rxyz(:, jat) - rxyz(:, iat)
                 where (per) dxyz = dxyz - astruct%cell_dim * nint(dxyz / astruct%cell_dim)
                 fact = (rcuts(astruct%iatype(iat)) + rcutH) / &
                      & sqrt(dxyz(1) * dxyz(1) + dxyz(2) * dxyz(2) + dxyz(3) * dxyz(3))
                 dxyz(:) = rxyz(:, iat) + dxyz(:) * fact
                 call add(hlist, list_new((/ .item. dxyz(1), .item. dxyz(2), .item. dxyz(3) /)), s)
              endif
           enddo
        end if
     enddo

     call f_free(rcuts)
  end if

  if (passivate .or. bufNeighbours > 0) then
     call deallocate_atomic_neighbours(nei)
  end if

  ! Start copying a subset of astruct into asub.
  asub%units = astruct%units
  asub%cell_dim = astruct%cell_dim
  asub%geocode = astruct%geocode
  asub%inputfile_format = astruct%inputfile_format

  ! Count the number of types in the subset.
  call dict_init(types)
  do iat = 1, astruct%nat
     if (mask_(iat) .and. .not. (trim(astruct%atomnames(astruct%iatype(iat))) .in. types)) &
          & call set(types // trim(astruct%atomnames(astruct%iatype(iat))), dict_size(types))
  end do
  if (dict_len(hlist) > 0 .and. .not. ("H" .in. types)) &
       & call set(types // "H", dict_size(types))
  call astruct_set_n_types(asub, dict_size(types))
  i = 1
  s => dict_iter(types)
  do while (associated(s))
     asub%atomnames(i) = trim(dict_key(s))
     i = i + 1
     s => dict_next(s)
  end do

  ! Count the number of atoms in the subset.
  nsub = 0
  do iat = 1, astruct%nat
     if (mask_(iat)) nsub = nsub + 1
  end do
  call astruct_set_n_atoms(asub, nsub + dict_len(hlist))
  i = 0
  do iat = 1, astruct%nat
     if (mask_(iat)) then
        i = i + 1
        asub%iatype(i) = types // trim(astruct%atomnames(astruct%iatype(iat)))
        asub%input_polarization(i) = astruct%input_polarization(iat)
        asub%rxyz(:, i) = astruct%rxyz(:, iat)
        if (associated(astruct%attributes(iat)%d)) then
           call dict_copy(asub%attributes(i)%d, astruct%attributes(iat)%d)
        end if
        if (.not. associated(asub%attributes(i)%d)) then
           call dict_init(asub%attributes(i)%d)
        end if
        call set(asub%attributes(i)%d // ASTRUCT_ATT_ORIG_ID, iat)
        if (len_trim(frozen) > 0 .and. .not. mask(iat)) then
           call set(asub%attributes(i)%d // ASTRUCT_ATT_FROZEN, frozen)
        end if
     end if
  end do
  s => dict_iter(hlist)
  do while(associated(s))
     i = i + 1
     asub%iatype(i) = types // "H"
     asub%rxyz(:, i) = s
     if (len_trim(frozen) > 0) then
        call dict_init(asub%attributes(i)%d)
        call set(asub%attributes(i)%d // ASTRUCT_ATT_FROZEN, frozen)
     end if
     s => dict_next(s)
  end do

  call dict_free(types)
  call dict_free(hlist)
  call f_free(mask_)

contains

  subroutine expand(nat, mask_, nmask)
    implicit none
    integer, intent(in) :: nat
    logical, dimension(nat), intent(inout) :: mask_
    logical, dimension(nat), intent(out) :: nmask
    integer :: iat, jat

    nmask = .false.
    do iat = 1, nat
       if (mask_(iat)) then
          call astruct_neighbours_iter(nei, iat)
          do while(astruct_neighbours_next(nei, jat))
             nmask(jat) = .true.
          end do
       end if
    end do
    do iat = 1, nat
       mask_(iat) = mask_(iat) .or. nmask(iat)
    end do
  end subroutine expand
END SUBROUTINE astruct_from_subset

subroutine astruct_set_cell(dict, cell)
  use dictionaries
  implicit none
  type(dictionary), pointer :: dict
  double precision, dimension(3), intent(in) :: cell

  call set(dict // "cell", cell)
end subroutine astruct_set_cell

subroutine astruct_add_atom(dict, xyz, symbol, slen)
  use dictionaries
  implicit none
  type(dictionary), pointer :: dict
  double precision, dimension(3), intent(in) :: xyz
  integer, intent(in) :: slen
  character(len = slen), intent(in) :: symbol

  call add(dict // "positions", dict_new(symbol .is. xyz))
end subroutine astruct_add_atom

subroutine astruct_get_types_dict(dict,types)
  use dictionaries
  use module_atoms, only: astruct_dict_get_types
  type(dictionary), pointer :: dict,types
  call astruct_dict_get_types(dict,types)
end subroutine astruct_get_types_dict
