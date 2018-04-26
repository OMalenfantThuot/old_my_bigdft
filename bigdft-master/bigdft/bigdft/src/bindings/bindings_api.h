/** @file
 * Bindings for the BigDFT package
 * @author
 * Copyright (C) 2013-2015 BigDFT group (DC)
 * This file is distributed under the terms of the
 * GNU General Public License, see ~/COPYING file
 * or http://www.gnu.org/copyleft/gpl.txt .
 * For the list of contributors, see ~/AUTHORS
**/
#ifndef BINDINGS_API_H
#define BINDINGS_API_H

#undef hz

/* atoms_get_amu src/init/atoms.f90:214 */
/* Fortran header:
subroutine atoms_get_amu(atoms, amu)
use module_defs, only: gp
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
real(gp), dimension(:), pointer :: amu
*/
void FC_FUNC_(atoms_get_amu, ATOMS_GET_AMU)(const f90_atoms_data *atoms, 
                                            f90_pointer_double *amu);
/* atoms_get_iatype src/init/atoms.f90:103 */
/* Fortran header:
subroutine atoms_get_iatype(atoms, iatype)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: iatype
*/
void FC_FUNC_(atoms_get_iatype, ATOMS_GET_IATYPE)(const f90_atoms_data *atoms, 
                                                  f90_pointer_int *iatype);
/* atoms_get_ifrztyp src/init/atoms.f90:135 */
/* Fortran header:
subroutine atoms_get_ifrztyp(atoms, ifrztyp)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: ifrztyp
*/
void FC_FUNC_(atoms_get_ifrztyp, ATOMS_GET_IFRZTYP)(const f90_atoms_data *atoms, 
                                                    f90_pointer_int *ifrztyp);
/* atoms_get_ixcpsp src/init/atoms.f90:204 */
/* Fortran header:
subroutine atoms_get_ixcpsp(atoms, ixcpsp)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: ixcpsp
*/
void FC_FUNC_(atoms_get_ixcpsp, ATOMS_GET_IXCPSP)(const f90_atoms_data *atoms, 
                                                  f90_pointer_int *ixcpsp);
/* atoms_get_natpol src/init/atoms.f90:125 */
/* Fortran header:
subroutine atoms_get_natpol(atoms, natpol)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: natpol
*/
void FC_FUNC_(atoms_get_natpol, ATOMS_GET_NATPOL)(const f90_atoms_data *atoms, 
                                                  f90_pointer_int *natpol);
/* atoms_get_nelpsp src/init/atoms.f90:156 */
/* Fortran header:
subroutine atoms_get_nelpsp(atoms, nelpsp)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: nelpsp
*/
void FC_FUNC_(atoms_get_nelpsp, ATOMS_GET_NELPSP)(const f90_atoms_data *atoms, 
                                                  f90_pointer_int *nelpsp);
/* atoms_get_nlcc_ngc src/init/atoms.f90:194 */
/* Fortran header:
subroutine atoms_get_nlcc_ngc(atoms, nlcc_ngc)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: nlcc_ngc
*/
void FC_FUNC_(atoms_get_nlcc_ngc, ATOMS_GET_NLCC_NGC)(const f90_atoms_data *atoms, 
                                                      f90_pointer_int *nlcc_ngc);
/* atoms_get_nlcc_ngv src/init/atoms.f90:186 */
/* Fortran header:
subroutine atoms_get_nlcc_ngv(atoms, nlcc_ngv)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: nlcc_ngv
*/
void FC_FUNC_(atoms_get_nlcc_ngv, ATOMS_GET_NLCC_NGV)(const f90_atoms_data *atoms, 
                                                      f90_pointer_int *nlcc_ngv);
/* atoms_get_nlccpar src/init/atoms.f90:257 */
/* Fortran header:
subroutine atoms_get_nlccpar(atoms, nlccpar)
use module_defs, only: gp
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
real(gp), dimension(:,:), pointer :: nlccpar
*/
void FC_FUNC_(atoms_get_nlccpar, ATOMS_GET_NLCCPAR)(const f90_atoms_data *atoms, 
                                                    f90_pointer_double_2D *nlccpar);
/* atoms_get_npspcode src/init/atoms.f90:166 */
/* Fortran header:
subroutine atoms_get_npspcode(atoms, npspcode)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: npspcode
*/
void FC_FUNC_(atoms_get_npspcode, ATOMS_GET_NPSPCODE)(const f90_atoms_data *atoms, 
                                                      f90_pointer_int *npspcode);
/* atoms_get_nzatom src/init/atoms.f90:176 */
/* Fortran header:
subroutine atoms_get_nzatom(atoms, nzatom)
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
integer, dimension(:), pointer :: nzatom
*/
void FC_FUNC_(atoms_get_nzatom, ATOMS_GET_NZATOM)(const f90_atoms_data *atoms, 
                                                  f90_pointer_int *nzatom);
/* atoms_get_psppar src/init/atoms.f90:247 */
/* Fortran header:
subroutine atoms_get_psppar(atoms, psppar)
use module_defs, only: gp
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
real(gp), dimension(:,:,:), pointer :: psppar
*/
void FC_FUNC_(atoms_get_psppar, ATOMS_GET_PSPPAR)(const f90_atoms_data *atoms, 
                                                  f90_pointer_double_3D *psppar);
/* atoms_get_radii_cf src/init/atoms.f90:236 */
/* Fortran header:
subroutine atoms_get_radii_cf(atoms, radii_cf)
use module_defs, only: gp
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
real(gp), dimension(:,:), pointer :: radii_cf
*/
void FC_FUNC_(atoms_get_radii_cf, ATOMS_GET_RADII_CF)(const f90_atoms_data *atoms, 
                                                      f90_pointer_double_2D *radii_cf);
/* atoms_get_rxyz src/init/atoms.f90:145 */
/* Fortran header:
subroutine atoms_get_rxyz(atoms, rxyz)
use module_defs, only: gp
use module_types
implicit none
type(atoms_data), intent(in) :: atoms
real(gp), dimension(:,:), pointer :: rxyz
*/
void FC_FUNC_(atoms_get_rxyz, ATOMS_GET_RXYZ)(const f90_atoms_data *atoms, 
                                              f90_pointer_double_2D *rxyz);
/* localfields_get_rhov src/bindings/bindingsf.f90:1031 */
/* Fortran header:
subroutine localfields_get_rhov(denspot, rhov)
use module_defs, only: dp
use module_types
implicit none
type(DFT_local_fields), intent(in) :: denspot
real(dp), dimension(:), pointer :: rhov
*/
void FC_FUNC_(localfields_get_rhov, LOCALFIELDS_GET_RHOV)(const f90_DFT_local_fields *denspot, 
                                                          f90_pointer_double *rhov);
/* localfields_get_v_ext src/bindings/bindingsf.f90:1042 */
/* Fortran header:
subroutine localfields_get_v_ext(denspot, v_ext)
use module_defs, only: wp
use module_types
implicit none
type(DFT_local_fields), intent(in) :: denspot
real(wp), dimension(:,:,:,:), pointer :: v_ext
*/
void FC_FUNC_(localfields_get_v_ext, LOCALFIELDS_GET_V_EXT)(const f90_DFT_local_fields *denspot, 
                                                            f90_pointer_double_4D *v_ext);
/* localfields_get_v_xc src/bindings/bindingsf.f90:1053 */
/* Fortran header:
subroutine localfields_get_v_xc(denspot, v_xc)
use module_defs, only: wp
use module_types
implicit none
type(DFT_local_fields), intent(in) :: denspot
real(wp), dimension(:,:,:,:), pointer :: v_xc
*/
void FC_FUNC_(localfields_get_v_xc, LOCALFIELDS_GET_V_XC)(const f90_DFT_local_fields *denspot, 
                                                          f90_pointer_double_4D *v_xc);
/* orbs_get_eval src/bindings/bindingsf.f90:820 */
/* Fortran header:
subroutine orbs_get_eval(orbs, eval)
use module_defs, only: wp
use module_types
implicit none
type(orbitals_data) :: orbs
real(wp), dimension(:), pointer :: eval
*/
void FC_FUNC_(orbs_get_eval, ORBS_GET_EVAL)(f90_orbitals_data *orbs, 
                                            f90_pointer_double *eval);
/* orbs_get_inwhichlocreg src/bindings/bindingsf.f90:864 */
/* Fortran header:
subroutine orbs_get_inwhichlocreg(orbs, locreg)
use module_types
implicit none
type(orbitals_data) :: orbs
integer, dimension(:), pointer :: locreg
*/
void FC_FUNC_(orbs_get_inwhichlocreg, ORBS_GET_INWHICHLOCREG)(f90_orbitals_data *orbs, 
                                                              f90_pointer_int *locreg);
/* orbs_get_kpts src/bindings/bindingsf.f90:842 */
/* Fortran header:
subroutine orbs_get_kpts(orbs, kpts)
use module_defs, only: gp
use module_types
implicit none
type(orbitals_data) :: orbs
real(gp), dimension(:,:), pointer :: kpts
*/
void FC_FUNC_(orbs_get_kpts, ORBS_GET_KPTS)(f90_orbitals_data *orbs, 
                                            f90_pointer_double_2D *kpts);
/* orbs_get_kwgts src/bindings/bindingsf.f90:853 */
/* Fortran header:
subroutine orbs_get_kwgts(orbs, kwgts)
use module_defs, only: gp
use module_types
implicit none
type(orbitals_data) :: orbs
real(gp), dimension(:), pointer :: kwgts
*/
void FC_FUNC_(orbs_get_kwgts, ORBS_GET_KWGTS)(f90_orbitals_data *orbs, 
                                              f90_pointer_double *kwgts);
/* orbs_get_occup src/bindings/bindingsf.f90:831 */
/* Fortran header:
subroutine orbs_get_occup(orbs, occup)
use module_defs, only: gp
use module_types
implicit none
type(orbitals_data) :: orbs
real(gp), dimension(:), pointer :: occup
*/
void FC_FUNC_(orbs_get_occup, ORBS_GET_OCCUP)(f90_orbitals_data *orbs, 
                                              f90_pointer_double *occup);
/* orbs_get_onwhichatom src/bindings/bindingsf.f90:874 */
/* Fortran header:
subroutine orbs_get_onwhichatom(orbs, atom)
use module_types
implicit none
type(orbitals_data) :: orbs
integer, dimension(:), pointer :: atom
*/
void FC_FUNC_(orbs_get_onwhichatom, ORBS_GET_ONWHICHATOM)(f90_orbitals_data *orbs, 
                                                          f90_pointer_int *atom);
/* orbs_get_onwhichmpi  */
/* allocate_atoms_nat src/modules/atoms_data.f90:1955 */
/* Fortran header:
subroutine allocate_atoms_nat(atoms)
use module_base
use module_atoms, only: atoms_data
use ao_inguess, only : aoig_data_null,lmax_ao
implicit none
type(atoms_data), intent(inout) :: atoms
integer :: iat
*/
void FC_FUNC_(allocate_atoms_nat, ALLOCATE_ATOMS_NAT)(f90_atoms_data *atoms);
/* allocate_atoms_ntypes src/modules/atoms_data.f90:1976 */
/* Fortran header:
subroutine allocate_atoms_ntypes(atoms)
use module_base
use module_atoms, only: atoms_data
implicit none
type(atoms_data), intent(inout) :: atoms

character(len = *), parameter :: subname='allocate_atoms_ntypes'
*/
void FC_FUNC_(allocate_atoms_ntypes, ALLOCATE_ATOMS_NTYPES)(f90_atoms_data *atoms);
/* allocaterhopot src/init/denspotd.f90:480 */
/* Fortran header:
subroutine allocateRhoPot(Glr,nspin,atoms,rxyz,denspot)
use module_base
use module_types
use module_interfaces, only: calculate_rhocore
implicit none
integer, intent(in) :: nspin
type(locreg_descriptors), intent(in) :: Glr
type(atoms_data), intent(in) :: atoms
real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
type(DFT_local_fields), intent(inout) :: denspot
*/
void FC_FUNC(allocaterhopot, ALLOCATERHOPOT)(const f90_locreg_descriptors *Glr, 
                                             const int *nspin, 
                                             const f90_atoms_data *atoms, 
                                             const double *rxyz, 
                                             f90_DFT_local_fields *denspot);
/* astruct_copy_alat src/init/atoms.f90:342 */
/* Fortran header:
subroutine astruct_copy_alat(astruct, alat)
use module_defs, only: gp
use module_types
implicit none
type(atomic_structure), intent(in) :: astruct
real(gp), intent(out) :: alat(3)
*/
void FC_FUNC_(astruct_copy_alat, ASTRUCT_COPY_ALAT)(const f90_atomic_structure *astruct, 
                                                    double *alat);
/* astruct_copy_geometry_data src/init/atoms.f90:289 */
/* Fortran header:
subroutine astruct_copy_geometry_data(astruct, geocode, format, units)
use module_types
implicit none
type(atomic_structure), intent(in) :: astruct
character(len = 1), intent(out) :: geocode 
character(len = 5), intent(out) :: format
character(len = 20), intent(out) :: units
*/
void FC_FUNC_(astruct_copy_geometry_data, ASTRUCT_COPY_GEOMETRY_DATA)(const f90_atomic_structure *astruct, 
                                                                      char *geocode, 
                                                                      char *format, 
                                                                      char *units, 
                                                                      int str_ln_1, 
                                                                      int str_ln_2, 
                                                                      int str_ln_3);
/* astruct_copy_name src/init/atoms.f90:315 */
/* Fortran header:
subroutine astruct_copy_name(astruct, ityp, name, ln)
use module_types
implicit none

type(atomic_structure), intent(in) :: astruct
integer, intent(in) :: ityp
character(len=1), dimension(20), intent(out) :: name

integer, intent(out) :: ln

integer :: i,lname
*/
void FC_FUNC_(astruct_copy_name, ASTRUCT_COPY_NAME)(const f90_atomic_structure *astruct, 
                                                    const int *ityp, 
                                                    char *name, 
                                                    int *ln, 
                                                    int str_ln_1);
/* astruct_copy_nat src/init/atoms.f90:83 */
/* Fortran header:
subroutine astruct_copy_nat(astruct, nat)
use module_types
implicit none
type(atomic_structure), intent(in) :: astruct
integer, intent(out) :: nat
*/
void FC_FUNC_(astruct_copy_nat, ASTRUCT_COPY_NAT)(const f90_atomic_structure *astruct, 
                                                  int *nat);
/* astruct_copy_ntypes src/init/atoms.f90:93 */
/* Fortran header:
subroutine astruct_copy_ntypes(astruct, ntypes)
use module_types
implicit none
type(atomic_structure), intent(in) :: astruct
integer, intent(out) :: ntypes
*/
void FC_FUNC_(astruct_copy_ntypes, ASTRUCT_COPY_NTYPES)(const f90_atomic_structure *astruct, 
                                                        int *ntypes);
/* astruct_merge_to_dict_binding src/init/atoms.f90:277 */
/* Fortran header:
subroutine astruct_merge_to_dict_binding(dict, astruct)
use module_atoms, only: wrapper => astruct_merge_to_dict
use module_types, only: atomic_structure
use dictionaries, only: dictionary
implicit none
type(dictionary), pointer :: dict
type(atomic_structure), intent(in) :: astruct
*/
void FC_FUNC_(astruct_merge_to_dict_binding, ASTRUCT_MERGE_TO_DICT_BINDING)(f90_dictionary_pointer *dict, 
                                                                            const f90_atomic_structure *astruct);
/* astruct_set_displacement src/modules/atoms_data.f90:2035 */
/* Fortran header:
subroutine astruct_set_displacement(astruct, randdis)
use module_defs, only: gp
use module_atoms, only: rxyz_inside_box,atomic_structure,move_this_coordinate
implicit none
type(atomic_structure), intent(inout) :: astruct
real(gp), intent(in) :: randdis 

integer :: iat,i
real :: tt
*/
void FC_FUNC_(astruct_set_displacement, ASTRUCT_SET_DISPLACEMENT)(f90_atomic_structure *astruct, 
                                                                  const double *randdis);
/* astruct_set_from_dict_binding src/bindings/bindingsf.f90:1997 */
/* Fortran header:
subroutine astruct_set_from_dict_binding(astruct, dict)
use dictionaries, only: dictionary
use module_atoms
implicit none
type(dictionary), pointer :: dict 

type(atomic_structure), intent(out) :: astruct
*/
void FC_FUNC_(astruct_set_from_dict_binding, ASTRUCT_SET_FROM_DICT_BINDING)(f90_atomic_structure *astruct, 
                                                                            f90_dictionary_pointer *dict);
/* astruct_set_from_file src/modules/atoms_data.f90:1848 */
/* Fortran header:
subroutine astruct_set_from_file(lstat, astruct, filename)
use module_base
use module_atoms, only: atomic_structure,read_atomic_file=>set_astruct_from_file
implicit none

logical, intent(out) :: lstat                     
type(atomic_structure), intent(inout) :: astruct  
character(len = *), intent(in) :: filename
*/
void FC_FUNC_(astruct_set_from_file, ASTRUCT_SET_FROM_FILE)(int *lstat, 
                                                            f90_atomic_structure *astruct, 
                                                            const char *filename, 
                                                            int str_ln_1);
/* astruct_set_geometry src/init/atoms.f90:54 */
/* Fortran header:
subroutine astruct_set_geometry(astruct, alat, geocode, format, units)
use module_defs, only: gp
use module_types
implicit none
type(atomic_structure), intent(inout) :: astruct
real(gp), intent(in) :: alat(3)
character(len=1), intent(in) :: geocode 
character, intent(in) :: format(5)
character, intent(in) :: units(20)
*/
void FC_FUNC_(astruct_set_geometry, ASTRUCT_SET_GEOMETRY)(f90_atomic_structure *astruct, 
                                                          const double *alat, 
                                                          const char *geocode, 
                                                          const char *format, 
                                                          const char *units, 
                                                          int str_ln_1, 
                                                          int str_ln_2, 
                                                          int str_ln_3);
/* astruct_set_n_atoms src/modules/atoms_data.f90:1784 */
/* Fortran header:
subroutine astruct_set_n_atoms(astruct, nat)
use module_base
use module_atoms, only: atomic_structure
implicit none
type(atomic_structure), intent(inout) :: astruct
integer, intent(in) :: nat

character(len=*), parameter :: subname='astruct_set_n_atoms' 
integer :: iat
*/
void FC_FUNC_(astruct_set_n_atoms, ASTRUCT_SET_N_ATOMS)(f90_atomic_structure *astruct, 
                                                        const int *nat);
/* astruct_set_n_types src/modules/atoms_data.f90:1821 */
/* Fortran header:
subroutine astruct_set_n_types(astruct, ntypes)
use module_base
use module_atoms, only: atomic_structure
implicit none
type(atomic_structure), intent(inout) :: astruct
integer, intent(in) :: ntypes


character(len=*), parameter :: subname='astruct_set_n_types'
*/
void FC_FUNC_(astruct_set_n_types, ASTRUCT_SET_N_TYPES)(f90_atomic_structure *astruct, 
                                                        const int *ntypes);
/* astruct_set_symmetries src/modules/atoms_data.f90:1866 */
/* Fortran header:
subroutine astruct_set_symmetries(astruct, disableSym, tol, elecfield, nspin)
use module_base
use module_atoms, only: atomic_structure,deallocate_symmetry_data

use m_ab6_symmetry
implicit none
type(atomic_structure), intent(inout) :: astruct
logical, intent(in) :: disableSym
real(gp), intent(in) :: tol
real(gp), intent(in) :: elecfield(3)
integer, intent(in) :: nspin

character(len=*), parameter :: subname='astruct_set_symmetries'
integer :: ierr
real(gp), dimension(3,3) :: rprimd
real(gp), dimension(:,:), allocatable :: xRed
integer, dimension(3, 3, AB6_MAX_SYMMETRIES) :: sym
integer, dimension(AB6_MAX_SYMMETRIES) :: symAfm
real(gp), dimension(3, AB6_MAX_SYMMETRIES) :: transNon
real(gp), dimension(3) :: genAfm
integer :: spaceGroupId, pointGroupMagn
*/
void FC_FUNC_(astruct_set_symmetries, ASTRUCT_SET_SYMMETRIES)(f90_atomic_structure *astruct, 
                                                              const int *disableSym, 
                                                              const double *tol, 
                                                              const double *elecfield, 
                                                              const int *nspin);
/* atoms_copy_alat  */
/* atoms_empty src/init/atoms.f90:32 */
/* Fortran header:
subroutine atoms_empty(atoms)
use module_atoms, only: atoms_data, deallocate_atoms_data
use f_refcnts, only: f_ref_new
implicit none
type(atoms_data), intent(inout) :: atoms
*/
void FC_FUNC_(atoms_empty, ATOMS_EMPTY)(f90_atoms_data *atoms);
/* atoms_free src/modules/atoms_data.f90:2018 */
/* Fortran header:
subroutine atoms_free(atoms)
use module_atoms, only: atoms_data,deallocate_atoms_data
use f_refcnts, only: f_ref_count, f_ref_new
implicit none
type(atoms_data), pointer :: atoms
*/
void FC_FUNC_(atoms_free, ATOMS_FREE)(f90_atoms_data_pointer *atoms);
/* atoms_get src/init/atoms.f90:71 */
/* Fortran header:
subroutine atoms_get(atoms, astruct, symObj)
use module_types
implicit none
type(atoms_data), intent(in), target :: atoms
type(atomic_structure), pointer :: astruct
type(symmetry_data), pointer :: symObj
*/
void FC_FUNC_(atoms_get, ATOMS_GET)(const f90_atoms_data *atoms, 
                                    f90_atomic_structure_pointer *astruct, 
                                    f90_symmetry_data_pointer *symObj);
/* atoms_new src/modules/interface2.f90:163 src/modules/atoms_data.f90:2004 */
/* Fortran header:
subroutine atoms_new(atoms)
use module_types
implicit none
type(atoms_data), pointer :: atoms
*/
void FC_FUNC_(atoms_new, ATOMS_NEW)(f90_atoms_data_pointer *atoms);
/* atoms_set_name src/init/atoms.f90:43 */
/* Fortran header:
subroutine atoms_set_name(atoms, ityp, name)
use module_types
implicit none
type(atoms_data), intent(inout) :: atoms
integer, intent(in) :: ityp
character(len=1), dimension(20), intent(in) :: name
*/
void FC_FUNC_(atoms_set_name, ATOMS_SET_NAME)(f90_atoms_data *atoms, 
                                              const int *ityp, 
                                              const char *name, 
                                              int str_ln_1);
/* atoms_write src/init/atoms.f90:9 */
/* Fortran header:
subroutine atoms_write(atoms, filename, forces, energy, comment)
use module_defs, only: gp
use module_types
use module_atoms, only: astruct_dump_to_file
implicit none
character(len = *), intent(in) :: comment
character(len = *), intent(in) :: filename
type(atoms_data), intent(in) :: atoms
real(gp), intent(in) :: energy
real(gp), dimension(:,:), pointer :: forces
*/
void FC_FUNC_(atoms_write, ATOMS_WRITE)(const f90_atoms_data *atoms, 
                                        const char *filename, 
                                        f90_pointer_double_2D *forces, 
                                        const double *energy, 
                                        const char *comment, 
                                        int str_ln_1, 
                                        int str_ln_2);
/* bigdft_exec src/bindings/bindingsf.f90:1674 */
/* Fortran header:
subroutine bigdft_exec(runObj,outs,infocode)
use bigdft_run, only: run_objects,state_properties,bigdft_state
implicit none
type(run_objects), intent(inout) :: runObj
type(state_properties), intent(inout) :: outs
integer, intent(inout) :: infocode
*/
void FC_FUNC_(bigdft_exec, BIGDFT_EXEC)(f90_run_objects *runObj, 
                                        f90_state_properties *outs, 
                                        int *infocode);
/* bigdft_finalize src/external.f90:123 */
/* Fortran header:
subroutine bigdft_finalize(ierr)
use BigDFT_API
implicit none
integer, intent(out) :: ierr
*/
void FC_FUNC_(bigdft_finalize, BIGDFT_FINALIZE)(int *ierr);
/* bigdft_init_mpi_env src/external.f90:73 */
/* Fortran header:
subroutine bigdft_init_mpi_env(mpi_info,mpi_groupsize, ierr)
use BigDFT_API
implicit none

integer, dimension(4), intent(out) :: mpi_info
integer, intent(in) :: mpi_groupsize
integer, intent(out) :: ierr

integer :: iproc,nproc,ngroup_size
*/
void FC_FUNC_(bigdft_init_mpi_env, BIGDFT_INIT_MPI_ENV)(int *mpi_info, 
                                                        const int *mpi_groupsize, 
                                                        int *ierr);
/* bigdft_init_mpi_force src/external.f90:113 */
/* Fortran header:
subroutine bigdft_init_mpi_force(igroup, ngroup)
use BigDFT_API
implicit none
integer, intent(in) :: igroup, ngroup
*/
void FC_FUNC_(bigdft_init_mpi_force, BIGDFT_INIT_MPI_FORCE)(const int *igroup, 
                                                            const int *ngroup);
/* bigdft_mpi_init src/external.f90:58 */
/* Fortran header:
subroutine bigdft_mpi_init(ierr)
use wrapper_mpi, only: wmpi_init_thread,MPI_SUCCESS
use module_types, only: bigdft_init_errors,bigdft_init_timing_categories
implicit none
integer, intent(out) :: ierr
*/
void FC_FUNC_(bigdft_mpi_init, BIGDFT_MPI_INIT)(int *ierr);
/* call_external_c_fromadd  */
/* check_linear_and_create_lzd src/linear/initAndUtils.f90:51 */
/* Fortran header:
subroutine check_linear_and_create_Lzd(iproc,nproc,linType,Lzd,atoms,orbs,nspin,rxyz)
use module_base
use module_types
use module_xc
use ao_inguess, only: atomic_info
use locregs, only: locreg_null,copy_locreg_descriptors
use public_enums
use locregs_init, only: determine_locreg_parallel, check_linear_inputguess
implicit none

integer, intent(in) :: iproc,nproc,nspin
type(local_zone_descriptors), intent(inout) :: Lzd
type(atoms_data), intent(in) :: atoms
type(orbitals_data), intent(inout) :: orbs
real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
integer, intent(in) :: linType


character(len=*), parameter :: subname='check_linear_and_create_Lzd'
logical :: linear
real(gp) :: rcov
integer :: iat,ityp,nspin_ig,ilr
real(gp), dimension(:), allocatable :: locrad
logical,dimension(:), allocatable :: calculateBounds
*/
void FC_FUNC_(check_linear_and_create_lzd, CHECK_LINEAR_AND_CREATE_LZD)(const int *iproc, 
                                                                        const int *nproc, 
                                                                        const int *linType, 
                                                                        f90_local_zone_descriptors *Lzd, 
                                                                        const f90_atoms_data *atoms, 
                                                                        f90_orbitals_data *orbs, 
                                                                        const int *nspin, 
                                                                        const double *rxyz);
/* close_file src/bindings/bindingsf.f90:149 */
/* Fortran header:
subroutine close_file(unitwf)
implicit none
integer, intent(in) :: unitwf
*/
void FC_FUNC_(close_file, CLOSE_FILE)(const int *unitwf);
/* createeffectiveionicpotential src/init/ionicpot.f90:1746 */
/* Fortran header:
subroutine createEffectiveIonicPotential(iproc, verb, input, atoms, rxyz, shift,  dpbox, pkernel, pot_ion, rho_ion, elecfield, psoffset)

use module_base
use module_dpbox, only: denspot_distribution
use module_types

implicit none


integer, intent(in) :: iproc
logical, intent(in) :: verb

real(gp), intent(in) :: psoffset
type(atoms_data), intent(in) :: atoms

type(input_variables), intent(in) :: input
type(denspot_distribution), intent(in) :: dpbox
real(gp), dimension(3), intent(in) :: elecfield
real(gp), dimension(3), intent(in) :: shift
real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
type(coulomb_operator), intent(inout) :: pkernel
real(wp), dimension(*), intent(inout) :: pot_ion
real(wp), dimension(*), intent(inout) :: rho_ion


logical :: counterions
real(dp), dimension(:), allocatable :: counter_ions
integer :: ncounter_ions
*/
void FC_FUNC(createeffectiveionicpotential, CREATEEFFECTIVEIONICPOTENTIAL)(const int *iproc, 
                                                                           const int *verb, 
                                                                           const f90_input_variables *input, 
                                                                           const f90_atoms_data *atoms, 
                                                                           const double *rxyz, 
                                                                           const double *shift, 
                                                                           const f90_denspot_distribution *dpbox, 
                                                                           f90_coulomb_operator *pkernel, 
                                                                           double *pot_ion, 
                                                                           double *rho_ion, 
                                                                           const double *elecfield, 
                                                                           const double *psoffset);
/* createprojectorsarrays src/init.f90:217 */
/* Fortran header:
subroutine createProjectorsArrays(lr,rxyz,at,ob,cpmult,fpmult,hx,hy,hz,dry_run,nl,init_projectors_completely)
use module_base
use psp_projectors_base, only: DFT_PSP_projectors_null, nonlocal_psp_descriptors_null, allocate_workarrays_projectors
use psp_projectors, only: set_nlpsp_to_wfd, bounds_to_plr_limits
use module_types
use gaussians, only: gaussian_basis, gaussian_basis_from_psp, gaussian_basis_from_paw
use public_enums, only: PSPCODE_PAW
use orbitalbasis
use ao_inguess, only: lmax_ao
implicit none
real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
type(locreg_descriptors),intent(in) :: lr
type(atoms_data), intent(in) :: at

type(orbital_basis), intent(in) :: ob
real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz

logical, intent(in) :: dry_run 
type(DFT_PSP_projectors), intent(out) :: nl
logical,intent(in) :: init_projectors_completely 

character(len=*), parameter :: subname='createProjectorsArrays'
integer :: n1,n2,n3,nl1,nl2,nl3,nu1,nu2,nu3,mseg,nbseg_dim,npack_dim,mproj_max
integer :: iat,iseg

integer, dimension(:), allocatable :: nbsegs_cf,keyg_lin
logical, dimension(:,:,:), allocatable :: logrid
*/
void FC_FUNC(createprojectorsarrays, CREATEPROJECTORSARRAYS)(const f90_locreg_descriptors *lr, 
                                                             const double *rxyz, 
                                                             const f90_atoms_data *at, 
                                                             const f90_orbital_basis *ob, 
                                                             const double *cpmult, 
                                                             const double *fpmult, 
                                                             const double *hx, 
                                                             const double *hy, 
                                                             const double *hz, 
                                                             const int *dry_run, 
                                                             f90_DFT_PSP_projectors *nl, 
                                                             const int *init_projectors_completely);
/* deallocate_double_1d src/bindings/bindingsf.f90:170 */
/* Fortran header:
subroutine deallocate_double_1D(array)
use dynamic_memory, only: f_free_ptr
implicit none

double precision, dimension(:), pointer :: array
*/
void FC_FUNC_(deallocate_double_1d, DEALLOCATE_DOUBLE_1D)(f90_pointer_double *array);
/* density_descriptors src/init/denspotd.f90:600 */
/* Fortran header:
subroutine density_descriptors(iproc,nproc,xc,nspin,crmult,frmult,atoms,dpbox,rho_commun,rxyz,rhodsc)
use module_base
use module_dpbox, only:  denspot_distribution
use module_types
use module_xc
implicit none
integer, intent(in) :: iproc,nproc,nspin
type(xc_info), intent(in) :: xc
real(gp), intent(in) :: crmult,frmult
type(atoms_data), intent(in) :: atoms
type(denspot_distribution), intent(in) :: dpbox
character(len=3), intent(in) :: rho_commun
real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz

type(rho_descriptors), intent(out) :: rhodsc
*/
void FC_FUNC_(density_descriptors, DENSITY_DESCRIPTORS)(const int *iproc, 
                                                        const int *nproc, 
                                                        const f90_xc_info *xc, 
                                                        const int *nspin, 
                                                        const double *crmult, 
                                                        const double *frmult, 
                                                        const f90_atoms_data *atoms, 
                                                        const f90_denspot_distribution *dpbox, 
                                                        const char *rho_commun, 
                                                        const double *rxyz, 
                                                        f90_rho_descriptors *rhodsc, 
                                                        int str_ln_1);
/* denspot_communications src/init/denspotd.f90:233 */
/* Fortran header:
subroutine denspot_communications(iproc,nproc,igpu,xc,nspin,geocode,SICapproach,dpbox)
use module_base
use module_dpbox, only: denspot_distribution
use module_types
use module_xc
implicit none
integer, intent(in) :: nspin,iproc,nproc,igpu
type(xc_info), intent(in) :: xc
character(len=1), intent(in) :: geocode 
character(len=4), intent(in) :: SICapproach
type(denspot_distribution), intent(inout) :: dpbox
*/
void FC_FUNC_(denspot_communications, DENSPOT_COMMUNICATIONS)(const int *iproc, 
                                                              const int *nproc, 
                                                              const int *igpu, 
                                                              const f90_xc_info *xc, 
                                                              const int *nspin, 
                                                              const char *geocode, 
                                                              const char *SICapproach, 
                                                              f90_denspot_distribution *dpbox, 
                                                              int str_ln_1, 
                                                              int str_ln_2);
/* denspot_full_density src/init/denspotd.f90:288 src/init/denspotd.f90:383 */
/* Fortran header:
subroutine denspot_full_density(denspot, rho_full, iproc, new)
use module_base
use module_types
use memory_profiling
implicit none
type(DFT_local_fields), intent(in) :: denspot
integer, intent(in) :: iproc
integer, intent(out) :: new
real(gp), dimension(:), pointer :: rho_full

character(len = *), parameter :: subname = "denspot_full_density"
integer :: nslice, ierr, irhodim, irhoxcsh
*/
void FC_FUNC_(denspot_full_density, DENSPOT_FULL_DENSITY)(const f90_DFT_local_fields *denspot, 
                                                          f90_pointer_double *rho_full, 
                                                          const int *iproc, 
                                                          int *new);
/* denspot_full_v_ext src/init/denspotd.f90:337 src/init/denspotd.f90:438 */
/* Fortran header:
subroutine denspot_full_v_ext(denspot, pot_full, iproc, new)
use module_base
use module_types
use memory_profiling
implicit none
type(DFT_local_fields), intent(in) :: denspot
integer, intent(in) :: iproc
integer, intent(out) :: new
real(gp), pointer :: pot_full(:)

character(len = *), parameter :: subname = "localfields_full_potential"
integer :: ierr
*/
void FC_FUNC_(denspot_full_v_ext, DENSPOT_FULL_V_EXT)(const f90_DFT_local_fields *denspot, 
                                                      f90_pointer_double *pot_full, 
                                                      const int *iproc, 
                                                      int *new);
/* dict_append src/bindings/bindingsf.f90:1808 */
/* Fortran header:
subroutine dict_append(dict)
use dictionaries, only: dictionary, operator(//), dict_len
implicit none
type(dictionary), pointer :: dict
*/
void FC_FUNC_(dict_append, DICT_APPEND)(f90_dictionary_pointer *dict);
/* dict_dump src/bindings/bindingsf.f90:1828 */
/* Fortran header:
subroutine dict_dump(dict, unit)
use dictionaries, only: dictionary
use yaml_output, only: yaml_dict_dump
implicit none
type(dictionary), pointer :: dict
integer, intent(in) :: unit
*/
void FC_FUNC_(dict_dump, DICT_DUMP)(f90_dictionary_pointer *dict, 
                                    const int *unit);
/* dict_dump_to_file src/bindings/bindingsf.f90:1843 */
/* Fortran header:
subroutine dict_dump_to_file(dict, path)
use dictionaries, only: dictionary
use yaml_output, only: yaml_dict_dump, yaml_set_stream, yaml_close_stream, yaml_get_default_stream
implicit none
type(dictionary), pointer :: dict
character(len = *), intent(in) :: path

integer :: unit
*/
void FC_FUNC_(dict_dump_to_file, DICT_DUMP_TO_FILE)(f90_dictionary_pointer *dict, 
                                                    const char *path, 
                                                    int str_ln_1);
/* dict_free_binding  */
/* dict_init_binding src/bindings/bindingsf.f90:1978 */
/* Fortran header:
subroutine dict_init_binding(dict)
use dictionaries, only: dictionary, wrapper => dict_init
implicit none
type(dictionary), pointer :: dict
*/
void FC_FUNC_(dict_init_binding, DICT_INIT_BINDING)(f90_dictionary_pointer *dict);
/* dict_insert src/bindings/bindingsf.f90:1797 */
/* Fortran header:
subroutine dict_insert(dict, key)
use dictionaries, only: dictionary, operator(//)
implicit none
type(dictionary), pointer :: dict
character(len = *), intent(in) :: key
*/
void FC_FUNC_(dict_insert, DICT_INSERT)(f90_dictionary_pointer *dict, 
                                        const char *key, 
                                        int str_ln_1);
/* dict_iter src/bindings/bindingsf.f90:1912 */
/* Fortran header:
subroutine dict_iter(dict, exists)
use dictionaries, only: dictionary, wrapper => dict_iter
implicit none
type(dictionary), pointer :: dict
logical, intent(out) :: exists

type(dictionary), pointer :: start
*/
void FC_FUNC_(dict_iter, DICT_ITER)(f90_dictionary_pointer *dict, 
                                    int *exists);
/* dict_key src/bindings/bindingsf.f90:1902 */
/* Fortran header:
subroutine dict_key(dict, buf)
use dictionaries, only: dictionary, max_field_length, wrapper => dict_key
implicit none
type(dictionary), pointer :: dict
character(len = max_field_length), intent(out) :: buf
*/
void FC_FUNC_(dict_key, DICT_KEY)(f90_dictionary_pointer *dict, 
                                  char *buf, 
                                  int str_ln_1);
/* dict_len src/bindings/bindingsf.f90:1940 */
/* Fortran header:
subroutine dict_len(dict, ln)
use dictionaries, only: dictionary, wrapper => dict_len
implicit none
type(dictionary), pointer :: dict
integer, intent(out) :: ln
*/
void FC_FUNC_(dict_len, DICT_LEN)(f90_dictionary_pointer *dict, 
                                  int *ln);
/* dict_move_to_item src/bindings/bindingsf.f90:1785 */
/* Fortran header:
subroutine dict_move_to_item(dict, exists, id)
use dictionaries, only: dictionary, operator(//), dict_len
implicit none
type(dictionary), pointer :: dict
logical, intent(out) :: exists
integer, intent(in) :: id
*/
void FC_FUNC_(dict_move_to_item, DICT_MOVE_TO_ITEM)(f90_dictionary_pointer *dict, 
                                                    int *exists, 
                                                    const int *id);
/* dict_move_to_key src/bindings/bindingsf.f90:1773 */
/* Fortran header:
subroutine dict_move_to_key(dict, exists, key)
use dictionaries, only: dictionary, operator(//), has_key
implicit none
type(dictionary), pointer :: dict
logical, intent(out) :: exists
character(len = *), intent(in) :: key
*/
void FC_FUNC_(dict_move_to_key, DICT_MOVE_TO_KEY)(f90_dictionary_pointer *dict, 
                                                  int *exists, 
                                                  const char *key, 
                                                  int str_ln_1);
/* dict_next src/bindings/bindingsf.f90:1926 */
/* Fortran header:
subroutine dict_next(dict, exists)
use dictionaries, only: dictionary, wrapper => dict_next
implicit none
type(dictionary), pointer :: dict
logical, intent(out) :: exists

type(dictionary), pointer :: next
*/
void FC_FUNC_(dict_next, DICT_NEXT)(f90_dictionary_pointer *dict, 
                                    int *exists);
/* dict_parse src/bindings/bindingsf.f90:1860 */
/* Fortran header:
subroutine dict_parse(dict, buf)
use dictionaries, only: dictionary, operator(//), dict_len,operator(.pop.),dict_free
use yaml_parse, only: yaml_parse_from_string
implicit none
type(dictionary), pointer :: dict
character(len = *), intent(in) :: buf
type(dictionary), pointer :: dict_load
*/
void FC_FUNC_(dict_parse, DICT_PARSE)(f90_dictionary_pointer *dict, 
                                      const char *buf, 
                                      int str_ln_1);
/* dict_pop src/bindings/bindingsf.f90:1877 */
/* Fortran header:
subroutine dict_pop(dict, exists, key)
use dictionaries, only: dictionary, has_key, dict_remove, dict_init
implicit none
type(dictionary), pointer :: dict
logical, intent(out) :: exists
character(len = *), intent(in) :: key
*/
void FC_FUNC_(dict_pop, DICT_POP)(f90_dictionary_pointer *dict, 
                                  int *exists, 
                                  const char *key, 
                                  int str_ln_1);
/* dict_put src/bindings/bindingsf.f90:1817 */
/* Fortran header:
subroutine dict_put(dict, val)
use dictionaries, only: dictionary, set
implicit none
type(dictionary), pointer :: dict
character(len = *), intent(in) :: val
*/
void FC_FUNC_(dict_put, DICT_PUT)(f90_dictionary_pointer *dict, 
                                  const char *val, 
                                  int str_ln_1);
/* dict_size src/bindings/bindingsf.f90:1950 */
/* Fortran header:
subroutine dict_size(dict, ln)
use dictionaries, only: dictionary, wrapper => dict_size
implicit none
type(dictionary), pointer :: dict
integer, intent(out) :: ln
*/
void FC_FUNC_(dict_size, DICT_SIZE)(f90_dictionary_pointer *dict, 
                                    int *ln);
/* dict_update_binding src/bindings/bindingsf.f90:1969 */
/* Fortran header:
subroutine dict_update_binding(dict, ref)
use dictionaries, only: dictionary, wrapper => dict_update
implicit none
type(dictionary), pointer :: dict, ref
*/
void FC_FUNC_(dict_update_binding, DICT_UPDATE_BINDING)(f90_dictionary_pointer *dict, 
                                                        f90_dictionary_pointer *ref);
/* dict_value_binding src/bindings/bindingsf.f90:1892 */
/* Fortran header:
subroutine dict_value_binding(dict, buf)
use dictionaries, only: dictionary, max_field_length, wrapper => dict_value
implicit none
type(dictionary), pointer :: dict
character(len = max_field_length), intent(out) :: buf
*/
void FC_FUNC_(dict_value_binding, DICT_VALUE_BINDING)(f90_dictionary_pointer *dict, 
                                                      char *buf, 
                                                      int str_ln_1);
/* dpbox_set_box src/init/denspotd.f90:126 */
/* Fortran header:
subroutine dpbox_set_box(dpbox,Lzd)
use module_base
use module_dpbox, only: denspot_distribution
use module_types
implicit none
type(local_zone_descriptors), intent(in) :: Lzd
type(denspot_distribution), intent(inout) :: dpbox
*/
void FC_FUNC_(dpbox_set_box, DPBOX_SET_BOX)(f90_denspot_distribution *dpbox, 
                                            const f90_local_zone_descriptors *Lzd);
/* energs_copy_data src/bindings/bindingsf.f90:1356 */
/* Fortran header:
subroutine energs_copy_data(energs, eh, exc, evxc, eion, edisp, ekin, epot,  eproj, eexctX, ebs, eKS, trH, evsum, evsic)
use module_defs, only: gp
use module_types
implicit none
type(energy_terms), intent(in) :: energs
real(gp), intent(out) :: eh, exc, evxc, eion, edisp, ekin, epot, eproj,  eexctX, ebs, eKS, trH, evsum, evsic
*/
void FC_FUNC_(energs_copy_data, ENERGS_COPY_DATA)(const f90_energy_terms *energs, 
                                                  double *eh, 
                                                  double *exc, 
                                                  double *evxc, 
                                                  double *eion, 
                                                  double *edisp, 
                                                  double *ekin, 
                                                  double *epot, 
                                                  double *eproj, 
                                                  double *eexctX, 
                                                  double *ebs, 
                                                  double *eKS, 
                                                  double *trH, 
                                                  double *evsum, 
                                                  double *evsic);
/* fill_logrid src/init/gridmanipulation.f90:510 */
/* Fortran header:
subroutine fill_logrid(geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,     ntypes,iatype,rxyz,radii,rmult,hx,hy,hz,logrid)
use module_base
implicit none

character(len=*), intent(in) :: geocode 
integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,ntypes
real(gp), intent(in) :: rmult,hx,hy,hz
integer, dimension(nat), intent(in) :: iatype
real(gp), dimension(ntypes), intent(in) :: radii
real(gp), dimension(3,nat), intent(in) :: rxyz
logical, dimension(0:n1,0:n2,0:n3), intent(out) :: logrid

real(kind=8), parameter :: eps_mach=1.d-12
integer :: i1,i2,i3,iat,ml1,ml2,ml3,mu1,mu2,mu3,j1,j2,j3,i1s,i1e,i2s,i2e,i3s,i3e
real(gp) :: dx,dy2,dz2,rad,dy2pdz2,radsq
*/
void FC_FUNC_(fill_logrid, FILL_LOGRID)(const char *geocode, 
                                        const int *n1, 
                                        const int *n2, 
                                        const int *n3, 
                                        const int *nl1, 
                                        const int *nu1, 
                                        const int *nl2, 
                                        const int *nu2, 
                                        const int *nl3, 
                                        const int *nu3, 
                                        const int *nbuf, 
                                        const int *nat, 
                                        const int *ntypes, 
                                        const int *iatype, 
                                        const double *rxyz, 
                                        const double *radii, 
                                        const double *rmult, 
                                        const double *hx, 
                                        const double *hy, 
                                        const double *hz, 
                                        int *logrid, 
                                        int str_ln_1);
/* f_lib_finalize  */
/* f_lib_initialize  */
/* free_wave_to_isf  */
/* geopt src/geopt/geometry.f90:45 */
/* Fortran header:
subroutine geopt(runObj,outs,nproc,iproc,ncount_bigdft)
use module_base
use bigdft_run
use yaml_output
use minpar
implicit none

type(run_objects), intent(inout) :: runObj
type(state_properties), intent(inout) :: outs
integer, intent(in) :: nproc,iproc
integer, intent(out) :: ncount_bigdft

integer, parameter :: ugeopt = 16
logical :: fail
integer :: ibfgs,ierr
character(len=6) :: outfile, fmt
character(len=5) :: fn4
character(len=40) :: comment
character(len=60) :: filename
*/
void FC_FUNC(geopt, GEOPT)(f90_run_objects *runObj, 
                           f90_state_properties *outs, 
                           const int *nproc, 
                           const int *iproc, 
                           int *ncount_bigdft);
/* glr_copy src/bindings/bindingsf.f90:199 */
/* Fortran header:
subroutine glr_copy(glr, d, wfd, from)
use module_types
use locregs, only: copy_locreg_descriptors
implicit none
type(locreg_descriptors), pointer :: glr
type(grid_dimensions), pointer :: d
type(wavefunctions_descriptors), pointer :: wfd
type(locreg_descriptors), intent(in) :: from
*/
void FC_FUNC_(glr_copy, GLR_COPY)(f90_locreg_descriptors_pointer *glr, 
                                  f90_grid_dimensions_pointer *d, 
                                  f90_wavefunctions_descriptors_pointer *wfd, 
                                  const f90_locreg_descriptors *from);
/* glr_empty src/bindings/bindingsf.f90:250 */
/* Fortran header:
subroutine glr_empty(glr)
use locregs
implicit none
type(locreg_descriptors), intent(inout) :: glr
*/
void FC_FUNC_(glr_empty, GLR_EMPTY)(f90_locreg_descriptors *glr);
/* glr_free src/bindings/bindingsf.f90:241 */
/* Fortran header:
subroutine glr_free(glr)
use module_types
implicit none
type(locreg_descriptors), pointer :: glr
*/
void FC_FUNC_(glr_free, GLR_FREE)(f90_locreg_descriptors_pointer *glr);
/* glr_get_data src/bindings/bindingsf.f90:229 */
/* Fortran header:
subroutine glr_get_data(glr, d, wfd)
use module_types
implicit none
type(locreg_descriptors), intent(inout), target :: glr
type(grid_dimensions), pointer :: d
type(wavefunctions_descriptors), pointer :: wfd
*/
void FC_FUNC_(glr_get_data, GLR_GET_DATA)(f90_locreg_descriptors *glr, 
                                          f90_grid_dimensions_pointer *d, 
                                          f90_wavefunctions_descriptors_pointer *wfd);
/* glr_get_dimensions src/bindings/bindingsf.f90:259 */
/* Fortran header:
subroutine glr_get_dimensions(glr , n, ni, ns, nsi, nfl, nfu, norb)
use module_types
implicit none
type(locreg_descriptors), intent(in) :: glr
integer, dimension(3), intent(out) :: n, ni, ns, nsi, nfl, nfu
integer, intent(out) :: norb
*/
void FC_FUNC_(glr_get_dimensions, GLR_GET_DIMENSIONS)(const f90_locreg_descriptors *glr, 
                                                      int *n, 
                                                      int *ni, 
                                                      int *ns, 
                                                      int *nsi, 
                                                      int *nfl, 
                                                      int *nfu, 
                                                      int *norb);
/* glr_get_locreg_data src/bindings/bindingsf.f90:320 */
/* Fortran header:
subroutine glr_get_locreg_data(glr, locrad, locregCenter)
use module_types
implicit none
type(locreg_descriptors), intent(in) :: glr
double precision, dimension(3), intent(out) :: locregCenter
double precision, intent(out) :: locrad
*/
void FC_FUNC_(glr_get_locreg_data, GLR_GET_LOCREG_DATA)(const f90_locreg_descriptors *glr, 
                                                        double *locrad, 
                                                        double *locregCenter);
/* glr_get_psi_size src/init/kswfn.f90:11 */
/* Fortran header:
subroutine glr_get_psi_size(glr, psisize)
use module_types
implicit none
type(locreg_descriptors), intent(in) :: glr
integer, intent(out) :: psisize
*/
void FC_FUNC_(glr_get_psi_size, GLR_GET_PSI_SIZE)(const f90_locreg_descriptors *glr, 
                                                  int *psisize);
/* glr_init src/bindings/bindingsf.f90:216 */
/* Fortran header:
subroutine glr_init(glr, d, wfd)
use module_types
implicit none
type(locreg_descriptors), intent(inout), target :: glr
type(grid_dimensions), pointer :: d
type(wavefunctions_descriptors), pointer :: wfd
*/
void FC_FUNC_(glr_init, GLR_INIT)(f90_locreg_descriptors *glr, 
                                  f90_grid_dimensions_pointer *d, 
                                  f90_wavefunctions_descriptors_pointer *wfd);
/* glr_new src/bindings/bindingsf.f90:190 */
/* Fortran header:
subroutine glr_new(glr)
use module_types
implicit none
type(locreg_descriptors), pointer :: glr
*/
void FC_FUNC_(glr_new, GLR_NEW)(f90_locreg_descriptors_pointer *glr);
/* glr_set_bounds src/bindings/bindingsf.f90:366 */
/* Fortran header:
subroutine glr_set_bounds(lr)
use module_types
use bounds, only: locreg_bounds
implicit none
type(locreg_descriptors), intent(inout) :: lr
*/
void FC_FUNC_(glr_set_bounds, GLR_SET_BOUNDS)(f90_locreg_descriptors *lr);
/* glr_set_dimensions src/bindings/bindingsf.f90:291 */
/* Fortran header:
subroutine glr_set_dimensions(glr, n, ni, ns, nsi, nfl, nfu)
use module_types
implicit none
type(locreg_descriptors), intent(inout) :: glr
integer, dimension(3), intent(in) :: n, ni, ns, nsi, nfl, nfu
*/
void FC_FUNC_(glr_set_dimensions, GLR_SET_DIMENSIONS)(f90_locreg_descriptors *glr, 
                                                      const int *n, 
                                                      const int *ni, 
                                                      const int *ns, 
                                                      const int *nsi, 
                                                      const int *nfl, 
                                                      const int *nfu);
/* glr_set_wave_descriptors src/bindings/bindingsf.f90:347 */
/* Fortran header:
subroutine glr_set_wave_descriptors(iproc,hx,hy,hz,atoms,rxyz,   crmult,frmult,Glr)
use module_base, only: gp
use module_types
use module_interfaces, only: createWavefunctionsDescriptors
implicit none

type(atoms_data), intent(in) :: atoms
integer, intent(in) :: iproc
real(gp), intent(in) :: hx,hy,hz,crmult,frmult
real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz

type(locreg_descriptors), intent(inout) :: Glr
*/
void FC_FUNC_(glr_set_wave_descriptors, GLR_SET_WAVE_DESCRIPTORS)(const int *iproc, 
                                                                  const double *hx, 
                                                                  const double *hy, 
                                                                  const double *hz, 
                                                                  const f90_atoms_data *atoms, 
                                                                  const double *rxyz, 
                                                                  const double *crmult, 
                                                                  const double *frmult, 
                                                                  f90_locreg_descriptors *Glr);
/* glr_set_wfd_dims src/bindings/bindingsf.f90:332 */
/* Fortran header:
subroutine glr_set_wfd_dims(glr, nseg_c, nseg_f, nvctr_c, nvctr_f)
use module_types
use locregs, only: allocate_wfd
implicit none
type(locreg_descriptors), intent(inout) :: glr
integer, intent(in) :: nseg_c, nseg_f, nvctr_c, nvctr_f
*/
void FC_FUNC_(glr_set_wfd_dims, GLR_SET_WFD_DIMS)(f90_locreg_descriptors *glr, 
                                                  const int *nseg_c, 
                                                  const int *nseg_f, 
                                                  const int *nvctr_c, 
                                                  const int *nvctr_f);
/* glr_wfd_get_data src/bindings/bindingsf.f90:378 */
/* Fortran header:
subroutine glr_wfd_get_data(wfd, nvctr_c, nvctr_f, nseg_c, nseg_f,  keyglob, keygloc, keyvglob, keyvloc)
use module_types
implicit none
type(wavefunctions_descriptors), intent(in) :: wfd
integer, intent(out) :: nvctr_c, nvctr_f, nseg_c, nseg_f
integer, dimension(:,:), pointer :: keyglob, keygloc
integer, dimension(:), pointer :: keyvglob, keyvloc
*/
void FC_FUNC_(glr_wfd_get_data, GLR_WFD_GET_DATA)(const f90_wavefunctions_descriptors *wfd, 
                                                  int *nvctr_c, 
                                                  int *nvctr_f, 
                                                  int *nseg_c, 
                                                  int *nseg_f, 
                                                  f90_pointer_int_2D *keyglob, 
                                                  f90_pointer_int_2D *keygloc, 
                                                  f90_pointer_int *keyvglob, 
                                                  f90_pointer_int *keyvloc);
/* gpu_free src/bindings/bindingsf.f90:1115 */
/* Fortran header:
subroutine gpu_free(GPU)
use module_types
implicit none
type(GPU_pointers), pointer :: GPU
*/
void FC_FUNC_(gpu_free, GPU_FREE)(f90_GPU_pointers_pointer *GPU);
/* gpu_new src/bindings/bindingsf.f90:1106 */
/* Fortran header:
subroutine gpu_new(GPU)
use module_types
implicit none
type(GPU_pointers), pointer :: GPU
*/
void FC_FUNC_(gpu_new, GPU_NEW)(f90_GPU_pointers_pointer *GPU);
/* initialize_dft_local_fields src/init/denspotd.f90:12 */
/* Fortran header:
subroutine initialize_DFT_local_fields(denspot, ixc, nspden)
use module_base
use module_dpbox, only: dpbox_null
use module_types
use module_xc
use public_enums
implicit none
type(DFT_local_fields), intent(inout) :: denspot
integer, intent(in) :: ixc, nspden
*/
void FC_FUNC_(initialize_dft_local_fields, INITIALIZE_DFT_LOCAL_FIELDS)(f90_DFT_local_fields *denspot, 
                                                                        const int *ixc, 
                                                                        const int *nspden);
/* init_orbitals_data_for_linear src/linear/initAndUtils.f90:350 */
/* Fortran header:
subroutine init_orbitals_data_for_linear(iproc, nproc, nspinor, input, astruct, rxyz, lorbs, norb_par_ref, norbu_par_ref, norbd_par_ref)
use module_base
use module_types
use module_interfaces, only: assignToLocreg2, orbitals_descriptors
use public_enums
implicit none


integer, intent(in) :: iproc, nproc, nspinor
type(input_variables), intent(in) :: input
type(atomic_structure), intent(in) :: astruct
real(kind=8),dimension(3,astruct%nat), intent(in) :: rxyz
type(orbitals_data), intent(out) :: lorbs
integer,dimension(0:nproc-1),intent(in),optional :: norb_par_ref, norbu_par_ref, norbd_par_ref


integer :: norb, norbu, norbd, ityp, iat, ilr, iorb, nlr, iiat, ispin
integer, dimension(:), allocatable :: norbsPerLocreg, norbsPerAtom
real(kind=8),dimension(:,:), allocatable :: locregCenter
character(len=*), parameter :: subname='init_orbitals_data_for_linear'
logical :: with_optional
logical,dimension(3) :: optional_present
*/
void FC_FUNC_(init_orbitals_data_for_linear, INIT_ORBITALS_DATA_FOR_LINEAR)(const int *iproc, 
                                                                            const int *nproc, 
                                                                            const int *nspinor, 
                                                                            const f90_input_variables *input, 
                                                                            const f90_atomic_structure *astruct, 
                                                                            const double *rxyz, 
                                                                            f90_orbitals_data *lorbs, 
                                                                            const int *norb_par_ref, 
                                                                            const int *norbu_par_ref, 
                                                                            const int *norbd_par_ref);
/* inputs_check_psi_id src/bindings/bindingsf.f90:685 */
/* Fortran header:
subroutine inputs_check_psi_id(inputpsi, input_wf_format, dir_output, ln, orbs, lorbs, iproc, nproc)
use module_types
use module_fragments
use module_interfaces, only: input_check_psi_id
use f_enums
implicit none
integer, intent(out) :: input_wf_format
type(f_enumerator), intent(inout) :: inputpsi
integer, intent(in) :: iproc, ln, nproc
character(len = ln), intent(in) :: dir_output
type(orbitals_data), intent(in) :: orbs, lorbs

type(system_fragment), dimension(:), pointer :: ref_frags  
character(len=100) :: frag_dir
*/
void FC_FUNC_(inputs_check_psi_id, INPUTS_CHECK_PSI_ID)(f90_f_enumerator *inputpsi, 
                                                        int *input_wf_format, 
                                                        const char *dir_output, 
                                                        const int *ln, 
                                                        const f90_orbitals_data *orbs, 
                                                        const f90_orbitals_data *lorbs, 
                                                        const int *iproc, 
                                                        const int *nproc, 
                                                        int str_ln_1);
/* inputs_free src/bindings/bindingsf.f90:520 */
/* Fortran header:
subroutine inputs_free(in)
use module_input_keys
implicit none
type(input_variables), pointer :: in
*/
void FC_FUNC_(inputs_free, INPUTS_FREE)(f90_input_variables_pointer *in);
/* inputs_get_dft src/bindings/bindingsf.f90:555 */
/* Fortran header:
subroutine inputs_get_dft(in, hx, hy, hz, crmult, frmult, ixc, qcharge, efield, nspin, mpol,  gnrm, itermax, nrepmax, ncong, idsx, dispcorr, inpsi, outpsi, outgrid,  rbuf, ncongt, davidson, nvirt, nplottedvirt, sym, last_run)
use module_defs, only: gp
use module_types
implicit none
type(input_variables), intent(in) :: in
real(gp), intent(out) :: hx, hy, hz, crmult, frmult, efield(3), gnrm, rbuf, qcharge
integer, intent(out) :: ixc, nspin, mpol, itermax, nrepmax, ncong, idsx,  dispcorr, inpsi, outpsi, outgrid, ncongt, davidson, nvirt, nplottedvirt,  sym, last_run
*/
void FC_FUNC_(inputs_get_dft, INPUTS_GET_DFT)(const f90_input_variables *in, 
                                              double *hx, 
                                              double *hy, 
                                              double *hz, 
                                              double *crmult, 
                                              double *frmult, 
                                              int *ixc, 
                                              double *qcharge, 
                                              double *efield, 
                                              int *nspin, 
                                              int *mpol, 
                                              double *gnrm, 
                                              int *itermax, 
                                              int *nrepmax, 
                                              int *ncong, 
                                              int *idsx, 
                                              int *dispcorr, 
                                              int *inpsi, 
                                              int *outpsi, 
                                              int *outgrid, 
                                              double *rbuf, 
                                              int *ncongt, 
                                              int *davidson, 
                                              int *nvirt, 
                                              int *nplottedvirt, 
                                              int *sym, 
                                              int *last_run);
/* inputs_get_geopt src/bindings/bindingsf.f90:624 */
/* Fortran header:
subroutine inputs_get_geopt(in, geopt_approach, ncount_cluster_x, frac_fluct, forcemax,  randdis, betax, history, ionmov, dtion, strtarget, qmass)
use module_defs, only: gp
use module_types
implicit none
type(input_variables), intent(in) :: in
character(len = 10), intent(out) :: geopt_approach
integer, intent(out) :: ncount_cluster_x, history, ionmov
real(gp), intent(out) :: frac_fluct, forcemax, randdis, betax, dtion, strtarget(6)
real(gp), pointer :: qmass(:)
*/
void FC_FUNC_(inputs_get_geopt, INPUTS_GET_GEOPT)(const f90_input_variables *in, 
                                                  char *geopt_approach, 
                                                  int *ncount_cluster_x, 
                                                  double *frac_fluct, 
                                                  double *forcemax, 
                                                  double *randdis, 
                                                  double *betax, 
                                                  int *history, 
                                                  int *ionmov, 
                                                  double *dtion, 
                                                  double *strtarget, 
                                                  f90_pointer_double *qmass, 
                                                  int str_ln_1);
/* inputs_get_linear src/bindings/bindingsf.f90:674 */
/* Fortran header:
subroutine inputs_get_linear(linear, inputPsiId)
use public_enums
implicit none
integer, intent(out) :: linear
integer, intent(in) :: inputPsiId
*/
void FC_FUNC_(inputs_get_linear, INPUTS_GET_LINEAR)(int *linear, 
                                                    const int *inputPsiId);
/* inputs_get_mix src/bindings/bindingsf.f90:600 */
/* Fortran header:
subroutine inputs_get_mix(in, iscf, itrpmax, norbsempty, occopt, alphamix, rpnrm_cv,  gnrm_startmix, Tel, alphadiis)
use module_defs, only: gp
use module_types
use module_input_keys, only: set_iscf

implicit none
type(input_variables), intent(in) :: in
integer, intent(out) :: iscf, itrpmax, norbsempty, occopt
real(gp), intent(out) :: alphamix, rpnrm_cv, gnrm_startmix, Tel, alphadiis
*/
void FC_FUNC_(inputs_get_mix, INPUTS_GET_MIX)(const f90_input_variables *in, 
                                              int *iscf, 
                                              int *itrpmax, 
                                              int *norbsempty, 
                                              int *occopt, 
                                              double *alphamix, 
                                              double *rpnrm_cv, 
                                              double *gnrm_startmix, 
                                              double *Tel, 
                                              double *alphadiis);
/* inputs_get_output src/bindings/bindingsf.f90:544 */
/* Fortran header:
subroutine inputs_get_output(in, dir_output)

use module_types
implicit none
type(input_variables), intent(in) :: in
character(len = 100), intent(out) :: dir_output
*/
void FC_FUNC_(inputs_get_output, INPUTS_GET_OUTPUT)(const f90_input_variables *in, 
                                                    char *dir_output, 
                                                    int str_ln_1);
/* inputs_get_perf src/bindings/bindingsf.f90:653 */
/* Fortran header:
subroutine inputs_get_perf(in, linear)
use module_types
implicit none
type(input_variables), intent(in) :: in
integer, intent(out) :: linear
*/
void FC_FUNC_(inputs_get_perf, INPUTS_GET_PERF)(const f90_input_variables *in, 
                                                int *linear);
/* inputs_set_dict src/bindings/bindingsf.f90:531 */
/* Fortran header:
subroutine inputs_set_dict(in, level, val)
use dictionaries
use module_input_keys, only: input_variables, input_set
implicit none
type(input_variables), intent(inout) :: in
character(len = *), intent(in) :: level
type(dictionary), pointer :: val
*/
void FC_FUNC_(inputs_set_dict, INPUTS_SET_DICT)(f90_input_variables *in, 
                                                const char *level, 
                                                f90_dictionary_pointer *val, 
                                                int str_ln_1);
/* input_wf src/cluster.f90:135 src/init.f90:2419 */
/* Fortran header:
subroutine input_wf(iproc,nproc,in,GPU,atoms,rxyz,denspot,denspot0,nlpsp,KSwfn,tmb,energs,inputpsi,input_wf_format,norbv,lzd_old,psi_old,rxyz_old,tmb_old,ref_frags,cdft,locregcenters)
use module_defs, only: gp,wp
use f_enums, only: f_enumerator
use module_types
use module_fragments
use constrained_dft
implicit none
integer, intent(in) :: iproc, nproc, input_wf_format
type(f_enumerator), intent(in) :: inputpsi
type(input_variables), intent(in) :: in
type(GPU_pointers), intent(inout) :: GPU
type(atoms_data), intent(inout) :: atoms
real(gp), dimension(3, atoms%astruct%nat), target, intent(in) :: rxyz
type(DFT_local_fields), intent(inout) :: denspot
type(DFT_wavefunction), intent(inout) :: KSwfn,tmb,tmb_old 
type(energy_terms), intent(inout) :: energs 
real(gp), dimension(*), intent(out) :: denspot0 
real(wp), dimension(:), pointer :: psi_old
integer, intent(out) :: norbv
type(DFT_PSP_projectors), intent(inout) :: nlpsp
real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz_old
type(local_zone_descriptors),intent(in):: lzd_old
type(system_fragment), dimension(:), pointer :: ref_frags
type(cdft_data), intent(out) :: cdft
real(kind=8),dimension(3,atoms%astruct%nat),intent(in),optional :: locregcenters
*/
void FC_FUNC_(input_wf, INPUT_WF)(const int *iproc, 
                                  const int *nproc, 
                                  const f90_input_variables *in, 
                                  f90_GPU_pointers *GPU, 
                                  f90_atoms_data *atoms, 
                                  const double *rxyz, 
                                  f90_DFT_local_fields *denspot, 
                                  double *denspot0, 
                                  f90_DFT_PSP_projectors *nlpsp, 
                                  f90_DFT_wavefunction *KSwfn, 
                                  f90_DFT_wavefunction *tmb, 
                                  f90_energy_terms *energs, 
                                  const f90_f_enumerator *inputpsi, 
                                  const int *input_wf_format, 
                                  int *norbv, 
                                  const f90_local_zone_descriptors *lzd_old, 
                                  f90_pointer_double *psi_old, 
                                  const double *rxyz_old, 
                                  f90_DFT_wavefunction *tmb_old, 
                                  f90_system_fragment_pointer *ref_frags, 
                                  f90_cdft_data *cdft, 
                                  const double *locregcenters);
/* ionicenergyandforces src/init/ionicpot.f90:12 */
/* Fortran header:
subroutine IonicEnergyandForces(iproc,nproc,dpbox,at,elecfield, rxyz,eion,fion,dispersion,edisp,fdisp,ewaldstr, pot_ion,pkernel,psoffset)
use module_base
use module_types
use Poisson_Solver, except_dp => dp, except_gp => gp
use gaussians, only: initialize_real_space_conversion, finalize_real_space_conversion,mp_exp
use module_atoms
use module_dpbox
use abi_interfaces_geometry, only: abi_metric
use abi_interfaces_common, only: abi_ewald, abi_ewald2
use m_paw_numeric, only: paw_splint
use abi_interfaces_numeric, only: abi_derf_ab
use vdwcorrection
use yaml_output
use public_enums, only: PSPCODE_PAW
use bounds, only: ext_buffers
implicit none

type(denspot_distribution), intent(in) :: dpbox
type(atoms_data), intent(in) :: at
integer, intent(in) :: iproc,nproc,dispersion
real(gp), dimension(3), intent(in) :: elecfield
real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
type(coulomb_operator), intent(inout) :: pkernel
real(gp), intent(out) :: eion,edisp,psoffset
real(dp), dimension(6),intent(out) :: ewaldstr
real(gp), dimension(:,:), pointer :: fion,fdisp
real(dp), dimension(*), intent(out) :: pot_ion

real(gp), parameter :: mp_tiny = 1.e-30_gp
logical :: slowion=.false.,use_iterator=.false.
logical :: perx,pery,perz,gox,goy,goz
integer ::  nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,n3i,n3pi,i3s
integer :: n1i,n2i,i,iat,ii,ityp,jat,jtyp
integer :: isx,iex,isy,iey,isz,iez,i1,i2,i3,j1,j2,j3,ind,ierr
real(gp) :: ucvol,rloc,rlocinv2sq,twopitothreehalf,atint,shortlength,charge,eself,rx,ry,rz
real(gp) :: fxion,fyion,fzion,dist,fxerf,fyerf,fzerf,cutoff
real(gp) :: hxh,hyh,hzh
real(gp) :: hxx,hxy,hxz,hyy,hyz,hzz,chgprod
real(gp) :: xp,Vel,prefactor,ehart,de,dv
real(gp) :: x,y,z,yp,zp,r2,arg

real(gp), dimension(3,3) :: gmet,rmet,rprimd,gprimd

real(gp), dimension(:,:), allocatable :: fewald,xred
real(gp), dimension(:), allocatable  :: mpx,mpy,mpz
real(gp), dimension(3) :: cc
real(gp), dimension(1) :: rr, vr
type(atoms_iterator) :: atit
type(dpbox_iterator) :: boxit
integer, dimension(2,3) :: nbox
*/
void FC_FUNC(ionicenergyandforces, IONICENERGYANDFORCES)(const int *iproc, 
                                                         const int *nproc, 
                                                         const f90_denspot_distribution *dpbox, 
                                                         const f90_atoms_data *at, 
                                                         const double *elecfield, 
                                                         const double *rxyz, 
                                                         double *eion, 
                                                         f90_pointer_double_2D *fion, 
                                                         const int *dispersion, 
                                                         double *edisp, 
                                                         f90_pointer_double_2D *fdisp, 
                                                         double *ewaldstr, 
                                                         double *pot_ion, 
                                                         f90_coulomb_operator *pkernel, 
                                                         double *psoffset);
/* kernel_get_comm src/bindings/bindingsf.f90:936 */
/* Fortran header:
subroutine kernel_get_comm(pkernel, igroup, ngroup, iproc_grp,  nproc_grp, mpi_comm)
use module_types
implicit none
type(coulomb_operator), intent(inout) :: pkernel
integer, intent(out) :: igroup, ngroup, iproc_grp, nproc_grp, mpi_comm
*/
void FC_FUNC_(kernel_get_comm, KERNEL_GET_COMM)(f90_coulomb_operator *pkernel, 
                                                int *igroup, 
                                                int *ngroup, 
                                                int *iproc_grp, 
                                                int *nproc_grp, 
                                                int *mpi_comm);
/* kswfn_init_comm src/init/kswfn.f90:118 */
/* Fortran header:
subroutine kswfn_init_comm(wfn, dpbox, iproc, nproc, nspin, imethod_overlap)
use module_dpbox, only: denspot_distribution
use module_types

use communications_base, only: comms_linear_null
use communications_init, only: init_comms_linear, init_comms_linear_sumrho, initialize_communication_potential
implicit none
integer, intent(in) :: iproc, nproc, nspin, imethod_overlap
type(DFT_wavefunction), intent(inout) :: wfn
type(denspot_distribution), intent(in) :: dpbox
*/
void FC_FUNC_(kswfn_init_comm, KSWFN_INIT_COMM)(f90_DFT_wavefunction *wfn, 
                                                const f90_denspot_distribution *dpbox, 
                                                const int *iproc, 
                                                const int *nproc, 
                                                const int *nspin, 
                                                const int *imethod_overlap);
/* kswfn_mpi_copy src/init/kswfn.f90:98 */
/* Fortran header:
subroutine kswfn_mpi_copy(psic, jproc, psiStart, psiSize)
use module_base
use module_types
implicit none
integer, intent(in) :: psiSize, jproc, psiStart
real(wp), intent(inout) :: psic(psiSize)

integer :: ierr
integer :: status(MPI_STATUS_SIZE)
*/
void FC_FUNC_(kswfn_mpi_copy, KSWFN_MPI_COPY)(double *psic, 
                                              const int *jproc, 
                                              const int *psiStart, 
                                              const int *psiSize);
/* kswfn_optimization_loop src/cluster.f90:1598 */
/* Fortran header:
subroutine kswfn_optimization_loop(iproc, nproc, opt,  alphamix, idsx, inputpsi, KSwfn, denspot, nlpsp, energs, atoms, GPU, xcstr,  in)
use module_base
use module_types
use module_interfaces, only: denspot_set_history, hpsitopsi, last_orthon, write_energies
use module_xc, only: XC_NO_HARTREE
use yaml_output
use public_enums
implicit none
real(dp), dimension(6), intent(out) :: xcstr
integer, intent(in) :: iproc, nproc, idsx
type(f_enumerator), intent(in) :: inputpsi
real(gp), intent(in) :: alphamix
type(DFT_optimization_loop), intent(inout) :: opt
type(DFT_wavefunction), intent(inout) :: KSwfn
type(DFT_local_fields), intent(inout) :: denspot
type(energy_terms), intent(inout) :: energs
type(atoms_data), intent(in) :: atoms
type(GPU_pointers), intent(inout) :: GPU
type(DFT_PSP_projectors), intent(inout) :: nlpsp
type(input_variables), intent(in) :: in 

character(len = *), parameter :: subname = "kswfn_optimization_loop"
logical :: endloop, scpot, endlooprp, lcs
integer :: ndiis_sd_sw, idsx_actual_before, linflag, ierr,iter_for_diis
integer :: ikpt_homo,ikpt_lumo,ispin_homo,ispin_lumo
real(gp) :: gnrm_zero,homo,lumo,occup_lumo,minres_gpe
character(len=5) :: final_out
*/
void FC_FUNC_(kswfn_optimization_loop, KSWFN_OPTIMIZATION_LOOP)(const int *iproc, 
                                                                const int *nproc, 
                                                                f90_DFT_optimization_loop *opt, 
                                                                const double *alphamix, 
                                                                const int *idsx, 
                                                                const f90_f_enumerator *inputpsi, 
                                                                f90_DFT_wavefunction *KSwfn, 
                                                                f90_DFT_local_fields *denspot, 
                                                                f90_DFT_PSP_projectors *nlpsp, 
                                                                f90_energy_terms *energs, 
                                                                const f90_atoms_data *atoms, 
                                                                f90_GPU_pointers *GPU, 
                                                                double *xcstr, 
                                                                const f90_input_variables *in);
/* kswfn_post_treatments src/cluster.f90:2016 */
/* Fortran header:
subroutine kswfn_post_treatments(iproc, nproc, KSwfn, tmb, linear,  fxyz, fnoise, fion, fdisp, fpulay,  strten, pressure, ewaldstr, xcstr,  GPU, denspot, atoms, rxyz, nlpsp,  output_denspot, dir_output, gridformat, refill_proj,  calculate_dipole, calculate_quadrupole, calculate_strten,nspin,  plot_pot_axes)
use module_base
use module_types
use module_interfaces, only: XC_potential, density_and_hpot
use Poisson_Solver, except_dp => dp, except_gp => gp
use yaml_output
use communications_base, only: deallocate_comms_linear, deallocate_p2pComms
use communications, only: synchronize_onesided_communication
use sparsematrix_base, only: deallocate_matrices, deallocate_sparse_matrix
use multipole, only: calculate_dipole_moment
use public_enums
use orbitalbasis
use io, only: plot_density
implicit none

type(DFT_wavefunction), intent(in) :: KSwfn
type(DFT_wavefunction), intent(inout) :: tmb
type(GPU_pointers), intent(inout) :: GPU
type(DFT_local_fields), intent(inout) :: denspot
type(atoms_data), intent(in) :: atoms
type(DFT_PSP_projectors), intent(inout) :: nlpsp
logical, intent(in) :: linear, refill_proj, calculate_dipole, calculate_quadrupole, calculate_strten
integer, intent(in) :: iproc, nproc, nspin
type(f_enumerator), intent(in) :: output_denspot
character(len = *), intent(in) :: dir_output
character(len = *), intent(in) :: gridformat
real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
real(gp), dimension(3, atoms%astruct%nat), intent(in) :: fdisp, fion, fpulay
real(dp), dimension(6), intent(in) :: ewaldstr
real(dp), dimension(6), intent(inout) :: xcstr
real(gp), intent(out) :: fnoise, pressure
real(gp), dimension(6), intent(out) :: strten
real(gp), dimension(3, atoms%astruct%nat), intent(out) :: fxyz
integer,dimension(3),intent(in) :: plot_pot_axes


character(len = *), parameter :: subname = "kswfn_post_treatments"
integer ::  jproc, nsize_psi, imode, i3xcsh_old

real(dp), dimension(6) :: hstrten
real(gp) :: ehart_fake, exc_fake, evxc_fake
type(orbital_basis) :: ob
*/
void FC_FUNC_(kswfn_post_treatments, KSWFN_POST_TREATMENTS)(const int *iproc, 
                                                            const int *nproc, 
                                                            const f90_DFT_wavefunction *KSwfn, 
                                                            f90_DFT_wavefunction *tmb, 
                                                            const int *linear, 
                                                            double *fxyz, 
                                                            double *fnoise, 
                                                            const double *fion, 
                                                            const double *fdisp, 
                                                            const double *fpulay, 
                                                            double *strten, 
                                                            double *pressure, 
                                                            const double *ewaldstr, 
                                                            double *xcstr, 
                                                            f90_GPU_pointers *GPU, 
                                                            f90_DFT_local_fields *denspot, 
                                                            const f90_atoms_data *atoms, 
                                                            const double *rxyz, 
                                                            f90_DFT_PSP_projectors *nlpsp, 
                                                            const f90_f_enumerator *output_denspot, 
                                                            const char *dir_output, 
                                                            const char *gridformat, 
                                                            const int *refill_proj, 
                                                            const int *calculate_dipole, 
                                                            const int *calculate_quadrupole, 
                                                            const int *calculate_strten, 
                                                            const int *nspin, 
                                                            const int *plot_pot_axes, 
                                                            int str_ln_1, 
                                                            int str_ln_2);
/* localfields_copy_metadata src/bindings/bindingsf.f90:1015 */
/* Fortran header:
subroutine localfields_copy_metadata(denspot, rhov_is, hgrid, ni, psoffset)
use module_defs, only: gp,dp
use module_types
implicit none
type(DFT_local_fields), intent(in) :: denspot
integer, intent(out) :: rhov_is, ni(3)
real(gp), intent(out) :: hgrid(3)
real(dp), intent(out) :: psoffset
*/
void FC_FUNC_(localfields_copy_metadata, LOCALFIELDS_COPY_METADATA)(const f90_DFT_local_fields *denspot, 
                                                                    int *rhov_is, 
                                                                    double *hgrid, 
                                                                    int *ni, 
                                                                    double *psoffset);
/* localfields_free src/bindings/bindingsf.f90:979 */
/* Fortran header:
subroutine localfields_free(denspotd, fion, fdisp)
use module_base
use module_dpbox, only: dpbox_free
use module_types
use Poisson_Solver, except_dp => dp, except_gp => gp
use memory_profiling
implicit none
type(DFT_local_fields), pointer :: denspotd
real(gp), dimension(:,:), pointer :: fion, fdisp

character(len = *), parameter :: subname = "localfields_free"
*/
void FC_FUNC_(localfields_free, LOCALFIELDS_FREE)(f90_DFT_local_fields_pointer *denspotd, 
                                                  f90_pointer_double_2D *fion, 
                                                  f90_pointer_double_2D *fdisp);
/* localfields_get_data src/bindings/bindingsf.f90:966 */
/* Fortran header:
subroutine localfields_get_data(denspotd, rhod, dpbox)
use module_dpbox
use module_types
implicit none
type(DFT_local_fields), intent(in), target :: denspotd
type(denspot_distribution), pointer :: dpbox
type(rho_descriptors), pointer :: rhod
*/
void FC_FUNC_(localfields_get_data, LOCALFIELDS_GET_DATA)(const f90_DFT_local_fields *denspotd, 
                                                          f90_rho_descriptors_pointer *rhod, 
                                                          f90_denspot_distribution_pointer *dpbox);
/* localfields_get_pkernel src/bindings/bindingsf.f90:1064 */
/* Fortran header:
subroutine localfields_get_pkernel(denspot, pkernel)
use module_types
implicit none
type(DFT_local_fields), intent(in), target :: denspot
type(coulomb_operator), pointer :: pkernel
*/
void FC_FUNC_(localfields_get_pkernel, LOCALFIELDS_GET_PKERNEL)(const f90_DFT_local_fields *denspot, 
                                                                f90_coulomb_operator_pointer *pkernel);
/* localfields_get_pkernelseq src/bindings/bindingsf.f90:1074 */
/* Fortran header:
subroutine localfields_get_pkernelseq(denspot, pkernelseq)
use module_types
implicit none
type(DFT_local_fields), intent(in), target :: denspot
type(coulomb_operator), pointer :: pkernelseq
*/
void FC_FUNC_(localfields_get_pkernelseq, LOCALFIELDS_GET_PKERNELSEQ)(const f90_DFT_local_fields *denspot, 
                                                                      f90_coulomb_operator_pointer *pkernelseq);
/* localfields_new src/bindings/bindingsf.f90:950 */
/* Fortran header:
subroutine localfields_new(self, denspotd, rhod, dpbox)
use module_dpbox
use module_types
implicit none
integer(kind = 8), intent(in) :: self
type(DFT_local_fields), pointer :: denspotd
type(denspot_distribution), pointer :: dpbox
type(rho_descriptors), pointer :: rhod
*/
void FC_FUNC_(localfields_new, LOCALFIELDS_NEW)(const long *self, 
                                                f90_DFT_local_fields_pointer *denspotd, 
                                                f90_rho_descriptors_pointer *rhod, 
                                                f90_denspot_distribution_pointer *dpbox);
/* lzd_copy_data src/bindings/bindingsf.f90:429 */
/* Fortran header:
subroutine lzd_copy_data(lzd, nlr)
use module_types
implicit none
type(local_zone_descriptors), intent(in) :: lzd
integer, intent(out) :: nlr
*/
void FC_FUNC_(lzd_copy_data, LZD_COPY_DATA)(const f90_local_zone_descriptors *lzd, 
                                            int *nlr);
/* lzd_empty src/bindings/bindingsf.f90:449 */
/* Fortran header:
subroutine lzd_empty(lzd)
use module_types
implicit none
type(local_zone_descriptors), intent(inout) :: lzd
*/
void FC_FUNC_(lzd_empty, LZD_EMPTY)(f90_local_zone_descriptors *lzd);
/* lzd_free src/bindings/bindingsf.f90:439 */
/* Fortran header:
subroutine lzd_free(lzd)
use module_types
implicit none
type(local_zone_descriptors), pointer :: lzd
*/
void FC_FUNC_(lzd_free, LZD_FREE)(f90_local_zone_descriptors_pointer *lzd);
/* lzd_get_data src/bindings/bindingsf.f90:419 */
/* Fortran header:
subroutine lzd_get_data(lzd, glr)
use module_types
implicit none
type(local_zone_descriptors), target, intent(inout) :: lzd
type(locreg_descriptors), pointer :: glr
*/
void FC_FUNC_(lzd_get_data, LZD_GET_DATA)(f90_local_zone_descriptors *lzd, 
                                          f90_locreg_descriptors_pointer *glr);
/* lzd_get_hgrids src/bindings/bindingsf.f90:484 */
/* Fortran header:
subroutine lzd_get_hgrids(Lzd, hgrids)
use module_base
use module_types
implicit none
type(local_zone_descriptors), intent(in) :: Lzd
real(gp), intent(out) :: hgrids(3)
*/
void FC_FUNC_(lzd_get_hgrids, LZD_GET_HGRIDS)(const f90_local_zone_descriptors *Lzd, 
                                              double *hgrids);
/* lzd_get_llr src/bindings/bindingsf.f90:495 */
/* Fortran header:
subroutine lzd_get_llr(Lzd, i, llr)
use module_base
use module_types
implicit none
type(local_zone_descriptors), intent(in) :: Lzd
integer, intent(in) :: i
type(locreg_descriptors), pointer :: llr
*/
void FC_FUNC_(lzd_get_llr, LZD_GET_LLR)(const f90_local_zone_descriptors *Lzd, 
                                        const int *i, 
                                        f90_locreg_descriptors_pointer *llr);
/* lzd_init src/bindings/bindingsf.f90:407 */
/* Fortran header:
subroutine lzd_init(lzd, glr)
use module_types
implicit none
type(local_zone_descriptors), target, intent(inout) :: lzd
type(locreg_descriptors), pointer :: glr
*/
void FC_FUNC_(lzd_init, LZD_INIT)(f90_local_zone_descriptors *lzd, 
                                  f90_locreg_descriptors_pointer *glr);
/* lzd_init_llr src/linear/initAndUtils.f90:481 */
/* Fortran header:
subroutine lzd_init_llr(iproc, nproc, input, astruct, rxyz, orbs, lzd)
use module_base
use module_types
use locregs, only: locreg_null
implicit none


integer, intent(in) :: iproc, nproc
type(input_variables), intent(in) :: input
type(atomic_structure), intent(in) :: astruct
real(kind=8),dimension(3,astruct%nat), intent(in) :: rxyz
type(orbitals_data), intent(in) :: orbs
type(local_zone_descriptors), intent(inout) :: lzd


integer :: iat, ityp, ilr, istat, iorb, iilr
real(kind=8),dimension(:,:), allocatable :: locregCenter
character(len=*), parameter :: subname='lzd_init_llr'
real(8):: t1, t2
*/
void FC_FUNC_(lzd_init_llr, LZD_INIT_LLR)(const int *iproc, 
                                          const int *nproc, 
                                          const f90_input_variables *input, 
                                          const f90_atomic_structure *astruct, 
                                          const double *rxyz, 
                                          const f90_orbitals_data *orbs, 
                                          f90_local_zone_descriptors *lzd);
/* lzd_new src/bindings/bindingsf.f90:398 */
/* Fortran header:
subroutine lzd_new(lzd)
use module_types
implicit none
type(local_zone_descriptors), pointer :: lzd
*/
void FC_FUNC_(lzd_new, LZD_NEW)(f90_local_zone_descriptors_pointer *lzd);
/* lzd_set_hgrids src/init/wavefunctions.f90:422 */
/* Fortran header:
subroutine lzd_set_hgrids(Lzd, hgrids)
use module_base
use module_types
implicit none
type(local_zone_descriptors), intent(inout) :: Lzd
real(gp), intent(in) :: hgrids(3)
*/
void FC_FUNC_(lzd_set_hgrids, LZD_SET_HGRIDS)(f90_local_zone_descriptors *Lzd, 
                                              const double *hgrids);
/* lzd_set_nlr src/bindings/bindingsf.f90:458 */
/* Fortran header:
subroutine lzd_set_nlr(lzd, nlr, geocode)
use locregs
use module_types, only: local_zone_descriptors
implicit none
type(local_zone_descriptors), intent(inout) :: lzd
integer, intent(in) :: nlr
character, intent(in) :: geocode 

integer :: i
*/
void FC_FUNC_(lzd_set_nlr, LZD_SET_NLR)(f90_local_zone_descriptors *lzd, 
                                        const int *nlr, 
                                        const char *geocode, 
                                        int str_ln_1);
/* mem_destroy src/bindings/bindingsf.f90:1736 */
/* Fortran header:
subroutine mem_destroy(mem)
use module_types, only: memory_estimation
implicit none
type(memory_estimation), pointer :: mem
*/
void FC_FUNC_(mem_destroy, MEM_DESTROY)(f90_memory_estimation_pointer *mem);
/* mem_new src/bindings/bindingsf.f90:1727 */
/* Fortran header:
subroutine mem_new(mem)
use module_types, only: memory_estimation
implicit none
type(memory_estimation), pointer :: mem
*/
void FC_FUNC_(mem_new, MEM_NEW)(f90_memory_estimation_pointer *mem);
/* memoryestimator src/profiling/memoryestimator.f90:12 */
/* Fortran header:
subroutine MemoryEstimator(nproc,idsx,lr,norb,nspinor,nkpt,nprojel,nspin,itrpmax,iscf,mem)

use module_base
use module_types
use Poisson_Solver
use locreg_operations, only: memspace_work_arrays_sumrho,memspace_work_arrays_locham
implicit none


integer, intent(in) :: nproc,idsx,norb,nspin,nprojel
integer, intent(in) :: nkpt,nspinor,itrpmax,iscf
type(locreg_descriptors), intent(in) :: lr
type(memory_estimation), intent(out) :: mem


real(kind=8), parameter :: eps_mach=1.d-12
integer :: norbp,nvctrp,n1,n2,n3
integer :: n01,n02,n03,m1,m2,m3,md1,md2,md3,nd1,nd2,nd3
integer(kind=8) :: mworkham, mworkrho
real(kind=8) :: omemwf,omemker,omemden,omempot,omemproj,nden,npotden,npotham,narr
real(kind=8) :: tt
*/
void FC_FUNC(memoryestimator, MEMORYESTIMATOR)(const int *nproc, 
                                               const int *idsx, 
                                               const f90_locreg_descriptors *lr, 
                                               const int *norb, 
                                               const int *nspinor, 
                                               const int *nkpt, 
                                               const int *nprojel, 
                                               const int *nspin, 
                                               const int *itrpmax, 
                                               const int *iscf, 
                                               f90_memory_estimation *mem);
/* mem_to_c src/bindings/bindingsf.f90:1746 */
/* Fortran header:
subroutine mem_to_c(mem, submat, ncomponents, norb, norbp, oneorb, allpsi_mpi,  psistorage, projarr, grid, workarr, kernel, density, psolver, ham, peak)
use module_types, only: memory_estimation
implicit none
type(memory_estimation), intent(in) :: mem
double precision, intent(out) :: submat, oneorb, allpsi_mpi,  psistorage, projarr, grid, workarr, kernel, density, psolver, ham, peak
integer, intent(out) :: ncomponents, norb, norbp
*/
void FC_FUNC_(mem_to_c, MEM_TO_C)(const f90_memory_estimation *mem, 
                                  double *submat, 
                                  int *ncomponents, 
                                  int *norb, 
                                  int *norbp, 
                                  double *oneorb, 
                                  double *allpsi_mpi, 
                                  double *psistorage, 
                                  double *projarr, 
                                  double *grid, 
                                  double *workarr, 
                                  double *kernel, 
                                  double *density, 
                                  double *psolver, 
                                  double *ham, 
                                  double *peak);
/* optloop_bcast src/bindings/bindingsf.f90:1510 */
/* Fortran header:
subroutine optloop_bcast(optloop, iproc)
use module_base
use module_types
use module_input_keys, only: set_iscf,set_scf_mode
implicit none
type(DFT_optimization_loop), intent(inout) :: optloop
integer, intent(in) :: iproc

integer, dimension(4) :: iData
real(gp), dimension(3) :: rData
*/
void FC_FUNC_(optloop_bcast, OPTLOOP_BCAST)(f90_DFT_optimization_loop *optloop, 
                                            const int *iproc);
/* optloop_copy_data src/bindings/bindingsf.f90:1401 */
/* Fortran header:
subroutine optloop_copy_data(optloop, gnrm_cv, rpnrm_cv, gnrm_startmix, gnrm, rpnrm,   itrpmax, nrepmax, itermax, itrp, itrep, iter, iscf, infocode)
use module_defs, only: gp
use module_types
use module_input_keys, only: set_iscf
implicit none
type(DFT_optimization_loop), intent(in) :: optloop
integer, intent(out) :: iscf, itrpmax, nrepmax, itermax, itrp, itrep, iter, infocode
real(gp), intent(out) :: gnrm, rpnrm, gnrm_cv, rpnrm_cv, gnrm_startmix
*/
void FC_FUNC_(optloop_copy_data, OPTLOOP_COPY_DATA)(const f90_DFT_optimization_loop *optloop, 
                                                    double *gnrm_cv, 
                                                    double *rpnrm_cv, 
                                                    double *gnrm_startmix, 
                                                    double *gnrm, 
                                                    double *rpnrm, 
                                                    int *itrpmax, 
                                                    int *nrepmax, 
                                                    int *itermax, 
                                                    int *itrp, 
                                                    int *itrep, 
                                                    int *iter, 
                                                    int *iscf, 
                                                    int *infocode);
/* optloop_free src/bindings/bindingsf.f90:1392 */
/* Fortran header:
subroutine optloop_free(optloop)
use module_types
implicit none
type(DFT_optimization_loop), pointer :: optloop
*/
void FC_FUNC_(optloop_free, OPTLOOP_FREE)(f90_DFT_optimization_loop_pointer *optloop);
/* optloop_new src/bindings/bindingsf.f90:1381 */
/* Fortran header:
subroutine optloop_new(self, optloop)
use module_types
implicit none
integer(kind = 8), intent(in) :: self
type(DFT_optimization_loop), pointer :: optloop
*/
void FC_FUNC_(optloop_new, OPTLOOP_NEW)(const long *self, 
                                        f90_DFT_optimization_loop_pointer *optloop);
/* optloop_sync_data src/bindings/bindingsf.f90:1429 */
/* Fortran header:
subroutine optloop_sync_data(optloop, gnrm_cv, rpnrm_cv, gnrm_startmix, gnrm, rpnrm,   itrpmax, nrepmax, itermax, itrp, itrep, iter, iscf, infocode)
use module_defs, only: gp
use module_types
use module_input_keys, only: set_scf_mode
implicit none
type(DFT_optimization_loop), intent(inout) :: optloop
integer, intent(in) :: iscf, itrpmax, nrepmax, itermax, itrp, itrep, iter, infocode
real(gp), intent(in) :: gnrm, rpnrm, gnrm_cv, rpnrm_cv, gnrm_startmix
*/
void FC_FUNC_(optloop_sync_data, OPTLOOP_SYNC_DATA)(f90_DFT_optimization_loop *optloop, 
                                                    const double *gnrm_cv, 
                                                    const double *rpnrm_cv, 
                                                    const double *gnrm_startmix, 
                                                    const double *gnrm, 
                                                    const double *rpnrm, 
                                                    const int *itrpmax, 
                                                    const int *nrepmax, 
                                                    const int *itermax, 
                                                    const int *itrp, 
                                                    const int *itrep, 
                                                    const int *iter, 
                                                    const int *iscf, 
                                                    const int *infocode);
/* orbs_comm_empty src/bindings/bindingsf.f90:784 */
/* Fortran header:
subroutine orbs_comm_empty(comms)
use module_base
use module_types
use communications_base, only: comms_cubic, deallocate_comms

implicit none
type(comms_cubic), intent(inout) :: comms
*/
void FC_FUNC_(orbs_comm_empty, ORBS_COMM_EMPTY)(f90_comms_cubic *comms);
/* orbs_comm_free src/bindings/bindingsf.f90:772 */
/* Fortran header:
subroutine orbs_comm_free(comms)
use module_base
use module_types
use communications_base, only: comms_cubic

implicit none
type(comms_cubic), pointer :: comms
*/
void FC_FUNC_(orbs_comm_free, ORBS_COMM_FREE)(f90_comms_cubic_pointer *comms);
/* orbs_comm_init src/bindings/bindingsf.f90:757 */
/* Fortran header:
subroutine orbs_comm_init(comms, orbs, lr, iproc, nproc)
use module_base
use module_types
use communications_base, only: comms_cubic
use communications_init, only: orbitals_communicators
implicit none
integer, intent(in) :: iproc,nproc
type(locreg_descriptors), intent(in) :: lr
type(orbitals_data), intent(inout) :: orbs
type(comms_cubic), intent(inout) :: comms
*/
void FC_FUNC_(orbs_comm_init, ORBS_COMM_INIT)(f90_comms_cubic *comms, 
                                              f90_orbitals_data *orbs, 
                                              const f90_locreg_descriptors *lr, 
                                              const int *iproc, 
                                              const int *nproc);
/* orbs_comm_new src/bindings/bindingsf.f90:744 */
/* Fortran header:
subroutine orbs_comm_new(comms)
use module_base
use module_types
use communications_base, only: comms_cubic

implicit none
type(comms_cubic), pointer :: comms
*/
void FC_FUNC_(orbs_comm_new, ORBS_COMM_NEW)(f90_comms_cubic_pointer *comms);
/* orbs_empty src/bindings/bindingsf.f90:735 */
/* Fortran header:
subroutine orbs_empty(orbs)
use module_types
implicit none
type(orbitals_data), intent(inout) :: orbs
*/
void FC_FUNC_(orbs_empty, ORBS_EMPTY)(f90_orbitals_data *orbs);
/* orbs_free src/bindings/bindingsf.f90:725 */
/* Fortran header:
subroutine orbs_free(orbs)
use module_types

implicit none
type(orbitals_data), pointer :: orbs
*/
void FC_FUNC_(orbs_free, ORBS_FREE)(f90_orbitals_data_pointer *orbs);
/* orbs_get_dimensions src/bindings/bindingsf.f90:798 */
/* Fortran header:
subroutine orbs_get_dimensions(orbs, norb, norbp, norbu, norbd, nspin, nspinor, npsidim,  nkpts, nkptsp, isorb, iskpts)
use module_types
implicit none
type(orbitals_data), intent(in) :: orbs
integer, intent(out) :: norb, norbp, norbu, norbd, nspin, nspinor, npsidim,  nkpts, nkptsp, isorb, iskpts
*/
void FC_FUNC_(orbs_get_dimensions, ORBS_GET_DIMENSIONS)(const f90_orbitals_data *orbs, 
                                                        int *norb, 
                                                        int *norbp, 
                                                        int *norbu, 
                                                        int *norbd, 
                                                        int *nspin, 
                                                        int *nspinor, 
                                                        int *npsidim, 
                                                        int *nkpts, 
                                                        int *nkptsp, 
                                                        int *isorb, 
                                                        int *iskpts);
/* orbs_get_iorbp src/bindings/bindingsf.f90:1253 */
/* Fortran header:
subroutine orbs_get_iorbp(orbs, iorbp, isorb, iproc, ikpt, iorb, ispin, ispinor)
use module_types
implicit none

integer, intent(out) :: iorbp, isorb, iproc
type(orbitals_data), intent(in) :: orbs
integer, intent(in) :: ikpt, iorb, ispin, ispinor
*/
void FC_FUNC_(orbs_get_iorbp, ORBS_GET_IORBP)(const f90_orbitals_data *orbs, 
                                              int *iorbp, 
                                              int *isorb, 
                                              int *iproc, 
                                              const int *ikpt, 
                                              const int *iorb, 
                                              const int *ispin, 
                                              const int *ispinor);
/* orbs_init src/bindings/bindingsf.f90:716 */
/* Fortran header:
subroutine orbs_init(orbs)
use module_types
implicit none
type(orbitals_data), intent(inout) :: orbs
*/
void FC_FUNC_(orbs_init, ORBS_INIT)(f90_orbitals_data *orbs);
/* orbs_new src/bindings/bindingsf.f90:707 */
/* Fortran header:
subroutine orbs_new(orbs)
use module_types
implicit none
type(orbitals_data), pointer :: orbs
*/
void FC_FUNC_(orbs_new, ORBS_NEW)(f90_orbitals_data_pointer *orbs);
/* orbs_open_file src/bindings/bindingsf.f90:884 */
/* Fortran header:
subroutine orbs_open_file(orbs, unitwf, name, ln, iformat, iorbp, ispinor)
use module_types
use public_enums
use module_interfaces, only: open_filename_of_iorb
implicit none
type(orbitals_data), intent(in) :: orbs
integer, intent(in) :: ln, iformat, iorbp, ispinor
character(len = 1), dimension(ln), intent(in) :: name
integer, intent(inout) :: unitwf

character(len = ln) :: filename
integer :: i, iorb_out
*/
void FC_FUNC_(orbs_open_file, ORBS_OPEN_FILE)(const f90_orbitals_data *orbs, 
                                              int *unitwf, 
                                              const char *name, 
                                              const int *ln, 
                                              const int *iformat, 
                                              const int *iorbp, 
                                              const int *ispinor, 
                                              int str_ln_1);
/* print_memory_estimation src/output.f90:1134 */
/* Fortran header:
subroutine print_memory_estimation(mem)
use module_types
use yaml_output
use yaml_strings
implicit none
type(memory_estimation), intent(in) :: mem
*/
void FC_FUNC_(print_memory_estimation, PRINT_MEMORY_ESTIMATION)(const f90_memory_estimation *mem);
/* proj_free src/bindings/bindingsf.f90:914 */
/* Fortran header:
subroutine proj_free(nlpspd)
use psp_projectors_base, only: free_DFT_PSP_projectors
use module_types
use memory_profiling
implicit none
type(DFT_PSP_projectors), pointer :: nlpspd
*/
void FC_FUNC_(proj_free, PROJ_FREE)(f90_DFT_PSP_projectors_pointer *nlpspd);
/* proj_get_dimensions src/bindings/bindingsf.f90:925 */
/* Fortran header:
subroutine proj_get_dimensions(nlpspd, nproj, nprojel)
use module_types
implicit none
type(DFT_PSP_projectors), intent(in) :: nlpspd
integer, intent(out) :: nproj, nprojel
*/
void FC_FUNC_(proj_get_dimensions, PROJ_GET_DIMENSIONS)(const f90_DFT_PSP_projectors *nlpspd, 
                                                        int *nproj, 
                                                        int *nprojel);
/* proj_new src/bindings/bindingsf.f90:905 */
/* Fortran header:
subroutine proj_new(nlpspd)
use module_types
implicit none
type(nonlocal_psp_descriptors), pointer :: nlpspd
*/
void FC_FUNC_(proj_new, PROJ_NEW)(f90_nonlocal_psp_descriptors_pointer *nlpspd);
/* read_orbital_variables  */
/* read_wave_descr src/restart.f90:674 */
/* Fortran header:
subroutine read_wave_descr(lstat, filename, ln,  norbu, norbd, iorbs, ispins, nkpt, ikpts, nspinor, ispinor)
use public_enums 
use module_input_keys
implicit none
integer, intent(in) :: ln
character, intent(in) :: filename(ln)
integer, intent(out) :: norbu, norbd, nkpt, nspinor
integer, intent(out) :: iorbs, ispins, ikpts, ispinor
logical, intent(out) :: lstat

character(len = 1024) :: filename_
integer :: iformat, i
character(len = 1024) :: testf
*/
void FC_FUNC_(read_wave_descr, READ_WAVE_DESCR)(int *lstat, 
                                                const char *filename, 
                                                const int *ln, 
                                                int *norbu, 
                                                int *norbd, 
                                                int *iorbs, 
                                                int *ispins, 
                                                int *nkpt, 
                                                int *ikpts, 
                                                int *nspinor, 
                                                int *ispinor, 
                                                int str_ln_1);
/* read_wave_to_isf src/restart.f90:622 */
/* Fortran header:
subroutine read_wave_to_isf(lstat, filename, ln, iorbp, hx, hy, hz,  n1, n2, n3, nspinor, psiscf)
use module_base
use module_types
use module_interfaces, only: readwavetoisf, readwavetoisf_etsf
use public_enums
use module_input_keys
implicit none

integer, intent(in) :: ln
character, intent(in) :: filename(ln)
integer, intent(in) :: iorbp
integer, intent(out) :: n1, n2, n3, nspinor
real(gp), intent(out) :: hx, hy, hz
real(wp), dimension(:,:,:,:), pointer :: psiscf
logical, intent(out) :: lstat

character(len = 1024) :: filename_
integer :: iformat, i
*/
void FC_FUNC_(read_wave_to_isf, READ_WAVE_TO_ISF)(int *lstat, 
                                                  const char *filename, 
                                                  const int *ln, 
                                                  const int *iorbp, 
                                                  double *hx, 
                                                  double *hy, 
                                                  double *hz, 
                                                  int *n1, 
                                                  int *n2, 
                                                  int *n3, 
                                                  int *nspinor, 
                                                  f90_pointer_double_4D *psiscf, 
                                                  int str_ln_1);
/* run_objects_c_obj src/bindings/bindingsf.f90:1602 */
/* Fortran header:
subroutine run_objects_c_obj(run, c_obj)
use bigdft_run
implicit none
type(run_objects), intent(in) :: run
integer(kind = 8), intent(out) :: c_obj
*/
void FC_FUNC_(run_objects_c_obj, RUN_OBJECTS_C_OBJ)(const f90_run_objects *run, 
                                                    long *c_obj);
/* run_objects_copy src/bindings/bindingsf.f90:1570 */
/* Fortran header:
subroutine run_objects_copy(run, from)
use bigdft_run
implicit none
type(run_objects), pointer :: run
type(run_objects), intent(in) :: from

type(run_objects), pointer :: intern
*/
void FC_FUNC_(run_objects_copy, RUN_OBJECTS_COPY)(f90_run_objects_pointer *run, 
                                                  const f90_run_objects *from);
/* run_objects_destroy src/bindings/bindingsf.f90:1584 */
/* Fortran header:
subroutine run_objects_destroy(runObj)
use bigdft_run
implicit none
type(run_objects), pointer :: runObj
*/
void FC_FUNC_(run_objects_destroy, RUN_OBJECTS_DESTROY)(f90_run_objects_pointer *runObj);
/* run_objects_dump_to_file src/bindings/bindingsf.f90:1628 */
/* Fortran header:
subroutine run_objects_dump_to_file(iostat, dict, fname, userOnly,ln)
use dictionaries, only: dictionary
use module_input_keys, only: input_keys_dump
use module_defs, only: UNINITIALIZED, gp
use yaml_output
use f_utils, only: f_get_free_unit
use yaml_strings, only: f_strcpy
implicit none
integer, intent(in) :: ln
integer, intent(out) :: iostat
type(dictionary), pointer :: dict
character, dimension(ln), intent(in) :: fname
logical, intent(in) :: userOnly

integer, parameter :: iunit_true = 145214 
integer :: iunit_def,iunit,iln
real(gp), dimension(3), parameter :: dummy = (/ 0._gp, 0._gp, 0._gp /)
character(len=256) :: filetmp
*/
void FC_FUNC_(run_objects_dump_to_file, RUN_OBJECTS_DUMP_TO_FILE)(int *iostat, 
                                                                  f90_dictionary_pointer *dict, 
                                                                  const char *fname, 
                                                                  const int *userOnly, 
                                                                  const int *ln, 
                                                                  int str_ln_1);
/* run_objects_get src/bindings/bindingsf.f90:1611 */
/* Fortran header:
subroutine run_objects_get(runObj, dict, inputs, atoms)
use bigdft_run
use module_types, only: input_variables
use module_atoms, only: atoms_data
use dictionaries
implicit none
type(run_objects), intent(in) :: runObj
type(dictionary), pointer :: dict
type(input_variables), pointer :: inputs
type(atoms_data), pointer :: atoms
*/
void FC_FUNC_(run_objects_get, RUN_OBJECTS_GET)(const f90_run_objects *runObj, 
                                                f90_dictionary_pointer *dict, 
                                                f90_input_variables_pointer *inputs, 
                                                f90_atoms_data_pointer *atoms);
/* run_objects_init_from_run_name src/modules/bigdft_run.f90:2480 */
/* Fortran header:
subroutine run_objects_init_from_run_name(runObj, radical, posinp)
use module_base
use bigdft_run
implicit none
type(run_objects), intent(out) :: runObj
character(len = *), intent(in) :: radical, posinp

type(dictionary), pointer :: run_dict
*/
void FC_FUNC_(run_objects_init_from_run_name, RUN_OBJECTS_INIT_FROM_RUN_NAME)(f90_run_objects *runObj, 
                                                                              const char *radical, 
                                                                              const char *posinp, 
                                                                              int str_ln_1, 
                                                                              int str_ln_2);
/* run_objects_new src/bindings/bindingsf.f90:1551 */
/* Fortran header:
subroutine run_objects_new(runObj)
use bigdft_run
implicit none
type(run_objects), pointer :: runObj

type(run_objects), pointer :: intern
*/
void FC_FUNC_(run_objects_new, RUN_OBJECTS_NEW)(f90_run_objects_pointer *runObj);
/* run_objects_nullify_dict src/bindings/bindingsf.f90:1697 */
/* Fortran header:
subroutine run_objects_nullify_dict(runObj)
use bigdft_run, only: run_objects
implicit none
type(run_objects), intent(inout) :: runObj
*/
void FC_FUNC_(run_objects_nullify_dict, RUN_OBJECTS_NULLIFY_DICT)(f90_run_objects *runObj);
/* run_objects_nullify_volatile src/bindings/bindingsf.f90:1706 */
/* Fortran header:
subroutine run_objects_nullify_volatile(runObj)
use f_enums
use public_enums
use bigdft_run, only: run_objects
use module_defs, only: verbose
use yaml_output, only: yaml_sequence_close
use module_base, only: bigdft_mpi
implicit none
type(run_objects), intent(inout) :: runObj
*/
void FC_FUNC_(run_objects_nullify_volatile, RUN_OBJECTS_NULLIFY_VOLATILE)(f90_run_objects *runObj);
/* run_objects_set_c_obj src/bindings/bindingsf.f90:1594 */
/* Fortran header:
subroutine run_objects_set_c_obj(run, c_obj)
use bigdft_run
implicit none
type(run_objects), intent(inout) :: run
integer(kind = 8), intent(in) :: c_obj
*/
void FC_FUNC_(run_objects_set_c_obj, RUN_OBJECTS_SET_C_OBJ)(f90_run_objects *run, 
                                                            const long *c_obj);
/* run_objects_system_setup src/modules/bigdft_run.f90:2538 */
/* Fortran header:
subroutine run_objects_system_setup(runObj, iproc, nproc, rxyz, shift, mem)
use module_base, only: gp,f_memcpy,f_enumerator,f_int
use bigdft_run
use module_types
use module_fragments
use module_interfaces, only: system_initialization
use psp_projectors_base, only: free_DFT_PSP_projectors
use communications_base, only: deallocate_comms
implicit none
type(run_objects), intent(inout) :: runObj
integer, intent(in) :: iproc, nproc
real(gp), dimension(3,runObj%atoms%astruct%nat), intent(out) :: rxyz
real(gp), dimension(3), intent(out) :: shift
type(memory_estimation), intent(out) :: mem

integer :: input_wf_format
type(DFT_PSP_projectors) :: nlpsp
type(f_enumerator) :: inputpsi
type(system_fragment), dimension(:), pointer :: ref_frags
character(len = *), parameter :: subname = "run_objects_estimate_memory"
*/
void FC_FUNC_(run_objects_system_setup, RUN_OBJECTS_SYSTEM_SETUP)(f90_run_objects *runObj, 
                                                                  const int *iproc, 
                                                                  const int *nproc, 
                                                                  double *rxyz, 
                                                                  double *shift, 
                                                                  f90_memory_estimation *mem);
/* run_objects_update src/modules/bigdft_run.f90:2500 */
/* Fortran header:
subroutine run_objects_update(runObj, dict)
use module_base, only: bigdft_mpi
use bigdft_run, only: run_objects,init_QM_restart_objects,init_MM_restart_objects,set_run_objects,bigdft_nat
use dictionaries
use yaml_output
use module_input_dicts, only: create_log_file
implicit none
type(run_objects), intent(inout) :: runObj
type(dictionary), pointer :: dict

type(dictionary), pointer :: item
logical :: dict_from_files
*/
void FC_FUNC_(run_objects_update, RUN_OBJECTS_UPDATE)(f90_run_objects *runObj, 
                                                      f90_dictionary_pointer *dict);
/* state_properties_alloc src/bindings/bindingsf.f90:1281 */
/* Fortran header:
subroutine state_properties_alloc(outs, nat)
use bigdft_run
implicit none
type(state_properties), pointer :: outs
integer, intent(in) :: nat

type(state_properties), pointer :: intern
*/
void FC_FUNC_(state_properties_alloc, STATE_PROPERTIES_ALLOC)(f90_state_properties_pointer *outs, 
                                                              const int *nat);
/* state_properties_c_obj src/bindings/bindingsf.f90:1346 */
/* Fortran header:
subroutine state_properties_c_obj(outs, c_obj)
use bigdft_run
implicit none
type(state_properties), intent(in) :: outs
integer(kind = 8), intent(out) :: c_obj
*/
void FC_FUNC_(state_properties_c_obj, STATE_PROPERTIES_C_OBJ)(const f90_state_properties *outs, 
                                                              long *c_obj);
/* state_properties_copy src/bindings/bindingsf.f90:1293 */
/* Fortran header:
subroutine state_properties_copy(outs, from)
use bigdft_run
implicit none
type(state_properties), pointer :: outs
type(state_properties), intent(in) :: from

type(state_properties), pointer :: intern
*/
void FC_FUNC_(state_properties_copy, STATE_PROPERTIES_COPY)(f90_state_properties_pointer *outs, 
                                                            const f90_state_properties *from);
/* state_properties_delete src/bindings/bindingsf.f90:1305 */
/* Fortran header:
subroutine state_properties_delete(outs)
use bigdft_run
implicit none
type(state_properties), pointer :: outs
*/
void FC_FUNC_(state_properties_delete, STATE_PROPERTIES_DELETE)(f90_state_properties_pointer *outs);
/* state_properties_get src/bindings/bindingsf.f90:1314 */
/* Fortran header:
subroutine state_properties_get(outs, energs, fxyz, fdim, fnoise, pressure, strten, etot)
use module_defs, only: gp
use module_types, only: energy_terms
use bigdft_run
implicit none
type(state_properties), intent(in), target :: outs
type(energy_terms), pointer :: energs
real(gp), dimension(:,:), pointer :: fxyz
integer, intent(out) :: fdim
real(gp), intent(out) :: fnoise, pressure
real(gp), dimension(6), intent(out) :: strten
real(gp), intent(out) :: etot
*/
void FC_FUNC_(state_properties_get, STATE_PROPERTIES_GET)(const f90_state_properties *outs, 
                                                          f90_energy_terms_pointer *energs, 
                                                          f90_pointer_double_2D *fxyz, 
                                                          int *fdim, 
                                                          double *fnoise, 
                                                          double *pressure, 
                                                          double *strten, 
                                                          double *etot);
/* state_properties_set_c_obj src/bindings/bindingsf.f90:1338 */
/* Fortran header:
subroutine state_properties_set_c_obj(outs, c_obj)
use bigdft_run
implicit none
type(state_properties), intent(inout) :: outs
integer(kind = 8), intent(in) :: c_obj
*/
void FC_FUNC_(state_properties_set_c_obj, STATE_PROPERTIES_SET_C_OBJ)(f90_state_properties *outs, 
                                                                      const long *c_obj);
/* symmetry_set_irreductible_zone  */
/* system_createkernels src/init/sysprop.f90:861 */
/* Fortran header:
subroutine system_createKernels(denspot, verb)
use module_base
use module_types
use Poisson_Solver, except_dp => dp, except_gp => gp
implicit none
logical, intent(in) :: verb
integer(kind=8)  :: iproc_node, nproc_node
type(DFT_local_fields), intent(inout) :: denspot
*/
void FC_FUNC_(system_createkernels, SYSTEM_CREATEKERNELS)(f90_DFT_local_fields *denspot, 
                                                          const int *verb);
/* system_initkernels src/init/sysprop.f90:824 */
/* Fortran header:
subroutine system_initKernels(verb, iproc, nproc, geocode, in, denspot)
use module_types
use module_xc
use Poisson_Solver, except_dp => dp, except_gp => gp
use module_base
implicit none
logical, intent(in) :: verb
integer, intent(in) :: iproc, nproc
character, intent(in) :: geocode 
type(input_variables), intent(in) :: in
type(DFT_local_fields), intent(inout) :: denspot

integer, parameter :: ndegree_ip = 16
*/
void FC_FUNC_(system_initkernels, SYSTEM_INITKERNELS)(const int *verb, 
                                                      const int *iproc, 
                                                      const int *nproc, 
                                                      const char *geocode, 
                                                      const f90_input_variables *in, 
                                                      f90_DFT_local_fields *denspot, 
                                                      int str_ln_1);
/* system_size src/init/gridmanipulation.f90:14 */
/* Fortran header:
subroutine system_size(atoms,rxyz,crmult,frmult,hx,hy,hz,OCLconv,Glr,shift)
use module_base
use module_types
use yaml_strings, only: yaml_toa
implicit none
type(atoms_data), intent(inout) :: atoms
real(gp), intent(in) :: crmult,frmult
real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz

real(gp), intent(inout) :: hx,hy,hz
logical, intent(in) :: OCLconv
type(locreg_descriptors), intent(out) :: Glr
real(gp), dimension(3), intent(out) :: shift


integer, parameter :: lupfil=14
real(gp), parameter :: eps_mach=1.e-12_gp
integer :: iat,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i
real(gp) :: ri,rad,cxmin,cxmax,cymin,cymax,czmin,czmax,alatrue1,alatrue2,alatrue3
*/
void FC_FUNC_(system_size, SYSTEM_SIZE)(f90_atoms_data *atoms, 
                                        double *rxyz, 
                                        const double *crmult, 
                                        const double *frmult, 
                                        double *hx, 
                                        double *hy, 
                                        double *hz, 
                                        const int *OCLconv, 
                                        f90_locreg_descriptors *Glr, 
                                        double *shift);
/* update_wavefunctions_size src/linear/initAndUtils.f90:815 */
/* Fortran header:
subroutine update_wavefunctions_size(lzd,npsidim_orbs,npsidim_comp,orbs,iproc,nproc)
use module_base
use module_types
implicit none


type(local_zone_descriptors), intent(in) :: lzd
type(orbitals_data), intent(in) :: orbs
integer, intent(in) :: iproc, nproc
integer, intent(out) :: npsidim_orbs, npsidim_comp


character(len = *), parameter :: subname = "update_wavefunctions_size"
integer :: npsidim, ilr, iorb
integer :: nvctr_tot,jproc
integer, allocatable, dimension(:) :: ncntt 
integer, allocatable, dimension(:,:) :: nvctr_par
*/
void FC_FUNC_(update_wavefunctions_size, UPDATE_WAVEFUNCTIONS_SIZE)(const f90_local_zone_descriptors *lzd, 
                                                                    int *npsidim_orbs, 
                                                                    int *npsidim_comp, 
                                                                    const f90_orbitals_data *orbs, 
                                                                    const int *iproc, 
                                                                    const int *nproc);
/* wf_empty src/bindings/bindingsf.f90:1172 */
/* Fortran header:
subroutine wf_empty(wf)
use module_base, only: f_free_ptr
use module_types
use memory_profiling
implicit none
type(DFT_wavefunction), intent(inout) :: wf
*/
void FC_FUNC_(wf_empty, WF_EMPTY)(f90_DFT_wavefunction *wf);
/* wf_free src/bindings/bindingsf.f90:1185 */
/* Fortran header:
subroutine wf_free(wf)
use module_types
use memory_profiling
implicit none
type(DFT_wavefunction), pointer :: wf
*/
void FC_FUNC_(wf_free, WF_FREE)(f90_DFT_wavefunction_pointer *wf);
/* wf_get_data src/bindings/bindingsf.f90:1157 */
/* Fortran header:
subroutine wf_get_data(wf, orbs, comm, lzd)
use module_types
use communications_base, only: comms_cubic
implicit none
type(DFT_wavefunction), target, intent(in) :: wf
type(orbitals_data), pointer :: orbs
type(comms_cubic), pointer :: comm
type(local_zone_descriptors), pointer :: lzd
*/
void FC_FUNC_(wf_get_data, WF_GET_DATA)(const f90_DFT_wavefunction *wf, 
                                        f90_orbitals_data_pointer *orbs, 
                                        f90_comms_cubic_pointer *comm, 
                                        f90_local_zone_descriptors_pointer *lzd);
/* wf_get_psi src/bindings/bindingsf.f90:1198 */
/* Fortran header:
subroutine wf_get_psi(wf, psi, hpsi)
use module_types
implicit none
type(DFT_wavefunction), intent(in) :: wf
integer(kind = 8), intent(out) :: psi
integer(kind = 8), intent(out) :: hpsi
*/
void FC_FUNC_(wf_get_psi, WF_GET_PSI)(const f90_DFT_wavefunction *wf, 
                                      long *psi, 
                                      long *hpsi);
/* wf_get_psi_size src/bindings/bindingsf.f90:1216 */
/* Fortran header:
subroutine wf_get_psi_size(psi, psiSize)
use module_defs, only: wp
use module_types
implicit none
real(wp), dimension(:), pointer :: psi
integer(kind = 8), intent(out) :: psiSize
*/
void FC_FUNC_(wf_get_psi_size, WF_GET_PSI_SIZE)(f90_pointer_double *psi, 
                                                long *psiSize);
/* wf_iorbp_to_psi src/bindings/bindingsf.f90:1227 */
/* Fortran header:
subroutine wf_iorbp_to_psi(psir, psi, lr)
use module_base, only: wp,f_zero
use locregs
use locreg_operations
implicit none
type(locreg_descriptors), intent(in) :: lr
real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i), intent(out) :: psir

character(len=*), parameter :: subname='wf_orb_to_psi'
type(workarr_sumrho) :: w
*/
void FC_FUNC_(wf_iorbp_to_psi, WF_IORBP_TO_PSI)(double *psir, 
                                                const double *psi, 
                                                const f90_locreg_descriptors *lr);
/* wf_new src/bindings/bindingsf.f90:1124 */
/* Fortran header:
subroutine wf_new(self, wf, orbs, comm, lzd)
use module_types
use communications_base, only: comms_cubic
implicit none
integer(kind = 8), intent(in) :: self
type(DFT_wavefunction), pointer :: wf
type(orbitals_data), pointer :: orbs
type(comms_cubic), pointer :: comm
type(local_zone_descriptors), pointer :: lzd
*/
void FC_FUNC_(wf_new, WF_NEW)(const long *self, 
                              f90_DFT_wavefunction_pointer *wf, 
                              f90_orbitals_data_pointer *orbs, 
                              f90_comms_cubic_pointer *comm, 
                              f90_local_zone_descriptors_pointer *lzd);
/* writeonewave src/wavelib/i-o.f90:632 */
/* Fortran header:
subroutine writeonewave(unitwf,useFormattedOutput,iorb,n1,n2,n3,hx,hy,hz,nat,rxyz,  nseg_c,nvctr_c,keyg_c,keyv_c,  nseg_f,nvctr_f,keyg_f,keyv_f, psi_c,psi_f,eval)
use module_base
use yaml_output
implicit none
logical, intent(in) :: useFormattedOutput
integer, intent(inout) :: unitwf,iorb,n1,n2,n3,nat,nseg_c,nvctr_c,nseg_f,nvctr_f
real(gp), intent(in) :: hx,hy,hz
real(wp), intent(in) :: eval
integer, dimension(nseg_c), intent(in) :: keyv_c
integer, dimension(nseg_f), intent(in) :: keyv_f
integer, dimension(2,nseg_c), intent(in) :: keyg_c
integer, dimension(2,nseg_f), intent(in) :: keyg_f
real(wp), dimension(nvctr_c), intent(in) :: psi_c
real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
real(gp), dimension(3,nat), intent(in) :: rxyz

integer :: iat,jj,j0,j1,ii,i0,i1,i2,i3,i,iseg,j
real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7
*/
void FC_FUNC(writeonewave, WRITEONEWAVE)(int *unitwf, 
                                         const int *useFormattedOutput, 
                                         int *iorb, 
                                         int *n1, 
                                         int *n2, 
                                         int *n3, 
                                         const double *hx, 
                                         const double *hy, 
                                         const double *hz, 
                                         int *nat, 
                                         const double *rxyz, 
                                         int *nseg_c, 
                                         int *nvctr_c, 
                                         const int *keyg_c, 
                                         const int *keyv_c, 
                                         int *nseg_f, 
                                         int *nvctr_f, 
                                         const int *keyg_f, 
                                         const int *keyv_f, 
                                         const double *psi_c, 
                                         const double *psi_f, 
                                         const double *eval);
/* writeonewave_linear src/modules/io.f90:901 */
/* Fortran header:
subroutine writeonewave_linear(unitwf,useFormattedOutput,iorb,n1,n2,n3,ns1,ns2,ns3,hx,hy,hz,locregCenter,locrad,confPotOrder,confPotprefac,nat,rxyz, nseg_c,nvctr_c,keyg_c,keyv_c,  nseg_f,nvctr_f,keyg_f,keyv_f, psi_c,psi_f,eval,onwhichatom)
use module_base
use yaml_output
implicit none
logical, intent(in) :: useFormattedOutput
integer, intent(in) :: unitwf,iorb,n1,n2,n3,ns1,ns2,ns3,nat,nseg_c,nvctr_c,nseg_f,nvctr_f,confPotOrder
real(gp), intent(in) :: hx,hy,hz,locrad,confPotprefac
real(wp), intent(in) :: eval
integer, dimension(nseg_c), intent(in) :: keyv_c
integer, dimension(nseg_f), intent(in) :: keyv_f
integer, dimension(2,nseg_c), intent(in) :: keyg_c
integer, dimension(2,nseg_f), intent(in) :: keyg_f
real(wp), dimension(nvctr_c), intent(in) :: psi_c
real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
real(gp), dimension(3,nat), intent(in) :: rxyz
real(gp), dimension(3), intent(in) :: locregCenter
integer, intent(in) :: onwhichatom

integer :: iat,jj,j0,j1,ii,i0,i1,i2,i3,i,iseg,j,np,n1p1
real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7
*/
void FC_FUNC_(writeonewave_linear, WRITEONEWAVE_LINEAR)(const int *unitwf, 
                                                        const int *useFormattedOutput, 
                                                        const int *iorb, 
                                                        const int *n1, 
                                                        const int *n2, 
                                                        const int *n3, 
                                                        const int *ns1, 
                                                        const int *ns2, 
                                                        const int *ns3, 
                                                        const double *hx, 
                                                        const double *hy, 
                                                        const double *hz, 
                                                        const double *locregCenter, 
                                                        const double *locrad, 
                                                        const int *confPotOrder, 
                                                        const double *confPotprefac, 
                                                        const int *nat, 
                                                        const double *rxyz, 
                                                        const int *nseg_c, 
                                                        const int *nvctr_c, 
                                                        const int *keyg_c, 
                                                        const int *keyv_c, 
                                                        const int *nseg_f, 
                                                        const int *nvctr_f, 
                                                        const int *keyg_f, 
                                                        const int *keyv_f, 
                                                        const double *psi_c, 
                                                        const double *psi_f, 
                                                        const double *eval, 
                                                        const int *onwhichatom);
#endif
