!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2015-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

module module_mhgps_state
    use module_base, only: gp !bigdft base module
    implicit none

    private

    public :: mhgps_state
    public :: init_mhgps_state
    public :: finalize_mhgps_state

    type mhgps_state
        !values that must not be changed during the execution
        character(len = 5) :: dirprefix
        character(len = 20) :: mhgps_version

        integer :: nid

        !state variables
        character(len=100), allocatable :: joblist(:,:)
        integer            :: njobsmax
        integer            :: njobs
        integer            :: ifolder
        integer            :: ijob
        character(len=8)   :: currDir !length is length(dirprefix)+
                                      !number of digits
        real(gp)           :: ef_counter
        integer            :: nsad
        integer            :: isad
        character(len=5)   :: isadc
        integer            :: nrestart
        integer            :: ntodo
        character(len=5)   :: ntodoc
        integer            :: isadprob
        character(len=5)   :: isadprobc
        integer            :: iproc
        integer            :: nproc
        integer            :: igroup
        integer            :: ngroups

        real(gp), allocatable :: attempted_connections(:,:,:,:)
        integer            :: nattemptedmax
        integer            :: nattempted
    end type

contains

subroutine init_mhgps_state(mhgpsst)
    use module_base
    implicit none
    !parameters
    type(mhgps_state), intent(inout) :: mhgpsst
    include 'mhgps_version_number-inc.f90'
 
    !attention: there exists many dependecies on
    !the exact length on direprefix in the whole code.
    !DO NOT simply change dirprefix 
    mhgpsst%dirprefix     = 'ioput'

    mhgpsst%njobsmax      = 999
    mhgpsst%njobs         = -1
    mhgpsst%ifolder       = -1
    mhgpsst%ijob          = -1
!   mhgpsst%joblist = f_malloc_str((/1.to.2, 1.to.999/),id='joblist') !how??
    allocate(mhgpsst%joblist(2,mhgpsst%njobsmax))
 
    write(mhgpsst%currDir,'(a,i3.3)')&
                     trim(adjustl(mhgpsst%dirprefix)),1
    mhgpsst%ef_counter      = 0.0_gp
    mhgpsst%nsad            = 0 
    mhgpsst%isad            = 0 
    mhgpsst%isadc           = ''
    mhgpsst%nrestart        = 0
    mhgpsst%ntodo           = 0
    mhgpsst%ntodoc          = ''
    mhgpsst%isadprob        = 0
    mhgpsst%isadprobc       = ''
    mhgpsst%iproc=bigdft_mpi%iproc
    mhgpsst%nproc=bigdft_mpi%nproc
    mhgpsst%igroup=bigdft_mpi%igroup
    !number of groups
    mhgpsst%ngroups=bigdft_mpi%ngroup!mpi_info(4)
    !actual value of iproc
    mhgpsst%iproc=mhgpsst%iproc+mhgpsst%igroup*mhgpsst%ngroups
    mhgpsst%nattemptedmax = 0
    mhgpsst%nattempted    = 0
end subroutine
subroutine finalize_mhgps_state(mhgpsst)
    use module_base
    implicit none
    !parameters
    type(mhgps_state), intent(inout) :: mhgpsst

    deallocate(mhgpsst%joblist)
    call f_free(mhgpsst%attempted_connections)
end subroutine

end module
