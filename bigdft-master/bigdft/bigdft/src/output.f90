!>  @file
!!  File where most relevant screen output are collected
!!  Routines which are present in this file should have *all* arguments as intent(in)
!!  Also, the master process only should acces these routines
!! @author
!!    Copyright (C) 2011-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Display the logo of BigDFT subroutine print_logo()
subroutine print_logo()
  use module_base
  use yaml_output
  implicit none
  integer :: namelen,ierr
  character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
  integer :: nthreads
  !integer, parameter :: ln = 1024
!$ integer :: omp_get_max_threads

!  call yaml_comment('Daubechies Wavelets for DFT Pseudopotential Calculations',hfill='=')


  call yaml_mapping_open('Code logo')
!logo of BigDFT, new version
  call yaml_scalar('"__________________________________ A fast and precise DFT wavelet code')
  call yaml_scalar('|     |     |     |     |     |                                        ')
  call yaml_scalar('|     |     |     |     |     |      BBBB         i       gggggg       ')
  call yaml_scalar('|_____|_____|_____|_____|_____|     B    B               g             ')
  call yaml_scalar('|     |  :  |  :  |     |     |    B     B        i     g              ')
  call yaml_scalar('|     |-0+--|-0+--|     |     |    B    B         i     g        g     ')
  call yaml_scalar('|_____|__:__|__:__|_____|_____|___ BBBBB          i     g         g    ')
  call yaml_scalar('|  :  |     |     |  :  |     |    B    B         i     g         g    ')
  call yaml_scalar('|--+0-|     |     |-0+--|     |    B     B     iiii     g         g    ')
  call yaml_scalar('|__:__|_____|_____|__:__|_____|    B     B        i      g        g    ')
  call yaml_scalar('|     |  :  |  :  |     |     |    B BBBB        i        g      g     ')
  call yaml_scalar('|     |-0+--|-0+--|     |     |    B        iiiii          gggggg      ')
  call yaml_scalar('|_____|__:__|__:__|_____|_____|__BBBBB                                 ')
  call yaml_scalar('|     |     |     |  :  |     |                           TTTTTTTTT    ')
  call yaml_scalar('|     |     |     |--+0-|     |  DDDDDD          FFFFF        T        ')
  call yaml_scalar('|_____|_____|_____|__:__|_____| D      D        F        TTTT T        ')
  call yaml_scalar('|     |     |     |  :  |     |D        D      F        T     T        ')
  call yaml_scalar('|     |     |     |--+0-|     |D         D     FFFF     T     T        ')
  call yaml_scalar('|_____|_____|_____|__:__|_____|D___      D     F         T    T        ')
  call yaml_scalar('|     |     |  :  |     |     |D         D     F          TTTTT        ')
  call yaml_scalar('|     |     |--+0-|     |     | D        D     F         T    T        ')
  call yaml_scalar('|_____|_____|__:__|_____|_____|          D     F        T     T        ')
  call yaml_scalar('|     |     |     |     |     |         D               T    T         ')
  call yaml_scalar('|     |     |     |     |     |   DDDDDD       F         TTTT          ')
  call yaml_scalar('|_____|_____|_____|_____|_____|______                    www.bigdft.org   "')

  call yaml_mapping_close()

  call yaml_map('Reference Paper','The Journal of Chemical Physics 129, 014109 (2008)')
  call yaml_map('Version Number',package_version)
  call yaml_map('Timestamp of this run',yaml_date_and_time_toa())

  call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)
  if (ierr ==0) call yaml_map('Root process Hostname',trim(nodename_local))
  call yaml_map('Number of MPI tasks',bigdft_mpi%nproc)

  nthreads = 0
!$  nthreads=omp_get_max_threads()
  call yaml_map('OpenMP parallelization',nthreads>0)
  if (nthreads > 0) then
     call yaml_map('Maximal OpenMP threads per MPI task',nthreads)
  endif

END SUBROUTINE print_logo


!> Display the options of the configure tool (autotools)
subroutine print_configure_options()
  use yaml_output
  implicit none
  integer, parameter :: ln = 1024
  character(len = ln), dimension(4) :: buf

  call yaml_comment('Code compiling options',hfill='-')
  call yaml_mapping_open("Compilation options")
  call bigdft_config_get_user_args(buf(1), ln)
  call yaml_map("Configure arguments", '"'//trim(buf(1))//'"')
  call bigdft_config_get_compilers(buf(1), buf(2), buf(3), ln)
  call yaml_map("Compilers (CC, FC, CXX)", buf(1:3))
  call bigdft_config_get_compiler_flags(buf(1), buf(2), buf(3), buf(4), ln)
  call yaml_mapping_open("Compiler flags")
  if (len_trim(buf(1))>0) call yaml_map("CFLAGS",   trim(buf(1)))
  if (len_trim(buf(2))>0) call yaml_map("FCFLAGS",  trim(buf(2)))
  if (len_trim(buf(3))>0) call yaml_map("CXXFLAGS", trim(buf(3)))
  if (len_trim(buf(4))>0) call yaml_map("CPPFLAGS", trim(buf(4)))
  call yaml_mapping_close()
!!$  call bigdft_config_get_linker(buf(1), buf(2), buf(3), buf(4), ln)
!!$  call yaml_mapping_open("Linker")
!!$   call yaml_map("LD",      trim(buf(1)))
!!$   call yaml_map("LDFLAGS", trim(buf(2)))
!!$   call yaml_map("LIBS",    trim(buf(3)))
!!$   call yaml_map("Full linking options", trim(buf(4)))
!!$  call yaml_mapping_close()

 call yaml_mapping_close()

end subroutine print_configure_options




!> Write the eigenvalues-related information
subroutine write_eigenvalues_data(etol,orbs,mom_vec)
  use module_base
  use module_types
  use yaml_output
  implicit none
  real(gp), intent(in) :: etol
  type(orbitals_data), intent(in) :: orbs
  real(gp), dimension(:,:,:), pointer :: mom_vec
  !local variables
  logical :: degup,degdw
  integer :: ikptw,iorb,ikpt,jorb,isorb,nwrtmsg,ndegen
  real(gp) :: spinsignw,mx,my,mz,mpol,tolerance
  character(len=64) :: message
  character(len=150) :: commentline
  real(wp), dimension(2) :: preval

  commentline=repeat(' ',len(commentline))

  if (etol > 1.0_gp) then
     tolerance=0.0_gp
  else
     tolerance=etol
  end if

  ! Calculate and print the magnetisation, no matter the verbosity
  if (orbs%nspin == 2) then
     mpol = 0._gp
     do ikpt=1,orbs%nkpts
        isorb = (ikpt - 1) * orbs%norb
        do iorb = 1, orbs%norbu
           mpol = mpol + orbs%occup(isorb + iorb) * orbs%kwgts(ikpt)
        end do
        do iorb = orbs%norbu + 1, orbs%norb, 1
           mpol = mpol - orbs%occup(isorb + iorb) * orbs%kwgts(ikpt)
        end do
     end do
     call yaml_map("Total magnetization",mpol,fmt='(f9.6)')
  end if

  if (get_verbose_level() > 1) then
     call yaml_comment('Eigenvalues and New Occupation Numbers')

     call yaml_sequence_open('Orbitals',flow=.true.)
     call yaml_newline()

     do ikpt=1,orbs%nkpts
        if (orbs%nkpts > 1 .and. orbs%nspinor >= 2) then
           write(commentline,"(1x,A,I4.4,A,3F12.6)") &
                &   "Kpt #", ikpt, " BZ coord. = ", orbs%kpts(:, ikpt)
           !write(*,'(a)')trim(commentline)
           call yaml_comment(trim(commentline))
           call yaml_newline()
           ikptw=ikpt
        else
           ikptw=UNINITIALIZED(1)
        end if
        preval=0.0_wp
        nwrtmsg=0
        ndegen=0
        isorb = (ikpt - 1) * orbs%norb
        if (orbs%nspin==1.or.orbs%nspinor==4) then
           spinsignw=UNINITIALIZED(1.0_gp)
           do iorb=1,orbs%norb
              if (orbs%nspinor ==4 .and. associated(mom_vec)) then
                 mx=(mom_vec(2,iorb,1)/mom_vec(1,iorb,1))
                 my=(mom_vec(3,iorb,1)/mom_vec(1,iorb,1))
                 mz=(mom_vec(4,iorb,1)/mom_vec(1,iorb,1))
              else
                 mx=UNINITIALIZED(1.0_gp)
                 my=UNINITIALIZED(1.0_gp)
                 mz=UNINITIALIZED(1.0_gp)
              end if
              degup = find_degeneracy_up(iorb+isorb)
              call yaml_sequence()
              call write_orbital_data(orbs%eval(isorb + iorb),orbs%occup(isorb+iorb),&
                   spinsignw,ikptw,mx,my,mz)
              !yaml output (carriage return)
              if (iorb == orbs%norb .and. ikpt == orbs%nkpts) then
                 call yaml_sequence_close(advance='no')
                 !print *,'there',nwrtmsg,message
              end if
              call yaml_comment(trim(yaml_toa(iorb,fmt='(i5.5)')),advance='no')
              if (nwrtmsg==1) then
                 call yaml_comment(adjustl(message))
              else
                 call yaml_newline()
                 !call yaml_stream_attributes()
              end if
           end do
        else
           mx=UNINITIALIZED(1.0_gp)
           my=UNINITIALIZED(1.0_gp)
           mz=UNINITIALIZED(1.0_gp)

           do iorb=1,min(orbs%norbu,orbs%norbd)
              jorb=orbs%norbu+iorb
              call yaml_sequence()
              call write_orbital_data(orbs%eval(isorb + iorb),orbs%occup(isorb+iorb),&
                   1.0_gp,ikptw,mx,my,mz)
              call yaml_sequence()
              call write_orbital_data(orbs%eval(isorb + jorb),orbs%occup(isorb+jorb),&
                   -1.0_gp,ikptw,mx,my,mz)
              !yaml output (carriage return)
              degup=find_degeneracy_up(iorb+isorb)
              degdw=find_degeneracy_up(jorb+isorb)
              nwrtmsg=0
              if (degup .or. degdw) nwrtmsg=1
              if (degup .and. degdw) message='  <-deg->  '
              if (iorb == orbs%norbu .and. orbs%norbu==orbs%norbd .and. ikpt == orbs%nkpts) then
                 call yaml_sequence_close(advance='no')
              end if
              call yaml_comment(trim(yaml_toa(iorb,fmt='(i5.5)')),advance='no')
              if (nwrtmsg==1) then
                 call yaml_comment(adjustl(message))
              else

                 call yaml_newline()
              end if

           end do
           if (orbs%norbu > orbs%norbd) then
              do iorb=orbs%norbd+1,orbs%norbu
                 call yaml_sequence()
                 call write_orbital_data(orbs%eval(isorb+iorb),orbs%occup(isorb+iorb),&
                      1.0_gp,ikptw,mx,my,mz)
                 !yaml output (carriage return)
                 degup = find_degeneracy_up(iorb+isorb)
                 if (iorb == orbs%norbu .and. ikpt == orbs%nkpts) then
                    call yaml_sequence_close(advance='no')
                 end if
                 call yaml_comment(trim(yaml_toa(iorb,fmt='(i5.5)')),advance='no')
                 if (nwrtmsg==1) then
                    call yaml_comment(adjustl(message))
                 else
                    call yaml_newline()
                 end if
              end do
           else if (orbs%norbd > orbs%norbu) then
              do iorb=2*orbs%norbu+1,orbs%norbu+orbs%norbd
                 call yaml_sequence()
                 call write_orbital_data(orbs%eval(isorb+iorb),orbs%occup(isorb+iorb),&
                      -1.0_gp,ikptw,mx,my,mz)
                 !yaml output (carriage return)
                 degdw = find_degeneracy_down(iorb+isorb)
                 if (iorb == orbs%norbu+orbs%norbd .and. ikpt == orbs%nkpts) then
                    call yaml_sequence_close(advance='no')
                 end if
                 call yaml_comment(trim(yaml_toa(iorb,fmt='(i5.5)')),advance='no')
                 if (nwrtmsg==1) then
                    call yaml_comment(adjustl(message))
                 else
                    call yaml_newline()
                 end if
              end do
           end if
        end if
     end do
     ! Close the map of Eigenvalues and New Occupations Numbers
     !call yaml_mapping_close()
  end if
  !find fermi level
  if (orbs%efermi /= uninitialized(orbs%efermi)) then
     call yaml_map('Fermi Energy',orbs%efermi,fmt='(1pe21.14)')
  end if


contains

  function find_degeneracy_up(iorb) result(wrt)
    implicit none
    integer, intent(in) :: iorb
!    logical find_degeneracy_up
    !local variables
    logical :: wrt

    if (nwrtmsg==0 .and. iorb+1 < orbs%norbu) then
       wrt = abs(orbs%eval(iorb)-orbs%eval(iorb+1)) <= tolerance
       if (wrt) preval(1)=orbs%eval(iorb)
    else
       wrt = abs(orbs%eval(iorb)-preval(1)) <= tolerance
    end if
    !print *,'etol',etol,orbs%eval(iorb),preval(1),wrt,iorb
    nwrtmsg=0

    !disable degeneracy finder, problems to be fixed
    wrt=.false.

    if (wrt) then
       nwrtmsg=1
       message='  <-deg    '
    end if

    if (.not. wrt) preval(1)=orbs%eval(iorb)
  end function find_degeneracy_up

  function find_degeneracy_down(iorb) result(wrt)
    implicit none
    integer, intent(in) :: iorb
!    logical find_degeneracy_down
    !local variables
    logical :: wrt

    if (nwrtmsg==0 .and. iorb+1 < orbs%norb) then
       wrt = abs(orbs%eval(iorb)-orbs%eval(iorb+1)) <= tolerance
       if (wrt) preval(2)=orbs%eval(iorb)
    else
       wrt = abs(orbs%eval(iorb)-preval(2)) <= tolerance
    end if

    !disable degeneracy finder, problems to be fixed
    wrt=.false.


    nwrtmsg=0
    if (wrt) then
       nwrtmsg=1
       message='    deg->  '
    end if

    if (.not. wrt) preval(2)=orbs%eval(iorb)

  end function find_degeneracy_down

END SUBROUTINE write_eigenvalues_data


!> Writing rules, control if the last eigenvector is degenerate
!! do this for each spin
!! for each spin it is supposed that only the last group is not completely passed
!! and also that the components of each of the group but the last are the same for up and 
!! down polarisation. Do not work properly in the other cases
subroutine write_ig_eigenvectors(etol,orbse,nspin,norb,norbu,norbd)
   use module_base
   use module_types
   use yaml_output
   implicit none
   integer, intent(in) :: nspin,norb,norbu,norbd
   real(gp), intent(in) :: etol
   type(orbitals_data), intent(in) :: orbse
   !local variables
   character(len=64) :: message
   character(len=25) :: gapstring
   integer :: iorb,ndegen,nwrtmsg,ikpt,iorbst,ikptw
   real(gp) :: HLIGgap,mx,my,mz,spinsignw
   real(wp), dimension(2) :: preval
  character(len=150) :: commentline

  commentline=repeat(' ',len(commentline))


   !loop over all the k-points of the IG
   iorbst=0 !starting orbital in the k-points distribution
   !check if norbu and norbd are equal
   if (nspin==2 .and. orbse%norbu /= orbse%norbd) then
      write(*,*)'ERROR (write_ig_eigenvectors): the IG orbs structure should have norbu=norbd',orbse%norbu,orbse%norbd
      stop
   end if

  call yaml_sequence_open('Input Guess Orbitals',flow=.true.)!,advance='no')
  call yaml_newline()

  !always without spinors in the IG
  mx=UNINITIALIZED(1.0_gp)
  my=UNINITIALIZED(1.0_gp)
  mz=UNINITIALIZED(1.0_gp)


   do ikpt=1,orbse%nkpts
      if (orbse%nkpts > 1 .and. orbse%nspinor >= 2) then
         write(commentline,"(1x,A,I4.4,A,3F12.6)") &
              &   "Kpt #", ikpt, " BZ coord. = ", orbse%kpts(:, ikpt)
         !write(*,'(a)')trim(commentline)
         call yaml_comment(trim(commentline))
         ikptw=ikpt
      else
         ikptw=UNINITIALIZED(1)
      end if

      preval=0.0_wp
      nwrtmsg=0
      ndegen=0
      do iorb=1,orbse%norbu
         if (nspin==1) then
            spinsignw=UNINITIALIZED(1.0_gp)
            if (nwrtmsg==1) then
               if (abs(orbse%eval(iorb+iorbst)-preval(1)) <= etol) then
                  !degeneracy found
                  message='  <- found degeneracy'
                  ndegen=ndegen+1
               else
                  nwrtmsg=0
               end if
            end if
            if (abs(iorb - norb) <= 5) then
               nwrtmsg=1
               message=' <- '
            end if
            if (iorb == norb) then
               !calculate the IG HOMO-LUMO gap
               if(norb<orbse%norbu) then
                  HLIGgap=orbse%eval(iorb+1+iorbst)-orbse%eval(iorb+iorbst)
                  write(gapstring,'(a,f8.4,a)') ', H-L IG gap: ',HLIGgap*Ha_eV,' eV'
               else
                  gapstring=''
               end if
               nwrtmsg=1
               message=' <- Last InputGuess eval'//gapstring
               preval(1)=orbse%eval(iorb+iorbst)
            end if
            if (iorb-1 == norb) then
               nwrtmsg=1
               message=' <- First virtual eval '
            end if
            call yaml_sequence()
            call write_orbital_data(orbse%eval(iorbst+iorb),&
                 orbse%occup(iorbst+iorb),spinsignw,ikptw,mx,my,mz)

            if (nwrtmsg == 1) then
               !write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
               !   &   'evale(',iorb,')=',orbse%eval(iorb+iorbst),trim(message)
            else
               !if ((iorb <= 5 .or. iorb >= orbse%norbu-5) .or. get_verbose_level() > 0) &
               !write(*,'(1x,a,i0,a,1x,1pe21.14)') &
               !   &   'evale(',iorb,')=',orbse%eval(iorb+iorbst)
            end if
         else
            if (nwrtmsg==1) then
               if (abs(orbse%eval(iorb+iorbst)-preval(1)) <= etol .and. &
                  &   abs(orbse%eval(iorb+orbse%norbu+iorbst)-preval(2)) <= etol) then
               !degeneracy found
               message='  <-deg->  '
               !ndegen=ndegen+1 removed, only for non magnetized cases
            else if (abs(orbse%eval(iorb+iorbst)-preval(1)) <= etol) then
               !degeneracy found
               message='  <-deg    '
            else if (abs(orbse%eval(iorb+orbse%norbu+iorbst)-preval(2)) <= etol) then
               !degeneracy found
               message='    deg->  '
            else
               nwrtmsg=0
            end if
         end if
         if (iorb == norbu .and. iorb == norbd) then
            nwrtmsg=1
            message='  <-Last-> '
            preval(1)=orbse%eval(iorb+iorbst)
            preval(2)=orbse%eval(iorb+orbse%norbu+iorbst)
         else if (iorb == norbu) then
            nwrtmsg=1
            message='  <-Last   '
            preval(1)=orbse%eval(iorb+iorbst)
         else if (iorb == norbd) then
            nwrtmsg=1
            message='    Last-> '
            preval(2)=orbse%eval(iorb+orbse%norbu+iorbst)
         end if
         if ((iorb <= 5 .or. iorb >= orbse%norbu-5) .or. get_verbose_level() > 0) then
            call yaml_sequence()
            call write_orbital_data(orbse%eval(iorb+iorbst),&
                 orbse%occup(iorb+iorbst),&
                 1.0_gp,ikptw,mx,my,mz)
            call yaml_sequence()
            call write_orbital_data(orbse%eval(iorb+iorbst+orbse%norbu),&
                 orbse%occup(iorb+iorbst+orbse%norbu),&
                 -1.0_gp,ikptw,mx,my,mz)
         end if

         if (nwrtmsg==1) then
            !write(*,'(1x,a,i4,a,1x,1pe21.14,a12,a,i4,a,1x,1pe21.14)') &
            !   &   'evale(',iorb,',u)=',orbse%eval(iorb+iorbst),message,&
            !   &   'evale(',iorb,',d)=',orbse%eval(iorb+orbse%norbu+iorbst)
         else
            !if ((iorb <= 5 .or. iorb >= orbse%norbu-5) .or. get_verbose_level() > 0) &
            !write(*,'(1x,a,i4,a,1x,1pe21.14,12x,a,i4,a,1x,1pe21.14)') &
            !   &   'evale(',iorb,',u)=',orbse%eval(iorb+iorbst),&
            !   &   'evale(',iorb,',d)=',orbse%eval(iorb+orbse%norbu+iorbst)
         end if
      end if
      if (iorb == orbse%norbu .and. ikpt == orbse%nkpts) then
         call yaml_sequence_close(advance='no')
      end if
      if (nwrtmsg==1) then
         call yaml_comment(adjustl(message))
      else
         call yaml_newline()
      end if
   end do
   !increment k-points shift
   iorbst=iorbst+orbse%norb
end do

if (orbse%efermi /= uninitialized(orbse%efermi)) then
   call yaml_map('Fermi Energy',orbse%efermi,fmt='(1pe21.11)')
end if

!call yaml_stream_attributes()

END SUBROUTINE write_ig_eigenvectors


!> Write orbital information with NO advance
subroutine write_orbital_data(eval,occup,spinsign,ikpt,mx,my,mz)
  use module_base
  use yaml_output
  implicit none
  integer, intent(in) :: ikpt !< k-point id
  real(gp), intent(in) :: eval !< orbital energy
  real(gp), intent(in) :: occup !< orbital occupation number
  real(gp), intent(in) :: spinsign !< orbital spin (collinear and averaged)
  real(gp), intent(in) :: mx,my,mz !< spin magnetisation directions
  !local variables
  logical :: smallfmt


  call yaml_mapping_open(flow=.true.)
  !change the format if the spin and the k-point are initialized at the same way
  smallfmt=ikpt /= UNINITIALIZED(ikpt) .and. spinsign /= UNINITIALIZED(spinsign)

  if (smallfmt) then
     call yaml_map('e',eval,fmt='(1pe13.6)')
  else
     call yaml_map('e',eval,fmt='(1pe19.12)')
  end if

  if (occup /= UNINITIALIZED(occup)) then
     if (smallfmt) then
        call yaml_map('f',occup,fmt='(f5.3)')
     else
        call yaml_map('f',occup,fmt='(f6.4)')
     end if
  end if
  if (spinsign /= UNINITIALIZED(spinsign)) call yaml_map('s',int(spinsign),fmt='(i2)')
  if (ikpt /= UNINITIALIZED(ikpt)) then
     !if (int(spinsign)==-1) then
     !   call yaml_stream_attributes()
     !end if
     call yaml_map('k',ikpt,fmt='(i5)')
  end if
  if (mx /= UNINITIALIZED(mx) .and. my /= UNINITIALIZED(my) .and. mz /= UNINITIALIZED(mz)) &
     call yaml_map('M',(/mx,my,mz/),fmt='(f8.5)')
  call yaml_mapping_close(advance='no')

END SUBROUTINE write_orbital_data


!> Write DIIS weights
subroutine write_diis_weights(ncplx,idsx,ngroup,nkpts,itdiis,rds)
  use module_base
  use yaml_output
  implicit none
  integer, intent(in) :: ncplx,idsx,ngroup,nkpts,itdiis
  real(tp), dimension(ncplx,idsx+1,ngroup,nkpts), intent(in) :: rds
  !local variables
  integer :: j,igroup,ikpt
  character(len=2) :: mesupdw
  if (get_verbose_level() < 10) then
     !we restrict the printing to the first k point only.
     if (ngroup==1) then
        if (get_verbose_level() >0) then
!!$           write(*,'(1x,a,2x,18(1x,1pe9.2))')&
!!$             'DIIS wgts:',reshape(rds(1:ncplx,1:itdiis+1,1,1),&
!!$             (/ncplx*(itdiis+1)/))!,&
           !yaml output
           call yaml_newline()
           call yaml_sequence_open('DIIS weights',flow=.true.)
           do j=1,itdiis+1
              call yaml_sequence(yaml_toa(rds(1,j,1,1),fmt='(1pe9.2)'))
           end do
           call yaml_sequence_close()
           !call yaml_map('DIIS weights',&
           !     (/(rds(1:ncplx,j,1,1),j=1,itdiis+1)/),fmt='(1pe9.2)')
        end if
!        write(70,'(1x,a,1pe9.2)',advance='no')'DIIS wgts: [ ',rds(1:ncplx,1,1,1)
        do j=2,itdiis+1
!           write(70,'(a,1pe9.2)',advance='no')', ',rds(1:ncplx,j,1,1)
        end do
!        write(70,'(a)')']'
        !'(',ttr,tti,')'
     else if (get_verbose_level() >0) then
        do igroup=1,ngroup
           if (igroup==1) mesupdw='up'
           if (igroup==2) mesupdw='dw'
           write(*,'(1x,a,2x,18(1x,1pe9.2))')'DIIS wgts'//mesupdw//':',&
                (rds(1:ncplx,j,igroup,1),j=1,itdiis+1)
        end do
     end if
  else if (get_verbose_level() >0) then
     do ikpt = 1, nkpts
        if (ngroup==1) then
           write(*,'(1x,a,I3.3,a,2x,9(1x,(1pe9.2)))')'DIIS wgts (kpt #', ikpt, &
                & ')',(rds(1:ncplx,j,1,ikpt),j=1,itdiis+1)
        else
           do igroup=1,ngroup
              if (igroup==1) mesupdw='up'
              if (igroup==2) mesupdw='dw'
              write(*,'(1x,a,I3.3,a,2x,9(1x,a,2(1pe9.2),a))')'DIIS wgts (kpt #', ikpt, &
                   & ')'//mesupdw//':',('(',rds(1:ncplx,j,igroup,ikpt),')',j=1,itdiis+1)
           end do
        end if
     end do
  end if
END SUBROUTINE write_diis_weights


!> Print gnrms (residue per orbital)
subroutine write_gnrms(nkpts,norb,gnrms)
  use module_base
  use yaml_output
  implicit none
  integer, intent(in) :: norb,nkpts
  real(wp), dimension(norb,nkpts), intent(in) :: gnrms
  !local variables
  integer :: ikpt,iorb

  call yaml_newline()
  call yaml_sequence_open('Residues per orbital',flow=.true.)
  call yaml_newline()

  do ikpt=1,nkpts
     if (nkpts > 1) call yaml_comment('Kpt #'//adjustl(trim(yaml_toa(ikpt,fmt='(i4.4)'))))
     do iorb=1,norb
        call yaml_sequence(trim(yaml_toa(gnrms(iorb,ikpt),fmt='(1pe19.12)')),advance='no')
        if (ikpt == nkpts .and. iorb == norb)   call yaml_sequence_close(advance='no')
        call yaml_comment(trim(yaml_toa(iorb,fmt='(i5.5)')))
     end do
  end do

END SUBROUTINE write_gnrms


!> Print the atomic forces
subroutine write_forces(astruct,fxyz)
   use module_base
   use module_atoms
   use yaml_output
   implicit none
   !Arguments
   type(atomic_structure), intent(in) :: astruct
   real(gp), dimension(3,astruct%nat), intent(in) :: fxyz !< Atomic forces
   !Local variables
   real(gp) :: sumx,sumy,sumz
   integer :: iat

   sumx=0.d0
   sumy=0.d0
   sumz=0.d0
   call yaml_comment('Atomic Forces',hfill='-')
   call yaml_sequence_open('Atomic Forces (Ha/Bohr)')
   do iat=1,astruct%nat
      call yaml_sequence(advance='no')
      call yaml_mapping_open(flow=.true.)
      call yaml_map(trim(astruct%atomnames(astruct%iatype(iat))),fxyz(1:3,iat),fmt='(1pe20.12)')
      !call yaml_map('AU',fxyz(1:3,iat),fmt='(1pe20.12)')
      !call yaml_map('eV/A',fxyz(1:3,iat)*Ha_eV/Bohr_Ang,fmt='(1pe9.2)')
      call yaml_mapping_close(advance='no')
      call yaml_comment(trim(yaml_toa(iat,fmt='(i4.4)')))
!      write(*,'(1x,i5,1x,a6,3(1x,1pe12.5))') &
!      iat,trim(astruct%atomnames(astruct%iatype(iat))),(fxyz(j,iat),j=1,3)
      sumx=sumx+fxyz(1,iat)
      sumy=sumy+fxyz(2,iat)
      sumz=sumz+fxyz(3,iat)
   end do
   call yaml_sequence_close()
END SUBROUTINE write_forces


!> Write stress tensor matrix
subroutine write_strten_info(fullinfo,strten,volume,pressure,message)
  use module_base
  use yaml_output
  implicit none
  logical, intent(in) :: fullinfo
  real(gp), intent(in) :: volume,pressure
  character(len=*), intent(in) :: message
  real(gp), dimension(6), intent(in) :: strten
  !local variables

  call yaml_sequence_open(trim(message)//' stress tensor matrix (Ha/Bohr^3)')
  call yaml_sequence(yaml_toa((/strten(1),strten(6),strten(5)/),fmt='(1pg20.12)'))
  call yaml_sequence(yaml_toa((/strten(6),strten(2),strten(4)/),fmt='(1pg20.12)'))
  call yaml_sequence(yaml_toa((/strten(5),strten(4),strten(3)/),fmt='(1pg20.12)'))
  call yaml_sequence_close()
  !write(*,'(1x,a)')'Stress Tensor, '//trim(message)//' contribution (Ha/Bohr^3):'
  !write(*,'(1x,t10,10x,a,t30,10x,a,t50,10x,a)')'x','y','z'
  !write(*,'(1x,a,t10,1pe20.12,t30,1pe20.12,t50,1pe20.12)')'x',strten(1),strten(6),strten(5)
  !write(*,'(1x,a,t30,1pe20.12,t50,1pe20.12)')'y',strten(2),strten(4)
  !write(*,'(1x,a,t50,1pe20.12)')'z',strten(3)

  if (fullinfo) then
     call yaml_mapping_open('Pressure')
     call yaml_map('Ha/Bohr^3',pressure,fmt='(1pg22.14)')
     call yaml_map('GPa',pressure*AU_GPa,fmt='(1pg14.6)')
     call yaml_map('PV (Ha)',pressure*volume,fmt='(1pg22.14)')
     call yaml_mapping_close()
     !write(*,'(1x,a,1pe22.14,a,1pe14.6,a,1pe22.14)')'Pressure:',pressure,&
     !     ' (',pressure*AU_GPa,' GPa), P V:',pressure*volume
  end if

END SUBROUTINE write_strten_info

subroutine write_atomic_density_matrix(nspin,astruct,nl)
  use psp_projectors_base, only: DFT_PSP_projectors
  use yaml_output
  use module_atoms
  use ao_inguess, only: lmax_ao,ishell_toa
  use yaml_strings
  implicit none
  integer, intent(in) :: nspin
  type(DFT_PSP_projectors), intent(in) :: nl
  type(atomic_structure), intent(in) :: astruct
  !local variables
  logical :: flowrite
  integer :: ispin,l
  integer, dimension(0:lmax_ao) :: igamma
  character(len=32) :: msg
  type(atoms_iterator) :: atit

  call yaml_stream_attributes(flowrite=flowrite)

  if (.not. associated(nl%iagamma)) return
  call yaml_newline()
  call yaml_sequence_open('Atomic density matrix in the PSP projectors')
  !iterate above atoms
  atit=atoms_iter(astruct)
  do while(atoms_iter_next(atit))
     igamma=nl%iagamma(:,atit%iat)
     if (all(igamma == 0)) cycle
     call yaml_newline()
     call yaml_sequence(advance='no')
     if (flowrite) call yaml_mapping_open()
     call yaml_map('Symbol',trim(atit%name),advance='no')
     call yaml_comment('Atom '//trim(yaml_toa(atit%iat)))
     do l=0,lmax_ao
        !for the moment no imaginary part printed out
        if (igamma(l) == 0) cycle
        call yaml_mapping_open('Channel '//ishell_toa(l))
        if (nspin==4) then
           call yaml_newline()
           call yaml_map('Matrix',&
                nl%gamma_mmp(:,1:2*l+1,1:2*l+1,igamma(l),1),fmt='(1pg12.2)')
        else
           do ispin=1,nspin
              call yaml_newline()
              if (nspin==1) then
                 call f_strcpy(src='Matrix',dest=msg)
              else if (ispin==1) then
                 call f_strcpy(src='Spin up',dest=msg)
              else if (ispin==2) then
                 call f_strcpy(src='Spin down',dest=msg)
              end if
              call yaml_map(trim(msg),&
                   nl%gamma_mmp(1,1:2*l+1,1:2*l+1,igamma(l),ispin),fmt='(1pg12.2)')
           end do
        end if
        call yaml_mapping_close()
     end do
     if (flowrite) call yaml_mapping_close()
  end do
  call yaml_sequence_close()
end subroutine write_atomic_density_matrix


!> Assign some of the physical system variables
!! Performs also some cross-checks with other variables
!! The pointers in atoms structure have to be associated or nullified.
subroutine print_atomic_variables(atoms, hmax, ixc)
  use module_base
  use module_types
  use public_enums, only: RADII_SOURCE,PSPCODE_HGH,PSPCODE_HGH_K,PSPCODE_HGH_K_NLCC,&
       PSPCODE_PAW,PSPCODE_GTH, PSPCODE_PSPIO
  use public_keys, only: COEFF_KEY
  use module_xc
  use yaml_output
  implicit none
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(in) :: hmax
  integer, intent(in) :: ixc
  !Local variables
  logical :: nonloc
  integer :: i,j,j0,l,ityp,iat,natyp,mproj,inlcc
  real(gp) :: minrad
  real(gp), dimension(3,3) :: hij
  real(gp), dimension(2,2,4) :: offdiagarr
  character(len=500) :: name_xc1, name_xc2

  !If no atoms...
  if (atoms%astruct%ntypes == 0) return

  !print the pseudopotential matrices
  do l=1,4
     do i=1,2
        do j=i+1,3
           offdiagarr(i,j-i,l)=0._gp
           if (l==1) then
              if (i==1) then
                 if (j==2)   offdiagarr(i,j-i,l)=-0.5_gp*sqrt(3._gp/5._gp)
                 if (j==3)   offdiagarr(i,j-i,l)=0.5_gp*sqrt(5._gp/21._gp)
              else
                 offdiagarr(i,j-i,l)=-0.5_gp*sqrt(100._gp/63._gp)
              end if
           else if (l==2) then
              if (i==1) then
                 if (j==2)   offdiagarr(i,j-i,l)=-0.5_gp*sqrt(5._gp/7._gp)
                 if (j==3)   offdiagarr(i,j-i,l)=1._gp/6._gp*sqrt(35._gp/11._gp)
              else
                 offdiagarr(i,j-i,l)=-7._gp/3._gp*sqrt(1._gp/11._gp)
              end if
           else if (l==3) then
              if (i==1) then
                 if (j==2)   offdiagarr(i,j-i,l)=-0.5_gp*sqrt(7._gp/9._gp)
                 if (j==3)   offdiagarr(i,j-i,l)=0.5_gp*sqrt(63._gp/143._gp)
              else
                 offdiagarr(i,j-i,l)=-9._gp*sqrt(1._gp/143._gp)
              end if
           else if (l==4) then
              if (i==1) then
                 if (j==2)   offdiagarr(i,j-i,l)=-0.5_gp*sqrt(9._gp/11._gp)
                 if (j==3)   offdiagarr(i,j-i,l)=0.5_gp*sqrt(33._gp/65._gp)
              else
                 offdiagarr(i,j-i,l)=-11._gp*sqrt(1._gp/195._gp)
              end if
           end if
        end do
     end do
  end do

!  write(*,'(1x,a)')&
  !       '------------------------------------ Pseudopotential coefficients (Upper Triangular)'
  inlcc=0
  call yaml_comment('System Properties',hfill='-')
  call yaml_sequence_open('Properties of atoms in the system')
  do ityp=1,atoms%astruct%ntypes
     call yaml_sequence(advance='no')
     call yaml_map('Symbol',trim(atoms%astruct%atomnames(ityp)),advance='no')
     call yaml_comment('Type No. '//trim(yaml_toa(ityp,fmt='(i2.2)')))
     call yaml_map('No. of Electrons',atoms%nelpsp(ityp))
     natyp=0
     do iat=1,atoms%astruct%nat
        if (atoms%astruct%iatype(iat) == ityp) natyp=natyp+1
     end do
     call yaml_map('No. of Atoms',natyp)

     call yaml_mapping_open('Radii of active regions (AU)')!,flow=.true.)
       call yaml_map('Coarse',atoms%radii_cf(ityp,1),fmt='(f8.5)')
       call yaml_map('Fine',atoms%radii_cf(ityp,2),fmt='(f8.5)')
       call yaml_map('Coarse PSP',atoms%radii_cf(ityp,3),fmt='(f8.5)')
       call yaml_map('Source',RADII_SOURCE(atoms%iradii_source(ityp)))
       !if (atoms%radii_cf(ityp, 1) == UNINITIALIZED(1.0_gp)) then
       !   call yaml_map('Source','Hard-Coded')
       !else
       !   call yaml_map('Source','PSP File')
       !end if
     call yaml_mapping_close()

     minrad=1.e10_gp
     do i=0,4
        if (atoms%psppar(i,0,ityp)/=0._gp) then
           minrad=min(minrad,atoms%psppar(i,0,ityp))
        end if
     end do
     if (minrad < 1e9_gp) then
        if (atoms%radii_cf(ityp,2) /=0.0_gp) then
           call yaml_map('Grid Spacing threshold (AU)',2.5_gp*minrad,fmt='(f5.2)')
        else
           call yaml_map('Grid Spacing threshold (AU)',1.25_gp*minrad,fmt='(f5.2)')
        end if
     end if
     !control whether the grid spacing is too high
     if (hmax > 2.5_gp*minrad) then
        call yaml_warning('Chosen Grid spacings seem too high for the '// &
           & trim(atoms%astruct%atomnames(ityp))//' atom type. At your own risk!')
     end if

     select case(atoms%npspcode(ityp))
     case(PSPCODE_GTH)
        call yaml_map('Pseudopotential type','GTH')
     case(PSPCODE_HGH)
        call yaml_map('Pseudopotential type','HGH')
     case(PSPCODE_HGH_K)
        call yaml_map('Pseudopotential type','HGH-K')
     case(PSPCODE_HGH_K_NLCC)
        call yaml_map('Pseudopotential type','HGH-K + NLCC')
     case(PSPCODE_PAW)
        call yaml_map('Pseudopotential type','PAW + HGH')
     case(PSPCODE_PSPIO)
        call yaml_map('Pseudopotential type','PSPIO')
     end select
     if (atoms%psppar(0,0,ityp)/=0) then
        call yaml_mapping_open('Local Pseudo Potential (HGH convention)')
          call yaml_map('Rloc',atoms%psppar(0,0,ityp),fmt='(f9.5)')
          if (atoms%npspcode(ityp) == PSPCODE_GTH .or. &
               & atoms%npspcode(ityp) == PSPCODE_HGH .or. &
               & atoms%npspcode(ityp) == PSPCODE_HGH_K .or. &
               & atoms%npspcode(ityp) == PSPCODE_HGH_K_NLCC .or. &
               & atoms%npspcode(ityp) == PSPCODE_PAW) & ! to be removed later
               & call yaml_map(COEFF_KEY,atoms%psppar(0,1:4,ityp),fmt='(f9.5)')
        call yaml_mapping_close()
     end if
     !nlcc term
     if (atoms%npspcode(ityp) == PSPCODE_HGH_K_NLCC) then
        inlcc=inlcc+1
        call yaml_mapping_open('Non Linear Core Correction term')
            call yaml_map('Rcore',atoms%nlccpar(0,inlcc),fmt='(f9.5)')
            call yaml_map('Core charge',atoms%nlccpar(1,inlcc),fmt='(f9.5)')
        call yaml_mapping_close()
     end if
     !see if nonlocal terms are present
     nonloc=.false.
     verify_nl: do l=1,4
        if (any(atoms%psppar(l,0:3,ityp) /= 0._gp)) then
           nonloc=.true.
           exit verify_nl
        end if
     end do verify_nl
     if (nonloc) then
        call yaml_sequence_open('NonLocal PSP Parameters')
        do l=1,4
           if (any(atoms%psppar(l,0:3,ityp) /= 0._gp)) then
              call yaml_sequence(advance='no')
              call yaml_map('Channel (l)',l-1)
              if (atoms%psppar(l,0,ityp) > 0._gp) &
                   & call yaml_map('Rloc',atoms%psppar(l,0,ityp),fmt='(f9.5)')
              hij=0._gp
              do i=1,3
                 hij(i,i)=atoms%psppar(l,i,ityp)
              end do
              if (atoms%npspcode(ityp) == PSPCODE_HGH) then !traditional HGH convention
                 hij(1,2)=offdiagarr(1,1,l)*atoms%psppar(l,2,ityp)
                 hij(1,3)=offdiagarr(1,2,l)*atoms%psppar(l,3,ityp)
                 hij(2,3)=offdiagarr(2,1,l)*atoms%psppar(l,3,ityp)
              else if (atoms%npspcode(ityp) == PSPCODE_HGH_K &
                   .or. atoms%npspcode(ityp) == PSPCODE_HGH_K_NLCC &
                   .or. atoms%npspcode(ityp) == PSPCODE_PSPIO) then !HGH-K convention
                 hij(1,2)=atoms%psppar(l,4,ityp)
                 hij(1,3)=atoms%psppar(l,5,ityp)
                 hij(2,3)=atoms%psppar(l,6,ityp)
              end if
              call yaml_sequence_open('h_ij matrix')
                call yaml_sequence(trim(yaml_toa(hij(1,1:3),fmt='(f9.5)')))
                call yaml_sequence(trim(yaml_toa((/hij(1,2),hij(2,2),hij(2,3)/),fmt='(f9.5)')))
                call yaml_sequence(trim(yaml_toa((/hij(1,3),hij(2,3),hij(3,3)/),fmt='(f9.5)')))
              call yaml_sequence_close()
           end if
        end do
        call yaml_sequence_close()
     end if
     ! PAW case.
     if (atoms%npspcode(ityp) == PSPCODE_PAW) then
        if (atoms%pawtab(ityp)%has_wvl > 0) then
           call yaml_sequence_open('NonLocal PSP Parameters (PAW)')
           j0 = 1
           l = 1
           do i = 1, atoms%pawtab(ityp)%basis_size
              call yaml_sequence(advance='no')
              call yaml_map('Channel (l,m,n)', atoms%pawtab(ityp)%indlmn(1:3, l))
              l = l + atoms%pawtab(ityp)%orbitals(i) * 2 + 1
              call yaml_sequence_open('gaussians')
              do j = j0, j0 + atoms%pawtab(ityp)%wvl%pngau(i) / 2 - 1
                 call yaml_sequence(advance='no')
                 call yaml_mapping_open(flow = .true.)
                 call yaml_map('factor', atoms%pawtab(ityp)%wvl%pfac(:,j),fmt='(f10.6)')
                 call yaml_map('exponent', atoms%pawtab(ityp)%wvl%parg(:,j),fmt='(f10.6)')
                 call yaml_mapping_close()
              end do
              j0 = j0 + atoms%pawtab(ityp)%wvl%pngau(i)
              call yaml_sequence_close()
           end do
           call yaml_sequence_close()
        end if
        mproj = 0
        do i = 1, atoms%pawtab(ityp)%basis_size
           mproj = mproj + 2 * atoms%pawtab(ityp)%orbitals(i) + 1
        end do
     else
        mproj = 0
        do l=1,4 
           do i=1,3 
              if (atoms%psppar(l,i,ityp) /= 0.0_gp) mproj=mproj+2*l-1
           enddo
        enddo
     end if
     !call numb_proj(ityp,atoms%astruct%ntypes,atoms%psppar,atoms%npspcode,mproj)
     call yaml_map('No. of projectors',mproj)

     !control if the PSP is calculated with the same XC value
     if (atoms%ixcpsp(ityp) < 0) then
        call xc_get_name(name_xc1, atoms%ixcpsp(ityp), XC_MIXED)
     else
        call xc_get_name(name_xc1, atoms%ixcpsp(ityp), XC_ABINIT)
     end if
     if (ixc < 0) then
        call xc_get_name(name_xc2, ixc, XC_MIXED)
     else
        call xc_get_name(name_xc2, ixc, XC_ABINIT)
     end if
     call yaml_map('PSP XC','"'//trim(name_xc1)//'"')
     if (trim(name_xc1) /= trim(name_xc2)) then
        call yaml_warning('PSP generated with a different XC. Input XC is "'//trim(name_xc2) // '"')
     end if
  end do
  call yaml_sequence_close()
!!!  tt=dble(norb)/dble(nproc)
!!!  norbp=int((1.d0-eps_mach*tt) + tt)
!!!  !if (verb.eq.0) write(*,'(1x,a,1x,i0)') 'norbp=',norbp

  ! if linear scaling applied with more then InputGuess, then go read input.lin for radii
  !  if (in%linear /= 'OFF' .and. in%linear /= 'LIG') then
  !     lin%nlr=atoms%astruct%nat
  !     call allocateBasicArrays(atoms, lin)
  !     call readLinearParameters(verb, nproc, lin, atoms, atomNames)
  !  end if
END SUBROUTINE print_atomic_variables


!> Display an estimation of the occupied memory
subroutine print_memory_estimation(mem)
  use module_types
  use yaml_output
  use yaml_strings
  implicit none
  type(memory_estimation), intent(in) :: mem

  call yaml_comment('Estimation of Memory Consumption',hfill='-')
  call yaml_mapping_open('Memory requirements for principal quantities (MiB.KiB)')
    call yaml_map('Subspace Matrix',trim(MibdotKib(mem%submat)),advance='no')
      call yaml_comment('(Number of Orbitals:'//trim(yaml_toa(mem%norb))//')',tabbing=50)
    call yaml_map('Single orbital',trim(MibdotKib(mem%oneorb)),advance='no')
      call yaml_comment('(Number of Components:'//trim(yaml_toa(mem%ncomponents))//')',tabbing=50)
    call yaml_map('All (distributed) orbitals',trim(MibdotKib(mem%allpsi_mpi)),advance='no')
      call yaml_comment('(Number of Orbitals per MPI task:'//trim(yaml_toa(mem%norbp))//')',tabbing=50)
    call yaml_map('Wavefunction storage size',trim(MibdotKib(mem%psistorage)),advance='no')
      call yaml_comment('(DIIS/SD workspaces included)',tabbing=50)
    call yaml_map('Nonlocal Pseudopotential Arrays',trim(MibdotKib(mem%projarr)))
    call yaml_map('Full Uncompressed (ISF) grid',trim(MibdotKib(mem%grid)))
    call yaml_map('Workspaces storage size',trim(MibdotKib(mem%workarr)))
  call yaml_mapping_close()

  call yaml_mapping_open('Accumulated memory requirements during principal run stages (MiB.KiB)')
     call yaml_map('Kernel calculation',trim(MibdotKib(mem%kernel)))
     call yaml_map('Density Construction',trim(MibdotKib(mem%density)))
     call yaml_map('Poisson Solver',trim(MibdotKib(mem%psolver)))
     call yaml_map('Hamiltonian application',trim(MibdotKib(mem%ham)))
     call yaml_map('Orbitals Orthonormalization',trim(MibdotKib(mem%ham+mem%submat)))
!           call yaml_comment('Wfn, Work, Den, Ker ',tabbing=50)
  call yaml_mapping_close()
  call yaml_map('Estimated Memory Peak (MB)',yaml_toa(mega(mem%peak)))

contains

  function mega(omemory)
    implicit none
    real(kind=8), intent(in) :: omemory
    integer(kind=8) :: mega
    mega=int(omemory/1048576.d0,kind=8)
  end function mega

  function kappa(omemory)
    implicit none
    real(kind=8), intent(in) :: omemory
    integer :: kappa
    kappa=ceiling((omemory-aint(omemory/1048576.d0)*1048576.d0)/1024.d0)
  end function kappa

  function MiBdotKiB(omemory)
    implicit none
    real(kind=8), intent(in) :: omemory
    character(len=50) MiBdotKiB

    MiBdotKiB=repeat(' ',len(MiBdotKiB))

    MiBdotKiB=trim(adjustl(yaml_toa(int(mega(omemory)))))//'.'//&
         trim(adjustl(yaml_toa(int(kappa(omemory)))))

  end function MiBdotKiB

END SUBROUTINE print_memory_estimation


!> Display information about the box and the grid
subroutine print_atoms_and_grid(Glr, atoms, rxyz, hx, hy, hz)
  use module_defs
  use numerics, only: Bohr_Ang
  use module_types
  use yaml_output
  use yaml_strings
  use locregs
  implicit none
  !Arguments
  type(atoms_data), intent(in) :: atoms
  type(locreg_descriptors), intent(in) :: Glr
  real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
  real(gp), intent(in) :: hx, hy, hz
  !Local variables
  integer :: iunit !, iat

  if (atoms%astruct%ntypes > 0) then
     call yaml_comment('Atom Positions (specified and grid units)',hfill='-')
     ! New version
     call yaml_mapping_open('Atomic structure')
     call yaml_get_default_stream(unit = iunit)
     call wtyaml(iunit, UNINITIALIZED(1.d0), rxyz, atoms%astruct, .false., rxyz, &
          .true., atoms%astruct%shift, (/ hx, hy, hz /))
     call yaml_mapping_close()
  end if
  call yaml_comment('Grid properties',hfill='-')
  call yaml_map('Box Grid spacings',(/hx,hy,hz/),fmt='(f7.4)')
  call yaml_mapping_open('Sizes of the simulation domain')
  call yaml_map('AU',(/atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),atoms%astruct%cell_dim(3)/),fmt='(1pg12.5)')
  call yaml_map('Angstroem',(/atoms%astruct%cell_dim(1)*Bohr_Ang,&
       atoms%astruct%cell_dim(2)*Bohr_Ang,atoms%astruct%cell_dim(3)*Bohr_Ang/),fmt='(1pg12.5)')
  call yaml_map('Grid Spacing Units',(/Glr%d%n1,Glr%d%n2,Glr%d%n3/),fmt='(i4)')
  call yaml_mapping_open('High resolution region boundaries (GU)',flow=.false.)
  call yaml_map('From',(/Glr%d%nfl1,Glr%d%nfl2,Glr%d%nfl3/),fmt='(i4)')
  call yaml_map('To',(/Glr%d%nfu1,Glr%d%nfu2,Glr%d%nfu3/),fmt='(i4)')
  call yaml_mapping_close()
  call yaml_mapping_close()
  call yaml_map('High Res. box is treated separately',Glr%hybrid_on)
END SUBROUTINE print_atoms_and_grid


!> Write atomic file in yaml format
subroutine wtyaml(iunit,energy,rxyz,astruct,wrtforces,forces, &
     & wrtlog, shift, hgrids)
  use module_base, only: f_err_throw,Bohr_Ang
  use module_defs, only: gp, UNINITIALIZED
  use yaml_output
  use yaml_strings
  use module_atoms, only: atomic_structure,frozen_itof
  use ao_inguess, only: charge_and_spol
  implicit none
  !Arguments
  logical, intent(in) :: wrtforces !< True if write the atomic forces
  logical, intent(in) :: wrtlog    !< Level of verbosity (false remove some parts)
  integer, intent(in) :: iunit     !< File unit
  type(atomic_structure), intent(in) :: astruct
  real(gp), intent(in) :: energy
  real(gp), dimension(3,astruct%nat), intent(in) :: rxyz,forces
  real(gp), dimension(3), intent(in) :: shift, hgrids
  !local variables
  logical :: reduced, perx, pery, perz
  character(len=4) :: frzchain
  real(gp), dimension(3) :: xred
  real(gp) :: factor
  integer :: iat,ichg,ispol

  reduced=.false.
  Units: select case(trim(astruct%units))
  case('angstroem','angstroemd0')
     call yaml_map('Units','angstroem', unit = iunit)
     factor=Bohr_Ang
  case('reduced')
     if (.not. wrtlog) then
        call yaml_map('Units','reduced', unit = iunit)
        reduced=.true.
     end if
     factor = 1.0_gp
  case('atomic','atomicd0','bohr','bohrd0')
     ! Default
     factor=1.0_gp
     !Important to display the default units (TD)
     if (wrtlog) call yaml_map('Units','bohr')
  case default
     call f_err_throw('Writing the atomic file. Error, unknown units ("'// trim(astruct%units)//'")', &
          & err_name='BIGDFT_RUNTIME_ERROR')
  end select Units

  !cell information
  perx = .false.
  pery = .false.
  perz = .false.
  BC :select case(astruct%geocode)
  case('S')
     call yaml_sequence_open('Cell', flow=.true., unit = iunit)
     call yaml_sequence(yaml_toa(astruct%cell_dim(1)*factor), unit = iunit) !x
     call yaml_sequence('.inf', unit = iunit)             !y
     call yaml_sequence(yaml_toa(astruct%cell_dim(3)*factor), unit = iunit) !z
     call yaml_sequence_close(unit = iunit)
     !angdeg to be added
     perx = .true.
     pery = .false.
     perz = .true.
  case('W')
     call yaml_sequence_open('Cell', flow=.true., unit = iunit)
     call yaml_sequence('.inf', unit = iunit)             !x
     call yaml_sequence('.inf', unit = iunit)             !y
     call yaml_sequence(yaml_toa(astruct%cell_dim(3)*factor), unit = iunit) !z
     call yaml_sequence_close(unit = iunit)
     perx = .false.
     pery = .false.
     perz = .true.
  case('P')
     call yaml_map('Cell',(/astruct%cell_dim(1)*factor, &
          & astruct%cell_dim(2)*factor, astruct%cell_dim(3)*factor/), unit = iunit)
     !angdeg to be added
     perx = .true.
     pery = .true.
     perz = .true.
  case('F')
     ! Default
     !call yaml_map('BC','free')
  end select BC

  !Write atomic positions
  call yaml_sequence_open('Positions', unit = iunit)
  do iat=1,astruct%nat
     !for very large systems, consider deactivating the printing, but should do this in a cleaner manner
     !if (astruct%nat < 500.or.(.not. wrtlog)) then
     call yaml_sequence(advance='no', unit = iunit)
     !end if
     if (extra_info(iat)) then
        call yaml_mapping_open(flow=.true., unit = iunit)
     end if
     xred=rxyz(:,iat)
     if (reduced) then
        if (perx) xred(1)=rxyz(1,iat)/astruct%cell_dim(1)
        if (pery) xred(2)=rxyz(2,iat)/astruct%cell_dim(2)
        if (perz) xred(3)=rxyz(3,iat)/astruct%cell_dim(3)
     else
        !Multiply by the factor to have the right units
        xred = xred*factor
     end if
     if (wrtlog) then
        !for very large systems, consider deactivating the printing, but should do this in a cleaner manner
        !if (astruct%nat < 500) then
        call print_one_atom(trim(astruct%atomnames(astruct%iatype(iat))),&
             xred,hgrids,iat)
        !end if
!!$        call yaml_map(trim(astruct%atomnames(astruct%iatype(iat))),&
!!$             & xred,fmt="(g18.10)", unit = iunit, advance = "no")
!!$        xred(1:3) = rxyz(1:3,iat) / hgrids
!!$        write(gu, "('[ 'F6.2', 'F6.2', 'F6.2'] 'I4.4)") xred, iat
!!$        call yaml_comment(gu, unit = iunit)
     else
        call yaml_map(trim(astruct%atomnames(astruct%iatype(iat))),&
             & xred,fmt="(1pg25.17)", unit = iunit)
     end if
     if (extra_info(iat)) then
        call charge_and_spol(astruct%input_polarization(iat),ichg,ispol)
        if (ispol /=0) call yaml_map('IGSpin',ispol, unit = iunit)
        if (ichg /=0) call yaml_map('IGChg',ichg, unit = iunit)
        if (astruct%ifrztyp(iat) /= 0) then
           call frozen_itof(astruct%ifrztyp(iat),frzchain)
           call yaml_map('Frozen',frzchain, unit = iunit)
        end if
        call yaml_mapping_close(unit = iunit)
     end if
  end do
  call yaml_sequence_close(unit = iunit) !positions

  !Write atomic forces
  if (wrtforces) then
     call yaml_sequence_open('Forces (Ha/Bohr)', unit = iunit)
     do iat=1,astruct%nat
        call yaml_sequence(advance='no', unit = iunit)
        call yaml_map(trim(astruct%atomnames(astruct%iatype(iat))),forces(:,iat),fmt='(1pg25.17)', unit = iunit)
     end do
     call yaml_sequence_close(unit = iunit) !values
  end if
  if (wrtlog) then
     call yaml_map('Rigid Shift Applied (AU)',(/-shift(1),-shift(2),-shift(3)/),fmt='(1pg12.5)')
  else
     call yaml_mapping_open('Properties', unit = iunit)
     call yaml_map('Timestamp',yaml_date_and_time_toa(), unit = iunit)
     if (energy /= 0.0_gp .and. energy /= UNINITIALIZED(energy)) then
        call yaml_map("Energy (Ha)", energy, unit = iunit)
     end if
     call yaml_mapping_close(unit = iunit) !properties
  end if

contains

  function extra_info(iat)
    implicit none
    integer, intent(in) :: iat
    logical extra_info
    extra_info=astruct%input_polarization(iat) /=100 .or. astruct%ifrztyp(iat)/=0
  end function extra_info

  subroutine print_one_atom(atomname,rxyz,hgrids,id)
    implicit none
    integer, intent(in) :: id
    character(len=*), intent(in) :: atomname
    double precision, dimension(3), intent(in) :: rxyz,hgrids
    !local variables
    character(len=*), parameter :: fmtat='(1pg18.10)',fmtg='(F7.2)',fmti='(i4.4)'
    integer :: i

    call yaml_sequence_open(atomname,flow=.true.)
    do i=1,3
       call yaml_sequence(yaml_toa(rxyz(i),fmt=fmtat))
    end do
    call yaml_sequence_close(advance='no')
    call yaml_comment(trim(yaml_toa(rxyz/hgrids/factor,fmt=fmtg))//trim(yaml_toa(id,fmt=fmti))) !we can also put tabbing=

  end subroutine print_one_atom

END SUBROUTINE wtyaml


!> Display the waefunctions descriptors (segments and points)
subroutine print_wfd(wfd)
  use yaml_output
  use compression
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd

  call yaml_mapping_open('Wavefunctions Descriptors, full simulation domain')
  !write(*,'(1x,a)')&
  !   &   '------------------------------------------------- Wavefunctions Descriptors Creation'

  !write(*,'(2(1x,a,i10))') &
  !     &   'Coarse resolution grid: Number of segments= ',Glr%wfd%nseg_c,'points=',Glr%wfd%nvctr_c
  call yaml_mapping_open('Coarse resolution grid')!,flow=.true.)
  call yaml_map('No. of segments',wfd%nseg_c)
  call yaml_map('No. of points',wfd%nvctr_c)
  call yaml_mapping_close()

  !write(*,'(2(1x,a,i10))')
  !'  Fine resolution grid: Number of segments= ',Glr%wfd%nseg_f,'points=',Glr%wfd%nvctr_f
  call yaml_mapping_open('Fine resolution grid')!,flow=.true.)
  call yaml_map('No. of segments',wfd%nseg_f)
  call yaml_map('No. of points',wfd%nvctr_f)
  call yaml_mapping_close()
  call yaml_mapping_close()

END SUBROUTINE print_wfd

subroutine print_nlpsp(nlpsp)
  use module_defs, only: gp
  use module_types, only: DFT_PSP_projectors
  use yaml_output
  implicit none
  type(DFT_PSP_projectors), intent(in) :: nlpsp
  !local variables
  integer :: iat,ilr,sizemask,maxmask,totmask,totpack

  call yaml_mapping_open('NonLocal PSP Projectors Descriptors')
  if (nlpsp%on_the_fly) then
     call yaml_map('Creation strategy','On-the-fly')
  else
     call yaml_map('Creation strategy','Once-and-for-all')
  end if
  call yaml_map('Total number of projectors',nlpsp%nproj)
  call yaml_map('Total number of components',nlpsp%nprojel)
  call yaml_map('Percent of zero components',nint(100.0_gp*nlpsp%zerovol))
  !calculate the amount of memory spent in the descriptor for the wavefunction
  maxmask=0
  totmask=0
  totpack=0
  do iat=1,nlpsp%natoms
     if (nlpsp%projs(iat)%mproj>0) then
        totpack=max(totpack,nlpsp%projs(iat)%region%plr%wfd%nvctr_c+&
             7*nlpsp%projs(iat)%region%plr%wfd%nvctr_f)
     end if
     sizemask=0
     if (associated(nlpsp%projs(iat)%region%tolr)) then
        !do ilr=1,nlpsp%pspd(iat)%nlr
        do ilr=1,size(nlpsp%projs(iat)%region%tolr)
           sizemask=sizemask+&
                nlpsp%projs(iat)%region%tolr(ilr)%nmseg_c+&
                nlpsp%projs(iat)%region%tolr(ilr)%nmseg_f
        end do
     end if
     maxmask=max(maxmask,sizemask)
     totmask=totmask+sizemask
  end do
  totpack=totpack*4
  if (associated(nlpsp%scpr)) totpack=totpack+size(nlpsp%scpr)
  if (associated(nlpsp%cproj)) totpack=totpack+size(nlpsp%cproj)*2
  if (totpack /=0) &
       call yaml_map('Size of workspaces',totpack)
  if (maxmask /=0) &
       call yaml_map('Maximum size of masking arrays for a projector',3*maxmask)
  if (totmask /=0) &
       call yaml_map('Cumulative size of masking arrays',3*totmask)

  call yaml_mapping_close()
END SUBROUTINE print_nlpsp


!> Display information about the electronic orbitals
subroutine print_orbitals(orbs, geocode)
  use module_types, only: orbitals_data
  use module_defs, only: gp
  use yaml_output
  use yaml_strings
  implicit none
  type(orbitals_data), intent(in) :: orbs
  character(len = 1), intent(in) :: geocode

  integer :: jproc, nproc, jpst, norbme, norbyou, nelec
  integer :: ikpts, iorb1, iorb
  real(gp) :: rocc

  nelec = int(sum(orbs%occup) + 1d-12) / orbs%nkpts

  call yaml_comment('Electronic Orbital Initialization',hfill='-')
  !call yaml_comment('Occupation Numbers',hfill='-')
  call yaml_map('Total Number of Electrons',nelec,fmt='(i8)')

  ! Number of orbitals
  if (orbs%nspin==1) then
     call yaml_map('Spin treatment','Averaged')
     if (mod(nelec,2).ne.0) then
        call yaml_warning('Odd number of electrons, no closed shell system')
        !write(*,'(1x,a)') 'WARNING: odd number of electrons, no closed shell system'
     end if
  else if(orbs%nspin==4) then
     call yaml_map('Spin treatment','Spinorial (non-collinearity possible)')
  else
     call yaml_map('Spin treatment','Collinear')
  end if

  !distribution of wavefunction arrays between processors
  !tuned for the moment only on the cubic distribution
  call yaml_mapping_open('Orbitals Repartition')
  jpst=0
  nproc = size(orbs%norb_par, 1)
  do jproc=0,nproc-1
     norbme=orbs%norb_par(jproc,0)
     norbyou=orbs%norb_par(min(jproc+1,nproc-1),0)
     if (norbme /= norbyou .or. jproc == nproc-1) then
        call yaml_map('MPI tasks '//trim(yaml_toa(jpst,fmt='(i0)'))//'-'//trim(yaml_toa(jproc,fmt='(i0)')),norbme,fmt='(i0)')
        !write(*,'(3(a,i0),a)')&
        !     ' Processes from ',jpst,' to ',jproc,' treat ',norbme,' orbitals '
        jpst=jproc+1
     end if
  end do
  !write(*,'(3(a,i0),a)')&
  !     ' Processes from ',jpst,' to ',nproc-1,' treat ',norbyou,' orbitals '
  call yaml_mapping_close()

  call yaml_map('Total Number of Orbitals',orbs%norb,fmt='(i8)')

  !No orbs finished
  if (orbs%norb == 0) return

  call yaml_sequence_open('Input Occupation Numbers')
  do ikpts=1,orbs%nkpts
     if (geocode /= 'F') then
        call yaml_comment('Kpt #' // adjustl(trim(yaml_toa(ikpts,fmt='(i4.4)'))) // ' BZ coord. = ' // &
        & trim(yaml_toa(orbs%kpts(:, ikpts),fmt='(f12.6)')))
     end if
     call yaml_sequence(advance='no')
     call yaml_mapping_open('Occupation Numbers',flow=.true.)
     !write(*,'(1x,a,t28,i8)') 'Total Number of Orbitals',norb
     iorb1=1
     rocc=orbs%occup(1+(ikpts-1)*orbs%norb)
     do iorb=1,orbs%norb
        if (orbs%occup(iorb+(ikpts-1)*orbs%norb) /= rocc) then
           if (iorb1 == iorb-1) then
              call yaml_map('Orbital No.'//trim(yaml_toa(iorb1)),rocc,fmt='(f6.4)')
              !write(*,'(1x,a,i0,a,f6.4)') 'occup(',iorb1,')= ',rocc
           else
           call yaml_map('Orbitals No.'//trim(yaml_toa(iorb1))//'-'//&
                adjustl(trim(yaml_toa(iorb-1))),rocc,fmt='(f6.4)')
           !write(*,'(1x,a,i0,a,i0,a,f6.4)') 'occup(',iorb1,':',iorb-1,')= ',rocc
           end if
           rocc=orbs%occup(iorb+(ikpts-1)*orbs%norb)
           iorb1=iorb
        end if
     enddo
     if (iorb1 == orbs%norb) then
        call yaml_map('Orbital No.'//trim(yaml_toa(orbs%norb)),orbs%occup(ikpts*orbs%norb),fmt='(f6.4)')
        !write(*,'(1x,a,i0,a,f6.4)') 'occup(',norb,')= ',occup(norb)
     else
        call yaml_map('Orbitals No.'//trim(yaml_toa(iorb1))//'-'//&
             adjustl(trim(yaml_toa(orbs%norb))),orbs%occup(ikpts*orbs%norb),fmt='(f6.4)')
        !write(*,'(1x,a,i0,a,i0,a,f6.4)') 'occup(',iorb1,':',norb,')= ',occup(norb)
     end if
     call yaml_mapping_close()
  end do
  call yaml_sequence_close()
END SUBROUTINE print_orbitals
