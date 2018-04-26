!> @file
!! Routines associated to the convolutions for the scalar products
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> Wrapper
subroutine wnrm_wrap(ncplx,mvctr_c,mvctr_f,psi,scpr)
  use module_defs, only: wp,dp
  implicit none
  integer, intent(in) :: mvctr_c,mvctr_f,ncplx
  real(wp), dimension((mvctr_c+7*mvctr_f)*ncplx), intent(in) :: psi
  real(dp), intent(out) :: scpr
  !local variables
  integer :: i_f
  real(dp) :: scalp

  i_f=min(mvctr_f,1)
 
  call wnrm(mvctr_c,mvctr_f,psi,psi(mvctr_c+i_f),scpr)

  if (ncplx ==2) then
     call wnrm(mvctr_c,mvctr_f,&
          psi(mvctr_c+7*mvctr_f+1),psi(mvctr_c+7*mvctr_f+mvctr_c+i_f),scalp)
     scpr=scpr+scalp
  end if
  
END SUBROUTINE wnrm_wrap


!> Calculates the norm SQUARED (scpr) of a wavefunction (in vector form)
!! given the distribution of the data also dnrm2 or ddot can be called
subroutine wnrm(mvctr_c,mvctr_f,psi_c,psi_f,scpr)
  use module_defs, only: wp,dp
  implicit none
  !Arguments
  integer, intent(in) :: mvctr_c,mvctr_f
  real(wp), dimension(mvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(in) :: psi_f
  real(dp), intent(out) :: scpr
  !local variables
  integer :: i
  real(dp) :: scpr0,scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7
!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!    real(gp) :: tel

!!!  !dee
!!!    open(unit=97,file='time_wnrm',status='unknown',position='append')
!!!    call system_clock(ncount0,ncount_rate,ncount_max)

    scpr=0.0_dp
    scpr0=0.0_dp
    scpr1=0.0_dp
    scpr2=0.0_dp
    scpr3=0.0_dp
    scpr4=0.0_dp
    scpr5=0.0_dp
    scpr6=0.0_dp
    scpr7=0.0_dp

!$omp parallel default(private)&
!$omp shared(mvctr_c,mvctr_f,psi_c,psi_f,scpr) firstprivate(scpr0,scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7)

!$omp do schedule(static)
 do i=1,mvctr_c
    scpr0=scpr0+real(psi_c(i),dp)**2
 enddo
!$omp end do nowait
!$omp do schedule(guided)
 do i=1,mvctr_f
    scpr1=scpr1+real(psi_f(1,i),dp)**2
    scpr2=scpr2+real(psi_f(2,i),dp)**2
    scpr3=scpr3+real(psi_f(3,i),dp)**2
    scpr4=scpr4+real(psi_f(4,i),dp)**2
    scpr5=scpr5+real(psi_f(5,i),dp)**2
    scpr6=scpr6+real(psi_f(6,i),dp)**2
    scpr7=scpr7+real(psi_f(7,i),dp)**2
enddo
!$omp end do
    scpr0=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

!$omp critical 
    scpr=scpr+scpr0
!$omp end critical

!$omp end parallel

!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!    write(97,*) 'wnrm:',tel
!!!    close(97)

END SUBROUTINE wnrm


!> Calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
!! The array mask is then constructed so that successive application of the projector on the same object can
!! be done without bitonic search
subroutine wpdot_keys(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,&
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,&
     scpr)
  use module_defs, only: wp
  implicit none

  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), intent(out) :: scpr
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,iboff,length,ja0,ja1,i,j
  real(wp) :: scpr1,scpr0,tt
  integer :: iaseg0
  !these arrays have to be allocatable
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  !Variables for OpenMP
  !$ integer :: ithread,nthread,nchunk
  !$ integer :: omp_get_thread_num,omp_get_num_threads
  !$ integer :: ibsegs, ibsege

  keyag_c_lin = keyag_c(1,:) !speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:) !speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_wp

  !$omp parallel default (none) &
  !$omp shared (maseg_c,keyav_c,keyag_c,keyag_c_lin,keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f)&
  !$omp shared (apsi_c,bpsi_c,bpsi_f,keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f)&
  !$omp shared (apsi_f,scpr) &
!!$  !$omp parallel default(shared) &
  !$omp private(jaj,iaoff,length,ja1,ja0,jb1,jb0,iboff,scpr0,scpr1) &
  !$omp private(jbj,ibseg,iaseg0,i,j,tt)!!!,ithread,nthread,ibsegs,ibsege,nchunk)

  scpr0=0.0_wp
  scpr1=0.0_wp

!!!!start of general region

  !alternative way of parallelizing the loop, to be tested to explore performances
  !LG  ibsegs=1
  !LG  ibsege=mbseg_c
  !LG  !$ ithread=omp_get_thread_num()
  !LG  !$ nthread=omp_get_num_threads() 

  iaseg0=1 

  !coarse part. Loop on the projectors segments
  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg_c/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg_c+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg_c)

  !$omp do schedule(static)
  do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     !     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb0=max(keybg_c(1,ibseg),keyag_c_lin(1))
     jb1=keybg_c(2,ibseg) !ending point of projector segment
     iboff = max(jb0-keybg_c(1,ibseg),0)

     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt_inline(keyag_c_lin,maseg_c,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)

        do i=0,length
           tt=apsi_c(jaj+iaoff+i)
           scpr0=scpr0+tt*bpsi_c(jbj+i+iboff)
        enddo

        !call op_c_inline(length,jaj+iaoff,jbj+iboff,scpr0)

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1
  enddo
  !stop
  !$omp end do nowait

  ! fine part
  !LG  ibsegs=1
  !LG  ibsege=mbseg_f

  iaseg0=1

  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg_f/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg_f+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg_f)

  !$omp do schedule(static)
  do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     !jb0=keybg_f(1,ibseg)
     jb0=max(keybg_f(1,ibseg),keyag_f_lin(1))
     jb1=keybg_f(2,ibseg)
     iboff = max(jb0-keybg_f(1,ibseg),0)
     call hunt_inline(keyag_f_lin,maseg_f,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside
        jaj=keyav_f(iaseg0)

        do i=0,length
           do j=1,7
              tt=apsi_f(j,jaj+iaoff+i)
              scpr1=scpr1+tt*bpsi_f(j,jbj+i+iboff)
           end do
        enddo
        !call op_f_inline(length,jaj+iaoff,jbj+iboff,scpr1)

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1
  enddo
  !$omp end do !implicit barrier 

!!!!end of general region

  scpr0=scpr0+scpr1

  !$omp critical 
  scpr=scpr+scpr0
  !$omp end critical

  !$omp end parallel

!  include 'wpdot-inc.f90'

contains

!!$  pure subroutine op_c_inline(length,jaj,jbj,scpr0)
!!$    implicit none
!!$    integer, intent(in) :: length,jaj,jbj
!!$    real(wp), intent(inout) :: scpr0
!!$    !local variables
!!$    integer :: i
!!$    real(wp) :: tt
!!$    do i=0,length
!!$       tt=apsi_c(jaj+i)!iaoff+i)
!!$       scpr0=scpr0+tt*bpsi_c(jbj+i)!+iboff)
!!$    enddo
!!$  end subroutine op_c_inline
!!$
!!$  pure subroutine op_f_inline(length,jaj,jbj,scpr1)
!!$    implicit none
!!$    integer, intent(in) :: length,jaj,jbj
!!$    real(wp), intent(inout) :: scpr1
!!$    !local variables
!!$    integer :: i,j
!!$    real(wp) :: tt
!!$    do i=0,length
!!$       do j=1,7
!!$          tt=apsi_f(j,jaj+i)!iaoff+i)
!!$          scpr1=scpr1+tt*bpsi_f(j,jbj+i)!+iboff)
!!$       end do
!!$    enddo
!!$  end subroutine op_f_inline

  include 'scalar_product-inc.f90'
  
END SUBROUTINE wpdot_keys

!> Calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
!! The array mask is then constructed so that successive application of the projector on the same object can
!! be done without bitonic search
!! The array of the wavefunction is also compressed
!! so that successive projector application can be performed also with linear algebra routines
subroutine wpdot_keys_pack(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,&
     apack_c,apack_f,scpr)
  use module_defs, only: wp
  implicit none
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), intent(out) :: scpr
  real(wp), dimension(mbvctr_c), intent(inout) :: apack_c
  real(wp), dimension(7,mbvctr_f), intent(inout) :: apack_f

  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,iboff,length,ja0,ja1,i,j
  real(wp) :: scpr1,scpr0,tt
  integer :: iaseg0
  !these arrays have to be allocatable
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  !Variables for OpenMP
  !$ integer :: ithread,nthread,nchunk
  !$ integer :: omp_get_thread_num,omp_get_num_threads
  !$ integer :: ibsegs,ibsege

  keyag_c_lin = keyag_c(1,:) !speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:) !speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_wp

  !$omp parallel default (none) &
  !$omp shared (maseg_c,keyav_c,keyag_c,keyag_c_lin,keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f)&
  !$omp shared (apsi_c,bpsi_c,bpsi_f,keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f)&
  !$omp shared (apsi_f,scpr,apack_c,apack_f) &
!!$  !$omp parallel default(shared) &
  !$omp private(jaj,iaoff,length,ja1,ja0,jb1,jb0,iboff,scpr0,scpr1) &
  !$omp private(jbj,ibseg,iaseg0,i,j,tt)!!!,ithread,nthread,ibsegs,ibsege,nchunk)

  scpr0=0.0_wp
  scpr1=0.0_wp

!!!!start of general region

  !alternative way of parallelizing the loop, to be tested to explore performances
  !LG  ibsegs=1
  !LG  ibsege=mbseg_c
  !LG  !$ ithread=omp_get_thread_num()
  !LG  !$ nthread=omp_get_num_threads() 

  iaseg0=1 

  !coarse part. Loop on the projectors segments
  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg_c/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg_c+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg_c)

  !$omp do schedule(static)
  do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     !     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb0=max(keybg_c(1,ibseg),keyag_c_lin(1))
     jb1=keybg_c(2,ibseg) !ending point of projector segment
     iboff = max(jb0-keybg_c(1,ibseg),0)

     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt_inline(keyag_c_lin,maseg_c,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)

        do i=0,length
           tt=apsi_c(jaj+iaoff+i)
           scpr0=scpr0+tt*bpsi_c(jbj+i+iboff)
           apack_c(jbj+i+iboff)=tt
           !apack_c(jbj+i)=tt
        enddo
        !call op_c_inline(length,jaj+iaoff,jbj+iboff,scpr0)

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1
  enddo
  !stop
  !$omp end do nowait

  ! fine part
  !LG  ibsegs=1
  !LG  ibsege=mbseg_f

  iaseg0=1

  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg_f/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg_f+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg_f)

  !$omp do schedule(static)
  do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     !jb0=keybg_f(1,ibseg)
     jb0=max(keybg_f(1,ibseg),keyag_f_lin(1))
     jb1=keybg_f(2,ibseg)
     iboff = max(jb0-keybg_f(1,ibseg),0)
     call hunt_inline(keyag_f_lin,maseg_f,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside
        jaj=keyav_f(iaseg0)

        do i=0,length
           do j=1,7
              tt=apsi_f(j,jaj+iaoff+i)
              scpr1=scpr1+tt*bpsi_f(j,jbj+i+iboff)
              apack_f(j,jbj+i+iboff)=tt
              !apack_f(j,jbj+i)=tt
           end do
        enddo
        !call op_f_inline(length,jaj+iaoff,jbj+iboff,scpr1)

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1
  enddo
  !$omp end do !implicit barrier 

!!!!end of general region

  scpr0=scpr0+scpr1

  !$omp critical 
  scpr=scpr+scpr0
  !$omp end critical

  !$omp end parallel

!  include 'wpdot-inc.f90'

contains

!!$  subroutine op_c_inline(length,jaj,jbj,scpr0)
!!$    implicit none
!!$    integer, intent(in) :: length,jaj,jbj
!!$    real(wp), intent(inout) :: scpr0
!!$    !local variables
!!$    integer :: i
!!$    real(wp) :: tt
!!$
!!$    do i=0,length
!!$       tt=apsi_c(jaj+i)!iaoff+i)
!!$       scpr0=scpr0+tt*bpsi_c(jbj+i)!+iboff)
!!$       !apack_c(jbj+i+iboff)=tt
!!$       apack_c(jbj+i)=tt
!!$    enddo
!!$  end subroutine op_c_inline
!!$
!!$  subroutine op_f_inline(length,jaj,jbj,scpr1)
!!$    implicit none
!!$    integer, intent(in) :: length,jaj,jbj
!!$    real(wp), intent(inout) :: scpr1
!!$    !local variables
!!$    integer :: i,j
!!$    real(wp) :: tt
!!$
!!$    do i=0,length
!!$       do j=1,7
!!$          tt=apsi_f(j,jaj+i)!iaoff+i)
!!$          scpr1=scpr1+tt*bpsi_f(j,jbj+i)!+iboff)
!!$          !apack_f(j,jbj+i+iboff)=tt
!!$          apack_f(j,jbj+i)=tt
!!$       end do
!!$    enddo
!!$  end subroutine op_f_inline

  include 'scalar_product-inc.f90'
  
END SUBROUTINE wpdot_keys_pack

!> Calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
!! The array mask is then constructed so that successive application of the projector on the same object can
!! be done without bitonic search
!! The array of the wavefunction is also compressed
!! so that successive projector application can be performed also with linear algebra routines
subroutine waxpy_keys_unpack(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,&
     apack_c,apack_f,scpr)
  use module_defs, only: wp
  implicit none
  real(wp), dimension(mavctr_c), intent(inout) :: apsi_c
  real(wp), dimension(7,mavctr_f), intent(inout) :: apsi_f
  real(wp), intent(in) :: scpr
  real(wp), dimension(mbvctr_c), intent(in) :: apack_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: apack_f
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,iboff,length,ja0,ja1,i,j
  integer :: iaseg0
  real(wp) :: tt
  !these arrays have to be allocatable
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  !Variables for OpenMP
  !$ integer :: ithread,nthread,nchunk
  !$ integer :: omp_get_thread_num,omp_get_num_threads
  !$ integer :: ibsegs,ibsege

  keyag_c_lin = keyag_c(1,:) !speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:) !speed up access in hunt subroutine by consecutive arrangement in memory

  !$omp parallel default (none) &
  !$omp shared (maseg_c,keyav_c,keyag_c,keyag_c_lin,keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f)&
  !$omp shared (apsi_c,bpsi_c,bpsi_f,keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f)&
  !$omp shared (apsi_f,scpr,apack_c,apack_f) &
!!$  !$omp parallel default(shared) &
  !$omp private(jaj,iaoff,length,ja1,ja0,jb1,jb0,iboff,i,j,tt) &
  !$omp private(jbj,ibseg,iaseg0)!!!,ithread,nthread,ibsegs,ibsege,nchunk)

  !alternative way of parallelizing the loop, to be tested to explore performances
  !LG  ibsegs=1
  !LG  ibsege=mbseg_c
  !LG  !$ ithread=omp_get_thread_num()
  !LG  !$ nthread=omp_get_num_threads() 

  iaseg0=1 

  !coarse part. Loop on the projectors segments
  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg_c/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg_c+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg_c)

  !$omp do schedule(static)
  do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     !     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb0=max(keybg_c(1,ibseg),keyag_c_lin(1))
     jb1=keybg_c(2,ibseg) !ending point of projector segment
     iboff = max(jb0-keybg_c(1,ibseg),0)

     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt_inline(keyag_c_lin,maseg_c,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)

        do i=0,length
           !tt=bpsi_c(jbj+i)
           tt=apack_c(jbj+i+iboff)+scpr*bpsi_c(jbj+i+iboff)
           apsi_c(jaj+i+iaoff)=apsi_c(jaj+i+iaoff)+tt
        enddo
!        call op_c_inline(length,jaj+iaoff,jbj+iboff,scpr)

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1
  enddo
  !stop
  !$omp end do nowait

  ! fine part
  !LG  ibsegs=1
  !LG  ibsege=mbseg_f

  iaseg0=1

  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg_f/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg_f+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg_f)

  !$omp do schedule(static)
  do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     !jb0=keybg_f(1,ibseg)
     jb0=max(keybg_f(1,ibseg),keyag_f_lin(1))
     jb1=keybg_f(2,ibseg)
     iboff = max(jb0-keybg_f(1,ibseg),0)
     call hunt_inline(keyag_f_lin,maseg_f,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside
        jaj=keyav_f(iaseg0)

        do i=0,length
           do j=1,7
              tt=apack_f(j,jbj+i+iboff)+scpr*bpsi_f(j,jbj+i+iboff)
              apsi_f(j,jaj+i+iaoff)=apsi_f(j,jaj+i+iaoff)+tt
           end do
        enddo
        !call op_f_inline(length,jaj+iaoff,jbj+iboff,scpr)

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1
  enddo
  !$omp end do !implicit barrier 
  !$omp end parallel


contains

!!$  subroutine op_c_inline(length,jaj,jbj,scpr)
!!$    implicit none
!!$    integer, intent(in) :: length,jaj,jbj
!!$    real(wp), intent(in) :: scpr
!!$    !local variables
!!$    integer :: i
!!$    real(wp) :: tt
!!$
!!$     do i=0,length
!!$        !tt=bpsi_c(jbj+i)
!!$        tt=apack_c(jbj+i)+scpr*bpsi_c(jbj+i)
!!$        apsi_c(jaj+i)=apsi_c(jaj+i)+tt
!!$     enddo
!!$
!!$  end subroutine op_c_inline
!!$
!!$  subroutine op_f_inline(length,jaj,jbj,scpr)
!!$    implicit none
!!$    integer, intent(in) :: length,jaj,jbj
!!$    real(wp), intent(in) :: scpr
!!$    !local variables
!!$    integer :: i,j
!!$    real(wp) :: tt
!!$
!!$     do i=0,length
!!$        do j=1,7
!!$           tt=apack_f(j,jbj+i)+scpr*bpsi_f(j,jbj+i)
!!$           apsi_f(j,jaj+i)=apsi_f(j,jaj+i)+tt
!!$        end do
!!$     enddo
!!$  end subroutine op_f_inline

  include 'scalar_product-inc.f90'
  
END SUBROUTINE waxpy_keys_unpack


!> Calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
!! The array mask is used so that application of the projector can
!! be done without bitonic search
!! The array of the wavefunction is also compressed
!! so that successive projector application can be performed also with linear algebra routines
subroutine wpdot_mask_pack(  &
     mavctr_c,mavctr_f,mseg_c,mseg_f,amask_c,amask_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,bpsi_c,bpsi_f,&
     apack_c,apack_f,scpr)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,mseg_c,mseg_f,mbvctr_c,mbvctr_f
  integer, dimension(3,mseg_c), intent(in) :: amask_c
  integer, dimension(3,mseg_f), intent(in) :: amask_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  real(wp), dimension(mbvctr_c), intent(inout) :: apack_c
  real(wp), dimension(7,mbvctr_f), intent(inout) :: apack_f
  real(wp), intent(out) :: scpr
  !local variables
  integer :: iaseg,jaj,jbj,length,i,j
  real(wp) :: scpr1,scpr0,tt

  scpr=0.0_wp

  !$omp parallel default(shared) &
  !$omp private(i,jaj,length,tt,j,scpr0,scpr1) &
  !$omp private(jbj,iaseg)

  scpr0=0.0_wp
  scpr1=0.0_wp

  !$omp do !schedule(static)
  do iaseg=1,mseg_c
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_c(1,iaseg) !number of elements to be copied
     jaj   =amask_c(2,iaseg) !starting point in original array
     jbj   =amask_c(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        tt=apsi_c(jaj+i)
        scpr0=scpr0+tt*bpsi_c(jbj+i)
        apack_c(jbj+i)=tt
     enddo
  end do
  !$omp end do nowait

  !$omp do !schedule(static)
  do iaseg=1,mseg_f
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_f(1,iaseg) !number of elements to be copied
     jaj   =amask_f(2,iaseg) !starting point in original array
     jbj   =amask_f(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        do j=1,7
           tt=apsi_f(j,jaj+i)
           scpr1=scpr1+tt*bpsi_f(j,jbj+i)
           apack_f(j,jbj+i)=tt
        end do
     enddo
  end do
  !$omp end do !implicit barrier 

  scpr0=scpr0+scpr1

  !$omp critical 
  scpr=scpr+scpr0
  !$omp end critical

  !$omp end parallel
  
END SUBROUTINE wpdot_mask_pack

!> Calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
!! The array mask is used so that application of the projector can
!! be done without bitonic search
subroutine wpdot_mask(  &
     mavctr_c,mavctr_f,mseg_c,mseg_f,amask_c,amask_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,bpsi_c,bpsi_f,&
     scpr)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,mseg_c,mseg_f,mbvctr_c,mbvctr_f
  integer, dimension(3,mseg_c), intent(in) :: amask_c
  integer, dimension(3,mseg_f), intent(in) :: amask_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  real(wp), intent(out) :: scpr
  !local variables
  integer :: iaseg,jaj,jbj,length,i,j
  real(wp) :: scpr1,scpr0,tt

  scpr=0.0_wp

  !$omp parallel default(shared) &
  !$omp private(i,jaj,length,tt,j,scpr0,scpr1) &
  !$omp private(jbj,iaseg)

  scpr0=0.0_wp
  scpr1=0.0_wp

  !$omp do !schedule(static)
  do iaseg=1,mseg_c
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_c(1,iaseg) !number of elements to be copied
     jaj   =amask_c(2,iaseg) !starting point in original array
     jbj   =amask_c(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        tt=apsi_c(jaj+i)
        scpr0=scpr0+tt*bpsi_c(jbj+i)
     enddo
  end do
  !$omp end do nowait

  !$omp do !schedule(static)
  do iaseg=1,mseg_f
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_f(1,iaseg) !number of elements to be copied
     jaj   =amask_f(2,iaseg) !starting point in original array
     jbj   =amask_f(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        do j=1,7
           tt=apsi_f(j,jaj+i)
           scpr1=scpr1+tt*bpsi_f(j,jbj+i)
        end do
     enddo
  end do
  !$omp end do !implicit barrier 

  scpr0=scpr0+scpr1

  !$omp critical 
  scpr=scpr+scpr0
  !$omp end critical

  !$omp end parallel
  
END SUBROUTINE wpdot_mask


!> Rank 1 update of wavefunction a with wavefunction b: apsi=apsi+scpr*bpsi
!! The update is only done in the localization region of apsi
!! The array mask is used so that application of the projector can
!! be done without bitonic search
subroutine waxpy_mask(  &
     mavctr_c,mavctr_f,mseg_c,mseg_f,amask_c,amask_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,bpsi_c,bpsi_f,scpr)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,mseg_c,mseg_f,mbvctr_c,mbvctr_f
  real(wp), intent(in) :: scpr
  integer, dimension(3,mseg_c), intent(in) :: amask_c
  integer, dimension(3,mseg_f), intent(in) :: amask_f
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  real(wp), dimension(mavctr_c), intent(inout) :: apsi_c
  real(wp), dimension(7,mavctr_f), intent(inout) :: apsi_f

  !local variables
  integer :: iaseg,jaj,jbj,length,i,j
  real(wp) :: tt

  !quick return if possible
  if (scpr==0.0_wp) return

  !$omp parallel default(shared) &
  !$omp private(i,jaj,length,tt,j) &
  !$omp private(jbj,iaseg)

  !$omp do !schedule(static)
  do iaseg=1,mseg_c
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_c(1,iaseg) !number of elements to be copied
     jaj   =amask_c(2,iaseg) !starting point in original array
     jbj   =amask_c(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        tt=bpsi_c(jbj+i)
        apsi_c(jaj+i)=apsi_c(jaj+i)+scpr*tt
     enddo
  end do
  !$omp end do nowait

  !$omp do !schedule(static)
  do iaseg=1,mseg_f
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_f(1,iaseg) !number of elements to be copied
     jaj   =amask_f(2,iaseg) !starting point in original array
     jbj   =amask_f(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        do j=1,7
           tt=bpsi_f(j,jbj+i)
           apsi_f(j,jaj+i)=apsi_f(j,jaj+i)+scpr*tt
        end do
     enddo
  end do
  !$omp end do !implicit barrier 

  !$omp end parallel
  
END SUBROUTINE waxpy_mask

!> Rank 1 update of wavefunction a with wavefunction b: apsi=apsi+scpr*bpsi
!! The update is only done in the localization region of apsi
!! The array mask is used so that application of the projector can
!! be done without bitonic search
subroutine waxpy_mask_unpack(  &
     mavctr_c,mavctr_f,mseg_c,mseg_f,amask_c,amask_f,apack_c,apack_f,&
     apsi_c,apsi_f, &
     mbvctr_c,mbvctr_f,bpsi_c,bpsi_f,scpr)
  use module_defs, only: wp
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,mseg_c,mseg_f,mbvctr_c,mbvctr_f
  real(wp), intent(in) :: scpr
  integer, dimension(3,mseg_c), intent(in) :: amask_c
  integer, dimension(3,mseg_f), intent(in) :: amask_f
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  real(wp), dimension(mbvctr_c), intent(in) :: apack_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: apack_f
  real(wp), dimension(mavctr_c), intent(inout) :: apsi_c
  real(wp), dimension(7,mavctr_f), intent(inout) :: apsi_f
  !local variables
  integer :: iaseg,jaj,jbj,length,i,j
  real(wp) :: tt

  !$omp parallel default(shared) &
  !$omp private(i,jaj,length,tt,j) &
  !$omp private(jbj,iaseg)

  !$omp do !schedule(static)
  do iaseg=1,mseg_c
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_c(1,iaseg) !number of elements to be copied
     jaj   =amask_c(2,iaseg) !starting point in original array
     jbj   =amask_c(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        !tt=bpsi_c(jbj+i)
        tt=apack_c(jbj+i)+scpr*bpsi_c(jbj+i)
        apsi_c(jaj+i)=apsi_c(jaj+i)+tt
     enddo
  end do
  !$omp end do nowait

  !$omp do !schedule(static)
  do iaseg=1,mseg_f
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_f(1,iaseg) !number of elements to be copied
     jaj   =amask_f(2,iaseg) !starting point in original array
     jbj   =amask_f(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        do j=1,7
           !tt=bpsi_f(j,jbj+i)
           tt=apack_f(j,jbj+i)+scpr*bpsi_f(j,jbj+i)
           apsi_f(j,jaj+i)=apsi_f(j,jaj+i)+tt
        end do
     enddo
  end do
  !$omp end do !implicit barrier 

  !$omp end parallel
  
END SUBROUTINE waxpy_mask_unpack

!> find the number of chunks which are needed to perform blas operations among two compressed wavefunctions
subroutine count_wblas_segs(maseg,mbseg,keyag_lin,keyag,keybg,nbsegs)
  implicit none
  integer, intent(in) :: maseg,mbseg
  integer, dimension(maseg), intent(in) :: keyag_lin !>linear version of second indices of keyag
  integer, dimension(2,maseg), intent(in) :: keyag !>values of the keys ordered for compression a
  integer, dimension(2,mbseg), intent(in) :: keybg !>values of the keys ordered for compression b
  integer, dimension(mbseg), intent(inout) :: nbsegs !>number of common segments for each segment of b
  !local variables
  integer :: ibseg,jb1,jb0,length,ja0,ja1,imask
  integer :: iaseg0
  !Variables for OpenMP
  !$ integer :: ibsegs,ibsege

  !$omp parallel default (none) &
  !$omp shared (maseg,keyag,keyag_lin,keybg,mbseg,nbsegs)&
!!$  !$omp parallel default(shared) &
  !$omp private(length,ja1,ja0,jb1,jb0,imask) &
  !$omp private(ibseg,iaseg0)!!!,ithread,nthread,ibsegs,ibsege,nchunk)

  !alternative way of parallelizing the loop, to be tested to explore performances
  !LG  ibsegs=1
  !LG  ibsege=mbseg
  !LG  !$ ithread=omp_get_thread_num()
  !LG  !$ nthread=omp_get_num_threads() 

  iaseg0=1 

  !coarse part. Loop on the projectors segments
  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg)

  !$omp do schedule(static)
  do ibseg=1,mbseg
     jb0=max(keybg(1,ibseg),keyag_lin(1)) !starting point of projector segment
     jb1=keybg(2,ibseg) !ending point of projector segment

     !first passage
     imask=0

     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt_inline(keyag_lin,maseg,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop: do while(iaseg0 <= maseg)
        ja0=keyag_lin(iaseg0)
        ja1=min(jb1,keyag(2,iaseg0)) 
        length = ja1-jb0

        !count the active segments for this ibseg
        if (length+1 > 0) then
           imask=imask+1
        end if
        
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg) exit nonconvex_loop !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg(1,ibseg),keyag_lin(iaseg0))
        if (keyag_lin(iaseg0)>jb1) exit nonconvex_loop !segment is not covered
     end do nonconvex_loop

     !first passage, fill the segments array
     nbsegs(ibseg)=imask

     !disable loop if the end is reached
     if (iaseg0 == maseg .and. keybg(1,ibseg)> keyag_lin(maseg)) iaseg0=iaseg0+1
  enddo
  !$omp end do 
  !$omp end parallel

  contains

    include 'scalar_product-inc.f90'

end subroutine count_wblas_segs

!> find the number of chunks which are needed to perform blas operations among two compressed wavefunctions
subroutine fill_wblas_segs(maseg,mbseg,mask_segs,isegs_offset,keyag_lin,keyag,keybg,keyav,keybv,amask)
  implicit none
  integer, intent(in) :: maseg,mbseg,mask_segs
  integer, dimension(maseg), intent(in) :: keyav !>position of the segments in compressed a storage
  integer, dimension(mbseg), intent(in) :: keybv !>position of the segments in compressed b storage
  integer, dimension(maseg), intent(in) :: keyag_lin !>linear version of second indices of keyag
  integer, dimension(mbseg), intent(in) :: isegs_offset !>displacement in common segments for each segment of b. 
                                                            !! Integral function of nbsegs vector
  integer, dimension(2,maseg), intent(in) :: keyag !>values of the keys ordered for compression a
  integer, dimension(2,mbseg), intent(in) :: keybg !>values of the keys ordered for compression b
  integer, dimension(3,mask_segs), intent(out) :: amask !>masking array and positions in compressed storage
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,iboff,length,ja0,ja1,imask
  integer :: iaseg0
  !Variables for OpenMP
  !$ integer :: ibsegs,ibsege

  !$omp parallel default (none) &
  !$omp shared (maseg,keyav,keyag,keyag_lin,keybg,mbseg,keybv,isegs_offset,amask)&
!!$  !$omp parallel default(shared) &
  !$omp private(jaj,iaoff,length,ja1,ja0,jb1,jb0,iboff,imask) &
  !$omp private(jbj,ibseg,iaseg0)!!!,ithread,nthread,ibsegs,ibsege,nchunk)

  !alternative way of parallelizing the loop, to be tested to explore performances
  !LG  ibsegs=1
  !LG  ibsege=mbseg
  !LG  !$ ithread=omp_get_thread_num()
  !LG  !$ nthread=omp_get_num_threads() 

  iaseg0=1 

  !coarse part. Loop on the projectors segments
  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg)

  !$omp do schedule(static)
  do ibseg=1,mbseg
     jbj=keybv(ibseg)
     jb0=max(keybg(1,ibseg),keyag_lin(1)) !starting point of projector segment
     jb1=keybg(2,ibseg) !ending point of projector segment
     iboff = max(jb0-keybg(1,ibseg),0)

     !second passage, retrieve starting point
     imask=isegs_offset(ibseg)

     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt_inline(keyag_lin,maseg,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop: do while(iaseg0 <= maseg)
        !length = jb1-jb0
        !iaoff = jb0-keyag_lin(iaseg0)!jb0-ja0

        ja0=keyag_lin(iaseg0)
        ja1=min(jb1,keyag(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav(iaseg0)
        
        !second passage: fill the masking elements as they have to be used
        if (length+1 > 0) then
           imask=imask+1
           !with the masking array there should be no need to perform the bitonic search anymore
           amask(1,imask)=length+1  !number of elements to be copied
           amask(2,imask)=jaj+iaoff !starting point in original array
           amask(3,imask)=jbj+iboff !starting point in packed array
        end if

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg) exit nonconvex_loop !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg(1,ibseg),keyag_lin(iaseg0))
        if (keyag_lin(iaseg0)>jb1) exit nonconvex_loop !segment is not covered
        jbj=jbj+max(jb0-keybg(1,ibseg),0)
     end do nonconvex_loop
     !disable loop if the end is reached
     if (iaseg0 == maseg .and. keybg(1,ibseg)> keyag_lin(maseg)) iaseg0=iaseg0+1
  enddo
  !$omp end do 
  !$omp end parallel
 contains

    include 'scalar_product-inc.f90'

  end subroutine fill_wblas_segs
