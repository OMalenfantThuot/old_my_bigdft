!> @file
!!  Routines to create the localisation region
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Determine a set of localisation regions from the centers and the radii.
!! cut in cubes the global reference system
!!subroutine determine_locreg(nlr,cxyz,locrad,hx,hy,hz,Glr,Llr)
!!  use module_base
!!  use module_types
!!  implicit none
!!  integer, intent(in) :: nlr
!!  real(gp), intent(in) :: hx,hy,hz
!!  type(locreg_descriptors), intent(in) :: Glr
!!  real(gp), dimension(nlr), intent(in) :: locrad
!!  real(gp), dimension(3,nlr), intent(in) :: cxyz
!!  type(locreg_descriptors), dimension(nlr), intent(out) :: Llr
!!  !local variables
!!  character(len=*), parameter :: subname='determine_locreg'
!!  logical :: perx,pery,perz
!!  integer :: ilr,isx,isy,isz,iex,iey,iez
!!  real(gp) :: rx,ry,rz,cutoff
!!
!!  !determine the limits of the different localisation regions
!!  do ilr=1,nlr
!!
!!     rx=cxyz(1,ilr)
!!     ry=cxyz(2,ilr)
!!     rz=cxyz(3,ilr)
!!
!!     cutoff=locrad(ilr)
!!
!!     isx=floor((rx-cutoff)/hx)
!!     isy=floor((ry-cutoff)/hy)
!!     isz=floor((rz-cutoff)/hz)
!!
!!     iex=ceiling((rx+cutoff)/hx)
!!     iey=ceiling((ry+cutoff)/hy)
!!     iez=ceiling((rz+cutoff)/hz)
!!
!!     !assign the geometric code to the localisation region
!!     select case(Glr%geocode)
!!     case('F')
!!        isx=max(isx,0)
!!        isy=max(isy,0)
!!        isz=max(isz,0)
!!        
!!        iex=min(iex,Glr%d%n1)
!!        iey=min(iey,Glr%d%n2)
!!        iez=min(iez,Glr%d%n3)
!!
!!        perx=.false.
!!        pery=.false.
!!        perz=.false.
!!        Llr(ilr)%geocode='F'
!!     case('S')
!!        if (iex - isx >= Glr%d%n1 - 14) then
!!           isx=0
!!           iex=Glr%d%n1
!!           perx=.true.
!!        else
!!           isx=modulo(isx,Glr%d%n1+1)
!!           iex=modulo(iex,Glr%d%n1+1)
!!           perx=.false.
!!        end if
!!
!!        isy=max(isy,0)
!!        iey=min(iey,Glr%d%n2)
!!        pery=.false.
!!
!!        !control the geometric code for the localisation region
!!        if (iez - isz >= Glr%d%n3 - 14) then
!!           isz=0
!!           iez=Glr%d%n3
!!           perz=.true.
!!        else
!!           isz=modulo(isz,Glr%d%n3+1)
!!           iez=modulo(iez,Glr%d%n3+1)
!!           perz=.false.
!!        end if
!!
!!        if (perx .and. perz) then
!!           Llr(ilr)%geocode='S'
!!        else if (.not.perx .and. .not.perz) then
!!           Llr(ilr)%geocode='F'
!!        else
!!           write(*,*)'ERROR: localisation region geometry not allowed'
!!           stop
!!        end if
!!
!!     case('P')
!!        if (iex - isx >= Glr%d%n1 - 14) then
!!           isx=0
!!           iex=Glr%d%n1
!!           perx=.true.
!!        else
!!           isx=modulo(isx,Glr%d%n1+1)
!!           iex=modulo(iex,Glr%d%n1+1)
!!           perx=.false.
!!        end if
!!        if (iey - isy >= Glr%d%n2 - 14) then
!!           isy=0
!!           iey=Glr%d%n2
!!           pery=.true.
!!        else
!!           isy=modulo(isy,Glr%d%n2+1)
!!           iey=modulo(iey,Glr%d%n2+1)
!!           pery=.false.
!!        end if
!!        if (iez - isz >= Glr%d%n3 - 14) then
!!           isz=0
!!           iez=Glr%d%n3
!!           perz=.true.
!!        else
!!           isz=modulo(isz,Glr%d%n3+1)
!!           iez=modulo(iez,Glr%d%n3+1)
!!           perz=.false.
!!        end if
!!        if (perx .and. perz .and. pery) then
!!           Llr(ilr)%geocode='P'
!!        else if (.not.perx .and. .not.perz .and. .not. pery) then
!!           Llr(ilr)%geocode='F'
!!        else if (perx .and. perz .and. .not. pery) then
!!           Llr(ilr)%geocode='S'
!!        else
!!           write(*,*)'ERROR: localisation region geometry not allowed'
!!           stop
!!        end if
!!
!!     end select
!!     
!!     !values for the starting point of the cube
!!     Llr(ilr)%ns1=isx
!!     Llr(ilr)%ns2=isy
!!     Llr(ilr)%ns3=isz
!!     !dimensions of the localisation region
!!     Llr(ilr)%d%n1=iex-isx
!!     Llr(ilr)%d%n2=iey-isy
!!     Llr(ilr)%d%n3=iez-isz
!!
!!     !dimensions of the fine grid inside the localisation region
!!     if (isx < iex) then
!!        Llr(ilr)%d%nfl1=max(isx,Glr%d%nfl1)-isx
!!        Llr(ilr)%d%nfu1=min(iex,Glr%d%nfu1)-isx
!!     else
!!        write(*,*)'Yet to be implemented (little effort?)'
!!        stop
!!     end if
!!
!!     if (isy < iey) then
!!        Llr(ilr)%d%nfl2=max(isy,Glr%d%nfl2)-isy
!!        Llr(ilr)%d%nfu2=min(iey,Glr%d%nfu2)-isy
!!     else
!!        write(*,*)'Yet to be implemented (little effort?)'
!!        stop
!!     end if
!!
!!     if (isz < iez) then
!!        Llr(ilr)%d%nfl3=max(isz,Glr%d%nfl3)-isz
!!        Llr(ilr)%d%nfu3=min(iez,Glr%d%nfu3)-isz
!!     else
!!        write(*,*)'Yet to be implemented (little effort?)'
!!        stop
!!     end if
!!     
!!     !dimensions of the interpolating scaling functions grid
!!     select case(Llr(ilr)%geocode)
!!     case('F')
!!        Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+31
!!        Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+31
!!        Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+31
!!     case('S')
!!        Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+2
!!        Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+31
!!        Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+2
!!     case('P')
!!        Llr(ilr)%d%n1i=2*Llr(ilr)%d%n1+2
!!        Llr(ilr)%d%n2i=2*Llr(ilr)%d%n2+2
!!        Llr(ilr)%d%n3i=2*Llr(ilr)%d%n3+2
!!     end select
!!
!!     
!!     !define the wavefunction descriptors inside the localisation region
!!     !calculate the number of point and segments for local localisation regions
!!     !coarse part
!!     call num_segkeys_loc(Glr%d%n1,Glr%d%n2,isx,iex,isy,iey,isz,iez,&
!!          Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
!!          Llr(ilr)%wfd%nseg_c,Llr(ilr)%wfd%nvctr_c)
!!     !fine part
!!     call num_segkeys_loc(Glr%d%n1,Glr%d%n2,isx,iex,isy,iey,isz,iez,&
!!          Glr%wfd%nseg_f,Glr%wfd%nvctr_f,&
!!          Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!          Llr(ilr)%wfd%nseg_f,Llr(ilr)%wfd%nvctr_f)
!!
!!     !allocate the wavefunction descriptors following the needs
!!     call allocate_wfd(Llr(ilr)%wfd,subname)
!!
!!     !fill such descriptors
!!     !coarse part
!!     call segkeys_loc(Glr%d%n1,Glr%d%n2,isx,iex,isy,iey,isz,iez,& !n(m)
!!          Glr%wfd%nseg_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1),&
!!          Llr(ilr)%wfd%nseg_c,Llr(ilr)%wfd%nvctr_c,&
!!          Llr(ilr)%wfd%keyg(1,1),Llr(ilr)%wfd%keyv(1))
!!     !fine part
!!     call segkeys_loc(Glr%d%n1,Glr%d%n2,isx,iex,isy,iey,isz,iez,& !n(m) 
!!          Glr%wfd%nseg_f,&
!!          Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!          Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),&
!!          Llr(ilr)%wfd%nseg_f,Llr(ilr)%wfd%nvctr_f,&
!!          Llr(ilr)%wfd%keyg(1,Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)),&
!!          Llr(ilr)%wfd%keyv(Llr(ilr)%wfd%nseg_c+min(1,Llr(ilr)%wfd%nseg_f)))
!!
!!     !if the localisation region is isolated build also the bounds
!!     if (Llr(ilr)%geocode=='F') then
!!        call locreg_bounds(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,&
!!             Llr(ilr)%d%nfl1,Llr(ilr)%d%nfu1,Llr(ilr)%d%nfl2,Llr(ilr)%d%nfu2,&
!!             Llr(ilr)%d%nfl3,Llr(ilr)%d%nfu3,Llr(ilr)%wfd,Llr(ilr)%bounds)
!!     end if
!!
!!  end do
!!
!!  !after all localisation regions are determined draw them
!!  !call draw_locregs(nlr,hx,hy,hz,Llr)
!!
!!END SUBROUTINE determine_locreg


subroutine draw_locregs(nlr,hx,hy,hz,Llr)
  use module_base
  use locregs
  use bounds, only: wfd_to_logrids
  implicit none
  integer, intent(in) :: nlr
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), dimension(nlr), intent(in) :: Llr
  !local variables
  character(len=*), parameter :: subname='draw_locregs'
  character(len=4) :: message
  integer :: i1,i2,i3,ilr,nvctr_tot
  logical, dimension(:,:,:), allocatable :: logrid_c,logrid_f

  !calculate total number
  nvctr_tot=0
  do ilr=1,nlr
     nvctr_tot=nvctr_tot+Llr(ilr)%wfd%nvctr_c+Llr(ilr)%wfd%nvctr_f
  end do

  !open file for writing
  open(unit=22,file='locregs.xyz',status='unknown')
  write(22,*) nvctr_tot,' atomic'
  write(22,*)'coarse and fine points of all the different localisation regions'

  do ilr=1,nlr
     !define logrids
     logrid_c = f_malloc((/ 0.to.Llr(ilr)%d%n1, 0.to.Llr(ilr)%d%n2, 0.to.Llr(ilr)%d%n3 /),id='logrid_c')
     logrid_f = f_malloc((/ 0.to.Llr(ilr)%d%n1, 0.to.Llr(ilr)%d%n2, 0.to.Llr(ilr)%d%n3 /),id='logrid_f')

     call wfd_to_logrids(Llr(ilr)%d%n1,Llr(ilr)%d%n2,Llr(ilr)%d%n3,Llr(ilr)%wfd,&
          logrid_c,logrid_f)

     write(message,'(1a,i0)')'g',ilr
     do i3=0,Llr(ilr)%d%n3  
        do i2=0,Llr(ilr)%d%n2  
           do i1=0,Llr(ilr)%d%n1
              if (logrid_c(i1,i2,i3))&
                   write(22,'(a4,2x,3(1x,e10.3))') &
                   message,real(i1+Llr(ilr)%ns1,gp)*hx,&
                   real(i2+Llr(ilr)%ns2,gp)*hy,real(i3+Llr(ilr)%ns3,gp)*hz
           enddo
        enddo
     end do
     write(message,'(1a,i0)')'G',ilr
     do i3=0,Llr(ilr)%d%n3 
        do i2=0,Llr(ilr)%d%n2 
           do i1=0,Llr(ilr)%d%n1
              if (logrid_f(i1,i2,i3))&
                   write(22,'(a4,2x,3(1x,e10.3))') &
                   message,real(i1+Llr(ilr)%ns1,gp)*hx,&
                   real(i2+Llr(ilr)%ns2,gp)*hy,real(i3+Llr(ilr)%ns3,gp)*hz
           enddo
        enddo
     enddo


     call f_free(logrid_c)
     call f_free(logrid_f)
  end do

  !close file for writing
  close(unit=22)  
END SUBROUTINE draw_locregs


