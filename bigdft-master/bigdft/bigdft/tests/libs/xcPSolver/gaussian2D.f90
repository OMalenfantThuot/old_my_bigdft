!> @file
!!  Use integral form for Poisson solver
!! @author
!!    Copyright (c) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Program testing new ideas like momentum-preserving gaussian integrals 2D
program MP_gaussian
  use module_base
  use gaussians
  use multipole_preserving
  use yaml_output
  implicit none
  integer, parameter :: iunit=16        !< File unit for the plot
  integer, parameter :: nmoms=16        !< Number of calculated moments
  integer, parameter :: nstep=10        !< Number of resolution to calculate the moments
  integer, parameter :: nsigma=1        !< Number of different gaussian functions
  integer, parameter :: npts=50         !< Arrays from -npts to npts
  real(gp), parameter :: hgrid = 1.0_gp !< Grid step
  real(gp), parameter :: sigma = 0.25_gp !< Sigma gaussian
  integer :: j
  integer :: imoms,pow,istep,isigma
  real(gp) :: pgauss,x0,y0,reference,max_fj
  real(gp), dimension(0:nmoms,2) :: moments
  real(gp), dimension(3,2,0:nmoms) :: avgmaxmin
  real(gp), dimension(:,:), allocatable :: fj_phi,fj_coll


  call f_lib_initialize()

  pow=0

  !pgauss=0.5_gp/((0.1_gp*hgrid)**2)!8.0e-3_dp*1.25_dp**(6*(8-1))
  !array where we have to write the value of the discretization
  fj_phi=f_malloc( (/ 0 .to. 2*npts, 0 .to. 2*npts/), id='fj_phi')
  fj_coll=f_malloc((/ 0 .to. 2*npts, 0 .to. 2*npts/), id='fj_coll')
  call initialize_real_space_conversion() !initialize the work arrays needed to integrate with isf

  ! Calculate for different nsigma sigma
  do isigma=1,nsigma
     pgauss=0.5_gp/((sigma+0.01_gp*(isigma-1)*hgrid)**2)
     call yaml_map('sigma/h',sqrt(0.5_gp/pgauss)/hgrid)
     !plot raw function (fort.iunit)
     do j=-npts,2*npts
        if (pow /= 0) then
           write(iunit,*) j,exp(-pgauss*(j*hgrid)**2)*((j*hgrid)**pow)
        else
           write(iunit,*) j,exp(-pgauss*(j*hgrid)**2)
        end if
     end do

     avgmaxmin=0.d0
     avgmaxmin(3,:,:)=1.d100
     max_fj=0.0_gp
     do istep=1,nstep
        x0=(-0.5_gp+real(istep-1,gp)/real(nstep,gp)+real(npts,gp))*hgrid
        y0=(-0.5_gp+real(istep-1,gp)/real(nstep,gp)+real(npts,gp))*hgrid
        !call yaml_map('x0',x0,advance='no')
        !call yaml_comment('Step No.'//trim(yaml_toa(istep)),tabbing=70)
        call evaluate_moments2D(nmoms,npts,hgrid,pgauss,pow,x0,y0,fj_phi,fj_coll,moments)
        max_fj=max(max_fj,maxval(abs(fj_coll-fj_phi)))
!!$  !print moments value
!!$  do imoms=0,nmoms
!!$     reference=gauint0(pgauss,imoms+pow)
!!$     if (reference /=0.0_gp) then
!!$        call yaml_map('Mom No.'//trim(yaml_toa(imoms)),&
!!$             (moments(imoms,:)-reference)/reference,fmt='(1pe22.14)',advance='no')
!!$     else
!!$        call yaml_map('Mom No.'//trim(yaml_toa(imoms)),&
!!$             moments(imoms,:),fmt='(1pe22.14)',advance='no')
!!$     end if
!!$     call yaml_comment('Ref: '//trim(yaml_toa(reference,fmt='(1pe22.14)')))
!!$  end do

        !calculate the average, maximum and minimum of each moment in function of the reference
        !j=1 use the elemental property of the mp_exp function with fj_phi
        !j=2 collocation array with fj_coll
        do j=1,2
           do imoms=0,nmoms
              reference=gauint0(pgauss,imoms+pow)**2
              print *,j,imoms,reference,moments(imoms,j)
              if (reference /= 0.0_gp) then
                 !x^even
                 moments(imoms,j) = abs((moments(imoms,j)-reference))!/reference)
              else
                 !x^odd
                 moments(imoms,j) = abs(moments(imoms,j))
              end if
              avgmaxmin(1,j,imoms) = avgmaxmin(1,j,imoms)+moments(imoms,j)/real(nstep,gp)
              avgmaxmin(2,j,imoms) = max(moments(imoms,j),avgmaxmin(2,j,imoms))
              avgmaxmin(3,j,imoms) = min(moments(imoms,j),avgmaxmin(3,j,imoms))
           end do
        end do
     end do

     !Plot fort.(iunit+1)
     write(iunit+1,'(104(1pe14.5))') sqrt(0.5_gp/pgauss)/hgrid,avgmaxmin
     call yaml_map('maxdiff' // trim(yaml_toa(isigma)), (/ sqrt(0.5_gp/pgauss)/hgrid, max_fj /) )
     !print *,'maxdiff',sqrt(0.5_gp/pgauss)/hgrid,max_fj
  end do

  call yaml_map('Results',reshape(avgmaxmin,(/6,nmoms+1/)),fmt='(1pe14.5)')

  call finalize_real_space_conversion()

  call f_free(fj_phi)
  call f_free(fj_coll)
  call f_lib_finalize()

end program MP_gaussian


!> Classify the quality of a multipole extraction in both cases
subroutine evaluate_moments2D(nmoms,npts,hgrid,pgauss,pow,x0,y0,fj_phi,fj_coll,moments)
  use module_base, only: gp
  use multipole_preserving, only: mp_exp
  implicit none
  !Arguments
  integer, intent(in) :: npts,pow,nmoms
  real(gp), intent(in) :: hgrid,pgauss,x0,y0
  real(gp), dimension(0:nmoms,2), intent(out) :: moments
  real(gp), dimension(0:2*npts,0:2*npts), intent(out) :: fj_phi,fj_coll
  !local variables
  integer :: j,jy

  !use the elemental property of the mp_exp function
  do jy=0,2*npts
     fj_phi(:,jy)=mp_exp(hgrid,x0,pgauss,(/(j,j=0,2*npts)/),pow,.true.)*mp_exp(hgrid,y0,pgauss,jy,pow,.true.)
  end do
  !scfdotf((/(j,j=-npts,npts)/),hgrid,pgauss,x0,pow)
  call moments_2d(2*npts+1,2*npts+1,fj_phi,x0,y0,hgrid,nmoms,moments(0,1))

  !collocation array
  do jy=0,2*npts
     fj_coll(:,jy)=mp_exp(hgrid,x0,pgauss,(/(j,j=0,2*npts)/),pow,.false.)*mp_exp(hgrid,y0,pgauss,jy,pow,.false.)
  end do
  !if (pow /=0) then
  !   fj_coll=(/(exp(-pgauss*(j*hgrid-x0)**2)*(j*hgrid-x0)**pow,j=-npts,npts)/)
  !else
  !   fj_coll=(/(exp(-pgauss*(j*hgrid-x0)**2),j=-npts,npts)/)
  !end if
  call moments_2d(2*npts+1,2*npts+1,fj_coll,x0,y0,hgrid,nmoms,moments(0,2))

end subroutine evaluate_moments2D


!> Calculate the moments of an array with respect to a reference point
subroutine moments_2d(nx,ny,array,x0,y0,h,nmoms,moments)
  use module_base, only: gp
  implicit none
  !Arguments
  integer, intent(in) :: nmoms,nx,ny
  real(gp), intent(in) :: x0,y0,h !< grid spacing
  real(gp), dimension(nx,ny), intent(in) :: array
  real(gp), dimension(0:nmoms), intent(out) :: moments
  !local variables
  integer :: j,kx,ky
  real(gp) :: x,y

  do j=0,nmoms
     moments(j)=0.0_gp
     do kx=1,nx
        x=real(kx-1,gp)*h-x0
        do ky =1,ny
           y=real(ky-1,gp)*h-y0
           moments(j)=moments(j)+(x**j*y**j)*array(kx,ky)
        end do
     end do
     moments(j)=moments(j)*h*h
  end do

end subroutine moments_2d
