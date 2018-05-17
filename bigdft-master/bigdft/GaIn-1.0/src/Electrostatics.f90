!> @file Electrostatics.f90
!!
!! Defines utility routines for computing coulomb potential and field
!! issued from gaussian basis elements.
!!
!! Author: I. Duchemin July 2015
!!

!> compute the potential generated at position \f$r_{test}\f$ by a Cubic Harmonic Orbital centered in r.
!!
!!  \f$
!!      \int dr' \, \frac{1}{|r_{test}-r'|} Y_{xyz}(r'-r)
!!  \f$
!! 
recursive function potential_from_C(r_test,alpha,r,nx,ny,nz) result(value)
  
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  
  implicit none
  
  ! input variables
  real(kind=8) :: r_test(3)  !< position where we compute the potential
  real(kind=8) :: alpha      !< exponent of the source charge distribution
  real(kind=8) :: r(3)       !< center of the source charge distribution
  integer      :: nx         !< x power for first cubic Harmonic
  integer      :: ny         !< y power for first cubic Harmonic
  integer      :: nz         !< z power for first cubic Harmonic
  
  ! return value
  real(kind=8) :: value
  
  ! parameters
  integer , parameter     :: l_test       =0
  integer , parameter     :: m_test       =0
  real(kind=8), parameter :: alpha_test   =1.0d16
  real(kind=8), parameter :: charge_norm  =2.0d0/(4.0d0*atan(1.0d0))*(alpha_test*sqrt(alpha_test))
  integer, parameter      :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! local variable
  integer       :: ir2
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: R_X_R
  
  ! decomposition of the S solid harmonic into cubic harmonic
  call R_from_Y(0,0,c1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx,ny,nz)
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx+ny+nz))=0.0d0
  c2(ir2)=1.0d0
  
  ! set value
  value = charge_norm * R_X_R(alpha_test,r_test,c1,0,alpha,r,c2,nx+ny+nz)
  
end function

!> compute the potential generated at position r_test by a Solid Harmonic Orbital centered in r.
!!
!!  \f$
!!      \int dr' \, \frac{1}{|r_{test}-r'|} Y_{lm}(r'-r)
!!  \f$
!!
recursive function potential_from_Y(r_test,alpha,r,l,m) result(value)
  
  implicit none
  
  ! input variables
  real(kind=8) :: r_test(3)  !< position where we compute the potential
  real(kind=8) :: alpha      !< exponent of the source charge distribution
  real(kind=8) :: r(3)       !< center of the source charge distribution
  integer      :: l          !< l momentum of source the charge distribution
  integer      :: m          !< m momentum of source the charge distribution
  
  ! return value
  real(kind=8) :: value
  
  ! parameters
  integer , parameter     :: l_test       =0
  integer , parameter     :: m_test       =0
  real(kind=8), parameter :: alpha_test   =1.0d16
  real(kind=8), parameter :: charge_norm  =2.0d0/(4.0d0*atan(1.0d0))*(alpha_test*sqrt(alpha_test))
  
  ! coulomb interaction function
  real(kind=8) :: Y_Coulomb_Y
  
  ! set value
  value = charge_norm * Y_Coulomb_Y(alpha_test, r_test, l_test, m_test, alpha, r, l, m)
  
end function
  
!> compute the three component of the electric field generated at position r_test by a Solid Harmonic Orbital centered in r
!!
!!  \f$
!!     \nabla \int dr' \, \frac{1}{|r_{test}-r'|} Y_{lm}(r'-r)
!!  \f$
!!
recursive subroutine field_from_Y(r_test,alpha,r,l,m,E_test)
  
  implicit none
  
  ! input variables
  real(kind=8) :: r_test(3)  !< position where we compute the potential
  real(kind=8) :: alpha      !< exponent of the source charge distribution
  real(kind=8) :: r(3)       !< center of the source charge distribution
  integer      :: l          !< l momentum of source the charge distribution
  integer      :: m          !< m momentum of source the charge distribution
  
  ! return value
  real(kind=8) :: E_test(3)  !< electric field at test position
  
  ! parameters
  real(kind=8), parameter :: alpha_test   =1.0d16
  real(kind=8), parameter :: charge_norm  =-4.0d0/(4.0d0*atan(1.0d0))/sqrt(3.0d0)*(alpha_test**2*sqrt(alpha_test))
  
  ! coulomb interaction function
  real(kind=8) :: Y_Coulomb_Y
  
  ! set value
  E_test(1) = charge_norm * Y_Coulomb_Y(alpha_test, r_test, 1, 1,alpha,r,l,m)
  E_test(2) = charge_norm * Y_Coulomb_Y(alpha_test, r_test, 1,-1,alpha,r,l,m)
  E_test(3) = charge_norm * Y_Coulomb_Y(alpha_test, r_test, 1, 0,alpha,r,l,m)
  
end subroutine  
  
!> compute the three component of the electric field generated at position r_test by a Cubic Harmonic Orbital centered in r
!!
!!  \f$
!!     \nabla \int dr' \, \frac{1}{|r_{test}-r'|} Y_{xyz}(r'-r)
!!  \f$
!!
recursive subroutine field_from_C(r_test,alpha,r,nx,ny,nz,E_test)
  
  use mod_R_from_Y
  use mod_XZY_power_to_ir
  
  implicit none
  
  ! input variables
  real(kind=8) :: r_test(3)  !< position where we compute the potential
  real(kind=8) :: alpha      !< exponent of the source charge distribution
  real(kind=8) :: r(3)       !< center of the source charge distribution
  integer      :: nx         !< x power for first cubic Harmonic
  integer      :: ny         !< y power for first cubic Harmonic
  integer      :: nz         !< z power for first cubic Harmonic
  
  ! return value
  real(kind=8) :: E_test(3)  !< electric field at test position
  
  ! parameters
  real(kind=8), parameter :: alpha_test   =1.0d16
  real(kind=8), parameter :: charge_norm  =-4.0d0/(4.0d0*atan(1.0d0))/sqrt(3.0d0)*(alpha_test**2*sqrt(alpha_test))
  integer, parameter      :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! local variable
  integer       :: ir2
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: R_X_R
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx,ny,nz)
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx+ny+nz))=0.0d0
  c2(ir2)=1.0d0
  
  ! decomposition of the Px solid harmonic into cubic harmonic
  call R_from_Y(1,1,c1)
  
  ! set value
  E_test(1) = charge_norm * R_X_R(alpha_test,r_test,c1,1,alpha,r,c2,nx+ny+nz)
  
  
  ! decomposition of the Py solid harmonic into cubic harmonic
  call R_from_Y(1,-1,c1)
  
  ! set value
  E_test(2) = charge_norm * R_X_R(alpha_test,r_test,c1,1,alpha,r,c2,nx+ny+nz)
  
  ! decomposition of the Pz solid harmonic into cubic harmonic
  call R_from_Y(1, 0,c1)
  
  ! set value
  E_test(3) = charge_norm * R_X_R(alpha_test,r_test,c1,1,alpha,r,c2,nx+ny+nz)
  
end subroutine


