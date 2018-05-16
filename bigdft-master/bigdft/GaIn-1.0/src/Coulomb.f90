
!> @file Coulomb.f90
!!
!! Defines utility routines for computing coulomb 
!! integrals between gaussian basis elements.
!!
!! Author: I. Duchemin July 2015
!!
recursive function R_X_R(alpha1,r1,c1,l1max,alpha2,r2,c2,l2max)
  
  use mod_CoulombUtils
  use mod_R_to_D_Conversion
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)   !< center for first cubic Harmonic
  real(kind=8) :: r2(3)   !< center for second cubic Harmonic
  real(kind=8) :: alpha1  !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2  !< exponent for second cubic Harmonic
  integer      :: l1max   !< max order for first cubic Harmonic
  integer      :: l2max   !< max order for second cubic Harmonic
  real(kind=8) :: c1(455) !< coefficients of first cubic Harmonic in the R basis
  real(kind=8) :: c2(455) !< coefficients of second cubic Harmonic in the R basis
  
  ! return value
  real(kind=8)  :: R_X_R
  
  ! local variables
  real(kind=8) :: coeffs1(455)
  real(kind=8) :: coeffs2(455)
  real(kind=8) :: coeffs_tmp(455)
  real(kind=8) :: x,y,z,q
  real(kind=8) :: p,d,erfpd,exppd
  integer      :: il,ir
  
  ! parameters
  real(kind=8), parameter :: PI=3.1415926535897932384626433832795d0
  integer     , parameter :: imax(0:12)=(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! form coeffs for derivatives on the left side
  coeffs1=0.0d0
  do il=1,imax(l1max)
    if ( c1(il).ne.0.0d0 ) then
      ! compute contribution from this index
      call R_to_D(alpha1,il,coeffs_tmp)
      ! add contrib
      coeffs1(1:il)=coeffs1(1:il)+c1(il)*coeffs_tmp(1:il)
    end if
  end do
  
  ! form coeffs for derivatives on the right side
  coeffs2=0.0d0
  do ir=1,imax(l2max)
    if ( c2(ir).ne.0.0d0 ) then
      ! compute contribution from this index
      call R_to_D(alpha2,ir,coeffs_tmp)
      ! add contrib
      coeffs2(1:ir)=coeffs2(1:ir)+c2(ir)*coeffs_tmp(1:ir)
    end if
  end do
  
  ! compute combined exponant
  q=alpha1*alpha2/(alpha1+alpha2)
  
  ! compute center distance
  x=r2(1)-r1(1)
  y=r2(2)-r1(2)
  z=r2(3)-r1(3)
  
  ! compute once and for all p=sqrt(q), d, erf(p*d) and exp(-p^2*d^2)
  p=sqrt(q)
  d=sqrt(x**2+y**2+z**2)
  erfpd=erf(p*d)
  exppd=exp(-p**2*d**2)
  
  ! loop on the components
  R_X_R=0.0d0
  do il=1,imax(l1max)
    if ( coeffs1(il).ne.0.0d0 ) then
      do ir=1,imax(l2max)
        if ( coeffs2(ir).ne.0.0d0 ) then
          R_X_R =R_X_R + coeffs1(il)*coeffs2(ir)*D_X_D(il,ir,x,y,z,q,p,d,erfpd,exppd)
        end if
      end do
    end if
  end do
  
  ! norm and return value
  R_X_R =R_X_R*sqrt(PI/alpha1)*(PI/alpha1)*sqrt(PI/alpha2)*(PI/alpha2)
  
end function
  
!> Two centers Repulsion integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{xyz}^{(1)}(r-R_1) \frac{1}{|r-r'|} Y_{xyz}^{(2)}(r'-R_2)
!!  \f$
!!
recursive function C_Coulomb_C(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2)
  
  use mod_XZY_power_to_ir
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2   !< exponent for second cubic Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  
  ! return value
  real(kind=8)  :: C_Coulomb_C
  
  ! local variable
  integer       :: ir1
  integer       :: ir2
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! compute coulomb integral
  C_Coulomb_C =R_X_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2)
  
end function
  
  
!> Two centers Repulsion integrals for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{lm}^{(1)}(r-R_1) \frac{1}{|r-r'|} Y_{lm}^{(2)}(r'-R_2)
!!  \f$
!!
recursive function Y_Coulomb_Y(alpha1,r1,l1,m1,alpha2,r2,l2,m2)
  
  use mod_R_from_Y
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  
  ! return value
  real(kind=8)  :: Y_Coulomb_Y
  
  ! local variable
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: R_X_R
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! compute coulomb integral
  Y_Coulomb_Y =R_X_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2)
  
end function
  
!> Two centers ionic integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, Y_{xyz}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) \frac{1}{|r-R_{ion}|} 
!!  \f$
!!
recursive function CC_Coulomb_Ion(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2,rion)
  
  use mod_R_1_norm
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: rion(3)  !< ion position
  real(kind=8) :: alpha1   !< exponent for first cubic Harmonic
  real(kind=8) :: alpha2   !< exponent for second cubic Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  
  ! return value
  real(kind=8)  :: CC_Coulomb_Ion
  
  ! local variable
  integer       :: ir1
  integer       :: ir2
  integer       :: l3
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: alpha3
  real(kind=8)  :: r3(3)
  real(kind=8)  :: R_X_R
  
  ! parameters
  real(kind=8)       :: alpha_ion=1.0e16
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2,alpha3,r3,c3,l3)
  
  ! set ion coeff in the R basis
  c4(1)=1.0d0/R_1_norm(alpha_ion,(/1.0d0/),0)
  
  ! compute coulomb integral
  CC_Coulomb_Ion =R_X_R(alpha3,r3,c3,l3,alpha_ion,rion,c4,0)
  
end function

!> Two centers ionic integrals for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, Y_{lm}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) \frac{1}{|r-R_{ion}|} 
!!  \f$
!!
recursive function YY_Coulomb_Ion(alpha1,r1,l1,m1,alpha2,r2,l2,m2,rion)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: rion(3)  !< ion position
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  
  ! return value
  real(kind=8)  :: YY_Coulomb_Ion
  
  ! local variable
  integer       :: l3
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: alpha3
  real(kind=8)  :: r3(3)
  real(kind=8)  :: R_X_R
  
  ! parameters
  real(kind=8)  :: alpha_ion=1.0e16
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,alpha3,r3,c3,l3)
  
  ! set ion coeff in the R basis
  c4(1)=1.0d0/R_1_norm(alpha_ion,(/1.0d0/),0)
  
  ! compute coulomb integral
  YY_Coulomb_Ion =R_X_R(alpha3,r3,c3,l3,alpha_ion,rion,c4,0)
  
end function


!> Three centers Repulsion integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{xyz}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) \frac{1}{|r-r'|} Y_{xyz}^{(3)}(r'-R_3)
!!  \f$
!!
recursive function CC_Coulomb_C(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2,alpha3,r3,nx3,ny3,nz3)
  
  use mod_R_1_norm
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: r3(3)    !< center for third cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  integer      :: nx3      !< x power for third cubic Harmonic
  integer      :: ny3      !< y power for third cubic Harmonic
  integer      :: nz3      !< z power for third cubic Harmonic
  
  ! return value
  real(kind=8)  :: CC_Coulomb_C
  
  ! local variable
  integer       :: ir1
  integer       :: ir2
  integer       :: ir3
  integer       :: l4
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: alpha4
  real(kind=8)  :: r4(3)
  real(kind=8)  :: R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! get 3rd member index in the R basis
  ir3 =ir_index(nx3,ny3,nz3)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! set 3rd member coeff in the R basis
  c3(1:imax(nx3+ny3+nz3))=0.0d0
  c3(ir3)=1.0d0
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2,alpha4,r4,c4,l4)
  
  ! compute coulomb integral
  CC_Coulomb_C =R_X_R(alpha3,r3,c3,nx3+ny3+nz3,alpha4,r4,c4,l4)
  
end function


!> Three centers Repulsion integrals for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{lm}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) \frac{1}{|r-r'|} Y_{lm}^{(3)}(r'-R_3)
!!  \f$
!!
recursive function YY_Coulomb_Y(alpha1,r1,l1,m1,alpha2,r2,l2,m2,alpha3,r3,l3,m3)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: r3(3)    !< center for third spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: l3       !< angular momentum for third spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  integer      :: m3       !< orbital momentum for third spherical Harmonic
  
  ! return value
  real(kind=8)  :: YY_Coulomb_Y
  
  ! local variable
  integer       :: l4
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: alpha4
  real(kind=8)  :: r4(3)
  real(kind=8)  :: R_X_R
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l3,m3,c3)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,alpha4,r4,c4,l4)
  
  ! compute coulomb integral
  YY_Coulomb_Y =R_X_R(alpha3,r3,c3,l3,alpha4,r4,c4,l4)
  
end function


!> Four centers Repulsion integrals for Cubic Spherical Harmonics (C)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{xyz}^{(1)}(r-R_1) Y_{xyz}^{(2)}(r-R_2) \frac{1}{|r-r'|} Y_{xyz}^{(3)}(r'-R_3) Y_{xyz}^{(4)}(r'-R_4)
!!  \f$
!!
recursive function CC_Coulomb_CC(alpha1,r1,nx1,ny1,nz1,alpha2,r2,nx2,ny2,nz2,alpha3,r3,nx3,ny3,nz3,alpha4,r4,nx4,ny4,nz4)
  
  use mod_R_1_norm
  use mod_XZY_power_to_ir
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first cubic Harmonic
  real(kind=8) :: r2(3)    !< center for second cubic Harmonic
  real(kind=8) :: r3(3)    !< center for third cubic Harmonic
  real(kind=8) :: r4(3)    !< center for fourth cubic Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
  real(kind=8) :: alpha4   !< exponent for fourth spherical Harmonic
  integer      :: nx1      !< x power for first cubic Harmonic
  integer      :: ny1      !< y power for first cubic Harmonic
  integer      :: nz1      !< z power for first cubic Harmonic
  integer      :: nx2      !< x power for second cubic Harmonic
  integer      :: ny2      !< y power for second cubic Harmonic
  integer      :: nz2      !< z power for second cubic Harmonic
  integer      :: nx3      !< x power for third cubic Harmonic
  integer      :: ny3      !< y power for third cubic Harmonic
  integer      :: nz3      !< z power for third cubic Harmonic
  integer      :: nx4      !< x power for fourth cubic Harmonic
  integer      :: ny4      !< y power for fourth cubic Harmonic
  integer      :: nz4      !< z power for fourth cubic Harmonic
  
  ! return value
  real(kind=8)  :: CC_Coulomb_CC
  
  ! local variable
  integer       :: ir1
  integer       :: ir2
  integer       :: ir3
  integer       :: ir4
  integer       :: l5
  integer       :: l6
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: c5(455)
  real(kind=8)  :: c6(455)
  real(kind=8)  :: alpha5
  real(kind=8)  :: alpha6
  real(kind=8)  :: r5(3)
  real(kind=8)  :: r6(3)
  real(kind=8)  :: R_X_R
  
  ! parameters
  integer, parameter :: imax(0:12) =(/1,4,10,20,35,56,84,120,165,220,286,364,455/) !< maximum decomposition coeff index for a given l
  
  ! get 1st member index in the R basis
  ir1 =ir_index(nx1,ny1,nz1)
  
  ! get 2nd member index in the R basis
  ir2 =ir_index(nx2,ny2,nz2)
  
  ! get 3rd member index in the R basis
  ir3 =ir_index(nx3,ny3,nz3)
  
  ! get 4th member index in the R basis
  ir4 =ir_index(nx4,ny4,nz4)
  
  ! set 1st member coeff in the R basis
  c1(1:imax(nx1+ny1+nz1))=0.0d0
  c1(ir1)=1.0d0
  
  ! set 2nd member coeff in the R basis
  c2(1:imax(nx2+ny2+nz2))=0.0d0
  c2(ir2)=1.0d0
  
  ! set 3rd member coeff in the R basis
  c3(1:imax(nx3+ny3+nz3))=0.0d0
  c3(ir3)=1.0d0
  
  ! set 4th member coeff in the R basis
  c4(1:imax(nx4+ny4+nz4))=0.0d0
  c4(ir4)=1.0d0
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,nx1+ny1+nz1,alpha2,r2,c2,nx2+ny2+nz2,alpha5,r5,c5,l5)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha3,r3,c3,nx3+ny3+nz3,alpha4,r4,c4,nx4+ny4+nz4,alpha6,r6,c6,l6)
  
  ! compute coulomb integral
  CC_Coulomb_CC =R_X_R(alpha5,r5,c5,l5,alpha6,r6,c6,l6)
  
end function

!> Four centers Repulsion integrals for Solid Spherical Harmonics (Y)
!!
!!  \f$
!!      \int dr \, dr' \, Y_{lm}^{(1)}(r-R_1) Y_{lm}^{(2)}(r-R_2) \frac{1}{|r-r'|} Y_{lm}^{(3)}(r'-R_3) Y_{lm}^{(4)}(r'-R_4)
!!  \f$
!!
recursive function YY_Coulomb_YY(alpha1,r1,l1,m1,alpha2,r2,l2,m2,alpha3,r3,l3,m3,alpha4,r4,l4,m4)
  
  use mod_R_1_norm
  use mod_R_from_Y
  use mod_CubicHarmonicsProduct
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: r3(3)    !< center for third spherical Harmonic
  real(kind=8) :: r4(3)    !< center for third spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
  real(kind=8) :: alpha4   !< exponent for third spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: l3       !< angular momentum for third spherical Harmonic
  integer      :: l4       !< angular momentum for third spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  integer      :: m3       !< orbital momentum for third spherical Harmonic
  integer      :: m4       !< orbital momentum for third spherical Harmonic
  
  ! return value
  real(kind=8)  :: YY_Coulomb_YY
  
  ! local variable
  integer       :: l5
  integer       :: l6
  real(kind=8)  :: c1(455)
  real(kind=8)  :: c2(455)
  real(kind=8)  :: c3(455)
  real(kind=8)  :: c4(455)
  real(kind=8)  :: c5(455)
  real(kind=8)  :: c6(455)
  real(kind=8)  :: alpha5
  real(kind=8)  :: alpha6
  real(kind=8)  :: r5(3)
  real(kind=8)  :: r6(3)
  real(kind=8)  :: R_X_R
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l3,m3,c3)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l4,m4,c4)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,alpha5,r5,c5,l5)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha3,r3,c3,l3,alpha4,r4,c4,l4,alpha6,r6,c6,l6)
  
  ! compute coulomb integral
  YY_Coulomb_YY =R_X_R(alpha5,r5,c5,l5,alpha6,r6,c6,l6)
  
end function


!> Two centers Modified Repulsion integrals for Solid Spherical Harmonics (Y)
!!
!!         /
!! compute | dr1 dr2 Ylm1(r1-R1) * Y00(r1-R2) * 1/|r1-r2| * Ylm2(r2-R2)
!!         /
!!
recursive function Y_ModCoulomb_Y(acut,alpha1,r1,l1,m1,alpha2,r2,l2,m2)
  
  use mod_CubicHarmonicsProduct
  use mod_R_from_Y
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: acut     !< exponent for cutoff spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  
  ! return value
  real(kind=8)  :: Y_ModCoulomb_Y
  
  ! local variable
  real(kind=8)  :: c1(455),c2(455),c3(455)
  real(kind=8)  :: ac
  integer       :: lmax
  real(kind=8)  :: rc(3)
  real(kind=8)  :: cc(455)
  real(kind=8)  :: R_X_R
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! decomposition of the cutoff solid harmonic into cubic harmonic
  c3   =0.0d0
  c3(1)=1.0d0
  
  ! get product of first solid harmonic times cutoff solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,acut,r2,c3,0,ac,rc,cc,lmax)
  
  ! compute coulomb integral
  Y_ModCoulomb_Y =R_X_R(ac,rc,cc,lmax,alpha2,r2,c2,l2)
  
end function


!> Three centers Electron Modified Repulsion integrals for Solid Spherical Harmonics (Y)
!!
!!         /
!! compute | dr1 dr2 Ylm1(r1-R1) * Ylm2(r1-R2) * Y00(r1-R3) * 1/|r1-r2| * Ylm3(r2-R3)
!!         /
!!
recursive function YY_ModCoulomb_Y(acut,alpha1,r1,l1,m1,alpha2,r2,l2,m2,alpha3,r3,l3,m3)
  
  use mod_CubicHarmonicsProduct
  use mod_R_from_Y
  
  implicit none
  
  ! inputs arguments
  real(kind=8) :: r1(3)    !< center for first spherical Harmonic
  real(kind=8) :: r2(3)    !< center for second spherical Harmonic
  real(kind=8) :: r3(3)    !< center for third spherical Harmonic
  real(kind=8) :: acut     !< exponent for cutoff spherical Harmonic
  real(kind=8) :: alpha1   !< exponent for first spherical Harmonic
  real(kind=8) :: alpha2   !< exponent for second spherical Harmonic
  real(kind=8) :: alpha3   !< exponent for third spherical Harmonic
  integer      :: l1       !< angular momentum for first spherical Harmonic
  integer      :: l2       !< angular momentum for second spherical Harmonic
  integer      :: l3       !< angular momentum for third spherical Harmonic
  integer      :: m1       !< orbital momentum for first spherical Harmonic
  integer      :: m2       !< orbital momentum for second spherical Harmonic
  integer      :: m3       !< orbital momentum for third spherical Harmonic
  
  ! return value
  real(kind=8)  :: YY_ModCoulomb_Y
  
  ! local variable
  real(kind=8)  :: c1(455),c2(455),c3(455)
  real(kind=8)  :: ac,atmp
  integer       :: lmax,ltmp
  real(kind=8)  :: rc(3),rtmp(3)
  real(kind=8)  :: R_X_R
  
  ! decomposition of the first solid harmonic into cubic harmonic
  call R_from_Y(l1,m1,c1)
  
  ! decomposition of the second solid harmonic into cubic harmonic
  call R_from_Y(l2,m2,c2)
  
  ! get product of first and second solid harmonic
  call RxR_to_R(alpha1,r1,c1,l1,alpha2,r2,c2,l2,atmp,rtmp,c3,ltmp)
  
  ! get product with cutoff solid harmonic
  c2(1)=1.0d0
  call RxR_to_R(atmp,rtmp,c3,ltmp,acut,r3,c2, 0,ac,rc,c1,lmax)
  
  ! decomposition of the third solid harmonic into cubic harmonic
  call R_from_Y(l3,m3,c3)
  
  ! compute coulomb integral
  YY_ModCoulomb_Y =R_X_R(ac,rc,c1,lmax,alpha3,r3,c3,l3)
  
end function

