
!> This module defines utility routines for computing the 
!! laplacian of x^n*y^p*z^q*exp(-a*(x^2+y^2+z^2)) functions 
!!
!! These routines are generated automatically through 
!! maxima scripting.
!!
!! Author: I. Duchemin Aout 2015
!!
module mod_R_Laplacian

  implicit none

contains

!> expand laplacian of x^n*y^p*z^q*exp(-a*(x^2+y^2+z^2))
!!
recursive subroutine R_Laplacian(a,c,lmax,t)
 
  implicit none
 
  real(kind=8), intent(in):: a!< exponant of the solid harmonic
  integer     , intent(in):: lmax!< maximum l momentum
  real(kind=8), intent(in) ,dimension(*):: c!< coefficients of the decomposition of the cubic harmonic in the R basis
  real(kind=8), intent(out),dimension(*):: t!< coefficients of the decomposition of the laplacian of the cubic harmonic in the R basis
  
  
  ! select lmax case
  select case(lmax)
    case(0)
      t(1)=0.0d0+c(1)*(-6.0d0*a)
      t(2)=0.0d0
      t(3)=0.0d0
      t(4)=0.0d0
      t(5)=0.0d0+c(1)*(4.0d0*a**2)
      t(6)=0.0d0
      t(7)=0.0d0
      t(8)=0.0d0+c(1)*(4.0d0*a**2)
      t(9)=0.0d0
      t(10)=0.0d0+c(1)*(4.0d0*a**2)
    case(1)
      t(1)=0.0d0+c(1)*(-6.0d0*a)
      t(2)=0.0d0+c(2)*(-10.0d0*a)
      t(3)=0.0d0+c(3)*(-10.0d0*a)
      t(4)=0.0d0+c(4)*(-10.0d0*a)
      t(5)=0.0d0+c(1)*(4.0d0*a**2)
      t(6)=0.0d0
      t(7)=0.0d0
      t(8)=0.0d0+c(1)*(4.0d0*a**2)
      t(9)=0.0d0
      t(10)=0.0d0+c(1)*(4.0d0*a**2)
      t(11)=0.0d0+c(2)*(4.0d0*a**2)
      t(12)=0.0d0+c(3)*(4.0d0*a**2)
      t(13)=0.0d0+c(4)*(4.0d0*a**2)
      t(14)=0.0d0+c(2)*(4.0d0*a**2)
      t(15)=0.0d0
      t(16)=0.0d0+c(2)*(4.0d0*a**2)
      t(17)=0.0d0+c(3)*(4.0d0*a**2)
      t(18)=0.0d0+c(4)*(4.0d0*a**2)
      t(19)=0.0d0+c(3)*(4.0d0*a**2)
      t(20)=0.0d0+c(4)*(4.0d0*a**2)
    case(2)
      t(1)=0.0d0+c(1)*(-6.0d0*a)+c(5)*(2.0d0)+c(8)*(2.0d0)+c(10)*(2.0d0)
      t(2)=0.0d0+c(2)*(-10.0d0*a)
      t(3)=0.0d0+c(3)*(-10.0d0*a)
      t(4)=0.0d0+c(4)*(-10.0d0*a)
      t(5)=0.0d0+c(1)*(4.0d0*a**2)+c(5)*(-14.0d0*a)
      t(6)=0.0d0+c(6)*(-14.0d0*a)
      t(7)=0.0d0+c(7)*(-14.0d0*a)
      t(8)=0.0d0+c(1)*(4.0d0*a**2)+c(8)*(-14.0d0*a)
      t(9)=0.0d0+c(9)*(-14.0d0*a)
      t(10)=0.0d0+c(1)*(4.0d0*a**2)+c(10)*(-14.0d0*a)
      t(11)=0.0d0+c(2)*(4.0d0*a**2)
      t(12)=0.0d0+c(3)*(4.0d0*a**2)
      t(13)=0.0d0+c(4)*(4.0d0*a**2)
      t(14)=0.0d0+c(2)*(4.0d0*a**2)
      t(15)=0.0d0
      t(16)=0.0d0+c(2)*(4.0d0*a**2)
      t(17)=0.0d0+c(3)*(4.0d0*a**2)
      t(18)=0.0d0+c(4)*(4.0d0*a**2)
      t(19)=0.0d0+c(3)*(4.0d0*a**2)
      t(20)=0.0d0+c(4)*(4.0d0*a**2)
      t(21)=0.0d0+c(5)*(4.0d0*a**2)
      t(22)=0.0d0+c(6)*(4.0d0*a**2)
      t(23)=0.0d0+c(7)*(4.0d0*a**2)
      t(24)=0.0d0+c(5)*(4.0d0*a**2)+c(8)*(4.0d0*a**2)
      t(25)=0.0d0+c(9)*(4.0d0*a**2)
      t(26)=0.0d0+c(5)*(4.0d0*a**2)+c(10)*(4.0d0*a**2)
      t(27)=0.0d0+c(6)*(4.0d0*a**2)
      t(28)=0.0d0+c(7)*(4.0d0*a**2)
      t(29)=0.0d0+c(6)*(4.0d0*a**2)
      t(30)=0.0d0+c(7)*(4.0d0*a**2)
      t(31)=0.0d0+c(8)*(4.0d0*a**2)
      t(32)=0.0d0+c(9)*(4.0d0*a**2)
      t(33)=0.0d0+c(8)*(4.0d0*a**2)+c(10)*(4.0d0*a**2)
      t(34)=0.0d0+c(9)*(4.0d0*a**2)
      t(35)=0.0d0+c(10)*(4.0d0*a**2)
    case(3)
      t(1)=0.0d0+c(1)*(-6.0d0*a)+c(5)*(2.0d0)+c(8)*(2.0d0)+c(10)*(2.0d0)
      t(2)=0.0d0+c(2)*(-10.0d0*a)+c(11)*(6.0d0)+c(14)*(2.0d0)+c(16)*(2.0d0)
      t(3)=0.0d0+c(3)*(-10.0d0*a)+c(12)*(2.0d0)+c(17)*(6.0d0)+c(19)*(2.0d0)
      t(4)=0.0d0+c(4)*(-10.0d0*a)+c(13)*(2.0d0)+c(18)*(2.0d0)+c(20)*(6.0d0)
      t(5)=0.0d0+c(1)*(4.0d0*a**2)+c(5)*(-14.0d0*a)
      t(6)=0.0d0+c(6)*(-14.0d0*a)
      t(7)=0.0d0+c(7)*(-14.0d0*a)
      t(8)=0.0d0+c(1)*(4.0d0*a**2)+c(8)*(-14.0d0*a)
      t(9)=0.0d0+c(9)*(-14.0d0*a)
      t(10)=0.0d0+c(1)*(4.0d0*a**2)+c(10)*(-14.0d0*a)
      t(11)=0.0d0+c(2)*(4.0d0*a**2)+c(11)*(-18.0d0*a)
      t(12)=0.0d0+c(3)*(4.0d0*a**2)+c(12)*(-18.0d0*a)
      t(13)=0.0d0+c(4)*(4.0d0*a**2)+c(13)*(-18.0d0*a)
      t(14)=0.0d0+c(2)*(4.0d0*a**2)+c(14)*(-18.0d0*a)
      t(15)=0.0d0+c(15)*(-18.0d0*a)
      t(16)=0.0d0+c(2)*(4.0d0*a**2)+c(16)*(-18.0d0*a)
      t(17)=0.0d0+c(3)*(4.0d0*a**2)+c(17)*(-18.0d0*a)
      t(18)=0.0d0+c(4)*(4.0d0*a**2)+c(18)*(-18.0d0*a)
      t(19)=0.0d0+c(3)*(4.0d0*a**2)+c(19)*(-18.0d0*a)
      t(20)=0.0d0+c(4)*(4.0d0*a**2)+c(20)*(-18.0d0*a)
      t(21)=0.0d0+c(5)*(4.0d0*a**2)
      t(22)=0.0d0+c(6)*(4.0d0*a**2)
      t(23)=0.0d0+c(7)*(4.0d0*a**2)
      t(24)=0.0d0+c(5)*(4.0d0*a**2)+c(8)*(4.0d0*a**2)
      t(25)=0.0d0+c(9)*(4.0d0*a**2)
      t(26)=0.0d0+c(5)*(4.0d0*a**2)+c(10)*(4.0d0*a**2)
      t(27)=0.0d0+c(6)*(4.0d0*a**2)
      t(28)=0.0d0+c(7)*(4.0d0*a**2)
      t(29)=0.0d0+c(6)*(4.0d0*a**2)
      t(30)=0.0d0+c(7)*(4.0d0*a**2)
      t(31)=0.0d0+c(8)*(4.0d0*a**2)
      t(32)=0.0d0+c(9)*(4.0d0*a**2)
      t(33)=0.0d0+c(8)*(4.0d0*a**2)+c(10)*(4.0d0*a**2)
      t(34)=0.0d0+c(9)*(4.0d0*a**2)
      t(35)=0.0d0+c(10)*(4.0d0*a**2)
      t(36)=0.0d0+c(11)*(4.0d0*a**2)
      t(37)=0.0d0+c(12)*(4.0d0*a**2)
      t(38)=0.0d0+c(13)*(4.0d0*a**2)
      t(39)=0.0d0+c(11)*(4.0d0*a**2)+c(14)*(4.0d0*a**2)
      t(40)=0.0d0+c(15)*(4.0d0*a**2)
      t(41)=0.0d0+c(11)*(4.0d0*a**2)+c(16)*(4.0d0*a**2)
      t(42)=0.0d0+c(12)*(4.0d0*a**2)+c(17)*(4.0d0*a**2)
      t(43)=0.0d0+c(13)*(4.0d0*a**2)+c(18)*(4.0d0*a**2)
      t(44)=0.0d0+c(12)*(4.0d0*a**2)+c(19)*(4.0d0*a**2)
      t(45)=0.0d0+c(13)*(4.0d0*a**2)+c(20)*(4.0d0*a**2)
      t(46)=0.0d0+c(14)*(4.0d0*a**2)
      t(47)=0.0d0+c(15)*(4.0d0*a**2)
      t(48)=0.0d0+c(14)*(4.0d0*a**2)+c(16)*(4.0d0*a**2)
      t(49)=0.0d0+c(15)*(4.0d0*a**2)
      t(50)=0.0d0+c(16)*(4.0d0*a**2)
      t(51)=0.0d0+c(17)*(4.0d0*a**2)
      t(52)=0.0d0+c(18)*(4.0d0*a**2)
      t(53)=0.0d0+c(17)*(4.0d0*a**2)+c(19)*(4.0d0*a**2)
      t(54)=0.0d0+c(18)*(4.0d0*a**2)+c(20)*(4.0d0*a**2)
      t(55)=0.0d0+c(19)*(4.0d0*a**2)
      t(56)=0.0d0+c(20)*(4.0d0*a**2)
    case(4)
      t(1)=0.0d0+c(1)*(-6.0d0*a)+c(5)*(2.0d0)+c(8)*(2.0d0)+c(10)*(2.0d0)
      t(2)=0.0d0+c(2)*(-10.0d0*a)+c(11)*(6.0d0)+c(14)*(2.0d0)+c(16)*(2.0d0)
      t(3)=0.0d0+c(3)*(-10.0d0*a)+c(12)*(2.0d0)+c(17)*(6.0d0)+c(19)*(2.0d0)
      t(4)=0.0d0+c(4)*(-10.0d0*a)+c(13)*(2.0d0)+c(18)*(2.0d0)+c(20)*(6.0d0)
      t(5)=0.0d0+c(1)*(4.0d0*a**2)+c(5)*(-14.0d0*a)+c(21)*(12.0d0)+c(24)*(2.0d0)+c(26)*(2.0d0)
      t(6)=0.0d0+c(6)*(-14.0d0*a)+c(22)*(6.0d0)+c(27)*(6.0d0)+c(29)*(2.0d0)
      t(7)=0.0d0+c(7)*(-14.0d0*a)+c(23)*(6.0d0)+c(28)*(2.0d0)+c(30)*(6.0d0)
      t(8)=0.0d0+c(1)*(4.0d0*a**2)+c(8)*(-14.0d0*a)+c(24)*(2.0d0)+c(31)*(12.0d0)+c(33)*(2.0d0)
      t(9)=0.0d0+c(9)*(-14.0d0*a)+c(25)*(2.0d0)+c(32)*(6.0d0)+c(34)*(6.0d0)
      t(10)=0.0d0+c(1)*(4.0d0*a**2)+c(10)*(-14.0d0*a)+c(26)*(2.0d0)+c(33)*(2.0d0)+c(35)*(12.0d0)
      t(11)=0.0d0+c(2)*(4.0d0*a**2)+c(11)*(-18.0d0*a)
      t(12)=0.0d0+c(3)*(4.0d0*a**2)+c(12)*(-18.0d0*a)
      t(13)=0.0d0+c(4)*(4.0d0*a**2)+c(13)*(-18.0d0*a)
      t(14)=0.0d0+c(2)*(4.0d0*a**2)+c(14)*(-18.0d0*a)
      t(15)=0.0d0+c(15)*(-18.0d0*a)
      t(16)=0.0d0+c(2)*(4.0d0*a**2)+c(16)*(-18.0d0*a)
      t(17)=0.0d0+c(3)*(4.0d0*a**2)+c(17)*(-18.0d0*a)
      t(18)=0.0d0+c(4)*(4.0d0*a**2)+c(18)*(-18.0d0*a)
      t(19)=0.0d0+c(3)*(4.0d0*a**2)+c(19)*(-18.0d0*a)
      t(20)=0.0d0+c(4)*(4.0d0*a**2)+c(20)*(-18.0d0*a)
      t(21)=0.0d0+c(5)*(4.0d0*a**2)+c(21)*(-22.0d0*a)
      t(22)=0.0d0+c(6)*(4.0d0*a**2)+c(22)*(-22.0d0*a)
      t(23)=0.0d0+c(7)*(4.0d0*a**2)+c(23)*(-22.0d0*a)
      t(24)=0.0d0+c(5)*(4.0d0*a**2)+c(8)*(4.0d0*a**2)+c(24)*(-22.0d0*a)
      t(25)=0.0d0+c(9)*(4.0d0*a**2)+c(25)*(-22.0d0*a)
      t(26)=0.0d0+c(5)*(4.0d0*a**2)+c(10)*(4.0d0*a**2)+c(26)*(-22.0d0*a)
      t(27)=0.0d0+c(6)*(4.0d0*a**2)+c(27)*(-22.0d0*a)
      t(28)=0.0d0+c(7)*(4.0d0*a**2)+c(28)*(-22.0d0*a)
      t(29)=0.0d0+c(6)*(4.0d0*a**2)+c(29)*(-22.0d0*a)
      t(30)=0.0d0+c(7)*(4.0d0*a**2)+c(30)*(-22.0d0*a)
      t(31)=0.0d0+c(8)*(4.0d0*a**2)+c(31)*(-22.0d0*a)
      t(32)=0.0d0+c(9)*(4.0d0*a**2)+c(32)*(-22.0d0*a)
      t(33)=0.0d0+c(8)*(4.0d0*a**2)+c(10)*(4.0d0*a**2)+c(33)*(-22.0d0*a)
      t(34)=0.0d0+c(9)*(4.0d0*a**2)+c(34)*(-22.0d0*a)
      t(35)=0.0d0+c(10)*(4.0d0*a**2)+c(35)*(-22.0d0*a)
      t(36)=0.0d0+c(11)*(4.0d0*a**2)
      t(37)=0.0d0+c(12)*(4.0d0*a**2)
      t(38)=0.0d0+c(13)*(4.0d0*a**2)
      t(39)=0.0d0+c(11)*(4.0d0*a**2)+c(14)*(4.0d0*a**2)
      t(40)=0.0d0+c(15)*(4.0d0*a**2)
      t(41)=0.0d0+c(11)*(4.0d0*a**2)+c(16)*(4.0d0*a**2)
      t(42)=0.0d0+c(12)*(4.0d0*a**2)+c(17)*(4.0d0*a**2)
      t(43)=0.0d0+c(13)*(4.0d0*a**2)+c(18)*(4.0d0*a**2)
      t(44)=0.0d0+c(12)*(4.0d0*a**2)+c(19)*(4.0d0*a**2)
      t(45)=0.0d0+c(13)*(4.0d0*a**2)+c(20)*(4.0d0*a**2)
      t(46)=0.0d0+c(14)*(4.0d0*a**2)
      t(47)=0.0d0+c(15)*(4.0d0*a**2)
      t(48)=0.0d0+c(14)*(4.0d0*a**2)+c(16)*(4.0d0*a**2)
      t(49)=0.0d0+c(15)*(4.0d0*a**2)
      t(50)=0.0d0+c(16)*(4.0d0*a**2)
      t(51)=0.0d0+c(17)*(4.0d0*a**2)
      t(52)=0.0d0+c(18)*(4.0d0*a**2)
      t(53)=0.0d0+c(17)*(4.0d0*a**2)+c(19)*(4.0d0*a**2)
      t(54)=0.0d0+c(18)*(4.0d0*a**2)+c(20)*(4.0d0*a**2)
      t(55)=0.0d0+c(19)*(4.0d0*a**2)
      t(56)=0.0d0+c(20)*(4.0d0*a**2)
      t(57)=0.0d0+c(21)*(4.0d0*a**2)
      t(58)=0.0d0+c(22)*(4.0d0*a**2)
      t(59)=0.0d0+c(23)*(4.0d0*a**2)
      t(60)=0.0d0+c(21)*(4.0d0*a**2)+c(24)*(4.0d0*a**2)
      t(61)=0.0d0+c(25)*(4.0d0*a**2)
      t(62)=0.0d0+c(21)*(4.0d0*a**2)+c(26)*(4.0d0*a**2)
      t(63)=0.0d0+c(22)*(4.0d0*a**2)+c(27)*(4.0d0*a**2)
      t(64)=0.0d0+c(23)*(4.0d0*a**2)+c(28)*(4.0d0*a**2)
      t(65)=0.0d0+c(22)*(4.0d0*a**2)+c(29)*(4.0d0*a**2)
      t(66)=0.0d0+c(23)*(4.0d0*a**2)+c(30)*(4.0d0*a**2)
      t(67)=0.0d0+c(24)*(4.0d0*a**2)+c(31)*(4.0d0*a**2)
      t(68)=0.0d0+c(25)*(4.0d0*a**2)+c(32)*(4.0d0*a**2)
      t(69)=0.0d0+c(24)*(4.0d0*a**2)+c(26)*(4.0d0*a**2)+c(33)*(4.0d0*a**2)
      t(70)=0.0d0+c(25)*(4.0d0*a**2)+c(34)*(4.0d0*a**2)
      t(71)=0.0d0+c(26)*(4.0d0*a**2)+c(35)*(4.0d0*a**2)
      t(72)=0.0d0+c(27)*(4.0d0*a**2)
      t(73)=0.0d0+c(28)*(4.0d0*a**2)
      t(74)=0.0d0+c(27)*(4.0d0*a**2)+c(29)*(4.0d0*a**2)
      t(75)=0.0d0+c(28)*(4.0d0*a**2)+c(30)*(4.0d0*a**2)
      t(76)=0.0d0+c(29)*(4.0d0*a**2)
      t(77)=0.0d0+c(30)*(4.0d0*a**2)
      t(78)=0.0d0+c(31)*(4.0d0*a**2)
      t(79)=0.0d0+c(32)*(4.0d0*a**2)
      t(80)=0.0d0+c(31)*(4.0d0*a**2)+c(33)*(4.0d0*a**2)
      t(81)=0.0d0+c(32)*(4.0d0*a**2)+c(34)*(4.0d0*a**2)
      t(82)=0.0d0+c(33)*(4.0d0*a**2)+c(35)*(4.0d0*a**2)
      t(83)=0.0d0+c(34)*(4.0d0*a**2)
      t(84)=0.0d0+c(35)*(4.0d0*a**2)
    case(5)
      t(1)=0.0d0+c(1)*(-6.0d0*a)+c(5)*(2.0d0)+c(8)*(2.0d0)+c(10)*(2.0d0)
      t(2)=0.0d0+c(2)*(-10.0d0*a)+c(11)*(6.0d0)+c(14)*(2.0d0)+c(16)*(2.0d0)
      t(3)=0.0d0+c(3)*(-10.0d0*a)+c(12)*(2.0d0)+c(17)*(6.0d0)+c(19)*(2.0d0)
      t(4)=0.0d0+c(4)*(-10.0d0*a)+c(13)*(2.0d0)+c(18)*(2.0d0)+c(20)*(6.0d0)
      t(5)=0.0d0+c(1)*(4.0d0*a**2)+c(5)*(-14.0d0*a)+c(21)*(12.0d0)+c(24)*(2.0d0)+c(26)*(2.0d0)
      t(6)=0.0d0+c(6)*(-14.0d0*a)+c(22)*(6.0d0)+c(27)*(6.0d0)+c(29)*(2.0d0)
      t(7)=0.0d0+c(7)*(-14.0d0*a)+c(23)*(6.0d0)+c(28)*(2.0d0)+c(30)*(6.0d0)
      t(8)=0.0d0+c(1)*(4.0d0*a**2)+c(8)*(-14.0d0*a)+c(24)*(2.0d0)+c(31)*(12.0d0)+c(33)*(2.0d0)
      t(9)=0.0d0+c(9)*(-14.0d0*a)+c(25)*(2.0d0)+c(32)*(6.0d0)+c(34)*(6.0d0)
      t(10)=0.0d0+c(1)*(4.0d0*a**2)+c(10)*(-14.0d0*a)+c(26)*(2.0d0)+c(33)*(2.0d0)+c(35)*(12.0d0)
      t(11)=0.0d0+c(2)*(4.0d0*a**2)+c(11)*(-18.0d0*a)+c(36)*(20.0d0)+c(39)*(2.0d0)+c(41)*(2.0d0)
      t(12)=0.0d0+c(3)*(4.0d0*a**2)+c(12)*(-18.0d0*a)+c(37)*(12.0d0)+c(42)*(6.0d0)+c(44)*(2.0d0)
      t(13)=0.0d0+c(4)*(4.0d0*a**2)+c(13)*(-18.0d0*a)+c(38)*(12.0d0)+c(43)*(2.0d0)+c(45)*(6.0d0)
      t(14)=0.0d0+c(2)*(4.0d0*a**2)+c(14)*(-18.0d0*a)+c(39)*(6.0d0)+c(46)*(12.0d0)+c(48)*(2.0d0)
      t(15)=0.0d0+c(15)*(-18.0d0*a)+c(40)*(6.0d0)+c(47)*(6.0d0)+c(49)*(6.0d0)
      t(16)=0.0d0+c(2)*(4.0d0*a**2)+c(16)*(-18.0d0*a)+c(41)*(6.0d0)+c(48)*(2.0d0)+c(50)*(12.0d0)
      t(17)=0.0d0+c(3)*(4.0d0*a**2)+c(17)*(-18.0d0*a)+c(42)*(2.0d0)+c(51)*(20.0d0)+c(53)*(2.0d0)
      t(18)=0.0d0+c(4)*(4.0d0*a**2)+c(18)*(-18.0d0*a)+c(43)*(2.0d0)+c(52)*(12.0d0)+c(54)*(6.0d0)
      t(19)=0.0d0+c(3)*(4.0d0*a**2)+c(19)*(-18.0d0*a)+c(44)*(2.0d0)+c(53)*(6.0d0)+c(55)*(12.0d0)
      t(20)=0.0d0+c(4)*(4.0d0*a**2)+c(20)*(-18.0d0*a)+c(45)*(2.0d0)+c(54)*(2.0d0)+c(56)*(20.0d0)
      t(21)=0.0d0+c(5)*(4.0d0*a**2)+c(21)*(-22.0d0*a)
      t(22)=0.0d0+c(6)*(4.0d0*a**2)+c(22)*(-22.0d0*a)
      t(23)=0.0d0+c(7)*(4.0d0*a**2)+c(23)*(-22.0d0*a)
      t(24)=0.0d0+c(5)*(4.0d0*a**2)+c(8)*(4.0d0*a**2)+c(24)*(-22.0d0*a)
      t(25)=0.0d0+c(9)*(4.0d0*a**2)+c(25)*(-22.0d0*a)
      t(26)=0.0d0+c(5)*(4.0d0*a**2)+c(10)*(4.0d0*a**2)+c(26)*(-22.0d0*a)
      t(27)=0.0d0+c(6)*(4.0d0*a**2)+c(27)*(-22.0d0*a)
      t(28)=0.0d0+c(7)*(4.0d0*a**2)+c(28)*(-22.0d0*a)
      t(29)=0.0d0+c(6)*(4.0d0*a**2)+c(29)*(-22.0d0*a)
      t(30)=0.0d0+c(7)*(4.0d0*a**2)+c(30)*(-22.0d0*a)
      t(31)=0.0d0+c(8)*(4.0d0*a**2)+c(31)*(-22.0d0*a)
      t(32)=0.0d0+c(9)*(4.0d0*a**2)+c(32)*(-22.0d0*a)
      t(33)=0.0d0+c(8)*(4.0d0*a**2)+c(10)*(4.0d0*a**2)+c(33)*(-22.0d0*a)
      t(34)=0.0d0+c(9)*(4.0d0*a**2)+c(34)*(-22.0d0*a)
      t(35)=0.0d0+c(10)*(4.0d0*a**2)+c(35)*(-22.0d0*a)
      t(36)=0.0d0+c(11)*(4.0d0*a**2)+c(36)*(-26.0d0*a)
      t(37)=0.0d0+c(12)*(4.0d0*a**2)+c(37)*(-26.0d0*a)
      t(38)=0.0d0+c(13)*(4.0d0*a**2)+c(38)*(-26.0d0*a)
      t(39)=0.0d0+c(11)*(4.0d0*a**2)+c(14)*(4.0d0*a**2)+c(39)*(-26.0d0*a)
      t(40)=0.0d0+c(15)*(4.0d0*a**2)+c(40)*(-26.0d0*a)
      t(41)=0.0d0+c(11)*(4.0d0*a**2)+c(16)*(4.0d0*a**2)+c(41)*(-26.0d0*a)
      t(42)=0.0d0+c(12)*(4.0d0*a**2)+c(17)*(4.0d0*a**2)+c(42)*(-26.0d0*a)
      t(43)=0.0d0+c(13)*(4.0d0*a**2)+c(18)*(4.0d0*a**2)+c(43)*(-26.0d0*a)
      t(44)=0.0d0+c(12)*(4.0d0*a**2)+c(19)*(4.0d0*a**2)+c(44)*(-26.0d0*a)
      t(45)=0.0d0+c(13)*(4.0d0*a**2)+c(20)*(4.0d0*a**2)+c(45)*(-26.0d0*a)
      t(46)=0.0d0+c(14)*(4.0d0*a**2)+c(46)*(-26.0d0*a)
      t(47)=0.0d0+c(15)*(4.0d0*a**2)+c(47)*(-26.0d0*a)
      t(48)=0.0d0+c(14)*(4.0d0*a**2)+c(16)*(4.0d0*a**2)+c(48)*(-26.0d0*a)
      t(49)=0.0d0+c(15)*(4.0d0*a**2)+c(49)*(-26.0d0*a)
      t(50)=0.0d0+c(16)*(4.0d0*a**2)+c(50)*(-26.0d0*a)
      t(51)=0.0d0+c(17)*(4.0d0*a**2)+c(51)*(-26.0d0*a)
      t(52)=0.0d0+c(18)*(4.0d0*a**2)+c(52)*(-26.0d0*a)
      t(53)=0.0d0+c(17)*(4.0d0*a**2)+c(19)*(4.0d0*a**2)+c(53)*(-26.0d0*a)
      t(54)=0.0d0+c(18)*(4.0d0*a**2)+c(20)*(4.0d0*a**2)+c(54)*(-26.0d0*a)
      t(55)=0.0d0+c(19)*(4.0d0*a**2)+c(55)*(-26.0d0*a)
      t(56)=0.0d0+c(20)*(4.0d0*a**2)+c(56)*(-26.0d0*a)
      t(57)=0.0d0+c(21)*(4.0d0*a**2)
      t(58)=0.0d0+c(22)*(4.0d0*a**2)
      t(59)=0.0d0+c(23)*(4.0d0*a**2)
      t(60)=0.0d0+c(21)*(4.0d0*a**2)+c(24)*(4.0d0*a**2)
      t(61)=0.0d0+c(25)*(4.0d0*a**2)
      t(62)=0.0d0+c(21)*(4.0d0*a**2)+c(26)*(4.0d0*a**2)
      t(63)=0.0d0+c(22)*(4.0d0*a**2)+c(27)*(4.0d0*a**2)
      t(64)=0.0d0+c(23)*(4.0d0*a**2)+c(28)*(4.0d0*a**2)
      t(65)=0.0d0+c(22)*(4.0d0*a**2)+c(29)*(4.0d0*a**2)
      t(66)=0.0d0+c(23)*(4.0d0*a**2)+c(30)*(4.0d0*a**2)
      t(67)=0.0d0+c(24)*(4.0d0*a**2)+c(31)*(4.0d0*a**2)
      t(68)=0.0d0+c(25)*(4.0d0*a**2)+c(32)*(4.0d0*a**2)
      t(69)=0.0d0+c(24)*(4.0d0*a**2)+c(26)*(4.0d0*a**2)+c(33)*(4.0d0*a**2)
      t(70)=0.0d0+c(25)*(4.0d0*a**2)+c(34)*(4.0d0*a**2)
      t(71)=0.0d0+c(26)*(4.0d0*a**2)+c(35)*(4.0d0*a**2)
      t(72)=0.0d0+c(27)*(4.0d0*a**2)
      t(73)=0.0d0+c(28)*(4.0d0*a**2)
      t(74)=0.0d0+c(27)*(4.0d0*a**2)+c(29)*(4.0d0*a**2)
      t(75)=0.0d0+c(28)*(4.0d0*a**2)+c(30)*(4.0d0*a**2)
      t(76)=0.0d0+c(29)*(4.0d0*a**2)
      t(77)=0.0d0+c(30)*(4.0d0*a**2)
      t(78)=0.0d0+c(31)*(4.0d0*a**2)
      t(79)=0.0d0+c(32)*(4.0d0*a**2)
      t(80)=0.0d0+c(31)*(4.0d0*a**2)+c(33)*(4.0d0*a**2)
      t(81)=0.0d0+c(32)*(4.0d0*a**2)+c(34)*(4.0d0*a**2)
      t(82)=0.0d0+c(33)*(4.0d0*a**2)+c(35)*(4.0d0*a**2)
      t(83)=0.0d0+c(34)*(4.0d0*a**2)
      t(84)=0.0d0+c(35)*(4.0d0*a**2)
      t(85)=0.0d0+c(36)*(4.0d0*a**2)
      t(86)=0.0d0+c(37)*(4.0d0*a**2)
      t(87)=0.0d0+c(38)*(4.0d0*a**2)
      t(88)=0.0d0+c(36)*(4.0d0*a**2)+c(39)*(4.0d0*a**2)
      t(89)=0.0d0+c(40)*(4.0d0*a**2)
      t(90)=0.0d0+c(36)*(4.0d0*a**2)+c(41)*(4.0d0*a**2)
      t(91)=0.0d0+c(37)*(4.0d0*a**2)+c(42)*(4.0d0*a**2)
      t(92)=0.0d0+c(38)*(4.0d0*a**2)+c(43)*(4.0d0*a**2)
      t(93)=0.0d0+c(37)*(4.0d0*a**2)+c(44)*(4.0d0*a**2)
      t(94)=0.0d0+c(38)*(4.0d0*a**2)+c(45)*(4.0d0*a**2)
      t(95)=0.0d0+c(39)*(4.0d0*a**2)+c(46)*(4.0d0*a**2)
      t(96)=0.0d0+c(40)*(4.0d0*a**2)+c(47)*(4.0d0*a**2)
      t(97)=0.0d0+c(39)*(4.0d0*a**2)+c(41)*(4.0d0*a**2)+c(48)*(4.0d0*a**2)
      t(98)=0.0d0+c(40)*(4.0d0*a**2)+c(49)*(4.0d0*a**2)
      t(99)=0.0d0+c(41)*(4.0d0*a**2)+c(50)*(4.0d0*a**2)
      t(100)=0.0d0+c(42)*(4.0d0*a**2)+c(51)*(4.0d0*a**2)
      t(101)=0.0d0+c(43)*(4.0d0*a**2)+c(52)*(4.0d0*a**2)
      t(102)=0.0d0+c(42)*(4.0d0*a**2)+c(44)*(4.0d0*a**2)+c(53)*(4.0d0*a**2)
      t(103)=0.0d0+c(43)*(4.0d0*a**2)+c(45)*(4.0d0*a**2)+c(54)*(4.0d0*a**2)
      t(104)=0.0d0+c(44)*(4.0d0*a**2)+c(55)*(4.0d0*a**2)
      t(105)=0.0d0+c(45)*(4.0d0*a**2)+c(56)*(4.0d0*a**2)
      t(106)=0.0d0+c(46)*(4.0d0*a**2)
      t(107)=0.0d0+c(47)*(4.0d0*a**2)
      t(108)=0.0d0+c(46)*(4.0d0*a**2)+c(48)*(4.0d0*a**2)
      t(109)=0.0d0+c(47)*(4.0d0*a**2)+c(49)*(4.0d0*a**2)
      t(110)=0.0d0+c(48)*(4.0d0*a**2)+c(50)*(4.0d0*a**2)
      t(111)=0.0d0+c(49)*(4.0d0*a**2)
      t(112)=0.0d0+c(50)*(4.0d0*a**2)
      t(113)=0.0d0+c(51)*(4.0d0*a**2)
      t(114)=0.0d0+c(52)*(4.0d0*a**2)
      t(115)=0.0d0+c(51)*(4.0d0*a**2)+c(53)*(4.0d0*a**2)
      t(116)=0.0d0+c(52)*(4.0d0*a**2)+c(54)*(4.0d0*a**2)
      t(117)=0.0d0+c(53)*(4.0d0*a**2)+c(55)*(4.0d0*a**2)
      t(118)=0.0d0+c(54)*(4.0d0*a**2)+c(56)*(4.0d0*a**2)
      t(119)=0.0d0+c(55)*(4.0d0*a**2)
      t(120)=0.0d0+c(56)*(4.0d0*a**2)
    case(6)
      t(1)=0.0d0+c(1)*(-6.0d0*a)+c(5)*(2.0d0)+c(8)*(2.0d0)+c(10)*(2.0d0)
      t(2)=0.0d0+c(2)*(-10.0d0*a)+c(11)*(6.0d0)+c(14)*(2.0d0)+c(16)*(2.0d0)
      t(3)=0.0d0+c(3)*(-10.0d0*a)+c(12)*(2.0d0)+c(17)*(6.0d0)+c(19)*(2.0d0)
      t(4)=0.0d0+c(4)*(-10.0d0*a)+c(13)*(2.0d0)+c(18)*(2.0d0)+c(20)*(6.0d0)
      t(5)=0.0d0+c(1)*(4.0d0*a**2)+c(5)*(-14.0d0*a)+c(21)*(12.0d0)+c(24)*(2.0d0)+c(26)*(2.0d0)
      t(6)=0.0d0+c(6)*(-14.0d0*a)+c(22)*(6.0d0)+c(27)*(6.0d0)+c(29)*(2.0d0)
      t(7)=0.0d0+c(7)*(-14.0d0*a)+c(23)*(6.0d0)+c(28)*(2.0d0)+c(30)*(6.0d0)
      t(8)=0.0d0+c(1)*(4.0d0*a**2)+c(8)*(-14.0d0*a)+c(24)*(2.0d0)+c(31)*(12.0d0)+c(33)*(2.0d0)
      t(9)=0.0d0+c(9)*(-14.0d0*a)+c(25)*(2.0d0)+c(32)*(6.0d0)+c(34)*(6.0d0)
      t(10)=0.0d0+c(1)*(4.0d0*a**2)+c(10)*(-14.0d0*a)+c(26)*(2.0d0)+c(33)*(2.0d0)+c(35)*(12.0d0)
      t(11)=0.0d0+c(2)*(4.0d0*a**2)+c(11)*(-18.0d0*a)+c(36)*(20.0d0)+c(39)*(2.0d0)+c(41)*(2.0d0)
      t(12)=0.0d0+c(3)*(4.0d0*a**2)+c(12)*(-18.0d0*a)+c(37)*(12.0d0)+c(42)*(6.0d0)+c(44)*(2.0d0)
      t(13)=0.0d0+c(4)*(4.0d0*a**2)+c(13)*(-18.0d0*a)+c(38)*(12.0d0)+c(43)*(2.0d0)+c(45)*(6.0d0)
      t(14)=0.0d0+c(2)*(4.0d0*a**2)+c(14)*(-18.0d0*a)+c(39)*(6.0d0)+c(46)*(12.0d0)+c(48)*(2.0d0)
      t(15)=0.0d0+c(15)*(-18.0d0*a)+c(40)*(6.0d0)+c(47)*(6.0d0)+c(49)*(6.0d0)
      t(16)=0.0d0+c(2)*(4.0d0*a**2)+c(16)*(-18.0d0*a)+c(41)*(6.0d0)+c(48)*(2.0d0)+c(50)*(12.0d0)
      t(17)=0.0d0+c(3)*(4.0d0*a**2)+c(17)*(-18.0d0*a)+c(42)*(2.0d0)+c(51)*(20.0d0)+c(53)*(2.0d0)
      t(18)=0.0d0+c(4)*(4.0d0*a**2)+c(18)*(-18.0d0*a)+c(43)*(2.0d0)+c(52)*(12.0d0)+c(54)*(6.0d0)
      t(19)=0.0d0+c(3)*(4.0d0*a**2)+c(19)*(-18.0d0*a)+c(44)*(2.0d0)+c(53)*(6.0d0)+c(55)*(12.0d0)
      t(20)=0.0d0+c(4)*(4.0d0*a**2)+c(20)*(-18.0d0*a)+c(45)*(2.0d0)+c(54)*(2.0d0)+c(56)*(20.0d0)
      t(21)=0.0d0+c(5)*(4.0d0*a**2)+c(21)*(-22.0d0*a)+c(57)*(30.0d0)+c(60)*(2.0d0)+c(62)*(2.0d0)
      t(22)=0.0d0+c(6)*(4.0d0*a**2)+c(22)*(-22.0d0*a)+c(58)*(20.0d0)+c(63)*(6.0d0)+c(65)*(2.0d0)
      t(23)=0.0d0+c(7)*(4.0d0*a**2)+c(23)*(-22.0d0*a)+c(59)*(20.0d0)+c(64)*(2.0d0)+c(66)*(6.0d0)
      t(24)=0.0d0+c(5)*(4.0d0*a**2)+c(8)*(4.0d0*a**2)+c(24)*(-22.0d0*a)+c(60)*(12.0d0)+c(67)*(12.0d0)+c(69)*(2.0d0)
      t(25)=0.0d0+c(9)*(4.0d0*a**2)+c(25)*(-22.0d0*a)+c(61)*(12.0d0)+c(68)*(6.0d0)+c(70)*(6.0d0)
      t(26)=0.0d0+c(5)*(4.0d0*a**2)+c(10)*(4.0d0*a**2)+c(26)*(-22.0d0*a)+c(62)*(12.0d0)+c(69)*(2.0d0)+c(71)*(12.0d0)
      t(27)=0.0d0+c(6)*(4.0d0*a**2)+c(27)*(-22.0d0*a)+c(63)*(6.0d0)+c(72)*(20.0d0)+c(74)*(2.0d0)
      t(28)=0.0d0+c(7)*(4.0d0*a**2)+c(28)*(-22.0d0*a)+c(64)*(6.0d0)+c(73)*(12.0d0)+c(75)*(6.0d0)
      t(29)=0.0d0+c(6)*(4.0d0*a**2)+c(29)*(-22.0d0*a)+c(65)*(6.0d0)+c(74)*(6.0d0)+c(76)*(12.0d0)
      t(30)=0.0d0+c(7)*(4.0d0*a**2)+c(30)*(-22.0d0*a)+c(66)*(6.0d0)+c(75)*(2.0d0)+c(77)*(20.0d0)
      t(31)=0.0d0+c(8)*(4.0d0*a**2)+c(31)*(-22.0d0*a)+c(67)*(2.0d0)+c(78)*(30.0d0)+c(80)*(2.0d0)
      t(32)=0.0d0+c(9)*(4.0d0*a**2)+c(32)*(-22.0d0*a)+c(68)*(2.0d0)+c(79)*(20.0d0)+c(81)*(6.0d0)
      t(33)=0.0d0+c(8)*(4.0d0*a**2)+c(10)*(4.0d0*a**2)+c(33)*(-22.0d0*a)+c(69)*(2.0d0)+c(80)*(12.0d0)+c(82)*(12.0d0)
      t(34)=0.0d0+c(9)*(4.0d0*a**2)+c(34)*(-22.0d0*a)+c(70)*(2.0d0)+c(81)*(6.0d0)+c(83)*(20.0d0)
      t(35)=0.0d0+c(10)*(4.0d0*a**2)+c(35)*(-22.0d0*a)+c(71)*(2.0d0)+c(82)*(2.0d0)+c(84)*(30.0d0)
      t(36)=0.0d0+c(11)*(4.0d0*a**2)+c(36)*(-26.0d0*a)
      t(37)=0.0d0+c(12)*(4.0d0*a**2)+c(37)*(-26.0d0*a)
      t(38)=0.0d0+c(13)*(4.0d0*a**2)+c(38)*(-26.0d0*a)
      t(39)=0.0d0+c(11)*(4.0d0*a**2)+c(14)*(4.0d0*a**2)+c(39)*(-26.0d0*a)
      t(40)=0.0d0+c(15)*(4.0d0*a**2)+c(40)*(-26.0d0*a)
      t(41)=0.0d0+c(11)*(4.0d0*a**2)+c(16)*(4.0d0*a**2)+c(41)*(-26.0d0*a)
      t(42)=0.0d0+c(12)*(4.0d0*a**2)+c(17)*(4.0d0*a**2)+c(42)*(-26.0d0*a)
      t(43)=0.0d0+c(13)*(4.0d0*a**2)+c(18)*(4.0d0*a**2)+c(43)*(-26.0d0*a)
      t(44)=0.0d0+c(12)*(4.0d0*a**2)+c(19)*(4.0d0*a**2)+c(44)*(-26.0d0*a)
      t(45)=0.0d0+c(13)*(4.0d0*a**2)+c(20)*(4.0d0*a**2)+c(45)*(-26.0d0*a)
      t(46)=0.0d0+c(14)*(4.0d0*a**2)+c(46)*(-26.0d0*a)
      t(47)=0.0d0+c(15)*(4.0d0*a**2)+c(47)*(-26.0d0*a)
      t(48)=0.0d0+c(14)*(4.0d0*a**2)+c(16)*(4.0d0*a**2)+c(48)*(-26.0d0*a)
      t(49)=0.0d0+c(15)*(4.0d0*a**2)+c(49)*(-26.0d0*a)
      t(50)=0.0d0+c(16)*(4.0d0*a**2)+c(50)*(-26.0d0*a)
      t(51)=0.0d0+c(17)*(4.0d0*a**2)+c(51)*(-26.0d0*a)
      t(52)=0.0d0+c(18)*(4.0d0*a**2)+c(52)*(-26.0d0*a)
      t(53)=0.0d0+c(17)*(4.0d0*a**2)+c(19)*(4.0d0*a**2)+c(53)*(-26.0d0*a)
      t(54)=0.0d0+c(18)*(4.0d0*a**2)+c(20)*(4.0d0*a**2)+c(54)*(-26.0d0*a)
      t(55)=0.0d0+c(19)*(4.0d0*a**2)+c(55)*(-26.0d0*a)
      t(56)=0.0d0+c(20)*(4.0d0*a**2)+c(56)*(-26.0d0*a)
      t(57)=0.0d0+c(21)*(4.0d0*a**2)+c(57)*(-30.0d0*a)
      t(58)=0.0d0+c(22)*(4.0d0*a**2)+c(58)*(-30.0d0*a)
      t(59)=0.0d0+c(23)*(4.0d0*a**2)+c(59)*(-30.0d0*a)
      t(60)=0.0d0+c(21)*(4.0d0*a**2)+c(24)*(4.0d0*a**2)+c(60)*(-30.0d0*a)
      t(61)=0.0d0+c(25)*(4.0d0*a**2)+c(61)*(-30.0d0*a)
      t(62)=0.0d0+c(21)*(4.0d0*a**2)+c(26)*(4.0d0*a**2)+c(62)*(-30.0d0*a)
      t(63)=0.0d0+c(22)*(4.0d0*a**2)+c(27)*(4.0d0*a**2)+c(63)*(-30.0d0*a)
      t(64)=0.0d0+c(23)*(4.0d0*a**2)+c(28)*(4.0d0*a**2)+c(64)*(-30.0d0*a)
      t(65)=0.0d0+c(22)*(4.0d0*a**2)+c(29)*(4.0d0*a**2)+c(65)*(-30.0d0*a)
      t(66)=0.0d0+c(23)*(4.0d0*a**2)+c(30)*(4.0d0*a**2)+c(66)*(-30.0d0*a)
      t(67)=0.0d0+c(24)*(4.0d0*a**2)+c(31)*(4.0d0*a**2)+c(67)*(-30.0d0*a)
      t(68)=0.0d0+c(25)*(4.0d0*a**2)+c(32)*(4.0d0*a**2)+c(68)*(-30.0d0*a)
      t(69)=0.0d0+c(24)*(4.0d0*a**2)+c(26)*(4.0d0*a**2)+c(33)*(4.0d0*a**2)+c(69)*(-30.0d0*a)
      t(70)=0.0d0+c(25)*(4.0d0*a**2)+c(34)*(4.0d0*a**2)+c(70)*(-30.0d0*a)
      t(71)=0.0d0+c(26)*(4.0d0*a**2)+c(35)*(4.0d0*a**2)+c(71)*(-30.0d0*a)
      t(72)=0.0d0+c(27)*(4.0d0*a**2)+c(72)*(-30.0d0*a)
      t(73)=0.0d0+c(28)*(4.0d0*a**2)+c(73)*(-30.0d0*a)
      t(74)=0.0d0+c(27)*(4.0d0*a**2)+c(29)*(4.0d0*a**2)+c(74)*(-30.0d0*a)
      t(75)=0.0d0+c(28)*(4.0d0*a**2)+c(30)*(4.0d0*a**2)+c(75)*(-30.0d0*a)
      t(76)=0.0d0+c(29)*(4.0d0*a**2)+c(76)*(-30.0d0*a)
      t(77)=0.0d0+c(30)*(4.0d0*a**2)+c(77)*(-30.0d0*a)
      t(78)=0.0d0+c(31)*(4.0d0*a**2)+c(78)*(-30.0d0*a)
      t(79)=0.0d0+c(32)*(4.0d0*a**2)+c(79)*(-30.0d0*a)
      t(80)=0.0d0+c(31)*(4.0d0*a**2)+c(33)*(4.0d0*a**2)+c(80)*(-30.0d0*a)
      t(81)=0.0d0+c(32)*(4.0d0*a**2)+c(34)*(4.0d0*a**2)+c(81)*(-30.0d0*a)
      t(82)=0.0d0+c(33)*(4.0d0*a**2)+c(35)*(4.0d0*a**2)+c(82)*(-30.0d0*a)
      t(83)=0.0d0+c(34)*(4.0d0*a**2)+c(83)*(-30.0d0*a)
      t(84)=0.0d0+c(35)*(4.0d0*a**2)+c(84)*(-30.0d0*a)
      t(85)=0.0d0+c(36)*(4.0d0*a**2)
      t(86)=0.0d0+c(37)*(4.0d0*a**2)
      t(87)=0.0d0+c(38)*(4.0d0*a**2)
      t(88)=0.0d0+c(36)*(4.0d0*a**2)+c(39)*(4.0d0*a**2)
      t(89)=0.0d0+c(40)*(4.0d0*a**2)
      t(90)=0.0d0+c(36)*(4.0d0*a**2)+c(41)*(4.0d0*a**2)
      t(91)=0.0d0+c(37)*(4.0d0*a**2)+c(42)*(4.0d0*a**2)
      t(92)=0.0d0+c(38)*(4.0d0*a**2)+c(43)*(4.0d0*a**2)
      t(93)=0.0d0+c(37)*(4.0d0*a**2)+c(44)*(4.0d0*a**2)
      t(94)=0.0d0+c(38)*(4.0d0*a**2)+c(45)*(4.0d0*a**2)
      t(95)=0.0d0+c(39)*(4.0d0*a**2)+c(46)*(4.0d0*a**2)
      t(96)=0.0d0+c(40)*(4.0d0*a**2)+c(47)*(4.0d0*a**2)
      t(97)=0.0d0+c(39)*(4.0d0*a**2)+c(41)*(4.0d0*a**2)+c(48)*(4.0d0*a**2)
      t(98)=0.0d0+c(40)*(4.0d0*a**2)+c(49)*(4.0d0*a**2)
      t(99)=0.0d0+c(41)*(4.0d0*a**2)+c(50)*(4.0d0*a**2)
      t(100)=0.0d0+c(42)*(4.0d0*a**2)+c(51)*(4.0d0*a**2)
      t(101)=0.0d0+c(43)*(4.0d0*a**2)+c(52)*(4.0d0*a**2)
      t(102)=0.0d0+c(42)*(4.0d0*a**2)+c(44)*(4.0d0*a**2)+c(53)*(4.0d0*a**2)
      t(103)=0.0d0+c(43)*(4.0d0*a**2)+c(45)*(4.0d0*a**2)+c(54)*(4.0d0*a**2)
      t(104)=0.0d0+c(44)*(4.0d0*a**2)+c(55)*(4.0d0*a**2)
      t(105)=0.0d0+c(45)*(4.0d0*a**2)+c(56)*(4.0d0*a**2)
      t(106)=0.0d0+c(46)*(4.0d0*a**2)
      t(107)=0.0d0+c(47)*(4.0d0*a**2)
      t(108)=0.0d0+c(46)*(4.0d0*a**2)+c(48)*(4.0d0*a**2)
      t(109)=0.0d0+c(47)*(4.0d0*a**2)+c(49)*(4.0d0*a**2)
      t(110)=0.0d0+c(48)*(4.0d0*a**2)+c(50)*(4.0d0*a**2)
      t(111)=0.0d0+c(49)*(4.0d0*a**2)
      t(112)=0.0d0+c(50)*(4.0d0*a**2)
      t(113)=0.0d0+c(51)*(4.0d0*a**2)
      t(114)=0.0d0+c(52)*(4.0d0*a**2)
      t(115)=0.0d0+c(51)*(4.0d0*a**2)+c(53)*(4.0d0*a**2)
      t(116)=0.0d0+c(52)*(4.0d0*a**2)+c(54)*(4.0d0*a**2)
      t(117)=0.0d0+c(53)*(4.0d0*a**2)+c(55)*(4.0d0*a**2)
      t(118)=0.0d0+c(54)*(4.0d0*a**2)+c(56)*(4.0d0*a**2)
      t(119)=0.0d0+c(55)*(4.0d0*a**2)
      t(120)=0.0d0+c(56)*(4.0d0*a**2)
      t(121)=0.0d0+c(57)*(4.0d0*a**2)
      t(122)=0.0d0+c(58)*(4.0d0*a**2)
      t(123)=0.0d0+c(59)*(4.0d0*a**2)
      t(124)=0.0d0+c(57)*(4.0d0*a**2)+c(60)*(4.0d0*a**2)
      t(125)=0.0d0+c(61)*(4.0d0*a**2)
      t(126)=0.0d0+c(57)*(4.0d0*a**2)+c(62)*(4.0d0*a**2)
      t(127)=0.0d0+c(58)*(4.0d0*a**2)+c(63)*(4.0d0*a**2)
      t(128)=0.0d0+c(59)*(4.0d0*a**2)+c(64)*(4.0d0*a**2)
      t(129)=0.0d0+c(58)*(4.0d0*a**2)+c(65)*(4.0d0*a**2)
      t(130)=0.0d0+c(59)*(4.0d0*a**2)+c(66)*(4.0d0*a**2)
      t(131)=0.0d0+c(60)*(4.0d0*a**2)+c(67)*(4.0d0*a**2)
      t(132)=0.0d0+c(61)*(4.0d0*a**2)+c(68)*(4.0d0*a**2)
      t(133)=0.0d0+c(60)*(4.0d0*a**2)+c(62)*(4.0d0*a**2)+c(69)*(4.0d0*a**2)
      t(134)=0.0d0+c(61)*(4.0d0*a**2)+c(70)*(4.0d0*a**2)
      t(135)=0.0d0+c(62)*(4.0d0*a**2)+c(71)*(4.0d0*a**2)
      t(136)=0.0d0+c(63)*(4.0d0*a**2)+c(72)*(4.0d0*a**2)
      t(137)=0.0d0+c(64)*(4.0d0*a**2)+c(73)*(4.0d0*a**2)
      t(138)=0.0d0+c(63)*(4.0d0*a**2)+c(65)*(4.0d0*a**2)+c(74)*(4.0d0*a**2)
      t(139)=0.0d0+c(64)*(4.0d0*a**2)+c(66)*(4.0d0*a**2)+c(75)*(4.0d0*a**2)
      t(140)=0.0d0+c(65)*(4.0d0*a**2)+c(76)*(4.0d0*a**2)
      t(141)=0.0d0+c(66)*(4.0d0*a**2)+c(77)*(4.0d0*a**2)
      t(142)=0.0d0+c(67)*(4.0d0*a**2)+c(78)*(4.0d0*a**2)
      t(143)=0.0d0+c(68)*(4.0d0*a**2)+c(79)*(4.0d0*a**2)
      t(144)=0.0d0+c(67)*(4.0d0*a**2)+c(69)*(4.0d0*a**2)+c(80)*(4.0d0*a**2)
      t(145)=0.0d0+c(68)*(4.0d0*a**2)+c(70)*(4.0d0*a**2)+c(81)*(4.0d0*a**2)
      t(146)=0.0d0+c(69)*(4.0d0*a**2)+c(71)*(4.0d0*a**2)+c(82)*(4.0d0*a**2)
      t(147)=0.0d0+c(70)*(4.0d0*a**2)+c(83)*(4.0d0*a**2)
      t(148)=0.0d0+c(71)*(4.0d0*a**2)+c(84)*(4.0d0*a**2)
      t(149)=0.0d0+c(72)*(4.0d0*a**2)
      t(150)=0.0d0+c(73)*(4.0d0*a**2)
      t(151)=0.0d0+c(72)*(4.0d0*a**2)+c(74)*(4.0d0*a**2)
      t(152)=0.0d0+c(73)*(4.0d0*a**2)+c(75)*(4.0d0*a**2)
      t(153)=0.0d0+c(74)*(4.0d0*a**2)+c(76)*(4.0d0*a**2)
      t(154)=0.0d0+c(75)*(4.0d0*a**2)+c(77)*(4.0d0*a**2)
      t(155)=0.0d0+c(76)*(4.0d0*a**2)
      t(156)=0.0d0+c(77)*(4.0d0*a**2)
      t(157)=0.0d0+c(78)*(4.0d0*a**2)
      t(158)=0.0d0+c(79)*(4.0d0*a**2)
      t(159)=0.0d0+c(78)*(4.0d0*a**2)+c(80)*(4.0d0*a**2)
      t(160)=0.0d0+c(79)*(4.0d0*a**2)+c(81)*(4.0d0*a**2)
      t(161)=0.0d0+c(80)*(4.0d0*a**2)+c(82)*(4.0d0*a**2)
      t(162)=0.0d0+c(81)*(4.0d0*a**2)+c(83)*(4.0d0*a**2)
      t(163)=0.0d0+c(82)*(4.0d0*a**2)+c(84)*(4.0d0*a**2)
      t(164)=0.0d0+c(83)*(4.0d0*a**2)
      t(165)=0.0d0+c(84)*(4.0d0*a**2)
    case default
      print*,'Error: R_Laplacian not implemented for lmax>',lmax
  end select
 
end subroutine

end module
