!> @file
!!  Old kinetic convolution routines
!! @deprecated
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


         subroutine Convolkinetic(n1,n2,n3, &
               nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
               cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_fc,x_f,y_c,y_f)
!   y = (kinetic energy operator)x + (cprec*I)x 
    implicit real(kind=8) (a-h,o-z)
    logical :: firstcall=.true. 
    integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
    dimension x_c(0:n1,0:n2,0:n3),y_c(0:n1,0:n2,0:n3)
    dimension x_fc(0:n1,0:n2,0:n3,3),x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
    dimension y_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
    dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
    dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)


    parameter(lowfil=-14,lupfil=14)
    dimension a(lowfil-1:lupfil+1),b(lowfil-1:lupfil+1),c(lowfil-3:lupfil+3),e(lowfil:lupfil)
    scale=-.5d0/hgrid**2
!---------------------------------------------------------------------------
! second derivative filters for Daubechies 16
!  <phi|D^2|phi_i>
    a(0)=   -3.5536922899131901941296809374d0*scale
    a(1)=    2.2191465938911163898794546405d0*scale
    a(2)=   -0.6156141465570069496314853949d0*scale
    a(3)=    0.2371780582153805636239247476d0*scale
    a(4)=   -0.0822663999742123340987663521d0*scale
    a(5)=    0.02207029188482255523789911295638968409d0*scale
    a(6)=   -0.409765689342633823899327051188315485d-2*scale
    a(7)=    0.45167920287502235349480037639758496d-3*scale
    a(8)=   -0.2398228524507599670405555359023135d-4*scale
    a(9)=    2.0904234952920365957922889447361d-6*scale
    a(10)=  -3.7230763047369275848791496973044d-7*scale
    a(11)=  -1.05857055496741470373494132287d-8*scale
    a(12)=  -5.813879830282540547959250667d-11*scale
    a(13)=   2.70800493626319438269856689037647576d-13*scale
    a(14)=  -6.924474940639200152025730585882d-18*scale
    do i=1,14
        a(-i)=a(i)
    enddo
    a( 15)=0.d0
    b( 15)=0.d0
    a(-15)=0.d0
    b(-15)=0.d0
!  <phi|D^2|psi_i>
    c(-14)=     -3.869102413147656535541850057188d-18*scale
    c(-13)=      1.5130616560866154733900029272077362d-13*scale
    c(-12)=     -3.2264702314010525539061647271983988409d-11*scale
    c(-11)=     -5.96264938781402337319841002642d-9*scale
    c(-10)=     -2.1656830629214041470164889350342d-7*scale
    c(-9 )=      8.7969704055286288323596890609625d-7*scale
    c(-8 )=     -0.00001133456724516819987751818232711775d0*scale
    c(-7 )=      0.00021710795484646138591610188464622454d0*scale
    c(-6 )=     -0.0021356291838797986414312219042358542d0*scale
    c(-5 )=      0.00713761218453631422925717625758502986d0*scale
    c(-4 )=     -0.0284696165863973422636410524436931061d0*scale
    c(-3 )=      0.14327329352510759457155821037742893841d0*scale
    c(-2 )=     -0.42498050943780130143385739554118569733d0*scale
    c(-1 )=      0.65703074007121357894896358254040272157d0*scale
    c( 0 )=     -0.42081655293724308770919536332797729898d0*scale
    c( 1 )=     -0.21716117505137104371463587747283267899d0*scale
    c( 2 )=      0.63457035267892488185929915286969303251d0*scale
    c( 3 )=     -0.53298223962800395684936080758073568406d0*scale
    c( 4 )=      0.23370490631751294307619384973520033236d0*scale
    c( 5 )=     -0.05657736973328755112051544344507997075d0*scale
    c( 6 )=      0.0080872029411844780634067667008050127d0*scale
    c( 7 )=     -0.00093423623304808664741804536808932984d0*scale
    c( 8 )=      0.00005075807947289728306309081261461095d0*scale
    c( 9 )=     -4.62561497463184262755416490048242d-6*scale
    c( 10)=      6.3919128513793415587294752371778d-7*scale
    c( 11)=      1.87909235155149902916133888931d-8*scale
    c( 12)=      1.04757345962781829480207861447155543883d-10*scale
    c( 13)=     -4.84665690596158959648731537084025836d-13*scale
    c( 14)=      1.2392629629188986192855777620877d-17*scale
    c(-15)=0.d0
    c(-16)=0.d0
    c(-17)=0.d0
    c( 15)=0.d0
    c( 16)=0.d0
    c( 17)=0.d0
!  <psi|D^2|phi_i>
    do i=-14,14
        b(i)=c(-i)
    enddo
    !<psi|D^2|psi_i>
    e(0)=   -24.875846029392331358907766562d0*scale
    e(1)=   -7.1440597663471719869313377994d0*scale
    e(2)=   -0.04251705323669172315864542163525830944d0*scale
    e(3)=   -0.26995931336279126953587091167128839196d0*scale
    e(4)=    0.08207454169225172612513390763444496516d0*scale
    e(5)=   -0.02207327034586634477996701627614752761d0*scale
    e(6)=    0.00409765642831595181639002667514310145d0*scale
    e(7)=   -0.00045167920287507774929432548999880117d0*scale
    e(8)=    0.00002398228524507599670405555359023135d0*scale
    e(9)=   -2.0904234952920365957922889447361d-6*scale
    e(10)=   3.7230763047369275848791496973044d-7*scale
    e(11)=   1.05857055496741470373494132287d-8*scale
    e(12)=   5.8138798302825405479592506674648873655d-11*scale
    e(13)=  -2.70800493626319438269856689037647576d-13*scale
    e(14)=   6.924474940639200152025730585882d-18*scale
    do i=1,14
        e(-i)=e(i)
    enddo


!  if (firstcall) then
!
!! (1/2) d^2/dx^2
!    mflop1=0
!    do i3=0,n3
!    do i2=0,n2
!        do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
!                do l=max(ibyz_c(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_c(2,i2,i3)-i1)
!                    mflop1=mflop1+4
!                enddo
!                    mflop1=mflop1+3
!        enddo
!
!        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
!                do l=max(ibyz_c(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_c(2,i2,i3)-i1)
!                    mflop1=mflop1+2
!                enddo
!        enddo
!    enddo
!    enddo
!
!! + (1/2) d^2/dy^2
!    mflop2=0
!    do i3=0,n3
!    do i1=0,n1
!        do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
!                do l=max(ibxz_c(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_c(2,i1,i3)-i2)
!                    mflop2=mflop2+4
!                enddo
!                    mflop2=mflop2+2
!        enddo
!
!        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
!                do l=max(ibxz_c(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_c(2,i1,i3)-i2)
!                    mflop2=mflop2+2
!                enddo
!        enddo
!    enddo
!    enddo
!
!
!! + (1/2) d^2/dz^2
!    mflop3=0
!    do i2=0,n2
!    do i1=0,n1
!        do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
!                do l=max(ibxy_c(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_c(2,i1,i2)-i3)
!                    mflop3=mflop3+4
!                enddo
!                    mflop3=mflop3+2
!        enddo
!
!        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
!                do l=max(ibxy_c(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_c(2,i1,i2)-i3)
!                    mflop3=mflop3+2
!                enddo
!        enddo
!    enddo
!    enddo
!
!! wavelet part
! ! (1/2) d^2/dx^2
!    nflop1=0
!    do i3=nfl3,nfu3
!        do i2=nfl2,nfu2
!            do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
!                do l=max(ibyz_f(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_f(2,i2,i3)-i1)
!                    nflop1=nflop1+26
!                enddo
!                    nflop1=nflop1+17
!            enddo
!        enddo
!    enddo
!
! ! + (1/2) d^2/dy^2
!    nflop2=0
!    do i3=nfl3,nfu3
!    do i1=nfl1,nfu1
!       do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
!                do l=max(ibxz_f(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_f(2,i1,i3)-i2)
!                    nflop2=nflop2+26
!                enddo
!                    nflop2=nflop2+7
!       enddo
!    enddo
!    enddo
!
! ! + (1/2) d^2/dz^2
!    nflop3=0
!    do i2=nfl2,nfu2
!    do i1=nfl1,nfu1
!       do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
!                do l=max(ibxy_f(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_f(2,i1,i2)-i3)
!                    nflop3=nflop3+26
!                enddo
!                    nflop3=nflop3+7
!       enddo
!    enddo
!    enddo
!    
!    firstcall=.false.
!    endif


!---------------------------------------------------------------------------
                                                            
! Scaling function part

!       call system_clock(ncount0,ncount_rate,ncount_max)

! (1/2) d^2/dx^2
    do i3=0,n3
    do i2=0,n2
        do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-1,2
                t111=0.d0 ; s111=0.d0
                u111=0.d0 ; v111=0.d0
                do ii=max(i1+lowfil,ibyz_c(1,i2,i3)),min(i1+1+lupfil,ibyz_c(2,i2,i3))
                    t111=t111 + x_c(ii,i2,i3)*a(ii-i1)
                    s111=s111 + x_fc(ii,i2,i3,1)*b(ii-i1)
                    u111=u111 + x_c(ii,i2,i3)*a(ii-i1-1)
                    v111=v111 + x_fc(ii,i2,i3,1)*b(ii-i1-1)
                enddo
                y_c(i1  ,i2,i3)=t111+s111+cprecr*x_c(i1  ,i2,i3)
                y_c(i1+1,i2,i3)=u111+v111+cprecr*x_c(i1+1,i2,i3)
        enddo
        do i1=i1,ibyz_c(2,i2,i3) ! remainder of unrolled loop
                t111=0.d0 ; s111=0.d0
                do ii=max(i1+lowfil,ibyz_c(1,i2,i3)),min(i1+lupfil,ibyz_c(2,i2,i3))
                    t111=t111 + x_c(ii,i2,i3)*a(ii-i1)
                    s111=s111 + x_fc(ii,i2,i3,1)*b(ii-i1)
                enddo
                y_c(i1  ,i2,i3)=t111+s111+cprecr*x_c(i1  ,i2,i3)
        enddo

        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-3,4
                t211=0.d0 ; s211=0.d0 ; u211=0.d0 ; v211=0.d0 
                do ii=max(i1+lowfil,ibyz_c(1,i2,i3)),min(i1+3+lupfil,ibyz_c(2,i2,i3))
                    t211=t211 + x_c(ii,i2,i3)*c(ii-i1)
                    s211=s211 + x_c(ii,i2,i3)*c(ii-i1-1)
                    u211=u211 + x_c(ii,i2,i3)*c(ii-i1-2)
                    v211=v211 + x_c(ii,i2,i3)*c(ii-i1-3)
                enddo
                y_f(1,i1  ,i2,i3)=t211
                y_f(1,i1+1,i2,i3)=s211
                y_f(1,i1+2,i2,i3)=u211
                y_f(1,i1+3,i2,i3)=v211
        enddo
        do i1=i1,ibyz_f(2,i2,i3) ! remainder of unrolled loop
                t211=0.d0 
                do ii=max(i1+lowfil,ibyz_c(1,i2,i3)),min(i1+lupfil,ibyz_c(2,i2,i3))
                    t211=t211 + x_c(ii,i2,i3)*c(ii-i1)
                enddo
                y_f(1,i1,i2,i3)=t211
        enddo
    enddo
    enddo

!       call system_clock(ncount1,ncount_rate,ncount_max)
!       tel=dble(ncount1-ncount0)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'FIRST PART:x',tel,1.d-6*mflop1/tel

! + (1/2) d^2/dy^2
    nb=16
    do i3=0,n3
    do ii2=0,n2,nb
    do i1=0,n1
        if (ii2.ge.ibxz_c(1,i1,i3) .and. ii2+(nb-1).le.ibxz_c(2,i1,i3)) then
        do i2=ii2,ii2+(nb-1)-1,2 ! middle
                t111=0.d0 ; s111=0.d0
                u111=0.d0 ; v111=0.d0
                do ii=max(i2+lowfil,ibxz_c(1,i1,i3)),min(i2+1+lupfil,ibxz_c(2,i1,i3))
                    t111=t111 + x_c(i1,ii,i3)*a(ii-i2)
                    s111=s111 + x_fc(i1,ii,i3,2)*b(ii-i2)
                    u111=u111 + x_c(i1,ii,i3)*a(ii-i2-1)
                    v111=v111 + x_fc(i1,ii,i3,2)*b(ii-i2-1)
                enddo
                y_c(i1,i2  ,i3)=y_c(i1,i2  ,i3)+t111+s111
                y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+u111+v111
        enddo
        else
        do i2=max(ii2,ibxz_c(1,i1,i3)),min(ii2+(nb-1),ibxz_c(2,i1,i3)) ! head and tail 
                t111=0.d0 ; s111=0.d0
                do ii=max(i2+lowfil,ibxz_c(1,i1,i3)),min(i2+lupfil,ibxz_c(2,i1,i3))
                    t111=t111 + x_c(i1,ii,i3)*a(ii-i2)
                    s111=s111 + x_fc(i1,ii,i3,2)*b(ii-i2)
                enddo
                y_c(i1,i2  ,i3)=y_c(i1,i2  ,i3)+t111+s111
        enddo
        endif

        if (ii2.ge.ibxz_f(1,i1,i3) .and. ii2+(nb-1).le.ibxz_f(2,i1,i3)) then 
        do i2=ii2,ii2+(nb-1)-3,4 ! middle
                t121=0.d0 ; s121=0.d0 ; u121=0.d0 ; v121=0.d0 
                do ii=max(i2+lowfil,ibxz_c(1,i1,i3)),min(i2+3+lupfil,ibxz_c(2,i1,i3))
                    t121=t121 + x_c(i1,ii,i3)*c(ii-i2)
                    s121=s121 + x_c(i1,ii,i3)*c(ii-i2-1)
                    u121=u121 + x_c(i1,ii,i3)*c(ii-i2-2)
                    v121=v121 + x_c(i1,ii,i3)*c(ii-i2-3)
                enddo
                y_f(2,i1,i2  ,i3)=t121
                y_f(2,i1,i2+1,i3)=s121
                y_f(2,i1,i2+2,i3)=u121
                y_f(2,i1,i2+3,i3)=v121
        enddo
        else
        do i2=max(ii2,ibxz_f(1,i1,i3)),min(ii2+(nb-1),ibxz_f(2,i1,i3)) ! head and tail 
                t121=0.d0 
                do ii=max(i2+lowfil,ibxz_c(1,i1,i3)),min(i2+lupfil,ibxz_c(2,i1,i3))
                    t121=t121 + x_c(i1,ii,i3)*c(ii-i2)
                enddo
                y_f(2,i1,i2,i3)=t121
        enddo
        endif
    enddo
    enddo
    enddo


!       call system_clock(ncount2,ncount_rate,ncount_max)
!       tel=dble(ncount2-ncount1)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'FIRST PART:y',tel,1.d-6*mflop2/tel

! + (1/2) d^2/dz^2
    nb=16
    do i2=0,n2
    do ii3=0,n3,nb
    do i1=0,n1
        if (ii3.ge.ibxy_c(1,i1,i2) .and. ii3+(nb-1).le.ibxy_c(2,i1,i2)) then
        do i3=ii3,ii3+(nb-1)-1,2 ! middle
                t111=0.d0 ; s111=0.d0
                u111=0.d0 ; v111=0.d0
                do ii=max(i3+lowfil,ibxy_c(1,i1,i2)),min(i3+1+lupfil,ibxy_c(2,i1,i2))
                    t111=t111 + x_c(i1,i2,ii)*a(ii-i3)
                    s111=s111 + x_fc(i1,i2,ii,3)*b(ii-i3)
                    u111=u111 + x_c(i1,i2,ii)*a(ii-i3-1)
                    v111=v111 + x_fc(i1,i2,ii,3)*b(ii-i3-1)
                enddo
                y_c(i1,i2,i3  )=y_c(i1,i2,i3  )+t111+s111
                y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+u111+v111
        enddo
        else
        do i3=max(ii3,ibxy_c(1,i1,i2)),min(ii3+(nb-1),ibxy_c(2,i1,i2)) ! head and tail 
                t111=0.d0 ; s111=0.d0
                do ii=max(i3+lowfil,ibxy_c(1,i1,i2)),min(i3+lupfil,ibxy_c(2,i1,i2))
                    t111=t111 + x_c(i1,i2,ii)*a(ii-i3)
                    s111=s111 + x_fc(i1,i2,ii,3)*b(ii-i3)
                enddo
                y_c(i1,i2,i3  )=y_c(i1,i2,i3)+t111+s111
        enddo
        endif

        if (ii3.ge.ibxy_f(1,i1,i2) .and. ii3+(nb-1).le.ibxy_f(2,i1,i2)) then
        do i3=ii3,ii3+(nb-1)-3,4 ! middle
                t112=0.d0 ; s112=0.d0 ; u112=0.d0 ; v112=0.d0 
                do ii=max(i3+lowfil,ibxy_c(1,i1,i2)),min(i3+3+lupfil,ibxy_c(2,i1,i2))
                    t112=t112 + x_c(i1,i2,ii)*c(ii-i3)
                    s112=s112 + x_c(i1,i2,ii)*c(ii-i3-1)
                    u112=u112 + x_c(i1,i2,ii)*c(ii-i3-2)
                    v112=v112 + x_c(i1,i2,ii)*c(ii-i3-3)
                enddo
                y_f(4,i1,i2,i3  )=t112
                y_f(4,i1,i2,i3+1)=s112
                y_f(4,i1,i2,i3+2)=u112
                y_f(4,i1,i2,i3+3)=v112
        enddo
        else
        do i3=max(ii3,ibxy_f(1,i1,i2)),min(ii3+(nb-1),ibxy_f(2,i1,i2)) ! head and tail 
                t112=0.d0 
                do ii=max(i3+lowfil,ibxy_c(1,i1,i2)),min(i3+lupfil,ibxy_c(2,i1,i2))
                    t112=t112 + x_c(i1,i2,ii)*c(ii-i3)
                enddo
                y_f(4,i1,i2,i3)=t112
        enddo
        endif
    enddo
    enddo
    enddo

!       call system_clock(ncount3,ncount_rate,ncount_max)
!       tel=dble(ncount3-ncount2)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'FIRST PART:z',tel,1.d-6*mflop3/tel

! wavelet part
 ! (1/2) d^2/dx^2
    do i3=nfl3,nfu3
        do i2=nfl2,nfu2
            do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
                t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
                do l=max(ibyz_f(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_f(2,i2,i3)-i1)
                    t112=t112 + x_f(4,i1+l,i2,i3)*a(l) + x_f(5,i1+l,i2,i3)*b(l)
                    t121=t121 + x_f(2,i1+l,i2,i3)*a(l) + x_f(3,i1+l,i2,i3)*b(l)
                    t122=t122 + x_f(6,i1+l,i2,i3)*a(l) + x_f(7,i1+l,i2,i3)*b(l)
                    t212=t212 + x_f(4,i1+l,i2,i3)*c(l) + x_f(5,i1+l,i2,i3)*e(l)
                    t221=t221 + x_f(2,i1+l,i2,i3)*c(l) + x_f(3,i1+l,i2,i3)*e(l)
                    t222=t222 + x_f(6,i1+l,i2,i3)*c(l) + x_f(7,i1+l,i2,i3)*e(l)
                    t211=t211 + x_f(1,i1+l,i2,i3)*e(l)
                enddo

                y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112+cprecr*x_f(4,i1,i2,i3)
                y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121+cprecr*x_f(2,i1,i2,i3)
                y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211+cprecr*x_f(1,i1,i2,i3)
                y_f(6,i1,i2,i3)=t122+cprecr*x_f(6,i1,i2,i3)
                y_f(5,i1,i2,i3)=t212+cprecr*x_f(5,i1,i2,i3)
                y_f(3,i1,i2,i3)=t221+cprecr*x_f(3,i1,i2,i3)
                y_f(7,i1,i2,i3)=t222+cprecr*x_f(7,i1,i2,i3)
            enddo
        enddo
    enddo

!       call system_clock(ncount4,ncount_rate,ncount_max)
!       tel=dble(ncount4-ncount3)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'SECND PART:x',tel,1.d-6*nflop1/tel


 ! + (1/2) d^2/dy^2
    nb=8
    do i3=nfl3,nfu3
    do ii2=nfl2,nfu2,nb
    do i1=nfl1,nfu1
       do i2=max(ii2,ibxz_f(1,i1,i3)),min(ii2+(nb-1),ibxz_f(2,i1,i3))
                t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
                do l=max(ibxz_f(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_f(2,i1,i3)-i2)
                    t112=t112 + x_f(4,i1,i2+l,i3)*a(l) + x_f(6,i1,i2+l,i3)*b(l)
                    t211=t211 + x_f(1,i1,i2+l,i3)*a(l) + x_f(3,i1,i2+l,i3)*b(l)
                    t122=t122 + x_f(4,i1,i2+l,i3)*c(l) + x_f(6,i1,i2+l,i3)*e(l)
                    t212=t212 + x_f(5,i1,i2+l,i3)*a(l) + x_f(7,i1,i2+l,i3)*b(l)
                    t221=t221 + x_f(1,i1,i2+l,i3)*c(l) + x_f(3,i1,i2+l,i3)*e(l)
                    t222=t222 + x_f(5,i1,i2+l,i3)*c(l) + x_f(7,i1,i2+l,i3)*e(l)
                    t121=t121 + x_f(2,i1,i2+l,i3)*e(l)
                enddo

                y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
                y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
                y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
                y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
                y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
                y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
                y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222
       enddo
    enddo
    enddo
    enddo

!       call system_clock(ncount5,ncount_rate,ncount_max)
!       tel=dble(ncount5-ncount4)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'SECND PART:y',tel,1.d-6*nflop2/tel

 ! + (1/2) d^2/dz^2
    nb=8
    do i2=nfl2,nfu2
    do ii3=nfl3,nfu3,nb
    do i1=nfl1,nfu1
       do i3=max(ii3,ibxy_f(1,i1,i2)),min(ii3+(nb-1),ibxy_f(2,i1,i2))
                t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
                do l=max(ibxy_f(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_f(2,i1,i2)-i3)
                    t121=t121 + x_f(2,i1,i2,i3+l)*a(l) + x_f(6,i1,i2,i3+l)*b(l)
                    t211=t211 + x_f(1,i1,i2,i3+l)*a(l) + x_f(5,i1,i2,i3+l)*b(l)
                    t122=t122 + x_f(2,i1,i2,i3+l)*c(l) + x_f(6,i1,i2,i3+l)*e(l)
                    t212=t212 + x_f(1,i1,i2,i3+l)*c(l) + x_f(5,i1,i2,i3+l)*e(l)
                    t221=t221 + x_f(3,i1,i2,i3+l)*a(l) + x_f(7,i1,i2,i3+l)*b(l)
                    t222=t222 + x_f(3,i1,i2,i3+l)*c(l) + x_f(7,i1,i2,i3+l)*e(l)
                    t112=t112 + x_f(4,i1,i2,i3+l)*e(l)
                enddo

                y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+t112
                y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
                y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211
                y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+t122
                y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+t212
                y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+t221
                y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+t222

       enddo
    enddo
    enddo
    enddo

!       call system_clock(ncount6,ncount_rate,ncount_max)
!       tel=dble(ncount6-ncount5)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'SECND PART:z',tel,1.d-6*nflop3/tel
!
!       tel=dble(ncount6-ncount0)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'ALL   PART',  & 
!            tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

    return
    end



         subroutine ConvolkineticP(n1,n2,n3, &
               nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
               hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x,y,ekin)
!   y = y + (kinetic energy operator)x 
    implicit real(kind=8) (a-h,o-z)
    logical :: firstcall=.true. 
    integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
    dimension x(0:n1,2,0:n2,2,0:n3,2),y(0:n1,2,0:n2,2,0:n3,2)
    dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
    dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)


    parameter(lowfil=-14,lupfil=14)
    dimension a(lowfil-1:lupfil+1),b(lowfil-1:lupfil+1),c(lowfil-3:lupfil+3),e(lowfil:lupfil)
    scale=-.5d0/hgrid**2

!---------------------------------------------------------------------------
! second derivative filters for Daubechies 16
!  <phi|D^2|phi_i>
    a(0)=   -3.5536922899131901941296809374d0*scale
    a(1)=    2.2191465938911163898794546405d0*scale
    a(2)=   -0.6156141465570069496314853949d0*scale
    a(3)=    0.2371780582153805636239247476d0*scale
    a(4)=   -0.0822663999742123340987663521d0*scale
    a(5)=    0.02207029188482255523789911295638968409d0*scale
    a(6)=   -0.409765689342633823899327051188315485d-2*scale
    a(7)=    0.45167920287502235349480037639758496d-3*scale
    a(8)=   -0.2398228524507599670405555359023135d-4*scale
    a(9)=    2.0904234952920365957922889447361d-6*scale
    a(10)=  -3.7230763047369275848791496973044d-7*scale
    a(11)=  -1.05857055496741470373494132287d-8*scale
    a(12)=  -5.813879830282540547959250667d-11*scale
    a(13)=   2.70800493626319438269856689037647576d-13*scale
    a(14)=  -6.924474940639200152025730585882d-18*scale
    do i=1,14
        a(-i)=a(i)
    enddo
    a( 15)=0.d0
    b( 15)=0.d0
    a(-15)=0.d0
    b(-15)=0.d0
!  <phi|D^2|psi_i>
    c(-14)=     -3.869102413147656535541850057188d-18*scale
    c(-13)=      1.5130616560866154733900029272077362d-13*scale
    c(-12)=     -3.2264702314010525539061647271983988409d-11*scale
    c(-11)=     -5.96264938781402337319841002642d-9*scale
    c(-10)=     -2.1656830629214041470164889350342d-7*scale
    c(-9 )=      8.7969704055286288323596890609625d-7*scale
    c(-8 )=     -0.00001133456724516819987751818232711775d0*scale
    c(-7 )=      0.00021710795484646138591610188464622454d0*scale
    c(-6 )=     -0.0021356291838797986414312219042358542d0*scale
    c(-5 )=      0.00713761218453631422925717625758502986d0*scale
    c(-4 )=     -0.0284696165863973422636410524436931061d0*scale
    c(-3 )=      0.14327329352510759457155821037742893841d0*scale
    c(-2 )=     -0.42498050943780130143385739554118569733d0*scale
    c(-1 )=      0.65703074007121357894896358254040272157d0*scale
    c( 0 )=     -0.42081655293724308770919536332797729898d0*scale
    c( 1 )=     -0.21716117505137104371463587747283267899d0*scale
    c( 2 )=      0.63457035267892488185929915286969303251d0*scale
    c( 3 )=     -0.53298223962800395684936080758073568406d0*scale
    c( 4 )=      0.23370490631751294307619384973520033236d0*scale
    c( 5 )=     -0.05657736973328755112051544344507997075d0*scale
    c( 6 )=      0.0080872029411844780634067667008050127d0*scale
    c( 7 )=     -0.00093423623304808664741804536808932984d0*scale
    c( 8 )=      0.00005075807947289728306309081261461095d0*scale
    c( 9 )=     -4.62561497463184262755416490048242d-6*scale
    c( 10)=      6.3919128513793415587294752371778d-7*scale
    c( 11)=      1.87909235155149902916133888931d-8*scale
    c( 12)=      1.04757345962781829480207861447155543883d-10*scale
    c( 13)=     -4.84665690596158959648731537084025836d-13*scale
    c( 14)=      1.2392629629188986192855777620877d-17*scale
    c(-15)=0.d0
    c(-16)=0.d0
    c(-17)=0.d0
    c( 15)=0.d0
    c( 16)=0.d0
    c( 17)=0.d0
!  <psi|D^2|phi_i>
    do i=-14,14
        b(i)=c(-i)
    enddo
    !<psi|D^2|psi_i>
    e(0)=   -24.875846029392331358907766562d0*scale
    e(1)=   -7.1440597663471719869313377994d0*scale
    e(2)=   -0.04251705323669172315864542163525830944d0*scale
    e(3)=   -0.26995931336279126953587091167128839196d0*scale
    e(4)=    0.08207454169225172612513390763444496516d0*scale
    e(5)=   -0.02207327034586634477996701627614752761d0*scale
    e(6)=    0.00409765642831595181639002667514310145d0*scale
    e(7)=   -0.00045167920287507774929432548999880117d0*scale
    e(8)=    0.00002398228524507599670405555359023135d0*scale
    e(9)=   -2.0904234952920365957922889447361d-6*scale
    e(10)=   3.7230763047369275848791496973044d-7*scale
    e(11)=   1.05857055496741470373494132287d-8*scale
    e(12)=   5.8138798302825405479592506674648873655d-11*scale
    e(13)=  -2.70800493626319438269856689037647576d-13*scale
    e(14)=   6.924474940639200152025730585882d-18*scale
    do i=1,14
        e(-i)=e(i)
    enddo


!  if (firstcall) then
!
!! (1/2) d^2/dx^2
!    mflop1=0
!    do i3=0,n3
!    do i2=0,n2
!        do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)
!                do l=max(ibyz_c(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_c(2,i2,i3)-i1)
!                    mflop1=mflop1+4
!                enddo
!                    mflop1=mflop1+4
!        enddo
!
!        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
!                do l=max(ibyz_c(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_c(2,i2,i3)-i1)
!                    mflop1=mflop1+3
!                enddo
!                    mflop1=mflop1+3
!        enddo
!    enddo
!    enddo
!
!! + (1/2) d^2/dy^2
!    mflop2=0
!    do i3=0,n3
!    do i1=0,n1
!        do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)
!                do l=max(ibxz_c(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_c(2,i1,i3)-i2)
!                    mflop2=mflop2+4
!                enddo
!                    mflop2=mflop2+4
!        enddo
!
!        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
!                do l=max(ibxz_c(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_c(2,i1,i3)-i2)
!                    mflop2=mflop2+3
!                enddo
!                    mflop2=mflop2+3
!        enddo
!    enddo
!    enddo
!
!
!! + (1/2) d^2/dz^2
!    mflop3=0
!    do i2=0,n2
!    do i1=0,n1
!        do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)
!                do l=max(ibxy_c(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_c(2,i1,i2)-i3)
!                    mflop3=mflop3+4
!                enddo
!                    mflop3=mflop3+4
!        enddo
!
!        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
!                do l=max(ibxy_c(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_c(2,i1,i2)-i3)
!                    mflop3=mflop3+3
!                enddo
!                    mflop3=mflop3+3
!        enddo
!    enddo
!    enddo
!
!! wavelet part
! ! (1/2) d^2/dx^2
!    nflop1=0
!    do i3=nfl3,nfu3
!        do i2=nfl2,nfu2
!            do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
!                do l=max(ibyz_f(1,i2,i3)-i1,lowfil),min(ibyz_f(2,i2,i3),nfu1-i1)
!                    nflop1=nflop1+26
!                enddo
!                    nflop1=nflop1+21
!            enddo
!        enddo
!    enddo
!
! ! + (1/2) d^2/dy^2
!    nflop2=0
!    do i3=nfl3,nfu3
!    do i1=nfl1,nfu1
!       do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
!                do l=max(ibxz_f(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_f(2,i1,i3)-i2)
!                    nflop2=nflop2+26
!                enddo
!                    nflop2=nflop2+21
!       enddo
!    enddo
!    enddo
!
! ! + (1/2) d^2/dz^2
!    nflop3=0
!    do i2=nfl2,nfu2
!    do i1=nfl1,nfu1
!       do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
!                do l=max(ibxy_f(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_f(2,i1,i2)-i3)
!                    nflop3=nflop3+26
!                enddo
!                    nflop3=nflop3+21
!       enddo
!    enddo
!    enddo
!    
!    firstcall=.false.
!    endif


!---------------------------------------------------------------------------

     ekin1=0.d0 ; ekin2=0.d0 ; ekin3=0.d0 ; ekin4=0.d0 ; ekin5=0.d0 ; ekin6=0.d0 ; ekin7=0.d0 
                                                            
! Scaling function part

!       call system_clock(ncount0,ncount_rate,ncount_max)


! (1/2) d^2/dx^2
    do i3=0,n3
    do i2=0,n2
        do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-1,2
                t111=0.d0 ; s111=0.d0
                u111=0.d0 ; v111=0.d0
                do ii=max(i1+lowfil,ibyz_c(1,i2,i3)),min(i1+1+lupfil,ibyz_c(2,i2,i3))
                    t111=t111 + x(ii,1,i2,1,i3,1)*a(ii-i1)
                    s111=s111 + x(ii,2,i2,1,i3,1)*b(ii-i1)
                    u111=u111 + x(ii,1,i2,1,i3,1)*a(ii-i1-1)
                    v111=v111 + x(ii,2,i2,1,i3,1)*b(ii-i1-1)
                enddo
                y(i1  ,1,i2,1,i3,1)=y(i1  ,1,i2,1,i3,1)+(t111+s111)
                y(i1+1,1,i2,1,i3,1)=y(i1+1,1,i2,1,i3,1)+(u111+v111)
                ekin1=ekin1+(t111+s111)*x(i1  ,1,i2,1,i3,1)
                ekin2=ekin2+(u111+v111)*x(i1+1,1,i2,1,i3,1)
        enddo
        do i1=i1,ibyz_c(2,i2,i3)  ! remainder of unrolled loop
                t111=0.d0 ; s111=0.d0
                do ii=max(i1+lowfil,ibyz_c(1,i2,i3)),min(i1+lupfil,ibyz_c(2,i2,i3))
                    t111=t111 + x(ii,1,i2,1,i3,1)*a(ii-i1)
                    s111=s111 + x(ii,2,i2,1,i3,1)*b(ii-i1)
                enddo
                y(i1,1,i2,1,i3,1)=y(i1,1,i2,1,i3,1)+(t111+s111)
                ekin3=ekin3+(t111+s111)*x(i1,1,i2,1,i3,1)
        enddo

        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-3,4
                t211=0.d0 ; s211=0.d0 ; u211=0.d0 ; v211=0.d0 
                do ii=max(i1+lowfil,ibyz_c(1,i2,i3)),min(i1+3+lupfil,ibyz_c(2,i2,i3))
                    t211=t211 + x(ii,1,i2,1,i3,1)*c(ii-i1)
                    s211=s211 + x(ii,1,i2,1,i3,1)*c(ii-i1-1)
                    u211=u211 + x(ii,1,i2,1,i3,1)*c(ii-i1-2)
                    v211=v211 + x(ii,1,i2,1,i3,1)*c(ii-i1-3)
                enddo
                y(i1  ,2,i2,1,i3,1)=y(i1  ,2,i2,1,i3,1)+t211
                y(i1+1,2,i2,1,i3,1)=y(i1+1,2,i2,1,i3,1)+s211
                y(i1+2,2,i2,1,i3,1)=y(i1+2,2,i2,1,i3,1)+u211
                y(i1+3,2,i2,1,i3,1)=y(i1+3,2,i2,1,i3,1)+v211
                ekin1=ekin1+t211*x(i1  ,2,i2,1,i3,1)
                ekin2=ekin2+s211*x(i1+1,2,i2,1,i3,1)
                ekin3=ekin3+u211*x(i1+2,2,i2,1,i3,1)
                ekin4=ekin4+v211*x(i1+3,2,i2,1,i3,1)
        enddo
        do i1=i1,ibyz_f(2,i2,i3) ! remainder of unrolled loop
                t211=0.d0 
                do ii=max(i1+lowfil,ibyz_c(1,i2,i3)),min(i1+lupfil,ibyz_c(2,i2,i3))
                    t211=t211 + x(ii,1,i2,1,i3,1)*c(ii-i1)
                enddo
                y(i1,2,i2,1,i3,1)=y(i1,2,i2,1,i3,1)+t211
                ekin5=ekin5+t211*x(i1,2,i2,1,i3,1)
        enddo
    enddo
    enddo

!       call system_clock(ncount1,ncount_rate,ncount_max)
!       tel=dble(ncount1-ncount0)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'P:FIRST PART:x',tel,1.d-6*mflop1/tel

! + (1/2) d^2/dy^2
    nb=16
    do i3=0,n3
    do ii2=0,n2,nb
    do i1=0,n1
        if (ii2.ge.ibxz_c(1,i1,i3) .and. ii2+(nb-1).le.ibxz_c(2,i1,i3)) then
        do i2=ii2,ii2+(nb-1)-1,2 ! middle
                t111=0.d0 ; s111=0.d0
                u111=0.d0 ; v111=0.d0
                do ii=max(i2+lowfil,ibxz_c(1,i1,i3)),min(i2+1+lupfil,ibxz_c(2,i1,i3))
                    t111=t111 + x(i1,1,ii,1,i3,1)*a(ii-i2)
                    s111=s111 + x(i1,1,ii,2,i3,1)*b(ii-i2)
                    u111=u111 + x(i1,1,ii,1,i3,1)*a(ii-i2-1)
                    v111=v111 + x(i1,1,ii,2,i3,1)*b(ii-i2-1)
                enddo
                y(i1,1,i2  ,1,i3,1)=y(i1,1,i2  ,1,i3,1)+(t111+s111)
                y(i1,1,i2+1,1,i3,1)=y(i1,1,i2+1,1,i3,1)+(u111+v111)
                ekin1=ekin1+(t111+s111)*x(i1,1,i2  ,1,i3,1)
                ekin2=ekin2+(u111+v111)*x(i1,1,i2+1,1,i3,1)
        enddo
        else
        do i2=max(ii2,ibxz_c(1,i1,i3)),min(ii2+(nb-1),ibxz_c(2,i1,i3)) ! head and tail 
                t111=0.d0 ; s111=0.d0
                do ii=max(i2+lowfil,ibxz_c(1,i1,i3)),min(i2+lupfil,ibxz_c(2,i1,i3))
                    t111=t111 + x(i1,1,ii,1,i3,1)*a(ii-i2)
                    s111=s111 + x(i1,1,ii,2,i3,1)*b(ii-i2)
                enddo
                y(i1,1,i2,1,i3,1)=y(i1,1,i2,1,i3,1)+(t111+s111)
                ekin3=ekin3+(t111+s111)*x(i1,1,i2,1,i3,1)
        enddo
        endif

        if (ii2.ge.ibxz_f(1,i1,i3) .and. ii2+(nb-1).le.ibxz_f(2,i1,i3)) then 
        do i2=ii2,ii2+(nb-1)-3,4 ! middle
                t121=0.d0 ; s121=0.d0 ; u121=0.d0 ; v121=0.d0 
                do ii=max(i2+lowfil,ibxz_c(1,i1,i3)),min(i2+3+lupfil,ibxz_c(2,i1,i3))
                    t121=t121 + x(i1,1,ii,1,i3,1)*c(ii-i2)
                    s121=s121 + x(i1,1,ii,1,i3,1)*c(ii-i2-1)
                    u121=u121 + x(i1,1,ii,1,i3,1)*c(ii-i2-2)
                    v121=v121 + x(i1,1,ii,1,i3,1)*c(ii-i2-3)
                enddo
                y(i1,1,i2  ,2,i3,1)=y(i1,1,i2  ,2,i3,1)+t121
                y(i1,1,i2+1,2,i3,1)=y(i1,1,i2+1,2,i3,1)+s121
                y(i1,1,i2+2,2,i3,1)=y(i1,1,i2+2,2,i3,1)+u121
                y(i1,1,i2+3,2,i3,1)=y(i1,1,i2+3,2,i3,1)+v121
                ekin1=ekin1+t121*x(i1,1,i2  ,2,i3,1)
                ekin2=ekin2+s121*x(i1,1,i2+1,2,i3,1)
                ekin3=ekin3+u121*x(i1,1,i2+2,2,i3,1)
                ekin4=ekin4+v121*x(i1,1,i2+3,2,i3,1)
        enddo
        else
        do i2=max(ii2,ibxz_f(1,i1,i3)),min(ii2+(nb-1),ibxz_f(2,i1,i3)) ! head and tail 
                t121=0.d0 
                do ii=max(i2+lowfil,ibxz_c(1,i1,i3)),min(i2+lupfil,ibxz_c(2,i1,i3))
                    t121=t121 + x(i1,1,ii,1,i3,1)*c(ii-i2)
                enddo
                y(i1,1,i2,2,i3,1)=y(i1,1,i2,2,i3,1)+t121
                ekin5=ekin5+t121*x(i1,1,i2,2,i3,1)
        enddo
        endif
    enddo
    enddo
    enddo

!       call system_clock(ncount2,ncount_rate,ncount_max)
!       tel=dble(ncount2-ncount1)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'P:FIRST PART:y',tel,1.d-6*mflop2/tel

! + (1/2) d^2/dz^2
    nb=16
    do i2=0,n2
    do ii3=0,n3,nb
    do i1=0,n1
        if (ii3.ge.ibxy_c(1,i1,i2) .and. ii3+(nb-1).le.ibxy_c(2,i1,i2)) then
        do i3=ii3,ii3+(nb-1)-1,2 ! middle
                t111=0.d0 ; s111=0.d0
                u111=0.d0 ; v111=0.d0
                do ii=max(i3+lowfil,ibxy_c(1,i1,i2)),min(i3+1+lupfil,ibxy_c(2,i1,i2))
                    t111=t111 + x(i1,1,i2,1,ii,1)*a(ii-i3)
                    s111=s111 + x(i1,1,i2,1,ii,2)*b(ii-i3)
                    u111=u111 + x(i1,1,i2,1,ii,1)*a(ii-i3-1)
                    v111=v111 + x(i1,1,i2,1,ii,2)*b(ii-i3-1)
                enddo
                y(i1,1,i2,1,i3  ,1)=y(i1,1,i2,1,i3  ,1)+(t111+s111)
                y(i1,1,i2,1,i3+1,1)=y(i1,1,i2,1,i3+1,1)+(u111+v111)
                ekin1=ekin1+(t111+s111)*x(i1,1,i2,1,i3  ,1)
                ekin2=ekin2+(u111+v111)*x(i1,1,i2,1,i3+1,1)
        enddo
        else
        do i3=max(ii3,ibxy_c(1,i1,i2)),min(ii3+(nb-1),ibxy_c(2,i1,i2)) ! head and tail 
                t111=0.d0 ; s111=0.d0
                do ii=max(i3+lowfil,ibxy_c(1,i1,i2)),min(i3+lupfil,ibxy_c(2,i1,i2))
                    t111=t111 + x(i1,1,i2,1,ii,1)*a(ii-i3)
                    s111=s111 + x(i1,1,i2,1,ii,2)*b(ii-i3)
                enddo
                y(i1,1,i2,1,i3,1)=y(i1,1,i2,1,i3,1)+(t111+s111)
                ekin3=ekin3+(t111+s111)*x(i1,1,i2,1,i3,1)
        enddo
        endif

        if (ii3.ge.ibxy_f(1,i1,i2) .and. ii3+(nb-1).le.ibxy_f(2,i1,i2)) then
        do i3=ii3,ii3+(nb-1)-3,4 ! middle
                t112=0.d0 ; s112=0.d0 ; u112=0.d0 ; v112=0.d0 
                do ii=max(i3+lowfil,ibxy_c(1,i1,i2)),min(i3+3+lupfil,ibxy_c(2,i1,i2))
                    t112=t112 + x(i1,1,i2,1,ii,1)*c(ii-i3)
                    s112=s112 + x(i1,1,i2,1,ii,1)*c(ii-i3-1)
                    u112=u112 + x(i1,1,i2,1,ii,1)*c(ii-i3-2)
                    v112=v112 + x(i1,1,i2,1,ii,1)*c(ii-i3-3)
                enddo
                y(i1,1,i2,1,i3  ,2)=y(i1,1,i2,1,i3  ,2)+t112
                y(i1,1,i2,1,i3+1,2)=y(i1,1,i2,1,i3+1,2)+s112
                y(i1,1,i2,1,i3+2,2)=y(i1,1,i2,1,i3+2,2)+u112
                y(i1,1,i2,1,i3+3,2)=y(i1,1,i2,1,i3+3,2)+v112
                ekin1=ekin1+t112*x(i1,1,i2,1,i3  ,2)
                ekin2=ekin2+s112*x(i1,1,i2,1,i3+1,2)
                ekin3=ekin3+u112*x(i1,1,i2,1,i3+2,2)
                ekin4=ekin4+v112*x(i1,1,i2,1,i3+3,2)
        enddo
        else
        do i3=max(ii3,ibxy_f(1,i1,i2)),min(ii3+(nb-1),ibxy_f(2,i1,i2)) ! head and tail 
                t112=0.d0 
                do ii=max(i3+lowfil,ibxy_c(1,i1,i2)),min(i3+lupfil,ibxy_c(2,i1,i2))
                    t112=t112 + x(i1,1,i2,1,ii,1)*c(ii-i3)
                enddo
                y(i1,1,i2,1,i3,2)=y(i1,1,i2,1,i3,2)+t112
                ekin5=ekin5+t112*x(i1,1,i2,1,i3,2)
        enddo
        endif
    enddo
    enddo
    enddo

!       call system_clock(ncount3,ncount_rate,ncount_max)
!       tel=dble(ncount3-ncount2)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'P:FIRST PART:z',tel,1.d-6*mflop3/tel

! wavelet part
 ! (1/2) d^2/dx^2
    do i3=nfl3,nfu3
        do i2=nfl2,nfu2
            do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
                t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
                do l=max(ibyz_f(1,i2,i3)-i1,lowfil),min(lupfil,ibyz_f(2,i2,i3)-i1)
                    t112=t112 + x(i1+l,1,i2,1,i3,2)*a(l) + x(i1+l,2,i2,1,i3,2)*b(l)
                    t121=t121 + x(i1+l,1,i2,2,i3,1)*a(l) + x(i1+l,2,i2,2,i3,1)*b(l)
                    t122=t122 + x(i1+l,1,i2,2,i3,2)*a(l) + x(i1+l,2,i2,2,i3,2)*b(l)
                    t212=t212 + x(i1+l,1,i2,1,i3,2)*c(l) + x(i1+l,2,i2,1,i3,2)*e(l)
                    t221=t221 + x(i1+l,1,i2,2,i3,1)*c(l) + x(i1+l,2,i2,2,i3,1)*e(l)
                    t222=t222 + x(i1+l,1,i2,2,i3,2)*c(l) + x(i1+l,2,i2,2,i3,2)*e(l)
                    t211=t211 + x(i1+l,2,i2,1,i3,1)*e(l)
                enddo

                y(i1,1,i2,1,i3,2)=y(i1,1,i2,1,i3,2)+t112
                y(i1,1,i2,2,i3,1)=y(i1,1,i2,2,i3,1)+t121
                y(i1,2,i2,1,i3,1)=y(i1,2,i2,1,i3,1)+t211
                y(i1,1,i2,2,i3,2)=y(i1,1,i2,2,i3,2)+t122
                y(i1,2,i2,1,i3,2)=y(i1,2,i2,1,i3,2)+t212
                y(i1,2,i2,2,i3,1)=y(i1,2,i2,2,i3,1)+t221
                y(i1,2,i2,2,i3,2)=y(i1,2,i2,2,i3,2)+t222
                ekin1=ekin1+x(i1,1,i2,1,i3,2)*t112
                ekin2=ekin2+x(i1,1,i2,2,i3,1)*t121
                ekin3=ekin3+x(i1,2,i2,1,i3,1)*t211
                ekin4=ekin4+x(i1,1,i2,2,i3,2)*t122
                ekin5=ekin5+x(i1,2,i2,1,i3,2)*t212
                ekin6=ekin6+x(i1,2,i2,2,i3,1)*t221
                ekin7=ekin7+x(i1,2,i2,2,i3,2)*t222
            enddo
        enddo
    enddo

!       call system_clock(ncount4,ncount_rate,ncount_max)
!       tel=dble(ncount4-ncount3)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'P:SECND PART:x',tel,1.d-6*nflop1/tel


 ! + (1/2) d^2/dy^2
    nb=8 
    do i3=nfl3,nfu3
    do ii2=nfl2,nfu2,nb
    do i1=nfl1,nfu1
       do i2=max(ii2,ibxz_f(1,i1,i3)),min(ii2+(nb-1),ibxz_f(2,i1,i3))
                t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
                do l=max(ibxz_f(1,i1,i3)-i2,lowfil),min(lupfil,ibxz_f(2,i1,i3)-i2)
                    t112=t112 + x(i1,1,i2+l,1,i3,2)*a(l) + x(i1,1,i2+l,2,i3,2)*b(l)
                    t211=t211 + x(i1,2,i2+l,1,i3,1)*a(l) + x(i1,2,i2+l,2,i3,1)*b(l)
                    t122=t122 + x(i1,1,i2+l,1,i3,2)*c(l) + x(i1,1,i2+l,2,i3,2)*e(l)
                    t212=t212 + x(i1,2,i2+l,1,i3,2)*a(l) + x(i1,2,i2+l,2,i3,2)*b(l)
                    t221=t221 + x(i1,2,i2+l,1,i3,1)*c(l) + x(i1,2,i2+l,2,i3,1)*e(l)
                    t222=t222 + x(i1,2,i2+l,1,i3,2)*c(l) + x(i1,2,i2+l,2,i3,2)*e(l)
                    t121=t121 + x(i1,1,i2+l,2,i3,1)*e(l)
                enddo

                y(i1,1,i2,1,i3,2)=y(i1,1,i2,1,i3,2)+t112
                y(i1,1,i2,2,i3,1)=y(i1,1,i2,2,i3,1)+t121
                y(i1,2,i2,1,i3,1)=y(i1,2,i2,1,i3,1)+t211
                y(i1,1,i2,2,i3,2)=y(i1,1,i2,2,i3,2)+t122
                y(i1,2,i2,1,i3,2)=y(i1,2,i2,1,i3,2)+t212
                y(i1,2,i2,2,i3,1)=y(i1,2,i2,2,i3,1)+t221
                y(i1,2,i2,2,i3,2)=y(i1,2,i2,2,i3,2)+t222
                ekin1=ekin1+x(i1,1,i2,1,i3,2)*t112
                ekin2=ekin2+x(i1,1,i2,2,i3,1)*t121
                ekin3=ekin3+x(i1,2,i2,1,i3,1)*t211
                ekin4=ekin4+x(i1,1,i2,2,i3,2)*t122
                ekin5=ekin5+x(i1,2,i2,1,i3,2)*t212
                ekin6=ekin6+x(i1,2,i2,2,i3,1)*t221
                ekin7=ekin7+x(i1,2,i2,2,i3,2)*t222
       enddo
    enddo
    enddo
    enddo

!       call system_clock(ncount5,ncount_rate,ncount_max)
!       tel=dble(ncount5-ncount4)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'P:SECND PART:y',tel,1.d-6*nflop2/tel

 ! + (1/2) d^2/dz^2
    nb=8 
    do i2=nfl2,nfu2
    do ii3=nfl3,nfu3,nb
    do i1=nfl1,nfu1
       do i3=max(ii3,ibxy_f(1,i1,i2)),min(ii3+(nb-1),ibxy_f(2,i1,i2))
                t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0 
                do l=max(ibxy_f(1,i1,i2)-i3,lowfil),min(lupfil,ibxy_f(2,i1,i2)-i3)
                    t121=t121 + x(i1,1,i2,2,i3+l,1)*a(l) + x(i1,1,i2,2,i3+l,2)*b(l)
                    t211=t211 + x(i1,2,i2,1,i3+l,1)*a(l) + x(i1,2,i2,1,i3+l,2)*b(l)
                    t122=t122 + x(i1,1,i2,2,i3+l,1)*c(l) + x(i1,1,i2,2,i3+l,2)*e(l)
                    t212=t212 + x(i1,2,i2,1,i3+l,1)*c(l) + x(i1,2,i2,1,i3+l,2)*e(l)
                    t221=t221 + x(i1,2,i2,2,i3+l,1)*a(l) + x(i1,2,i2,2,i3+l,2)*b(l)
                    t222=t222 + x(i1,2,i2,2,i3+l,1)*c(l) + x(i1,2,i2,2,i3+l,2)*e(l)
                    t112=t112 + x(i1,1,i2,1,i3+l,2)*e(l)
                enddo

                y(i1,1,i2,1,i3,2)=y(i1,1,i2,1,i3,2)+t112
                y(i1,1,i2,2,i3,1)=y(i1,1,i2,2,i3,1)+t121
                y(i1,2,i2,1,i3,1)=y(i1,2,i2,1,i3,1)+t211
                y(i1,1,i2,2,i3,2)=y(i1,1,i2,2,i3,2)+t122
                y(i1,2,i2,1,i3,2)=y(i1,2,i2,1,i3,2)+t212
                y(i1,2,i2,2,i3,1)=y(i1,2,i2,2,i3,1)+t221
                y(i1,2,i2,2,i3,2)=y(i1,2,i2,2,i3,2)+t222
                ekin1=ekin1+x(i1,1,i2,1,i3,2)*t112
                ekin2=ekin2+x(i1,1,i2,2,i3,1)*t121
                ekin3=ekin3+x(i1,2,i2,1,i3,1)*t211
                ekin4=ekin4+x(i1,1,i2,2,i3,2)*t122
                ekin5=ekin5+x(i1,2,i2,1,i3,2)*t212
                ekin6=ekin6+x(i1,2,i2,2,i3,1)*t221
                ekin7=ekin7+x(i1,2,i2,2,i3,2)*t222

       enddo
    enddo
    enddo
    enddo

       ekin=ekin1+ekin2+ekin3+ekin4+ekin5+ekin6+ekin7

!       call system_clock(ncount6,ncount_rate,ncount_max)
!       tel=dble(ncount6-ncount5)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'P:SECND PART:z',tel,1.d-6*nflop3/tel
!
!       tel=dble(ncount6-ncount0)/dble(ncount_rate)
!       write(99,'(a40,2(1x,e10.3))') 'P:ALL   PART',  & 
!            tel,1.d-6*(mflop1+mflop2+mflop3+nflop1+nflop2+nflop3)/tel

    return
    end
