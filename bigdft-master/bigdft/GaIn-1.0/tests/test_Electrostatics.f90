program test_Electrostatics
  
  use mod_R_1_norm
  
  implicit none
  
  ! local variables
  integer, parameter      :: nmax = 3
  integer, parameter      :: ng   = 127
  real(kind=8), parameter :: a    = 1.0d16
  real(kind=8), parameter :: b    = 3.5d0
  real(kind=8), parameter :: d    = 2.0d0
  real(kind=8), parameter :: lg   = 7.0d0
  real(kind=8), parameter :: x0   = -12.0d0
  real(kind=8), parameter :: y0   = -13.5d0
  real(kind=8), parameter :: z0   =  15.0d0
  integer, parameter      :: nx2  = 2
  integer, parameter      :: ny2  = 0
  integer, parameter      :: nz2  = 1
  integer, parameter      :: l2   = 2
  integer, parameter      :: m2   =-1
  integer, parameter      :: lmax = 6
  real(kind=8), parameter :: r1(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
  real(kind=8), parameter :: r2(3) = (/ 0.5d0,-0.3d0, 0.2d0 /)
  integer      :: i
  integer      :: j
  integer      :: k
  integer      :: l
  integer      :: m
  integer      :: nx
  integer      :: ny
  integer      :: nz
  real(kind=8) :: Ex
  real(kind=8) :: Ey
  real(kind=8) :: Ez
  real(kind=8) :: x
  real(kind=8) :: y
  real(kind=8) :: z
  real(kind=8) :: tmp
  real(kind=8) :: V
  real(kind=8) :: S1
  real(kind=8) :: S2
  real(kind=8) :: S3
  real(kind=8) :: S1_
  real(kind=8) :: E_(3)
  real(kind=8) :: potential_from_C
  real(kind=8) :: potential_from_Y
  real(kind=8) :: Y_Value
  
  !!
  !! check for C routines
  !!
  do nx=0,nmax
    do ny=0,nmax
      do nz=0,nmax
        
        ! explicit calculation
        S1=0.0d0
        
        do k=0,ng
          z=-lg+(2.0d0*lg*k)/ng
          do j=0,ng
            y=-lg+(2.0d0*lg*j)/ng
            do i=0,ng
              x=-lg+(2.0d0*lg*i)/ng
              V  = 1.0d0/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
              tmp= x**nx * y**ny * z**nz * exp( -b * (x**2+y**2+z**2) )
              S1 = S1 + tmp * V
            end do
          end do
        end do
        S1=S1*8.0d0*lg**3/ng**3
        
        ! get potential
        S1_= potential_from_C((/x0,y0,z0/),b,r1,nx,ny,nz)
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit interaction=',S1,'  potential_from_C=',S1_
        
        ! test result
        if ( abs((S1-S1_)/S1) > 1.0d-6 ) then
          print *,'Error: potential_from_C differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Electrostatics failed'
          print *,''
          stop 1
        end if
        
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Electrostatics passed up to nx=',nmax,' ny=',nmax,' nz=',nmax
  print *,''
  
  !!
  !! check for Y routine
  !!
  do l=0,lmax
    do m=-l,l
      ! explicit calculation
      S1=0.0d0
      do i=0,ng
        x=-lg+(2.0d0*lg*i)/ng
        do j=0,ng
          y=-lg+(2.0d0*lg*j)/ng
          do k=0,ng
            z=-lg+(2.0d0*lg*k)/ng
            V  = 1.0d0/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
            tmp=Y_Value(b,l,m,(/x,y,z/))
            S1 = S1 + tmp * V
          end do
        end do
      end do
      S1=S1*8.0d0*lg**3/ng**3
      
      ! get corresponding interaction
      S1_= potential_from_Y((/x0,y0,z0/),b,r1,l,m)
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  explicit interaction=',S1,'  potential_from_Y=',S1_
        
      ! test result
      if ( abs((S1-S1_)/S1) > 1.0d-5 ) then
        print *,'Error: potential_from_Y differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Electrostatics failed'
        print *,''
        stop 1
      end if
      
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Electrostatics passed up to lmax=',lmax
  print *,''
  
  !!
  !! check for C routines
  !!
  do nx=0,nmax
    do ny=0,nmax
      do nz=0,nmax
        
        ! explicit calculation
        S1=0.0d0
        S2=0.0d0
        S3=0.0d0
        do k=0,ng
          z=-lg+(2.0d0*lg*k)/ng
          do j=0,ng
            y=-lg+(2.0d0*lg*j)/ng
            do i=0,ng
              x=-lg+(2.0d0*lg*i)/ng
              Ex  = -(x-x0)/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)**3
              Ey  = -(y-y0)/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)**3
              Ez  = -(z-z0)/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)**3
              tmp= x**nx * y**ny * z**nz * exp( -b * (x**2+y**2+z**2) )
              S1 = S1 + tmp * Ex
              S2 = S2 + tmp * Ey
              S3 = S3 + tmp * Ez
            end do
          end do
        end do
        S1=S1*8.0d0*lg**3/ng**3
        S2=S2*8.0d0*lg**3/ng**3
        S3=S3*8.0d0*lg**3/ng**3
        
        ! get field
        call field_from_C((/x0,y0,z0/),b,r1,nx,ny,nz,E_)
        
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit interaction=',S1,'  field_from_C(1)=',E_(1)
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit interaction=',S2,'  field_from_C(2)=',E_(2)
        write (*,'(A,I2,A,I2,A,I2,A,E22.16,A,E22.16)')  & 
          'nx=',nx,'  ny=',ny,'  nz=',nz,'  explicit interaction=',S3,'  field_from_C(3)=',E_(3)
        
        ! test result
        if ( abs((S1-E_(1))/S1) > 1.0d-6 ) then
          print *,'Error: field_from_C differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Electrostatics failed'
          print *,''
          stop 1
        end if
        if ( abs((S2-E_(2))/S2) > 1.0d-6 ) then
          print *,'Error: field_from_C differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Electrostatics failed'
          print *,''
          stop 1
        end if
        if ( abs((S3-E_(3))/S3) > 1.0d-6 ) then
          print *,'Error: field_from_C differs from explicit calculation for case'
          print *,'       nx =',nx
          print *,'       ny =',ny
          print *,'       nz =',nz
          print *,''
          print *,'test_Electrostatics failed'
          print *,''
          stop 1
        end if
        
      end do
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Electrostatics passed up to nx=',nmax,' ny=',nmax,' nz=',nmax
  print *,''
  
  !!
  !! check for Y routine
  !!
  do l=0,lmax
    do m=-l,l
      ! explicit calculation
      S1=0.0d0
      S2=0.0d0
      S3=0.0d0
      do i=0,ng
        x=-lg+(2.0d0*lg*i)/ng
        do j=0,ng
          y=-lg+(2.0d0*lg*j)/ng
          do k=0,ng
            z=-lg+(2.0d0*lg*k)/ng
            Ex  = -(x-x0)/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)**3
            Ey  = -(y-y0)/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)**3
            Ez  = -(z-z0)/sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)**3
            tmp=Y_Value(b,l,m,(/x,y,z/))
            S1 = S1 + tmp * Ex
            S2 = S2 + tmp * Ey
            S3 = S3 + tmp * Ez
          end do
        end do
      end do
      S1=S1*8.0d0*lg**3/ng**3
      S2=S2*8.0d0*lg**3/ng**3
      S3=S3*8.0d0*lg**3/ng**3
      
      ! get field
      call field_from_Y((/x0,y0,z0/),b,r1,l,m,E_)
      
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  explicit interaction=',S1,'  field_from_Y(1)=',E_(1)
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  explicit interaction=',S2,'  field_from_Y(2)=',E_(2)
      write (*,'(A,I2,A,I2,A,E22.16,A,E22.16)')  & 
        'l=',l,'  m=',m,'  explicit interaction=',S3,'  field_from_Y(3)=',E_(3)
      
      ! test result
      if ( abs((S1-E_(1))/S1) > 1.0d-5 ) then
        print *,'Error: field_from_Y differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Electrostatics failed'
        print *,''
        stop 1
      end if
      if ( abs((S2-E_(2))/S2) > 1.0d-5 ) then
        print *,'Error: field_from_Y differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Electrostatics failed'
        print *,''
        stop 1
      end if
      if ( abs((S3-E_(3))/S3) > 1.0d-5 ) then
        print *,'Error: field_from_Y differs from explicit calculation for case'
        print *,'       l =',l
        print *,'       m =',m
        print *,''
        print *,'test_Electrostatics failed'
        print *,''
        stop 1
      end if
      
    end do
  end do
  
  print *,''
  write (*,'(A,I1,A,I1,A,I1)') 'test_Electrostatics passed up to lmax=',lmax
  print *,''
  
end program

