!!
!! Solver for the Generalized Poisson Equation. 
!! @author Giuseppe Fisicaro

program GPS_3D

   use wrapper_mpi
   use Poisson_Solver
   use PSbox
   use yaml_output
   use dynamic_memory
   use dictionaries, dict_set => set
   use time_profiling
   use f_utils
   use yaml_strings
   implicit none
   
   !now these parameters have to be specified from the command line
!!$   integer, parameter :: n01 = 200
!!$   integer, parameter :: n02 = 200
!!$   integer, parameter :: n03 = 200
   integer, parameter :: nspden = 1
   character(len=4) :: PSol !, parameter :: PSol='PI'!'PCG' !    PCG = Preconditioned Conjugate Gradient Method.
                                             !    PI  = Polarization Iterative Method.
   ! To set 1 for analytical epsilon, 2 for analytical electron dependence,
   ! 3 for real electron density from cube file, 4 for rigid cavity.
   integer :: SetEps! = 1!3 
   logical :: usegpu
   logical :: lin_PB
   logical, parameter :: PCGstart = .false. !.true.
   logical, parameter :: PIstart = .false. !.true.
   logical, parameter :: CFgrid = .false. !.true.
   logical, parameter :: Fgrid = .false. !.true.

   real(kind=8), parameter :: acell = 10.d0 !10.d0
   real(kind=8), parameter :: rad_cav = 1.7d0 !1.7d0 ! Radius of the dielectric rigid cavity = rad_cav*acell (with nat=1).
   real(kind=8), parameter :: multp = 1.d0
   integer :: nat = 1 ! Number of atoms to build rigid cavity with nat=1.
   real(kind=8) :: erfL  ! To set 1 for Vacuum and correct analitic comparison with gaussian potential.
   real(kind=8) :: erfR  
   real(kind=8), parameter :: sigmaeps = 0.05d0*acell
   !> To set 1 for gaussian rho, 2 for gaussian potential. Beaware to use 2 or 3 when SetEps is setted to 1!!! 
   integer :: Setrho = 1 
!!$giu   character(len=2), parameter :: geocode = 'F'
!   character(len=2), parameter :: geocode = 'P'
   character(len=2) :: geocode
   character(len=2) :: geocodeprova
   character(len=2), parameter :: datacode = 'G'
   !> Order of accuracy for derivatives into ApplyLaplace subroutine 
   !!= Total number of points at left and right of the x0 where we want to calculate the derivative.
   integer, parameter :: nord = 16
   integer, dimension(3) :: ndims,ndimsc,ndimsf
   real(8), dimension(3) :: hgrids,hgridsc,hgridsf
   real(kind=8), parameter :: a_gauss = 1.0d0,a2 = a_gauss**2
   !integer :: m1,m2,m3,md1,md2,md3,nd1,nd2,nd3,n1,n2,n3,
   integer :: itype_scf,i_all,i_stat,n_cell,iproc,nproc,ixc,n01,n02,n03,iat,n_iter
   real(kind=8) :: hx,hy,hz,freq,fz,fz1,fz2,pi,curr,average,CondNum,wcurr,ave1,ave2,rhores2,En1,En2,dVnorm,hgrid,sume,delta
   real(kind=8) :: Adiag,ersqrt,ercurr,factor,r,r2,max_diff,max_diffpot,fact,x1,x2,x3,derf_tt,diffcurr,diffcurrS,divprod,einit
   real(kind=8) :: ehartree,offset,epol
   real(kind=8), dimension(:,:,:,:), allocatable :: density,rhopot,rvApp,rhoele,rhoion,potsol,rhopotf,densityf

   logical :: logyes
   ! Now start modification for check.
   real(kind=8), dimension(:,:,:,:,:), allocatable :: dens_check,pot_check
   integer :: i_check,unt,igpu
   ! To set 1 for normal run, 3 for check V[\rho,\epsilon] + V[\rho_ion,epsilon] is = to V[\rho + \rho_ion, epsilon]
   integer, parameter :: n_check = 1 

   real(kind=8), dimension(:,:), allocatable :: rxyz
   real(kind=8), dimension(:), allocatable :: radii
   real(kind=8), pointer :: kernel(:)
   type(coulomb_operator) :: pkernel
!   type(mpi_environment), intent(in), optional :: mpi_env
   real(kind=8), dimension(:,:,:), allocatable :: eps,potential,pot_ion,epsf,potentialf
   integer :: i1,i2,i3,isp,whichone,i,ii,j,info,icurr,ip,isd,i1_max,i2_max,i3_max,n3d,n3p,n3pi,i3xcsh,i3s,n3pr2,n3pr1,ierr
!   type(mpi_environment) :: bigdft_mpi
  type(dictionary), pointer :: options,dict_input

  real(kind=8), dimension(:,:,:,:), allocatable :: dlogeps
  !> inverse of epsilon. Needed for PI method.
  !! if absent, it will be calculated from the array of epsilon
  real(kind=8), dimension(:,:,:), allocatable :: oneoeps
  !> inverse square root of epsilon. Needed for PCG method.
  !! if absent, it will be calculated from the array of epsilon
  real(kind=8), dimension(:,:,:), allocatable :: oneosqrteps
  !> correction term of the Generalized Laplacian
  !! if absent, it will be calculated from the array of epsilon
  real(kind=8), dimension(:,:,:), allocatable :: corr

   call f_lib_initialize()

   !read command line
   call PS_Check_command_line_options(options)

   call f_zero(PSol)
   PSol=options .get. 'method'
   if (len_trim(PSol)==0) call f_strcpy(src='VAC',dest=PSol)
   ndims=options // 'ndim'
!!$giu
   geocode=options//'geocode'
   SetEps =options//'seteps'
   usegpu = options // 'accel'
   logyes= options // 'logfile'
   delta=0.3d0
   delta= options .get. 'deltacav'
   lin_PB = SetEps == 17

   call dict_init(dict_input)

   if ('input' .in. options) &
        call dict_copy(dest=dict_input,src=options//'input')


   call dict_free(options)


   igpu=0
   if (usegpu) igpu=1

   n01=ndims(1)
   n02=ndims(2)
   n03=ndims(3)
   ndimsc(1) = n01/2
   ndimsc(2) = n02/2
   ndimsc(3) = n03/2
   ndimsf(1) = n01*2
   ndimsf(2) = n02*2
   ndimsf(3) = n03*2

   if (SetEps.eq.3) then
    call get_size_from_cube(n01,n02,n03,hx,hy,hz,nat)
    ndims(1) = n01
    ndims(2) = n02
    ndims(3) = n03
   else
    hx=acell/real(n01-1,kind=8)
    hy=acell/real(n02-1,kind=8)
    hz=acell/real(n03-1,kind=8)
   end if
   hgridsc(1)=hx*2.d0 
   hgridsc(2)=hy*2.d0 
   hgridsc(3)=hz*2.d0 
   !hgridsc(1)=acell/real(n01/2-1,kind=8)
   !hgridsc(2)=acell/real(n02/2-1,kind=8)
   !hgridsc(3)=acell/real(n03/2-1,kind=8)
   hgridsf(1)=hx/2.d0 
   hgridsf(2)=hy/2.d0 
   hgridsf(3)=hz/2.d0 

   call mpiinit()
   iproc=mpirank()
   nproc=mpisize()

   !control memory profiling
   call f_malloc_set_status(memory_limit=0.e0,iproc=iproc)
   if (iproc ==0) then
        if (logyes) then
         call yaml_set_stream(record_length=92,tabbing=30,unit=70,filename='log.yaml',position='rewind')
        else
         call yaml_set_stream(record_length=92,tabbing=30)
        end if
      call yaml_new_document()
   end if

   erfL = 78.36d0 
   erfR = 1.0d0  
   if (iproc ==0) then
    call yaml_map('rad_cav',rad_cav)
    call yaml_map('multp',multp)
    if ((SetEps.eq.5).and.( trim(PSol)=='VAC')) then
     write(*,*)'Running a Generalized Poisson calculation'
    else if ((SetEps.eq.6).and.( trim(PSol)=='VAC')) then
     write(*,*)'Running a Poisson-Boltzmann calculation'
    else if ((SetEps.eq.7).and.( trim(PSol)=='VAC')) then
     write(*,*)'Running a Generalized Poisson calculation with restart'
    end if
   end if
!   call MPI_INIT(ierr)
!   call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
!   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

!   bigdft_mpi%mpi_comm=MPI_COMM_WORLD !workaround to be removed

   !allocate arrays to avoid stack overflow

   density=f_malloc([n01,n02,n03,nspden],id='density')
   rhopot =f_malloc([n01,n02,n03,nspden],id='rhopot')
   potsol =f_malloc([n01,n02,n03,nspden],id='potsol')
   rvApp  =f_malloc([n01,n02,n03,nspden],id='rvApp')
   rhoele =f_malloc([n01,n02,n03,nspden],id='rhoele')
   rhoion =f_malloc([n01,n02,n03,nspden],id='rhoion')
   dens_check =f_malloc([n01,n02,n03,nspden,n_check],id='dens_check')
   pot_check =f_malloc([n01,n02,n03,nspden,n_check],id='pot_check')
   rxyz   =f_malloc([3,nat],id='rxyz')
   radii   =f_malloc([nat],id='radii')
   rhopotf =f_malloc([ndimsf(1),ndimsf(2),ndimsf(3),nspden],id='rhopotf')
   densityf =f_malloc([ndimsf(1),ndimsf(2),ndimsf(3),nspden],id='densityf')
   
   eps=f_malloc([n01,n02,n03],id='eps')
   dlogeps=f_malloc([3,n01,n02,n03],id='dlogeps')
   oneoeps=f_malloc([n01,n02,n03],id='oneoeps')
   oneosqrteps=f_malloc([n01,n02,n03],id='oneosqrteps')
   corr=f_malloc([n01,n02,n03],id='corr')
   potential=f_malloc([n01,n02,n03],id='potential')
   pot_ion=f_malloc([n01,n02,n03],id='pot_ion')
   epsf=f_malloc([ndimsf(1),ndimsf(2),ndimsf(3)],id='epsf')
   potentialf=f_malloc([ndimsf(1),ndimsf(2),ndimsf(3)],id='potentialf')

   n_cell = max(n01,n02,n03)
   hgrid=max(hx,hy,hz)
   rxyz=0.d0
   ixc=0
   itype_scf=16
!!$   iproc=0
!!$   nproc=1
!-------------------------------------------------------------------------

   ! Create the Kernel.
   ! Calculate the kernel in parallel for each processor.

   hgrids=(/hx,hy,hz/)

   ehartree=0.0
!------------------------------------------------------------------------

!------------------------------------------------------------------------
   ! Set environment, namely permittivity epsilon, rhoele = electron charge density and rhoion = ion charge density .

   if (SetEps.eq.3) then
      call get_rho(n01,n02,n03,nspden,nat,acell,a_gauss,hx,hy,hz,rhoele,rhoion,sume,rxyz,iproc)
   else
    if (nat.eq.1) then
     !delta=0.30d0 !6.d0*max(hx,hy,hz)
     rxyz(1,1) = hx*real((n01-1)/2,kind=8)
     rxyz(2,1) = hy*real((n02-1)/2,kind=8)
     rxyz(3,1) = hz*real((n03-1)/2,kind=8)
!     rxyz(1,1) = hx*real(10,kind=8)
!     rxyz(2,1) = hy*real(10,kind=8)
!     rxyz(3,1) = hz*real(10,kind=8)
     radii(1)=rad_cav!*1.5d0/0.52917721092d0
    else if (nat.eq.2) then
     delta=0.3d0 !6.d0*max(hx,hy,hz)
     rxyz(1,1) = hx*real(20,kind=8)
     rxyz(2,1) = hy*real(10,kind=8)
     rxyz(3,1) = hz*real(20,kind=8)
     rxyz(1,2) = hx*real(10,kind=8)
     rxyz(2,2) = hy*real(10,kind=8)
     rxyz(3,2) = hz*real(10,kind=8)
     radii(1)=rad_cav!*1.5d0/0.52917721092d0
     radii(2)=rad_cav!*1.5d0/0.52917721092d0
    else if (nat.eq.3) then
      rxyz(1:3,1)=[9.300000d0, 9.300337d0, 9.243250d0]!-[2.30d0,2.85d0,3.7d0]
      rxyz(1:3,2)=[9.300000d0, 8.415319d0, 8.700265d0]!-[2.30d0,2.85d0,3.7d0]
      rxyz(1:3,3)=[9.300000d0, 7.299663d0, 9.156750d0]!-[2.30d0,2.85d0,3.7d0]
!      rxyz(1:3,4)=[7.300000d0, 7.299663d0,9.156750d0]-[2.30d0,2.85d0,3.7d0]
!      rxyz(1:3,5)=[7.300000d0, 7.299663d0,8.156750d0]-[2.30d0,2.85d0,3.7d0]
!      radii=[1.5d0,2.0d0,1.5d0]
!      radii=[1.5d0,2.0d0,1.5d0,1.5d0]
!      radii=[1.5d0,1.2d0,1.2d0,1.2d0,1.2d0]
!      rxyz(1:3,2)=[5.00000d0, 5.00000d0 + 0.2d0, 5.000000d0]

     delta=0.3d0
!     delta=2.0d0
!     delta=delta*0.25d0 ! Divided by 4 because both rigid cavities are 4*delta widespread
     radii=[1.4d0,1.0d0,1.0d0]
!     do iat=1,nat
!      radii(iat) = rad_cav*radii(iat)/0.52917721092d0 + 1.22d0*delta
!     end do
    end if

    rhoele(:,:,:,:) = 0.d0
    rhoion(:,:,:,:) = 0.d0
   end if

!   SetEps=4
     call SetEpsilon(n01,n02,n03,nspden,nord,nat,iproc,acell,a_gauss,hx,hy,hz,erfL,erfR,sigmaeps,&
          4,geocode,PSol,eps,dlogeps,oneoeps,oneosqrteps,corr,rhoele,rad_cav,rxyz,radii,delta)

    if (SetEps.lt.5) then
     if ( trim(PSol)=='VAC') then
      eps=1.d0
!      SetEps=1
      erfL=1.d0
      dlogeps=0.d0
      oneoeps=1.d0
      oneosqrteps=1.d0
      corr=0.d0
     end if
    end if

!    SetEps=1
!    Setrho=1
    call print_PB_function(n01,n02,n03,iproc,hx,hy,hz,nord,acell)
!------------------------------------------------------------------------
! Vacuum
!   eps=1.d0
!   dlogeps=0.d0
!   oneoeps=1.d0
!   oneosqrteps=1.d0
!   corr=0.d0
! Water
!   eps=78.36d0
!   dlogeps=0.d0
!   oneoeps=1.d0/78.36d0
!   oneosqrteps=1.d0/sqrt(78.36d0)
!   corr=0.d0

   ! Set initial density, and the associated analitical potential for the Standard Poisson Equation.
   call SetInitDensPot(n01,n02,n03,nspden,iproc,nat,eps,dlogeps,sigmaeps,SetEps,&
        erfL,erfR,acell,a_gauss,a2,hx,hy,hz,Setrho,density,potential,geocode,offset,einit,multp,rxyz,lin_PB)
!   call SetInitDensPot(n01,n02,n03,nspden,iproc,eps,dlogeps,sigmaeps,1,erfL,erfR,&
!   acell,a_gauss,a2,hx,hy,hz,1,density,potential,geocode,offset,einit,multp)

!   eps=1.d0
!   corr=0.d0
!   oneosqrteps=1.d0
!------------------------------------------------------------------------

!------------------------------------------------------------------------
   geocodeprova='F'
   ! Calculate the charge starting from the potential applying the proper Laplace operator.
   call ApplyLaplace(geocodeprova,n01,n02,n03,nspden,hx,hy,hz,potential,rvApp,acell,eps,nord,.false.,multp)

  if (iproc==0) then
   call yaml_mapping_open('Comparison between Generalized Poisson operator and analytical density')
   call writeroutinePot(n01,n02,n03,1,density,0,rvApp)
   call yaml_mapping_close()
  end if

   geocodeprova='F'
   ! Calculate the charge starting from the potential applying the proper Laplace operator.
   call ApplyLaplace2(geocodeprova,n01,n02,n03,nspden,hx,hy,hz,potential,rvApp,acell,eps,dlogeps,nord,.false.,multp)

  if (iproc==0) then
     call yaml_mapping_open('Comparison between Generalized Poisson operator 2 and analytical density')
     call writeroutinePot(n01,n02,n03,1,density,1,rvApp)
     call yaml_mapping_close()
  end if

!   geocodeprova='F'
!   ! Calculate the charge starting from the potential applying the improved ISF Laplace operator.
!   call ApplyLaplace_ISF(geocodeprova,n01,n02,n03,nspden,hx,hy,hz,potential,rvApp,acell,eps,nord,.false.,multp)
!
!  if (iproc==0) then
!   call yaml_comment('Comparison between improved ISF Generalized Poisson operator and analytical density')
!   call writeroutinePot(n01,n02,n03,1,density,0,rvApp)
!  end if

   geocodeprova='F'
   ! Calculate the charge starting from the potential applying the improved ISF Laplace operator PCG-style.
   call ApplyLaplace_corr(geocodeprova,n01,n02,n03,nspden,hx,hy,hz,potential,rvApp,acell,eps,corr,oneosqrteps,nord,.false.,multp)

  if (iproc==0) then
     call yaml_mapping_open('Comparison between Generalized Poisson operator PCG-style and analytical density')
     call writeroutinePot(n01,n02,n03,1,density,0,rvApp)
     call yaml_mapping_close()
  end if

!   geocodeprova='F'
!   ! Calculate the charge starting from the potential applying the improved ISF Laplace operator PCG-style.
!   call ApplyLaplace_ISF_corr(geocodeprova,n01,n02,n03,nspden,hx,hy,hz,potential,rvApp,acell,eps,corr,oneosqrteps,nord,.false.,multp)
!
!  if (iproc==0) then
!   call yaml_comment('Comparison between improved ISF Generalized Poisson operator PCG-style and analytical density')
!   call writeroutinePot(n01,n02,n03,1,density,0,rvApp)
!  end if
!------------------------------------------------------------------------

   geocodeprova='F'
   ! Calculate the charge starting from the potential applying the proper Laplace operator.
   call ApplyLaplace(geocodeprova,n01,n02,n03,nspden,hx,hy,hz,potential,rvApp,acell,eps,nord,lin_PB,multp)

  if (iproc==0) then
     call yaml_comment('Comparison between Poisson-Boltzmann operator and analytical density')
     call writeroutinePot(n01,n02,n03,1,density,0,rvApp)
  end if

!!$   max_diffpot = 0.d0
!!$   i1_max = 1
!!$   i2_max = 1
!!$   i3_max = 1
!!$   do i3=1,n03
!!$    do i2=1,n02
!!$     do i1=1,n01
!!$     fact=abs(density(i1,i2,i3,1) - rvApp(i1,i2,i3,1))
!!$      if (max_diffpot < fact) then
!!$       max_diffpot = fact
!!$       i1_max = i1
!!$       i2_max = i2
!!$       i3_max = i3
!!$      end if
!!$     end do
!!$    end do
!!$   end do
!!$   write(*,*)'Max abs difference between analytic density - ApplyLaplace to potential'
!!$   write(*,'(3(1x,I4),2(1x,e14.7))')i1_max,i2_max,i3_max,max_diffpot,abs(density(n01/2,n02/2,n03/2,1)-rvApp(n01/2,n02/2,n03/2,1))
  if (n_check.eq.1) then
!   dens_check(:,:,:,:,1) = rhoion(:,:,:,:) - rhoele(:,:,:,:)
   dens_check(:,:,:,:,1) = density(:,:,:,:)
  else if (n_check.eq.3) then
   dens_check(:,:,:,:,1) = rhoion(:,:,:,:)
   dens_check(:,:,:,:,2) = -rhoele(:,:,:,:)
   dens_check(:,:,:,:,3) = rhoion(:,:,:,:) - rhoele(:,:,:,:)
  end if

 do i_check=1,n_check
 
  if (n_check.eq.3 .and. SetEps.eq.3 .and. iproc==0) then
     !write(*,'(a)')'--------------------------------------------------------------------------------------------'
     call yaml_comment('Calculation for ion+ele potential difference '//i_check**'(i4)',hfill='-')
     !write(*,'(a,i4)')'Calculation for ion+ele potential difference',i_check
  end if

!  if (SetEps.eq.3) then
   rhopot(:,:,:,:) = dens_check(:,:,:,:,i_check)
!  else
!   rhopot(:,:,:,:) = density(:,:,:,:)
!  end if


   if (usegpu) call dict_set(dict_input//'setup'//'accel','CUDA')
   call dict_set(dict_input//'environment'//'delta',delta)
   if (trim(PSol) /= 'VAC') then
      call dict_set(dict_input//'environment'//'cavity','rigid')
      call dict_set(dict_input//'environment'//'gps_algorithm',PSol)
   end if

  !new method
   if (CFgrid) then
    pkernel=pkernel_init(iproc,nproc,dict_input,geocode,ndimsc,hgridsc)
   else if (Fgrid) then
    pkernel=pkernel_init(iproc,nproc,dict_input,geocode,ndimsf,hgridsf)
   else
    pkernel=pkernel_init(iproc,nproc,dict_input,geocode,ndims,hgrids)
   end if

!!$  pkernel=pkernel_init(.true.,iproc,nproc,igpu,geocode,ndims,hgrids,itype_scf,alg=PSol)
  call dict_free(dict_input)
  call pkernel_set(pkernel,verbose=.true.)

  !pkernel%minres=1.0d-12

  if ( trim(PSol)=='PCG') then
     
!!$     !set the coulomb operator of this system
!!$   if (SetEps.eq.2 .or. SetEps.eq.3) then
!!$     call pkernel_set(pkernel,verbose=.true.,eps=eps,oneosqrteps=oneosqrteps,corr=corr)
!!$   else
!!$     call pkernel_set(pkernel,verbose=.true.,eps=eps)
!!$   end if
     if (Fgrid) then
        !pkernel%minres=1.0d-12
        call back_trans_ISF_3D(geocode,8,ndimsf(1),ndimsf(2),ndimsf(3),eps,epsf)
        call back_trans_ISF_3D(geocode,8,ndimsf(1),ndimsf(2),ndimsf(3),rhopot,rhopotf)
        call pkernel_set_epsilon(pkernel,eps=epsf)
     else
      if (any(SetEps == [2,3,4,17]))  then
        call pkernel_set_epsilon(pkernel,oneosqrteps=oneosqrteps,corr=corr,eps=eps)
      else
        call pkernel_set_epsilon(pkernel,eps=eps)
      end if
     end if
!!$     call H_potential('G',pkernel,rhopot,rhopot,ehartree,0.d0,.false.)
!!$     if (iproc==0) call writeroutinePot(n01,n02,n03,nspden,rhopot,pkernel%max_iter,&
!!$          potential)

!!$     !old method, outside of PSolver
!!$     pkernel=pkernel_init(.true.,iproc,nproc,0,geocode,ndims,hgrids,itype_scf,method='VAC')
!!$     !set the coulomb operator of this system
!!$     call pkernel_set(pkernel,verbose=.true.)
!!$     call Prec_conjugate_gradient(n01,n02,n03,nspden,hx,hy,hz,rhopot,acell,eps,nord,pkernel,potential)

  else if (trim(PSol)=='PI') then

     !set the coulomb operator of this system
!!$   if (SetEps.eq.2 .or. SetEps.eq.3) then
!!$     call pkernel_set(pkernel,verbose=.true.,eps=eps,dlogeps=dlogeps,oneoeps=oneoeps)
!!$   else
!!$     call pkernel_set(pkernel,verbose=.true.,eps=eps)
!!$   end if
   if (any(SetEps == [2,3,4,17])) then
      call pkernel_set_epsilon(pkernel,eps=eps,oneoeps=oneoeps,dlogeps=dlogeps)
   else
      call pkernel_set_epsilon(pkernel,eps=eps)
   end if

!!$     call H_potential('G',pkernel,rhopot,rhopot,ehartree,0.d0,.false.)
!!$     if (iproc==0) then
!!$        call writeroutinePot(n01,n02,n03,nspden,rhopot,pkernel%max_iter,potential)
!!$     end if
     

!!$     !old method, outside of PSolver
!!$     pkernel=pkernel_init(.true.,iproc,nproc,0,geocode,ndims,hgrids,itype_scf,method='VAC')
!!$     call PolarizationIteration(n01,n02,n03,nspden,rhopot,acell,eps,nord,pkernel,potential)

!!$  else
!!$     call f_err_throw('Unrecognized method (provided "'//trim(PSol)//'")')
  end if

  select case(SetEps)
     !if (any(SetEps == [2,3,4])) then
  case(2,3,4)
     if (Fgrid) then
      call H_potential('D',pkernel,rhopotf(1,1,pkernel%grid%istart+1,1),&
           rhopotf(1,1,pkernel%grid%istart+1,1),ehartree,offset,.false.)
      call PS_gather(src=rhopotf,kernel=pkernel)
      call for_trans_ISF_3D(geocode,8,ndimsf(1),ndimsf(2),ndimsf(3),rhopotf,rhopot)
     else
      call H_potential('D',pkernel,rhopot(1,1,pkernel%grid%istart+1,1),rhopot(1,1,pkernel%grid%istart+1,1),&
           ehartree,offset,.false.)
      call PS_gather(src=rhopot,kernel=pkernel)
!!$      call H_potential('G',pkernel,rhopot,rhopot,ehartree,offset,.false.)
     end if
  case(5)
  !else if (any(SetEps == [5])) then
     call Prec_conjugate_gradient(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,density,acell,&
          eps,SetEps,nord,pkernel,potential,corr,oneosqrteps,dlogeps,multp,offset,geocode,lin_PB,.false.,CFgrid)
   if (PIstart) call PolarizationIteration(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,density,acell,&
        eps,nord,pkernel,potential,oneoeps,dlogeps,multp,offset,geocode,PIstart)
  !     call PolarizationIteration(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,acell,&
   ! eps,nord,pkernel,potential,oneoeps,dlogeps,multp,offset,geocode)
  case(6)
!   call Poisson_Boltzmann(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,acell,eps,SetEps,nord,pkernel,potential,corr,oneosqrteps,multp)
!   call Poisson_Boltzmann_improved(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,acell,&
    ! eps,SetEps,nord,pkernel,potential,corr,oneosqrteps,multp)
   call Poisson_Boltzmann_improved2(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,acell,eps,SetEps,nord,pkernel,&
        potential,corr,oneosqrteps,multp)
  case(7)
  ! Simple PCG with restart
     call Prec_conjugate_gradient_restart(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,acell,&
          eps,SetEps,nord,pkernel,potential,corr,oneosqrteps,multp,offset,geocode,lin_PB)
  case(8)
       call PolarizationIteration(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,density,acell,&
            eps,nord,pkernel,potential,oneoeps,dlogeps,multp,offset,geocode,.false.)
   if (PCGstart) call Prec_conjugate_gradient(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,density,acell,&
          eps,SetEps,nord,pkernel,potential,corr,oneosqrteps,dlogeps,multp,offset,geocode,lin_PB,PCGstart,CFgrid)
  case(9)
       call Prec_Steepest_Descent(n01,n02,n03,nspden,hx,hy,hz,rhopot,acell,eps,dlogeps,nord,pkernel,potential,geocode)
  case(10)
     call Prec_conjugate_gradient_PSD(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,acell,&
          eps,SetEps,nord,pkernel,potential,corr,oneosqrteps,oneoeps,dlogeps,multp,offset,geocode,lin_PB)
  case(11)
     call Prec_conjugate_gradient_Inputguess(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,density,acell,&
          eps,SetEps,nord,pkernel,potential,corr,oneosqrteps,dlogeps,multp,offset,geocode,lin_PB)
  case(12)
       call PolarizationIteration_Inputguess(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,density,acell,&
            eps,nord,pkernel,potential,oneoeps,dlogeps,multp,offset,geocode)
  case(13)
       call Prec_Steepest_Descent_Inputguess(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,density,acell,eps,&
            dlogeps,nord,pkernel,potential,geocode,multp)
  case(14)
       call Prec_conjugate_gradient(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,density,acell,&
            eps,SetEps,nord,pkernel,potential,corr,oneosqrteps,dlogeps,multp,offset,geocode,lin_PB,.false.,CFgrid)
       call PolarizationIteration_Inputguess(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,density,acell,&
            eps,nord,pkernel,potential,oneoeps,dlogeps,multp,offset,geocode)
  case(15)
       call Prec_conjugate_gradient(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,density,acell,&
            eps,SetEps,nord,pkernel,potential,corr,oneosqrteps,dlogeps,multp,offset,geocode,lin_PB,.false.,CFgrid)
       call Prec_Steepest_Descent_Inputguess(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,density,acell,eps,&
            dlogeps,nord,pkernel,potential,geocode,multp)
  case(16)
       call PolarizationIteration(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,density,acell,&
            eps,nord,pkernel,potential,oneoeps,dlogeps,multp,offset,geocode,.false.)
       call Prec_conjugate_gradient_Inputguess(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,density,acell,&
            eps,SetEps,nord,pkernel,potential,corr,oneosqrteps,dlogeps,multp,offset,geocode,lin_PB)
  case(17)
       call Poisson_Boltzmann_good(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,acell,eps,SetEps,nord,&
            pkernel,potential,density,oneoeps,dlogeps,corr,oneosqrteps,multp,geocode)
  end select

  pot_check(:,:,:,:,i_check) = rhopot(:,:,:,:)

  if (iproc==0) then
     call writeroutinePot(n01,n02,n03,nspden,rhopot,pkernel%max_iter,potential)
     call yaml_map('Expected hartree energy',einit)
     call yaml_map('Computed Hartree energy',ehartree)
     call yaml_map('Diff of expected-computed Hartree energy',einit-ehartree)
  end if

 end do ! End do check

   geocodeprova='F'
   ! Calculate the charge starting from the potential applying the proper Laplace operator.
   call ApplyLaplace(geocodeprova,n01,n02,n03,nspden,hx,hy,hz,rhopot(:,:,:,1),rvApp,acell,eps,nord,.false.,multp)

  if (iproc==0) then
   call yaml_mapping_open('Comparison between numerical and starting analytical density with old GPoperator')
   call writeroutinePot(n01,n02,n03,1,density,0,rvApp)
   call yaml_mapping_close()
  end if

   ! Calculate the charge starting from the potential applying the proper Laplace operator.
   call ApplyLaplace2(geocodeprova,n01,n02,n03,nspden,hx,hy,hz,rhopot(:,:,:,1),rvApp,acell,eps,dlogeps,nord,.false.,multp)

  if (iproc==0) then
   call yaml_mapping_open('Comparison between numerical and starting analytical density witn new 2 GPoperator')
   call writeroutinePot(n01,n02,n03,1,density,0,rvApp)
   call yaml_mapping_close()
  end if

!   ! Calculate the charge starting from the potential applying the improved ISF Laplace operator.
!   call ApplyLaplace_ISF(geocodeprova,n01,n02,n03,nspden,hx,hy,hz,rhopot(:,:,:,1),rvApp,acell,eps,nord,.false.,multp)
!
!  if (iproc==0) then
!   call yaml_comment('Comparison between improved ISF Generalized Poisson operator and analytical density',hfill='-')
!   call writeroutinePot(n01,n02,n03,1,density,0,rvApp)
!  end if

   ! Calculate the charge starting from the potential applying the improved ISF Laplace operator PCG-style.
   call ApplyLaplace_corr(geocodeprova,n01,n02,n03,nspden,hx,hy,hz,&
        rhopot(:,:,:,1),rvApp,acell,eps,corr,oneosqrteps,nord,.false.,multp)

  if (iproc==0) then
   call yaml_mapping_open('Comparison between Generalized Poisson operator PCG-style and analytical density')
   call writeroutinePot(n01,n02,n03,1,density,0,rvApp)
   call yaml_mapping_close()
  end if

!   ! Calculate the charge starting from the potential applying the improved ISF Laplace operator PCG-style.
!   call ApplyLaplace_ISF_corr(geocodeprova,n01,n02,n03,nspden,hx,hy,hz,rhopot(:,:,:,1),rvApp,acell,eps,corr,oneosqrteps,nord,.false.,multp)
!
!  if (iproc==0) then
!   call yaml_comment('Comparison between improved ISF Generalized Poisson operator PCG-style and analytical density',hfill='-')
!   call writeroutinePot(n01,n02,n03,1,density,0,rvApp)
!  end if
!
! call Polarization_charge(n01,n02,n03,nspden,hx,hy,hz,rhopot,rvApp,acell,eps,nord)

!------------------------------------------------------------------
 if (n_check.eq.3 .and. SetEps.eq.3) then

  write(*,'(a)')'--------------------------------------------------------------------------------------------'
  write(*,'(a)')'Difference between V[\rho,\epsilon] + V[\rho_ion,epsilon] and V[\rho + \rho_ion, epsilon]'
  potential(:,:,:)= pot_check(:,:,:,1,1) + pot_check(:,:,:,1,2)
  call writeroutinePot(n01,n02,n03,nspden,rhopot,pkernel%max_iter,potential)

      unt=f_get_free_unit(21)
      call f_open_file(unt,file='final_ion_ele.dat')
      i1=(n01-1)/2+1
      do i2=1,n02
         do i3=1,n03
            write(unt,'(2(1x,I4),2(1x,e14.7))')i2,i3,pot_check(i1,i2,i3,1,1),pot_check(i1,i2,i3,1,2)
         end do
      end do
      call f_close(unt)

      unt=f_get_free_unit(22)
      call f_open_file(unt,file='final_ion_ele_line.dat')
      i1=(n01-1)/2+1
      i3=(n03-1)/2+1
      do i2=1,n02
       write(unt,'(1x,I8,2(1x,e22.15))')i2,pot_check(i1,i2,i3,1,1),pot_check(i1,i2,i3,1,2)
      end do
      call f_close(unt)

 end if

!------------------------------------------------------------------

!!$  call pkernel_free(pkernel)
!!$
!!$  PSol='VAC'
!!$  eps=1.d0
!!$  erfL=1.d0
!!$  pkernel=pkernel_init(.true.,iproc,nproc,0,geocode,ndims,hgrids,itype_scf,alg=PSol)
!!$  call pkernel_set(pkernel,verbose=.true.)
!!$
!!$  potsol(:,:,:,:)=dens_check(:,:,:,:,1)
!!$
!!$  call H_potential('G',pkernel,potsol,potsol,ehartree,offset,.false.)
!!$
!!$  unt=f_get_free_unit(23)
!!$  call f_open_file(unt,file='finalpot_vacuum_line.dat')
!!$  do i2=1,n02
!!$   write(unt,'(1x,I8,2(1x,e22.15))')i2,potsol(n01/2,i2,n03/2,1),potsol(1,i2,n03/2,1)
!!$  end do
!!$  call f_close(unt)
!!$
!!$
!!$!-------------------------------------------------
!!$    !calculate polarization energy.
!!$    epol=0.d0
!!$    isp=1
!!$    do i3=1,n03
!!$     do i2=1,n02
!!$      do i1=1,n01
!!$       epol= epol + potsol(i1,i2,i3,isp)*rvApp(i1,i2,i3,isp)
!!$      end do
!!$     end do
!!$    end do
!!$    epol=0.5*hx*hy*hz*epol
!!$  if (iproc==0) then
!!$    call yaml_map('Vacuum Hartree energy',ehartree)
!!$    call yaml_map('Computed polarization energy in hartree',epol)
!!$     epol=epol*627.509469d0
!!$     call yaml_map('Computed polarization energy in kcal/mol',epol)
!!$
!!$  end if

!-------------------------------------------------

  call pkernel_free(pkernel)
  call f_free(density)
  call f_free(densityf)
  call f_free(rhopot)
  call f_free(rhopotf)
  call f_free(potsol)
  call f_free(rhoele)
  call f_free(rhoion)
  call f_free(dens_check)
  call f_free(pot_check)
  call f_free(rxyz)
  call f_free(radii)
  call f_free(rvApp)
  call f_free(eps)
  call f_free(epsf)
  call f_free(dlogeps)
  call f_free(oneoeps)
  call f_free(oneosqrteps)
  call f_free(corr)
  call f_free(potential)
  call f_free(potentialf)
  call f_free(pot_ion)

!!$  if (iproc ==0) then
!!$     call yaml_release_document()
!!$     call yaml_close_all_streams()
!!$  end if
  call mpifinalize()
  call f_lib_finalize()
  
contains

  !>identify the options from command line
  !! and write the result in options dict
  subroutine PS_Check_command_line_options(options)
    use yaml_parse
    use dictionaries
    implicit none
    !> dictionary of the options of the run
    !! on entry, it contains the options for initializing
    !! on exit, it contains in the key "BigDFT", a list of the 
    !! dictionaries of each of the run that the local instance of BigDFT
    !! code has to execute
    type(dictionary), pointer :: options
    !local variables
    type(yaml_cl_parse) :: parser !< command line parser

    !define command-line options
    parser=yaml_cl_parse_null()
    !between these lines, for another executable using BigDFT as a blackbox,
    !other command line options can be specified
    !then the bigdft options can be specified
    call PS_check_options(parser)
    !parse command line, and retrieve arguments
    call yaml_cl_parse_cmd_line(parser,args=options)
    !free command line parser information
    call yaml_cl_parse_free(parser)

  end subroutine PS_Check_command_line_options

end program GPS_3D

subroutine PS_Check_options(parser)
  use yaml_parse
  use dictionaries, only: dict_new,operator(.is.)
  implicit none
  type(yaml_cl_parse), intent(inout) :: parser

  call yaml_cl_parse_option(parser,'ndim','None',&
       'Domain Sizes','n',&
       dict_new('Usage' .is. &
       'Sizes of the simulation domain of the check',&
       'Allowed values' .is. &
       'Yaml list of integers. If a scalar integer is given, all the dimensions will have this size.'),first_option=.true.)

  call yaml_cl_parse_option(parser,'geocode','F',&
       'Boundary conditions','g',&
       dict_new('Usage' .is. &
       'Set the boundary conditions of the run',&
       'Allowed values' .is. &
       'String scalar. "F","S","W","P" boundary conditions are allowed'))

  call yaml_cl_parse_option(parser,'method','None',&
       'Embedding method','m',&
       dict_new('Usage' .is. &
       'Set the embedding method used. A non present value implies vacuum treatment.',&
       'Allowed values' .is. &
       dict_new("PI" .is. 'Polarization iteration Method',&
                "PCG" .is. 'Preconditioned Conjugate Gradient')))

  call yaml_cl_parse_option(parser,'seteps','4',&
       'Epsilon determination method','e',&
       dict_new('Usage' .is. &
       'Set the dielectric constant determination method.',&
       'Allowed values' .is. &
       dict_new('1' .is. 'Analytical epsilon' ,&
                '2' .is. 'analytical electron dependence',&
                '3' .is. 'real electron density from cube file (need electroninc_density.cube)',&
                '4' .is. 'calculate the cavity and dump it on disk',&
                '5' .is. 'Solves GPe with PCG customized (should be identical to 4 + PCG)',&
                '6' .is. 'Modified Poisson Botzmann Equation solver',&
                '7' .is. 'Solves GPe with PCG customized and a restart is implemented',&
                '8' .is. 'Solves GPe with PI customized (should be identical to 4 + PI)',&
                '9' .is. 'Solves GPe with PSD',&
                '10' .is. 'Solves GPe with PCG customized and is coupled with PSD')))

  call yaml_cl_parse_option(parser,'accel','No',&
       'GPU Acceleration','a',&
       dict_new('Usage' .is. &
       'Boolean, set the GPU acceleration'))

  call yaml_cl_parse_option(parser,'logfile','Yes',&
       'Write logfile','l',&
       dict_new('Usage' .is. &
       'Boolean, set the logfile as log.yaml'))

  call yaml_cl_parse_option(parser,'deltacav','None',&
       'Delta rigid cavity','d',&
       dict_new('Usage' .is. &
       'Sizes of the delta for error function',&
       'Allowed values' .is. &
       'Real value'))

  call yaml_cl_parse_option(parser,'input','None',&
       'Inputfile of Poisson solver','i',&
       dict_new('Usage' .is. &
       'Put the dictionary of PS inputfile to preload some parameters',&
       'Allowed values' .is. &
       'yaml Dictionary'))


end subroutine PS_Check_options

subroutine print_PB_function(n01,n02,n03,iproc,hx,hy,hz,nord,acell)

  use f_utils

  implicit none

  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: iproc
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell

  !local variables
  integer :: i,unt,n
  real(kind=8) :: PB_charge,v,dv
  integer, parameter :: n_points = 20001
  real(kind=8), parameter :: length = 1.d0 ! One side length.
  real(kind=8), dimension(n_points) :: func

  unt=f_get_free_unit(25)
  call f_open_file(unt,file='PB_function.dat')

  n=(n_points-1)/2
  dv=length/real(n,kind=8)
  do i=-n,n
   v=dv*real(i,kind=8)
   func(n+i+1)=PB_charge(v)
   write(unt,'(1x,I8,2(1x,e22.15))')i,v,func(n+i+1)
  end do

  call f_close(unt)

end subroutine print_PB_function

!> Calcultion of the Poisson-Boltzmann function.
!! following the definitions given in J. J. López-García, J. Horno, C. Grosse Langmuir 27, 13970-13974 (2011).
pure function PB_charge(x) result(ions_conc)

 use numerics, only: safe_exp

  implicit none

  !> argument of the Poisson-Botzmann function
  real(8), intent(in) :: x
  real(8) :: ions_conc
  ! Values to be given
  integer, parameter :: n_ions = 2 !< number of ionic species in the dielectric liquid system
  real(8), dimension(n_ions) :: z_ions !< valence of ionic species
  real(8), dimension(n_ions) :: c_ions !< bulk concentations of ionic species [mol/m^3]
  real(8), dimension(n_ions) :: c_max  !< maximum local concentration that ionic species can attain [mol/m^3]
  real(8), dimension(n_ions) :: r_ions !< effective ionic radius of ionic species [m]
  real(8), parameter :: Temp = 300 ! Temperature of the liquid system [K]
  !> packing coefficient p = 1 for perfect packing, p = pi_greek/(3(2)^{1/2}) ≈ 0.74 for close packing,
  real(8), parameter :: p = 0.74d0 
                                   !! p ≈ 0.64 for random close packing, and p = pi_greek/6 ≈ 0.52 for simple cubic packing.
  ! Nedeed constant
  real(8), parameter :: n_avo = 6.0221412927d23 ! Avogadro's number [1/mol]
  real(8), parameter :: k_b = 3.166811429d-6 ! Boltzmann constant in atomic unit [E_{H}/K]
  real(8), parameter :: bohr = 5.291772109217d-11 ! m
  !local variables
  integer :: i,j
  real(8) :: pi,fact,vol_bohr,K_bT,t,fact1,sumc,y,h,l
  real(8), dimension(n_ions) :: c_ratio  !< c_ions/c_max
  integer, parameter :: PBeq=3 ! Set 1 for linear, 2 for standard, 3 for
                               ! modified Poisson-Boltzmann equation.

  pi = 4.d0*datan(1.d0)
  k_bT = k_b*Temp
  vol_bohr=bohr*bohr*bohr
  fact=n_avo*vol_bohr
  fact1=(4.d0/3.d0)*pi*n_avo
  l=0.d0

  z_ions(1)=1.0d0
  z_ions(2)=-1.0d0
  c_ions(1)=100.0d0
  c_ions(2)=100.0d0
  r_ions(1)=3.0d-10
  r_ions(2)=3.0d-10
  sumc=0.d0
  do i=1,n_ions
   c_max(i)=p/(fact1*(r_ions(1)**3))
   c_ratio(i)=c_ions(i)/c_max(i)
   sumc=sumc+c_ratio(i)
  end do

  !--------------------------------------------------------
  if (PBeq.eq.1) then
   ! Linear Poisson-Boltzmann Equation.
   ions_conc = 0.d0
    do i=1,n_ions
     t = -z_ions(i)*x/k_bT !*0.01d0
     ions_conc = ions_conc + z_ions(i)*c_ions(i)*t
    end do
    ions_conc = ions_conc*fact!*1.d3  
  else if (PBeq.eq.2) then
   ! Standard Poisson-Boltzmann Equation.
   ions_conc = 0.d0
    do i=1,n_ions
     t = -z_ions(i)*x/k_bT!*0.01d0
     t=safe_exp(t) ! Comment this line for linear Poisson-Boltzmann Equation.
     ions_conc = ions_conc + z_ions(i)*c_ions(i)*t
    end do
    ions_conc = ions_conc*fact!*1.d3  
  else if (PBeq.eq.3) then  
   ! Modified Poisson-Boltzmann Equation.
   ions_conc = 0.d0
   do i=1,n_ions
    y=x/k_bT ! correct one
    !y=x/k_bT*0.05d0
    t = -z_ions(i)*y 
    h=0.d0
    do j=1,n_ions
     h=h+c_ratio(j)*safe_exp((z_ions(i)-z_ions(j))*y)
    end do
    l=safe_exp(z_ions(i)*y)*(1.d0-sumc)+h
    t=1.d0/l
    ions_conc = ions_conc + z_ions(i)*c_ions(i)*t 
   end do
   ions_conc = ions_conc*fact ! correct one
   !ions_conc = ions_conc*fact*5.0d2
  end if
  !--------------------------------------------------------

end function PB_charge

subroutine PolarizationIteration(n01,n02,n03,nspden,iproc,&
     hx,hy,hz,b,bb,acell,eps,nord,pkernel,potential,oneoeps,dlogeps,multp,offset,geocode,PIstart)
  use yaml_output
  use Poisson_Solver
  use wrapper_linalg
  use dynamic_memory
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden,iproc
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp,offset
  type(coulomb_operator), intent(inout) :: pkernel
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps,potential,oneoeps
  real(kind=8), dimension(3,n01,n02,n03), intent(in) :: dlogeps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: b
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: bb
  character(len=2), intent(in) :: geocode
  logical, intent(in) :: PIstart

  real(kind=8), dimension(n01,n02,n03)  :: pot_ion
  real(kind=8), parameter :: eta = 0.6d0 !1.0d0 ! Polarization Iterative Method parameter.
  real(kind=8), parameter :: taupol = 1.0d-8 ! Polarization Iterative Method parameter.
  integer, parameter :: maxiterpol=50
  !real(kind=8), dimension(n01,n02,n03,nspden,3) :: dlv
  real(kind=8), dimension(:,:,:,:), allocatable :: dlv,deps,rhosol,rhopol,rhotot
  real(kind=8), dimension(:,:,:,:), allocatable :: rhopolnew,rhopolold,rhores,lv
  integer :: i1,i2,i3,i,j,ip,isp
  real(kind=8) :: divprod,rhores2,diffcurr,pi,ehartree,res,rho,normr,rpoints
  
  rhosol=f_malloc([n01,n02,n03,nspden],id='rhosol')
  rhopol=f_malloc([n01,n02,n03,nspden],id='rhopol')
  rhotot=f_malloc([n01,n02,n03,nspden],id='rhotot')
  rhopolnew=f_malloc([n01,n02,n03,nspden],id='rhopolnew')
  rhopolold=f_malloc([n01,n02,n03,nspden],id='rhopolold')
  rhores=f_malloc([n01,n02,n03,nspden],id='rhores')
  lv=f_malloc([n01,n02,n03,nspden],id='lv')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  dlv=f_malloc([3,n01,n02,n03],id='dlv')

  pi = 4.d0*datan(1.d0)
  rpoints=product(real([n01,n02,n03],kind=8))

  open(unit=18,file='PI_normr_'//trim(geocode)//'.dat',status='unknown')
  open(unit=38,file='PI_accuracy_'//trim(geocode)//'.dat',status='unknown')

  if (iproc ==0) then
   write(18,'(1x,a)')'iter normr rhores2'
   write(38,'(1x,a)')'iter i1_max i2_max i3_max max_val max_center'
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting Polarization Iteration '
  end if

  if (iproc ==0) then
   call yaml_map('rpoints',rpoints)
   call yaml_sequence_open('Embedded PSolver, Polarization Iteration Method')
  end if

!  call fssnord3DmatNabla3var(n01,n02,n03,nspden,hx,hy,hz,eps,deps,nord,acell)

  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     rhopol(i1,i2,i3,isp) = 0.d0
     rhosol(i1,i2,i3,isp) = bb(i1,i2,i3,isp)

!     !switch and create the logarithmic derivative of epsilon
!     dlv(1,i1,i2,i3)=deps(i1,i2,i3,1)/eps(i1,i2,i3)
!     dlv(2,i1,i2,i3)=deps(i1,i2,i3,2)/eps(i1,i2,i3)
!     dlv(3,i1,i2,i3)=deps(i1,i2,i3,3)/eps(i1,i2,i3)

    end do
   end do
  end do

  if (PIstart) then
   call fssnord3DmatNabla_LG2(n01,n02,n03,hx,hy,hz,b,nord,acell,1.d0,dlogeps,eps,rhopol,rhores2)
   normr=sqrt(rhores2/rpoints)
  end if


!    call writeroutinePot(n01,n02,n03,nspden,potential,0,potential)
!    call writeroutine(n01,n02,n03,nspden,rhosol,0)

    do ip=1,maxiterpol

     if (iproc ==0) then
      write(*,'(a)')'--------------------------------------------------------------------------------------------'
      write(*,*)'Starting PI iteration ',ip
     end if

     isp=1
     do i3=1,n03
      do i2=1,n02
       do i1=1,n01
          !rhotot(i1,i2,i3,isp)=(1/eps(i1,i2,i3))*rhosol(i1,i2,i3,isp)+rhopol(i1,i2,i3,isp)
          !rhotot(i1,i2,i3,isp)=pkernel%oneoeps(i1,i2,i3)*rhosol(i1,i2,i3,isp)+rhopol(i1,i2,i3,isp)
          rhotot(i1,i2,i3,isp)=oneoeps(i1,i2,i3)*rhosol(i1,i2,i3,isp)+rhopol(i1,i2,i3,isp)
          lv(i1,i2,i3,isp) = rhotot(i1,i2,i3,isp)
       end do
      end do
     end do
     
     call yaml_sequence(advance='no')
     call H_potential('G',pkernel,lv,pot_ion,ehartree,offset,.false.)

     if (iproc ==0) then
      call writeroutinePot(n01,n02,n03,nspden,lv,ip,potential)
     end if

     !call fssnord3DmatNabla_LG2(n01,n02,n03,hx,hy,hz,lv,nord,acell,eta,pkernel%dlogeps,rhopol,rhores2)
     call fssnord3DmatNabla_LG2(n01,n02,n03,hx,hy,hz,lv,nord,acell,eta,dlogeps,eps,rhopol,rhores2)
     !call fssnord3DmatNabla_LG2(n01,n02,n03,hx,hy,hz,lv,nord,acell,eta,dlv,rhopol,rhores2)

     normr=sqrt(rhores2/rpoints)
     !!!call fssnord3DmatNabla(n01,n02,n03,nspden,hx,hy,hz,lv,dlv,nord,acell)
     !!!
     !!!isp=1
     !!!rhores2=0.d0
     !!!do i3=1,n03
     !!! do i2=1,n02
     !!!  do i1=1,n01
     !!!   divprod = 0.d0
     !!!   do j=1,3
     !!!      divprod = divprod + deps(i1,i2,i3,j)*dlv(i1,i2,i3,isp,j)
     !!!   end do
     !!!   !rhopolnew(i1,i2,i3,isp)=(1/(4.d0*pi))*(1/eps(i1,i2,i3))*divprod
     !!!   res=(1/(4.d0*pi))*(1/eps(i1,i2,i3))*divprod
     !!!   rho=rhopol(i1,i2,i3,isp)
     !!!   res=res-rho
     !!!   res=eta*res
     !!!   rhores2=rhores2+res*res
     !!!   rhopol(i1,i2,i3,isp)=res+rho
!!$  !!!      rhopolold(i1,i2,i3,isp)=rhopol(i1,i2,i3,isp)
!!$  !!!      rhopol(i1,i2,i3,isp)=eta*rhopolnew(i1,i2,i3,isp) + (1.d0-eta)*rhopolold(i1,i2,i3,isp)
!!$  !!!      rhores(i1,i2,i3,isp) = rhopol(i1,i2,i3,isp) - rhopolold(i1,i2,i3,isp)
     !!!  end do
     !!! end do
     !!!end do

!!$     call axpy(n01*n02*n03,1.d0,rhopol(1,1,1,1),1,rhores(1,1,1,1),1)
!!$     call axpy(n01*n02*n03,eta,rhores(1,1,1,1),1,rhopol(1,1,1,1),1)
!!$     rhores2=dot(n01*n02*n03,rhores(1,1,1,1),1,rhores(1,1,1,1),1)

!!$     rhores2 = 0.d0
!!$     isp=1
!!$     do i3=1,n03
!!$      do i2=1,n02
!!$       do i1=1,n01
!!$        rhores2 = rhores2 + rhores(i1,i2,i3,isp)*rhores(i1,i2,i3,isp)
!!$       end do
!!$      end do
!!$     end do

    if (iproc ==0) then
     write(18,'(1x,I8,2(1x,e14.7))')ip,normr,rhores2
    end if
     !write(*,'(1x,I8,1x,e14.7)')ip,rhores2

     call EPS_iter_output_LG(ip,0.0_dp,normr,0.0_dp,0.0_dp,0.0_dp)
     if (normr.lt.taupol) exit

!     call writeroutine(n01,n02,n03,nspden,rhores,ip)

  end do

    isp=1
    do i3=1,n03
     do i2=1,n02
      do i1=1,n01
       b(i1,i2,i3,isp) = lv(i1,i2,i3,isp)
      end do
     end do
    end do

    call yaml_sequence_close()

    close(unit=18)
    close(unit=38)
    !write(*,*)
    !write(*,'(1x,a,1x,i8)')'Polarization iterations =',ip
    !write(*,'(1x,a,1x,e14.7)')'rhores polarization square =',rhores2
    !write(*,*)
    !write(*,'(a)')'Max abs difference between analytic potential and the computed one'
    if (iproc ==0) then
     call writeroutinePot(n01,n02,n03,nspden,b,ip,potential)
    end if
    !write(*,*)
    !write(*,'(a)')'Termination of Polarization Iteration'
    !write(*,'(a)')'--------------------------------------------------------------------------------------------'

    call f_free(rhosol)
    call f_free(rhopol)
    call f_free(rhotot)
    call f_free(rhopolnew)
    call f_free(rhopolold)
    call f_free(rhores)
    call f_free(lv)
    call f_free(deps)
    call f_free(dlv)

end subroutine PolarizationIteration

subroutine PolarizationIteration_Inputguess(n01,n02,n03,nspden,iproc,&
     hx,hy,hz,b,bb,acell,eps,nord,pkernel,potential,oneoeps,dlogeps,multp,offset,geocode)
  use yaml_output
  use Poisson_Solver
  use wrapper_linalg
  use dynamic_memory
  use f_utils
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden,iproc
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp,offset
  type(coulomb_operator), intent(inout) :: pkernel
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps,potential,oneoeps
  real(kind=8), dimension(3,n01,n02,n03), intent(in) :: dlogeps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: b
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: bb
  character(len=2), intent(in) :: geocode
  logical :: PIig

  real(kind=8), dimension(n01,n02,n03)  :: pot_ion
  real(kind=8), parameter :: eta = 0.6d0 !1.0d0 ! Polarization Iterative Method parameter.
  real(kind=8), parameter :: taupol = 1.0d-12 ! Polarization Iterative Method parameter.
  integer, parameter :: maxiterpol=50
  integer, parameter :: max_scf=1
  !real(kind=8), dimension(n01,n02,n03,nspden,3) :: dlv
  real(kind=8), dimension(:,:,:,:), allocatable :: dlv,deps,rhosol,rhopol,rhotot
  real(kind=8), dimension(:,:,:,:), allocatable :: rhopolnew,rhopolold,rhores,lv,rvApp
  integer :: i1,i2,i3,i,j,ip,isp,i_scf
  real(kind=8) :: divprod,rhores2,diffcurr,pi,ehartree,res,rho,normr,rpoints,errorcur
  
  rhosol=f_malloc([n01,n02,n03,nspden],id='rhosol')
  rhopol=f_malloc([n01,n02,n03,nspden],id='rhopol')
  rhotot=f_malloc([n01,n02,n03,nspden],id='rhotot')
  rhopolnew=f_malloc([n01,n02,n03,nspden],id='rhopolnew')
  rhopolold=f_malloc([n01,n02,n03,nspden],id='rhopolold')
  rhores=f_malloc([n01,n02,n03,nspden],id='rhores')
  lv=f_malloc([n01,n02,n03,nspden],id='lv')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  dlv=f_malloc([3,n01,n02,n03],id='dlv')
  rvApp=f_malloc([n01,n02,n03,nspden],id='rvApp')

  pi = 4.d0*datan(1.d0)
  rpoints=product(real([n01,n02,n03],kind=8))

  open(unit=18,file='PI_normr_'//trim(geocode)//'.dat',status='unknown')
  open(unit=38,file='PI_accuracy_'//trim(geocode)//'.dat',status='unknown')

  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting an Inputguess PCG calculation'
  end if

  do i_scf=1,max_scf ! scf loop external

  if (i_scf.eq.1) then
   errorcur=taupol
   PIig=.true.
  else
   !PIig=.false.
   PIig=.true.
   errorcur=taupol
   !errorcur=errorcur
  end if

  call f_zero(rvApp)
  if (PIig) then

   call f_memcpy(src=b,dest=lv)

   ! Calculate the charge starting from the potential applying the proper
   ! Laplace operator.
   call ApplyLaplace(geocode,n01,n02,n03,nspden,hx,hy,hz,lv,rvApp,acell,eps,nord,.false.,multp)

   if (iproc==0) then
    call yaml_mapping_open('Comparison between Generalized Poisson operator and analytical density')
    call writeroutinePot(n01,n02,n03,1,bb,0,rvApp)
    call yaml_mapping_close()
   end if

!   ! Calculate the charge starting from the potential applying the proper
!   Laplace operator.
   call ApplyLaplace2(geocode,n01,n02,n03,nspden,hx,hy,hz,lv,rvApp,acell,eps,dlogeps,nord,.false.,multp)

  if (iproc==0) then
   call yaml_mapping_open('Comparison between Generalized Poisson operator 2 and analytical density')
   call writeroutinePot(n01,n02,n03,1,bb,1,rvApp)
   call yaml_mapping_close()
  end if

  end if

  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------!'
   write(*,*)'Starting Inputguess iteration ',i_scf
  end if


  if (iproc ==0) then
   write(18,'(1x,a)')'iter normr rhores2'
   write(38,'(1x,a)')'iter i1_max i2_max i3_max max_val max_center'
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting Polarization Iteration '
  end if

  if (iproc ==0) then
   call yaml_map('rpoints',rpoints)
   call yaml_sequence_open('Embedded PSolver, Polarization Iteration Method')
  end if

!  call fssnord3DmatNabla3var(n01,n02,n03,nspden,hx,hy,hz,eps,deps,nord,acell)

  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     rhopol(i1,i2,i3,isp) = 0.d0
     rhosol(i1,i2,i3,isp) = bb(i1,i2,i3,isp)

!     !switch and create the logarithmic derivative of epsilon
!     dlv(1,i1,i2,i3)=deps(i1,i2,i3,1)/eps(i1,i2,i3)
!     dlv(2,i1,i2,i3)=deps(i1,i2,i3,2)/eps(i1,i2,i3)
!     dlv(3,i1,i2,i3)=deps(i1,i2,i3,3)/eps(i1,i2,i3)

    end do
   end do
  end do

  if (PIig) then
   call fssnord3DmatNabla_LG2(n01,n02,n03,hx,hy,hz,lv,nord,acell,1.d0,dlogeps,eps,rhopol,rhores2)
   normr=sqrt(rhores2/rpoints)
  end if


!    call writeroutinePot(n01,n02,n03,nspden,potential,0,potential)
!    call writeroutine(n01,n02,n03,nspden,rhosol,0)

    do ip=1,maxiterpol

     if (iproc ==0) then
      write(*,'(a)')'--------------------------------------------------------------------------------------------'
      write(*,*)'Starting PI iteration ',ip
     end if

     isp=1
     do i3=1,n03
      do i2=1,n02
       do i1=1,n01
          !rhotot(i1,i2,i3,isp)=(1/eps(i1,i2,i3))*rhosol(i1,i2,i3,isp)+rhopol(i1,i2,i3,isp)
          !rhotot(i1,i2,i3,isp)=pkernel%oneoeps(i1,i2,i3)*rhosol(i1,i2,i3,isp)+rhopol(i1,i2,i3,isp)
          rhotot(i1,i2,i3,isp)=oneoeps(i1,i2,i3)*rhosol(i1,i2,i3,isp)+rhopol(i1,i2,i3,isp)
          lv(i1,i2,i3,isp) = rhotot(i1,i2,i3,isp)
       end do
      end do
     end do
     
     call yaml_sequence(advance='no')
     call H_potential('G',pkernel,lv,pot_ion,ehartree,offset,.false.)

     if (iproc ==0) then
      call writeroutinePot(n01,n02,n03,nspden,lv,ip,potential)
     end if

     !call fssnord3DmatNabla_LG2(n01,n02,n03,hx,hy,hz,lv,nord,acell,eta,pkernel%dlogeps,rhopol,rhores2)
     call fssnord3DmatNabla_LG2(n01,n02,n03,hx,hy,hz,lv,nord,acell,eta,dlogeps,eps,rhopol,rhores2)
     !call fssnord3DmatNabla_LG2(n01,n02,n03,hx,hy,hz,lv,nord,acell,eta,dlv,rhopol,rhores2)

     normr=sqrt(rhores2/rpoints)

    if (iproc ==0) then
     write(18,'(1x,I8,2(1x,e14.7))')ip,normr,rhores2
    end if
     !write(*,'(1x,I8,1x,e14.7)')ip,rhores2

     call EPS_iter_output_LG(ip,0.0_dp,normr,0.0_dp,0.0_dp,0.0_dp)
     if (normr.lt.errorcur) exit

!     call writeroutine(n01,n02,n03,nspden,rhores,ip)

  end do ! PI loop

   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'End Polarization Iteration '

  end do ! scf external loop

    isp=1
    do i3=1,n03
     do i2=1,n02
      do i1=1,n01
       b(i1,i2,i3,isp) = lv(i1,i2,i3,isp)
      end do
     end do
    end do

    call yaml_sequence_close()

    close(unit=18)
    close(unit=38)
    !write(*,*)
    !write(*,'(1x,a,1x,i8)')'Polarization iterations =',ip
    !write(*,'(1x,a,1x,e14.7)')'rhores polarization square =',rhores2
    !write(*,*)
    !write(*,'(a)')'Max abs difference between analytic potential and the computed one'
    if (iproc ==0) then
     call writeroutinePot(n01,n02,n03,nspden,b,ip,potential)
    end if
    write(*,*)
    write(*,'(a)')'End external scf loop'
    write(*,'(a)')'--------------------------------------------------------------------------------------------'

    call f_free(rhosol)
    call f_free(rhopol)
    call f_free(rhotot)
    call f_free(rhopolnew)
    call f_free(rhopolold)
    call f_free(rhores)
    call f_free(rvApp)
    call f_free(lv)
    call f_free(deps)
    call f_free(dlv)

end subroutine PolarizationIteration_Inputguess

subroutine Prec_conjugate_gradient(n01,n02,n03,nspden,iproc,hx,hy,hz,b,bb,&
     acell,eps,SetEps,nord,pkernel,potential,corr3,oneosqrteps,dlogeps,multp,&
     offset,geocode,lin_PB,PCGstart,CFgrid)

  use Poisson_Solver
  use yaml_output
  use f_utils
  use dynamic_memory
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden,iproc
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp,offset
  type(coulomb_operator), intent(inout) :: pkernel
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  integer, intent(in) :: SetEps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,corr3,oneosqrteps
  real(kind=8), dimension(3,n01,n02,n03), intent(in) :: dlogeps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: b
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: bb
  character(len=2), intent(in) :: geocode
  logical, intent(in) :: lin_PB,PCGstart,CFgrid

  real(kind=8), dimension(:,:,:,:), allocatable :: x,r,z,p,q,qold,lv,corr,deps,rvApp,lvc
  !real(kind=8), dimension(n01,n02,n03,3) :: deps
  real(kind=8), dimension(:,:,:), allocatable :: de2,ddeps,zcorr
  integer, parameter :: max_iter = 100
  real(kind=8), parameter :: max_ratioex = 1.0d10
  real(kind=8) :: alpha,beta,beta0,betanew,normb,normr,ratio,k,epsc,zeta,pval,qval,rval,pbval,multvar
  real(kind=8) :: eps1,eps2,eps3
  integer :: i,ii,j,i1,i2,i3,isp
  real(kind=8), parameter :: error = 1.0d-13
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8), dimension(n01,n02,n03) ::pot_ion
  real(kind=8) :: ehartree,pi,switch,rpoints,a
  real(kind=8) :: PB_charge,shift,offsetnew,prod

  !allocate heap arrays
  x=f_malloc([n01,n02,n03,nspden],id='x')
  r=f_malloc([n01,n02,n03,nspden],id='r')
  rvApp=f_malloc([n01,n02,n03,nspden],id='rvApp')
  z=f_malloc([n01,n02,n03,nspden],id='z')
  zcorr=f_malloc([n01,n02,n03],id='zcorr')
  p=f_malloc([n01,n02,n03,nspden],id='p')
  q=f_malloc([n01,n02,n03,nspden],id='q')
  qold=f_malloc([n01,n02,n03,nspden],id='qold')
  lv=f_malloc([n01,n02,n03,nspden],id='lv')
  lvc=f_malloc([n01/2,n02/2,n03/2,nspden],id='lvc')
  corr=f_malloc([n01,n02,n03,nspden],id='corr')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  ddeps=f_malloc([n01,n02,n03],id='ddeps')
  de2=f_malloc([n01,n02,n03],id='de2')

  pi = 4.d0*datan(1.d0)   
  rpoints=product(real([n01,n02,n03],kind=8))

  if (lin_PB) then
   open(unit=18,file='LinPB_PCG_normr_'//trim(geocode)//'.dat',status='unknown')
   open(unit=38,file='LinPB_PCG_accuracy_'//trim(geocode)//'.dat',status='unknown')
  else
   open(unit=18,file='PCG_normr_'//trim(geocode)//'.dat',status='unknown')
   open(unit=38,file='PCG_accuracy_'//trim(geocode)//'.dat',status='unknown')
  end if


  if (iproc ==0) then
   write(18,'(1x,a)')'iter normr ratio beta'
   write(38,'(1x,a)')'iter i1_max i2_max i3_max max_val max_center'
   call yaml_map('rpoints',rpoints)
   call yaml_sequence_open('Embedded PSolver, Preconditioned Conjugate Gradient Method')
  end if

  switch=0.0d0
  if (lin_PB) then
   switch=1.0d0
  end if

  call f_zero(rvApp)
  if (PCGstart) then

   ! Calculate the charge starting from the potential applying the proper Laplace operator.
   call ApplyLaplace(geocode,n01,n02,n03,nspden,hx,hy,hz,b,rvApp,acell,eps,nord,.false.,multp)

   if (iproc==0) then
    call yaml_mapping_open('Comparison between Generalized Poisson operator and analytical density')
    call writeroutinePot(n01,n02,n03,1,bb,0,rvApp)
    call yaml_mapping_close()
   end if

   ! Calculate the charge starting from the potential applying the proper Laplace operator.
   call ApplyLaplace2(geocode,n01,n02,n03,nspden,hx,hy,hz,b,rvApp,acell,eps,dlogeps,nord,.false.,multp)

  if (iproc==0) then
   call yaml_mapping_open('Comparison between Generalized Poisson operator 2 and analytical density')
   call writeroutinePot(n01,n02,n03,1,bb,1,rvApp)
   call yaml_mapping_close()
  end if

  end if

  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting Preconditioned Conjugate Gradient'
  end if

!------------------------------------------------------------------------------------
! Set the correction vector for the Generalized Laplace operator

!  call fssnordEpsilonDerivative(n01,n02,n03,nspden,hx,hy,hz,eps,de2,ddeps,nord,acell)

!  call fssnord3DmatNabla3varde2(n01,n02,n03,nspden,hx,hy,hx,eps,deps,de2,nord,acell)
!  call fssnord3DmatDiv3var(n01,n02,n03,nspden,hx,hy,hz,deps,ddeps,nord,acell)

!  isp=1
!  do i3=1,n03
!   do i2=1,n02
!    do i1=1,n01
!     corr(i1,i2,i3,isp)=(-0.125d0/pi)*(0.5d0*de2(i1,i2,i3)/eps(i1,i2,i3)-ddeps(i1,i2,i3))
!    end do
!   end do
!  end do


!------------------------------------------------------------------------------------
! Apply the Preconditioner

  normb=0.d0
  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     normb=normb+b(i1,i2,i3,isp)*b(i1,i2,i3,isp)
     !!lv(i1,i2,i3,isp) = b(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
     if (PCGstart) then
      x(i1,i2,i3,isp) = b(i1,i2,i3,isp)
      b(i1,i2,i3,isp) = bb(i1,i2,i3,isp)-rvApp(i1,i2,i3,isp)
     end if
    end do
   end do
  end do
!!  normb=dsqrt(normb)
  normb=sqrt(normb/rpoints)

!!$  call yaml_sequence(advance='no')
!!$  call H_potential('G',pkernel,lv,pot_ion,ehartree,offset,.false.)
!!$
!!$  isp=1
!!$  do i3=1,n03
!!$   do i2=1,n02
!!$    do i1=1,n01
!!$     p(i1,i2,i3,isp) = lv(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
!!$    end do
!!$   end do
!!$  end do
!!$
!!$!------------------------------------------------------------------------------------
!!$! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential correction
!!$
!!$  beta=0.d0
!!$  k=0.d0
!!$  isp=1
!!$  do i3=1,n03
!!$   do i2=1,n02
!!$    do i1=1,n01
!!$     q(i1,i2,i3,isp)=b(i1,i2,i3,isp)+p(i1,i2,i3,isp)*corr(i1,i2,i3,isp)
!!$     qold(i1,i2,i3,isp)=q(i1,i2,i3,isp)
!!$     beta=beta+b(i1,i2,i3,isp)*p(i1,i2,i3,isp)
!!$     k=k+p(i1,i2,i3,isp)*q(i1,i2,i3,isp)
!!$    end do
!!$   end do
!!$  end do
!!$
!!$!------------------------------------------------------------------------------------
!!$
!!$  alpha = beta/k
!!$  !write(*,*)alpha
!!$  normr=0.d0
!!$  isp=1
!!$  do i3=1,n03
!!$   do i2=1,n02
!!$    do i1=1,n01
!!$     x(i1,i2,i3,isp) = alpha*p(i1,i2,i3,isp)
!!$     r(i1,i2,i3,isp) = b(i1,i2,i3,isp) - alpha*q(i1,i2,i3,isp)
!!$     normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
!!$     lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
!!$    end do
!!$   end do
!!$  end do
!!$  normr=dsqrt(normr)
!!$
!!$  ratio=normr/normb
!!$
!!$  call EPS_iter_output(1,normb,normr,ratio,alpha,beta)
!!$!  call writeroutine(n01,n02,n03,nspden,r,1)
!!$  call writeroutinePot(n01,n02,n03,nspden,potential,0,potential)
!!$  call writeroutinePot(n01,n02,n03,nspden,x,1,potential)
!!$
!!$  write(18,'(1x,I8,2(1x,e14.7))')1,ratio,beta
  !write(*,'(1x,I8,2(1x,e14.7))')1,ratio,beta
  !initialization of the components
  call f_memcpy(src=b,dest=r)
  if (.not.PCGstart) call f_zero(x)
  call f_zero(q)
  call f_zero(p)
  beta=1.d0
  ratio=1.d0
  normr=normb
  shift=0.d0
  offsetnew=0.d0 

  ! Calculate the product in a finer grid by means of ISF.
!  call ISF_product(geocode,8,n01,n02,n03,r,oneosqrteps,zcorr,a)
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
           !lv(i1,i2,i3,isp) = pkernel%oneoeps(i1,i2,i3)*r(i1,i2,i3,isp)
           lv(i1,i2,i3,isp) = oneosqrteps(i1,i2,i3)*r(i1,i2,i3,isp)
!           lv(i1,i2,i3,isp) = zcorr(i1,i2,i3)
        end do
     end do
  end do

!  multvar=1.d0
!  if (iproc ==0) then
!   call yaml_map('iter',i)
!   call yaml_map('multvar',multvar)
!  end if

  do i=1,max_iter

   if (normr.lt.error) exit
   if (normr.gt.max_ratioex) exit

!   if ((modulo(i,20).eq.0) .and.(i.lt.61)) then
!    multvar=multvar*10.d0
!    if (iproc ==0) then
!     call yaml_map('iter',i)
!     call yaml_map('multvar',multvar)
!    end if
!   end if

   if (iproc==0) then
    write(*,'(a)')'--------------------------------------------------------------------------------------------!'
    write(*,*)'Starting PCG iteration ',i
   end if

!  Apply the Preconditioner

   if (iproc ==0) then
    call yaml_sequence(advance='no')
   end if
   
   if (CFgrid) then
    call for_trans_ISF_3D(geocode,8,n01,n02,n03,lv,lvc)
    call H_potential('G',pkernel,lvc,pot_ion,ehartree,0.0_dp,.false.)
    call back_trans_ISF_3D(geocode,8,n01,n02,n03,lvc,lv)
   else
    call H_potential('G',pkernel,lv,pot_ion,ehartree,0.0_dp,.false.)
   end if

   beta0 = beta
   beta=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
        z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
        beta=beta+r(i1,i2,i3,isp)*z(i1,i2,i3,isp)
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential correction
      !q(i1,i2,i3,isp)=r(i1,i2,i3,isp)+z(i1,i2,i3,isp)*corr(i1,i2,i3,isp)
     end do
    end do
   end do
!  call ISF_product_new(geocode,8,n01,n02,n03,hx,hy,hz,nord,acell,r,z,zcorr,beta)

  ! Calculate the product in a finer grid by means of ISF.
  !call ISF_product(geocode,8,n01,n02,n03,z,corr3,zcorr)
!  call ISF_product_new(geocode,8,n01,n02,n03,hx,hy,hz,nord,acell,z,corr3,zcorr,a)

  k=0.d0
  isp=1
  do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        zeta=z(i1,i2,i3,isp)
        !epsc=corr(i1,i2,i3,isp)
        !epsc=pkernel%corr(i1,i2,i3)
        epsc=corr3(i1,i2,i3)
        pval=p(i1,i2,i3,isp)
        qval=q(i1,i2,i3,isp)
        rval=r(i1,i2,i3,isp)
        pval = zeta+(beta/beta0)*pval
        ! Additional contribution to the Generalized Poisson operator
        pbval=-switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(zeta) 
                                                                          ! for the Poisson-Boltzmann solution.
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dsinh(multp*zeta)
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*zeta*dcosh(multp*x(i1,i2,i3,isp))
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dtanh(multp*zeta)
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*zeta
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*(zeta**2)
        !qval = zeta*epsc+rval+pbval+(beta/beta0)*qval
        prod=zeta*epsc
        !prod=zcorr(i1,i2,i3)
        qval = prod+rval+pbval+(beta/beta0)*qval
        k = k + pval*qval
        p(i1,i2,i3,isp) = pval
        q(i1,i2,i3,isp) = qval
        !p(i1,i2,i3,isp) = z(i1,i2,i3,isp)+(beta/beta0)*p(i1,i2,i3,isp)
        !q(i1,i2,i3,isp) = q(i1,i2,i3,isp)+(beta/beta0)*qold(i1,i2,i3,isp)
        !qold(i1,i2,i3,isp)=q(i1,i2,i3,isp)
        !k=k+p(i1,i2,i3,isp)*q(i1,i2,i3,isp)
     end do
    end do
   end do
!   call ISF_product_new(geocode,8,n01,n02,n03,hx,hy,hz,nord,acell,p,q,zcorr,k)

   alpha = beta/k
   !write(*,*)alpha

   offsetnew=0.d0
   normr=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      x(i1,i2,i3,isp) = x(i1,i2,i3,isp) + alpha*p(i1,i2,i3,isp)
      offsetnew=offsetnew+alpha*p(i1,i2,i3,isp)
      r(i1,i2,i3,isp) = r(i1,i2,i3,isp) - alpha*q(i1,i2,i3,isp)
      normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
      lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
     end do
    end do
   end do

  ! Calculate the product in a finer grid by means of ISF.
!  call ISF_product(geocode,8,n01,n02,n03,r,oneosqrteps,zcorr,a)
!   do i3=1,n03
!    do i2=1,n02
!     do i1=1,n01
!      lv(i1,i2,i3,isp) = zcorr(i1,i2,i3)
!     end do
!    end do
!   end do

!   normr=dsqrt(normr)
   normr=sqrt(normr/rpoints)

   if (trim(geocode) == 'P') then

    if (i.eq.1) then
     shift=offsetnew/(n01*n02*n03)-offset/(acell**3)
    else if (i.gt.1) then
     shift=offsetnew/(n01*n02*n03)
    end if

    isp=1
    do i3=1,n03
     do i2=1,n02
      do i1=1,n01
       x(i1,i2,i3,isp) = x(i1,i2,i3,isp) - shift
      end do
     end do
    end do

   end if

   ratio=normr/normb
   if (iproc ==0) then
   write(18,'(1x,I8,3(1x,e14.7))')i,normr,ratio,beta
   !write(*,'(1x,I8,2(1x,e14.7))')i,ratio,beta
   call EPS_iter_output_LG(i,normb,normr,ratio,alpha,beta)
!   call writeroutine(n01,n02,n03,nspden,r,i)
   call writeroutinePot(n01,n02,n03,nspden,x,i,potential)
   end if

  end do

   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      b(i1,i2,i3,isp) = x(i1,i2,i3,isp)
     end do
    end do
   end do

  call yaml_sequence_close()
   !write(*,*)
   !write(*,'(1x,a,1x,I8)')'PCG iterations =',i-1
   !write(*,'(1x,a,1x,e14.7)')'PCG error =',ratio
   !write(*,*)
   !write(*,*)'Max abs difference between analytic potential and the computed one'
!  if (iproc==0) then
!   call writeroutinePot(n01,n02,n03,nspden,b,i-1,potential)
!   write(*,*)
!  end if

  close(unit=18)
  close(unit=38)

  if (iproc==0) then
   write(*,'(a)')'Termination of Preconditioned Conjugate Gradient'
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
  end if

  call f_free(x)
  call f_free(r)
  call f_free(rvApp)
  call f_free(z)
  call f_free(zcorr)
  call f_free(p)
  call f_free(q)
  call f_free(qold)
  call f_free(lv)
  call f_free(corr)
  call f_free(deps)
  call f_free(ddeps)
  call f_free(de2)
  call f_free(lvc)

end subroutine  Prec_conjugate_gradient

subroutine Prec_conjugate_gradient_Inputguess(n01,n02,n03,nspden,iproc,hx,hy,hz,b,bb,&
     acell,eps,SetEps,nord,pkernel,potential,corr3,oneosqrteps,dlogeps,multp,offset,geocode,lin_PB)

  use Poisson_Solver
  use yaml_output
  use f_utils
  use dynamic_memory
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden,iproc
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp,offset
  type(coulomb_operator), intent(inout) :: pkernel
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  integer, intent(in) :: SetEps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,corr3,oneosqrteps
  real(kind=8), dimension(3,n01,n02,n03), intent(in) :: dlogeps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: b
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: bb
  character(len=2), intent(in) :: geocode
  logical, intent(in) :: lin_PB

  real(kind=8), dimension(:,:,:,:), allocatable :: x,r,z,p,q,qold,lv,corr,deps,rvApp,r0
  !real(kind=8), dimension(n01,n02,n03,3) :: deps
  real(kind=8), dimension(:,:,:), allocatable :: de2,ddeps,zcorr
  integer, parameter :: max_iter = 100
  integer, parameter :: max_scf = 1
  real(kind=8), parameter :: max_ratioex = 1.0d10
  real(kind=8) :: alpha,beta,beta0,betanew,normb,normr,ratio,k,epsc,zeta,pval,qval,rval,pbval,multvar
  real(kind=8) :: eps1,eps2,eps3
  integer :: i,ii,j,i1,i2,i3,isp,i_scf,inn
  real(kind=8), parameter :: error = 1.0d-13
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8), dimension(n01,n02,n03) ::pot_ion
  real(kind=8) :: ehartree,pi,switch,rpoints,a,errorcur
  real(kind=8) :: PB_charge,shift,offsetnew,prod
  logical :: PCGig

  !allocate heap arrays
  x=f_malloc([n01,n02,n03,nspden],id='x')
  r=f_malloc([n01,n02,n03,nspden],id='r')
  r0=f_malloc([n01,n02,n03,nspden],id='r0')
  rvApp=f_malloc([n01,n02,n03,nspden],id='rvApp')
  z=f_malloc([n01,n02,n03,nspden],id='z')
  zcorr=f_malloc([n01,n02,n03],id='zcorr')
  p=f_malloc([n01,n02,n03,nspden],id='p')
  q=f_malloc([n01,n02,n03,nspden],id='q')
  qold=f_malloc([n01,n02,n03,nspden],id='qold')
  lv=f_malloc([n01,n02,n03,nspden],id='lv')
  corr=f_malloc([n01,n02,n03,nspden],id='corr')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  ddeps=f_malloc([n01,n02,n03],id='ddeps')
  de2=f_malloc([n01,n02,n03],id='de2')

  pi = 4.d0*datan(1.d0)   
  rpoints=product(real([n01,n02,n03],kind=8))

  if (lin_PB) then
   open(unit=18,file='LinPB_PCG_normr_'//trim(geocode)//'.dat',status='unknown')
   open(unit=38,file='LinPB_PCG_accuracy_'//trim(geocode)//'.dat',status='unknown')
  else
   open(unit=18,file='PCG_normr_'//trim(geocode)//'.dat',status='unknown')
   open(unit=38,file='PCG_accuracy_'//trim(geocode)//'.dat',status='unknown')
  end if


  if (iproc ==0) then
   write(18,'(1x,a)')'iter normr ratio beta'
   write(38,'(1x,a)')'iter i1_max i2_max i3_max max_val max_center'
   call yaml_map('rpoints',rpoints)
   call yaml_sequence_open('Embedded PSolver, Preconditioned Conjugate Gradient Method')
  end if

  switch=0.0d0
  if (lin_PB) then
   switch=1.0d0
  end if


  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting an Inputguess PCG calculation'
  end if

  do i_scf=1,max_scf ! scf loop external

  if (i_scf.eq.1) then
   errorcur=error
   PCGig=.true.
  else
   !PCGig=.false.
   PCGig=.true.
   !errorcur=error
   errorcur=error
  end if

  call f_zero(rvApp)
  if (PCGig) then

   call f_memcpy(src=b,dest=x)

   ! Calculate the charge starting from the potential applying the proper Laplace operator.
   call ApplyLaplace(geocode,n01,n02,n03,nspden,hx,hy,hz,x,rvApp,acell,eps,nord,.true.,multp)

   if (iproc==0) then
    call yaml_mapping_open('Comparison between Generalized Poisson operator and analytical density')
    call writeroutinePot(n01,n02,n03,1,bb,0,rvApp)
    call yaml_mapping_close()
   end if

!   ! Calculate the charge starting from the potential applying the proper Laplace operator.
!   call ApplyLaplace2(geocode,n01,n02,n03,nspden,hx,hy,hz,x,rvApp,acell,eps,dlogeps,nord,.false.,multp)
!
!  if (iproc==0) then
!   call yaml_comment('Comparison between Generalized Poisson operator 2 and analytical density',hfill='-')
!   call writeroutinePot(n01,n02,n03,1,bb,1,rvApp)
!  end if

  end if

  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------!'
   write(*,*)'Starting Inputguess iteration ',i_scf
  end if

  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting Preconditioned Conjugate Gradient'
  end if

  normb=0.d0
  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     !!lv(i1,i2,i3,isp) = b(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
     if (PCGig) then
      !b(i1,i2,i3,isp) = bb(i1,i2,i3,isp)-rvApp(i1,i2,i3,isp)
      r0(i1,i2,i3,isp) = x(i1,i2,i3,isp)*corr3(i1,i2,i3)
      b(i1,i2,i3,isp) = bb(i1,i2,i3,isp) - r0(i1,i2,i3,isp)
     end if
     normb=normb+b(i1,i2,i3,isp)*b(i1,i2,i3,isp)
    end do
   end do
  end do
!!  normb=dsqrt(normb)
  normb=sqrt(normb/rpoints)

  call f_memcpy(src=b,dest=r)
  if (.not.PCGig) call f_zero(x)
  if (.not.PCGig) call f_zero(q)
  if (.not.PCGig) call f_zero(p)
  call f_zero(x)
  call f_zero(q)
  call f_zero(p)
  !call f_memcpy(src=b,dest=p)

  beta=1.d0
  ratio=1.d0
  normr=normb
  shift=0.d0
  offsetnew=0.d0 

  ! Calculate the product in a finer grid by means of ISF.
!  call ISF_product(geocode,8,n01,n02,n03,r,oneosqrteps,zcorr,a)
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
           !lv(i1,i2,i3,isp) = pkernel%oneoeps(i1,i2,i3)*r(i1,i2,i3,isp)
           lv(i1,i2,i3,isp) = oneosqrteps(i1,i2,i3)*r(i1,i2,i3,isp)
!           lv(i1,i2,i3,isp) = zcorr(i1,i2,i3)
        end do
     end do
  end do

  inn=1
  if (PCGig) then
   i=1
   inn=2

   call H_potential('G',pkernel,lv,pot_ion,ehartree,0.0_dp,.false.)

   beta0 = beta
   beta=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
        z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
        beta=beta+r(i1,i2,i3,isp)*z(i1,i2,i3,isp)
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential correction
      !q(i1,i2,i3,isp)=r(i1,i2,i3,isp)+z(i1,i2,i3,isp)*corr(i1,i2,i3,isp)
     end do
    end do
   end do

  k=0.d0
  isp=1
  do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        zeta=z(i1,i2,i3,isp)
        !epsc=corr(i1,i2,i3,isp)
        !epsc=pkernel%corr(i1,i2,i3)
        epsc=corr3(i1,i2,i3)
        pval=p(i1,i2,i3,isp)
        qval=q(i1,i2,i3,isp)
        rval=r(i1,i2,i3,isp)
        pval = zeta+(beta/beta0)*pval
        ! Additional contribution to the Generalized Poisson operator
        !pbval=-switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(zeta) 
        prod=zeta*epsc
        qval = prod !+rval+pbval+(beta/beta0)*qval
        !k = k + pval*(qval+rval)
        k = k + pval*(qval+rval)
        p(i1,i2,i3,isp) = pval
        q(i1,i2,i3,isp) = qval
     end do
    end do
   end do

   alpha = beta/k
   !write(*,*)alpha
   !alpha=1.d0
!   if ((PCGig).and.(i.eq.1)) alpha = 1.d0

   offsetnew=0.d0
   normr=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      !x(i1,i2,i3,isp) = x(i1,i2,i3,isp) + alpha*p(i1,i2,i3,isp)
      x(i1,i2,i3,isp) = p(i1,i2,i3,isp)
      offsetnew=offsetnew+alpha*p(i1,i2,i3,isp)
      !r(i1,i2,i3,isp) = alpha*(r0(i1,i2,i3,isp) - q(i1,i2,i3,isp))
      r(i1,i2,i3,isp) = r0(i1,i2,i3,isp) - q(i1,i2,i3,isp)
      !r(i1,i2,i3,isp) = r(i1,i2,i3,isp) - alpha*q(i1,i2,i3,isp)
      normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
      lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
      p(i1,i2,i3,isp) = 0.d0
      q(i1,i2,i3,isp) = 0.d0
     end do
    end do
   end do
  beta=1.d0
  ratio=1.d0

!   normr=dsqrt(normr)
   normr=sqrt(normr/rpoints)

   if (trim(geocode) == 'P') then

    if (i.eq.1) then
     shift=offsetnew/(n01*n02*n03)-offset/(acell**3)
    else if (i.gt.1) then
     shift=offsetnew/(n01*n02*n03)
    end if

    isp=1
    do i3=1,n03
     do i2=1,n02
      do i1=1,n01
       x(i1,i2,i3,isp) = x(i1,i2,i3,isp) - shift
      end do
     end do
    end do

   end if

   ratio=normr/normb
   if (iproc ==0) then
   write(18,'(1x,I8,3(1x,e14.7))')i,normr,ratio,beta
   !write(*,'(1x,I8,2(1x,e14.7))')i,ratio,beta
   call EPS_iter_output_LG(i,normb,normr,ratio,alpha,beta)
!   call writeroutine(n01,n02,n03,nspden,r,i)
   call writeroutinePot(n01,n02,n03,nspden,x,i,potential)
   end if

  end if
!  if (PCGig) stop
  do i=inn,max_iter

   if (normr.lt.errorcur) exit
   if (normr.gt.max_ratioex) exit

   if (iproc==0) then
    write(*,'(a)')'--------------------------------------------------------------------------------------------!'
    write(*,*)'Starting PCG iteration ',i
   end if


   if (iproc ==0) then
    call yaml_sequence(advance='no')
   end if
   
    call H_potential('G',pkernel,lv,pot_ion,ehartree,0.0_dp,.false.)

   beta0 = beta
   beta=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
        z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
        beta=beta+r(i1,i2,i3,isp)*z(i1,i2,i3,isp)
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential correction
      !q(i1,i2,i3,isp)=r(i1,i2,i3,isp)+z(i1,i2,i3,isp)*corr(i1,i2,i3,isp)
     end do
    end do
   end do

  k=0.d0
  isp=1
  do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        zeta=z(i1,i2,i3,isp)
        !epsc=corr(i1,i2,i3,isp)
        !epsc=pkernel%corr(i1,i2,i3)
        epsc=corr3(i1,i2,i3)
        pval=p(i1,i2,i3,isp)
        qval=q(i1,i2,i3,isp)
        rval=r(i1,i2,i3,isp)
        pval = zeta+(beta/beta0)*pval
        ! Additional contribution to the Generalized Poisson operator
        pbval=-switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(zeta) 
        pbval=0.d0
        prod=zeta*epsc
        qval = prod+rval+pbval+(beta/beta0)*qval
        k = k + pval*qval
        p(i1,i2,i3,isp) = pval
        q(i1,i2,i3,isp) = qval
     end do
    end do
   end do

   alpha = beta/k
   !write(*,*)alpha

!   if ((PCGig).and.(i.eq.1)) alpha = 1.d0

   offsetnew=0.d0
   normr=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      !x(i1,i2,i3,isp) = x(i1,i2,i3,isp) + alpha*p(i1,i2,i3,isp)
      x(i1,i2,i3,isp) = x(i1,i2,i3,isp) + alpha*p(i1,i2,i3,isp)
      offsetnew=offsetnew+alpha*p(i1,i2,i3,isp)
      r(i1,i2,i3,isp) = r(i1,i2,i3,isp) - alpha*q(i1,i2,i3,isp)
      normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
      lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
     end do
    end do
   end do

!   normr=dsqrt(normr)
   normr=sqrt(normr/rpoints)

   if (trim(geocode) == 'P') then

    if (i.eq.1) then
     shift=offsetnew/(n01*n02*n03)-offset/(acell**3)
    else if (i.gt.1) then
     shift=offsetnew/(n01*n02*n03)
    end if

    isp=1
    do i3=1,n03
     do i2=1,n02
      do i1=1,n01
       x(i1,i2,i3,isp) = x(i1,i2,i3,isp) - shift
      end do
     end do
    end do

   end if

   ratio=normr/normb
   if (iproc ==0) then
   write(18,'(1x,I8,3(1x,e14.7))')i,normr,ratio,beta
   !write(*,'(1x,I8,2(1x,e14.7))')i,ratio,beta
   call EPS_iter_output_LG(i,normb,normr,ratio,alpha,beta)
!   call writeroutine(n01,n02,n03,nspden,r,i)
   call writeroutinePot(n01,n02,n03,nspden,x,i,potential)
   end if

  end do ! PCG loop

  if (iproc==0) then
   write(*,'(a)')'Termination of Preconditioned Conjugate Gradient'
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
  end if

  end do !scf loop

  call yaml_sequence_close()
   !write(*,*)
   !write(*,'(1x,a,1x,I8)')'PCG iterations =',i-1
   !write(*,'(1x,a,1x,e14.7)')'PCG error =',ratio
   !write(*,*)
   !write(*,*)'Max abs difference between analytic potential and the computed one'
!  if (iproc==0) then
!   call writeroutinePot(n01,n02,n03,nspden,b,i-1,potential)
!   write(*,*)
!  end if

  close(unit=18)
  close(unit=38)
  
  call f_memcpy(src=x,dest=b)

  if (iproc==0) then
   write(*,'(a)')'Termination of restart loop'
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
  end if

  call f_free(x)
  call f_free(r)
  call f_free(r0)
  call f_free(rvApp)
  call f_free(z)
  call f_free(zcorr)
  call f_free(p)
  call f_free(q)
  call f_free(qold)
  call f_free(lv)
  call f_free(corr)
  call f_free(deps)
  call f_free(ddeps)
  call f_free(de2)

end subroutine  Prec_conjugate_gradient_Inputguess

subroutine Prec_conjugate_gradient_PSD(n01,n02,n03,nspden,iproc,hx,hy,hz,b,&
     acell,eps,SetEps,nord,pkernel,potential,corr3,oneosqrteps,oneoeps,dlogeps,multp,offset,geocode,lin_PB)

  use Poisson_Solver
  use yaml_output
  use f_utils
  use dynamic_memory
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden,iproc
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp,offset
  type(coulomb_operator), intent(inout) :: pkernel
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  integer, intent(in) :: SetEps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,corr3,oneosqrteps,oneoeps
  real(kind=8), dimension(3,n01,n02,n03), intent(in) :: dlogeps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: b
  character(len=2), intent(in) :: geocode
  logical, intent(in) :: lin_PB

  real(kind=8), dimension(:,:,:,:), allocatable :: x,r,z,p,q,qold,lv,corr,deps
  real(kind=8), dimension(:,:,:,:,:), allocatable :: dx
  !real(kind=8), dimension(n01,n02,n03,3) :: deps
  real(kind=8), dimension(:,:,:), allocatable :: de2,ddeps
  integer, parameter :: max_iter = 100
  real(kind=8), parameter :: max_ratioex = 1.0d10
  real(kind=8), parameter :: eta = 0.6d0 !1.0d0 !Like eta in the Polarization Iterative Method.
  real(kind=8) :: alpha,beta,beta0,betanew,normb,normr,ratio,k,epsc,zeta,pval,qval,rval,pbval,multvar
  real(kind=8) :: betak,divprod
  integer :: i,ii,j,i1,i2,i3,isp
  real(kind=8), parameter :: error = 1.0d-12
  real(kind=8), parameter :: PSDstart = 1.0d0
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8), dimension(n01,n02,n03) ::pot_ion
  real(kind=8) :: ehartree,pi,switch,rpoints
  real(kind=8) :: PB_charge,shift,offsetnew
  logical  :: PSDloop = .true.

  !allocate heap arrays
  x=f_malloc([n01,n02,n03,nspden],id='x')
  dx=f_malloc([n01,n02,n03,nspden,3],id='dx')
  r=f_malloc([n01,n02,n03,nspden],id='r')
  z=f_malloc([n01,n02,n03,nspden],id='z')
  p=f_malloc([n01,n02,n03,nspden],id='p')
  q=f_malloc([n01,n02,n03,nspden],id='q')
  qold=f_malloc([n01,n02,n03,nspden],id='qold')
  lv=f_malloc([n01,n02,n03,nspden],id='lv')
  corr=f_malloc([n01,n02,n03,nspden],id='corr')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  ddeps=f_malloc([n01,n02,n03],id='ddeps')
  de2=f_malloc([n01,n02,n03],id='de2')
 write(*,*)'im here'

  pi = 4.d0*datan(1.d0)   
  rpoints=product(real([n01,n02,n03],kind=8))

  if (lin_PB) then
   open(unit=18,file='LinPB_PCG_normr_'//trim(geocode)//'.dat',status='unknown')
   open(unit=38,file='LinPB_PCG_accuracy_'//trim(geocode)//'.dat',status='unknown')
  else
   open(unit=18,file='PCG_normr_'//trim(geocode)//'.dat',status='unknown')
   open(unit=38,file='PCG_accuracy_'//trim(geocode)//'.dat',status='unknown')
  end if


  if (iproc ==0) then
   write(18,'(1x,a)')'iter normr ratio beta'
   write(38,'(1x,a)')'iter i1_max i2_max i3_max max_val max_center'
   call yaml_map('rpoints',rpoints)
   call yaml_sequence_open('Embedded PSolver, Preconditioned Conjugate Gradient Method')
  end if

  switch=0.0d0
  if (lin_PB) then
   switch=1.0d0
  end if

  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting Preconditioned Conjugate Gradient'
  end if

  normb=0.d0
  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     normb=normb+b(i1,i2,i3,isp)*b(i1,i2,i3,isp)
     !!lv(i1,i2,i3,isp) = b(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
    end do
   end do
  end do
!!  normb=dsqrt(normb)
  normb=sqrt(normb/rpoints)

  call f_memcpy(src=b,dest=r)
  call f_zero(x)
  call f_zero(q)
  call f_zero(p)
  beta=1.d0
  betak=1.d0
  ratio=1.d0
  normr=normb
  shift=0.d0
  offsetnew=0.d0 

  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
           !lv(i1,i2,i3,isp) = pkernel%oneoeps(i1,i2,i3)*r(i1,i2,i3,isp)
           if (PSDloop) then
            lv(i1,i2,i3,isp) = oneoeps(i1,i2,i3)*r(i1,i2,i3,isp)
           else
            lv(i1,i2,i3,isp) = oneosqrteps(i1,i2,i3)*r(i1,i2,i3,isp)
           end if
        end do
     end do
  end do

!  multvar=1.d0
!  if (iproc ==0) then
!   call yaml_map('iter',i)
!   call yaml_map('multvar',multvar)
!  end if

  do i=1,max_iter

   if (normr.lt.error) exit
   if (normr.gt.max_ratioex) exit

!   if ((modulo(i,20).eq.0) .and.(i.lt.61)) then
!    multvar=multvar*10.d0
!    if (iproc ==0) then
!     call yaml_map('iter',i)
!     call yaml_map('multvar',multvar)
!    end if
!   end if

   if (iproc==0) then
    write(*,'(a)')'--------------------------------------------------------------------------------------------!'
    if (PSDloop) then
     write(*,*)'Starting PSD iteration ',i
    else
     write(*,*)'Starting PCG iteration ',i
    end if
   end if

!  Apply the Preconditioner

   if (iproc ==0) then
    call yaml_sequence(advance='no')
   end if

    call H_potential('G',pkernel,lv,pot_ion,ehartree,0.0_dp,.false.)

    if (PSDloop) then
     call fssnord3DmatNabla(geocode,n01,n02,n03,nspden,hx,hy,hz,lv,dx,nord,acell)
    end if

   beta0 = beta
   beta=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
           if (PSDloop) then
            z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)!*oneoeps(i1,i2,i3)
           else
            z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
           end if
        beta=beta+r(i1,i2,i3,isp)*z(i1,i2,i3,isp)
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential correction
      !q(i1,i2,i3,isp)=r(i1,i2,i3,isp)+z(i1,i2,i3,isp)*corr(i1,i2,i3,isp)
     end do
    end do
   end do


!   if (PSDloop) betak = 0.d0
   k=0.d0
   isp=1

   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        zeta=z(i1,i2,i3,isp)
        !epsc=corr(i1,i2,i3,isp)
        !epsc=pkernel%corr(i1,i2,i3)
        epsc=corr3(i1,i2,i3)
        pval=p(i1,i2,i3,isp)
        qval=q(i1,i2,i3,isp)
        rval=r(i1,i2,i3,isp)
        if (PSDloop) then
         pval = zeta
         divprod = 0.d0
         do j=1,3
         divprod = divprod + dlogeps(j,i1,i2,i3)*dx(i1,i2,i3,isp,j)
         end do
         qval = -0.25d0*divprod*eps(i1,i2,i3)/pi + rval
        else
         pval = zeta+betak*(beta/beta0)*pval
         ! Additional contribution to the Generalized Poisson operator
         pbval=-switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(zeta) 
                                                                          ! for the Poisson-Boltzmann solution.
         qval = zeta*epsc+rval+pbval+betak*(beta/beta0)*qval
        end if
!        qval = zeta*epsc+rval+(beta/beta0)*qval
        k = k + pval*qval
        p(i1,i2,i3,isp) = pval
        q(i1,i2,i3,isp) = qval
        !p(i1,i2,i3,isp) = z(i1,i2,i3,isp)+(beta/beta0)*p(i1,i2,i3,isp)
        !q(i1,i2,i3,isp) = q(i1,i2,i3,isp)+(beta/beta0)*qold(i1,i2,i3,isp)
        !qold(i1,i2,i3,isp)=q(i1,i2,i3,isp)
        !k=k+p(i1,i2,i3,isp)*q(i1,i2,i3,isp)
     end do
    end do
   end do

   alpha = beta/k
   if (PSDloop) alpha = eta
   write(*,*)alpha

   offsetnew=0.d0
   normr=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      x(i1,i2,i3,isp) = x(i1,i2,i3,isp) + alpha*p(i1,i2,i3,isp)
      offsetnew=offsetnew+alpha*p(i1,i2,i3,isp)
      r(i1,i2,i3,isp) = r(i1,i2,i3,isp) - alpha*q(i1,i2,i3,isp)
      normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
      if (.not.PSDloop) lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
     end do
    end do
   end do
!   normr=dsqrt(normr)
   normr=sqrt(normr/rpoints)

   if (normr.lt.PSDstart) then
    PSDloop=.true.
    isp=1
    do i3=1,n03
     do i2=1,n02
      do i1=1,n01
        lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*oneoeps(i1,i2,i3)
      end do
     end do
    end do
   end if

   if (trim(geocode) == 'P') then

    if (i.eq.1) then
     shift=offsetnew/(n01*n02*n03)-offset/(acell**3)
    else if (i.gt.1) then
     shift=offsetnew/(n01*n02*n03)
    end if

    isp=1
    do i3=1,n03
     do i2=1,n02
      do i1=1,n01
       x(i1,i2,i3,isp) = x(i1,i2,i3,isp) - shift
      end do
     end do
    end do

   end if

   ratio=normr/normb
   if (iproc ==0) then
   write(18,'(1x,I8,3(1x,e14.7))')i,normr,ratio,beta
   !write(*,'(1x,I8,2(1x,e14.7))')i,ratio,beta
   call EPS_iter_output_LG(i,normb,normr,ratio,alpha,beta)
!   call writeroutine(n01,n02,n03,nspden,r,i)
   call writeroutinePot(n01,n02,n03,nspden,x,i,potential)
   end if

  end do

   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      b(i1,i2,i3,isp) = x(i1,i2,i3,isp)
     end do
    end do
   end do

  call yaml_sequence_close()
   !write(*,*)
   !write(*,'(1x,a,1x,I8)')'PCG iterations =',i-1
   !write(*,'(1x,a,1x,e14.7)')'PCG error =',ratio
   !write(*,*)
   !write(*,*)'Max abs difference between analytic potential and the computed one'
!  if (iproc==0) then
!   call writeroutinePot(n01,n02,n03,nspden,b,i-1,potential)
!   write(*,*)
!  end if

  close(unit=18)
  close(unit=38)

  if (iproc==0) then
   write(*,'(a)')'Termination of Preconditioned Conjugate Gradient'
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
  end if
  call f_free(x)
  call f_free(dx)
  call f_free(r)
  call f_free(z)
  call f_free(p)
  call f_free(q)
  call f_free(qold)
  call f_free(lv)
  call f_free(corr)
  call f_free(deps)
  call f_free(ddeps)
  call f_free(de2)

end subroutine  Prec_conjugate_gradient_PSD

subroutine Prec_conjugate_gradient_restart(n01,n02,n03,nspden,iproc,hx,hy,hz,b,&
     acell,eps,SetEps,nord,pkernel,potential,corr3,oneosqrteps,multp,offset,geocode,lin_PB)

  use Poisson_Solver
  use yaml_output
  use f_utils
  use dynamic_memory
  use wrapper_linalg
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden,iproc
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp,offset
  type(coulomb_operator), intent(inout) :: pkernel
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  integer, intent(in) :: SetEps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,corr3,oneosqrteps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: b
  character(len=2), intent(in) :: geocode
  logical, intent(in) :: lin_PB

  real(kind=8), dimension(:,:,:,:), allocatable :: x,r,z,p,q,qold,lv,corr,deps,b_old
  !real(kind=8), dimension(n01,n02,n03,3) :: deps
  real(kind=8), dimension(:,:,:), allocatable :: de2,ddeps
  integer, parameter :: max_iter = 100
  integer, parameter :: max_iter_scf = 2
  real(kind=8), parameter :: max_ratioex = 1.0d10
  real(kind=8) :: alpha,beta,beta0,betanew,normb,normr,ratio,k,epsc,zeta,pval,qval,rval,pbval,multvar
  integer :: i,ii,j,i1,i2,i3,isp,i_scf
  real(kind=8), parameter :: error = 1.0d-8
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8) :: scal 
  real(kind=8), dimension(n01,n02,n03) ::pot_ion
  real(kind=8) :: ehartree,pi,switch,rpoints
  real(kind=8) :: PB_charge,shift,offsetnew

  !allocate heap arrays
  x=f_malloc([n01,n02,n03,nspden],id='x')
  r=f_malloc([n01,n02,n03,nspden],id='r')
  z=f_malloc([n01,n02,n03,nspden],id='z')
  p=f_malloc([n01,n02,n03,nspden],id='p')
  q=f_malloc([n01,n02,n03,nspden],id='q')
  qold=f_malloc([n01,n02,n03,nspden],id='qold')
  lv=f_malloc([n01,n02,n03,nspden],id='lv')
  corr=f_malloc([n01,n02,n03,nspden],id='corr')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  ddeps=f_malloc([n01,n02,n03],id='ddeps')
  de2=f_malloc([n01,n02,n03],id='de2')
  b_old=f_malloc([n01,n02,n03,nspden],id='b_old')

  pi = 4.d0*datan(1.d0)   
  rpoints=product(real([n01,n02,n03],kind=8))

  if (lin_PB) then
   open(unit=18,file='LinPB_PCG_normr_'//trim(geocode)//'.dat',status='unknown')
   open(unit=38,file='LinPB_PCG_accuracy_'//trim(geocode)//'.dat',status='unknown')
  else
   open(unit=18,file='PCG_normr_'//trim(geocode)//'.dat',status='unknown')
   open(unit=38,file='PCG_accuracy_'//trim(geocode)//'.dat',status='unknown')
  end if


  if (iproc ==0) then
   write(18,'(1x,a)')'iter normr ratio beta'
   write(38,'(1x,a)')'iter i1_max i2_max i3_max max_val max_center'
   call yaml_map('rpoints',rpoints)
   call yaml_sequence_open('Embedded PSolver, Preconditioned Conjugate Gradient Method')
  end if

  switch=0.0d0
  if (lin_PB) then
   switch=1.0d0
  end if
  
  normb=0.d0
  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     normb=normb+b(i1,i2,i3,isp)*b(i1,i2,i3,isp)
     !!lv(i1,i2,i3,isp) = b(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
    end do
   end do
  end do
!!  normb=dsqrt(normb)
  normb=sqrt(normb/rpoints)

  do i_scf=1,max_iter_scf ! restart loop.

   if (iproc==0) then
    write(*,'(a)')'--------------------------------------------------------------------------------------------!'
    write(*,*)'Starting GPS ',i_scf
   end if

  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting Preconditioned Conjugate Gradient'
  end if

!------------------------------------------------------------------------------------
! Set the correction vector for the Generalized Laplace operator

!  call fssnordEpsilonDerivative(n01,n02,n03,nspden,hx,hy,hz,eps,de2,ddeps,nord,acell)

!  call fssnord3DmatNabla3varde2(n01,n02,n03,nspden,hx,hy,hx,eps,deps,de2,nord,acell)
!  call fssnord3DmatDiv3var(n01,n02,n03,nspden,hx,hy,hz,deps,ddeps,nord,acell)

!  isp=1
!  do i3=1,n03
!   do i2=1,n02
!    do i1=1,n01
!     corr(i1,i2,i3,isp)=(-0.125d0/pi)*(0.5d0*de2(i1,i2,i3)/eps(i1,i2,i3)-ddeps(i1,i2,i3))
!    end do
!   end do
!  end do


!------------------------------------------------------------------------------------
! Apply the Preconditioner


  if (i_scf==1) then
   scal=1.d0
   call f_memcpy(src=b,dest=r)
   call f_zero(x)
!   call f_zero(q)
!   call f_zero(p)
  else
   scal=1.0d0
    b=b+30.d0
   call vscal(size(r),scal,b(1,1,1,1),1)
   call f_memcpy(src=b,dest=r)
   call f_zero(x)
!   call f_zero(r)
!   call axpy(size(r),-1.0_gp,b_old(1,1,1,1),1,r(1,1,1,1),1)
!   call axpy(size(r),1.0_gp,b(1,1,1,1),1,r(1,1,1,1),1)
  end if
  call f_memcpy(src=b,dest=b_old)
  
  call f_zero(q)
  call f_zero(p)
  beta=1.d0
  ratio=1.d0
  shift=0.d0
  offsetnew=0.d0 

  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
           !lv(i1,i2,i3,isp) = pkernel%oneoeps(i1,i2,i3)*r(i1,i2,i3,isp)
           lv(i1,i2,i3,isp) = oneosqrteps(i1,i2,i3)*r(i1,i2,i3,isp)
        end do
     end do
  end do

!  multvar=1.d0
!  if (iproc ==0) then
!   call yaml_map('iter',i)
!   call yaml_map('multvar',multvar)
!  end if

  do i=1,max_iter

!   if (normr.lt.error) exit
!   if (normr.gt.max_ratioex) exit

!   if ((modulo(i,20).eq.0) .and.(i.lt.61)) then
!    multvar=multvar*10.d0
!    if (iproc ==0) then
!     call yaml_map('iter',i)
!     call yaml_map('multvar',multvar)
!    end if
!   end if

   if (iproc==0) then
    write(*,'(a)')'--------------------------------------------------------------------------------------------!'
    write(*,*)'Starting PCG iteration ',i
   end if

!  Apply the Preconditioner

   if (iproc ==0) then
    call yaml_sequence(advance='no')
   end if

    call H_potential('G',pkernel,lv,pot_ion,ehartree,0.0_dp,.false.)

   beta0 = beta
   beta=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
        z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
        beta=beta+r(i1,i2,i3,isp)*z(i1,i2,i3,isp)
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential correction
      !q(i1,i2,i3,isp)=r(i1,i2,i3,isp)+z(i1,i2,i3,isp)*corr(i1,i2,i3,isp)
     end do
    end do
   end do


   k=0.d0
   isp=1

  do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        zeta=z(i1,i2,i3,isp)
        !epsc=corr(i1,i2,i3,isp)
        !epsc=pkernel%corr(i1,i2,i3)
        epsc=corr3(i1,i2,i3)
        pval=p(i1,i2,i3,isp)
        qval=q(i1,i2,i3,isp)
        rval=r(i1,i2,i3,isp)
        pval = zeta+(beta/beta0)*pval
        ! Additional contribution to the Generalized Poisson operator
        pbval=-switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(zeta) 
                                                                          ! for the Poisson-Boltzmann solution.
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dsinh(multp*zeta)
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*zeta*dcosh(multp*x(i1,i2,i3,isp))
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dtanh(multp*zeta)
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*zeta
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*(zeta**2)
        qval = zeta*epsc+rval+pbval+(beta/beta0)*qval
!        qval = zeta*epsc+rval+(beta/beta0)*qval
        k = k + pval*qval
        p(i1,i2,i3,isp) = pval
        q(i1,i2,i3,isp) = qval
        !p(i1,i2,i3,isp) = z(i1,i2,i3,isp)+(beta/beta0)*p(i1,i2,i3,isp)
        !q(i1,i2,i3,isp) = q(i1,i2,i3,isp)+(beta/beta0)*qold(i1,i2,i3,isp)
        !qold(i1,i2,i3,isp)=q(i1,i2,i3,isp)
        !k=k+p(i1,i2,i3,isp)*q(i1,i2,i3,isp)
     end do
    end do
   end do

   alpha = beta/k
   !write(*,*)alpha

   offsetnew=0.d0
   normr=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      x(i1,i2,i3,isp) = x(i1,i2,i3,isp) + alpha*p(i1,i2,i3,isp)
      offsetnew=offsetnew+alpha*p(i1,i2,i3,isp)
      r(i1,i2,i3,isp) = r(i1,i2,i3,isp) - alpha*q(i1,i2,i3,isp)
      normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
      lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
     end do
    end do
   end do
!   normr=dsqrt(normr)
   normr=sqrt(normr/rpoints)

   if (trim(geocode) == 'P') then

    if (i.eq.1) then
     shift=offsetnew/(n01*n02*n03)-offset/(acell**3)
    else if (i.gt.1) then
     shift=offsetnew/(n01*n02*n03)
    end if

    isp=1
    do i3=1,n03
     do i2=1,n02
      do i1=1,n01
       x(i1,i2,i3,isp) = x(i1,i2,i3,isp) - shift
      end do
     end do
    end do

   end if

   ratio=normr/normb
   if (iproc ==0) then
   write(18,'(1x,I8,3(1x,e14.7))')i,normr,ratio,beta
   !write(*,'(1x,I8,2(1x,e14.7))')i,ratio,beta
   call EPS_iter_output_LG(i,normb,normr,ratio,alpha,beta)
!   call writeroutine(n01,n02,n03,nspden,r,i)
   call writeroutinePot(n01,n02,n03,nspden,x,i,scal*potential)
   end if

!   if (i_scf==1.and.i==3) exit   
   if (normr.lt.error) exit
   if (normr.gt.max_ratioex) exit

  end do ! end PCG loop

  if (iproc==0) then
   write(*,'(a)')'Termination of Preconditioned Conjugate Gradient'
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
  end if

  end do ! end restart loop

   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      b(i1,i2,i3,isp) = x(i1,i2,i3,isp)
     end do
    end do
   end do

  call yaml_sequence_close()
   !write(*,*)
   !write(*,'(1x,a,1x,I8)')'PCG iterations =',i-1
   !write(*,'(1x,a,1x,e14.7)')'PCG error =',ratio
   !write(*,*)
   !write(*,*)'Max abs difference between analytic potential and the computed one'
!  if (iproc==0) then
!   call writeroutinePot(n01,n02,n03,nspden,b,i-1,potential)
!   write(*,*)
!  end if

  close(unit=18)
  close(unit=38)

  if (iproc==0) then
   write(*,'(a)')'Termination of restart loop'
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
  end if
  call f_free(x)
  call f_free(r)
  call f_free(z)
  call f_free(p)
  call f_free(q)
  call f_free(qold)
  call f_free(lv)
  call f_free(corr)
  call f_free(deps)
  call f_free(ddeps)
  call f_free(de2)
  call f_free(b_old)

end subroutine  Prec_conjugate_gradient_restart

subroutine Poisson_Boltzmann(n01,n02,n03,nspden,iproc,hx,hy,hz,b,acell,eps,SetEps,nord,pkernel,potential,corr3,oneosqrteps,multp)

  use Poisson_Solver
  use yaml_output
  use f_utils
  use dynamic_memory
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden,iproc
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp
  type(coulomb_operator), intent(inout) :: pkernel
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  integer, intent(in) :: SetEps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,corr3,oneosqrteps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: b

  real(kind=8), dimension(:,:,:,:), allocatable :: x,r,z,p,q,qold,lv,corr,deps,r_PB
  !real(kind=8), dimension(n01,n02,n03,3) :: deps
  real(kind=8), dimension(:,:,:), allocatable :: de2,ddeps
  integer, parameter :: max_iter = 50
  integer, parameter :: max_iter_PB = 150
  real(kind=8), parameter :: max_ratioex = 1.0d10
  real(kind=8), parameter :: max_ratioex_PB = 1.0d10
  real(kind=8) :: alpha,beta,beta0,betanew,normb,normr,ratio,k,epsc,zeta,pval,qval,rval,pbval,multvar
  integer :: i,ii,j,i1,i2,i3,isp,i_PB
  real(kind=8), parameter :: error = 1.0d-13
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8), parameter :: eta = 1.0d0 ! Mixing parameter for the Poisson-Boltzmann ionic charge.
  real(kind=8), parameter :: tauPB = 1.0d-13 ! Polarization Iterative Method parameter.
  real(kind=8), dimension(n01,n02,n03) ::pot_ion
  real(kind=8) :: ehartree,offset,pi,switch,rpoints,res,rho,rhores2,normrPB
  real(kind=8) :: PB_charge

  !allocate heap arrays
  x=f_malloc([n01,n02,n03,nspden],id='x')
  r=f_malloc([n01,n02,n03,nspden],id='r')
  z=f_malloc([n01,n02,n03,nspden],id='z')
  p=f_malloc([n01,n02,n03,nspden],id='p')
  q=f_malloc([n01,n02,n03,nspden],id='q')
  qold=f_malloc([n01,n02,n03,nspden],id='qold')
  lv=f_malloc([n01,n02,n03,nspden],id='lv')
  corr=f_malloc([n01,n02,n03,nspden],id='corr')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  ddeps=f_malloc([n01,n02,n03],id='ddeps')
  de2=f_malloc([n01,n02,n03],id='de2')
  r_PB=f_malloc([n01,n02,n03,nspden],id='r_PB')

  pi = 4.d0*datan(1.d0)   
  rpoints=product(real([n01,n02,n03],kind=8))

  open(unit=18,file='PB_PCG_normr.dat',status='unknown')
  open(unit=38,file='PB_PCG_accuracy.dat',status='unknown')

  if (iproc ==0) then
   write(18,'(1x,a)')'iter_PB normrPB rhores2'
   write(38,'(1x,a)')'iter i1_max i2_max i3_max max_val max_center'
   call yaml_map('rpoints',rpoints)
   call yaml_sequence_open('Embedded PSolver, Preconditioned Conjugate Gradient Method')
  end if

  switch=0.0d0
  if (SetEps.eq.6) then
   switch=1.0d0
  end if

!--------------------------------------------------------------------------------------------
! Set the correction vector for the Generalized Laplace operator

!  call fssnordEpsilonDerivative(n01,n02,n03,nspden,hx,hy,hz,eps,de2,ddeps,nord,acell)

!  call fssnord3DmatNabla3varde2(n01,n02,n03,nspden,hx,hy,hx,eps,deps,de2,nord,acell)
!  call fssnord3DmatDiv3var(n01,n02,n03,nspden,hx,hy,hz,deps,ddeps,nord,acell)

!  isp=1
!  do i3=1,n03
!   do i2=1,n02
!    do i1=1,n01
!     corr(i1,i2,i3,isp)=(-0.125d0/pi)*(0.5d0*de2(i1,i2,i3)/eps(i1,i2,i3)-ddeps(i1,i2,i3))
!    end do
!   end do
!  end do
!--------------------------------------------------------------------------------------------

  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting a Poisson-Bolzmann run'
  end if
  
  call f_zero(x)
  call f_zero(r_PB)
  call f_memcpy(src=b,dest=r)

  beta=1.d0
  ratio=1.d0

  do i_PB=1,max_iter_PB ! Poisson-Boltzmann loop.

   if (iproc==0) then
    write(*,'(a)')'--------------------------------------------------------------------------------------------!'
    write(*,*)'Starting Poisson-Boltzmann iteration ',i_PB
   end if


  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting Preconditioned Conjugate Gradient'
  end if

  normb=0.d0
  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     normb=normb+b(i1,i2,i3,isp)*b(i1,i2,i3,isp)
     !!lv(i1,i2,i3,isp) = b(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
    end do
   end do
  end do
!!  normb=dsqrt(normb)
  normb=sqrt(normb/rpoints)

!  call f_memcpy(src=b,dest=r)
  call f_zero(x)
  call f_zero(q)
  call f_zero(p)
  beta=1.d0
  ratio=1.d0
  normr=normb

  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
           !lv(i1,i2,i3,isp) = pkernel%oneoeps(i1,i2,i3)*r(i1,i2,i3,isp)
           lv(i1,i2,i3,isp) = oneosqrteps(i1,i2,i3)*r(i1,i2,i3,isp)
        end do
     end do
  end do

  do i=1,max_iter

   if (normr.lt.error) exit
   if (normr.gt.max_ratioex) exit

   if (iproc==0) then
    write(*,'(a)')'--------------------------------------------------------------------------------------------!'
    write(*,*)'Starting PCG iteration ',i
   end if

!  Apply the Preconditioner

   if (iproc ==0) then
    call yaml_sequence(advance='no')
   end if
   call H_potential('G',pkernel,lv,pot_ion,ehartree,offset,.false.)

   beta0 = beta
   beta=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
        z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
        beta=beta+r(i1,i2,i3,isp)*z(i1,i2,i3,isp)
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential correction
      !q(i1,i2,i3,isp)=r(i1,i2,i3,isp)+z(i1,i2,i3,isp)*corr(i1,i2,i3,isp)
     end do
    end do
   end do


   k=0.d0
   isp=1

  do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        zeta=z(i1,i2,i3,isp)
        !epsc=corr(i1,i2,i3,isp)
        !epsc=pkernel%corr(i1,i2,i3)
        epsc=corr3(i1,i2,i3)
        pval=p(i1,i2,i3,isp)
        qval=q(i1,i2,i3,isp)
        rval=r(i1,i2,i3,isp)
        pval = zeta+(beta/beta0)*pval
        pbval=0.d0
!        pbval=-switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(zeta) ! Additional contribution to the Generalized Poisson operator
!                                                                           ! for the Poisson-Boltzmann solution.
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dsinh(multp*zeta)
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*zeta*dcosh(multp*x(i1,i2,i3,isp))
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dtanh(multp*zeta)
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*zeta
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*(zeta**2)
        qval = zeta*epsc+rval+pbval+(beta/beta0)*qval
        k = k + pval*qval
        p(i1,i2,i3,isp) = pval
        q(i1,i2,i3,isp) = qval
        !p(i1,i2,i3,isp) = z(i1,i2,i3,isp)+(beta/beta0)*p(i1,i2,i3,isp)
        !q(i1,i2,i3,isp) = q(i1,i2,i3,isp)+(beta/beta0)*qold(i1,i2,i3,isp)
        !qold(i1,i2,i3,isp)=q(i1,i2,i3,isp)
        !k=k+p(i1,i2,i3,isp)*q(i1,i2,i3,isp)
     end do
    end do
   end do

   alpha = beta/k
   !write(*,*)alpha

   normr=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      x(i1,i2,i3,isp) = x(i1,i2,i3,isp) + alpha*p(i1,i2,i3,isp)
      r(i1,i2,i3,isp) = r(i1,i2,i3,isp) - alpha*q(i1,i2,i3,isp)
      normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
      lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
     end do
    end do
   end do
!   normr=dsqrt(normr)
   normr=sqrt(normr/rpoints)

   ratio=normr/normb
   if (iproc ==0) then
!   write(18,'(1x,I8,3(1x,e14.7))')i,normr,ratio,beta
   !write(*,'(1x,I8,2(1x,e14.7))')i,ratio,beta
   call EPS_iter_output_LG(i,normb,normr,ratio,alpha,beta)
!   call writeroutine(n01,n02,n03,nspden,r,i)
   call writeroutinePot(n01,n02,n03,nspden,x,i,potential)
   end if

  end do ! PCG loop

   rhores2=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      zeta=x(i1,i2,i3,isp)
      res=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(zeta) ! Additional contribution to the Generalized Poisson operator
                                                                      ! for the Poisson-Boltzmann equation.
      rho=r_PB(i1,i2,i3,isp)
      res=res-rho
      res=eta*res
      rhores2=rhores2+res*res
      r_PB(i1,i2,i3,isp)=res+rho
      r(i1,i2,i3,isp) = b(i1,i2,i3,isp) + r_PB(i1,i2,i3,isp)
     end do
    end do
   end do

  normrPB=sqrt(rhores2/rpoints)

   if (iproc==0) then
    write(*,'(a)')'--------------------------------------------------------------------------------------------!'
    write(*,*)'End Poisson-Boltzmann iteration ',i_PB
   end if

  if (iproc ==0) then
   call yaml_map('iter PB',i_PB)
   call yaml_map('normrPB',normrPB)
   call yaml_map('rhores2',rhores2)
   write(18,'(1x,I8,3(1x,e14.7))')i_PB,normrPB,rhores2
   call writeroutinePot(n01,n02,n03,nspden,x,i_PB,potential)
  end if

   if (normrPB.lt.tauPB) exit
   if (normrPB.gt.max_ratioex_PB) exit

 end do ! Poisson-Boltzmann loop

   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      b(i1,i2,i3,isp) = x(i1,i2,i3,isp)
     end do
    end do
   end do

  call yaml_sequence_close()
   !write(*,*)
   !write(*,'(1x,a,1x,I8)')'PCG iterations =',i-1
   !write(*,'(1x,a,1x,e14.7)')'PCG error =',ratio
   !write(*,*)
   !write(*,*)'Max abs difference between analytic potential and the computed one'
!  if (iproc==0) then
!   call writeroutinePot(n01,n02,n03,nspden,b,i-1,potential)
!   write(*,*)
!  end if

  close(unit=18)
  close(unit=38)

  if (iproc==0) then
   write(*,'(a)')'Termination of Preconditioned Conjugate Gradient'
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
  end if

  call f_free(x)
  call f_free(r)
  call f_free(z)
  call f_free(p)
  call f_free(q)
  call f_free(qold)
  call f_free(lv)
  call f_free(corr)
  call f_free(deps)
  call f_free(ddeps)
  call f_free(de2)
  call f_free(r_PB)

end subroutine Poisson_Boltzmann

subroutine Poisson_Boltzmann_good(n01,n02,n03,nspden,iproc,hx,hy,hz,b,acell,eps,SetEps,nord,pkernel,potential,&
  density,oneoeps,dlogeps,corr3,oneosqrteps,multp,geocode)

  use Poisson_Solver
  use yaml_output
  use f_utils
  use dynamic_memory
  use PSbox, only: PS_gather
  use f_enums
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden,iproc
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp
  type(coulomb_operator), intent(inout) :: pkernel
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  integer, intent(in) :: SetEps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,oneoeps,dlogeps,corr3,oneosqrteps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: b
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: density
  character(len=2), intent(in) :: geocode

  real(kind=8), dimension(:,:,:,:), allocatable :: x,r,z,p,q,qold,lv,corr,deps,r_PB
  !real(kind=8), dimension(n01,n02,n03,3) :: deps
  real(kind=8), dimension(:,:,:), allocatable :: de2,ddeps
  integer, parameter :: max_iter = 50
  integer, parameter :: max_iter_PB = 150
  real(kind=8), parameter :: max_ratioex = 1.0d10
  real(kind=8), parameter :: max_ratioex_PB = 1.0d10
  real(kind=8) :: alpha,beta,beta0,betanew,normb,normr,ratio,k,epsc,zeta,pval,qval,rval,pbval,multvar
  integer :: i,ii,j,i1,i2,i3,isp,i_PB
  real(kind=8), parameter :: error = 1.0d-13
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8), parameter :: eta = 1.0d0 ! Mixing parameter for the Poisson-Boltzmann ionic charge.
  real(kind=8), parameter :: tauPB = 1.0d-13 ! Polarization Iterative Method parameter.
  real(kind=8), dimension(n01,n02,n03) ::pot_ion
  real(kind=8) :: ehartree,offset,pi,switch,rpoints,res,rho,rhores2,normrPB
  real(kind=8) :: PB_charge

  !allocate heap arrays
  x=f_malloc([n01,n02,n03,nspden],id='x')
  r=f_malloc([n01,n02,n03,nspden],id='r')
  z=f_malloc([n01,n02,n03,nspden],id='z')
  p=f_malloc([n01,n02,n03,nspden],id='p')
  q=f_malloc([n01,n02,n03,nspden],id='q')
  qold=f_malloc([n01,n02,n03,nspden],id='qold')
  lv=f_malloc([n01,n02,n03,nspden],id='lv')
  corr=f_malloc([n01,n02,n03,nspden],id='corr')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  ddeps=f_malloc([n01,n02,n03],id='ddeps')
  de2=f_malloc([n01,n02,n03],id='de2')
  r_PB=f_malloc([n01,n02,n03,nspden],id='r_PB')

  pi = 4.d0*datan(1.d0)
  rpoints=product(real([n01,n02,n03],kind=8))

  open(unit=18,file='PB_PCG_normr.dat',status='unknown')
  open(unit=38,file='PB_PCG_accuracy.dat',status='unknown')

  if (iproc ==0) then
     write(18,'(1x,a)')'iter_PB normrPB rhores2'
     write(38,'(1x,a)')'iter i1_max i2_max i3_max max_val max_center'
  end if

  !switch=0.0d0
  !if (SetEps.eq.6) then
  switch=1.0d0
  !end if

  if (.not. (pkernel%method .hasattr. 'PB_NONE')) then
     call H_potential('D',pkernel,b(1,1,pkernel%grid%istart+1,1),b(1,1,pkernel%grid%istart+1,1),&
          ehartree,offset,.false.)
  else
     call yaml_map('rpoints',rpoints)
     call yaml_sequence_open('Embedded PSolver, Preconditioned Conjugate Gradient Method')
     if (iproc==0) then
        write(*,'(a)')'--------------------------------------------------------------------------------------------'
        write(*,'(a)')'Starting a Poisson-Bolzmann run'
     end if

     call f_zero(x)
     call f_zero(r_PB)
     call f_memcpy(src=b,dest=r)

     beta=1.d0
     ratio=1.d0

     do i_PB=1,max_iter_PB ! Poisson-Boltzmann loop.

        if (iproc==0) then
           write(*,'(a)')'--------------------------------------------------------------------------------------------!'
           write(*,*)'Starting Poisson-Boltzmann iteration ',i_PB
        end if


        if (iproc==0) then
           write(*,'(a)')'--------------------------------------------------------------------------------------------'
           write(*,'(a)')'Starting Preconditioned Conjugate Gradient'
        end if

        !-------------------------------------------------------------------------------------
!!$  if (i_PB.eq.1) then
!!$   call Prec_conjugate_gradient(n01,n02,n03,nspden,iproc,hx,hy,hz,r,density,acell,&
!!$        eps,SetEps,nord,pkernel,potential,corr3,oneosqrteps,dlogeps,multp,offset,geocode,.false.,.false.,.false.)
!!$  else
!!$   call Prec_conjugate_gradient_Inputguess(n01,n02,n03,nspden,iproc,hx,hy,hz,x,r,acell,&
!!$        eps,SetEps,nord,pkernel,potential,corr3,oneosqrteps,dlogeps,multp,offset,geocode,.false.)
!!$  end if
!!$  !call PolarizationIteration_Inputguess(n01,n02,n03,nspden,iproc,hx,hy,hz,rhopot,density,acell,&
!!$  !     eps,nord,pkernel,potential,oneoeps,dlogeps,multp,offset,geocode)
        !pkernel%minres=1.0d-13
        !pkernel%opt%use_input_guess=.true.
        call H_potential('D',pkernel,r(1,1,pkernel%grid%istart+1,1),r(1,1,pkernel%grid%istart+1,1),&
             ehartree,offset,.false.)
        call PS_gather(r,pkernel)
        call f_memcpy(src=r,dest=x)

        !-------------------------------------------------------------------------------------

        if (i_PB.eq.1) call f_memcpy(src=r,dest=x)
        rhores2=0.d0
        isp=1
        do i3=1,n03
           do i2=1,n02
              do i1=1,n01
                 zeta=x(i1,i2,i3,isp)
                 !x(i1,i2,i3,isp)=zeta
                 res=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(zeta) ! Additional contribution to the Generalized Poisson operator
                 ! for the Poisson-Boltzmann equation.
                 rho=r_PB(i1,i2,i3,isp)
                 res=res-rho
                 res=eta*res
                 rhores2=rhores2+res*res
                 r_PB(i1,i2,i3,isp)=res+rho
                 r(i1,i2,i3,isp) = b(i1,i2,i3,isp) + r_PB(i1,i2,i3,isp)
              end do
           end do
        end do

        normrPB=sqrt(rhores2/rpoints)

        if (iproc==0) then
           write(*,'(a)')'--------------------------------------------------------------------------------------------!'
           write(*,*)'End Poisson-Boltzmann iteration ',i_PB
        end if

        if (iproc ==0) then
           call yaml_map('iter PB',i_PB)
           call yaml_map('normrPB',normrPB)
           call yaml_map('rhores2',rhores2)
           write(18,'(1x,I8,3(1x,e14.7))')i_PB,normrPB,rhores2
           call writeroutinePot(n01,n02,n03,nspden,x,i_PB,potential)
        end if

        if (normrPB.lt.tauPB) exit
        if (normrPB.gt.max_ratioex_PB) exit

     end do ! Poisson-Boltzmann loop

     isp=1
     do i3=1,n03
        do i2=1,n02
           do i1=1,n01
              b(i1,i2,i3,isp) = x(i1,i2,i3,isp)
           end do
        end do
     end do

     call yaml_sequence_close()
     !write(*,*)
     !write(*,'(1x,a,1x,I8)')'PCG iterations =',i-1
     !write(*,'(1x,a,1x,e14.7)')'PCG error =',ratio
     !write(*,*)
     !write(*,*)'Max abs difference between analytic potential and the computed one'
     !  if (iproc==0) then
     !   call writeroutinePot(n01,n02,n03,nspden,b,i-1,potential)
     !   write(*,*)
     !  end if

  end if

  close(unit=18)
  close(unit=38)

  if (iproc==0) then
     write(*,'(a)')'Termination of Preconditioned Conjugate Gradient'
     write(*,'(a)')'--------------------------------------------------------------------------------------------'
  end if

  call f_free(x)
  call f_free(r)
  call f_free(z)
  call f_free(p)
  call f_free(q)
  call f_free(qold)
  call f_free(lv)
  call f_free(corr)
  call f_free(deps)
  call f_free(ddeps)
  call f_free(de2)
  call f_free(r_PB)

end subroutine Poisson_Boltzmann_good

subroutine Poisson_Boltzmann_improved(n01,n02,n03,nspden,iproc,hx,hy,hz,b,&
     acell,eps,SetEps,nord,pkernel,potential,corr3,oneosqrteps,multp)

  use Poisson_Solver
  use yaml_output
  use f_utils
  use dynamic_memory
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden,iproc
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp
  type(coulomb_operator), intent(inout) :: pkernel
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  integer, intent(in) :: SetEps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,corr3,oneosqrteps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: b

  real(kind=8), dimension(:,:,:,:), allocatable :: x,r,z,p,q,qold,lv,corr,deps,r_PB,x_PB,x_check
  !real(kind=8), dimension(n01,n02,n03,3) :: deps
  real(kind=8), dimension(:,:,:), allocatable :: de2,ddeps
  integer, parameter :: max_iter = 50
  integer, parameter :: max_iter_PB = 150
  real(kind=8), parameter :: max_ratioex = 1.0d10
  real(kind=8), parameter :: max_ratioex_PB = 1.0d10
  real(kind=8) :: alpha,beta,beta0,betanew,normb,normr,ratio,k,epsc,zeta,pval,qval,rval,pbval,multvar
  integer :: i,ii,j,i1,i2,i3,isp,i_PB
  real(kind=8), parameter :: error = 1.0d-13
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8), parameter :: eta = 1.0d0 ! Mixing parameter for the Poisson-Boltzmann ionic charge.
  real(kind=8), parameter :: tauPB = 1.0d-13 ! Polarization Iterative Method parameter.
  real(kind=8), dimension(n01,n02,n03) ::pot_ion
  real(kind=8) :: ehartree,offset,pi,switch,rpoints,res,rho,rhores2,normrPB
  real(kind=8) :: PB_charge

  !allocate heap arrays
  x=f_malloc([n01,n02,n03,nspden],id='x')
  r=f_malloc([n01,n02,n03,nspden],id='r')
  z=f_malloc([n01,n02,n03,nspden],id='z')
  p=f_malloc([n01,n02,n03,nspden],id='p')
  q=f_malloc([n01,n02,n03,nspden],id='q')
  qold=f_malloc([n01,n02,n03,nspden],id='qold')
  lv=f_malloc([n01,n02,n03,nspden],id='lv')
  corr=f_malloc([n01,n02,n03,nspden],id='corr')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  ddeps=f_malloc([n01,n02,n03],id='ddeps')
  de2=f_malloc([n01,n02,n03],id='de2')
  r_PB=f_malloc([n01,n02,n03,nspden],id='r_PB')
  x_PB=f_malloc([n01,n02,n03,nspden],id='x_PB')
  x_check=f_malloc([n01,n02,n03,nspden],id='x_check')

  pi = 4.d0*datan(1.d0)   
  rpoints=product(real([n01,n02,n03],kind=8))

  open(unit=18,file='PBimpro_PCG_normr.dat',status='unknown')
  open(unit=38,file='PCGimpro_accuracy.dat',status='unknown')

  if (iproc ==0) then
   write(18,'(1x,a)')'iter_PB normrPB rhores2'
   write(38,'(1x,a)')'iter i1_max i2_max i3_max max_val max_center'
   call yaml_map('rpoints',rpoints)
   call yaml_sequence_open('Embedded PSolver, Preconditioned Conjugate Gradient Method')
  end if

  switch=0.0d0
  if (SetEps.eq.6) then
   switch=1.0d0
  end if

!--------------------------------------------------------------------------------------------
! Set the correction vector for the Generalized Laplace operator

!  call fssnordEpsilonDerivative(n01,n02,n03,nspden,hx,hy,hz,eps,de2,ddeps,nord,acell)

!  call fssnord3DmatNabla3varde2(n01,n02,n03,nspden,hx,hy,hx,eps,deps,de2,nord,acell)
!  call fssnord3DmatDiv3var(n01,n02,n03,nspden,hx,hy,hz,deps,ddeps,nord,acell)

!  isp=1
!  do i3=1,n03
!   do i2=1,n02
!    do i1=1,n01
!     corr(i1,i2,i3,isp)=(-0.125d0/pi)*(0.5d0*de2(i1,i2,i3)/eps(i1,i2,i3)-ddeps(i1,i2,i3))
!    end do
!   end do
!  end do
!--------------------------------------------------------------------------------------------

  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting a Poisson-Bolzmann run'
  end if
  
  call f_zero(x)
  call f_zero(r_PB)
  call f_zero(x_PB)
  call f_memcpy(src=b,dest=r)

  beta=1.d0
  ratio=1.d0

  do i_PB=1,max_iter_PB ! Poisson-Boltzmann loop.

   if (iproc==0) then
    write(*,'(a)')'--------------------------------------------------------------------------------------------!'
    write(*,*)'Starting Poisson-Boltzmann iteration ',i_PB
   end if


  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting Preconditioned Conjugate Gradient'
  end if

  normb=0.d0
  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     normb=normb+b(i1,i2,i3,isp)*b(i1,i2,i3,isp)
     !!lv(i1,i2,i3,isp) = b(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
    end do
   end do
  end do
!!  normb=dsqrt(normb)
  normb=sqrt(normb/rpoints)

!  call f_memcpy(src=b,dest=r)
  call f_zero(x)
  call f_zero(q)
  call f_zero(p)
  beta=1.d0
  ratio=1.d0
  normr=normb

  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
           !lv(i1,i2,i3,isp) = pkernel%oneoeps(i1,i2,i3)*r(i1,i2,i3,isp)
           lv(i1,i2,i3,isp) = oneosqrteps(i1,i2,i3)*r(i1,i2,i3,isp)
        end do
     end do
  end do

  do i=1,max_iter

   if (normr.lt.error) exit
   if (normr.gt.max_ratioex) exit

   if (iproc==0) then
    write(*,'(a)')'--------------------------------------------------------------------------------------------!'
    write(*,*)'Starting PCG iteration ',i
   end if

!  Apply the Preconditioner

   if (iproc ==0) then
    call yaml_sequence(advance='no')
   end if
   call H_potential('G',pkernel,lv,pot_ion,ehartree,offset,.false.)

   beta0 = beta
   beta=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
        z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
        beta=beta+r(i1,i2,i3,isp)*z(i1,i2,i3,isp)
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential correction
      !q(i1,i2,i3,isp)=r(i1,i2,i3,isp)+z(i1,i2,i3,isp)*corr(i1,i2,i3,isp)
     end do
    end do
   end do


   k=0.d0
   isp=1

  do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        zeta=z(i1,i2,i3,isp)
        !epsc=corr(i1,i2,i3,isp)
        !epsc=pkernel%corr(i1,i2,i3)
        epsc=corr3(i1,i2,i3)
        pval=p(i1,i2,i3,isp)
        qval=q(i1,i2,i3,isp)
        rval=r(i1,i2,i3,isp)
        pval = zeta+(beta/beta0)*pval
        pbval=0.d0
 ! Additional contribution to the Generalized Poisson operator
!        pbval=-switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(zeta)
!                                                                           ! for the Poisson-Boltzmann solution.
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dsinh(multp*zeta)
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*zeta*dcosh(multp*x(i1,i2,i3,isp))
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dtanh(multp*zeta)
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*zeta
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*(zeta**2)
        qval = zeta*epsc+rval+pbval+(beta/beta0)*qval
        k = k + pval*qval
        p(i1,i2,i3,isp) = pval
        q(i1,i2,i3,isp) = qval
        !p(i1,i2,i3,isp) = z(i1,i2,i3,isp)+(beta/beta0)*p(i1,i2,i3,isp)
        !q(i1,i2,i3,isp) = q(i1,i2,i3,isp)+(beta/beta0)*qold(i1,i2,i3,isp)
        !qold(i1,i2,i3,isp)=q(i1,i2,i3,isp)
        !k=k+p(i1,i2,i3,isp)*q(i1,i2,i3,isp)
     end do
    end do
   end do

   alpha = beta/k
   !write(*,*)alpha

   normr=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      x(i1,i2,i3,isp) = x(i1,i2,i3,isp) + alpha*p(i1,i2,i3,isp)
      r(i1,i2,i3,isp) = r(i1,i2,i3,isp) - alpha*q(i1,i2,i3,isp)
      normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
      lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
      x_check(i1,i2,i3,isp) = x_PB(i1,i2,i3,isp) + x(i1,i2,i3,isp)
     end do
    end do
   end do
!   normr=dsqrt(normr)
   normr=sqrt(normr/rpoints)

   ratio=normr/normb
   if (iproc ==0) then
!   write(18,'(1x,I8,3(1x,e14.7))')i,normr,ratio,beta
   !write(*,'(1x,I8,2(1x,e14.7))')i,ratio,beta
   call EPS_iter_output_LG(i,normb,normr,ratio,alpha,beta)
!   call writeroutine(n01,n02,n03,nspden,r,i)
   call writeroutinePot(n01,n02,n03,nspden,x_check,i,potential)
   end if

  end do ! PCG loop

   rhores2=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      zeta=x(i1,i2,i3,isp)
      res=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(zeta) ! Additional contribution to the Generalized Poisson operator
                                                                      ! for the Poisson-Boltzmann equation.
      rho=r_PB(i1,i2,i3,isp)
      res=res-rho
      res=eta*res
      rhores2=rhores2+res*res
      r_PB(i1,i2,i3,isp)=res+rho
      x_PB(i1,i2,i3,isp) = x_PB(i1,i2,i3,isp) + x(i1,i2,i3,isp)
!      r(i1,i2,i3,isp) = b(i1,i2,i3,isp) + r_PB(i1,i2,i3,isp)
      r(i1,i2,i3,isp) = r(i1,i2,i3,isp) + r_PB(i1,i2,i3,isp)
     end do
    end do
   end do

  normrPB=sqrt(rhores2/rpoints)

   if (iproc==0) then
    write(*,'(a)')'--------------------------------------------------------------------------------------------!'
    write(*,*)'End Poisson-Boltzmann iteration ',i_PB
   end if

  if (iproc ==0) then
   call yaml_map('iter PB',i_PB)
   call yaml_map('normrPB',normrPB)
   call yaml_map('rhores2',rhores2)
   write(18,'(1x,I8,3(1x,e14.7))')i_PB,normrPB,rhores2
   call writeroutinePot(n01,n02,n03,nspden,x_PB,i_PB,potential)
  end if

   if (normrPB.lt.tauPB) exit
   if (normrPB.gt.max_ratioex_PB) exit

 end do ! Poisson-Boltzmann loop

   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      b(i1,i2,i3,isp) = x_PB(i1,i2,i3,isp)
     end do
    end do
   end do

  call yaml_sequence_close()
   !write(*,*)
   !write(*,'(1x,a,1x,I8)')'PCG iterations =',i-1
   !write(*,'(1x,a,1x,e14.7)')'PCG error =',ratio
   !write(*,*)
   !write(*,*)'Max abs difference between analytic potential and the computed one'
!  if (iproc==0) then
!   call writeroutinePot(n01,n02,n03,nspden,b,i-1,potential)
!   write(*,*)
!  end if

  close(unit=18)
  close(unit=38)

  if (iproc==0) then
   write(*,'(a)')'Termination of Preconditioned Conjugate Gradient'
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
  end if

  call f_free(x)
  call f_free(r)
  call f_free(z)
  call f_free(p)
  call f_free(q)
  call f_free(qold)
  call f_free(lv)
  call f_free(corr)
  call f_free(deps)
  call f_free(ddeps)
  call f_free(de2)
  call f_free(r_PB)
  call f_free(x_PB)
  call f_free(x_check)

end subroutine Poisson_Boltzmann_improved

subroutine Poisson_Boltzmann_improved2(n01,n02,n03,nspden,iproc,hx,hy,hz,b,&
     acell,eps,SetEps,nord,pkernel,potential,corr3,oneosqrteps,multp)

  use Poisson_Solver
  use yaml_output
  use f_utils
  use dynamic_memory
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden,iproc
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp
  type(coulomb_operator), intent(inout) :: pkernel
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  integer, intent(in) :: SetEps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,corr3,oneosqrteps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: b

  real(kind=8), dimension(:,:,:,:), allocatable :: x,r,z,p,q,qold,lv,corr,deps,r_PB,x_PB,x_check
  !real(kind=8), dimension(n01,n02,n03,3) :: deps
  real(kind=8), dimension(:,:,:), allocatable :: de2,ddeps
  integer, parameter :: max_iter = 50
  integer, parameter :: max_iter_PB = 150
  real(kind=8), parameter :: max_ratioex = 1.0d10
  real(kind=8), parameter :: max_ratioex_PB = 1.0d10
  real(kind=8) :: alpha,beta,beta0,betanew,normb,normr,ratio,k,epsc,zeta,pval,qval,rval,pbval,multvar
  integer :: i,ii,j,i1,i2,i3,isp,i_PB
  real(kind=8), parameter :: error = 1.0d-13 !1.0d-13
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8), parameter :: eta = 1.0d0 ! Mixing parameter for the Poisson-Boltzmann ionic charge.
  real(kind=8), parameter :: tauPB = 1.0d-13 ! Exit of Poisson-Boltzmann loop, to be = error for GPE.
  real(kind=8), dimension(n01,n02,n03) ::pot_ion
  real(kind=8) :: ehartree,offset,pi,switch,rpoints,res,rho,rhores2,normrPB,errorvar
  real(kind=8) :: PB_charge

  !allocate heap arrays
  x=f_malloc([n01,n02,n03,nspden],id='x')
  r=f_malloc([n01,n02,n03,nspden],id='r')
  z=f_malloc([n01,n02,n03,nspden],id='z')
  p=f_malloc([n01,n02,n03,nspden],id='p')
  q=f_malloc([n01,n02,n03,nspden],id='q')
  qold=f_malloc([n01,n02,n03,nspden],id='qold')
  lv=f_malloc([n01,n02,n03,nspden],id='lv')
  corr=f_malloc([n01,n02,n03,nspden],id='corr')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  ddeps=f_malloc([n01,n02,n03],id='ddeps')
  de2=f_malloc([n01,n02,n03],id='de2')
  r_PB=f_malloc([n01,n02,n03,nspden],id='r_PB')
  x_PB=f_malloc([n01,n02,n03,nspden],id='x_PB')
  x_check=f_malloc([n01,n02,n03,nspden],id='x_check')

  pi = 4.d0*datan(1.d0)   
  rpoints=product(real([n01,n02,n03],kind=8))

  open(unit=18,file='PBimpro_PCG_normr.dat',status='unknown')
  open(unit=38,file='PCGimpro_accuracy.dat',status='unknown')

  if (iproc ==0) then
   write(18,'(1x,a)')'iter_PB normrPB rhores2'
   write(38,'(1x,a)')'iter i1_max i2_max i3_max max_val max_center'
   call yaml_map('rpoints',rpoints)
   call yaml_sequence_open('Embedded PSolver, Preconditioned Conjugate Gradient Method')
  end if

  switch=0.0d0
  if (SetEps.eq.6) then
   switch=1.0d0
  end if

!--------------------------------------------------------------------------------------------
! Set the correction vector for the Generalized Laplace operator

!  call fssnordEpsilonDerivative(n01,n02,n03,nspden,hx,hy,hz,eps,de2,ddeps,nord,acell)

!  call fssnord3DmatNabla3varde2(n01,n02,n03,nspden,hx,hy,hx,eps,deps,de2,nord,acell)
!  call fssnord3DmatDiv3var(n01,n02,n03,nspden,hx,hy,hz,deps,ddeps,nord,acell)

!  isp=1
!  do i3=1,n03
!   do i2=1,n02
!    do i1=1,n01
!     corr(i1,i2,i3,isp)=(-0.125d0/pi)*(0.5d0*de2(i1,i2,i3)/eps(i1,i2,i3)-ddeps(i1,i2,i3))
!    end do
!   end do
!  end do
!--------------------------------------------------------------------------------------------

  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting a Poisson-Bolzmann run'
  end if
  
  call f_zero(x)
  call f_zero(r_PB)
  call f_zero(x_PB)
  call f_memcpy(src=b,dest=r)

  beta=1.d0
  ratio=1.d0

  do i_PB=1,max_iter_PB ! Poisson-Boltzmann loop.

   if (iproc==0) then
    write(*,'(a)')'--------------------------------------------------------------------------------------------!'
    write(*,*)'Starting Poisson-Boltzmann iteration ',i_PB
   end if


  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting Preconditioned Conjugate Gradient'
  end if

  normb=0.d0
  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     normb=normb+b(i1,i2,i3,isp)*b(i1,i2,i3,isp)
     !!lv(i1,i2,i3,isp) = b(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
    end do
   end do
  end do
!!  normb=dsqrt(normb)
  normb=sqrt(normb/rpoints)

!  call f_memcpy(src=b,dest=r)
!  call f_zero(x)
  call f_zero(q)
  call f_zero(p)
  beta=1.d0
  ratio=1.d0
  normr=normb

  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
           !lv(i1,i2,i3,isp) = pkernel%oneoeps(i1,i2,i3)*r(i1,i2,i3,isp)
           lv(i1,i2,i3,isp) = oneosqrteps(i1,i2,i3)*r(i1,i2,i3,isp)
        end do
     end do
  end do

  if (i_PB.eq.1) then
   errorvar=1.0d-6
   errorvar=error
   errorvar=max(error,errorvar)
  else if (i_PB.eq.2) then
   errorvar=1.0d-8
   errorvar=error
   errorvar=max(error,errorvar)
  else
   errorvar=error
  end if
  if (iproc ==0) then
   call yaml_map('errorvar',errorvar)
  end if

  do i=1,max_iter

   if (normr.lt.errorvar) exit
   if (normr.gt.max_ratioex) exit

   if (iproc==0) then
    write(*,'(a)')'--------------------------------------------------------------------------------------------!'
    write(*,*)'Starting PCG iteration ',i
   end if

!  Apply the Preconditioner

   if (iproc ==0) then
    call yaml_sequence(advance='no')
   end if
   call H_potential('G',pkernel,lv,pot_ion,ehartree,offset,.false.)

   beta0 = beta
   beta=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
        !z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
        z(i1,i2,i3,isp) = lv(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
        beta=beta+r(i1,i2,i3,isp)*z(i1,i2,i3,isp)
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential correction
      !q(i1,i2,i3,isp)=r(i1,i2,i3,isp)+z(i1,i2,i3,isp)*corr(i1,i2,i3,isp)
     end do
    end do
   end do


   k=0.d0
   isp=1

  do i3=1,n03
    do i2=1,n02
     do i1=1,n01
        zeta=z(i1,i2,i3,isp)
        !epsc=corr(i1,i2,i3,isp)
        !epsc=pkernel%corr(i1,i2,i3)
        epsc=corr3(i1,i2,i3)
        pval=p(i1,i2,i3,isp)
        qval=q(i1,i2,i3,isp)
        rval=r(i1,i2,i3,isp)
        pval = zeta+(beta/beta0)*pval
        pbval=0.d0
!        pbval=-switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(zeta) ! Additional contribution to the Generalized Poisson operator
!                                                                           ! for the Poisson-Boltzmann solution.
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dsinh(multp*zeta)
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*zeta*dcosh(multp*x(i1,i2,i3,isp))
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dtanh(multp*zeta)
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*zeta
!        pbval=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*(zeta**2)
        qval = zeta*epsc+rval+pbval+(beta/beta0)*qval
        k = k + pval*qval
        p(i1,i2,i3,isp) = pval
        q(i1,i2,i3,isp) = qval
        !p(i1,i2,i3,isp) = z(i1,i2,i3,isp)+(beta/beta0)*p(i1,i2,i3,isp)
        !q(i1,i2,i3,isp) = q(i1,i2,i3,isp)+(beta/beta0)*qold(i1,i2,i3,isp)
        !qold(i1,i2,i3,isp)=q(i1,i2,i3,isp)
        !k=k+p(i1,i2,i3,isp)*q(i1,i2,i3,isp)
     end do
    end do
   end do

   alpha = beta/k
   !write(*,*)alpha

   normr=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      x(i1,i2,i3,isp) = x(i1,i2,i3,isp) + alpha*p(i1,i2,i3,isp)
      r(i1,i2,i3,isp) = r(i1,i2,i3,isp) - alpha*q(i1,i2,i3,isp)
      normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/dsqrt(eps(i1,i2,i3))
      !lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*pkernel%oneoeps(i1,i2,i3)
      lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)*oneosqrteps(i1,i2,i3)
      !x_check(i1,i2,i3,isp) = x_PB(i1,i2,i3,isp) + x(i1,i2,i3,isp)
!      x_check(i1,i2,i3,isp) = x(i1,i2,i3,isp)
     end do
    end do
   end do
!   normr=dsqrt(normr)
   normr=sqrt(normr/rpoints)

   ratio=normr/normb
   if (iproc ==0) then
!   write(18,'(1x,I8,3(1x,e14.7))')i,normr,ratio,beta
   !write(*,'(1x,I8,2(1x,e14.7))')i,ratio,beta
   call EPS_iter_output_LG(i,normb,normr,ratio,alpha,beta)
!   call writeroutine(n01,n02,n03,nspden,r,i)
   call writeroutinePot(n01,n02,n03,nspden,x,i,potential)
   end if

  end do ! PCG loop

   rhores2=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      zeta=x(i1,i2,i3,isp)
      res=switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(zeta) ! Additional contribution to the Generalized Poisson operator
                                                                      ! for the Poisson-Boltzmann equation.
      rho=r_PB(i1,i2,i3,isp)
      res=res-rho
      res=eta*res
      rhores2=rhores2+res*res
      r_PB(i1,i2,i3,isp)=res+rho
!      x_PB(i1,i2,i3,isp) = x_PB(i1,i2,i3,isp) + x(i1,i2,i3,isp)
!      r(i1,i2,i3,isp) = b(i1,i2,i3,isp) + r_PB(i1,i2,i3,isp)
      r(i1,i2,i3,isp) = r(i1,i2,i3,isp) + r_PB(i1,i2,i3,isp) - rho
     end do
    end do
   end do

  normrPB=sqrt(rhores2/rpoints)

   if (iproc==0) then
    write(*,'(a)')'--------------------------------------------------------------------------------------------!'
    write(*,*)'End Poisson-Boltzmann iteration ',i_PB
   end if

  if (iproc ==0) then
   call yaml_map('iter PB',i_PB)
   call yaml_map('normrPB',normrPB)
   call yaml_map('rhores2',rhores2)
   write(18,'(1x,I8,3(1x,e14.7))')i_PB,normrPB,rhores2
   call writeroutinePot(n01,n02,n03,nspden,x,i_PB,potential)
  end if

   if (normrPB.lt.tauPB) exit
   if (normrPB.gt.max_ratioex_PB) exit

 end do ! Poisson-Boltzmann loop

   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      b(i1,i2,i3,isp) = x(i1,i2,i3,isp)
     end do
    end do
   end do

  call yaml_sequence_close()
   !write(*,*)
   !write(*,'(1x,a,1x,I8)')'PCG iterations =',i-1
   !write(*,'(1x,a,1x,e14.7)')'PCG error =',ratio
   !write(*,*)
   !write(*,*)'Max abs difference between analytic potential and the computed one'
!  if (iproc==0) then
!   call writeroutinePot(n01,n02,n03,nspden,b,i-1,potential)
!   write(*,*)
!  end if

  close(unit=18)
  close(unit=38)

  if (iproc==0) then
   write(*,'(a)')'Termination of Preconditioned Conjugate Gradient'
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
  end if

  call f_free(x)
  call f_free(r)
  call f_free(z)
  call f_free(p)
  call f_free(q)
  call f_free(qold)
  call f_free(lv)
  call f_free(corr)
  call f_free(deps)
  call f_free(ddeps)
  call f_free(de2)
  call f_free(r_PB)
  call f_free(x_PB)
  call f_free(x_check)

end subroutine Poisson_Boltzmann_improved2

subroutine Prec_Steepest_Descent(n01,n02,n03,nspden,hx,hy,hz,b,acell,eps,dlogeps,nord,pkernel,potential,geocode)

  use Poisson_Solver

  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell
  type(coulomb_operator), intent(inout) :: pkernel
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(3,n01,n02,n03), intent(in) :: dlogeps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: b
  character(len=2), intent(in) :: geocode

  real(kind=8), dimension(n01,n02,n03,nspden) :: x,r,p,lv
  real(kind=8), dimension(n01,n02,n03,nspden,3) :: dx
  real(kind=8), dimension(n01,n02,n03,3) :: deps
  real(kind=8), parameter :: eta = 0.6d0 !1.0d0 !Like eta in the Polarization Iterative Method.
  integer, parameter :: max_iter = 50
  real(kind=8), parameter :: max_ratioex = 1.0d10
  real(kind=8) :: alpha,beta,normb,normr,ratio,k,rpoints
  integer :: i,ii,j,i1,i2,i3,isp
  real(kind=8), parameter :: error = 1.0d-12
  real(kind=8), dimension(n01,n02,n03) ::pot_ion
  real(kind=8) :: ehartree,offset,pi,divprod

  pi = 4.d0*datan(1.d0)
  rpoints=product(real([n01,n02,n03],kind=8))

  open(unit=18,file='PSDConvergence.dat',status='unknown')
  open(unit=38,file='MaxAnalysisPSD.dat',status='unknown')

  write(*,'(a)')'--------------------------------------------------------------------------------------------'
  write(*,'(a)')'Starting Preconditioned Steepest Descent'
  write(*,'(a)')'Starting PSD iteration 1'

!------------------------------------------------------------------------------------
! Set the correction vector for the Generalized Laplace operator

!  call fssnord3DmatNabla3var(n01,n02,n03,nspden,hx,hy,hz,eps,deps,nord,acell)

  normb=0.d0
  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     normb=normb+b(i1,i2,i3,isp)*b(i1,i2,i3,isp)
     lv(i1,i2,i3,isp) = b(i1,i2,i3,isp)/eps(i1,i2,i3)
    end do
   end do
  end do
  normb=dsqrt(normb/rpoints)

  !call H_potential('G',pkernel,lv,pot_ion,ehartree,0.0_dp,.false.)
  call H_potential('G',pkernel,lv,pot_ion,ehartree,offset,.false.)

!------------------------------------------------------------------------------------
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential
! correction

  call fssnord3DmatNabla(geocode,n01,n02,n03,nspden,hx,hy,hz,lv,dx,nord,acell)

  normr=0.d0
  beta=0.d0
  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     x(i1,i2,i3,isp) = lv(i1,i2,i3,isp)
     beta=beta+b(i1,i2,i3,isp)*lv(i1,i2,i3,isp)
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential
! correction
     divprod = 0.d0
     do j=1,3
      !divprod = divprod + deps(i1,i2,i3,j)*dx(i1,i2,i3,isp,j)
      divprod = divprod + dlogeps(j,i1,i2,i3)*dx(i1,i2,i3,isp,j)
     end do
     r(i1,i2,i3,isp)=eta*0.25d0*divprod*eps(i1,i2,i3)/pi
     normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
     !normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)/(eps(i1,i2,i3)**2)
     lv(i1,i2,i3,isp) = eta*0.25d0*divprod/pi
    end do
   end do
  end do
  normr=dsqrt(normr/rpoints)

!------------------------------------------------------------------------------------

  ratio=normr/normb

!  call writeroutine(n01,n02,n03,nspden,r,1)
  call writeroutinePot(n01,n02,n03,nspden,potential,0,potential)
  call writeroutinePot(n01,n02,n03,nspden,x,1,potential)

  write(18,'(1x,I8,2(1x,e14.7))')1,normr,beta
  write(*,'(1x,I8,2(1x,e14.7))')1,normr,beta


  do i=2,max_iter

   if (normr.lt.error) exit
   if (ratio.gt.max_ratioex) exit

   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,*)'Starting PSD iteration ',i

!  Apply the Preconditioner

   call H_potential('G',pkernel,lv,pot_ion,ehartree,offset,.false.)

   call fssnord3DmatNabla(geocode,n01,n02,n03,nspden,hx,hy,hz,lv,dx,nord,acell)
   !call fssnord3DmatNabla(n01,n02,n03,nspden,lv,dx,nord,acell)

   beta=0.d0
   normr=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      x(i1,i2,i3,isp) = x(i1,i2,i3,isp) + lv(i1,i2,i3,isp)
      beta=beta+r(i1,i2,i3,isp)*lv(i1,i2,i3,isp)
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential
! correction
      divprod = 0.d0
      do j=1,3
       !divprod = divprod + deps(i1,i2,i3,j)*dx(i1,i2,i3,isp,j)
       divprod = divprod + dlogeps(j,i1,i2,i3)*dx(i1,i2,i3,isp,j)
      end do
      r(i1,i2,i3,isp)=eta*0.25d0*divprod*eps(i1,i2,i3)/pi +(1.d0-eta)*r(i1,i2,i3,isp)
      normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
      !normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)/(eps(i1,i2,i3)**2)
      lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/eps(i1,i2,i3)
     end do
    end do
   end do
   normr=dsqrt(normr/rpoints)

   ratio=normr/normb
   write(18,'(1x,I8,2(1x,e14.7))')i,normr,beta
   write(*,'(1x,I8,2(1x,e14.7))')i,normr,beta
!   call writeroutine(n01,n02,n03,nspden,r,i)
   call writeroutinePot(n01,n02,n03,nspden,x,i,potential)
   call EPS_iter_output_LG(i,0.0_dp,normr,0.0_dp,0.0_dp,0.0_dp)

  end do

   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      b(i1,i2,i3,isp) = x(i1,i2,i3,isp)
     end do
    end do
   end do

  write(*,*)
  write(*,'(1x,a,1x,I8)')'PSD iterations =',i-1
  write(*,'(1x,a,1x,e14.7)')'PSD error =',normr
  write(*,*)
  write(*,*)'Max abs difference between analytic potential and the computed one'
  call writeroutinePot(n01,n02,n03,nspden,b,i-1,potential)
  write(*,*)

  close(unit=18)
  close(unit=38)

  write(*,'(a)')'Termination of Preconditioned Steepest Descent'
  write(*,'(a)')'--------------------------------------------------------------------------------------------'

end subroutine  Prec_Steepest_Descent

subroutine Prec_Steepest_Descent_Inputguess(n01,n02,n03,nspden,iproc,&
     hx,hy,hz,b,bb,acell,eps,dlogeps,nord,pkernel,potential,geocode,multp)
  use yaml_output
  use Poisson_Solver
  use f_utils
  use dynamic_memory

  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden,iproc
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp
  type(coulomb_operator), intent(inout) :: pkernel
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(3,n01,n02,n03), intent(in) :: dlogeps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential
  real(kind=8), dimension(n01,n02,n03,nspden), intent(inout) :: b
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: bb
  character(len=2), intent(in) :: geocode

  real(kind=8), dimension(n01,n02,n03,nspden) :: x,r,p,lv,rvApp,r0
  real(kind=8), dimension(n01,n02,n03,nspden,3) :: dx
  real(kind=8), dimension(n01,n02,n03,3) :: deps
  real(kind=8), parameter :: eta = 0.6d0 !1.0d0 !Like eta in the Polarization Iterative Method.
  integer, parameter :: max_iter = 50
  integer, parameter :: max_scf = 1
  real(kind=8), parameter :: max_ratioex = 1.0d10
  real(kind=8) :: alpha,beta,normb,normr,ratio,k,rpoints
  integer :: i,ii,j,i1,i2,i3,isp,i1_max,i2_max,i3_max,i_scf
  real(kind=8), parameter :: error = 1.0d-12
  real(kind=8), dimension(n01,n02,n03) ::pot_ion,re
  real(kind=8) :: ehartree,offset,pi,divprod,fact,max_val,errorcur,rr,kk
  logical :: PSDig

  pi = 4.d0*datan(1.d0)
  rpoints=product(real([n01,n02,n03],kind=8))

  open(unit=18,file='PSDConvergence.dat',status='unknown')
  open(unit=38,file='MaxAnalysisPSD.dat',status='unknown')

  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,'(a)')'Starting an Inputguess PSD calculation'
  end if

  do i_scf=1,max_scf ! scf loop external

  if (i_scf.eq.1) then
   errorcur=error
   PSDig=.true.
  else
   !PSDig=.false.
   PSDig=.true.
   !errorcur=error
   errorcur=error
  end if

  if (PSDig) then

   call f_memcpy(src=b,dest=x)

   ! Calculate the charge starting from the potential applying the proper
   ! Laplace operator.
   call ApplyLaplace(geocode,n01,n02,n03,nspden,hx,hy,hz,x,rvApp,acell,eps,nord,.false.,multp)

   if (iproc==0) then
    call yaml_mapping_open('Comparison between Generalized Poisson operator and analytical density')
    call writeroutinePot(n01,n02,n03,1,bb,0,rvApp)
    call yaml_mapping_close()
   end if

!   ! Calculate the charge starting from the potential applying the proper
!   Laplace operator.
!   call
!   ApplyLaplace2(geocode,n01,n02,n03,nspden,hx,hy,hz,x,rvApp,acell,eps,dlogeps,nord,.false.,multp)
!
!  if (iproc==0) then
!   call yaml_comment('Comparison between Generalized Poisson operator 2 and
!   analytical density',hfill='-')
!   call writeroutinePot(n01,n02,n03,1,bb,1,rvApp)
!  end if

  end if

  if (iproc==0) then
   write(*,'(a)')'--------------------------------------------------------------------------------------------!'
   write(*,*)'Starting Inputguess iteration ',i_scf
  end if


  write(*,'(a)')'--------------------------------------------------------------------------------------------'
  write(*,'(a)')'Starting Preconditioned Steepest Descent'
  write(*,'(a)')'Starting PSD iteration 1'

!------------------------------------------------------------------------------------
! Set the correction vector for the Generalized Laplace operator

!  call fssnord3DmatNabla3var(n01,n02,n03,nspden,hx,hy,hz,eps,deps,nord,acell)
  if (PSDig) call fssnord3DmatNabla(geocode,n01,n02,n03,nspden,hx,hy,hz,x,dx,nord,acell)

  normb=0.d0
  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     r0(i1,i2,i3,isp)=0.d0
     kk=0.d0
     if (PSDig) then
      kk=1.d0
      divprod = 0.d0
      do j=1,3
       !divprod = divprod + deps(i1,i2,i3,j)*dx(i1,i2,i3,isp,j)
       divprod = divprod + dlogeps(j,i1,i2,i3)*dx(i1,i2,i3,isp,j)
      end do
      r0(i1,i2,i3,isp)=0.25d0*divprod*eps(i1,i2,i3)/pi
     end if
     r(i1,i2,i3,isp)=bb(i1,i2,i3,isp)+r0(i1,i2,i3,isp)
     normb=normb+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
     lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/eps(i1,i2,i3)
    end do
   end do
  end do
  normb=dsqrt(normb/rpoints)

  !call H_potential('G',pkernel,lv,pot_ion,ehartree,0.0_dp,.false.)
  call H_potential('G',pkernel,lv,pot_ion,ehartree,offset,.false.)

!------------------------------------------------------------------------------------
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential
! correction

  call fssnord3DmatNabla(geocode,n01,n02,n03,nspden,hx,hy,hz,lv,dx,nord,acell)

  normr=0.d0
  beta=0.d0
  isp=1
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     x(i1,i2,i3,isp) = lv(i1,i2,i3,isp)
     beta=beta+b(i1,i2,i3,isp)*lv(i1,i2,i3,isp)
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential
! correction
     divprod = 0.d0
     do j=1,3
      !divprod = divprod + deps(i1,i2,i3,j)*dx(i1,i2,i3,isp,j)
      divprod = divprod + dlogeps(j,i1,i2,i3)*dx(i1,i2,i3,isp,j)
     end do
     !r(i1,i2,i3,isp)=eta*0.25d0*divprod*eps(i1,i2,i3)/pi +kk*(1.d0-eta)*r0(i1,i2,i3,isp)*eps(i1,i2,i3)
     !r(i1,i2,i3,isp)=eta*0.25d0*divprod*eps(i1,i2,i3)/pi +kk*(1.d0-eta)*r(i1,i2,i3,isp)
     r(i1,i2,i3,isp)=eta*(0.25d0*divprod*eps(i1,i2,i3)/pi - kk*r0(i1,i2,i3,isp))
     normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
     !normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)/(eps(i1,i2,i3)**2)
     !lv(i1,i2,i3,isp) = eta*0.25d0*divprod/pi
     lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/eps(i1,i2,i3)
    end do
   end do
  end do
  normr=dsqrt(normr/rpoints)

!------------------------------------------------------------------------------------

  ratio=normr/normb

!  call writeroutine(n01,n02,n03,nspden,r,1)
  call writeroutinePot(n01,n02,n03,nspden,potential,0,potential)
  call writeroutinePot(n01,n02,n03,nspden,x,1,potential)

  write(18,'(1x,I8,2(1x,e14.7))')1,normr,beta
  !write(*,'(1x,I8,2(1x,e14.7))')1,normr,beta
  call EPS_iter_output_LG(1,0.0_dp,normr,0.0_dp,0.0_dp,0.0_dp)
  !if (PSDig) stop

  do i=2,max_iter

   if (normr.lt.errorcur) exit
   if (ratio.gt.max_ratioex) exit

   write(*,'(a)')'--------------------------------------------------------------------------------------------'
   write(*,*)'Starting PSD iteration ',i

!  Apply the Preconditioner

   call H_potential('G',pkernel,lv,pot_ion,ehartree,offset,.false.)

   call fssnord3DmatNabla(geocode,n01,n02,n03,nspden,hx,hy,hz,lv,dx,nord,acell)
   !call fssnord3DmatNabla(n01,n02,n03,nspden,lv,dx,nord,acell)

   beta=0.d0
   normr=0.d0
   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      x(i1,i2,i3,isp) = x(i1,i2,i3,isp) + lv(i1,i2,i3,isp)
      beta=beta+r(i1,i2,i3,isp)*lv(i1,i2,i3,isp)
! Apply the Generalized Laplace operator nabla(eps*nabla) to the potential
! correction
      divprod = 0.d0
      do j=1,3
       !divprod = divprod + deps(i1,i2,i3,j)*dx(i1,i2,i3,isp,j)
       divprod = divprod + dlogeps(j,i1,i2,i3)*dx(i1,i2,i3,isp,j)
      end do
      r(i1,i2,i3,isp)=eta*0.25d0*divprod*eps(i1,i2,i3)/pi +(1.d0-eta)*r(i1,i2,i3,isp)
      normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)
      !normr=normr+r(i1,i2,i3,isp)*r(i1,i2,i3,isp)/(eps(i1,i2,i3)**2)
      lv(i1,i2,i3,isp) = r(i1,i2,i3,isp)/eps(i1,i2,i3)
     end do
    end do
   end do
   normr=dsqrt(normr/rpoints)

   ratio=normr/normb
   write(18,'(1x,I8,2(1x,e14.7))')i,normr,beta
   write(*,'(1x,I8,2(1x,e14.7))')i,normr,beta
!   call writeroutine(n01,n02,n03,nspden,r,i)
   call writeroutinePot(n01,n02,n03,nspden,x,i,potential)
   call EPS_iter_output_LG(i,0.0_dp,normr,0.0_dp,0.0_dp,0.0_dp)

  end do ! Loop PSD.

  end do ! Scf external loop.

   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01
      b(i1,i2,i3,isp) = x(i1,i2,i3,isp)
     end do
    end do
   end do

  write(*,*)
  write(*,'(1x,a,1x,I8)')'PSD iterations =',i-1
  write(*,'(1x,a,1x,e14.7)')'PSD error =',normr
  write(*,*)
  write(*,*)'Max abs difference between analytic potential and the computed one'
  call writeroutinePot(n01,n02,n03,nspden,b,i-1,potential)
  write(*,*)

  close(unit=18)
  close(unit=38)

  write(*,'(a)')'Termination of Preconditioned Steepest Descent'
  write(*,'(a)')'--------------------------------------------------------------------------------------------'

end subroutine  Prec_Steepest_Descent_Inputguess

subroutine EPS_iter_output_LG(iter,normb,normr,ratio,alpha,beta)
  !use module_defs, only: dp
  use yaml_output
  use yaml_strings
  implicit none
  integer, intent(in) :: iter
  integer, parameter :: dp=8
  real(dp), intent(in) :: normb,normr,ratio,beta,alpha

  call yaml_mapping_open('Iteration quality',flow=.true.)
  call yaml_comment('Iteration '//iter,hfill='_')
  !write the PCG iteration
  call yaml_map('iter',iter,fmt='(i4)')
  !call yaml_map('rho_norm',normb)
  if (normr/=0.0_dp) call yaml_map('normr',normr,fmt='(1pe16.4)')
  if (ratio /= 0.0_dp) call yaml_map('ratio',ratio,fmt='(1pe16.4)')
  if (alpha /= 0.0_dp) call yaml_map('alpha',alpha,fmt='(1pe16.4)')
  if (beta /= 0.0_dp) call yaml_map('beta',beta,fmt='(1pe16.4)')

  call yaml_mapping_close()
end subroutine EPS_iter_output_LG


subroutine writeroutine(n01,n02,n03,nspden,r,i)

  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  integer, intent(in) :: i
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: r
  integer :: i1,i2,i3,j,i1_max,i2_max,i3_max,jj
  real(kind=8) :: max_val,fact

     j=i+40
     jj=i+100

    if (i.le.50) then

     i3=n03/2
     do i2=1,n02
      do i1=1,n01
       write(j,'(2(1x,I4),1x,e14.7)') i1,i2,r(i1,i2,i3,1)
      end do
      write(j,*) 
     end do

      do i1=1,n01
!       write(jj,'(1x,I4,2(1x,e22.15))') i1,r(i1,n02/2,i3,1),r(i1,1,i3,1)
       write(jj,'(1x,I4,2(1x,e22.15))') i1,r(i1,i1,i1,1),r(i1,n02/2,i3,1)
      end do

    end if

      max_val = 0.d0
      i1_max = 1
      i2_max = 1
      i3_max = 1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               fact=abs(r(i1,i2,i3,1))
               if (max_val < fact) then
                  max_val = fact
                  i1_max = i1
                  i2_max = i2
                  i3_max = i3
               end if
            end do
         end do
      end do
      write(39,'(4(1x,I4),4(1x,e22.15))')i,i1_max,i2_max,i3_max,max_val,&
           r(n01/2,n02/2,n03/2,1),r(2,n02/2,n03/2,1),r(10,n02/2,n03/2,1)

end subroutine writeroutine

subroutine writeroutinePot(n01,n02,n03,nspden,ri,i,potential)
  use yaml_output
  use dynamic_memory
  use f_utils
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  integer, intent(in) :: i
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: ri
  real(kind=8), dimension(n01,n02,n03),intent(in) :: potential
  !automatic array, to be check is stack poses problem
  real(kind=8), dimension(:,:,:,:), allocatable :: re
  integer :: i1,i2,i3,j,i1_max,i2_max,i3_max,jj,unt
  real(kind=8) :: max_val,fact
  character(len=20) :: str
  re=f_malloc([n01,n02,n03,nspden],id='re')
  
      max_val = 0.d0
      i1_max = 1
      i2_max = 1
      i3_max = 1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               re(i1,i2,i3,1) = ri(i1,i2,i3,1) - potential(i1,i2,i3)
               fact=abs(re(i1,i2,i3,1))
               if (max_val < fact) then
                  max_val = fact
                  i1_max = i1
                  i2_max = i2
                  i3_max = i3
               end if
            end do
         end do
      end do
      write(38,'(4(1x,I4),2(1x,e22.15))')i,i1_max,i2_max,i3_max,max_val,&
           re(n01/2,n02/2,n03/2,1)
!!$      write(38,'(4(1x,I4),4(1x,e22.15))')i,i1_max,i2_max,i3_max,max_val,&
!!$           re(n01/2,n02/2,n03/2,1),re(2,n02/2,n03/2,1),re(10,n02/2,n03/2,1)
      !write(*,'(4(1x,I4),4(1x,e22.15))')i,i1_max,i2_max,i3_max,max_val,&
      !     re(n01/2,n02/2,n03/2,1),re(2,n02/2,n03/2,1),re(10,n02/2,n03/2,1)
      if (max_val == 0.d0) then
         call yaml_map('Inf. Norm difference with reference',0.d0)
      else
         call yaml_mapping_open('Inf. Norm difference with reference')
         call yaml_map('Value',max_val,fmt='(1pe22.15)')
         call yaml_map('Point',[i1_max,i2_max,i3_max],fmt='(i4)')
         call yaml_map('Some values',[re(n01/2,n02/2,n03/2,1),re(2,n02/2,n03/2,1),re(10,n02/2,n03/2,1)],&
              fmt='(1pe22.15)')
         call yaml_mapping_close()
      end if
      
      unt=f_get_free_unit(21)
      call f_open_file(unt,file='final.dat')
      !i1=n01/2
      i1=(n01-1)/2+1
      do i2=1,n02
         do i3=1,n03
            write(unt,'(2(1x,I4),2(1x,e14.7))')i2,i3,ri(i1,i2,i3,1),potential(i1,i2,i3)
         end do
      end do
      call f_close(unt)

      unt=f_get_free_unit(22)
      i1=(n01-1)/2+1
      i3=(n03-1)/2+1
      call f_open_file(unt,file='final_line.dat')
      !write (str, *) i
      !str = adjustl(str)
      !call f_open_file(unt,file='final_line'//trim(str)//'.dat')
      do i2=1,n02
       write(unt,'(1x,I8,3(1x,e22.15))') i2,ri(i1,i2,i3,1),potential(i1,i2,i3),re(i1,i2,i3,1)
      end do
      call f_close(unt)

      call f_free(re)
end subroutine writeroutinePot

subroutine FluxSurface(n01,n02,n03,nspden,hx,hy,hz,x,acell,eps,nord)
  use dynamic_memory
  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell
  real(kind=8), dimension(n01,n02,n03), intent(in) :: x
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(:,:,:,:,:), allocatable :: dx
  real(kind=8) :: pi,flux
  integer :: i1,i2,i3,isp,i
  
  dx=f_malloc([n01,n02,n03,nspden,3],id='dx')

  pi = 4.d0*datan(1.d0)

   call fssnord3DmatNabla("S ",n01,n02,n03,nspden,hx,hy,hz,x,dx,nord,acell)

     flux=0.d0
      isp=1
      do i3=1,n03
       do i2=1,n02
        do i1=1,n01

         if (i1.eq.1) then
          flux=flux-eps(i1,i2,i3)*dx(i1,i2,i3,isp,1)*hy*hz
         end if
         if (i1.eq.n01) then
          flux=flux+eps(i1,i2,i3)*dx(i1,i2,i3,isp,1)*hy*hz
         end if

         if (i2.eq.1) then
          flux=flux-eps(i1,i2,i3)*dx(i1,i2,i3,isp,2)*hx*hz
         end if
         if (i2.eq.n02) then
          flux=flux+eps(i1,i2,i3)*dx(i1,i2,i3,isp,2)*hx*hz
         end if

         if (i3.eq.1) then
          flux=flux-eps(i1,i2,i3)*dx(i1,i2,i3,isp,3)*hx*hy
         end if
         if (i3.eq.n03) then
          flux=flux+eps(i1,i2,i3)*dx(i1,i2,i3,isp,3)*hx*hy
         end if

        end do
       end do
      end do

    write(*,'(1x,a,1x,e14.7)')'Surface flux is',flux

    call f_free(dx)

end subroutine FluxSurface 

subroutine  ISF_product(geocode,itype,n01,n02,n03,x,c,y)
  use yaml_output
  use dynamic_memory
  implicit none
  character(len=2), intent(in) :: geocode
  integer, intent(in) :: itype
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  real(kind=8), dimension(n01,n02,n03), intent(in) :: x
  real(kind=8), dimension(n01,n02,n03), intent(in) :: c
  real(kind=8), dimension(n01,n02,n03), intent(out) :: y

  ! Local variables.
  real(kind=8), dimension(:,:,:), allocatable :: xx,cc,yy,ycur,re
  real(kind=8), dimension(:), allocatable :: coarse1,finer1,coarse2,finer2
  integer :: i,i1,i2,i3,n01d,n02d,n03d,nd,nt
  integer :: i1_max,i2_max,i3_max
  real(kind=8) :: max_val,fact

  n01d=2*n01
  n02d=2*n02
  n03d=2*n03
  xx=f_malloc([n01d,n02d,n03d],id='xx')
  cc=f_malloc([n01d,n02d,n03d],id='cc')
  yy=f_malloc([n01d,n02d,n03d],id='yy')
  ycur=f_malloc([n01,n02,n03],id='ycur')
  re=f_malloc([n01,n02,n03],id='re')

  !Input grid product.
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     y(i1,i2,i3)=x(i1,i2,i3)*c(i1,i2,i3)
     ycur(i1,i2,i3)=x(i1,i2,i3)*c(i1,i2,i3)
    end do
   end do
  end do

!---------------------------------------------------------
  ! 3D-backward transformation.
  ! x-axis
  nd=n01d
  nt=n01d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  coarse2=f_malloc(nd,id='coarse2')
  finer2=f_malloc(nd,id='finer2')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
   coarse2(i)=0.d0
   finer2(i)=0.d0
  end do
  do i3=1,n03
   do i2=1,n02
    coarse1(1:n01)=x(1:n01,i2,i3)
    coarse2(1:n01)=c(1:n01,i2,i3)
    call back_trans_ISC(itype,nd,nt,coarse1,finer1,coarse2,finer2)
    xx(1:n01d,i2,i3)=finer1(1:n01d)
    cc(1:n01d,i2,i3)=finer2(1:n01d)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
  call f_free(coarse2)
  call f_free(finer2)
  ! y-axis
  nd=n02d
  nt=n02d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  coarse2=f_malloc(nd,id='coarse2')
  finer2=f_malloc(nd,id='finer2')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
   coarse2(i)=0.d0
   finer2(i)=0.d0
  end do
  do i3=1,n03
   do i1=1,n01
    coarse1(1:n02)=x(i1,1:n02,i3)
    coarse2(1:n02)=c(i1,1:n02,i3)
    call back_trans_ISC(itype,nd,nt,coarse1,finer1,coarse2,finer2)
    xx(i1,1:n02d,i3)=finer1(1:n02d)
    cc(i1,1:n02d,i3)=finer2(1:n02d)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
  call f_free(coarse2)
  call f_free(finer2)
  ! z-axis
  nd=n03d
  nt=n03d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  coarse2=f_malloc(nd,id='coarse2')
  finer2=f_malloc(nd,id='finer2')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
   coarse2(i)=0.d0
   finer2(i)=0.d0
  end do
  do i2=1,n02
   do i1=1,n01
    coarse1(1:n03)=x(i1,i2,1:n03)
    coarse2(1:n03)=c(i1,i2,1:n03)
    call back_trans_ISC(itype,nd,nt,coarse1,finer1,coarse2,finer2)
    xx(i1,i2,1:n03d)=finer1(1:n03d)
    cc(i1,i2,1:n03d)=finer2(1:n03d)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
  call f_free(coarse2)
  call f_free(finer2)
!---------------------------------------------------------

  !Finer grid product.
  do i3=1,n03d
   do i2=1,n02d
    do i1=1,n01d
     yy(i1,i2,i3)=xx(i1,i2,i3)*cc(i1,i2,i3)
    end do
   end do
  end do

!---------------------------------------------------------
  ! 3D-forward transformation.
  ! x-axis
  nd=n01d
  nt=n01d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
  end do
  do i3=1,n03
   do i2=1,n02
    finer1(1:nd)=yy(1:nd,i2,i3)
    call for_trans_ISC(itype,nd,nt,finer1,coarse1)
    y(1:n01,i2,i3)=coarse1(1:n01)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
  ! y-axis
  nd=n02d
  nt=n02d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
  end do
  do i3=1,n03
   do i1=1,n01
    finer1(1:nd)=yy(i1,1:nd,i3)
    call for_trans_ISC(itype,nd,nt,finer1,coarse1)
    y(i1,1:n02,i3)=coarse1(1:n02)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
  ! z-axis
  nd=n03d
  nt=n03d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
  end do
  do i2=1,n02
   do i1=1,n01
    finer1(1:nd)=yy(i1,i2,1:nd)
    call for_trans_ISC(itype,nd,nt,finer1,coarse1)
    y(i1,i2,1:n03)=coarse1(1:n03)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
!---------------------------------------------------------

!---------------------------------------------------------
!Checking difference between coarse/finer grid product.

  max_val = 0.d0
  i1_max = 1
  i2_max = 1
  i3_max = 1
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           re(i1,i2,i3) = ycur(i1,i2,i3) - y(i1,i2,i3)
           fact=abs(re(i1,i2,i3))
           if (max_val < fact) then
              max_val = fact
              i1_max = i1
              i2_max = i2
              i3_max = i3
           end if
        end do
     end do
  end do
  call yaml_map('Diff c-f product',max_val,fmt='(1pe22.15)')
!---------------------------------------------------------


  call f_free(xx)
  call f_free(cc)
  call f_free(yy)
  call f_free(ycur)
  call f_free(re)

end subroutine ISF_product

subroutine ISF_product_new(geocode,itype,n01,n02,n03,hx,hy,hz,nord,acell,x,c,y,k)
  use yaml_output
  use dynamic_memory
  implicit none
  character(len=2), intent(in) :: geocode
  integer, intent(in) :: itype
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  real(kind=8), intent(in) :: hx,hy,hz,acell
  integer, intent(in) :: nord
  real(kind=8), dimension(n01,n02,n03), intent(in) :: x
  real(kind=8), dimension(n01,n02,n03), intent(in) :: c
  real(kind=8), dimension(n01,n02,n03), intent(out) :: y
  real(kind=8), intent(out) :: k

  ! Local variables.
  real(kind=8), dimension(:,:,:), allocatable :: xx,cc,yy,ycur,re
  real(kind=8), dimension(:,:,:,:), allocatable :: nxx,ncc
  real(kind=8), dimension(3) :: hgrids
  real(kind=8), dimension(:), allocatable :: coarse1,finer1,coarse2,finer2
  integer :: i,i1,i2,i3,i11,i22,i33,n01d,n02d,n03d,nd,nt
  integer :: i1_max,i2_max,i3_max
  real(kind=8) :: max_val,fact,xxx,ccc

  n01d=2*n01
  n02d=2*n02
  n03d=2*n03
  xx=f_malloc([n01d,n02d,n03d],id='xx')
  cc=f_malloc([n01d,n02d,n03d],id='cc')
  yy=f_malloc([n01d,n02d,n03d],id='yy')
  ycur=f_malloc([n01,n02,n03],id='ycur')
  re=f_malloc([n01,n02,n03],id='re')
  nxx=f_malloc([n01d,n02d,n03d,3],id='nxx')
  ncc=f_malloc([n01d,n02d,n03d,3],id='ncc')

  hgrids(1)=hx
  hgrids(2)=hy
  hgrids(3)=hz

  !Input grid product.
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     ycur(i1,i2,i3)=x(i1,i2,i3)*c(i1,i2,i3)
     xx(i1,i2,i3)=x(i1,i2,i3)
     cc(i1,i2,i3)=c(i1,i2,i3)
     i11=i1+n01
     i22=i2+n02
     i33=i3+n03
     xx(i11,i22,i33)=0.d0
     cc(i11,i22,i33)=0.d0
    end do
   end do
  end do

!---------------------------------------------------------
  ! 3D-backward transformation.
  ! z-axis
  nd=n03d
  nt=n03d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  coarse2=f_malloc(nd,id='coarse2')
  finer2=f_malloc(nd,id='finer2')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
   coarse2(i)=0.d0
   finer2(i)=0.d0
  end do
  do i2=1,n02
   do i1=1,n01
    coarse1(1:n03)=xx(i1,i2,1:n03)
    coarse2(1:n03)=cc(i1,i2,1:n03)
    call back_trans_ISC(itype,nd,nt,coarse1,finer1,coarse2,finer2)
    xx(i1,i2,1:n03d)=finer1(1:n03d)
    cc(i1,i2,1:n03d)=finer2(1:n03d)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
  call f_free(coarse2)
  call f_free(finer2)
  ! y-axis
  nd=n02d
  nt=n02d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  coarse2=f_malloc(nd,id='coarse2')
  finer2=f_malloc(nd,id='finer2')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
   coarse2(i)=0.d0
   finer2(i)=0.d0
  end do
  do i3=1,n03d
   do i1=1,n01
    coarse1(1:n02)=xx(i1,1:n02,i3)
    coarse2(1:n02)=cc(i1,1:n02,i3)
    call back_trans_ISC(itype,nd,nt,coarse1,finer1,coarse2,finer2)
    xx(i1,1:n02d,i3)=finer1(1:n02d)
    cc(i1,1:n02d,i3)=finer2(1:n02d)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
  call f_free(coarse2)
  call f_free(finer2)
  ! x-axis
  nd=n01d
  nt=n01d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  coarse2=f_malloc(nd,id='coarse2')
  finer2=f_malloc(nd,id='finer2')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
   coarse2(i)=0.d0
   finer2(i)=0.d0
  end do
  do i3=1,n03d
   do i2=1,n02d
    coarse1(1:n01)=xx(1:n01,i2,i3)
    coarse2(1:n01)=cc(1:n01,i2,i3)
    call back_trans_ISC(itype,nd,nt,coarse1,finer1,coarse2,finer2)
    xx(1:n01d,i2,i3)=finer1(1:n01d)
    cc(1:n01d,i2,i3)=finer2(1:n01d)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
  call f_free(coarse2)
  call f_free(finer2)
!---------------------------------------------------------

!  call fssnord3DmatNabla(geocode,n01d,n02d,n03d,1,hx,hy,hz,xx,nxx,nord,acell)
!  call fssnord3DmatNabla(geocode,n01d,n02d,n03d,1,hx,hy,hz,cc,ncc,nord,acell)

  !Finer grid product.
  k=0.d0
  do i3=1,n03d
   do i2=1,n02d
    do i1=1,n01d
!     yy(i1,i2,i3)=xx(i1,i2,i3)*cc(i1,i2,i3)
!     xxx=0.d0
!     ccc=0.d0
!     do i=1,3
!      xxx=xxx+nxx(i1,i2,i3,i)*hgrids(i)
!      ccc=ccc+ncc(i1,i2,i3,i)*hgrids(i)
!     end do
!     xx(i1,i2,i3)=xx(i1,i2,i3)+0.5d0*xxx
!     cc(i1,i2,i3)=cc(i1,i2,i3)+0.5d0*ccc
     yy(i1,i2,i3)=xx(i1,i2,i3)*cc(i1,i2,i3)
     k=k+yy(i1,i2,i3)
    end do
   end do
  end do

!---------------------------------------------------------
  ! 3D-forward transformation.
  ! x-axis
  nd=n01d
  nt=n01d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
  end do
  do i3=1,n03d
   do i2=1,n02d
    finer1(1:nd)=yy(1:nd,i2,i3)
    call for_trans_ISC(itype,nd,nt,finer1,coarse1)
    yy(1:n01,i2,i3)=coarse1(1:n01)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
  ! y-axis
  nd=n02d
  nt=n02d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
  end do
  do i3=1,n03d
   do i1=1,n01
    finer1(1:nd)=yy(i1,1:nd,i3)
    call for_trans_ISC(itype,nd,nt,finer1,coarse1)
    yy(i1,1:n02,i3)=coarse1(1:n02)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
  ! z-axis
  nd=n03d
  nt=n03d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
  end do
  do i2=1,n02
   do i1=1,n01
    finer1(1:nd)=yy(i1,i2,1:nd)
    call for_trans_ISC(itype,nd,nt,finer1,coarse1)
    yy(i1,i2,1:n03)=coarse1(1:n03)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
!---------------------------------------------------------

  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     y(i1,i2,i3)=yy(i1,i2,i3)
    end do
   end do
  end do

!---------------------------------------------------------
!Checking difference between coarse/finer grid product.

  max_val = 0.d0
  i1_max = 1
  i2_max = 1
  i3_max = 1
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           re(i1,i2,i3) = ycur(i1,i2,i3) - y(i1,i2,i3)
           fact=abs(re(i1,i2,i3))
           if (max_val < fact) then
              max_val = fact
              i1_max = i1
              i2_max = i2
              i3_max = i3
           end if
        end do
     end do
  end do
  call yaml_map('Diff c-f product',max_val,fmt='(1pe22.15)')
!---------------------------------------------------------

  call f_free(xx)
  call f_free(cc)
  call f_free(yy)
  call f_free(ycur)
  call f_free(re)
  call f_free(nxx)
  call f_free(ncc)

end subroutine ISF_product_new

subroutine for_trans_ISF_3D(geocode,itype,n01c,n02c,n03c,x,y)
  use yaml_output
  use dynamic_memory
  implicit none
  character(len=2), intent(in) :: geocode
  integer, intent(in) :: itype
  integer, intent(in) :: n01c
  integer, intent(in) :: n02c
  integer, intent(in) :: n03c
  real(kind=8), dimension(n01c,n02c,n03c), intent(in) :: x
  real(kind=8), dimension(n01c/2,n02c/2,n03c/2), intent(out) :: y

  ! Local variables.
  real(kind=8), dimension(:,:,:), allocatable :: xx
  real(kind=8), dimension(:), allocatable :: coarse1,finer1
  integer :: i,i1,i2,i3,i11,i22,i33,n01d,n02d,n03d,nd,nt
  integer :: i1_max,i2_max,i3_max,n01,n02,n03
  real(kind=8) :: max_val,fact,xxx,ccc

  n01=n01c/2
  n02=n02c/2
  n03=n03c/2
  n01d=n01c
  n02d=n02c
  n03d=n03c
  
  xx=f_malloc([n01d,n02d,n03d],id='xx')

  do i3=1,n03d
   do i2=1,n02d
    do i1=1,n01d
     xx(i1,i2,i3)=x(i1,i2,i3)
    end do
   end do
  end do

!---------------------------------------------------------
  ! 3D-forward transformation.
  ! x-axis
  nd=n01d
  nt=n01d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
  end do
  do i3=1,n03d
   do i2=1,n02d
    finer1(1:nd)=xx(1:nd,i2,i3)
    call for_trans_ISC(itype,nd,nt,finer1,coarse1)
    xx(1:n01,i2,i3)=coarse1(1:n01)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
  ! y-axis
  nd=n02d
  nt=n02d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
  end do
  do i3=1,n03d
   do i1=1,n01
    finer1(1:nd)=xx(i1,1:nd,i3)
    call for_trans_ISC(itype,nd,nt,finer1,coarse1)
    xx(i1,1:n02,i3)=coarse1(1:n02)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
  ! z-axis
  nd=n03d
  nt=n03d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
  end do
  do i2=1,n02
   do i1=1,n01
    finer1(1:nd)=xx(i1,i2,1:nd)
    call for_trans_ISC(itype,nd,nt,finer1,coarse1)
    xx(i1,i2,1:n03)=coarse1(1:n03)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
!---------------------------------------------------------

  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     y(i1,i2,i3)=xx(i1,i2,i3)
    end do
   end do
  end do

!---------------------------------------------------------
  
  call f_free(xx)

end subroutine for_trans_ISF_3D

subroutine back_trans_ISF_3D(geocode,itype,n01c,n02c,n03c,x,y)
  use yaml_output
  use dynamic_memory
  implicit none
  character(len=2), intent(in) :: geocode
  integer, intent(in) :: itype
  integer, intent(in) :: n01c
  integer, intent(in) :: n02c
  integer, intent(in) :: n03c
  real(kind=8), dimension(n01c/2,n02c/2,n03c/2), intent(in) :: x
  real(kind=8), dimension(n01c,n02c,n03c), intent(out) :: y

  ! Local variables.
  real(kind=8), dimension(:,:,:), allocatable :: xx
  real(kind=8), dimension(:), allocatable :: coarse1,finer1
  integer :: i,i1,i2,i3,i11,i22,i33,n01d,n02d,n03d,nd,nt
  integer :: i1_max,i2_max,i3_max,n01,n02,n03
  real(kind=8) :: max_val,fact,xxx,ccc

  n01=n01c/2
  n02=n02c/2
  n03=n03c/2
  n01d=n01c
  n02d=n02c
  n03d=n03c
  
  xx=f_malloc([n01d,n02d,n03d],id='xx')
  
  do i3=1,n03
   do i2=1,n02
    do i1=1,n01
     xx(i1,i2,i3)=x(i1,i2,i3)
     i11=i1+n01
     i22=i2+n02
     i33=i3+n03
     xx(i11,i22,i33)=0.d0
    end do
   end do
  end do

!---------------------------------------------------------
  ! 3D-backward transformation.
  ! z-axis
  nd=n03d
  nt=n03d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
  end do
  do i2=1,n02
   do i1=1,n01
    coarse1(1:n03)=xx(i1,i2,1:n03)
    call back_trans_ISC1vect(itype,nd,nt,coarse1,finer1)
    xx(i1,i2,1:n03d)=finer1(1:n03d)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
  ! y-axis
  nd=n02d
  nt=n02d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
  end do
  do i3=1,n03d
   do i1=1,n01
    coarse1(1:n02)=xx(i1,1:n02,i3)
    call back_trans_ISC1vect(itype,nd,nt,coarse1,finer1)
    xx(i1,1:n02d,i3)=finer1(1:n02d)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
  ! x-axis
  nd=n01d
  nt=n01d
  coarse1=f_malloc(nd,id='coarse1')
  finer1=f_malloc(nd,id='finer1')
  do i=1,nd
   coarse1(i)=0.d0
   finer1(i)=0.d0
  end do
  do i3=1,n03d
   do i2=1,n02d
    coarse1(1:n01)=xx(1:n01,i2,i3)
    call back_trans_ISC1vect(itype,nd,nt,coarse1,finer1)
    xx(1:n01d,i2,i3)=finer1(1:n01d)
   end do
  end do
  call f_free(coarse1)
  call f_free(finer1)
!---------------------------------------------------------

  do i3=1,n03d
   do i2=1,n02d
    do i1=1,n01d
     y(i1,i2,i3)=xx(i1,i2,i3)
    end do
   end do
  end do

  call f_free(xx)

end subroutine back_trans_ISF_3D

subroutine back_trans_ISC1vect(itype,nd,nt,x1,y1)
! backward wavelet transform
! nd: length of data set
! nt length of data in data set to be transformed
! m filter length (m has to be even!)
! x input data, y output data
  implicit real*8 (a-h,o-z)
  integer, intent(in) :: itype
  dimension x1(0:nd-1),y1(0:nd-1)

  ! Allowed ISF 2,4,6,8,14,16,20,24,30,40,50,60,100
  !include 'lazy_8_lifted.inc'
  !include 'lazy_8.inc'
  include 'lazy_16.inc'
  !include 'lazy_30.inc'

        do 100,i=0,nt/2-1
        y1(2*i+0)=0.d0
        y1(2*i+1)=0.d0

        do 50,j=-m/2,m/2-1

! periodically wrap index if necessary
        ind=i-j
99      continue
        if (ind.lt.0) then
        ind=ind+nt/2
        goto 99
        endif
        if (ind.ge.nt/2) then
        ind=ind-nt/2
        goto 99
        endif

        y1(2*i+0)=y1(2*i+0) + ch(2*j-0)*x1(ind)+cg(2*j-0)*x1(ind+nt/2)
        y1(2*i+1)=y1(2*i+1) + ch(2*j+1)*x1(ind)+cg(2*j+1)*x1(ind+nt/2)

50      continue

100     continue

        return
end subroutine back_trans_ISC1vect

subroutine back_trans_ISC(itype,nd,nt,x1,y1,x2,y2)
! backward wavelet transform
! nd: length of data set
! nt length of data in data set to be transformed
! m filter length (m has to be even!)
! x input data, y output data
  implicit real*8 (a-h,o-z)
  integer, intent(in) :: itype
  dimension x1(0:nd-1),y1(0:nd-1)
  dimension x2(0:nd-1),y2(0:nd-1)

  ! Allowed ISF 2,4,6,8,14,16,20,24,30,40,50,60,100
  !include 'lazy_8_lifted.inc'
  !include 'lazy_8.inc'
  include 'lazy_16.inc'
  !include 'lazy_30.inc'

        do 100,i=0,nt/2-1
        y1(2*i+0)=0.d0
        y1(2*i+1)=0.d0

        y2(2*i+0)=0.d0
        y2(2*i+1)=0.d0

        do 50,j=-m/2,m/2-1

! periodically wrap index if necessary
        ind=i-j
99      continue
        if (ind.lt.0) then
        ind=ind+nt/2
        goto 99
        endif
        if (ind.ge.nt/2) then
        ind=ind-nt/2
        goto 99
        endif

        y1(2*i+0)=y1(2*i+0) + ch(2*j-0)*x1(ind)+cg(2*j-0)*x1(ind+nt/2)
        y1(2*i+1)=y1(2*i+1) + ch(2*j+1)*x1(ind)+cg(2*j+1)*x1(ind+nt/2)

        y2(2*i+0)=y2(2*i+0) + ch(2*j-0)*x2(ind)+cg(2*j-0)*x2(ind+nt/2)
        y2(2*i+1)=y2(2*i+1) + ch(2*j+1)*x2(ind)+cg(2*j+1)*x2(ind+nt/2)
50      continue

100     continue

        return
end subroutine back_trans_ISC

subroutine for_trans_ISC(itype,nd,nt,x,y)
! forward wavelet transform
! nd: length of data set
! nt length of data in data set to be transformed
! m filter length (m has to be even!)
! x input data, y output data
  implicit real*8 (a-h,o-z)
  integer, intent(in) :: itype
  dimension x(0:nd-1),y(0:nd-1)
  
  ! Allowed ISF 2,4,6,8,14,16,20,24,30,40,50,60,100
  !include 'lazy_8_lifted.inc'
  !include 'lazy_8.inc'
  include 'lazy_16.inc'
  !include 'lazy_30.inc'

        do 100,i=0,nt/2-1
        y(     i)=0.d0
        y(nt/2+i)=0.d0

        do 50,j=-m+1,m

! periodically wrap index if necessary
        ind=j+2*i
99      continue
        if (ind.lt.0) then
        ind=ind+nt
        goto 99
        endif
        if (ind.ge.nt) then
        ind=ind-nt
        goto 99
        endif

        y(     i)=y(     i)+cht(j)*x(ind)
        y(nt/2+i)=y(nt/2+i)+cgt(j)*x(ind)
50      continue

100     continue

        return
end subroutine for_trans_ISC

subroutine ApplyLaplace(geocode,n01,n02,n03,nspden,hx,hy,hz,x,y,acell,eps,nord,lin_PB,multp)
  use dynamic_memory
  implicit none
  character(len=2), intent(in) :: geocode
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: x
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: y
  logical, intent(in) :: lin_PB

  ! Local variables.
  real(kind=8), dimension(:,:,:,:), allocatable :: ddx
  real(kind=8), dimension(:,:,:,:,:), allocatable :: dx
  real(kind=8), dimension(:,:,:,:), allocatable :: deps
  real(kind=8) :: pi,switch
  integer :: i1,i2,i3,isp,i
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8) :: PB_charge

  pi = 4.d0*datan(1.d0)   

  ddx=f_malloc([n01,n02,n03,nspden],id='ddx')
  dx=f_malloc([n01,n02,n03,nspden,3],id='dx')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  !write(*,*)geocode
  call fssnord3DmatNabla(geocode,n01,n02,n03,nspden,hx,hy,hz,x,dx,nord,acell)

      isp=1
      do i3=1,n03
       do i2=1,n02
        do i1=1,n01
         do i=1,3
          dx(i1,i2,i3,isp,i)=eps(i1,i2,i3)*dx(i1,i2,i3,isp,i)
         end do
        end do
       end do
      end do

   call fssnord3DmatDiv(geocode,n01,n02,n03,nspden,hx,hy,hz,dx,y,nord,acell)

   y(:,:,:,:)=-y(:,:,:,:)/(4.d0*pi)

   if (lin_PB) then
    isp=1
    do i3=1,n03
     do i2=1,n02
      do i1=1,n01
       y(i1,i2,i3,isp) = y(i1,i2,i3,isp) - ((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(x(i1,i2,i3,isp))
!       y(i1,i2,i3,isp) = y(i1,i2,i3,isp) + ((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dsinh(multp*x(i1,i2,i3,isp))
!       y(i1,i2,i3,isp) = y(i1,i2,i3,isp) + ((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dtanh(multp*x(i1,i2,i3,isp))
!       y(i1,i2,i3,isp) = y(i1,i2,i3,isp) + ((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*x(i1,i2,i3,isp)
!       y(i1,i2,i3,isp) = y(i1,i2,i3,isp) + ((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*(x(i1,i2,i3,isp)**2)
      end do
     end do
    end do
   end if

   call f_free(deps)
   call f_free(ddx)
   call f_free(dx)

end subroutine ApplyLaplace

subroutine ApplyLaplace_ISF(geocode,n01,n02,n03,nspden,hx,hy,hz,x,y,acell,eps,nord,lin_PB,multp)
  use yaml_output
  use dynamic_memory
  implicit none
  character(len=2), intent(in) :: geocode
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: x
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: y
  logical, intent(in) :: lin_PB

  ! Local variables.
  real(kind=8), dimension(:,:,:,:), allocatable :: ddx,xx,yy,ycur
  real(kind=8), dimension(:,:,:,:,:), allocatable :: dx,dxx
  real(kind=8), dimension(:,:,:,:), allocatable :: deps
  real(kind=8), dimension(:,:,:), allocatable :: eeps,re
  real(kind=8) :: pi,switch,hxd,hyd,hzd
  integer :: i1,i2,i3,isp,i,n01d,n02d,n03d,i1_max,i2_max,i3_max
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8) :: PB_charge,max_val,fact

  pi = 4.d0*datan(1.d0)   
  n01d=n01*2
  n02d=n02*2
  n03d=n03*2
  hxd=hx/2.d0
  hyd=hy/2.d0
  hzd=hz/2.d0
  !hxd=acell/real(n01d-1,kind=8)
  !hyd=acell/real(n02d-1,kind=8)
  !hzd=acell/real(n03d-1,kind=8)

  ddx=f_malloc([n01,n02,n03,nspden],id='ddx')
  ycur=f_malloc([n01,n02,n03,nspden],id='ycur')
  xx=f_malloc([n01d,n02d,n03d,nspden],id='xx')
  yy=f_malloc([n01d,n02d,n03d,nspden],id='yy')
  dx=f_malloc([n01,n02,n03,nspden,3],id='dx')
  dxx=f_malloc([n01d,n02d,n03d,nspden,3],id='dxx')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  eeps=f_malloc([n01d,n02d,n03d],id='eeps')
  re=f_malloc([n01,n02,n03],id='re')
  !write(*,*)geocode
  call fssnord3DmatNabla(geocode,n01,n02,n03,nspden,hx,hy,hz,x,dx,nord,acell)

      isp=1
      do i3=1,n03
       do i2=1,n02
        do i1=1,n01
         do i=1,3
          dx(i1,i2,i3,isp,i)=eps(i1,i2,i3)*dx(i1,i2,i3,isp,i)
         end do
        end do
       end do
      end do

   call fssnord3DmatDiv(geocode,n01,n02,n03,nspden,hx,hy,hz,dx,y,nord,acell)

   y(:,:,:,:)=-y(:,:,:,:)/(4.d0*pi)
   ycur(:,:,:,:)=y(:,:,:,:)

!------------------------------------------------------------------------------------
   !Apply the improved ISF Generalized Poisson operator.

   !call H_potential('G',pkernel,lvc,pot_ion,ehartree,0.0_dp,.false.)
   call back_trans_ISF_3D(geocode,16,n01d,n02d,n03d,x,xx)
   call back_trans_ISF_3D(geocode,16,n01d,n02d,n03d,eps,eeps)

   call fssnord3DmatNabla(geocode,n01d,n02d,n03d,nspden,hxd,hyd,hzd,xx,dxx,nord,acell)

      isp=1
      do i3=1,n03d
       do i2=1,n02d
        do i1=1,n01d
         do i=1,3
          dxx(i1,i2,i3,isp,i)=eeps(i1,i2,i3)*dxx(i1,i2,i3,isp,i)
         end do
        end do
       end do
      end do

   call fssnord3DmatDiv(geocode,n01d,n02d,n03d,nspden,hxd,hyd,hzd,dxx,yy,nord,acell)

   call for_trans_ISF_3D(geocode,16,n01d,n02d,n03d,yy,y)

   y(:,:,:,:)=-y(:,:,:,:)/(4.d0*pi)
   yy(:,:,:,:)=-yy(:,:,:,:)/(4.d0*pi)

!---------------------------------------------------------
!Checking difference between coarse/finer grid product.

  max_val = 0.d0
  i1_max = 1
  i2_max = 1
  i3_max = 1
  isp=1
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           !re(i1,i2,i3) = ycur(i1,i2,i3,1) - y(i1,i2,i3,1)
           re(i1,i2,i3) = ycur(i1,i2,i3,1) - yy(2*i1-1,2*i2-1,2*i3-1,1)
           fact=abs(re(i1,i2,i3))
           if (max_val < fact) then
              max_val = fact
              i1_max = i1
              i2_max = i2
              i3_max = i3
           end if
        end do
     end do
  end do
  call yaml_map('Diff ISF operator c-f ISF grid',max_val,fmt='(1pe22.15)')
  write(*,*)i1_max,i2_max,i3_max
!---------------------------------------------------------

!------------------------------------------------------------------------------------
   !Apply linear PB operator.
   if (lin_PB) then
    isp=1
    do i3=1,n03
     do i2=1,n02
      do i1=1,n01
       y(i1,i2,i3,isp) = y(i1,i2,i3,isp) - ((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(x(i1,i2,i3,isp))
      end do
     end do
    end do
   end if

   call f_free(ddx)
   call f_free(ycur)
   call f_free(xx)
   call f_free(yy)
   call f_free(dx)
   call f_free(dxx)
   call f_free(deps)
   call f_free(eeps)
   call f_free(re)

end subroutine ApplyLaplace_ISF

subroutine ApplyLaplace_corr(geocode,n01,n02,n03,nspden,hx,hy,hz,x,y,acell,eps,corr,oneosqrteps,nord,lin_PB,multp)
  use yaml_output
  use dynamic_memory
  implicit none
  character(len=2), intent(in) :: geocode
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: x
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: corr
  real(kind=8), dimension(n01,n02,n03), intent(in) :: oneosqrteps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: y
  logical, intent(in) :: lin_PB

  ! Local variables.
  real(kind=8), dimension(:,:,:,:), allocatable :: ddx,xx,yy,ycur,xp,dddx,xxp
  real(kind=8), dimension(:,:,:,:,:), allocatable :: dx,dxx
  real(kind=8), dimension(:,:,:,:), allocatable :: deps
  real(kind=8), dimension(:,:,:), allocatable :: re,ccorr,ooneosqrteps
  real(kind=8) :: pi,switch,hxd,hyd,hzd
  integer :: i1,i2,i3,isp,i,n01d,n02d,n03d,i1_max,i2_max,i3_max
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8) :: PB_charge,max_val,fact

  pi = 4.d0*datan(1.d0)   
  n01d=n01*2
  n02d=n02*2
  n03d=n03*2
  hxd=hx/2.d0
  hyd=hy/2.d0
  hzd=hz/2.d0
  !hxd=acell/real(n01d-1,kind=8)
  !hyd=acell/real(n02d-1,kind=8)
  !hzd=acell/real(n03d-1,kind=8)

  ddx=f_malloc([n01,n02,n03,nspden],id='ddx')
  dddx=f_malloc([n01d,n02d,n03d,nspden],id='ddx')
  xp=f_malloc([n01,n02,n03,nspden],id='xp')
  ycur=f_malloc([n01,n02,n03,nspden],id='ycur')
  xx=f_malloc([n01d,n02d,n03d,nspden],id='xx')
  xxp=f_malloc([n01d,n02d,n03d,nspden],id='xxp')
  yy=f_malloc([n01d,n02d,n03d,nspden],id='yy')
  dx=f_malloc([n01,n02,n03,nspden,3],id='dx')
  dxx=f_malloc([n01d,n02d,n03d,nspden,3],id='dxx')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  ccorr=f_malloc([n01d,n02d,n03d],id='eeps')
  ooneosqrteps=f_malloc([n01d,n02d,n03d],id='eeps')
  re=f_malloc([n01,n02,n03],id='re')
  !write(*,*)geocode

      isp=1
      do i3=1,n03
       do i2=1,n02
        do i1=1,n01
          xp(i1,i2,i3,isp)=x(i1,i2,i3,isp)/oneosqrteps(i1,i2,i3)
        end do
       end do
      end do

!  call fssnord3DmatNabla(geocode,n01,n02,n03,nspden,hx,hy,hz,xp,dx,nord,acell)
!  call fssnord3DmatDiv(geocode,n01,n02,n03,nspden,hx,hy,hz,dx,y,nord,acell)
  call fssnordEpsilonDerivative(n01,n02,n03,nspden,hx,hy,hz,xp,ddx,y,nord,acell)

      isp=1
      do i3=1,n03
       do i2=1,n02
        do i1=1,n01
          y(i1,i2,i3,isp)=-y(i1,i2,i3,isp)/oneosqrteps(i1,i2,i3)*0.25d0/pi+x(i1,i2,i3,isp)*corr(i1,i2,i3)
        end do
       end do
      end do


!   y(:,:,:,:)=-y(:,:,:,:)/(4.d0*pi)

!------------------------------------------------------------------------------------
   !Apply linear PB operator.
   if (lin_PB) then
    isp=1
    do i3=1,n03
     do i2=1,n02
      do i1=1,n01
       y(i1,i2,i3,isp) = y(i1,i2,i3,isp) - ((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(x(i1,i2,i3,isp))
      end do
     end do
    end do
   end if

   call f_free(ddx)
   call f_free(dddx)
   call f_free(ycur)
   call f_free(xx)
   call f_free(xxp)
   call f_free(xp)
   call f_free(yy)
   call f_free(dx)
   call f_free(dxx)
   call f_free(deps)
   call f_free(ccorr)
   call f_free(ooneosqrteps)
   call f_free(re)

end subroutine ApplyLaplace_corr

subroutine ApplyLaplace_ISF_corr(geocode,n01,n02,n03,nspden,hx,hy,hz,x,y,acell,eps,corr,oneosqrteps,nord,lin_PB,multp)
  use yaml_output
  use dynamic_memory
  implicit none
  character(len=2), intent(in) :: geocode
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: x
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: corr
  real(kind=8), dimension(n01,n02,n03), intent(in) :: oneosqrteps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: y
  logical, intent(in) :: lin_PB

  ! Local variables.
  real(kind=8), dimension(:,:,:,:), allocatable :: ddx,xx,yy,ycur,xp,dddx,xxp
  real(kind=8), dimension(:,:,:,:,:), allocatable :: dx,dxx
  real(kind=8), dimension(:,:,:,:), allocatable :: deps
  real(kind=8), dimension(:,:,:), allocatable :: re,ccorr,ooneosqrteps
  real(kind=8) :: pi,switch,hxd,hyd,hzd
  integer :: i1,i2,i3,isp,i,n01d,n02d,n03d,i1_max,i2_max,i3_max
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8) :: PB_charge,max_val,fact

  pi = 4.d0*datan(1.d0)   
  n01d=n01*2
  n02d=n02*2
  n03d=n03*2
  hxd=hx/2.d0
  hyd=hy/2.d0
  hzd=hz/2.d0
  !hxd=acell/real(n01d-1,kind=8)
  !hyd=acell/real(n02d-1,kind=8)
  !hzd=acell/real(n03d-1,kind=8)

  ddx=f_malloc([n01,n02,n03,nspden],id='ddx')
  dddx=f_malloc([n01d,n02d,n03d,nspden],id='ddx')
  xp=f_malloc([n01,n02,n03,nspden],id='xp')
  ycur=f_malloc([n01,n02,n03,nspden],id='ycur')
  xx=f_malloc([n01d,n02d,n03d,nspden],id='xx')
  xxp=f_malloc([n01d,n02d,n03d,nspden],id='xxp')
  yy=f_malloc([n01d,n02d,n03d,nspden],id='yy')
  dx=f_malloc([n01,n02,n03,nspden,3],id='dx')
  dxx=f_malloc([n01d,n02d,n03d,nspden,3],id='dxx')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  ccorr=f_malloc([n01d,n02d,n03d],id='eeps')
  ooneosqrteps=f_malloc([n01d,n02d,n03d],id='eeps')
  re=f_malloc([n01,n02,n03],id='re')
  !write(*,*)geocode

      isp=1
      do i3=1,n03
       do i2=1,n02
        do i1=1,n01
          xp(i1,i2,i3,isp)=x(i1,i2,i3,isp)/oneosqrteps(i1,i2,i3)
        end do
       end do
      end do

!  call fssnord3DmatNabla(geocode,n01,n02,n03,nspden,hx,hy,hz,xp,dx,nord,acell)
!  call fssnord3DmatDiv(geocode,n01,n02,n03,nspden,hx,hy,hz,dx,y,nord,acell)
  call fssnordEpsilonDerivative(n01,n02,n03,nspden,hx,hy,hz,xp,ddx,y,nord,acell)

      isp=1
      do i3=1,n03
       do i2=1,n02
        do i1=1,n01
          y(i1,i2,i3,isp)=-y(i1,i2,i3,isp)/oneosqrteps(i1,i2,i3)*0.25d0/pi+x(i1,i2,i3,isp)*corr(i1,i2,i3)
        end do
       end do
      end do


!   y(:,:,:,:)=-y(:,:,:,:)/(4.d0*pi)
   ycur(:,:,:,:)=y(:,:,:,:)
!
!!------------------------------------------------------------------------------------
!   !Apply the improved ISF Generalized Poisson operator.
!
   call back_trans_ISF_3D(geocode,16,n01d,n02d,n03d,x,xx)
   call back_trans_ISF_3D(geocode,16,n01d,n02d,n03d,oneosqrteps,ooneosqrteps)
   call back_trans_ISF_3D(geocode,16,n01d,n02d,n03d,corr,ccorr)

      isp=1
      do i3=1,n03d
       do i2=1,n02d
        do i1=1,n01d
          xxp(i1,i2,i3,isp)=xx(i1,i2,i3,isp)/ooneosqrteps(i1,i2,i3)
        end do
       end do
      end do
!  call fssnord3DmatNabla(geocode,n01,n02,n03,nspden,hx,hy,hz,xp,dx,nord,acell)
!  call fssnord3DmatDiv(geocode,n01,n02,n03,nspden,hx,hy,hz,dx,y,nord,acell)
  call fssnordEpsilonDerivative(n01d,n02d,n03d,nspden,hxd,hyd,hzd,xxp,dddx,yy,nord,acell)
      
      isp=1
      do i3=1,n03d
       do i2=1,n02d
        do i1=1,n01d
          yy(i1,i2,i3,isp)=-yy(i1,i2,i3,isp)/ooneosqrteps(i1,i2,i3)*0.25d0/pi+xx(i1,i2,i3,isp)*ccorr(i1,i2,i3)
        end do
       end do
      end do


   call for_trans_ISF_3D(geocode,16,n01d,n02d,n03d,yy,y)
!
!   y(:,:,:,:)=-y(:,:,:,:)/(4.d0*pi)
!   yy(:,:,:,:)=-yy(:,:,:,:)/(4.d0*pi)

!---------------------------------------------------------
!Checking difference between coarse/finer grid product.

  max_val = 0.d0
  i1_max = 1
  i2_max = 1
  i3_max = 1
  isp=1
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           !re(i1,i2,i3) = ycur(i1,i2,i3,1) - y(i1,i2,i3,1)
           re(i1,i2,i3) = ycur(i1,i2,i3,1) - yy(2*i1-1,2*i2-1,2*i3-1,1)
           fact=abs(re(i1,i2,i3))
           if (max_val < fact) then
              max_val = fact
              i1_max = i1
              i2_max = i2
              i3_max = i3
           end if
        end do
     end do
  end do
  call yaml_map('Diff ISF operator c-f ISF grid',max_val,fmt='(1pe22.15)')
  write(*,*)i1_max,i2_max,i3_max
!---------------------------------------------------------

!------------------------------------------------------------------------------------
   !Apply linear PB operator.
   if (lin_PB) then
    isp=1
    do i3=1,n03
     do i2=1,n02
      do i1=1,n01
       y(i1,i2,i3,isp) = y(i1,i2,i3,isp) - ((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(x(i1,i2,i3,isp))
      end do
     end do
    end do
   end if

   call f_free(ddx)
   call f_free(dddx)
   call f_free(ycur)
   call f_free(xx)
   call f_free(xxp)
   call f_free(xp)
   call f_free(yy)
   call f_free(dx)
   call f_free(dxx)
   call f_free(deps)
   call f_free(ccorr)
   call f_free(ooneosqrteps)
   call f_free(re)

end subroutine ApplyLaplace_ISF_corr

subroutine ApplyLaplace2(geocode,n01,n02,n03,nspden,hx,hy,hz,x,y,acell,eps,dlogeps,nord,lin_PB,multp)
  use dynamic_memory
  implicit none
  character(len=2), intent(in) :: geocode
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,multp
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: x
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(3,n01,n02,n03), intent(in) :: dlogeps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: y
  logical, intent(in) :: lin_PB

  ! Local variables.
  real(kind=8), dimension(:,:,:,:), allocatable :: ddx
  real(kind=8), dimension(:,:,:,:,:), allocatable :: dx
  real(kind=8), dimension(:,:,:,:), allocatable :: deps
  real(kind=8) :: pi,switch,prod
  integer :: i1,i2,i3,isp,i
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8) :: PB_charge

  pi = 4.d0*datan(1.d0)   

  ddx=f_malloc([n01,n02,n03,nspden],id='ddx')
  dx=f_malloc([n01,n02,n03,nspden,3],id='dx')
  deps=f_malloc([n01,n02,n03,3],id='deps')
  !write(*,*)geocode
  call fssnord3DmatNabla(geocode,n01,n02,n03,nspden,hx,hy,hz,x,dx,nord,acell)
  call fssnord3DmatDiv(geocode,n01,n02,n03,nspden,hx,hy,hz,dx,y,nord,acell)

      isp=1
      do i3=1,n03
       do i2=1,n02
        do i1=1,n01
         prod=0.d0
         do i=1,3
          prod = prod + dlogeps(i,i1,i2,i3)*dx(i1,i2,i3,isp,i)
         end do
          y(i1,i2,i3,isp) = -0.25d0*eps(i1,i2,i3)*(prod+y(i1,i2,i3,isp))/pi
        end do
       end do
      end do


!   y(:,:,:,:)=-y(:,:,:,:)/(4.d0*pi)

   if (lin_PB) then
    isp=1
    do i3=1,n03
     do i2=1,n02
      do i1=1,n01
       y(i1,i2,i3,isp) = y(i1,i2,i3,isp) - ((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(x(i1,i2,i3,isp))
!       y(i1,i2,i3,isp) = y(i1,i2,i3,isp) + ((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dsinh(multp*x(i1,i2,i3,isp))
!       y(i1,i2,i3,isp) = y(i1,i2,i3,isp) + ((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dtanh(multp*x(i1,i2,i3,isp))
!       y(i1,i2,i3,isp) = y(i1,i2,i3,isp) + ((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*x(i1,i2,i3,isp)
!       y(i1,i2,i3,isp) = y(i1,i2,i3,isp) + ((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*(x(i1,i2,i3,isp)**2)
      end do
     end do
    end do
   end if

   call f_free(deps)
   call f_free(ddx)
   call f_free(dx)

end subroutine ApplyLaplace2

subroutine Polarization_charge(n01,n02,n03,nspden,hx,hy,hz,x,y,acell,eps,nord)

  use dynamic_memory

  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  real(kind=8), intent(in) :: hx,hy,hz
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell
  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: x
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: y
  real(kind=8), dimension(:,:,:,:), allocatable :: ddx
  real(kind=8), dimension(:,:,:,:,:), allocatable :: dx
  real(kind=8), dimension(:,:,:,:), allocatable :: deps
  real(kind=8) :: pi
  integer :: i1,i2,i3,isp,i

  pi = 4.d0*datan(1.d0)

  open(unit=23,file='Pol_charge.dat',status='unknown')
  open(unit=24,file='Pol_charge_line.dat',status='unknown')

  ddx=f_malloc([n01,n02,n03,nspden],id='ddx')
  dx=f_malloc([n01,n02,n03,nspden,3],id='dx')
  deps=f_malloc([n01,n02,n03,3],id='deps')

  call fssnord3DmatNabla("F ",n01,n02,n03,nspden,hx,hy,hz,x,dx,nord,acell)

      isp=1
      do i3=1,n03
       do i2=1,n02
        do i1=1,n01
         do i=1,3
          dx(i1,i2,i3,isp,i)=(eps(i1,i2,i3)-1)*dx(i1,i2,i3,isp,i)/(4.d0*pi)
         end do
        end do
       end do
      end do

   call fssnord3DmatDiv("F ",n01,n02,n03,nspden,hx,hy,hz,dx,y,nord,acell)

     i3=1!n03/2
     do i2=1,n02
      do i1=1,n01
       write(23,'(2(1x,I4),2(1x,e14.7))')i1,i2,y(i1,i2,i3,1),y(i1,i2,n03/2,1)
      end do
      write(23,*)
     end do

     do i1=1,n01
      write(24,'(1x,I8,2(1x,e22.15))') i1,y(i1,n02/2,n03/2,1)!,y(n01/2,i1,n03/2,1)
     end do

  close(unit=23)
  close(unit=24)

   call f_free(deps)
   call f_free(ddx)
   call f_free(dx)

end subroutine Polarization_charge

subroutine fssnord3DmatNabla(geocode,n01,n02,n03,nspden,hx,hy,hz,u,du,nord,acell)
      implicit none

!c..this routine computes 'nord' order accurate first derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.

!c..input:
!c..ngrid       = number of points in the grid, 
!c..u(ngrid)    = function values at the grid points

!c..output:
!c..du(ngrid)   = first derivative values at the grid points

!c..declare the pass
      character(len=2), intent(in) :: geocode
      integer, intent(in) :: n01,n02,n03,nspden,nord
      real(kind=8), intent(in) :: hx,hy,hz
      real(kind=8), intent(in) :: acell
      real(kind=8), dimension(n01,n02,n03,nspden) :: u
      real(kind=8), dimension(n01,n02,n03,nspden,3) :: du

!c..local variables
      integer :: n,m,n_cell
      integer :: i,j,ib,i1,i2,i3,isp,i1_max,i2_max,ii
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D,c1DF
      real(kind=8) :: max_diff,fact
      logical :: perx,pery,perz

      n = nord+1
      m = nord/2
      n_cell = max(n01,n02,n03)

      !buffers associated to the geocode
      !conditions for periodicity in the three directions
      perx=(geocode /= 'F')
      pery=(geocode == 'P')
      perz=(geocode /= 'F')

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
        c1DF(i,j)=0.d0
       end do
      end do

       include 'FiniteDiffCorff.inc'

      isp=1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

             du(i1,i2,i3,isp,1) = 0.0d0

             if (i1.le.m) then
              if (perx) then
               do j=-m,m
                ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                du(i1,i2,i3,isp,1) = du(i1,i2,i3,isp,1) + c1D(j,0)*u(ii,i2,i3,isp)!/hx
               end do
              else
               do j=-m,m
                du(i1,i2,i3,isp,1) = du(i1,i2,i3,isp,1) + c1D(j,i1-m-1)*u(j+m+1,i2,i3,isp)!/hx
               end do
              end if
             else if (i1.gt.n01-m) then
              if (perx) then
               do j=-m,m
                ii=modulo(i1 + j - 1, n01 ) + 1
                du(i1,i2,i3,isp,1) = du(i1,i2,i3,isp,1) + c1D(j,0)*u(ii,i2,i3,isp)!/hx
               end do
              else
               do j=-m,m
                du(i1,i2,i3,isp,1) = du(i1,i2,i3,isp,1) + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3,isp)!/hx
               end do
              end if
             else
              do j=-m,m
               du(i1,i2,i3,isp,1) = du(i1,i2,i3,isp,1) + c1D(j,0)*u(i1 + j,i2,i3,isp)!/hx
              end do
             end if
             du(i1,i2,i3,isp,1)=du(i1,i2,i3,isp,1)/hx

             du(i1,i2,i3,isp,2) = 0.0d0

             if (i2.le.m) then
              if (pery) then
               do j=-m,m
                ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                du(i1,i2,i3,isp,2) = du(i1,i2,i3,isp,2) + c1D(j,0)*u(i1,ii,i3,isp)!/hy
               end do
              else
               do j=-m,m
                du(i1,i2,i3,isp,2) = du(i1,i2,i3,isp,2) + c1D(j,i2-m-1)*u(i1,j+m+1,i3,isp)!/hy
               end do 
              end if
             else if (i2.gt.n02-m) then
              if (pery) then
               do j=-m,m
                ii=modulo(i2 + j - 1, n02 ) + 1
                du(i1,i2,i3,isp,2) = du(i1,i2,i3,isp,2) + c1D(j,0)*u(i1,ii,i3,isp)!/hy
               end do
              else
               do j=-m,m
                du(i1,i2,i3,isp,2) = du(i1,i2,i3,isp,2) + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3,isp)!/hy
               end do
              end if
             else
              do j=-m,m
               du(i1,i2,i3,isp,2) = du(i1,i2,i3,isp,2) + c1D(j,0)*u(i1,i2 + j,i3,isp)!/hy
              end do
             end if
              du(i1,i2,i3,isp,2)=du(i1,i2,i3,isp,2)/hy

             du(i1,i2,i3,isp,3) = 0.0d0

             if (i3.le.m) then
              if (perz) then
               do j=-m,m
                ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                du(i1,i2,i3,isp,3) = du(i1,i2,i3,isp,3) + c1D(j,0)*u(i1,i2,ii,isp)!/hx
               end do
              else
               do j=-m,m
                du(i1,i2,i3,isp,3) = du(i1,i2,i3,isp,3) + c1D(j,i3-m-1)*u(i1,i2,j+m+1,isp)!/hz
               end do
              end if
             else if (i3.gt.n03-m) then
              if (perz) then
               do j=-m,m
                ii=modulo(i3 + j - 1, n03 ) + 1
                du(i1,i2,i3,isp,3) = du(i1,i2,i3,isp,3) + c1D(j,0)*u(i1,i2,ii,isp)!/hx
               end do
              else
               do j=-m,m
                du(i1,i2,i3,isp,3) = du(i1,i2,i3,isp,3) + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m,isp)!/hz
               end do
              end if
             else
              do j=-m,m
               du(i1,i2,i3,isp,3) = du(i1,i2,i3,isp,3) + c1D(j,0)*u(i1,i2,i3 + j,isp)!/hz
              end do
             end if
              du(i1,i2,i3,isp,3)=du(i1,i2,i3,isp,3)/hz

            end do
         end do
      end do

end subroutine fssnord3DmatNabla

!> Like fssnord3DmatNabla but corrected such that the index goes at the beginning
!! Multiplies also times (nabla epsilon)/(4pi*epsilon)= nabla (log(epsilon))/(4*pi)
subroutine fssnord3DmatNabla_LG2(n01,n02,n03,hx,hy,hz,u,nord,acell,eta,dlogeps,eps,rhopol,rhores2)
  !use module_defs, only: pi_param
  implicit none

  !c..this routine computes 'nord' order accurate first derivatives 
  !c..on a equally spaced grid with coefficients from 'Matematica' program.

  !c..input:
  !c..ngrid       = number of points in the grid, 
  !c..u(ngrid)    = function values at the grid points

  !c..output:
  !c..du(ngrid)   = first derivative values at the grid points

  !c..declare the pass

  integer, intent(in) :: n01,n02,n03,nord
  real(kind=8), intent(in) :: hx,hy,hz
  real(kind=8), intent(in) :: acell,eta
  real(kind=8), dimension(n01,n02,n03), intent(in) :: u
  real(kind=8), dimension(3,n01,n02,n03), intent(in) :: dlogeps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(n01,n02,n03), intent(inout) :: rhopol
  real(kind=8), intent(out) :: rhores2

  !c..local variables
  integer :: n,m,n_cell
  integer :: i,j,ib,i1,i2,i3,isp,i1_max,i2_max
  !real(kind=8), parameter :: oneo4pi=0.25d0/pi_param
  real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D,c1DF
  real(kind=8) :: max_diff,fact,dx,dy,dz,res,rho
  real(kind=8) :: oneo4pi

  oneo4pi=1.0d0/(16.d0*atan(1.d0))

  n = nord+1
  m = nord/2
  n_cell = max(n01,n02,n03)

  ! Beware that n_cell has to be > than n.
  if (n_cell.lt.n) then
     write(*,*)'ngrid in has to be setted > than n=nord + 1'
     stop
  end if

  ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
  !Only nord=2,4,6,8,16
  if (all(nord /=[2,4,6,8,16])) then
     write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
     stop
  end if

  do i=-m,m
     do j=-m,m
        c1D(i,j)=0.d0
        c1DF(i,j)=0.d0
     end do
  end do

  include 'FiniteDiffCorff.inc'

  rhores2=0.d0
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01

           dx=0.d0

           if (i1.le.m) then
              do j=-m,m
                 dx = dx + c1D(j,i1-m-1)*u(j+m+1,i2,i3)
              end do
           else if (i1.gt.n01-m) then
              do j=-m,m
                 dx = dx + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3)
              end do
           else
              do j=-m,m
                 dx = dx + c1D(j,0)*u(i1 + j,i2,i3)
              end do
           end if
           dx=dx/hx

           dy = 0.0d0
           if (i2.le.m) then
              do j=-m,m
                 dy = dy + c1D(j,i2-m-1)*u(i1,j+m+1,i3)
              end do
           else if (i2.gt.n02-m) then
              do j=-m,m
                 dy = dy + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3)
              end do
           else
              do j=-m,m
                 dy = dy + c1D(j,0)*u(i1,i2 + j,i3)
              end do
           end if
           dy=dy/hy

           dz = 0.0d0
           if (i3.le.m) then
              do j=-m,m
                 dz = dz + c1D(j,i3-m-1)*u(i1,i2,j+m+1)
              end do
           else if (i3.gt.n03-m) then
              do j=-m,m
                 dz = dz + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m)
              end do
           else
              do j=-m,m
                 dz = dz + c1D(j,0)*u(i1,i2,i3 + j)
              end do
           end if
           dz=dz/hz

           !retrieve the previous treatment
           res = dlogeps(1,i1,i2,i3)*dx + &
                dlogeps(2,i1,i2,i3)*dy + dlogeps(3,i1,i2,i3)*dz
           res = res*oneo4pi
           rho=rhopol(i1,i2,i3)
           res=res-rho
           res=eta*res
           rhores2=rhores2+res*res*eps(i1,i2,i3)**2
           rhopol(i1,i2,i3)=res+rho

        end do
     end do
  end do

end subroutine fssnord3DmatNabla_LG2


subroutine fssnord3DmatNabla3var(n01,n02,n03,nspden,hx,hy,hz,u,du,nord,acell)
      implicit none

!c..this routine computes 'nord' order accurate first derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.

!c..input:
!c..ngrid       = number of points in the grid, 
!c..u(ngrid)    = function values at the grid points

!c..output:
!c..du(ngrid)   = first derivative values at the grid points

!c..declare the pass
      integer, intent(in) :: n01,n02,n03,nspden,nord
      real(kind=8), intent(in) :: hx,hy,hz
      real(kind=8), intent(in) :: acell
      real(kind=8), dimension(n01,n02,n03) :: u
      real(kind=8), dimension(n01,n02,n03,3) :: du

!c..local variables
      integer :: n,m,n_cell
      integer :: i,j,ib,i1,i2,i3
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D

      n = nord+1
      m = nord/2
      n_cell = max(n01,n02,n03)

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
       end do
      end do

       include 'FiniteDiffCorff.inc'

      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

             du(i1,i2,i3,1) = 0.0d0

             if (i1.le.m) then
              do j=-m,m
               du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,i1-m-1)*u(j+m+1,i2,i3)/hx
              end do
             else if (i1.gt.n01-m) then
              do j=-m,m
               du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3)/hx
              end do
             else
              do j=-m,m
               du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,0)*u(i1 + j,i2,i3)/hx
              end do
             end if

             du(i1,i2,i3,2) = 0.0d0

             if (i2.le.m) then
             do j=-m,m
              du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,i2-m-1)*u(i1,j+m+1,i3)/hy
             end do 
             else if (i2.gt.n02-m) then
              do j=-m,m
               du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3)/hy
              end do
             else
              do j=-m,m
               du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,0)*u(i1,i2 + j,i3)/hy
              end do
             end if

             du(i1,i2,i3,3) = 0.0d0

             if (i3.le.m) then
             do j=-m,m
              du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,i3-m-1)*u(i1,i2,j+m+1)/hz
             end do
             else if (i3.gt.n03-m) then
              do j=-m,m
               du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m)/hz
              end do
             else
              do j=-m,m
               du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,0)*u(i1,i2,i3 + j)/hz
              end do
             end if

            end do
         end do
      end do

end subroutine fssnord3DmatNabla3var

subroutine fssnord3DmatNabla3varde2(n01,n02,n03,nspden,hx,hy,hz,u,du,du2,nord,acell)
      implicit none

!c..this routine computes 'nord' order accurate first derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.

!c..input:
!c..ngrid       = number of points in the grid, 
!c..u(ngrid)    = function values at the grid points

!c..output:
!c..du(ngrid)   = first derivative values at the grid points

!c..declare the pass
      integer, intent(in) :: n01,n02,n03,nspden,nord
      real(kind=8), intent(in) :: hx,hy,hz
      real(kind=8), intent(in) :: acell
      real(kind=8), dimension(n01,n02,n03) :: u
      real(kind=8), dimension(n01,n02,n03,3) :: du
      real(kind=8), dimension(n01,n02,n03) :: du2

!c..local variables
      integer :: n,m,n_cell
      integer :: i,j,ib,i1,i2,i3
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D

      n = nord+1
      m = nord/2
      n_cell = max(n01,n02,n03)

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
       end do
      end do

       include 'FiniteDiffCorff.inc'

      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

             du(i1,i2,i3,1) = 0.0d0
             du2(i1,i2,i3) = 0.0d0

             if (i1.le.m) then
              do j=-m,m
               du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,i1-m-1)*u(j+m+1,i2,i3)/hx
              end do
             else if (i1.gt.n01-m) then
              do j=-m,m
               du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3)/hx
              end do
             else
              do j=-m,m
               du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,0)*u(i1 + j,i2,i3)/hx
              end do
             end if

             du2(i1,i2,i3) = du(i1,i2,i3,1)*du(i1,i2,i3,1)
             du(i1,i2,i3,2) = 0.0d0

             if (i2.le.m) then
             do j=-m,m
              du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,i2-m-1)*u(i1,j+m+1,i3)/hy
             end do 
             else if (i2.gt.n02-m) then
              do j=-m,m
               du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3)/hy
              end do
             else
              do j=-m,m
               du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,0)*u(i1,i2 + j,i3)/hy
              end do
             end if

             du2(i1,i2,i3) = du2(i1,i2,i3) + du(i1,i2,i3,2)*du(i1,i2,i3,2)

             du(i1,i2,i3,3) = 0.0d0

             if (i3.le.m) then
             do j=-m,m
              du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,i3-m-1)*u(i1,i2,j+m+1)/hz
             end do
             else if (i3.gt.n03-m) then
              do j=-m,m
               du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m)/hz
              end do
             else
              do j=-m,m
               du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,0)*u(i1,i2,i3 + j)/hz
              end do
             end if

             du2(i1,i2,i3) = du2(i1,i2,i3) + du(i1,i2,i3,3)*du(i1,i2,i3,3)

            end do
         end do
      end do

end subroutine fssnord3DmatNabla3varde2

subroutine fssnordEpsilonDerivative(n01,n02,n03,nspden,hx,hy,hz,u,du2,ddu,nord,acell)
      implicit none

!c..this routine computes 'nord' order accurate first and second derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.


      integer, intent(in) :: n01,n02,n03,nspden,nord
      real(kind=8), intent(in) :: hx,hy,hz
      real(kind=8), intent(in) :: acell
      real(kind=8), dimension(n01,n02,n03), intent(in) :: u
      real(kind=8), dimension(n01,n02,n03), intent(out) :: du2,ddu

      integer :: n,m,n_cell
      integer :: i,j,ib,i1,i2,i3,i1_max,i2_max
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D,c2D
      real(kind=8) :: d,dd

      n = nord+1
      m = nord/2
      n_cell = min(n01,n02,n03)

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first and second derivative coefficients from 'Matematica'.
      !Only nord=2,4,6,8,16 admitted.

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
        c2D(i,j)=0.d0
       end do
      end do


       include 'FiniteDiffCorff.inc'
       include 'FiniteDiffCorff_2der.inc'

      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

             du2(i1,i2,i3) = 0.0d0
             ddu(i1,i2,i3) = 0.0d0
             d=0.d0
             dd=0.d0

             if (i1.le.m) then
              do j=-m,m
               d = d + c1D(j,i1-m-1)*u(j+m+1,i2,i3)
               dd = dd + c2D(j,i1-m-1)*u(j+m+1,i2,i3)
              end do
             else if (i1.gt.n01-m) then
              do j=-m,m
               d = d + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3)
               dd = dd + c2D(j,i1-n01+m)*u(n01 + j - m,i2,i3)
              end do
             else
              do j=-m,m
               d = d + c1D(j,0)*u(i1 + j,i2,i3)
               dd = dd + c2D(j,0)*u(i1 + j,i2,i3)
              end do
             end if

             du2(i1,i2,i3) = d*d/hx/hx
             ddu(i1,i2,i3) = dd/hx/hx

             d=0.d0
             dd=0.d0

             if (i2.le.m) then
             do j=-m,m
              d = d + c1D(j,i2-m-1)*u(i1,j+m+1,i3)
              dd = dd + c2D(j,i2-m-1)*u(i1,j+m+1,i3)
             end do 
             else if (i2.gt.n02-m) then
              do j=-m,m
               d = d + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3)
               dd = dd + c2D(j,i2-n02+m)*u(i1,n02 + j - m,i3)
              end do
             else
              do j=-m,m
               d = d + c1D(j,0)*u(i1,i2 + j,i3)
               dd = dd + c2D(j,0)*u(i1,i2 + j,i3)
              end do
             end if

             du2(i1,i2,i3) = du2(i1,i2,i3) + d*d/hy/hy
             ddu(i1,i2,i3) = ddu(i1,i2,i3) + dd/hy/hy

             d=0.d0
             dd=0.d0

             if (i3.le.m) then
             do j=-m,m
              d = d + c1D(j,i3-m-1)*u(i1,i2,j+m+1)
              dd = dd + c2D(j,i3-m-1)*u(i1,i2,j+m+1)
             end do
             else if (i3.gt.n03-m) then
              do j=-m,m
               d = d + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m)
               dd = dd + c2D(j,i3-n03+m)*u(i1,i2,n03 + j - m)
              end do
             else
              do j=-m,m
               d = d + c1D(j,0)*u(i1,i2,i3 + j)
               dd = dd + c2D(j,0)*u(i1,i2,i3 + j)
              end do
             end if

             du2(i1,i2,i3) = du2(i1,i2,i3) + d*d/hz/hz
             ddu(i1,i2,i3) = ddu(i1,i2,i3) + dd/hz/hz

            end do
         end do
      end do

end subroutine fssnordEpsilonDerivative

subroutine fssnord3DmatDiv(geocode,n01,n02,n03,nspden,hx,hy,hz,u,du,nord,acell)
      implicit none

!c..this routine computes 'nord' order accurate first derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.

!c..input:
!c..ngrid       = number of points in the grid, 
!c..u(ngrid)    = function values at the grid points

!c..output:
!c..du(ngrid)   = first derivative values at the grid points

!c..declare the pass
      character(len=2), intent(in) :: geocode
      integer, intent(in) :: n01,n02,n03,nspden,nord
      real(kind=8), intent(in) :: hx,hy,hz
      real(kind=8), intent(in) :: acell
      real(kind=8), dimension(n01,n02,n03,nspden,3) :: u
      real(kind=8), dimension(n01,n02,n03,nspden) :: du

!c..local variables
      integer :: n,m,n_cell
      integer :: i,j,ib,i1,i2,i3,isp,ii
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D
      real(kind=8) :: d1,d2,d3
      real(kind=8), parameter :: zero = 0.d0! 1.0d-11
      logical :: perx,pery,perz

      n = nord+1
      m = nord/2
      n_cell = max(n01,n02,n03)

      !buffers associated to the geocode
      !conditions for periodicity in the three directions
      perx=(geocode /= 'F')
      pery=(geocode == 'P')
      perz=(geocode /= 'F')

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
       end do
      end do

       include 'FiniteDiffCorff.inc'

      isp=1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

             du(i1,i2,i3,isp) = 0.0d0

             d1 = 0.d0
             if (i1.le.m) then
              if (perx) then
               do j=-m,m
                ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                d1 = d1 + c1D(j,0)*u(ii,i2,i3,isp,1)!/hx
               end do
              else
               do j=-m,m
                d1 = d1 + c1D(j,i1-m-1)*u(j+m+1,i2,i3,isp,1)!/hx
               end do
              end if
             else if (i1.gt.n01-m) then
              if (perx) then
               do j=-m,m
                ii=modulo(i1 + j - 1, n01 ) + 1
                d1 = d1 + c1D(j,0)*u(ii,i2,i3,isp,1)!/hx
               end do
              else
               do j=-m,m
                d1 = d1 + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3,isp,1)!/hx
               end do
              end if
             else
              do j=-m,m
               d1 = d1 + c1D(j,0)*u(i1 + j,i2,i3,isp,1)!/hx
              end do
             end if
              d1=d1/hx

             d2 = 0.d0
             if (i2.le.m) then
              if (pery) then
               do j=-m,m
                ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                d2 = d2 + c1D(j,0)*u(i1,ii,i3,isp,2)!/hy
               end do
              else
               do j=-m,m
                d2 = d2 + c1D(j,i2-m-1)*u(i1,j+m+1,i3,isp,2)!/hy
               end do 
              end if
             else if (i2.gt.n02-m) then
              if (pery) then
               do j=-m,m
                ii=modulo(i2 + j - 1, n02 ) + 1
                d2 = d2 + c1D(j,0)*u(i1,ii,i3,isp,2)!/hy
               end do
              else
               do j=-m,m
                d2 = d2 + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3,isp,2)!/hy
               end do
              end if
             else
              do j=-m,m
               d2 = d2 + c1D(j,0)*u(i1,i2 + j,i3,isp,2)!/hy
              end do
             end if
              d2=d2/hy

             d3 = 0.d0
             if (i3.le.m) then
              if (perz) then
               do j=-m,m
                ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                d3 = d3 + c1D(j,0)*u(i1,i2,ii,isp,3)!/hz
               end do
              else
               do j=-m,m
                d3 = d3 + c1D(j,i3-m-1)*u(i1,i2,j+m+1,isp,3)!/hz
               end do
              end if
             else if (i3.gt.n03-m) then
              if (perz) then
               do j=-m,m
                ii=modulo(i3 + j - 1, n03 ) + 1
                d3 = d3 + c1D(j,0)*u(i1,i2,ii,isp,3)!/hz
               end do
              else
               do j=-m,m
                d3 = d3 + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m,isp,3)!/hz
               end do
              end if
             else
              do j=-m,m
               d3 = d3 + c1D(j,0)*u(i1,i2,i3 + j,isp,3)!/hz
              end do
             end if
              d3=d3/hz

             du(i1,i2,i3,isp) = d1+d2+d3

            end do
         end do
      end do

end subroutine fssnord3DmatDiv

subroutine fssnord3DmatDiv3var(n01,n02,n03,nspden,hx,hy,hz,u,du,nord,acell)
      implicit none

!c..this routine computes 'nord' order accurate first derivatives 
!c..on a equally spaced grid with coefficients from 'Matematica' program.

!c..input:
!c..ngrid       = number of points in the grid, 
!c..u(ngrid)    = function values at the grid points

!c..output:
!c..du(ngrid)   = first derivative values at the grid points

!c..declare the pass
      integer, intent(in) :: n01,n02,n03,nspden,nord
      real(kind=8), intent(in) :: hx,hy,hz
      real(kind=8), intent(in) :: acell
      real(kind=8), dimension(n01,n02,n03,3) :: u
      real(kind=8), dimension(n01,n02,n03) :: du

!c..local variables
      integer :: n,m,n_cell
      integer :: i,j,ib,i1,i2,i3,isp
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D
      real(kind=8) :: d1,d2,d3
      real(kind=8), parameter :: zero = 0.d0! 1.0d-11

      n = nord+1
      m = nord/2
      n_cell = max(n01,n02,n03)

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
       write(*,*)'ngrid in has to be setted > than n=nord + 1'
       stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
       !O.K.
      case default
       write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
       stop
      end select

      do i=-m,m
       do j=-m,m
        c1D(i,j)=0.d0
       end do
      end do

       include 'FiniteDiffCorff.inc'

      isp=1
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

             du(i1,i2,i3) = 0.0d0

             d1 = 0.d0
             if (i1.le.m) then
              do j=-m,m
               d1 = d1 + c1D(j,i1-m-1)*u(j+m+1,i2,i3,1)!/hx
              end do
             else if (i1.gt.n01-m) then
              do j=-m,m
               d1 = d1 + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3,1)!/hx
              end do
             else
              do j=-m,m
               d1 = d1 + c1D(j,0)*u(i1 + j,i2,i3,1)!/hx
              end do
             end if
              d1=d1/hx

             d2 = 0.d0
             if (i2.le.m) then
             do j=-m,m
              d2 = d2 + c1D(j,i2-m-1)*u(i1,j+m+1,i3,2)!/hy
             end do 
             else if (i2.gt.n02-m) then
              do j=-m,m
               d2 = d2 + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3,2)!/hy
              end do
             else
              do j=-m,m
               d2 = d2 + c1D(j,0)*u(i1,i2 + j,i3,2)!/hy
              end do
             end if
              d2=d2/hy

             d3 = 0.d0
             if (i3.le.m) then
             do j=-m,m
              d3 = d3 + c1D(j,i3-m-1)*u(i1,i2,j+m+1,3)!/hz
             end do
             else if (i3.gt.n03-m) then
              do j=-m,m
               d3 = d3 + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m,3)!/hz
              end do
             else
              do j=-m,m
               d3 = d3 + c1D(j,0)*u(i1,i2,i3 + j,3)!/hz
              end do
             end if
              d3=d3/hz

             du(i1,i2,i3) = d1+d2+d3

            end do
         end do
      end do

end subroutine fssnord3DmatDiv3var

subroutine SetInitDensPot(n01,n02,n03,nspden,iproc,natreal,eps,dlogeps,sigmaeps,SetEps,erfL,erfR,&
     acell,a_gauss,a2,hx,hy,hz,Setrho,density,potential,geocode,offset,einit,multp,rxyzreal,lin_PB)
  use dynamic_memory
  use yaml_output
  use f_utils
  use box
  use numerics, only: safe_exp,oneofourpi

  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden,iproc,natreal
  integer, intent(in) :: Setrho
  real(kind=8), intent(in) :: acell,a_gauss,a2,hx,hy,hz,sigmaeps,erfL,erfR
  integer, intent(in) :: SetEps
  real(kind=8), dimension(n01,n02,n03), intent(in) :: eps
  real(kind=8), dimension(3,n01,n02,n03), intent(in) :: dlogeps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: density
  real(kind=8), dimension(n01,n02,n03), intent(out) :: potential
  character(len=2), intent(in) :: geocode
  real(kind=8), intent(out) :: offset,einit
  real(kind=8), intent(in) :: multp
  logical, intent(in) :: lin_PB
  real(kind=8), dimension(n01,n02,n03) :: potential1
  real(kind=8), dimension(3,natreal), intent(inout) :: rxyzreal
  real(kind=8), dimension(:,:,:,:), allocatable :: density1,density2
  real(kind=8), dimension(:,:,:), allocatable :: potential2
  integer :: i,i1,i2,i3,ifx,ify,ifz,unt,iat,jat,ii,nat,j,k,l,px,py,pz,imax
  real(kind=8) :: sigma,sigma1,sigma2,pi,sumd,sump,tt1,tt2,x0,r12
  real(kind=8) :: x1,x2,x3,r,r2,r1,r22,derf_tt1,derf_tt2,factor,factor1,factor2,value,d
  real(kind=8) :: length,denval,derf_tt,k1,k2,switch,dens
  real(kind=8) :: x,y,fx,fx2,fy,fy2,fz,fz2,a,ax,ay,az,bx,by,bz,tt,fx1,fy1,fz1,valuemax
  real(kind=8), dimension(3) :: r_v
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8) :: PB_charge
  !local variables
  logical, parameter :: old=.false.
  real(kind=8), parameter :: valuebc=1.d-15
  real(kind=8), dimension(3) :: rv,shift,sh
  real(kind=8), dimension(6) :: plandist
  integer, dimension(6) :: ba
  real(kind=8), dimension(3,27*natreal) :: rxyztot
  real(kind=8), dimension(:,:), allocatable :: rxyz
  logical :: perx,pery,perz
  integer :: nbl1,nbl2,nbl3,nbr1,nbr2,nbr3,nn
  type(cell) :: mesh
  type(box_iterator) :: bit

  density1=f_malloc([n01,n02,n03,nspden],id='density1')
  density2=f_malloc([n01,n02,n03,nspden],id='density2')
  potential2=f_malloc([n01,n02,n03],id='potential2')

  density=0.d0
  potential=0.d0

  pi = 4.d0*datan(1.d0)
  offset=0.d0

 if (SetEps.eq.1) then

  if (trim(geocode) == 'F') then

  if (Setrho.eq.1) then
! Set initial density as gaussian (or double gaussian with zero total charge) and potential as error function. It works only
! in a vacuum environment.

   sigma1 = 0.033d0*acell
   sigma2 = 2.0d0*sigma1
   x0 = 0.d0 ! hx*real(25-n01/2,kind=8)

         !Normalization
         factor1 = 1.d0/((sigma1**3)*sqrt((2.d0*pi)**3))
!         factor2 = 1.d0/((sigma2**3)*sqrt((2.d0*pi)**3))
         factor2 = 1.d0/((sigma2**3)*sqrt((2.d0*pi)**3))
         !gaussian function for the density.
         sumd=0.d0
         sump=0.d0
         do i3=1,n03
            x3 = hz*real(i3-n03/2,kind=8)
            do i2=1,n02
               x2 = hy*real(i2-n02/2,kind=8)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2,kind=8)
                  r12=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
                  r22 = x1*x1+x2*x2+x3*x3
                  do i=1,nspden
                     density(i1,i2,i3,i) = 1.d0/real(nspden,kind=8)*max(factor1*exp(-0.5d0*r12/(sigma1**2)),1d-24)&
                                         - 1.d0/real(nspden,kind=8)*max(factor2*exp(-0.5d0*r22/(sigma2**2)),1d-24)
                  sumd=sumd+density(i1,i2,i3,i)
                  end do
                  r1 = sqrt(r12)
                  r2 = sqrt(r22)
                  !Potential from a gaussian
                  if (r1 == 0.d0) then
                     potential(i1,i2,i3) = 2.d0/(sqrt(2.d0*pi)*sigma1)
                  else
                     call derf_local(derf_tt1,r1/(sqrt(2.d0)*sigma1))
                     potential(i1,i2,i3) = derf_tt1/r1
                  end if
                  if (r2 == 0.d0) then
                     potential(i1,i2,i3) = potential(i1,i2,i3) - 2.d0/(sqrt(2.d0*pi)*sigma2)
                  else
                     call derf_local(derf_tt2,r2/(sqrt(2.d0)*sigma2))
                     potential(i1,i2,i3) = potential(i1,i2,i3) - derf_tt2/r2
                  end if
                  sump=sump+potential(i1,i2,i3)
               end do
            end do
         end do

  else if (Setrho.eq.2) then

! Set initial potential as gaussian and density as the correct Generalized Laplace operator. It works with a gaussian epsilon.

   sigma = 0.03d0*acell
   x0 = 0.d0 ! hx*real(25-n01/2,kind=8)
!print *,'we should be here for vacuum'
         !Normalization
         factor = 1.d0/((sigma**3)*sqrt((2.d0*pi)**3))
         !gaussian function for the potential.
         sump=0.d0
         do i3=1,n03
            x3 = hz*real(i3-n03/2,kind=8)
            do i2=1,n02
               x2 = hy*real(i2-n02/2,kind=8)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2,kind=8)
!                  r2=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
                  r2 = x1*x1+x2*x2+x3*x3
                  potential(i1,i2,i3) = 1.d0/real(nspden,kind=8)*max(factor*exp(-0.5d0*r2/(sigma**2)),1d-24)
                   sump=sump+potential(i1,i2,i3)
               end do
            end do
          end do

         sumd=0.d0
         !analitic density calculation.
         do i3=1,n03
            x3 = hz*real(i3-n03/2,kind=8)
            do i2=1,n02
               x2 = hy*real(i2-n02/2,kind=8)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2,kind=8)
!                  r2=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
                  r2 = x1*x1+x2*x2+x3*x3
                  do i=1,nspden
                   density(i1,i2,i3,i) =(-1.d0/(4.d0*pi))*potential(i1,i2,i3)*(1.d0/(sigma**2))*&
                                        (r2*(1.d0/(sigmaeps**2))*(eps(i1,i2,i3)-erfR)+eps(i1,i2,i3)*(r2*(1.d0/(sigma**2))-3.d0))
                   sumd=sumd+density(i1,i2,i3,i)
                  end do
               end do
            end do
          end do

  else if (Setrho.eq.3) then

! Set initial potential as double gaussian and density as the correct 
!!Generalized Laplace operator. It works with a gaussian epsilon.

!   sigma1 = 0.033d0*acell
!   sigma2 = 2.d0*sigma1
   sigma2 = 0.033d0*acell
   sigma1 = 2.d0*sigma2
   x0 = 0.d0 ! hx*real(25-n01/2,kind=8)

         !Normalization
         factor1 = 1.d0/((sigma1**3)*sqrt((2.d0*pi)**3))
         factor2 = 1.d0/((sigma2**3)*sqrt((2.d0*pi)**3))
         !gaussian function for the potential.
         sump=0.d0
         do i3=1,n03
            x3 = hz*real(i3-n03/2,kind=8)
            do i2=1,n02
               x2 = hy*real(i2-n02/2,kind=8)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2,kind=8)
!                  r2=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
                  r2 = x1*x1+x2*x2+x3*x3
                  potential1(i1,i2,i3) = 1.d0/real(nspden,kind=8)*max(factor1*exp(-0.5d0*r2/(sigma1**2)),1d-24)
                  potential2(i1,i2,i3) = 1.d0/real(nspden,kind=8)*max(factor2*exp(-0.5d0*r2/(sigma2**2)),1d-24)
                  potential(i1,i2,i3) = potential1(i1,i2,i3) - potential2(i1,i2,i3)
                  sump=sump+potential(i1,i2,i3)
               end do
            end do
          end do

         sumd=0.d0
         !analitic density calculation.
         do i3=1,n03
            x3 = hz*real(i3-n03/2,kind=8)
            do i2=1,n02
               x2 = hy*real(i2-n02/2,kind=8)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2,kind=8)
!                  r2=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
                  r2 = x1*x1+x2*x2+x3*x3
                  do i=1,nspden
                   density1(i1,i2,i3,i) = (-1.d0/(4.d0*pi))*potential1(i1,i2,i3)*(1.d0/(sigma1**2))*&
                                        (r2*(1.d0/(sigmaeps**2))*(eps(i1,i2,i3)-erfR)+eps(i1,i2,i3)*(r2*(1.d0/(sigma1**2))-3.d0))
                   density2(i1,i2,i3,i) = (-1.d0/(4.d0*pi))*potential2(i1,i2,i3)*(1.d0/(sigma2**2))*&
                                        (r2*(1.d0/(sigmaeps**2))*(eps(i1,i2,i3)-erfR)+eps(i1,i2,i3)*(r2*(1.d0/(sigma2**2))-3.d0))
                   density(i1,i2,i3,i) = density1(i1,i2,i3,i) - density2(i1,i2,i3,i)

                   sumd=sumd+density(i1,i2,i3,i)
                  end do
               end do
            end do
          end do


  end if

   denval=0.d0

  else if (trim(geocode) == 'P') then

         !parameters for the test functions
         length=acell
         a=0.5d0/a_gauss**2
         !test functions in the three directions
         ifx=5
         ify=5
         ifz=5
         !parameters of the test functions
         ax=length
         ay=length
         az=length
         bx=2.d0!real(nu,kind=8)
         by=2.d0!real(nu,kind=8)
         bz=2.d0

              !plot of the functions used
              do i1=1,n03
                 x = hx*real(i1,kind=8)!valid if hy=hz
                 y = hz*real(i1,kind=8)
                 call functions(x,ax,bx,fx,fx1,fx2,ifx)
                 call functions(y,az,bz,fz,fz1,fz2,ifz)
                 write(20,'(1x,I8,4(1x,e22.15))')i1,fx,fx2,fz,fz2
              end do

         !Initialization of density and potential
         sumd=0.d0
         sump=0.d0
         do i3=1,n03
            x3 = hz*real(i3-n03/2-1,kind=8)
            call functions(x3,az,bz,fz,fz1,fz2,ifz)
            do i2=1,n02
               x2 = hy*real(i2-n02/2-1,kind=8)
               call functions(x2,ay,by,fy,fy1,fy2,ify)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2-1,kind=8)
                  call functions(x1,ax,bx,fx,fx1,fx2,ifx)
                  factor = (erfL-erfR)/(dexp(1.d0)**3-dexp(-1.d0)**3)
                  do i=1,nspden
                     density(i1,i2,i3,i) = 1.d0/real(nspden,kind=8)*((fx2*fy*fz+fx*fy2*fz+fx*fy*fz2)*eps(i1,i2,i3)+&
                                           factor*((fx1*fy*fz)**2+(fx*fy1*fz)**2+(fx*fy*fz1)**2))
                   sumd=sumd+density(i1,i2,i3,i)
                  end do
                  potential(i1,i2,i3) = -4.d0*pi*fx*fy*fz
                  sump=sump+potential(i1,i2,i3)
               end do
            end do
         end do

         denval=0.d0

      offset=0.d0
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               offset=offset+potential(i1,i2,i3)
            end do
         end do
      end do

      offset=offset*hx*hy*hz

!      write(*,*)'offset',offset

  else if (trim(geocode) == 'S') then

         !parameters for the test functions
         length=acell
         !parameters of the test functions
         ax=length
         az=length
         bx=2.d0!real(nu,kind=8)
         bz=2.d0!real(nu,kind=8)
         !non-periodic dimension
         ay=length
         by=a

         !Initialisation of density and potential
         sumd=0.d0
         sump=0.d0
         do i3=1,n03
            x3 = hz*real(i3-n03/2-1,kind=8)
            call functions(x3,az,bz,fz,fz1,fz2,ifz)
            do i2=1,n02
               x2 = hy*real(i2-n02/2-1,kind=8)
               call functions(x2,ay,by,fy,fy1,fy2,ify)
               do i1=1,n01
                  x1 = hx*real(i1-n02/2-1,kind=8)
                  call functions(x1,ax,bx,fx,fx1,fx2,ifx)
                  factor = (erfL-erfR)/(dexp(1.d0)**2)
                  do i=1,nspden
                     density(i1,i2,i3,i) = 1.d0/real(nspden,kind=8)/(4.d0*pi)*((fx2*fy*fz+fx*fy2*fz+fx*fy*fz2)*eps(i1,i2,i3)+&
                                           factor*((fx1*fy*fz)**2+(fx*fy1*fz)**2+(fx*fy*fz1)**2))
                   sumd=sumd+density(i1,i2,i3,i)
                  end do
                  potential(i1,i2,i3) = -fx*fy*fz
                  sump=sump+potential(i1,i2,i3)
               end do
            end do
         end do

         denval=0.d0

  end if

 else if (any(SetEps == [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])) then

! Set initial potential as gaussian and density as the correct Generalized
! Laplace operator. It works with a gaussian epsilon.
  sigma = 0.05d0*acell
  factor = 1.d0/((sigma**3)*sqrt((2.d0*pi)**3))

  shift(1)=hx*real(n01,kind=8)
  shift(2)=hy*real(n02,kind=8)
  shift(3)=hz*real(n03,kind=8)

  !buffers associated to the geocode
  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

!!$  call ext_buffers(perx,nbl1,nbr1)
!!$  call ext_buffers(pery,nbl2,nbr2)
!!$  call ext_buffers(perz,nbl3,nbr3)

  nbl1=0
  nbl2=0
  nbl3=0
  nbr1=0
  nbr2=0
  nbr3=0


!------------------------------------------------------------------------------------------------------
! Depending of Free, Periodic or Surface bc, image atoms are or not included.

  if (iproc==0) then
   call yaml_map('Starting initial potential calculation','Yes')
   call yaml_map('nbl1',nbl1)
   call yaml_map('nbl2',nbl2)
   call yaml_map('nbl3',nbl3)
   do iat=1,natreal
    call yaml_map('real input atoms',iat)
    call yaml_map('rxyz',rxyzreal(:,iat))
   end do
  end if

  px=0
  py=0
  pz=0
  if (perx) px=1
  if (pery) py=1
  if (perz) pz=1

  rxyztot(:,:)=0.d0

  i=0
  do iat=1,natreal

   ba(1:6)=0
   ! checking what are the image atoms to include in the calculation of the
   ! cavity.
   rv(1:3)=rxyzreal(1:3,iat)
   plandist(1)=dabs(rv(1))
   plandist(2)=dabs(shift(1)-rv(1))
   plandist(3)=dabs(rv(2))
   plandist(4)=dabs(shift(2)-rv(2))
   plandist(5)=dabs(rv(3))
   plandist(6)=dabs(shift(3)-rv(3))       !Normalization
   do ii=1,6
    valuemax=0.d0
    d=plandist(ii)
    r2=d*d
    value=factor*safe_exp(-0.5d0*r2/(sigma**2))
    if (value.gt.valuebc) then ! valuebc is the value to check on the box border to accept or refuse an image atom.
     if (abs(value).gt.valuemax) then
      valuemax=abs(value)
      imax=ii
     end if
     select case(ii)
     case (1)
      ba(1)=1*px
     case (2)
      ba(2)=-1*px
     case (3)
      ba(3)=1*py
     case (4)
      ba(4)=-1*py
     case (5)
      ba(5)=1*pz
     case (6)
      ba(6)=-1*pz
     end select
    end if
   end do

   do j=ba(6),ba(5)
    sh(3)=real(j,kind=8)*shift(3)
    do k=ba(4),ba(3)
     sh(2)=real(k,kind=8)*shift(2)
     do l=ba(2),ba(1)
      sh(1)=real(l,kind=8)*shift(1)
      rv(1:3)=rxyzreal(1:3,iat) + sh(1:3)
      i=i+1
      rxyztot(1:3,i)=rv(1:3)
     end do
    end do
   end do

  end do

  nat=i

  rxyz=f_malloc([3,nat],id='rxyz')
  rxyz(1:3,1:nat)=rxyztot(1:3,1:nat)

  if (iproc==0) then
     !write(*,*)plandist
     !write(*,'(1x,a,1x,e14.7,1x,a,1x,i4)')'Value max =',valuemax,'at bc side',imax
   call yaml_map('nat',nat)
   do iat=1,nat
    call yaml_map('atom',iat)
    call yaml_map('rxyz',rxyz(:,iat))
   end do
  end if

!------------------------------------------------------------------------------------------------------
! Starting the potential and density building for rxyztot atoms=real+image atoms (total
! natcurr) for periodic and surface boundary conditions or atoms=real for free bc.

  density(:,:,:,:)=0.d0
  potential(:,:,:)=0.d0
  potential1(:,:,:)=0.d0
  offset=0.d0
  
  do iat=1,nat
     nn=n01/2+1
     !gaussian function for the potential.
     sump=0.d0
     if (old) then
        do i3=1,n03
           x3=hz*(i3-1-nbl3)
           !x3 = hz*real(i3-n03/2,kind=8)
           do i2=1,n02
              x2 = hy*(i2-1-nbl2)
              !x2 = hy*real(i2-n02/2,kind=8)
              do i1=1,n01
                 x1 = hx*(i1-1-nbl1)
                 !x1 = hx*real(i1-n01/2,kind=8)
                 r2=(x1-rxyz(1,iat))**2+(x2-rxyz(2,iat))**2+(x3-rxyz(3,iat))**2
                 !r2 = x1*x1+x2*x2+x3*x3
                 !      potential(i1,i2,i3) = 1.d0/real(nspden,kind=8)*max(factor*exp(-0.5d0*r2/(sigma**2)),1d-24)
                 potential1(i1,i2,i3) = factor*safe_exp(-0.5d0*r2/(sigma**2))
                 potential(i1,i2,i3) = potential(i1,i2,i3) + potential1(i1,i2,i3)
                 sump=sump+potential(i1,i2,i3)
                 offset=offset+potential(i1,i2,i3)
                 !      if (i1.eq.nn.and.i2.eq.nn.and.i3.eq.nn) then
                 !       write(*,*)i1,i2,i3
                 !       write(*,*)x1,x2,x3,rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),potential(i1,i2,i3)
                 !      end if
              end do
           end do
        end do
        offset=offset*hx*hy*hz
        if (iproc == 0) call yaml_map('offset',offset)
        switch=0.0d0
        if (lin_PB) then
           switch=1.0d0
        end if

        sumd=0.d0
        !analitic density calculation.
        do i3=1,n03
           x3=hz*(i3-1-nbl3)
           !x3 = hz*real(i3-n03/2,kind=8)
           r_v(3)=x3-rxyz(3,iat)
           do i2=1,n02
              x2 = hy*(i2-1-nbl2)
              !x2 = hy*real(i2-n02/2,kind=8)
              r_v(2)=x2-rxyz(2,iat)
              do i1=1,n01
                 x1 = hx*(i1-1-nbl1)
                 !x1 = hx*real(i1-n01/2,kind=8)
                 r_v(1)=x1-rxyz(1,iat)
                 r2=(x1-rxyz(1,iat))**2+(x2-rxyz(2,iat))**2+(x3-rxyz(3,iat))**2
                 !r2 = x1*x1+x2*x2+x3*x3
                 k1=0.d0
                 do i=1,3
                    k1 = k1 + dlogeps(i,i1,i2,i3)*potential1(i1,i2,i3)*(-r_v(i)/(sigma**2))
                 end do
                 k2 = potential1(i1,i2,i3)*(r2/(sigma**2)-3.d0)/(sigma**2)
                 do i=1,nspden
                    !        density(i1,i2,i3,i) =(-1.d0/(4.d0*pi))*eps(i1,i2,i3)*(k1+k2)&
                    !                             -switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(potential(i1,i2,i3))
                    dens = (-1.d0/(4.d0*pi))*eps(i1,i2,i3)*(k1+k2)&
                         -switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(potential1(i1,i2,i3))
                    !                             +switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dsinh(multp*potential(i1,i2,i3))
                    !                             +switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dtanh(multp*potential(i1,i2,i3))
                    !                             +switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*potential(i1,i2,i3)
                    !                             +switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*(potential(i1,i2,i3)**2)
                    density(i1,i2,i3,i) = density(i1,i2,i3,i) + dens
                    sumd=sumd+density(i1,i2,i3,i)
                 end do
              end do
           end do
        end do
     else
        switch=0.0d0
        if (lin_PB) switch=1.0d0

        mesh=cell_new(geocode,[n01,n02,n03],[hx,hy,hz])

        sumd=0.d0
        bit=box_iter(mesh,origin=rxyz(:,iat))
        do while(box_next_point(bit))
           r2=square(mesh,bit%rxyz)
           potential1(bit%i,bit%j,bit%k) = factor*safe_exp(-0.5d0*r2/(sigma**2))
           potential(bit%i,bit%j,bit%k) = potential(bit%i,bit%j,bit%k) + potential1(bit%i,bit%j,bit%k)
           sump=sump+potential(bit%i,bit%j,bit%k)
           offset=offset+potential(bit%i,bit%j,bit%k)
           k1=dotp(mesh,dlogeps(:,bit%i,bit%j,bit%k),bit%rxyz)
           k1=-potential1(bit%i,bit%j,bit%k)*k1/(sigma**2)
           k2 = potential1(bit%i,bit%j,bit%k)*(r2/(sigma**2)-3.d0)/(sigma**2)
           dens = (-oneofourpi)*eps(bit%i,bit%j,bit%k)*(k1+k2)&
                -switch*((eps(bit%i,bit%j,bit%k)-1.0d0)/(eps0-1.0d0))*PB_charge(potential1(bit%i,bit%j,bit%k))
           do i=1,nspden
              density(bit%i,bit%j,bit%k,i) = density(bit%i,bit%j,bit%k,i) + dens
              sumd=sumd+density(bit%i,bit%j,bit%k,i)
           end do
        end do
        offset=offset*mesh%volume_element!hx*hy*hz
        if (iproc == 0) call yaml_map('offset',offset)
     end if

  end do ! loop iat over atoms

  call f_free(rxyz)

else if (SetEps.eq.18) then

   ! Set initial potential as gaussian and density as the correct Generalized
   ! Laplace operator. It works with a gaussian epsilon.

  offset=0.d0
  sigma = 0.05d0*acell
  x0 = 0.d0 ! hx*real(25-n01/2,kind=8)
  print *,'we should be here for old analytical functions'
  !Normalization
  factor = 1.d0/((sigma**3)*sqrt((2.d0*pi)**3))
  !gaussian function for the potential.
  sump=0.d0
  do i3=1,n03
     x3 = hz*real(i3-n03/2,kind=8)
     do i2=1,n02
        x2 = hy*real(i2-n02/2,kind=8)
        do i1=1,n01
           x1 = hx*real(i1-n01/2,kind=8)
           !                  r2=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
           r2 = x1*x1+x2*x2+x3*x3
           potential(i1,i2,i3) = 1.d0/real(nspden,kind=8)*max(factor*exp(-0.5d0*r2/(sigma**2)),1d-24)
           sump=sump+potential(i1,i2,i3)
           offset=offset+potential(i1,i2,i3)
        end do
     end do
  end do
  offset=offset*hx*hy*hz

  !      write(*,*)'offset',offset

  switch=0.0d0
  if (lin_PB) then
     switch=1.0d0
  end if

  sumd=0.d0
  !analitic density calculation.
  do i3=1,n03
     x3 = hz*real(i3-n03/2,kind=8)
     r_v(3)=x3
     do i2=1,n02
        x2 = hy*real(i2-n02/2,kind=8)
        r_v(2)=x2
        do i1=1,n01
           x1 = hx*real(i1-n01/2,kind=8)
           r_v(1)=x1
           !                  r2=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
           r2 = x1*x1+x2*x2+x3*x3
           k1=0.d0
           do i=1,3
              k1 = k1 + dlogeps(i,i1,i2,i3)*potential(i1,i2,i3)*(-r_v(i)/(sigma**2))
           end do
           k2 = potential(i1,i2,i3)*(r2/(sigma**2)-3.d0)/(sigma**2)
           do i=1,nspden
              density(i1,i2,i3,i) =(-1.d0/(4.d0*pi))*eps(i1,i2,i3)*(k1+k2)&
                   -switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*PB_charge(potential(i1,i2,i3))
              !                                 +switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dsinh(multp*potential(i1,i2,i3))
              !                                 +switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*dtanh(multp*potential(i1,i2,i3))
              !                                 +switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*potential(i1,i2,i3)
              !                                 +switch*((eps(i1,i2,i3)-1.0d0)/(eps0-1.0d0))*multp*(potential(i1,i2,i3)**2)
              sumd=sumd+density(i1,i2,i3,i)
           end do
        end do
     end do
  end do

end if

!plot of the starting conditions
 unt=f_get_free_unit(21)
 call f_open_file(unt,file='initial.dat')
 !i1=n01/2
 i1=(n01-1)/2+1 
 !i1=1
 do i2=1,n02
    do i3=1,n03
       write(unt,'(2(1x,I4),3(1x,e14.7))')i2,i3,density(i1,i2,i3,1),potential(i1,i2,i3),eps(i1,i2,i3)
    end do
 end do
 call f_close(unt)

 unt=f_get_free_unit(22)
 call f_open_file(unt,file='initial_line.dat')
 i1=(n01-1)/2+1
 !i1=1
 i3=(n03-1)/2+1
! i3=1  
 do i2=1,n02
!  write(unt,'(1x,I8,3(1x,e22.15))') i2,density(n01/2,i2,n03/2,1),potential(n01/2,i2,n03/2),eps(n01/2,i2,n03/2)
  y=hy*real(i2-1,kind=8)
  write(unt,'(1x,I8,4(1x,e22.15))') i2,y,density(i1,i2,i3,1),potential(i1,i2,i3),eps(i1,i2,i3)
 end do
 call f_close(unt)

 !calculate hartree energy
 einit=0.d0
 if (.not. old) then
 bit=box_iter(mesh)
 do while(box_next_point(bit))
    einit= einit+density(bit%i,bit%j,bit%k,1)*potential(bit%i,bit%j,bit%k) 
 end do
 einit=0.5*mesh%volume_element*einit
 if (iproc ==0) then
    call yaml_map('Total Charge',sumd*mesh%volume_element)
    call yaml_map('Potential monopole',sump*mesh%volume_element)
    call yaml_map('Potential at the boundary 1 n02/2 1',&
         potential(1,n02/2,1))
    call yaml_map('Density at the boundary 1 n02/2 1',density(1,n02/2,1,1))
 end if

 else
    do i3=1,n03
       do i2=1,n02
          do i1=1,n01
             einit= einit + density(i1,i2,i3,1)*potential(i1,i2,i3) 
          end do
       end do
    end do
    einit=0.5*hx*hy*hz*einit
  if (iproc ==0) then
     call yaml_map('Total Charge',sumd*hx*hy*hz)
     call yaml_map('Potential monopole',sump*hx*hy*hz)
     call yaml_map('Potential at the boundary 1 n02/2 1',&
          potential(1,n02/2,1))
     call yaml_map('Density at the boundary 1 n02/2 1',density(1,n02/2,1,1))
   end if

 end if
  !write(*,*) 'charge sumd',sumd*hx*hy*hz,'potential sum',sump*hx*hy*hz
  !write(*,'(1x,a,1x,e14.7)')'Potential at the boundary 1 n02/2 1',poteantial(1,n02/2,1)
  !write(*,'(1x,a,1x,e14.7)')'Density at the boundary 1 n02/2 1',density(1,n02/2,1,1)

  call f_free(density1)
  call f_free(density2)
  call f_free(potential2)

end subroutine SetInitDensPot

subroutine functions(x,a,b,f,f1,f2,whichone)
      implicit none
      integer, intent(in) :: whichone
      real(kind=8), intent(in) :: x,a,b
      real(kind=8), intent(out) :: f,f1,f2
      !local variables
      real(kind=8) :: r,r2,y,yp,ys,factor,pi,g,h,g1,g2,h1,h2
      real(kind=8) :: length,frequency,nu,sigma,agauss

      pi = 4.d0*datan(1.d0)
      select case(whichone)
      case(1)
         !constant
         f=1.d0
         f2=0.d0
      case(2)
         !gaussian of sigma s.t. a=1/(2*sigma^2)
         r2=a*x**2
         f=dexp(-r2)
         f2=(-2.d0*a+4.d0*a*r2)*dexp(-r2)
      case(3)
         !gaussian "shrinked" with a=length of the system
         length=a
         r=pi*x/length
         y=dtan(r)
         yp=pi/length*1.d0/(dcos(r))**2
         ys=2.d0*pi/length*y*yp
         factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
         f2=factor*dexp(-y**2)
         f=dexp(-y**2)
      case(4)
         !cosine with a=length, b=frequency
         length=a
         frequency=b
         r=frequency*pi*x/length
         f=dcos(r)
         f2=-(frequency*pi/length)**2*dcos(r)
      case(5)
         !exp of a cosine, a=length
         nu=2.d0
         r=pi*nu/a*x
         y=dcos(r)
         yp=dsin(r)
         f=dexp(y)
         factor=(pi*nu/a)**2*(-y+yp**2)
         f1=-f*yp*(pi*nu/a)
         f2=factor*f
      case(6)
         !gaussian times "shrinked" gaussian, sigma=length/10
         length=a
         r=pi*x/length
         y=dtan(r)
         yp=pi/length*1.d0/(dcos(r))**2
         ys=2.d0*pi/length*y*yp
         factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
         g=dexp(-y**2)
         g1=-2.d0*y*yp*g
         g2=factor*dexp(-y**2)

         sigma=length/10
         agauss=0.5d0/sigma**2
         r2=agauss*x**2
         h=max(dexp(-r2),1.0d-24)
         h1=-2.d0*agauss*x*h
         h2=(-2.d0*agauss+4.d0*agauss*r2)*h
         f=max(g,1.d-24)*h
         f1=g1*h+max(g,1.d-24)*h1
         f2=g2*h+max(g,1.d-24)*h2+2.d0*g1*h1
      case(7)
         !sine with a=length, b=frequency
         length=a
         frequency=b
         r=frequency*pi*x/length
         f=dsin(r)
         f2=-(frequency*pi/length)**2*dsin(r)
      end select

end subroutine functions

subroutine SetEpsilon(n01,n02,n03,nspden,nord,nat,iproc,acell,a_gauss,hx,hy,hz,&
     erfL,erfR,sigmaeps,SetEps,geocode,PSol,eps,dlogeps,oneoeps,oneosqrteps,corr,&
     rhoele,rad_cav,rxyz,radii,delta)

  use dynamic_memory
  use yaml_output

  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  integer, intent(in) :: nord
  integer, intent(in) :: nat,iproc
  real(kind=8), intent(in) :: acell,a_gauss,hx,hy,hz,erfL,erfR,sigmaeps,rad_cav,delta
  integer, intent(in) :: SetEps
  character(len=2), intent(in) :: geocode
  character(len=4), intent(in) :: PSol

  !> dielectric function. Needed for non VAC methods, given in full dimensions
  real(kind=8), dimension(n01,n02,n03), intent(out) :: eps
  !> logarithmic derivative of epsilon. Needed for PCG method.
  !! if absent, it will be calculated from the array of epsilon
  real(kind=8), dimension(3,n01,n02,n03), intent(out) :: dlogeps
  !> inverse of epsilon. Needed for PI method.
  !! if absent, it will be calculated from the array of epsilon
  real(kind=8), dimension(n01,n02,n03), intent(out) :: oneoeps
  !> inverse square root of epsilon. Needed for PCG method.
  !! if absent, it will be calculated from the array of epsilon
  real(kind=8), dimension(n01,n02,n03), intent(out) :: oneosqrteps
  !> correction term of the Generalized Laplacian
  !! if absent, it will be calculated from the array of epsilon
  real(kind=8), dimension(n01,n02,n03), intent(out) :: corr

  real(kind=8), dimension(n01,n02,n03,nspden), intent(in) :: rhoele
  real(kind=8), dimension(3,nat), intent(inout) :: rxyz
  real(kind=8), dimension(nat), intent(inout) :: radii

  ! local variables.
  real(kind=8), dimension(n01,n02,n03,nspden) :: edens
  real(kind=8), dimension(n01,n02,n03,nspden,3) :: nabla_edens ! Nabla of the electron density.
  real(kind=8), dimension(n01,n02,n03,nspden) :: ddt_edens ! Laplacian of the electron density.
  real(kind=8), dimension(n01,n02,n03,3) :: deps ! Nabla of the electron density.
  integer :: i,i1,i2,i3,ifx,ify,ifz,isp,iat
  real(kind=8) :: edensmax = 0.0035d0!!!!
  real(kind=8) :: edensmin = 0.0001d0
  real(kind=8), parameter :: eps0 = 78.36d0
  real(kind=8) :: x1,x2,x3,r,t,pi,r2,sigma,x0,factor,length,oneoeps0,oneosqrteps0
  real(kind=8) :: x,y,fx,fx2,fy,fy2,fz,fz2,a,ax,ay,az,bx,by,bz,tt,fx1,fy1,fz1
  real(kind=8) :: fact1,fact2,fact3,dtx,d2,dd,coeff,coeff1

  open(unit=21,file='epsilon.dat',status='unknown')
  open(unit=22,file='epsilon_line.dat',status='unknown')

  nabla_edens=0.d0
  ddt_edens=0.d0

 if (SetEps.eq.1) then

  if (trim(geocode) == 'F') then

   sigma = sigmaeps
   x0 = 0.d0 ! hx*real(25-n01/2,kind=8)

         !Normalization
         factor = erfL - erfR
         !factor = 1.d0/((sigma**3)*sqrt((2.d0*pi)**3))
         !gaussian function
         do i3=1,n03
            x3 = hz*real(i3-n03/2,kind=8)
            do i2=1,n02
               x2 = hy*real(i2-n02/2,kind=8)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2,kind=8)
                  r2=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
                  !r2 = x1*x1+x2*x2+x3*x3
                     eps(i1,i2,i3) = max(factor*exp(-0.5d0*r2/(sigma**2)),1d-24) + erfR
               end do
            end do
         end do

  else if (trim(geocode) == 'P') then

         !parameters for the test functions
         length=acell
         a=0.5d0/a_gauss**2
         !test functions in the three directions
         ifx=5
         ify=5
         ifz=5
         !parameters of the test functions
         ax=length
         ay=length
         az=length
         bx=2.d0!real(nu,kind=8)
         by=2.d0!real(nu,kind=8)
         bz=2.d0

              !plot of the functions used
              do i1=1,n03
                 x = hx*real(i1-n01/2-1,kind=8)!valid if hy=hz
                 y = hy*real(i1-n01/2-1,kind=8)
                 call functions(x,ax,bx,fx,fx1,fx2,ifx)
                 call functions(y,ay,by,fy,fy1,fy2,ify)
                 write(36,'(1x,I8,6(1x,e22.15))')i1,fx,fx1,fx2,fy,fy1,fy2
              end do

         !Initialization of dielectric constant for periodic boundary
         !conditions.
         do i3=1,n03
            x3 = hz*real(i3-n03/2-1,kind=8)
            call functions(x3,az,bz,fz,fz1,fz2,ifz)
            do i2=1,n02
               x2 = hy*real(i2-n02/2-1,kind=8)
               call functions(x2,ay,by,fy,fy1,fy2,ify)
               do i1=1,n01
                  x1 = hx*real(i1-n01/2-1,kind=8)
                  call functions(x1,ax,bx,fx,fx1,fx2,ifx)
                  factor = (erfL-erfR)/(dexp(1.d0)**3-dexp(-1.d0)**3)
                  eps(i1,i2,i3) = factor*(fx*fy*fz-dexp(-1.d0)**3) + erfR
               end do
            end do
         end do


  else if (trim(geocode) == 'S') then

         !parameters for the test functions
         length=acell
         a=0.5d0/a_gauss**2
         !test functions in the three directions
         ifx=5
         ifz=5
         !non-periodic dimension
         ify=6
         !parameters of the test functions
         ax=length
         az=length
         bx=2.d0!real(nu,kind=8)
         bz=2.d0!real(nu,kind=8)
         !non-periodic dimension
         ay=length
         by=a

              !plot of the functions used
              do i1=1,n03
                 x = hx*real(i1-n01/2-1,kind=8)!valid if hy=hz
                 y = hy*real(i1-n01/2-1,kind=8)
                 call functions(x,ax,bx,fx,fx1,fx2,ifx)
                 call functions(y,ay,by,fy,fy1,fy2,ify)
                 write(36,'(1x,I8,6(1x,e22.15))')i1,fx,fx1,fx2,fy,fy1,fy2
              end do

         !Initialisation of density and potential
         do i3=1,n03
            x3 = hz*real(i3-n03/2-1,kind=8)
            call functions(x3,az,bz,fz,fz1,fz2,ifz)
            do i2=1,n02
               x2 = hy*real(i2-n02/2-1,kind=8)
               call functions(x2,ay,by,fy,fy1,fy2,ify)
               do i1=1,n01
                  x1 = hx*real(i1-n02/2-1,kind=8)
                  call functions(x1,ax,bx,fx,fx1,fx2,ifx)
                  factor = (erfL-erfR)/(dexp(1.d0)**2)
                  eps(i1,i2,i3) = factor*(fx*fy*fz) + erfR
               end do
            end do
         end do

  end if

 else if (SetEps.eq.2 .or. SetEps.eq.3) then

   if (SetEps.eq.2) then
    call SetEledens(n01,n02,n03,nspden,nord,acell,a_gauss,hx,hy,hz,SetEps,edens,nabla_edens,ddt_edens)
   else if (SetEps.eq.3) then
    edens(:,:,:,:) = rhoele(:,:,:,:)
    call fssnord3DmatNabla(geocode,n01,n02,n03,nspden,hx,hy,hz,edens,nabla_edens,nord,acell)
    call fssnord3DmatDiv(geocode,n01,n02,n03,nspden,hx,hy,hz,nabla_edens,ddt_edens,nord,acell)
   end if

!   r2=(rad_cav/0.52917721092d0)**2
!   edensmax = max(exp(-0.5d0*r2/(0.16d0)),1d-24)
!   delta=12.d0*max(hx,hy,hz)
!   r2=((rad_cav+delta)/0.52917721092d0)**2
!   edensmin = max(exp(-0.5d0*r2/(0.16d0)),1d-24)

   pi = 4.d0*datan(1.d0)
   r=0.d0
   t=0.d0
   oneoeps0=1.d0/eps0
   oneosqrteps0=1.d0/dsqrt(eps0)
   fact1=2.d0*pi/(dlog(edensmax)-dlog(edensmin))
   fact2=(dlog(eps0))/(2.d0*pi)
   fact3=(dlog(eps0))/(dlog(edensmax)-dlog(edensmin))

   if ( trim(PSol)=='PCG') then

   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01

      if (dabs(edens(i1,i2,i3,isp)).gt.edensmax) then
       eps(i1,i2,i3)=1.d0
       oneosqrteps(i1,i2,i3)=1.d0
       do i=1,3
        dlogeps(i,i1,i2,i3)=0.d0
       end do
       corr(i1,i2,i3)=0.d0
      else if (dabs(edens(i1,i2,i3,isp)).lt.edensmin) then
       eps(i1,i2,i3)=eps0
       oneosqrteps(i1,i2,i3)=oneosqrteps0
       do i=1,3
        dlogeps(i,i1,i2,i3)=0.d0
       end do
       corr(i1,i2,i3)=0.d0
      else
       r=fact1*(dlog(edensmax)-dlog(dabs(edens(i1,i2,i3,isp))))
       t=fact2*(r-dsin(r))
       eps(i1,i2,i3)=dexp(t)
       oneosqrteps(i1,i2,i3)=dexp(-0.5d0*t)
       coeff=fact3*(1.d0-dcos(r))
       dtx=-coeff/dabs(edens(i1,i2,i3,isp))
       d2=0.d0
       do i=1,3
        dlogeps(i,i1,i2,i3)=dtx*nabla_edens(i1,i2,i3,isp,i)
        d2 = d2+nabla_edens(i1,i2,i3,isp,i)**2
       end do
        dd = ddt_edens(i1,i2,i3,isp)
       coeff1=(0.5d0*(coeff**2)+fact3*fact1*dsin(r)+coeff)/((edens(i1,i2,i3,isp))**2)
       corr(i1,i2,i3)=(0.125d0/pi)*dexp(t)*(coeff1*d2+dtx*dd)
      end if

     end do
    end do
   end do

   else if ( trim(PSol)=='PI') then

   isp=1
   do i3=1,n03
    do i2=1,n02
     do i1=1,n01

      if (dabs(edens(i1,i2,i3,isp)).gt.edensmax) then
       eps(i1,i2,i3)=1.d0
       oneoeps(i1,i2,i3)=1.d0
       do i=1,3
        dlogeps(i,i1,i2,i3)=0.d0
       end do
      else if (dabs(edens(i1,i2,i3,isp)).lt.edensmin) then
       eps(i1,i2,i3)=eps0
       oneoeps(i1,i2,i3)=oneoeps0
       do i=1,3
        dlogeps(i,i1,i2,i3)=0.d0
       end do
      else
       r=fact1*(dlog(edensmax)-dlog(dabs(edens(i1,i2,i3,isp))))
       t=fact2*(r-dsin(r))
       eps(i1,i2,i3)=dexp(t)
       oneoeps(i1,i2,i3)=dexp(-t)
       coeff=fact3*(1.d0-dcos(r))
       dtx=-coeff/dabs(edens(i1,i2,i3,isp))
       do i=1,3
        dlogeps(i,i1,i2,i3)=dtx*nabla_edens(i1,i2,i3,isp,i)
       end do
      end if

     end do
    end do
   end do

   end if

else if (SetEps ==4) then

   do iat=1,nat
    radii(iat) = radii(iat) + delta
   end do

      if (iproc==0) then
       call yaml_map('Delta cavity',delta)
       call yaml_map('radii',radii)
       do iat=1,nat
        call yaml_map('atom',iat)
        call yaml_map('rxyz',rxyz(:,iat))
       end do
      end if


!      call Eps_rigid_cavity([n01,n02,n03],nspden,nord,acell,[hx,hy,hz],nat,rxyz,radii,eps0,delta,eps,dlogeps,oneoeps,oneosqrteps,corr)

      call Eps_rigid_cavity_multiatoms([n01,n02,n03],nspden,nord,acell,[hx,hy,hz],&
           nat,rxyz,radii,eps0,delta,eps,dlogeps,oneoeps,oneosqrteps,corr,geocode,iproc)

!      call Eps_rigid_cavity_new([n01,n02,n03],nspden,nord,acell,[hx,hy,hz],nat,rxyz,radii,eps0,delta,eps,dlogeps,oneoeps,oneosqrteps,corr)
!      call Eps_rigid_cavity_new_multiatoms([n01,n02,n03],nspden,nord,acell,[hx,hy,hz],nat,rxyz,radii,eps0,delta,eps,dlogeps,oneoeps,oneosqrteps,corr)
!      call Eps_rigid_cavity_new2([n01,n02,n03],nspden,nord,acell,[hx,hy,hz],nat,rxyz,radii,eps0,delta,eps,dlogeps,oneoeps,oneosqrteps,corr)
!      call Eps_rigid_cavity_new2_multiatoms([n01,n02,n03],nspden,nord,acell,[hx,hy,hz],nat,rxyz,radii,eps0,delta,eps,dlogeps,oneoeps,oneosqrteps,corr)
!!$
!!$print *,'New dlogeps calculation'
!!$      call fssnord3DmatNabla3var(n01,n02,n03,nspden,hx,hy,hz,eps,deps,nord,acell)
!!$      isp=1
!!$      do i3=1,n03
!!$       do i2=1,n02
!!$        do i1=1,n01
!!$         do i=1,3
!!$          dlogeps(i,i1,i2,i3)=deps(i1,i2,i3,i)/eps(i1,i2,i3)
!!$         end do
!!$       end do
!!$      end do
!!$     end do

end if


!!$     i3=1!n03/2
!!$     do i2=1,n02
!!$      do i1=1,n01
!!$       write(21,'(2(1x,I4),2(1x,e14.7))')i1,i2,eps(i1,i2,i3),eps(i1,i2,n03/2)
!!$      end do
!!$      write(21,*)
!!$     end do

     !i1=n01/2+1
     i1=(n01-1)/2+1
     !i1=1
     do i2=1,n02
      do i3=1,n03
       write(21,'(2(1x,I4),2(1x,e14.7))')i2,i3,eps(i2,i1,i3),eps(i1,i2,i3)
      end do
      write(21,*)
     end do


     i1=(n01-1)/2+1
     i3=(n03-1)/2+1
     do i2=1,n02
      y=hy*real(i2-1,kind=8)
      write(22,'(1x,I8,6(1x,e22.15))') i2,y,eps(i1,i2,i3),dlogeps(1,i1,i2,i3),dlogeps(2,i1,i2,i3),dlogeps(3,i1,i2,i3),corr(i1,i2,i3)
     end do

  close(unit=21)
  close(unit=22)

end subroutine SetEpsilon


!> Calculates the value of the dielectric function for a smoothed cavity 
!! given a set of centres and radii.
!! Need the epsilon0 as well as the radius of the cavit and its smoothness
subroutine Eps_rigid_cavity(ndims,nspden,nord,acell,hgrids,nat,rxyz,radii,epsilon0,delta,eps,dlogeps,oneoeps,oneosqrteps,corr)
  implicit none
  integer, intent(in) :: nat !< number of centres defining the cavity
  integer, intent(in) :: nspden
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell
  real(kind=8), intent(in) :: epsilon0 !< dielectric constant of th solvent
  real(kind=8), intent(in) :: delta !< smoothness factor of the cavity
  integer, dimension(3), intent(in) :: ndims   !< dimensions of the simulation box
  real(kind=8), dimension(3), intent(in) :: hgrids !< grid spacings
  real(kind=8), dimension(nat), intent(in) :: radii !< radii of each of the atoms
  !> position of all the atoms in the grid coordinates
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: eps !< dielectric function
  real(kind=8), dimension(3,ndims(1),ndims(2),ndims(3)), intent(out) :: dlogeps !< dlogeps
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneoeps !< inverse of epsilon. Needed for PI method.
  !> inverse square root of epsilon. Needed for PCG method.
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneosqrteps 
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: corr !< correction term of the Generalized Laplacian.
  !local variables
  integer :: i,i1,i2,i3,iat
  real(kind=8) :: r2,x,y2,z2,d2,y,z,eps_min,eps1,pi,de2,dde,d1,oneod,h,coeff
  real(kind=8), dimension(3) :: deps,ddeps,v
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)) :: v_de2
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)) :: v_ddeps

  pi = 4.d0*datan(1.d0)


  do i3=1,ndims(3)
     z=hgrids(3)*i3 !(i3-1) for 0 axis start.
     z2=z*z
     v(3)=z
     do i2=1,ndims(2)
        y=hgrids(2)*i2 !*(i2-1) for 0 axis start.
        y2=y*y
        v(2)=y
        do i1=1,ndims(1)
           x=hgrids(1)*i1 !(i1-1) for 0 axis start.
           v(1)=x
           r2=x*x+y2+z2
           !choose the closest atom
           eps_min=1.d100
           do iat=1,nat
            d2=(x-rxyz(1,iat))**2+(y-rxyz(2,iat))**2+(z-rxyz(3,iat))**2
            if (d2.eq.0.d0) then
               d2=1.0d-30
               eps1=epsl(sqrt(d2),radii(iat),delta,epsilon0)
               d1=d1eps(sqrt(d2),radii(iat),delta,epsilon0)
               oneod=1.d0/sqrt(d2)
!               coeff=-2.d0*((sqrt(d2)-radii(iat))/(delta**2))
!               coeff=oneod+2.d0*((sqrt(d2)-radii(iat))/(delta**2))
               do i=1,3
                h=0.d0
                deps(i) =0.d0
                ddeps(i)=0.d0
               end do
               eps_min=eps1
             exit
            else
               oneod=1.d0/sqrt(d2)
               eps1=epsl(sqrt(d2),radii(iat),delta,epsilon0)
              if (eps1< eps_min) then
                 d1=d1eps(sqrt(d2),radii(iat),delta,epsilon0)
                 coeff=oneod+2.d0*((sqrt(d2)-radii(iat))/(delta**2))
                 do i=1,3
                  h=(v(i)-rxyz(i,iat))*oneod
                  deps(i) =d1*h
                  ddeps(i)=d1*(oneod-coeff*(h**2))
                 end do
                 eps_min=eps1
              end if
              if (abs(eps_min-1.d0) < epsilon(1.d0)) exit
            end if
           end do

           if (nat==0) then
              eps_min=1.d0
              deps=0.d0
              ddeps=0.d0
           end if
           eps(i1,i2,i3)=eps_min
           oneoeps(i1,i2,i3)=1.d0/eps_min
           oneosqrteps(i1,i2,i3)=1.d0/sqrt(eps_min)
           de2=0.d0
           dde=0.d0
           do i=1,3
            dlogeps(i,i1,i2,i3)=deps(i)/eps_min
            de2 = de2 + deps(i)**2
            dde = dde + ddeps(i)
           end do
!            de2 = d1**2
!            dde = d1*(3.d0*oneod-coeff*(h**2))
            corr(i1,i2,i3)=-(0.125d0/pi)*(0.5d0/eps_min*de2-dde)
        end do
     end do
  end do

!!$  call fssnordEpsilonDerivative(ndims(1),ndims(2),ndims(3),1,hgrids(1),hgrids(2),hgrids(3),eps,v_de2,v_ddeps,nord,acell)
!!$
!!$  do i3=1,ndims(3)
!!$     z=hgrids(3)*i3 !(i3-1) for 0 axis start.
!!$     z2=z*z
!!$     v(3)=z
!!$     do i2=1,ndims(2)
!!$        y=hgrids(2)*i2 !*(i2-1) for 0 axis start.
!!$        y2=y*y
!!$        v(2)=y
!!$        do i1=1,ndims(1)
!!$           x=hgrids(1)*i1 !(i1-1) for 0 axis start.
!!$           v(1)=x
!!$           r2=x*x+y2+z2
!!$           !choose the closest atom
!!$           eps_min=1.d100
!!$           do iat=1,nat
!!$            d2=(x-rxyz(1,iat))**2+(y-rxyz(2,iat))**2+(z-rxyz(3,iat))**2
!!$            if (d2.eq.0.d0) then
!!$               d2=1.0d-30
!!$               eps1=epsl(sqrt(d2),radii(iat),delta,epsilon0)
!!$               d1=d1eps(sqrt(d2),radii(iat),delta,epsilon0)
!!$               oneod=1.d0/sqrt(d2)
!!$!               coeff=-2.d0*((sqrt(d2)-radii(iat))/(delta**2))
!!$!               coeff=oneod+2.d0*((sqrt(d2)-radii(iat))/(delta**2))
!!$               do i=1,3
!!$                h=0.d0
!!$                deps(i) =0.d0
!!$                ddeps(i)=0.d0
!!$               end do
!!$               eps_min=eps1
!!$             exit
!!$            else
!!$               oneod=1.d0/sqrt(d2)
!!$               eps1=epsl(sqrt(d2),radii(iat),delta,epsilon0)
!!$              if (eps1< eps_min) then
!!$                 d1=d1eps(sqrt(d2),radii(iat),delta,epsilon0)
!!$                 coeff=oneod+2.d0*((sqrt(d2)-radii(iat))/(delta**2))
!!$                 do i=1,3
!!$                  h=(v(i)-rxyz(i,iat))*oneod
!!$                  deps(i) =d1*h
!!$                  ddeps(i)=d1*(oneod-coeff*(h**2))
!!$                 end do
!!$                 eps_min=eps1
!!$              end if
!!$              if (abs(eps_min-1.d0) < epsilon(1.d0)) exit
!!$            end if
!!$           end do
!!$
!!$           if (nat==0) then
!!$              eps_min=1.d0
!!$              deps=0.d0
!!$              ddeps=0.d0
!!$           end if
!!$           eps(i1,i2,i3)=eps_min
!!$           oneoeps(i1,i2,i3)=1.d0/eps_min
!!$           oneosqrteps(i1,i2,i3)=1.d0/sqrt(eps_min)
!!$           dde=v_ddeps(i1,i2,i3)
!!$!           de2=v_de2(i1,i2,i3)
!!$           de2=0.d0
!!$!           dde=0.d0
!!$           do i=1,3
!!$            dlogeps(i,i1,i2,i3)=deps(i)/eps_min
!!$            de2 = de2 + deps(i)**2
!!$!            dde = dde + ddeps(i)
!!$           end do
!!$!            de2 = d1**2
!!$!            dde = d1*(3.d0*oneod-coeff*(h**2))
!!$            corr(i1,i2,i3)=-(0.125d0/pi)*(0.5d0/eps_min*de2-dde)
!!$        end do
!!$     end do
!!$  end do

  contains

    pure function epsl(r,rc,delta,epsilon0)
      implicit none
      real(kind=8), intent(in) :: r,rc,delta,epsilon0
      real(kind=8) :: epsl
      !local variables
      real(kind=8) :: d

      d=(r-rc)/delta
      epsl=0.5d0*(epsilon0-1.d0)*(erf(d)+1.d0)+1.d0
    end function epsl

    pure function d1eps(r,rc,delta,epsilon0)
      implicit none
      real(kind=8), intent(in) :: r,rc,delta,epsilon0
      real(kind=8) :: d1eps
      !local variables
      real(kind=8) :: d

      d=(r-rc)/delta
      d1eps=((epsilon0-1.d0)/(delta*sqrt(pi)))*max(exp(-d**2),1.0d-24)
    end function d1eps

end subroutine Eps_rigid_cavity

!> calculates the value of the dielectric function for a smoothed cavity 
!! given a set of centres and radii. Based on the Andreussi epsilon function
!! with a gaussian \rho^{elec}.
!! Need the epsilon0 as well as the radius of the cavit and its smoothness
subroutine Eps_rigid_cavity_new(ndims,nspden,nord,acell,hgrids,nat,rxyz,radii,epsilon0,delta,eps,dlogeps,oneoeps,oneosqrteps,corr)
  implicit none
  integer, intent(in) :: nat !< number of centres defining the cavity
  integer, intent(in) :: nspden
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell
  real(kind=8), intent(in) :: epsilon0 !< dielectric constant of th solvent
  real(kind=8), intent(in) :: delta !< smoothness factor of the cavity
  integer, dimension(3), intent(in) :: ndims   !< dimensions of the simulation box
  real(kind=8), dimension(3), intent(in) :: hgrids !< grid spacings
  real(kind=8), dimension(nat), intent(in) :: radii !< radii of each of the atoms
  !> position of all the atoms in the grid coordinates
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: eps !< dielectric function
  real(kind=8), dimension(3,ndims(1),ndims(2),ndims(3)), intent(out) :: dlogeps !< dlogeps
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneoeps !< inverse of epsilon. Needed for PI method.
  !> inverse square root of epsilon. Needed for PCG method.
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneosqrteps 
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: corr !< correction term of the Generalized Laplacian.
  !local variables
  integer :: i,i1,i2,i3,iat
  real(kind=8) :: r2,x,y2,z2,d,d2,d12,y,z,eps_min,eps1,pi,de2,dde,d1,oneod,h,coeff,dmin,dmax,oneoeps0,oneosqrteps0
  real(kind=8) :: r,t,fact1,fact2,fact3,dd,dtx
  real(kind=8), dimension(3) :: deps,ddeps,v
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)) :: v_de2
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)) :: v_ddeps

   pi = 4.d0*datan(1.d0)
   r=0.d0
   t=0.d0
   oneoeps0=1.d0/epsilon0
   oneosqrteps0=1.d0/dsqrt(epsilon0)

  do i3=1,ndims(3)
   z=hgrids(3)*i3 !(i3-1) for 0 axis start.
   v(3)=z
   do i2=1,ndims(2)
    y=hgrids(2)*i2 !*(i2-1) for 0 axis start.
    v(2)=y
    do i1=1,ndims(1)
     x=hgrids(1)*i1 !(i1-1) for 0 axis start.
     v(1)=x
     !choose the closest atom
     do iat=1,nat
      dmax = radii(iat) - 3.d0*delta
      dmin = radii(iat) + 1.d0*delta
      fact1=2.d0*pi/(-(dmax**2) + dmin**2)
      fact2=(dlog(epsilon0))/(2.d0*pi)
      fact3=(dlog(epsilon0))/(-(dmax**2) + dmin**2)
      d2=(x-rxyz(1,iat))**2+(y-rxyz(2,iat))**2+(z-rxyz(3,iat))**2
      if (d2.eq.0.d0) d2=1.0d-30
      d=dsqrt(d2)

      if (d.lt.dmax) then
       eps(i1,i2,i3)=1.d0
       oneoeps(i1,i2,i3)=1.d0
       oneosqrteps(i1,i2,i3)=1.d0
       do i=1,3
        dlogeps(i,i1,i2,i3)=0.d0
       end do
       corr(i1,i2,i3)=0.d0
      else if (d.gt.dmin) then
       eps(i1,i2,i3)=epsilon0
       oneoeps(i1,i2,i3)=oneoeps0
       oneosqrteps(i1,i2,i3)=oneosqrteps0
       do i=1,3
        dlogeps(i,i1,i2,i3)=0.d0
       end do
       corr(i1,i2,i3)=0.d0
      else
       write(40,*)d,dmin,dmax
       r=fact1*(-(dmax**2) + d2)
       t=fact2*(r-dsin(r))
       eps(i1,i2,i3)=dexp(t)
       oneoeps(i1,i2,i3)=dexp(-t)
       oneosqrteps(i1,i2,i3)=dexp(-0.5d0*t)
!       coeff=fact3*(1.d0-dcos(r))
       dtx=fact3*(1.d0-dcos(r))
!       dtx=-coeff
!       d12=0.d0
       do i=1,3
        dlogeps(i,i1,i2,i3)=dtx*2.d0*(v(i)-rxyz(i,iat))
!        d12 = d12+(dtx*2.d0*(v(i)-rxyz(i,iat)))**2
       end do
       d12 = 4.d0*(dtx**2)*d2
       dd =  4.d0*fact1*fact3*dsin(r)*d2 + 6.d0*dtx
       corr(i1,i2,i3)=(0.125d0/pi)*dexp(t)*(0.5d0*d12+dd)
      end if

     end do
    end do
   end do
  end do

end subroutine Eps_rigid_cavity_new

!> calculates the value of the dielectric function for a smoothed cavity 
!! given a set of centres and radii. Based on the Andreussi paper (Eq. 40) epsilon function
!! with a gaussian \rho^{elec}.
!! Need the epsilon0 as well as the radius of the cavit and its smoothness
subroutine Eps_rigid_cavity_new2(ndims,nspden,nord,acell,hgrids,nat,rxyz,radii,epsilon0,delta,eps,dlogeps,oneoeps,oneosqrteps,corr)
  implicit none
  integer, intent(in) :: nat !< number of centres defining the cavity
  integer, intent(in) :: nspden
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell
  real(kind=8), intent(in) :: epsilon0 !< dielectric constant of th solvent
  real(kind=8), intent(in) :: delta !< smoothness factor of the cavity
  integer, dimension(3), intent(in) :: ndims   !< dimensions of the simulation box
  real(kind=8), dimension(3), intent(in) :: hgrids !< grid spacings
  real(kind=8), dimension(nat), intent(in) :: radii !< radii of each of the atoms
  !> position of all the atoms in the grid coordinates
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: eps !< dielectric function
  real(kind=8), dimension(3,ndims(1),ndims(2),ndims(3)), intent(out) :: dlogeps !< dlogeps
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneoeps !< inverse of epsilon. Needed for PI method.
  !> inverse square root of epsilon. Needed for PCG method.
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneosqrteps 
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: corr !< correction term of the Generalized Laplacian.
  !local variables
  integer :: i,i1,i2,i3,iat
  real(kind=8) :: r2,x,y2,z2,d,d2,d12,y,z,eps_min,eps1,pi,de2,dde,d1,oneod,h,coeff,dmin,dmax,oneoeps0,oneosqrteps0
  real(kind=8) :: r,t,fact1,fact2,fact3,dd,dtx
  real(kind=8), dimension(3) :: deps,ddeps,v
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)) :: v_de2
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)) :: v_ddeps

   pi = 4.d0*datan(1.d0)
   r=0.d0
   t=0.d0
   oneoeps0=1.d0/epsilon0
   oneosqrteps0=1.d0/dsqrt(epsilon0)

  do i3=1,ndims(3)
   z=hgrids(3)*i3 !(i3-1) for 0 axis start.
   v(3)=z
   do i2=1,ndims(2)
    y=hgrids(2)*i2 !*(i2-1) for 0 axis start.
    v(2)=y
    do i1=1,ndims(1)
     x=hgrids(1)*i1 !(i1-1) for 0 axis start.
     v(1)=x
     !choose the closest atom
     do iat=1,nat
      dmax = radii(iat) - 2.d0*delta
      dmin = radii(iat) + 2.d0*delta
      fact1=2.d0*pi/(-(dmax**2) + dmin**2)
      fact2=(epsilon0-1.d0)/(2.d0*pi)
      fact3=(epsilon0-1.d0)/(-(dmax**2) + dmin**2)
      d2=(x-rxyz(1,iat))**2+(y-rxyz(2,iat))**2+(z-rxyz(3,iat))**2
      if (d2.eq.0.d0) d2=1.0d-30
      d=dsqrt(d2)

      if (d.lt.dmax) then
       eps(i1,i2,i3)=1.d0
       oneoeps(i1,i2,i3)=1.d0
       oneosqrteps(i1,i2,i3)=1.d0
       do i=1,3
        dlogeps(i,i1,i2,i3)=0.d0
       end do
       corr(i1,i2,i3)=0.d0
      else if (d.gt.dmin) then
       eps(i1,i2,i3)=epsilon0
       oneoeps(i1,i2,i3)=oneoeps0
       oneosqrteps(i1,i2,i3)=oneosqrteps0
       do i=1,3
        dlogeps(i,i1,i2,i3)=0.d0
       end do
       corr(i1,i2,i3)=0.d0
      else
       write(40,*)d,dmin,dmax
       r=fact1*(-(dmax**2) + d2)
       t=fact2*(r-dsin(r))
       eps(i1,i2,i3)=1.d0 + t
       oneoeps(i1,i2,i3)=1.d0/eps(i1,i2,i3)
       oneosqrteps(i1,i2,i3)=1.d0/dsqrt(eps(i1,i2,i3))
       dtx=fact3*(1.d0-dcos(r))
       do i=1,3
        dlogeps(i,i1,i2,i3)=dtx*2.d0*(v(i)-rxyz(i,iat))/eps(i1,i2,i3)
       end do
       d12 = 4.d0*(dtx**2)*d2
       dd =  4.d0*fact1*fact3*dsin(r)*d2 + 6.d0*dtx
       corr(i1,i2,i3)=(-0.125d0/pi)*(0.5d0*d12/eps(i1,i2,i3)-dd)
      end if

     end do
    end do
   end do
  end do

end subroutine Eps_rigid_cavity_new2

!!$! New value!! Modified nl=14 and nr=15 after else.
!!$subroutine ext_buffers(periodic,nl,nr)
!!$  implicit none
!!$  logical, intent(in) :: periodic
!!$  integer, intent(out) :: nl,nr
!!$
!!$  if (periodic) then
!!$     nl=0
!!$     nr=0
!!$  else
!!$     nl=0
!!$     nr=0
!!$  end if
!!$END SUBROUTINE ext_buffers

!> calculates the value of the dielectric function for a smoothed cavity 
!! given a set of centres and radii. Based on error function.
!! Need the epsilon0 as well as the radius of the cavit and its smoothness
subroutine Eps_rigid_cavity_multiatoms(ndims,nspden,nord,acell,hgrids,natreal,rxyzreal,&
     radiireal,epsilon0,delta,eps,dlogeps,oneoeps,oneosqrteps,corr,geocode,iproc)
  use dynamic_memory
  use yaml_output
  use f_utils
  use f_enums
  use box
  implicit none
  integer, intent(in) :: natreal !< number of centres defining the cavity
  integer, intent(in) :: nspden
  integer, intent(in) :: nord,iproc
  real(kind=8), intent(in) :: acell
  real(kind=8), intent(in) :: epsilon0 !< dielectric constant of th solvent
  real(kind=8), intent(in) :: delta !< smoothness factor of the cavity
  integer, dimension(3), intent(in) :: ndims   !< dimensions of the simulation box
  real(kind=8), dimension(3), intent(in) :: hgrids !< grid spacings
  real(kind=8), dimension(natreal), intent(in) :: radiireal !< radii of each of the atoms
  character(len=2), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  !> position of all the atoms in the grid coordinates
  real(kind=8), dimension(3,natreal), intent(in) :: rxyzreal
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: eps !< dielectric function
  real(kind=8), dimension(3,ndims(1),ndims(2),ndims(3)), intent(out) :: dlogeps !< dlogeps
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneoeps !< inverse of epsilon. Needed for PI method.
  !> inverse square root of epsilon. Needed for PCG method.
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneosqrteps 
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: corr !< correction term of the Generalized Laplacian.
  !local variables
  logical, parameter :: old=.false.
  integer :: i,i1,i2,i3,iat,jat,ii,nat,j,k,l,px,py,pz
  real(kind=8) :: r2,x,y2,z2,d,d2,d12,y,z,eps_min,eps1,pi,de2,dde,d1,oneod,h,coeff,dmin,dmax,oneoeps0,oneosqrteps0
  real(kind=8) :: r,t,fact1,fact2,fact3,dd,dtx,curr,value,valuemin
  real(kind=8), parameter :: valuebc=1.d0
  real(kind=8), dimension(3) :: deps,ddeps,v,rv,shift,sh
  real(kind=8), dimension(6) :: plandist
  integer, dimension(6) :: ba
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)) :: v_de2
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)) :: v_ddeps
  real(kind=8), dimension(:), allocatable :: ep,ddep
  real(kind=8), dimension(:,:), allocatable :: dep
  real(kind=8), dimension(3,27*natreal) :: rxyztot
  real(kind=8), dimension(27*natreal) :: radiitot
  real(kind=8), dimension(:), allocatable :: radii
  real(kind=8), dimension(:,:), allocatable :: rxyz
  logical :: perx,pery,perz
  integer :: nbl1,nbl2,nbl3,nbr1,nbr2,nbr3,imin
  type(cell) :: mesh
  type(box_iterator) :: bit

  pi = 4.d0*datan(1.d0)
  r=0.d0
  t=0.d0
  oneoeps0=1.d0/epsilon0
  oneosqrteps0=1.d0/dsqrt(epsilon0)

  shift(1)=hgrids(1)*ndims(1)
  shift(2)=hgrids(2)*ndims(2)
  shift(3)=hgrids(3)*ndims(3)

  !buffers associated to the geocode
  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

!!$  call ext_buffers(perx,nbl1,nbr1)
!!$  call ext_buffers(pery,nbl2,nbr2)
!!$  call ext_buffers(perz,nbl3,nbr3)
  nbl1=0
  nbl2=0
  nbl3=0
  nbr1=0
  nbr2=0
  nbr3=0


!------------------------------------------------------------------------------------------------------
! Depending of Free, Periodic or Surface bc, image atoms are or not included.

  if (iproc==0) then
   call yaml_map('nbl1',nbl1)
   call yaml_map('nbl2',nbl2)
   call yaml_map('nbl3',nbl3)
   do iat=1,natreal
    call yaml_map('real input atoms',iat)
    call yaml_map('radii',radiireal(iat))
    call yaml_map('rxyz',rxyzreal(:,iat))
   end do
  end if

  px=0 
  py=0
  pz=0
  if (perx) px=1
  if (pery) py=1
  if (perz) pz=1

  rxyztot(:,:)=0.d0

  i=0
  do iat=1,natreal

   ba(1:6)=0
   ! checking what are the image atoms to include in the calculation of the
   ! cavity.
   rv(1:3)=rxyzreal(1:3,iat)
   plandist(1)=dabs(rv(1))
   plandist(2)=dabs(shift(1)-rv(1))
   plandist(3)=dabs(rv(2))
   plandist(4)=dabs(shift(2)-rv(2))
   plandist(5)=dabs(rv(3))
   plandist(6)=dabs(shift(3)-rv(3))
   do ii=1,6
    valuemin=1.d0
    d=plandist(ii)
    value=epsl(d,radiireal(iat),delta)
    if (value.lt.valuebc) then ! valuebc is the value to check on the box border to accept or refuse an image atom.
     if (abs(value).lt.valuemin) then
      valuemin=abs(value)
      imin=ii
     end if
     select case(ii)
     case (1)
      ba(1)=1*px
     case (2)
      ba(2)=-1*px
     case (3)
      ba(3)=1*py
     case (4)
      ba(4)=-1*py
     case (5)
      ba(5)=1*pz
     case (6)
      ba(6)=-1*pz
     end select
    end if
   end do

   do j=ba(6),ba(5)
    sh(3)=real(j,kind=8)*shift(3)
    do k=ba(4),ba(3)
     sh(2)=real(k,kind=8)*shift(2)
     do l=ba(2),ba(1)
      sh(1)=real(l,kind=8)*shift(1)
      rv(1:3)=rxyzreal(1:3,iat) + sh(1:3)
      i=i+1
      rxyztot(1:3,i)=rv(1:3)
      radiitot(i)=radiireal(iat)
     end do
    end do
   end do

  end do

   nat=i

   ep=f_malloc(nat,id='ep')
   ddep=f_malloc(nat,id='ddep')
   dep=f_malloc([3,nat],id='dep')
   rxyz=f_malloc([3,nat],id='rxyz')
   radii=f_malloc(nat,id='radii')

   rxyz(1:3,1:nat)=rxyztot(1:3,1:nat)
   radii(1:nat)=radiitot(1:nat)

   if (iproc==0) then
      !write(*,*)plandist
      !write(*,'(1x,a,1x,e14.7,1x,a,1x,i4)')'Value min =',valuemin,'at bc side',imin
    call yaml_map('nat',nat)
    do iat=1,nat
     call yaml_map('atom',iat)
     call yaml_map('radii',radii(iat))
     call yaml_map('rxyz',rxyz(:,iat))
    end do
 end if

 !------------------------------------------------------------------------------------------------------
 ! Starting the cavity building for rxyztot atoms=real+image atoms (total natcurr) for periodic
 ! and surface boundary conditions or atoms=real for free bc.
 if (old) then

    do i3=1,ndims(3)
       z=hgrids(3)*(i3-1-nbl3)
       !z=hgrids(3)*i3 !(i3-1) for 0 axis start.
       v(3)=z
       do i2=1,ndims(2)
          y=hgrids(2)*(i2-1-nbl2)
          !y=hgrids(2)*i2 !*(i2-1) for 0 axis start.
          v(2)=y
          do i1=1,ndims(1)
             x=hgrids(1)*(i1-1-nbl1)
             !x=hgrids(1)*i1 !(i1-1) for 0 axis start.
             v(1)=x
             !choose the closest atom
             do iat=1,nat
                d2=(x-rxyz(1,iat))**2+(y-rxyz(2,iat))**2+(z-rxyz(3,iat))**2
!print *,'oldi,j,k',i1,i2,i3,d2,iat,v
                d=dsqrt(d2)
                !------------------------------------------------------------------
                ! Loop over the atoms, we neeed to fill the vector
                ! ep(iat) -> atom centered function from 0 (inside the cavity) to 1 outside
                ! dep(i,iat) -> nabla of ep(iat)
                ! ddep(iat) -> laplacian of ep(iat)
                ! Using error function as ep(iat).
                if (d2.eq.0.d0) then
                   d2=1.0d-30
                   ep(iat)=epsl(d,radii(iat),delta)
                   do i=1,3
                      dep(i,iat)=0.d0
                   end do
                   ddep(iat)=2.0d0/(delta**2)*d1eps(0.d0,radii(iat),delta)
                else
                   oneod=1.d0/d
                   ep(iat)=epsl(d,radii(iat),delta)
                   d1=d1eps(d,radii(iat),delta)
                   coeff=2.d0*((sqrt(d2)-radii(iat))/(delta**2))
                   do i=1,3
                      h=(v(i)-rxyz(i,iat))*oneod
                      dep(i,iat) =d1*h
                   end do
                   ddep(iat)=d1*(2.d0*oneod-coeff)
                end if

                ! Using gaussian function as ep(iat).
                !      ep(iat)=epgauss(d2,radii(iat),delta)
                !      d1=d1epgauss(d2,radii(iat),delta)
                !      do i=1,3
                !        h=v(i)-rxyz(i,iat)
                !       dep(i,iat)=d1*h
                !      end do
                !      ddep(iat)=d1*(3.d0-d2/(delta**2))

                !------------------------------------------------------------------

             end do

             eps(i1,i2,i3)=(epsilon0-1.d0)*product(ep)+1.d0
             oneoeps(i1,i2,i3)=1.d0/eps(i1,i2,i3)
             oneosqrteps(i1,i2,i3)=1.d0/dsqrt(eps(i1,i2,i3))

             do i=1,3
                deps(i)=0.d0
                do jat=0,nat-1
                   curr=dep(i,jat+1)
                   do iat=1,nat-1
                      curr=curr*ep(modulo(iat+jat,nat)+1)
                   end do
                   deps(i) = deps(i) + curr
                end do
                deps(i) = deps(i)*(epsilon0-1.d0)
             end do

             d12=0.d0
             do i=1,3
                dlogeps(i,i1,i2,i3)=deps(i)/eps(i1,i2,i3)
                d12 = d12 + deps(i)**2
             end do

             dd=0.d0
             do jat=1,nat
                curr=ddep(jat)
                do iat=1,nat-1
                   curr=curr*ep(modulo(iat+jat-1,nat)+1)
                end do
                dd = dd + curr
             end do

             do i=1,3
                do iat=1,nat-1
                   do jat=iat+1,nat
                      curr=dep(i,iat)*dep(i,jat)
                      do ii=1,nat
                         if ((ii.eq.iat).or.(ii.eq.jat)) then
                         else
                            curr=curr*ep(ii)
                         end if
                      end do
                      curr=curr*2.d0
                      dd = dd + curr
                   end do
                end do
             end do

             dd=dd*(epsilon0-1.d0)
             corr(i1,i2,i3)=(-0.125d0/pi)*(0.5d0*d12/eps(i1,i2,i3)-dd)

          end do
       end do
    end do
 else

    mesh=cell_new(geocode,ndims,hgrids) 
    bit=box_iter(mesh)
    do while(box_next_point(bit))
       !choose the closest atom
       do iat=1,nat
          d2=square(mesh,bit%rxyz-rxyz(:,iat))
!print *,'i,j,k',bit%i,bit%j,bit%k,d2,iat,bit%rxyz
          d=dsqrt(d2)
          !------------------------------------------------------------------
          ! Loop over the atoms, we neeed to fill the vector
          ! ep(iat) -> atom centered function from 0 (inside the cavity) to 1 outside
          ! dep(i,iat) -> nabla of ep(iat)
          ! ddep(iat) -> laplacian of ep(iat)
          ! Using error function as ep(iat).
          if (d2.eq.0.d0) then
             d2=1.0d-30
             ep(iat)=epsl(d,radii(iat),delta)
             do i=1,3
                dep(i,iat)=0.d0
             end do
             ddep(iat)=2.0d0/(delta**2)*d1eps(0.d0,radii(iat),delta)
          else
             oneod=1.d0/d
             ep(iat)=epsl(d,radii(iat),delta)
             d1=d1eps(d,radii(iat),delta)
             coeff=2.d0*((d-radii(iat))/(delta**2))            
             do i=1,3
                h=(bit%rxyz(i)-rxyz(i,iat))*oneod
                dep(i,iat) =d1*h
             end do
             ddep(iat)=d1*(2.d0*oneod-coeff)
          end if

          ! Using gaussian function as ep(iat).
          !      ep(iat)=epgauss(d2,radii(iat),delta)
          !      d1=d1epgauss(d2,radii(iat),delta)
          !      do i=1,3
          !        h=v(i)-rxyz(i,iat)
          !       dep(i,iat)=d1*h
          !      end do
          !      ddep(iat)=d1*(3.d0-d2/(delta**2))

          !------------------------------------------------------------------

       end do
       eps(bit%i,bit%j,bit%k)=(epsilon0-1.d0)*product(ep)+1.d0
       oneoeps(bit%i,bit%j,bit%k)=1.d0/eps(bit%i,bit%j,bit%k)
       oneosqrteps(bit%i,bit%j,bit%k)=1.d0/dsqrt(eps(bit%i,bit%j,bit%k))

       do i=1,3
          deps(i)=0.d0
          do jat=0,nat-1
             curr=dep(i,jat+1)
             do iat=1,nat-1
                curr=curr*ep(modulo(iat+jat,nat)+1)
             end do
             deps(i) = deps(i) + curr
          end do
          deps(i) = deps(i)*(epsilon0-1.d0)
       end do

!       d12=0.d0
       do i=1,3
          dlogeps(i,bit%i,bit%j,bit%k)=deps(i)/eps(bit%i,bit%j,bit%k)
!          d12 = d12 + deps(i)**2
       end do

       d12=square(mesh,deps)

       dd=0.d0
       do jat=1,nat
          curr=ddep(jat)
          do iat=1,nat-1
             curr=curr*ep(modulo(iat+jat-1,nat)+1)
          end do
          dd = dd + curr
       end do

       do i=1,3
          do iat=1,nat-1
             do jat=iat+1,nat
                curr=dep(i,iat)*dep(i,jat)
                do ii=1,nat
                   if ((ii.eq.iat).or.(ii.eq.jat)) then
                   else
                      curr=curr*ep(ii)
                   end if
                end do
                curr=curr*2.d0
                dd = dd + curr
             end do
          end do
       end do

       dd=dd*(epsilon0-1.d0)
       corr(bit%i,bit%j,bit%k)=(-0.125d0/pi)*(0.5d0*d12/eps(bit%i,bit%j,bit%k)-dd)

    end do

 end if
 call f_free(ep)
 call f_free(ddep)
  call f_free(dep)
  call f_free(rxyz)
  call f_free(radii)

  contains

    pure function epsl(r,rc,delta)
      implicit none
      real(kind=8), intent(in) :: r,rc,delta
      real(kind=8) :: epsl
      !local variables
      real(kind=8) :: d

      d=(r-rc)/delta
      epsl=0.5d0*(erf(d)+1.d0)
    end function epsl

    pure function d1eps(r,rc,delta)
      use numerics, only: safe_exp
      implicit none
      real(kind=8), intent(in) :: r,rc,delta
      real(kind=8) :: d1eps
      !local variables
      real(kind=8) :: d

      d=(r-rc)/delta
      d1eps=(1.d0/(delta*sqrt(pi)))*max(safe_exp(-d**2),1.0d-24) !max function IMPORTANT for
      !big r (about 10) to avoid oscillations far away.
      !d1eps=(1.d0/(delta*sqrt(pi)))*safe_exp(-d**2)
    end function d1eps
    
    pure function epgauss(r2,rc,delta)
      use numerics, only: safe_exp
      implicit none
      real(kind=8), intent(in) :: r2,rc,delta
      real(kind=8) :: epgauss
      !local variables
      real(kind=8) :: d

      d=r2/(2.d0*delta**2)
      epgauss=1.d0 - safe_exp(-d)
    end function epgauss
    
    pure function d1epgauss(r2,rc,delta)
      use numerics, only: safe_exp
      implicit none
      real(kind=8), intent(in) :: r2,rc,delta
      real(kind=8) :: d1epgauss
      !local variables
      real(kind=8) :: d

      d=r2/(2.d0*delta**2)
      d1epgauss=safe_exp(-d)/(delta**2)
    end function d1epgauss

end subroutine Eps_rigid_cavity_multiatoms

!> calculates the value of the dielectric function for a smoothed cavity 
!! given a set of centres and radii. Based on the Andreussi epsilon function
!! with a gaussian \rho^{elec}.
!! Need the epsilon0 as well as the radius of the cavit and its smoothness
subroutine Eps_rigid_cavity_new_multiatoms(ndims,nspden,nord,acell,hgrids,nat,&
     rxyz,radii,epsilon0,delta,eps,dlogeps,oneoeps,oneosqrteps,corr)
  implicit none
  integer, intent(in) :: nat !< number of centres defining the cavity
  integer, intent(in) :: nspden
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell
  real(kind=8), intent(in) :: epsilon0 !< dielectric constant of th solvent
  real(kind=8), intent(in) :: delta !< smoothness factor of the cavity
  integer, dimension(3), intent(in) :: ndims   !< dimensions of the simulation box
  real(kind=8), dimension(3), intent(in) :: hgrids !< grid spacings
  real(kind=8), dimension(nat), intent(in) :: radii !< radii of each of the atoms
  !> position of all the atoms in the grid coordinates
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: eps !< dielectric function
  real(kind=8), dimension(3,ndims(1),ndims(2),ndims(3)), intent(out) :: dlogeps !< dlogeps
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneoeps !< inverse of epsilon. Needed for PI method.
  !> inverse square root of epsilon. Needed for PCG method.
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneosqrteps 
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: corr !< correction term of the Generalized Laplacian.
  !local variables
  integer :: i,i1,i2,i3,iat,jat,ii
  real(kind=8) :: r2,x,y2,z2,d,d2,d12,y,z,eps_min,eps1,pi,de2,dde,d1,oneod,h,coeff,dmin,dmax,oneoeps0,oneosqrteps0
  real(kind=8) :: r,t,fact1,fact2,fact3,dd,dtx,curr
  real(kind=8), dimension(3) :: deps,ddeps,v
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)) :: v_de2
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)) :: v_ddeps
  real(kind=8), dimension(nat) :: ep,ddep
  real(kind=8), dimension(3,nat) :: dep

   pi = 4.d0*datan(1.d0)
   r=0.d0
   t=0.d0
   oneoeps0=1.d0/epsilon0
   oneosqrteps0=1.d0/dsqrt(epsilon0)

  do i3=1,ndims(3)
   z=hgrids(3)*i3 !(i3-1) for 0 axis start.
   v(3)=z
   do i2=1,ndims(2)
    y=hgrids(2)*i2 !*(i2-1) for 0 axis start.
    v(2)=y
    do i1=1,ndims(1)
     x=hgrids(1)*i1 !(i1-1) for 0 axis start.
     v(1)=x
     !choose the closest atom
     do iat=1,nat
      dmax = radii(iat) - 2.36d0*delta
      dmin = radii(iat) + 1.64d0*delta
      fact1=2.d0*pi/(-(dmax**2) + dmin**2)
      fact2=(dlog(2.d0))/(2.d0*pi)
      fact3=(dlog(2.d0))/(-(dmax**2) + dmin**2)
      d2=(x-rxyz(1,iat))**2+(y-rxyz(2,iat))**2+(z-rxyz(3,iat))**2
      if (d2.eq.0.d0) d2=1.0d-30
      d=dsqrt(d2)

      if (d.lt.dmax) then
       ep(iat)=0.d0
       do i=1,3
        dep(i,iat)=0.d0
       end do
       ddep(iat)=0.d0
      else if (d.gt.dmin) then
       ep(iat)=1.d0
       do i=1,3
        dep(i,iat)=0.d0
       end do
       ddep(iat)=0.d0
      else
       r=fact1*(-(dmax**2) + d2)
       t=fact2*(r-dsin(r))
       ep(iat)=dexp(t)-1.d0
       dtx=fact3*(1.d0-dcos(r))
       do i=1,3
        dep(i,iat)=dexp(t)*dtx*2.d0*(v(i)-rxyz(i,iat))
       end do
       ddep(iat) = dexp(t)*(4.d0*(dtx**2)*d2 + 4.d0*fact1*fact3*dsin(r)*d2 + 6.d0*dtx)
      end if
     end do

     eps(i1,i2,i3)=(epsilon0-1.d0)*product(ep)+1.d0
     oneoeps(i1,i2,i3)=1.d0/eps(i1,i2,i3)
     oneosqrteps(i1,i2,i3)=1.d0/dsqrt(eps(i1,i2,i3))

     do i=1,3
      deps(i)=0.d0
      do jat=0,nat-1
       curr=dep(i,jat+1)
       do iat=1,nat-1
        curr=curr*ep(modulo(iat+jat,nat)+1)
       end do
        deps(i) = deps(i) + curr
      end do
      deps(i) = deps(i)*(epsilon0-1.d0)
     end do

     d12=0.d0
     do i=1,3
      dlogeps(i,i1,i2,i3)=deps(i)/eps(i1,i2,i3)
      d12 = d12 + deps(i)**2
     end do

     dd=0.d0
     do jat=1,nat
      curr=ddep(jat)
      do iat=1,nat-1
       curr=curr*ep(modulo(iat+jat-1,nat)+1)
      end do
      dd = dd + curr
     end do

      do i=1,3
       do iat=1,nat-1
        do jat=iat+1,nat
         curr=dep(i,iat)*dep(i,jat)
         do ii=1,nat
          if ((ii.eq.iat).or.(ii.eq.jat)) then
          else
           curr=curr*ep(ii)
          end if
         end do
         curr=curr*2.d0
         dd = dd + curr
        end do
       end do
      end do

     dd=dd*(epsilon0-1.d0)
     corr(i1,i2,i3)=(-0.125d0/pi)*(0.5d0*d12/eps(i1,i2,i3)-dd)

    end do
   end do
  end do

end subroutine Eps_rigid_cavity_new_multiatoms

!> calculates the value of the dielectric function for a smoothed cavity 
!! given a set of centres and radii. Based on the Andreussi paper (Eq. 40) epsilon function
!! with a gaussian \rho^{elec}.
!! Need the epsilon0 as well as the radius of the cavit and its smoothness
subroutine Eps_rigid_cavity_new2_multiatoms(ndims,nspden,nord,acell,hgrids,nat,&
     rxyz,radii,epsilon0,delta,eps,dlogeps,oneoeps,oneosqrteps,corr)
  implicit none
  integer, intent(in) :: nat !< number of centres defining the cavity
  integer, intent(in) :: nspden
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell
  real(kind=8), intent(in) :: epsilon0 !< dielectric constant of th solvent
  real(kind=8), intent(in) :: delta !< smoothness factor of the cavity
  integer, dimension(3), intent(in) :: ndims   !< dimensions of the simulation box
  real(kind=8), dimension(3), intent(in) :: hgrids !< grid spacings
  real(kind=8), dimension(nat), intent(in) :: radii !< radii of each of the atoms
  !> position of all the atoms in the grid coordinates
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: eps !< dielectric function
  real(kind=8), dimension(3,ndims(1),ndims(2),ndims(3)), intent(out) :: dlogeps !< dlogeps
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneoeps !< inverse of epsilon. Needed for PI method.
  !> inverse square root of epsilon. Needed for PCG method.
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneosqrteps
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: corr !< correction term of the Generalized Laplacian.
  !local variables
  integer :: i,i1,i2,i3,iat,jat,ii
  real(kind=8) :: r2,x,y2,z2,d,d2,d12,y,z,eps_min,eps1,pi,de2,dde,d1,oneod,h,coeff,dmin,dmax,oneoeps0,oneosqrteps0
  real(kind=8) :: r,t,fact1,fact2,fact3,dd,dtx,curr
  real(kind=8), dimension(3) :: deps,ddeps,v
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)) :: v_de2
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)) :: v_ddeps
  real(kind=8), dimension(nat) :: ep,ddep
  real(kind=8), dimension(3,nat) :: dep

   pi = 4.d0*datan(1.d0)
   r=0.d0
   t=0.d0
   oneoeps0=1.d0/epsilon0
   oneosqrteps0=1.d0/dsqrt(epsilon0)

  do i3=1,ndims(3)
   z=hgrids(3)*i3 !(i3-1) for 0 axis start.
   v(3)=z
   do i2=1,ndims(2)
    y=hgrids(2)*i2 !*(i2-1) for 0 axis start.
    v(2)=y
    do i1=1,ndims(1)
     x=hgrids(1)*i1 !(i1-1) for 0 axis start.
     v(1)=x
     !choose the closest atom
     do iat=1,nat
      dmax = radii(iat) - 2.d0*delta
      dmin = radii(iat) + 2.d0*delta
      fact1=2.d0*pi/(-(dmax**2) + dmin**2)
      fact2=1.d0/(2.d0*pi)
      fact3=1.d0/(-(dmax**2) + dmin**2)
      d2=(x-rxyz(1,iat))**2+(y-rxyz(2,iat))**2+(z-rxyz(3,iat))**2
      if (d2.eq.0.d0) d2=1.0d-30
      d=dsqrt(d2)

      if (d.lt.dmax) then
       ep(iat)=0.d0
       do i=1,3
        dep(i,iat)=0.d0
       end do
       ddep(iat)=0.d0
      else if (d.gt.dmin) then
       ep(iat)=1.d0
       do i=1,3
        dep(i,iat)=0.d0
       end do
       ddep(iat)=0.d0
      else
       write(40,*)d,dmin,dmax
       r=fact1*(-(dmax**2) + d2)
       t=fact2*(r-dsin(r))
       ep(iat)=t
       dtx=fact3*(1.d0-dcos(r))
       do i=1,3
        dep(i,iat)=dtx*2.d0*(v(i)-rxyz(i,iat))
       end do
       ddep(iat) =  4.d0*fact1*fact3*dsin(r)*d2 + 6.d0*dtx
      end if
     end do

     eps(i1,i2,i3)=(epsilon0-1.d0)*product(ep)+1.d0
     oneoeps(i1,i2,i3)=1.d0/eps(i1,i2,i3)
     oneosqrteps(i1,i2,i3)=1.d0/dsqrt(eps(i1,i2,i3))

     do i=1,3
      deps(i)=0.d0
      do jat=0,nat-1
       curr=dep(i,jat+1)
       do iat=1,nat-1
        curr=curr*ep(modulo(iat+jat,nat)+1)
       end do
        deps(i) = deps(i) + curr
      end do
      deps(i) = deps(i)*(epsilon0-1.d0)
     end do

     d12=0.d0
     do i=1,3
      dlogeps(i,i1,i2,i3)=deps(i)/eps(i1,i2,i3)
      d12 = d12 + deps(i)**2
     end do

     dd=0.d0
     do jat=1,nat
      curr=ddep(jat)
      do iat=1,nat-1
       curr=curr*ep(modulo(iat+jat-1,nat)+1)
      end do
      dd = dd + curr
     end do

      do i=1,3
       do iat=1,nat-1
        do jat=iat+1,nat
         curr=dep(i,iat)*dep(i,jat)
         do ii=1,nat
!          if ((ii.ne.iat).or.(ii.ne.jat)) then
          if ((ii.eq.iat).or.(ii.eq.jat)) then
          else
           curr=curr*ep(ii)
          end if
         end do
         curr=curr*2.d0
         dd = dd + curr
        end do
       end do
      end do

     dd=dd*(epsilon0-1.d0)
     corr(i1,i2,i3)=(-0.125d0/pi)*(0.5d0*d12/eps(i1,i2,i3)-dd)

    end do
   end do
  end do

end subroutine Eps_rigid_cavity_new2_multiatoms

subroutine SetEledens(n01,n02,n03,nspden,nord,acell,a_gauss,hx,hy,hz,SetEps,edens,nabla,ddt)

  use yaml_output

  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  integer, intent(in) :: nord
  real(kind=8), intent(in) :: acell,a_gauss,hx,hy,hz
  integer, intent(in) :: SetEps
  real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: edens
  real(kind=8), dimension(n01,n02,n03,nspden,3), intent(out) :: nabla
  real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: ddt

  real(kind=8) :: x1,x2,x3,pi,factor,r2,fac,sigma,x0,sume,dd
  real(kind=8), dimension(3) :: x
  integer :: i,i1,i2,i3,isp

  open(unit=18,file='ElectronDensity.out',status='unknown')
  open(unit=19,file='ElectronDensity_line.out',status='unknown')

  if (SetEps.eq.2) then

   pi = 4.d0*atan(1.d0)
   sigma = 0.04d0*acell
   x0 = 0.d0 ! hx*real(25-n01/2,kind=8)

         !Normalization
!         factor = 1.d0/((sigma**3)*sqrt((2.d0*pi)**3))
         factor = 1.d0

         !gaussian function for the potential.
         sume=0.d0
         isp=1
         do i3=1,n03
            x3 = hz*real(i3-n03/2,kind=8)
            x(3) = x3
            do i2=1,n02
               x2 = hy*real(i2-n02/2,kind=8)
               x(2) = x2
               do i1=1,n01
                  x1 = hx*real(i1-n01/2,kind=8)
                  x(1) = x1
!                  r2=(x1-x0)*(x1-x0)+(x2-x0)*(x2-x0)+(x3)*(x3)
                  r2 = x1*x1+x2*x2+x3*x3
                  edens(i1,i2,i3,isp) = max(factor*exp(-0.5d0*r2/(sigma**2)),1d-24)
                  sume = sume + edens(i1,i2,i3,isp)
                  fac=-edens(i1,i2,i3,isp)/(sigma**2)
                  do i=1,3
                   nabla(i1,i2,i3,isp,i)= fac*x(i)
                  end do      
                  ddt(i1,i2,i3,isp)= fac*(3.d0-(r2/(sigma**2)))
               end do
            end do
          end do

          sume = sume*hx*hy*hz

  else
   write(*,*)'Error in the Set Eledens subroutine'
   stop
  end if

    call yaml_map('Total electron Charge',sume)

     i3=1!n03/2
     do i2=1,n02
      do i1=1,n01
       write(18,'(2(1x,I4),2(1x,e14.7))')i1,i2,edens(i1,i2,i3,isp),edens(i1,i2,n03/2,isp)
      end do
      write(18,*)
     end do

     do i1=1,n01
      write(19,'(1x,I8,2(1x,e22.15))') i1,edens(i1,n02/2,n03/2,isp),edens(n01/2,i1,n03/2,isp)
     end do

  close(unit=18)
  close(unit=19)

end subroutine SetEledens

subroutine get_size_from_cube(n01,n02,n03,hx,hy,hz,nat)

  use yaml_output

  implicit none
  integer, intent(out) :: n01
  integer, intent(out) :: n02
  integer, intent(out) :: n03
  real(kind=8), intent(out) :: hx
  real(kind=8), intent(out) :: hy
  real(kind=8), intent(out) :: hz
  integer, intent(out) :: nat

  real(8):: z_incr, z1_incr,x1_incr, x2_incr, y_incr,y2_incr
  real(8):: origin_x, origin_y, origin_z
  character (100) :: line1,line2

  open (52,file='electronic_density.cube')

  read(52,*)line1
  read(52,*)line2
  read(52,*)nat,origin_x,origin_y,origin_z
  read(52,*)n01, hx, x1_incr, x2_incr
  read(52,*)n02, y_incr, hy, y2_incr
  read(52,*)n03, z_incr, z1_incr, hz

  close(unit=52)

end subroutine get_size_from_cube

subroutine get_rho(n01,n02,n03,nspden,natoms,acell,a_gauss,hx,hy,hz,rhoele,rhoion,sume,rxyz,iproc)

  use yaml_output
  use f_utils

  implicit none
  integer, intent(in) :: n01
  integer, intent(in) :: n02
  integer, intent(in) :: n03
  integer, intent(in) :: nspden
  integer, intent(in) :: natoms,iproc
  real(kind=8), intent(in) :: acell,a_gauss,hx,hy,hz
  real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: rhoele
  real(kind=8), dimension(n01,n02,n03,nspden), intent(out) :: rhoion
  real(kind=8), intent(out) :: sume
  real(kind=8), dimension(3,natoms), intent(out) :: rxyz

  real(kind=8) :: x1,x2,x3,pi,factor,r2,fac,sigma,x0,sumi
  real(kind=8), dimension(3) :: vec
  integer :: i,isp

  integer:: ngp,igp,I1,I2,I3,grid_x, grid_y, grid_z
  integer:: N_line,ilmp,iat,nat,no_line,unt
  real(8), allocatable:: ATOM_info(:,:), chg(:,:,:), sigat(:)
  real(8):: x,y,z,z_incr, z1_incr,x1_incr, x2_incr, y_incr,y2_incr
  real(8):: origin_x, origin_y, origin_z
  real(8):: total,xcent,ycent,zcent,hgridx_h,hgridy_h,hgridz_h,charge
  character(14) filename
  character(3) fn
  character (100) :: line1,line2

  open (50,file='electronic_density.cube')
  open (51,file='sigat.dat')

  read(50,*)line1
  read(50,*)line2
  read(50,*)nat,origin_x,origin_y,origin_z
  read(50,*)grid_x, hgridx_h, x1_incr, x2_incr
  read(50,*)grid_y, y_incr, hgridy_h, y2_incr
  read(50,*)grid_z, z_incr, z1_incr, hgridz_h

  if ((n01.ne.grid_x).or.(n02.ne.grid_y).or.(n03.ne.grid_z).or.(nat.ne.natoms)) then
   write(*,*)'Error get_rho subroutine'
   stop
  end if

  N_line=NINT(((grid_y*grid_x*grid_z)/6)*1.d0)
!  write(*,*)N_line
  no_line=((grid_y*grid_x*grid_z))

  allocate(chg(0:grid_x-1,0:grid_y-1,0:grid_z-1),ATOM_info(5,nat),sigat(nat))

  do I1= 1,nat
   read(50,*) (ATOM_info(I3,I1),I3=1,5)
   read(51,*) sigat(I1)
  end do

  do I1=0,grid_x-1
   do I2=0,grid_y-1
    read(50,*)(chg(I1,I2,I3),I3=0,grid_z-1)
   end do
  end do

!  hgridx=.4d0
!  hgridy=.4d0
!  hgridz=.4d0

!  hgridx_h=.5d0*hgridx
!  hgridy_h=.5d0*hgridy
!  hgridz_h=.5d0*hgridz

  total=0.d0
  xcent=0.d0
  ycent=0.d0
  zcent=0.d0
!  write(*,*)hgridz_h,hgridy_h,hgridx_h

  do I3=0,grid_z-1
   z=(I3)*hgridz_h
   do I2=0,grid_y-1
    y=(I2)*hgridy_h
    do I1=0,grid_x-1
     x=(I1)*hgridx_h

     charge=chg(I1,I2,I3)
     total=total+charge
     xcent=xcent+charge*x
     ycent=ycent+charge*y
     zcent=zcent+charge*z

    end do
   end do
  end do

  total=-total*hgridx_h*hgridy_h*hgridz_h
  xcent=-xcent*hgridx_h*hgridy_h*hgridz_h
  ycent=-ycent*hgridx_h*hgridy_h*hgridz_h
  zcent=-zcent*hgridx_h*hgridy_h*hgridz_h
  sume = total

  if (iproc ==0) then
   write(*,*) 'total, centers  Without ions'
   write(*,*)  total,xcent,ycent,zcent
  end if

  do iat= 1,nat
   if (iproc ==0) then
    write(*,*) ATOM_info(:,iat)
   end if
   total=total+ATOM_info(2,iat)
   xcent=xcent+ATOM_info(2,iat)*ATOM_info(3,iat)
   ycent=ycent+ATOM_info(2,iat)*ATOM_info(4,iat)
   zcent=zcent+ATOM_info(2,iat)*ATOM_info(5,iat)
 end do

   if (iproc ==0) then
    write(*,*) 'total, centers  with    ions'
    write(*,*)  total,xcent,ycent,zcent
   end if

!-----------------------------------------------------------------------
! Compute ion charge density

  pi = 4.d0*atan(1.d0)

  isp=1
  do I3=0,grid_z-1
   do I2=0,grid_y-1
    do I1=0,grid_x-1
     rhoele(I1+1,I2+1,I3+1,isp) = chg(I1,I2,I3)
     rhoion(I1+1,I2+1,I3+1,isp) = 0.d0
    end do
   end do
  end do

  do iat= 1,nat

   do i=1,3
    vec(i) = ATOM_info(i+2,iat)
    rxyz(i,iat) = ATOM_info(i+2,iat)
   end do

   sigma = sigat(iat) !*acell
   x0 = 0.d0 ! hx*real(25-n01/2,kind=8)

   !Normalization
   factor = ATOM_info(2,iat)/((sigma**3)*sqrt((2.d0*pi)**3))

   !gaussian function for the potential.
   sumi=0.d0
   isp=1

  do I3=0,grid_z-1
   x3=(I3)*hgridz_h
   do I2=0,grid_y-1
    x2=(I2)*hgridy_h
    do I1=0,grid_x-1
     x1=(I1)*hgridx_h
      r2=(x1-vec(1))**2+(x2-vec(2))**2+(x3-vec(3))**2
!                  r2 = x1*x1+x2*x2+x3*x3
      rhoion(I1+1,I2+1,I3+1,isp) = rhoion(I1+1,I2+1,I3+1,isp)  +  max(factor*exp(-0.5d0*r2/(sigma**2)),1d-24)
      sumi = sumi + rhoion(I1+1,I2+1,I3+1,isp)
     end do
    end do
   end do

  end do

     unt=f_get_free_unit(40)
     call f_open_file(unt,file='ElectronDensity_fromc.dat')
     i3=1!n03/2
     do i2=1,n02
      do i1=1,n01
       write(unt,'(2(1x,I4),2(1x,e14.7))')i1,i2,rhoele(i1,i2,i3,isp),rhoele(i1,i2,n03/2,isp)
      end do
      write(unt,*)
     end do
     call f_close(unt)

     unt=f_get_free_unit(41)
     call f_open_file(unt,file='ElectronDensity_fromc_line.dat')
     do i2=1,n02
      write(unt,'(1x,I8,3(1x,e22.15))')i2,rhoele(i2,n02/2,n03/2,isp),rhoele(n01/2,i2,n03/2,isp),rhoele(n01/2,n02/2,i2,isp)
     end do
     call f_close(unt)

     unt=f_get_free_unit(42)
     call f_open_file(unt,file='IonicDensity.dat')
     i3=1!n03/2
     do i2=1,n02
      do i1=1,n01
       write(unt,'(2(1x,I4),2(1x,e14.7))')i1,i2,rhoion(i1,i2,i3,isp),rhoion(i1,i2,n03/2,isp)
      end do
      write(unt,*)
     end do
     call f_close(unt)

     unt=f_get_free_unit(43)
     call f_open_file(unt,file='IonicDensity_line.dat')
     do i2=1,n02
      write(unt,'(1x,I8,1(1x,e22.15))')i2,rhoion(n01/2,i2,n03/2,isp)
     end do
     call f_close(unt)


   if (iproc ==0) then
    write(*,*)'Total ions charge'
    write(*,*)sumi*hgridx_h*hgridy_h*hgridz_h
   end if

  deallocate(chg,ATOM_info,sigat)

  close(unit=50)
  close(unit=51)

end subroutine get_rho

subroutine derf_local(derf_yy,yy)

 implicit none
 integer, parameter :: dp=8
 real(dp),intent(in) :: yy
 real(dp),intent(out) :: derf_yy
 integer          ::  done,ii,isw
 real(dp), parameter :: &
       ! coefficients for 0.0 <= yy < .477
       &  pp(5)=(/ 113.8641541510502e0_dp, 377.4852376853020e0_dp,  &
       &           3209.377589138469e0_dp, .1857777061846032e0_dp,  &
       &           3.161123743870566e0_dp /)
  real(dp), parameter :: &
       &  qq(4)=(/ 244.0246379344442e0_dp, 1282.616526077372e0_dp,  &
       &           2844.236833439171e0_dp, 23.60129095234412e0_dp/)
  ! coefficients for .477 <= yy <= 4.0
  real(dp), parameter :: &
       &  p1(9)=(/ 8.883149794388376e0_dp, 66.11919063714163e0_dp,  &
       &           298.6351381974001e0_dp, 881.9522212417691e0_dp,  &
       &           1712.047612634071e0_dp, 2051.078377826071e0_dp,  &
       &           1230.339354797997e0_dp, 2.153115354744038e-8_dp, &
       &           .5641884969886701e0_dp /)
  real(dp), parameter :: &
       &  q1(8)=(/ 117.6939508913125e0_dp, 537.1811018620099e0_dp,  &
       &           1621.389574566690e0_dp, 3290.799235733460e0_dp,  &
       &           4362.619090143247e0_dp, 3439.367674143722e0_dp,  &
       &           1230.339354803749e0_dp, 15.74492611070983e0_dp/)
  ! coefficients for 4.0 < y,
  real(dp), parameter :: &
       &  p2(6)=(/ -3.603448999498044e-01_dp, -1.257817261112292e-01_dp,   &
       &           -1.608378514874228e-02_dp, -6.587491615298378e-04_dp,   &
       &           -1.631538713730210e-02_dp, -3.053266349612323e-01_dp/)
  real(dp), parameter :: &
       &  q2(5)=(/ 1.872952849923460e0_dp   , 5.279051029514284e-01_dp,    &
       &           6.051834131244132e-02_dp , 2.335204976268692e-03_dp,    &
       &           2.568520192289822e0_dp /)
  real(dp), parameter :: &
       &  sqrpi=.5641895835477563e0_dp, xbig=13.3e0_dp, xlarge=6.375e0_dp, xmin=1.0e-10_dp
  real(dp) ::  res,xden,xi,xnum,xsq,xx

 xx = yy
 isw = 1
!Here change the sign of xx, and keep track of it thanks to isw
 if (xx<0.0e0_dp) then
  isw = -1
  xx = -xx
 end if

 done=0

!Residual value, if yy < -6.375e0_dp
 res=-1.0e0_dp

!abs(yy) < .477, evaluate approximation for erfc
 if (xx<0.477e0_dp) then
! xmin is a very small number
  if (xx<xmin) then
   res = xx*pp(3)/qq(3)
  else
   xsq = xx*xx
   xnum = pp(4)*xsq+pp(5)
   xden = xsq+qq(4)
   do ii = 1,3
    xnum = xnum*xsq+pp(ii)
    xden = xden*xsq+qq(ii)
   end do
   res = xx*xnum/xden
  end if
  if (isw==-1) res = -res
  done=1
 end if

!.477 < abs(yy) < 4.0 , evaluate approximation for erfc
!.477 < abs(yy) < 4.0 , evaluate approximation for erfc
 if (xx<=4.0e0_dp .and. done==0 ) then
  xsq = xx*xx
  xnum = p1(8)*xx+p1(9)
  xden = xx+q1(8)
  do ii=1,7
   xnum = xnum*xx+p1(ii)
   xden = xden*xx+q1(ii)
  end do
  res = xnum/xden
  res = res* exp(-xsq)
  if (isw.eq.-1) then
     res = res-1.0e0_dp
  else
     res=1.0e0_dp-res
  end if
  done=1
 end if

!y > 13.3e0_dp
 if (isw > 0 .and. xx > xbig .and. done==0 ) then
  res = 1.0e0_dp
  done=1
 end if

!4.0 < yy < 13.3e0_dp  .or. -6.375e0_dp < yy < -4.0
!evaluate minimax approximation for erfc
 if ( ( isw > 0 .or. xx < xlarge ) .and. done==0 ) then
  xsq = xx*xx
  xi = 1.0e0_dp/xsq
  xnum= p2(5)*xi+p2(6)
  xden = xi+q2(5)
  do ii = 1,4
   xnum = xnum*xi+p2(ii)
   xden = xden*xi+q2(ii)
  end do
  res = (sqrpi+xi*xnum/xden)/xx
  res = res* exp(-xsq)
  if (isw.eq.-1) then
     res = res-1.0e0_dp
  else
     res=1.0e0_dp-res
  end if
 end if

!All cases have been investigated
 derf_yy = res

end subroutine derf_local
