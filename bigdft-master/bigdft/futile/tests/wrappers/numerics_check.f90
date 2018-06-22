!> @file
!!  Test of some functionalities of the numeric groups
!! @author
!!    Copyright (C) 2016-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
program numeric_check
  use futile
  use f_harmonics
  character(len=*), parameter :: input1=&
       "  {name: ndim, shortname: n, default: 30,"//&
       "  help_string: Size of the array for multipoles,"//&
       "  help_dict: {Allowed values: integer}}"
  character(len=*), parameter :: input2=&
       "  {name: boldify, shortname: b, default: None,"//&
       "  help_string: Boldify the string as a test,"//&
       "  help_dict: {Allowed values: string scalar}}"
  character(len=*), parameter :: input3=&
       "  {name: blinkify, shortname: l, default: None,"//&
       "  help_string: Make the string blinking,"//&
       "  help_dict: {Allowed values: string scalar}}"

  character(len=*), parameter :: inputs=&
       '-'//input1//f_cr//&
       '-'//input2//f_cr//&
       '-'//input3


  integer :: n,i
  type(f_multipoles) :: mp
  type(dictionary), pointer :: options
  real(f_double), dimension(3) :: rxyz
  real(f_double), dimension(:), allocatable :: density
  
  call f_lib_initialize()
  call yaml_new_document()
  call yaml_argparse(options,inputs)
  n=options//'ndim'

  density=f_malloc(n,id='density')
  call f_random_number(density)

  rxyz=1.0_f_double

  !create random multipoles
  call f_multipoles_create(mp,2)

  do i=1,n
     call f_multipoles_accumulate(mp,rxyz,density(i))
  end do

  !here we may print the results of the multipole calculations
  call yaml_mapping_open('Multipoles of the array')
  call yaml_map('q0',sum(density))
  call yaml_mapping_close()

  call yaml_mapping_open('Calculated multipoles')
  call yaml_map('monopole',get_monopole(mp))
  call yaml_map('dipole',get_dipole(mp))
  call yaml_map('quadrupole',get_quadrupole(mp))
  call yaml_map('quadrupole intensity',get_quadrupole_intensities(mp))
  call yaml_mapping_close()

  call f_multipoles_release(mp)

  call f_free(density)
!  call dict_free(options)

!!$  !test of the multipole preserving routine
!!$  !initialize the work arrays needed to integrate with isf
!!$  !names of the routines to be redefined
!!$  call initialize_real_space_conversion(isf_m=mp_isf_order)
!!$
!!$  boxit = box_iter(mesh,origin=rxyz,cutoff=cutoff)
!!$  call finalize_real_space_conversion()

  call test_f_functions()

  call dict_free(options)
  !here some tests about the box usage
  call test_box_functions()

  call f_lib_finalize()

end program numeric_check

subroutine test_f_functions()
  use futile, gp => f_double
  use f_functions
  use numerics
  implicit none
  !local variables
  type(f_function) :: func1,func2!,func3
  type(f_grid_1d) :: grid
  integer :: unit

  !start with the simple evaluation
  grid=f_grid_1d_new(UNIFORM_GRID,[-1.0_gp,1.0_gp],npts=1000)
  
  !func1=f_function_new(F_GAUSSIAN,exponent=onehalf/0.01_gp)
  func1=f_function_new(F_POLYNOMIAL,coefficients=[0.0_gp,one])
  func2=f_function_new(F_POLYNOMIAL,coefficients=[0.0_gp,one])
!!$
!!$  call f_function_product(func1,func2,func3)
!!$
!!$  func3=func1*func2
!!$
  unit=12
  call f_open_file(unit=unit,file='testfunction.txt')
  call f_function_dump(unit,func1,grid)
  call f_close(unit=unit)

!!$  call f_open_file(unit=unit,file='testproduct.txt')
!!$  call f_function_dump(unit,func3,grid)
!!$  call f_close(unit=unit)
!!$
!!$  call f_free

end subroutine test_f_functions

subroutine test_box_functions()
  use futile, gp=>f_double
  use box
  use numerics, only: pi
  implicit none
  !local variables
  integer(f_long) :: tomp,tseq
  type(cell) :: mesh_ortho,mesh_noortho
  integer, dimension(3) :: ndims
  real(gp), dimension(:,:,:,:), allocatable :: v1,v2
  real(gp), dimension(3) :: angrad

  ndims=[300,300,300]

  mesh_ortho=cell_new('S',ndims,[1.0_gp,1.0_gp,1.0_gp])

  v1=f_malloc([3,ndims(1),ndims(2),ndims(3)],id='v1')
  v2=f_malloc([3,ndims(1),ndims(2),ndims(3)],id='v2')

! to be uncommented
!tbu  call loop_dotp('SEQ',mesh_ortho,v1,v2,tseq)   
!tbu  call yaml_map('Normal loop, seq (ns)',tseq)
!tbu  
!tbu  call loop_dotp('OMP',mesh_ortho,v1,v2,tomp)
!tbu  call yaml_map('Normal loop, omp (ns)',tomp)
!tbu
!tbu  call loop_dotp('ITR',mesh_ortho,v1,v2,tseq)
!tbu  call yaml_map('Normal loop, itr (ns)',tseq)
!tbu
!tbu  call loop_dotp('IOM',mesh_ortho,v1,v2,tseq)
!tbu  call yaml_map('Normal loop, iom (ns)',tseq)
!tbu
!tbu  call loop_dotp('ITM',mesh_ortho,v1,v2,tseq)
!tbu  call yaml_map('Normal loop, mpi (ns)',tseq)

  ndims=100

  mesh_ortho=cell_null()
  mesh_ortho=cell_new('F',ndims,[1.0_gp,1.0_gp,1.0_gp])
!tbu  call loop_box_function('distance',mesh_ortho)
!  call loop_box_function('box_cutoff',mesh_ortho)
   
  mesh_ortho=cell_null()
  mesh_ortho=cell_new('S',ndims,[1.0_gp,1.0_gp,1.0_gp])
!tbu  call loop_box_function('distance',mesh_ortho)
!  call loop_box_function('box_cutoff',mesh_ortho)

  mesh_ortho=cell_null()
  mesh_ortho=cell_new('P',ndims,[1.0_gp,1.0_gp,1.0_gp])
!tbu  call loop_box_function('distance',mesh_ortho)
!  call loop_box_function('box_cutoff',mesh_ortho)

  angrad(1) = 90.0_gp/180.0_gp*pi
  angrad(2) = 70.0_gp/180.0_gp*pi
  angrad(3) = 90.0_gp/180.0_gp*pi
  
  mesh_noortho=cell_new('S',ndims,[1.0_gp,1.0_gp,1.0_gp],alpha_bc=angrad(1),beta_ac=angrad(2),gamma_ab=angrad(3)) 
!tbu  call loop_box_function('distance',mesh_noortho)
!  call loop_box_function('box_cutoff',mesh_noortho)

  angrad(1) = 80.0_gp/180.0_gp*pi
  angrad(2) = 80.0_gp/180.0_gp*pi
  angrad(3) = 80.0_gp/180.0_gp*pi
  
  mesh_noortho=cell_null()
  mesh_noortho=cell_new('P',ndims,[1.0_gp,1.0_gp,1.0_gp],alpha_bc=angrad(1),beta_ac=angrad(2),gamma_ab=angrad(3)) 
  call loop_box_function('distance',mesh_noortho)
  call loop_box_function('box_cutoff',mesh_noortho)

  angrad(1) = 80.0_gp/180.0_gp*pi
  angrad(2) = 76.0_gp/180.0_gp*pi
  angrad(3) = 101.0_gp/180.0_gp*pi
  
  mesh_noortho=cell_null()
  mesh_noortho=cell_new('P',ndims,[1.0_gp,1.0_gp,1.0_gp],alpha_bc=angrad(1),beta_ac=angrad(2),gamma_ab=angrad(3)) 
  call loop_box_function('distance',mesh_noortho)
  call loop_box_function('box_cutoff',mesh_noortho)
  

  ndims=200
  angrad(1) = 20.0_gp/180.0_gp*pi
  angrad(2) = 25.0_gp/180.0_gp*pi
  angrad(3) = 30.0_gp/180.0_gp*pi
  
  mesh_noortho=cell_null()
  mesh_noortho=cell_new('P',ndims,[1.0_gp,1.0_gp,1.0_gp],alpha_bc=angrad(1),beta_ac=angrad(2),gamma_ab=angrad(3)) 
  call loop_box_function('distance',mesh_noortho)
  call loop_box_function('box_cutoff',mesh_noortho)

  call f_free(v1)
  call f_free(v2)

end subroutine test_box_functions

subroutine loop_box_function(fcheck,mesh)
  use futile
  use f_precisions
  use box
  use f_utils
  use yaml_strings
  use wrapper_MPI
  use numerics, only:pi
  implicit none
  character(len=*), intent(in) :: fcheck
  type(cell), intent(in) :: mesh
  !local variables
  integer :: i,ii 
  real(f_double) :: totvolS,totvolS1,totvolS2,totvolC,r,IntaS,IntaC,cen,errorS,errorC,d2
  real(f_double) :: diff,diff_old,dist1,dist2,cutoff,totvol_Bcutoff
  real(f_double), dimension(3) :: rxyz0,rd,rv, angdeg
  type(box_iterator) :: bit
  integer, dimension(2,3) :: nbox

  angdeg = mesh%angrad*180.0_f_double/pi
  select case(trim(fcheck))
  case('distance')
     bit=box_iter(mesh)
     r=20.0_f_double
     ! Full list of functions in box.f90 to be checked:
     ! rxyz_ortho, distance, r_wrap, closest_r, 
     ! square_gu, square_gd, dotp_gu, dotp_gd.
     call yaml_mapping_open('Check of functions distance, closest_r, rxyz_ortho')
     call yaml_map('Cell orthorhombic',mesh%orthorhombic)
     call yaml_map('Cell ndims',mesh%ndims)
     call yaml_map('Cell hgrids',mesh%hgrids)
     call yaml_map('Cell angles deg',angdeg)
     call yaml_map('Cell angles rad',mesh%angrad)
     call yaml_map('Cell periodity (FREE=0,PERIODIC=1)',mesh%bc)
     call yaml_map('Volume element',mesh%volume_element)
     call yaml_map('Contravariant matrix',mesh%gu)
     call yaml_map('Covariant matrix',mesh%gd)
     call yaml_map('Product of the two',matmul(mesh%gu,mesh%gd))
     call yaml_map('uabc matrix',mesh%uabc)
     call yaml_map('Sphere radius or cube side',r)
     do i=1,3

        totvolS=0.0_f_double
        totvolS1=0.0_f_double
        totvolS2=0.0_f_double
        totvolC=0.0_f_double
        diff=0.0_f_double
        if (i==1) cen=0.0_f_double
        if (i==2) cen=mesh%ndims(1)*0.5_f_double
        if (i==3) cen=mesh%ndims(1)*1.5_f_double
        rxyz0=[cen,cen,cen]
        if (mesh%bc(1)==0) rxyz0(1)=mesh%ndims(1)*0.5_f_double
        if (mesh%bc(2)==0) rxyz0(2)=mesh%ndims(2)*0.5_f_double
        if (mesh%bc(3)==0) rxyz0(3)=mesh%ndims(3)*0.5_f_double
        do while(box_next_point(bit))
           ! Sphere volume integral with distance
           if (distance(bit%mesh,bit%rxyz,rxyz0) .le. r) then
              totvolS=totvolS+1.0_f_double
           end if
           ! Sphere volume integral with rxyz_ortho
           rd=closest_r(bit%mesh,bit%rxyz,rxyz0)
           rv=rxyz_ortho(bit%mesh,rd)
           d2=0.0_f_double
           do ii=1,3
              d2=d2+rv(ii)**2
           end do
           dist1=sqrt(d2)
           if (dist1 .le. r) then
              totvolS1=totvolS1+1.0_f_double
           end if
           d2=square_gd(mesh,rd)
           dist2=sqrt(d2)
           ! Sphere volume integral with rxyz_ortho
           if (dist2 .le. r) then
              totvolS2=totvolS2+1.0_f_double
           end if
           diff_old=abs(dist2-dist1)
           if (diff_old.gt.diff) then
              diff=diff_old
!              call yaml_map('position x',bit%i)
!              call yaml_map('position y',bit%j)
!              call yaml_map('position z',bit%k)
!              call yaml_map('bit%rxyz',bit%rxyz)
!              call yaml_map('rxyz0',rxyz0)
!              call yaml_map('closest_r -> rd',rd)
!              call yaml_map('rxyz_ortho -> rv',rv)
!              call yaml_map('square_gd(mesh,rd) -> dist2',dist2)
!              call yaml_map('ortho square of rv -> dist1',dist1)
!              call yaml_map('diff',diff)
           end if
           ! Cube volume integral
           if ((rv(1).ge.-r .and. rv(1).lt.r) .and.&
               (rv(2).ge.-r .and. rv(2).lt.r) .and.&
               (rv(3).ge.-r .and. rv(3).lt.r)) then
              totvolC=totvolC+1.0_f_double
           end if
        end do
        totvolS=totvolS*mesh%volume_element
        totvolS1=totvolS1*mesh%volume_element
        totvolS2=totvolS2*mesh%volume_element
        totvolC=totvolC*mesh%volume_element
        IntaS=4.0_f_double/3.0_f_double*pi*r**3
!        if(mesh%bc(2)==0 .and. i==1) IntaS=IntaS*0.5_f_double
        !IntaC=(2.0_f_double*r+1.0_f_double)**3
        IntaC=(2.0_f_double*r)**3
!        if(mesh%bc(2)==0 .and. i==1) IntaC=IntaC*0.5_f_double
        errorS=abs((totvolS-IntaS)/IntaS)
        errorC=abs((totvolC-IntaC)/IntaC)
        call yaml_mapping_open('center')
        call yaml_map('Sphere or cube center',rxyz0)
        call yaml_map('Numerical sphere integral with distance',totvolS)
        call yaml_map('Numerical sphere integral with closest_r and rxyz_ortho',totvolS1)
        call yaml_map('Numerical sphere integral with closest_r and square_gd',totvolS2)
        call yaml_map('Analytical sphere integral',IntaS)
        call yaml_map('Sphere integral error',errorS)
        call yaml_map('Numerical cube integral',totvolC)
        call yaml_map('Analytical cube integral',IntaC)
        call yaml_map('Cube integral error',errorC)
        call yaml_map('Maximum difference between closest_r and square_gd',diff)
        call yaml_mapping_close()
     end do
     call yaml_mapping_close()
  case('box_cutoff')
     r=20.0_f_double
     cutoff = r
     ! To check the functions box_nbox_from_cutoff and 
     ! cell_cutoff_extremao of box.f90.
     call yaml_mapping_open('Check of box cutoff')
     call yaml_map('Cell orthorhombic',mesh%orthorhombic)
     call yaml_map('Cell ndims',mesh%ndims)
     call yaml_map('Cell hgrids',mesh%hgrids)
     call yaml_map('Cell angles deg',angdeg)
     call yaml_map('Cell angles rad',mesh%angrad)
     call yaml_map('Cell periodity (FREE=0,PERIODIC=1)',mesh%bc)
     call yaml_map('Volume element',mesh%volume_element)
     call yaml_map('Contravariant matrix',mesh%gu)
     call yaml_map('Covariant matrix',mesh%gd)
     call yaml_map('Product of the two',matmul(mesh%gu,mesh%gd))
     call yaml_map('uabc matrix',mesh%uabc)
     call yaml_map('Sphere radius or cube side',r)
     call yaml_map('Box cube cutoff',cutoff)
!     do i=1,3
        i=2
        totvolC=0.0_f_double
        totvol_Bcutoff=0.0_f_double
        diff=0.0_f_double
        if (i==1) cen=0.0_f_double
        if (i==2) cen=mesh%ndims(1)*0.5_f_double
        if (i==3) cen=mesh%ndims(1)*1.5_f_double
        rxyz0=[cen,cen,cen]
        if (mesh%bc(1)==0) rxyz0(1)=mesh%ndims(1)*0.5_f_double
        if (mesh%bc(2)==0) rxyz0(2)=mesh%ndims(2)*0.5_f_double
        if (mesh%bc(3)==0) rxyz0(3)=mesh%ndims(3)*0.5_f_double
        call yaml_map('Sphere or cube center',rxyz0)
        bit=box_iter(mesh)
        call yaml_map('bit%nbox whole box',bit%nbox)
        do while(box_next_point(bit))
           rd=closest_r(bit%mesh,bit%rxyz,rxyz0)
           rv=rxyz_ortho(bit%mesh,rd)
           d2=0.0_f_double
           do ii=1,3
              d2=d2+rv(ii)**2
           end do
           dist1=sqrt(d2)
           d2=square_gd(mesh,rd)
           dist2=sqrt(d2)
           diff_old=abs(dist2-dist1)
           if (diff_old.gt.diff) then
              diff=diff_old
           end if
           ! Cube volume integral
           if ((rv(1).ge.-r .and. rv(1).le.r) .and.&
               (rv(2).ge.-r .and. rv(2).le.r) .and.&
               (rv(3).ge.-r .and. rv(3).le.r)) then
              totvolC=totvolC+1.0_f_double
           end if
        end do
        totvolC=totvolC*mesh%volume_element
        IntaC=(2.0_f_double*r)**3
        errorC=abs((totvolC-IntaC)/IntaC)
        
        ! Using the internal nbox
        !nbox = box_nbox_from_cutoff(mesh,rxyz0,cutoff)
        !call yaml_map('nbox',nbox)
        !bit=box_iter(mesh,nbox)
        bit=box_iter(mesh,origin=rxyz0,cutoff=cutoff)
        call yaml_map('bit%nbox reduced box',bit%nbox)
        do while(box_next_point(bit))
           totvol_Bcutoff=totvol_Bcutoff+1.0_f_double
        end do
        totvol_Bcutoff=totvol_Bcutoff*mesh%volume_element

        call yaml_mapping_open('center')
        call yaml_map('Sphere or cube center',rxyz0)
        call yaml_map('Numerical cube integral',totvolC)
        call yaml_map('Analytical cube integral',IntaC)
        call yaml_map('Cube integral error',errorC)
        call yaml_map('Maximum difference between closest_r and square_gd',diff)
        call yaml_map('Numerical box cutoff integral',totvol_Bcutoff)
        call yaml_mapping_close()
!     end do
     call yaml_mapping_close()
  case('other')
  end select

end subroutine loop_box_function

subroutine loop_dotp(strategy,mesh,v1,v2,time)
  use f_precisions
  use box
  use f_utils
  use yaml_strings
  use wrapper_MPI
  implicit none
  character(len=*), intent(in) :: strategy
  type(cell), intent(in) :: mesh
  real(f_double), dimension(3,mesh%ndims(1),mesh%ndims(2),mesh%ndims(3)), intent(inout) :: v1
  real(f_double), dimension(3,mesh%ndims(1),mesh%ndims(2),mesh%ndims(3)), intent(inout) :: v2
  integer(f_long), intent(out) :: time
  !local variables
  integer :: i1,i2,i3,n3p,i3s,ithread,nthread
  integer(f_long) :: t0,t1
  real(f_double) :: totdot,res
  type(box_iterator) :: bit
  integer, dimension(2,3) :: nbox
  !$ integer, external ::  omp_get_thread_num,omp_get_num_threads

  !initialization
  do i3=1,mesh%ndims(3)
     do i2=1,mesh%ndims(2)
        do i1=1,mesh%ndims(1)
           !the scalar product of these objects is 20.0
           v1(:,i1,i2,i3)=[1.0_f_double,2.0_f_double,3.0_f_double]
           v2(:,i1,i2,i3)=[2.0_f_double,3.0_f_double,4.0_f_double]
        end do
     end do
  end do
  t0=0
  t1=0

  select case(trim(strategy))
  case('SEQ')
     totdot=0.0_f_double
     t0=f_time()
     do i3=1,mesh%ndims(3)
        do i2=1,mesh%ndims(2)
           do i1=1,mesh%ndims(1)
              res=dotp_gd(mesh,v1(1,i1,i2,i3),v2(:,i1,i2,i3))
              res=res/20.0_f_double
              totdot=totdot+res
              v2(:,i1,i2,i3)=res
           end do
        end do
     end do
     t1=f_time()
  case('OMP')
     totdot=0.0_f_double
     t0=f_time()
     !$omp parallel do default(shared) &
     !$omp private(i1,i2,i3,res)&
     !$omp reduction(+:totdot)
     do i3=1,mesh%ndims(3)
        do i2=1,mesh%ndims(2)
           do i1=1,mesh%ndims(1)
              res=dotp_gd(mesh,v1(1,i1,i2,i3),v2(:,i1,i2,i3))
              res=res/20.0_f_double
              totdot=totdot+res
              v2(:,i1,i2,i3)=res
           end do
        end do
     end do
     !$omp end parallel do
     t1=f_time()
  case('ITR') !iterator case
     bit=box_iter(mesh)
     totdot=0.0_f_double
     t0=f_time()
     do while(box_next_point(bit))
        res=dotp_gd(bit%mesh,v1(1,bit%i,bit%j,bit%k),v2(:,bit%i,bit%j,bit%k))
        res=res/20.0_f_double
        totdot=totdot+res
        v2(:,bit%i,bit%j,bit%k)=res
     end do
     t1=f_time()
  case('IOM') !iterator with omp
     bit=box_iter(mesh)
     totdot=0.0_f_double
     t0=f_time()
     nthread=1
     !$omp parallel default(shared) &
     !$omp private(res,ithread)&
     !$omp firstprivate(bit)&
     !$omp reduction(+:totdot)
     ithread=0
     !$ ithread=omp_get_thread_num()
     !$ nthread=omp_get_num_threads()
     call box_iter_split(bit,nthread,ithread)
     do while(box_next_point(bit))
        res=dotp_gd(bit%mesh,v1(1,bit%i,bit%j,bit%k),v2(:,bit%i,bit%j,bit%k))
        res=res/20.0_f_double
        totdot=totdot+res
        v2(:,bit%i,bit%j,bit%k)=res
     end do
     call box_iter_merge(bit)
     !$omp end parallel
     t1=f_time()
  case('ITM') !iterator with mpi
     call mpiinit()
     nbox(1,:)=1
     nbox(2,:)=mesh%ndims
     call distribute_on_tasks(mesh%ndims(3),mpirank(),mpisize(),n3p,i3s)
     nbox(1,3)=i3s+1
     nbox(2,3)=i3s+n3p
     !print *,'here',mpisize(),mpirank(),i3s,n3p,mesh%ndims(3)
     bit=box_iter(mesh,i3s=i3s+1,n3p=n3p)
     totdot=0.0_f_double
     t0=f_time()
     do while(box_next_point(bit))
        res=dotp_gd(bit%mesh,v1(1,bit%i,bit%j,bit%k),v2(:,bit%i,bit%j,bit%k))
        res=res/20.0_f_double
        totdot=totdot+res
        v2(:,bit%i,bit%j,bit%k)=res
     end do
     call fmpi_allreduce(totdot,1,op=FMPI_SUM)
     !call mpigather
     t1=f_time()
     call mpifinalize()
  end select
  
  !totdot should be the size of the array
  call f_assert(int(totdot,f_long) == mesh%ndim,&
       'Wrong reduction, found "'+totdot+'" instead of "'+mesh%ndim+'"')
  !call yaml_map('TotDot',totdot)
  !totsum should be the size of the array multipliled by three
  call f_assert(sum(v2) == f_size(v2),&
       'Wrong array writing, found "'+sum(v2)+'" instead of "'+f_size(v2)+'"')
  !call yaml_map('TotSum',sum(v2))

  time=t1-t0

end subroutine loop_dotp

! Parallelization a number n over nproc nasks
subroutine distribute_on_tasks(n, iproc, nproc, np, is)
  implicit none
  ! Calling arguments
  integer,intent(in) :: n, iproc, nproc
  integer,intent(out) :: np, is

  ! Local variables
  integer :: ii

  ! First distribute evenly... (LG: if n is, say, 34 and nproc is 7 - thus 8 MPI processes)
  np = n/nproc                !(LG: we here have np=4) 
  is = iproc*np               !(LG: is=iproc*4 : 0,4,8,12,16,20,24,28)
  ! ... and now distribute the remaining objects.
  ii = n-nproc*np             !(LG: ii=34-28=6)
  if (iproc<ii) np = np + 1   !(LG: the first 6 tasks (iproc<=5) will have np=5)
  is = is + min(iproc,ii)     !(LG: update is, so (iproc,np,is): (0,5,0),(1,5,5),(2,5,10),(3,5,15),(4,5,20),(5,5,25),(6,4,30),(7,4,34))

end subroutine distribute_on_tasks
