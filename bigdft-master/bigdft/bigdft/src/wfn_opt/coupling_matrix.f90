!> @file
!!  Routines related to coupling matrix (TD-DFT Casida's formalism)
!! @author
!!    Copyright (C) 2009-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
subroutine center_of_charge(at,rxyz,cc)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(3), intent(out) :: cc
  !local variables
  integer :: iat,ityp
  real(gp) :: zatom,rx,ry,rz,qtot

  cc(1)=0.0_gp
  cc(2)=0.0_gp
  cc(3)=0.0_gp
  qtot=0.0_gp
  do iat=1,at%astruct%nat
     ityp=at%astruct%iatype(iat)
     zatom=real(at%nelpsp(ityp),gp)
     qtot=qtot+zatom
     !coordinates of the center
     rx=rxyz(1,iat)
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)
     cc(1)=cc(1)+rx*zatom
     cc(2)=cc(2)+ry*zatom
     cc(3)=cc(3)+rz*zatom
  end do
  if (qtot /= 0.0_gp) then
     cc(1)=cc(1)/qtot
     cc(2)=cc(2)/qtot
     cc(3)=cc(3)/qtot
  end if
END SUBROUTINE center_of_charge


!> Calculate the coupling matrix needed for Casida's TDDFT approach
subroutine calculate_coupling_matrix(iproc,nproc,dir_output,boxit,tddft_approach,nspin,ndimp,orbsocc,orbsvirt,&
     center_of_charge,pkernel,dvxcdrho,psirocc,psivirtr)
  use module_base
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use yaml_output
  use box
  use f_harmonics
  implicit none
  character(len=4), intent(in) :: tddft_approach
  character(len=*), intent(in) :: dir_output
  integer, intent(in) :: iproc,nproc,nspin,ndimp
  real(gp), dimension(3) :: center_of_charge
  type(orbitals_data), intent(in) :: orbsocc,orbsvirt
  real(wp), dimension(ndimp,orbsocc%norb), intent(in) :: psirocc
  real(wp), dimension(ndimp,orbsvirt%norb), intent(in) :: psivirtr
  real(wp), dimension(ndimp,max((nspin*(nspin+1))/2,2)), intent(in) :: dvxcdrho
  type(coulomb_operator) :: pkernel
  type(box_iterator) :: boxit
  !local variables
  integer, parameter :: ALPHA_=2,P_=1,SPIN_=3
  integer :: imulti,jmulti,spinindex,ialpha,ip,ibeta,iq,ispin,jspin,ntda
  integer :: nmulti,ndipoles,iap,ibq,nalphap,istep
  real(wp) :: eap,ebq,krpa,kfxc,q,rho_bq,kfxc_od
  type(f_progress_bar) :: bar
  type(f_multipoles) :: mp
  integer, dimension(:,:), allocatable :: transitions
  real(wp), dimension(:,:), allocatable :: Kaux!,dipoles
  real(wp), dimension(:,:), allocatable :: transition_quantities
  real(wp), dimension(:,:), pointer :: Kbig,K
  real(wp), dimension(:), allocatable :: v_ias
  real(wp), dimension(:), allocatable :: rho_ias
  !$ integer, external :: omp_get_max_threads
  !$ integer :: nthreads

  call f_routine('calculate_coupling_matrix')

  if(iproc==0) call yaml_comment('Linear-Response TDDFT calculations',hfill='-')

  call allocate_transitions_lookup()
  nmulti=nalphap
  ndipoles=nmulti
  if (nspin == 1) ndipoles=2*nmulti

  !Allocate partial densities and potentials
  !rho_ias = f_malloc([ndimp,nmulti],id='rho_ias')
  rho_ias = f_malloc(ndimp,id='rho_ias')
  v_ias = f_malloc(ndimp,id='v_ias')

  !Allocation of dipoles (computed in order to get the oscillator strength)
  !dipoles = f_malloc0([3, ndipoles],id='dipoles')
  transition_quantities = f_malloc0([10, ndipoles],id='transition_quantities')
  K = f_malloc0_ptr([nmulti, nmulti],id='K')
  !For nspin=1, define an auxiliary matrix for spin-off-diagonal terms.
  if (nspin==1) Kaux = f_malloc0([nmulti, nmulti],id='Kaux')

  if (iproc==0) bar=f_progress_bar_new(nstep=((nalphap+1)*nalphap)/2)

  !$ nthreads=omp_get_max_threads()

  call PS_set_options(pkernel,verbose=.false.)
  istep=0
  do iap=1,nalphap
       ialpha=transitions(ALPHA_,iap)
       ip=transitions(P_,iap)
       ispin=transitions(SPIN_,iap)
       eap=orbsvirt%eval(ialpha)-orbsocc%eval(ip)

       call f_multipoles_create(mp,2,center=center_of_charge)
       !extraction of the coefficients d_i
       do while(box_next_point(boxit))
         !fill the rho_ias array
         rho_ias(boxit%ind)=psirocc(boxit%ind,ip)*psivirtr(boxit%ind,ialpha)/boxit%mesh%volume_element
         q=rho_ias(boxit%ind)
         boxit%tmp=boxit%rxyz-center_of_charge
         !dipoles(:,iap)=dipoles(:,iap)+boxit%tmp*q
         call f_multipoles_accumulate(mp,boxit%tmp,q*boxit%mesh%volume_element)
       end do
       transition_quantities(:,iap)=mp%monomials
       call f_multipoles_release(mp)
       !dipoles(:,iap)=sqrt(eap)*dipoles(:,iap)*boxit%mesh%volume_element

       !for every rho iap  calculate the corresponding potential
       !copy the transition  density in the inout structure
       call f_memcpy(src=rho_ias,dest=v_ias)
       call Electrostatic_Solver(pkernel,v_ias)

       !now we have to calculate the corresponding element of the RPA part of the coupling matrix
       !$omp parallel do if (iap>=6*nthreads .and. nthreads>1) default(none) &
       !$omp shared(transitions,orbsvirt,orbsocc,ispin,iap,psirocc,psivirtr)&
       !$omp shared(v_ias,rho_ias,dvxcdrho,K,Kaux,nspin,tddft_approach,eap)&
       !$omp private(ibq,ibeta,iq,jspin,ebq,spinindex,krpa,kfxc,kfxc_od,rho_bq)&
       !$omp firstprivate(boxit)
       do ibq=1,iap
         ibeta=transitions(ALPHA_,ibq)
         iq=transitions(P_,ibq)
         jspin=transitions(SPIN_,ibq)
         ebq=orbsvirt%eval(ibeta)-orbsocc%eval(iq)

         spinindex=ispin+jspin-1

         !and of the rpa part
         !reconstruct the density on the fly
         krpa=0.0_dp
         kfxc=0.0_wp
         kfxc_od=0.0_wp
         do while(box_next_point(boxit))
            rho_bq=psirocc(boxit%ind,iq)*psivirtr(boxit%ind,ibeta) !volume element for the integration is considered inside
            krpa=krpa+v_ias(boxit%ind)*rho_bq
            kfxc=kfxc+rho_ias(boxit%ind)*rho_bq*dvxcdrho(boxit%ind,spinindex)
            !in the nspin=1 case we also have to calculate the off-diagonal term
            if (nspin==1) kfxc_od=kfxc_od+rho_ias(boxit%ind)*rho_bq*dvxcdrho(boxit%ind,2)
         end do
!!$         krpa=dot(ndimp,rho_ias(1,ibq),1,v_ias(1),1)*boxit%mesh%volume_element

         K(iap,ibq)=krpa+kfxc
         if (nspin==1) Kaux(iap,ibq)=krpa+kfxc_od

!!$         K(iap,ibq)=krpa
!!$         if (nspin==1) Kaux(iap,ibq)=krpa
!!$         !calculate the fxc term
!!$         kfxc=0.0_wp
!!$         do while(box_next_point(boxit))
!!$           kfxc=kfxc+rho_ias(boxit%ind,iap)*rho_ias(boxit%ind,ibq)*dvxcdrho(boxit%ind,spinindex)*boxit%mesh%volume_element
!!$         end do
!!$         K(iap,ibq)=K(iap,ibq)+kfxc
!!$         !in the nspin=1 case we also have to calculate the off-diagonal term
!!$         if (nspin==1) then
!!$           kfxc=0.0_wp
!!$           do while(box_next_point(boxit))
!!$             kfxc=kfxc+rho_ias(boxit%ind,iap)*rho_ias(boxit%ind,ibq)*dvxcdrho(boxit%ind,2)*boxit%mesh%volume_element
!!$           end do
!!$           Kaux(iap,ibq)=Kaux(iap,ibq)+kfxc
!!$         end if

         !add factors from energy occupation numbers (for non-tda case)
         !If full TDDFT, then multiply the coupling matrix element by the "2*sqrt(\omega_q*\omega_{q'})" coefficient of eq. 2.
         if (tddft_approach=='full') then
           K(iap,ibq)=K(iap,ibq)*2.0_wp*sqrt(eap)*sqrt(ebq)
           if (nspin ==1) Kaux(iap,ibq)=Kaux(iap,ibq)*2.0_wp*sqrt(eap)*sqrt(ebq)
         end if
         !istep=istep+1
      end do
      !$omp end parallel do
       if (iproc==0) call dump_progress_bar(bar,step=((iap+1)*iap)/2)
    end do
  !If more than one processor, then perform the MPI_all_reduce of K (and of Kaux if nspin=1) and of dipoles.
  if (nproc > 1) then
     call fmpi_allreduce(K,FMPI_SUM,comm=bigdft_mpi%mpi_comm)
     if (nspin ==1) call fmpi_allreduce(Kaux,FMPI_SUM,comm=bigdft_mpi%mpi_comm)
     !call fmpi_allreduce(dipoles(1,1),3*nmulti,FMPI_SUM,comm=bigdft_mpi%mpi_comm)
     call fmpi_allreduce(transition_quantities(1,1),10*nmulti,FMPI_SUM,comm=bigdft_mpi%mpi_comm)
  end if


  !Build the upper triangular part of K (and Kaux if nspin=1)
  do imulti=1,nmulti
     do jmulti=imulti+1,nmulti
        K(imulti,jmulti)=K(jmulti,imulti)
     end do
  end do
  if (nspin==1) then
     do imulti=1,nmulti
        do jmulti=imulti+1,nmulti
           Kaux(imulti,jmulti)=Kaux(jmulti,imulti)
        end do
     end do
  end if

  !Exponent used to differentiate full TDDFT and TDA
  ntda=1
  if (tddft_approach=='TDA') then
     ntda=0
  else if (tddft_approach=='full') then
     ntda=1
  end if

  !Add the diagonal part: {\omega_i}^{1+ntda} \delta_{i,j}
  do iap=1,nalphap
     ialpha=transitions(ALPHA_,iap)
     ip=transitions(P_,iap)
     eap=orbsvirt%eval(ialpha)-orbsocc%eval(ip)
     K(iap,iap)=K(iap,iap)+eap**(ntda+1)
     !put the sqrt of the transition energies in the first place of the array
     !as the monopole of a codensity is always 0
     transition_quantities(1,iap)=sqrt(eap)
  end do
  !Construction of the whole coupling matrix Kbig for nspin=1
  if (nspin == 1) then
     Kbig = f_malloc0_ptr((/ 2*nmulti, 2*nmulti /),id='Kbig')
     do iap=1,nalphap
        do ibq=1,nalphap
           Kbig(iap,ibq)=K(iap,ibq)
           Kbig(iap+nmulti,ibq+nmulti)=K(iap,ibq)
           Kbig(iap+nmulti,ibq)=Kaux(iap,ibq)
           Kbig(iap,ibq+nmulti)=Kaux(iap,ibq)
        end do
     end do
     !Copy the values of the dipoles in the second part of the array
     !call f_memcpy(n=3*nmulti,src=dipoles(1,1),dest=dipoles(1,nmulti+1))
     call f_memcpy(n=10*nmulti,src=transition_quantities(1,1),&
          dest=transition_quantities(1,nmulti+1))
  else
    Kbig=>K
  end if

  !here the matrix can be written on disk, together with the transition dipoles
  if (iproc==0) then
     call f_savetxt(dir_output//'coupling_matrix.txt',Kbig)
     !call f_savetxt(dir_output//'transition_dipoles.txt',dipoles)
     call f_savetxt(dir_output//'transition_quantities.txt',transition_quantities)
  end if

  !Deallocations
  call f_free_ptr(K)
  !call f_free(dipoles)
  call f_free(transition_quantities)
  call f_free(transitions)
  if (nspin ==1 ) then
     call f_free(Kaux)
     call f_free_ptr(Kbig)
  end if
  call f_free(rho_ias)
  call f_free(v_ias)

  call f_release_routine()

contains

  !> internal routine to define the lookup arrays
  subroutine allocate_transitions_lookup()
    use dictionaries
    implicit none
    integer :: imulti,iorbi,iorba,i
    integer, dimension(3) :: tmp
    type(dictionary), pointer :: lookup,local

    call dict_init(lookup)
    loop_i: do imulti=1,orbsvirt%norb*orbsocc%norb
      iorbi=(imulti-1)/orbsvirt%norb+1 !occ. state index
      iorba=imulti-(iorbi-1)*orbsvirt%norb !virt. state index

      !Check the spin of the occ. and virt. orbitals considered
      if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
        !Check the sign of the spin (used later for the choice of exchange correlation kernel).
        if (orbsocc%spinsgn(iorbi) == 1.0_gp) then
          ispin=1
        else if (orbsocc%spinsgn(iorbi) == -1.0_gp) then
          ispin=2
        end if
      else
        !If the orbitals do not have the same spin, then cycle
        cycle loop_i
      end if
      tmp(P_)=iorbi
      tmp(ALPHA_)=iorba
      tmp(SPIN_)=ispin
      local=>list_new(.item. tmp)
      call add(lookup,local)
   end do loop_i

   nalphap=dict_len(lookup)

   !now allocate lookup array accordingly
   transitions=f_malloc([3,nalphap],id='transitions')
   do iap=0,nalphap-1
      !here we write the transitions information
      tmp=lookup//iap
      do i=1,3
         transitions(i,iap+1)=tmp(i)
      end do
   end do
   call dict_free(lookup)

 end subroutine allocate_transitions_lookup
end subroutine calculate_coupling_matrix

!> Calculate the coupling matrix needed for Casida's TDDFT approach
subroutine coupling_matrix_prelim(iproc,nproc,geocode,tddft_approach,nspin,lr,orbsocc,orbsvirt,i3s,n3p,&
     hxh,hyh,hzh,chargec,pkernel,dvxcdrho,psirocc,psivirtr,exc_fac)
  use module_base
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use yaml_output
  use bounds, only: ext_buffers
  use locregs
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  character(len=4), intent(in) :: tddft_approach
  integer, intent(in) :: iproc,nproc,n3p,nspin,i3s
  real(gp), intent(in) :: hxh,hyh,hzh
  real(gp), dimension(3) :: chargec
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbsocc,orbsvirt
  real(wp), dimension(lr%d%n1i*lr%d%n2i*n3p*orbsocc%norb), intent(in) :: psirocc
  real(wp), dimension(lr%d%n1i*lr%d%n2i*n3p*orbsvirt%norb), intent(in) :: psivirtr
  real(wp), dimension(lr%d%n1i,lr%d%n2i,n3p,max((nspin*(nspin+1))/2,2)), intent(in) :: dvxcdrho
  type(coulomb_operator) :: pkernel
  real(gp), intent(in) :: exc_fac
  !local variables
  character(len=*), parameter :: subname='coupling_matrix_prelim'
  !logical :: tda=.true.
  logical :: dofH=.true.,dofxc=.true.,dodiag=.true.,perx,pery,perz
  integer :: imulti,jmulti,jorba,jorbi,spinindex
  integer :: i1,i2,i3p,iorbi,iorba,indi,indj,inda,indb,ind2,ind3,ntda,ispin,jspin
  integer :: ik,jk,nmulti,lwork,info,nbl1,nbl2,nbl3,nbr3,nbr2,nbr1,ndipoles
  real(gp) :: ehart,hfac,x,y,z,fsumrule_test
  real(wp), dimension(:), allocatable :: omega,work
  real(wp), dimension(:,:), allocatable :: K,Kbig,Kaux,dipoles,fi
  real(wp), dimension(:,:,:), allocatable :: v_ias
  real(wp), dimension(:,:,:,:), allocatable :: rho_ias
  real(wp), dimension(:,:,:), allocatable :: rhov_ij
  real(wp), dimension(:,:,:), allocatable :: rho_ab
  real(wp) :: K_u_HF

  !Preliminary comments :
  !
  ! In Linear Response TDDFT, under an adiabatic approximation, one has to solve Casida's equation:
  ! ( A   B ) (X)          ( 1   0 ) (X)
  ! (       ) ( ) = \omega (       ) ( )                                                                        (eq. 1)
  ! ( B*  A*) (Y)          ( 0  -1 ) (Y)
  ! where A and B are matrices defined by:
  ! A_{q,q'} = \omega_q \delta_{q,q'} + K_{q,q'},
  ! B_{q,q'} = K_{q,q'}.
  !
  ! q is a shortcut for a,i,\sigma, where a represents a Kohn-Sham (KS) virtual state, i a KS occupied state,
  ! and sigma the spin of both states (which have to be equal).
  ! \omega_q = E_a - A_i, with E_a (respectively E_i) the energy of the virtual (resp. occupied) state.
  !
  ! The important element is the coupling matrix K, its elements being defined by:
  ! K_{q,q'} = \int dr \int dr' \rho_q(r) f_{Hxc}(r,r') \rho(r')
  ! f_{Hxc} is the sum of Hartree (H) and exchange correlation (xc) kernels: f_{Hxc} = f_H + f_{xc}.
  ! \rho_q(r) = \psi_a(r) \psi_i(r) is the partial density of the KS transition.
  !
  ! One can reduce without any approximation the size of the problem of eq. 1 by recasting it into:
  ! \Omega_{q,q'} F = {\omega_q}^2 \delta_{q,q'} + 2 sqrt(\omega_q \omega_{q'}) K_{q,q'} F = \omega^2 F         (eq. 2)
  ! This will be called full TDDFT in the following (when the input parameter tddft_approach is set to 'full')
  !
  ! The Tamm-Dancoff Approximation can also be applied, and one needs only to solve:
  ! A X = \omega_q \delta_{q,q'} + K_{q,q'} X = \omega X                                                        (eq. 3)
  ! This will be icalled TDA in the following (when the input parameter tddft_approach is set to 'TDA').
  !
  ! In the implementation, we will use the fact that both equations 2 and 3 have a similar form:
  ! {\omega_{q}}^{1+n} \delta_{q,q^\prime} + ( 2 \sqrt{\omega_{q} \omega_{q^\prime}} )^{n} K_{q,q^\prime}
  !                     = \omega^{1+n} F^{n} X^{1-n}                                                            (eq. 4)
  ! The exponent n will allow the differentiation between full TDDFT or TDA:
  ! - if n=0, then TDA calculation is performed
  ! - if n=1, then full TDDFT calculation is performed
  ! In the implementation, n is named 'ntda'.
  !
  !
  ! By looking at the kernel f_{Hxc} definition, one sees that there is no spin dependance for the Hartree part,
  ! while there is one for the exchange-correlation part.
  ! This is of importance in the implementation, since spin-averaged or spin-polarized calculation can be done.
  ! In the spin averaged case (nspin=1), each orbital has an occupation number of 0 or 2, while
  ! in the spin polarized case (nspin=2), each orbital has an occupation number of 0 or 1.
  ! This means that, in the spin averaged case, one can divide the coupling matrix into four sub-matrices,
  ! two using one definition of f_{xc}, the other two using another definition.
  ! The whole coupling matrix still needs to be computed in the spin-polarized case.
  !
  !
  ! To summarize, this means that the implementation needs to take into account two parameters:
  !1- Is the full-TDDFT problem solved, or is the TDA applied ?
  !2- What about the spin: spin-averaged or spin-polarized ?
  call f_routine('coupling_matrix_prelim')

  if(iproc==0) call yaml_comment('Linear-Response TDDFT calculations',hfill='-')
  !if(iproc==0) write(*,'(1x,a)')"=========================================================="
  !if(iproc==0) write(*,'(1x,a)')" Linear-Response TDDFT calculations"



  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)


  !Calculate nmulti, the number of allowed transitions between an occupied KS state and a virtual KS state.
  !Both states must have the same spin for the transition to occur.
  !For that purpose, we use a multi-index imulti: imulti=a+(i-1)*nvirt, where:
  !a is the virtual state index, i is the occupied state index, and nvirt is the total number of virtual states.
  !The length of this multi-index is the number of occupied states times the number of virtual states (nocc*nvirt)
  nmulti=0 !initialize the counter of allowed transitions
  do imulti=1,orbsvirt%norb*orbsocc%norb
     !calculate the orbital index
     iorbi=(imulti-1)/orbsvirt%norb+1 !index of the occupied state considered
     iorba=imulti-(iorbi-1)*orbsvirt%norb !index of the virtual state considered
     !check if the spin of both orbitals is the same
     if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
        !if same spin, then increment the counter of allowed transitions
        nmulti=nmulti+1
     end if
  end do

  !Allocate partial densities and potentials
  rho_ias = f_malloc((/ lr%d%n1i, lr%d%n2i, max(n3p,1), nmulti /),id='rho_ias')
  v_ias = f_malloc((/ lr%d%n1i, lr%d%n2i, max(n3p,1) /),id='v_ias')
  if (exc_fac .ne. 0.0 .and. ispin==jspin) then
     rhov_ij = f_malloc((/ lr%d%n1i, lr%d%n2i, max(n3p,1) /),id='rhov_ij')
     rho_ab  = f_malloc((/ lr%d%n1i, lr%d%n2i, max(n3p,1) /),id='rho_ab')
  end if

  !Preliminary comments on the spin :
  !1- In the spin-averaged case (nspin=1), the whole coupling matrix Kbig is made of 4 symmetric sub-matrices of size nmulti*nmulti.
  !   This is due to the fact that each orbital has an occupation number of 0 or 2.
  !   The K_H (Hartree) terms are the same in each of these 4 sub-matrices.
  !   However, the exchange correlation part depends on the spin,
  !   so that for the diagonal sub-matrices K, one definition of f_xc is used,
  !   while another one is for the off-diagonal sub-matrices Kaux.
  !2- In the spin-polarized case (nspin=2), the whole coupling matrix K (and not Kbig) is a matrix of size nmulti*nmulti.
  !   This is due to the fact that each orbital has an occupation number of 0 or 1.
  !   There is no spin dependence on Hartree part of the coupling matrix K_H
  !   However, the exchange correlation part depends on the spin,
  !   so that three definitions of the exchange correlation are used :
  !   one for f_xc^{up,up}, one for f_xc^{down,down} and one for f_xc^{down,up}==f_xc^{up,down}.

  !We define the size of the whole coupling matrix (and the number of dipoles to be considered)
  if (nspin == 1) then
     !If spin-averaged calculation, then we have to take into account the fact that each
     !orbital represents two atoms, so the dimension of the coupling matrix is doubled.
     ndipoles=2*nmulti
  else
     !otherwise, the size of the coupling matrix is nmulti.
     ndipoles=nmulti
  end if

  !Allocation of dipoles (computed in order to get the oscillator strength)
  dipoles = f_malloc0((/ 3, ndipoles /),id='dipoles')
  !call to_zero(3*ndipoles,dipoles)

  !Allocate the coupling matrix elements.
  !Recall that K has two different meanings:
  !- if nspin=1 then it represents one sub-matrix of the whole coupling matrix
  !- else it represents the whole coupling matrix.
  K = f_malloc0((/ nmulti, nmulti /),id='K')

  !For nspin=1, define an auxiliary matrix for spin-off-diagonal terms.
  if (nspin==1) then
     Kaux = f_malloc0((/ nmulti, nmulti /),id='Kaux')
  end if

  !Exponent used to differentiate full TDDFT and TDA
  if (tddft_approach=='TDA') then
     ntda=0
  else if (tddft_approach=='full') then
     ntda=1
  end if

  hfac=1.0_gp/(hxh*hyh*hzh)

  !Now we can start to build the partial densities and the corresponding partial potentials.
  !loop_i is a loop over a multiindex, the same as the one used above to find the number of allowed transitions.
  !The multindex is organised as I=a+(i-1)*norbv.
  !!$$if (iproc==0) call yaml_comment('Before computation of all the different coupling matrix elements')
  ik=0
  !!$$if (iproc==0) call yaml_sequence_open('Sequential building of the coupling matrix element')
  loop_i: do imulti=1,orbsvirt%norb*orbsocc%norb
     !!$$if (iproc==0) call yaml_comment(trim(yaml_toa(imulti,fmt='(i8.8)')))
     !Calculate the orbital index
     iorbi=(imulti-1)/orbsvirt%norb+1 !occ. state index
     iorba=imulti-(iorbi-1)*orbsvirt%norb !virt. state index

     !Check the spin of the occ. and virt. orbitals considered
     if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
        !If they have the same spin, then increment the counter of transitions
        ik=ik+1
        !Check the sign of the spin (used later for the choice of exchange correlation kernel).
        if (orbsocc%spinsgn(iorbi) == 1.0_gp) then
           ispin=1
        else if (orbsocc%spinsgn(iorbi) == -1.0_gp) then
           ispin=2
        end if
     else
        !If the orbitals do not have the same spin, then cycle
        cycle loop_i
     end if


     !Partial density definition (psi_occ*psi_virt) and dipole moments calculation.
     !We use three loops, one for each direction of space.
     !Another multi-index is used here for the position of the grid point.
     do i3p=1,n3p !loop over z
        ind3=(i3p-1)*lr%d%n1i*lr%d%n2i !z-index
        z=real(i3p+i3s-nbl3-1,gp)*hzh-chargec(3) !z (origin at the center of charge)

        do i2=1,lr%d%n2i !loop over y
           ind2=(i2-1)*lr%d%n1i+ind3 !y-index
           y=real(i2-nbl2-1,gp)*hyh-chargec(2) !y (origin at the center of charge)

           do i1=1,lr%d%n1i !loop over x
              x=real(i1-nbl1-1,gp)*hxh-chargec(1) !x (origin at the center of charge)

              indi=i1+ind2+(iorbi-1)*lr%d%n1i*lr%d%n2i*n3p !multiindex giving the right index for the occupied state wavefunction
              inda=i1+ind2+(iorba-1)*lr%d%n1i*lr%d%n2i*n3p !multiindex giving the right index for the virtual state wavefunction
              rho_ias(i1,i2,i3p,ik)=hfac*psivirtr(inda)*psirocc(indi) !partial density, multiplied by the integration factor hfac

              !calculate dipole moments
              dipoles(1,ik)=dipoles(1,ik)+x*rho_ias(i1,i2,i3p,ik) !dipole moment along x-axis
              dipoles(2,ik)=dipoles(2,ik)+y*rho_ias(i1,i2,i3p,ik) !dipole moment along y-axis
              dipoles(3,ik)=dipoles(3,ik)+z*rho_ias(i1,i2,i3p,ik) !dipole moment along z-axis

           end do !loop over x
        end do !loop over y
     end do !loop over z
     !!$$if (iproc==0) call yaml_comment('dipole moment computed')

     !Calculate the Hartree potential corresponding to the partial density
     if (dofH) then
        !Copy the partial density onto the partial potential space to pass it to PSolver
        if (n3p > 0) call vcopy(lr%d%n1i*lr%d%n2i*n3p,rho_ias(1,1,1,ik),1,v_ias(1,1,1),1)
        !!$$if (iproc==0) call yaml_comment('vcopy done')
        !Partial potential term for each partial density
!        if (iproc == 0 .and. verbose > 1) then
!           write(*,*)'Poisson Solver application: orbitals (virt,occ):',iorba,iorbi
!        end if

        !!! For large hgrids, the code freezes during the H_potential
        call H_potential('D',pkernel,v_ias(1,1,1),rho_ias,ehart,0.0_dp,.false.,&
             quiet='YES')
!        if (iproc ==0) write(*,*) 'ehart',ehart*2.0_gp
     end if
     !!$$if (iproc==0) call yaml_comment('Hartree part computed for imulti')

     !After the Poisson Solver we can calculate the lower triangular part of K (and Kaux if nspin=1)
     jk=0
     loop_j: do jmulti=1,imulti
        !!$$if (iproc==0) call yaml_comment(trim(yaml_toa(jmulti,fmt='(i4.4)')))
        !Calculate the orbital indexes for the second transition of the coupling matrix element.
        !We use the same multi index as above for the first transtion.
        jorbi=(jmulti-1)/orbsvirt%norb+1 !index of the occupied state
        jorba=jmulti-(jorbi-1)*orbsvirt%norb !index of the virtual state

        !Check the spin of the orbitals considered (needed for the exchange-correlation part).
        !The same convention as above is used
        if (orbsocc%spinsgn(jorbi) == orbsvirt%spinsgn(jorba)) then
           jk=jk+1
           if (orbsocc%spinsgn(jorbi) == 1.0_gp) then
              jspin=1
           else if (orbsocc%spinsgn(jorbi) == -1.0_gp) then
              jspin=2
           end if
        else
           cycle loop_j
        end if

        !Calculation of the Hartree part of the coupling matrix K_H(ik,jk)
        !K_H(ik,jk) = \int V_Hartree(x,y,z) \rho_bj(x,y,z) dx dy dz
        !Add a self-exchange term in the case of TDHF under TDA.
        if (dofH) then

           !Multiplication of the RPA part
           K(ik,jk)=hxh*hyh*hzh*&
                dot(lr%d%n1i*lr%d%n2i*n3p,rho_ias(1,1,1,jk),1,v_ias(1,1,1),1)

           !In the non spin-pol case, the RPA part is the same for spin off diagonal
           if (nspin ==1) then
              Kaux(ik,jk)=K(ik,jk)
           end if


           !Add the self-exchange term (needed when one is doing TDHF with TDA)
           if (exc_fac .ne. 0.0 .and. ispin==jspin) then

              do i3p=1,n3p !loop over z
                 ind3=(i3p-1)*lr%d%n1i*lr%d%n2i !z-index
                 z=real(i3p+i3s-nbl3-1,gp)*hzh-chargec(3) !z (origin at the center of charge)

                 do i2=1,lr%d%n2i !loop over y
                    ind2=(i2-1)*lr%d%n1i+ind3 !y-index
                    y=real(i2-nbl2-1,gp)*hyh-chargec(2) !y (origin at the center of charge)

                    do i1=1,lr%d%n1i !loop over x
                       x=real(i1-nbl1-1,gp)*hxh-chargec(1) !x (origin at the center of charge)

                       indi=i1+ind2+(iorbi-1)*lr%d%n1i*lr%d%n2i*n3p !multiindex giving the right index for the occupied state wavefunction
                       indj=i1+ind2+(jorbi-1)*lr%d%n1i*lr%d%n2i*n3p !multiindex giving the right index for the other occupied state wavefunction
                       inda=i1+ind2+(iorba-1)*lr%d%n1i*lr%d%n2i*n3p !multiindex giving the right index for the virtual state wavefunction
                       indb=i1+ind2+(jorba-1)*lr%d%n1i*lr%d%n2i*n3p !multiindex giving the right index for the other virtual state wavefunction
                       rhov_ij(i1,i2,i3p)=hfac*psirocc(indj)*psirocc(indi) !partial density, multiplied by the integration factor hfac
                       rho_ab(i1,i2,i3p)=hfac*psivirtr(inda)*psivirtr(indb) !partial density, multiplied by the integration factor hfac

                    end do !loop over x
                 end do !loop over y
              end do !loop over z

              call H_potential('D',pkernel,rhov_ij(1,1,1),rhov_ij,ehart,0.0_dp,.false.,&
                   quiet='YES')

              K_u_HF=hxh*hyh*hzh*&
                   dot(lr%d%n1i*lr%d%n2i*n3p,rho_ab(1,1,1),1,rhov_ij(1,1,1),1) !exact-exchange term of the coupling matrix

              K(ik,jk) = K(ik,jk) - exc_fac * K_u_HF

           end if


        end if

        !Add the XC contribution
        !Map the spin couples to the index of dvxcdrho, in the abinit convention
        if (dofxc) then
           spinindex=ispin+jspin-1

           do i3p=1,n3p
              do i2=1,lr%d%n2i
                 do i1=1,lr%d%n1i
                    K(ik,jk)=K(ik,jk)+hxh*hyh*hzh*&
                         rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                         dvxcdrho(i1,i2,i3p,spinindex)
                    !!begin test
                    !if(iproc==0 .and. ik==2 .and. jk==1 .and. i1==lr%d%n1i/2 .and. i2==lr%d%n2i/2) then
                    !   write(*,'(5(i4),1x,e16.9e2)') iproc, spinindex, i1, i2, i3p, dvxcdrho(i1,i2,i3p,spinindex)
                    !   write(*,'(5(i4),1x,e16.9e2)') iproc, 2, i1, i2, i3p, dvxcdrho(i1,i2,i3p,2)
                    !end if
                    !if(iproc==2 .and. ik==2 .and. jk==1 .and. i1==lr%d%n1i/2 .and. i2==lr%d%n2i/2) then
                    !   write(*,'(5(i4),1x,e16.9e2)') iproc, spinindex, i1, i2, i3p, dvxcdrho(i1,i2,i3p,spinindex)
                    !   write(*,'(5(i4),1x,e16.9e2)') iproc, 2, i1, i2, i3p, dvxcdrho(i1,i2,i3p,2)
                    !end if
                    !!end test
                 end do
              end do
           end do

           !Calculate the spin off-diagonal term for nspin==1
           if (nspin ==1) then
              spinindex=2
              do i3p=1,n3p
                 do i2=1,lr%d%n2i
                    do i1=1,lr%d%n1i
                       Kaux(ik,jk)=Kaux(ik,jk)+hxh*hyh*hzh*&
                            rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                            dvxcdrho(i1,i2,i3p,spinindex)
                       !!begin test
                       !if(iproc==2 .and. ik==2 .and. jk==1 .and. i1==lr%d%n1i/2 .and. i2==lr%d%n2i/2)&
                       !   write(*,*) spinindex, i1, i2, i3p, dvxcdrho(i1,i2,i3p,spinindex)
                       !!end test
                    end do
                 end do
              end do
           end if

        end if

        !add factors from energy occupation numbers (for non-tda case)
        !If full TDDFT, then multiply the coupling matrix element by the "2*sqrt(\omega_q*\omega_{q'})" coefficient of eq. 2.
        if (tddft_approach=='full') then
           K(ik,jk)=K(ik,jk)*2.0_wp*sqrt(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))*&
                sqrt(orbsvirt%eval(jorba)-orbsocc%eval(jorbi))
           if (nspin ==1) then
              Kaux(ik,jk)=Kaux(ik,jk)*2.0_wp*sqrt(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))*&
                   sqrt(orbsvirt%eval(jorba)-orbsocc%eval(jorbi))
           end if
        end if

     end do loop_j
  end do loop_i
  !!$$if (iproc==0) call yaml_sequence_close()

  !!$$if (iproc==0) call yaml_comment('Computation of all the different coupling matrix elements done')

  !Deallocations
  if (exc_fac .ne. 0.0 .and. ispin==jspin) then
     call f_free(rhov_ij)
     call f_free(rho_ab)
  end if
  call f_free(v_ias)
  call f_free(rho_ias)

  !If more than one processor, then perform the MPI_all_reduce of K (and of Kaux if nspin=1) and of dipoles.
  if (nproc > 1) then
     call fmpi_allreduce(K,FMPI_SUM,comm=bigdft_mpi%mpi_comm)
     if (nspin ==1) call fmpi_allreduce(Kaux,FMPI_SUM,comm=bigdft_mpi%mpi_comm)
     call fmpi_allreduce(dipoles(1,1),3*nmulti,FMPI_SUM,comm=bigdft_mpi%mpi_comm)
  end if

  !Copy the values of the dipoles in the second part of the array, given the occupation number of each orbital.
  if (nspin==1) then
     call vcopy(3*nmulti,dipoles(1,1),1,dipoles(1,nmulti+1),1)
  end if

 
  call dscal(3*ndipoles,hxh*hyh*hzh,dipoles(1,1),1)

  !Build the upper triangular part of K (and Kaux if nspin=1)
  do imulti=1,nmulti
     do jmulti=imulti+1,nmulti
        K(imulti,jmulti)=K(jmulti,imulti)
     end do
  end do
  if (nspin==1) then
     do imulti=1,nmulti
        do jmulti=imulti+1,nmulti
           Kaux(imulti,jmulti)=Kaux(jmulti,imulti)
        end do
     end do
  end if

  !Add the diagonal part: {\omega_i}^{1+ntda} \delta_{i,j}
  if (dodiag) then
     !Loop over the transitions
     ik=0 !Initialise the counter of KS transitions
     loop_i2: do imulti=1,orbsvirt%norb*orbsocc%norb
        !Calculate the orbital index
        iorbi=(imulti-1)/orbsvirt%norb+1 !occ. state index
        iorba=imulti-(iorbi-1)*orbsvirt%norb !virt. state index
        !Check the spin of the considered orbitals to increment the counter
        if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
           ik=ik+1
        else
           cycle loop_i2
        end if
        !If both states have the same spin, then add the diagonal term (with the right power 1+ntda)
        K(ik,ik)=K(ik,ik)+(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))**(ntda+1)
     end do loop_i2
  end if

!  !print out the result
!  ik=0
!  loop_i3: do imulti=1,orbsvirt%norb*orbsocc%norb
!     !calculate the orbital index
!     iorbi=(imulti-1)/orbsvirt%norb+1
!     iorba=imulti-(iorbi-1)*orbsvirt%norb
!     if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
!        ik=ik+1
!     else
!        cycle loop_i3
!     end if
!     if (iproc == 0) write(*,*) 'iorba,iorbi,K',iorba,iorbi,K(ik,ik)
!  end do loop_i3


  !Construction of the whole coupling matrix Kbig for nspin=1
  if (nspin == 1) then
     Kbig = f_malloc0((/ 2*nmulti, 2*nmulti /),id='Kbig')

     do ik=1,nmulti
        do jk=1,nmulti
           Kbig(ik,jk)=K(ik,jk)
           Kbig(ik+nmulti,jk+nmulti)=K(ik,jk)
           Kbig(ik+nmulti,jk)=Kaux(ik,jk)
           Kbig(ik,jk+nmulti)=Kaux(ik,jk)
        end do
     end do
  end if

  !Allocate the oscillator strengths
  fi = f_malloc((/ 3, ndipoles /),id='fi')
  !Allocate the excitation energy vector (actually \omega^2 when full TDDFT)
  omega = f_malloc(ndipoles,id='omega')

  lwork = 3*ndipoles !safe value
  work = f_malloc(lwork,id='work')

  !begin test: print out the matrix elements of K
  if (nspin==1) then
     do jmulti = 1, ndipoles
        fsumrule_test=0.0
        do imulti = 1, ndipoles
           !if (iproc==0) write(*,*) jmulti, imulti, Kbig(jmulti, imulti)
           fsumrule_test=fsumrule_test+Kbig(jmulti, imulti)**2
        end do
        !if (iproc==0) write(*,*) jmulti, fsumrule_test
     end do
  else
     do jmulti = 1, ndipoles
        fsumrule_test=0.0
        do imulti = 1, ndipoles
           !if (iproc==0) write(*,*) jmulti, imulti, K(jmulti, imulti)
           fsumrule_test=fsumrule_test+K(jmulti, imulti)**2
        end do
        !if (iproc==0) write(*,*) jmulti, fsumrule_test
     end do
  end if
  !end test

  call f_savetxt('old_coupling_matrix.txt',Kbig) ! we are in the nspin=1 case

  !Find the excitattion energy.
  if (nspin==1) then
     call DSYEV('V','U',ndipoles,Kbig,ndipoles,omega,work,lwork,info)
  else
     call DSYEV('V','U',ndipoles,K   ,ndipoles,omega,work,lwork,info)
  end if
  if (info /= 0) then
     call yaml_warning('Error, DSYEV' // trim(yaml_toa(info)))
  end if

  !If the problem is about full TDDFT, then rewrite the values so that it is the actual excitation energy
  !(and not the excitation energy squared).
  if (tddft_approach=='full') then
     do imulti=1,ndipoles
        omega(imulti)=sqrt(omega(imulti))
     end do
  end if

  !Calculation of the transition dipoles
  !Note that, after DSYEV, K represents the coefficients of the KS transition expansion of the true excitations.
  if (nspin==1) then
     call gemm('N','N',3,ndipoles,ndipoles,1.0_wp,dipoles(1,1),3,&
          Kbig(1,1),ndipoles,0.0_wp,fi(1,1),3)
  else
     call gemm('N','N',3,ndipoles,ndipoles,1.0_wp,dipoles(1,1),3,&
          K(1,1),ndipoles,0.0_wp,fi(1,1),3)
  end if

  !Summary of the results and pretty printing
  if (iproc == 0) then

     if (tddft_approach=='TDA') call yaml_comment('TAMM-DANCOFF APPROXIMATION',hfill='-')
     if (tddft_approach=='full') call yaml_comment('FULL TDDFT',hfill='-')

     !Print the excitation energies and Oscillator strength
     !A yaml output is done as well as a 'tddft_spectrum.txt' file is written
     call yaml_sequence_open('Excitation Energy and Oscillator Strength')
     open(unit=9, file='tddft_spectrum.txt')

     write(9,*) ndipoles

     do imulti = 1, ndipoles
        call yaml_sequence(trim(yaml_toa((/ Ha_eV*omega(imulti),&
              omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2) /),&
              fmt='(1pe10.3)')),advance='no')
        call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))
        write(9,'(f16.12,5x,e16.9e2)') Ha_eV*omega(imulti),&
             omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
     end do

     call yaml_sequence_close()
     close(unit=9)
!!        do imulti = 1, ndipoles
!!           write(9,'(f9.4,5x,1pe10.3)') Ha_eV*sqrt(omega(imulti)),&
!!                sqrt(omega(imulti))*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
!!        end do
!!     end if
!     close(unit=9)


     !Print the excitation energies, and the information about each KS transition sufficiently involved in it:
     !the index of the occupied and virtual state, and the absolute value of the coefficient associated to it.
     !A yaml output is done as well as a 'transitions.txt' file is written
     call yaml_sequence_open('Transition energies (eV)')
     open(unit=10, file='transitions.txt') !This will be used to plot further outputs.

     !Count the number of KS transitions needed to account for all true excitations
     !This number is written in the first line of transitions.txt.
     ik=0 !initialization of the counter of transitions
     !if (nspin==1) then
     !   do imulti = 1, nmulti
     !      !jk=0 !initialization of the counter of allowed transitions
     !      do jmulti = 1, nmulti!orbsocc%norb*orbsvirt%norb
     !         !jorbi=(jmulti-1)/orbsvirt%norb+1
     !         !jorba=jmulti-(jorbi-1)*orbsvirt%norb
     !         !if (orbsocc%spinsgn(jorbi) == orbsvirt%spinsgn(jorba)) then
     !         !   jk=jk+1
     !            if ( abs(Kbig(imulti       ,jmulti       )) > 5.d-02 ) ik=ik+1 !We chose a minimal value for the transition to be taken into account
     !            if ( abs(Kbig(imulti+nmulti,jmulti       )) > 5.d-02 ) ik=ik+1 !We chose a minimal value for the transition to be taken into account
     !            if ( abs(Kbig(imulti       ,jmulti+nmulti)) > 5.d-02 ) ik=ik+1 !We chose a minimal value for the transition to be taken into account
     !            if ( abs(Kbig(imulti+nmulti,jmulti+nmulti)) > 5.d-02 ) ik=ik+1 !We chose a minimal value for the transition to be taken into account
     !         !else
     !         !   cycle
     !         !end if
     !      end do
     !   end do
     !else
     !   do imulti = 1, nmulti
     !      !jk=0 !initialization of the counter of allowed transitions
     !      do jmulti = 1, nmulti !orbsocc%norb*orbsvirt%norb
     !         !jorbi=(jmulti-1)/orbsvirt%norb+1
     !         !jorba=jmulti-(jorbi-1)*orbsvirt%norb
     !         !if (orbsocc%spinsgn(jorbi) == orbsvirt%spinsgn(jorba)) then
     !         !   jk=jk+1
     !            if ( abs(K(imulti,jmulti)) > 5.d-02 ) ik=ik+1
     !         !else
     !         !   cycle
     !         !end if
     !      end do
     !   end do
     !end if
     if (nspin==1) then
        do imulti = 1, ndipoles
           do jmulti = 1, ndipoles
              if ( abs(Kbig(imulti,jmulti)) > 5.d-02 ) ik=ik+1 !We chose a minimal value for the transition to be taken into account
           end do
        end do
     else
        do imulti = 1, ndipoles
           do jmulti = 1, ndipoles
              if ( abs(K(imulti,jmulti)) > 5.d-02 ) ik=ik+1 !We chose a minimal value for the transition to be taken into account
           end do
        end do
     end if

     write(10,*) ik !number of transitions written

     !Then we loop ever the matrix elements of K
     !(K is now containing the coefficients of the KS transitions reproducing the true excitations.)
     if (nspin==1) then

        do imulti=1,nmulti
           call yaml_sequence(advance='no')
           call yaml_sequence_open(advance='no',flow=.true.)
           call yaml_map('Energy',trim(yaml_toa(Ha_eV*omega(imulti),fmt='(f10.5)')))
           ik=0 !transition counter when nspin=1
           do jorbi = 1, orbsocc%norb
              do jorba = 1, orbsvirt%norb
                 jmulti =  (jorbi-1)*orbsvirt%norb+ jorba
                 if (abs(Kbig(jmulti,imulti)) > 5.d-02) then
                    if (ik /= 0) call yaml_newline()
                    ik = ik + 1
                    call yaml_mapping_open(flow=.true.)
                       call yaml_map('Transition',trim(yaml_toa((/ jorbi, jorba /))))
                       call yaml_map('Coeff',trim(yaml_toa(abs(Kbig(jmulti,imulti))**2,fmt='(1pe10.3)')))
                    call yaml_mapping_close()
                    write(10,'(f16.12,2(2x,i4,2x,E16.9E2),2(2x,E16.9E2))') Ha_eV*omega(imulti), jorbi,&
                            orbsocc%eval(jorbi), jorba, orbsvirt%eval(jorba), abs(Kbig(jmulti,imulti)),&
                            omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
                 end if
                 if (abs(Kbig(jmulti+nmulti,imulti)) > 5.d-02) then
                    if (ik /= 0) call yaml_newline()
                    ik = ik + 1
                    call yaml_mapping_open(flow=.true.)
                       call yaml_map('Transition',trim(yaml_toa((/ jorbi+orbsocc%norb, jorba+orbsvirt%norb /))))
                       call yaml_map('Coeff',trim(yaml_toa(abs(Kbig(jmulti+nmulti,imulti))**2,fmt='(1pe10.3)')))
                    call yaml_mapping_close()
                    write(10,'(f16.12,2(2x,i4,2x,E16.9E2),2(2x,E16.9E2))') Ha_eV*omega(imulti), jorbi,&
                            orbsocc%eval(jorbi), jorba, orbsvirt%eval(jorba), abs(Kbig(jmulti+nmulti,imulti)),&
                            omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
                 end if
              end do
           end do
           call yaml_sequence_close(advance='no')
           call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))
        end do
        do imulti=1,nmulti
           call yaml_sequence(advance='no')
           call yaml_sequence_open(advance='no',flow=.true.)
           call yaml_map('Energy',trim(yaml_toa(Ha_eV*omega(imulti+nmulti),fmt='(f10.5)')))
           do jorbi = 1, orbsocc%norb
              do jorba = 1, orbsvirt%norb
                 jmulti =  (jorbi-1)*orbsvirt%norb+ jorba
                 if (abs(Kbig(jmulti,imulti+nmulti)) > 5.d-02) then
                    if (ik /= 0) call yaml_newline()
                    ik = ik + 1
                    call yaml_mapping_open(flow=.true.)
                       call yaml_map('Transition',trim(yaml_toa((/ jorbi, jorba /))))
                       call yaml_map('Coeff',trim(yaml_toa(abs(Kbig(jmulti,imulti+nmulti))**2,fmt='(1pe10.3)')))
                    call yaml_mapping_close()
                    write(10,'(f16.12,2(2x,i4,2x,E16.9E2),2(2x,E16.9E2))') Ha_eV*omega(imulti+nmulti), jorbi,&
                            orbsocc%eval(jorbi), jorba, orbsvirt%eval(jorba), abs(Kbig(jmulti,imulti+nmulti)),&
                            omega(imulti+nmulti)*(2.0_gp/3.0_gp)*&
                            (fi(1,imulti+nmulti)**2+fi(2,imulti+nmulti)**2+fi(3,imulti+nmulti)**2)
                 end if
                 if (abs(Kbig(jmulti+nmulti,imulti+nmulti)) > 5.d-02) then
                    if (ik /= 0) call yaml_newline()
                    ik = ik + 1
                    call yaml_mapping_open(flow=.true.)
                       call yaml_map('Transition',trim(yaml_toa((/ jorbi+orbsocc%norb, jorba+orbsvirt%norb /))))
                       call yaml_map('Coeff',trim(yaml_toa(abs(Kbig(jmulti+nmulti,imulti+nmulti))**2,fmt='(1pe10.3)')))
                    call yaml_mapping_close()
                    write(10,'(f16.12,2(2x,i4,2x,E16.9E2),2(2x,E16.9E2))') Ha_eV*omega(imulti+nmulti), jorbi,&
                            orbsocc%eval(jorbi), jorba, orbsvirt%eval(jorba), abs(Kbig(jmulti+nmulti,imulti+nmulti)),&
                            omega(imulti+nmulti)*(2.0_gp/3.0_gp)*&
                            (fi(1,imulti+nmulti)**2+fi(2,imulti+nmulti)**2+fi(3,imulti+nmulti)**2)
                 end if
              end do
           end do
           call yaml_sequence_close(advance='no')
           call yaml_comment(trim(yaml_toa(imulti+nmulti,fmt='(i4.4)')))
        end do

     else

        do imulti=1,ndipoles
           call yaml_sequence(advance='no')
           call yaml_sequence_open(advance='no',flow=.true.)
           call yaml_map('Energy',trim(yaml_toa(Ha_eV*omega(imulti),fmt='(f10.5)')))
           ik=0 !identifier of the need of a new line when nspin=2
           !We have to check the spin sign, and therefore use a new transition counter.
           jk=0 !transition counter
           do jmulti = 1, orbsocc%norb*orbsvirt%norb
              jorbi=(jmulti-1)/orbsvirt%norb+1
              jorba=jmulti-(jorbi-1)*orbsvirt%norb
              if (orbsocc%spinsgn(jorbi) == orbsvirt%spinsgn(jorba)) then
                 jk=jk+1
                 if (abs(K(jk,imulti)) > 5.d-02) then
                    if (ik /= 0) call yaml_newline()
                    ik = ik + 1
                    call yaml_mapping_open(flow=.true.)
                       call yaml_map('Transition',trim(yaml_toa((/ jorbi, jorba /))))
                       call yaml_map('Coeff',trim(yaml_toa(abs(K(jk,imulti))**2,fmt='(1pe10.3)')))
                    call yaml_mapping_close()
                    write(10,'(f16.12,2(2x,i4,2x,E16.9E2),2(2x,E16.9E2))') Ha_eV*omega(imulti), jorbi,&
                            orbsocc%eval(jorbi), jorba, orbsvirt%eval(jorba), abs(K(jk,imulti)),&
                           &omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
                 end if
              !else
              !   cycle
              end if
           end do
           call yaml_sequence_close(advance='no')
           call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))
        end do

     end if


!     do imulti = 1,nmulti
!        call yaml_sequence(advance='no')
!        call yaml_sequence_open(advance='no',flow=.true.)
!!        if (tddft_approach=='TDA') call yaml_map('Energy',trim(yaml_toa(Ha_eV*omega(imulti),fmt='(f10.5)')))
!!        if (tddft_approach=='full') call yaml_map('Energy',trim(yaml_toa(Ha_eV*sqrt(omega(imulti)),fmt='(f10.5)')))
!        call yaml_map('Energy',trim(yaml_toa(Ha_eV*omega(imulti),fmt='(f10.5)')))
!        ik=0 !transition counter when nspin=1, identifier of the need of a new line when nspin=2
!
!        if (nspin==1) then
!
!           !No need to check the spin sign
!           do iorbi = 1, orbsocc%norb
!              do iorba = 1, orbsvirt%norb
!                 jmulti =  (iorbi-1)*orbsvirt%norb+ iorba
!                 if (abs(Kbig(jmulti,imulti)) > 5.d-02) then
!                    if (ik /= 0) call yaml_newline()
!                    ik = ik + 1
!                    call yaml_mapping_open(flow=.true.)
!                       call yaml_map('Transition',trim(yaml_toa((/ iorbi, iorba /))))
!                       call yaml_map('Coeff',trim(yaml_toa(abs(Kbig(jmulti,imulti))**2,fmt='(1pe10.3)')))
!                    call yaml_mapping_close()
!                    write(10,'(f16.12,2(2x,i4,2x,E16.9E2),2(2x,E16.9E2))') Ha_eV*omega(imulti), iorbi,&
!                            orbsocc%eval(iorbi), iorba, orbsvirt%eval(iorba), abs(Kbig(jmulti,imulti)),&
!                           &omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
!                 end if
!                 !We must consider the fact that ndipoles=2*nmulti, so we explicitly have to search in the second half of the Kbig-matrix
!                 if (abs(Kbig(jmulti+nmulti,imulti)) > 5.d-02) then
!                    if (ik /= 0) call yaml_newline()
!                    ik = ik + 1
!                    call yaml_mapping_open(flow=.true.)
!                       call yaml_map('Transition',trim(yaml_toa((/ iorbi, iorba /))))
!                       call yaml_map('Coeff',trim(yaml_toa(abs(Kbig(jmulti+nmulti,imulti))**2,fmt='(1pe10.3)')))
!                    call yaml_mapping_close()
!                    write(10,'(f16.12,2(2x,i4,2x,E16.9E2),2(2x,E16.9E2))') Ha_eV*omega(imulti), iorbi,&
!                            orbsocc%eval(iorbi), iorba, orbsvirt%eval(iorba), abs(Kbig(jmulti+nmulti,imulti)),&
!                           &omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
!                 end if
!                 if (abs(Kbig(jmulti,imulti+nmulti)) > 5.d-02) then
!                    if (ik /= 0) call yaml_newline()
!                    ik = ik + 1
!                    call yaml_mapping_open(flow=.true.)
!                       call yaml_map('Transition',trim(yaml_toa((/ iorbi, iorba /))))
!                       call yaml_map('Coeff',trim(yaml_toa(abs(Kbig(jmulti,imulti+nmulti))**2,fmt='(1pe10.3)')))
!                    call yaml_mapping_close()
!                    write(10,'(f16.12,2(2x,i4,2x,E16.9E2),2(2x,E16.9E2))') Ha_eV*omega(imulti+nmulti), iorbi,&
!                            orbsocc%eval(iorbi), iorba, orbsvirt%eval(iorba), abs(Kbig(jmulti,imulti+nmulti)),&
!                           &omega(imulti+nmulti)*(2.0_gp/3.0_gp)*(fi(1,imulti+nmulti)**2+fi(2,imulti+nmulti)**2+&
!                           &fi(3,imulti+nmulti)**2)
!                 end if
!                 if (abs(Kbig(jmulti+nmulti,imulti+nmulti)) > 5.d-02) then
!                    if (ik /= 0) call yaml_newline()
!                    ik = ik + 1
!                    call yaml_mapping_open(flow=.true.)
!                       call yaml_map('Transition',trim(yaml_toa((/ iorbi, iorba /))))
!                       call yaml_map('Coeff',trim(yaml_toa(abs(Kbig(jmulti+nmulti,imulti+nmulti))**2,fmt='(1pe10.3)')))
!                    call yaml_mapping_close()
!                    write(10,'(f16.12,2(2x,i4,2x,E16.9E2),2(2x,E16.9E2))') Ha_eV*omega(imulti+nmulti), iorbi,&
!                            orbsocc%eval(iorbi), iorba, orbsvirt%eval(iorba), abs(Kbig(jmulti+nmulti,imulti+nmulti)),&
!                           &omega(imulti+nmulti)*(2.0_gp/3.0_gp)*(fi(1,imulti+nmulti)**2+fi(2,imulti+nmulti)**2+&
!                           &fi(3,imulti+nmulti)**2)
!                 end if
!              end do
!           end do
!
!        else
!
!           !We have to check the spin sign, and therefore use a new transition counter.
!           jk=0 !transition counter
!           do jmulti = 1, orbsocc%norb*orbsvirt%norb
!              jorbi=(jmulti-1)/orbsvirt%norb+1
!              jorba=jmulti-(jorbi-1)*orbsvirt%norb
!              if (orbsocc%spinsgn(jorbi) == orbsvirt%spinsgn(jorba)) then
!                 jk=jk+1
!                 if (abs(K(jk,imulti)) > 5.d-02) then
!                    if (ik /= 0) call yaml_newline()
!                    ik = ik + 1
!                    call yaml_mapping_open(flow=.true.)
!                       call yaml_map('Transition',trim(yaml_toa((/ jorbi, jorba /))))
!                       call yaml_map('Coeff',trim(yaml_toa(abs(K(jk,imulti))**2,fmt='(1pe10.3)')))
!                    call yaml_mapping_close()
!                    write(10,'(f16.12,2(2x,i4,2x,E16.9E2),2(2x,E16.9E2))') Ha_eV*omega(imulti), jorbi,&
!                            orbsocc%eval(jorbi), jorba, orbsvirt%eval(jorba), abs(K(jk,imulti)),&
!                           &omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
!                 end if
!              else
!                 cycle
!              end if
!           end do
!
!        end if
!
!        call yaml_sequence_close(advance='no')
!        call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))
!
!     end do !imulti

     call yaml_sequence_close()
     close(unit=10)

  end if !iproc=0

  !Deallocations
  call f_free(omega)
  call f_free(work)
  call f_free(fi)
  call f_free(K)
  call f_free(dipoles)
  if (nspin ==1 ) then
     call f_free(Kaux)
     call f_free(Kbig)
  end if

  call f_release_routine()

END SUBROUTINE coupling_matrix_prelim
