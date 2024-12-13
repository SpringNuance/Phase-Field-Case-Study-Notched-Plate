!====================================================================
!          Program for phase field damage coupled with 
!          mechanical loading and hydrogen diffusion 
!          Mechanical model: standard Hooke's law elasticity 
!                            isotropic Von Mises plasticity
!          Hydrogen diffusion model: Fick's second law
!          Damage model: phase field damage model
!          by Nguyen Xuan Binh
!          binh.nguyen@aalto.fi
!          July 2024, Abaqus 2023
!          Model extended from the work of Emilio Martinez-Paneda
!          Paper: A phase field formulation for hydrogen assisted cracking
!          Computer Methods in Applied Mechanics and Engineering 342: 742-761 
!          (2018) doi: 10.1016/j.cma.2018.07.021
!          DO NOT DISTRIBUTE WITHOUT AUTHOR'S PERMISSION
!====================================================================

!     State variables  
    ! 1 to 6: S11, S22, S33, S12, S13, S23
    ! 7 to 12: E11, E22, E33, E12, E13, E23
    ! "1, AR1_sig11, AR1_sig11   ",
    ! "2, AR2_sig22, AR2_sig22   ",
    ! "3, AR3_sig33, AR3_sig33   ",
    ! "4, AR4_sig12, AR4_sig12   ",
    ! "5, AR5_sig13, AR5_sig13   ",
    ! "6, AR6_sig23, AR6_sig23   ",
    ! "7, AR7_eps11, AR7_eps11   ",
    ! "8, AR8_eps22, AR8_eps22   ",
    ! "9, AR9_eps33, AR9_eps33   ",
    ! "10, AR10_eps12, AR10_eps12   ",
    ! "11, AR11_eps13, AR11_eps13   ",
    ! "12, AR12_eps23, AR12_eps23   ",     
    ! "13, AR13_eqplas, AR13_eqplas   ",
    ! "14, AR14_sig_H, AR14_sig_H   ",
    ! "15, AR15_sig_vonMises, AR15_sig_vonMises   ",
    ! "16, AR16_triax, AR16_triax   ",
    ! "17, AR17_lode, AR17_lode   ",
    ! "18, AR18_phi, AR18_phi   ",
    ! "19, AR19_history, AR19_history   ",
    ! "20, AR20_Gc, AR20_Gc   ",
    ! "21, AR21_theta_coverage, AR21_theta_coverage   ",
    ! "22, AR22_C_mol, AR22_C_mol   ",	
    ! "23, AR23_CL_mol, AR23_CL_mol   ", 
    ! "24, AR24_CT_mol, AR24_CT_mol   ", 

!***********************************************************************

module precision
    use iso_fortran_env
    integer, parameter :: dp = real64
end module precision

!***********************************************************************

module common_block
    use precision
    implicit none

    real(kind=dp), parameter :: molar_mass_H = 1.00784d0 ! g/mol
    real(kind=dp), parameter :: molar_mass_Fe = 55.845d0 ! g/mol
    real(kind=dp), parameter :: ratio_molar_mass_Fe_H = 55.415d0
    real(kind=dp), parameter :: density_metal = 7900.0d0 ! kg/m^3
    real(kind=dp), parameter :: avogadro = 6.022e23 ! 1/mol

    integer, parameter :: start_mech_props_index = 1 ! Index of the first mechanical property in props
    integer, parameter :: start_phase_props_idx = 9 ! Index of the first phase field property in props
    integer, parameter :: start_hydro_props_idx = 17 ! Index of the first hydrogen diffusion property in props
    integer, parameter :: start_flow_props_idx = 41 ! Index of the first flow curve data in props
    
    integer, parameter :: eqplas_idx = 13 ! Index of the equivalent plastic strain in statev
    integer, parameter :: sig_H_idx = 14 ! Index of the hydrogen concentration in statev
    integer, parameter :: sig_vonMises_idx = 15 ! Index of the equivalent von Mises stress in statev
    integer, parameter :: triax_idx = 16 ! Index of the triaxiality in statev
    integer, parameter :: lode_idx = 17 ! Index of the Lode parameter in statev
    integer, parameter :: phi_idx = 18 ! Index of the phase field parameter in statev
    integer, parameter :: history_idx = 19 ! Index of the history variable in statev
    integer, parameter :: Gc_idx = 20 ! Index of the fracture energy in statev
    integer, parameter :: theta_coverage_idx = 21 ! Index of the hydrogen coverage in statev
    integer, parameter :: C_mol_idx = 22 ! Index of the hydrogen concentration in the matrix in statev
    integer, parameter :: CL_mol_idx = 23 ! Index of the hydrogen concentration in the lattice in statev
    integer, parameter :: CT_mol_idx = 24 ! Index of the hydrogen concentration in the trap in statev
    
    ! First dim: maximum number of elements to accomodate varying number of elements when remeshed
    ! Second dim: number of solution state dependent variables (nsvars in UEL and nstatev in UMAT)
    ! Third dim: number of integration points

    real(kind=dp) :: user_vars(100000, 24, 8)
    ! Run it with 1 cpu, and run it with 4 cpus, and compare the results
    ! 
    integer :: adjacency_matrix !
    integer :: nelem ! Storing the actual number of elements
    
    save
    ! The save command is very important. 
    ! It allows the values to be stored and shared between subroutines 
    ! without resetting them to zero every time the subroutine is called

end module

!***********************************************************************

subroutine UEXTERNALDB(lop,lrestart,time,dtime,kstep,kinc)
    use precision
    use common_block
    include 'aba_param.inc' 
    dimension time(2)
    
    ! LOP=0 indicates that the subroutine is being called at the start of the analysis.
    if (lop == 0) then 
        user_vars = 0.0d0
        ! Create that adjacency matrix. 
    end if

return
end

!*****************************************************************
!  8-node     8---------------7
!  brick     /|              /|       zeta (positive)
!           / |  x 7   x 8  / |       
!          5---------------6  |       |     eta (positive)
!          |  | x 5   x 6  |  |       |   /
!          |  |            |  |       |  /
!          |  4------------|--3       | /
!          | /   x 3   x 4 | /        |/
!          |/   x 1   x 2  |/         O--------- xi (positive)
!          1---------------2           origin at cube center
!          
!          Outer number is nodal points
!         Inner number marked with x is intergration (int) points
!
!*****************************************************************

module iso_module

    use precision
    real(kind=dp), parameter :: coord_inter = 1.0d0
    real(kind=dp), parameter :: int_inter = 1.0d0 / sqrt(3.0d0)
    real(kind=dp), parameter :: coord_extra = sqrt(3.0d0)
    real(kind=dp), parameter :: int_extra = 1.0d0

    ! weight is the integration point weight for their shape function contribution
    real(kind=dp), parameter :: weight(8) = (/1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0/)
    
    ! Interpolating coordinates (nodal to int)
    ! Isoparametric coordinates for nodal points in hexahedral 3D element
    real(kind=dp), parameter :: xi_nodal_inter(8)   = (/ -coord_inter,  coord_inter,  coord_inter, -coord_inter, &
                                                         -coord_inter,  coord_inter,  coord_inter, -coord_inter /)
    real(kind=dp), parameter :: eta_nodal_inter(8)  = (/ -coord_inter, -coord_inter,  coord_inter,  coord_inter, &
                                                         -coord_inter, -coord_inter,  coord_inter,  coord_inter /)
    real(kind=dp), parameter :: zeta_nodal_inter(8) = (/ -coord_inter, -coord_inter, -coord_inter, -coord_inter, &
                                                          coord_inter,  coord_inter,  coord_inter,  coord_inter /)

    ! Isoparametric coordinates for integration points in hexahedral 3D element
    real(kind=dp), parameter :: xi_int_inter(8)   = (/ -int_inter,  int_inter, -int_inter,  int_inter, &
                                                         -int_inter,  int_inter, -int_inter,  int_inter /)
    real(kind=dp), parameter :: eta_int_inter(8)  = (/ -int_inter, -int_inter,  int_inter,  int_inter, &
                                                         -int_inter, -int_inter,  int_inter,  int_inter /)
    real(kind=dp), parameter :: zeta_int_inter(8) = (/ -int_inter, -int_inter, -int_inter, -int_inter, &
                                                          int_inter,  int_inter,  int_inter,  int_inter /)


    ! Extrapolating coordinates (int to nodal)
    real(kind=dp), parameter :: xi_nodal_extra(8)   = (/ -coord_extra,  coord_extra,  coord_extra, -coord_extra, &
                                                         -coord_extra,  coord_extra,  coord_extra, -coord_extra /)
    real(kind=dp), parameter :: eta_nodal_extra(8)  = (/ -coord_extra, -coord_extra,  coord_extra,  coord_extra, &
                                                         -coord_extra, -coord_extra,  coord_extra,  coord_extra /)
    real(kind=dp), parameter :: zeta_nodal_extra(8) = (/ -coord_extra, -coord_extra, -coord_extra, -coord_extra, &
                                                          coord_extra,  coord_extra,  coord_extra,  coord_extra /)

    real(kind=dp), parameter :: xi_int_extra(8)   = (/ -int_extra,  int_extra, -int_extra,  int_extra, &
                                                         -int_extra,  int_extra, -int_extra,  int_extra /)
    real(kind=dp), parameter :: eta_int_extra(8)  = (/ -int_extra, -int_extra,  int_extra,  int_extra, &
                                                         -int_extra, -int_extra,  int_extra,  int_extra /)
    real(kind=dp), parameter :: zeta_int_extra(8) = (/ -int_extra, -int_extra, -int_extra, -int_extra, &
                                                          int_extra,  int_extra,  int_extra,  int_extra /)

end module iso_module

subroutine kshapefcn(kintk,ninpt,nnode,ndim,shape_node_to_int,B_deriv_local)
!   
    use iso_module
    include 'aba_param.inc'
!
    dimension shape_node_to_int(1, nnode),B_deriv_local(ndim,nnode)

!   determine (g,h,r)
    f = xi_int_inter(kintk)
    g = eta_int_inter(kintk)
    h = zeta_int_inter(kintk)

!   shape functions
    shape_node_to_int(1,1)=0.125d0 * (1.d0-f) * (1.d0-g) * (1.d0-h)
    shape_node_to_int(1,2)=0.125d0 * (1.d0+f) * (1.d0-g) * (1.d0-h)
    shape_node_to_int(1,3)=0.125d0 * (1.d0+f) * (1.d0+g) * (1.d0-h)
    shape_node_to_int(1,4)=0.125d0 * (1.d0-f) * (1.d0+g) * (1.d0-h)
    shape_node_to_int(1,5)=0.125d0 * (1.d0-f) * (1.d0-g) * (1.d0+h)
    shape_node_to_int(1,6)=0.125d0 * (1.d0+f) * (1.d0-g) * (1.d0+h)
    shape_node_to_int(1,7)=0.125d0 * (1.d0+f) * (1.d0+g) * (1.d0+h)
    shape_node_to_int(1,8)=0.125d0 * (1.d0-f) * (1.d0+g) * (1.d0+h)

!   derivative d(Ni)/d(f)
    B_deriv_local(1,1)=-0.125d0 * (1.d0-g) * (1.d0-h)
    B_deriv_local(1,2)= 0.125d0 * (1.d0-g) * (1.d0-h)
    B_deriv_local(1,3)= 0.125d0 * (1.d0+g) * (1.d0-h)
    B_deriv_local(1,4)=-0.125d0 * (1.d0+g) * (1.d0-h)
    B_deriv_local(1,5)=-0.125d0 * (1.d0-g) * (1.d0+h)
    B_deriv_local(1,6)= 0.125d0 * (1.d0-g) * (1.d0+h)
    B_deriv_local(1,7)= 0.125d0 * (1.d0+g) * (1.d0+h)
    B_deriv_local(1,8)=-0.125d0 * (1.d0+g) * (1.d0+h)

!     derivative d(Ni)/d(g)
    B_deriv_local(2,1)=-0.125d0 * (1.d0-f) * (1.d0-h)
    B_deriv_local(2,2)=-0.125d0 * (1.d0+f) * (1.d0-h)
    B_deriv_local(2,3)= 0.125d0 * (1.d0+f) * (1.d0-h)
    B_deriv_local(2,4)= 0.125d0 * (1.d0-f) * (1.d0-h)
    B_deriv_local(2,5)=-0.125d0 * (1.d0-f) * (1.d0+h)
    B_deriv_local(2,6)=-0.125d0 * (1.d0+f) * (1.d0+h)
    B_deriv_local(2,7)= 0.125d0 * (1.d0+f) * (1.d0+h)
    B_deriv_local(2,8)= 0.125d0 * (1.d0-f) * (1.d0+h)

!     derivative d(Ni)/d(h)
    B_deriv_local(3,1)=-0.125d0 * (1.d0-f) * (1.d0-g)
    B_deriv_local(3,2)=-0.125d0 * (1.d0+f) * (1.d0-g)
    B_deriv_local(3,3)=-0.125d0 * (1.d0+f) * (1.d0+g)
    B_deriv_local(3,4)=-0.125d0 * (1.d0-f) * (1.d0+g)
    B_deriv_local(3,5)= 0.125d0 * (1.d0-f) * (1.d0-g)
    B_deriv_local(3,6)= 0.125d0 * (1.d0+f) * (1.d0-g)
    B_deriv_local(3,7)= 0.125d0 * (1.d0+f) * (1.d0+g)
    B_deriv_local(3,8)= 0.125d0 * (1.d0-f) * (1.d0+g)

return
end


subroutine kjacobian(jelem,ndim,nnode,coords,B_deriv_local,djac,B_deriv_global,mcrd)

!     Notation: djac - Jac determinant; xjaci - inverse of Jac matrix
!     B_deriv_global - shape functions derivatives w.r.t. global coordinates

    include 'aba_param.inc'

    dimension xjac(ndim,ndim),xjaci(ndim,ndim),coords(mcrd,nnode), &
        B_deriv_local(ndim,nnode),B_deriv_global(ndim,nnode)

    xjac=0.d0

    do inode=1,nnode
        do idim=1,ndim
            do jdim=1,ndim
                xjac(jdim,idim)=xjac(jdim,idim)+ &
                    B_deriv_local(jdim,inode) * coords(idim,inode)
            end do
        end do
    end do

    if (ndim==3) then

        djac = xjac(1,1) * xjac(2,2) * xjac(3,3)+xjac(2,1) * xjac(3,2) * xjac(1,3) &
              +xjac(3,1) * xjac(2,3) * xjac(1,2)-xjac(3,1) * xjac(2,2) * xjac(1,3) &
              -xjac(2,1) * xjac(1,2) * xjac(3,3)-xjac(1,1) * xjac(2,3) * xjac(3,2)
        !if (djac>0.d0) then ! jacobian is positive - o.k.
        xjaci(1,1)=(xjac(2,2) * xjac(3,3)-xjac(2,3) * xjac(3,2))/djac
        xjaci(1,2)=(xjac(1,3) * xjac(3,2)-xjac(1,2) * xjac(3,3))/djac
        xjaci(1,3)=(xjac(1,2) * xjac(2,3)-xjac(1,3) * xjac(2,2))/djac
        xjaci(2,1)=(xjac(2,3) * xjac(3,1)-xjac(2,1) * xjac(3,3))/djac
        xjaci(2,2)=(xjac(1,1) * xjac(3,3)-xjac(1,3) * xjac(3,1))/djac
        xjaci(2,3)=(xjac(1,3) * xjac(2,1)-xjac(1,1) * xjac(2,3))/djac
        xjaci(3,1)=(xjac(2,1) * xjac(3,2)-xjac(2,2) * xjac(3,1))/djac
        xjaci(3,2)=(xjac(1,2) * xjac(3,1)-xjac(1,1) * xjac(3,2))/djac
        xjaci(3,3)=(xjac(1,1) * xjac(2,2)-xjac(1,2) * xjac(2,1))/djac
        !else ! negative or zero jacobian
        ! write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
        !endif
    else if (ndim==2) then
        djac=xjac(1,1) * xjac(2,2)-xjac(1,2) * xjac(2,1)
        xjaci(1,1)=xjac(2,2)/djac
        xjaci(2,2)=xjac(1,1)/djac
        xjaci(1,2)=-xjac(1,2)/djac
        xjaci(2,1)=-xjac(2,1)/djac
    endif

    B_deriv_global=matmul(xjaci,B_deriv_local)

return
end

subroutine kbmatrix(B_deriv_global,ntens,nnode,ndim,Bu_matrix)

    !   Notation, strain tensor: e11, e22, e33, e12, e13, e23
    include 'aba_param.inc'

    dimension B_deriv_global(ndim,nnode),Bu_matrix(ntens,nnode * ndim)

    Bu_matrix=0.d0
    do inode=1,nnode
        Bu_matrix(1,ndim * inode-ndim+1)=B_deriv_global(1,inode)
        Bu_matrix(2,ndim * inode-ndim+2)=B_deriv_global(2,inode)
        Bu_matrix(4,ndim * inode-ndim+1)=B_deriv_global(2,inode)
        Bu_matrix(4,ndim * inode-ndim+2)=B_deriv_global(1,inode)
        if (ndim==3) then
            Bu_matrix(3,ndim * inode)=B_deriv_global(3,inode)
            Bu_matrix(5,ndim * inode-2)=B_deriv_global(3,inode)
            Bu_matrix(5,ndim * inode)=B_deriv_global(1,inode)
            Bu_matrix(6,ndim * inode-1)=B_deriv_global(3,inode)
            Bu_matrix(6,ndim * inode)=B_deriv_global(2,inode)
        endif
    end do

return
end

subroutine kstatevar(npt,nsvint,statev,statev_ip,icopy)

!   Transfer data to/from element-level state variable array from/to
!   material-point level state variable array.

    include 'aba_param.inc'

    dimension statev(*),statev_ip(*)

    isvinc=(npt-1) * nsvint     ! integration point increment

    if (icopy==1) then ! Prepare arrays for entry into umat
        do i=1,nsvint
            statev_ip(i)=statev(i+isvinc)
        enddo
    else ! Update element state variables upon return from umat
        do i=1,nsvint
            statev(i+isvinc)=statev_ip(i)
        enddo
    end if

return
end

subroutine UMAT_elastic(props,nprops,ddsdde,stress,dstran,ntens,ndi,nshr,statev)
    
    use precision
    use common_block
    include 'aba_param.inc'
      
!   Subroutine with the material model
      
    dimension props(*),ddsdde(ntens,ntens),stress(ntens),statev(*),dstran(ntens)
    real(kind=dp) :: E, nu, mu, lambda
    
    E = props(1)           ! Young's modulus 
    nu = props(2)          ! Poisson's ratio 

    ! Lame's parameters
    mu = E/(2.0d0 * (1.0d0 + nu))  ! Shear modulus
    lambda = E * nu/((1.0d0 + nu) * (1.0d0 - 2.0d0 * nu)) ! Lame's first constant

    ! print *, 'UMAT_elastic: mu = ', mu
    ! print *, 'UMAT_elastic: lambda = ', lambda
    ! initialize as 0
    ddsdde = 0.0 ! Their unit is Pa
    
    do i = 1, ndi
        do j = 1, ndi
            ddsdde(j, i) = lambda
        end do 
        ddsdde(i,i) = lambda + 2.0d0 * mu
    end do 

    ! Shear contribution
    do i = ndi + 1, ntens
        ddsdde(i,i) = mu
    end do 

    ! Stress increment evaluation
    stress = stress + matmul(ddsdde,dstran)      

return
end

! This is isotropic von Mises plasticity model

subroutine UMAT_von_Mises(props,nprops,ddsdde,stress,dstran,ntens,ndi,nshr,statev)
    
    use precision
    use common_block
    include 'aba_param.inc'
      
    dimension props(*),ddsdde(ntens,ntens),stress(ntens),statev(*), dstran(ntens)

    real(kind=dp) :: E, nu, lambda, mu, eqplas, deqpl, syield, syiel0, sig_Mises, shydro, rhs 
    real(kind=dp) :: effective_mu, effective_lambda, effective_hard    

    dimension eelas(ntens), eplas(ntens), flow(ntens), hard(3)

    real(kind=8), parameter :: toler = 1e-12
    real(kind=8), parameter :: newton = 100
    integer :: nvalue, k_newton

    ! LOCAL ARRAYS
    ! ----------------------------------------------------------------
    ! EELAS - ELASTIC STRAINS
    ! EPLAS - PLASTIC STRAINS
    ! FLOW - DIRECTION OF PLASTIC FLOW
    ! ----------------------------------------------------------------
    
    ! ----------------------------------------------------------------
    ! UMAT FOR ISOTROPIC ELASTICITY AND ISOTROPIC MISES PLASTICITY
    ! CANNOT BE USED FOR PLANE STRESS
    ! ----------------------------------------------------------------
    ! PROPS(1) - E
    ! PROPS(2) - NU
    ! PROPS(3:nprops) - SYIELD AN HARDENING DATA
    ! props(3) - syiel0, props(4) - eqpl0, props(5) - syiel1, props(6) - eqpl1, ...
    ! and props(nprops-1) - SYIELD_N, props(nprops) - EQPLAS_N
    ! CALLS UHARD FOR CURVE OF YIELD STRESS VS. PLASTIC STRAIN
    ! ----------------------------------------------------------------

    ! material properties

    E = props(1)           ! Young's modulus 
    nu = props(2)          ! Poisson's ratio 
    ! props(1:8) is for mechanical properties
    ! props(9:16) is for phase field properties
    ! props(17:40) is for hydrogen diffusion properties
    ! props(17:40) is for hydrogen diffusion properties
    

    ! Lame's parameters
    mu = E/(2.0d0 * (1.0 + nu))  ! Shear modulus
    lambda = E*nu/((1.0 + nu) * (1.0 - 2.0 * nu)) ! Lame's first constant

    ! initialize as 0
    ddsdde = 0.0 ! Their unit is Pa
    
    do i = 1, ndi
        do j = 1, ndi
            ddsdde(j, i) = lambda
        end do 
        ddsdde(i,i) = lambda + 2.0 * mu
    end do 

    ! Shear contribution
    do i = ndi + 1, ntens
        ddsdde(i,i) = mu
    end do 

    ! print hello here
    ! write(7,*) 'Hello from UMAT'
    ! recover elastic and plastic strains and rotate forward
    ! also recover equivalent plastic strain

    call rotsig(statev(      1), drot, eelas, 2, ndi, nshr)
    call rotsig(statev(ntens+1), drot, eplas, 2, ndi, nshr)
    
    eqplas = statev(1+2*ntens)      ! Equivalent plastic strain

    ! calculate predictor stress and elastic strain

    do i = 1, ntens
        do j = 1, ntens
            stress(j) = stress(j) + ddsdde(j, i) * dstran(i)
        end do
        eelas(i) = eelas(i) + dstran(i)
    end do

    ! Calculate equivalent von Mises stress
    
    sig_Mises = (stress(1) - stress(2))**2 + &
                (stress(2) - stress(3))**2 + &
                (stress(3) - stress(1))**2
    
    do i = ndi + 1, ntens
        sig_Mises = sig_Mises + 6.d0 * stress(i)**2
    end do
    sig_Mises = sqrt(sig_Mises/2.d0) ! Unit is Pa
    
    ! get yield stress from the specified hardening curve
    ! nvalue equal to number of points on the hardening curve
    
    nvalue = (nprops - start_flow_props_idx + 1) / 2

    ! print *, 'nvalue = ', nvalue ! 100
    ! print *, 'start_flow_props_idx = ', start_flow_props_idx ! 41
    
    call UHARD_von_Mises(syiel0, hard, eqplas, &
                                statev, nvalue, props(start_flow_props_idx))
    
    ! Determine if active yielding

    if (sig_Mises > (1.d0 + toler) * syiel0) then

        ! actively yielding
        ! separate the hydrostatic from the deviatoric stress
        ! calculate the flow direction

        shydro = (stress(1) + stress(2) + stress(3))/3.d0
        do i=1,ndi
            flow(i) = (stress(i) - shydro)/sig_Mises
        end do
        do i=ndi+1,ntens
            flow(i) = stress(i)/sig_Mises
        end do
        
        ! solve for equivalent von Mises stress and equivalent plastic strain increment 
        ! using Newton-Raphson iteration

        syield = syiel0
        deqpl = 0.d0
        do k_newton = 1, newton
            rhs = sig_Mises - (3.d0 * mu * deqpl) - syield
            deqpl = deqpl + rhs / ((3.d0 * mu) + hard(1))

            call UHARD_von_Mises(syield, hard, eqplas + deqpl, &
                                statev, nvalue, props(start_flow_props_idx))
                                
            if (abs(rhs) < toler * syiel0) exit
        end do

        if (k_newton == newton) write(7,*) 'WARNING: plasticity loop failed'

        ! Update stresses, elastic and plastic strains
        do i = 1, ndi
            stress(i) = flow(i) * syield + shydro
            eplas(i) = eplas(i) + 3.d0/2.d0 * flow(i) * deqpl
            eelas(i) = eelas(i) - 3.d0/2.d0 * flow(i) * deqpl
        end do
        
        do i = ndi + 1, ntens
            stress(i) = flow(i) * syield
            eplas(i) = eplas(i) + 3.d0 * flow(i) * deqpl
            eelas(i) = eelas(i) - 3.d0 * flow(i) * deqpl
        end do

        ! Finally, we update the equivalent plastic strain
        eqplas = eqplas + deqpl

        ! Calculate the plastic strain energy density
        spd = deqpl * (syiel0 + syield) / 2.d0
       
        ! Formulate the jacobian (material tangent)   

        ! effective shear modulus
        effective_mu = mu * syield / sig_Mises 

        ! effective Lame's constant
        effective_lambda = (E/(1.d0 - 2.d0 * nu) - 2.d0 * effective_mu)/3.d0 

        ! effective hardening modulus
        effective_hard = 3.d0 * mu * hard(1)/(3.d0 * mu + hard(1)) - 3.d0 * effective_mu 

        do i = 1, ndi
            do j = 1, ndi
                ddsdde(j,i) = effective_lambda
            end do
            ddsdde(i,i) = 2.d0 * effective_mu + effective_lambda
        end do

        do i = ndi + 1, ntens
            ddsdde(i,i) = effective_mu
        end do

        do i = 1, ntens
            do j = 1, ntens
                ddsdde(j,i) = ddsdde(j,i) + effective_hard * flow(j) * flow(i)
            end do
        end do
    endif

return
end

!***********************************************************************

subroutine UHARD_von_Mises(syield, hard, eqplas, statev, nvalue, table)

    include 'aba_param.inc'

    character*80 cmname
    dimension hard(3),statev(*),table(2, nvalue)
    
    ! Variables to be defined
    ! syield: Yield stress for isotropic plasticity. Yield surface size for combined hardening.
    ! hard(1): Variation of SYIELD with respect to the equivalent plastic strain
    ! hard(2): Variation of SYIELD with respect to the equivalent plastic strain rate
    ! hard(3): Variation of SYIELD with respect to the temperature (only in coupled thermal-mechanical analysis)
    ! statev(nstatv): State variables
     
    ! Variables passed in for information
    ! eqplas: Equivalent plastic strain
    ! nvalue: Number of hardening properties entered for this user-defined hardening definition.
    ! table: Hardening properties entered for this user-defined hardening definition.
    
    ! Note: Fortran srt params are passed by reference, not by values
    ! nvalue here refer to the number of properties entered for the hardening definition
    ! table(N) is the address of the first element of the array of hardening properties
    
    ! table in UMAT is 1D array, and now we view it as a 2D array
    ! UHARD interpret every pair of elements in table as a column pair of (yield stress, eqplas)
    
    ! set yield stress to last value of table, hardening to zero
    
    syield = table(1, nvalue)
    hard(1) = 0.d0
    
    ! we print the first values of the table to see if it is correct
    ! print *, 'table(1, 1): ', table(1, 1)
    ! print *, 'table(2, 1): ', table(2, 1)
    ! print *, 'table(1, 2): ', table(1, 2)
    ! print *, 'table(2, 2): ', table(2, 2)

    ! Now print the last
    ! print *, 'table(1, nvalue): ', table(1, nvalue)
    ! print *, 'table(2, nvalue): ', table(2, nvalue)

    ! if more than one entry, search table
    
    if (nvalue > 1) then
        do k1 = 1, nvalue - 1
            eqpl1 = table(2, k1 + 1)
            if (eqplas < eqpl1) then
                eqpl0 = table(2, k1)
                if (eqpl1 <= eqpl0) then
                    write(7,*) 'error - plastic strain must be entered in ascending order'
                end if

                ! current yield stress and hardening

                deqpl = eqpl1 - eqpl0
                syiel0 = table(1, k1)
                syiel1 = table(1, k1 + 1)
                dsyiel = syiel1 - syiel0
                hard(1) = dsyiel/deqpl
                syield = syiel0 + (eqplas - eqpl0) * hard(1)
                exit
            endif
        end do
    endif

return
end


!***********************************************************************

subroutine UEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars, &
    props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime, &
    kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf, &
    lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njpro,period)

    use precision
    use common_block
    use iso_module
    include 'aba_param.inc' !implicit real(a-h o-z)
      
    dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),svars(*), &
        energy(8),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel), &
        a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*), &
        ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*)

    integer, parameter :: ndim = 3 ! Number of spatial dimensions
    integer, parameter :: ntens = 6 ! Number of stress-strain components
    integer, parameter :: ndi = 3 ! Number of direct stress-strain components
    integer, parameter :: nshr = 3 ! Number of shear stress-strain components
    integer, parameter :: ninpt = 8 ! Number of integration points
    integer, parameter :: nsvint = 24 ! Number of state variables at integration points
    integer, parameter :: ndof = 5 ! Number of degrees of freedom per node (3 displacements + damage + hydrogen concentration)
    
    integer, parameter :: start_u_idx = 1
    integer, parameter :: end_u_idx = 24 ! = nnode * ndim
    integer, parameter :: start_phi_idx = 25
    integer, parameter :: end_phi_idx = 32
    integer, parameter :: start_CL_idx = 33
    integer, parameter :: end_CL_idx = 40

    ! We scale up amatrx and rhs to avoid numerical issues of very small numbers
    ! In other words, we avoid warnings of 
    ! There is zero force everywhere 
    ! There is zero flux everywhere
    ! There is zero moment everywhere
    ! and so on
    ! Example: The solution is the same between these two equations
    ! 5 * x = 7
    ! 5e-23 * x = 7e-23
    ! The solution is x = 7/5 = 1.4, but the second equation will give x = 0 due to
    ! very small rhs that Abaqus treats rhs as zero
    real(kind=dp), parameter :: scale_up = 1.0d6

    
    ! The following data is not part of UEL, defined by the user    
    real(kind=dp), dimension(ndim * nnode) :: u_prev, du_prev
    real(kind=dp), dimension(nnode) :: phi_prev, dphi_prev
    real(kind=dp), dimension(nnode) :: CL_prev, dCL_prev
    
    !real(kind=dp), dimension(ninpt) :: weight                      ! Weights for integration points
    real(kind=dp), dimension(1,nnode) :: shape_node_to_int    ! Shape function that interpolates from nodal points to integration points
                                                                    ! The extra 1 dimension is for matrix multiplication, otherwise it would be a vector
    real(kind=dp), dimension(ninpt) :: shape_int_to_node         ! Shape function that extrapolates from integration points to nodal points
    real(kind=dp), dimension(ndim,nnode) :: B_deriv_local           ! Derivatives of N_shape_nodal_to_int with respect to isoparametric coordinates
    real(kind=dp), dimension(ndim,nnode) :: B_deriv_global          ! Derivatives of N_shape_nodal_to_int with respect to global coordinates
                                                                    ! This is the collection of vectors with spatial derivatives of N_shape_nodal_to_int
                                                                    ! Each column is the vector B_i in Emilio et al. 
    real(kind=dp), dimension(ntens,nnode * ndim) :: Bu_matrix         ! Strain-displacement matrix (B matrix)
    real(kind=dp), dimension(nnode, nnode) :: BB                    ! B matrix transpose times B matrix
    real(kind=dp), dimension(ntens,ntens) :: ddsdde                 ! Tangent stiffness matrix 
    real(kind=dp), dimension(ntens) :: stress                       ! Stress vector of the current element jelem
    real(kind=dp), dimension(ntens) :: stran                        ! Strain vector of the current element jelem
    real(kind=dp), dimension(ntens) :: dstran                       ! Incremental strain vector of the current element jelem
    real(kind=dp), dimension(nnode,nnode)  :: M_conc_capacity       ! Concentration capacity matrix (for hydrogen diffusion)
    real(kind=dp), dimension(nnode,nnode) :: K_diffusitivity        ! Diffusitivity matrix (for hydrogen diffusion)
    real(kind=dp), dimension(nnode,1) :: sig_H_node                 ! Hydrostatic stress at each node
    real(kind=dp), dimension(nnode,1) :: softened_sig_H_node       ! Softened hydrostatic stress at each node
    real(kind=dp), dimension(nnode,nnode) :: softened_sig_H_int         ! Softened hydrostatic stress at each integration point
    real(kind=dp), dimension(ndim,nnode) :: softened_grad_sig_H_int        ! Softened hydrostatic stress gradient at each node
    real(kind=dp), dimension(nsvint) :: statev_element              ! Local state variables of the current element jelem
    real(kind=dp), dimension(ndim, ndim) :: jac, inv_jac            ! Jacobian and its inverse

    ! Declaring all props as real(kind=dp) for consistency
    real(kind=dp) :: length_scale, Gc0, xkap, chi_DFT, delta_g_b0, R, T, VH, DL, KL
    real(kind=dp) :: CL_mol, CL_wtppm, CL_molfrac

!   initialising
    do k1=1,ndofel
        rhs(k1,1)=0.d0
    end do
    amatrx=0.d0
    
!   find total number of elements and stored it to nelem       
    if (dtime==0.d0) then
        if (jelem==1) then
            nelem=jelem
        else
            if (jelem>nelem) then
                nelem=jelem 
            end if
        endif 
    endif      

    ! Extract from the variable u and du
    u_prev(1:ndim * nnode)  = u(start_u_idx:end_u_idx)
    phi_prev(1:nnode)       = u(start_phi_idx:end_phi_idx)
    CL_prev(1:nnode)        = u(start_CL_idx:end_CL_idx)

    du_prev(1:ndim * nnode) = du(start_u_idx:end_u_idx, 1)
    dphi_prev(1:nnode)      = du(start_phi_idx:end_phi_idx, 1)
    dCL_prev(1:nnode)       = du(start_CL_idx:end_CL_idx, 1)

!   reading parameters
    length_scale=props(start_phase_props_idx)
    Gc0 = props(start_phase_props_idx+1) ! Critical energy release rate in the absence of hydrogen
    xkap = props(start_phase_props_idx+2) ! Well-conditioning parameter
    chi_DFT = props(start_phase_props_idx+3) ! Fitting slope to the DFT data
    delta_g_b0 = props(start_phase_props_idx+4) ! Gibbs free energy difference between the decohering interface and the surrounding material
    
    R = props(start_hydro_props_idx) ! Universal gas constant R (N*m)/(mol*K))
    T = props(start_hydro_props_idx+1) ! Temperature (K)
    VH = props(start_hydro_props_idx+2) ! Molar volume of H (m^3/mol)
    DL = props(start_hydro_props_idx+3) ! Diffusion coefficient of lattice species hydrogen

!   compute the hydrostatic stress
    sig_H_node = 0.d0    
    
    do inode = 1,nnode
        xi_node = xi_nodal_extra(inode)
        eta_node = eta_nodal_extra(inode)
        zeta_node = zeta_nodal_extra(inode)
        ! Change them to xi_node, eta_node, zeta_node
        shape_int_to_node(1) = 0.125d0 * (1.d0 - xi_node) * (1.d0 - eta_node) * (1.d0 - zeta_node)
        shape_int_to_node(2) = 0.125d0 * (1.d0 + xi_node) * (1.d0 - eta_node) * (1.d0 - zeta_node)
        shape_int_to_node(3) = 0.125d0 * (1.d0 - xi_node) * (1.d0 + eta_node) * (1.d0 - zeta_node)
        shape_int_to_node(4) = 0.125d0 * (1.d0 + xi_node) * (1.d0 + eta_node) * (1.d0 - zeta_node)
        shape_int_to_node(5) = 0.125d0 * (1.d0 - xi_node) * (1.d0 - eta_node) * (1.d0 + zeta_node)
        shape_int_to_node(6) = 0.125d0 * (1.d0 + xi_node) * (1.d0 - eta_node) * (1.d0 + zeta_node)
        shape_int_to_node(7) = 0.125d0 * (1.d0 - xi_node) * (1.d0 + eta_node) * (1.d0 + zeta_node)
        shape_int_to_node(8) = 0.125d0 * (1.d0 + xi_node) * (1.d0 + eta_node) * (1.d0 + zeta_node)
        
        ! print *, 'shape_int_to_node', shape_int_to_node
        do i = 1,ninpt
            isvinc = (i-1) * nsvint
            sig_H_node(inode,1) = sig_H_node(inode,1) &
                                + shape_int_to_node(i) * svars(isvinc + sig_H_idx)
        end do
    end do      
      
    do kintk=1,ninpt
        !   Compute shape_node_to_int and B_deriv_local
        call kshapefcn(kintk,ninpt,nnode,ndim,shape_node_to_int,B_deriv_local)      
        !   Compute djac and B_deriv_global
        call kjacobian(jelem,ndim,nnode,coords,B_deriv_local,djac,B_deriv_global,mcrd)
        !   Calculate strain displacement B-matrix
        call kbmatrix(B_deriv_global,ntens,nnode,ndim,Bu_matrix)     
        !   Differential volume at the integration point
        dvol = weight(kintk) * djac        
       
        ! print*, 'okay'
        !   compute from nodal values
        phi = 0.d0
        CL_mol = 0.d0
        do inode = 1,nnode
            phi = phi + shape_node_to_int(1,inode) * phi_prev(inode)
            CL_mol = CL_mol + shape_node_to_int(1,inode) * CL_prev(inode)
        end do   
        if (phi > 1.d0) then 
            phi = 1.d0
        else if (phi < 0.d0) then
            phi = 0.d0
        endif
       
        !   hydrogen contribution  
        KL = dexp(-delta_g_b0/(R * T))
        ! Convert CL from mol to wtppm
        CL_wtppm = (CL_mol * molar_mass_H) / (density_metal * 1.d-03)
        ! Conversion formula, if CL is in mol/m^3 and Cbar_L is in wtppm
        ! CL (wtppm) = [ Cbar_L (mol/m^3) * molar_mass_H (g/mol) ] / [ density_metal (kg/m^3) * 1e-06 (1g/1000kg) * 1000 (g/kg) ]

        ! Convert CL from wtppm to molfrac
        CL_molfrac = (CL_wtppm * 1.d-6) * (ratio_molar_mass_Fe_H)
        theta_coverage = CL_molfrac/(CL_molfrac + KL)
        Gc = Gc0 * (1.d0-chi_DFT * theta_coverage)
           
        !   compute the increment of strain and recover history variables
        dstran = matmul(Bu_matrix,du_prev)

        !   Transfer data from svars to statev_element for current element
        call kstatevar(kintk,nsvint,svars,statev_element,1)

        ! Extract stress, strain, phi, history from statev_element
        stress = statev_element(1:ntens)
        stran(1:ntens) = statev_element((ntens+1):(2 * ntens))
        phi_n = statev_element(phi_idx)
        history_n = statev_element(history_idx)
        
        if (dtime == 0.d0) then
            phi_n = phi
        endif
       
        !   compute strain energy density from the previous increment       
        psi = 0.d0
        do k1 = 1,ntens
            psi = psi + stress(k1) * stran(k1) * 0.5d0
        end do
       
        !   call umat to obtain stresses and constitutive matrix 
        ! call UMAT_elastic(props,nprops,ddsdde,stress,dstran,ntens,ndi,nshr,statev_element)
        
        call UMAT_von_Mises(props,nprops,ddsdde,stress,dstran,ntens,ndi,nshr,statev_element)
        
        stran = stran + dstran
       
        !   enforcing Karush-Kuhn-Tucker conditions
        if (psi > history_n) then
            history = psi
        else
            history = history_n
        endif
        
        !C_wtppm = (C_mol * molar_mass_H) / (density_metal * 1.0d-6 * 1000.d0)
        !CL_wtppm = (CL_mol * molar_mass_H) / (density_metal * 1.0d-6 * 1000.d0)
        !CT_wtppm = (CT_mol * molar_mass_H) / (density_metal * 1.0d-6 * 1000.d0)

        ! Update the state variables
        statev_element(1:ntens) = stress(1:ntens)
        statev_element((ntens+1):(2 * ntens)) = stran(1:ntens)
        statev_element(eqplas_idx) = 0.0
        statev_element(sig_H_idx) = (stress(1)+stress(2)+stress(3))/3.d0
        statev_element(sig_vonMises_idx) = 0.0
        statev_element(triax_idx) = 0.0
        statev_element(lode_idx) = 0.0
        statev_element(phi_idx) = phi
        statev_element(history_idx) = history
        statev_element(Gc_idx) = Gc
        statev_element(theta_coverage_idx) = theta_coverage
        statev_element(C_mol_idx) = 0.0
        statev_element(CL_mol_idx) = CL_mol
        statev_element(CT_mol_idx) = 0.0
        
        !   Transfer data from statev_element to svars
        !   This stage basically updates the state variables for the current elemennt in UEL
        call kstatevar(kintk,nsvint,svars,statev_element,0)
        
        softened_factor = (1.d0 - phi_n)**2 + xkap

        if (softened_factor > 1.d0) then
            softened_factor = 1.d0
        else if (softened_factor < 0.d0) then
            softened_factor = 0.d0
        endif
        ! ********************************************!
        ! DISPLACEMENT CONTRIBUTION TO amatrx AND rhs !
        ! ********************************************!

        ! 3D case
        ! 8 nodes x 3 displacement dofs ux, uy, uz = 24

        amatrx(start_u_idx:end_u_idx,start_u_idx:end_u_idx) = &
            amatrx(start_u_idx:end_u_idx,start_u_idx:end_u_idx) + dvol * (softened_factor &
                                * matmul(matmul(transpose(Bu_matrix),ddsdde),Bu_matrix))
            
        rhs(start_u_idx:end_u_idx,1) = rhs(start_u_idx:end_u_idx,1) - &
            dvol * (matmul(transpose(Bu_matrix),stress) * softened_factor)       

        ! *******************************************!
        ! PHASE FIELD CONTRIBUTION TO amatrx AND rhs !
        ! *******************************************!

        ! 3D case
        ! 8 nodes x 1 phase field dof = 8

        amatrx(start_phi_idx:end_phi_idx,start_phi_idx:end_phi_idx) = &
            amatrx(start_phi_idx:end_phi_idx,start_phi_idx:end_phi_idx) &
            + dvol * (matmul(transpose(B_deriv_global),B_deriv_global) * Gc * length_scale &
            + matmul(transpose(shape_node_to_int),shape_node_to_int) * (Gc/length_scale+2.d0 * history))   

        rhs(start_phi_idx:end_phi_idx,1) = rhs(start_phi_idx:end_phi_idx,1) &
            -dvol * (matmul(transpose(B_deriv_global),matmul(B_deriv_global,phi_prev)) &
             * Gc * length_scale+shape_node_to_int(1,:) &
             * ((Gc/length_scale+2.d0 * history) * phi-2.d0 * history))
            
        ! **************************************************!
        ! HYDROGEN DIFFUSION CONTRIBUTION TO amatrx AND rhs !
        ! **************************************************!

        ! 3D case
        ! 8 nodes x 1 hydrogen concentration dof = 8

        M_conc_capacity = matmul(transpose(shape_node_to_int),shape_node_to_int)/DL
        BB = matmul(transpose(B_deriv_global),B_deriv_global)

        K_diffusitivity = BB - VH/(R*T) * matmul(BB,matmul((sig_H_node * softened_factor),shape_node_to_int))

        amatrx(start_CL_idx:end_CL_idx,start_CL_idx:end_CL_idx) = &
            amatrx(start_CL_idx:end_CL_idx,start_CL_idx:end_CL_idx) &
            +dvol * (M_conc_capacity/dtime+K_diffusitivity)
            
        rhs(start_CL_idx:end_CL_idx,1) = rhs(start_CL_idx:end_CL_idx,1) - &
            dvol * (matmul(K_diffusitivity,CL_prev)+ &
            matmul(M_conc_capacity,dCL_prev)/dtime)

        !   Transfer data from statev_element to dummy mesh for visualization
        user_vars(jelem,1:ntens,kintk) = statev_element(1:6) * softened_factor
        user_vars(jelem,ntens+1:ntens * 2,kintk) = statev_element((ntens+1):(ntens * 2))
        user_vars(jelem,ntens*2:nsvint,kintk) = statev_element(ntens*2:nsvint)
        
    end do       ! end loop on material integration points
    
    ! **************************!
    ! Scaling up amatrx AND rhs !
    ! **************************!

    
    ! rhs(start_u_idx:end_u_idx,1) = rhs(start_u_idx:end_u_idx,1) * scale_up
    rhs(start_phi_idx:end_phi_idx,1) = rhs(start_phi_idx:end_phi_idx,1) * scale_up
    rhs(start_CL_idx:end_CL_idx,1) = rhs(start_CL_idx:end_CL_idx,1) * scale_up
    
    ! amatrx(start_u_idx:end_u_idx,start_u_idx:end_u_idx) = &
    !     amatrx(start_u_idx:end_u_idx,start_u_idx:end_u_idx) * scale_up
    amatrx(start_phi_idx:end_phi_idx,start_phi_idx:end_phi_idx) = &
        amatrx(start_phi_idx:end_phi_idx,start_phi_idx:end_phi_idx) * scale_up
    amatrx(start_CL_idx:end_CL_idx,start_CL_idx:end_CL_idx) = &
        amatrx(start_CL_idx:end_CL_idx,start_CL_idx:end_CL_idx) * scale_up

    ! print *, 'okay'
return
end
    

!***********************************************************************

subroutine UMAT(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt, &
    drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred, &
    cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt, &
    celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

    use common_block
    include 'aba_param.inc' 

    character*8 cmname
    dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens), &
        ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens), &
        time(2),predef(1),dpred(1),props(nprops),coords(3),drot(3,3), &
        dfgrd0(3,3),dfgrd1(3,3),jstep(4)

    ddsdde = 0.0d0
    noffset = noel - nelem    
    ! nelem: number of elements of UEL: [1, nelem]
    ! noel: number of elements of UMAT: [nelem + 1, 2 * nelem]
    ! => noffset: number of elements of UMAT offset by nelem: [nelem + 1, 2 * nelem] - nelem = [1, nelem]
    statev(1:nstatv) = user_vars(noffset,1:nstatv,npt)

return
end