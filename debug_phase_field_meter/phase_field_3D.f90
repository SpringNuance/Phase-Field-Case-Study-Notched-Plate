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
!     statev(1:6): stress tensor 11, 22, 33, 12, 13, 23        
!     statev(7:12) : strain tensor 11, 22, 33, 12, 13, 23  
!     statev(13) : equivalent plastic strain PEEQ (eqplas)
!     statev(14) : hydrostatic stress (σm)
!     statev(15) : crack phase field (ϕ)
!     statev(16) : history variable field (H)
!     statev(17) : dislocation density rho_d (ρd)
!     statev(18) : total hydrogen concentration (C)
!     statev(19) : hydrogen concentration in lattice sites (CL)
!     statev(20) : hydrogen concentration in trap sites (CT)

!     Debug version
!     statev(1:6): stress tensor 11, 22, 33, 12, 13, 23        
!     statev(7:12) : strain tensor 11, 22, 33, 12, 13, 23  
!     statev(13) : crack phase field (ϕ)
!     statev(14) : hydrostatic stress (σm)
!     statev(15) : hydrogen concentration in lattice sites (CL)

!***********************************************************************

module precision
    use iso_fortran_env
    integer, parameter :: dp = real64
end module precision

!***********************************************************************

module common_block
    use precision
    implicit none
    ! First dim: maximum number of elements to accomodate varying number of elements when remeshed
    ! Second dim: number of solution state dependent variables (nsvars in UEL and nstatev in UMAT)
    ! Third dim: number of integration points

    integer, parameter :: start_mech_index = 1 ! Index of the first mechanical property in props
    integer, parameter :: start_phase_index = 9 ! Index of the first phase field property in props
    integer, parameter :: start_hydro_index = 17 ! Index of the first hydrogen diffusion property in props
    integer, parameter :: start_UHARD_index = 41 ! Index of the first hardening curve in props
    real(kind=dp) :: user_vars(100000,15,8)
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
    
    nvalue = (nprops - start_UHARD_index + 1) / 2

    ! print *, 'nvalue = ', nvalue ! 100
    ! print *, 'start_UHARD_index = ', start_UHARD_index ! 41
    
    call UHARD_von_Mises(syiel0, hard, eqplas, &
                                statev, nvalue, props(start_UHARD_index))
    
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
                                statev, nvalue, props(start_UHARD_index))
                                
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

    ! Generally, any variable that contains the letter "j" means its datatype is integer
    ! For example, "jelem" is the element number and "jtype" is the element type2, all in integer

    ! rhs: An array containing the contributions of this element to the 
    !      right-hand-side vectors of the overall system of equations

    ! amatrx: An array containing the contribution of this element to the 
    !         Jacobian (stiffness) or other matrix of the overall system of equations

    ! svars: An array containing the values of the solution-dependent state variables 
    !        associated with this element. The number of such variables is nsvars

    ! props: A floating point array containing the nprops real property values 
    !        defined for use with this element. nprops is the user-specified number 
    !        of real property values. We usually use this only

    ! jprops: An integer array containing the njprop integer property values 
    !         defined for use with this element. njprop is the user-specified number 
    !         of integer property values. We don't use this

    ! coords: An array containing the original coordinates of the nodes of the element. 
    !         coords(k1,k2) is the k1-th coordinate of the k2-th node of the element.

    ! u, du, v, a: Arrays containing the current estimates of the basic solution 
    !             variables (displacements, rotations, temperatures, 
    !              depending on the degree of freedom) 
    ! u(k1)	Total values of the variables
    ! du(k1, 1)	Incremental values of the variables
    ! v(k1)		Time rate of change of the variables (velocities, rate of rotations, etc)
    ! a(k1)		Acceleartion of the variables (accelerations, angular accelerations, etc)
    
    ! time(1): Current value of step time or frequency.
    ! time(2): Current value of total time.
    ! dtime: time increment
    ! ndofel: Number of degrees of freedom in the element.
    !         this is equal to nnode x number of degree of freedom per node
    !         For user element that resembles C3D8, ndofel = 8 nodes x 5 dof = 40
    !         where 5 dof are ux, uy, yz (displacement), phi (crack phase field), CL (lattice hydrogen concentration)
    
    ! mlvarx: Dimensioning parameter used when several displacement or right-hand-side vectors are used.
    ! nrhs: Number of load vectors. NRHS is 1 in most nonlinear problems: it is 2 for the modified Riks static procedure7
    ! mcrd: It is the maximum of the user-defined maximum number of coordinates required at any node point 
    !       and the value of the largest active degree of freedom of the user element 
    !       that is less than or equal to 3. For example, if you specify that the 
    !       maximum number of coordinates is 1 and the active degrees of freedom of 
    !       the user element are 2, 3, and 6, MCRD will be 3. If you specify 
    !       that the maximum number of coordinates is 2 and the active degrees of 
    !       freedom of the user element are 11 and 12, MCRD will be 2

    ! nnode: User-defined number of nodes on the element (defined in User Element section in input file in node flag)
    ! jtype: Integer defining the element type. This is the user-defined integer value n in element type Un 
    !        (defined in User Element section in input file in type=Un flag)
    
    ! jelem: Current user element number.

    !     State variables  
    !     nsvars(1:nsvint) - statev for first integration point
    !     nsvars(nsvint+1:2*nsvint) - statev for second integration point
    !     ...
    !     nsvars((ninpt-1)*nsvint+1:ninpt*nsvint) - statev for last integration point
    !     where 0, nsvint, 2*nsvint, ..., (ninpt-1)*nsvint are the starting indices of the state variables for each integration point
    !     defined as nsvinc=(i-1)*nsvint where i is the integration point number (1 to ninpt)

    !     Now for each integration point, the state variables are defined as follows: 
    !     nsvars(1:6): stress tensor 11, 22, 33, 12, 13, 23  for first integration point     
    !     nsvars(7:12) : strain tensor 11, 22, 33, 12, 13, 23  for first integration point
    !     nsvars(13) : equivalent plastic strain PEEQ (eqplas) for first integration point
    !     nsvars(14) : delta equivalent plastic strain PEEQ (deqplas) for first integration point
    !     nsvars(15) : Crack Phase Field (ϕ) for first integration point
    !     nsvars(16) : hydrostatic stress (σm) for first integration point
    !     nsvars(17) : History variable field (H) for first integration point
    !     nsvars(18) : Total hydrogen concentration (C) for first integration point
    !     nsvars(19) : Hydrogen concentration in lattice sites (CL) for first integration point
    !     nsvars(20) : Hydrogen concentration in trap sites (CT) for first integration point

    !     nsvars(21:26): stress tensor 11, 22, 33, 12, 13, 23  for second integration point
    !     ...
    !     nsvars(39): Hydrogen concentration in trap sites (CT) for last integration point
    !     and so on until nsvars(nsvint*ninpt) for the last integration point

    parameter(ndim=3, ntens=6, ndi=3, nshr=3, ninpt=8, nsvint=15, ndof=5)
      
    real(kind=dp), parameter :: molar_mass_H = 1.00784d0 ! g/mol
    real(kind=dp), parameter :: molar_mass_Fe = 55.845d0 ! g/mol
    real(kind=dp), parameter :: ratio_molar_mass_Fe_H = 55.415d0

    integer, parameter :: start_u_idx = 1
    integer, parameter :: end_u_idx = 24 ! = nnode * ndim
    integer, parameter :: start_phi_idx = 25
    integer, parameter :: end_phi_idx = 32
    integer, parameter :: start_CL_idx = 33
    integer, parameter :: end_CL_idx = 40
    
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
    real(kind=dp), dimension(nnode,nnode) :: BB
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
    real(kind=dp) :: CL_wtppm, CL_molfrac

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
    length_scale=props(9)
    Gc0 = props(10) ! Critical energy release rate in the absence of hydrogen
    xkap = props(11) ! Well-conditioning parameter
    chi_DFT = props(12) ! Fitting slope to the DFT data
    delta_g_b0 = props(13) ! Gibbs free energy difference between the decohering interface and the surrounding material
    R = props(17) ! Universal gas constant R (N*m)/(mol*K))
    T = props(18) ! Temperature (K)
    VH = props(19) ! Molar volume of H (m^3/mol)
    DL = props(20) ! Diffusion coefficient of lattice species hydrogen

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
            sig_H_node(inode,1) = sig_H_node(inode,1) + shape_int_to_node(i) * svars(isvinc+14)
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
        CL_wtppm = 0.d0
        do inode = 1,nnode
            phi = phi+shape_node_to_int(1,inode) * phi_prev(inode)
            CL_wtppm = CL_wtppm+shape_node_to_int(1,inode) * CL_prev(inode)
        end do   
        if (phi>1.d0) then 
            phi = 1.d0
        else if (phi<0.d0) then
            phi = 0.d0
        endif
       
        !   hydrogen contribution  
        KL = dexp(-delta_g_b0/(R * T))
        CL_molfrac = (CL_wtppm * 1.d-6) * (ratio_molar_mass_Fe_H)
        thetaL = CL_molfrac/(CL_molfrac + KL)
        Gc = Gc0 * (1.d0-chi_DFT * thetaL)
           
        !   compute the increment of strain and recover history variables
        dstran = matmul(Bu_matrix,du_prev)
        call kstatevar(kintk,nsvint,svars,statev_element,1)
        stress = statev_element(1:ntens)
        stran(1:ntens) = statev_element((ntens+1):(2 * ntens))
        history_n = statev_element(2 * ntens+3)
        phi_n = statev_element(2 * ntens+1)
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
       
        statev_element(1:ntens) = stress(1:ntens)
        statev_element((ntens+1):(2 * ntens)) = stran(1:ntens)
        statev_element(2 * ntens+1) = phi
        statev_element(2 * ntens+2) = (stress(1)+stress(2)+stress(3))/3.d0 !SH
        statev_element(2 * ntens+3) = history
        
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
        softened_sig_H_node = sig_H_node * softened_factor ! shape (nnode, 1)
        softened_sig_H_int = matmul(softened_sig_H_node, shape_node_to_int)
        ! shape (nnode, 1) x (1, nnode) = (nnode, nnode)
        softened_grad_sig_H_int = matmul(B_deriv_global, softened_sig_H_int) 
        ! shape (ndim, nnode) x (nnode, nnode) = (ndim, nnode)
        ! softened_sig_H_nodal = sig_H * softened_factor
        K_diffusitivity = matmul(transpose(B_deriv_global),B_deriv_global) &
            - VH/(R*T) * matmul(BB,matmul((sig_H_node * softened_factor),shape_node_to_int))

        amatrx(start_CL_idx:end_CL_idx,start_CL_idx:end_CL_idx) = &
            amatrx(start_CL_idx:end_CL_idx,start_CL_idx:end_CL_idx) &
            +dvol * (M_conc_capacity/dtime+K_diffusitivity)
            
        rhs(start_CL_idx:end_CL_idx,1) = rhs(start_CL_idx:end_CL_idx,1) - &
            dvol * (matmul(K_diffusitivity,CL_prev)+ &
            matmul(M_conc_capacity,dCL_prev)/dtime)

! output
        user_vars(jelem,1:ntens,kintk) = statev_element(1:6) * softened_factor
        user_vars(jelem,ntens+1:ntens * 2,kintk) = statev_element((ntens+1):(ntens * 2))
        user_vars(jelem,(2 * ntens+1),kintk) = statev_element(2 * ntens+1)
        user_vars(jelem,(2 * ntens+2),kintk) = statev_element(2 * ntens+2)
        user_vars(jelem,(2 * ntens+3),kintk) = CL_wtppm
        
    end do       ! end loop on material integration points
    
    ! print *, 'okay'
return
end
    
!***********************************************************************

subroutine UVARM(uvar,direct,t,time,dtime,cmname,orname, &
    nuvarm,noel,npt,layer,kspt,kstep,kinc,ndi,nshr,coord, &
    jmac,jmatyp,matlayo,laccfla)
    
    use precision
    include 'aba_param.inc'
!
    character*80 cmname,orname
    character*3 flgray(15)
    dimension uvar(nuvarm),direct(3,3),t(3,3),time(2)
    dimension array(15),jarray(15),jmac(*),jmatyp(*),coord(*)

!     the dimensions of the variables flgray, array and jarray
!     must be set equal to or greater than 15.

    ! Variables to Be Defined
    ! uvar(nuvarm)
    ! An array containing the user-defined output variables. 
    ! These are passed in as the values at the beginning of the increment 
    ! and must be returned as the values at the end of the increment.
    
    ! volume of the Thermal Desorption Spectroscopy (TDS) sample is 
    ! real(kind=dp), parameter :: volume = 6000.0d-11 ! 0.015 m * 0.004 m * 0.001 m = 6000e-11 m^3
    real(kind=dp), parameter :: density_metal = 7900.0d0 ! kg/m^3
    real(kind=dp), parameter :: molar_mass_H = 1.00784 ! g/mol
    real(kind=dp), parameter :: molar_mass_Fe = 55.845 ! g/mol
    real(kind=dp), parameter :: avogadro = 6.022e23 ! 1/mol
    
    ! mass_metal = volume * density_metal ! kg

    ! Conversion formula, if CL is in wtppm and Cbar_L is in mol/m^3
    ! Cbar_L (mol/m^3) = [ CL (wtppm) * 1e-06 (1g/1000kg) * density_metal (kg/m^3) * 1000 (g/kg) ] / molar_mass_H (g/mol)
    ! We use this case

    ! Conversion formula, if CL is in mol/m^3 and Cbar_L is in wtppm
    ! CL (wtppm) = [ Cbar_L (mol/m^3) * molar_mass_H (g/mol) ] / [ density_metal (kg/m^3) * 1e-06 (1g/1000kg) * 1000 (g/kg) ]

    ! Convert H+ concentration to wt.ppm
    
    call GETVRM('SDV',array,jarray,flgray,jcrd,jmac,jmatyp,matlayo,laccfla)
    C_mol = array(3) ! Unit is mol/m^3
    CL_mol = array(4) ! Unit is mol/m^3
    CT_mol = array(5) ! Unit is mol/m^3
    CT_dis_mol = array(6) ! Unit is mol/m^3
    CT_gb_mol = array(7) ! Unit is mol/m^3
    CT_carb_mol = array(8) ! Unit is mol/m^3

    ! Convert H+ concentration to wt.ppm
    C_wtppm = (C_mol * molar_mass_H) / (density_metal * 1.0d-6 * 1000.d0)
    CL_wtppm = (CL_mol * molar_mass_H) / (density_metal * 1.0d-6 * 1000.d0)
    CT_wtppm = (CT_mol * molar_mass_H) / (density_metal * 1.0d-6 * 1000.d0)
    CT_dis_wtppm = (CT_dis_mol * molar_mass_H) / (density_metal * 1.0d-6 * 1000.d0)
    CT_gb_wtppm = (CT_gb_mol * molar_mass_H) / (density_metal * 1.0d-6 * 1000.d0)
    CT_carb_wtppm = (CT_carb_mol * molar_mass_H) / (density_metal * 1.0d-6 * 1000.d0)

    ! Convert H+ concentration to mole fraction
    C_mole_fraction = (C_wtppm * 1.0d-6) * (molar_mass_Fe/ molar_mass_H)
    CL_mole_fraction = (CL_wtppm * 1.0d-6) * (molar_mass_Fe/ molar_mass_H)
    CT_mole_fraction = (CT_wtppm * 1.0d-6) * (molar_mass_Fe/ molar_mass_H)
    CT_dis_mole_fraction = (CT_dis_wtppm * 1.0d-6) * (molar_mass_Fe/ molar_mass_H)
    CT_gb_mole_fraction = (CT_gb_wtppm * 1.0d-6) * (molar_mass_Fe/ molar_mass_H)
    CT_carb_mole_fraction = (CT_carb_wtppm * 1.0d-6) * (molar_mass_Fe/ molar_mass_H)

    ! Convert H+ concentration to 1/m^3
    C = C_mol * avogadro ! (1/m^3)
    CL = CL_mol * avogadro ! (1/m^3)
    CT = CT_mol * avogadro ! (1/m^3)
    CT_dis = CT_dis_mol * avogadro ! (1/m^3)
    CT_gb = CT_gb_mol * avogadro ! (1/m^3)
    CT_carb = CT_carb_mol * avogadro ! (1/m^3)

    uvar(1) = C_wtppm
    uvar(2) = CL_wtppm
    uvar(3) = CT_wtppm
    uvar(4) = CT_dis_wtppm
    uvar(5) = CT_gb_wtppm
    uvar(6) = CT_carb_wtppm

    uvar(7) = C_molfrac
    uvar(8) = CL_molfrac
    uvar(9) = CT_molfrac
    uvar(10) = CT_dis_molfrac
    uvar(11) = CT_gb_molfrac
    uvar(12) = CT_carb_molfrac

    uvar(13) = C
    uvar(14) = CL
    uvar(15) = CT
    uvar(16) = CT_dis
    uvar(17) = CT_gb
    uvar(18) = CT_carb

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