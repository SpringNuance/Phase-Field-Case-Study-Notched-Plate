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

subroutine kshapefcn(kintk,ninpt,nnode,ndim,deriv_shape_node_to_int,B_deriv_local)
!   
    use iso_module
    include 'aba_param.inc'
!
    dimension deriv_shape_node_to_int(1, nnode),B_deriv_local(ndim,nnode)

!   determine (g,h,r)
    f = xi_int_inter(kintk)
    g = eta_int_inter(kintk)
    h = zeta_int_inter(kintk)

!   shape functions
    deriv_shape_node_to_int(1,1)=0.125d0*(1.d0-f)*(1.d0-g)*(1.d0-h)
    deriv_shape_node_to_int(1,2)=0.125d0*(1.d0+f)*(1.d0-g)*(1.d0-h)
    deriv_shape_node_to_int(1,3)=0.125d0*(1.d0+f)*(1.d0+g)*(1.d0-h)
    deriv_shape_node_to_int(1,4)=0.125d0*(1.d0-f)*(1.d0+g)*(1.d0-h)
    deriv_shape_node_to_int(1,5)=0.125d0*(1.d0-f)*(1.d0-g)*(1.d0+h)
    deriv_shape_node_to_int(1,6)=0.125d0*(1.d0+f)*(1.d0-g)*(1.d0+h)
    deriv_shape_node_to_int(1,7)=0.125d0*(1.d0+f)*(1.d0+g)*(1.d0+h)
    deriv_shape_node_to_int(1,8)=0.125d0*(1.d0-f)*(1.d0+g)*(1.d0+h)

!   derivative d(Ni)/d(f)
    B_deriv_local(1,1)=-0.125d0*(1.d0-g)*(1.d0-h)
    B_deriv_local(1,2)= 0.125d0*(1.d0-g)*(1.d0-h)
    B_deriv_local(1,3)= 0.125d0*(1.d0+g)*(1.d0-h)
    B_deriv_local(1,4)=-0.125d0*(1.d0+g)*(1.d0-h)
    B_deriv_local(1,5)=-0.125d0*(1.d0-g)*(1.d0+h)
    B_deriv_local(1,6)= 0.125d0*(1.d0-g)*(1.d0+h)
    B_deriv_local(1,7)= 0.125d0*(1.d0+g)*(1.d0+h)
    B_deriv_local(1,8)=-0.125d0*(1.d0+g)*(1.d0+h)

!     derivative d(Ni)/d(g)
    B_deriv_local(2,1)=-0.125d0*(1.d0-f)*(1.d0-h)
    B_deriv_local(2,2)=-0.125d0*(1.d0+f)*(1.d0-h)
    B_deriv_local(2,3)= 0.125d0*(1.d0+f)*(1.d0-h)
    B_deriv_local(2,4)= 0.125d0*(1.d0-f)*(1.d0-h)
    B_deriv_local(2,5)=-0.125d0*(1.d0-f)*(1.d0+h)
    B_deriv_local(2,6)=-0.125d0*(1.d0+f)*(1.d0+h)
    B_deriv_local(2,7)= 0.125d0*(1.d0+f)*(1.d0+h)
    B_deriv_local(2,8)= 0.125d0*(1.d0-f)*(1.d0+h)

!     derivative d(Ni)/d(h)
    B_deriv_local(3,1)=-0.125d0*(1.d0-f)*(1.d0-g)
    B_deriv_local(3,2)=-0.125d0*(1.d0+f)*(1.d0-g)
    B_deriv_local(3,3)=-0.125d0*(1.d0+f)*(1.d0+g)
    B_deriv_local(3,4)=-0.125d0*(1.d0-f)*(1.d0+g)
    B_deriv_local(3,5)= 0.125d0*(1.d0-f)*(1.d0-g)
    B_deriv_local(3,6)= 0.125d0*(1.d0+f)*(1.d0-g)
    B_deriv_local(3,7)= 0.125d0*(1.d0+f)*(1.d0+g)
    B_deriv_local(3,8)= 0.125d0*(1.d0-f)*(1.d0+g)

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
                    B_deriv_local(jdim,inode)*coords(idim,inode)
            end do
        end do
    end do

    if (ndim==3) then

        djac = xjac(1,1)*xjac(2,2)*xjac(3,3)+xjac(2,1)*xjac(3,2)*xjac(1,3) &
              +xjac(3,1)*xjac(2,3)*xjac(1,2)-xjac(3,1)*xjac(2,2)*xjac(1,3) &
              -xjac(2,1)*xjac(1,2)*xjac(3,3)-xjac(1,1)*xjac(2,3)*xjac(3,2)
        !if (djac>0.d0) then ! jacobian is positive - o.k.
        xjaci(1,1)=(xjac(2,2)*xjac(3,3)-xjac(2,3)*xjac(3,2))/djac
        xjaci(1,2)=(xjac(1,3)*xjac(3,2)-xjac(1,2)*xjac(3,3))/djac
        xjaci(1,3)=(xjac(1,2)*xjac(2,3)-xjac(1,3)*xjac(2,2))/djac
        xjaci(2,1)=(xjac(2,3)*xjac(3,1)-xjac(2,1)*xjac(3,3))/djac
        xjaci(2,2)=(xjac(1,1)*xjac(3,3)-xjac(1,3)*xjac(3,1))/djac
        xjaci(2,3)=(xjac(1,3)*xjac(2,1)-xjac(1,1)*xjac(2,3))/djac
        xjaci(3,1)=(xjac(2,1)*xjac(3,2)-xjac(2,2)*xjac(3,1))/djac
        xjaci(3,2)=(xjac(1,2)*xjac(3,1)-xjac(1,1)*xjac(3,2))/djac
        xjaci(3,3)=(xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1))/djac
        !else ! negative or zero jacobian
        ! write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
        !endif
    else if (ndim==2) then
        djac=xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1)
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

    dimension B_deriv_global(ndim,nnode),Bu_matrix(ntens,nnode*ndim)

    Bu_matrix=0.d0
    do inode=1,nnode
        Bu_matrix(1,ndim*inode-ndim+1)=B_deriv_global(1,inode)
        Bu_matrix(2,ndim*inode-ndim+2)=B_deriv_global(2,inode)
        Bu_matrix(4,ndim*inode-ndim+1)=B_deriv_global(2,inode)
        Bu_matrix(4,ndim*inode-ndim+2)=B_deriv_global(1,inode)
        if (ndim==3) then
            Bu_matrix(3,ndim*inode)=B_deriv_global(3,inode)
            Bu_matrix(5,ndim*inode-2)=B_deriv_global(3,inode)
            Bu_matrix(5,ndim*inode)=B_deriv_global(1,inode)
            Bu_matrix(6,ndim*inode-1)=B_deriv_global(3,inode)
            Bu_matrix(6,ndim*inode)=B_deriv_global(2,inode)
        endif
    end do

return
end

subroutine kstatevar(npt,nsvint,statev,statev_ip,icopy)

!   Transfer data to/from element-level state variable array from/to
!   material-point level state variable array.

    include 'aba_param.inc'

    dimension statev(*),statev_ip(*)

    isvinc=(npt-1)*nsvint     ! integration point increment

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

subroutine UMAT_elastic(props,ddsdde,stress,dstran,ntens,ndi,nshr,statev)
    
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
    lambda = E*nu/((1.0d0 + nu) * (1.0d0 - 2.0d0 * nu)) ! Lame's first constant

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

    parameter(ndim=3,ntens=6,ndi=3,nshr=3,ninpt=8,nsvint=15)
      
    integer, parameter :: start_u_idx = 1
    integer, parameter :: end_u_idx = 24 ! = nnode*ndim
    integer, parameter :: start_phi_idx = 25
    integer, parameter :: end_phi_idx = 32
    integer, parameter :: start_CL_idx = 33
    integer, parameter :: end_CL_idx = 40
    
    ! The following data is not part of UEL, defined by the user    
    real(kind=dp), dimension(ndim*nnode) :: u_prev, du_prev
    real(kind=dp), dimension(nnode) :: phi_prev, dphi_prev
    real(kind=dp), dimension(nnode) :: CL_prev, dCL_prev
    
    !real(kind=dp), dimension(ninpt) :: weight                      ! Weights for integration points
    real(kind=dp), dimension(1,nnode) :: deriv_shape_node_to_int    ! Shape function that interpolates from nodal points to integration points
                                                                    ! The extra 1 dimension is for matrix multiplication, otherwise it would be a vector
    real(kind=dp), dimension(ninpt) :: deriv_shape_int_to_node         ! Shape function that extrapolates from integration points to nodal points
    real(kind=dp), dimension(ndim,nnode) :: B_deriv_local           ! Derivatives of N_shape_nodal_to_int with respect to isoparametric coordinates
    real(kind=dp), dimension(ndim,nnode) :: B_deriv_global          ! Derivatives of N_shape_nodal_to_int with respect to global coordinates
                                                                    ! This is the collection of vectors with spatial derivatives of N_shape_nodal_to_int
                                                                    ! Each column is the vector B_i in Emilio et al. 
    real(kind=dp), dimension(ntens,nnode*ndim) :: Bu_matrix         ! Strain-displacement matrix (B matrix)
    real(kind=dp), dimension(ndim,nnode) :: bC
    real(kind=dp), dimension(nnode,nnode) :: BB
    real(kind=dp), dimension(ntens,ntens) :: ddsdde                 ! Tangent stiffness matrix 
    real(kind=dp), dimension(ntens) :: stress                       ! Stress vector of the current element jelem
    real(kind=dp), dimension(ntens) :: stran                        ! Strain vector of the current element jelem
    real(kind=dp), dimension(ntens) :: dstran                       ! Incremental strain vector of the current element jelem
    real(kind=dp), dimension(nnode,nnode)  :: M_conc_capacity       ! Concentration capacity matrix (for hydrogen diffusion)
    real(kind=dp), dimension(nnode,nnode) :: K_diffusitivity        ! Diffusitivity matrix (for hydrogen diffusion)
    real(kind=dp), dimension(nnode,1) :: sig_H                      ! Hydrostatic stress at each node
    real(kind=dp), dimension(nnode,1) :: softened_sig_H_nodal       ! Softened hydrostatic stress at each node
    real(kind=dp), dimension(nnode,1) :: softened_sig_H_int         ! Softened hydrostatic stress at each integration point
    real(kind=dp), dimension(ndim, 1) :: softened_grad_sig_H        ! Softened hydrostatic stress gradient at each node
    real(kind=dp), dimension(nsvint) :: statev_element              ! Local state variables of the current element jelem
    real(kind=dp), dimension(ndim, ndim) :: jac, inv_jac            ! Jacobian and its inverse

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
    u_prev(1:ndim*nnode)  = u(start_u_idx:end_u_idx)
    phi_prev(1:nnode)     = u(start_phi_idx:end_phi_idx)
    CL_prev(1:nnode)      = u(start_CL_idx:end_CL_idx)

    du_prev(1:ndim*nnode) = du(start_u_idx:end_u_idx, 1)
    dphi_prev(1:nnode)    = du(start_phi_idx:end_phi_idx, 1)
    dCL_prev(1:nnode)     = du(start_CL_idx:end_CL_idx, 1)

!   reading parameters
    xlc=props(9)
    Gc0 = props(10) ! Critical energy release rate in the absence of hydrogen
    xkap = props(11) ! Well-conditioning parameter
    R = props(17) ! Universal gas constant R (N*m)/(mol*K))
    T = props(18) ! Temperature (K)
    VH = props(19) ! Molar volume of H (m^3/mol)
    DL = props(20) ! Diffusion coefficient of lattice species hydrogen

!   compute the hydrostatic stress
    sig_H=0.d0    
    
    do inode=1,nnode
        xi_node = xi_nodal_extra(inode)
        eta_node = eta_nodal_extra(inode)
        zeta_node = zeta_nodal_extra(inode)
        ! Change them to xi_node, eta_node, zeta_node
        deriv_shape_int_to_node(1) = 0.125d0 * (1.d0 - xi_node) * (1.d0 - eta_node) * (1.d0 - zeta_node)
        deriv_shape_int_to_node(2) = 0.125d0 * (1.d0 + xi_node) * (1.d0 - eta_node) * (1.d0 - zeta_node)
        deriv_shape_int_to_node(3) = 0.125d0 * (1.d0 - xi_node) * (1.d0 + eta_node) * (1.d0 - zeta_node)
        deriv_shape_int_to_node(4) = 0.125d0 * (1.d0 + xi_node) * (1.d0 + eta_node) * (1.d0 - zeta_node)
        deriv_shape_int_to_node(5) = 0.125d0 * (1.d0 - xi_node) * (1.d0 - eta_node) * (1.d0 + zeta_node)
        deriv_shape_int_to_node(6) = 0.125d0 * (1.d0 + xi_node) * (1.d0 - eta_node) * (1.d0 + zeta_node)
        deriv_shape_int_to_node(7) = 0.125d0 * (1.d0 - xi_node) * (1.d0 + eta_node) * (1.d0 + zeta_node)
        deriv_shape_int_to_node(8) = 0.125d0 * (1.d0 + xi_node) * (1.d0 + eta_node) * (1.d0 + zeta_node)
        
        ! print *, 'deriv_shape_int_to_node', deriv_shape_int_to_node
        do i=1,ninpt
            isvinc=(i-1)*nsvint
            sig_H(inode,1)=sig_H(inode,1)+deriv_shape_int_to_node(i)*svars(isvinc+14)
        end do
    end do      
      
    do kintk=1,ninpt
!   evaluate shape functions and derivatives
        call kshapefcn(kintk,ninpt,nnode,ndim,deriv_shape_node_to_int,B_deriv_local)      
        ! kshapefcn(kintk,ninpt,nnode,ndim,deriv_shape_node_to_int,B_deriv_local)
        call kjacobian(jelem,ndim,nnode,coords,B_deriv_local,djac,B_deriv_global,mcrd)
        ! subroutine kjacobian(jelem,ndim,nnode,coords,B_deriv_local,djac,B_deriv_global,mcrd)
    
!   form B-matrix
        call kbmatrix(B_deriv_global,ntens,nnode,ndim,Bu_matrix)     
        dvol=weight(kintk)*djac        
       
        ! print*, 'okay'
!   compute from nodal values
        phi=0.d0
        CL_wtppm=0.d0
        do inode=1,nnode
            phi=phi+deriv_shape_node_to_int(1,inode)*phi_prev(inode)
            CL_wtppm=CL_wtppm+deriv_shape_node_to_int(1,inode)*CL_prev(inode)
        end do   
        if (phi>1.d0) phi=1.d0
       
!   hydrogen contribution  
        thetaL = CL_wtppm*5.5d-05/(CL_wtppm*5.5d-05+dexp(-3.d7/(R*T)))
        Gc = Gc0*(1.d0-0.89d0*thetaL)
           
!   compute the increment of strain and recover history variables
        dstran=matmul(Bu_matrix,du_prev)
        call kstatevar(kintk,nsvint,svars,statev_element,1)
        stress=statev_element(1:ntens)
        stran(1:ntens)=statev_element((ntens+1):(2*ntens))
        history_n=statev_element(2*ntens+3)
        phi_n=statev_element(2*ntens+1)
        if (dtime==0.d0) phi_n=phi
       
!   compute strain energy density from the previous increment       
        psi=0.d0
        do k1=1,ntens
            psi=psi+stress(k1)*stran(k1)*0.5d0
        end do
       
!   call umat to obtain stresses and constitutive matrix 
        call UMAT_elastic(props,ddsdde,stress,dstran,ntens,ndi,nshr,statev_element)
        ! subroutine kumat(props,ddsdde,stress,dstran,ntens,statev)
        stran=stran+dstran
       
!   enforcing Karush-Kuhn-Tucker conditions
        if (psi>history_n) then
            history=psi
        else
            history=history_n
        endif
       
        statev_element(1:ntens)=stress(1:ntens)
        statev_element((ntens+1):(2*ntens))=stran(1:ntens)
        statev_element(2*ntens+1)=phi
        statev_element(2*ntens+2)=(stress(1)+stress(2)+stress(3))/3.d0 !SH
        statev_element(2*ntens+3)=history
        
        call kstatevar(kintk,nsvint,svars,statev_element,0)
        
        softened_factor = (1.d0 - phi_n)**2 + xkap
        ! ********************************************!
        ! DISPLACEMENT CONTRIBUTION TO amatrx AND rhs !
        ! ********************************************!

        ! 3D case
        ! 8 nodes x 3 displacement dofs ux, uy, uz = 24

        amatrx(start_u_idx:end_u_idx,start_u_idx:end_u_idx) = &
            amatrx(start_u_idx:end_u_idx,start_u_idx:end_u_idx)+dvol*(softened_factor &
                                * matmul(matmul(transpose(Bu_matrix),ddsdde),Bu_matrix))
            
        rhs(start_u_idx:end_u_idx,1)=rhs(start_u_idx:end_u_idx,1) - &
            dvol*(matmul(transpose(Bu_matrix),stress) * softened_factor)       

        ! *******************************************!
        ! PHASE FIELD CONTRIBUTION TO amatrx AND rhs !
        ! *******************************************!

        ! 3D case
        ! 8 nodes x 1 phase field dof = 8

        amatrx(start_phi_idx:end_phi_idx,start_phi_idx:end_phi_idx) = &
            amatrx(start_phi_idx:end_phi_idx,start_phi_idx:end_phi_idx) &
            + dvol*(matmul(transpose(B_deriv_global),B_deriv_global)*Gc*xlc &
            + matmul(transpose(deriv_shape_node_to_int),deriv_shape_node_to_int)*(Gc/xlc+2.d0*history))   

        rhs(start_phi_idx:end_phi_idx,1) = rhs(start_phi_idx:end_phi_idx,1) &
            -dvol*(matmul(transpose(B_deriv_global),matmul(B_deriv_global,phi_prev)) &
            *Gc*xlc+deriv_shape_node_to_int(1,:)*((Gc/xlc+2.d0*history)*phi-2.d0*history))
            
        ! **************************************************!
        ! HYDROGEN DIFFUSION CONTRIBUTION TO amatrx AND rhs !
        ! **************************************************!

        ! 3D case
        ! 8 nodes x 1 hydrogen concentration dof = 8

        M_conc_capacity = matmul(transpose(deriv_shape_node_to_int),deriv_shape_node_to_int)/DL
        BB = matmul(transpose(B_deriv_global),B_deriv_global)
        ! softened_sig_H_nodal = sig_H * softened_factor
        K_diffusitivity = matmul(transpose(B_deriv_global),B_deriv_global) &
            - VH/(R*T)*matmul(BB,matmul((sig_H * softened_factor),deriv_shape_node_to_int))

        ! M_conc_capacity = matmul(transpose(N_shape_nodal_to_int),N_shape_nodal_to_int)/DL
        ! ! B_deriv_global shape is (ndim, nnode)
        ! ! sig_H has shape (nnode, 1)
        ! softened_sig_H = sig_H * softened_factor ! shape (nnode, 1)
        ! softened_grad_sig_H = matmul(B_deriv_global, softened_sig_H) 
        ! ! shape (ndim, nnode) x (nnode, 1) = (ndim, 1)
        ! ! Equation 45
        ! K_diffusitivity = matmul(transpose(B_deriv_global),B_deriv_global) &
        !                 - VH/(R*T) * matmul(B_deriv_global, matmul(softened_grad_sig_H, N_shape_nodal_to_int))
        
        amatrx(start_CL_idx:end_CL_idx,start_CL_idx:end_CL_idx)= &
            amatrx(start_CL_idx:end_CL_idx,start_CL_idx:end_CL_idx) &
            +dvol*(M_conc_capacity/dtime+K_diffusitivity)
            
        rhs(start_CL_idx:end_CL_idx,1)=rhs(start_CL_idx:end_CL_idx,1) - &
            dvol*(matmul(K_diffusitivity,CL_prev)+ &
            matmul(M_conc_capacity,dCL_prev)/dtime)

! output
        user_vars(jelem,1:ntens,kintk) = statev_element(1:6) * softened_factor
        user_vars(jelem,ntens+1:ntens*2,kintk) = statev_element((ntens+1):(ntens*2))
        user_vars(jelem,(2*ntens+1),kintk) = statev_element(2*ntens+1)
        user_vars(jelem,(2*ntens+2),kintk) = statev_element(2*ntens+2)
        user_vars(jelem,(2*ntens+3),kintk) = CL_wtppm
        
    end do       ! end loop on material integration points
    
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