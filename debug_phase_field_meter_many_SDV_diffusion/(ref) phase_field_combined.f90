
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
!         Inner number marked with x is intergration (gauss) points
!
!*****************************************************************

module iso_module

    use precision
    real(kind=dp), parameter :: coord_inter = 1.0d0
    real(kind=dp), parameter :: gauss_inter = 1.0d0 / sqrt(3.0d0)
    real(kind=dp), parameter :: coord_extra = sqrt(3.0d0)
    real(kind=dp), parameter :: gauss_extra = 1.0d0

    ! weight is the integration point weight for their shape function contribution
    real(kind=dp), parameter :: weight(8) = (/1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0/)
    
    ! Interpolating coordinates (nodal to gauss)
    ! Isoparametric coordinates for nodal points in hexahedral 3D element
    real(kind=dp), parameter :: xi_nodal_inter(8)   = (/ -coord_inter,  coord_inter,  coord_inter, -coord_inter, &
                                                         -coord_inter,  coord_inter,  coord_inter, -coord_inter /)
    real(kind=dp), parameter :: eta_nodal_inter(8)  = (/ -coord_inter, -coord_inter,  coord_inter,  coord_inter, &
                                                         -coord_inter, -coord_inter,  coord_inter,  coord_inter /)
    real(kind=dp), parameter :: zeta_nodal_inter(8) = (/ -coord_inter, -coord_inter, -coord_inter, -coord_inter, &
                                                          coord_inter,  coord_inter,  coord_inter,  coord_inter /)

    ! Isoparametric coordinates for integration points in hexahedral 3D element
    real(kind=dp), parameter :: xi_gauss_inter(8)   = (/ -gauss_inter,  gauss_inter, -gauss_inter,  gauss_inter, &
                                                         -gauss_inter,  gauss_inter, -gauss_inter,  gauss_inter /)
    real(kind=dp), parameter :: eta_gauss_inter(8)  = (/ -gauss_inter, -gauss_inter,  gauss_inter,  gauss_inter, &
                                                         -gauss_inter, -gauss_inter,  gauss_inter,  gauss_inter /)
    real(kind=dp), parameter :: zeta_gauss_inter(8) = (/ -gauss_inter, -gauss_inter, -gauss_inter, -gauss_inter, &
                                                          gauss_inter,  gauss_inter,  gauss_inter,  gauss_inter /)


    ! Extrapolating coordinates (gauss to nodal)
    real(kind=dp), parameter :: xi_nodal_extra(8)   = (/ -coord_extra,  coord_extra,  coord_extra, -coord_extra, &
                                                         -coord_extra,  coord_extra,  coord_extra, -coord_extra /)
    real(kind=dp), parameter :: eta_nodal_extra(8)  = (/ -coord_extra, -coord_extra,  coord_extra,  coord_extra, &
                                                         -coord_extra, -coord_extra,  coord_extra,  coord_extra /)
    real(kind=dp), parameter :: zeta_nodal_extra(8) = (/ -coord_extra, -coord_extra, -coord_extra, -coord_extra, &
                                                          coord_extra,  coord_extra,  coord_extra,  coord_extra /)

    real(kind=dp), parameter :: xi_gauss_extra(8)   = (/ -gauss_extra,  gauss_extra, -gauss_extra,  gauss_extra, &
                                                         -gauss_extra,  gauss_extra, -gauss_extra,  gauss_extra /)
    real(kind=dp), parameter :: eta_gauss_extra(8)  = (/ -gauss_extra, -gauss_extra,  gauss_extra,  gauss_extra, &
                                                         -gauss_extra, -gauss_extra,  gauss_extra,  gauss_extra /)
    real(kind=dp), parameter :: zeta_gauss_extra(8) = (/ -gauss_extra, -gauss_extra, -gauss_extra, -gauss_extra, &
                                                          gauss_extra,  gauss_extra,  gauss_extra,  gauss_extra /)

end module iso_module



!***********************************************************************
! subroutine kshapefunc(kintk,ninpt,nnode,ndim,N_shape_nodal_to_gauss,B_deriv_local)

!     use precision
!     include 'aba_param.inc'

!     parameter (gaussCoord=0.577350269d0)
!     dimension N_shape_nodal_to_gauss(1, nnode),B_deriv_local(ndim,nnode),coord38(3,8)
!     real(kind=dp) :: f, g, h
!     integer :: kintk
!     data  coord38 /-1.d0, -1.d0, -1.d0, &
!                     1.d0, -1.d0, -1.d0, &
!                    -1.d0,  1.d0, -1.d0, &
!                     1.d0,  1.d0, -1.d0, &
!                    -1.d0, -1.d0,  1.d0, &
!                     1.d0, -1.d0,  1.d0, &
!                    -1.d0,  1.d0,  1.d0, &
!                     1.d0,  1.d0,  1.d0/ 
!     ! Ensure kintk is within the valid range
!     ! print *, 'hello'
! !   determine (f,g,h)
!     f = coord38(1,kintk)*gaussCoord
!     g = coord38(2,kintk)*gaussCoord
!     h = coord38(3,kintk)*gaussCoord

    
!     print *, 'f = ', f
!     print *, 'g = ', g
!     print *, 'h = ', h

! !   shape functions
!     N_shape_nodal_to_gauss(1,1) = 0.125d0*(1.d0-f)*(1.d0-g)*(1.d0-h)
!     N_shape_nodal_to_gauss(1,2) = 0.125d0*(1.d0+f)*(1.d0-g)*(1.d0-h)
!     N_shape_nodal_to_gauss(1,3) = 0.125d0*(1.d0+f)*(1.d0+g)*(1.d0-h)
!     N_shape_nodal_to_gauss(1,4) = 0.125d0*(1.d0-f)*(1.d0+g)*(1.d0-h)
!     N_shape_nodal_to_gauss(1,5) = 0.125d0*(1.d0-f)*(1.d0-g)*(1.d0+h)
!     N_shape_nodal_to_gauss(1,6) = 0.125d0*(1.d0+f)*(1.d0-g)*(1.d0+h)
!     N_shape_nodal_to_gauss(1,7) = 0.125d0*(1.d0+f)*(1.d0+g)*(1.d0+h)
!     N_shape_nodal_to_gauss(1,8) = 0.125d0*(1.d0-f)*(1.d0+g)*(1.d0+h)

! !     derivative d(Ni)/d(f)
!     B_deriv_local(1,1) = -0.125d0*(1.d0-g)*(1.d0-h)
!     B_deriv_local(1,2) =  0.125d0*(1.d0-g)*(1.d0-h)
!     B_deriv_local(1,3) =  0.125d0*(1.d0+g)*(1.d0-h)
!     B_deriv_local(1,4) = -0.125d0*(1.d0+g)*(1.d0-h)
!     B_deriv_local(1,5) = -0.125d0*(1.d0-g)*(1.d0+h)
!     B_deriv_local(1,6) =  0.125d0*(1.d0-g)*(1.d0+h)
!     B_deriv_local(1,7) =  0.125d0*(1.d0+g)*(1.d0+h)
!     B_deriv_local(1,8) = -0.125d0*(1.d0+g)*(1.d0+h)

! !     derivative d(Ni)/d(g)
!     B_deriv_local(2,1) = -0.125d0*(1.d0-f)*(1.d0-h)
!     B_deriv_local(2,2) = -0.125d0*(1.d0+f)*(1.d0-h)
!     B_deriv_local(2,3) =  0.125d0*(1.d0+f)*(1.d0-h)
!     B_deriv_local(2,4) =  0.125d0*(1.d0-f)*(1.d0-h)
!     B_deriv_local(2,5) = -0.125d0*(1.d0-f)*(1.d0+h)
!     B_deriv_local(2,6) = -0.125d0*(1.d0+f)*(1.d0+h)
!     B_deriv_local(2,7) =  0.125d0*(1.d0+f)*(1.d0+h)
!     B_deriv_local(2,8) =  0.125d0*(1.d0-f)*(1.d0+h)

! !   derivative d(Ni)/d(h)
!     B_deriv_local(3,1) = -0.125d0*(1.d0-f)*(1.d0-g)
!     B_deriv_local(3,2) = -0.125d0*(1.d0+f)*(1.d0-g)
!     B_deriv_local(3,3) = -0.125d0*(1.d0+f)*(1.d0+g)
!     B_deriv_local(3,4) = -0.125d0*(1.d0-f)*(1.d0+g)
!     B_deriv_local(3,5) =  0.125d0*(1.d0-f)*(1.d0-g)
!     B_deriv_local(3,6) =  0.125d0*(1.d0+f)*(1.d0-g)
!     B_deriv_local(3,7) =  0.125d0*(1.d0+f)*(1.d0+g)
!     B_deriv_local(3,8) =  0.125d0*(1.d0-f)*(1.d0+g)
! return
! end

!In our case, we have nnode = 8 and ndim = 3
      

! subroutine kshapefunc(kintk, ninpt, nnode, ndim, &
!                       N_shape_nodal_to_gauss, B_deriv_local)
  
! !   kintk: current integration point number
! !   ninpt: total number of integration points
! !   nnode: number of nodes per element (8 for a hexahedral element)
! !   ndim: number of spatial dimensions (3 in this case)
! !   N_shape_nodal_to_gauss: shape functions
! !   B_deriv_local: derivatives of shape functions with respect to isoparametric coordinates (xi, eta, zeta)

!     use precision
!     use iso_module
!     include 'aba_param.inc'

!     ! Output parameters
!     dimension N_shape_nodal_to_gauss(1, nnode), B_deriv_local(ndim, nnode)

!     ! Local variables
!     real(kind=dp) :: g, h, r

! !   First-order brick: 3D 8-nodes

! !   determine (xi, eta, zeta)
!     g = xi_gauss_inter(kintk)
!     h = eta_gauss_inter(kintk)
!     r = zeta_gauss_inter(kintk)

!     print *, 'kintk = ', kintk
!     print *, 'g = ', g
!     print *, 'h = ', h
!     print *, 'r = ', r

! ! https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/stm/default.htm?startat=ch03s02ath62.html
! !   shape functions 
!     ! Channge xi to g, eta to h and zeta to r
!     N_shape_nodal_to_gauss(1, 1) = 0.125d0 * (1.d0 - g) * (1.d0 - h) * (1.d0 - r)
!     N_shape_nodal_to_gauss(1, 2) = 0.125d0 * (1.d0 + g) * (1.d0 - h) * (1.d0 - r)
!     N_shape_nodal_to_gauss(1, 3) = 0.125d0 * (1.d0 + g) * (1.d0 + h) * (1.d0 - r)
!     N_shape_nodal_to_gauss(1, 4) = 0.125d0 * (1.d0 - g) * (1.d0 + h) * (1.d0 - r)
!     N_shape_nodal_to_gauss(1, 5) = 0.125d0 * (1.d0 - g) * (1.d0 - h) * (1.d0 + r)
!     N_shape_nodal_to_gauss(1, 6) = 0.125d0 * (1.d0 + g) * (1.d0 - h) * (1.d0 + r)
!     N_shape_nodal_to_gauss(1, 7) = 0.125d0 * (1.d0 + g) * (1.d0 + h) * (1.d0 + r)
!     N_shape_nodal_to_gauss(1, 8) = 0.125d0 * (1.d0 - g) * (1.d0 + h) * (1.d0 + r)

!     ! print *, 'N_shape_nodal_to_gauss = ', N_shape_nodal_to_gauss

! !   derivative d(Ni)/d(xi)
!     B_deriv_local(1, 1) = -0.125d0 * (1.d0 - h) * (1.d0 - r)
!     B_deriv_local(1, 2) =  0.125d0 * (1.d0 - h) * (1.d0 - r)
!     B_deriv_local(1, 3) =  0.125d0 * (1.d0 + h) * (1.d0 - r)
!     B_deriv_local(1, 4) = -0.125d0 * (1.d0 + h) * (1.d0 - r)
!     B_deriv_local(1, 5) = -0.125d0 * (1.d0 - h) * (1.d0 + r)
!     B_deriv_local(1, 6) =  0.125d0 * (1.d0 - h) * (1.d0 + r)
!     B_deriv_local(1, 7) =  0.125d0 * (1.d0 + h) * (1.d0 + r)
!     B_deriv_local(1, 8) = -0.125d0 * (1.d0 + h) * (1.d0 + r)


! !   derivative d(Ni)/d(xi)
!     B_deriv_local(2, 1) = -0.125d0 * (1.d0 - g) * (1.d0 - r)
!     B_deriv_local(2, 2) = -0.125d0 * (1.d0 + g) * (1.d0 - r)
!     B_deriv_local(2, 3) =  0.125d0 * (1.d0 + g) * (1.d0 - r)
!     B_deriv_local(2, 4) =  0.125d0 * (1.d0 - g) * (1.d0 - r)
!     B_deriv_local(2, 5) = -0.125d0 * (1.d0 - g) * (1.d0 + r)
!     B_deriv_local(2, 6) = -0.125d0 * (1.d0 + g) * (1.d0 + r)
!     B_deriv_local(2, 7) =  0.125d0 * (1.d0 + g) * (1.d0 + r)
!     B_deriv_local(2, 8) =  0.125d0 * (1.d0 - g) * (1.d0 + r)

! !   derivative d(Ni)/d(zeta)
!     B_deriv_local(3, 1) = -0.125d0 * (1.d0 - g) * (1.d0 - h)
!     B_deriv_local(3, 2) = -0.125d0 * (1.d0 + g) * (1.d0 - h)
!     B_deriv_local(3, 3) = -0.125d0 * (1.d0 + g) * (1.d0 + h)
!     B_deriv_local(3, 4) = -0.125d0 * (1.d0 - g) * (1.d0 + h)
!     B_deriv_local(3, 5) =  0.125d0 * (1.d0 - g) * (1.d0 - h)
!     B_deriv_local(3, 6) =  0.125d0 * (1.d0 + g) * (1.d0 - h)
!     B_deriv_local(3, 7) =  0.125d0 * (1.d0 + g) * (1.d0 + h)
!     B_deriv_local(3, 8) =  0.125d0 * (1.d0 - g) * (1.d0 + h)

! return
! end


! In our case, we have nnode = 8 and ndim = 3

! subroutine kjacobian(jelem, ndim, nnode, mcrd, coords, B_deriv_local, djac, B_deriv_global)

! !   jelem: current user element number
! !   ndim: number of spatial dimensions (3 in this case)
! !   nnode: number of nodes per element (8 for a hexahedral element)
! !   coords: nodal coordinates in the global coordinate system
! !   B_deriv_local: derivatives of shape functions with respect to isoparametric coordinates (xi, eta, zeta)
! !   djac: determinant of the Jacobian matrix
! !   B_deriv_global: derivatives of shape functions with respect to global coordinates (x, y, z)
! !   jac: Jacobian matrix
! !   inv_jac: inverse of the Jacobian matrix

!     use precision
!     use iso_module
!     include 'aba_param.inc'

!     dimension coords(mcrd, nnode), B_deriv_local(ndim, nnode), B_deriv_global(ndim, nnode), &
!                 jac(ndim, ndim), inv_jac(ndim, ndim)
!     ! real(kind=dp), dimension(mcrd, nnode) :: coords          ! Nodal coordinates in the global coordinate system
!     ! real(kind=dp), dimension(ndim, nnode) :: B_deriv_local    ! Derivatives of shape functions w.r.t. isoparametric coordinates
!     ! real(kind=dp), dimension(ndim, ndim) :: jac              ! Jacobian matrix
!     ! real(kind=dp), dimension(ndim, ndim) :: inv_jac          ! Inverse of the Jacobian matrix
!     ! real(kind=dp), dimension(ndim, nnode) :: B_deriv_global  ! Derivatives of shape functions w.r.t. global coordinates

!     ! integer :: inode, idim, jdim

!     jac = 0.d0

! !   Compute the Jacobian matrix (jac)
!     do inode = 1, nnode
!         do idim = 1, ndim
!             do jdim = 1, ndim
!                 jac(jdim, idim) = jac(jdim, idim) + B_deriv_local(jdim, inode) * coords(idim, inode)
!             end do
!         end do
!     end do

!     if (ndim == 3) then
!     !   Calculate the determinant of the Jacobian matrix (djac)
!         djac =  jac(1,1)*jac(2,2)*jac(3,3)+jac(2,1)*jac(3,2)*jac(1,3) &
!               + jac(3,1)*jac(2,3)*jac(1,2)-jac(3,1)*jac(2,2)*jac(1,3) &
!               - jac(2,1)*jac(1,2)*jac(3,3)-jac(1,1)*jac(2,3)*jac(3,2)
        
!         inv_jac(1,1)=(jac(2,2)*jac(3,3)-jac(2,3)*jac(3,2))/djac
!         inv_jac(1,2)=(jac(1,3)*jac(3,2)-jac(1,2)*jac(3,3))/djac
!         inv_jac(1,3)=(jac(1,2)*jac(2,3)-jac(1,3)*jac(2,2))/djac
!         inv_jac(2,1)=(jac(2,3)*jac(3,1)-jac(2,1)*jac(3,3))/djac
!         inv_jac(2,2)=(jac(1,1)*jac(3,3)-jac(1,3)*jac(3,1))/djac
!         inv_jac(2,3)=(jac(1,3)*jac(2,1)-jac(1,1)*jac(2,3))/djac
!         inv_jac(3,1)=(jac(2,1)*jac(3,2)-jac(2,2)*jac(3,1))/djac
!         inv_jac(3,2)=(jac(1,2)*jac(3,1)-jac(1,1)*jac(3,2))/djac
!         inv_jac(3,3)=(jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1))/djac
!         if (djac < 0.d0) then ! negative or zero jacobian
!             write(7,*) 'WARNING: element', jelem, 'has neg. Jacobian'
!         endif

!     else if (ndim == 2) then
!         djac = jac(1,1) * jac(2,2) - jac(1,2) * jac(2,1)
!         if (djac > 0.d0) then ! jacobian is positive - o.k.
!             inv_jac(1,1) =  jac(2,2)/djac
!             inv_jac(2,2) =  jac(1,1)/djac
!             inv_jac(1,2) = -jac(1,2)/djac
!             inv_jac(2,1) = -jac(2,1)/djac
!         else ! negative or zero jacobian
!             write(7,*) 'WARNING: element', jelem, 'has neg. Jacobian'
!         endif

!     endif

! !   Compute the derivatives of shape functions with respect to global coordinates (B_deriv_global)
!     B_deriv_global = matmul(inv_jac, B_deriv_local)

! return
! end

!***********************************************************************

! subroutine kbmatrix(B_deriv_global,ntens,nnode,ndim,Bu_matrix)

! !   Notation, strain tensor: e11, e22, e33, e12, e13, e23
!     include 'aba_param.inc'

!     dimension B_deriv_global(ndim,nnode),Bu_matrix(ntens,nnode*ndim)

!     Bu_matrix=0.d0
!     do inode=1,nnode
!         Bu_matrix(1,ndim*inode-ndim+1) = B_deriv_global(1,inode)
!         Bu_matrix(2,ndim*inode-ndim+2) = B_deriv_global(2,inode)
!         Bu_matrix(4,ndim*inode-ndim+1) = B_deriv_global(2,inode)
!         Bu_matrix(4,ndim*inode-ndim+2) = B_deriv_global(1,inode)
!         if (ndim == 3) then
!             Bu_matrix(3,ndim*inode) = B_deriv_global(3,inode)
!             Bu_matrix(5,ndim*inode-2) = B_deriv_global(3,inode)
!             Bu_matrix(5,ndim*inode) = B_deriv_global(1,inode)
!             Bu_matrix(6,ndim*inode-1) = B_deriv_global(3,inode)
!             Bu_matrix(6,ndim*inode) = B_deriv_global(2,inode)
!         endif
!     end do

! return
! end


!   *****************************************************************

subroutine kstatevar(kintk,nsvint,svars,statev_local,icopy)

!   Transfer data to/from element-level state variable array from/to
!   material-point level state variable array.

    include 'aba_param.inc'

    dimension svars(*),statev_local(*)


    isvinc = (kintk-1)*nsvint     ! integration point increment

    if (icopy == 1) then ! Prepare arrays for entry into umat
        do i=1, nsvint
            statev_local(i) = svars(i+isvinc)
        end do
    else ! Update element state variables upon return from umat
        do i=1, nsvint
            svars(i+isvinc) = statev_local(i)
        end do
    end if

return
end

!***********************************************************************

! Important: Since UEL is called for each element, the variables defined in the UEL subroutine are local to each element.
! Therefore, we do not need the common block to store coordinates and hydrostatic stress gradient like in UMAT + UMATHT approach

subroutine UEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars, &
    props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime, &
    kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf, &
    lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period)

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

    real(kind=dp) :: E, nu, length_scale, Gc0, xkap, chi_DFT, R, T, VH, DL, K, &
                     thetaL, Gc_theta, psi, history, delta_g_b0, CL_mole_fraction, softened_factor
    real(kind=dp) :: xi_node, eta_node, zeta_node, xi_gauss, eta_gauss, zeta_gauss, djac
    
    integer, parameter :: ndim=3, ntens=6, ndi=3, nshr=3, ninpt=8, nsvint=15, ndof=5 ! C3D8 element
    ! integer, parameter :: eqplas_index = 13, sig_H_index = 14, phi_index = 15,  &
    !                       history_index = 16, rho_d_index = 17, C_index = 18, CL_index = 19, CT_index = 20
    integer, parameter :: phi_index = 13, sig_H_index = 14, CL_index = 15
    real(kind=dp), parameter :: jitter = 0.0 
    ! integer :: kintk, inode, i, j, k1, k2, isvinc

    real(kind=dp), parameter :: molar_mass_H = 1.00784d0 ! g/mol
    real(kind=dp), parameter :: molar_mass_Fe = 55.845d0 ! g/mol
    ! mcrd: Number of spatial dimensions (maximum number of coordinates) (3 for 3D)
    ! mcrd is usually equal to ndim
    ! ntens: Number of stress components
    ! ndi: Number of principal commponents
    ! nshr: Number of shear components
    ! ninpt: Number of integration points
    ! nsvint: Number of state variables per integration point
    ! nsvinc: current state variable index for the integration point defined as (i-1)*nsvint 
    !        where i is the integration point number (1 to ninpt)
    
    integer :: start_u_idx, start_phi_idx, start_CL_idx
    integer :: end_u_idx, end_phi_idx, end_CL_idx

    ! The following data is not part of UEL, defined by the user    
    real(kind=dp), dimension(ndim*nnode) :: u_prev, du_prev
    real(kind=dp), dimension(nnode) :: phi_prev, dphi_prev
    real(kind=dp), dimension(nnode) :: CL_prev, dCL_prev
    
    !real(kind=dp), dimension(ninpt) :: weight                       ! Weights for integration points
    real(kind=dp), dimension(1,nnode) :: N_shape_nodal_to_gauss     ! Shape function that interpolates from nodal points to integration points
                                                                    ! The extra 1 dimension is for matrix multiplication, otherwise it would be a vector
    real(kind=dp), dimension(ninpt) :: N_shape_gauss_to_nodal       ! Shape function that extrapolates from integration points to nodal points
    real(kind=dp), dimension(ndim,nnode) :: B_deriv_local           ! Derivatives of N_shape_nodal_to_gauss with respect to isoparametric coordinates
    real(kind=dp), dimension(ndim,nnode) :: B_deriv_global          ! Derivatives of N_shape_nodal_to_gauss with respect to global coordinates
                                                                    ! This is the collection of vectors with spatial derivatives of N_shape_nodal_to_gauss
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
    real(kind=dp), dimension(nnode,1) :: softened_sig_H             ! Softened hydrostatic stress at each node
    real(kind=dp), dimension(ndim, 1) :: softened_grad_sig_H    ! Softened hydrostatic stress gradient at each node
    real(kind=dp), dimension(nsvint) :: statev_element              ! Local state variables of the current element jelem
    real(kind=dp), dimension(ndim, ndim) :: jac, inv_jac            ! Jacobian and its inverse

    logging = 0

    if (logging == 1) then
        print *, 'UEL: jelem: ', jelem
        print *, 'UEL: nnode: ', nnode
        print *, 'UEL: ndim: ', ndim
        print *, 'UEL: ntens: ', ntens
        print *, 'UEL: ndi: ', ndi
        print *, 'UEL: nshr: ', nshr
        print *, 'UEL: ninpt: ', ninpt
        print *, 'UEL: nsvint: ', nsvint
        print *, 'UEL: ndof: ', ndof
    end if


    ! Define the starting and ending indices of each dof in u, du, v, a
    start_u_idx = 1 ! 1
    end_u_idx = ndim * nnode ! 3*8 = 24

    start_phi_idx = ndim * nnode + 1 ! 3*8 + 1 = 25
    end_phi_idx = (ndim+1) * nnode ! 4*8 = 32

    start_CL_idx = (ndim+1) * nnode + 1 ! 4*8 + 1 = 33
    end_CL_idx = (ndim+2) * nnode ! 5*8 = 40

    if (logging == 1) then
        print *, 'start_u_idx: ', start_u_idx
        print *, 'start_phi_idx: ', start_phi_idx
        print *, 'start_CL_idx: ', start_CL_idx
        print *, 'end_u_idx: ', end_u_idx
        print *, 'end_phi_idx: ', end_phi_idx
        print *, 'end_CL_idx: ', end_CL_idx
    end if
    ! Extract from the variable u and du
    u_prev(1:ndim*nnode)  = u(start_u_idx:end_u_idx)
    phi_prev(1:nnode)     = u(start_phi_idx:end_phi_idx)
    CL_prev(1:nnode)      = u(start_CL_idx:end_CL_idx)

    if (logging == 1) then
        print *,  'UEL: u_prev: ', u_prev
        print *,  'UEL: phi_prev: ', phi_prev
        print *,  'UEL: CL_prev: ', CL_prev
    end if

    du_prev(1:ndim*nnode) = du(start_u_idx:end_u_idx, 1)
    dphi_prev(1:nnode)    = du(start_phi_idx:end_phi_idx, 1)
    dCL_prev(1:nnode)     = du(start_CL_idx:end_CL_idx, 1)

    if (logging == 1) then
        print *,  'UEL: du_prev: ', du_prev
        print *,  'UEL: dphi_prev: ', dphi_prev
        print *,  'UEL: dCL_prev: ', dCL_prev
    end if
        
!   initialising
    do k1=1, ndofel
        rhs(k1,1) = 0.d0
    end do
    amatrx = 0.d0

!   nelem is from the common_block module
!   find total number of elements and stored it to nelem       
    if (dtime == 0.d0) then
        if (jelem == 1) then
            nelem = jelem
        else
            if (jelem >= nelem) then
                nelem = jelem 
            end if
        end if 
    end if      
      

! Mechanical parameters are from 1 to 8
! Phase field parameters are from 9 to 16
! Hydrogen diffusion parameters are from 17 to 40
! flow curve parameters are from 41 to 40 + (num points on flow curve) * 2
! If flow curve has 100 points, then flow curve parameters are from 41 to 240

    E = props(1) ! Young's modulus
    nu = props(2) ! Poisson's ratio

    length_scale = props(9)
    Gc0 = props(10) ! Critical energy release rate in the absence of hydrogen
    xkap = props(11) ! Well-conditioning parameter
    chi_DFT = props(12) ! Fitting slope to the DFT data
    delta_g_b0 = props(13) ! Gibbs free energy difference between the decohering interface and the surrounding material
! Hydrogen diffusion parameters are from 17 to 31

    R = props(17) ! Universal gas constant R (N*m)/(mol*K))
    T = props(18) ! Temperature (K)
    VH = props(19) ! Molar volume of H (m^3/mol)
    DL = props(20) ! Diffusion coefficient of lattice species hydrogen

    if (logging == 1) then
        print *, 'E: ', E
        print *, 'nu: ', nu
        print *, 'lc: ', length_scale
        print *, 'Gc0: ', Gc0
        print *, 'xkap: ', xkap
        print *, 'chi_DFT: ', chi_DFT
        print *, 'delta_g_b0: ', delta_g_b0
        print *, 'R: ', R
        print *, 'T: ', T
        print *, 'VH: ', VH
        print *, 'DL: ', DL
    end if

    ! The idea is to compute the hydrostatic stress (sig_H) at each nodal point by integrating 
    ! over the Gauss points and extrapolate them to nodal points. 
    ! We calculate the shape functions at each Gauss point 
    ! using the nodal coordinates and then use these to update the hydrostatic stress.
    
    ! VERIFIED

    !   compute the hydrostatic stress
    sig_H = 0.d0
    
    do inode=1, nnode
        !print *, 'inode: ', inode
        xi_node = xi_nodal_extra(inode)
        eta_node = eta_nodal_extra(inode)
        zeta_node = zeta_nodal_extra(inode)
        
        if (logging == 1) then
            print *, 'xi_node: ', xi_node
            print *, 'eta_node: ', eta_node
            print *, 'zeta_node: ', zeta_node
        end if

        ! Change them to xi_node, eta_node, zeta_node
        N_shape_gauss_to_nodal(1) = 0.125d0 * (1.d0 - xi_node) * (1.d0 - eta_node) * (1.d0 - zeta_node)
        N_shape_gauss_to_nodal(2) = 0.125d0 * (1.d0 + xi_node) * (1.d0 - eta_node) * (1.d0 - zeta_node)
        N_shape_gauss_to_nodal(3) = 0.125d0 * (1.d0 - xi_node) * (1.d0 + eta_node) * (1.d0 - zeta_node)
        N_shape_gauss_to_nodal(4) = 0.125d0 * (1.d0 + xi_node) * (1.d0 + eta_node) * (1.d0 - zeta_node)
        N_shape_gauss_to_nodal(5) = 0.125d0 * (1.d0 - xi_node) * (1.d0 - eta_node) * (1.d0 + zeta_node)
        N_shape_gauss_to_nodal(6) = 0.125d0 * (1.d0 + xi_node) * (1.d0 - eta_node) * (1.d0 + zeta_node)
        N_shape_gauss_to_nodal(7) = 0.125d0 * (1.d0 - xi_node) * (1.d0 + eta_node) * (1.d0 + zeta_node)
        N_shape_gauss_to_nodal(8) = 0.125d0 * (1.d0 + xi_node) * (1.d0 + eta_node) * (1.d0 + zeta_node)
        
        if (logging == 1) then
            print *, 'N_shape_gauss_to_nodal: ', N_shape_gauss_to_nodal
        end if

        do i = 1, ninpt
            isvinc = (i-1) * nsvint
            sig_H(inode,1) = sig_H(inode,1) + &
                            N_shape_gauss_to_nodal(i) * svars(isvinc + sig_H_index)
        end do
    end do    

    if (logging == 1) then
        print *, 'sig_H: ', sig_H
    end if

      
    ! ******************************************!
    ! Start looping over all integration points !
    ! ******************************************!

    ! print *, 'ninpt: ', ninpt
    
    do kintk=1, ninpt ! (8 integration points)
        !print *, 'kintk: ', kintk
    !   evaluate shape functions and derivatives
    !   Update N_shape_nodal_to_gauss and B_deriv_local
        !print *, 'hello'
        ! call kshapefunc(kintk, ninpt, nnode, ndim, N_shape_nodal_to_gauss, B_deriv_local)    
    !        kshapefunc(kintk,ninpt,nnode,ndim,N_shape_nodal_to_gauss,B_deriv_local)

        xi_gauss = xi_gauss_inter(kintk)
        eta_gauss = eta_gauss_inter(kintk)
        zeta_gauss = zeta_gauss_inter(kintk)

        if (logging == 1) then
            print *, 'xi_gauss: ', xi_gauss
            print *, 'eta_gauss: ', eta_gauss
            print *, 'zeta_gauss: ', zeta_gauss
        end if

    !   shape functions

        N_shape_nodal_to_gauss(1,1) = 0.125d0*(1.d0-xi_gauss)*(1.d0-eta_gauss)*(1.d0-zeta_gauss)
        N_shape_nodal_to_gauss(1,2) = 0.125d0*(1.d0+xi_gauss)*(1.d0-eta_gauss)*(1.d0-zeta_gauss)
        N_shape_nodal_to_gauss(1,3) = 0.125d0*(1.d0+xi_gauss)*(1.d0+eta_gauss)*(1.d0-zeta_gauss)
        N_shape_nodal_to_gauss(1,4) = 0.125d0*(1.d0-xi_gauss)*(1.d0+eta_gauss)*(1.d0-zeta_gauss)
        N_shape_nodal_to_gauss(1,5) = 0.125d0*(1.d0-xi_gauss)*(1.d0-eta_gauss)*(1.d0+zeta_gauss)
        N_shape_nodal_to_gauss(1,6) = 0.125d0*(1.d0+xi_gauss)*(1.d0-eta_gauss)*(1.d0+zeta_gauss)
        N_shape_nodal_to_gauss(1,7) = 0.125d0*(1.d0+xi_gauss)*(1.d0+eta_gauss)*(1.d0+zeta_gauss)
        N_shape_nodal_to_gauss(1,8) = 0.125d0*(1.d0-xi_gauss)*(1.d0+eta_gauss)*(1.d0+zeta_gauss)

        if (logging == 1) then
            print *, 'N_shape_nodal_to_gauss: ', N_shape_nodal_to_gauss
        end if
       
    !     derivative d(Ni)/d(f)

        B_deriv_local(1,1) = -0.125d0*(1.d0-eta_gauss)*(1.d0-zeta_gauss)
        B_deriv_local(1,2) =  0.125d0*(1.d0-eta_gauss)*(1.d0-zeta_gauss)
        B_deriv_local(1,3) =  0.125d0*(1.d0+eta_gauss)*(1.d0-zeta_gauss)
        B_deriv_local(1,4) = -0.125d0*(1.d0+eta_gauss)*(1.d0-zeta_gauss)
        B_deriv_local(1,5) = -0.125d0*(1.d0-eta_gauss)*(1.d0+zeta_gauss)
        B_deriv_local(1,6) =  0.125d0*(1.d0-eta_gauss)*(1.d0+zeta_gauss)
        B_deriv_local(1,7) =  0.125d0*(1.d0+eta_gauss)*(1.d0+zeta_gauss)
        B_deriv_local(1,8) = -0.125d0*(1.d0+eta_gauss)*(1.d0+zeta_gauss)

    !     derivative d(Ni)/d(g)

        B_deriv_local(2,1) = -0.125d0*(1.d0-xi_gauss)*(1.d0-zeta_gauss)
        B_deriv_local(2,2) = -0.125d0*(1.d0+xi_gauss)*(1.d0-zeta_gauss)
        B_deriv_local(2,3) =  0.125d0*(1.d0+xi_gauss)*(1.d0-zeta_gauss)
        B_deriv_local(2,4) =  0.125d0*(1.d0-xi_gauss)*(1.d0-zeta_gauss)
        B_deriv_local(2,5) = -0.125d0*(1.d0-xi_gauss)*(1.d0+zeta_gauss)
        B_deriv_local(2,6) = -0.125d0*(1.d0+xi_gauss)*(1.d0+zeta_gauss)
        B_deriv_local(2,7) =  0.125d0*(1.d0+xi_gauss)*(1.d0+zeta_gauss)
        B_deriv_local(2,8) =  0.125d0*(1.d0-xi_gauss)*(1.d0+zeta_gauss)

    !   derivative d(Ni)/d(h)

        B_deriv_local(3,1) = -0.125d0*(1.d0-xi_gauss)*(1.d0-eta_gauss)
        B_deriv_local(3,2) = -0.125d0*(1.d0+xi_gauss)*(1.d0-eta_gauss)
        B_deriv_local(3,3) = -0.125d0*(1.d0+xi_gauss)*(1.d0+eta_gauss)
        B_deriv_local(3,4) = -0.125d0*(1.d0-xi_gauss)*(1.d0+eta_gauss)
        B_deriv_local(3,5) =  0.125d0*(1.d0-xi_gauss)*(1.d0-eta_gauss)
        B_deriv_local(3,6) =  0.125d0*(1.d0+xi_gauss)*(1.d0-eta_gauss)
        B_deriv_local(3,7) =  0.125d0*(1.d0+xi_gauss)*(1.d0+eta_gauss)
        B_deriv_local(3,8) =  0.125d0*(1.d0-xi_gauss)*(1.d0+eta_gauss)

        if (logging == 1) then
            print *, 'B_deriv_local: ', B_deriv_local
        end if

    
    !   Update djac and B_deriv_global
        ! call kjacobian(jelem, ndim, nnode, mcrd, coords, B_deriv_local, &
        !               djac, B_deriv_global)
        
        jac = 0.d0

    !   Compute the Jacobian matrix (jac)
        do inode = 1, nnode
            do idim = 1, ndim
                do jdim = 1, ndim
                    jac(jdim, idim) = jac(jdim, idim) + B_deriv_local(jdim, inode) * coords(idim, inode)
                end do
            end do
        end do

        !   Calculate the determinant of the Jacobian matrix (djac)
        djac =  jac(1,1)*jac(2,2)*jac(3,3)+jac(2,1)*jac(3,2)*jac(1,3) &
            + jac(3,1)*jac(2,3)*jac(1,2)-jac(3,1)*jac(2,2)*jac(1,3) &
            - jac(2,1)*jac(1,2)*jac(3,3)-jac(1,1)*jac(2,3)*jac(3,2)
        
        inv_jac(1,1)=(jac(2,2)*jac(3,3)-jac(2,3)*jac(3,2))/djac
        inv_jac(1,2)=(jac(1,3)*jac(3,2)-jac(1,2)*jac(3,3))/djac
        inv_jac(1,3)=(jac(1,2)*jac(2,3)-jac(1,3)*jac(2,2))/djac
        inv_jac(2,1)=(jac(2,3)*jac(3,1)-jac(2,1)*jac(3,3))/djac
        inv_jac(2,2)=(jac(1,1)*jac(3,3)-jac(1,3)*jac(3,1))/djac
        inv_jac(2,3)=(jac(1,3)*jac(2,1)-jac(1,1)*jac(2,3))/djac
        inv_jac(3,1)=(jac(2,1)*jac(3,2)-jac(2,2)*jac(3,1))/djac
        inv_jac(3,2)=(jac(1,2)*jac(3,1)-jac(1,1)*jac(3,2))/djac
        inv_jac(3,3)=(jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1))/djac
        if (djac < 0.d0) then ! negative or zero jacobian
            write(7,*) 'WARNING: element', jelem, 'has neg. Jacobian'
        endif

    !   Compute the derivatives of shape functions with respect to global coordinates (B_deriv_global)
        B_deriv_global = matmul(inv_jac, B_deriv_local)


        if (logging == 1) then
            print *, 'B_deriv_global: ', B_deriv_global
            print *, 'djac: ', djac
        end if
        ! kjacobian(jelem, ndim, nnode, mcrd, coords, B_deriv_local, djac, B_deriv_global)

        dvol = weight(kintk) * djac

        Bu_matrix=0.d0
        do inode=1,nnode
            bC(1, inode) = B_deriv_global(1, inode) ! ∂N/∂x
            bC(2, inode) = B_deriv_global(2, inode) ! ∂N/∂y
            bC(3, inode) = B_deriv_global(3, inode) ! ∂N/∂z
            Bu_matrix(1,ndim*inode-ndim+1) = B_deriv_global(1,inode)
            Bu_matrix(2,ndim*inode-ndim+2) = B_deriv_global(2,inode)
            Bu_matrix(3,ndim*inode) = B_deriv_global(3,inode)
            Bu_matrix(4,ndim*inode-ndim+1) = B_deriv_global(2,inode)
            Bu_matrix(4,ndim*inode-ndim+2) = B_deriv_global(1,inode)
            Bu_matrix(5,ndim*inode-2) = B_deriv_global(3,inode)
            Bu_matrix(5,ndim*inode) = B_deriv_global(1,inode)
            Bu_matrix(6,ndim*inode-1) = B_deriv_global(3,inode)
            Bu_matrix(6,ndim*inode) = B_deriv_global(2,inode)
        end do
        
        if (logging == 1) then
            print *, 'Bu_matrix: ', Bu_matrix
        end if
    
        !   compute from nodal values
        phi = 0.d0
        CL = 0.d0
        do inode = 1,nnode
            phi = phi + N_shape_nodal_to_gauss(1,inode) * phi_prev(inode)
            CL = CL + N_shape_nodal_to_gauss(1,inode) * CL_prev(inode)
        end do   
        
        if (logging == 1) then
            print *, 'phi: ', phi
            print *, 'CL: ', CL
        end if

        ! Phase field should be between 0 and 1
        if (phi > 1.d0) then 
            phi = 1.d0
        end if

        ! ! think about this code
        ! call UMATHT_H2_diffusion(u,dudt,dudg,flux,dfdt,dfdg, &
        !         statev,temp,dtemp,dtemdx,time,dtime,predef,dpred, &
        !         cmname,ntgrd,nstatv,props,nprops,coords,pnewdt, &
        !         noel,npt,layer,kspt,kstep,kinc)
       
        ! Hydrogen contribution  

        ! Equilibrium constant
        K = exp(-delta_g_b0/(R*T))
        
        ! CL is now in wtppm (weight parts per million)
        ! with CL_mole_fraction given in units of impurity mole fraction
        CL_mole_fraction = (CL / (1000000.d0)) * (molar_mass_Fe/ molar_mass_H) ! CL in wtppm to mole fraction
        
        ! Langmuir–McLean isotherm to compute the surface coverage from the bulk hydrogen concentration CL
        ! Equation 15
        thetaL = CL_mole_fraction/(CL_mole_fraction + K)
        

        ! print *, 'thetaL: ', thetaL
        ! Equivalently, one can define the critical energy release rate dependence on the hydrogen coverage
        ! Based on fitting the DFT data
        ! Equation 9 
        Gc_theta = Gc0 * (1.d0 - chi_DFT * thetaL)

        if (logging == 1) then
            print *, 'K: ', K
            print *, 'CL_mole_fraction: ', CL_mole_fraction
            print *, 'thetaL: ', thetaL
            print *, 'Gc_theta: ', Gc_theta
        end if

           
        ! Compute the increment of strain 
        dstran = matmul(Bu_matrix, du_prev) ! Dimension is (ntens, nnode x ndim) x (nnode x ndim) = (ntens)
        
        if (logging == 1) then            
            print *, 'dstran: ', dstran
        end if
        
        ! Moving svars to statev_element
        call kstatevar(kintk, nsvint, svars, statev_element, 1)
        
        ! Obtaining the stress and strain
        stress = statev_element(1:ntens)
        stran(1:ntens) = statev_element((ntens+1):(2*ntens))

        if (logging == 1) then
            print *, 'stress: ', stress
            print *, 'stran: ', stran
        end if

        ! Recover history variables
        history_n = statev_element(CL_index)
        phi_n = statev_element(phi_index)

        if (dtime == 0.d0) then
            phi_n = phi
        end if

        if (logging == 1) then
            print *, 'history_n: ', history_n
            print *, 'phi_n: ', phi_n
        end if
    !   compute strain energy density from the previous increment   
    !   Equation 26    
        psi = 0.d0
        do k1=1, ntens
            psi = psi + stress(k1) * stran(k1) * 0.5d0
        end do

        if (logging == 1) then
            print *, 'psi: ', psi
        end if
    !   call umat to obtain stresses and constitutive matrix 
        call UMAT_elastic(props,ddsdde,stress,dstran,ntens,ndi,nshr,statev_element)
            !UMAT_elastic(props,ddsdde,stress,dstran,ntens,ndi,nshr,statev)

        stran = stran + dstran
       
        if (logging == 1) then
            print *, 'ddsdde: ', ddsdde
            print *, 'stress: ', stress
            print *, 'dstran: ', dstran
            print *, 'stran: ', stran
        end if


    !   Enforcing Karush-Kuhn-Tucker conditions for history variable H
    !   Equation 38
        if (psi >= history_n) then
            history = psi
        else
            history = history_n
        endif

        if (logging == 1) then
            print *, 'history_n: ', history_n
            print *, 'psi: ', psi
            print *, 'history: ', history
        end if
        
        ! statev_element(1:ntens) = stress(1:ntens)
        ! statev_element((ntens+1):(2*ntens)) = stran(1:ntens)
        ! statev_element(eqplas_index) = 0.d0
        ! statev_element(rho_d_index) = 0.d0
        ! statev_element(phi_index) = phi
        ! statev_element(sig_H_index) = (stress(1)+stress(2)+stress(3))/3.d0
        ! statev_element(history_index) = history
        ! statev_element(C_index) = 0.d0
        ! statev_element(CL_index) = CL
        ! statev_element(CT_index) = 0.d0

        statev_element(1:ntens) = stress(1:ntens)
        statev_element((ntens+1):(2*ntens)) = stran(1:ntens)
        statev_element(phi_index) = phi
        statev_element(sig_H_index) = (stress(1)+stress(2)+stress(3))/3.d0
        statev_element(CL_index) = history

        ! Moving statev_element to svars
        call kstatevar(kintk, nsvint, svars, statev_element, 0)
        ! subroutine kstatevar(kintk,nsvint,svars,statev_local,icopy)
        
        ! Equation 5
        softened_factor = (1.d0 - phi_n)**2 + xkap

        if (softened_factor > 1.d0) then
            softened_factor = 1.d0
        end if

        if (logging == 1) then
            print *, 'phi_n: ', phi_n
            print *, 'softened_factor: ', softened_factor
        end if
        ! ********************************************!
        ! DISPLACEMENT CONTRIBUTION TO amatrx AND rhs !
        ! ********************************************!

        ! 3D case
        ! 8 nodes x 3 displacement dofs ux, uy, uz = 24
        
        ! print *, 'dvol: ', dvol
        ! Equation 39 
        amatrx(start_u_idx:end_u_idx, start_u_idx:end_u_idx) = &
            amatrx(start_u_idx:end_u_idx, start_u_idx:end_u_idx) & 
                + dvol* softened_factor &
                * matmul(matmul(transpose(Bu_matrix),ddsdde), Bu_matrix)
        
        if (logging == 1) then
            print *, 'amatrx for uxyz: ', amatrx(1:24,1:24)
        end if
        
        ! 3D case
        
        ! Equation 36
        rhs(start_u_idx:end_u_idx, 1) = rhs(start_u_idx:end_u_idx, 1) &
            - dvol * softened_factor * matmul(transpose(Bu_matrix),stress)      

        if (logging == 1) then
            print *, 'rhs for uxyz: ', rhs(1:24,1) 
        end if
        
        ! *******************************************!
        ! PHASE FIELD CONTRIBUTION TO amatrx AND rhs !
        ! *******************************************!

        ! 3D case
        ! 8 nodes x 1 phase field dof = 8

        ! Equation 40
        amatrx(start_phi_idx:end_phi_idx, start_phi_idx:end_phi_idx) = &
                amatrx(start_phi_idx:end_phi_idx, start_phi_idx:end_phi_idx) &
            + dvol * (matmul(transpose(B_deriv_global),B_deriv_global) * Gc_theta * length_scale &
            + matmul(transpose(N_shape_nodal_to_gauss),N_shape_nodal_to_gauss) &
            * (Gc_theta/length_scale + 2.d0*history)) 
        
        if (logging == 1) then
            print *, 'amatrx for phi: ', amatrx(25:32,25:32)
        end if

        ! 3D case

        ! Equation 37
        rhs(start_phi_idx:end_phi_idx,1) = rhs(start_phi_idx:end_phi_idx,1) &
            - dvol*(matmul(transpose(B_deriv_global),matmul(B_deriv_global,phi_prev)) &
            * Gc_theta * length_scale + N_shape_nodal_to_gauss(1,:) &
            *((Gc_theta/length_scale+2.d0*history)*phi-2.d0*history))

        if (logging == 1) then
            print *, 'rhs for phi: ', rhs(25:32,1)
        end if
            
        ! **************************************************!
        ! HYDROGEN DIFFUSION CONTRIBUTION TO amatrx AND rhs !
        ! **************************************************!

        ! ! Diffusitivity matrix M
        ! ! Equation 46
        ! M_conc_capacity = matmul(transpose(N_shape_nodal_to_gauss),N_shape_nodal_to_gauss)/DL
        ! ! B_deriv_global shape is (ndim, nnode)
        ! ! sig_H has shape (nnode, 1)
        ! softened_sig_H = sig_H * softened_factor ! shape (nnode, 1)
        ! softened_grad_sig_H = matmul(B_deriv_global, softened_sig_H) 
        ! ! shape (ndim, nnode) x (nnode, 1) = (ndim, 1)
        ! ! Equation 45
        ! K_diffusitivity = matmul(transpose(B_deriv_global),B_deriv_global) &
        !                 - VH/(R*T) * matmul(B_deriv_global, matmul(softened_grad_sig_H, N_shape_nodal_to_gauss))
        
        M_conc_capacity = matmul(transpose(N_shape_nodal_to_gauss),N_shape_nodal_to_gauss)/DL
        BB = matmul(transpose(bC),bC)
        K_diffusitivity = BB - VH/(R*T) * matmul(BC, matmul(sig_H * softened_factor, N_shape_nodal_to_gauss))
        
        if (logging == 1) then
            print *, 'M_conc_capacity: ', M_conc_capacity
            ! print *, 'softened_sig_H: ', softened_sig_H
            print *, 'K_diffusitivity: ', K_diffusitivity
        end if        


        ! 3D case
        ! 8 nodes x 1 hydrogen concentration dof = 8

        if (logging == 1) then
            print *, 'M_conc_capacity: ', M_conc_capacity
            print *, 'K_diffusitivity: ', K_diffusitivity
        end if

        ! Equation 48
        amatrx(start_CL_idx:end_CL_idx,start_CL_idx:end_CL_idx) = &
            amatrx(start_CL_idx:end_CL_idx,start_CL_idx:end_CL_idx) & 
            + dvol * (M_conc_capacity/(dtime + jitter) + K_diffusitivity)
        
        if (logging == 1) then
            print *, 'amatrx for CL: ', amatrx(33:40,33:40) 
        end if            
        
        ! 3D case

        ! Equation 44
        rhs(start_CL_idx:end_CL_idx,1) = rhs(start_CL_idx:end_CL_idx,1) &
            - dvol * (matmul(K_diffusitivity,CL_prev) + &
            matmul(M_conc_capacity, dCL_prev)/(dtime + jitter))

        if (logging == 1) then
            print *, 'rhs for CL: ', rhs(33:40,1)
        end if
        ! Saving statev_element

        ! user_vars(jelem,1:ntens, kintk) = statev_element(1:ntens) * softened_factor
        ! user_vars(jelem,ntens+1:2*ntens,kintk) = statev_element((ntens+1):(2*ntens))
        ! user_vars(jelem,eqplas_index,kintk) = statev_element(eqplas_index)
        ! user_vars(jelem,rho_d_index,kintk) = statev_element(rho_d_index)
        ! user_vars(jelem,phi_index,kintk) = statev_element(phi_index)
        ! user_vars(jelem,sig_H_index,kintk) = statev_element(sig_H_index)
        ! user_vars(jelem,history_index,kintk) = statev_element(history_index)
        ! user_vars(jelem,C_index,kintk) = statev_element(C_index)
        ! user_vars(jelem,CL_index,kintk) = statev_element(CL_index)
        ! user_vars(jelem,CT_index,kintk) = statev_element(CT_index)

        ! if (kinc == 20) then
        !     print *, 'statev_element(1:ntens): ', statev_element(1:ntens)
        !     print *, 'statev_element((ntens+1):(2*ntens)): ', statev_element((ntens+1):(2*ntens))
        !     print *, 'statev_element(phi_index): ', statev_element(phi_index)
        !     print *, 'statev_element(sig_H_index): ', statev_element(sig_H_index)
        !     print *, 'statev_element(CL_index): ', statev_element(CL_index)
        ! end if
        user_vars(jelem,1:ntens, kintk) = statev_element(1:ntens) * softened_factor
        user_vars(jelem,ntens+1:2*ntens,kintk) = statev_element((ntens+1):(2*ntens))
        user_vars(jelem,phi_index,kintk) = statev_element(phi_index)
        user_vars(jelem,sig_H_index,kintk) = statev_element(sig_H_index)
        user_vars(jelem,CL_index,kintk) = CL

        ! print *, "Okay here"
    end do       ! end loop on material integration points
    
    ! print *, "Okay here outside loop"
return
end

! This is isotropic von Mises plasticity model

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

!*****************************************************************

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

    !write(*,*) 'UMAT: noel = ', noel, ' npt = ', npt
    ddsdde = 0.0d0
    noffset = noel - nelem    
    ! nelem: number of elements of UEL: [1, nelem]
    ! noel: number of elements of UMAT: [nelem + 1, 2 * nelem]
    ! => noffset: number of elements of UMAT offset by nelem: [nelem + 1, 2 * nelem] - nelem = [1, nelem]
    ! print *, 'UMAT: noel = ', noel, ' nelem = ', nelem
    ! print *, 'UMAT: noffset = ', noffset
    statev(1:nstatv) = user_vars(noffset,1:nstatv,npt)
    !  write(*,*) 'UMAT: statev = ', statev
return
end

