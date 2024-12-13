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