!***********************************************************************

! ! In our case, we have nnode = 8 and ndim = 3

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

!     real(kind=dp), parameter :: gauss_coord = 1.d0/sqrt(3.d0)
!     dimension N_shape_nodal_to_gauss(1, nnode), B_deriv_local(ndim, nnode)
      
     
! !   First-order brick: 3D 8-nodes

! !   determine (xi, eta, zeta)
!     g = xi_gauss_inter(kintk)
!     h = eta_gauss_inter(kintk)
!     r = zeta_gauss_inter(kintk)

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