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