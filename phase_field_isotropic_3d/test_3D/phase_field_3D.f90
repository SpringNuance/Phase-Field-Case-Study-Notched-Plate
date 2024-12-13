! User element subroutine for the phase field model for fracture
! including coupling with hydrogen diffusion (mass transport)     
! Quadratic quadrilateral elements, semi-implicit integration (2 fields)
! The code is distributed under a BSD license     
      
! If using this code for research or industrial purposes, please cite:
! E. Mart�nez-Pa�eda, A. Golahmar and C.F. Niordson. 
! A phase field formulation for hydrogen assisted cracking
! Computer Methods in Applied Mechanics and Engineering 342: 742-761 
! (2018) doi: 10.1016/j.cma.2018.07.021
      
! Emilio Mart�nez-Pa�eda (mail@empaneda.com)
! University of Cambridge

module precision
    use iso_fortran_env
    integer, parameter :: dp = real64
end module precision

module visualization
    use precision
    implicit none
    real(kind=dp) :: user_var(100000,15,8)
    integer nelem
    save
end module
      
subroutine UEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars, &
    props,nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime, &
    kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf, &
    lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njpro,period)

    use precision
    use visualization
    include 'aba_param.inc' !implicit real(a-h o-z)
      
    dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),svars(*), &
        energy(8),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel), &
        a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*), &
        ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*)

    parameter(ndim=3,ntens=6,ndi=3,nshr=3,ninpt=8,nsvint=15)
      
    dimension dN(1,nnode),dNdz(ndim,nnode),dNS(ninpt), &
        dNdx(ndim,nnode),b(ntens,nnode*ndim),ddsdde(ntens,ntens), &
        stress(ntens),stran(ntens),bC(ndim,nnode),xm(nnode,nnode), &
        xk(nnode,nnode),BB(nnode,nnode),SHa(nnode,1),coord38(ndim,nnode), &
        dstran(ntens),statev_local(nsvint)
      

    real(kind=dp), parameter :: wght(8) = (/1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0/)
    real(kind=dp) :: xi_node, eta_node, zeta_node
    real(kind=dp), parameter :: coord_inter = 1.0d0
    real(kind=dp), parameter :: gauss_inter = 1.0d0 / sqrt(3.0d0)
    real(kind=dp), parameter :: coord_extra = sqrt(3.0d0)
    real(kind=dp), parameter :: gauss_extra = 1.0d0

    real(kind=dp), parameter :: xi_nodal_extra(8)   = (/ -coord_extra,  coord_extra,  coord_extra, -coord_extra, &
                                                         -coord_extra,  coord_extra,  coord_extra, -coord_extra /)
    real(kind=dp), parameter :: eta_nodal_extra(8)  = (/ -coord_extra, -coord_extra,  coord_extra,  coord_extra, &
                                                         -coord_extra, -coord_extra,  coord_extra,  coord_extra /)
    real(kind=dp), parameter :: zeta_nodal_extra(8) = (/ -coord_extra, -coord_extra, -coord_extra, -coord_extra, &
                                                          coord_extra,  coord_extra,  coord_extra,  coord_extra /)

    logging = 0
    if (logging == 2) then
        print*, 'nsvint', nsvint ! 15
        print*, 'ndofel', ndofel ! 32 = 8 nodes x 4 dofs
        print*, 'nrhs', nrhs ! 1
        print*, 'nprops', nprops ! 5
        print*, 'mcrd', mcrd ! 3
        print*, 'nnode', nnode ! 8
        print*, 'mlvarx', mlvarx ! 36
        print*, 'nsvars', nsvars ! 44
    end if
    
    
!   initialising
    do k1=1,ndofel
        rhs(k1,1)=0.d0
    end do
    amatrx=0.d0
    
    if (logging == 1) then
        print *, 'u_disp', u(1:24)
        print *, 'du_disp', du(1:24,1)
        print *, 'u_phi', u(25:32)
        print *, 'du_phi', du(25:32,1)
        print *, 'u_CL', u(33:40)
        print *, 'du_CL', du(33:40,1)
    endif

!   find total number of elements and stored it to nelem       
    if (dtime.eq.0.d0) then
        if (jelem.eq.1) then
            nelem=jelem
        else
            if (jelem.gt.nelem) then
                nelem=jelem 
            end if
        endif 
    endif      
    

!   reading parameters
    xlc=props(9)
    Gc0=props(10)
    D=props(20) ! Diffusion coefficient
    !print *, 'xlc', xlc
    !print *, 'Gc0', Gc0
    !print *, 'D', D
    xkap=1.d-7 ! well-conditioning parameter
    Vh=2000.d0 ! Molar volume of H
    T=300.d0 ! Temperature
    R=8314.5d0 ! Gas constant

!   compute the hydrostatic stress
    SHa=0.d0
    coord38=0.d0
    
    
    do inod=1,nnode
        xi_node = xi_nodal_extra(inod)
        eta_node = eta_nodal_extra(inod)
        zeta_node = zeta_nodal_extra(inod)
        ! Change them to xi_node, eta_node, zeta_node
        dNS(1) = 0.125d0 * (1.d0 - xi_node) * (1.d0 - eta_node) * (1.d0 - zeta_node)
        dNS(2) = 0.125d0 * (1.d0 + xi_node) * (1.d0 - eta_node) * (1.d0 - zeta_node)
        dNS(3) = 0.125d0 * (1.d0 - xi_node) * (1.d0 + eta_node) * (1.d0 - zeta_node)
        dNS(4) = 0.125d0 * (1.d0 + xi_node) * (1.d0 + eta_node) * (1.d0 - zeta_node)
        dNS(5) = 0.125d0 * (1.d0 - xi_node) * (1.d0 - eta_node) * (1.d0 + zeta_node)
        dNS(6) = 0.125d0 * (1.d0 + xi_node) * (1.d0 - eta_node) * (1.d0 + zeta_node)
        dNS(7) = 0.125d0 * (1.d0 - xi_node) * (1.d0 + eta_node) * (1.d0 + zeta_node)
        dNS(8) = 0.125d0 * (1.d0 + xi_node) * (1.d0 + eta_node) * (1.d0 + zeta_node)
        
        ! print *, 'dNS', dNS
        do i=1,ninpt
            isvinc=(i-1)*nsvint
            SHa(inod,1)=SHa(inod,1)+dNS(i)*svars(isvinc+14)
        end do
    end do      
      
    

    do kintk=1,ninpt
!   evaluate shape functions and derivatives
        call kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)      
        ! kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
        call kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
        ! subroutine kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
        dvol=wght(kintk)*djac

        ! print *, 'dvol', dvol
       
!   form B-matrix
      b=0.d0
      do inod=1,nnode
        bC(1,inod)=dNdx(1,inod)
        bC(2,inod)=dNdx(2,inod)
        bC(3,inod)=dNdx(3,inod)
       b(1,ndim*inod-ndim+1)=dNdx(1,inod)
       b(2,ndim*inod-ndim+2)=dNdx(2,inod)
       b(4,ndim*inod-ndim+1)=dNdx(2,inod)
       b(4,ndim*inod-ndim+2)=dNdx(1,inod)
       if (ndim.eq.3) then
        b(3,ndim*inod)=dNdx(3,inod)
        b(5,ndim*inod-2)=dNdx(3,inod)
        b(5,ndim*inod)=dNdx(1,inod)
        b(6,ndim*inod-1)=dNdx(3,inod)
        b(6,ndim*inod)=dNdx(2,inod)
       endif
      end do                
       
        ! print*, 'okay'
!   compute from nodal values
        phi=0.d0
        cL=0.d0
        do inod=1,nnode
            phi=phi+dN(1,inod)*u(ndim*nnode+inod)
            cL=cL+dN(1,inod)*u((ndim+1)*nnode+inod)
        end do   
        if (phi.gt.1.d0) phi=1.d0
       
!   hydrogen contribution  
        Theta=cL*5.5d-05/(cL*5.5d-05+dexp(-3.d7/(R*T)))
        Gc=Gc0*(1.d0-0.89d0*Theta)
           
!   compute the increment of strain and recover history variables
        dstran=matmul(b,du(1:ndim*nnode,1))
        call kstatevar(kintk,nsvint,svars,statev_local,1)
        stress=statev_local(1:ntens)
        stran(1:ntens)=statev_local((ntens+1):(2*ntens))
        Hn=statev_local(2*ntens+3)
        phin=statev_local(2*ntens+1)
        if (dtime.eq.0.d0) phin=phi
       
!   compute strain energy density from the previous increment       
        Psi=0.d0
        do k1=1,ntens
            Psi=Psi+stress(k1)*stran(k1)*0.5d0
        end do
       
!   call umat to obtain stresses and constitutive matrix 
        call kumat(props,ddsdde,stress,dstran,ntens,statev_local)
        ! subroutine kumat(props,ddsdde,stress,dstran,ntens,statev)
        stran=stran+dstran
       
!   enforcing Karush-Kuhn-Tucker conditions
        if (Psi.gt.Hn) then
            H=Psi
        else
            H=Hn
        endif
       
        statev_local(1:ntens)=stress(1:ntens)
        statev_local((ntens+1):(2*ntens))=stran(1:ntens)
        statev_local(2*ntens+1)=phi
        statev_local(2*ntens+2)=(stress(1)+stress(2)+stress(3))/3.d0 !SH
        statev_local(2*ntens+3)=H
        
        call kstatevar(kintk,nsvint,svars,statev_local,0)
        
        amatrx(1:24,1:24)=amatrx(1:24,1:24)+dvol*(((1.d0-phin)**2+xkap) &
                                * matmul(matmul(transpose(b),ddsdde),b))
            
        rhs(1:24,1)=rhs(1:24,1)- dvol*(matmul(transpose(b),stress)*((1.d0-phin)**2+xkap))       
        
        if (logging == 1) then
            print *, 'amatrix disp', amatrx(1:24,1:24)
            print *, 'rhs disp', rhs(1:24,1)
        endif

        amatrx(25:32,25:32)=amatrx(25:32,25:32) &
            +dvol*(matmul(transpose(dNdx),dNdx)*Gc*xlc &
            +matmul(transpose(dN),dN)*(Gc/xlc+2.d0*H))   

        rhs(25:32,1)=rhs(25:32,1) &
            -dvol*(matmul(transpose(dNdx),matmul(dNdx,u(25:32))) &
            *Gc*xlc+dN(1,:)*((Gc/xlc+2.d0*H)*phi-2.d0*H))
        
        if (logging == 1) then
            print *, 'amatrix phi', amatrx(25:32,25:32)
            print *, 'rhs phi', rhs(25:32,1)
        endif
            
        xm=matmul(transpose(dN),dN)/D
        BB=matmul(transpose(bC),bC)
        xk=BB-Vh/(R*T)*matmul(BB,matmul((SHa*((1.d0-phin)**2+xkap)),dN))

        amatrx(33:40,33:40)=amatrx(33:40,33:40)+dvol*(xm/dtime+xk)
            
        rhs(33:40,1)=rhs(33:40,1)-dvol*(matmul(xk,u(33:40))+ &
            matmul(xm,du(33:40,1))/dtime)

        if (logging == 1) then
            print *, 'amatrix CL', amatrx(33:40,33:40)
            print *, 'rhs CL', rhs(33:40,1)
        endif

! output
        user_var(jelem,1:6,kintk)=statev_local(1:6)*((1.d0-phin)**2+xkap)
        user_var(jelem,7:14,kintk)=statev_local((ntens+1):(2*ntens+2))
        user_var(jelem,(2*ntens+3),kintk)=cL
        
    end do       ! end loop on material integration points
    
    ! print *, 'okay'
return
end
    

!***********************************************************************
      subroutine kshapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
!
      include 'aba_param.inc'
!
      parameter (gaussCoord=0.577350269d0)
      dimension dN(1, nnode),dNdz(ndim,nnode),coord38(3,8)

      data  coord38 /-1.d0, -1.d0, -1.d0, &
                      1.d0, -1.d0, -1.d0, &
                     -1.d0,  1.d0, -1.d0, &
                      1.d0,  1.d0, -1.d0, &
                     -1.d0, -1.d0,  1.d0, &
                      1.d0, -1.d0,  1.d0, &
                     -1.d0,  1.d0,  1.d0, &
                      1.d0,  1.d0,  1.d0/ 

!     determine (g,h,r)
       f=coord38(1,kintk)*gaussCoord
       g=coord38(2,kintk)*gaussCoord
       h=coord38(3,kintk)*gaussCoord

!     shape functions
       dN(1,1)=0.125d0*(1.d0-f)*(1.d0-g)*(1.d0-h)
       dN(1,2)=0.125d0*(1.d0+f)*(1.d0-g)*(1.d0-h)
       dN(1,3)=0.125d0*(1.d0+f)*(1.d0+g)*(1.d0-h)
       dN(1,4)=0.125d0*(1.d0-f)*(1.d0+g)*(1.d0-h)
       dN(1,5)=0.125d0*(1.d0-f)*(1.d0-g)*(1.d0+h)
       dN(1,6)=0.125d0*(1.d0+f)*(1.d0-g)*(1.d0+h)
       dN(1,7)=0.125d0*(1.d0+f)*(1.d0+g)*(1.d0+h)
       dN(1,8)=0.125d0*(1.d0-f)*(1.d0+g)*(1.d0+h)

!     derivative d(Ni)/d(f)
       dNdz(1,1)=-0.125d0*(1.d0-g)*(1.d0-h)
       dNdz(1,2)= 0.125d0*(1.d0-g)*(1.d0-h)
       dNdz(1,3)= 0.125d0*(1.d0+g)*(1.d0-h)
       dNdz(1,4)=-0.125d0*(1.d0+g)*(1.d0-h)
       dNdz(1,5)=-0.125d0*(1.d0-g)*(1.d0+h)
       dNdz(1,6)= 0.125d0*(1.d0-g)*(1.d0+h)
       dNdz(1,7)= 0.125d0*(1.d0+g)*(1.d0+h)
       dNdz(1,8)=-0.125d0*(1.d0+g)*(1.d0+h)

!     derivative d(Ni)/d(g)
       dNdz(2,1)=-0.125d0*(1.d0-f)*(1.d0-h)
       dNdz(2,2)=-0.125d0*(1.d0+f)*(1.d0-h)
       dNdz(2,3)= 0.125d0*(1.d0+f)*(1.d0-h)
       dNdz(2,4)= 0.125d0*(1.d0-f)*(1.d0-h)
       dNdz(2,5)=-0.125d0*(1.d0-f)*(1.d0+h)
       dNdz(2,6)=-0.125d0*(1.d0+f)*(1.d0+h)
       dNdz(2,7)= 0.125d0*(1.d0+f)*(1.d0+h)
       dNdz(2,8)= 0.125d0*(1.d0-f)*(1.d0+h)

!     derivative d(Ni)/d(h)
       dNdz(3,1)=-0.125d0*(1.d0-f)*(1.d0-g)
       dNdz(3,2)=-0.125d0*(1.d0+f)*(1.d0-g)
       dNdz(3,3)=-0.125d0*(1.d0+f)*(1.d0+g)
       dNdz(3,4)=-0.125d0*(1.d0-f)*(1.d0+g)
       dNdz(3,5)= 0.125d0*(1.d0-f)*(1.d0-g)
       dNdz(3,6)= 0.125d0*(1.d0+f)*(1.d0-g)
       dNdz(3,7)= 0.125d0*(1.d0+f)*(1.d0+g)
       dNdz(3,8)= 0.125d0*(1.d0-f)*(1.d0+g)


      return
      end

!***********************************************************************
      subroutine kjacobian(jelem,ndim,nnode,coords,dNdz,djac,dNdx,mcrd)
!     Notation: djac - Jac determinant; xjaci - inverse of Jac matrix
!     dNdx - shape functions derivatives w.r.t. global coordinates
      include 'aba_param.inc'

      dimension xjac(ndim,ndim),xjaci(ndim,ndim),coords(mcrd,nnode), &
        dNdz(ndim,nnode),dNdx(ndim,nnode)

      xjac=0.d0

      do inod=1,nnode
       do idim=1,ndim
        do jdim=1,ndim
         xjac(jdim,idim)=xjac(jdim,idim)+ &
             dNdz(jdim,inod)*coords(idim,inod)
        end do
       end do
      end do

      if (ndim.eq.3) then

       djac=xjac(1,1)*xjac(2,2)*xjac(3,3)+xjac(2,1)*xjac(3,2)*xjac(1,3) &
       +xjac(3,1)*xjac(2,3)*xjac(1,2)-xjac(3,1)*xjac(2,2)*xjac(1,3) &
       -xjac(2,1)*xjac(1,2)*xjac(3,3)-xjac(1,1)*xjac(2,3)*xjac(3,2)
       !if (djac.gt.0.d0) then ! jacobian is positive - o.k.
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

      else if (ndim.eq.2) then

       djac=xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1)
       !if (djac.gt.0.d0) then ! jacobian is positive - o.k.
        xjaci(1,1)=xjac(2,2)/djac
        xjaci(2,2)=xjac(1,1)/djac
        xjaci(1,2)=-xjac(1,2)/djac
        xjaci(2,1)=-xjac(2,1)/djac
       !else ! negative or zero jacobian
       ! write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
       !endif

      endif

      dNdx=matmul(xjaci,dNdz)

      return
      end

!***********************************************************************
      subroutine kbmatrix(dNdx,ntens,nnode,ndim,b)
!     Notation, strain tensor: e11, e22, e33, e12, e13, e23
      include 'aba_param.inc'

      dimension dNdx(ndim,nnode),b(ntens,nnode*ndim)

      b=0.d0
      do inod=1,nnode
       b(1,ndim*inod-ndim+1)=dNdx(1,inod)
       b(2,ndim*inod-ndim+2)=dNdx(2,inod)
       b(4,ndim*inod-ndim+1)=dNdx(2,inod)
       b(4,ndim*inod-ndim+2)=dNdx(1,inod)
       if (ndim.eq.3) then
        b(3,ndim*inod)=dNdx(3,inod)
        b(5,ndim*inod-2)=dNdx(3,inod)
        b(5,ndim*inod)=dNdx(1,inod)
        b(6,ndim*inod-1)=dNdx(3,inod)
        b(6,ndim*inod)=dNdx(2,inod)
       endif
      end do

      return
      end

!***********************************************************************

      subroutine kstatevar(npt,nsvint,statev,statev_ip,icopy)
!
!     Transfer data to/from element-level state variable array from/to
!     material-point level state variable array.
!
      include 'aba_param.inc'

      dimension statev(*),statev_ip(*)

      isvinc=(npt-1)*nsvint     ! integration point increment

      if (icopy.eq.1) then ! Prepare arrays for entry into umat
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

!***********************************************************************
      subroutine kumat(props,ddsdde,stress,dstran,ntens,statev)
!
!     Subroutine with the material model
!
      include 'aba_param.inc' !implicit real(a-h o-z)

      dimension props(*),ddsdde(ntens,ntens),stress(ntens),statev(*),dstran(ntens)

!     Initialization
      ddsdde=0.d0
      E=props(1) ! Young's modulus
      xnu=props(2) ! Poisson's ratio

!     Build stiffness matrix
      eg2=E/(1.d0+xnu)
      elam=(E/(1.d0-2.d0*xnu)-eg2)/3.d0

!     Update stresses
      do k1=1,3
       do k2=1,3
        ddsdde(k2,k1)=elam
       end do
       ddsdde(k1,k1)=eg2+elam
      end do
      do k1=4,ntens
       ddsdde(k1,k1)=eg2/2.d0
      end do

      stress=stress+matmul(ddsdde,dstran)

      return
      end

!***********************************************************************
subroutine UMAT(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt, &
    drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred, &
    cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt, &
    celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

    use visualization
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
    statev(1:nstatv) = user_var(noffset,1:nstatv,npt)
    !  write(*,*) 'UMAT: statev = ', statev
return
end