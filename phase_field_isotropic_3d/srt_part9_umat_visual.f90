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