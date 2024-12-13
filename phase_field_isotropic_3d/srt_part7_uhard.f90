!***********************************************************************

subroutine UHARD(syield, hard, eqplas, eqplasrt, time, dtime, temp, &
                dtemp, noel, npt, layer, kspt, kstep, kinc, cmname, nstatv, &
                statev, numfieldv, predef, dpred, nvalue, table)

    include 'aba_param.inc'

    character*80 cmname
    dimension hard(3),statev(nstatv),time(*), &
          predef(numfieldv),dpred(*)
    
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

    dimension table(2, nvalue)
    
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
