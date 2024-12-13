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