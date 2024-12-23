C -----------------------------------------------------------C
C Subroutine for Abaqus/Explicit for isotropic elasticity    C
C and isotropic plasticity with pressure and lode            C
C dependence, according to plasticity model by BAI, Y. and   C
C WIERZBICKI, T. (published in "A new model of metal         C
C plasticity and fracture with pressure and Lode             C
C dependence", International Journal of Plasticity 24 (2008) C
C pages 1071-1096).                                          C
C															 C
C - Damage initiation according to 3D fracture locus         C
C   (equivalent plastic strain to failure initiation being   C
C   function of both hydrostatic pressure and lode angle)    C
C - Damage evolution based on energy dissipated during       C
C   damage process                                           C
C															 C
C - Temperature softening									 C
C - Strain rate hardening									 C
C                                                            C
C Originally based on the von Mises' J2 plasticity theory.   C
C Elastic predictor, radial corrector algorithm used. With   C
C explicit forward Euler scheme for integration of flow rule.C
C                                                            C
C Not suitable for solving problem settings including plane  C
C stresses. Not suitable for 2D models.                      C
C -----------------------------------------------------------C
C Solution dependent variables (SDVs):						 C
C															 C
C SDV1: 	Equivalent plastic strain						 C
C SDV2: 	Damage variable									 C
C SDV3: 	Yield stress at damage initiation				 C
C SDV4: 	Flag (=0 element not damaged,					 C
C		      =1 element experienced brittle fracture,	 C
C		      =2 element experienced ductile damage,     C
C                     =3 element experienced ductile fracture)	 C
C SDV5-10:	Total strain tensor								 C
C SDV11-16:	Plastic strain tensor							 C
C SDV17:	Equivalent plastic strain increment				 C
C SDV18:	Temperature softening correction function		 C
C SDV19:	Mutliplicative strain rate hardening correction	 C
C			function		 								 C
C SDV20:	Equivalent plastic strain rate					 C
C SDV21:	Equivalent plastic strain rate (computed over	 C
C           a user-defined no. of increments)				 C
C SDV22:	Increment counter (for eq. plas. str. rate       C
C           computation)									 C
C SDV23:	Sum of equivalent plastic strain increments      C
C           before the computation of eq. plas. str. rate    C
C SDV24:	Sum of time increments before computation of eq. C
C           plas. str. rate									 C
C SDV25:	Additive strain rate hardening correction		 C
C			function										 C
C SDV26:	Maximum principal stress						 C
C SDV27:	Equivalent plastic strain (averaged over a		 C
C			user-defined no. of increments)					 C
C SDV28:	Element deletion flag							 C
C -----------------------------------------------------------C
C Contact info:                                              C
C															 C
C Mohamed Sharaf                                             C
C RWTH Aachen                                                C
C Institut fuer Eisenhuettenkunde                            C
C Intzestrasse 1                                             C
C 52072 Aachen                                               C
C Germany                                                    C
C mohamed.sharaf@iehk.rwth-aachen.de                         C
C                                                            C
C Aachen, 6 March 2012                                       C
C------------------------------------------------------------C
C SDV and parameters modified by B. Wu, J. He and F. Shen    C
C                                                            C
C SDV29:	Eta         								     C
C SDV30:	Thetabar										 C
C SDV31:	Ductile damage initiation indicator				 C
C SDV32:	Ductile failure indicator						 C
C SDV33:	Damage initiation strain at the current step	 C
C SDV34:	Dcrit at the current step						 C
C SDV35:	Cleavage initiation strain at the current step	 C
C SDV36:	Cleavage failure indicator based on strain		 C
C SDV37:	Damage fracture strain at the current step	     C
C SDV38:	Equivalent plastic strain at damage initiation	 C
C SDV39:        Hydrostatic stress
C -----------------------------------------------------------C
C Additional parameters props(28)-props(36)                  C
C                                                            C
C Additional parameter for cut-off value is defined			 C
C Additional parameters Ddf1-Ddf4 are defined for Ddf-locus  C
C Additional parameters Cc1-Cc4 are defined for Cleavage 	 C
C -----------------------------------------------------------C
C Modifications by B. Wu, J. He and F. Shen                  C
C                                                            C
C Damage indicator: Damage initiates when SDV31>1            C
C Failure indicator: Failure happens when SDV32>1            C
C Cleavage failure indicator: Cleavage happens when SDV36>1  C
C Non-proportional Loading paths with varying stress states  C
C are considered are considered with weighting scheme        C
C Cut-off value concept is considered                        C
C -----------------------------------------------------------C
        subroutine vumat(
C Read only
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     1  stepTime, totalTime, dt, cmname, coordMp, charLength,
     1  props, density, strainInc, relSpinInc,
     1  tempOld, stretchOld, defgradOld, fieldOld,
     1  stressOld, stateOld, enerInternOld, enerInelasOld,
     1  tempNew, stretchNew, defgradNew, fieldNew,
C Write only
     1  stressNew, stateNew, enerInternNew, enerInelasNew)
C
        include 'vaba_param.inc'
C
        dimension props(nprops), density(nblock),
     1  coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     1  relSpinInc(*), tempOld(nblock),
     1  stretchOld(*), defgradOld(*),
     1  fieldOld(*), stressOld(nblock,ndir+nshr),
     1  stateOld(nblock,nstatev), enerInternOld(nblock),
     1  enerInelasOld(nblock), tempNew(nblock),
     1  stretchNew(*), defgradNew(*), fieldNew(*),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     1  enerInternNew(nblock), enerInelasNew(nblock)
C
        character*80 cmname
C
C Defining numerical constants
        parameter( zero = 0., one = 1., two = 2., three = 3.,
     1  third = one/three, half = .5, twoThirds = two/three,
     1  threeHalves = 1.5, thousand = 1000., tol=1.D-6 )
C        
       parameter(zNL=5.1D20, zmuR_L=-19.576D6, zmu0_L=28.6D6,
     1 zT0=300., zRgas=8314.46262, zEB=60D6, zVH=2000.   )
C
C
C
C Defining material properties constants
        e      = props(1)
        xnu    = props(2)
        ceta   = props(3)
        eta0   = props(4)
        cthetas= props(5)
        cthetat= props(6)
        cthetac= props(7)
        om     = props(8)
C Ductile damage initiation locus function parameters
        Ddi1   = props(9)
        Ddi2   = props(10)
        Ddi3   = props(11)
        Ddi4   = props(12)
        Ddi5   = props(13)
        Ddi6   = props(14)
C Fracture energy dissipation (energy dissipated per 
C element unit area after ductile damage dissipation)
        gf     = props(15)
C Maximum value for ductile damage variable
C       Dcrit   = props(16)
C 5 parameters for temperature softening correction
        cT1    = props(17)
        cT2    = props(18)
        cT3    = props(19)
        eta2   = props(20)
        cp     = props(21)
        t0     = props(22)
C 3 parameters for strain rate hardening correction
        cE1    = props(23)
        cE2    = props(24)
        cE3    = props(25)
C Number of increments to update plastic strain rate
        strrInc= props(26)
C Brittle damage initiation stress value
        sigdmg = props(27)
C Cut-off value
        Cut    = props(28)
C Ductile frature strain locus function parameters
        Ddf1   = props(29)
        Ddf2   = props(30)
        Ddf3   = props(31)
        Ddf4   = props(32)
C Cleavage initiation locus
        Cc1    = props(33)
        Cc2    = props(34)
        Cc3    = props(35)
        Cc4    = props(36)
C 	PROPS(37)=Switch temerapture function on(1)/off(0)
C  	PROPS(38)=Switch strain function on(1)/off(0)
        ntempSwitch=PROPS(37)
        nstrainSwitch=PROPS(38)
        
        nvalue = (nprops-40)/2
C
C Checking for valid entries
        if (om.lt.zero) then
          write(6,5)
 5        format(//,30X,'***ERROR - m MUST BE A NON-NEGATIVE INTEGER')
        endif
C
C Defining constants
        twomu  = e / ( one + xnu )
        thremu = threeHalves * twomu
        sixmu  = three * twomu
        alamda = twomu * ( e - twomu ) / ( sixmu - two * e )
        akappa = e * twomu * half / (three * (threeHalves * twomu - e
     1           ))
        term   = one / ( twomu * ( one + hard/thremu ) )
        con1   = sqrt( twoThirds )
        pi     = 4. * atan(one)
        con2   = (sqrt(three)/two) / (one - sqrt(three)/two)
C
C Computation per material point starts here
      do 100 i = 1,nblock
C
       Sigma_Hydro_old = stateOld(i,39) 
       stateNew(i,40) = Sigma_Hydro_old

       PEEQ_old = stateOld(i,1)
       stateNew(i,41) = PEEQ_old
C If brittle fracture has previously occured, stress is zeroed
        if (stateOld(i,4).gt.half .and. stateOld(i,4).lt.(one+half))
     1    then
          stressNew(i,1) = zero
          stressNew(i,2) = zero
          stressNew(i,3) = zero
          stressNew(i,4) = zero
          stressNew(i,5) = zero
          stressNew(i,6) = zero
          stateNew(i,2) = one
          stateNew(i,4) = one
          stateNew(i,28) = zero
          goto 100
        endif
C If ductile fracture has previously occured, stress is zeroed
        if ( (stateOld(i,4).gt.(two+half)) .and. 
     1     ( stateOld(i,4).lt.(three+half)) ) then
          stressNew(i,1) = zero
          stressNew(i,2) = zero
          stressNew(i,3) = zero
          stressNew(i,4) = zero
          stressNew(i,5) = zero
          stressNew(i,6) = zero
          stateNew(i,2) = one
          stateNew(i,4) = three
          stateNew(i,32) = one
          stateNew(i,28) = zero
          goto 100
        endif
C Updating total strain (states 5-10)
        stateNew(i,5)  = stateOld(i,5)  + strainInc(i,1)
        stateNew(i,6)  = stateOld(i,6)  + strainInc(i,2)
        stateNew(i,7)  = stateOld(i,7)  + strainInc(i,3)
        stateNew(i,8)  = stateOld(i,8)  + strainInc(i,4)
        stateNew(i,9)  = stateOld(i,9)  + strainInc(i,5)
        stateNew(i,10) = stateOld(i,10) + strainInc(i,6)
C
C Trace of total strain tensor
        epsilontrace = stateNew(i,5) + stateNew(i,6) + stateNew(i,7)
C
C Trace of strain increment tensor
        trace = strainInc(i,1) + strainInc(i,2) + strainInc(i,3)
C
C Calculating softening correction term due to temperature rise 
          facT = one
          stateNew(i,18)= facT
C
C Calculating hardening correction terms due to straining rate
        if (nstrainSwitch.NE.ZERO) then
          facStrrateM = cE1*log(stateOld(i,21)) + cE2
          facStrrateA = cE3*stateOld(i,21)
        else
          facStrrateM = one
          facStrrateA = zero
        endif
C
C Preventing negative values for straining rate
        if (facStrrateM .lt. one) facStrrateM = one
        stateNew(i,19)= facStrrateM
C
C Restoring previous softening correction term due to ductile damage
        facDctlDmgDev = one - stateOld(i,2)
        facDctlDmgVol = facDctlDmgDev
        if (stateOld(i,29).lt.zero) then
          facDctlDmgVol = one
        endif
C
C Trial stress
        sig1 = stressOld(i,1) + facDctlDmgVol*alamda*trace
     1                        + facDctlDmgDev*twomu*strainInc(i,1)
        sig2 = stressOld(i,2) + facDctlDmgVol*alamda*trace
     1                        + facDctlDmgDev*twomu*strainInc(i,2)
        sig3 = stressOld(i,3) + facDctlDmgVol*alamda*trace
     1                        + facDctlDmgDev*twomu*strainInc(i,3)
        sig4 = stressOld(i,4) + facDctlDmgDev*twomu*strainInc(i,4)
        sig5 = stressOld(i,5) + facDctlDmgDev*twomu*strainInc(i,5)
        sig6 = stressOld(i,6) + facDctlDmgDev*twomu*strainInc(i,6)
C
C Following equation numbering of publication mentioned in code title text
C Equ. (4) - Calculating the deviatoric part of trial stress
        sigmean = third * ( sig1 + sig2 + sig3 )
        ds1 = sig1 - sigmean
        ds2 = sig2 - sigmean
        ds3 = sig3 - sigmean
C
C Calculating the magnitude of the deviatoric trial stress tensor
        dsmag = sqrt( ds1**2 + ds2**2 + ds3**2 + two*sig4**2 + two
     1   *sig5
     1   **2 + two*sig6**2 )
C
C Preventing a divide by zero when dsmag is zero. When dsmag is zero, computation is still in the elastic zone
        if (dsmag.lt.tol .and. dsmag.ge.zero) dsmag = tol
        if (dsmag.gt.(zero-tol) .and. dsmag.le.zero) dsmag = zero-tol
C
C Following equation numbering of publication mentioned in code title text
C Eq. (1) - Calculating the 1st invariant of the stress tensor
        p   = zero - sigmean
C
C Eq. (2) - Calculating the 2nd invariant of the stress tensor
        q   = dsmag / con1
C
C Eq. (3) - Calculating the 3rd invariant of the stress tensor
        r   = ( three/twoThirds * (ds1*(ds1**2 + sig4**2 + sig6**2)
     1  + 2.
     1  *sig4*(ds1*sig4 + ds2*sig4 + sig6*sig5) + 2.*sig6*(ds1*sig6
     1  + sig4
     1  *sig5 + ds3*sig6) + ds2*(sig4**2 + ds2**2 + sig5**2) + 2.*sig5
     1  *(sig4*sig6 + ds2*sig5 + ds3*sig5) + ds3*(sig6**2 + sig5**2
     1  + ds3**2)) )**third
C
C Eq. (5) - Calculating the dimensionless hydrostatic pressure eta
        eta = sigmean / q
        stateNew(i,29)= eta
C
C Eq. (6) - Calculating the Lode angle theta
        cosine=(r/q)**3
C       Assuring that -1<cosine<1
        if (cosine.gt.one) cosine=one
        if (cosine.lt.(zero-one)) cosine=zero-one
        theta = third*acos(cosine)
C
C Eq. (7) - Calculating the normalized Lode angle thetabar
        thetabar = one - 6.*theta/pi
        stateNew(i,30)= thetabar
C
C Eq. (12) - Calculating the Lode angle dependant parameter gamma
        gamma = con2 * (one/cos(theta-pi/6.) - one)
C
C Eq. (13) - Determining cthetaax, either tension case or compression case
        if( thetabar .ge. zero ) then
          cthetaax = cthetat
        else
          cthetaax = cthetac
        endif
C
C Fetching yield stress from flow curve
        call ahard(sigmayield,hard,stateOld(i,1),props(41),nvalue,ntempSwitch,nstrainSwitch)
C
C Eq. (11) - Calculating radius of elastic zone
        radius = (sigmayield * facStrrateM + facStrrateA) * ( one
     1  - ceta
     1  *(eta-eta0) ) * (cthetas + ((cthetaax-cthetas)* (gamma- ((gamma
     1  **(om+one))/(om+one)) ) ) ) * facT * facDctlDmgVol
C
C Eq. (15) - Checking if yielding. The yielding factor facyld is set to be zero
C in the elastic zone, one if yielding
        facyld = zero
        if( q - radius .ge. zero ) facyld = one
C
C Eq. (20) - Calculating the tensor of the normal direction of plastic flow
C       Avoiding a divide by zero by avoiding that theta = zero
        if (theta.eq.zero) theta = theta + tol
C
        fac1 = con1 * three/(two*q)

        fac2 = facT * facDctlDmgVol * con1 * (sigmayield * facStrrateM
     1  + facStrrateA)*ceta*(cthetas+((cthetaax-cthetas)*(gamma
     1  -( (gamma**(om+one)) / (om+one) )))) * three*eta/(two*(q**2))
     
        fac3 = facT * facDctlDmgVol * con1 * (sigmayield * facStrrateM
     1  + facStrrateA) * (one-ceta*(eta-eta0)) * (cthetaax-cthetas)
     1  * (one -(gamma**om)) * (three*sqrt(three)/(two-sqrt(three)))
     1  * (tan(theta-pi/6.)/cos(theta-pi/6.)) *(one/(q*sin(three
     1  *theta)))
C
        on1 = fac1*ds1
     1  - fac2 * ds1
     1  - fac3 * (one/three+cos(three*theta)*ds1/(two*q)-three*((ds1
     1  **2)+(sig4**2)+(sig6**2))/(two*(q**2)))
C
        on2 = fac1*ds2
     1  - fac2 * ds2
     1  - fac3 * (one/three+cos(three*theta)*ds2/(two*q)-three*((sig4
     1  **2)+(ds2**2)+(sig5**2))/(two*(q**2)))
C
        on3 = fac1*ds3
     1  - fac2 * ds3
     1  - fac3 * (one/three+cos(three*theta)*ds3/(two*q)-three*((sig6
     1  **2)+(sig5**2)+(ds3**2))/(two*(q**2)))
C
        on4 = fac1*sig4
     1  - fac2 * sig4
     1  - fac3 * (cos(three*theta)*sig4/(two*q)-three*((ds1*sig4)
     1  +(ds2*sig4)+(sig6*sig5))/(two*(q**2)))
C
        on5 = fac1*sig5
     1  - fac2 * sig5
     1  - fac3 * (cos(three*theta)*sig5/(two*q)-three*((sig4*sig6)
     1  +(ds2*sig5)+(ds3*sig5))/(two*(q**2)))
C
        on6 = fac1*sig6
     1  - fac2 * sig6
     1  - fac3 * (cos(three*theta)*sig6/(two*q)-three*((ds1*sig6)
     1  +(sig4*sig5)+(ds3*sig6))/(two*(q**2)))
C
C Calculating trial yield function
        phitrial   = facyld * con1 * (q - radius) / facDctlDmgVol
C
C Calculating equivalent plastic strain
        deqps  = con1 * phitrial / (twomu + twoThirds * hard)
        if (deqps.lt.zero) then
          stateNew(i,1) = stateOld(i,1)
        goto 100
        endif
          stateNew(i,17) = deqps
C
C Updating equivalent plastic strain
        stateNew(i,1) = stateOld(i,1) + deqps
C
C Updating plastic strain (states 11-16)
        stateNew(i,11) = stateOld(i,11) + deqps * on1 / con1
        stateNew(i,12) = stateOld(i,12) + deqps * on2 / con1
        stateNew(i,13) = stateOld(i,13) + deqps * on3 / con1
        stateNew(i,14) = stateOld(i,14) + deqps * on4 / con1
        stateNew(i,15) = stateOld(i,15) + deqps * on5 / con1
        stateNew(i,16) = stateOld(i,16) + deqps * on6 / con1
C
C Updating equivalent plastic strain rate
        stateNew(i,22) = stateOld(i,22) + one
        stateNew(i,23) = stateOld(i,23) + deqps
        stateNew(i,24) = stateOld(i,24) + dt
        if (strrInc.gt.half .and. strrInc.lt.(one+half)) then
          stateNew(i,21) = deqps/dt
          stateNew(i,27) = deqps
        else
          stateNew(i,21) = (stateNew(i,23)/stateNew(i,24)
     1                   +  stateOld(i,21))/two
          stateNew(i,27) = (stateNew(i,23)/stateNew(i,22)
     1                   +  stateOld(i,27))/two
          if (stateNew(i,22).gt.strrInc-half .and. stateNew(i,22)
     1    .lt.strrInc+half) then
            stateNew(i,22) = zero
            stateNew(i,23) = zero
            stateNew(i,24) = zero
          endif
        endif
        stateNew(i,20) = deqps / dt
C
C Updating temperature
C
C
C Calculating equivalent plastic strain to damage initiation
        epsilondi = (Ddi1*exp(zero-Ddi2*eta) - 
     1                    Ddi3*exp(zero-Ddi4*eta))* thetabar**2
     1           +  Ddi3*exp(zero-Ddi4*eta)
          stateNew(i,33)= epsilondi
c  update damage initiation indicator, damage initition is suppressed when 
c  the stress triaxiality is less than the cut-off value
        if (eta.lt. Cut .and. stateNew(i,31).lt.1.) then
          stateNew(i,31) = stateOld(i,31)
        else
        stateNew(i,31) = stateOld(i,31)+deqps/epsilondi
        endif
        if (stateNew(i,31).ge.1.) then
        stateNew(i,31) = 1.0
        endif
c
C Calculating equivalent plastic strain to cleavage failure
       epsilonc=  ( Cc1*exp(zero-Cc2*eta) - 
     1                    Cc3*exp(zero-Cc4*eta))* thetabar**2
     1           +  Cc3*exp(zero-Cc4*eta)
       stateNew(i,35) = epsilonc
       stateNew(i,36) = stateOld(i,36)+ deqps/epsilonc
        if (stateNew(i,36).ge.1.) then
        stateNew(i,36) = 1.0
        endif
c	   
C Calculating equivalent plastic strain to ductile failure
       epsilondf=  (Ddf1*exp(zero-Ddf2*eta) - 
     1                    Ddf3*exp(zero-Ddf4*eta))* thetabar**2
     1           +  Ddf3*exp(zero-Ddf4*eta)
       stateNew(i,37) = epsilondf   
C
C Ductile damage
        if (stateNew(i,31).ge.1.0 .and. stateNew(i,1).gt.tol) then
          if (stateOld(i,3).le.tol .and. stateOld(i,3).ge.-tol) then
C           Registering the yield stress at the onset of damage for the first time
            stateNew(i,3) = radius
            stateNew(i,38) = stateNew(i,1)
          else
C           Keeping the yield stress at the onset of damage for next step
            stateNew(i,3) = stateOld(i,3)
            stateNew(i,38) = stateOld(i,38)
          endif
C         Updating failure indicator and damage fraction variable
          Dcrit=(stateNew(i,3))*(stateNew(i,37)-stateNew(i,33))/(gf)
          stateNew(i,34)= Dcrit
          if (eta.gt. Cut) then
            if (stateOld(i,32).lt.one) then
             stateNew(i,32) = stateOld(i,32) + deqps/epsilondf
             stateNew(i,2) = stateOld(i,2) + Dcrit*(deqps/epsilondf)
            else
             stateNew(i,2) = one
             stateNew(i,32) = one
             stateNew(i,4) = three
             stateNew(i,28) = zero
            endif
          endif
C       Elements under hydrostatic compression don't experience spherical damage
        if (eta.lt. Cut .and. stateNew(i,32).lt.one) then
          stateNew(i,2) = stateOld(i,2)
          stateNew(i,32) = stateOld(i,32)
        endif
        
        if(stateNew(i,31).ge.1.0 .and. stateNew(i,32).lt.one) then
              stateNew(i,4) = two
        else
              stateNew(i,4) = three
        endif
        else
C
C         In case no damage is reached, yield stress at damage onset
C         and damage variable are kept for next step
          stateNew(i,3) = stateOld(i,3)
          stateNew(i,38) = stateOld(i,38)
          stateNew(i,2) = stateOld(i,2)
          stateNew(i,32) = stateOld(i,32)
        endif
C
C       Previous ductile damage
        if (stateOld(i,4).gt.(one+half)
     1  .and. stateOld(i,4).lt.(two+half)) then
          stateNew(i,4) = two
        endif
        if (stateOld(i,4).gt.(two+half)
     1  .and. stateOld(i,4).lt.(three+half)) then
          stateNew(i,4) = three
        endif
C       Calculating new softening correction terms due to ductile damage
        stateNew(i,28)= one
        if (stateNew(i,32).ge.one) then
          facDctlDmgDev = zero
          facDctlDmgVol = zero
          stateNew(i,32) = one
          stateNew(i,2) = one
          stateNew(i,4) = three
          stateNew(i,28)= zero
          goto 100
        else
          facDctlDmgDev = one - stateNew(i,2)
          facDctlDmgVol = facDctlDmgDev
        endif
C       Elements under hydrostatic compression don't experience spherical damage
        fac1 = facDctlDmgVol * akappa * epsilontrace
        fac2 = facDctlDmgDev * twomu * deqps / con1
        sig1 = facDctlDmgDev * twomu * (stateNew(i,5)
     1       - epsilontrace/three - stateOld(i,11))
        sig2 = facDctlDmgDev * twomu * (stateNew(i,6)
     1       - epsilontrace/three -stateOld(i,12))
        sig3 = facDctlDmgDev * twomu * (stateNew(i,7)
     1       - epsilontrace/three -stateOld(i,13))
        sig4 = facDctlDmgDev * twomu * (stateNew(i,8) -stateOld(i,14))
        sig5 = facDctlDmgDev * twomu * (stateNew(i,9) -stateOld(i,15))
        sig6 = facDctlDmgDev * twomu * (stateNew(i,10)-stateOld(i,16))
C Update stress		
        stressNew(i,1) = sig1 + fac1 - fac2 * on1
        stressNew(i,2) = sig2 + fac1 - fac2 * on2
        stressNew(i,3) = sig3 + fac1 - fac2 * on3
        stressNew(i,4) = sig4        - fac2 * on4
        stressNew(i,5) = sig5        - fac2 * on5
        stressNew(i,6) = sig6        - fac2 * on6

        stateNew(i,39) = (stressNew(i,1) + stressNew(i,2) +
     1                    stressNew(i,3))/three 

C       Calculating invariants of stress tensor
        sig1 = stressNew(i,1)
        sig2 = stressNew(i,2)
        sig3 = stressNew(i,3)
        sig4 = stressNew(i,4)
        sig5 = stressNew(i,5)
        sig6 = stressNew(i,6)
        SI1 =   sig1 + sig2 + sig3
        SI2 =   sig1*sig2-sig4*sig4
     1        + sig1*sig3-sig6*sig6
     1        + sig2*sig3-sig5*sig5
        SI3 =   sig1*(sig2*sig3-sig5*sig5)
     1        - sig4*(sig4*sig3-sig5*sig6)
     1        + sig6*(sig4*sig5-sig2*sig6)
C       Preparing subvalues for calculating the principal stresses values
        cosine2 = (two*SI1*SI1*SI1-three*three*SI1*SI2+three*three
     1            *three*SI3)
     1            /(two*(SI1*SI1-three*SI2)**(threehalves))
C       Assuring that -1<cosine2<1
        if (cosine2.gt.one) cosine2=one
        if (cosine2.lt.(zero-one)) cosine2=zero-one
        alpha2 = acos(cosine2)
C       Calculating the principal stress values
        SP1 = SI1/three + twoThirds*sqrt(SI1*SI1-three*SI2)
     1                   *cos(alpha2/three)
        SP2 = SI1/three + twoThirds*sqrt(SI1*SI1-three*SI2)
     1                   *cos(alpha2/three+twoThirds*pi)
        SP3 = SI1/three + twoThirds*sqrt(SI1*SI1-three*SI2)
     1                   *cos(twoThirds*pi-alpha2/three)
C       Fetching the highest of the principal stress values
        sigmamax = max(abs(SP1),abs(SP2),abs(SP3))
        stateNew(i,26) = sigmamax
c
C       Brittle damage
      if (stateNew(i,36).ge.one .and. stateNew(i,1).gt.tol) then
        if (sigmamax.ge.sigdmg .and. stressNew(i,1).gt.zero) then
          stressNew(i,1) = zero
          stressNew(i,2) = zero
          stressNew(i,3) = zero
          stressNew(i,4) = zero
          stressNew(i,5) = zero
          stressNew(i,6) = zero
          stateNew(i,4)  = one
          stateNew(i,2)  = one
          stateNew(i,28) = zero
           write(*,*) "Coordinate:    ", coordMp(i,1), coordMp(i,2), coordMp(i,3)
           write(*,*) "stateNew(i,1): ", stateNew(i,1), epsilonc
           write(*,*) "eta thetabar: ", eta, thetabar
           write(*,*) "maximum principal stress:", sigmamax
        endif
      endif
C      Diaz2024_Eq29      
       zNT = (10.0d0)**(23.26 - 2.55*EXP( -5.5*stateNew(i,1) )) / 1D9
C
       zdNTdeps = 5876497539218220.0d0*
     1    EXP(-5.87159*EXP(-5.5*stateNew(i,1)) - 5.5*stateNew(i,1))

       zCL_r = zNL*EXP( (zmuR_L - zmu0_L)/(zRgas*zT0) )

       zKT = EXP(zEB/(zRgas*zT0))
C     Diaz2024_Eq21
        zCL_bar = (zNL/zCL_r)*
     1  EXP((tempNew(i)*zmuR_L-zmu0_L + zVH*stateNew(i,39))/(zRgas*zT0))
C     Diaz2024_Eq28
         zDstern = one + 
     1    (zKT*zNT/zNL)/(one + zKT*zCL_bar*zCL_r/zNL)**2 
C     Diaz2024_Eq30
       zThetaT = (zKT*zCL_bar*zCL_r/zNL)/(one + zKT*zCL_bar*zCL_r/zNL)         
C     Diaz2024_Eq48
       IF(totalTime .EQ. 0.0d0) THEN
         DeltaInelasEnergy = 0.0d0
       ELSE
         DeltaInelasEnergy = -zDstern*zCL_bar*zVH/zmuR_L*
     1    (stateNew(i,39) - stateOld(i,39))/density(i) - 
     1    zThetaT*zRgas*zT0/(zCL_r*zmuR_L) * zdNTdeps *
     1    (stateNew(i,1) - stateOld(i,1))/density(i)
       ENDIF
          
        stateNew(i,42) = DeltaInelasEnergy
C     Diaz2024_Eq47
       enerInelasNew(i) = enerInelasOld(i) + DeltaInelasEnergy

       stateNew(i,43) = enerInelasNew(i)
C     Test
       stateNew(i,65) = stateNew(i,39) - stateOld(i,39)
       stateNew(i,64) = stateNew(i,39) - stateNew(i,40)

       stateNew(i,63) = stateNew(i,1) - stateOld(i,1)
       stateNew(i,62) = stateNew(i,1) - stateNew(i,41)
C
  100 continue
C
      return
      end
C
C
      subroutine ahard(sigmayield,hard,eqplas,table,nvalue)
C
      include 'vaba_param.inc'
      dimension table(2,nvalue)
C
C     Set yield stress to last value of table, hardening to zero
      sigmayield=table(1,nvalue)
      hard=0.0
C
C     If more than one entry, search table
      if(nvalue.gt.1) then
        do 10 k1=1,nvalue-1
          eqpl1=table(2,k1+1)
          if(eqplas.lt.eqpl1) then
            eqpl0=table(2,k1)
            if(eqpl1.le.eqpl0) then
              write(6,7)
 7            format(//,30X,'***ERROR - PLASTIC STRAIN MUST BE ',
     1               'ENTERED IN ASCENDING ORDER,')
C
C             Subroutine XIT terminates execution and closes all files
              call XPLB_EXIT
            endif
            deqpl=eqpl1-eqpl0
            sigmayield0=table(1,k1)
            sigmayield1=table(1,k1+1)
            dsigmayield=sigmayield1-sigmayield0
            hard=dsigmayield/deqpl
            sigmayield=sigmayield0+(eqplas-eqpl0)*hard
            goto 20
          endif
 10     continue
 20     continue
        if(eqplas.gt.table(2,nvalue)) then
          hard=(table(1,nvalue)-table(1,nvalue-1))
     1        /(table(2,nvalue)-table(2,nvalue-1))
          sigmayield=table(1,nvalue)+(eqplas-table(2,nvalue))*hard
        endif
      endif
      return
C
C Iteration ends here
      end
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
      subroutine vumatht (
C Read only (unmodifiable) variables -
     *     nblock, nElem, nIntPt, nLayer, nSectPt, 
     *     ntgrad, nstatev, nfieldv, nprops,  
     *     cmname, stepTime, totalTime, dt,  
     *     coordMp, density, props,  
     *     tempOld, fieldOld, stateOld, enerThermOld, 
     *     tempNew, tempgradNew, fieldNew, 
C Write only (modifiable) variables -
     *     stateNew, fluxNew, enerThermNew, dEnerThDTemp, condEff )
C
      include 'vaba_param.inc'
C
      dimension nElem(nblock)
C
      dimension coordMp(nblock,*), density(nblock), props(nprops),
     *    tempOld(nblock), fieldOld(nblock, nfieldv), 
     *    stateOld(nblock, nstatev), enerThermOld(nblock),
     *    tempNew(nblock), tempgradNew(nblock, ntgrad), 
     *    fieldNew(nblock, nfieldv), stateNew(nblock, nstatev),
     *    fluxNew(nblock, ntgrad), enerThermNew(nblock), 
     *    dEnerThDTemp(nblock,2), condEff(nblock)
C
      character*80 cmname
C
      real*8 :: mu0_L,CL0,muR_L,xNL,EB,T0,Rgas,DL,VH 
      real*8 :: CL_r,CL_bar, KT,xNT,Dstern,dEnerTh,xThetaT,xdNTdeps
      real*8 :: xAvg,xBeta,C_trap_mol,C_total,CL_mol,xtime 

      mu0_L = props(1) !28.6D6 mJ/mol
      CL0   = props(2) !2.084D12 atoms/mm3
      muR_L = props(3) !-19.576D6 mJ/mol
      xNL   = props(4) !5.1D20 sites/mm3
      EB    = props(5) !60D6 mJ/mol 
      T0    = props(6) !300K
      Rgas  = props(7) !8314.46262 mJ/(molK)
      VH    = props(8) !2000.0 mm3/mol
      DL    = props(9) !1.27D-2 mm2/s
      xAvg  = props(10) !6.023D23
      xBeta = props(11) !6 
      xtime = props(12)
C      
      do km = 1, nblock
C     Diaz2024_Eq19                              
       CL_r = xNL*EXP( (muR_L - mu0_L)/(Rgas*T0) )
       stateNew(km,50) = CL_r 

       KT = EXP(EB/(Rgas*T0))
       stateNew(km,51) = KT
C     Diaz2024_Eq21
       IF(totalTime .EQ. 0.0d0) THEN
        CL_bar = (xNL/CL_r)*EXP( (muR_L - mu0_L)/(Rgas*T0) )
       ELSE
        CL_bar = (xNL/CL_r)*
     1  EXP((tempNew(km)*muR_L - mu0_L + VH*stateNew(km,39))/(Rgas*T0))
       ENDIF  

        stateNew(km,44) = CL_bar 

C        update heat flux vector (Diaz2024_Eq26) scaled with real testing time 300sec.
         do i = 1, ntgrad
            fluxNew(km, i) = -DL*xtime*CL_bar*tempgradNew(km,i)
         end do
C     Diaz2024_Eq24 (scaled with real testing time 300sec.)
         condEff(km) = DL*xtime*CL_bar 

C         IF(stateNew(km,1) .LE. 0d0) THEN
C            stateNew(km,1) = 0.0d0 
C         ENDIF  
C
C      Diaz2024_Eq29       
         xNT = (10.0d0)**(23.26 - 2.55*EXP( -5.5*stateNew(km,1) )) / 1D9

         stateNew(km,45) = xNT

         xdNTdeps = 5876497539218220.0d0*
     1    EXP(-5.87159*EXP(-5.5*stateNew(km,1)) - 5.5*stateNew(km,1))

         stateNew(km,53) = xdNTdeps
C      Diaz2024_Eq28    
         Dstern = 1.0d0 + 
     1    (KT*xNT/xNL)/(1.0d0 + KT*CL_bar*CL_r/xNL)**2

         stateNew(km,46) = Dstern
C      Diaz2024_Eq25
         if (totalTime .eq. 0.0d0) then
            dEnerTh = Dstern*CL_bar*(tempNew(km))/density(km)
         else
            dEnerTh = Dstern*CL_bar*
     1                (tempNew(km) - tempOld(km))/density(km)
         end if
         enerThermNew(km) = enerThermOld(km) + dEnerTh
         dEnerThDTemp(km,1) = Dstern*CL_bar/density(km)
         dEnerThDTemp(km,2) = 0.0d0    

         stateNew(km,47) = dEnerTh
         stateNew(km,48) = dEnerThDTemp(km,1)     
C      Diaz2024_Eq30         
        xThetaT = (KT*CL_bar*CL_r/xNL)/(1.0d0 + KT*CL_bar*CL_r/xNL)
        stateNew(km,49) = xThetaT

        stateNew(km,54) = Dstern*CL_bar*VH/muR_L

        stateNew(km,55) = xThetaT*Rgas*T0/(CL_r*muR_L) * xdNTdeps
C     Hydrogen concentation in metal lattice Barrera2016_Eq5
        CL_mol = CL_bar*CL_r/xAvg

        stateNew(km,57) = CL_mol
C     Hydrogen concentration in dislocations Barrera2016_Eq16
        C_trap_mol = (xNT/xAvg)*(KT*CL_mol)/
     1               (xBeta*(xNL/xAvg) + KT*CL_mol)

        stateNew(km,58) = C_trap_mol
C      Total hydrogen concentration
        C_total = CL_mol + C_trap_mol

        stateNew(km,56) = C_total

      end do
C
      return
      end