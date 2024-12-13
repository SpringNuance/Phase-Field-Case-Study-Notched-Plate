1

   Abaqus 2023.HF4                                  Date 30-Aug-2024   Time 22:17:33
   For use by AALTO UNIVERSITY under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    1     INCREMENT     1     STEP TIME    0.00    


                        S T E P       1     C O U P L E D - T E M P E R A T U R E - D I S P L A C E M E N T   A N A L Y S I S


                                                                                          

     AUTOMATIC TIME CONTROL WITH -
          A SUGGESTED INITIAL TIME INCREMENT OF                1.000E+03
          AND A TOTAL TIME PERIOD OF                           1.000E+07
          THE MINIMUM TIME INCREMENT ALLOWED IS                 100.    
          THE MAXIMUM TIME INCREMENT ALLOWED IS                1.000E+03

          CREEP AND SWELLING EFFECTS ARE OMITTED IN THIS STEP

     THE SIZE OF THE TIME INCREMENT IS CONTROLLED BY -
          THE TEMPERATURE CHANGE PER INCREMENT NOT EXCEEDING            10.0    
          THE RATE DEPENDENCY WILL BE CONSIDERED IN THIS STEP

          CREEP AND SWELLING EFFECTS WILL BE IGNORED.

 CONVERGENCE TOLERANCE PARAMETERS FOR FORCE    
     CRITERION FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR DISP.    CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INITIAL VALUE OF TIME AVERAGE FORCE                               1.000E-02
     AVERAGE FORCE     IS TIME AVERAGE FORCE    
     ALTERNATE CRIT. FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO FORCE     RELATIVE TO TIME AVRG. FORCE         1.000E-05
     CRITERION FOR RESIDUAL FORCE     WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR DISP.    CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL FORCE     FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO FORCE     REL. TO TIME AVRG. MAX. FORCE        1.000E-05
     CRITERION FOR ZERO DISP.    RELATIVE TO CHARACTERISTIC LENGTH     1.000E-08

 CONVERGENCE TOLERANCE PARAMETERS FOR HEAT FLUX
     CRITERION FOR RESIDUAL HEAT FLUX FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR TEMP.    CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INITIAL VALUE OF TIME AVERAGE HEAT FLUX                           1.000E-02
     AVERAGE HEAT FLUX IS TIME AVERAGE HEAT FLUX
     ALTERNATE CRIT. FOR RESIDUAL HEAT FLUX FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO HEAT FLUX RELATIVE TO TIME AVRG. HEAT FLUX     1.000E-05
     CRITERION FOR RESIDUAL HEAT FLUX WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR TEMP.    CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL HEAT FLUX FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO HEAT FLUX REL. TO TIME AVRG. MAX. HEAT FLUX    1.000E-05

     VOLUMETRIC STRAIN COMPATIBILITY TOLERANCE FOR HYBRID SOLIDS       1.000E-05
     AXIAL STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS             1.000E-05
     TRANS. SHEAR STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS      1.000E-05
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P>P0          5.000E-03
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P=0.0         0.100    
     CONTACT FORCE ERROR TOLERANCE FOR CONVERT SDI=YES                 1.00    
     DISPLACEMENT COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS           1.000E-05
     ROTATION COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS               1.000E-05

 EQUILIBRIUM WILL BE CHECKED FOR SEVERE DISCONTINUITY ITERATIONS

 TIME INCREMENTATION CONTROL PARAMETERS:
     FIRST EQUILIBRIUM ITERATION FOR CONSECUTIVE DIVERGENCE CHECK              4
     EQUILIBRIUM ITERATION AT WHICH LOG. CONVERGENCE RATE CHECK BEGINS         8
     EQUILIBRIUM ITERATION AFTER WHICH ALTERNATE RESIDUAL IS USED              9
     MAXIMUM EQUILIBRIUM ITERATIONS ALLOWED                                   16
     EQUILIBRIUM ITERATION COUNT FOR CUT-BACK IN NEXT INCREMENT               10
     MAXIMUM EQUILIB. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      4
     MAXIMUM ITERATIONS FOR SEVERE DISCONTINUITIES                            50
     MAXIMUM ATTEMPTS ALLOWED IN AN INCREMENT                                  5
     MAXIMUM DISCON. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      50
     MAXIMUM CONTACT AUGMENTATIONS FOR *SURFACE BEHAVIOR,AUGMENTED LAGRANGE   50
     CONSECUTIVE INCS MEETING TIME INTEG. TOL. FOR TIME INC. INCREASE          3
     CUT-BACK FACTOR AFTER DIVERGENCE                                 0.2500    
     CUT-BACK FACTOR FOR TOO SLOW CONVERGENCE                         0.5000    
     CUT-BACK FACTOR AFTER TOO MANY EQUILIBRIUM ITERATIONS            0.7500    
     CUT-BACK FACTOR AFTER TOO MANY SEVERE DISCONTINUITY ITERATIONS   0.2500    
     CUT-BACK FACTOR AFTER PROBLEMS IN ELEMENT ASSEMBLY               0.2500    
     INCREASE FACTOR AFTER TWO INCREMENTS THAT CONVERGE QUICKLY        1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED                       1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DYNAMICS)            1.250    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DIFFUSION)           2.000    
     MINIMUM TIME INCREMENT RATIO FOR EXTRAPOLATION TO OCCUR          0.1000    
     CUT-BACK FACTOR WHEN TIME INTEGRATION ACCURACY TOL. EXCEEDED     0.8500    
     MAX. RATIO OF TIME INTEG. MEASURE TO TOL. FOR TIME INC. INCREASE 0.7500    
     INCREASE FACTOR FOR TIME INCREMENT WHEN POSSIBLE                 0.8000    
     MIN. TIME INC. INCREASE FACTOR IN LINEAR TRANSIENT PROBLEMS      0.9500    
     MAX. RATIO OF TIME INCREMENT TO STABILITY LIMIT                   1.000    
     FRACTION OF STABILITY LIMIT FOR NEW TIME INCREMENT               0.9500    
     TIME INCREMENT INCREASE FACTOR BEFORE A TIME POINT                1.000    
     GLOBAL STABILIZATION CONTROL IS NOT USED

          PRINT OF INCREMENT NUMBER, TIME, ETC., EVERY    1  INCREMENTS

     THE MAXIMUM NUMBER OF INCREMENTS IN THIS STEP IS                 1000000

          LARGE DISPLACEMENT THEORY WILL BE USED

          UNSYMMETRIC MATRIX STORAGE AND SOLUTION WILL BE USED

     LINEAR EXTRAPOLATION WILL BE USED

     CHARACTERISTIC ELEMENT LENGTH     8.915E-06

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     ELEMENT OPERATIONS WILL BE CARRIED OUT IN PARALLEL USING   2 THREADS ON 1 DOMAIN

     PARALLEL CONTACT TRACKING ENABLED


     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+03
 
     INITIAL ELEMENT LOAD BALANCE
            DOMAIN    1     THREAD    NUMBER OF ELEMENTS
                              1            22920
                              2            18272
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.0s
 

 ***WARNING: THE SYSTEM MATRIX HAS 15626 NEGATIVE EIGENVALUES.

     NEGATIVE EIGENVALUES MEAN THAT THE SYSTEM MATRIX IS NOT 
     POSITIVE DEFINITE:
     FOR EXAMPLE, A BIFURCATION (BUCKLING) LOAD MAY HAVE BEEN EXCEEDED.

 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 4883 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 4853 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 4885 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 4886 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18518 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18482 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18163 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18369 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18393 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18610 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18477 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18877 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18874 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18870 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19171 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18971 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19178 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19081 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19076 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19415 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19449 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19437 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19446 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 15941 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 4819 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18167 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 8121 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18270 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18172 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18504 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18862 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18859 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18863 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18620 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18489 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 18633 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19277 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19186 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19349 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19276 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19278 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19173 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19645 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 20056 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19044 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19976 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19143 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19611 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: SOLVER PROBLEM.  ZERO PIVOT WHEN PROCESSING NODE 19193 INSTANCE 
             PLATE-1 D.O.F. 11
  Typically a zero pivot is associated with an overconstraint in the model. 
 However, Abaqus was not able to identify any chains of constraints that might 
 lead to an overconstraint at this node. This suggests that the model might be 
 in fact insufficiently constrained in this particular degree of freedom. 
 Please check the constraints used in this model. 


 ***WARNING: Solver problem. Zero pivot when processing D.O.F. 11 of 49 nodes. 
             The nodes have been identified in node set 
             WarnNodeSolvProbZeroPiv_11_1_1_1_1.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.223E-07   TIME AVG. FORCE       4.223E-07
 LARGEST SCALED RESIDUAL FORCE     -5.183E-17   AT NODE       1540   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.183E-17
 LARGEST INCREMENT OF DISP.         1.064E-09   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.064E-09   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  8.709E+13   TIME AVG. HEAT FLUX   8.709E+13
 LARGEST RESIDUAL HEAT FLUX         5.290E+17   AT NODE       3552   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.098E+11   AT NODE       2335   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.098E+11   AT NODE       2335   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.1s
 

 ***WARNING: THE SYSTEM MATRIX HAS 17102 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.223E-07   TIME AVG. FORCE       4.223E-07
 LARGEST SCALED RESIDUAL FORCE     -6.685E-18   AT NODE       1071   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.685E-18
 LARGEST INCREMENT OF DISP.         1.064E-09   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.799E-21   AT NODE       1488   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  5.184E+16   TIME AVG. HEAT FLUX   5.184E+16
 LARGEST RESIDUAL HEAT FLUX         5.461E+20   AT NODE       2022   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.561E+14   AT NODE      18785   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.561E+14   AT NODE      18785   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.1s
 

 ***WARNING: THE SYSTEM MATRIX HAS 16588 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.223E-07   TIME AVG. FORCE       4.223E-07
 LARGEST SCALED RESIDUAL FORCE     -5.713E-18   AT NODE        715   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.713E-18
 LARGEST INCREMENT OF DISP.         1.064E-09   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -2.746E-24   AT NODE         83   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.270E+18   TIME AVG. HEAT FLUX   2.270E+18
 LARGEST RESIDUAL HEAT FLUX        -2.421E+22   AT NODE       1851   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.019E+17   AT NODE       6717   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.019E+17   AT NODE       6717   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.0s
 

 ***WARNING: THE SYSTEM MATRIX HAS 16454 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.223E-07   TIME AVG. FORCE       4.223E-07
 LARGEST SCALED RESIDUAL FORCE     -6.363E-18   AT NODE       1535   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.363E-18
 LARGEST INCREMENT OF DISP.         1.064E-09   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.043E-24   AT NODE        959   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.328E+15   TIME AVG. HEAT FLUX   1.328E+15
 LARGEST RESIDUAL HEAT FLUX         7.392E+18   AT NODE       2557   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.280E+13   AT NODE       6193   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.019E+17   AT NODE       6717   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.0s
 

 ***WARNING: THE SYSTEM MATRIX HAS 16603 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      4.223E-07   TIME AVG. FORCE       4.223E-07
 LARGEST SCALED RESIDUAL FORCE      5.570E-18   AT NODE       1072   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.570E-18
 LARGEST INCREMENT OF DISP.         1.064E-09   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.350E-24   AT NODE         50   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.766E+15   TIME AVG. HEAT FLUX   2.766E+15
 LARGEST RESIDUAL HEAT FLUX         4.478E+18   AT NODE       4911   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.590E+13   AT NODE       2096   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.300E+13   AT NODE       2096   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.2s
 

 ***WARNING: THE SYSTEM MATRIX HAS 16618 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      4.223E-07   TIME AVG. FORCE       4.223E-07
 LARGEST SCALED RESIDUAL FORCE      6.855E-18   AT NODE       1896   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.855E-18
 LARGEST INCREMENT OF DISP.         1.064E-09   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.760E-24   AT NODE         53   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.079E+15   TIME AVG. HEAT FLUX   1.079E+15
 LARGEST RESIDUAL HEAT FLUX        -2.494E+18   AT NODE       1469   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -6.629E+13   AT NODE       5607   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.629E+13   AT NODE       5607   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.2s
 

 ***WARNING: THE SYSTEM MATRIX HAS 16589 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      4.223E-07   TIME AVG. FORCE       4.223E-07
 LARGEST SCALED RESIDUAL FORCE     -5.929E-18   AT NODE       1419   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.929E-18
 LARGEST INCREMENT OF DISP.         1.064E-09   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -3.547E-24   AT NODE          4   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  6.340E+12   TIME AVG. HEAT FLUX   6.340E+12
 LARGEST RESIDUAL HEAT FLUX         9.251E+15   AT NODE       5191   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.321E+11   AT NODE      20049   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.629E+13   AT NODE       5607   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.1s
 

 ***WARNING: THE SYSTEM MATRIX HAS 17505 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      4.223E-07   TIME AVG. FORCE       4.223E-07
 LARGEST SCALED RESIDUAL FORCE     -6.565E-18   AT NODE       1578   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.565E-18
 LARGEST INCREMENT OF DISP.         1.064E-09   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.139E-24   AT NODE        996   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.892E+13   TIME AVG. HEAT FLUX   2.892E+13
 LARGEST RESIDUAL HEAT FLUX         1.025E+17   AT NODE      15995   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.130E+12   AT NODE       2096   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.126E+12   AT NODE       2096   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.1s
 

 ***WARNING: THE SYSTEM MATRIX HAS 18011 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                      4.223E-07   TIME AVG. FORCE       4.223E-07
 LARGEST SCALED RESIDUAL FORCE     -7.527E-18   AT NODE       1101   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.527E-18
 LARGEST INCREMENT OF DISP.         1.064E-09   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -2.583E-24   AT NODE       1488   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.914E+13   TIME AVG. HEAT FLUX   4.914E+13
 LARGEST RESIDUAL HEAT FLUX        -1.365E+17   AT NODE       3451   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.435E+11   AT NODE       5547   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.138E+12   AT NODE       2096   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT     1 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT   250.    
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.0s
 

 ***WARNING: THE SYSTEM MATRIX HAS 17009 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.056E-07   TIME AVG. FORCE       1.056E-07
 LARGEST SCALED RESIDUAL FORCE     -1.296E-17   AT NODE       1540   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.296E-17
 LARGEST INCREMENT OF DISP.         2.661E-10   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.661E-10   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.793E+13   TIME AVG. HEAT FLUX   1.793E+13
 LARGEST RESIDUAL HEAT FLUX        -9.668E+16   AT NODE       2139   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.000E+12   AT NODE       5582   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.000E+12   AT NODE       5582   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.0s
 

 ***WARNING: THE SYSTEM MATRIX HAS 16854 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.056E-07   TIME AVG. FORCE       1.056E-07
 LARGEST SCALED RESIDUAL FORCE     -1.671E-18   AT NODE       1071   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.671E-18
 LARGEST INCREMENT OF DISP.         2.661E-10   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.450E-21   AT NODE       1488   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  7.111E+12   TIME AVG. HEAT FLUX   7.111E+12
 LARGEST RESIDUAL HEAT FLUX         1.798E+16   AT NODE      18972   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.145E+11   AT NODE       6780   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.000E+12   AT NODE       5582   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.0s
 

 ***WARNING: THE SYSTEM MATRIX HAS 17353 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.056E-07   TIME AVG. FORCE       1.056E-07
 LARGEST SCALED RESIDUAL FORCE     -1.428E-18   AT NODE        715   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.428E-18
 LARGEST INCREMENT OF DISP.         2.661E-10   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.866E-25   AT NODE         83   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.533E+15   TIME AVG. HEAT FLUX   1.533E+15
 LARGEST RESIDUAL HEAT FLUX        -1.073E+19   AT NODE       2022   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.017E+13   AT NODE      24289   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.012E+13   AT NODE      24289   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.0s
 

 ***WARNING: THE SYSTEM MATRIX HAS 16680 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.056E-07   TIME AVG. FORCE       1.056E-07
 LARGEST SCALED RESIDUAL FORCE     -1.591E-18   AT NODE       1535   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.591E-18
 LARGEST INCREMENT OF DISP.         2.661E-10   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.607E-25   AT NODE        959   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.247E+16   TIME AVG. HEAT FLUX   1.247E+16
 LARGEST RESIDUAL HEAT FLUX        -1.212E+20   AT NODE       2096   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.797E+13   AT NODE      18307   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.786E+13   AT NODE      18307   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT     1 STARTS. ATTEMPT NUMBER  3, TIME INCREMENT   100.    
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  2.9s
 

 ***WARNING: THE SYSTEM MATRIX HAS 17189 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.223E-08   TIME AVG. FORCE       1.000E-02
 LARGEST SCALED RESIDUAL FORCE     -3.258E-18   AT NODE      21957   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.258E-18
 LARGEST INCREMENT OF DISP.         1.064E-10   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.064E-10   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE HEAT FLUX                  5.298E+12   TIME AVG. HEAT FLUX   5.298E+12
 LARGEST RESIDUAL HEAT FLUX        -2.623E+16   AT NODE       2647   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -6.476E+10   AT NODE       5748   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.476E+10   AT NODE       5748   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  2.9s
 

 ***WARNING: THE SYSTEM MATRIX HAS 17529 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.223E-08   TIME AVG. FORCE       1.000E-02
 LARGEST SCALED RESIDUAL FORCE      8.760E-19   AT NODE        715   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      8.760E-19
 LARGEST INCREMENT OF DISP.         1.064E-10   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.648E-22   AT NODE      15985   DOF  3
   INSTANCE: PLATE-1                                                                         
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE HEAT FLUX                  6.581E+12   TIME AVG. HEAT FLUX   6.581E+12
 LARGEST RESIDUAL HEAT FLUX         4.459E+16   AT NODE       6486   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.070E+11   AT NODE       4758   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.070E+11   AT NODE       4758   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.0s
 

 ***WARNING: THE SYSTEM MATRIX HAS 17601 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.223E-08   TIME AVG. FORCE       1.000E-02
 LARGEST SCALED RESIDUAL FORCE      6.234E-19   AT NODE       1062   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.234E-19
 LARGEST INCREMENT OF DISP.         1.064E-10   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.257E-25   AT NODE      16039   DOF  3
   INSTANCE: PLATE-1                                                                         
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE HEAT FLUX                  6.408E+13   TIME AVG. HEAT FLUX   6.408E+13
 LARGEST RESIDUAL HEAT FLUX        -6.290E+17   AT NODE       2351   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.277E+11   AT NODE       1843   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.070E+11   AT NODE       4758   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  131150
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      2.469e+11
        SOLVER ELAPSED TIME:  3.0s
 

 ***WARNING: THE SYSTEM MATRIX HAS 17050 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.223E-08   TIME AVG. FORCE       1.000E-02
 LARGEST SCALED RESIDUAL FORCE      4.368E-19   AT NODE       1575   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.368E-19
 LARGEST INCREMENT OF DISP.         1.064E-10   AT NODE        960   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.417E-25   AT NODE         84   DOF  3
   INSTANCE: PLATE-1                                                                         
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE HEAT FLUX                  2.915E+16   TIME AVG. HEAT FLUX   2.915E+16
 LARGEST RESIDUAL HEAT FLUX        -4.285E+20   AT NODE       5319   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.357E+15   AT NODE       5548   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.357E+15   AT NODE       5548   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.
 

 ***ERROR: TIME INCREMENT REQUIRED IS LESS THAN THE MINIMUM SPECIFIED



     ANALYSIS SUMMARY:
     TOTAL OF          1  INCREMENTS
                       2  CUTBACKS IN AUTOMATIC INCREMENTATION
                      17  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                      17  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                      17  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       3  WARNING MESSAGES DURING USER INPUT PROCESSING
                    1166  WARNING MESSAGES DURING ANALYSIS
                    1145  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                      17  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       1  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =      56.    
       SYSTEM TIME (SEC)    =      2.5    
       TOTAL CPU TIME (SEC) =      58.    
       WALLCLOCK TIME (SEC) =           83
       MEMORY PEAK (GB)     =            0
