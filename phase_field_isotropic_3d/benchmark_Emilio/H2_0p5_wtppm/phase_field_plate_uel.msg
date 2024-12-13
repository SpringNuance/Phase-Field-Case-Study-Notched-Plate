1

   Abaqus 2023.HF4                                  Date 12-Aug-2024   Time 21:12:25
   For use by AALTO UNIVERSITY under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    1     INCREMENT     1     STEP TIME    0.00    


                        S T E P       1     C O U P L E D - T E M P E R A T U R E - D I S P L A C E M E N T   A N A L Y S I S


                                                                                          

          FIXED TIME INCREMENTS
          TIME INCREMENT IS                                    1.000E+04
          TIME PERIOD IS                                       1.000E+07

          CREEP AND SWELLING EFFECTS ARE OMITTED IN THIS STEP
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

 ITERATION CONTROL PARAMETERS:
     FIRST EQUILIBRIUM ITERATION FOR CONSECUTIVE DIVERGENCE CHECK              4
     EQUILIBRIUM ITERATION AFTER WHICH ALTERNATE RESIDUAL IS USED              9
     MAXIMUM EQUILIBRIUM ITERATIONS ALLOWED                                   16
     MAXIMUM ITERATIONS FOR SEVERE DISCONTINUITIES                            12
     GLOBAL STABILIZATION CONTROL IS NOT USED

          PRINT OF INCREMENT NUMBER, TIME, ETC., EVERY    1  INCREMENTS

     THE MAXIMUM NUMBER OF INCREMENTS IN THIS STEP IS                 1000000

          UNSYMMETRIC MATRIX STORAGE AND SOLUTION WILL BE USED

     LINEAR EXTRAPOLATION WILL BE USED

     CHARACTERISTIC ELEMENT LENGTH     1.066E-02

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     ELEMENT OPERATIONS WILL BE CARRIED OUT IN PARALLEL USING   2 THREADS ON 1 DOMAIN

     PARALLEL CONTACT TRACKING ENABLED


     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
     INITIAL ELEMENT LOAD BALANCE
            DOMAIN    1     THREAD    NUMBER OF ELEMENTS
                              1            5352
                              2            5456
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  197ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.271E-03   TIME AVG. FORCE       1.271E-03
 LARGEST SCALED RESIDUAL FORCE     -1.395E-14   AT NODE       6388   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.395E-14
 LARGEST INCREMENT OF DISP.         1.065E-05   AT NODE         38   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.065E-05   AT NODE         38   DOF  2
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.209E-16   TIME AVG. HEAT FLUX   1.000E-02
 LARGEST RESIDUAL HEAT FLUX         1.082E-15   AT NODE       3583   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         4.881E-15   AT NODE       1335   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.881E-15   AT NODE       1335   DOF 11
   INSTANCE: PLATE-1                                                                         
          THERE IS ZERO HEAT FLUX EVERYWHERE
 

 ***WARNING: THERE IS ZERO HEAT FLUX EVERYWHERE IN THE MODEL BASED ON THE 
             DEFAULT CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE HEAT FLUX 
             DURING THE CURRENT ITERATION TO VERIFY THAT THE HEAT FLUX IS SMALL 
             ENOUGH TO BE TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION 
             CONTROLS TO RESET THE CRITERION FOR ZERO HEAT FLUX.

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  1.000E-03
 STEP TIME COMPLETED       1.000E+04,  TOTAL TIME COMPLETED        1.000E+04


  INCREMENT     2 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  17ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.389E-03   TIME AVG. FORCE       3.830E-03
 LARGEST SCALED RESIDUAL FORCE      6.172E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.172E-11
 LARGEST INCREMENT OF DISP.         1.065E-05   AT NODE         38   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.159E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.505E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.980E-06   TIME AVG. HEAT FLUX   1.980E-06
 LARGEST RESIDUAL HEAT FLUX         6.475E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.569E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.569E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  190ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.389E-03   TIME AVG. FORCE       3.830E-03
 LARGEST SCALED RESIDUAL FORCE      1.999E-13   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.999E-13
 LARGEST INCREMENT OF DISP.         1.065E-05   AT NODE         38   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.639E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.979E-06   TIME AVG. HEAT FLUX   1.979E-06
 LARGEST RESIDUAL HEAT FLUX        -1.448E-15   AT NODE      11986   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.599E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.018E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  2.000E-03
 STEP TIME COMPLETED       2.000E+04,  TOTAL TIME COMPLETED        2.000E+04


  INCREMENT     3 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  191ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.583E-03   TIME AVG. FORCE       5.748E-03
 LARGEST SCALED RESIDUAL FORCE      1.595E-12   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.595E-12
 LARGEST INCREMENT OF DISP.         1.250E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.334E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.979E-06   TIME AVG. HEAT FLUX   2.979E-06
 LARGEST RESIDUAL HEAT FLUX         1.235E-15   AT NODE       2347   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.660E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.076E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  3.000E-03
 STEP TIME COMPLETED       3.000E+04,  TOTAL TIME COMPLETED        3.000E+04


  INCREMENT     4 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  16ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.278E-02   TIME AVG. FORCE       7.505E-03
 LARGEST SCALED RESIDUAL FORCE      1.275E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.275E-10
 LARGEST INCREMENT OF DISP.         2.084E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.345E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.649E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.973E-06   TIME AVG. HEAT FLUX   3.977E-06
 LARGEST RESIDUAL HEAT FLUX         5.604E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.721E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.090E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  180ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.278E-02   TIME AVG. FORCE       7.505E-03
 LARGEST SCALED RESIDUAL FORCE      8.414E-15   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      8.414E-15
 LARGEST INCREMENT OF DISP.         2.085E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.279E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  5.974E-06   TIME AVG. HEAT FLUX   3.977E-06
 LARGEST RESIDUAL HEAT FLUX         1.381E-15   AT NODE       3878   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.721E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.127E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  4.000E-03
 STEP TIME COMPLETED       4.000E+04,  TOTAL TIME COMPLETED        4.000E+04


  INCREMENT     5 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  193ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.597E-02   TIME AVG. FORCE       9.198E-03
 LARGEST SCALED RESIDUAL FORCE      6.286E-12   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.286E-12
 LARGEST INCREMENT OF DISP.         2.921E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.366E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  7.972E-06   TIME AVG. HEAT FLUX   4.976E-06
 LARGEST RESIDUAL HEAT FLUX        -1.540E-15   AT NODE       3805   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.782E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.161E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  5.000E-03
 STEP TIME COMPLETED       5.000E+04,  TOTAL TIME COMPLETED        5.000E+04


  INCREMENT     6 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  15ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.916E-02   TIME AVG. FORCE       1.086E-02
 LARGEST SCALED RESIDUAL FORCE      1.342E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.342E-10
 LARGEST INCREMENT OF DISP.         3.759E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.377E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.921E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  9.974E-06   TIME AVG. HEAT FLUX   5.976E-06
 LARGEST RESIDUAL HEAT FLUX         5.670E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.844E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.187E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  171ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.916E-02   TIME AVG. FORCE       1.086E-02
 LARGEST SCALED RESIDUAL FORCE      2.311E-14   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.311E-14
 LARGEST INCREMENT OF DISP.         3.759E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.336E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  9.974E-06   TIME AVG. HEAT FLUX   5.976E-06
 LARGEST RESIDUAL HEAT FLUX        -1.513E-15   AT NODE       5434   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.844E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.157E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  6.000E-03
 STEP TIME COMPLETED       6.000E+04,  TOTAL TIME COMPLETED        6.000E+04


  INCREMENT     7 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  197ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.236E-02   TIME AVG. FORCE       1.250E-02
 LARGEST SCALED RESIDUAL FORCE      1.399E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.399E-11
 LARGEST INCREMENT OF DISP.         4.599E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.398E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.198E-05   TIME AVG. HEAT FLUX   6.976E-06
 LARGEST RESIDUAL HEAT FLUX        -1.417E-15   AT NODE       4717   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.907E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.257E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  7.000E-03
 STEP TIME COMPLETED       7.000E+04,  TOTAL TIME COMPLETED        7.000E+04


  INCREMENT     8 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  17ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.555E-02   TIME AVG. FORCE       1.413E-02
 LARGEST SCALED RESIDUAL FORCE      1.538E-10   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.538E-10
 LARGEST INCREMENT OF DISP.         5.440E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.410E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.899E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.399E-05   TIME AVG. HEAT FLUX   7.978E-06
 LARGEST RESIDUAL HEAT FLUX         5.732E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.970E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.283E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  211ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.555E-02   TIME AVG. FORCE       1.413E-02
 LARGEST SCALED RESIDUAL FORCE      4.473E-14   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.473E-14
 LARGEST INCREMENT OF DISP.         5.441E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.458E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.399E-05   TIME AVG. HEAT FLUX   7.978E-06
 LARGEST RESIDUAL HEAT FLUX        -1.302E-15   AT NODE       4503   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.970E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.186E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  8.000E-03
 STEP TIME COMPLETED       8.000E+04,  TOTAL TIME COMPLETED        8.000E+04


  INCREMENT     9 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  212ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.874E-02   TIME AVG. FORCE       1.576E-02
 LARGEST SCALED RESIDUAL FORCE      2.461E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.461E-11
 LARGEST INCREMENT OF DISP.         6.284E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.430E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.600E-05   TIME AVG. HEAT FLUX   8.980E-06
 LARGEST RESIDUAL HEAT FLUX        -1.321E-15   AT NODE       2776   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.034E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.352E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.000E-03
 STEP TIME COMPLETED       9.000E+04,  TOTAL TIME COMPLETED        9.000E+04


  INCREMENT    10 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  14ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.194E-02   TIME AVG. FORCE       1.737E-02
 LARGEST SCALED RESIDUAL FORCE      2.837E-10   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.837E-10
 LARGEST INCREMENT OF DISP.         7.128E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.443E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    8.076E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.801E-05   TIME AVG. HEAT FLUX   9.984E-06
 LARGEST RESIDUAL HEAT FLUX         5.793E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.097E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.378E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  187ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.194E-02   TIME AVG. FORCE       1.737E-02
 LARGEST SCALED RESIDUAL FORCE      7.298E-14   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      7.298E-14
 LARGEST INCREMENT OF DISP.         7.129E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.644E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.802E-05   TIME AVG. HEAT FLUX   9.984E-06
 LARGEST RESIDUAL HEAT FLUX         1.507E-15   AT NODE       4061   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.098E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.214E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  1.000E-02
 STEP TIME COMPLETED       1.000E+05,  TOTAL TIME COMPLETED        1.000E+05


  INCREMENT    11 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  376ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.513E-02   TIME AVG. FORCE       1.899E-02
 LARGEST SCALED RESIDUAL FORCE      3.802E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.802E-11
 LARGEST INCREMENT OF DISP.         7.975E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.463E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.003E-05   TIME AVG. HEAT FLUX   1.099E-05
 LARGEST RESIDUAL HEAT FLUX         1.608E-15   AT NODE        697   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.162E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.447E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  1.100E-02
 STEP TIME COMPLETED       1.100E+05,  TOTAL TIME COMPLETED        1.100E+05


  INCREMENT    12 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.832E-02   TIME AVG. FORCE       2.060E-02
 LARGEST SCALED RESIDUAL FORCE      4.531E-10   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.531E-10
 LARGEST INCREMENT OF DISP.         8.823E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.478E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.027E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.206E-05   TIME AVG. HEAT FLUX   1.200E-05
 LARGEST RESIDUAL HEAT FLUX         5.853E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.227E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.472E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.832E-02   TIME AVG. FORCE       2.060E-02
 LARGEST SCALED RESIDUAL FORCE      1.076E-13   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.076E-13
 LARGEST INCREMENT OF DISP.         8.823E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.893E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.206E-05   TIME AVG. HEAT FLUX   1.200E-05
 LARGEST RESIDUAL HEAT FLUX         1.410E-15   AT NODE       4314   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.227E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.241E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  1.200E-02
 STEP TIME COMPLETED       1.200E+05,  TOTAL TIME COMPLETED        1.200E+05


  INCREMENT    13 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.151E-02   TIME AVG. FORCE       2.221E-02
 LARGEST SCALED RESIDUAL FORCE      5.411E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.411E-11
 LARGEST INCREMENT OF DISP.         9.673E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.497E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.408E-05   TIME AVG. HEAT FLUX   1.300E-05
 LARGEST RESIDUAL HEAT FLUX        -1.372E-15   AT NODE       2237   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.292E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.542E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  1.300E-02
 STEP TIME COMPLETED       1.300E+05,  TOTAL TIME COMPLETED        1.300E+05


  INCREMENT    14 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.470E-02   TIME AVG. FORCE       2.381E-02
 LARGEST SCALED RESIDUAL FORCE      6.620E-10   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.620E-10
 LARGEST INCREMENT OF DISP.         1.052E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.513E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.248E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.611E-05   TIME AVG. HEAT FLUX   1.401E-05
 LARGEST RESIDUAL HEAT FLUX         5.910E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.358E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.566E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.470E-02   TIME AVG. FORCE       2.381E-02
 LARGEST SCALED RESIDUAL FORCE      1.481E-13   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.481E-13
 LARGEST INCREMENT OF DISP.         1.052E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.201E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.611E-05   TIME AVG. HEAT FLUX   1.401E-05
 LARGEST RESIDUAL HEAT FLUX        -1.366E-15   AT NODE       9876   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.358E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.268E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  1.400E-02
 STEP TIME COMPLETED       1.400E+05,  TOTAL TIME COMPLETED        1.400E+05


  INCREMENT    15 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.789E-02   TIME AVG. FORCE       2.542E-02
 LARGEST SCALED RESIDUAL FORCE      7.278E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      7.278E-11
 LARGEST INCREMENT OF DISP.         1.138E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.531E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.815E-05   TIME AVG. HEAT FLUX   1.502E-05
 LARGEST RESIDUAL HEAT FLUX        -1.189E-15   AT NODE        826   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.425E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.636E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  1.500E-02
 STEP TIME COMPLETED       1.500E+05,  TOTAL TIME COMPLETED        1.500E+05


  INCREMENT    16 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.108E-02   TIME AVG. FORCE       2.702E-02
 LARGEST SCALED RESIDUAL FORCE      9.105E-10   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      9.105E-10
 LARGEST INCREMENT OF DISP.         1.223E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.548E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.471E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.018E-05   TIME AVG. HEAT FLUX   1.603E-05
 LARGEST RESIDUAL HEAT FLUX         5.966E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.491E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.660E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.108E-02   TIME AVG. FORCE       2.702E-02
 LARGEST SCALED RESIDUAL FORCE      1.944E-13   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.944E-13
 LARGEST INCREMENT OF DISP.         1.223E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.568E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.018E-05   TIME AVG. HEAT FLUX   1.603E-05
 LARGEST RESIDUAL HEAT FLUX        -1.349E-15   AT NODE       2535   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.492E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.295E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  1.600E-02
 STEP TIME COMPLETED       1.600E+05,  TOTAL TIME COMPLETED        1.600E+05


  INCREMENT    17 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.427E-02   TIME AVG. FORCE       2.863E-02
 LARGEST SCALED RESIDUAL FORCE      9.389E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      9.389E-11
 LARGEST INCREMENT OF DISP.         1.309E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.566E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.222E-05   TIME AVG. HEAT FLUX   1.704E-05
 LARGEST RESIDUAL HEAT FLUX        -1.267E-15   AT NODE      14649   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.559E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.729E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  1.700E-02
 STEP TIME COMPLETED       1.700E+05,  TOTAL TIME COMPLETED        1.700E+05


  INCREMENT    18 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.745E-02   TIME AVG. FORCE       3.023E-02
 LARGEST SCALED RESIDUAL FORCE      1.199E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.199E-09
 LARGEST INCREMENT OF DISP.         1.395E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.585E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.697E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.426E-05   TIME AVG. HEAT FLUX   1.806E-05
 LARGEST RESIDUAL HEAT FLUX         6.019E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.626E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.752E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.745E-02   TIME AVG. FORCE       3.023E-02
 LARGEST SCALED RESIDUAL FORCE      2.473E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.473E-13
 LARGEST INCREMENT OF DISP.         1.395E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.991E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.426E-05   TIME AVG. HEAT FLUX   1.806E-05
 LARGEST RESIDUAL HEAT FLUX         1.325E-15   AT NODE        262   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.627E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.320E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  1.800E-02
 STEP TIME COMPLETED       1.800E+05,  TOTAL TIME COMPLETED        1.800E+05


  INCREMENT    19 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.064E-02   TIME AVG. FORCE       3.183E-02
 LARGEST SCALED RESIDUAL FORCE      1.173E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.173E-10
 LARGEST INCREMENT OF DISP.         1.481E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.601E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.631E-05   TIME AVG. HEAT FLUX   1.907E-05
 LARGEST RESIDUAL HEAT FLUX        -1.413E-15   AT NODE       1903   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.695E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.822E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  1.900E-02
 STEP TIME COMPLETED       1.900E+05,  TOTAL TIME COMPLETED        1.900E+05


  INCREMENT    20 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  24ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.382E-02   TIME AVG. FORCE       3.343E-02
 LARGEST SCALED RESIDUAL FORCE      1.527E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.527E-09
 LARGEST INCREMENT OF DISP.         1.567E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.622E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.924E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.835E-05   TIME AVG. HEAT FLUX   2.009E-05
 LARGEST RESIDUAL HEAT FLUX         6.071E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.763E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.844E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.382E-02   TIME AVG. FORCE       3.343E-02
 LARGEST SCALED RESIDUAL FORCE      3.069E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.069E-13
 LARGEST INCREMENT OF DISP.         1.567E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.467E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.836E-05   TIME AVG. HEAT FLUX   2.009E-05
 LARGEST RESIDUAL HEAT FLUX         1.599E-15   AT NODE       3667   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.763E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.345E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  2.000E-02
 STEP TIME COMPLETED       2.000E+05,  TOTAL TIME COMPLETED        2.000E+05


  INCREMENT    21 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.700E-02   TIME AVG. FORCE       3.503E-02
 LARGEST SCALED RESIDUAL FORCE      1.430E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.430E-10
 LARGEST INCREMENT OF DISP.         1.654E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.637E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.041E-05   TIME AVG. HEAT FLUX   2.110E-05
 LARGEST RESIDUAL HEAT FLUX        -1.397E-15   AT NODE       4600   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.833E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.913E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  2.100E-02
 STEP TIME COMPLETED       2.100E+05,  TOTAL TIME COMPLETED        2.100E+05


  INCREMENT    22 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.018E-02   TIME AVG. FORCE       3.662E-02
 LARGEST SCALED RESIDUAL FORCE      1.895E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.895E-09
 LARGEST INCREMENT OF DISP.         1.740E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.660E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.153E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.247E-05   TIME AVG. HEAT FLUX   2.212E-05
 LARGEST RESIDUAL HEAT FLUX         6.121E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.902E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.935E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.018E-02   TIME AVG. FORCE       3.662E-02
 LARGEST SCALED RESIDUAL FORCE      3.723E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.723E-13
 LARGEST INCREMENT OF DISP.         1.740E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.995E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.247E-05   TIME AVG. HEAT FLUX   2.212E-05
 LARGEST RESIDUAL HEAT FLUX         1.411E-15   AT NODE       4771   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.902E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.369E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  2.200E-02
 STEP TIME COMPLETED       2.200E+05,  TOTAL TIME COMPLETED        2.200E+05


  INCREMENT    23 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  353ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.337E-02   TIME AVG. FORCE       3.822E-02
 LARGEST SCALED RESIDUAL FORCE      1.707E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.707E-10
 LARGEST INCREMENT OF DISP.         1.827E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.673E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.453E-05   TIME AVG. HEAT FLUX   2.314E-05
 LARGEST RESIDUAL HEAT FLUX         1.421E-15   AT NODE       3131   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.972E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.004E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  2.300E-02
 STEP TIME COMPLETED       2.300E+05,  TOTAL TIME COMPLETED        2.300E+05


  INCREMENT    24 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.655E-02   TIME AVG. FORCE       3.982E-02
 LARGEST SCALED RESIDUAL FORCE      2.302E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.302E-09
 LARGEST INCREMENT OF DISP.         1.914E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.699E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.384E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.659E-05   TIME AVG. HEAT FLUX   2.416E-05
 LARGEST RESIDUAL HEAT FLUX         6.168E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.042E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.025E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.655E-02   TIME AVG. FORCE       3.982E-02
 LARGEST SCALED RESIDUAL FORCE      4.432E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.432E-13
 LARGEST INCREMENT OF DISP.         1.914E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.573E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.659E-05   TIME AVG. HEAT FLUX   2.416E-05
 LARGEST RESIDUAL HEAT FLUX        -1.274E-15   AT NODE       1286   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.043E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.392E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  2.400E-02
 STEP TIME COMPLETED       2.400E+05,  TOTAL TIME COMPLETED        2.400E+05


  INCREMENT    25 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.972E-02   TIME AVG. FORCE       4.142E-02
 LARGEST SCALED RESIDUAL FORCE      2.004E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.004E-10
 LARGEST INCREMENT OF DISP.         2.001E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.710E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.865E-05   TIME AVG. HEAT FLUX   2.518E-05
 LARGEST RESIDUAL HEAT FLUX         1.368E-15   AT NODE        504   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.114E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.093E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  2.500E-02
 STEP TIME COMPLETED       2.500E+05,  TOTAL TIME COMPLETED        2.500E+05


  INCREMENT    26 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  24ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.290E-02   TIME AVG. FORCE       4.301E-02
 LARGEST SCALED RESIDUAL FORCE      2.751E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.751E-09
 LARGEST INCREMENT OF DISP.         2.089E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.738E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.618E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.073E-05   TIME AVG. HEAT FLUX   2.620E-05
 LARGEST RESIDUAL HEAT FLUX         6.213E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.185E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.114E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      8.290E-02   TIME AVG. FORCE       4.301E-02
 LARGEST SCALED RESIDUAL FORCE      5.195E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.195E-13
 LARGEST INCREMENT OF DISP.         2.089E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.197E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  5.073E-05   TIME AVG. HEAT FLUX   2.620E-05
 LARGEST RESIDUAL HEAT FLUX         1.506E-15   AT NODE       1724   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.185E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.415E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  2.600E-02
 STEP TIME COMPLETED       2.600E+05,  TOTAL TIME COMPLETED        2.600E+05


  INCREMENT    27 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.607E-02   TIME AVG. FORCE       4.461E-02
 LARGEST SCALED RESIDUAL FORCE      2.318E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.318E-10
 LARGEST INCREMENT OF DISP.         2.176E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.748E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  5.281E-05   TIME AVG. HEAT FLUX   2.722E-05
 LARGEST RESIDUAL HEAT FLUX         1.333E-15   AT NODE       4501   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.257E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.182E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  2.700E-02
 STEP TIME COMPLETED       2.700E+05,  TOTAL TIME COMPLETED        2.700E+05


  INCREMENT    28 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.924E-02   TIME AVG. FORCE       4.620E-02
 LARGEST SCALED RESIDUAL FORCE      3.239E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.239E-09
 LARGEST INCREMENT OF DISP.         2.264E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.779E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.854E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.489E-05   TIME AVG. HEAT FLUX   2.825E-05
 LARGEST RESIDUAL HEAT FLUX         6.255E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.329E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.201E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      8.924E-02   TIME AVG. FORCE       4.620E-02
 LARGEST SCALED RESIDUAL FORCE      6.008E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.008E-13
 LARGEST INCREMENT OF DISP.         2.264E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.865E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  5.489E-05   TIME AVG. HEAT FLUX   2.825E-05
 LARGEST RESIDUAL HEAT FLUX        -1.280E-15   AT NODE       3162   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.329E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.436E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  2.800E-02
 STEP TIME COMPLETED       2.800E+05,  TOTAL TIME COMPLETED        2.800E+05


  INCREMENT    29 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.242E-02   TIME AVG. FORCE       4.779E-02
 LARGEST SCALED RESIDUAL FORCE      2.650E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.650E-10
 LARGEST INCREMENT OF DISP.         2.352E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.786E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  5.697E-05   TIME AVG. HEAT FLUX   2.927E-05
 LARGEST RESIDUAL HEAT FLUX        -1.568E-15   AT NODE       5487   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.402E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.268E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  2.900E-02
 STEP TIME COMPLETED       2.900E+05,  TOTAL TIME COMPLETED        2.900E+05


  INCREMENT    30 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  24ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.559E-02   TIME AVG. FORCE       4.939E-02
 LARGEST SCALED RESIDUAL FORCE      3.768E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.768E-09
 LARGEST INCREMENT OF DISP.         2.440E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.820E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.093E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.905E-05   TIME AVG. HEAT FLUX   3.030E-05
 LARGEST RESIDUAL HEAT FLUX         6.294E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.475E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.287E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      9.559E-02   TIME AVG. FORCE       4.939E-02
 LARGEST SCALED RESIDUAL FORCE      6.869E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.869E-13
 LARGEST INCREMENT OF DISP.         2.440E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.057E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  5.905E-05   TIME AVG. HEAT FLUX   3.030E-05
 LARGEST RESIDUAL HEAT FLUX         1.474E-15   AT NODE       1775   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.475E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.456E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  3.000E-02
 STEP TIME COMPLETED       3.000E+05,  TOTAL TIME COMPLETED        3.000E+05


  INCREMENT    31 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  368ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.876E-02   TIME AVG. FORCE       5.098E-02
 LARGEST SCALED RESIDUAL FORCE      2.997E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.997E-10
 LARGEST INCREMENT OF DISP.         2.528E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.825E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  6.115E-05   TIME AVG. HEAT FLUX   3.133E-05
 LARGEST RESIDUAL HEAT FLUX         1.296E-15   AT NODE       3681   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.548E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.353E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  3.100E-02
 STEP TIME COMPLETED       3.100E+05,  TOTAL TIME COMPLETED        3.100E+05


  INCREMENT    32 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  24ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.102       TIME AVG. FORCE       5.257E-02
 LARGEST SCALED RESIDUAL FORCE      4.339E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.339E-09
 LARGEST INCREMENT OF DISP.         2.617E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.862E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.334E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.325E-05   TIME AVG. HEAT FLUX   3.236E-05
 LARGEST RESIDUAL HEAT FLUX         6.331E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.622E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.372E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.102       TIME AVG. FORCE       5.257E-02
 LARGEST SCALED RESIDUAL FORCE      7.775E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      7.775E-13
 LARGEST INCREMENT OF DISP.         2.617E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.132E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  6.325E-05   TIME AVG. HEAT FLUX   3.236E-05
 LARGEST RESIDUAL HEAT FLUX        -1.187E-15   AT NODE       3223   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.622E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.476E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  3.200E-02
 STEP TIME COMPLETED       3.200E+05,  TOTAL TIME COMPLETED        3.200E+05


  INCREMENT    33 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.105       TIME AVG. FORCE       5.416E-02
 LARGEST SCALED RESIDUAL FORCE      3.357E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.357E-10
 LARGEST INCREMENT OF DISP.         2.706E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.865E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  6.534E-05   TIME AVG. HEAT FLUX   3.339E-05
 LARGEST RESIDUAL HEAT FLUX        -1.840E-15   AT NODE        682   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.697E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.437E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  3.300E-02
 STEP TIME COMPLETED       3.300E+05,  TOTAL TIME COMPLETED        3.300E+05


  INCREMENT    34 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.108       TIME AVG. FORCE       5.575E-02
 LARGEST SCALED RESIDUAL FORCE      4.950E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.950E-09
 LARGEST INCREMENT OF DISP.         2.795E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.904E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.577E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.744E-05   TIME AVG. HEAT FLUX   3.442E-05
 LARGEST RESIDUAL HEAT FLUX         6.365E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.771E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.454E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.108       TIME AVG. FORCE       5.575E-02
 LARGEST SCALED RESIDUAL FORCE      8.723E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      8.723E-13
 LARGEST INCREMENT OF DISP.         2.795E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.211E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  6.745E-05   TIME AVG. HEAT FLUX   3.442E-05
 LARGEST RESIDUAL HEAT FLUX        -1.503E-15   AT NODE       2181   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.772E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.494E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  3.400E-02
 STEP TIME COMPLETED       3.400E+05,  TOTAL TIME COMPLETED        3.400E+05


  INCREMENT    35 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.111       TIME AVG. FORCE       5.734E-02
 LARGEST SCALED RESIDUAL FORCE      3.731E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.731E-10
 LARGEST INCREMENT OF DISP.         2.884E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.905E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  6.956E-05   TIME AVG. HEAT FLUX   3.546E-05
 LARGEST RESIDUAL HEAT FLUX         1.382E-15   AT NODE       2938   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.847E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.518E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  3.500E-02
 STEP TIME COMPLETED       3.500E+05,  TOTAL TIME COMPLETED        3.500E+05


  INCREMENT    36 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.115       TIME AVG. FORCE       5.893E-02
 LARGEST SCALED RESIDUAL FORCE      5.603E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.603E-09
 LARGEST INCREMENT OF DISP.         2.974E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.948E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.823E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.166E-05   TIME AVG. HEAT FLUX   3.649E-05
 LARGEST RESIDUAL HEAT FLUX         6.395E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.922E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.534E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  371ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.115       TIME AVG. FORCE       5.893E-02
 LARGEST SCALED RESIDUAL FORCE      9.708E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      9.708E-13
 LARGEST INCREMENT OF DISP.         2.974E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.292E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  7.166E-05   TIME AVG. HEAT FLUX   3.649E-05
 LARGEST RESIDUAL HEAT FLUX        -1.486E-15   AT NODE       4090   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.922E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.511E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  3.600E-02
 STEP TIME COMPLETED       3.600E+05,  TOTAL TIME COMPLETED        3.600E+05


  INCREMENT    37 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  373ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.118       TIME AVG. FORCE       6.052E-02
 LARGEST SCALED RESIDUAL FORCE      4.115E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.115E-10
 LARGEST INCREMENT OF DISP.         3.063E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.947E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  7.377E-05   TIME AVG. HEAT FLUX   3.753E-05
 LARGEST RESIDUAL HEAT FLUX         1.196E-15   AT NODE       4859   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.998E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.597E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  3.700E-02
 STEP TIME COMPLETED       3.700E+05,  TOTAL TIME COMPLETED        3.700E+05


  INCREMENT    38 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  27ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.121       TIME AVG. FORCE       6.211E-02
 LARGEST SCALED RESIDUAL FORCE      6.297E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.297E-09
 LARGEST INCREMENT OF DISP.         3.153E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.992E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.072E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.588E-05   TIME AVG. HEAT FLUX   3.856E-05
 LARGEST RESIDUAL HEAT FLUX         6.423E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.074E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.613E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.121       TIME AVG. FORCE       6.211E-02
 LARGEST SCALED RESIDUAL FORCE      1.073E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.073E-12
 LARGEST INCREMENT OF DISP.         3.153E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.377E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  7.589E-05   TIME AVG. HEAT FLUX   3.856E-05
 LARGEST RESIDUAL HEAT FLUX        -1.265E-15   AT NODE       4315   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.075E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.526E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  3.800E-02
 STEP TIME COMPLETED       3.800E+05,  TOTAL TIME COMPLETED        3.800E+05


  INCREMENT    39 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  377ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.124       TIME AVG. FORCE       6.370E-02
 LARGEST SCALED RESIDUAL FORCE      4.508E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.508E-10
 LARGEST INCREMENT OF DISP.         3.243E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.988E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  7.802E-05   TIME AVG. HEAT FLUX   3.960E-05
 LARGEST RESIDUAL HEAT FLUX         1.324E-15   AT NODE       2618   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.151E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.674E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  3.900E-02
 STEP TIME COMPLETED       3.900E+05,  TOTAL TIME COMPLETED        3.900E+05


  INCREMENT    40 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  24ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.127       TIME AVG. FORCE       6.528E-02
 LARGEST SCALED RESIDUAL FORCE      7.034E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      7.034E-09
 LARGEST INCREMENT OF DISP.         3.333E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.038E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.324E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.015E-05   TIME AVG. HEAT FLUX   4.064E-05
 LARGEST RESIDUAL HEAT FLUX         6.447E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.228E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.689E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  376ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.127       TIME AVG. FORCE       6.528E-02
 LARGEST SCALED RESIDUAL FORCE      1.178E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.178E-12
 LARGEST INCREMENT OF DISP.         3.334E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.465E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  8.015E-05   TIME AVG. HEAT FLUX   4.064E-05
 LARGEST RESIDUAL HEAT FLUX        -1.438E-15   AT NODE       4362   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.229E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.540E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  4.000E-02
 STEP TIME COMPLETED       4.000E+05,  TOTAL TIME COMPLETED        4.000E+05


  INCREMENT    41 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.130       TIME AVG. FORCE       6.687E-02
 LARGEST SCALED RESIDUAL FORCE      4.910E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.910E-10
 LARGEST INCREMENT OF DISP.         3.424E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.031E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  8.228E-05   TIME AVG. HEAT FLUX   4.168E-05
 LARGEST RESIDUAL HEAT FLUX         1.408E-15   AT NODE       1201   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.306E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.749E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  4.100E-02
 STEP TIME COMPLETED       4.100E+05,  TOTAL TIME COMPLETED        4.100E+05


  INCREMENT    42 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  24ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       6.845E-02
 LARGEST SCALED RESIDUAL FORCE      7.813E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      7.813E-09
 LARGEST INCREMENT OF DISP.         3.515E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.084E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.578E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.439E-05   TIME AVG. HEAT FLUX   4.272E-05
 LARGEST RESIDUAL HEAT FLUX         6.467E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.384E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.762E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  382ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       6.845E-02
 LARGEST SCALED RESIDUAL FORCE      1.286E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.286E-12
 LARGEST INCREMENT OF DISP.         3.515E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.555E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  8.440E-05   TIME AVG. HEAT FLUX   4.272E-05
 LARGEST RESIDUAL HEAT FLUX        -1.289E-15   AT NODE        136   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.384E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.553E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  4.200E-02
 STEP TIME COMPLETED       4.200E+05,  TOTAL TIME COMPLETED        4.200E+05


  INCREMENT    43 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  375ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.137       TIME AVG. FORCE       7.004E-02
 LARGEST SCALED RESIDUAL FORCE      5.319E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.319E-10
 LARGEST INCREMENT OF DISP.         3.606E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.075E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  8.654E-05   TIME AVG. HEAT FLUX   4.377E-05
 LARGEST RESIDUAL HEAT FLUX         1.243E-15   AT NODE         76   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.462E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.821E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  4.300E-02
 STEP TIME COMPLETED       4.300E+05,  TOTAL TIME COMPLETED        4.300E+05


  INCREMENT    44 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  24ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.140       TIME AVG. FORCE       7.162E-02
 LARGEST SCALED RESIDUAL FORCE      8.635E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      8.635E-09
 LARGEST INCREMENT OF DISP.         3.697E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.131E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.836E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.867E-05   TIME AVG. HEAT FLUX   4.481E-05
 LARGEST RESIDUAL HEAT FLUX         6.484E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.540E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.833E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  384ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.140       TIME AVG. FORCE       7.162E-02
 LARGEST SCALED RESIDUAL FORCE      1.397E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.397E-12
 LARGEST INCREMENT OF DISP.         3.697E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.646E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  8.867E-05   TIME AVG. HEAT FLUX   4.481E-05
 LARGEST RESIDUAL HEAT FLUX         1.297E-15   AT NODE       4989   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.541E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.564E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  4.400E-02
 STEP TIME COMPLETED       4.400E+05,  TOTAL TIME COMPLETED        4.400E+05


  INCREMENT    45 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  383ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.143       TIME AVG. FORCE       7.321E-02
 LARGEST SCALED RESIDUAL FORCE      5.732E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.732E-10
 LARGEST INCREMENT OF DISP.         3.788E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.119E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  9.080E-05   TIME AVG. HEAT FLUX   4.586E-05
 LARGEST RESIDUAL HEAT FLUX         1.445E-15   AT NODE       1448   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.620E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.890E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  4.500E-02
 STEP TIME COMPLETED       4.500E+05,  TOTAL TIME COMPLETED        4.500E+05


  INCREMENT    46 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.146       TIME AVG. FORCE       7.479E-02
 LARGEST SCALED RESIDUAL FORCE      9.501E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      9.501E-09
 LARGEST INCREMENT OF DISP.         3.880E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.179E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.096E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  9.294E-05   TIME AVG. HEAT FLUX   4.690E-05
 LARGEST RESIDUAL HEAT FLUX         6.498E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.699E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.901E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.146       TIME AVG. FORCE       7.479E-02
 LARGEST SCALED RESIDUAL FORCE      1.509E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.509E-12
 LARGEST INCREMENT OF DISP.         3.880E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.740E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  9.295E-05   TIME AVG. HEAT FLUX   4.690E-05
 LARGEST RESIDUAL HEAT FLUX         1.349E-15   AT NODE       4059   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.699E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.574E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  4.600E-02
 STEP TIME COMPLETED       4.600E+05,  TOTAL TIME COMPLETED        4.600E+05


  INCREMENT    47 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.149       TIME AVG. FORCE       7.637E-02
 LARGEST SCALED RESIDUAL FORCE      6.149E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.149E-10
 LARGEST INCREMENT OF DISP.         3.972E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.164E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  9.511E-05   TIME AVG. HEAT FLUX   4.795E-05
 LARGEST RESIDUAL HEAT FLUX        -1.427E-15   AT NODE       5412   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.778E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.956E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  4.700E-02
 STEP TIME COMPLETED       4.700E+05,  TOTAL TIME COMPLETED        4.700E+05


  INCREMENT    48 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.152       TIME AVG. FORCE       7.795E-02
 LARGEST SCALED RESIDUAL FORCE      1.041E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.041E-08
 LARGEST INCREMENT OF DISP.         4.064E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.228E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.360E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  9.726E-05   TIME AVG. HEAT FLUX   4.900E-05
 LARGEST RESIDUAL HEAT FLUX         6.507E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.858E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.966E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.152       TIME AVG. FORCE       7.795E-02
 LARGEST SCALED RESIDUAL FORCE      1.623E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.623E-12
 LARGEST INCREMENT OF DISP.         4.064E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.835E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  9.727E-05   TIME AVG. HEAT FLUX   4.900E-05
 LARGEST RESIDUAL HEAT FLUX        -1.452E-15   AT NODE       1074   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.858E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.582E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  4.800E-02
 STEP TIME COMPLETED       4.800E+05,  TOTAL TIME COMPLETED        4.800E+05


  INCREMENT    49 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  392ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.155       TIME AVG. FORCE       7.953E-02
 LARGEST SCALED RESIDUAL FORCE      6.568E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.568E-10
 LARGEST INCREMENT OF DISP.         4.156E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.210E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  9.945E-05   TIME AVG. HEAT FLUX   5.005E-05
 LARGEST RESIDUAL HEAT FLUX         1.472E-15   AT NODE      11759   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.939E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.018E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  4.900E-02
 STEP TIME COMPLETED       4.900E+05,  TOTAL TIME COMPLETED        4.900E+05


  INCREMENT    50 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.158       TIME AVG. FORCE       8.111E-02
 LARGEST SCALED RESIDUAL FORCE      1.136E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.136E-08
 LARGEST INCREMENT OF DISP.         4.249E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.278E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.627E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.016E-04   TIME AVG. HEAT FLUX   5.110E-05
 LARGEST RESIDUAL HEAT FLUX         6.512E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.019E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.028E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.158       TIME AVG. FORCE       8.111E-02
 LARGEST SCALED RESIDUAL FORCE      1.738E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.738E-12
 LARGEST INCREMENT OF DISP.         4.249E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.931E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.016E-04   TIME AVG. HEAT FLUX   5.110E-05
 LARGEST RESIDUAL HEAT FLUX        -1.397E-15   AT NODE        632   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.019E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.588E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  5.000E-02
 STEP TIME COMPLETED       5.000E+05,  TOTAL TIME COMPLETED        5.000E+05


  INCREMENT    51 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.162       TIME AVG. FORCE       8.268E-02
 LARGEST SCALED RESIDUAL FORCE      6.987E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.987E-10
 LARGEST INCREMENT OF DISP.         4.342E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.257E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.038E-04   TIME AVG. HEAT FLUX   5.216E-05
 LARGEST RESIDUAL HEAT FLUX         1.287E-15   AT NODE       4962   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.100E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.078E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  5.100E-02
 STEP TIME COMPLETED       5.100E+05,  TOTAL TIME COMPLETED        5.100E+05


  INCREMENT    52 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.165       TIME AVG. FORCE       8.426E-02
 LARGEST SCALED RESIDUAL FORCE      1.236E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.236E-08
 LARGEST INCREMENT OF DISP.         4.435E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.328E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.897E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.059E-04   TIME AVG. HEAT FLUX   5.321E-05
 LARGEST RESIDUAL HEAT FLUX         6.513E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.181E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.086E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  380ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.165       TIME AVG. FORCE       8.426E-02
 LARGEST SCALED RESIDUAL FORCE      1.854E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.854E-12
 LARGEST INCREMENT OF DISP.         4.435E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.028E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.059E-04   TIME AVG. HEAT FLUX   5.321E-05
 LARGEST RESIDUAL HEAT FLUX         1.400E-15   AT NODE       8733   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.181E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.593E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  5.200E-02
 STEP TIME COMPLETED       5.200E+05,  TOTAL TIME COMPLETED        5.200E+05


  INCREMENT    53 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  380ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.168       TIME AVG. FORCE       8.584E-02
 LARGEST SCALED RESIDUAL FORCE      7.405E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      7.405E-10
 LARGEST INCREMENT OF DISP.         4.529E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.305E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.081E-04   TIME AVG. HEAT FLUX   5.427E-05
 LARGEST RESIDUAL HEAT FLUX         1.293E-15   AT NODE        239   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.262E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.133E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  5.300E-02
 STEP TIME COMPLETED       5.300E+05,  TOTAL TIME COMPLETED        5.300E+05


  INCREMENT    54 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.171       TIME AVG. FORCE       8.741E-02
 LARGEST SCALED RESIDUAL FORCE      1.340E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.340E-08
 LARGEST INCREMENT OF DISP.         4.622E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.380E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.171E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.103E-04   TIME AVG. HEAT FLUX   5.532E-05
 LARGEST RESIDUAL HEAT FLUX         6.510E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.344E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.140E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  388ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.171       TIME AVG. FORCE       8.741E-02
 LARGEST SCALED RESIDUAL FORCE      1.970E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.970E-12
 LARGEST INCREMENT OF DISP.         4.623E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.125E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.103E-04   TIME AVG. HEAT FLUX   5.532E-05
 LARGEST RESIDUAL HEAT FLUX        -1.662E-15   AT NODE        594   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.344E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.596E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  5.400E-02
 STEP TIME COMPLETED       5.400E+05,  TOTAL TIME COMPLETED        5.400E+05


  INCREMENT    55 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.174       TIME AVG. FORCE       8.899E-02
 LARGEST SCALED RESIDUAL FORCE      7.820E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      7.820E-10
 LARGEST INCREMENT OF DISP.         4.716E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.354E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.125E-04   TIME AVG. HEAT FLUX   5.638E-05
 LARGEST RESIDUAL HEAT FLUX         1.667E-15   AT NODE        215   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.426E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.185E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  5.500E-02
 STEP TIME COMPLETED       5.500E+05,  TOTAL TIME COMPLETED        5.500E+05


  INCREMENT    56 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.177       TIME AVG. FORCE       9.056E-02
 LARGEST SCALED RESIDUAL FORCE      1.449E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.449E-08
 LARGEST INCREMENT OF DISP.         4.810E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.433E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.448E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.147E-04   TIME AVG. HEAT FLUX   5.744E-05
 LARGEST RESIDUAL HEAT FLUX         6.502E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.508E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.191E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  373ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.177       TIME AVG. FORCE       9.056E-02
 LARGEST SCALED RESIDUAL FORCE      2.086E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.086E-12
 LARGEST INCREMENT OF DISP.         4.811E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.222E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.147E-04   TIME AVG. HEAT FLUX   5.744E-05
 LARGEST RESIDUAL HEAT FLUX         1.322E-15   AT NODE       2290   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.508E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.596E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  5.600E-02
 STEP TIME COMPLETED       5.600E+05,  TOTAL TIME COMPLETED        5.600E+05


  INCREMENT    57 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.180       TIME AVG. FORCE       9.213E-02
 LARGEST SCALED RESIDUAL FORCE      8.231E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      8.231E-10
 LARGEST INCREMENT OF DISP.         4.905E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.404E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.169E-04   TIME AVG. HEAT FLUX   5.850E-05
 LARGEST RESIDUAL HEAT FLUX         1.301E-15   AT NODE       3803   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.590E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.233E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  5.700E-02
 STEP TIME COMPLETED       5.700E+05,  TOTAL TIME COMPLETED        5.700E+05


  INCREMENT    58 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.183       TIME AVG. FORCE       9.370E-02
 LARGEST SCALED RESIDUAL FORCE      1.563E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.563E-08
 LARGEST INCREMENT OF DISP.         5.000E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.487E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.729E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.190E-04   TIME AVG. HEAT FLUX   5.956E-05
 LARGEST RESIDUAL HEAT FLUX         6.490E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.673E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.237E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  350ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.183       TIME AVG. FORCE       9.370E-02
 LARGEST SCALED RESIDUAL FORCE      2.201E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.201E-12
 LARGEST INCREMENT OF DISP.         5.000E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.319E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.190E-04   TIME AVG. HEAT FLUX   5.956E-05
 LARGEST RESIDUAL HEAT FLUX         1.533E-15   AT NODE        629   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.673E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.595E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  5.800E-02
 STEP TIME COMPLETED       5.800E+05,  TOTAL TIME COMPLETED        5.800E+05


  INCREMENT    59 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.186       TIME AVG. FORCE       9.527E-02
 LARGEST SCALED RESIDUAL FORCE      8.636E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      8.636E-10
 LARGEST INCREMENT OF DISP.         5.094E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.455E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.212E-04   TIME AVG. HEAT FLUX   6.063E-05
 LARGEST RESIDUAL HEAT FLUX        -1.341E-15   AT NODE       1973   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.756E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.276E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  5.900E-02
 STEP TIME COMPLETED       5.900E+05,  TOTAL TIME COMPLETED        5.900E+05


  INCREMENT    60 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.189       TIME AVG. FORCE       9.684E-02
 LARGEST SCALED RESIDUAL FORCE      1.681E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.681E-08
 LARGEST INCREMENT OF DISP.         5.190E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.542E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.014E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.234E-04   TIME AVG. HEAT FLUX   6.169E-05
 LARGEST RESIDUAL HEAT FLUX         6.473E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.839E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.279E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  377ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.189       TIME AVG. FORCE       9.684E-02
 LARGEST SCALED RESIDUAL FORCE      2.315E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.315E-12
 LARGEST INCREMENT OF DISP.         5.190E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.415E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.234E-04   TIME AVG. HEAT FLUX   6.169E-05
 LARGEST RESIDUAL HEAT FLUX         1.438E-15   AT NODE       5199   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.839E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.592E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  6.000E-02
 STEP TIME COMPLETED       6.000E+05,  TOTAL TIME COMPLETED        6.000E+05


  INCREMENT    61 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  385ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.192       TIME AVG. FORCE       9.841E-02
 LARGEST SCALED RESIDUAL FORCE      9.033E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      9.033E-10
 LARGEST INCREMENT OF DISP.         5.285E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.507E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.257E-04   TIME AVG. HEAT FLUX   6.276E-05
 LARGEST RESIDUAL HEAT FLUX        -1.474E-15   AT NODE         55   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.922E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.315E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  6.100E-02
 STEP TIME COMPLETED       6.100E+05,  TOTAL TIME COMPLETED        6.100E+05


  INCREMENT    62 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.196       TIME AVG. FORCE       0.100    
 LARGEST SCALED RESIDUAL FORCE      1.804E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.804E-08
 LARGEST INCREMENT OF DISP.         5.381E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.598E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.303E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.279E-04   TIME AVG. HEAT FLUX   6.382E-05
 LARGEST RESIDUAL HEAT FLUX         6.451E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.001E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.317E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  380ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.196       TIME AVG. FORCE       0.100    
 LARGEST SCALED RESIDUAL FORCE      2.428E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.428E-12
 LARGEST INCREMENT OF DISP.         5.381E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.509E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.279E-04   TIME AVG. HEAT FLUX   6.382E-05
 LARGEST RESIDUAL HEAT FLUX         1.355E-15   AT NODE       2583   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.001E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.587E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  6.200E-02
 STEP TIME COMPLETED       6.200E+05,  TOTAL TIME COMPLETED        6.200E+05


  INCREMENT    63 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.199       TIME AVG. FORCE       0.102    
 LARGEST SCALED RESIDUAL FORCE      9.422E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      9.422E-10
 LARGEST INCREMENT OF DISP.         5.477E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.560E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.301E-04   TIME AVG. HEAT FLUX   6.489E-05
 LARGEST RESIDUAL HEAT FLUX        -1.359E-15   AT NODE       6195   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.009E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.348E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  6.300E-02
 STEP TIME COMPLETED       6.300E+05,  TOTAL TIME COMPLETED        6.300E+05


  INCREMENT    64 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.202       TIME AVG. FORCE       0.103    
 LARGEST SCALED RESIDUAL FORCE      1.932E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.932E-08
 LARGEST INCREMENT OF DISP.         5.573E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.655E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.596E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.323E-04   TIME AVG. HEAT FLUX   6.596E-05
 LARGEST RESIDUAL HEAT FLUX         6.424E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.017E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.349E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.202       TIME AVG. FORCE       0.103    
 LARGEST SCALED RESIDUAL FORCE      2.538E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.538E-12
 LARGEST INCREMENT OF DISP.         5.574E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.603E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.323E-04   TIME AVG. HEAT FLUX   6.596E-05
 LARGEST RESIDUAL HEAT FLUX        -1.356E-15   AT NODE       1799   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.017E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.579E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  6.400E-02
 STEP TIME COMPLETED       6.400E+05,  TOTAL TIME COMPLETED        6.400E+05


  INCREMENT    65 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.205       TIME AVG. FORCE       0.105    
 LARGEST SCALED RESIDUAL FORCE      9.799E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      9.799E-10
 LARGEST INCREMENT OF DISP.         5.670E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.614E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.345E-04   TIME AVG. HEAT FLUX   6.704E-05
 LARGEST RESIDUAL HEAT FLUX        -1.649E-15   AT NODE       2317   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.026E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.377E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  6.500E-02
 STEP TIME COMPLETED       6.500E+05,  TOTAL TIME COMPLETED        6.500E+05


  INCREMENT    66 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.208       TIME AVG. FORCE       0.106    
 LARGEST SCALED RESIDUAL FORCE      2.065E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.065E-08
 LARGEST INCREMENT OF DISP.         5.767E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.714E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.893E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.368E-04   TIME AVG. HEAT FLUX   6.811E-05
 LARGEST RESIDUAL HEAT FLUX         6.392E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.034E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.376E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  377ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.208       TIME AVG. FORCE       0.106    
 LARGEST SCALED RESIDUAL FORCE      2.646E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.646E-12
 LARGEST INCREMENT OF DISP.         5.767E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.694E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.368E-04   TIME AVG. HEAT FLUX   6.811E-05
 LARGEST RESIDUAL HEAT FLUX        -1.290E-15   AT NODE       1208   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.034E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.569E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  6.600E-02
 STEP TIME COMPLETED       6.600E+05,  TOTAL TIME COMPLETED        6.600E+05


  INCREMENT    67 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  380ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.211       TIME AVG. FORCE       0.108    
 LARGEST SCALED RESIDUAL FORCE      1.016E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.016E-09
 LARGEST INCREMENT OF DISP.         5.864E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.670E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.390E-04   TIME AVG. HEAT FLUX   6.918E-05
 LARGEST RESIDUAL HEAT FLUX        -1.441E-15   AT NODE       1157   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.042E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.400E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  6.700E-02
 STEP TIME COMPLETED       6.700E+05,  TOTAL TIME COMPLETED        6.700E+05


  INCREMENT    68 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.214       TIME AVG. FORCE       0.109    
 LARGEST SCALED RESIDUAL FORCE      2.202E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.202E-08
 LARGEST INCREMENT OF DISP.         5.962E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.773E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    8.194E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.413E-04   TIME AVG. HEAT FLUX   7.026E-05
 LARGEST RESIDUAL HEAT FLUX         6.355E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.051E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.398E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  383ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.214       TIME AVG. FORCE       0.109    
 LARGEST SCALED RESIDUAL FORCE      2.751E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.751E-12
 LARGEST INCREMENT OF DISP.         5.962E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.782E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.413E-04   TIME AVG. HEAT FLUX   7.026E-05
 LARGEST RESIDUAL HEAT FLUX         1.345E-15   AT NODE        691   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.051E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.557E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  6.800E-02
 STEP TIME COMPLETED       6.800E+05,  TOTAL TIME COMPLETED        6.800E+05


  INCREMENT    69 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  385ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.217       TIME AVG. FORCE       0.111    
 LARGEST SCALED RESIDUAL FORCE      1.052E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.052E-09
 LARGEST INCREMENT OF DISP.         6.059E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.727E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.435E-04   TIME AVG. HEAT FLUX   7.134E-05
 LARGEST RESIDUAL HEAT FLUX        -1.277E-15   AT NODE       4840   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.059E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.418E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  6.900E-02
 STEP TIME COMPLETED       6.900E+05,  TOTAL TIME COMPLETED        6.900E+05


  INCREMENT    70 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.220       TIME AVG. FORCE       0.112    
 LARGEST SCALED RESIDUAL FORCE      2.345E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.345E-08
 LARGEST INCREMENT OF DISP.         6.158E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.834E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    8.500E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.458E-04   TIME AVG. HEAT FLUX   7.241E-05
 LARGEST RESIDUAL HEAT FLUX         6.312E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.068E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.414E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.220       TIME AVG. FORCE       0.112    
 LARGEST SCALED RESIDUAL FORCE      2.852E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.852E-12
 LARGEST INCREMENT OF DISP.         6.158E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.868E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.458E-04   TIME AVG. HEAT FLUX   7.241E-05
 LARGEST RESIDUAL HEAT FLUX         1.489E-15   AT NODE       4386   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.068E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.543E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  7.000E-02
 STEP TIME COMPLETED       7.000E+05,  TOTAL TIME COMPLETED        7.000E+05


  INCREMENT    71 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  375ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.223       TIME AVG. FORCE       0.114    
 LARGEST SCALED RESIDUAL FORCE      1.085E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.085E-09
 LARGEST INCREMENT OF DISP.         6.256E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.785E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.480E-04   TIME AVG. HEAT FLUX   7.350E-05
 LARGEST RESIDUAL HEAT FLUX         1.355E-15   AT NODE        478   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.076E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.430E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  7.100E-02
 STEP TIME COMPLETED       7.100E+05,  TOTAL TIME COMPLETED        7.100E+05


  INCREMENT    72 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.226       TIME AVG. FORCE       0.116    
 LARGEST SCALED RESIDUAL FORCE      2.493E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.493E-08
 LARGEST INCREMENT OF DISP.         6.355E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.897E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    8.810E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.503E-04   TIME AVG. HEAT FLUX   7.458E-05
 LARGEST RESIDUAL HEAT FLUX         6.264E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.085E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.424E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.226       TIME AVG. FORCE       0.116    
 LARGEST SCALED RESIDUAL FORCE      2.949E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.949E-12
 LARGEST INCREMENT OF DISP.         6.355E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.951E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.503E-04   TIME AVG. HEAT FLUX   7.458E-05
 LARGEST RESIDUAL HEAT FLUX         1.333E-15   AT NODE       5260   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.085E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.526E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  7.200E-02
 STEP TIME COMPLETED       7.200E+05,  TOTAL TIME COMPLETED        7.200E+05


  INCREMENT    73 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.229       TIME AVG. FORCE       0.117    
 LARGEST SCALED RESIDUAL FORCE      1.117E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.117E-09
 LARGEST INCREMENT OF DISP.         6.453E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.845E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.526E-04   TIME AVG. HEAT FLUX   7.566E-05
 LARGEST RESIDUAL HEAT FLUX        -1.445E-15   AT NODE      10679   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.093E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.435E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  7.300E-02
 STEP TIME COMPLETED       7.300E+05,  TOTAL TIME COMPLETED        7.300E+05


  INCREMENT    74 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.232       TIME AVG. FORCE       0.119    
 LARGEST SCALED RESIDUAL FORCE      2.646E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.646E-08
 LARGEST INCREMENT OF DISP.         6.553E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.960E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.126E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.548E-04   TIME AVG. HEAT FLUX   7.674E-05
 LARGEST RESIDUAL HEAT FLUX         6.211E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.101E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.428E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  382ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.232       TIME AVG. FORCE       0.119    
 LARGEST SCALED RESIDUAL FORCE      3.041E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.041E-12
 LARGEST INCREMENT OF DISP.         6.553E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.030E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.548E-04   TIME AVG. HEAT FLUX   7.674E-05
 LARGEST RESIDUAL HEAT FLUX        -1.334E-15   AT NODE       1462   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.101E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.506E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  7.400E-02
 STEP TIME COMPLETED       7.400E+05,  TOTAL TIME COMPLETED        7.400E+05


  INCREMENT    75 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  379ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.235       TIME AVG. FORCE       0.120    
 LARGEST SCALED RESIDUAL FORCE      1.147E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.147E-09
 LARGEST INCREMENT OF DISP.         6.652E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.906E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.571E-04   TIME AVG. HEAT FLUX   7.783E-05
 LARGEST RESIDUAL HEAT FLUX        -1.519E-15   AT NODE        124   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.110E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.434E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  7.500E-02
 STEP TIME COMPLETED       7.500E+05,  TOTAL TIME COMPLETED        7.500E+05


  INCREMENT    76 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.238       TIME AVG. FORCE       0.122    
 LARGEST SCALED RESIDUAL FORCE      2.804E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.804E-08
 LARGEST INCREMENT OF DISP.         6.753E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.003E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.446E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.594E-04   TIME AVG. HEAT FLUX   7.892E-05
 LARGEST RESIDUAL HEAT FLUX         6.151E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.118E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.426E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  373ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.238       TIME AVG. FORCE       0.122    
 LARGEST SCALED RESIDUAL FORCE      3.128E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.128E-12
 LARGEST INCREMENT OF DISP.         6.753E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.105E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.594E-04   TIME AVG. HEAT FLUX   7.892E-05
 LARGEST RESIDUAL HEAT FLUX         1.343E-15   AT NODE       2443   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.118E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.484E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  7.600E-02
 STEP TIME COMPLETED       7.600E+05,  TOTAL TIME COMPLETED        7.600E+05


  INCREMENT    77 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  388ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.241       TIME AVG. FORCE       0.123    
 LARGEST SCALED RESIDUAL FORCE      1.175E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.175E-09
 LARGEST INCREMENT OF DISP.         6.853E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.969E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.617E-04   TIME AVG. HEAT FLUX   8.001E-05
 LARGEST RESIDUAL HEAT FLUX         1.468E-15   AT NODE       5001   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.127E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.427E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  7.700E-02
 STEP TIME COMPLETED       7.700E+05,  TOTAL TIME COMPLETED        7.700E+05


  INCREMENT    78 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.244       TIME AVG. FORCE       0.125    
 LARGEST SCALED RESIDUAL FORCE      2.968E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.968E-08
 LARGEST INCREMENT OF DISP.         6.954E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.009E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.771E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.640E-04   TIME AVG. HEAT FLUX   8.110E-05
 LARGEST RESIDUAL HEAT FLUX         6.086E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.135E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.416E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.244       TIME AVG. FORCE       0.125    
 LARGEST SCALED RESIDUAL FORCE      3.210E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.210E-12
 LARGEST INCREMENT OF DISP.         6.954E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.176E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.640E-04   TIME AVG. HEAT FLUX   8.110E-05
 LARGEST RESIDUAL HEAT FLUX         1.291E-15   AT NODE       5419   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.135E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.459E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  7.800E-02
 STEP TIME COMPLETED       7.800E+05,  TOTAL TIME COMPLETED        7.800E+05


  INCREMENT    79 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.247       TIME AVG. FORCE       0.126    
 LARGEST SCALED RESIDUAL FORCE      1.201E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.201E-09
 LARGEST INCREMENT OF DISP.         7.054E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.003E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.662E-04   TIME AVG. HEAT FLUX   8.219E-05
 LARGEST RESIDUAL HEAT FLUX        -1.443E-15   AT NODE       2417   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.144E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.412E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  7.900E-02
 STEP TIME COMPLETED       7.900E+05,  TOTAL TIME COMPLETED        7.900E+05


  INCREMENT    80 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.250       TIME AVG. FORCE       0.128    
 LARGEST SCALED RESIDUAL FORCE      3.137E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.137E-08
 LARGEST INCREMENT OF DISP.         7.156E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.016E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.010E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.685E-04   TIME AVG. HEAT FLUX   8.328E-05
 LARGEST RESIDUAL HEAT FLUX         6.015E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.152E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.400E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  379ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.250       TIME AVG. FORCE       0.128    
 LARGEST SCALED RESIDUAL FORCE      3.285E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.285E-12
 LARGEST INCREMENT OF DISP.         7.156E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.241E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.685E-04   TIME AVG. HEAT FLUX   8.328E-05
 LARGEST RESIDUAL HEAT FLUX         1.308E-15   AT NODE       1386   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.152E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.432E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  8.000E-02
 STEP TIME COMPLETED       8.000E+05,  TOTAL TIME COMPLETED        8.000E+05


  INCREMENT    81 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  372ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.253       TIME AVG. FORCE       0.129    
 LARGEST SCALED RESIDUAL FORCE      1.224E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.224E-09
 LARGEST INCREMENT OF DISP.         7.257E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.010E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.708E-04   TIME AVG. HEAT FLUX   8.438E-05
 LARGEST RESIDUAL HEAT FLUX        -1.338E-15   AT NODE       2567   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.160E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.390E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  8.100E-02
 STEP TIME COMPLETED       8.100E+05,  TOTAL TIME COMPLETED        8.100E+05


  INCREMENT    82 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.256       TIME AVG. FORCE       0.131    
 LARGEST SCALED RESIDUAL FORCE      3.312E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.312E-08
 LARGEST INCREMENT OF DISP.         7.359E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.023E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.044E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.731E-04   TIME AVG. HEAT FLUX   8.547E-05
 LARGEST RESIDUAL HEAT FLUX         5.938E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.169E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.376E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  391ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.256       TIME AVG. FORCE       0.131    
 LARGEST SCALED RESIDUAL FORCE      3.354E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.354E-12
 LARGEST INCREMENT OF DISP.         7.360E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.302E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.731E-04   TIME AVG. HEAT FLUX   8.547E-05
 LARGEST RESIDUAL HEAT FLUX         1.440E-15   AT NODE       4837   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.169E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.402E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  8.200E-02
 STEP TIME COMPLETED       8.200E+05,  TOTAL TIME COMPLETED        8.200E+05


  INCREMENT    83 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  389ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.259       TIME AVG. FORCE       0.133    
 LARGEST SCALED RESIDUAL FORCE      1.245E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.245E-09
 LARGEST INCREMENT OF DISP.         7.461E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.017E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.754E-04   TIME AVG. HEAT FLUX   8.657E-05
 LARGEST RESIDUAL HEAT FLUX         1.447E-15   AT NODE       1416   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.177E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.360E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  8.300E-02
 STEP TIME COMPLETED       8.300E+05,  TOTAL TIME COMPLETED        8.300E+05


  INCREMENT    84 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.262       TIME AVG. FORCE       0.134    
 LARGEST SCALED RESIDUAL FORCE      3.492E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.492E-08
 LARGEST INCREMENT OF DISP.         7.564E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.030E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.078E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.777E-04   TIME AVG. HEAT FLUX   8.767E-05
 LARGEST RESIDUAL HEAT FLUX         5.855E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.186E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.344E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.262       TIME AVG. FORCE       0.134    
 LARGEST SCALED RESIDUAL FORCE      3.415E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.415E-12
 LARGEST INCREMENT OF DISP.         7.565E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.356E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.777E-04   TIME AVG. HEAT FLUX   8.767E-05
 LARGEST RESIDUAL HEAT FLUX        -1.603E-15   AT NODE        550   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.186E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.369E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  8.400E-02
 STEP TIME COMPLETED       8.400E+05,  TOTAL TIME COMPLETED        8.400E+05


  INCREMENT    85 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.265       TIME AVG. FORCE       0.136    
 LARGEST SCALED RESIDUAL FORCE      1.263E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.263E-09
 LARGEST INCREMENT OF DISP.         7.667E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.024E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.801E-04   TIME AVG. HEAT FLUX   8.877E-05
 LARGEST RESIDUAL HEAT FLUX        -1.443E-15   AT NODE       4383   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.194E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.322E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  8.500E-02
 STEP TIME COMPLETED       8.500E+05,  TOTAL TIME COMPLETED        8.500E+05


  INCREMENT    86 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.267       TIME AVG. FORCE       0.137    
 LARGEST SCALED RESIDUAL FORCE     -3.727E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.727E-08
 LARGEST INCREMENT OF DISP.         7.771E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.037E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.127E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.824E-04   TIME AVG. HEAT FLUX   8.987E-05
 LARGEST RESIDUAL HEAT FLUX         5.766E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.202E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.304E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  381ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.267       TIME AVG. FORCE       0.137    
 LARGEST SCALED RESIDUAL FORCE      3.468E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.468E-12
 LARGEST INCREMENT OF DISP.         7.771E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.404E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.824E-04   TIME AVG. HEAT FLUX   8.987E-05
 LARGEST RESIDUAL HEAT FLUX        -1.239E-15   AT NODE       1305   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.202E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.333E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  8.600E-02
 STEP TIME COMPLETED       8.600E+05,  TOTAL TIME COMPLETED        8.600E+05


  INCREMENT    87 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  356ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.270       TIME AVG. FORCE       0.139    
 LARGEST SCALED RESIDUAL FORCE      1.279E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.279E-09
 LARGEST INCREMENT OF DISP.         7.874E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.031E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.848E-04   TIME AVG. HEAT FLUX   9.097E-05
 LARGEST RESIDUAL HEAT FLUX         1.283E-15   AT NODE       1134   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.211E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.276E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  8.700E-02
 STEP TIME COMPLETED       8.700E+05,  TOTAL TIME COMPLETED        8.700E+05


  INCREMENT    88 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  27ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.273       TIME AVG. FORCE       0.140    
 LARGEST SCALED RESIDUAL FORCE     -3.975E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.975E-08
 LARGEST INCREMENT OF DISP.         7.979E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.045E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.179E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.872E-04   TIME AVG. HEAT FLUX   9.208E-05
 LARGEST RESIDUAL HEAT FLUX         5.670E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.219E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.256E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.273       TIME AVG. FORCE       0.140    
 LARGEST SCALED RESIDUAL FORCE      3.513E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.513E-12
 LARGEST INCREMENT OF DISP.         7.979E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.446E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.872E-04   TIME AVG. HEAT FLUX   9.208E-05
 LARGEST RESIDUAL HEAT FLUX         1.322E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.219E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.294E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  8.800E-02
 STEP TIME COMPLETED       8.800E+05,  TOTAL TIME COMPLETED        8.800E+05


  INCREMENT    89 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  372ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.276       TIME AVG. FORCE       0.142    
 LARGEST SCALED RESIDUAL FORCE      1.291E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.291E-09
 LARGEST INCREMENT OF DISP.         8.083E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.038E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.895E-04   TIME AVG. HEAT FLUX   9.319E-05
 LARGEST RESIDUAL HEAT FLUX        -1.609E-15   AT NODE       2286   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.227E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.221E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  8.900E-02
 STEP TIME COMPLETED       8.900E+05,  TOTAL TIME COMPLETED        8.900E+05


  INCREMENT    90 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.279       TIME AVG. FORCE       0.143    
 LARGEST SCALED RESIDUAL FORCE     -4.236E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.236E-08
 LARGEST INCREMENT OF DISP.         8.188E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.052E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.233E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.919E-04   TIME AVG. HEAT FLUX   9.429E-05
 LARGEST RESIDUAL HEAT FLUX         5.568E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.235E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.199E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  373ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.279       TIME AVG. FORCE       0.143    
 LARGEST SCALED RESIDUAL FORCE      3.550E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.550E-12
 LARGEST INCREMENT OF DISP.         8.188E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.480E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.919E-04   TIME AVG. HEAT FLUX   9.429E-05
 LARGEST RESIDUAL HEAT FLUX         1.752E-15   AT NODE       3847   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.235E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.252E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.000E-02
 STEP TIME COMPLETED       9.000E+05,  TOTAL TIME COMPLETED        9.000E+05


  INCREMENT    91 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.282       TIME AVG. FORCE       0.145    
 LARGEST SCALED RESIDUAL FORCE      1.310E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.310E-09
 LARGEST INCREMENT OF DISP.         8.293E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.045E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.942E-04   TIME AVG. HEAT FLUX   9.541E-05
 LARGEST RESIDUAL HEAT FLUX        -2.470E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.243E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.158E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.100E-02
 STEP TIME COMPLETED       9.100E+05,  TOTAL TIME COMPLETED        9.100E+05


  INCREMENT    92 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.285       TIME AVG. FORCE       0.146    
 LARGEST SCALED RESIDUAL FORCE     -4.511E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.511E-08
 LARGEST INCREMENT OF DISP.         8.399E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.060E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.289E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.966E-04   TIME AVG. HEAT FLUX   9.652E-05
 LARGEST RESIDUAL HEAT FLUX         5.460E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.252E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.133E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.285       TIME AVG. FORCE       0.146    
 LARGEST SCALED RESIDUAL FORCE      3.577E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.577E-12
 LARGEST INCREMENT OF DISP.         8.399E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.507E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.966E-04   TIME AVG. HEAT FLUX   9.652E-05
 LARGEST RESIDUAL HEAT FLUX         1.388E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.252E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.207E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.200E-02
 STEP TIME COMPLETED       9.200E+05,  TOTAL TIME COMPLETED        9.200E+05


  INCREMENT    93 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.288       TIME AVG. FORCE       0.148    
 LARGEST SCALED RESIDUAL FORCE      1.336E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.336E-09
 LARGEST INCREMENT OF DISP.         8.505E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.053E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.990E-04   TIME AVG. HEAT FLUX   9.763E-05
 LARGEST RESIDUAL HEAT FLUX         1.808E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.260E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.085E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.300E-02
 STEP TIME COMPLETED       9.300E+05,  TOTAL TIME COMPLETED        9.300E+05


  INCREMENT    94 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.291       TIME AVG. FORCE       0.149    
 LARGEST SCALED RESIDUAL FORCE     -4.800E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.800E-08
 LARGEST INCREMENT OF DISP.         8.611E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.068E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.347E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.014E-04   TIME AVG. HEAT FLUX   9.875E-05
 LARGEST RESIDUAL HEAT FLUX         5.345E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.268E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.058E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.291       TIME AVG. FORCE       0.149    
 LARGEST SCALED RESIDUAL FORCE      3.594E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.594E-12
 LARGEST INCREMENT OF DISP.         8.612E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.525E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.014E-04   TIME AVG. HEAT FLUX   9.875E-05
 LARGEST RESIDUAL HEAT FLUX         1.470E-15   AT NODE       2777   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.268E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.159E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.400E-02
 STEP TIME COMPLETED       9.400E+05,  TOTAL TIME COMPLETED        9.400E+05


  INCREMENT    95 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.293       TIME AVG. FORCE       0.151    
 LARGEST SCALED RESIDUAL FORCE      1.360E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.360E-09
 LARGEST INCREMENT OF DISP.         8.718E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.061E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.038E-04   TIME AVG. HEAT FLUX   9.986E-05
 LARGEST RESIDUAL HEAT FLUX         1.471E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.276E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.002E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.500E-02
 STEP TIME COMPLETED       9.500E+05,  TOTAL TIME COMPLETED        9.500E+05


  INCREMENT    96 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  20ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.296       TIME AVG. FORCE       0.152    
 LARGEST SCALED RESIDUAL FORCE     -5.105E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.105E-08
 LARGEST INCREMENT OF DISP.         8.825E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.076E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.408E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.062E-04   TIME AVG. HEAT FLUX   1.010E-04
 LARGEST RESIDUAL HEAT FLUX         5.224E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.284E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.973E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  371ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.296       TIME AVG. FORCE       0.152    
 LARGEST SCALED RESIDUAL FORCE      3.601E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.601E-12
 LARGEST INCREMENT OF DISP.         8.826E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.535E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.062E-04   TIME AVG. HEAT FLUX   1.010E-04
 LARGEST RESIDUAL HEAT FLUX         1.262E-15   AT NODE       2554   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.284E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.109E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.600E-02
 STEP TIME COMPLETED       9.600E+05,  TOTAL TIME COMPLETED        9.600E+05


  INCREMENT    97 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.299       TIME AVG. FORCE       0.154    
 LARGEST SCALED RESIDUAL FORCE      1.382E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.382E-09
 LARGEST INCREMENT OF DISP.         8.933E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.069E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.086E-04   TIME AVG. HEAT FLUX   1.021E-04
 LARGEST RESIDUAL HEAT FLUX         1.279E-15   AT NODE       4198   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.292E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.910E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.700E-02
 STEP TIME COMPLETED       9.700E+05,  TOTAL TIME COMPLETED        9.700E+05


  INCREMENT    98 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.302       TIME AVG. FORCE       0.155    
 LARGEST SCALED RESIDUAL FORCE     -5.425E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.425E-08
 LARGEST INCREMENT OF DISP.         9.041E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.085E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.471E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.110E-04   TIME AVG. HEAT FLUX   1.032E-04
 LARGEST RESIDUAL HEAT FLUX         5.097E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.300E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.878E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.800E-02
 STEP TIME COMPLETED       9.800E+05,  TOTAL TIME COMPLETED        9.800E+05


  INCREMENT    99 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  386ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.305       TIME AVG. FORCE       0.157    
 LARGEST SCALED RESIDUAL FORCE      1.408E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.408E-09
 LARGEST INCREMENT OF DISP.         9.150E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.085E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.134E-04   TIME AVG. HEAT FLUX   1.044E-04
 LARGEST RESIDUAL HEAT FLUX        -1.382E-15   AT NODE       4560   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.307E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.848E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.900E-02
 STEP TIME COMPLETED       9.900E+05,  TOTAL TIME COMPLETED        9.900E+05


  INCREMENT   100 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.308       TIME AVG. FORCE       0.158    
 LARGEST SCALED RESIDUAL FORCE     -5.787E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.787E-08
 LARGEST INCREMENT OF DISP.         9.259E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.090E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.541E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.159E-04   TIME AVG. HEAT FLUX   1.055E-04
 LARGEST RESIDUAL HEAT FLUX         4.948E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.315E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.747E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.100    
 STEP TIME COMPLETED       1.000E+06,  TOTAL TIME COMPLETED        1.000E+06


  INCREMENT   101 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.311       TIME AVG. FORCE       0.160    
 LARGEST SCALED RESIDUAL FORCE      1.425E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.425E-09
 LARGEST INCREMENT OF DISP.         9.368E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.093E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.182E-04   TIME AVG. HEAT FLUX   1.066E-04
 LARGEST RESIDUAL HEAT FLUX        -1.476E-15   AT NODE       2140   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.323E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.741E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.101    
 STEP TIME COMPLETED       1.010E+06,  TOTAL TIME COMPLETED        1.010E+06


  INCREMENT   102 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.313       TIME AVG. FORCE       0.161    
 LARGEST SCALED RESIDUAL FORCE     -6.142E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.142E-08
 LARGEST INCREMENT OF DISP.         9.478E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.099E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.610E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.206E-04   TIME AVG. HEAT FLUX   1.077E-04
 LARGEST RESIDUAL HEAT FLUX         4.805E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.330E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.628E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.102    
 STEP TIME COMPLETED       1.020E+06,  TOTAL TIME COMPLETED        1.020E+06


  INCREMENT   103 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  378ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.316       TIME AVG. FORCE       0.163    
 LARGEST SCALED RESIDUAL FORCE      1.439E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.439E-09
 LARGEST INCREMENT OF DISP.         9.588E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.102E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.230E-04   TIME AVG. HEAT FLUX   1.089E-04
 LARGEST RESIDUAL HEAT FLUX        -1.293E-15   AT NODE        492   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.338E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.615E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.103    
 STEP TIME COMPLETED       1.030E+06,  TOTAL TIME COMPLETED        1.030E+06


  INCREMENT   104 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.319       TIME AVG. FORCE       0.164    
 LARGEST SCALED RESIDUAL FORCE     -6.514E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.514E-08
 LARGEST INCREMENT OF DISP.         9.699E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.108E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.682E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.254E-04   TIME AVG. HEAT FLUX   1.100E-04
 LARGEST RESIDUAL HEAT FLUX         4.659E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.346E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.501E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.104    
 STEP TIME COMPLETED       1.040E+06,  TOTAL TIME COMPLETED        1.040E+06


  INCREMENT   105 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.322       TIME AVG. FORCE       0.166    
 LARGEST SCALED RESIDUAL FORCE      1.449E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.449E-09
 LARGEST INCREMENT OF DISP.         9.810E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.111E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.279E-04   TIME AVG. HEAT FLUX   1.111E-04
 LARGEST RESIDUAL HEAT FLUX         1.416E-15   AT NODE       2407   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.353E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.478E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.105    
 STEP TIME COMPLETED       1.050E+06,  TOTAL TIME COMPLETED        1.050E+06


  INCREMENT   106 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.325       TIME AVG. FORCE       0.167    
 LARGEST SCALED RESIDUAL FORCE     -6.907E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.907E-08
 LARGEST INCREMENT OF DISP.         9.922E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.117E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.756E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.303E-04   TIME AVG. HEAT FLUX   1.123E-04
 LARGEST RESIDUAL HEAT FLUX         4.506E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.360E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.363E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.106    
 STEP TIME COMPLETED       1.060E+06,  TOTAL TIME COMPLETED        1.060E+06


  INCREMENT   107 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  382ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.327       TIME AVG. FORCE       0.169    
 LARGEST SCALED RESIDUAL FORCE      1.456E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.456E-09
 LARGEST INCREMENT OF DISP.         1.003E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.120E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.327E-04   TIME AVG. HEAT FLUX   1.134E-04
 LARGEST RESIDUAL HEAT FLUX        -1.565E-15   AT NODE      12450   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.368E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.329E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.107    
 STEP TIME COMPLETED       1.070E+06,  TOTAL TIME COMPLETED        1.070E+06


  INCREMENT   108 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.330       TIME AVG. FORCE       0.170    
 LARGEST SCALED RESIDUAL FORCE     -7.319E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.319E-08
 LARGEST INCREMENT OF DISP.         1.015E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.127E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.835E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.351E-04   TIME AVG. HEAT FLUX   1.145E-04
 LARGEST RESIDUAL HEAT FLUX         4.347E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.375E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.213E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.108    
 STEP TIME COMPLETED       1.080E+06,  TOTAL TIME COMPLETED        1.080E+06


  INCREMENT   109 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.333       TIME AVG. FORCE       0.172    
 LARGEST SCALED RESIDUAL FORCE      1.460E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.460E-09
 LARGEST INCREMENT OF DISP.         1.026E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.130E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.375E-04   TIME AVG. HEAT FLUX   1.157E-04
 LARGEST RESIDUAL HEAT FLUX         1.411E-15   AT NODE       3413   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.382E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.167E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.109    
 STEP TIME COMPLETED       1.090E+06,  TOTAL TIME COMPLETED        1.090E+06


  INCREMENT   110 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.336       TIME AVG. FORCE       0.173    
 LARGEST SCALED RESIDUAL FORCE     -7.753E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.753E-08
 LARGEST INCREMENT OF DISP.         1.037E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.137E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.916E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.400E-04   TIME AVG. HEAT FLUX   1.168E-04
 LARGEST RESIDUAL HEAT FLUX         4.182E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.389E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.051E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.110    
 STEP TIME COMPLETED       1.100E+06,  TOTAL TIME COMPLETED        1.100E+06


  INCREMENT   111 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  389ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.338       TIME AVG. FORCE       0.175    
 LARGEST SCALED RESIDUAL FORCE      1.460E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.460E-09
 LARGEST INCREMENT OF DISP.         1.049E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.140E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.425E-04   TIME AVG. HEAT FLUX   1.180E-04
 LARGEST RESIDUAL HEAT FLUX         1.370E-15   AT NODE       5274   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.396E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.993E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.111    
 STEP TIME COMPLETED       1.110E+06,  TOTAL TIME COMPLETED        1.110E+06


  INCREMENT   112 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.341       TIME AVG. FORCE       0.176    
 LARGEST SCALED RESIDUAL FORCE     -8.209E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.209E-08
 LARGEST INCREMENT OF DISP.         1.060E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.147E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.002E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.449E-04   TIME AVG. HEAT FLUX   1.191E-04
 LARGEST RESIDUAL HEAT FLUX         4.010E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.403E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.877E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.112    
 STEP TIME COMPLETED       1.120E+06,  TOTAL TIME COMPLETED        1.120E+06


  INCREMENT   113 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  383ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.344       TIME AVG. FORCE       0.178    
 LARGEST SCALED RESIDUAL FORCE      1.456E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.456E-09
 LARGEST INCREMENT OF DISP.         1.072E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.150E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.473E-04   TIME AVG. HEAT FLUX   1.202E-04
 LARGEST RESIDUAL HEAT FLUX         1.428E-15   AT NODE       5381   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.410E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.806E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.113    
 STEP TIME COMPLETED       1.130E+06,  TOTAL TIME COMPLETED        1.130E+06


  INCREMENT   114 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.347       TIME AVG. FORCE       0.179    
 LARGEST SCALED RESIDUAL FORCE     -8.690E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.690E-08
 LARGEST INCREMENT OF DISP.         1.083E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.158E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.091E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.498E-04   TIME AVG. HEAT FLUX   1.214E-04
 LARGEST RESIDUAL HEAT FLUX         3.833E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.417E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.689E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.114    
 STEP TIME COMPLETED       1.140E+06,  TOTAL TIME COMPLETED        1.140E+06


  INCREMENT   115 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  386ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.349       TIME AVG. FORCE       0.181    
 LARGEST SCALED RESIDUAL FORCE      1.449E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.449E-09
 LARGEST INCREMENT OF DISP.         1.095E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.161E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.523E-04   TIME AVG. HEAT FLUX   1.225E-04
 LARGEST RESIDUAL HEAT FLUX         1.597E-15   AT NODE       2361   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.423E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.606E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.115    
 STEP TIME COMPLETED       1.150E+06,  TOTAL TIME COMPLETED        1.150E+06


  INCREMENT   116 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.352       TIME AVG. FORCE       0.182    
 LARGEST SCALED RESIDUAL FORCE     -9.195E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.195E-08
 LARGEST INCREMENT OF DISP.         1.107E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.169E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.185E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.547E-04   TIME AVG. HEAT FLUX   1.237E-04
 LARGEST RESIDUAL HEAT FLUX         3.649E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.430E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.489E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.116    
 STEP TIME COMPLETED       1.160E+06,  TOTAL TIME COMPLETED        1.160E+06


  INCREMENT   117 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  394ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.355       TIME AVG. FORCE       0.184    
 LARGEST SCALED RESIDUAL FORCE      1.436E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.436E-09
 LARGEST INCREMENT OF DISP.         1.118E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.171E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.572E-04   TIME AVG. HEAT FLUX   1.248E-04
 LARGEST RESIDUAL HEAT FLUX         1.402E-15   AT NODE      13977   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.436E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.391E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.117    
 STEP TIME COMPLETED       1.170E+06,  TOTAL TIME COMPLETED        1.170E+06


  INCREMENT   118 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.357       TIME AVG. FORCE       0.185    
 LARGEST SCALED RESIDUAL FORCE     -9.727E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.727E-08
 LARGEST INCREMENT OF DISP.         1.130E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.180E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.283E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.597E-04   TIME AVG. HEAT FLUX   1.260E-04
 LARGEST RESIDUAL HEAT FLUX         3.460E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.442E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.275E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.118    
 STEP TIME COMPLETED       1.180E+06,  TOTAL TIME COMPLETED        1.180E+06


  INCREMENT   119 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.360       TIME AVG. FORCE       0.187    
 LARGEST SCALED RESIDUAL FORCE      1.420E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.420E-09
 LARGEST INCREMENT OF DISP.         1.142E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.183E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.621E-04   TIME AVG. HEAT FLUX   1.271E-04
 LARGEST RESIDUAL HEAT FLUX         1.316E-15   AT NODE       1930   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.448E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.163E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.119    
 STEP TIME COMPLETED       1.190E+06,  TOTAL TIME COMPLETED        1.190E+06


  INCREMENT   120 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.362       TIME AVG. FORCE       0.188    
 LARGEST SCALED RESIDUAL FORCE     -1.029E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.029E-07
 LARGEST INCREMENT OF DISP.         1.154E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.191E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.386E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.646E-04   TIME AVG. HEAT FLUX   1.283E-04
 LARGEST RESIDUAL HEAT FLUX         3.264E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.455E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.047E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.120    
 STEP TIME COMPLETED       1.200E+06,  TOTAL TIME COMPLETED        1.200E+06


  INCREMENT   121 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  377ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.365       TIME AVG. FORCE       0.190    
 LARGEST SCALED RESIDUAL FORCE      1.399E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.399E-09
 LARGEST INCREMENT OF DISP.         1.166E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.194E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.671E-04   TIME AVG. HEAT FLUX   1.295E-04
 LARGEST RESIDUAL HEAT FLUX        -1.349E-15   AT NODE        401   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.460E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.920E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.121    
 STEP TIME COMPLETED       1.210E+06,  TOTAL TIME COMPLETED        1.210E+06


  INCREMENT   122 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.368       TIME AVG. FORCE       0.191    
 LARGEST SCALED RESIDUAL FORCE     -1.088E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.088E-07
 LARGEST INCREMENT OF DISP.         1.178E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.203E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.493E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.696E-04   TIME AVG. HEAT FLUX   1.306E-04
 LARGEST RESIDUAL HEAT FLUX         3.063E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.466E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.805E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.122    
 STEP TIME COMPLETED       1.220E+06,  TOTAL TIME COMPLETED        1.220E+06


  INCREMENT   123 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  386ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.370       TIME AVG. FORCE       0.193    
 LARGEST SCALED RESIDUAL FORCE      1.373E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.373E-09
 LARGEST INCREMENT OF DISP.         1.190E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.206E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.722E-04   TIME AVG. HEAT FLUX   1.318E-04
 LARGEST RESIDUAL HEAT FLUX         1.478E-15   AT NODE       3587   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.472E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.662E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.123    
 STEP TIME COMPLETED       1.230E+06,  TOTAL TIME COMPLETED        1.230E+06


  INCREMENT   124 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.373       TIME AVG. FORCE       0.194    
 LARGEST SCALED RESIDUAL FORCE     -1.150E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.150E-07
 LARGEST INCREMENT OF DISP.         1.202E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.216E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.607E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.746E-04   TIME AVG. HEAT FLUX   1.329E-04
 LARGEST RESIDUAL HEAT FLUX         2.857E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.477E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.548E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.124    
 STEP TIME COMPLETED       1.240E+06,  TOTAL TIME COMPLETED        1.240E+06


  INCREMENT   125 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.376       TIME AVG. FORCE       0.195    
 LARGEST SCALED RESIDUAL FORCE      1.341E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.341E-09
 LARGEST INCREMENT OF DISP.         1.214E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.218E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.771E-04   TIME AVG. HEAT FLUX   1.341E-04
 LARGEST RESIDUAL HEAT FLUX         1.316E-15   AT NODE       4135   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.483E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.389E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.125    
 STEP TIME COMPLETED       1.250E+06,  TOTAL TIME COMPLETED        1.250E+06


  INCREMENT   126 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.378       TIME AVG. FORCE       0.197    
 LARGEST SCALED RESIDUAL FORCE     -1.215E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.215E-07
 LARGEST INCREMENT OF DISP.         1.226E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.228E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.726E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.796E-04   TIME AVG. HEAT FLUX   1.353E-04
 LARGEST RESIDUAL HEAT FLUX         2.645E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.488E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.276E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.126    
 STEP TIME COMPLETED       1.260E+06,  TOTAL TIME COMPLETED        1.260E+06


  INCREMENT   127 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.381       TIME AVG. FORCE       0.198    
 LARGEST SCALED RESIDUAL FORCE      1.305E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.305E-09
 LARGEST INCREMENT OF DISP.         1.239E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.231E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.821E-04   TIME AVG. HEAT FLUX   1.364E-04
 LARGEST RESIDUAL HEAT FLUX        -1.380E-15   AT NODE       1189   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.493E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.100E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.127    
 STEP TIME COMPLETED       1.270E+06,  TOTAL TIME COMPLETED        1.270E+06


  INCREMENT   128 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.383       TIME AVG. FORCE       0.200    
 LARGEST SCALED RESIDUAL FORCE     -1.284E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.284E-07
 LARGEST INCREMENT OF DISP.         1.251E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.242E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.846E-04   TIME AVG. HEAT FLUX   1.376E-04
 LARGEST RESIDUAL HEAT FLUX         2.429E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.498E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.989E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.128    
 STEP TIME COMPLETED       1.280E+06,  TOTAL TIME COMPLETED        1.280E+06


  INCREMENT   129 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.386       TIME AVG. FORCE       0.201    
 LARGEST SCALED RESIDUAL FORCE      1.262E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.262E-09
 LARGEST INCREMENT OF DISP.         1.264E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.244E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.871E-04   TIME AVG. HEAT FLUX   1.388E-04
 LARGEST RESIDUAL HEAT FLUX        -1.243E-15   AT NODE       1266   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.503E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.795E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.129    
 STEP TIME COMPLETED       1.290E+06,  TOTAL TIME COMPLETED        1.290E+06


  INCREMENT   130 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.388       TIME AVG. FORCE       0.203    
 LARGEST SCALED RESIDUAL FORCE     -1.357E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.357E-07
 LARGEST INCREMENT OF DISP.         1.276E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.255E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.896E-04   TIME AVG. HEAT FLUX   1.399E-04
 LARGEST RESIDUAL HEAT FLUX         2.207E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.508E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.686E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.130    
 STEP TIME COMPLETED       1.300E+06,  TOTAL TIME COMPLETED        1.300E+06


  INCREMENT   131 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.391       TIME AVG. FORCE       0.204    
 LARGEST SCALED RESIDUAL FORCE      1.214E-09   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.214E-09
 LARGEST INCREMENT OF DISP.         1.289E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.258E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.920E-04   TIME AVG. HEAT FLUX   1.411E-04
 LARGEST RESIDUAL HEAT FLUX         1.375E-15   AT NODE       1064   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.512E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.474E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.131    
 STEP TIME COMPLETED       1.310E+06,  TOTAL TIME COMPLETED        1.310E+06


  INCREMENT   132 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.394       TIME AVG. FORCE       0.206    
 LARGEST SCALED RESIDUAL FORCE     -1.434E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.434E-07
 LARGEST INCREMENT OF DISP.         1.301E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.269E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.945E-04   TIME AVG. HEAT FLUX   1.423E-04
 LARGEST RESIDUAL HEAT FLUX         1.981E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.517E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.367E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.132    
 STEP TIME COMPLETED       1.320E+06,  TOTAL TIME COMPLETED        1.320E+06


  INCREMENT   133 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.396       TIME AVG. FORCE       0.207    
 LARGEST SCALED RESIDUAL FORCE     -1.176E-09   AT NODE       6990   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.176E-09
 LARGEST INCREMENT OF DISP.         1.314E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.272E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.971E-04   TIME AVG. HEAT FLUX   1.434E-04
 LARGEST RESIDUAL HEAT FLUX         1.450E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.521E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.136E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.133    
 STEP TIME COMPLETED       1.330E+06,  TOTAL TIME COMPLETED        1.330E+06


  INCREMENT   134 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.399       TIME AVG. FORCE       0.208    
 LARGEST SCALED RESIDUAL FORCE     -1.514E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.514E-07
 LARGEST INCREMENT OF DISP.         1.327E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.284E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.995E-04   TIME AVG. HEAT FLUX   1.446E-04
 LARGEST RESIDUAL HEAT FLUX         1.751E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.525E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.031E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.134    
 STEP TIME COMPLETED       1.340E+06,  TOTAL TIME COMPLETED        1.340E+06


  INCREMENT   135 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  371ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.401       TIME AVG. FORCE       0.210    
 LARGEST SCALED RESIDUAL FORCE     -1.152E-09   AT NODE       6990   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.152E-09
 LARGEST INCREMENT OF DISP.         1.340E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.286E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.020E-04   TIME AVG. HEAT FLUX   1.458E-04
 LARGEST RESIDUAL HEAT FLUX        -1.298E-15   AT NODE       3919   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.528E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.780E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.135    
 STEP TIME COMPLETED       1.350E+06,  TOTAL TIME COMPLETED        1.350E+06


  INCREMENT   136 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.404       TIME AVG. FORCE       0.211    
 LARGEST SCALED RESIDUAL FORCE     -1.600E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.600E-07
 LARGEST INCREMENT OF DISP.         1.353E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.299E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.046E-04   TIME AVG. HEAT FLUX   1.470E-04
 LARGEST RESIDUAL HEAT FLUX         1.517E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.532E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.679E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.136    
 STEP TIME COMPLETED       1.360E+06,  TOTAL TIME COMPLETED        1.360E+06


  INCREMENT   137 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  352ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.406       TIME AVG. FORCE       0.213    
 LARGEST SCALED RESIDUAL FORCE     -1.123E-09   AT NODE       6990   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.123E-09
 LARGEST INCREMENT OF DISP.         1.366E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.301E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.071E-04   TIME AVG. HEAT FLUX   1.482E-04
 LARGEST RESIDUAL HEAT FLUX        -1.373E-15   AT NODE       4628   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.536E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.407E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.137    
 STEP TIME COMPLETED       1.370E+06,  TOTAL TIME COMPLETED        1.370E+06


  INCREMENT   138 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.409       TIME AVG. FORCE       0.214    
 LARGEST SCALED RESIDUAL FORCE     -1.690E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.690E-07
 LARGEST INCREMENT OF DISP.         1.379E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.315E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.096E-04   TIME AVG. HEAT FLUX   1.493E-04
 LARGEST RESIDUAL HEAT FLUX         1.280E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.539E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.309E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.138    
 STEP TIME COMPLETED       1.380E+06,  TOTAL TIME COMPLETED        1.380E+06


  INCREMENT   139 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.411       TIME AVG. FORCE       0.216    
 LARGEST SCALED RESIDUAL FORCE     -1.089E-09   AT NODE       6990   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.089E-09
 LARGEST INCREMENT OF DISP.         1.392E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.317E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.122E-04   TIME AVG. HEAT FLUX   1.505E-04
 LARGEST RESIDUAL HEAT FLUX        -1.381E-15   AT NODE       2333   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.542E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.016E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.139    
 STEP TIME COMPLETED       1.390E+06,  TOTAL TIME COMPLETED        1.390E+06


  INCREMENT   140 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.414       TIME AVG. FORCE       0.217    
 LARGEST SCALED RESIDUAL FORCE     -1.785E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.785E-07
 LARGEST INCREMENT OF DISP.         1.406E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.331E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.147E-04   TIME AVG. HEAT FLUX   1.517E-04
 LARGEST RESIDUAL HEAT FLUX         1.039E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.545E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.922E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.140    
 STEP TIME COMPLETED       1.400E+06,  TOTAL TIME COMPLETED        1.400E+06


  INCREMENT   141 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  345ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.416       TIME AVG. FORCE       0.218    
 LARGEST SCALED RESIDUAL FORCE     -1.050E-09   AT NODE       6990   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.050E-09
 LARGEST INCREMENT OF DISP.         1.419E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.333E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.172E-04   TIME AVG. HEAT FLUX   1.529E-04
 LARGEST RESIDUAL HEAT FLUX        -1.627E-15   AT NODE       5212   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.547E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.607E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.141    
 STEP TIME COMPLETED       1.410E+06,  TOTAL TIME COMPLETED        1.410E+06


  INCREMENT   142 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.418       TIME AVG. FORCE       0.220    
 LARGEST SCALED RESIDUAL FORCE     -1.886E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.886E-07
 LARGEST INCREMENT OF DISP.         1.432E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.348E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.196E-04   TIME AVG. HEAT FLUX   1.541E-04
 LARGEST RESIDUAL HEAT FLUX         7.956E-08   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.550E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.517E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.142    
 STEP TIME COMPLETED       1.420E+06,  TOTAL TIME COMPLETED        1.420E+06


  INCREMENT   143 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.421       TIME AVG. FORCE       0.221    
 LARGEST SCALED RESIDUAL FORCE     -1.006E-09   AT NODE       6990   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.006E-09
 LARGEST INCREMENT OF DISP.         1.446E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.350E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.221E-04   TIME AVG. HEAT FLUX   1.552E-04
 LARGEST RESIDUAL HEAT FLUX         1.307E-15   AT NODE       4618   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.552E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.179E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.143    
 STEP TIME COMPLETED       1.430E+06,  TOTAL TIME COMPLETED        1.430E+06


  INCREMENT   144 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.423       TIME AVG. FORCE       0.223    
 LARGEST SCALED RESIDUAL FORCE     -1.993E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.993E-07
 LARGEST INCREMENT OF DISP.         1.459E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.365E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.245E-04   TIME AVG. HEAT FLUX   1.564E-04
 LARGEST RESIDUAL HEAT FLUX         5.497E-08   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.554E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.094E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.144    
 STEP TIME COMPLETED       1.440E+06,  TOTAL TIME COMPLETED        1.440E+06


  INCREMENT   145 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.426       TIME AVG. FORCE       0.224    
 LARGEST SCALED RESIDUAL FORCE     -1.001E-09   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.001E-09
 LARGEST INCREMENT OF DISP.         1.473E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.367E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.270E-04   TIME AVG. HEAT FLUX   1.576E-04
 LARGEST RESIDUAL HEAT FLUX        -1.480E-15   AT NODE       1662   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.556E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.732E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.145    
 STEP TIME COMPLETED       1.450E+06,  TOTAL TIME COMPLETED        1.450E+06


  INCREMENT   146 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.428       TIME AVG. FORCE       0.225    
 LARGEST SCALED RESIDUAL FORCE     -2.106E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.106E-07
 LARGEST INCREMENT OF DISP.         1.487E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.384E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.295E-04   TIME AVG. HEAT FLUX   1.588E-04
 LARGEST RESIDUAL HEAT FLUX         3.020E-08   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.558E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.652E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.146    
 STEP TIME COMPLETED       1.460E+06,  TOTAL TIME COMPLETED        1.460E+06


  INCREMENT   147 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.431       TIME AVG. FORCE       0.227    
 LARGEST SCALED RESIDUAL FORCE     -1.061E-09   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.061E-09
 LARGEST INCREMENT OF DISP.         1.501E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.385E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.320E-04   TIME AVG. HEAT FLUX   1.600E-04
 LARGEST RESIDUAL HEAT FLUX         1.401E-15   AT NODE       4106   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.559E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.265E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.147    
 STEP TIME COMPLETED       1.470E+06,  TOTAL TIME COMPLETED        1.470E+06


  INCREMENT   148 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.433       TIME AVG. FORCE       0.228    
 LARGEST SCALED RESIDUAL FORCE     -2.225E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.225E-07
 LARGEST INCREMENT OF DISP.         1.515E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.403E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.344E-04   TIME AVG. HEAT FLUX   1.612E-04
 LARGEST RESIDUAL HEAT FLUX        -2.753E-08   AT NODE       6706   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.560E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.191E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.148    
 STEP TIME COMPLETED       1.480E+06,  TOTAL TIME COMPLETED        1.480E+06


  INCREMENT   149 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.435       TIME AVG. FORCE       0.230    
 LARGEST SCALED RESIDUAL FORCE     -1.122E-09   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.122E-09
 LARGEST INCREMENT OF DISP.         1.529E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.404E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.369E-04   TIME AVG. HEAT FLUX   1.624E-04
 LARGEST RESIDUAL HEAT FLUX         1.461E-15   AT NODE       5003   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.561E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.779E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.149    
 STEP TIME COMPLETED       1.490E+06,  TOTAL TIME COMPLETED        1.490E+06


  INCREMENT   150 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.438       TIME AVG. FORCE       0.231    
 LARGEST SCALED RESIDUAL FORCE     -2.352E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.352E-07
 LARGEST INCREMENT OF DISP.         1.543E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.422E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.394E-04   TIME AVG. HEAT FLUX   1.635E-04
 LARGEST RESIDUAL HEAT FLUX        -2.759E-08   AT NODE       6706   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.562E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.112E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.150    
 STEP TIME COMPLETED       1.500E+06,  TOTAL TIME COMPLETED        1.500E+06


  INCREMENT   151 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  331ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.440       TIME AVG. FORCE       0.232    
 LARGEST SCALED RESIDUAL FORCE     -1.186E-09   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.186E-09
 LARGEST INCREMENT OF DISP.         1.557E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.424E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.419E-04   TIME AVG. HEAT FLUX   1.647E-04
 LARGEST RESIDUAL HEAT FLUX         2.470E-15   AT NODE       6705   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.562E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.302E-06   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.151    
 STEP TIME COMPLETED       1.510E+06,  TOTAL TIME COMPLETED        1.510E+06


  INCREMENT   152 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.442       TIME AVG. FORCE       0.234    
 LARGEST SCALED RESIDUAL FORCE     -2.486E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.486E-07
 LARGEST INCREMENT OF DISP.         1.572E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.443E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.446E-04   TIME AVG. HEAT FLUX   1.659E-04
 LARGEST RESIDUAL HEAT FLUX        -4.476E-08   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.562E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.405E-06   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.152    
 STEP TIME COMPLETED       1.520E+06,  TOTAL TIME COMPLETED        1.520E+06


  INCREMENT   153 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  371ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.445       TIME AVG. FORCE       0.235    
 LARGEST SCALED RESIDUAL FORCE     -1.251E-09   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.251E-09
 LARGEST INCREMENT OF DISP.         1.586E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.444E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.470E-04   TIME AVG. HEAT FLUX   1.671E-04
 LARGEST RESIDUAL HEAT FLUX         1.325E-15   AT NODE        737   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.562E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.419E-06   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.153    
 STEP TIME COMPLETED       1.530E+06,  TOTAL TIME COMPLETED        1.530E+06


  INCREMENT   154 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.447       TIME AVG. FORCE       0.236    
 LARGEST SCALED RESIDUAL FORCE     -2.628E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.628E-07
 LARGEST INCREMENT OF DISP.         1.601E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.464E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.495E-04   TIME AVG. HEAT FLUX   1.683E-04
 LARGEST RESIDUAL HEAT FLUX        -6.976E-08   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.561E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.504E-06   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.154    
 STEP TIME COMPLETED       1.540E+06,  TOTAL TIME COMPLETED        1.540E+06


  INCREMENT   155 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  338ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.449       TIME AVG. FORCE       0.238    
 LARGEST SCALED RESIDUAL FORCE     -1.319E-09   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.319E-09
 LARGEST INCREMENT OF DISP.         1.616E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.465E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.519E-04   TIME AVG. HEAT FLUX   1.695E-04
 LARGEST RESIDUAL HEAT FLUX        -1.368E-15   AT NODE       4372   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.561E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.570E-06   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.155    
 STEP TIME COMPLETED       1.550E+06,  TOTAL TIME COMPLETED        1.550E+06


  INCREMENT   156 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  16ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.452       TIME AVG. FORCE       0.239    
 LARGEST SCALED RESIDUAL FORCE     -2.779E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.779E-07
 LARGEST INCREMENT OF DISP.         1.630E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.486E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.543E-04   TIME AVG. HEAT FLUX   1.707E-04
 LARGEST RESIDUAL HEAT FLUX        -9.468E-08   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.560E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.635E-06   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.156    
 STEP TIME COMPLETED       1.560E+06,  TOTAL TIME COMPLETED        1.560E+06


  INCREMENT   157 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.454       TIME AVG. FORCE       0.241    
 LARGEST SCALED RESIDUAL FORCE     -1.389E-09   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.389E-09
 LARGEST INCREMENT OF DISP.         1.645E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.487E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.569E-04   TIME AVG. HEAT FLUX   1.719E-04
 LARGEST RESIDUAL HEAT FLUX        -2.033E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.558E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.376E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.157    
 STEP TIME COMPLETED       1.570E+06,  TOTAL TIME COMPLETED        1.570E+06


  INCREMENT   158 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.456       TIME AVG. FORCE       0.242    
 LARGEST SCALED RESIDUAL FORCE     -2.940E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.940E-07
 LARGEST INCREMENT OF DISP.         1.660E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.509E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.594E-04   TIME AVG. HEAT FLUX   1.731E-04
 LARGEST RESIDUAL HEAT FLUX        -1.194E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.557E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.408E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.158    
 STEP TIME COMPLETED       1.580E+06,  TOTAL TIME COMPLETED        1.580E+06


  INCREMENT   159 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.459       TIME AVG. FORCE       0.243    
 LARGEST SCALED RESIDUAL FORCE     -1.471E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.471E-09
 LARGEST INCREMENT OF DISP.         1.675E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.510E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.619E-04   TIME AVG. HEAT FLUX   1.743E-04
 LARGEST RESIDUAL HEAT FLUX         1.998E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.555E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.969E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.159    
 STEP TIME COMPLETED       1.590E+06,  TOTAL TIME COMPLETED        1.590E+06


  INCREMENT   160 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  20ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.461       TIME AVG. FORCE       0.245    
 LARGEST SCALED RESIDUAL FORCE     -3.110E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.110E-07
 LARGEST INCREMENT OF DISP.         1.691E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.533E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.644E-04   TIME AVG. HEAT FLUX   1.755E-04
 LARGEST RESIDUAL HEAT FLUX        -1.440E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.553E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.989E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.160    
 STEP TIME COMPLETED       1.600E+06,  TOTAL TIME COMPLETED        1.600E+06


  INCREMENT   161 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.463       TIME AVG. FORCE       0.246    
 LARGEST SCALED RESIDUAL FORCE     -1.688E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.688E-09
 LARGEST INCREMENT OF DISP.         1.706E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.534E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.669E-04   TIME AVG. HEAT FLUX   1.767E-04
 LARGEST RESIDUAL HEAT FLUX         2.200E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.550E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.583E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.161    
 STEP TIME COMPLETED       1.610E+06,  TOTAL TIME COMPLETED        1.610E+06


  INCREMENT   162 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.465       TIME AVG. FORCE       0.247    
 LARGEST SCALED RESIDUAL FORCE     -3.292E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.292E-07
 LARGEST INCREMENT OF DISP.         1.722E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.559E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.694E-04   TIME AVG. HEAT FLUX   1.779E-04
 LARGEST RESIDUAL HEAT FLUX        -1.682E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.548E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.591E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.162    
 STEP TIME COMPLETED       1.620E+06,  TOTAL TIME COMPLETED        1.620E+06


  INCREMENT   163 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.467       TIME AVG. FORCE       0.249    
 LARGEST SCALED RESIDUAL FORCE     -1.917E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.917E-09
 LARGEST INCREMENT OF DISP.         1.737E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.559E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.719E-04   TIME AVG. HEAT FLUX   1.791E-04
 LARGEST RESIDUAL HEAT FLUX        -1.776E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.545E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.221E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.163    
 STEP TIME COMPLETED       1.630E+06,  TOTAL TIME COMPLETED        1.630E+06


  INCREMENT   164 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  20ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.470       TIME AVG. FORCE       0.250    
 LARGEST SCALED RESIDUAL FORCE     -3.485E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.485E-07
 LARGEST INCREMENT OF DISP.         1.753E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.585E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.744E-04   TIME AVG. HEAT FLUX   1.803E-04
 LARGEST RESIDUAL HEAT FLUX        -1.920E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.541E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.215E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.164    
 STEP TIME COMPLETED       1.640E+06,  TOTAL TIME COMPLETED        1.640E+06


  INCREMENT   165 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  348ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.472       TIME AVG. FORCE       0.251    
 LARGEST SCALED RESIDUAL FORCE     -2.157E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.157E-09
 LARGEST INCREMENT OF DISP.         1.769E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.585E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.768E-04   TIME AVG. HEAT FLUX   1.815E-04
 LARGEST RESIDUAL HEAT FLUX        -1.848E-15   AT NODE       2575   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.538E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.881E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.165    
 STEP TIME COMPLETED       1.650E+06,  TOTAL TIME COMPLETED        1.650E+06


  INCREMENT   166 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.474       TIME AVG. FORCE       0.253    
 LARGEST SCALED RESIDUAL FORCE     -3.691E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.691E-07
 LARGEST INCREMENT OF DISP.         1.785E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.613E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.792E-04   TIME AVG. HEAT FLUX   1.827E-04
 LARGEST RESIDUAL HEAT FLUX        -2.154E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.534E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.861E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.166    
 STEP TIME COMPLETED       1.660E+06,  TOTAL TIME COMPLETED        1.660E+06


  INCREMENT   167 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.476       TIME AVG. FORCE       0.254    
 LARGEST SCALED RESIDUAL FORCE     -2.410E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.410E-09
 LARGEST INCREMENT OF DISP.         1.801E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.613E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.817E-04   TIME AVG. HEAT FLUX   1.839E-04
 LARGEST RESIDUAL HEAT FLUX        -1.360E-15   AT NODE       5212   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.529E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.565E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.167    
 STEP TIME COMPLETED       1.670E+06,  TOTAL TIME COMPLETED        1.670E+06


  INCREMENT   168 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.478       TIME AVG. FORCE       0.255    
 LARGEST SCALED RESIDUAL FORCE     -3.911E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.911E-07
 LARGEST INCREMENT OF DISP.         1.818E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.641E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.841E-04   TIME AVG. HEAT FLUX   1.851E-04
 LARGEST RESIDUAL HEAT FLUX        -2.382E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.525E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.529E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.168    
 STEP TIME COMPLETED       1.680E+06,  TOTAL TIME COMPLETED        1.680E+06


  INCREMENT   169 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  352ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.480       TIME AVG. FORCE       0.257    
 LARGEST SCALED RESIDUAL FORCE     -2.677E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.677E-09
 LARGEST INCREMENT OF DISP.         1.834E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.641E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.865E-04   TIME AVG. HEAT FLUX   1.863E-04
 LARGEST RESIDUAL HEAT FLUX         1.804E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.519E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.273E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.169    
 STEP TIME COMPLETED       1.690E+06,  TOTAL TIME COMPLETED        1.690E+06


  INCREMENT   170 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  343ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.483       TIME AVG. FORCE       0.258    
 LARGEST SCALED RESIDUAL FORCE     -2.610E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.610E-09
 LARGEST INCREMENT OF DISP.         1.851E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.670E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.889E-04   TIME AVG. HEAT FLUX   1.875E-04
 LARGEST RESIDUAL HEAT FLUX        -1.610E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.514E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.233E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.170    
 STEP TIME COMPLETED       1.700E+06,  TOTAL TIME COMPLETED        1.700E+06


  INCREMENT   171 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  344ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.485       TIME AVG. FORCE       0.259    
 LARGEST SCALED RESIDUAL FORCE     -2.914E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.914E-09
 LARGEST INCREMENT OF DISP.         1.868E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.677E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.913E-04   TIME AVG. HEAT FLUX   1.887E-04
 LARGEST RESIDUAL HEAT FLUX        -1.371E-15   AT NODE       3379   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.508E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.913E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.171    
 STEP TIME COMPLETED       1.710E+06,  TOTAL TIME COMPLETED        1.710E+06


  INCREMENT   172 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  347ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.487       TIME AVG. FORCE       0.261    
 LARGEST SCALED RESIDUAL FORCE     -2.965E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.965E-09
 LARGEST INCREMENT OF DISP.         1.885E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.695E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.938E-04   TIME AVG. HEAT FLUX   1.899E-04
 LARGEST RESIDUAL HEAT FLUX         1.479E-15   AT NODE       4816   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.502E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.112E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.172    
 STEP TIME COMPLETED       1.720E+06,  TOTAL TIME COMPLETED        1.720E+06


  INCREMENT   173 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  329ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.489       TIME AVG. FORCE       0.262    
 LARGEST SCALED RESIDUAL FORCE     -3.141E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.141E-09
 LARGEST INCREMENT OF DISP.         1.902E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.711E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.962E-04   TIME AVG. HEAT FLUX   1.911E-04
 LARGEST RESIDUAL HEAT FLUX         3.303E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.496E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.544E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.173    
 STEP TIME COMPLETED       1.730E+06,  TOTAL TIME COMPLETED        1.730E+06


  INCREMENT   174 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.491       TIME AVG. FORCE       0.263    
 LARGEST SCALED RESIDUAL FORCE     -3.279E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.279E-09
 LARGEST INCREMENT OF DISP.         1.919E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.727E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.985E-04   TIME AVG. HEAT FLUX   1.922E-04
 LARGEST RESIDUAL HEAT FLUX         1.679E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.489E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.906E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.174    
 STEP TIME COMPLETED       1.740E+06,  TOTAL TIME COMPLETED        1.740E+06


  INCREMENT   175 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.493       TIME AVG. FORCE       0.265    
 LARGEST SCALED RESIDUAL FORCE     -3.432E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.432E-09
 LARGEST INCREMENT OF DISP.         1.936E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.744E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.010E-04   TIME AVG. HEAT FLUX   1.934E-04
 LARGEST RESIDUAL HEAT FLUX        -2.054E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.481E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.295E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.175    
 STEP TIME COMPLETED       1.750E+06,  TOTAL TIME COMPLETED        1.750E+06


  INCREMENT   176 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  356ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.495       TIME AVG. FORCE       0.266    
 LARGEST SCALED RESIDUAL FORCE     -3.587E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.587E-09
 LARGEST INCREMENT OF DISP.         1.954E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.762E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.033E-04   TIME AVG. HEAT FLUX   1.946E-04
 LARGEST RESIDUAL HEAT FLUX         2.415E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.474E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.686E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.176    
 STEP TIME COMPLETED       1.760E+06,  TOTAL TIME COMPLETED        1.760E+06


  INCREMENT   177 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.497       TIME AVG. FORCE       0.267    
 LARGEST SCALED RESIDUAL FORCE     -3.746E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.746E-09
 LARGEST INCREMENT OF DISP.         1.972E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.779E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.058E-04   TIME AVG. HEAT FLUX   1.958E-04
 LARGEST RESIDUAL HEAT FLUX        -1.729E-15   AT NODE        201   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.466E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.084E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.177    
 STEP TIME COMPLETED       1.770E+06,  TOTAL TIME COMPLETED        1.770E+06


  INCREMENT   178 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.499       TIME AVG. FORCE       0.269    
 LARGEST SCALED RESIDUAL FORCE     -3.910E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.910E-09
 LARGEST INCREMENT OF DISP.         1.990E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.797E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.082E-04   TIME AVG. HEAT FLUX   1.970E-04
 LARGEST RESIDUAL HEAT FLUX        -1.556E-15   AT NODE        747   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.457E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.488E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.178    
 STEP TIME COMPLETED       1.780E+06,  TOTAL TIME COMPLETED        1.780E+06


  INCREMENT   179 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.501       TIME AVG. FORCE       0.270    
 LARGEST SCALED RESIDUAL FORCE     -4.077E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.077E-09
 LARGEST INCREMENT OF DISP.         2.008E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.816E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.105E-04   TIME AVG. HEAT FLUX   1.982E-04
 LARGEST RESIDUAL HEAT FLUX        -2.193E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.448E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.898E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.179    
 STEP TIME COMPLETED       1.790E+06,  TOTAL TIME COMPLETED        1.790E+06


  INCREMENT   180 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  368ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.503       TIME AVG. FORCE       0.271    
 LARGEST SCALED RESIDUAL FORCE     -4.248E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.248E-09
 LARGEST INCREMENT OF DISP.         2.026E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.835E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.128E-04   TIME AVG. HEAT FLUX   1.994E-04
 LARGEST RESIDUAL HEAT FLUX        -1.509E-15   AT NODE       4569   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.439E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.315E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.180    
 STEP TIME COMPLETED       1.800E+06,  TOTAL TIME COMPLETED        1.800E+06


  INCREMENT   181 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.505       TIME AVG. FORCE       0.273    
 LARGEST SCALED RESIDUAL FORCE     -4.423E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.423E-09
 LARGEST INCREMENT OF DISP.         2.045E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.855E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.153E-04   TIME AVG. HEAT FLUX   2.006E-04
 LARGEST RESIDUAL HEAT FLUX        -1.485E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.429E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.738E-05   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.181    
 STEP TIME COMPLETED       1.810E+06,  TOTAL TIME COMPLETED        1.810E+06


  INCREMENT   182 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.507       TIME AVG. FORCE       0.274    
 LARGEST SCALED RESIDUAL FORCE     -4.603E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.603E-09
 LARGEST INCREMENT OF DISP.         2.064E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.875E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.176E-04   TIME AVG. HEAT FLUX   2.018E-04
 LARGEST RESIDUAL HEAT FLUX         1.440E-15   AT NODE       3730   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.419E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.017E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.182    
 STEP TIME COMPLETED       1.820E+06,  TOTAL TIME COMPLETED        1.820E+06


  INCREMENT   183 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  377ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.509       TIME AVG. FORCE       0.275    
 LARGEST SCALED RESIDUAL FORCE     -4.788E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.788E-09
 LARGEST INCREMENT OF DISP.         2.082E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.895E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.199E-04   TIME AVG. HEAT FLUX   2.030E-04
 LARGEST RESIDUAL HEAT FLUX        -1.499E-15   AT NODE      11986   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.408E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.060E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.183    
 STEP TIME COMPLETED       1.830E+06,  TOTAL TIME COMPLETED        1.830E+06


  INCREMENT   184 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  373ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.511       TIME AVG. FORCE       0.276    
 LARGEST SCALED RESIDUAL FORCE     -4.977E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.977E-09
 LARGEST INCREMENT OF DISP.         2.102E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.916E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.223E-04   TIME AVG. HEAT FLUX   2.042E-04
 LARGEST RESIDUAL HEAT FLUX        -1.443E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.397E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.105E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.184    
 STEP TIME COMPLETED       1.840E+06,  TOTAL TIME COMPLETED        1.840E+06


  INCREMENT   185 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  385ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.513       TIME AVG. FORCE       0.278    
 LARGEST SCALED RESIDUAL FORCE     -5.171E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.171E-09
 LARGEST INCREMENT OF DISP.         2.121E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.938E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.247E-04   TIME AVG. HEAT FLUX   2.054E-04
 LARGEST RESIDUAL HEAT FLUX        -1.307E-15   AT NODE       4296   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.386E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.150E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.185    
 STEP TIME COMPLETED       1.850E+06,  TOTAL TIME COMPLETED        1.850E+06


  INCREMENT   186 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  379ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.515       TIME AVG. FORCE       0.279    
 LARGEST SCALED RESIDUAL FORCE     -5.370E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.370E-09
 LARGEST INCREMENT OF DISP.         2.141E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.960E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.270E-04   TIME AVG. HEAT FLUX   2.066E-04
 LARGEST RESIDUAL HEAT FLUX         1.323E-15   AT NODE       3261   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.374E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.195E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.186    
 STEP TIME COMPLETED       1.860E+06,  TOTAL TIME COMPLETED        1.860E+06


  INCREMENT   187 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.517       TIME AVG. FORCE       0.280    
 LARGEST SCALED RESIDUAL FORCE     -5.574E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.574E-09
 LARGEST INCREMENT OF DISP.         2.160E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.983E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.293E-04   TIME AVG. HEAT FLUX   2.078E-04
 LARGEST RESIDUAL HEAT FLUX         1.804E-15   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.361E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.241E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.187    
 STEP TIME COMPLETED       1.870E+06,  TOTAL TIME COMPLETED        1.870E+06


  INCREMENT   188 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  379ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.519       TIME AVG. FORCE       0.281    
 LARGEST SCALED RESIDUAL FORCE     -5.784E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.784E-09
 LARGEST INCREMENT OF DISP.         2.181E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.006E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.315E-04   TIME AVG. HEAT FLUX   2.090E-04
 LARGEST RESIDUAL HEAT FLUX        -1.494E-15   AT NODE        366   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.349E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.288E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.188    
 STEP TIME COMPLETED       1.880E+06,  TOTAL TIME COMPLETED        1.880E+06


  INCREMENT   189 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  380ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.521       TIME AVG. FORCE       0.283    
 LARGEST SCALED RESIDUAL FORCE     -5.998E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.998E-09
 LARGEST INCREMENT OF DISP.         2.201E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.030E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.336E-04   TIME AVG. HEAT FLUX   2.102E-04
 LARGEST RESIDUAL HEAT FLUX         1.401E-15   AT NODE        483   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.335E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.336E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.189    
 STEP TIME COMPLETED       1.890E+06,  TOTAL TIME COMPLETED        1.890E+06


  INCREMENT   190 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  408ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.523       TIME AVG. FORCE       0.284    
 LARGEST SCALED RESIDUAL FORCE     -6.219E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.219E-09
 LARGEST INCREMENT OF DISP.         2.221E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.055E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.360E-04   TIME AVG. HEAT FLUX   2.114E-04
 LARGEST RESIDUAL HEAT FLUX         1.411E-15   AT NODE       2575   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.321E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.384E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.190    
 STEP TIME COMPLETED       1.900E+06,  TOTAL TIME COMPLETED        1.900E+06


  INCREMENT   191 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  383ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.525       TIME AVG. FORCE       0.285    
 LARGEST SCALED RESIDUAL FORCE     -6.445E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.445E-09
 LARGEST INCREMENT OF DISP.         2.242E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.080E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.383E-04   TIME AVG. HEAT FLUX   2.126E-04
 LARGEST RESIDUAL HEAT FLUX        -1.846E-15   AT NODE       6434   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.307E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.433E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.191    
 STEP TIME COMPLETED       1.910E+06,  TOTAL TIME COMPLETED        1.910E+06


  INCREMENT   192 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  379ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.527       TIME AVG. FORCE       0.287    
 LARGEST SCALED RESIDUAL FORCE     -6.677E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.677E-09
 LARGEST INCREMENT OF DISP.         2.263E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.107E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.406E-04   TIME AVG. HEAT FLUX   2.138E-04
 LARGEST RESIDUAL HEAT FLUX         1.707E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.292E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.483E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.192    
 STEP TIME COMPLETED       1.920E+06,  TOTAL TIME COMPLETED        1.920E+06


  INCREMENT   193 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.528       TIME AVG. FORCE       0.288    
 LARGEST SCALED RESIDUAL FORCE     -6.915E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.915E-09
 LARGEST INCREMENT OF DISP.         2.285E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.133E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.426E-04   TIME AVG. HEAT FLUX   2.150E-04
 LARGEST RESIDUAL HEAT FLUX         2.942E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.277E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.534E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.193    
 STEP TIME COMPLETED       1.930E+06,  TOTAL TIME COMPLETED        1.930E+06


  INCREMENT   194 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  380ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.530       TIME AVG. FORCE       0.289    
 LARGEST SCALED RESIDUAL FORCE     -7.160E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.160E-09
 LARGEST INCREMENT OF DISP.         2.306E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.161E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.447E-04   TIME AVG. HEAT FLUX   2.162E-04
 LARGEST RESIDUAL HEAT FLUX        -1.589E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.261E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.585E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.194    
 STEP TIME COMPLETED       1.940E+06,  TOTAL TIME COMPLETED        1.940E+06


  INCREMENT   195 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.532       TIME AVG. FORCE       0.290    
 LARGEST SCALED RESIDUAL FORCE     -7.411E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.411E-09
 LARGEST INCREMENT OF DISP.         2.328E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.189E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.469E-04   TIME AVG. HEAT FLUX   2.174E-04
 LARGEST RESIDUAL HEAT FLUX         1.832E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.245E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.637E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.195    
 STEP TIME COMPLETED       1.950E+06,  TOTAL TIME COMPLETED        1.950E+06


  INCREMENT   196 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  356ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.534       TIME AVG. FORCE       0.292    
 LARGEST SCALED RESIDUAL FORCE     -7.669E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.669E-09
 LARGEST INCREMENT OF DISP.         2.350E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.219E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.490E-04   TIME AVG. HEAT FLUX   2.186E-04
 LARGEST RESIDUAL HEAT FLUX        -1.887E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.228E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.690E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.196    
 STEP TIME COMPLETED       1.960E+06,  TOTAL TIME COMPLETED        1.960E+06


  INCREMENT   197 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  355ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.536       TIME AVG. FORCE       0.293    
 LARGEST SCALED RESIDUAL FORCE     -7.934E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.934E-09
 LARGEST INCREMENT OF DISP.         2.373E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.249E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.511E-04   TIME AVG. HEAT FLUX   2.198E-04
 LARGEST RESIDUAL HEAT FLUX         1.700E-15   AT NODE       5298   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.210E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.743E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.197    
 STEP TIME COMPLETED       1.970E+06,  TOTAL TIME COMPLETED        1.970E+06


  INCREMENT   198 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.538       TIME AVG. FORCE       0.294    
 LARGEST SCALED RESIDUAL FORCE     -8.206E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.206E-09
 LARGEST INCREMENT OF DISP.         2.396E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.280E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.534E-04   TIME AVG. HEAT FLUX   2.209E-04
 LARGEST RESIDUAL HEAT FLUX        -1.382E-15   AT NODE       5184   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.192E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.797E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.198    
 STEP TIME COMPLETED       1.980E+06,  TOTAL TIME COMPLETED        1.980E+06


  INCREMENT   199 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.539       TIME AVG. FORCE       0.295    
 LARGEST SCALED RESIDUAL FORCE     -8.486E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.486E-09
 LARGEST INCREMENT OF DISP.         2.419E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.311E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.555E-04   TIME AVG. HEAT FLUX   2.221E-04
 LARGEST RESIDUAL HEAT FLUX        -2.658E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.174E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.852E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.199    
 STEP TIME COMPLETED       1.990E+06,  TOTAL TIME COMPLETED        1.990E+06


  INCREMENT   200 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.541       TIME AVG. FORCE       0.296    
 LARGEST SCALED RESIDUAL FORCE     -8.774E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.774E-09
 LARGEST INCREMENT OF DISP.         2.442E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.344E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.576E-04   TIME AVG. HEAT FLUX   2.233E-04
 LARGEST RESIDUAL HEAT FLUX         2.144E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.155E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.908E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.200    
 STEP TIME COMPLETED       2.000E+06,  TOTAL TIME COMPLETED        2.000E+06


  INCREMENT   201 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  345ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.543       TIME AVG. FORCE       0.298    
 LARGEST SCALED RESIDUAL FORCE     -9.070E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.070E-09
 LARGEST INCREMENT OF DISP.         2.466E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.378E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.597E-04   TIME AVG. HEAT FLUX   2.245E-04
 LARGEST RESIDUAL HEAT FLUX        -1.543E-15   AT NODE       1996   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.135E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.964E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.201    
 STEP TIME COMPLETED       2.010E+06,  TOTAL TIME COMPLETED        2.010E+06


  INCREMENT   202 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  353ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.545       TIME AVG. FORCE       0.299    
 LARGEST SCALED RESIDUAL FORCE     -9.375E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.375E-09
 LARGEST INCREMENT OF DISP.         2.490E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.413E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.618E-04   TIME AVG. HEAT FLUX   2.257E-04
 LARGEST RESIDUAL HEAT FLUX         1.700E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.115E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.022E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.202    
 STEP TIME COMPLETED       2.020E+06,  TOTAL TIME COMPLETED        2.020E+06


  INCREMENT   203 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  349ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.546       TIME AVG. FORCE       0.300    
 LARGEST SCALED RESIDUAL FORCE     -9.688E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.688E-09
 LARGEST INCREMENT OF DISP.         2.515E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.449E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.637E-04   TIME AVG. HEAT FLUX   2.268E-04
 LARGEST RESIDUAL HEAT FLUX        -1.825E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.094E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.080E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.203    
 STEP TIME COMPLETED       2.030E+06,  TOTAL TIME COMPLETED        2.030E+06


  INCREMENT   204 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.548       TIME AVG. FORCE       0.301    
 LARGEST SCALED RESIDUAL FORCE     -1.001E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.001E-08
 LARGEST INCREMENT OF DISP.         2.539E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.486E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.658E-04   TIME AVG. HEAT FLUX   2.280E-04
 LARGEST RESIDUAL HEAT FLUX         1.441E-15   AT NODE        505   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.073E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.139E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.204    
 STEP TIME COMPLETED       2.040E+06,  TOTAL TIME COMPLETED        2.040E+06


  INCREMENT   205 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.550       TIME AVG. FORCE       0.303    
 LARGEST SCALED RESIDUAL FORCE     -1.034E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.034E-08
 LARGEST INCREMENT OF DISP.         2.565E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.524E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.679E-04   TIME AVG. HEAT FLUX   2.292E-04
 LARGEST RESIDUAL HEAT FLUX         1.665E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.051E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.199E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.205    
 STEP TIME COMPLETED       2.050E+06,  TOTAL TIME COMPLETED        2.050E+06


  INCREMENT   206 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.552       TIME AVG. FORCE       0.304    
 LARGEST SCALED RESIDUAL FORCE     -1.068E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.068E-08
 LARGEST INCREMENT OF DISP.         2.590E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.564E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.699E-04   TIME AVG. HEAT FLUX   2.304E-04
 LARGEST RESIDUAL HEAT FLUX        -1.368E-15   AT NODE       1714   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.028E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.260E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.206    
 STEP TIME COMPLETED       2.060E+06,  TOTAL TIME COMPLETED        2.060E+06


  INCREMENT   207 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.553       TIME AVG. FORCE       0.305    
 LARGEST SCALED RESIDUAL FORCE     -1.104E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.104E-08
 LARGEST INCREMENT OF DISP.         2.616E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.604E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.719E-04   TIME AVG. HEAT FLUX   2.315E-04
 LARGEST RESIDUAL HEAT FLUX        -1.558E-15   AT NODE       2093   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.005E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.321E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.207    
 STEP TIME COMPLETED       2.070E+06,  TOTAL TIME COMPLETED        2.070E+06


  INCREMENT   208 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.555       TIME AVG. FORCE       0.306    
 LARGEST SCALED RESIDUAL FORCE     -1.140E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.140E-08
 LARGEST INCREMENT OF DISP.         2.643E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.647E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.666E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.739E-04   TIME AVG. HEAT FLUX   2.327E-04
 LARGEST RESIDUAL HEAT FLUX         1.776E-15   AT NODE       6434   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.810E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.384E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.208    
 STEP TIME COMPLETED       2.080E+06,  TOTAL TIME COMPLETED        2.080E+06


  INCREMENT   209 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.557       TIME AVG. FORCE       0.307    
 LARGEST SCALED RESIDUAL FORCE     -1.177E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.177E-08
 LARGEST INCREMENT OF DISP.         2.670E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.690E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.705E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.759E-04   TIME AVG. HEAT FLUX   2.339E-04
 LARGEST RESIDUAL HEAT FLUX        -1.517E-15   AT NODE       5109   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.566E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.447E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.209    
 STEP TIME COMPLETED       2.090E+06,  TOTAL TIME COMPLETED        2.090E+06


  INCREMENT   210 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.558       TIME AVG. FORCE       0.309    
 LARGEST SCALED RESIDUAL FORCE     -1.216E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.216E-08
 LARGEST INCREMENT OF DISP.         2.697E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.735E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.745E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.776E-04   TIME AVG. HEAT FLUX   2.351E-04
 LARGEST RESIDUAL HEAT FLUX        -1.438E-15   AT NODE        806   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.315E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.511E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.210    
 STEP TIME COMPLETED       2.100E+06,  TOTAL TIME COMPLETED        2.100E+06


  INCREMENT   211 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.560       TIME AVG. FORCE       0.310    
 LARGEST SCALED RESIDUAL FORCE     -1.256E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.256E-08
 LARGEST INCREMENT OF DISP.         2.725E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.781E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.786E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.796E-04   TIME AVG. HEAT FLUX   2.362E-04
 LARGEST RESIDUAL HEAT FLUX         1.547E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.057E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.577E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.211    
 STEP TIME COMPLETED       2.110E+06,  TOTAL TIME COMPLETED        2.110E+06


  INCREMENT   212 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.561       TIME AVG. FORCE       0.311    
 LARGEST SCALED RESIDUAL FORCE     -1.297E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.297E-08
 LARGEST INCREMENT OF DISP.         2.753E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.830E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.828E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.814E-04   TIME AVG. HEAT FLUX   2.374E-04
 LARGEST RESIDUAL HEAT FLUX         1.603E-15   AT NODE        179   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.793E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.643E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.212    
 STEP TIME COMPLETED       2.120E+06,  TOTAL TIME COMPLETED        2.120E+06


  INCREMENT   213 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.563       TIME AVG. FORCE       0.312    
 LARGEST SCALED RESIDUAL FORCE     -1.339E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.339E-08
 LARGEST INCREMENT OF DISP.         2.782E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.879E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.871E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.832E-04   TIME AVG. HEAT FLUX   2.385E-04
 LARGEST RESIDUAL HEAT FLUX         2.755E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.522E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.710E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.213    
 STEP TIME COMPLETED       2.130E+06,  TOTAL TIME COMPLETED        2.130E+06


  INCREMENT   214 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.564       TIME AVG. FORCE       0.313    
 LARGEST SCALED RESIDUAL FORCE     -1.383E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.383E-08
 LARGEST INCREMENT OF DISP.         2.811E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.931E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.915E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.850E-04   TIME AVG. HEAT FLUX   2.397E-04
 LARGEST RESIDUAL HEAT FLUX         1.922E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.244E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.779E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.214    
 STEP TIME COMPLETED       2.140E+06,  TOTAL TIME COMPLETED        2.140E+06


  INCREMENT   215 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  355ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.566       TIME AVG. FORCE       0.314    
 LARGEST SCALED RESIDUAL FORCE     -1.428E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.428E-08
 LARGEST INCREMENT OF DISP.         2.841E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.984E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.960E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.870E-04   TIME AVG. HEAT FLUX   2.408E-04
 LARGEST RESIDUAL HEAT FLUX        -1.293E-15   AT NODE        360   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.959E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.848E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.215    
 STEP TIME COMPLETED       2.150E+06,  TOTAL TIME COMPLETED        2.150E+06


  INCREMENT   216 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.567       TIME AVG. FORCE       0.316    
 LARGEST SCALED RESIDUAL FORCE     -1.475E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.475E-08
 LARGEST INCREMENT OF DISP.         2.871E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.039E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.006E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.887E-04   TIME AVG. HEAT FLUX   2.420E-04
 LARGEST RESIDUAL HEAT FLUX        -2.269E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.667E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.918E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.216    
 STEP TIME COMPLETED       2.160E+06,  TOTAL TIME COMPLETED        2.160E+06


  INCREMENT   217 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  355ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.569       TIME AVG. FORCE       0.317    
 LARGEST SCALED RESIDUAL FORCE     -1.524E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.524E-08
 LARGEST INCREMENT OF DISP.         2.902E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.096E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.054E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.903E-04   TIME AVG. HEAT FLUX   2.432E-04
 LARGEST RESIDUAL HEAT FLUX        -1.628E-15   AT NODE       1168   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.368E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.990E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.217    
 STEP TIME COMPLETED       2.170E+06,  TOTAL TIME COMPLETED        2.170E+06


  INCREMENT   218 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  353ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.570       TIME AVG. FORCE       0.318    
 LARGEST SCALED RESIDUAL FORCE     -1.574E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.574E-08
 LARGEST INCREMENT OF DISP.         2.934E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.156E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.102E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.920E-04   TIME AVG. HEAT FLUX   2.443E-04
 LARGEST RESIDUAL HEAT FLUX        -1.303E-15   AT NODE        341   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.062E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.062E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.218    
 STEP TIME COMPLETED       2.180E+06,  TOTAL TIME COMPLETED        2.180E+06


  INCREMENT   219 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.572       TIME AVG. FORCE       0.319    
 LARGEST SCALED RESIDUAL FORCE     -1.626E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.626E-08
 LARGEST INCREMENT OF DISP.         2.966E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.220E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.154E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.940E-04   TIME AVG. HEAT FLUX   2.454E-04
 LARGEST RESIDUAL HEAT FLUX        -1.998E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.748E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.136E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.219    
 STEP TIME COMPLETED       2.190E+06,  TOTAL TIME COMPLETED        2.190E+06


  INCREMENT   220 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  346ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.573       TIME AVG. FORCE       0.320    
 LARGEST SCALED RESIDUAL FORCE     -1.682E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.682E-08
 LARGEST INCREMENT OF DISP.         2.999E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.283E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.207E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.956E-04   TIME AVG. HEAT FLUX   2.466E-04
 LARGEST RESIDUAL HEAT FLUX         2.498E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.427E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.215E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.220    
 STEP TIME COMPLETED       2.200E+06,  TOTAL TIME COMPLETED        2.200E+06


  INCREMENT   221 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.575       TIME AVG. FORCE       0.321    
 LARGEST SCALED RESIDUAL FORCE     -1.735E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.735E-08
 LARGEST INCREMENT OF DISP.         3.032E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.349E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.256E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.971E-04   TIME AVG. HEAT FLUX   2.477E-04
 LARGEST RESIDUAL HEAT FLUX         1.374E-15   AT NODE       2716   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.098E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.288E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.221    
 STEP TIME COMPLETED       2.210E+06,  TOTAL TIME COMPLETED        2.210E+06


  INCREMENT   222 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  355ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.576       TIME AVG. FORCE       0.323    
 LARGEST SCALED RESIDUAL FORCE     -1.793E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.793E-08
 LARGEST INCREMENT OF DISP.         3.067E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.417E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.311E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.988E-04   TIME AVG. HEAT FLUX   2.489E-04
 LARGEST RESIDUAL HEAT FLUX        -2.609E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.762E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.365E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.222    
 STEP TIME COMPLETED       2.220E+06,  TOTAL TIME COMPLETED        2.220E+06


  INCREMENT   223 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.577       TIME AVG. FORCE       0.324    
 LARGEST SCALED RESIDUAL FORCE     -1.853E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.853E-08
 LARGEST INCREMENT OF DISP.         3.102E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.489E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.367E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.006E-04   TIME AVG. HEAT FLUX   2.500E-04
 LARGEST RESIDUAL HEAT FLUX        -1.478E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.417E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.444E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.223    
 STEP TIME COMPLETED       2.230E+06,  TOTAL TIME COMPLETED        2.230E+06


  INCREMENT   224 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.579       TIME AVG. FORCE       0.325    
 LARGEST SCALED RESIDUAL FORCE     -1.915E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.915E-08
 LARGEST INCREMENT OF DISP.         3.137E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.563E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.424E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.023E-04   TIME AVG. HEAT FLUX   2.511E-04
 LARGEST RESIDUAL HEAT FLUX        -1.525E-15   AT NODE       5498   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.065E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.524E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.224    
 STEP TIME COMPLETED       2.240E+06,  TOTAL TIME COMPLETED        2.240E+06


  INCREMENT   225 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.580       TIME AVG. FORCE       0.326    
 LARGEST SCALED RESIDUAL FORCE     -1.980E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.980E-08
 LARGEST INCREMENT OF DISP.         3.174E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.641E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.484E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.036E-04   TIME AVG. HEAT FLUX   2.523E-04
 LARGEST RESIDUAL HEAT FLUX         1.420E-15   AT NODE       3910   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         4.704E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.606E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.225    
 STEP TIME COMPLETED       2.250E+06,  TOTAL TIME COMPLETED        2.250E+06


  INCREMENT   226 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.581       TIME AVG. FORCE       0.327    
 LARGEST SCALED RESIDUAL FORCE     -2.047E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.047E-08
 LARGEST INCREMENT OF DISP.         3.211E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.721E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.545E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.050E-04   TIME AVG. HEAT FLUX   2.534E-04
 LARGEST RESIDUAL HEAT FLUX        -1.443E-15   AT NODE      11984   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         4.335E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.689E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.226    
 STEP TIME COMPLETED       2.260E+06,  TOTAL TIME COMPLETED        2.260E+06


  INCREMENT   227 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.583       TIME AVG. FORCE       0.328    
 LARGEST SCALED RESIDUAL FORCE     -2.117E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.117E-08
 LARGEST INCREMENT OF DISP.         3.249E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.805E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.608E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.067E-04   TIME AVG. HEAT FLUX   2.545E-04
 LARGEST RESIDUAL HEAT FLUX        -1.561E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         3.958E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.773E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.227    
 STEP TIME COMPLETED       2.270E+06,  TOTAL TIME COMPLETED        2.270E+06


  INCREMENT   228 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.584       TIME AVG. FORCE       0.329    
 LARGEST SCALED RESIDUAL FORCE     -2.190E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.190E-08
 LARGEST INCREMENT OF DISP.         3.288E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.892E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.673E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.081E-04   TIME AVG. HEAT FLUX   2.556E-04
 LARGEST RESIDUAL HEAT FLUX         1.418E-15   AT NODE       3395   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         3.572E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.859E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.228    
 STEP TIME COMPLETED       2.280E+06,  TOTAL TIME COMPLETED        2.280E+06


  INCREMENT   229 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  371ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.585       TIME AVG. FORCE       0.330    
 LARGEST SCALED RESIDUAL FORCE     -2.266E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.266E-08
 LARGEST INCREMENT OF DISP.         3.328E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.984E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.741E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.096E-04   TIME AVG. HEAT FLUX   2.567E-04
 LARGEST RESIDUAL HEAT FLUX        -1.887E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         3.204E-03   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.947E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.229    
 STEP TIME COMPLETED       2.290E+06,  TOTAL TIME COMPLETED        2.290E+06


  INCREMENT   230 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.586       TIME AVG. FORCE       0.332    
 LARGEST SCALED RESIDUAL FORCE     -2.345E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.345E-08
 LARGEST INCREMENT OF DISP.         3.368E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.079E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.810E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.107E-04   TIME AVG. HEAT FLUX   2.578E-04
 LARGEST RESIDUAL HEAT FLUX         1.457E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         3.026E-03   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.036E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.230    
 STEP TIME COMPLETED       2.300E+06,  TOTAL TIME COMPLETED        2.300E+06


  INCREMENT   231 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.588       TIME AVG. FORCE       0.333    
 LARGEST SCALED RESIDUAL FORCE     -2.427E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.427E-08
 LARGEST INCREMENT OF DISP.         3.410E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.178E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.883E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.121E-04   TIME AVG. HEAT FLUX   2.589E-04
 LARGEST RESIDUAL HEAT FLUX        -1.943E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.843E-03   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.128E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.231    
 STEP TIME COMPLETED       2.310E+06,  TOTAL TIME COMPLETED        2.310E+06


  INCREMENT   232 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.589       TIME AVG. FORCE       0.334    
 LARGEST SCALED RESIDUAL FORCE     -2.513E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.513E-08
 LARGEST INCREMENT OF DISP.         3.453E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.281E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.958E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.134E-04   TIME AVG. HEAT FLUX   2.600E-04
 LARGEST RESIDUAL HEAT FLUX        -2.449E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.682E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.221E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.232    
 STEP TIME COMPLETED       2.320E+06,  TOTAL TIME COMPLETED        2.320E+06


  INCREMENT   233 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  380ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.590       TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -2.603E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.603E-08
 LARGEST INCREMENT OF DISP.         3.497E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.390E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.035E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.146E-04   TIME AVG. HEAT FLUX   2.611E-04
 LARGEST RESIDUAL HEAT FLUX         2.755E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.653E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.316E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.233    
 STEP TIME COMPLETED       2.330E+06,  TOTAL TIME COMPLETED        2.330E+06


  INCREMENT   234 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  377ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.591       TIME AVG. FORCE       0.336    
 LARGEST SCALED RESIDUAL FORCE     -2.697E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.697E-08
 LARGEST INCREMENT OF DISP.         3.542E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.503E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.116E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.160E-04   TIME AVG. HEAT FLUX   2.622E-04
 LARGEST RESIDUAL HEAT FLUX         1.554E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.623E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.413E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.234    
 STEP TIME COMPLETED       2.340E+06,  TOTAL TIME COMPLETED        2.340E+06


  INCREMENT   235 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.592       TIME AVG. FORCE       0.337    
 LARGEST SCALED RESIDUAL FORCE     -2.795E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.795E-08
 LARGEST INCREMENT OF DISP.         3.588E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.621E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.200E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.170E-04   TIME AVG. HEAT FLUX   2.633E-04
 LARGEST RESIDUAL HEAT FLUX        -1.943E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.593E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.512E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.235    
 STEP TIME COMPLETED       2.350E+06,  TOTAL TIME COMPLETED        2.350E+06


  INCREMENT   236 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.593       TIME AVG. FORCE       0.338    
 LARGEST SCALED RESIDUAL FORCE     -2.898E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.898E-08
 LARGEST INCREMENT OF DISP.         3.636E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.746E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.287E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.183E-04   TIME AVG. HEAT FLUX   2.644E-04
 LARGEST RESIDUAL HEAT FLUX        -3.081E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.562E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.613E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.236    
 STEP TIME COMPLETED       2.360E+06,  TOTAL TIME COMPLETED        2.360E+06


  INCREMENT   237 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.595       TIME AVG. FORCE       0.339    
 LARGEST SCALED RESIDUAL FORCE     -3.006E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.006E-08
 LARGEST INCREMENT OF DISP.         3.684E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.876E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.378E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.195E-04   TIME AVG. HEAT FLUX   2.655E-04
 LARGEST RESIDUAL HEAT FLUX        -1.475E-15   AT NODE        177   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.530E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.717E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.237    
 STEP TIME COMPLETED       2.370E+06,  TOTAL TIME COMPLETED        2.370E+06


  INCREMENT   238 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.596       TIME AVG. FORCE       0.340    
 LARGEST SCALED RESIDUAL FORCE     -3.119E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.119E-08
 LARGEST INCREMENT OF DISP.         3.734E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.012E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.472E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.204E-04   TIME AVG. HEAT FLUX   2.666E-04
 LARGEST RESIDUAL HEAT FLUX        -2.387E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.497E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.823E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.238    
 STEP TIME COMPLETED       2.380E+06,  TOTAL TIME COMPLETED        2.380E+06


  INCREMENT   239 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.597       TIME AVG. FORCE       0.341    
 LARGEST SCALED RESIDUAL FORCE     -3.237E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.237E-08
 LARGEST INCREMENT OF DISP.         3.786E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.155E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.571E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.215E-04   TIME AVG. HEAT FLUX   2.676E-04
 LARGEST RESIDUAL HEAT FLUX        -1.853E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.463E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.932E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.239    
 STEP TIME COMPLETED       2.390E+06,  TOTAL TIME COMPLETED        2.390E+06


  INCREMENT   240 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.598       TIME AVG. FORCE       0.342    
 LARGEST SCALED RESIDUAL FORCE     -3.362E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.362E-08
 LARGEST INCREMENT OF DISP.         3.839E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.305E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.674E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.223E-04   TIME AVG. HEAT FLUX   2.687E-04
 LARGEST RESIDUAL HEAT FLUX        -2.331E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.428E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.044E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.240    
 STEP TIME COMPLETED       2.400E+06,  TOTAL TIME COMPLETED        2.400E+06


  INCREMENT   241 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  412ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.599       TIME AVG. FORCE       0.344    
 LARGEST SCALED RESIDUAL FORCE     -3.494E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.494E-08
 LARGEST INCREMENT OF DISP.         3.894E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.463E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.781E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.235E-04   TIME AVG. HEAT FLUX   2.698E-04
 LARGEST RESIDUAL HEAT FLUX         1.665E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.393E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.158E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.241    
 STEP TIME COMPLETED       2.410E+06,  TOTAL TIME COMPLETED        2.410E+06


  INCREMENT   242 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.600       TIME AVG. FORCE       0.345    
 LARGEST SCALED RESIDUAL FORCE     -3.632E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.632E-08
 LARGEST INCREMENT OF DISP.         3.950E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.629E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.894E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.246E-04   TIME AVG. HEAT FLUX   2.708E-04
 LARGEST RESIDUAL HEAT FLUX         1.444E-15   AT NODE       2682   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.841E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.276E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.242    
 STEP TIME COMPLETED       2.420E+06,  TOTAL TIME COMPLETED        2.420E+06


  INCREMENT   243 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  352ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.601       TIME AVG. FORCE       0.346    
 LARGEST SCALED RESIDUAL FORCE     -3.778E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.778E-08
 LARGEST INCREMENT OF DISP.         4.008E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.804E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.011E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.257E-04   TIME AVG. HEAT FLUX   2.719E-04
 LARGEST RESIDUAL HEAT FLUX         2.887E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.381E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.397E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.243    
 STEP TIME COMPLETED       2.430E+06,  TOTAL TIME COMPLETED        2.430E+06


  INCREMENT   244 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  346ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.602       TIME AVG. FORCE       0.347    
 LARGEST SCALED RESIDUAL FORCE     -3.932E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.932E-08
 LARGEST INCREMENT OF DISP.         4.068E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.989E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.135E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.263E-04   TIME AVG. HEAT FLUX   2.729E-04
 LARGEST RESIDUAL HEAT FLUX        -2.387E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.933E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.521E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.244    
 STEP TIME COMPLETED       2.440E+06,  TOTAL TIME COMPLETED        2.440E+06


  INCREMENT   245 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  348ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.603       TIME AVG. FORCE       0.348    
 LARGEST SCALED RESIDUAL FORCE     -4.094E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.094E-08
 LARGEST INCREMENT OF DISP.         4.130E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.183E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.265E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.270E-04   TIME AVG. HEAT FLUX   2.740E-04
 LARGEST RESIDUAL HEAT FLUX         2.276E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.498E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.649E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.245    
 STEP TIME COMPLETED       2.450E+06,  TOTAL TIME COMPLETED        2.450E+06


  INCREMENT   246 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  350ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.604       TIME AVG. FORCE       0.349    
 LARGEST SCALED RESIDUAL FORCE     -4.266E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.266E-08
 LARGEST INCREMENT OF DISP.         4.194E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.388E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.401E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.278E-04   TIME AVG. HEAT FLUX   2.750E-04
 LARGEST RESIDUAL HEAT FLUX        -1.499E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.076E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.781E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.246    
 STEP TIME COMPLETED       2.460E+06,  TOTAL TIME COMPLETED        2.460E+06


  INCREMENT   247 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.605       TIME AVG. FORCE       0.350    
 LARGEST SCALED RESIDUAL FORCE     -4.449E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.449E-08
 LARGEST INCREMENT OF DISP.         4.260E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.605E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.544E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.287E-04   TIME AVG. HEAT FLUX   2.760E-04
 LARGEST RESIDUAL HEAT FLUX        -1.597E-15   AT NODE       1105   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.668E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.917E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.247    
 STEP TIME COMPLETED       2.470E+06,  TOTAL TIME COMPLETED        2.470E+06


  INCREMENT   248 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  351ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.606       TIME AVG. FORCE       0.351    
 LARGEST SCALED RESIDUAL FORCE     -4.642E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.642E-08
 LARGEST INCREMENT OF DISP.         4.328E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.834E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.695E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.295E-04   TIME AVG. HEAT FLUX   2.771E-04
 LARGEST RESIDUAL HEAT FLUX        -2.144E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -6.273E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.057E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.248    
 STEP TIME COMPLETED       2.480E+06,  TOTAL TIME COMPLETED        2.480E+06


  INCREMENT   249 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.607       TIME AVG. FORCE       0.352    
 LARGEST SCALED RESIDUAL FORCE     -4.848E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.848E-08
 LARGEST INCREMENT OF DISP.         4.399E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.077E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.854E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.298E-04   TIME AVG. HEAT FLUX   2.781E-04
 LARGEST RESIDUAL HEAT FLUX         2.623E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -6.894E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.202E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.249    
 STEP TIME COMPLETED       2.490E+06,  TOTAL TIME COMPLETED        2.490E+06


  INCREMENT   250 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  348ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.608       TIME AVG. FORCE       0.353    
 LARGEST SCALED RESIDUAL FORCE     -5.066E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.066E-08
 LARGEST INCREMENT OF DISP.         4.472E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.334E-05   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.022E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.307E-04   TIME AVG. HEAT FLUX   2.791E-04
 LARGEST RESIDUAL HEAT FLUX        -1.628E-15   AT NODE       2554   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -7.529E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.352E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.250    
 STEP TIME COMPLETED       2.500E+06,  TOTAL TIME COMPLETED        2.500E+06


  INCREMENT   251 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  356ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.609       TIME AVG. FORCE       0.354    
 LARGEST SCALED RESIDUAL FORCE     -5.299E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.299E-08
 LARGEST INCREMENT OF DISP.         4.548E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.610E-05   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.202E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.312E-04   TIME AVG. HEAT FLUX   2.801E-04
 LARGEST RESIDUAL HEAT FLUX         1.690E-15   AT NODE       6434   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -8.180E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.508E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.251    
 STEP TIME COMPLETED       2.510E+06,  TOTAL TIME COMPLETED        2.510E+06


  INCREMENT   252 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  355ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.610       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -5.548E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.548E-08
 LARGEST INCREMENT OF DISP.         4.627E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.903E-05   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.393E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.317E-04   TIME AVG. HEAT FLUX   2.811E-04
 LARGEST RESIDUAL HEAT FLUX        -1.380E-15   AT NODE       2566   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -8.847E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.669E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.252    
 STEP TIME COMPLETED       2.520E+06,  TOTAL TIME COMPLETED        2.520E+06


  INCREMENT   253 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  355ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.610       TIME AVG. FORCE       0.356    
 LARGEST SCALED RESIDUAL FORCE     -5.814E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.814E-08
 LARGEST INCREMENT OF DISP.         4.709E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.215E-05   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.595E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.320E-04   TIME AVG. HEAT FLUX   2.821E-04
 LARGEST RESIDUAL HEAT FLUX         3.941E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -9.530E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.837E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.253    
 STEP TIME COMPLETED       2.530E+06,  TOTAL TIME COMPLETED        2.530E+06


  INCREMENT   254 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.611       TIME AVG. FORCE       0.357    
 LARGEST SCALED RESIDUAL FORCE     -6.098E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.098E-08
 LARGEST INCREMENT OF DISP.         4.795E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.547E-05   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.811E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.324E-04   TIME AVG. HEAT FLUX   2.831E-04
 LARGEST RESIDUAL HEAT FLUX         1.665E-15   AT NODE       6991   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.023E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.012E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.254    
 STEP TIME COMPLETED       2.540E+06,  TOTAL TIME COMPLETED        2.540E+06


  INCREMENT   255 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.612       TIME AVG. FORCE       0.358    
 LARGEST SCALED RESIDUAL FORCE     -6.404E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.404E-08
 LARGEST INCREMENT OF DISP.         4.883E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.901E-05   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.040E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.323E-04   TIME AVG. HEAT FLUX   2.841E-04
 LARGEST RESIDUAL HEAT FLUX         1.452E-15   AT NODE       1581   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.095E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.194E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.255    
 STEP TIME COMPLETED       2.550E+06,  TOTAL TIME COMPLETED        2.550E+06


  INCREMENT   256 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.613       TIME AVG. FORCE       0.359    
 LARGEST SCALED RESIDUAL FORCE     -6.732E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.732E-08
 LARGEST INCREMENT OF DISP.         4.976E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.279E-05   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.285E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.324E-04   TIME AVG. HEAT FLUX   2.850E-04
 LARGEST RESIDUAL HEAT FLUX        -1.693E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.169E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.383E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.256    
 STEP TIME COMPLETED       2.560E+06,  TOTAL TIME COMPLETED        2.560E+06


  INCREMENT   257 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.613       TIME AVG. FORCE       0.360    
 LARGEST SCALED RESIDUAL FORCE     -7.086E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.086E-08
 LARGEST INCREMENT OF DISP.         5.072E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.686E-05   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.548E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.329E-04   TIME AVG. HEAT FLUX   2.860E-04
 LARGEST RESIDUAL HEAT FLUX         1.568E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.245E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.582E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.257    
 STEP TIME COMPLETED       2.570E+06,  TOTAL TIME COMPLETED        2.570E+06


  INCREMENT   258 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.614       TIME AVG. FORCE       0.361    
 LARGEST SCALED RESIDUAL FORCE     -7.468E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.468E-08
 LARGEST INCREMENT OF DISP.         5.173E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.012E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.831E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.330E-04   TIME AVG. HEAT FLUX   2.870E-04
 LARGEST RESIDUAL HEAT FLUX         1.804E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.323E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.790E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.258    
 STEP TIME COMPLETED       2.580E+06,  TOTAL TIME COMPLETED        2.580E+06


  INCREMENT   259 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.615       TIME AVG. FORCE       0.362    
 LARGEST SCALED RESIDUAL FORCE     -7.882E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.882E-08
 LARGEST INCREMENT OF DISP.         5.279E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.059E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.136E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.331E-04   TIME AVG. HEAT FLUX   2.879E-04
 LARGEST RESIDUAL HEAT FLUX        -3.102E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.403E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.008E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.259    
 STEP TIME COMPLETED       2.590E+06,  TOTAL TIME COMPLETED        2.590E+06


  INCREMENT   260 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.615       TIME AVG. FORCE       0.363    
 LARGEST SCALED RESIDUAL FORCE     -8.331E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.331E-08
 LARGEST INCREMENT OF DISP.         5.389E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.109E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.464E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.334E-04   TIME AVG. HEAT FLUX   2.889E-04
 LARGEST RESIDUAL HEAT FLUX         1.498E-15   AT NODE       4062   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.485E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.238E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.260    
 STEP TIME COMPLETED       2.600E+06,  TOTAL TIME COMPLETED        2.600E+06


  INCREMENT   261 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  327ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.616       TIME AVG. FORCE       0.364    
 LARGEST SCALED RESIDUAL FORCE     -8.820E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.820E-08
 LARGEST INCREMENT OF DISP.         5.505E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.163E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.818E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.333E-04   TIME AVG. HEAT FLUX   2.898E-04
 LARGEST RESIDUAL HEAT FLUX        -1.451E-15   AT NODE       2776   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.570E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.480E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.261    
 STEP TIME COMPLETED       2.610E+06,  TOTAL TIME COMPLETED        2.610E+06


  INCREMENT   262 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.617       TIME AVG. FORCE       0.365    
 LARGEST SCALED RESIDUAL FORCE     -9.354E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.354E-08
 LARGEST INCREMENT OF DISP.         5.626E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.222E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    8.203E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.331E-04   TIME AVG. HEAT FLUX   2.907E-04
 LARGEST RESIDUAL HEAT FLUX         1.756E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.657E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.736E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.262    
 STEP TIME COMPLETED       2.620E+06,  TOTAL TIME COMPLETED        2.620E+06


  INCREMENT   263 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  356ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.617       TIME AVG. FORCE       0.366    
 LARGEST SCALED RESIDUAL FORCE     -9.938E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.938E-08
 LARGEST INCREMENT OF DISP.         5.753E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.285E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    8.620E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.328E-04   TIME AVG. HEAT FLUX   2.917E-04
 LARGEST RESIDUAL HEAT FLUX         1.512E-15   AT NODE       2777   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.747E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.007E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.263    
 STEP TIME COMPLETED       2.630E+06,  TOTAL TIME COMPLETED        2.630E+06


  INCREMENT   264 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  346ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.618       TIME AVG. FORCE       0.367    
 LARGEST SCALED RESIDUAL FORCE     -1.058E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.058E-07
 LARGEST INCREMENT OF DISP.         5.887E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.353E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.074E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.327E-04   TIME AVG. HEAT FLUX   2.926E-04
 LARGEST RESIDUAL HEAT FLUX         1.392E-15   AT NODE       1591   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.840E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.293E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.264    
 STEP TIME COMPLETED       2.640E+06,  TOTAL TIME COMPLETED        2.640E+06


  INCREMENT   265 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  343ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.618       TIME AVG. FORCE       0.368    
 LARGEST SCALED RESIDUAL FORCE     -1.129E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.129E-07
 LARGEST INCREMENT OF DISP.         6.029E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.426E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.572E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.325E-04   TIME AVG. HEAT FLUX   2.935E-04
 LARGEST RESIDUAL HEAT FLUX         3.171E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.936E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.600E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.265    
 STEP TIME COMPLETED       2.650E+06,  TOTAL TIME COMPLETED        2.650E+06


  INCREMENT   266 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.619       TIME AVG. FORCE       0.369    
 LARGEST SCALED RESIDUAL FORCE     -1.207E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.207E-07
 LARGEST INCREMENT OF DISP.         6.178E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.506E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.011E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.322E-04   TIME AVG. HEAT FLUX   2.944E-04
 LARGEST RESIDUAL HEAT FLUX        -2.220E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.035E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.925E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.266    
 STEP TIME COMPLETED       2.660E+06,  TOTAL TIME COMPLETED        2.660E+06


  INCREMENT   267 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.619       TIME AVG. FORCE       0.370    
 LARGEST SCALED RESIDUAL FORCE     -1.293E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.293E-07
 LARGEST INCREMENT OF DISP.         6.335E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.592E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.071E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.317E-04   TIME AVG. HEAT FLUX   2.953E-04
 LARGEST RESIDUAL HEAT FLUX         3.719E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.138E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.027E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.267    
 STEP TIME COMPLETED       2.670E+06,  TOTAL TIME COMPLETED        2.670E+06


  INCREMENT   268 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.619       TIME AVG. FORCE       0.371    
 LARGEST SCALED RESIDUAL FORCE     -1.389E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.389E-07
 LARGEST INCREMENT OF DISP.         6.501E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.686E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.136E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.313E-04   TIME AVG. HEAT FLUX   2.962E-04
 LARGEST RESIDUAL HEAT FLUX        -1.499E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.245E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.064E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.268    
 STEP TIME COMPLETED       2.680E+06,  TOTAL TIME COMPLETED        2.680E+06


  INCREMENT   269 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.620       TIME AVG. FORCE       0.371    
 LARGEST SCALED RESIDUAL FORCE     -1.495E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.495E-07
 LARGEST INCREMENT OF DISP.         6.678E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.788E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.208E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.312E-04   TIME AVG. HEAT FLUX   2.970E-04
 LARGEST RESIDUAL HEAT FLUX        -1.936E-15   AT NODE      11987   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.355E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.103E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.269    
 STEP TIME COMPLETED       2.690E+06,  TOTAL TIME COMPLETED        2.690E+06


  INCREMENT   270 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.620       TIME AVG. FORCE       0.372    
 LARGEST SCALED RESIDUAL FORCE     -1.614E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.614E-07
 LARGEST INCREMENT OF DISP.         6.865E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.899E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.287E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.305E-04   TIME AVG. HEAT FLUX   2.979E-04
 LARGEST RESIDUAL HEAT FLUX         1.569E-15   AT NODE        568   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.469E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.145E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.270    
 STEP TIME COMPLETED       2.700E+06,  TOTAL TIME COMPLETED        2.700E+06


  INCREMENT   271 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.620       TIME AVG. FORCE       0.373    
 LARGEST SCALED RESIDUAL FORCE     -1.747E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.747E-07
 LARGEST INCREMENT OF DISP.         7.063E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.021E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.375E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.303E-04   TIME AVG. HEAT FLUX   2.988E-04
 LARGEST RESIDUAL HEAT FLUX         1.586E-15   AT NODE       2472   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.588E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.190E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.271    
 STEP TIME COMPLETED       2.710E+06,  TOTAL TIME COMPLETED        2.710E+06


  INCREMENT   272 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  353ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.620       TIME AVG. FORCE       0.374    
 LARGEST SCALED RESIDUAL FORCE     -1.896E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.896E-07
 LARGEST INCREMENT OF DISP.         7.275E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.154E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.472E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.299E-04   TIME AVG. HEAT FLUX   2.996E-04
 LARGEST RESIDUAL HEAT FLUX        -1.433E-15   AT NODE       3278   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.712E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.238E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.272    
 STEP TIME COMPLETED       2.720E+06,  TOTAL TIME COMPLETED        2.720E+06


  INCREMENT   273 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  332ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.620       TIME AVG. FORCE       0.375    
 LARGEST SCALED RESIDUAL FORCE     -2.064E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.064E-07
 LARGEST INCREMENT OF DISP.         7.501E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.300E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.579E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.297E-04   TIME AVG. HEAT FLUX   3.005E-04
 LARGEST RESIDUAL HEAT FLUX        -1.402E-15   AT NODE       3632   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.841E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.289E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.273    
 STEP TIME COMPLETED       2.730E+06,  TOTAL TIME COMPLETED        2.730E+06


  INCREMENT   274 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  353ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.620       TIME AVG. FORCE       0.376    
 LARGEST SCALED RESIDUAL FORCE     -2.254E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.254E-07
 LARGEST INCREMENT OF DISP.         7.742E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.461E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.700E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.289E-04   TIME AVG. HEAT FLUX   3.013E-04
 LARGEST RESIDUAL HEAT FLUX        -1.311E-15   AT NODE       6446   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.975E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.344E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.274    
 STEP TIME COMPLETED       2.740E+06,  TOTAL TIME COMPLETED        2.740E+06


  INCREMENT   275 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  338ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.620       TIME AVG. FORCE       0.377    
 LARGEST SCALED RESIDUAL FORCE     -2.469E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.469E-07
 LARGEST INCREMENT OF DISP.         8.000E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.639E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.834E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.284E-04   TIME AVG. HEAT FLUX   3.021E-04
 LARGEST RESIDUAL HEAT FLUX         2.859E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.116E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.403E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.275    
 STEP TIME COMPLETED       2.750E+06,  TOTAL TIME COMPLETED        2.750E+06


  INCREMENT   276 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  347ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.620       TIME AVG. FORCE       0.378    
 LARGEST SCALED RESIDUAL FORCE     -2.714E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.714E-07
 LARGEST INCREMENT OF DISP.         8.277E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.835E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.984E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.278E-04   TIME AVG. HEAT FLUX   3.030E-04
 LARGEST RESIDUAL HEAT FLUX         1.582E-15   AT NODE       6705   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.262E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.465E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.276    
 STEP TIME COMPLETED       2.760E+06,  TOTAL TIME COMPLETED        2.760E+06


  INCREMENT   277 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  345ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.620       TIME AVG. FORCE       0.379    
 LARGEST SCALED RESIDUAL FORCE     -2.994E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.994E-07
 LARGEST INCREMENT OF DISP.         8.574E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.054E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.154E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.261E-04   TIME AVG. HEAT FLUX   3.038E-04
 LARGEST RESIDUAL HEAT FLUX        -1.635E-15   AT NODE       1223   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.415E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.533E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.277    
 STEP TIME COMPLETED       2.770E+06,  TOTAL TIME COMPLETED        2.770E+06


  INCREMENT   278 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  347ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.620       TIME AVG. FORCE       0.379    
 LARGEST SCALED RESIDUAL FORCE     -3.316E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.316E-07
 LARGEST INCREMENT OF DISP.         8.894E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.296E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.347E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.247E-04   TIME AVG. HEAT FLUX   3.046E-04
 LARGEST RESIDUAL HEAT FLUX         1.749E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.576E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.606E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.278    
 STEP TIME COMPLETED       2.780E+06,  TOTAL TIME COMPLETED        2.780E+06


  INCREMENT   279 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  352ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.619       TIME AVG. FORCE       0.380    
 LARGEST SCALED RESIDUAL FORCE     -3.686E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.686E-07
 LARGEST INCREMENT OF DISP.         9.240E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.566E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.563E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.230E-04   TIME AVG. HEAT FLUX   3.054E-04
 LARGEST RESIDUAL HEAT FLUX         1.998E-15   AT NODE       6436   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.744E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.684E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.279    
 STEP TIME COMPLETED       2.790E+06,  TOTAL TIME COMPLETED        2.790E+06


  INCREMENT   280 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  352ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.619       TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -4.114E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.114E-07
 LARGEST INCREMENT OF DISP.         9.618E-03   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.868E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.810E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.213E-04   TIME AVG. HEAT FLUX   3.061E-04
 LARGEST RESIDUAL HEAT FLUX         1.846E-15   AT NODE       2548   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.921E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.767E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.280    
 STEP TIME COMPLETED       2.800E+06,  TOTAL TIME COMPLETED        2.800E+06


  INCREMENT   281 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  339ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.619       TIME AVG. FORCE       0.382    
 LARGEST SCALED RESIDUAL FORCE     -4.612E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.612E-07
 LARGEST INCREMENT OF DISP.         1.003E-02   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.206E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.090E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.201E-04   TIME AVG. HEAT FLUX   3.069E-04
 LARGEST RESIDUAL HEAT FLUX        -1.616E-15   AT NODE        814   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.107E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.857E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.281    
 STEP TIME COMPLETED       2.810E+06,  TOTAL TIME COMPLETED        2.810E+06


  INCREMENT   282 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.618       TIME AVG. FORCE       0.383    
 LARGEST SCALED RESIDUAL FORCE     -5.193E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.193E-07
 LARGEST INCREMENT OF DISP.         1.048E-02   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.584E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.410E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.183E-04   TIME AVG. HEAT FLUX   3.076E-04
 LARGEST RESIDUAL HEAT FLUX         1.488E-15   AT NODE       1450   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.302E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.953E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.282    
 STEP TIME COMPLETED       2.820E+06,  TOTAL TIME COMPLETED        2.820E+06


  INCREMENT   283 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  392ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.618       TIME AVG. FORCE       0.384    
 LARGEST SCALED RESIDUAL FORCE     -5.876E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.876E-07
 LARGEST INCREMENT OF DISP.         1.098E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.010E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.777E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.170E-04   TIME AVG. HEAT FLUX   3.084E-04
 LARGEST RESIDUAL HEAT FLUX         1.526E-15   AT NODE       5011   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.508E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.055E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.283    
 STEP TIME COMPLETED       2.830E+06,  TOTAL TIME COMPLETED        2.830E+06


  INCREMENT   284 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  353ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.617       TIME AVG. FORCE       0.385    
 LARGEST SCALED RESIDUAL FORCE     -6.685E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.685E-07
 LARGEST INCREMENT OF DISP.         1.153E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.489E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.134E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.148E-04   TIME AVG. HEAT FLUX   3.091E-04
 LARGEST RESIDUAL HEAT FLUX        -1.749E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.724E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.165E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.284    
 STEP TIME COMPLETED       2.840E+06,  TOTAL TIME COMPLETED        2.840E+06


  INCREMENT   285 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.616       TIME AVG. FORCE       0.385    
 LARGEST SCALED RESIDUAL FORCE     -7.640E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.640E-07
 LARGEST INCREMENT OF DISP.         1.213E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.026E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.518E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.133E-04   TIME AVG. HEAT FLUX   3.098E-04
 LARGEST RESIDUAL HEAT FLUX         1.915E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.952E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.279E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.285    
 STEP TIME COMPLETED       2.850E+06,  TOTAL TIME COMPLETED        2.850E+06


  INCREMENT   286 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.615       TIME AVG. FORCE       0.386    
 LARGEST SCALED RESIDUAL FORCE     -8.772E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.772E-07
 LARGEST INCREMENT OF DISP.         1.279E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.637E-04   AT NODE       6446   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.947E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.121E-04   TIME AVG. HEAT FLUX   3.105E-04
 LARGEST RESIDUAL HEAT FLUX         2.089E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.192E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.398E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.286    
 STEP TIME COMPLETED       2.860E+06,  TOTAL TIME COMPLETED        2.860E+06


  INCREMENT   287 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.614       TIME AVG. FORCE       0.387    
 LARGEST SCALED RESIDUAL FORCE     -1.013E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.013E-06
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.328E-04   AT NODE       6446   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.425E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.109E-04   TIME AVG. HEAT FLUX   3.112E-04
 LARGEST RESIDUAL HEAT FLUX         1.573E-15   AT NODE       1770   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.444E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.520E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.287    
 STEP TIME COMPLETED       2.870E+06,  TOTAL TIME COMPLETED        2.870E+06


  INCREMENT   288 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.612       TIME AVG. FORCE       0.388    
 LARGEST SCALED RESIDUAL FORCE     -1.176E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.176E-06
 LARGEST INCREMENT OF DISP.         1.432E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.101E-04   AT NODE       6446   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.955E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.103E-04   TIME AVG. HEAT FLUX   3.119E-04
 LARGEST RESIDUAL HEAT FLUX        -1.405E-15   AT NODE        177   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.708E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.645E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.288    
 STEP TIME COMPLETED       2.880E+06,  TOTAL TIME COMPLETED        2.880E+06


  INCREMENT   289 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  371ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.611       TIME AVG. FORCE       0.388    
 LARGEST SCALED RESIDUAL FORCE     -1.371E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.371E-06
 LARGEST INCREMENT OF DISP.         1.520E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.970E-04   AT NODE        178   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.521E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.103E-04   TIME AVG. HEAT FLUX   3.126E-04
 LARGEST RESIDUAL HEAT FLUX         1.398E-15   AT NODE       5306   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.984E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.761E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.289    
 STEP TIME COMPLETED       2.890E+06,  TOTAL TIME COMPLETED        2.890E+06


  INCREMENT   290 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  389ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.609       TIME AVG. FORCE       0.389    
 LARGEST SCALED RESIDUAL FORCE     -1.608E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.608E-06
 LARGEST INCREMENT OF DISP.         1.618E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.939E-04   AT NODE        178   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.150E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.093E-04   TIME AVG. HEAT FLUX   3.133E-04
 LARGEST RESIDUAL HEAT FLUX         1.570E-15   AT NODE       3076   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -6.272E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.875E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.290    
 STEP TIME COMPLETED       2.900E+06,  TOTAL TIME COMPLETED        2.900E+06


  INCREMENT   291 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  389ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.607       TIME AVG. FORCE       0.390    
 LARGEST SCALED RESIDUAL FORCE     -1.886E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.886E-06
 LARGEST INCREMENT OF DISP.         1.725E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.101E-03   AT NODE        178   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.761E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.087E-04   TIME AVG. HEAT FLUX   3.140E-04
 LARGEST RESIDUAL HEAT FLUX        -1.447E-15   AT NODE       3142   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -6.567E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.953E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.291    
 STEP TIME COMPLETED       2.910E+06,  TOTAL TIME COMPLETED        2.910E+06


  INCREMENT   292 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  404ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.605       TIME AVG. FORCE       0.391    
 LARGEST SCALED RESIDUAL FORCE     -2.228E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.228E-06
 LARGEST INCREMENT OF DISP.         1.844E-02   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.217E-03   AT NODE       6438   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    8.419E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.091E-04   TIME AVG. HEAT FLUX   3.147E-04
 LARGEST RESIDUAL HEAT FLUX         1.344E-15   AT NODE        982   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -6.868E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.012E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.292    
 STEP TIME COMPLETED       2.920E+06,  TOTAL TIME COMPLETED        2.920E+06


  INCREMENT   293 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  395ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.602       TIME AVG. FORCE       0.391    
 LARGEST SCALED RESIDUAL FORCE     -2.636E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.636E-06
 LARGEST INCREMENT OF DISP.         1.973E-02   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.341E-03   AT NODE       6438   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.056E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.100E-04   TIME AVG. HEAT FLUX   3.153E-04
 LARGEST RESIDUAL HEAT FLUX         1.547E-15   AT NODE       2175   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -7.171E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.022E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.293    
 STEP TIME COMPLETED       2.930E+06,  TOTAL TIME COMPLETED        2.930E+06


  INCREMENT   294 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  407ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.599       TIME AVG. FORCE       0.392    
 LARGEST SCALED RESIDUAL FORCE     -3.106E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.106E-06
 LARGEST INCREMENT OF DISP.         2.113E-02   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.473E-03   AT NODE        179   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.602E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.123E-04   TIME AVG. HEAT FLUX   3.160E-04
 LARGEST RESIDUAL HEAT FLUX         1.957E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -7.464E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.938E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.294    
 STEP TIME COMPLETED       2.940E+06,  TOTAL TIME COMPLETED        2.940E+06


  INCREMENT   295 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  368ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.596       TIME AVG. FORCE       0.393    
 LARGEST SCALED RESIDUAL FORCE     -3.684E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.684E-06
 LARGEST INCREMENT OF DISP.         2.264E-02   AT NODE       6446   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.607E-03   AT NODE        180   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.014E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.129E-04   TIME AVG. HEAT FLUX   3.167E-04
 LARGEST RESIDUAL HEAT FLUX        -2.245E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -7.744E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.039E-03   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.295    
 STEP TIME COMPLETED       2.950E+06,  TOTAL TIME COMPLETED        2.950E+06


  INCREMENT   296 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  389ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.593       TIME AVG. FORCE       0.394    
 LARGEST SCALED RESIDUAL FORCE     -4.289E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.289E-06
 LARGEST INCREMENT OF DISP.         2.420E-02   AT NODE       6446   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.740E-03   AT NODE        181   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.039E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.162E-04   TIME AVG. HEAT FLUX   3.173E-04
 LARGEST RESIDUAL HEAT FLUX         2.776E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -7.989E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.364E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.296    
 STEP TIME COMPLETED       2.960E+06,  TOTAL TIME COMPLETED        2.960E+06


  INCREMENT   297 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.589       TIME AVG. FORCE       0.394    
 LARGEST SCALED RESIDUAL FORCE     -4.871E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.871E-06
 LARGEST INCREMENT OF DISP.         2.578E-02   AT NODE        178   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.873E-03   AT NODE        182   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.033E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.214E-04   TIME AVG. HEAT FLUX   3.180E-04
 LARGEST RESIDUAL HEAT FLUX        -1.293E-15   AT NODE        189   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -8.176E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.049E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.297    
 STEP TIME COMPLETED       2.970E+06,  TOTAL TIME COMPLETED        2.970E+06


  INCREMENT   298 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  384ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.584       TIME AVG. FORCE       0.395    
 LARGEST SCALED RESIDUAL FORCE     -5.332E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.332E-06
 LARGEST INCREMENT OF DISP.         2.727E-02   AT NODE       6438   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.005E-03   AT NODE       6457   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.928E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.293E-04   TIME AVG. HEAT FLUX   3.187E-04
 LARGEST RESIDUAL HEAT FLUX        -1.349E-15   AT NODE       3014   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -8.275E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.791E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.298    
 STEP TIME COMPLETED       2.980E+06,  TOTAL TIME COMPLETED        2.980E+06


  INCREMENT   299 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  388ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.579       TIME AVG. FORCE       0.395    
 LARGEST SCALED RESIDUAL FORCE     -5.660E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.660E-06
 LARGEST INCREMENT OF DISP.         2.863E-02   AT NODE       6442   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.135E-03   AT NODE       5723   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.391E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.402E-04   TIME AVG. HEAT FLUX   3.195E-04
 LARGEST RESIDUAL HEAT FLUX         1.363E-15   AT NODE       3487   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -8.268E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.583E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.299    
 STEP TIME COMPLETED       2.990E+06,  TOTAL TIME COMPLETED        2.990E+06


  INCREMENT   300 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  388ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.573       TIME AVG. FORCE       0.396    
 LARGEST SCALED RESIDUAL FORCE     -6.365E-06   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.365E-06
 LARGEST INCREMENT OF DISP.         2.987E-02   AT NODE       6452   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.255E-03   AT NODE        186   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.710E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.539E-04   TIME AVG. HEAT FLUX   3.203E-04
 LARGEST RESIDUAL HEAT FLUX        -1.429E-15   AT NODE        366   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -8.120E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.181E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.300    
 STEP TIME COMPLETED       3.000E+06,  TOTAL TIME COMPLETED        3.000E+06


  INCREMENT   301 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  424ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.567       TIME AVG. FORCE       0.397    
 LARGEST SCALED RESIDUAL FORCE     -7.575E-06   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.575E-06
 LARGEST INCREMENT OF DISP.         3.100E-02   AT NODE       6455   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.352E-03   AT NODE       6458   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.084E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.721E-04   TIME AVG. HEAT FLUX   3.211E-04
 LARGEST RESIDUAL HEAT FLUX        -2.109E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -7.789E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.882E-03   AT NODE       6446   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.301    
 STEP TIME COMPLETED       3.010E+06,  TOTAL TIME COMPLETED        3.010E+06


  INCREMENT   302 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  398ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.560       TIME AVG. FORCE       0.397    
 LARGEST SCALED RESIDUAL FORCE     -7.788E-06   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.788E-06
 LARGEST INCREMENT OF DISP.         3.206E-02   AT NODE        183   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -2.939E-03   AT NODE       2548   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.382E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.940E-04   TIME AVG. HEAT FLUX   3.220E-04
 LARGEST RESIDUAL HEAT FLUX        -1.462E-15   AT NODE       2263   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -7.230E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.311E-03   AT NODE        178   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.302    
 STEP TIME COMPLETED       3.020E+06,  TOTAL TIME COMPLETED        3.020E+06


  INCREMENT   303 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  397ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.553       TIME AVG. FORCE       0.398    
 LARGEST SCALED RESIDUAL FORCE      8.235E-06   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      8.235E-06
 LARGEST INCREMENT OF DISP.         3.304E-02   AT NODE        184   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.265E-03   AT NODE      11988   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.492E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.204E-04   TIME AVG. HEAT FLUX   3.230E-04
 LARGEST RESIDUAL HEAT FLUX        -1.560E-15   AT NODE       1113   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -6.494E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.380E-03   AT NODE       6438   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.303    
 STEP TIME COMPLETED       3.030E+06,  TOTAL TIME COMPLETED        3.030E+06


  INCREMENT   304 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  394ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.545       TIME AVG. FORCE       0.398    
 LARGEST SCALED RESIDUAL FORCE      9.116E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      9.116E-06
 LARGEST INCREMENT OF DISP.         3.384E-02   AT NODE       6902   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.114E-03   AT NODE       6446   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -4.297E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.495E-04   TIME AVG. HEAT FLUX   3.241E-04
 LARGEST RESIDUAL HEAT FLUX        -1.411E-15   AT NODE       5415   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.649E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.449E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.304    
 STEP TIME COMPLETED       3.040E+06,  TOTAL TIME COMPLETED        3.040E+06


  INCREMENT   305 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  398ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.538       TIME AVG. FORCE       0.399    
 LARGEST SCALED RESIDUAL FORCE      1.567E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.567E-05
 LARGEST INCREMENT OF DISP.         3.440E-02   AT NODE       6907   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.502E-03   AT NODE       6441   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.093E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.775E-04   TIME AVG. HEAT FLUX   3.252E-04
 LARGEST RESIDUAL HEAT FLUX         1.447E-15   AT NODE       2562   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.099E-02   AT NODE       2548   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        9.141E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.305    
 STEP TIME COMPLETED       3.050E+06,  TOTAL TIME COMPLETED        3.050E+06


  INCREMENT   306 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  386ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.531       TIME AVG. FORCE       0.399    
 LARGEST SCALED RESIDUAL FORCE      2.287E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.287E-05
 LARGEST INCREMENT OF DISP.         3.463E-02   AT NODE       6458   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.270E-03   AT NODE       6439   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.016E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.026E-04   TIME AVG. HEAT FLUX   3.265E-04
 LARGEST RESIDUAL HEAT FLUX         1.613E-15   AT NODE       1987   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.671E-02   AT NODE        364   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        9.364E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.306    
 STEP TIME COMPLETED       3.060E+06,  TOTAL TIME COMPLETED        3.060E+06


  INCREMENT   307 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  385ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.523       TIME AVG. FORCE       0.399    
 LARGEST SCALED RESIDUAL FORCE      2.694E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.694E-05
 LARGEST INCREMENT OF DISP.         3.467E-02   AT NODE        188   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.114E-03   AT NODE       6451   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.496E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.229E-04   TIME AVG. HEAT FLUX   3.278E-04
 LARGEST RESIDUAL HEAT FLUX         1.687E-15   AT NODE       3633   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.018E-02   AT NODE        364   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.271E-02   AT NODE       3126   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.307    
 STEP TIME COMPLETED       3.070E+06,  TOTAL TIME COMPLETED        3.070E+06


  INCREMENT   308 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  385ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.515       TIME AVG. FORCE       0.400    
 LARGEST SCALED RESIDUAL FORCE      2.889E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.889E-05
 LARGEST INCREMENT OF DISP.         3.458E-02   AT NODE        189   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.971E-03   AT NODE       6454   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.402E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.390E-04   TIME AVG. HEAT FLUX   3.291E-04
 LARGEST RESIDUAL HEAT FLUX         1.573E-15   AT NODE       1252   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.399E-02   AT NODE        365   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.173E-02   AT NODE       3126   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.308    
 STEP TIME COMPLETED       3.080E+06,  TOTAL TIME COMPLETED        3.080E+06


  INCREMENT   309 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  352ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.508       TIME AVG. FORCE       0.400    
 LARGEST SCALED RESIDUAL FORCE      2.630E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.630E-05
 LARGEST INCREMENT OF DISP.         3.445E-02   AT NODE        190   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.812E-03   AT NODE       3050   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.222E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.532E-04   TIME AVG. HEAT FLUX   3.305E-04
 LARGEST RESIDUAL HEAT FLUX        -1.454E-15   AT NODE       1152   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.971E-02   AT NODE        367   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.006E-02   AT NODE       3127   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.309    
 STEP TIME COMPLETED       3.090E+06,  TOTAL TIME COMPLETED        3.090E+06


  INCREMENT   310 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  379ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.501       TIME AVG. FORCE       0.400    
 LARGEST SCALED RESIDUAL FORCE      2.152E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.152E-05
 LARGEST INCREMENT OF DISP.         3.430E-02   AT NODE       6471   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.751E-03   AT NODE        377   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.714E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.647E-04   TIME AVG. HEAT FLUX   3.319E-04
 LARGEST RESIDUAL HEAT FLUX        -1.469E-15   AT NODE       3258   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.696E-02   AT NODE        368   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.128E-03   AT NODE       3125   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.310    
 STEP TIME COMPLETED       3.100E+06,  TOTAL TIME COMPLETED        3.100E+06


  INCREMENT   311 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  389ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.493       TIME AVG. FORCE       0.401    
 LARGEST SCALED RESIDUAL FORCE      1.623E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.623E-05
 LARGEST INCREMENT OF DISP.         3.417E-02   AT NODE       6883   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.685E-03   AT NODE       3051   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.122E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.753E-04   TIME AVG. HEAT FLUX   3.333E-04
 LARGEST RESIDUAL HEAT FLUX        -1.384E-15   AT NODE        473   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.518E-02   AT NODE        376   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.773E-03   AT NODE       3128   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.311    
 STEP TIME COMPLETED       3.110E+06,  TOTAL TIME COMPLETED        3.110E+06


  INCREMENT   312 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  406ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.486       TIME AVG. FORCE       0.401    
 LARGEST SCALED RESIDUAL FORCE      1.346E-05   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.346E-05
 LARGEST INCREMENT OF DISP.         3.404E-02   AT NODE       6473   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.630E-03   AT NODE        379   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.788E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.842E-04   TIME AVG. HEAT FLUX   3.348E-04
 LARGEST RESIDUAL HEAT FLUX         1.536E-15   AT NODE        886   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.377E-02   AT NODE        377   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.863E-03   AT NODE       3129   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.312    
 STEP TIME COMPLETED       3.120E+06,  TOTAL TIME COMPLETED        3.120E+06


  INCREMENT   313 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  379ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.479       TIME AVG. FORCE       0.401    
 LARGEST SCALED RESIDUAL FORCE      9.032E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      9.032E-06
 LARGEST INCREMENT OF DISP.         3.395E-02   AT NODE       6477   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.595E-03   AT NODE       3133   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.241E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.926E-04   TIME AVG. HEAT FLUX   3.362E-04
 LARGEST RESIDUAL HEAT FLUX        -1.399E-15   AT NODE       4855   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.279E-02   AT NODE        378   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.532E-03   AT NODE       3130   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.313    
 STEP TIME COMPLETED       3.130E+06,  TOTAL TIME COMPLETED        3.130E+06


  INCREMENT   314 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  394ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.472       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE      1.227E-05   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.227E-05
 LARGEST INCREMENT OF DISP.         3.384E-02   AT NODE       6888   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.555E-03   AT NODE        381   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.710E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.997E-04   TIME AVG. HEAT FLUX   3.377E-04
 LARGEST RESIDUAL HEAT FLUX         1.396E-15   AT NODE        650   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.197E-02   AT NODE        379   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.079E-03   AT NODE       3050   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.314    
 STEP TIME COMPLETED       3.140E+06,  TOTAL TIME COMPLETED        3.140E+06


  INCREMENT   315 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  378ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.465       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE     -1.288E-05   AT NODE       6709   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.288E-05
 LARGEST INCREMENT OF DISP.         3.380E-02   AT NODE        197   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.538E-03   AT NODE       3053   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.801E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.068E-04   TIME AVG. HEAT FLUX   3.392E-04
 LARGEST RESIDUAL HEAT FLUX        -1.478E-15   AT NODE        617   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.141E-02   AT NODE        380   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.008E-03   AT NODE       3132   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.315    
 STEP TIME COMPLETED       3.150E+06,  TOTAL TIME COMPLETED        3.150E+06


  INCREMENT   316 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  371ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.458       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE      1.402E-05   AT NODE       6709   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.402E-05
 LARGEST INCREMENT OF DISP.         3.373E-02   AT NODE        198   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.520E-03   AT NODE        369   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.006E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.119E-04   TIME AVG. HEAT FLUX   3.407E-04
 LARGEST RESIDUAL HEAT FLUX         1.547E-15   AT NODE       4883   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.093E-02   AT NODE       6463   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.174E-03   AT NODE        297   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.316    
 STEP TIME COMPLETED       3.160E+06,  TOTAL TIME COMPLETED        3.160E+06


  INCREMENT   317 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  384ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.451       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE     -1.725E-05   AT NODE       6709   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.725E-05
 LARGEST INCREMENT OF DISP.         3.376E-02   AT NODE          2   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.507E-03   AT NODE       3055   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.451E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.183E-04   TIME AVG. HEAT FLUX   3.422E-04
 LARGEST RESIDUAL HEAT FLUX        -1.308E-15   AT NODE       3731   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.058E-02   AT NODE       6143   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.441E-03   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.317    
 STEP TIME COMPLETED       3.170E+06,  TOTAL TIME COMPLETED        3.170E+06


  INCREMENT   318 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  384ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.445       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE     -3.372E-05   AT NODE       6458   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.372E-05
 LARGEST INCREMENT OF DISP.         3.371E-02   AT NODE        199   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.508E-03   AT NODE        388   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.876E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.228E-04   TIME AVG. HEAT FLUX   3.437E-04
 LARGEST RESIDUAL HEAT FLUX         1.375E-15   AT NODE       1580   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.025E-02   AT NODE       6468   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.424E-03   AT NODE        297   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.318    
 STEP TIME COMPLETED       3.180E+06,  TOTAL TIME COMPLETED        3.180E+06


  INCREMENT   319 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.438       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE     -7.120E-05   AT NODE       6459   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.120E-05
 LARGEST INCREMENT OF DISP.         3.378E-02   AT NODE        200   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.514E-03   AT NODE      12015   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.031E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.295E-04   TIME AVG. HEAT FLUX   3.453E-04
 LARGEST RESIDUAL HEAT FLUX         1.517E-15   AT NODE        770   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.001E-02   AT NODE       6466   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.762E-03   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.319    
 STEP TIME COMPLETED       3.190E+06,  TOTAL TIME COMPLETED        3.190E+06


  INCREMENT   320 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.431       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE     -1.333E-04   AT NODE       6144   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.333E-04
 LARGEST INCREMENT OF DISP.         3.376E-02   AT NODE       6493   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.545E-03   AT NODE       6891   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.975E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.322E-04   TIME AVG. HEAT FLUX   3.468E-04
 LARGEST RESIDUAL HEAT FLUX        -1.494E-15   AT NODE       4762   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.978E-02   AT NODE       6470   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.028E-03   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.320    
 STEP TIME COMPLETED       3.200E+06,  TOTAL TIME COMPLETED        3.200E+06


  INCREMENT   321 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  378ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.425       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE     -2.092E-04   AT NODE       6146   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.092E-04
 LARGEST INCREMENT OF DISP.         3.388E-02   AT NODE       6496   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.537E-03   AT NODE      12022   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.083E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.391E-04   TIME AVG. HEAT FLUX   3.483E-04
 LARGEST RESIDUAL HEAT FLUX        -1.432E-15   AT NODE       2828   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.962E-02   AT NODE       6884   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.457E-03   AT NODE        297   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.321    
 STEP TIME COMPLETED       3.210E+06,  TOTAL TIME COMPLETED        3.210E+06


  INCREMENT   322 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  392ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.418       TIME AVG. FORCE       0.403    
 LARGEST SCALED RESIDUAL FORCE     -2.742E-04   AT NODE       6467   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.742E-04
 LARGEST INCREMENT OF DISP.         3.387E-02   AT NODE       6772   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.601E-03   AT NODE       6893   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.204E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.404E-04   TIME AVG. HEAT FLUX   3.499E-04
 LARGEST RESIDUAL HEAT FLUX        -1.455E-15   AT NODE       4410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.945E-02   AT NODE       6475   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.898E-03   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.322    
 STEP TIME COMPLETED       3.220E+06,  TOTAL TIME COMPLETED        3.220E+06


  INCREMENT   323 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.412       TIME AVG. FORCE       0.403    
 LARGEST SCALED RESIDUAL FORCE     -3.715E-04   AT NODE       6465   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.715E-04
 LARGEST INCREMENT OF DISP.         3.402E-02   AT NODE        205   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.569E-03   AT NODE      12024   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.094E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.475E-04   TIME AVG. HEAT FLUX   3.514E-04
 LARGEST RESIDUAL HEAT FLUX         1.507E-15   AT NODE        395   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.933E-02   AT NODE       6479   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.465E-03   AT NODE        297   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.323    
 STEP TIME COMPLETED       3.230E+06,  TOTAL TIME COMPLETED        3.230E+06


  INCREMENT   324 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  388ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.405       TIME AVG. FORCE       0.403    
 LARGEST SCALED RESIDUAL FORCE     -4.674E-04   AT NODE       6472   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.674E-04
 LARGEST INCREMENT OF DISP.         3.402E-02   AT NODE        206   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.668E-03   AT NODE       6899   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -6.654E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.471E-04   TIME AVG. HEAT FLUX   3.530E-04
 LARGEST RESIDUAL HEAT FLUX         1.471E-15   AT NODE       3467   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.922E-02   AT NODE       3060   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.100E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.324    
 STEP TIME COMPLETED       3.240E+06,  TOTAL TIME COMPLETED        3.240E+06


  INCREMENT   325 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  399ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.399       TIME AVG. FORCE       0.403    
 LARGEST SCALED RESIDUAL FORCE     -6.047E-04   AT NODE       6472   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.047E-04
 LARGEST INCREMENT OF DISP.         3.421E-02   AT NODE        207   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.832E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -7.774E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.546E-04   TIME AVG. HEAT FLUX   3.545E-04
 LARGEST RESIDUAL HEAT FLUX         1.285E-15   AT NODE        927   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.914E-02   AT NODE        374   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.169E-02   AT NODE        297   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.325    
 STEP TIME COMPLETED       3.250E+06,  TOTAL TIME COMPLETED        3.250E+06


  INCREMENT   326 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  391ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.392       TIME AVG. FORCE       0.403    
 LARGEST SCALED RESIDUAL FORCE     -6.706E-04   AT NODE       6476   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.706E-04
 LARGEST INCREMENT OF DISP.         3.420E-02   AT NODE       5713   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.597E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.039E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.531E-04   TIME AVG. HEAT FLUX   3.560E-04
 LARGEST RESIDUAL HEAT FLUX         1.477E-15   AT NODE       1040   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.909E-02   AT NODE       3061   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.329E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.326    
 STEP TIME COMPLETED       3.260E+06,  TOTAL TIME COMPLETED        3.260E+06


  INCREMENT   327 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  425ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.386       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE     -8.174E-04   AT NODE       6474   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.174E-04
 LARGEST INCREMENT OF DISP.         3.442E-02   AT NODE       6152   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.837E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.162E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.611E-04   TIME AVG. HEAT FLUX   3.576E-04
 LARGEST RESIDUAL HEAT FLUX        -1.428E-15   AT NODE       3642   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.907E-02   AT NODE       6486   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.402E-02   AT NODE        297   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.327    
 STEP TIME COMPLETED       3.270E+06,  TOTAL TIME COMPLETED        3.270E+06


  INCREMENT   328 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  385ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.379       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE     -8.640E-04   AT NODE       6478   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.640E-04
 LARGEST INCREMENT OF DISP.         3.441E-02   AT NODE       6511   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.607E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.451E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.579E-04   TIME AVG. HEAT FLUX   3.591E-04
 LARGEST RESIDUAL HEAT FLUX        -1.401E-15   AT NODE        205   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.906E-02   AT NODE       6775   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.569E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.328    
 STEP TIME COMPLETED       3.280E+06,  TOTAL TIME COMPLETED        3.280E+06


  INCREMENT   329 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  382ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.373       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE     -8.863E-04   AT NODE       6483   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.863E-04
 LARGEST INCREMENT OF DISP.         3.464E-02   AT NODE        211   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.885E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.392E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.665E-04   TIME AVG. HEAT FLUX   3.607E-04
 LARGEST RESIDUAL HEAT FLUX         1.575E-15   AT NODE       1333   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.901E-02   AT NODE       6488   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.635E-02   AT NODE        297   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.329    
 STEP TIME COMPLETED       3.290E+06,  TOTAL TIME COMPLETED        3.290E+06


  INCREMENT   330 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  378ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.366       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE     -8.788E-04   AT NODE       6481   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.788E-04
 LARGEST INCREMENT OF DISP.         3.464E-02   AT NODE        212   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -7.629E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.587E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.617E-04   TIME AVG. HEAT FLUX   3.622E-04
 LARGEST RESIDUAL HEAT FLUX         1.725E-15   AT NODE       5455   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.897E-02   AT NODE       6493   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.817E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.330    
 STEP TIME COMPLETED       3.300E+06,  TOTAL TIME COMPLETED        3.300E+06


  INCREMENT   331 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  368ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.359       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE     -8.433E-04   AT NODE       5718   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.433E-04
 LARGEST INCREMENT OF DISP.         3.488E-02   AT NODE        213   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.936E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.455E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.721E-04   TIME AVG. HEAT FLUX   3.637E-04
 LARGEST RESIDUAL HEAT FLUX         1.296E-15   AT NODE       5378   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.895E-02   AT NODE        323   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.888E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.331    
 STEP TIME COMPLETED       3.310E+06,  TOTAL TIME COMPLETED        3.310E+06


  INCREMENT   332 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.353       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE     -6.181E-04   AT NODE       6487   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.181E-04
 LARGEST INCREMENT OF DISP.         3.487E-02   AT NODE       6530   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -8.536E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.178E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.661E-04   TIME AVG. HEAT FLUX   3.652E-04
 LARGEST RESIDUAL HEAT FLUX         1.250E-15   AT NODE       2782   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.895E-02   AT NODE       6502   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.065E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.332    
 STEP TIME COMPLETED       3.320E+06,  TOTAL TIME COMPLETED        3.320E+06


  INCREMENT   333 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  371ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.346       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE     -4.798E-04   AT NODE       6491   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.798E-04
 LARGEST INCREMENT OF DISP.         3.513E-02   AT NODE       6526   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.762E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -8.961E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.776E-04   TIME AVG. HEAT FLUX   3.668E-04
 LARGEST RESIDUAL HEAT FLUX         1.348E-15   AT NODE       4644   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.901E-02   AT NODE       6500   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.135E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.333    
 STEP TIME COMPLETED       3.330E+06,  TOTAL TIME COMPLETED        3.330E+06


  INCREMENT   334 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  368ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.340       TIME AVG. FORCE       0.402    
 LARGEST SCALED RESIDUAL FORCE     -3.230E-04   AT NODE       6494   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.230E-04
 LARGEST INCREMENT OF DISP.         3.511E-02   AT NODE        217   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -9.095E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -6.472E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.702E-04   TIME AVG. HEAT FLUX   3.683E-04
 LARGEST RESIDUAL HEAT FLUX         1.453E-15   AT NODE        623   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.902E-02   AT NODE       6776   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.274E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.334    
 STEP TIME COMPLETED       3.340E+06,  TOTAL TIME COMPLETED        3.340E+06


  INCREMENT   335 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  344ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.334       TIME AVG. FORCE       0.401    
 LARGEST SCALED RESIDUAL FORCE     -1.836E-04   AT NODE       6497   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.836E-04
 LARGEST INCREMENT OF DISP.         3.538E-02   AT NODE        218   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.075E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.629E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.819E-04   TIME AVG. HEAT FLUX   3.698E-04
 LARGEST RESIDUAL HEAT FLUX        -1.552E-15   AT NODE       2681   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.903E-02   AT NODE       6506   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.308E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.335    
 STEP TIME COMPLETED       3.350E+06,  TOTAL TIME COMPLETED        3.350E+06


  INCREMENT   336 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.327       TIME AVG. FORCE       0.401    
 LARGEST SCALED RESIDUAL FORCE      1.163E-04   AT NODE       6709   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.163E-04
 LARGEST INCREMENT OF DISP.         3.535E-02   AT NODE       6541   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -9.252E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.519E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.726E-04   TIME AVG. HEAT FLUX   3.713E-04
 LARGEST RESIDUAL HEAT FLUX         1.389E-15   AT NODE        799   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.899E-02   AT NODE        325   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.372E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.336    
 STEP TIME COMPLETED       3.360E+06,  TOTAL TIME COMPLETED        3.360E+06


  INCREMENT   337 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.321       TIME AVG. FORCE       0.401    
 LARGEST SCALED RESIDUAL FORCE     -1.149E-04   AT NODE       6709   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.149E-04
 LARGEST INCREMENT OF DISP.         3.564E-02   AT NODE       6544   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.070E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.316E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.840E-04   TIME AVG. HEAT FLUX   3.729E-04
 LARGEST RESIDUAL HEAT FLUX         1.280E-15   AT NODE       1377   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.904E-02   AT NODE        326   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.356E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.337    
 STEP TIME COMPLETED       3.370E+06,  TOTAL TIME COMPLETED        3.370E+06


  INCREMENT   338 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  275ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.315       TIME AVG. FORCE       0.401    
 LARGEST SCALED RESIDUAL FORCE      1.161E-04   AT NODE       6709   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.161E-04
 LARGEST INCREMENT OF DISP.         3.560E-02   AT NODE          4   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -9.071E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.562E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.730E-04   TIME AVG. HEAT FLUX   3.743E-04
 LARGEST RESIDUAL HEAT FLUX         1.447E-15   AT NODE        884   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.910E-02   AT NODE       6518   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.382E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.338    
 STEP TIME COMPLETED       3.380E+06,  TOTAL TIME COMPLETED        3.380E+06


  INCREMENT   339 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  376ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.309       TIME AVG. FORCE       0.400    
 LARGEST SCALED RESIDUAL FORCE     -1.133E-04   AT NODE       6709   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.133E-04
 LARGEST INCREMENT OF DISP.         3.589E-02   AT NODE        222   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.804E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.317E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.840E-04   TIME AVG. HEAT FLUX   3.759E-04
 LARGEST RESIDUAL HEAT FLUX        -1.419E-15   AT NODE       3267   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.914E-02   AT NODE       6516   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.338E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.339    
 STEP TIME COMPLETED       3.390E+06,  TOTAL TIME COMPLETED        3.390E+06


  INCREMENT   340 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  347ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.303       TIME AVG. FORCE       0.400    
 LARGEST SCALED RESIDUAL FORCE      1.130E-04   AT NODE       6709   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.130E-04
 LARGEST INCREMENT OF DISP.         3.586E-02   AT NODE       6550   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -8.712E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.512E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.728E-04   TIME AVG. HEAT FLUX   3.773E-04
 LARGEST RESIDUAL HEAT FLUX         1.312E-15   AT NODE       4307   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.913E-02   AT NODE       6524   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.334E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.340    
 STEP TIME COMPLETED       3.400E+06,  TOTAL TIME COMPLETED        3.400E+06


  INCREMENT   341 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.297       TIME AVG. FORCE       0.400    
 LARGEST SCALED RESIDUAL FORCE     -1.095E-04   AT NODE       6709   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.095E-04
 LARGEST INCREMENT OF DISP.         3.616E-02   AT NODE       6553   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.410E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.267E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.835E-04   TIME AVG. HEAT FLUX   3.788E-04
 LARGEST RESIDUAL HEAT FLUX        -1.445E-15   AT NODE       3832   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.917E-02   AT NODE        398   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.276E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.341    
 STEP TIME COMPLETED       3.410E+06,  TOTAL TIME COMPLETED        3.410E+06


  INCREMENT   342 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  315ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.291       TIME AVG. FORCE       0.399    
 LARGEST SCALED RESIDUAL FORCE      1.084E-04   AT NODE       6709   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.084E-04
 LARGEST INCREMENT OF DISP.         3.613E-02   AT NODE        226   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -8.261E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.417E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.720E-04   TIME AVG. HEAT FLUX   3.803E-04
 LARGEST RESIDUAL HEAT FLUX         1.509E-15   AT NODE       3629   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.924E-02   AT NODE       6526   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.255E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.342    
 STEP TIME COMPLETED       3.420E+06,  TOTAL TIME COMPLETED        3.420E+06


  INCREMENT   343 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.286       TIME AVG. FORCE       0.399    
 LARGEST SCALED RESIDUAL FORCE     -1.412E-04   AT NODE      12054   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.412E-04
 LARGEST INCREMENT OF DISP.         3.643E-02   AT NODE        227   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.959E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.950E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.822E-04   TIME AVG. HEAT FLUX   3.817E-04
 LARGEST RESIDUAL HEAT FLUX         1.478E-15   AT NODE       4459   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.931E-02   AT NODE       6536   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.192E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.343    
 STEP TIME COMPLETED       3.430E+06,  TOTAL TIME COMPLETED        3.430E+06


  INCREMENT   344 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  339ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.279       TIME AVG. FORCE       0.399    
 LARGEST SCALED RESIDUAL FORCE     -3.158E-04   AT NODE      12056   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.158E-04
 LARGEST INCREMENT OF DISP.         3.642E-02   AT NODE       5697   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -7.786E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.042E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.700E-04   TIME AVG. HEAT FLUX   3.831E-04
 LARGEST RESIDUAL HEAT FLUX        -1.377E-15   AT NODE        221   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.933E-02   AT NODE       6532   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.162E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.344    
 STEP TIME COMPLETED       3.440E+06,  TOTAL TIME COMPLETED        3.440E+06


  INCREMENT   345 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  341ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.273       TIME AVG. FORCE       0.398    
 LARGEST SCALED RESIDUAL FORCE     -5.543E-04   AT NODE      12059   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.543E-04
 LARGEST INCREMENT OF DISP.         3.673E-02   AT NODE        230   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.500E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.165E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.788E-04   TIME AVG. HEAT FLUX   3.846E-04
 LARGEST RESIDUAL HEAT FLUX        -1.393E-15   AT NODE       3212   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.937E-02   AT NODE        219   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.099E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.345    
 STEP TIME COMPLETED       3.450E+06,  TOTAL TIME COMPLETED        3.450E+06


  INCREMENT   346 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  292ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.267       TIME AVG. FORCE       0.398    
 LARGEST SCALED RESIDUAL FORCE     -7.755E-04   AT NODE      12061   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.755E-04
 LARGEST INCREMENT OF DISP.         3.675E-02   AT NODE        231   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -7.321E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.727E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.660E-04   TIME AVG. HEAT FLUX   3.860E-04
 LARGEST RESIDUAL HEAT FLUX         1.561E-15   AT NODE       1299   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.949E-02   AT NODE       6544   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.064E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.346    
 STEP TIME COMPLETED       3.460E+06,  TOTAL TIME COMPLETED        3.460E+06


  INCREMENT   347 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  340ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.260       TIME AVG. FORCE       0.398    
 LARGEST SCALED RESIDUAL FORCE     -8.516E-04   AT NODE      12061   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.516E-04
 LARGEST INCREMENT OF DISP.         3.707E-02   AT NODE       5636   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.059E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.797E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.727E-04   TIME AVG. HEAT FLUX   3.874E-04
 LARGEST RESIDUAL HEAT FLUX         1.301E-15   AT NODE          5   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.957E-02   AT NODE       6964   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.005E-02   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.347    
 STEP TIME COMPLETED       3.470E+06,  TOTAL TIME COMPLETED        3.470E+06


  INCREMENT   348 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  346ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.253       TIME AVG. FORCE       0.397    
 LARGEST SCALED RESIDUAL FORCE     -9.132E-04   AT NODE      12063   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.132E-04
 LARGEST INCREMENT OF DISP.         3.712E-02   AT NODE      12091   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.887E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.950E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.581E-04   TIME AVG. HEAT FLUX   3.887E-04
 LARGEST RESIDUAL HEAT FLUX         1.333E-15   AT NODE       3802   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.960E-02   AT NODE       6546   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.969E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.348    
 STEP TIME COMPLETED       3.480E+06,  TOTAL TIME COMPLETED        3.480E+06


  INCREMENT   349 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  334ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.247       TIME AVG. FORCE       0.397    
 LARGEST SCALED RESIDUAL FORCE     -9.609E-04   AT NODE      12065   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.609E-04
 LARGEST INCREMENT OF DISP.         3.747E-02   AT NODE        234   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.649E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.032E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.645E-04   TIME AVG. HEAT FLUX   3.901E-04
 LARGEST RESIDUAL HEAT FLUX         1.132E-15   AT NODE        401   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.969E-02   AT NODE        406   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.917E-02   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.349    
 STEP TIME COMPLETED       3.490E+06,  TOTAL TIME COMPLETED        3.490E+06


  INCREMENT   350 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.240       TIME AVG. FORCE       0.396    
 LARGEST SCALED RESIDUAL FORCE     -1.022E-03   AT NODE      12066   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.022E-03
 LARGEST INCREMENT OF DISP.         3.757E-02   AT NODE       5640   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.490E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.030E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.495E-04   TIME AVG. HEAT FLUX   3.914E-04
 LARGEST RESIDUAL HEAT FLUX        -1.463E-15   AT NODE        804   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.984E-02   AT NODE       6553   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.879E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.350    
 STEP TIME COMPLETED       3.500E+06,  TOTAL TIME COMPLETED        3.500E+06


  INCREMENT   351 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.233       TIME AVG. FORCE       0.396    
 LARGEST SCALED RESIDUAL FORCE     -1.081E-03   AT NODE      11280   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.081E-03
 LARGEST INCREMENT OF DISP.         3.795E-02   AT NODE        237   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.277E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.132E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.546E-04   TIME AVG. HEAT FLUX   3.928E-04
 LARGEST RESIDUAL HEAT FLUX         1.302E-15   AT NODE       4439   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.991E-02   AT NODE       6768   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.834E-02   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.351    
 STEP TIME COMPLETED       3.510E+06,  TOTAL TIME COMPLETED        3.510E+06


  INCREMENT   352 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.226       TIME AVG. FORCE       0.395    
 LARGEST SCALED RESIDUAL FORCE     -1.128E-03   AT NODE      11279   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.128E-03
 LARGEST INCREMENT OF DISP.         3.815E-02   AT NODE      12101   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.131E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.117E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.392E-04   TIME AVG. HEAT FLUX   3.940E-04
 LARGEST RESIDUAL HEAT FLUX         1.501E-15   AT NODE       1162   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.001E-02   AT NODE        408   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.797E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.352    
 STEP TIME COMPLETED       3.520E+06,  TOTAL TIME COMPLETED        3.520E+06


  INCREMENT   353 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  342ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.219       TIME AVG. FORCE       0.395    
 LARGEST SCALED RESIDUAL FORCE     -1.164E-03   AT NODE      12068   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.164E-03
 LARGEST INCREMENT OF DISP.         3.859E-02   AT NODE       6598   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.943E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.159E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.425E-04   TIME AVG. HEAT FLUX   3.953E-04
 LARGEST RESIDUAL HEAT FLUX        -1.261E-15   AT NODE        917   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.017E-02   AT NODE       5697   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.758E-02   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.353    
 STEP TIME COMPLETED       3.530E+06,  TOTAL TIME COMPLETED        3.530E+06


  INCREMENT   354 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  326ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.211       TIME AVG. FORCE       0.394    
 LARGEST SCALED RESIDUAL FORCE     -1.207E-03   AT NODE      12070   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.207E-03
 LARGEST INCREMENT OF DISP.         3.890E-02   AT NODE      12106   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.810E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.114E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.265E-04   TIME AVG. HEAT FLUX   3.965E-04
 LARGEST RESIDUAL HEAT FLUX         1.330E-15   AT NODE        553   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.030E-02   AT NODE       6561   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.723E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.354    
 STEP TIME COMPLETED       3.540E+06,  TOTAL TIME COMPLETED        3.540E+06


  INCREMENT   355 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  325ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.204       TIME AVG. FORCE       0.394    
 LARGEST SCALED RESIDUAL FORCE     -1.260E-03   AT NODE      12071   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.260E-03
 LARGEST INCREMENT OF DISP.         3.944E-02   AT NODE       6600   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.643E-03   AT NODE        297   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.207E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.284E-04   TIME AVG. HEAT FLUX   3.977E-04
 LARGEST RESIDUAL HEAT FLUX         1.338E-15   AT NODE       3614   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.044E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.690E-02   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.355    
 STEP TIME COMPLETED       3.550E+06,  TOTAL TIME COMPLETED        3.550E+06


  INCREMENT   356 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  287ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.196       TIME AVG. FORCE       0.393    
 LARGEST SCALED RESIDUAL FORCE     -1.295E-03   AT NODE      12073   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.295E-03
 LARGEST INCREMENT OF DISP.         3.990E-02   AT NODE      11770   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.744E-03   AT NODE       6799   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.225E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.106E-04   TIME AVG. HEAT FLUX   3.989E-04
 LARGEST RESIDUAL HEAT FLUX         1.187E-15   AT NODE       3921   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.068E-02   AT NODE       5636   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.657E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.356    
 STEP TIME COMPLETED       3.560E+06,  TOTAL TIME COMPLETED        3.560E+06


  INCREMENT   357 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  305ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.188       TIME AVG. FORCE       0.393    
 LARGEST SCALED RESIDUAL FORCE     -1.322E-03   AT NODE      12074   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.322E-03
 LARGEST INCREMENT OF DISP.         4.060E-02   AT NODE       6816   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.722E-03   AT NODE      12354   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.283E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.104E-04   TIME AVG. HEAT FLUX   4.001E-04
 LARGEST RESIDUAL HEAT FLUX         1.186E-15   AT NODE        909   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.082E-02   AT NODE       6572   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.630E-02   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.357    
 STEP TIME COMPLETED       3.570E+06,  TOTAL TIME COMPLETED        3.570E+06


  INCREMENT   358 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  339ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.180       TIME AVG. FORCE       0.392    
 LARGEST SCALED RESIDUAL FORCE     -1.366E-03   AT NODE      12079   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.366E-03
 LARGEST INCREMENT OF DISP.         4.130E-02   AT NODE      12114   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.892E-03   AT NODE       6580   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.367E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.908E-04   TIME AVG. HEAT FLUX   4.012E-04
 LARGEST RESIDUAL HEAT FLUX         1.488E-15   AT NODE        691   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.110E-02   AT NODE       6577   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.600E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.358    
 STEP TIME COMPLETED       3.580E+06,  TOTAL TIME COMPLETED        3.580E+06


  INCREMENT   359 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  339ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.172       TIME AVG. FORCE       0.392    
 LARGEST SCALED RESIDUAL FORCE     -1.402E-03   AT NODE      12081   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.402E-03
 LARGEST INCREMENT OF DISP.         4.229E-02   AT NODE      12118   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.915E-03   AT NODE       3218   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.464E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.878E-04   TIME AVG. HEAT FLUX   4.022E-04
 LARGEST RESIDUAL HEAT FLUX        -1.366E-15   AT NODE       1597   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.135E-02   AT NODE       5640   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.577E-02   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.359    
 STEP TIME COMPLETED       3.590E+06,  TOTAL TIME COMPLETED        3.590E+06


  INCREMENT   360 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.163       TIME AVG. FORCE       0.391    
 LARGEST SCALED RESIDUAL FORCE     -1.431E-03   AT NODE      12083   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.431E-03
 LARGEST INCREMENT OF DISP.         4.331E-02   AT NODE      12122   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.100E-03   AT NODE       6597   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.500E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.652E-04   TIME AVG. HEAT FLUX   4.032E-04
 LARGEST RESIDUAL HEAT FLUX         1.194E-15   AT NODE       4019   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.162E-02   AT NODE       6585   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.550E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.360    
 STEP TIME COMPLETED       3.600E+06,  TOTAL TIME COMPLETED        3.600E+06


  INCREMENT   361 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  390ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.154       TIME AVG. FORCE       0.390    
 LARGEST SCALED RESIDUAL FORCE     -1.462E-03   AT NODE      12089   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.462E-03
 LARGEST INCREMENT OF DISP.         4.477E-02   AT NODE      12123   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.187E-03   AT NODE       3273   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.621E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.588E-04   TIME AVG. HEAT FLUX   4.042E-04
 LARGEST RESIDUAL HEAT FLUX        -1.347E-15   AT NODE       3803   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.195E-02   AT NODE       6583   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.534E-02   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.361    
 STEP TIME COMPLETED       3.610E+06,  TOTAL TIME COMPLETED        3.610E+06


  INCREMENT   362 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  391ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.144       TIME AVG. FORCE       0.390    
 LARGEST SCALED RESIDUAL FORCE     -1.500E-03   AT NODE      11789   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.500E-03
 LARGEST INCREMENT OF DISP.         4.642E-02   AT NODE      12128   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.424E-03   AT NODE       6599   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.698E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.330E-04   TIME AVG. HEAT FLUX   4.051E-04
 LARGEST RESIDUAL HEAT FLUX         1.186E-15   AT NODE       3682   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.228E-02   AT NODE       6587   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.508E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.362    
 STEP TIME COMPLETED       3.620E+06,  TOTAL TIME COMPLETED        3.620E+06


  INCREMENT   363 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.134       TIME AVG. FORCE       0.389    
 LARGEST SCALED RESIDUAL FORCE     -1.537E-03   AT NODE      11788   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.537E-03
 LARGEST INCREMENT OF DISP.         4.863E-02   AT NODE       6631   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.592E-03   AT NODE       3221   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.777E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.218E-04   TIME AVG. HEAT FLUX   4.060E-04
 LARGEST RESIDUAL HEAT FLUX         1.456E-15   AT NODE        691   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.260E-02   AT NODE       6591   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.500E-02   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.363    
 STEP TIME COMPLETED       3.630E+06,  TOTAL TIME COMPLETED        3.630E+06


  INCREMENT   364 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.123       TIME AVG. FORCE       0.388    
 LARGEST SCALED RESIDUAL FORCE     -1.575E-03   AT NODE      12091   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.575E-03
 LARGEST INCREMENT OF DISP.         5.126E-02   AT NODE      12134   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.940E-03   AT NODE       6812   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.943E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.915E-04   TIME AVG. HEAT FLUX   4.068E-04
 LARGEST RESIDUAL HEAT FLUX        -1.248E-15   AT NODE        814   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.296E-02   AT NODE       6807   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.476E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.364    
 STEP TIME COMPLETED       3.640E+06,  TOTAL TIME COMPLETED        3.640E+06


  INCREMENT   365 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  371ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.111       TIME AVG. FORCE       0.387    
 LARGEST SCALED RESIDUAL FORCE     -1.603E-03   AT NODE      12094   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.603E-03
 LARGEST INCREMENT OF DISP.         5.467E-02   AT NODE      12138   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -7.248E-03   AT NODE      12371   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.007E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.756E-04   TIME AVG. HEAT FLUX   4.075E-04
 LARGEST RESIDUAL HEAT FLUX        -1.218E-15   AT NODE        347   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.318E-02   AT NODE       6816   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.478E-02   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.365    
 STEP TIME COMPLETED       3.650E+06,  TOTAL TIME COMPLETED        3.650E+06


  INCREMENT   366 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  350ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.780E-02   TIME AVG. FORCE       0.387    
 LARGEST SCALED RESIDUAL FORCE     -1.635E-03   AT NODE      12099   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.635E-03
 LARGEST INCREMENT OF DISP.         5.860E-02   AT NODE       6842   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -7.797E-03   AT NODE       6611   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.151E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.385E-04   TIME AVG. HEAT FLUX   4.082E-04
 LARGEST RESIDUAL HEAT FLUX         1.239E-15   AT NODE       4037   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.328E-02   AT NODE        246   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.455E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.366    
 STEP TIME COMPLETED       3.660E+06,  TOTAL TIME COMPLETED        3.660E+06


  INCREMENT   367 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  340ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.363E-02   TIME AVG. FORCE       0.386    
 LARGEST SCALED RESIDUAL FORCE     -1.675E-03   AT NODE      12101   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.675E-03
 LARGEST INCREMENT OF DISP.         6.302E-02   AT NODE        264   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -8.417E-03   AT NODE       6613   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.273E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.162E-04   TIME AVG. HEAT FLUX   4.087E-04
 LARGEST RESIDUAL HEAT FLUX        -1.337E-15   AT NODE      12450   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.309E-02   AT NODE        248   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.477E-02   AT NODE       6706   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.367    
 STEP TIME COMPLETED       3.670E+06,  TOTAL TIME COMPLETED        3.670E+06


  INCREMENT   368 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  332ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.819E-02   TIME AVG. FORCE       0.385    
 LARGEST SCALED RESIDUAL FORCE     -1.681E-03   AT NODE      12104   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.681E-03
 LARGEST INCREMENT OF DISP.         6.665E-02   AT NODE        266   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -9.465E-03   AT NODE        253   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -3.140E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.731E-04   TIME AVG. HEAT FLUX   4.092E-04
 LARGEST RESIDUAL HEAT FLUX        -1.261E-15   AT NODE      13297   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.246E-02   AT NODE       6613   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.451E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.368    
 STEP TIME COMPLETED       3.680E+06,  TOTAL TIME COMPLETED        3.680E+06


  INCREMENT   369 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  338ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.176E-02   TIME AVG. FORCE       0.384    
 LARGEST SCALED RESIDUAL FORCE     -1.664E-03   AT NODE      12107   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.664E-03
 LARGEST INCREMENT OF DISP.         6.793E-02   AT NODE       6699   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.078E-02   AT NODE      11100   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -3.119E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.435E-04   TIME AVG. HEAT FLUX   4.096E-04
 LARGEST RESIDUAL HEAT FLUX        -1.095E-15   AT NODE       2461   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.104E-02   AT NODE       6824   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.487E-02   AT NODE       6706   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.369    
 STEP TIME COMPLETED       3.690E+06,  TOTAL TIME COMPLETED        3.690E+06


  INCREMENT   370 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  318ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.512E-02   TIME AVG. FORCE       0.383    
 LARGEST SCALED RESIDUAL FORCE     -1.660E-03   AT NODE      12110   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.660E-03
 LARGEST INCREMENT OF DISP.         6.399E-02   AT NODE       6658   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.289E-02   AT NODE        258   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -3.112E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.916E-04   TIME AVG. HEAT FLUX   4.098E-04
 LARGEST RESIDUAL HEAT FLUX        -1.507E-15   AT NODE        456   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.851E-02   AT NODE        255   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.458E-02   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.370    
 STEP TIME COMPLETED       3.700E+06,  TOTAL TIME COMPLETED        3.700E+06


  INCREMENT   371 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  335ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.032E-02   TIME AVG. FORCE       0.382    
 LARGEST SCALED RESIDUAL FORCE     -1.652E-03   AT NODE      11769   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.652E-03
 LARGEST INCREMENT OF DISP.         5.333E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.553E-02   AT NODE      12148   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -3.535E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.537E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX         1.162E-15   AT NODE       4337   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.490E-02   AT NODE       6635   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.405E-02   AT NODE       6706   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  0.371    
 STEP TIME COMPLETED       3.710E+06,  TOTAL TIME COMPLETED        3.710E+06


  INCREMENT   372 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  350ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.707E-03   AT NODE      12120   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.707E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.760E-02   AT NODE       7098   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX        -1.229E-15   AT NODE       4354   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.458E-02   AT NODE       6990   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  303ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.659E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.659E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.968E-03   AT NODE      12138   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX        -1.283E-15   AT NODE       5362   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.672E-15   AT NODE       8642   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  347ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.628E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.628E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -3.199E-03   AT NODE      12134   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX        -1.276E-15   AT NODE       5193   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.818E-15   AT NODE       2186   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  331ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.606E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.606E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -2.309E-03   AT NODE      12132   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX         1.381E-15   AT NODE       4629   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.779E-15   AT NODE      11985   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  351ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.590E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.590E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.800E-03   AT NODE      11100   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX        -1.549E-15   AT NODE       4041   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.222E-14   AT NODE       2847   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  347ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.577E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.577E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.422E-03   AT NODE       6631   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX         1.413E-15   AT NODE       3143   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.647E-15   AT NODE       3818   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  323ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.566E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.566E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.221E-03   AT NODE      11101   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX        -1.450E-15   AT NODE      13537   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.497E-15   AT NODE        192   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  347ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.557E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.557E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.045E-03   AT NODE       6830   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX        -1.164E-15   AT NODE      11986   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.487E-15   AT NODE        608   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  341ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.548E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.548E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -8.891E-04   AT NODE       6830   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX        -1.136E-15   AT NODE       4611   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.210E-14   AT NODE       6706   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  315ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    10


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.539E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.539E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -7.683E-04   AT NODE        254   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX         1.146E-15   AT NODE       5302   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.512E-14   AT NODE       6706   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    11


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.530E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.530E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -7.237E-04   AT NODE        254   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX         1.402E-15   AT NODE       2930   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.081E-14   AT NODE       4208   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  377ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    12


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.521E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.521E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.640E-04   AT NODE      12130   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX        -1.302E-15   AT NODE       3727   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.786E-15   AT NODE       2793   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  385ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    13


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.513E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.513E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.759E-04   AT NODE        254   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX         1.137E-15   AT NODE       4643   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.606E-15   AT NODE       3594   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    14


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.504E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.504E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.302E-04   AT NODE        254   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX        -1.308E-15   AT NODE       3555   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.238E-15   AT NODE       6706   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    15


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.495E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.495E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.125E-04   AT NODE       6627   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX         1.207E-15   AT NODE       3555   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.345E-15   AT NODE      14485   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  394ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    16


 AVERAGE FORCE                      4.462E-03   TIME AVG. FORCE       0.381    
 LARGEST SCALED RESIDUAL FORCE     -1.485E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.485E-03
 LARGEST INCREMENT OF DISP.         3.732E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.954E-04   AT NODE       6627   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  3.969E-04   TIME AVG. HEAT FLUX   4.099E-04
 LARGEST RESIDUAL HEAT FLUX         1.071E-15   AT NODE       5228   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.147E-02   AT NODE        336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.978E-15   AT NODE       4573   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 

 ***NOTE: THE SOLUTION FAILS TO CONVERGE IN THE MAXIMUM NUMBER OF EQUILIBRIUM 
          ITERATIONS ALLOWED.
 

 ***ERROR: FIXED TIME INCREMENT IS TOO LARGE

 ***ERROR: THE ANALYSIS HAS BEEN TERMINATED DUE TO PREVIOUS ERRORS. ALL OUTPUT 
           REQUESTS HAVE BEEN WRITTEN FOR THE LAST CONVERGED INCREMENT.



     ANALYSIS SUMMARY:
     TOTAL OF        372  INCREMENTS
                       1  CUTBACKS IN AUTOMATIC INCREMENTATION
                     435  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                     435  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                     351  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       5  WARNING MESSAGES DURING USER INPUT PROCESSING
                       1  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                       0  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       2  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =     4.15E+02
       SYSTEM TIME (SEC)    =      22.    
       TOTAL CPU TIME (SEC) =     4.36E+02
       WALLCLOCK TIME (SEC) =          368
       MEMORY PEAK (GB)     =            0
