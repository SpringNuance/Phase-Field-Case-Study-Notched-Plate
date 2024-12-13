1

   Abaqus 2023.HF4                                  Date 31-Aug-2024   Time 18:34:30
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

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  350ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  29ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  366ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  362ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  25ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  357ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  339ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  26ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  329ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  332ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  27ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  375ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  332ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  26ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  328ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  338ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  26ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  340ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  338ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  26ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  330ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  356ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  25ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  340ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  347ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  27ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  378ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  351ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  25ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  348ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  340ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  24ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  345ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  355ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  30ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  357ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  337ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  347ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  360ms

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
 	1 HOST x 1 MPI RANK PER HOST x 1 THREAD PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  26ms

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
