1

   Abaqus 2023.HF4                                  Date 12-Aug-2024   Time 21:12:29
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
                              1            5280
                              2            5528
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  186ms

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

 AVERAGE HEAT FLUX                  2.417E-16   TIME AVG. HEAT FLUX   1.000E-02
 LARGEST RESIDUAL HEAT FLUX         2.165E-15   AT NODE       3583   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.762E-15   AT NODE       1335   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        9.762E-15   AT NODE       1335   DOF 11
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
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.389E-03   TIME AVG. FORCE       3.830E-03
 LARGEST SCALED RESIDUAL FORCE      5.065E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.065E-11
 LARGEST INCREMENT OF DISP.         1.065E-05   AT NODE         38   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.671E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.160E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.960E-06   TIME AVG. HEAT FLUX   3.960E-06
 LARGEST RESIDUAL HEAT FLUX         1.295E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.114E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.114E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  204ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.389E-03   TIME AVG. FORCE       3.830E-03
 LARGEST SCALED RESIDUAL FORCE      1.640E-13   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.640E-13
 LARGEST INCREMENT OF DISP.         1.065E-05   AT NODE         38   DOF  2
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.945E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.958E-06   TIME AVG. HEAT FLUX   3.958E-06
 LARGEST RESIDUAL HEAT FLUX        -2.897E-15   AT NODE      11986   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.120E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.035E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  163ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.583E-03   TIME AVG. FORCE       5.748E-03
 LARGEST SCALED RESIDUAL FORCE      1.306E-12   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.306E-12
 LARGEST INCREMENT OF DISP.         1.703E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.136E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  7.957E-06   TIME AVG. HEAT FLUX   5.958E-06
 LARGEST RESIDUAL HEAT FLUX        -2.527E-15   AT NODE       1023   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.132E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.215E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  19ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.278E-02   TIME AVG. FORCE       7.505E-03
 LARGEST SCALED RESIDUAL FORCE      1.037E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.037E-10
 LARGEST INCREMENT OF DISP.         2.840E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.137E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.341E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.195E-05   TIME AVG. HEAT FLUX   7.954E-06
 LARGEST RESIDUAL HEAT FLUX         1.118E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.144E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.216E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  175ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.278E-02   TIME AVG. FORCE       7.505E-03
 LARGEST SCALED RESIDUAL FORCE      6.853E-15   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.853E-15
 LARGEST INCREMENT OF DISP.         2.841E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.793E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.195E-05   TIME AVG. HEAT FLUX   7.954E-06
 LARGEST RESIDUAL HEAT FLUX        -2.503E-15   AT NODE       5411   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.144E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.243E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  217ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.597E-02   TIME AVG. FORCE       9.198E-03
 LARGEST SCALED RESIDUAL FORCE      5.111E-12   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.111E-12
 LARGEST INCREMENT OF DISP.         3.980E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.139E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.594E-05   TIME AVG. HEAT FLUX   9.952E-06
 LARGEST RESIDUAL HEAT FLUX         2.723E-15   AT NODE       3730   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.156E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.229E-04   AT NODE          8   DOF 11
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
 LARGEST SCALED RESIDUAL FORCE      1.178E-10   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.178E-10
 LARGEST INCREMENT OF DISP.         5.121E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.141E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.943E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.995E-05   TIME AVG. HEAT FLUX   1.195E-05
 LARGEST RESIDUAL HEAT FLUX         1.129E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.169E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.233E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  186ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.916E-02   TIME AVG. FORCE       1.086E-02
 LARGEST SCALED RESIDUAL FORCE      1.869E-14   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.869E-14
 LARGEST INCREMENT OF DISP.         5.122E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.753E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.995E-05   TIME AVG. HEAT FLUX   1.195E-05
 LARGEST RESIDUAL HEAT FLUX         2.882E-15   AT NODE        322   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.169E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.292E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  196ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.236E-02   TIME AVG. FORCE       1.250E-02
 LARGEST SCALED RESIDUAL FORCE      1.131E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.131E-11
 LARGEST INCREMENT OF DISP.         6.264E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.143E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.396E-05   TIME AVG. HEAT FLUX   1.395E-05
 LARGEST RESIDUAL HEAT FLUX        -3.012E-15   AT NODE        747   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.181E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.246E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  14ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.555E-02   TIME AVG. FORCE       1.413E-02
 LARGEST SCALED RESIDUAL FORCE      2.860E-10   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.860E-10
 LARGEST INCREMENT OF DISP.         7.409E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.144E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.094E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.797E-05   TIME AVG. HEAT FLUX   1.595E-05
 LARGEST RESIDUAL HEAT FLUX         1.138E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.194E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.249E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  195ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.555E-02   TIME AVG. FORCE       1.413E-02
 LARGEST SCALED RESIDUAL FORCE      3.592E-14   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.592E-14
 LARGEST INCREMENT OF DISP.         7.409E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.782E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.797E-05   TIME AVG. HEAT FLUX   1.595E-05
 LARGEST RESIDUAL HEAT FLUX         2.558E-15   AT NODE       5400   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.194E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.336E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  211ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.874E-02   TIME AVG. FORCE       1.576E-02
 LARGEST SCALED RESIDUAL FORCE      1.975E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.975E-11
 LARGEST INCREMENT OF DISP.         8.556E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.147E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.200E-05   TIME AVG. HEAT FLUX   1.796E-05
 LARGEST RESIDUAL HEAT FLUX         2.994E-15   AT NODE       5035   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.206E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.262E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  15ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.193E-02   TIME AVG. FORCE       1.737E-02
 LARGEST SCALED RESIDUAL FORCE      5.275E-10   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.275E-10
 LARGEST INCREMENT OF DISP.         9.704E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.148E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.497E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.602E-05   TIME AVG. HEAT FLUX   1.997E-05
 LARGEST RESIDUAL HEAT FLUX         1.147E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.219E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.265E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  187ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.193E-02   TIME AVG. FORCE       1.737E-02
 LARGEST SCALED RESIDUAL FORCE      5.818E-14   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.818E-14
 LARGEST INCREMENT OF DISP.         9.705E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.878E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.603E-05   TIME AVG. HEAT FLUX   1.997E-05
 LARGEST RESIDUAL HEAT FLUX        -2.697E-15   AT NODE       3383   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.219E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.378E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.513E-02   TIME AVG. FORCE       1.899E-02
 LARGEST SCALED RESIDUAL FORCE      3.032E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.032E-11
 LARGEST INCREMENT OF DISP.         1.086E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.150E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.007E-05   TIME AVG. HEAT FLUX   2.198E-05
 LARGEST RESIDUAL HEAT FLUX         2.828E-15   AT NODE       2258   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.232E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.278E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.832E-02   TIME AVG. FORCE       2.060E-02
 LARGEST SCALED RESIDUAL FORCE      8.426E-10   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      8.426E-10
 LARGEST INCREMENT OF DISP.         1.201E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.152E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.903E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.411E-05   TIME AVG. HEAT FLUX   2.399E-05
 LARGEST RESIDUAL HEAT FLUX         1.156E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.245E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.281E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.832E-02   TIME AVG. FORCE       2.060E-02
 LARGEST SCALED RESIDUAL FORCE      8.509E-14   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      8.509E-14
 LARGEST INCREMENT OF DISP.         1.201E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.035E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  4.411E-05   TIME AVG. HEAT FLUX   2.399E-05
 LARGEST RESIDUAL HEAT FLUX         2.719E-15   AT NODE       2697   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.245E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.418E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.151E-02   TIME AVG. FORCE       2.221E-02
 LARGEST SCALED RESIDUAL FORCE      4.287E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.287E-11
 LARGEST INCREMENT OF DISP.         1.316E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.154E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.816E-05   TIME AVG. HEAT FLUX   2.600E-05
 LARGEST RESIDUAL HEAT FLUX        -3.140E-15   AT NODE       3036   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.258E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.293E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.469E-02   TIME AVG. FORCE       2.381E-02
 LARGEST SCALED RESIDUAL FORCE      1.231E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.231E-09
 LARGEST INCREMENT OF DISP.         1.432E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.156E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.311E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.222E-05   TIME AVG. HEAT FLUX   2.802E-05
 LARGEST RESIDUAL HEAT FLUX         1.163E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.271E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.297E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  356ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.469E-02   TIME AVG. FORCE       2.381E-02
 LARGEST SCALED RESIDUAL FORCE      1.163E-13   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.163E-13
 LARGEST INCREMENT OF DISP.         1.432E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.251E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  5.222E-05   TIME AVG. HEAT FLUX   2.802E-05
 LARGEST RESIDUAL HEAT FLUX         2.728E-15   AT NODE       2512   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.271E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.456E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.788E-02   TIME AVG. FORCE       2.542E-02
 LARGEST SCALED RESIDUAL FORCE      5.726E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.726E-11
 LARGEST INCREMENT OF DISP.         1.548E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.158E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  5.628E-05   TIME AVG. HEAT FLUX   3.004E-05
 LARGEST RESIDUAL HEAT FLUX        -3.151E-15   AT NODE       4517   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.284E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.309E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  24ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.107E-02   TIME AVG. FORCE       2.702E-02
 LARGEST SCALED RESIDUAL FORCE      1.694E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.694E-09
 LARGEST INCREMENT OF DISP.         1.664E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.160E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.723E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.035E-05   TIME AVG. HEAT FLUX   3.206E-05
 LARGEST RESIDUAL HEAT FLUX         1.170E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.297E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.312E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.107E-02   TIME AVG. FORCE       2.702E-02
 LARGEST SCALED RESIDUAL FORCE      1.519E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.519E-13
 LARGEST INCREMENT OF DISP.         1.664E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.522E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  6.035E-05   TIME AVG. HEAT FLUX   3.206E-05
 LARGEST RESIDUAL HEAT FLUX         3.329E-15   AT NODE       2945   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.297E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.491E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.425E-02   TIME AVG. FORCE       2.862E-02
 LARGEST SCALED RESIDUAL FORCE      7.336E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      7.336E-11
 LARGEST INCREMENT OF DISP.         1.780E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.162E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  6.442E-05   TIME AVG. HEAT FLUX   3.408E-05
 LARGEST RESIDUAL HEAT FLUX        -3.040E-15   AT NODE       5212   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.310E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.324E-04   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      5.744E-02   TIME AVG. FORCE       3.022E-02
 LARGEST SCALED RESIDUAL FORCE      2.231E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.231E-09
 LARGEST INCREMENT OF DISP.         1.896E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.164E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.138E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.851E-05   TIME AVG. HEAT FLUX   3.611E-05
 LARGEST RESIDUAL HEAT FLUX         1.177E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.323E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.326E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.744E-02   TIME AVG. FORCE       3.022E-02
 LARGEST SCALED RESIDUAL FORCE      1.925E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.925E-13
 LARGEST INCREMENT OF DISP.         1.897E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.843E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  6.851E-05   TIME AVG. HEAT FLUX   3.611E-05
 LARGEST RESIDUAL HEAT FLUX         3.009E-15   AT NODE       4430   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.323E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.524E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  372ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.062E-02   TIME AVG. FORCE       3.182E-02
 LARGEST SCALED RESIDUAL FORCE      9.103E-11   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      9.103E-11
 LARGEST INCREMENT OF DISP.         2.013E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.166E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  7.260E-05   TIME AVG. HEAT FLUX   3.813E-05
 LARGEST RESIDUAL HEAT FLUX         3.160E-15   AT NODE       4474   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.337E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.338E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.380E-02   TIME AVG. FORCE       3.342E-02
 LARGEST SCALED RESIDUAL FORCE      2.842E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.842E-09
 LARGEST INCREMENT OF DISP.         2.130E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.169E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.556E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.668E-05   TIME AVG. HEAT FLUX   4.016E-05
 LARGEST RESIDUAL HEAT FLUX         1.183E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.350E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.341E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.380E-02   TIME AVG. FORCE       3.342E-02
 LARGEST SCALED RESIDUAL FORCE      2.371E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.371E-13
 LARGEST INCREMENT OF DISP.         2.130E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.212E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  7.668E-05   TIME AVG. HEAT FLUX   4.016E-05
 LARGEST RESIDUAL HEAT FLUX         2.518E-15   AT NODE       1831   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.350E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.554E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.698E-02   TIME AVG. FORCE       3.502E-02
 LARGEST SCALED RESIDUAL FORCE      1.101E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.101E-10
 LARGEST INCREMENT OF DISP.         2.247E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.171E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  8.079E-05   TIME AVG. HEAT FLUX   4.219E-05
 LARGEST RESIDUAL HEAT FLUX         2.751E-15   AT NODE        989   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.364E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.352E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  24ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.015E-02   TIME AVG. FORCE       3.662E-02
 LARGEST SCALED RESIDUAL FORCE      3.527E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.527E-09
 LARGEST INCREMENT OF DISP.         2.364E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.173E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.978E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.490E-05   TIME AVG. HEAT FLUX   4.423E-05
 LARGEST RESIDUAL HEAT FLUX         1.188E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.377E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.354E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.015E-02   TIME AVG. FORCE       3.662E-02
 LARGEST SCALED RESIDUAL FORCE      2.855E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.855E-13
 LARGEST INCREMENT OF DISP.         2.365E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.624E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  8.490E-05   TIME AVG. HEAT FLUX   4.423E-05
 LARGEST RESIDUAL HEAT FLUX        -2.488E-15   AT NODE       4598   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.377E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.581E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.333E-02   TIME AVG. FORCE       3.821E-02
 LARGEST SCALED RESIDUAL FORCE      1.305E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.305E-10
 LARGEST INCREMENT OF DISP.         2.482E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.175E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  8.900E-05   TIME AVG. HEAT FLUX   4.626E-05
 LARGEST RESIDUAL HEAT FLUX         3.103E-15   AT NODE       5298   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.391E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.366E-04   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      7.650E-02   TIME AVG. FORCE       3.981E-02
 LARGEST SCALED RESIDUAL FORCE      4.289E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.289E-09
 LARGEST INCREMENT OF DISP.         2.600E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.178E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.404E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  9.311E-05   TIME AVG. HEAT FLUX   4.830E-05
 LARGEST RESIDUAL HEAT FLUX         1.192E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.405E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.368E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.650E-02   TIME AVG. FORCE       3.981E-02
 LARGEST SCALED RESIDUAL FORCE      3.372E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.372E-13
 LARGEST INCREMENT OF DISP.         2.600E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.075E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  9.311E-05   TIME AVG. HEAT FLUX   4.830E-05
 LARGEST RESIDUAL HEAT FLUX        -2.883E-15   AT NODE       5478   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.405E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.606E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.967E-02   TIME AVG. FORCE       4.140E-02
 LARGEST SCALED RESIDUAL FORCE      1.520E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.520E-10
 LARGEST INCREMENT OF DISP.         2.718E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.179E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  9.725E-05   TIME AVG. HEAT FLUX   5.034E-05
 LARGEST RESIDUAL HEAT FLUX         2.891E-15   AT NODE       1500   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.419E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.379E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.285E-02   TIME AVG. FORCE       4.300E-02
 LARGEST SCALED RESIDUAL FORCE      5.125E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.125E-09
 LARGEST INCREMENT OF DISP.         2.836E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.183E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.833E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.014E-04   TIME AVG. HEAT FLUX   5.238E-05
 LARGEST RESIDUAL HEAT FLUX         1.195E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.432E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.380E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  380ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      8.285E-02   TIME AVG. FORCE       4.300E-02
 LARGEST SCALED RESIDUAL FORCE      3.921E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.921E-13
 LARGEST INCREMENT OF DISP.         2.836E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.561E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.014E-04   TIME AVG. HEAT FLUX   5.238E-05
 LARGEST RESIDUAL HEAT FLUX         2.705E-15   AT NODE       8403   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.432E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.626E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.601E-02   TIME AVG. FORCE       4.459E-02
 LARGEST SCALED RESIDUAL FORCE      1.746E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.746E-10
 LARGEST INCREMENT OF DISP.         2.955E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.184E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.055E-04   TIME AVG. HEAT FLUX   5.443E-05
 LARGEST RESIDUAL HEAT FLUX        -2.602E-15   AT NODE       5299   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.446E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.391E-04   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      8.918E-02   TIME AVG. FORCE       4.618E-02
 LARGEST SCALED RESIDUAL FORCE      6.038E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.038E-09
 LARGEST INCREMENT OF DISP.         3.073E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.188E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.267E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.097E-04   TIME AVG. HEAT FLUX   5.647E-05
 LARGEST RESIDUAL HEAT FLUX         1.198E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.460E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.392E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  375ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      8.918E-02   TIME AVG. FORCE       4.618E-02
 LARGEST SCALED RESIDUAL FORCE      4.498E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.498E-13
 LARGEST INCREMENT OF DISP.         3.073E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.078E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.097E-04   TIME AVG. HEAT FLUX   5.647E-05
 LARGEST RESIDUAL HEAT FLUX         2.877E-15   AT NODE       3786   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.460E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.644E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.234E-02   TIME AVG. FORCE       4.777E-02
 LARGEST SCALED RESIDUAL FORCE      1.980E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.980E-10
 LARGEST INCREMENT OF DISP.         3.192E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.189E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.139E-04   TIME AVG. HEAT FLUX   5.852E-05
 LARGEST RESIDUAL HEAT FLUX         3.168E-15   AT NODE       2020   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.474E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.403E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.550E-02   TIME AVG. FORCE       4.936E-02
 LARGEST SCALED RESIDUAL FORCE      7.028E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      7.028E-09
 LARGEST INCREMENT OF DISP.         3.312E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.193E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.705E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.180E-04   TIME AVG. HEAT FLUX   6.057E-05
 LARGEST RESIDUAL HEAT FLUX         1.200E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.488E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.404E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      9.550E-02   TIME AVG. FORCE       4.936E-02
 LARGEST SCALED RESIDUAL FORCE      5.099E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.099E-13
 LARGEST INCREMENT OF DISP.         3.312E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.621E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.180E-04   TIME AVG. HEAT FLUX   6.057E-05
 LARGEST RESIDUAL HEAT FLUX         2.771E-15   AT NODE       4342   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.488E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.658E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.866E-02   TIME AVG. FORCE       5.095E-02
 LARGEST SCALED RESIDUAL FORCE      2.221E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.221E-10
 LARGEST INCREMENT OF DISP.         3.431E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.193E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.222E-04   TIME AVG. HEAT FLUX   6.263E-05
 LARGEST RESIDUAL HEAT FLUX        -2.252E-15   AT NODE       2085   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.503E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.414E-04   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.102       TIME AVG. FORCE       5.254E-02
 LARGEST SCALED RESIDUAL FORCE      8.096E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      8.096E-09
 LARGEST INCREMENT OF DISP.         3.551E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.198E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.147E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.264E-04   TIME AVG. HEAT FLUX   6.468E-05
 LARGEST RESIDUAL HEAT FLUX         1.201E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.517E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.415E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.102       TIME AVG. FORCE       5.254E-02
 LARGEST SCALED RESIDUAL FORCE      5.721E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.721E-13
 LARGEST INCREMENT OF DISP.         3.551E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.019E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.264E-04   TIME AVG. HEAT FLUX   6.468E-05
 LARGEST RESIDUAL HEAT FLUX        -2.763E-15   AT NODE        603   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.517E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.669E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  368ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.105       TIME AVG. FORCE       5.413E-02
 LARGEST SCALED RESIDUAL FORCE      2.467E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.467E-10
 LARGEST INCREMENT OF DISP.         3.671E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.198E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.305E-04   TIME AVG. HEAT FLUX   6.674E-05
 LARGEST RESIDUAL HEAT FLUX         2.686E-15   AT NODE       4978   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.531E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.424E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.108       TIME AVG. FORCE       5.572E-02
 LARGEST SCALED RESIDUAL FORCE      9.241E-09   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      9.241E-09
 LARGEST INCREMENT OF DISP.         3.791E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.203E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.593E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.347E-04   TIME AVG. HEAT FLUX   6.880E-05
 LARGEST RESIDUAL HEAT FLUX         1.201E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.545E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.425E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  393ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.108       TIME AVG. FORCE       5.572E-02
 LARGEST SCALED RESIDUAL FORCE      6.362E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.362E-13
 LARGEST INCREMENT OF DISP.         3.791E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.077E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.347E-04   TIME AVG. HEAT FLUX   6.880E-05
 LARGEST RESIDUAL HEAT FLUX         2.675E-15   AT NODE       3482   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.545E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.675E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.111       TIME AVG. FORCE       5.731E-02
 LARGEST SCALED RESIDUAL FORCE      2.718E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.718E-10
 LARGEST INCREMENT OF DISP.         3.911E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.203E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.389E-04   TIME AVG. HEAT FLUX   7.087E-05
 LARGEST RESIDUAL HEAT FLUX         3.168E-15   AT NODE       2626   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.560E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.434E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.114       TIME AVG. FORCE       5.889E-02
 LARGEST SCALED RESIDUAL FORCE      1.047E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.047E-08
 LARGEST INCREMENT OF DISP.         4.032E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.208E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.045E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.431E-04   TIME AVG. HEAT FLUX   7.293E-05
 LARGEST RESIDUAL HEAT FLUX         1.200E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.574E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.434E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.114       TIME AVG. FORCE       5.889E-02
 LARGEST SCALED RESIDUAL FORCE      7.016E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      7.016E-13
 LARGEST INCREMENT OF DISP.         4.032E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.137E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.432E-04   TIME AVG. HEAT FLUX   7.293E-05
 LARGEST RESIDUAL HEAT FLUX         2.528E-15   AT NODE        386   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.574E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.678E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.118       TIME AVG. FORCE       6.048E-02
 LARGEST SCALED RESIDUAL FORCE      2.972E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.972E-10
 LARGEST INCREMENT OF DISP.         4.153E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.209E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.473E-04   TIME AVG. HEAT FLUX   7.500E-05
 LARGEST RESIDUAL HEAT FLUX         2.644E-15   AT NODE       3786   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.588E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.443E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.121       TIME AVG. FORCE       6.206E-02
 LARGEST SCALED RESIDUAL FORCE      1.177E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.177E-08
 LARGEST INCREMENT OF DISP.         4.275E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.214E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.501E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.516E-04   TIME AVG. HEAT FLUX   7.707E-05
 LARGEST RESIDUAL HEAT FLUX         1.198E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.603E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.442E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  383ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.121       TIME AVG. FORCE       6.206E-02
 LARGEST SCALED RESIDUAL FORCE      7.680E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      7.680E-13
 LARGEST INCREMENT OF DISP.         4.275E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.197E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.516E-04   TIME AVG. HEAT FLUX   7.707E-05
 LARGEST RESIDUAL HEAT FLUX         2.886E-15   AT NODE       3838   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.603E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.677E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.124       TIME AVG. FORCE       6.364E-02
 LARGEST SCALED RESIDUAL FORCE      3.227E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.227E-10
 LARGEST INCREMENT OF DISP.         4.396E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.214E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.558E-04   TIME AVG. HEAT FLUX   7.914E-05
 LARGEST RESIDUAL HEAT FLUX         2.661E-15   AT NODE       3921   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.617E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.450E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.127       TIME AVG. FORCE       6.523E-02
 LARGEST SCALED RESIDUAL FORCE      1.316E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.316E-08
 LARGEST INCREMENT OF DISP.         4.518E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.220E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.963E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.600E-04   TIME AVG. HEAT FLUX   8.121E-05
 LARGEST RESIDUAL HEAT FLUX         1.195E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.632E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.450E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  383ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.127       TIME AVG. FORCE       6.523E-02
 LARGEST SCALED RESIDUAL FORCE      8.351E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      8.351E-13
 LARGEST INCREMENT OF DISP.         4.518E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.258E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.600E-04   TIME AVG. HEAT FLUX   8.121E-05
 LARGEST RESIDUAL HEAT FLUX         3.015E-15   AT NODE       3964   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.632E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.671E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  397ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.130       TIME AVG. FORCE       6.681E-02
 LARGEST SCALED RESIDUAL FORCE      3.482E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.482E-10
 LARGEST INCREMENT OF DISP.         4.640E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.219E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.642E-04   TIME AVG. HEAT FLUX   8.329E-05
 LARGEST RESIDUAL HEAT FLUX         2.835E-15   AT NODE       4423   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.647E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.457E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       6.839E-02
 LARGEST SCALED RESIDUAL FORCE      1.462E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.462E-08
 LARGEST INCREMENT OF DISP.         4.763E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.225E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    8.430E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.684E-04   TIME AVG. HEAT FLUX   8.537E-05
 LARGEST RESIDUAL HEAT FLUX         1.191E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.661E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.457E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  386ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.133       TIME AVG. FORCE       6.839E-02
 LARGEST SCALED RESIDUAL FORCE      9.024E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      9.024E-13
 LARGEST INCREMENT OF DISP.         4.763E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.319E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.685E-04   TIME AVG. HEAT FLUX   8.537E-05
 LARGEST RESIDUAL HEAT FLUX         3.033E-15   AT NODE       2600   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.661E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.661E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.136       TIME AVG. FORCE       6.997E-02
 LARGEST SCALED RESIDUAL FORCE      3.736E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.736E-10
 LARGEST INCREMENT OF DISP.         4.885E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.225E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.727E-04   TIME AVG. HEAT FLUX   8.744E-05
 LARGEST RESIDUAL HEAT FLUX         3.091E-15   AT NODE       2931   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.676E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.464E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.139       TIME AVG. FORCE       7.155E-02
 LARGEST SCALED RESIDUAL FORCE      1.617E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.617E-08
 LARGEST INCREMENT OF DISP.         5.009E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.231E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    8.902E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.769E-04   TIME AVG. HEAT FLUX   8.953E-05
 LARGEST RESIDUAL HEAT FLUX         1.186E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.690E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.462E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  355ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.139       TIME AVG. FORCE       7.155E-02
 LARGEST SCALED RESIDUAL FORCE      9.696E-13   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      9.696E-13
 LARGEST INCREMENT OF DISP.         5.009E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.379E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.770E-04   TIME AVG. HEAT FLUX   8.953E-05
 LARGEST RESIDUAL HEAT FLUX         2.830E-15   AT NODE       4536   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.690E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.646E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  368ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.143       TIME AVG. FORCE       7.312E-02
 LARGEST SCALED RESIDUAL FORCE      3.987E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.987E-10
 LARGEST INCREMENT OF DISP.         5.132E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.231E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.812E-04   TIME AVG. HEAT FLUX   9.161E-05
 LARGEST RESIDUAL HEAT FLUX        -2.831E-15   AT NODE       1424   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.705E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.469E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.146       TIME AVG. FORCE       7.470E-02
 LARGEST SCALED RESIDUAL FORCE      1.781E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.781E-08
 LARGEST INCREMENT OF DISP.         5.256E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.238E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.381E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.855E-04   TIME AVG. HEAT FLUX   9.370E-05
 LARGEST RESIDUAL HEAT FLUX         1.180E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.720E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.467E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  376ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.146       TIME AVG. FORCE       7.470E-02
 LARGEST SCALED RESIDUAL FORCE      1.036E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.036E-12
 LARGEST INCREMENT OF DISP.         5.256E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.439E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.855E-04   TIME AVG. HEAT FLUX   9.370E-05
 LARGEST RESIDUAL HEAT FLUX        -2.866E-15   AT NODE       3461   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.720E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.626E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  382ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.149       TIME AVG. FORCE       7.628E-02
 LARGEST SCALED RESIDUAL FORCE      4.233E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.233E-10
 LARGEST INCREMENT OF DISP.         5.379E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.237E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.897E-04   TIME AVG. HEAT FLUX   9.578E-05
 LARGEST RESIDUAL HEAT FLUX        -2.453E-15   AT NODE        179   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.735E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.472E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  24ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.152       TIME AVG. FORCE       7.785E-02
 LARGEST SCALED RESIDUAL FORCE      1.953E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.953E-08
 LARGEST INCREMENT OF DISP.         5.504E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.244E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.865E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.940E-04   TIME AVG. HEAT FLUX   9.787E-05
 LARGEST RESIDUAL HEAT FLUX         1.172E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.749E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.470E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  384ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.152       TIME AVG. FORCE       7.785E-02
 LARGEST SCALED RESIDUAL FORCE      1.102E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.102E-12
 LARGEST INCREMENT OF DISP.         5.504E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.497E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  1.940E-04   TIME AVG. HEAT FLUX   9.787E-05
 LARGEST RESIDUAL HEAT FLUX         2.583E-15   AT NODE       4318   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.749E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.602E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  381ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.155       TIME AVG. FORCE       7.942E-02
 LARGEST SCALED RESIDUAL FORCE      4.474E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.474E-10
 LARGEST INCREMENT OF DISP.         5.628E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.243E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  1.983E-04   TIME AVG. HEAT FLUX   9.997E-05
 LARGEST RESIDUAL HEAT FLUX         2.808E-15   AT NODE       4618   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.764E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.475E-04   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.158       TIME AVG. FORCE       8.099E-02
 LARGEST SCALED RESIDUAL FORCE      2.133E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.133E-08
 LARGEST INCREMENT OF DISP.         5.753E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.250E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.036E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.026E-04   TIME AVG. HEAT FLUX   1.021E-04
 LARGEST RESIDUAL HEAT FLUX         1.164E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.779E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.473E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  368ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.158       TIME AVG. FORCE       8.099E-02
 LARGEST SCALED RESIDUAL FORCE      1.167E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.167E-12
 LARGEST INCREMENT OF DISP.         5.753E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.553E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.026E-04   TIME AVG. HEAT FLUX   1.021E-04
 LARGEST RESIDUAL HEAT FLUX         2.726E-15   AT NODE       1587   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.779E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.573E-07   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.161       TIME AVG. FORCE       8.256E-02
 LARGEST SCALED RESIDUAL FORCE      4.707E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.707E-10
 LARGEST INCREMENT OF DISP.         5.878E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.249E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.069E-04   TIME AVG. HEAT FLUX   1.042E-04
 LARGEST RESIDUAL HEAT FLUX        -2.701E-15   AT NODE       2511   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.794E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.477E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  24ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.164       TIME AVG. FORCE       8.413E-02
 LARGEST SCALED RESIDUAL FORCE      2.323E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.323E-08
 LARGEST INCREMENT OF DISP.         6.004E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.257E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.085E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.112E-04   TIME AVG. HEAT FLUX   1.063E-04
 LARGEST RESIDUAL HEAT FLUX         1.154E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.808E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.474E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  399ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.164       TIME AVG. FORCE       8.413E-02
 LARGEST SCALED RESIDUAL FORCE      1.229E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.229E-12
 LARGEST INCREMENT OF DISP.         6.004E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.607E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.112E-04   TIME AVG. HEAT FLUX   1.063E-04
 LARGEST RESIDUAL HEAT FLUX        -2.822E-15   AT NODE       2715   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.808E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.538E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  378ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.167       TIME AVG. FORCE       8.570E-02
 LARGEST SCALED RESIDUAL FORCE      4.932E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.932E-10
 LARGEST INCREMENT OF DISP.         6.130E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.256E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.155E-04   TIME AVG. HEAT FLUX   1.084E-04
 LARGEST RESIDUAL HEAT FLUX        -2.932E-15   AT NODE       5039   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.823E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.477E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.170       TIME AVG. FORCE       8.727E-02
 LARGEST SCALED RESIDUAL FORCE      2.521E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.521E-08
 LARGEST INCREMENT OF DISP.         6.256E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.264E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.136E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.199E-04   TIME AVG. HEAT FLUX   1.105E-04
 LARGEST RESIDUAL HEAT FLUX         1.143E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.838E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.474E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  381ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.170       TIME AVG. FORCE       8.727E-02
 LARGEST SCALED RESIDUAL FORCE      1.290E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.290E-12
 LARGEST INCREMENT OF DISP.         6.256E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.657E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.199E-04   TIME AVG. HEAT FLUX   1.105E-04
 LARGEST RESIDUAL HEAT FLUX        -2.637E-15   AT NODE       3098   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.838E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.498E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  385ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.173       TIME AVG. FORCE       8.884E-02
 LARGEST SCALED RESIDUAL FORCE      5.146E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.146E-10
 LARGEST INCREMENT OF DISP.         6.382E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.262E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.241E-04   TIME AVG. HEAT FLUX   1.126E-04
 LARGEST RESIDUAL HEAT FLUX         3.442E-15   AT NODE       1743   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.853E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.476E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.176       TIME AVG. FORCE       9.040E-02
 LARGEST SCALED RESIDUAL FORCE      2.728E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.728E-08
 LARGEST INCREMENT OF DISP.         6.510E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.271E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.187E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.284E-04   TIME AVG. HEAT FLUX   1.147E-04
 LARGEST RESIDUAL HEAT FLUX         1.131E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.867E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.472E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  386ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.176       TIME AVG. FORCE       9.040E-02
 LARGEST SCALED RESIDUAL FORCE      1.348E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.348E-12
 LARGEST INCREMENT OF DISP.         6.510E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.705E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.284E-04   TIME AVG. HEAT FLUX   1.147E-04
 LARGEST RESIDUAL HEAT FLUX         3.156E-15   AT NODE       3409   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.867E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.452E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  379ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.180       TIME AVG. FORCE       9.196E-02
 LARGEST SCALED RESIDUAL FORCE      5.350E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.350E-10
 LARGEST INCREMENT OF DISP.         6.637E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.269E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.328E-04   TIME AVG. HEAT FLUX   1.168E-04
 LARGEST RESIDUAL HEAT FLUX         2.864E-15   AT NODE       4341   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.882E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.474E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.183       TIME AVG. FORCE       9.352E-02
 LARGEST SCALED RESIDUAL FORCE      2.945E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      2.945E-08
 LARGEST INCREMENT OF DISP.         6.764E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.278E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.239E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.371E-04   TIME AVG. HEAT FLUX   1.189E-04
 LARGEST RESIDUAL HEAT FLUX         1.118E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.897E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.470E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.183       TIME AVG. FORCE       9.352E-02
 LARGEST SCALED RESIDUAL FORCE      1.403E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.403E-12
 LARGEST INCREMENT OF DISP.         6.765E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.748E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.371E-04   TIME AVG. HEAT FLUX   1.189E-04
 LARGEST RESIDUAL HEAT FLUX        -2.588E-15   AT NODE       4937   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.897E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.401E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  372ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.186       TIME AVG. FORCE       9.509E-02
 LARGEST SCALED RESIDUAL FORCE      5.540E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.540E-10
 LARGEST INCREMENT OF DISP.         6.892E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.276E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.415E-04   TIME AVG. HEAT FLUX   1.210E-04
 LARGEST RESIDUAL HEAT FLUX        -2.682E-15   AT NODE       5433   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.912E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.470E-04   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.189       TIME AVG. FORCE       9.664E-02
 LARGEST SCALED RESIDUAL FORCE      3.171E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.171E-08
 LARGEST INCREMENT OF DISP.         7.021E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.285E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.291E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.458E-04   TIME AVG. HEAT FLUX   1.231E-04
 LARGEST RESIDUAL HEAT FLUX         1.104E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.926E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.465E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  377ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.189       TIME AVG. FORCE       9.664E-02
 LARGEST SCALED RESIDUAL FORCE      1.455E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.455E-12
 LARGEST INCREMENT OF DISP.         7.021E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.788E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.458E-04   TIME AVG. HEAT FLUX   1.231E-04
 LARGEST RESIDUAL HEAT FLUX         2.682E-15   AT NODE       3837   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.926E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.345E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.192       TIME AVG. FORCE       9.820E-02
 LARGEST SCALED RESIDUAL FORCE      5.717E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.717E-10
 LARGEST INCREMENT OF DISP.         7.149E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.283E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.501E-04   TIME AVG. HEAT FLUX   1.252E-04
 LARGEST RESIDUAL HEAT FLUX         2.428E-15   AT NODE       1755   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.941E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.465E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.195       TIME AVG. FORCE       0.100    
 LARGEST SCALED RESIDUAL FORCE      3.407E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.407E-08
 LARGEST INCREMENT OF DISP.         7.279E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.293E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.345E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.545E-04   TIME AVG. HEAT FLUX   1.274E-04
 LARGEST RESIDUAL HEAT FLUX         1.088E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.956E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.460E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.195       TIME AVG. FORCE       0.100    
 LARGEST SCALED RESIDUAL FORCE      1.503E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.503E-12
 LARGEST INCREMENT OF DISP.         7.279E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.822E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.545E-04   TIME AVG. HEAT FLUX   1.274E-04
 LARGEST RESIDUAL HEAT FLUX        -2.842E-15   AT NODE       2204   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.956E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.282E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  378ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.198       TIME AVG. FORCE       0.101    
 LARGEST SCALED RESIDUAL FORCE      5.878E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.878E-10
 LARGEST INCREMENT OF DISP.         7.408E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.291E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.589E-04   TIME AVG. HEAT FLUX   1.295E-04
 LARGEST RESIDUAL HEAT FLUX         3.521E-15   AT NODE       4937   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.970E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.459E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.201       TIME AVG. FORCE       0.103    
 LARGEST SCALED RESIDUAL FORCE      3.652E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.652E-08
 LARGEST INCREMENT OF DISP.         7.538E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.301E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.399E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.633E-04   TIME AVG. HEAT FLUX   1.316E-04
 LARGEST RESIDUAL HEAT FLUX         1.070E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.985E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.452E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.201       TIME AVG. FORCE       0.103    
 LARGEST SCALED RESIDUAL FORCE      1.546E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.546E-12
 LARGEST INCREMENT OF DISP.         7.538E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.851E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.633E-04   TIME AVG. HEAT FLUX   1.316E-04
 LARGEST RESIDUAL HEAT FLUX        -2.646E-15   AT NODE       4746   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.985E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.214E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.204       TIME AVG. FORCE       0.104    
 LARGEST SCALED RESIDUAL FORCE      6.022E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.022E-10
 LARGEST INCREMENT OF DISP.         7.668E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.299E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.677E-04   TIME AVG. HEAT FLUX   1.337E-04
 LARGEST RESIDUAL HEAT FLUX        -2.718E-15   AT NODE       4508   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.999E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.450E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.207       TIME AVG. FORCE       0.106    
 LARGEST SCALED RESIDUAL FORCE      3.907E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      3.907E-08
 LARGEST INCREMENT OF DISP.         7.799E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.309E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.454E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.721E-04   TIME AVG. HEAT FLUX   1.359E-04
 LARGEST RESIDUAL HEAT FLUX         1.052E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.014E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.443E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.207       TIME AVG. FORCE       0.106    
 LARGEST SCALED RESIDUAL FORCE      1.585E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.585E-12
 LARGEST INCREMENT OF DISP.         7.799E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.874E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.721E-04   TIME AVG. HEAT FLUX   1.359E-04
 LARGEST RESIDUAL HEAT FLUX         2.890E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.014E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.139E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  379ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.210       TIME AVG. FORCE       0.108    
 LARGEST SCALED RESIDUAL FORCE      6.148E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.148E-10
 LARGEST INCREMENT OF DISP.         7.930E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.307E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.765E-04   TIME AVG. HEAT FLUX   1.380E-04
 LARGEST RESIDUAL HEAT FLUX         2.849E-15   AT NODE       3939   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.028E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.440E-04   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.213       TIME AVG. FORCE       0.109    
 LARGEST SCALED RESIDUAL FORCE      4.173E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.173E-08
 LARGEST INCREMENT OF DISP.         8.061E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.317E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.510E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.809E-04   TIME AVG. HEAT FLUX   1.401E-04
 LARGEST RESIDUAL HEAT FLUX         1.032E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.043E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.433E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  380ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.213       TIME AVG. FORCE       0.109    
 LARGEST SCALED RESIDUAL FORCE      1.620E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.620E-12
 LARGEST INCREMENT OF DISP.         8.062E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.891E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.809E-04   TIME AVG. HEAT FLUX   1.401E-04
 LARGEST RESIDUAL HEAT FLUX         2.763E-15   AT NODE       3199   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.043E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.058E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  380ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.216       TIME AVG. FORCE       0.111    
 LARGEST SCALED RESIDUAL FORCE      6.254E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.254E-10
 LARGEST INCREMENT OF DISP.         8.193E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.315E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.854E-04   TIME AVG. HEAT FLUX   1.423E-04
 LARGEST RESIDUAL HEAT FLUX        -2.786E-15   AT NODE       4580   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.057E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.428E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  24ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.219       TIME AVG. FORCE       0.112    
 LARGEST SCALED RESIDUAL FORCE      4.448E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.448E-08
 LARGEST INCREMENT OF DISP.         8.326E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.326E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.568E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.898E-04   TIME AVG. HEAT FLUX   1.444E-04
 LARGEST RESIDUAL HEAT FLUX         1.011E-06   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.071E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.420E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  375ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.219       TIME AVG. FORCE       0.112    
 LARGEST SCALED RESIDUAL FORCE      1.648E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.648E-12
 LARGEST INCREMENT OF DISP.         8.326E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.901E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.898E-04   TIME AVG. HEAT FLUX   1.444E-04
 LARGEST RESIDUAL HEAT FLUX        -2.567E-15   AT NODE       5212   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.071E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.971E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  381ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.222       TIME AVG. FORCE       0.114    
 LARGEST SCALED RESIDUAL FORCE      6.340E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.340E-10
 LARGEST INCREMENT OF DISP.         8.458E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.323E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  2.942E-04   TIME AVG. HEAT FLUX   1.465E-04
 LARGEST RESIDUAL HEAT FLUX        -3.034E-15   AT NODE       2976   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.085E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.415E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  24ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.225       TIME AVG. FORCE       0.115    
 LARGEST SCALED RESIDUAL FORCE      4.735E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.735E-08
 LARGEST INCREMENT OF DISP.         8.592E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.335E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.626E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  2.986E-04   TIME AVG. HEAT FLUX   1.487E-04
 LARGEST RESIDUAL HEAT FLUX         9.882E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.099E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.406E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.225       TIME AVG. FORCE       0.115    
 LARGEST SCALED RESIDUAL FORCE      1.671E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.671E-12
 LARGEST INCREMENT OF DISP.         8.592E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.904E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  2.986E-04   TIME AVG. HEAT FLUX   1.487E-04
 LARGEST RESIDUAL HEAT FLUX        -2.781E-15   AT NODE       4023   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.099E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.878E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.228       TIME AVG. FORCE       0.117    
 LARGEST SCALED RESIDUAL FORCE      6.404E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.404E-10
 LARGEST INCREMENT OF DISP.         8.725E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.332E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.031E-04   TIME AVG. HEAT FLUX   1.508E-04
 LARGEST RESIDUAL HEAT FLUX        -2.801E-15   AT NODE       5328   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.113E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.399E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.230       TIME AVG. FORCE       0.118    
 LARGEST SCALED RESIDUAL FORCE      5.032E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.032E-08
 LARGEST INCREMENT OF DISP.         8.859E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.344E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.685E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.075E-04   TIME AVG. HEAT FLUX   1.530E-04
 LARGEST RESIDUAL HEAT FLUX         9.642E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.127E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.390E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.230       TIME AVG. FORCE       0.118    
 LARGEST SCALED RESIDUAL FORCE      1.688E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.688E-12
 LARGEST INCREMENT OF DISP.         8.860E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.899E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.075E-04   TIME AVG. HEAT FLUX   1.530E-04
 LARGEST RESIDUAL HEAT FLUX        -2.995E-15   AT NODE       2368   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.127E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.778E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.233       TIME AVG. FORCE       0.120    
 LARGEST SCALED RESIDUAL FORCE      6.444E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.444E-10
 LARGEST INCREMENT OF DISP.         8.994E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.341E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.119E-04   TIME AVG. HEAT FLUX   1.551E-04
 LARGEST RESIDUAL HEAT FLUX        -2.768E-15   AT NODE        260   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.141E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.382E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.236       TIME AVG. FORCE       0.121    
 LARGEST SCALED RESIDUAL FORCE      5.340E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.340E-08
 LARGEST INCREMENT OF DISP.         9.129E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.353E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.745E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.165E-04   TIME AVG. HEAT FLUX   1.573E-04
 LARGEST RESIDUAL HEAT FLUX         9.388E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.155E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.371E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  376ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.236       TIME AVG. FORCE       0.121    
 LARGEST SCALED RESIDUAL FORCE      1.698E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.698E-12
 LARGEST INCREMENT OF DISP.         9.129E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.885E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.165E-04   TIME AVG. HEAT FLUX   1.573E-04
 LARGEST RESIDUAL HEAT FLUX         2.536E-15   AT NODE       1930   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.155E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.672E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  390ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.239       TIME AVG. FORCE       0.123    
 LARGEST SCALED RESIDUAL FORCE      6.459E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.459E-10
 LARGEST INCREMENT OF DISP.         9.264E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.351E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.210E-04   TIME AVG. HEAT FLUX   1.594E-04
 LARGEST RESIDUAL HEAT FLUX         3.131E-15   AT NODE       3184   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.168E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.362E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.242       TIME AVG. FORCE       0.124    
 LARGEST SCALED RESIDUAL FORCE      5.660E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.660E-08
 LARGEST INCREMENT OF DISP.         9.401E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.363E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.807E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.255E-04   TIME AVG. HEAT FLUX   1.616E-04
 LARGEST RESIDUAL HEAT FLUX         9.120E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.182E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.351E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  371ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.242       TIME AVG. FORCE       0.124    
 LARGEST SCALED RESIDUAL FORCE      1.701E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.701E-12
 LARGEST INCREMENT OF DISP.         9.401E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.863E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.255E-04   TIME AVG. HEAT FLUX   1.616E-04
 LARGEST RESIDUAL HEAT FLUX        -2.703E-15   AT NODE        840   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.182E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.559E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  381ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.245       TIME AVG. FORCE       0.126    
 LARGEST SCALED RESIDUAL FORCE      6.556E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.556E-10
 LARGEST INCREMENT OF DISP.         9.537E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.360E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.299E-04   TIME AVG. HEAT FLUX   1.637E-04
 LARGEST RESIDUAL HEAT FLUX         2.533E-15   AT NODE       1650   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.195E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.341E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.248       TIME AVG. FORCE       0.127    
 LARGEST SCALED RESIDUAL FORCE      5.991E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.991E-08
 LARGEST INCREMENT OF DISP.         9.674E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.373E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.869E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.343E-04   TIME AVG. HEAT FLUX   1.659E-04
 LARGEST RESIDUAL HEAT FLUX         8.838E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.209E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.329E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  375ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.248       TIME AVG. FORCE       0.127    
 LARGEST SCALED RESIDUAL FORCE      1.697E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.697E-12
 LARGEST INCREMENT OF DISP.         9.674E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.832E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.343E-04   TIME AVG. HEAT FLUX   1.659E-04
 LARGEST RESIDUAL HEAT FLUX        -3.320E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.209E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.440E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  390ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.251       TIME AVG. FORCE       0.129    
 LARGEST SCALED RESIDUAL FORCE      6.661E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.661E-10
 LARGEST INCREMENT OF DISP.         9.811E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.370E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.388E-04   TIME AVG. HEAT FLUX   1.681E-04
 LARGEST RESIDUAL HEAT FLUX        -2.914E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.222E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.317E-04   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.254       TIME AVG. FORCE       0.131    
 LARGEST SCALED RESIDUAL FORCE      6.334E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.334E-08
 LARGEST INCREMENT OF DISP.         9.950E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.383E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.933E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.432E-04   TIME AVG. HEAT FLUX   1.702E-04
 LARGEST RESIDUAL HEAT FLUX         8.543E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.235E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.304E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          HEAT FLUX EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  381ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.254       TIME AVG. FORCE       0.131    
 LARGEST SCALED RESIDUAL FORCE      1.685E-12   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.685E-12
 LARGEST INCREMENT OF DISP.         9.950E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.790E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  3.432E-04   TIME AVG. HEAT FLUX   1.702E-04
 LARGEST RESIDUAL HEAT FLUX        -2.701E-15   AT NODE       4426   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.235E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.314E-07   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  397ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.257       TIME AVG. FORCE       0.132    
 LARGEST SCALED RESIDUAL FORCE      6.749E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.749E-10
 LARGEST INCREMENT OF DISP.         1.009E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.381E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.476E-04   TIME AVG. HEAT FLUX   1.724E-04
 LARGEST RESIDUAL HEAT FLUX         2.515E-15   AT NODE       2452   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.248E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.291E-04   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.259       TIME AVG. FORCE       0.134    
 LARGEST SCALED RESIDUAL FORCE      6.690E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.690E-08
 LARGEST INCREMENT OF DISP.         1.023E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.394E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.999E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.522E-04   TIME AVG. HEAT FLUX   1.745E-04
 LARGEST RESIDUAL HEAT FLUX         8.233E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.261E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.277E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  8.400E-02
 STEP TIME COMPLETED       8.400E+05,  TOTAL TIME COMPLETED        8.400E+05


  INCREMENT    85 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  379ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.262       TIME AVG. FORCE       0.135    
 LARGEST SCALED RESIDUAL FORCE      6.848E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.848E-10
 LARGEST INCREMENT OF DISP.         1.037E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.395E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.567E-04   TIME AVG. HEAT FLUX   1.767E-04
 LARGEST RESIDUAL HEAT FLUX         2.650E-15   AT NODE        223   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.273E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.268E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.265       TIME AVG. FORCE       0.137    
 LARGEST SCALED RESIDUAL FORCE      7.046E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      7.046E-08
 LARGEST INCREMENT OF DISP.         1.051E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.403E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.064E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.613E-04   TIME AVG. HEAT FLUX   1.789E-04
 LARGEST RESIDUAL HEAT FLUX         7.903E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.286E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.246E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  8.600E-02
 STEP TIME COMPLETED       8.600E+05,  TOTAL TIME COMPLETED        8.600E+05


  INCREMENT    87 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  376ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.268       TIME AVG. FORCE       0.138    
 LARGEST SCALED RESIDUAL FORCE      6.893E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.893E-10
 LARGEST INCREMENT OF DISP.         1.065E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.406E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.658E-04   TIME AVG. HEAT FLUX   1.811E-04
 LARGEST RESIDUAL HEAT FLUX        -3.069E-15   AT NODE       3924   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.298E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.236E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.271       TIME AVG. FORCE       0.140    
 LARGEST SCALED RESIDUAL FORCE      7.425E-08   AT NODE      11094   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      7.425E-08
 LARGEST INCREMENT OF DISP.         1.079E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.415E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.132E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.704E-04   TIME AVG. HEAT FLUX   1.832E-04
 LARGEST RESIDUAL HEAT FLUX         7.568E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.310E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.215E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  8.800E-02
 STEP TIME COMPLETED       8.800E+05,  TOTAL TIME COMPLETED        8.800E+05


  INCREMENT    89 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.274       TIME AVG. FORCE       0.141    
 LARGEST SCALED RESIDUAL FORCE      6.918E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.918E-10
 LARGEST INCREMENT OF DISP.         1.093E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.417E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.750E-04   TIME AVG. HEAT FLUX   1.854E-04
 LARGEST RESIDUAL HEAT FLUX        -3.117E-15   AT NODE       4309   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.322E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.203E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.276       TIME AVG. FORCE       0.143    
 LARGEST SCALED RESIDUAL FORCE     -7.890E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.890E-08
 LARGEST INCREMENT OF DISP.         1.107E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.427E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.222E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.795E-04   TIME AVG. HEAT FLUX   1.876E-04
 LARGEST RESIDUAL HEAT FLUX         7.221E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.334E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.181E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.000E-02
 STEP TIME COMPLETED       9.000E+05,  TOTAL TIME COMPLETED        9.000E+05


  INCREMENT    91 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.279       TIME AVG. FORCE       0.144    
 LARGEST SCALED RESIDUAL FORCE      6.920E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.920E-10
 LARGEST INCREMENT OF DISP.         1.122E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.429E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.841E-04   TIME AVG. HEAT FLUX   1.898E-04
 LARGEST RESIDUAL HEAT FLUX        -3.442E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.346E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.166E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.282       TIME AVG. FORCE       0.146    
 LARGEST SCALED RESIDUAL FORCE     -8.421E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.421E-08
 LARGEST INCREMENT OF DISP.         1.136E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.439E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.328E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.887E-04   TIME AVG. HEAT FLUX   1.920E-04
 LARGEST RESIDUAL HEAT FLUX         6.861E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.357E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.144E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.200E-02
 STEP TIME COMPLETED       9.200E+05,  TOTAL TIME COMPLETED        9.200E+05


  INCREMENT    93 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.285       TIME AVG. FORCE       0.147    
 LARGEST SCALED RESIDUAL FORCE      6.898E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.898E-10
 LARGEST INCREMENT OF DISP.         1.150E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.441E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  3.932E-04   TIME AVG. HEAT FLUX   1.942E-04
 LARGEST RESIDUAL HEAT FLUX         2.786E-15   AT NODE       1987   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.368E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.127E-04   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.288       TIME AVG. FORCE       0.149    
 LARGEST SCALED RESIDUAL FORCE     -8.983E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.983E-08
 LARGEST INCREMENT OF DISP.         1.165E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.452E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.439E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  3.978E-04   TIME AVG. HEAT FLUX   1.963E-04
 LARGEST RESIDUAL HEAT FLUX         6.488E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.380E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.105E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.400E-02
 STEP TIME COMPLETED       9.400E+05,  TOTAL TIME COMPLETED        9.400E+05


  INCREMENT    95 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  372ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.290       TIME AVG. FORCE       0.150    
 LARGEST SCALED RESIDUAL FORCE      6.849E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.849E-10
 LARGEST INCREMENT OF DISP.         1.179E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.454E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.024E-04   TIME AVG. HEAT FLUX   1.985E-04
 LARGEST RESIDUAL HEAT FLUX         2.750E-15   AT NODE       1939   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.390E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.086E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.293       TIME AVG. FORCE       0.152    
 LARGEST SCALED RESIDUAL FORCE     -9.577E-08   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.577E-08
 LARGEST INCREMENT OF DISP.         1.194E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.465E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.555E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.070E-04   TIME AVG. HEAT FLUX   2.007E-04
 LARGEST RESIDUAL HEAT FLUX         6.102E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.401E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.064E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.000E+04,  FRACTION OF STEP COMPLETED  9.600E-02
 STEP TIME COMPLETED       9.600E+05,  TOTAL TIME COMPLETED        9.600E+05


  INCREMENT    97 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.000E+04
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.296       TIME AVG. FORCE       0.153    
 LARGEST SCALED RESIDUAL FORCE      6.774E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.774E-10
 LARGEST INCREMENT OF DISP.         1.209E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.467E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.117E-04   TIME AVG. HEAT FLUX   2.029E-04
 LARGEST RESIDUAL HEAT FLUX        -2.876E-15   AT NODE       2828   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.411E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.041E-04   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.299       TIME AVG. FORCE       0.155    
 LARGEST SCALED RESIDUAL FORCE     -1.020E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.020E-07
 LARGEST INCREMENT OF DISP.         1.224E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.478E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.677E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.162E-04   TIME AVG. HEAT FLUX   2.051E-04
 LARGEST RESIDUAL HEAT FLUX         5.705E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.422E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.019E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.301       TIME AVG. FORCE       0.156    
 LARGEST SCALED RESIDUAL FORCE      6.670E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.670E-10
 LARGEST INCREMENT OF DISP.         1.238E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.480E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.208E-04   TIME AVG. HEAT FLUX   2.073E-04
 LARGEST RESIDUAL HEAT FLUX         3.014E-15   AT NODE        932   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.432E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        9.936E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.304       TIME AVG. FORCE       0.158    
 LARGEST SCALED RESIDUAL FORCE     -1.087E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.087E-07
 LARGEST INCREMENT OF DISP.         1.253E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.493E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.805E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.253E-04   TIME AVG. HEAT FLUX   2.095E-04
 LARGEST RESIDUAL HEAT FLUX         5.295E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.441E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        9.719E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  384ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.307       TIME AVG. FORCE       0.159    
 LARGEST SCALED RESIDUAL FORCE      6.536E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.536E-10
 LARGEST INCREMENT OF DISP.         1.268E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.494E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.300E-04   TIME AVG. HEAT FLUX   2.117E-04
 LARGEST RESIDUAL HEAT FLUX         3.202E-15   AT NODE       1417   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.451E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        9.431E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  20ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.310       TIME AVG. FORCE       0.160    
 LARGEST SCALED RESIDUAL FORCE     -1.157E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.157E-07
 LARGEST INCREMENT OF DISP.         1.283E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.507E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.940E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.346E-04   TIME AVG. HEAT FLUX   2.139E-04
 LARGEST RESIDUAL HEAT FLUX         4.875E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.460E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        9.216E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.312       TIME AVG. FORCE       0.162    
 LARGEST SCALED RESIDUAL FORCE      6.370E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.370E-10
 LARGEST INCREMENT OF DISP.         1.298E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.509E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.392E-04   TIME AVG. HEAT FLUX   2.161E-04
 LARGEST RESIDUAL HEAT FLUX        -2.831E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.469E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.896E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  20ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.315       TIME AVG. FORCE       0.163    
 LARGEST SCALED RESIDUAL FORCE     -1.232E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.232E-07
 LARGEST INCREMENT OF DISP.         1.314E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.523E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -3.082E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.438E-04   TIME AVG. HEAT FLUX   2.184E-04
 LARGEST RESIDUAL HEAT FLUX         4.443E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.478E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.684E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.318       TIME AVG. FORCE       0.165    
 LARGEST SCALED RESIDUAL FORCE      6.170E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      6.170E-10
 LARGEST INCREMENT OF DISP.         1.329E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.524E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.484E-04   TIME AVG. HEAT FLUX   2.206E-04
 LARGEST RESIDUAL HEAT FLUX         2.898E-15   AT NODE       4615   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.486E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.330E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.320       TIME AVG. FORCE       0.166    
 LARGEST SCALED RESIDUAL FORCE     -1.311E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.311E-07
 LARGEST INCREMENT OF DISP.         1.344E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.539E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -3.232E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.530E-04   TIME AVG. HEAT FLUX   2.228E-04
 LARGEST RESIDUAL HEAT FLUX         4.002E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.494E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.121E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.323       TIME AVG. FORCE       0.168    
 LARGEST SCALED RESIDUAL FORCE      5.936E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.936E-10
 LARGEST INCREMENT OF DISP.         1.360E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.540E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.577E-04   TIME AVG. HEAT FLUX   2.250E-04
 LARGEST RESIDUAL HEAT FLUX         2.873E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.502E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.732E-05   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.326       TIME AVG. FORCE       0.169    
 LARGEST SCALED RESIDUAL FORCE     -1.394E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.394E-07
 LARGEST INCREMENT OF DISP.         1.375E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.555E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -3.390E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.623E-04   TIME AVG. HEAT FLUX   2.272E-04
 LARGEST RESIDUAL HEAT FLUX         3.551E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.509E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.527E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  394ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.328       TIME AVG. FORCE       0.171    
 LARGEST SCALED RESIDUAL FORCE      5.666E-10   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      5.666E-10
 LARGEST INCREMENT OF DISP.         1.391E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.556E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.670E-04   TIME AVG. HEAT FLUX   2.294E-04
 LARGEST RESIDUAL HEAT FLUX        -2.692E-15   AT NODE       2666   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.516E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.101E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.331       TIME AVG. FORCE       0.172    
 LARGEST SCALED RESIDUAL FORCE     -1.483E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.483E-07
 LARGEST INCREMENT OF DISP.         1.406E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.572E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -3.557E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.716E-04   TIME AVG. HEAT FLUX   2.317E-04
 LARGEST RESIDUAL HEAT FLUX         3.090E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.523E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.901E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  371ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.333       TIME AVG. FORCE       0.174    
 LARGEST SCALED RESIDUAL FORCE     -5.423E-10   AT NODE       6990   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.423E-10
 LARGEST INCREMENT OF DISP.         1.422E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.573E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.763E-04   TIME AVG. HEAT FLUX   2.339E-04
 LARGEST RESIDUAL HEAT FLUX         2.778E-15   AT NODE      15745   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.530E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.436E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  21ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.336       TIME AVG. FORCE       0.175    
 LARGEST SCALED RESIDUAL FORCE     -1.577E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.577E-07
 LARGEST INCREMENT OF DISP.         1.438E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.590E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -3.733E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.809E-04   TIME AVG. HEAT FLUX   2.361E-04
 LARGEST RESIDUAL HEAT FLUX         2.622E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.536E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.242E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  376ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.339       TIME AVG. FORCE       0.177    
 LARGEST SCALED RESIDUAL FORCE     -5.285E-10   AT NODE       6990   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.285E-10
 LARGEST INCREMENT OF DISP.         1.454E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.591E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.856E-04   TIME AVG. HEAT FLUX   2.383E-04
 LARGEST RESIDUAL HEAT FLUX        -2.731E-15   AT NODE       5047   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.542E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.737E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.341       TIME AVG. FORCE       0.178    
 LARGEST SCALED RESIDUAL FORCE     -1.676E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.676E-07
 LARGEST INCREMENT OF DISP.         1.470E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.609E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -3.920E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.903E-04   TIME AVG. HEAT FLUX   2.406E-04
 LARGEST RESIDUAL HEAT FLUX         2.146E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.547E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.549E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  373ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.344       TIME AVG. FORCE       0.179    
 LARGEST SCALED RESIDUAL FORCE     -5.121E-10   AT NODE       6990   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.121E-10
 LARGEST INCREMENT OF DISP.         1.486E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.610E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  4.948E-04   TIME AVG. HEAT FLUX   2.428E-04
 LARGEST RESIDUAL HEAT FLUX         4.781E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.552E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.003E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.346       TIME AVG. FORCE       0.181    
 LARGEST SCALED RESIDUAL FORCE     -1.782E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.782E-07
 LARGEST INCREMENT OF DISP.         1.503E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.628E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -4.117E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  4.993E-04   TIME AVG. HEAT FLUX   2.450E-04
 LARGEST RESIDUAL HEAT FLUX         1.663E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.557E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.823E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  387ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.349       TIME AVG. FORCE       0.182    
 LARGEST SCALED RESIDUAL FORCE     -4.928E-10   AT NODE       6990   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.928E-10
 LARGEST INCREMENT OF DISP.         1.519E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.629E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  5.040E-04   TIME AVG. HEAT FLUX   2.473E-04
 LARGEST RESIDUAL HEAT FLUX         2.649E-15   AT NODE       4658   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.561E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.232E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.351       TIME AVG. FORCE       0.184    
 LARGEST SCALED RESIDUAL FORCE     -1.895E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.895E-07
 LARGEST INCREMENT OF DISP.         1.535E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.649E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -4.327E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.086E-04   TIME AVG. HEAT FLUX   2.495E-04
 LARGEST RESIDUAL HEAT FLUX         1.174E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.565E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        4.061E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  394ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.354       TIME AVG. FORCE       0.185    
 LARGEST SCALED RESIDUAL FORCE     -5.194E-10   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.194E-10
 LARGEST INCREMENT OF DISP.         1.552E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.649E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  5.132E-04   TIME AVG. HEAT FLUX   2.517E-04
 LARGEST RESIDUAL HEAT FLUX         3.483E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.569E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.425E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.356       TIME AVG. FORCE       0.187    
 LARGEST SCALED RESIDUAL FORCE     -2.014E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.014E-07
 LARGEST INCREMENT OF DISP.         1.568E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.670E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -4.550E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.179E-04   TIME AVG. HEAT FLUX   2.540E-04
 LARGEST RESIDUAL HEAT FLUX         6.808E-08   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.572E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        3.264E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  378ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.359       TIME AVG. FORCE       0.188    
 LARGEST SCALED RESIDUAL FORCE     -5.588E-10   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.588E-10
 LARGEST INCREMENT OF DISP.         1.585E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.670E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  5.225E-04   TIME AVG. HEAT FLUX   2.562E-04
 LARGEST RESIDUAL HEAT FLUX        -2.693E-15   AT NODE       4397   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.575E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.581E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.361       TIME AVG. FORCE       0.189    
 LARGEST SCALED RESIDUAL FORCE     -2.141E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.141E-07
 LARGEST INCREMENT OF DISP.         1.602E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.692E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -4.787E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.272E-04   TIME AVG. HEAT FLUX   2.584E-04
 LARGEST RESIDUAL HEAT FLUX        -5.558E-08   AT NODE       6706   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.577E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.430E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  376ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.364       TIME AVG. FORCE       0.191    
 LARGEST SCALED RESIDUAL FORCE     -5.999E-10   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.999E-10
 LARGEST INCREMENT OF DISP.         1.619E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.692E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  5.317E-04   TIME AVG. HEAT FLUX   2.607E-04
 LARGEST RESIDUAL HEAT FLUX         2.845E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.579E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.698E-05   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.366       TIME AVG. FORCE       0.192    
 LARGEST SCALED RESIDUAL FORCE     -2.276E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.276E-07
 LARGEST INCREMENT OF DISP.         1.636E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.715E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -5.040E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.362E-04   TIME AVG. HEAT FLUX   2.629E-04
 LARGEST RESIDUAL HEAT FLUX        -5.630E-08   AT NODE       6706   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.580E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.560E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.369       TIME AVG. FORCE       0.194    
 LARGEST SCALED RESIDUAL FORCE     -6.429E-10   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.429E-10
 LARGEST INCREMENT OF DISP.         1.653E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.715E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  5.408E-04   TIME AVG. HEAT FLUX   2.652E-04
 LARGEST RESIDUAL HEAT FLUX         2.856E-15   AT NODE       3615   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.581E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.289E-05   AT NODE       6990   DOF 11
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
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.371       TIME AVG. FORCE       0.195    
 LARGEST SCALED RESIDUAL FORCE     -2.420E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.420E-07
 LARGEST INCREMENT OF DISP.         1.671E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.739E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -5.309E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.454E-04   TIME AVG. HEAT FLUX   2.674E-04
 LARGEST RESIDUAL HEAT FLUX        -8.173E-08   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.582E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.312E-05   AT NODE       6990   DOF 11
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
        SOLVER ELAPSED TIME:  381ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.373       TIME AVG. FORCE       0.197    
 LARGEST SCALED RESIDUAL FORCE     -6.877E-10   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.877E-10
 LARGEST INCREMENT OF DISP.         1.688E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.739E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  5.499E-04   TIME AVG. HEAT FLUX   2.696E-04
 LARGEST RESIDUAL HEAT FLUX         2.925E-15   AT NODE       2435   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.582E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.507E-05   AT NODE       6990   DOF 11
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


 AVERAGE FORCE                      0.376       TIME AVG. FORCE       0.198    
 LARGEST SCALED RESIDUAL FORCE     -2.573E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.573E-07
 LARGEST INCREMENT OF DISP.         1.706E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.764E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -5.596E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.545E-04   TIME AVG. HEAT FLUX   2.719E-04
 LARGEST RESIDUAL HEAT FLUX        -1.319E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.581E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.527E-05   AT NODE       6990   DOF 11
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
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.378       TIME AVG. FORCE       0.199    
 LARGEST SCALED RESIDUAL FORCE     -7.345E-10   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.345E-10
 LARGEST INCREMENT OF DISP.         1.723E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.764E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  5.589E-04   TIME AVG. HEAT FLUX   2.741E-04
 LARGEST RESIDUAL HEAT FLUX        -3.775E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.580E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.733E-05   AT NODE       6990   DOF 11
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


 AVERAGE FORCE                      0.381       TIME AVG. FORCE       0.201    
 LARGEST SCALED RESIDUAL FORCE     -2.737E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.737E-07
 LARGEST INCREMENT OF DISP.         1.741E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.791E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -5.903E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.634E-04   TIME AVG. HEAT FLUX   2.764E-04
 LARGEST RESIDUAL HEAT FLUX        -1.819E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.579E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.748E-05   AT NODE       6990   DOF 11
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
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.383       TIME AVG. FORCE       0.202    
 LARGEST SCALED RESIDUAL FORCE     -7.832E-10   AT NODE       6706   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.832E-10
 LARGEST INCREMENT OF DISP.         1.759E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.790E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  5.680E-04   TIME AVG. HEAT FLUX   2.786E-04
 LARGEST RESIDUAL HEAT FLUX         3.502E-15   AT NODE       2767   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.577E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.232E-05   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.385       TIME AVG. FORCE       0.203    
 LARGEST SCALED RESIDUAL FORCE     -2.911E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.911E-07
 LARGEST INCREMENT OF DISP.         1.777E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.819E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -6.233E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.725E-04   TIME AVG. HEAT FLUX   2.809E-04
 LARGEST RESIDUAL HEAT FLUX        -2.316E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.574E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.303E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.388       TIME AVG. FORCE       0.205    
 LARGEST SCALED RESIDUAL FORCE     -9.034E-10   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.034E-10
 LARGEST INCREMENT OF DISP.         1.796E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.817E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  5.771E-04   TIME AVG. HEAT FLUX   2.831E-04
 LARGEST RESIDUAL HEAT FLUX        -2.668E-15   AT NODE       5357   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.571E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.318E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.390       TIME AVG. FORCE       0.206    
 LARGEST SCALED RESIDUAL FORCE     -3.097E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.097E-07
 LARGEST INCREMENT OF DISP.         1.814E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.848E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -6.585E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.818E-04   TIME AVG. HEAT FLUX   2.853E-04
 LARGEST RESIDUAL HEAT FLUX        -2.810E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.568E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.368E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.392       TIME AVG. FORCE       0.208    
 LARGEST SCALED RESIDUAL FORCE     -1.040E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.040E-09
 LARGEST INCREMENT OF DISP.         1.832E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.846E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  5.862E-04   TIME AVG. HEAT FLUX   2.876E-04
 LARGEST RESIDUAL HEAT FLUX         3.176E-15   AT NODE      12450   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.563E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.447E-05   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.395       TIME AVG. FORCE       0.209    
 LARGEST SCALED RESIDUAL FORCE     -3.297E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.297E-07
 LARGEST INCREMENT OF DISP.         1.851E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.878E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -6.964E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.904E-04   TIME AVG. HEAT FLUX   2.898E-04
 LARGEST RESIDUAL HEAT FLUX        -3.297E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.559E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.473E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.397       TIME AVG. FORCE       0.210    
 LARGEST SCALED RESIDUAL FORCE     -1.185E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.185E-09
 LARGEST INCREMENT OF DISP.         1.870E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.876E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  5.950E-04   TIME AVG. HEAT FLUX   2.921E-04
 LARGEST RESIDUAL HEAT FLUX        -2.611E-15   AT NODE       3325   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.553E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.620E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  23ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.399       TIME AVG. FORCE       0.212    
 LARGEST SCALED RESIDUAL FORCE     -3.510E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.510E-07
 LARGEST INCREMENT OF DISP.         1.889E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.910E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -7.372E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  5.997E-04   TIME AVG. HEAT FLUX   2.943E-04
 LARGEST RESIDUAL HEAT FLUX        -3.777E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.547E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.619E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  373ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.402       TIME AVG. FORCE       0.213    
 LARGEST SCALED RESIDUAL FORCE     -1.339E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.339E-09
 LARGEST INCREMENT OF DISP.         1.908E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.908E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  6.044E-04   TIME AVG. HEAT FLUX   2.966E-04
 LARGEST RESIDUAL HEAT FLUX         3.220E-15   AT NODE       5006   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.541E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.837E-05   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.404       TIME AVG. FORCE       0.214    
 LARGEST SCALED RESIDUAL FORCE     -3.738E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.738E-07
 LARGEST INCREMENT OF DISP.         1.928E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.943E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -7.811E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.087E-04   TIME AVG. HEAT FLUX   2.988E-04
 LARGEST RESIDUAL HEAT FLUX        -4.247E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.534E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.808E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.406       TIME AVG. FORCE       0.216    
 LARGEST SCALED RESIDUAL FORCE     -1.501E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.501E-09
 LARGEST INCREMENT OF DISP.         1.947E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.941E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  6.130E-04   TIME AVG. HEAT FLUX   3.011E-04
 LARGEST RESIDUAL HEAT FLUX        -3.065E-15   AT NODE       1151   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.526E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.100E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.408       TIME AVG. FORCE       0.217    
 LARGEST SCALED RESIDUAL FORCE     -3.982E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.982E-07
 LARGEST INCREMENT OF DISP.         1.967E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.979E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -8.284E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.176E-04   TIME AVG. HEAT FLUX   3.033E-04
 LARGEST RESIDUAL HEAT FLUX        -4.706E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.518E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.038E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  356ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.411       TIME AVG. FORCE       0.219    
 LARGEST SCALED RESIDUAL FORCE     -1.672E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.672E-09
 LARGEST INCREMENT OF DISP.         1.987E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.975E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  6.219E-04   TIME AVG. HEAT FLUX   3.056E-04
 LARGEST RESIDUAL HEAT FLUX         2.712E-15   AT NODE       1735   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.508E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.409E-05   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.413       TIME AVG. FORCE       0.220    
 LARGEST SCALED RESIDUAL FORCE     -4.245E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.245E-07
 LARGEST INCREMENT OF DISP.         2.007E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.016E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -8.796E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.263E-04   TIME AVG. HEAT FLUX   3.078E-04
 LARGEST RESIDUAL HEAT FLUX        -5.151E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.499E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.312E-05   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.415       TIME AVG. FORCE       0.221    
 LARGEST SCALED RESIDUAL FORCE     -1.853E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.853E-09
 LARGEST INCREMENT OF DISP.         2.027E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.012E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  6.309E-04   TIME AVG. HEAT FLUX   3.100E-04
 LARGEST RESIDUAL HEAT FLUX         2.956E-15   AT NODE       3136   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.488E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.077E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  22ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.417       TIME AVG. FORCE       0.223    
 LARGEST SCALED RESIDUAL FORCE     -4.527E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.527E-07
 LARGEST INCREMENT OF DISP.         2.047E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.055E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -9.351E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.352E-04   TIME AVG. HEAT FLUX   3.123E-04
 LARGEST RESIDUAL HEAT FLUX        -5.581E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.477E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.063E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.419       TIME AVG. FORCE       0.224    
 LARGEST SCALED RESIDUAL FORCE     -2.044E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.044E-09
 LARGEST INCREMENT OF DISP.         2.068E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.050E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  6.395E-04   TIME AVG. HEAT FLUX   3.145E-04
 LARGEST RESIDUAL HEAT FLUX        -2.873E-15   AT NODE       3278   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.465E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.217E-04   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.422       TIME AVG. FORCE       0.225    
 LARGEST SCALED RESIDUAL FORCE     -4.830E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.830E-07
 LARGEST INCREMENT OF DISP.         2.089E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.096E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -9.952E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.438E-04   TIME AVG. HEAT FLUX   3.168E-04
 LARGEST RESIDUAL HEAT FLUX        -5.992E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.453E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.199E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  349ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.424       TIME AVG. FORCE       0.227    
 LARGEST SCALED RESIDUAL FORCE     -2.246E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.246E-09
 LARGEST INCREMENT OF DISP.         2.110E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.091E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 AVERAGE HEAT FLUX                  6.483E-04   TIME AVG. HEAT FLUX   3.190E-04
 LARGEST RESIDUAL HEAT FLUX         2.977E-15   AT NODE       3227   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.440E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.362E-04   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.426       TIME AVG. FORCE       0.228    
 LARGEST SCALED RESIDUAL FORCE     -5.157E-07   AT NODE      11986   DOF  2
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.157E-07
 LARGEST INCREMENT OF DISP.         2.131E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.139E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.061E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.525E-04   TIME AVG. HEAT FLUX   3.212E-04
 LARGEST RESIDUAL HEAT FLUX        -6.383E-07   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.426E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.340E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.428       TIME AVG. FORCE       0.229    
 LARGEST SCALED RESIDUAL FORCE     -2.460E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.460E-09
 LARGEST INCREMENT OF DISP.         2.153E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.134E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  6.568E-04   TIME AVG. HEAT FLUX   3.235E-04
 LARGEST RESIDUAL HEAT FLUX         3.136E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.411E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.513E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.430       TIME AVG. FORCE       0.231    
 LARGEST SCALED RESIDUAL FORCE     -2.393E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.393E-09
 LARGEST INCREMENT OF DISP.         2.174E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.182E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -4.909E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.609E-04   TIME AVG. HEAT FLUX   3.257E-04
 LARGEST RESIDUAL HEAT FLUX        -2.719E-15   AT NODE       2613   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.396E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.488E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

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
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.432       TIME AVG. FORCE       0.232    
 LARGEST SCALED RESIDUAL FORCE     -2.650E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.650E-09
 LARGEST INCREMENT OF DISP.         2.196E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.189E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE HEAT FLUX                  6.651E-04   TIME AVG. HEAT FLUX   3.279E-04
 LARGEST RESIDUAL HEAT FLUX        -2.706E-15   AT NODE       3599   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.380E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.646E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.434       TIME AVG. FORCE       0.233    
 LARGEST SCALED RESIDUAL FORCE     -2.680E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.680E-09
 LARGEST INCREMENT OF DISP.         2.218E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.219E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -5.447E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.694E-04   TIME AVG. HEAT FLUX   3.302E-04
 LARGEST RESIDUAL HEAT FLUX        -2.812E-15   AT NODE       4141   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.363E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.678E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

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
        SOLVER ELAPSED TIME:  352ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.436       TIME AVG. FORCE       0.234    
 LARGEST SCALED RESIDUAL FORCE     -2.830E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.830E-09
 LARGEST INCREMENT OF DISP.         2.241E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.242E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -5.743E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.737E-04   TIME AVG. HEAT FLUX   3.324E-04
 LARGEST RESIDUAL HEAT FLUX        -2.960E-15   AT NODE        643   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.345E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.775E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.438       TIME AVG. FORCE       0.236    
 LARGEST SCALED RESIDUAL FORCE     -2.936E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.936E-09
 LARGEST INCREMENT OF DISP.         2.263E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.268E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -5.962E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.779E-04   TIME AVG. HEAT FLUX   3.346E-04
 LARGEST RESIDUAL HEAT FLUX        -3.086E-15   AT NODE       3921   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.327E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.848E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

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
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.440       TIME AVG. FORCE       0.237    
 LARGEST SCALED RESIDUAL FORCE     -3.063E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.063E-09
 LARGEST INCREMENT OF DISP.         2.286E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.294E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -6.221E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.820E-04   TIME AVG. HEAT FLUX   3.369E-04
 LARGEST RESIDUAL HEAT FLUX        -2.871E-15   AT NODE       5283   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.307E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.932E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  338ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.442       TIME AVG. FORCE       0.238    
 LARGEST SCALED RESIDUAL FORCE     -3.187E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.187E-09
 LARGEST INCREMENT OF DISP.         2.310E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.320E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -6.477E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.862E-04   TIME AVG. HEAT FLUX   3.391E-04
 LARGEST RESIDUAL HEAT FLUX         3.165E-15   AT NODE       2775   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.287E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.014E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

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
        SOLVER ELAPSED TIME:  351ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.444       TIME AVG. FORCE       0.240    
 LARGEST SCALED RESIDUAL FORCE     -3.317E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.317E-09
 LARGEST INCREMENT OF DISP.         2.333E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.348E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -6.746E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.905E-04   TIME AVG. HEAT FLUX   3.413E-04
 LARGEST RESIDUAL HEAT FLUX         2.698E-15   AT NODE       3045   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.266E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.098E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  346ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.446       TIME AVG. FORCE       0.241    
 LARGEST SCALED RESIDUAL FORCE     -3.450E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.450E-09
 LARGEST INCREMENT OF DISP.         2.357E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.376E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -7.022E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.948E-04   TIME AVG. HEAT FLUX   3.435E-04
 LARGEST RESIDUAL HEAT FLUX         2.632E-15   AT NODE       2533   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.245E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.183E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

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
        SOLVER ELAPSED TIME:  347ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.448       TIME AVG. FORCE       0.242    
 LARGEST SCALED RESIDUAL FORCE     -3.587E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.587E-09
 LARGEST INCREMENT OF DISP.         2.381E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.406E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -7.309E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.989E-04   TIME AVG. HEAT FLUX   3.458E-04
 LARGEST RESIDUAL HEAT FLUX        -4.385E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.222E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.270E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  349ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.450       TIME AVG. FORCE       0.244    
 LARGEST SCALED RESIDUAL FORCE     -3.727E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.727E-09
 LARGEST INCREMENT OF DISP.         2.405E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.436E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -7.605E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.030E-04   TIME AVG. HEAT FLUX   3.480E-04
 LARGEST RESIDUAL HEAT FLUX         2.685E-15   AT NODE       1832   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.198E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.357E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

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
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.452       TIME AVG. FORCE       0.245    
 LARGEST SCALED RESIDUAL FORCE     -3.872E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.872E-09
 LARGEST INCREMENT OF DISP.         2.430E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.467E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -7.912E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.071E-04   TIME AVG. HEAT FLUX   3.502E-04
 LARGEST RESIDUAL HEAT FLUX        -3.287E-15   AT NODE       2238   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.174E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.447E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.454       TIME AVG. FORCE       0.246    
 LARGEST SCALED RESIDUAL FORCE     -4.020E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.020E-09
 LARGEST INCREMENT OF DISP.         2.455E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.498E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -8.230E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.109E-04   TIME AVG. HEAT FLUX   3.524E-04
 LARGEST RESIDUAL HEAT FLUX        -3.173E-15   AT NODE       3903   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.148E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.537E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

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
        SOLVER ELAPSED TIME:  355ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.456       TIME AVG. FORCE       0.247    
 LARGEST SCALED RESIDUAL FORCE     -4.173E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.173E-09
 LARGEST INCREMENT OF DISP.         2.480E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.531E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -8.560E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.151E-04   TIME AVG. HEAT FLUX   3.546E-04
 LARGEST RESIDUAL HEAT FLUX         7.272E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.122E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.629E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.458       TIME AVG. FORCE       0.249    
 LARGEST SCALED RESIDUAL FORCE     -4.330E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.330E-09
 LARGEST INCREMENT OF DISP.         2.506E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.565E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -8.902E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.190E-04   TIME AVG. HEAT FLUX   3.568E-04
 LARGEST RESIDUAL HEAT FLUX         2.772E-15   AT NODE       2185   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.095E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.722E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

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
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.460       TIME AVG. FORCE       0.250    
 LARGEST SCALED RESIDUAL FORCE     -4.492E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.492E-09
 LARGEST INCREMENT OF DISP.         2.532E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.600E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -9.258E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.229E-04   TIME AVG. HEAT FLUX   3.590E-04
 LARGEST RESIDUAL HEAT FLUX        -3.074E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.067E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.817E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  346ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.462       TIME AVG. FORCE       0.251    
 LARGEST SCALED RESIDUAL FORCE     -4.659E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.659E-09
 LARGEST INCREMENT OF DISP.         2.558E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.636E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.548E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.270E-04   TIME AVG. HEAT FLUX   3.612E-04
 LARGEST RESIDUAL HEAT FLUX         2.761E-15   AT NODE       9623   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.038E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.913E-04   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

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
        SOLVER ELAPSED TIME:  334ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.464       TIME AVG. FORCE       0.252    
 LARGEST SCALED RESIDUAL FORCE     -4.830E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.830E-09
 LARGEST INCREMENT OF DISP.         2.585E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.673E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.796E-10
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.307E-04   TIME AVG. HEAT FLUX   3.634E-04
 LARGEST RESIDUAL HEAT FLUX        -2.826E-15   AT NODE       3446   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.008E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.011E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  355ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.465       TIME AVG. FORCE       0.254    
 LARGEST SCALED RESIDUAL FORCE     -5.006E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.006E-09
 LARGEST INCREMENT OF DISP.         2.612E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.712E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.005E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.346E-04   TIME AVG. HEAT FLUX   3.656E-04
 LARGEST RESIDUAL HEAT FLUX        -3.400E-15   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.976E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.110E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.467       TIME AVG. FORCE       0.255    
 LARGEST SCALED RESIDUAL FORCE     -5.188E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.188E-09
 LARGEST INCREMENT OF DISP.         2.640E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.751E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.030E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.386E-04   TIME AVG. HEAT FLUX   3.678E-04
 LARGEST RESIDUAL HEAT FLUX         2.772E-15   AT NODE       1964   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.944E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.211E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.469       TIME AVG. FORCE       0.256    
 LARGEST SCALED RESIDUAL FORCE     -5.375E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.375E-09
 LARGEST INCREMENT OF DISP.         2.668E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.792E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.057E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.424E-04   TIME AVG. HEAT FLUX   3.700E-04
 LARGEST RESIDUAL HEAT FLUX        -2.584E-15   AT NODE       4635   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.911E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.313E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  345ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.471       TIME AVG. FORCE       0.257    
 LARGEST SCALED RESIDUAL FORCE     -5.568E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.568E-09
 LARGEST INCREMENT OF DISP.         2.696E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.835E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.083E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.464E-04   TIME AVG. HEAT FLUX   3.722E-04
 LARGEST RESIDUAL HEAT FLUX        -2.628E-15   AT NODE        886   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.877E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.416E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  342ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.473       TIME AVG. FORCE       0.259    
 LARGEST SCALED RESIDUAL FORCE     -5.766E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.766E-09
 LARGEST INCREMENT OF DISP.         2.725E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.878E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.110E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.500E-04   TIME AVG. HEAT FLUX   3.744E-04
 LARGEST RESIDUAL HEAT FLUX         2.948E-15   AT NODE       2777   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.842E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.521E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.474       TIME AVG. FORCE       0.260    
 LARGEST SCALED RESIDUAL FORCE     -5.971E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.971E-09
 LARGEST INCREMENT OF DISP.         2.754E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.923E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.138E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.538E-04   TIME AVG. HEAT FLUX   3.766E-04
 LARGEST RESIDUAL HEAT FLUX         3.456E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.806E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.628E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.476       TIME AVG. FORCE       0.261    
 LARGEST SCALED RESIDUAL FORCE     -6.183E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.183E-09
 LARGEST INCREMENT OF DISP.         2.784E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.970E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.166E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.574E-04   TIME AVG. HEAT FLUX   3.787E-04
 LARGEST RESIDUAL HEAT FLUX         2.611E-15   AT NODE       6545   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.768E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.737E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.478       TIME AVG. FORCE       0.262    
 LARGEST SCALED RESIDUAL FORCE     -6.401E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.401E-09
 LARGEST INCREMENT OF DISP.         2.814E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.018E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.194E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.612E-04   TIME AVG. HEAT FLUX   3.809E-04
 LARGEST RESIDUAL HEAT FLUX         2.936E-15   AT NODE        682   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.730E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.847E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  350ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.479       TIME AVG. FORCE       0.264    
 LARGEST SCALED RESIDUAL FORCE     -6.626E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.626E-09
 LARGEST INCREMENT OF DISP.         2.844E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.068E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.224E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.648E-04   TIME AVG. HEAT FLUX   3.831E-04
 LARGEST RESIDUAL HEAT FLUX         3.553E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.690E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.958E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  179ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.481       TIME AVG. FORCE       0.265    
 LARGEST SCALED RESIDUAL FORCE     -6.858E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.858E-09
 LARGEST INCREMENT OF DISP.         2.876E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.120E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.254E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.685E-04   TIME AVG. HEAT FLUX   3.852E-04
 LARGEST RESIDUAL HEAT FLUX        -3.553E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.649E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.072E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  190ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.483       TIME AVG. FORCE       0.266    
 LARGEST SCALED RESIDUAL FORCE     -7.098E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.098E-09
 LARGEST INCREMENT OF DISP.         2.907E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.174E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.284E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.724E-04   TIME AVG. HEAT FLUX   3.874E-04
 LARGEST RESIDUAL HEAT FLUX        -2.831E-15   AT NODE      12450   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.608E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.187E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  169ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.485       TIME AVG. FORCE       0.267    
 LARGEST SCALED RESIDUAL FORCE     -7.346E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.346E-09
 LARGEST INCREMENT OF DISP.         2.940E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.229E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.315E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.761E-04   TIME AVG. HEAT FLUX   3.896E-04
 LARGEST RESIDUAL HEAT FLUX        -2.580E-15   AT NODE       3828   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.564E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.304E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  169ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.486       TIME AVG. FORCE       0.268    
 LARGEST SCALED RESIDUAL FORCE     -7.602E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.602E-09
 LARGEST INCREMENT OF DISP.         2.973E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.287E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.347E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.793E-04   TIME AVG. HEAT FLUX   3.917E-04
 LARGEST RESIDUAL HEAT FLUX        -3.143E-15   AT NODE       3834   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.520E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.423E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  180ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.488       TIME AVG. FORCE       0.270    
 LARGEST SCALED RESIDUAL FORCE     -7.867E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.867E-09
 LARGEST INCREMENT OF DISP.         3.006E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.346E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.380E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.829E-04   TIME AVG. HEAT FLUX   3.939E-04
 LARGEST RESIDUAL HEAT FLUX         2.794E-15   AT NODE       2814   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.475E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.544E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  167ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.489       TIME AVG. FORCE       0.271    
 LARGEST SCALED RESIDUAL FORCE     -8.141E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.141E-09
 LARGEST INCREMENT OF DISP.         3.040E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.408E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.413E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.860E-04   TIME AVG. HEAT FLUX   3.960E-04
 LARGEST RESIDUAL HEAT FLUX         3.308E-15   AT NODE       2165   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.428E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.666E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  195ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.491       TIME AVG. FORCE       0.272    
 LARGEST SCALED RESIDUAL FORCE     -8.424E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.424E-09
 LARGEST INCREMENT OF DISP.         3.075E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.473E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.448E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.893E-04   TIME AVG. HEAT FLUX   3.982E-04
 LARGEST RESIDUAL HEAT FLUX        -3.007E-15   AT NODE       3418   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.380E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.791E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  177ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.493       TIME AVG. FORCE       0.273    
 LARGEST SCALED RESIDUAL FORCE     -8.718E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.718E-09
 LARGEST INCREMENT OF DISP.         3.110E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.539E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.483E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.929E-04   TIME AVG. HEAT FLUX   4.003E-04
 LARGEST RESIDUAL HEAT FLUX         3.097E-15   AT NODE       3674   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.331E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.918E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  182ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.494       TIME AVG. FORCE       0.274    
 LARGEST SCALED RESIDUAL FORCE     -9.023E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.023E-09
 LARGEST INCREMENT OF DISP.         3.146E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.609E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.519E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.962E-04   TIME AVG. HEAT FLUX   4.024E-04
 LARGEST RESIDUAL HEAT FLUX        -2.640E-15   AT NODE       4201   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.281E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.046E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  193ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.496       TIME AVG. FORCE       0.276    
 LARGEST SCALED RESIDUAL FORCE     -9.338E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.338E-09
 LARGEST INCREMENT OF DISP.         3.183E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.681E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.556E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.994E-04   TIME AVG. HEAT FLUX   4.045E-04
 LARGEST RESIDUAL HEAT FLUX        -3.281E-15   AT NODE        248   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.229E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.177E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  187ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.497       TIME AVG. FORCE       0.277    
 LARGEST SCALED RESIDUAL FORCE     -9.665E-09   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.665E-09
 LARGEST INCREMENT OF DISP.         3.221E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.755E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.594E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.027E-04   TIME AVG. HEAT FLUX   4.067E-04
 LARGEST RESIDUAL HEAT FLUX        -3.941E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.176E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.310E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  174ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.499       TIME AVG. FORCE       0.278    
 LARGEST SCALED RESIDUAL FORCE     -1.000E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.000E-08
 LARGEST INCREMENT OF DISP.         3.259E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.836E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.634E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.058E-04   TIME AVG. HEAT FLUX   4.088E-04
 LARGEST RESIDUAL HEAT FLUX         3.222E-15   AT NODE       2636   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.121E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.446E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  184ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.500       TIME AVG. FORCE       0.279    
 LARGEST SCALED RESIDUAL FORCE     -1.037E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.037E-08
 LARGEST INCREMENT OF DISP.         3.298E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.915E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.675E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.090E-04   TIME AVG. HEAT FLUX   4.109E-04
 LARGEST RESIDUAL HEAT FLUX         3.456E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.065E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.590E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  207ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.502       TIME AVG. FORCE       0.280    
 LARGEST SCALED RESIDUAL FORCE     -1.072E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.072E-08
 LARGEST INCREMENT OF DISP.         3.338E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.000E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.714E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.121E-04   TIME AVG. HEAT FLUX   4.130E-04
 LARGEST RESIDUAL HEAT FLUX         3.545E-15   AT NODE       4827   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         1.008E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.722E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  176ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.503       TIME AVG. FORCE       0.281    
 LARGEST SCALED RESIDUAL FORCE     -1.111E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.111E-08
 LARGEST INCREMENT OF DISP.         3.379E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.087E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.758E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.151E-04   TIME AVG. HEAT FLUX   4.151E-04
 LARGEST RESIDUAL HEAT FLUX        -3.103E-15   AT NODE        987   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         9.495E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.868E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  167ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.505       TIME AVG. FORCE       0.283    
 LARGEST SCALED RESIDUAL FORCE     -1.150E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.150E-08
 LARGEST INCREMENT OF DISP.         3.421E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.179E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.801E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.183E-04   TIME AVG. HEAT FLUX   4.172E-04
 LARGEST RESIDUAL HEAT FLUX         3.497E-15   AT NODE        177   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.893E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.012E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  188ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.506       TIME AVG. FORCE       0.284    
 LARGEST SCALED RESIDUAL FORCE     -1.192E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.192E-08
 LARGEST INCREMENT OF DISP.         3.464E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.275E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.847E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.212E-04   TIME AVG. HEAT FLUX   4.192E-04
 LARGEST RESIDUAL HEAT FLUX         2.700E-15   AT NODE       4315   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         8.277E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.161E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  166ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.508       TIME AVG. FORCE       0.285    
 LARGEST SCALED RESIDUAL FORCE     -1.234E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.234E-08
 LARGEST INCREMENT OF DISP.         3.507E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.375E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.894E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.238E-04   TIME AVG. HEAT FLUX   4.213E-04
 LARGEST RESIDUAL HEAT FLUX        -3.256E-15   AT NODE      12450   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         7.646E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.312E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  183ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.509       TIME AVG. FORCE       0.286    
 LARGEST SCALED RESIDUAL FORCE     -1.279E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.279E-08
 LARGEST INCREMENT OF DISP.         3.552E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.479E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.943E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.270E-04   TIME AVG. HEAT FLUX   4.234E-04
 LARGEST RESIDUAL HEAT FLUX         3.025E-15   AT NODE       4677   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.999E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.466E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  178ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.510       TIME AVG. FORCE       0.287    
 LARGEST SCALED RESIDUAL FORCE     -1.326E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.326E-08
 LARGEST INCREMENT OF DISP.         3.598E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.588E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.993E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.296E-04   TIME AVG. HEAT FLUX   4.254E-04
 LARGEST RESIDUAL HEAT FLUX         3.255E-15   AT NODE        668   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         6.337E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.623E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  181ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.512       TIME AVG. FORCE       0.288    
 LARGEST SCALED RESIDUAL FORCE     -1.375E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.375E-08
 LARGEST INCREMENT OF DISP.         3.645E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.702E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.045E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.326E-04   TIME AVG. HEAT FLUX   4.275E-04
 LARGEST RESIDUAL HEAT FLUX        -3.469E-15   AT NODE       6434   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.659E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.783E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  165ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.513       TIME AVG. FORCE       0.289    
 LARGEST SCALED RESIDUAL FORCE     -1.426E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.426E-08
 LARGEST INCREMENT OF DISP.         3.693E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.821E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.099E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.355E-04   TIME AVG. HEAT FLUX   4.296E-04
 LARGEST RESIDUAL HEAT FLUX        -3.136E-15   AT NODE       3142   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         5.097E-03   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.947E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  169ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.514       TIME AVG. FORCE       0.290    
 LARGEST SCALED RESIDUAL FORCE     -1.479E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.479E-08
 LARGEST INCREMENT OF DISP.         3.743E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.945E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.154E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.382E-04   TIME AVG. HEAT FLUX   4.316E-04
 LARGEST RESIDUAL HEAT FLUX        -2.725E-15   AT NODE       5030   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         4.764E-03   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.114E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  190ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.516       TIME AVG. FORCE       0.292    
 LARGEST SCALED RESIDUAL FORCE     -1.535E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.535E-08
 LARGEST INCREMENT OF DISP.         3.793E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.076E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.212E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.402E-04   TIME AVG. HEAT FLUX   4.336E-04
 LARGEST RESIDUAL HEAT FLUX         3.432E-15   AT NODE       1640   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         4.429E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.285E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  172ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.517       TIME AVG. FORCE       0.293    
 LARGEST SCALED RESIDUAL FORCE     -1.593E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.593E-08
 LARGEST INCREMENT OF DISP.         3.846E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.213E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.273E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.427E-04   TIME AVG. HEAT FLUX   4.357E-04
 LARGEST RESIDUAL HEAT FLUX        -2.483E-15   AT NODE        923   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         4.364E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.460E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.518       TIME AVG. FORCE       0.294    
 LARGEST SCALED RESIDUAL FORCE     -1.655E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.655E-08
 LARGEST INCREMENT OF DISP.         3.899E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.357E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.335E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.450E-04   TIME AVG. HEAT FLUX   4.377E-04
 LARGEST RESIDUAL HEAT FLUX         3.159E-15   AT NODE       1764   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         4.296E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.639E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.519       TIME AVG. FORCE       0.295    
 LARGEST SCALED RESIDUAL FORCE     -1.719E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.719E-08
 LARGEST INCREMENT OF DISP.         3.954E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.508E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.400E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.476E-04   TIME AVG. HEAT FLUX   4.397E-04
 LARGEST RESIDUAL HEAT FLUX        -2.780E-15   AT NODE       1048   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         4.227E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.822E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.520       TIME AVG. FORCE       0.296    
 LARGEST SCALED RESIDUAL FORCE     -1.786E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.786E-08
 LARGEST INCREMENT OF DISP.         4.011E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.666E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.468E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.506E-04   TIME AVG. HEAT FLUX   4.417E-04
 LARGEST RESIDUAL HEAT FLUX         3.830E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         4.156E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.009E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.522       TIME AVG. FORCE       0.297    
 LARGEST SCALED RESIDUAL FORCE     -1.857E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.857E-08
 LARGEST INCREMENT OF DISP.         4.069E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.833E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.539E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.524E-04   TIME AVG. HEAT FLUX   4.437E-04
 LARGEST RESIDUAL HEAT FLUX         3.088E-15   AT NODE       4536   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         4.083E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.202E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.523       TIME AVG. FORCE       0.298    
 LARGEST SCALED RESIDUAL FORCE     -1.932E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.932E-08
 LARGEST INCREMENT OF DISP.         4.129E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.008E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.613E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.546E-04   TIME AVG. HEAT FLUX   4.457E-04
 LARGEST RESIDUAL HEAT FLUX        -2.731E-15   AT NODE       1544   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         4.008E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.399E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.524       TIME AVG. FORCE       0.299    
 LARGEST SCALED RESIDUAL FORCE     -2.011E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.011E-08
 LARGEST INCREMENT OF DISP.         4.191E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.193E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.691E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.566E-04   TIME AVG. HEAT FLUX   4.476E-04
 LARGEST RESIDUAL HEAT FLUX         3.941E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         3.931E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.602E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.525       TIME AVG. FORCE       0.300    
 LARGEST SCALED RESIDUAL FORCE     -2.094E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.094E-08
 LARGEST INCREMENT OF DISP.         4.255E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.388E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.772E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.585E-04   TIME AVG. HEAT FLUX   4.496E-04
 LARGEST RESIDUAL HEAT FLUX        -6.106E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         3.853E-03   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.811E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.526       TIME AVG. FORCE       0.301    
 LARGEST SCALED RESIDUAL FORCE     -2.182E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.182E-08
 LARGEST INCREMENT OF DISP.         4.321E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.594E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.857E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.607E-04   TIME AVG. HEAT FLUX   4.516E-04
 LARGEST RESIDUAL HEAT FLUX        -2.949E-15   AT NODE       6434   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.873E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.025E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  353ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.527       TIME AVG. FORCE       0.302    
 LARGEST SCALED RESIDUAL FORCE     -2.275E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.275E-08
 LARGEST INCREMENT OF DISP.         4.389E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.811E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.946E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.628E-04   TIME AVG. HEAT FLUX   4.535E-04
 LARGEST RESIDUAL HEAT FLUX         3.442E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.797E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.246E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  345ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.528       TIME AVG. FORCE       0.303    
 LARGEST SCALED RESIDUAL FORCE     -2.373E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.373E-08
 LARGEST INCREMENT OF DISP.         4.460E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.041E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.040E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.652E-04   TIME AVG. HEAT FLUX   4.555E-04
 LARGEST RESIDUAL HEAT FLUX        -5.065E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.745E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.474E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.529       TIME AVG. FORCE       0.305    
 LARGEST SCALED RESIDUAL FORCE     -2.478E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.478E-08
 LARGEST INCREMENT OF DISP.         4.532E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.284E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.138E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.666E-04   TIME AVG. HEAT FLUX   4.574E-04
 LARGEST RESIDUAL HEAT FLUX         3.303E-15   AT NODE      11986   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -6.715E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.709E-04   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.530       TIME AVG. FORCE       0.306    
 LARGEST SCALED RESIDUAL FORCE     -2.588E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.588E-08
 LARGEST INCREMENT OF DISP.         4.608E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.541E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.242E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.689E-04   TIME AVG. HEAT FLUX   4.593E-04
 LARGEST RESIDUAL HEAT FLUX         2.921E-15   AT NODE       4565   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -7.711E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.952E-04   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  352ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.531       TIME AVG. FORCE       0.307    
 LARGEST SCALED RESIDUAL FORCE     -2.706E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.706E-08
 LARGEST INCREMENT OF DISP.         4.686E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.814E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.352E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.700E-04   TIME AVG. HEAT FLUX   4.612E-04
 LARGEST RESIDUAL HEAT FLUX        -2.877E-15   AT NODE       2777   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -8.731E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.020E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  348ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.532       TIME AVG. FORCE       0.308    
 LARGEST SCALED RESIDUAL FORCE     -2.832E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.832E-08
 LARGEST INCREMENT OF DISP.         4.767E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.104E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.467E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.721E-04   TIME AVG. HEAT FLUX   4.631E-04
 LARGEST RESIDUAL HEAT FLUX        -3.275E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -9.777E-03   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.046E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.533       TIME AVG. FORCE       0.309    
 LARGEST SCALED RESIDUAL FORCE     -2.966E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.966E-08
 LARGEST INCREMENT OF DISP.         4.851E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.412E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.590E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.732E-04   TIME AVG. HEAT FLUX   4.650E-04
 LARGEST RESIDUAL HEAT FLUX         2.570E-15   AT NODE       3186   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.085E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.073E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  352ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.534       TIME AVG. FORCE       0.310    
 LARGEST SCALED RESIDUAL FORCE     -3.109E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.109E-08
 LARGEST INCREMENT OF DISP.         4.938E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.739E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.719E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.749E-04   TIME AVG. HEAT FLUX   4.669E-04
 LARGEST RESIDUAL HEAT FLUX         2.749E-15   AT NODE       1724   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.195E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.101E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  352ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.534       TIME AVG. FORCE       0.311    
 LARGEST SCALED RESIDUAL FORCE     -3.262E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.262E-08
 LARGEST INCREMENT OF DISP.         5.029E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.088E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.857E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.756E-04   TIME AVG. HEAT FLUX   4.688E-04
 LARGEST RESIDUAL HEAT FLUX         2.670E-15   AT NODE       4208   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.308E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.130E-03   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.535       TIME AVG. FORCE       0.312    
 LARGEST SCALED RESIDUAL FORCE     -3.426E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.426E-08
 LARGEST INCREMENT OF DISP.         5.124E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.461E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.003E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.770E-04   TIME AVG. HEAT FLUX   4.706E-04
 LARGEST RESIDUAL HEAT FLUX         3.664E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.424E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.160E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.536       TIME AVG. FORCE       0.313    
 LARGEST SCALED RESIDUAL FORCE     -3.602E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.602E-08
 LARGEST INCREMENT OF DISP.         5.223E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.858E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.159E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.780E-04   TIME AVG. HEAT FLUX   4.725E-04
 LARGEST RESIDUAL HEAT FLUX         3.032E-15   AT NODE       1879   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.543E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.192E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  372ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.537       TIME AVG. FORCE       0.314    
 LARGEST SCALED RESIDUAL FORCE     -3.792E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.792E-08
 LARGEST INCREMENT OF DISP.         5.325E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.028E-04   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.325E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.784E-04   TIME AVG. HEAT FLUX   4.743E-04
 LARGEST RESIDUAL HEAT FLUX        -5.163E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.666E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.225E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  368ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.538       TIME AVG. FORCE       0.315    
 LARGEST SCALED RESIDUAL FORCE     -3.997E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.997E-08
 LARGEST INCREMENT OF DISP.         5.433E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.074E-04   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.504E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.794E-04   TIME AVG. HEAT FLUX   4.761E-04
 LARGEST RESIDUAL HEAT FLUX        -2.606E-15   AT NODE       3110   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.792E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.259E-03   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.538       TIME AVG. FORCE       0.316    
 LARGEST SCALED RESIDUAL FORCE     -4.218E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.218E-08
 LARGEST INCREMENT OF DISP.         5.545E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.124E-04   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.696E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.816E-04   TIME AVG. HEAT FLUX   4.779E-04
 LARGEST RESIDUAL HEAT FLUX         3.639E-15   AT NODE       2697   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -1.921E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.296E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.539       TIME AVG. FORCE       0.317    
 LARGEST SCALED RESIDUAL FORCE     -4.458E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.458E-08
 LARGEST INCREMENT OF DISP.         5.663E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.177E-04   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.903E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.825E-04   TIME AVG. HEAT FLUX   4.797E-04
 LARGEST RESIDUAL HEAT FLUX         4.718E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.055E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.333E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.540       TIME AVG. FORCE       0.318    
 LARGEST SCALED RESIDUAL FORCE     -4.719E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.719E-08
 LARGEST INCREMENT OF DISP.         5.786E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.234E-04   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.126E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.830E-04   TIME AVG. HEAT FLUX   4.815E-04
 LARGEST RESIDUAL HEAT FLUX         2.868E-15   AT NODE       1116   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.192E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.373E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.540       TIME AVG. FORCE       0.319    
 LARGEST SCALED RESIDUAL FORCE     -5.004E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.004E-08
 LARGEST INCREMENT OF DISP.         5.915E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.295E-04   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.367E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.832E-04   TIME AVG. HEAT FLUX   4.833E-04
 LARGEST RESIDUAL HEAT FLUX         3.747E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.334E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.415E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  351ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.541       TIME AVG. FORCE       0.320    
 LARGEST SCALED RESIDUAL FORCE     -5.315E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.315E-08
 LARGEST INCREMENT OF DISP.         6.051E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.362E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.629E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.836E-04   TIME AVG. HEAT FLUX   4.850E-04
 LARGEST RESIDUAL HEAT FLUX         4.038E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.480E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.460E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.541       TIME AVG. FORCE       0.321    
 LARGEST SCALED RESIDUAL FORCE     -5.656E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.656E-08
 LARGEST INCREMENT OF DISP.         6.194E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.434E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.914E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.846E-04   TIME AVG. HEAT FLUX   4.868E-04
 LARGEST RESIDUAL HEAT FLUX         3.150E-15   AT NODE        414   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.630E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.507E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  351ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.542       TIME AVG. FORCE       0.322    
 LARGEST SCALED RESIDUAL FORCE     -6.032E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.032E-08
 LARGEST INCREMENT OF DISP.         6.345E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.513E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.226E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.852E-04   TIME AVG. HEAT FLUX   4.885E-04
 LARGEST RESIDUAL HEAT FLUX         3.078E-15   AT NODE       1899   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.786E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.557E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  355ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.543       TIME AVG. FORCE       0.323    
 LARGEST SCALED RESIDUAL FORCE     -6.447E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.447E-08
 LARGEST INCREMENT OF DISP.         6.504E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.598E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.567E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.857E-04   TIME AVG. HEAT FLUX   4.902E-04
 LARGEST RESIDUAL HEAT FLUX        -4.201E-15   AT NODE       4061   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.947E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.611E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.543       TIME AVG. FORCE       0.324    
 LARGEST SCALED RESIDUAL FORCE     -6.907E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.907E-08
 LARGEST INCREMENT OF DISP.         6.672E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.691E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.941E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.859E-04   TIME AVG. HEAT FLUX   4.919E-04
 LARGEST RESIDUAL HEAT FLUX         3.185E-15   AT NODE       6434   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.114E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.668E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  355ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.543       TIME AVG. FORCE       0.324    
 LARGEST SCALED RESIDUAL FORCE     -7.419E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.419E-08
 LARGEST INCREMENT OF DISP.         6.849E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.791E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.352E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.865E-04   TIME AVG. HEAT FLUX   4.936E-04
 LARGEST RESIDUAL HEAT FLUX        -3.744E-15   AT NODE       3277   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.287E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.729E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  359ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.544       TIME AVG. FORCE       0.325    
 LARGEST SCALED RESIDUAL FORCE     -7.988E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.988E-08
 LARGEST INCREMENT OF DISP.         7.038E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.901E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.803E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.864E-04   TIME AVG. HEAT FLUX   4.953E-04
 LARGEST RESIDUAL HEAT FLUX         3.719E-15   AT NODE       6434   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.466E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.794E-03   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.544       TIME AVG. FORCE       0.326    
 LARGEST SCALED RESIDUAL FORCE     -8.625E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.625E-08
 LARGEST INCREMENT OF DISP.         7.238E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.021E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    8.302E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.869E-04   TIME AVG. HEAT FLUX   4.970E-04
 LARGEST RESIDUAL HEAT FLUX         2.465E-15   AT NODE       2012   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.652E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.863E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.544       TIME AVG. FORCE       0.327    
 LARGEST SCALED RESIDUAL FORCE     -9.339E-08   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.339E-08
 LARGEST INCREMENT OF DISP.         7.451E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.152E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    8.854E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.868E-04   TIME AVG. HEAT FLUX   4.986E-04
 LARGEST RESIDUAL HEAT FLUX        -4.163E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.846E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.938E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.545       TIME AVG. FORCE       0.328    
 LARGEST SCALED RESIDUAL FORCE     -1.014E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.014E-07
 LARGEST INCREMENT OF DISP.         7.679E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.296E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    9.466E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.872E-04   TIME AVG. HEAT FLUX   5.003E-04
 LARGEST RESIDUAL HEAT FLUX         2.960E-15   AT NODE       1247   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.048E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.018E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.545       TIME AVG. FORCE       0.329    
 LARGEST SCALED RESIDUAL FORCE     -1.105E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.105E-07
 LARGEST INCREMENT OF DISP.         7.921E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.454E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.014E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.874E-04   TIME AVG. HEAT FLUX   5.019E-04
 LARGEST RESIDUAL HEAT FLUX         4.607E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.258E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.103E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.545       TIME AVG. FORCE       0.330    
 LARGEST SCALED RESIDUAL FORCE     -1.208E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.208E-07
 LARGEST INCREMENT OF DISP.         8.180E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.628E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.090E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.876E-04   TIME AVG. HEAT FLUX   5.035E-04
 LARGEST RESIDUAL HEAT FLUX         5.856E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.478E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.196E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.545       TIME AVG. FORCE       0.331    
 LARGEST SCALED RESIDUAL FORCE     -1.324E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.324E-07
 LARGEST INCREMENT OF DISP.         8.458E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.820E-04   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.175E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.874E-04   TIME AVG. HEAT FLUX   5.051E-04
 LARGEST RESIDUAL HEAT FLUX        -3.858E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.707E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.294E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  371ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.545       TIME AVG. FORCE       0.332    
 LARGEST SCALED RESIDUAL FORCE     -1.457E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.457E-07
 LARGEST INCREMENT OF DISP.         8.756E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.034E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.269E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.864E-04   TIME AVG. HEAT FLUX   5.067E-04
 LARGEST RESIDUAL HEAT FLUX         3.164E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.947E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.400E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  376ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.545       TIME AVG. FORCE       0.333    
 LARGEST SCALED RESIDUAL FORCE     -1.610E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.610E-07
 LARGEST INCREMENT OF DISP.         9.077E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.270E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.377E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.843E-04   TIME AVG. HEAT FLUX   5.082E-04
 LARGEST RESIDUAL HEAT FLUX        -3.567E-15   AT NODE        177   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.199E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.514E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.545       TIME AVG. FORCE       0.333    
 LARGEST SCALED RESIDUAL FORCE     -1.784E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.784E-07
 LARGEST INCREMENT OF DISP.         9.423E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.534E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.497E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.828E-04   TIME AVG. HEAT FLUX   5.098E-04
 LARGEST RESIDUAL HEAT FLUX        -2.899E-15   AT NODE       2085   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.462E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.635E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  371ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.545       TIME AVG. FORCE       0.334    
 LARGEST SCALED RESIDUAL FORCE     -1.987E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.987E-07
 LARGEST INCREMENT OF DISP.         9.797E-03   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        3.828E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.634E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.812E-04   TIME AVG. HEAT FLUX   5.113E-04
 LARGEST RESIDUAL HEAT FLUX        -5.662E-15   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.739E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.768E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  375ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.544       TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -2.222E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.222E-07
 LARGEST INCREMENT OF DISP.         1.020E-02   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.158E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.789E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.802E-04   TIME AVG. HEAT FLUX   5.128E-04
 LARGEST RESIDUAL HEAT FLUX        -3.127E-15   AT NODE       2058   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -6.030E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.909E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.544       TIME AVG. FORCE       0.336    
 LARGEST SCALED RESIDUAL FORCE     -2.498E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.498E-07
 LARGEST INCREMENT OF DISP.         1.064E-02   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.526E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.966E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.795E-04   TIME AVG. HEAT FLUX   5.143E-04
 LARGEST RESIDUAL HEAT FLUX        -4.441E-15   AT NODE        176   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -6.336E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.063E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.543       TIME AVG. FORCE       0.337    
 LARGEST SCALED RESIDUAL FORCE     -2.818E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.818E-07
 LARGEST INCREMENT OF DISP.         1.112E-02   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        4.940E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.167E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.779E-04   TIME AVG. HEAT FLUX   5.158E-04
 LARGEST RESIDUAL HEAT FLUX         4.108E-15   AT NODE        550   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -6.659E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.224E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.543       TIME AVG. FORCE       0.338    
 LARGEST SCALED RESIDUAL FORCE     -3.198E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.198E-07
 LARGEST INCREMENT OF DISP.         1.165E-02   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.406E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.398E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.762E-04   TIME AVG. HEAT FLUX   5.172E-04
 LARGEST RESIDUAL HEAT FLUX         2.755E-15   AT NODE       3251   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -6.998E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.398E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  419ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.542       TIME AVG. FORCE       0.339    
 LARGEST SCALED RESIDUAL FORCE     -3.646E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.646E-07
 LARGEST INCREMENT OF DISP.         1.224E-02   AT NODE        176   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        5.931E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.663E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.749E-04   TIME AVG. HEAT FLUX   5.187E-04
 LARGEST RESIDUAL HEAT FLUX        -3.090E-15   AT NODE      15932   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -7.357E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.583E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  368ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.541       TIME AVG. FORCE       0.339    
 LARGEST SCALED RESIDUAL FORCE     -4.183E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.183E-07
 LARGEST INCREMENT OF DISP.         1.288E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        6.521E-04   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.971E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.747E-04   TIME AVG. HEAT FLUX   5.201E-04
 LARGEST RESIDUAL HEAT FLUX         2.681E-15   AT NODE       3633   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -7.735E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.781E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.540       TIME AVG. FORCE       0.340    
 LARGEST SCALED RESIDUAL FORCE     -4.822E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.822E-07
 LARGEST INCREMENT OF DISP.         1.360E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.192E-04   AT NODE       6446   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.325E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.746E-04   TIME AVG. HEAT FLUX   5.215E-04
 LARGEST RESIDUAL HEAT FLUX         2.882E-15   AT NODE       2259   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -8.133E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -3.984E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.539       TIME AVG. FORCE       0.341    
 LARGEST SCALED RESIDUAL FORCE     -5.589E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.589E-07
 LARGEST INCREMENT OF DISP.         1.439E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        7.949E-04   AT NODE       6446   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.659E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.754E-04   TIME AVG. HEAT FLUX   5.229E-04
 LARGEST RESIDUAL HEAT FLUX        -2.658E-15   AT NODE       4921   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -8.553E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.193E-03   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.538       TIME AVG. FORCE       0.342    
 LARGEST SCALED RESIDUAL FORCE     -6.511E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.511E-07
 LARGEST INCREMENT OF DISP.         1.526E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        8.799E-04   AT NODE       6446   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.006E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.762E-04   TIME AVG. HEAT FLUX   5.243E-04
 LARGEST RESIDUAL HEAT FLUX         3.198E-15   AT NODE       2293   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -8.993E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.401E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.536       TIME AVG. FORCE       0.342    
 LARGEST SCALED RESIDUAL FORCE     -7.634E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.634E-07
 LARGEST INCREMENT OF DISP.         1.622E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        9.754E-04   AT NODE        178   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.388E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.755E-04   TIME AVG. HEAT FLUX   5.257E-04
 LARGEST RESIDUAL HEAT FLUX        -3.038E-15   AT NODE        177   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -9.453E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.606E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.535       TIME AVG. FORCE       0.343    
 LARGEST SCALED RESIDUAL FORCE     -8.980E-07   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.980E-07
 LARGEST INCREMENT OF DISP.         1.728E-02   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.082E-03   AT NODE        178   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    4.790E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.762E-04   TIME AVG. HEAT FLUX   5.271E-04
 LARGEST RESIDUAL HEAT FLUX         2.850E-15   AT NODE       4211   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -9.932E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.783E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.533       TIME AVG. FORCE       0.344    
 LARGEST SCALED RESIDUAL FORCE     -1.063E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.063E-06
 LARGEST INCREMENT OF DISP.         1.844E-02   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.196E-03   AT NODE        178   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.200E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.781E-04   TIME AVG. HEAT FLUX   5.284E-04
 LARGEST RESIDUAL HEAT FLUX        -4.226E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -0.104       AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.936E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  367ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.531       TIME AVG. FORCE       0.345    
 LARGEST SCALED RESIDUAL FORCE     -1.252E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.252E-06
 LARGEST INCREMENT OF DISP.         1.973E-02   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.325E-03   AT NODE       6438   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.577E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.804E-04   TIME AVG. HEAT FLUX   5.298E-04
 LARGEST RESIDUAL HEAT FLUX        -2.798E-15   AT NODE       3797   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -0.109       AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.978E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  352ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.528       TIME AVG. FORCE       0.345    
 LARGEST SCALED RESIDUAL FORCE     -1.499E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.499E-06
 LARGEST INCREMENT OF DISP.         2.112E-02   AT NODE        177   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.457E-03   AT NODE        179   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.008E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.851E-04   TIME AVG. HEAT FLUX   5.312E-04
 LARGEST RESIDUAL HEAT FLUX        -3.187E-15   AT NODE       4101   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -0.114       AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.010E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  368ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.526       TIME AVG. FORCE       0.346    
 LARGEST SCALED RESIDUAL FORCE     -1.763E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.763E-06
 LARGEST INCREMENT OF DISP.         2.262E-02   AT NODE       6446   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.597E-03   AT NODE        180   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.283E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.863E-04   TIME AVG. HEAT FLUX   5.325E-04
 LARGEST RESIDUAL HEAT FLUX         3.338E-15   AT NODE       3126   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -0.119       AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -4.791E-03   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  353ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.522       TIME AVG. FORCE       0.347    
 LARGEST SCALED RESIDUAL FORCE     -2.095E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.095E-06
 LARGEST INCREMENT OF DISP.         2.421E-02   AT NODE       6446   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.739E-03   AT NODE       6452   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.579E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.921E-04   TIME AVG. HEAT FLUX   5.339E-04
 LARGEST RESIDUAL HEAT FLUX         3.920E-15   AT NODE       6707   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -0.124       AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.026E-03   AT NODE       6707   DOF 11
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
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.519       TIME AVG. FORCE       0.347    
 LARGEST SCALED RESIDUAL FORCE     -2.439E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.439E-06
 LARGEST INCREMENT OF DISP.         2.584E-02   AT NODE        178   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        1.880E-03   AT NODE        182   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.674E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  9.002E-04   TIME AVG. HEAT FLUX   5.353E-04
 LARGEST RESIDUAL HEAT FLUX         2.829E-15   AT NODE      15477   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -0.127       AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.786E-03   AT NODE        176   DOF 11
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


 AVERAGE FORCE                      0.515       TIME AVG. FORCE       0.348    
 LARGEST SCALED RESIDUAL FORCE     -2.741E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.741E-06
 LARGEST INCREMENT OF DISP.         2.742E-02   AT NODE       6438   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.022E-03   AT NODE        183   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.521E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  9.108E-04   TIME AVG. HEAT FLUX   5.368E-04
 LARGEST RESIDUAL HEAT FLUX         3.084E-15   AT NODE        794   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -0.130       AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.918E-03   AT NODE        176   DOF 11
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
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.510       TIME AVG. FORCE       0.349    
 LARGEST SCALED RESIDUAL FORCE     -2.958E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -2.958E-06
 LARGEST INCREMENT OF DISP.         2.890E-02   AT NODE        179   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.162E-03   AT NODE       5723   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.150E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  9.266E-04   TIME AVG. HEAT FLUX   5.382E-04
 LARGEST RESIDUAL HEAT FLUX        -2.650E-15   AT NODE        900   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -0.131       AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -8.130E-03   AT NODE        176   DOF 11
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
        SOLVER ELAPSED TIME:  356ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.505       TIME AVG. FORCE       0.349    
 LARGEST SCALED RESIDUAL FORCE     -3.081E-06   AT NODE      11984   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.081E-06
 LARGEST INCREMENT OF DISP.         3.024E-02   AT NODE        180   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.295E-03   AT NODE        186   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.729E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  9.480E-04   TIME AVG. HEAT FLUX   5.398E-04
 LARGEST RESIDUAL HEAT FLUX        -2.754E-15   AT NODE      15021   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -0.130       AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -9.345E-03   AT NODE        176   DOF 11
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
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.499       TIME AVG. FORCE       0.350    
 LARGEST SCALED RESIDUAL FORCE     -3.617E-06   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.617E-06
 LARGEST INCREMENT OF DISP.         3.146E-02   AT NODE       6455   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.        2.415E-03   AT NODE       6458   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    6.259E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  9.754E-04   TIME AVG. HEAT FLUX   5.414E-04
 LARGEST RESIDUAL HEAT FLUX         2.784E-15   AT NODE       1263   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -0.126       AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.033E-02   AT NODE        177   DOF 11
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
        SOLVER ELAPSED TIME:  378ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.493       TIME AVG. FORCE       0.350    
 LARGEST SCALED RESIDUAL FORCE     -4.094E-06   AT NODE       6434   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.094E-06
 LARGEST INCREMENT OF DISP.         3.261E-02   AT NODE       6447   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -2.637E-03   AT NODE       3126   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    7.206E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.009E-03   TIME AVG. HEAT FLUX   5.432E-04
 LARGEST RESIDUAL HEAT FLUX         3.156E-15   AT NODE        388   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -0.119       AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.165E-02   AT NODE       6446   DOF 11
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
        SOLVER ELAPSED TIME:  344ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.486       TIME AVG. FORCE       0.351    
 LARGEST SCALED RESIDUAL FORCE      4.365E-06   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.365E-06
 LARGEST INCREMENT OF DISP.         3.370E-02   AT NODE        184   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -3.924E-03   AT NODE       6435   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.222E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.051E-03   TIME AVG. HEAT FLUX   5.451E-04
 LARGEST RESIDUAL HEAT FLUX        -2.777E-15   AT NODE       3415   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -0.109       AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.197E-02   AT NODE       6438   DOF 11
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
        SOLVER ELAPSED TIME:  352ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.479       TIME AVG. FORCE       0.351    
 LARGEST SCALED RESIDUAL FORCE      4.489E-06   AT NODE          8   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      4.489E-06
 LARGEST INCREMENT OF DISP.         3.463E-02   AT NODE       6902   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.251E-03   AT NODE       6435   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.174E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.098E-03   TIME AVG. HEAT FLUX   5.472E-04
 LARGEST RESIDUAL HEAT FLUX        -2.829E-15   AT NODE       3884   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -9.568E-02   AT NODE          8   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.284E-02   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  342ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.472       TIME AVG. FORCE       0.352    
 LARGEST SCALED RESIDUAL FORCE      7.031E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      7.031E-06
 LARGEST INCREMENT OF DISP.         3.536E-02   AT NODE       6907   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.757E-03   AT NODE       6441   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -5.525E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.147E-03   TIME AVG. HEAT FLUX   5.494E-04
 LARGEST RESIDUAL HEAT FLUX         2.663E-15   AT NODE       3894   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -8.419E-02   AT NODE       2548   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.446E-02   AT NODE          8   DOF 11
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
        SOLVER ELAPSED TIME:  349ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.465       TIME AVG. FORCE       0.352    
 LARGEST SCALED RESIDUAL FORCE      1.114E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.114E-05
 LARGEST INCREMENT OF DISP.         3.577E-02   AT NODE          1   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.667E-03   AT NODE       6439   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.190E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.193E-03   TIME AVG. HEAT FLUX   5.518E-04
 LARGEST RESIDUAL HEAT FLUX        -3.218E-15   AT NODE       4565   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -7.959E-02   AT NODE       2548   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.547E-02   AT NODE          8   DOF 11
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


 AVERAGE FORCE                      0.458       TIME AVG. FORCE       0.353    
 LARGEST SCALED RESIDUAL FORCE      1.424E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.424E-05
 LARGEST INCREMENT OF DISP.         3.591E-02   AT NODE        188   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.491E-03   AT NODE       6451   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.537E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.231E-03   TIME AVG. HEAT FLUX   5.543E-04
 LARGEST RESIDUAL HEAT FLUX        -2.838E-15   AT NODE        202   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -7.046E-02   AT NODE        364   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.802E-02   AT NODE       3126   DOF 11
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
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.451       TIME AVG. FORCE       0.353    
 LARGEST SCALED RESIDUAL FORCE      1.512E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.512E-05
 LARGEST INCREMENT OF DISP.         3.587E-02   AT NODE       6469   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.299E-03   AT NODE      12273   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.084E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.261E-03   TIME AVG. HEAT FLUX   5.569E-04
 LARGEST RESIDUAL HEAT FLUX        -2.779E-15   AT NODE       2935   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.935E-02   AT NODE        365   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.174E-02   AT NODE       3126   DOF 11
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
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.444       TIME AVG. FORCE       0.353    
 LARGEST SCALED RESIDUAL FORCE      1.427E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.427E-05
 LARGEST INCREMENT OF DISP.         3.579E-02   AT NODE       6464   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.165E-03   AT NODE        376   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.257E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.286E-03   TIME AVG. HEAT FLUX   5.595E-04
 LARGEST RESIDUAL HEAT FLUX         3.219E-15   AT NODE       3658   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -5.068E-02   AT NODE       6451   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.724E-02   AT NODE       3127   DOF 11
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
        SOLVER ELAPSED TIME:  356ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.437       TIME AVG. FORCE       0.354    
 LARGEST SCALED RESIDUAL FORCE      1.163E-05   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      1.163E-05
 LARGEST INCREMENT OF DISP.         3.566E-02   AT NODE        192   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.094E-03   AT NODE       3132   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.883E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.306E-03   TIME AVG. HEAT FLUX   5.623E-04
 LARGEST RESIDUAL HEAT FLUX        -3.185E-15   AT NODE       3278   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.589E-02   AT NODE        368   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.360E-02   AT NODE       3127   DOF 11
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
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.431       TIME AVG. FORCE       0.354    
 LARGEST SCALED RESIDUAL FORCE      8.483E-06   AT NODE       6707   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE      8.483E-06
 LARGEST INCREMENT OF DISP.         3.555E-02   AT NODE        193   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.015E-03   AT NODE        378   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.391E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.323E-03   TIME AVG. HEAT FLUX   5.650E-04
 LARGEST RESIDUAL HEAT FLUX         3.063E-15   AT NODE       2214   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.260E-02   AT NODE        376   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        9.799E-03   AT NODE       3128   DOF 11
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
        SOLVER ELAPSED TIME:  373ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.424       TIME AVG. FORCE       0.354    
 LARGEST SCALED RESIDUAL FORCE     -1.739E-05   AT NODE       6450   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.739E-05
 LARGEST INCREMENT OF DISP.         3.544E-02   AT NODE        194   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.982E-03   AT NODE      11215   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.890E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.339E-03   TIME AVG. HEAT FLUX   5.678E-04
 LARGEST RESIDUAL HEAT FLUX        -2.610E-15   AT NODE       4285   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.040E-02   AT NODE       5721   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.577E-03   AT NODE       3129   DOF 11
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
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.417       TIME AVG. FORCE       0.354    
 LARGEST SCALED RESIDUAL FORCE     -6.440E-05   AT NODE       6448   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.440E-05
 LARGEST INCREMENT OF DISP.         3.535E-02   AT NODE       6888   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.918E-03   AT NODE       6460   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.087E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.353E-03   TIME AVG. HEAT FLUX   5.707E-04
 LARGEST RESIDUAL HEAT FLUX        -2.925E-15   AT NODE       4042   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.866E-02   AT NODE       6904   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.601E-03   AT NODE       3130   DOF 11
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
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.411       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -1.455E-04   AT NODE       6457   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.455E-04
 LARGEST INCREMENT OF DISP.         3.529E-02   AT NODE       6480   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.915E-03   AT NODE      12277   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.498E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.365E-03   TIME AVG. HEAT FLUX   5.735E-04
 LARGEST RESIDUAL HEAT FLUX         3.067E-15   AT NODE        609   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.746E-02   AT NODE       6875   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.257E-03   AT NODE       3050   DOF 11
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
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.404       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -3.135E-04   AT NODE       5722   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.135E-04
 LARGEST INCREMENT OF DISP.         3.524E-02   AT NODE        198   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.876E-03   AT NODE       6927   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    5.458E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.376E-03   TIME AVG. HEAT FLUX   5.764E-04
 LARGEST RESIDUAL HEAT FLUX        -3.095E-15   AT NODE       3350   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.644E-02   AT NODE       6460   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.764E-03   AT NODE       3132   DOF 11
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
        SOLVER ELAPSED TIME:  365ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.398       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -4.824E-04   AT NODE       5724   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.824E-04
 LARGEST INCREMENT OF DISP.         3.524E-02   AT NODE          2   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.864E-03   AT NODE      12008   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    8.553E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.386E-03   TIME AVG. HEAT FLUX   5.793E-04
 LARGEST RESIDUAL HEAT FLUX         3.132E-15   AT NODE       3608   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.569E-02   AT NODE       6463   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.657E-03   AT NODE       3051   DOF 11
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
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.392       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -6.556E-04   AT NODE       5724   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.556E-04
 LARGEST INCREMENT OF DISP.         3.522E-02   AT NODE        199   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.880E-03   AT NODE       6887   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.177E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.395E-03   TIME AVG. HEAT FLUX   5.822E-04
 LARGEST RESIDUAL HEAT FLUX         2.971E-15   AT NODE       4942   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.506E-02   AT NODE        369   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.347E-03   AT NODE       3052   DOF 11
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
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.385       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -7.013E-04   AT NODE       6903   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.013E-04
 LARGEST INCREMENT OF DISP.         3.527E-02   AT NODE       6488   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.864E-03   AT NODE       3134   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.257E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.404E-03   TIME AVG. HEAT FLUX   5.851E-04
 LARGEST RESIDUAL HEAT FLUX        -2.699E-15   AT NODE       1343   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.457E-02   AT NODE        370   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.369E-03   AT NODE       3133   DOF 11
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
        SOLVER ELAPSED TIME:  355ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.379       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -7.114E-04   AT NODE       6461   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.114E-04
 LARGEST INCREMENT OF DISP.         3.529E-02   AT NODE       6493   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.842E-03   AT NODE       6475   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.276E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.411E-03   TIME AVG. HEAT FLUX   5.880E-04
 LARGEST RESIDUAL HEAT FLUX         2.993E-15   AT NODE        758   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.422E-02   AT NODE       6471   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.162E-03   AT NODE       3131   DOF 11
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
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.372       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -7.198E-04   AT NODE       6459   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.198E-04
 LARGEST INCREMENT OF DISP.         3.537E-02   AT NODE        203   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.889E-03   AT NODE      12022   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.307E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.421E-03   TIME AVG. HEAT FLUX   5.910E-04
 LARGEST RESIDUAL HEAT FLUX         2.744E-15   AT NODE       3529   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.393E-02   AT NODE       6883   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.159E-03   AT NODE       3053   DOF 11
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


 AVERAGE FORCE                      0.366       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -7.239E-04   AT NODE       6459   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.239E-04
 LARGEST INCREMENT OF DISP.         3.542E-02   AT NODE        204   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.922E-03   AT NODE       6893   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.347E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.425E-03   TIME AVG. HEAT FLUX   5.939E-04
 LARGEST RESIDUAL HEAT FLUX        -2.743E-15   AT NODE       2101   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.363E-02   AT NODE       6473   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.958E-03   AT NODE       3054   DOF 11
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
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.360       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -7.271E-04   AT NODE       6146   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.271E-04
 LARGEST INCREMENT OF DISP.         3.552E-02   AT NODE       6500   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.908E-03   AT NODE       3135   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.358E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.433E-03   TIME AVG. HEAT FLUX   5.968E-04
 LARGEST RESIDUAL HEAT FLUX        -3.180E-15   AT NODE       3340   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.333E-02   AT NODE       6477   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.954E-03   AT NODE       3055   DOF 11
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
        SOLVER ELAPSED TIME:  360ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.354       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -7.352E-04   AT NODE       6467   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.352E-04
 LARGEST INCREMENT OF DISP.         3.560E-02   AT NODE       6776   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.959E-03   AT NODE       6899   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.386E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.437E-03   TIME AVG. HEAT FLUX   5.998E-04
 LARGEST RESIDUAL HEAT FLUX         2.924E-15   AT NODE       1743   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.316E-02   AT NODE        196   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.880E-03   AT NODE        414   DOF 11
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
        SOLVER ELAPSED TIME:  358ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.348       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -7.380E-04   AT NODE       6465   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.380E-04
 LARGEST INCREMENT OF DISP.         3.572E-02   AT NODE        208   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -4.959E-03   AT NODE      12026   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.381E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.444E-03   TIME AVG. HEAT FLUX   6.027E-04
 LARGEST RESIDUAL HEAT FLUX        -2.603E-15   AT NODE       2812   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.309E-02   AT NODE       6484   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.852E-03   AT NODE       3057   DOF 11
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
        SOLVER ELAPSED TIME:  363ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.342       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -7.443E-04   AT NODE       6472   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.443E-04
 LARGEST INCREMENT OF DISP.         3.581E-02   AT NODE        209   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.006E-03   AT NODE       6901   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.420E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.446E-03   TIME AVG. HEAT FLUX   6.056E-04
 LARGEST RESIDUAL HEAT FLUX         3.165E-15   AT NODE       3075   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.303E-02   AT NODE       5719   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.011E-03   AT NODE        414   DOF 11
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
        SOLVER ELAPSED TIME:  386ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.336       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -7.533E-04   AT NODE       6476   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.533E-04
 LARGEST INCREMENT OF DISP.         3.595E-02   AT NODE       6511   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.016E-03   AT NODE       3063   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.450E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.452E-03   TIME AVG. HEAT FLUX   6.085E-04
 LARGEST RESIDUAL HEAT FLUX         2.736E-15   AT NODE       3329   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.291E-02   AT NODE       6486   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.878E-03   AT NODE       3059   DOF 11
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
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.329       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -7.536E-04   AT NODE       6474   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.536E-04
 LARGEST INCREMENT OF DISP.         3.604E-02   AT NODE       6518   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.068E-03   AT NODE       6940   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.465E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.454E-03   TIME AVG. HEAT FLUX   6.114E-04
 LARGEST RESIDUAL HEAT FLUX        -3.034E-15   AT NODE       4529   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.280E-02   AT NODE        200   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.898E-03   AT NODE        414   DOF 11
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
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.323       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -7.518E-04   AT NODE       6478   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.518E-04
 LARGEST INCREMENT OF DISP.         3.620E-02   AT NODE        212   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.067E-03   AT NODE      12036   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.451E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.459E-03   TIME AVG. HEAT FLUX   6.143E-04
 LARGEST RESIDUAL HEAT FLUX        -2.639E-15   AT NODE        714   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.281E-02   AT NODE       6493   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.746E-03   AT NODE      11209   DOF 11
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
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.317       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -7.537E-04   AT NODE       6483   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.537E-04
 LARGEST INCREMENT OF DISP.         3.630E-02   AT NODE       6520   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.136E-03   AT NODE       6943   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.481E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.460E-03   TIME AVG. HEAT FLUX   6.172E-04
 LARGEST RESIDUAL HEAT FLUX         3.093E-15   AT NODE        888   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.284E-02   AT NODE       6496   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.658E-03   AT NODE       3045   DOF 11
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
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.311       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -7.510E-04   AT NODE       6481   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.510E-04
 LARGEST INCREMENT OF DISP.         3.646E-02   AT NODE       6530   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.139E-03   AT NODE      12042   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.481E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.466E-03   TIME AVG. HEAT FLUX   6.201E-04
 LARGEST RESIDUAL HEAT FLUX        -2.628E-15   AT NODE       4253   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.278E-02   AT NODE       6772   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.750E-03   AT NODE       3136   DOF 11
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
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.305       TIME AVG. FORCE       0.355    
 LARGEST SCALED RESIDUAL FORCE     -7.425E-04   AT NODE       5718   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.425E-04
 LARGEST INCREMENT OF DISP.         3.658E-02   AT NODE        216   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.190E-03   AT NODE       6950   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.471E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.468E-03   TIME AVG. HEAT FLUX   6.230E-04
 LARGEST RESIDUAL HEAT FLUX         3.384E-15   AT NODE        228   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.275E-02   AT NODE        205   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.761E-03   AT NODE       3062   DOF 11
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
        SOLVER ELAPSED TIME:  354ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.298       TIME AVG. FORCE       0.354    
 LARGEST SCALED RESIDUAL FORCE     -7.318E-04   AT NODE       5720   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.318E-04
 LARGEST INCREMENT OF DISP.         3.673E-02   AT NODE       6532   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.193E-03   AT NODE      12048   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    1.453E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.473E-03   TIME AVG. HEAT FLUX   6.258E-04
 LARGEST RESIDUAL HEAT FLUX         2.736E-15   AT NODE       5489   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.284E-02   AT NODE       6776   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.723E-03   AT NODE       3063   DOF 11
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
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.292       TIME AVG. FORCE       0.354    
 LARGEST SCALED RESIDUAL FORCE     -7.192E-04   AT NODE       6491   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.192E-04
 LARGEST INCREMENT OF DISP.         3.685E-02   AT NODE       5709   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.289E-03   AT NODE       6510   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.456E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.473E-03   TIME AVG. HEAT FLUX   6.287E-04
 LARGEST RESIDUAL HEAT FLUX        -2.608E-15   AT NODE       3859   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.288E-02   AT NODE       6505   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.607E-03   AT NODE      12035   DOF 11
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
        SOLVER ELAPSED TIME:  378ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.287       TIME AVG. FORCE       0.354    
 LARGEST SCALED RESIDUAL FORCE     -7.056E-04   AT NODE       6489   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.056E-04
 LARGEST INCREMENT OF DISP.         3.702E-02   AT NODE        220   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.267E-03   AT NODE      12261   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.417E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.475E-03   TIME AVG. HEAT FLUX   6.315E-04
 LARGEST RESIDUAL HEAT FLUX        -2.981E-15   AT NODE       3291   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.282E-02   AT NODE        209   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.635E-03   AT NODE       3066   DOF 11
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
        SOLVER ELAPSED TIME:  380ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.281       TIME AVG. FORCE       0.354    
 LARGEST SCALED RESIDUAL FORCE     -6.824E-04   AT NODE       6494   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.824E-04
 LARGEST INCREMENT OF DISP.         3.713E-02   AT NODE       6964   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.314E-03   AT NODE       6514   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.370E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.473E-03   TIME AVG. HEAT FLUX   6.344E-04
 LARGEST RESIDUAL HEAT FLUX         2.751E-15   AT NODE       4603   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.293E-02   AT NODE       6511   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.595E-03   AT NODE       6530   DOF 11
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
        SOLVER ELAPSED TIME:  410ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.275       TIME AVG. FORCE       0.353    
 LARGEST SCALED RESIDUAL FORCE     -6.525E-04   AT NODE       6497   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.525E-04
 LARGEST INCREMENT OF DISP.         3.731E-02   AT NODE       6546   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.351E-03   AT NODE      12313   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.332E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.474E-03   TIME AVG. HEAT FLUX   6.372E-04
 LARGEST RESIDUAL HEAT FLUX         2.708E-15   AT NODE       2217   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.301E-02   AT NODE       6518   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.572E-03   AT NODE       6526   DOF 11
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
        SOLVER ELAPSED TIME:  395ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.269       TIME AVG. FORCE       0.353    
 LARGEST SCALED RESIDUAL FORCE     -6.189E-04   AT NODE       6504   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.189E-04
 LARGEST INCREMENT OF DISP.         3.744E-02   AT NODE        223   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.391E-03   AT NODE        398   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.262E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.471E-03   TIME AVG. HEAT FLUX   6.399E-04
 LARGEST RESIDUAL HEAT FLUX        -2.590E-15   AT NODE        826   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.299E-02   AT NODE        212   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.634E-03   AT NODE        217   DOF 11
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
        SOLVER ELAPSED TIME:  404ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.263       TIME AVG. FORCE       0.353    
 LARGEST SCALED RESIDUAL FORCE     -5.799E-04   AT NODE       6499   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.799E-04
 LARGEST INCREMENT OF DISP.         3.761E-02   AT NODE       6553   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.391E-03   AT NODE      12319   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.178E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.471E-03   TIME AVG. HEAT FLUX   6.427E-04
 LARGEST RESIDUAL HEAT FLUX        -2.917E-15   AT NODE       4336   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.310E-02   AT NODE       6520   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.631E-03   AT NODE       5709   DOF 11
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
        SOLVER ELAPSED TIME:  376ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.257       TIME AVG. FORCE       0.352    
 LARGEST SCALED RESIDUAL FORCE     -3.493E-04   AT NODE       6518   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -3.493E-04
 LARGEST INCREMENT OF DISP.         3.775E-02   AT NODE        226   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.482E-03   AT NODE       6786   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -7.194E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.468E-03   TIME AVG. HEAT FLUX   6.454E-04
 LARGEST RESIDUAL HEAT FLUX         2.587E-15   AT NODE       3353   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.322E-02   AT NODE       6530   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.706E-03   AT NODE       6541   DOF 11
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
        SOLVER ELAPSED TIME:  380ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.252       TIME AVG. FORCE       0.352    
 LARGEST SCALED RESIDUAL FORCE     -4.778E-04   AT NODE      12051   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -4.778E-04
 LARGEST INCREMENT OF DISP.         3.794E-02   AT NODE        227   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.484E-03   AT NODE      12325   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -9.948E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.467E-03   TIME AVG. HEAT FLUX   6.481E-04
 LARGEST RESIDUAL HEAT FLUX        -3.013E-15   AT NODE       4210   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.324E-02   AT NODE        216   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.682E-03   AT NODE       6964   DOF 11
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
        SOLVER ELAPSED TIME:  366ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.246       TIME AVG. FORCE       0.352    
 LARGEST SCALED RESIDUAL FORCE     -5.606E-04   AT NODE      12052   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -5.606E-04
 LARGEST INCREMENT OF DISP.         3.810E-02   AT NODE       5697   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.526E-03   AT NODE       6920   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.154E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.462E-03   TIME AVG. HEAT FLUX   6.508E-04
 LARGEST RESIDUAL HEAT FLUX         2.543E-15   AT NODE       4466   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.341E-02   AT NODE       6532   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.773E-03   AT NODE       6546   DOF 11
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
        SOLVER ELAPSED TIME:  389ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.240       TIME AVG. FORCE       0.351    
 LARGEST SCALED RESIDUAL FORCE     -6.142E-04   AT NODE      12054   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.142E-04
 LARGEST INCREMENT OF DISP.         3.831E-02   AT NODE        230   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.551E-03   AT NODE      11248   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.271E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.460E-03   TIME AVG. HEAT FLUX   6.534E-04
 LARGEST RESIDUAL HEAT FLUX         2.896E-15   AT NODE       3288   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.353E-02   AT NODE       5709   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.752E-03   AT NODE        223   DOF 11
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
        SOLVER ELAPSED TIME:  400ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.234       TIME AVG. FORCE       0.351    
 LARGEST SCALED RESIDUAL FORCE     -6.652E-04   AT NODE      12056   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -6.652E-04
 LARGEST INCREMENT OF DISP.         3.850E-02   AT NODE       6566   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.595E-03   AT NODE       6547   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.398E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.453E-03   TIME AVG. HEAT FLUX   6.560E-04
 LARGEST RESIDUAL HEAT FLUX         2.279E-15   AT NODE        673   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.362E-02   AT NODE        220   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.853E-03   AT NODE       6553   DOF 11
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
        SOLVER ELAPSED TIME:  389ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.227       TIME AVG. FORCE       0.351    
 LARGEST SCALED RESIDUAL FORCE     -7.116E-04   AT NODE      12059   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.116E-04
 LARGEST INCREMENT OF DISP.         3.874E-02   AT NODE        233   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.625E-03   AT NODE      12328   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.492E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.446E-03   TIME AVG. HEAT FLUX   6.586E-04
 LARGEST RESIDUAL HEAT FLUX        -3.246E-15   AT NODE       1632   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.384E-02   AT NODE       6964   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        5.992E-03   AT NODE        329   DOF 11
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
        SOLVER ELAPSED TIME:  387ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.221       TIME AVG. FORCE       0.350    
 LARGEST SCALED RESIDUAL FORCE     -7.689E-04   AT NODE      12061   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -7.689E-04
 LARGEST INCREMENT OF DISP.         3.896E-02   AT NODE       6792   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.680E-03   AT NODE        331   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.615E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.436E-03   TIME AVG. HEAT FLUX   6.611E-04
 LARGEST RESIDUAL HEAT FLUX        -3.043E-15   AT NODE        402   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.393E-02   AT NODE       6546   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -5.911E-03   AT NODE       6560   DOF 11
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
        SOLVER ELAPSED TIME:  390ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.214       TIME AVG. FORCE       0.350    
 LARGEST SCALED RESIDUAL FORCE     -8.225E-04   AT NODE      12063   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.225E-04
 LARGEST INCREMENT OF DISP.         3.925E-02   AT NODE        235   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.691E-03   AT NODE      12335   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.734E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.428E-03   TIME AVG. HEAT FLUX   6.636E-04
 LARGEST RESIDUAL HEAT FLUX         2.632E-15   AT NODE       3835   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.414E-02   AT NODE       6550   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.129E-03   AT NODE        400   DOF 11
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
        SOLVER ELAPSED TIME:  382ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.208       TIME AVG. FORCE       0.349    
 LARGEST SCALED RESIDUAL FORCE     -8.619E-04   AT NODE      12065   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -8.619E-04
 LARGEST INCREMENT OF DISP.         3.954E-02   AT NODE       6578   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.775E-03   AT NODE        408   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.823E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.417E-03   TIME AVG. HEAT FLUX   6.660E-04
 LARGEST RESIDUAL HEAT FLUX         3.259E-15   AT NODE        831   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.436E-02   AT NODE       6553   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.002E-03   AT NODE        230   DOF 11
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
        SOLVER ELAPSED TIME:  373ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.201       TIME AVG. FORCE       0.349    
 LARGEST SCALED RESIDUAL FORCE     -9.044E-04   AT NODE      12066   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.044E-04
 LARGEST INCREMENT OF DISP.         3.991E-02   AT NODE      12101   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.780E-03   AT NODE       5699   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -1.927E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.406E-03   TIME AVG. HEAT FLUX   6.684E-04
 LARGEST RESIDUAL HEAT FLUX         2.646E-15   AT NODE        750   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.451E-02   AT NODE        226   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.181E-03   AT NODE        403   DOF 11
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
        SOLVER ELAPSED TIME:  380ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.194       TIME AVG. FORCE       0.348    
 LARGEST SCALED RESIDUAL FORCE     -9.542E-04   AT NODE      11280   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.542E-04
 LARGEST INCREMENT OF DISP.         4.029E-02   AT NODE       6598   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.884E-03   AT NODE       6983   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.050E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.394E-03   TIME AVG. HEAT FLUX   6.707E-04
 LARGEST RESIDUAL HEAT FLUX        -2.717E-15   AT NODE        242   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.485E-02   AT NODE       6560   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.126E-03   AT NODE        233   DOF 11
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
        SOLVER ELAPSED TIME:  392ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.187       TIME AVG. FORCE       0.348    
 LARGEST SCALED RESIDUAL FORCE     -9.849E-04   AT NODE      11279   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -9.849E-04
 LARGEST INCREMENT OF DISP.         4.078E-02   AT NODE      12106   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.903E-03   AT NODE      12349   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.112E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.382E-03   TIME AVG. HEAT FLUX   6.730E-04
 LARGEST RESIDUAL HEAT FLUX         2.469E-15   AT NODE        528   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.504E-02   AT NODE        229   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        6.200E-03   AT NODE        331   DOF 11
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
        SOLVER ELAPSED TIME:  397ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.180       TIME AVG. FORCE       0.347    
 LARGEST SCALED RESIDUAL FORCE     -1.014E-03   AT NODE      12070   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.014E-03
 LARGEST INCREMENT OF DISP.         4.131E-02   AT NODE       6600   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.009E-03   AT NODE       6796   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.229E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.366E-03   TIME AVG. HEAT FLUX   6.752E-04
 LARGEST RESIDUAL HEAT FLUX        -2.682E-15   AT NODE       4320   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.543E-02   AT NODE       6570   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.281E-03   AT NODE        235   DOF 11
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
        SOLVER ELAPSED TIME:  399ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.172       TIME AVG. FORCE       0.347    
 LARGEST SCALED RESIDUAL FORCE     -1.036E-03   AT NODE      12071   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.036E-03
 LARGEST INCREMENT OF DISP.         4.198E-02   AT NODE      11770   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.068E-03   AT NODE      12271   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.268E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.350E-03   TIME AVG. HEAT FLUX   6.773E-04
 LARGEST RESIDUAL HEAT FLUX         2.719E-15   AT NODE       1489   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.570E-02   AT NODE        232   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.367E-03   AT NODE       6578   DOF 11
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
        SOLVER ELAPSED TIME:  393ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.165       TIME AVG. FORCE       0.346    
 LARGEST SCALED RESIDUAL FORCE     -1.054E-03   AT NODE      12074   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.054E-03
 LARGEST INCREMENT OF DISP.         4.276E-02   AT NODE      12112   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.177E-03   AT NODE       6797   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.342E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.331E-03   TIME AVG. HEAT FLUX   6.794E-04
 LARGEST RESIDUAL HEAT FLUX        -2.564E-15   AT NODE       1171   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.618E-02   AT NODE       6572   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.492E-03   AT NODE        238   DOF 11
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
        SOLVER ELAPSED TIME:  394ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.157       TIME AVG. FORCE       0.346    
 LARGEST SCALED RESIDUAL FORCE     -1.072E-03   AT NODE      12077   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.072E-03
 LARGEST INCREMENT OF DISP.         4.370E-02   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.281E-03   AT NODE      12358   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.390E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.311E-03   TIME AVG. HEAT FLUX   6.813E-04
 LARGEST RESIDUAL HEAT FLUX         2.769E-15   AT NODE        459   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.657E-02   AT NODE        234   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.610E-03   AT NODE       6598   DOF 11
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
        SOLVER ELAPSED TIME:  547ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.148       TIME AVG. FORCE       0.345    
 LARGEST SCALED RESIDUAL FORCE     -1.088E-03   AT NODE      12081   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.088E-03
 LARGEST INCREMENT OF DISP.         4.484E-02   AT NODE       6611   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.422E-03   AT NODE        405   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.433E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.287E-03   TIME AVG. HEAT FLUX   6.832E-04
 LARGEST RESIDUAL HEAT FLUX        -2.516E-15   AT NODE       3423   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.713E-02   AT NODE       5640   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.780E-03   AT NODE       6591   DOF 11
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
        SOLVER ELAPSED TIME:  384ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.140       TIME AVG. FORCE       0.344    
 LARGEST SCALED RESIDUAL FORCE     -1.104E-03   AT NODE      12083   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.104E-03
 LARGEST INCREMENT OF DISP.         4.627E-02   AT NODE       6613   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.581E-03   AT NODE       6598   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.512E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.260E-03   TIME AVG. HEAT FLUX   6.850E-04
 LARGEST RESIDUAL HEAT FLUX        -2.673E-15   AT NODE       4545   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.765E-02   AT NODE       6585   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -6.903E-03   AT NODE       6600   DOF 11
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
        SOLVER ELAPSED TIME:  393ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.131       TIME AVG. FORCE       0.344    
 LARGEST SCALED RESIDUAL FORCE     -1.119E-03   AT NODE      12089   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.119E-03
 LARGEST INCREMENT OF DISP.         4.803E-02   AT NODE       6824   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.792E-03   AT NODE        241   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.557E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.229E-03   TIME AVG. HEAT FLUX   6.867E-04
 LARGEST RESIDUAL HEAT FLUX        -2.219E-15   AT NODE        343   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.820E-02   AT NODE        239   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.131E-03   AT NODE       6604   DOF 11
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
        SOLVER ELAPSED TIME:  385ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.121       TIME AVG. FORCE       0.343    
 LARGEST SCALED RESIDUAL FORCE     -1.144E-03   AT NODE      11789   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.144E-03
 LARGEST INCREMENT OF DISP.         5.029E-02   AT NODE      11101   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -7.039E-03   AT NODE      12110   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.632E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.194E-03   TIME AVG. HEAT FLUX   6.883E-04
 LARGEST RESIDUAL HEAT FLUX        -2.527E-15   AT NODE       3348   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.887E-02   AT NODE       6587   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.296E-03   AT NODE       6608   DOF 11
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
        SOLVER ELAPSED TIME:  381ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.111       TIME AVG. FORCE       0.342    
 LARGEST SCALED RESIDUAL FORCE     -1.141E-03   AT NODE      12091   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.141E-03
 LARGEST INCREMENT OF DISP.         5.309E-02   AT NODE       6635   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -7.378E-03   AT NODE       6604   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.702E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.155E-03   TIME AVG. HEAT FLUX   6.897E-04
 LARGEST RESIDUAL HEAT FLUX        -2.478E-15   AT NODE       5340   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.944E-02   AT NODE       6600   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -7.501E-03   AT NODE       6820   DOF 11
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
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.940E-02   TIME AVG. FORCE       0.342    
 LARGEST SCALED RESIDUAL FORCE     -1.146E-03   AT NODE      12093   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.146E-03
 LARGEST INCREMENT OF DISP.         5.664E-02   AT NODE      12138   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -7.793E-03   AT NODE       6608   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.776E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.110E-03   TIME AVG. HEAT FLUX   6.910E-04
 LARGEST RESIDUAL HEAT FLUX         2.656E-15   AT NODE       5418   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.986E-02   AT NODE       6604   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.836E-03   AT NODE        405   DOF 11
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
        SOLVER ELAPSED TIME:  379ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.727E-02   TIME AVG. FORCE       0.341    
 LARGEST SCALED RESIDUAL FORCE     -1.144E-03   AT NODE      12097   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.144E-03
 LARGEST INCREMENT OF DISP.         6.085E-02   AT NODE       6842   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -8.357E-03   AT NODE       6820   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.810E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.058E-03   TIME AVG. HEAT FLUX   6.922E-04
 LARGEST RESIDUAL HEAT FLUX         2.035E-15   AT NODE       3399   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -4.004E-02   AT NODE       6608   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.318E-03   AT NODE        339   DOF 11
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
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.394E-02   TIME AVG. FORCE       0.340    
 LARGEST SCALED RESIDUAL FORCE     -1.146E-03   AT NODE      12099   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.146E-03
 LARGEST INCREMENT OF DISP.         6.538E-02   AT NODE        264   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -9.117E-03   AT NODE      12122   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    2.602E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  1.000E-03   TIME AVG. HEAT FLUX   6.931E-04
 LARGEST RESIDUAL HEAT FLUX         2.234E-15   AT NODE       3285   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.968E-02   AT NODE       6820   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        9.230E-03   AT NODE       6593   DOF 11
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
        SOLVER ELAPSED TIME:  357ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.952E-02   TIME AVG. FORCE       0.339    
 LARGEST SCALED RESIDUAL FORCE     -1.148E-03   AT NODE      12102   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.148E-03
 LARGEST INCREMENT OF DISP.         6.905E-02   AT NODE        266   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.018E-02   AT NODE       6824   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.539E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  9.360E-04   TIME AVG. HEAT FLUX   6.939E-04
 LARGEST RESIDUAL HEAT FLUX         2.460E-15   AT NODE       5063   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.837E-02   AT NODE       6621   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.074E-02   AT NODE        341   DOF 11
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
        SOLVER ELAPSED TIME:  364ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.432E-02   TIME AVG. FORCE       0.338    
 LARGEST SCALED RESIDUAL FORCE     -1.149E-03   AT NODE      12106   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.149E-03
 LARGEST INCREMENT OF DISP.         6.972E-02   AT NODE       6699   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.171E-02   AT NODE       6631   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.525E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  8.635E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX        -2.197E-15   AT NODE       5175   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.558E-02   AT NODE       6824   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.337E-02   AT NODE       6601   DOF 11
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
        SOLVER ELAPSED TIME:  361ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.917E-02   TIME AVG. FORCE       0.337    
 LARGEST SCALED RESIDUAL FORCE     -1.183E-03   AT NODE      12110   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.183E-03
 LARGEST INCREMENT OF DISP.         6.453E-02   AT NODE       6658   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.395E-02   AT NODE       5505   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION   -2.660E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  7.818E-04   TIME AVG. HEAT FLUX   6.946E-04
 LARGEST RESIDUAL HEAT FLUX        -2.683E-15   AT NODE       3019   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -3.080E-02   AT NODE       6830   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.550E-02   AT NODE       6607   DOF 11
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
        SOLVER ELAPSED TIME:  369ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.574E-03   TIME AVG. FORCE       0.336    
 LARGEST SCALED RESIDUAL FORCE     -1.225E-03   AT NODE      12112   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.225E-03
 LARGEST INCREMENT OF DISP.         5.205E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.681E-02   AT NODE       6648   DOF  3
   INSTANCE: PLATE-1                                                                         
 ESTIMATE OF DISP.    CORRECTION    3.212E-04
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE HEAT FLUX                  6.943E-04   TIME AVG. HEAT FLUX   6.946E-04
 LARGEST RESIDUAL HEAT FLUX        -2.062E-15   AT NODE       4056   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.        -2.409E-02   AT NODE       6635   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.711E-02   AT NODE       5509   DOF 11
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
        SOLVER ELAPSED TIME:  381ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.330E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.330E-03
 LARGEST INCREMENT OF DISP.         3.490E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.863E-02   AT NODE      12221   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX         2.292E-15   AT NODE       3586   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.627E-02   AT NODE       6616   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  386ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.298E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.298E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.656E-03   AT NODE      12139   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX        -2.456E-15   AT NODE       1377   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        8.796E-15   AT NODE       3688   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  388ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.275E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.275E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -3.602E-03   AT NODE        258   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX         2.779E-15   AT NODE       4270   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.926E-14   AT NODE       3240   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  392ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.259E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.259E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -2.654E-03   AT NODE       6635   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX        -2.495E-15   AT NODE       5132   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        2.467E-14   AT NODE       4211   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  377ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.246E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.246E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -2.107E-03   AT NODE          7   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX         2.455E-15   AT NODE       5132   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.278E-14   AT NODE       4277   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  424ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.236E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.236E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.704E-03   AT NODE       6838   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX         2.499E-15   AT NODE       4061   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.410E-14   AT NODE        620   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  378ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.228E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.228E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.474E-03   AT NODE        256   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX         2.540E-15   AT NODE       5413   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.116E-14   AT NODE       4764   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  382ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.220E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.220E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.311E-03   AT NODE       6631   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX        -2.664E-15   AT NODE      14907   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        7.931E-15   AT NODE       1853   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  392ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.212E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.212E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.098E-03   AT NODE        255   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX        -2.299E-15   AT NODE       3509   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.757E-14   AT NODE       3769   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  383ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    10


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.205E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.205E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -1.017E-03   AT NODE        255   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX         2.222E-15   AT NODE       5306   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.554E-14   AT NODE       4087   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  374ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    11


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.198E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.198E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -9.024E-04   AT NODE        255   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX        -2.543E-15   AT NODE       3399   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.294E-14   AT NODE       4093   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  378ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    12


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.191E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.191E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -8.206E-04   AT NODE       6830   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX         2.484E-15   AT NODE       3399   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.919E-14   AT NODE       6706   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  370ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    13


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.184E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.184E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -7.773E-04   AT NODE       6830   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX        -2.440E-15   AT NODE        682   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -1.191E-14   AT NODE       2912   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  373ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    14


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.177E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.177E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.716E-04   AT NODE       6830   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX         2.391E-15   AT NODE       4888   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.       -2.308E-14   AT NODE       1660   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  368ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    15


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.169E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.169E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -6.100E-04   AT NODE        254   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX         2.562E-15   AT NODE       3443   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.877E-14   AT NODE       9544   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 
	UNSYMMETRIC NEW DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST x 1 MPI RANK PER HOST x 2 THREADS PER RANK
        NUMBER OF EQUATIONS:  65612
        NUMBER OF RHS:        1
        NUMBER OF FLOPS:      9.693e+09
        SOLVER ELAPSED TIME:  362ms

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    16


 AVERAGE FORCE                      3.214E-03   TIME AVG. FORCE       0.335    
 LARGEST SCALED RESIDUAL FORCE     -1.162E-03   AT NODE      12116   DOF  3
   INSTANCE: PLATE-1                                                                         
  CORRESPONDING RESIDUAL FORCE     -1.162E-03
 LARGEST INCREMENT OF DISP.         3.484E-02   AT NODE          9   DOF  3
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO DISP.       -5.685E-04   AT NODE        254   DOF  3
   INSTANCE: PLATE-1                                                                         
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE HEAT FLUX                  6.068E-04   TIME AVG. HEAT FLUX   6.944E-04
 LARGEST RESIDUAL HEAT FLUX        -2.845E-15   AT NODE       3828   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST INCREMENT OF TEMP.         2.032E-02   AT NODE        410   DOF 11
   INSTANCE: PLATE-1                                                                         
 LARGEST CORRECTION TO TEMP.        1.381E-14   AT NODE       6885   DOF 11
   INSTANCE: PLATE-1                                                                         
          THE HEAT FLUX EQUILIBRIUM EQUATIONS HAVE CONVERGED
 

 ***NOTE: THE SOLUTION FAILS TO CONVERGE IN THE MAXIMUM NUMBER OF EQUILIBRIUM 
          ITERATIONS ALLOWED.
 

 ***ERROR: FIXED TIME INCREMENT IS TOO LARGE

 ***ERROR: THE ANALYSIS HAS BEEN TERMINATED DUE TO PREVIOUS ERRORS. ALL OUTPUT 
           REQUESTS HAVE BEEN WRITTEN FOR THE LAST CONVERGED INCREMENT.



     ANALYSIS SUMMARY:
     TOTAL OF        332  INCREMENTS
                       1  CUTBACKS IN AUTOMATIC INCREMENTATION
                     388  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                     388  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                     313  INVOLVE MATRIX DECOMPOSITION, INCLUDING
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
       USER TIME (SEC)      =     3.59E+02
       SYSTEM TIME (SEC)    =      20.    
       TOTAL CPU TIME (SEC) =     3.79E+02
       WALLCLOCK TIME (SEC) =          319
       MEMORY PEAK (GB)     =            0
