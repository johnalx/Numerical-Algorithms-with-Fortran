      PROGRAM TEST 
!                                        (Thomas Meuser) 88/07/13 , 94/0
!***********************************************************************
!                                                                       
!     Test program for the subroutine BVP                               
!     Solves a boundary value problem for a first degree system of ordin
!     differential equations of the form:                               
!                                                                       
!                   Y1' = F1(X,Y1,Y2,...,YN)                            
!                   Y2' = F2(X,Y1,Y2,...,YN)                            
!                             ...                                       
!                   YN' = FN(X,Y1,Y2,...,YN)                            
!                                                                       
!     It uses the shooting method and gives an approximate value for an 
!     approximation to the starting value  Y(A) for the system.         
!                                                                       
!     The output for the given test data should be:                     
!                                                                       
![                                                                      
![                                                                      
![                                                                      
![ TEST EXAMPLE (WITH IVP PROGRAM NO 1):                                
![ =====================================                                
![ GIVEN SET OF DIFFERENTIAL EQUATIONS:                                 
![ ------------------------------------                                 
![          Y'( 1) = Y(2)                                               
![          Y'( 2) = -Y(1)**3                                           
![ WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)                            
![   AT POINT  A =   .000D+00  :  .000D+00                              
![   AT POINT  B =   .100D+01  :  .000D+00                              
![                                                                      
![ AT  A =   .000D+00 THE FOLLOWING FUNCTION VALUES ARE PROVIDED:       
![                                                                      
![   Y( 1) =   .000D+00                                                 
![   Y( 2) =   .120D+02                                                 
![                                                                      
![ REQUIRED PARAMETER:                                                  
![ -------------------                                                  
![ -PRECISION FOR THE IVP =     .1000000000D-06                         
![ -PRECISION FOR THE BVP =     .1000000000D-07                         
![ -STEP WIDTH FOR THE IVP=     .1000000000D-01                         
![ -MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE IVP'S =  10000     
![ -MAX. NUMBER OF NEWTON-ITERATION STEPS                  =   1000     
![                                                                      
![ SOLUTION OF THE PROBLEM:                                             
![ ------------------------                                             
![ THE STARTING VALUE OF A SOLUTION OF THE BVP IS APPROXIMATED          
![ AT A =  .000D+00 AS FOLLOWS:                                         
![                                                                      
![   Y( 1) =      .0000000000D+00                                       
![   Y( 2) =      .9722981025D+01                                       
![                                                                      
![ THIS REQUIRES     3 NEWTON-ITERATIONS!                               
![                                                                      
![ DETERMINED AS DESIRED!                                               
![                                                                      
![                                                                      
![                                                                      
![ TEST EXAMPLE (WITH IVP PROGRAM NO 2):                                
![ =====================================                                
![ GIVEN SET OF DIFFERENTIAL EQUATIONS:                                 
![ ------------------------------------                                 
![          Y'( 1) = Y(2)                                               
![          Y'( 2) = -Y(1)**3                                           
![ WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)                            
![   AT POINT  A =   .000D+00  :  .000D+00                              
![   AT POINT  B =   .100D+01  :  .000D+00                              
![                                                                      
![ AT  A =   .000D+00 THE FOLLOWING FUNCTION VALUES ARE PROVIDED:       
![                                                                      
![   Y( 1) =   .000D+00                                                 
![   Y( 2) =   .120D+02                                                 
![                                                                      
![ REQUIRED PARAMETER:                                                  
![ -------------------                                                  
![ -PRECISION FOR THE IVP =     .1000000000D-06                         
![ -PRECISION FOR THE BVP =     .1000000000D-07                         
![ -STEP WIDTH FOR THE IVP=     .1000000000D-01                         
![ -MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE IVP'S =  10000     
![ -MAX. NUMBER OF NEWTON-ITERATION STEPS                  =   1000     
![                                                                      
![ SOLUTION OF THE PROBLEM:                                             
![ ------------------------                                             
![ THE STARTING VALUE OF A SOLUTION OF THE BVP IS APPROXIMATED          
![ AT A =  .000D+00 AS FOLLOWS:                                         
![                                                                      
![   Y( 1) =      .0000000000D+00                                       
![   Y( 2) =      .9771778775D+01                                       
![                                                                      
![ THIS REQUIRES     3 NEWTON-ITERATIONS!                               
![                                                                      
![ DETERMINED AS DESIRED!                                               
![                                                                      
![                                                                      
![                                                                      
![ TEST EXAMPLE (WITH IVP PROGRAM NO 3):                                
![ =====================================                                
![ GIVEN SET OF DIFFERENTIAL EQUATIONS:                                 
![ ------------------------------------                                 
![          Y'( 1) = Y(2)                                               
![          Y'( 2) = -Y(1)**3                                           
![ WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)                            
![   AT POINT  A =   .000D+00  :  .000D+00                              
![   AT POINT  B =   .100D+01  :  .000D+00                              
![                                                                      
![ AT  A =   .000D+00 THE FOLLOWING FUNCTION VALUES ARE PROVIDED:       
![                                                                      
![   Y( 1) =   .000D+00                                                 
![   Y( 2) =   .120D+02                                                 
![                                                                      
![ REQUIRED PARAMETER:                                                  
![ -------------------                                                  
![ -PRECISION FOR THE IVP =     .1000000000D-06                         
![ -PRECISION FOR THE BVP =     .1000000000D-07                         
![ -STEP WIDTH FOR THE IVP=     .1000000000D-01                         
![ -MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE IVP'S =  10000     
![ -MAX. NUMBER OF NEWTON-ITERATION STEPS                  =   1000     
![                                                                      
![ SOLUTION OF THE PROBLEM:                                             
![ ------------------------                                             
![ THE STARTING VALUE OF A SOLUTION OF THE BVP IS APPROXIMATED          
![ AT A =  .000D+00 AS FOLLOWS:                                         
![                                                                      
![   Y( 1) =      .0000000000D+00                                       
![   Y( 2) =      .9722980349D+01                                       
![                                                                      
![ THIS REQUIRES     4 NEWTON-ITERATIONS!                               
![                                                                      
![ DETERMINED AS DESIRED!                                               
![                                                                      
![                                                                      
![                                                                      
![ TEST EXAMPLE (WITH IVP PROGRAM NO 4):                                
![ =====================================                                
![ GIVEN SET OF DIFFERENTIAL EQUATIONS:                                 
![ ------------------------------------                                 
![          Y'( 1) = Y(2)                                               
![          Y'( 2) = -Y(1)**3                                           
![ WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)                            
![   AT POINT  A =   .000D+00  :  .000D+00                              
![   AT POINT  B =   .100D+01  :  .000D+00                              
![                                                                      
![ AT  A =   .000D+00 THE FOLLOWING FUNCTION VALUES ARE PROVIDED:       
![                                                                      
![   Y( 1) =   .000D+00                                                 
![   Y( 2) =   .120D+02                                                 
![                                                                      
![ REQUIRED PARAMETER:                                                  
![ -------------------                                                  
![ -PRECISION FOR THE IVP =     .1000000000D-06                         
![ -PRECISION FOR THE BVP =     .1000000000D-07                         
![ -STEP WIDTH FOR THE IVP=     .1000000000D-01                         
![ -MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE IVP'S =  10000     
![ -MAX. NUMBER OF NEWTON-ITERATION STEPS                  =   1000     
![                                                                      
![ SOLUTION OF THE PROBLEM:                                             
![ ------------------------                                             
![ THE STARTING VALUE OF A SOLUTION OF THE BVP IS APPROXIMATED          
![ AT A =  .000D+00 AS FOLLOWS:                                         
![                                                                      
![   Y( 1) =      .0000000000D+00                                       
![   Y( 2) =      .9327367735D+01                                       
![                                                                      
![ THIS REQUIRES     5 NEWTON-ITERATIONS!                               
![                                                                      
![ DETERMINED AS DESIRED!                                               
!***********************************************************************
!                                                                       
!      (    Depending on the starting value for  Y(2)   ( Y(1) = 0 is   
!      (    given . )  in each of the examples one of the following     
!      (    values for  Y(2) should be reached approximately at  A :    
!      (  0 ; +- 9.7229810... ; +- 38.8919241... ; +- 87.5068292... ; ..
!                                                                       
!         Further test runs can be made by modifying the DATA values    
!         and the EXTERNAL SUBROUTINES that define the DEs.             
!                                                                       
!***********************************************************************
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (N = 2) 
      DIMENSION YANF (N), AUXF (N, 12), WORK1 (5 * N + 78) 
      DIMENSION WORK21 (N, N), WORK22 (N + 9, 9), WORK23 (2 * N + 10,   &
      10)                                                               
      CHARACTER(80) YSTRI (N) 
      EXTERNAL DGL, RAND 
!                                                                       
!     Preassign test values. Exchange test examples while adjusting     
!     the  EXTERNAL SUBROUTINES and the parameter N.                    
!                                                                       
      DATA A, B, H, WA, WB / 0.D0, 1.D0, .01D0, 0.D0, 0.D0 / 
      DATA IFMAX, ITMAX / 10000, 1000 / 
      EPSAWP = 1.D-7 
      EPSRB = 1.D-8 
!                                                                       
!     Set up output                                                     
!                                                                       
      DATA (YSTRI (I) , I = 1, N) / 'Y(2)', '-Y(1)**3' / 
!                                                                       
!     Set up test example output for various examples                   
!                                                                       
      DO 100 IV = 1, 4 
!                                                                       
         YANF (1) = 0.0D0 
         YANF (2) = 12.0D0 
         WRITE ( * , 2060) ' ' 
         WRITE ( *, 2000) IV 
         WRITE ( *, 2010) (I, YSTRI (I), I = 1, N) 
         WRITE ( *, 2015) A, WA, B, WB 
         WRITE ( *, 2020) A, (I, YANF (I), I = 1, N) 
         WRITE ( *, 2030) EPSAWP, EPSRB, H, IFMAX, ITMAX 
!                                                                       
!                                                                       
         CALL BVP (A, B, H, YANF, N, DGL, RAND, IV, EPSAWP, EPSRB,      &
         IFMAX, ITMAX, ITER, IFEHL, AUXF, N, WORK1, WORK21, WORK22,     &
         WORK23)                                                        
!                                                                       
!     Put out 'break-off' notice according to IFEHL (0/1/2/3/4/5/6/7/8) 
!                                                                       
         IF (IFEHL.EQ.1) THEN 
            WRITE ( * , 2060) 'ERRORBOUND(S) TOO SMALL' 
         ELSEIF (IFEHL.EQ.2) THEN 
            WRITE ( * , 2060) 'ERROR: B <= A' 
         ELSEIF (IFEHL.EQ.3) THEN 
            WRITE ( * , 2060) 'ERROR: STEP SIZE H <= 0' 
         ELSEIF (IFEHL.EQ.4) THEN 
            WRITE ( * , 2060) 'ERROR: NUMBER N OF DEQ''S NOT CORRECT' 
         ELSEIF (IFEHL.EQ.5) THEN 
            WRITE ( * , 2060) 'ERROR: WRONG IVP-PROGRAM CHOICE' 
         ELSEIF (IFEHL.EQ.6) THEN 
      WRITE ( * , 2060) 'IFMAX TO SMALL FOR SOLUTION OF THE IVP-PROGRAM' 
         ELSEIF (IFEHL.EQ.7) THEN 
      WRITE ( * , 2060) 'AFTER ITMAX STEPS PRECISION WAS NOT REACHED!' 
         ELSEIF (IFEHL.EQ.8) THEN 
            WRITE ( * , 2060) 'ERROR: JACOBI-MATRIX IS SINGULAR!' 
         ELSE 
!                                                                       
!     Put out solution                                                  
!                                                                       
            WRITE ( *, 2040) A, (I, YANF (I), I = 1, N) 
            WRITE ( *, 2050) ITER 
            WRITE ( * , 2060) 'DETERMINED AS DESIRED!' 
         ENDIF 
  100 END DO 
      STOP 
!                                                                       
!                                                                       
 2000 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',1X,'TEST EXAMPLE (WITH IVP PROGRAM NO ',I1,'):',  &
     &        T78,']*',/,1X,'C[',1X,37('='),T78,']*',/,                 &
     &        1X,'C[',1X,'GIVEN SET OF DIFFERENTIAL EQUATIONS:',        &
     &        T78,']*',/,1X,'C[',1X,36('-'),T78,']*')                   
 2010 FORMAT ((1X,'C[',10X,'Y''(',I2,') = ',A30,T78,']*')) 
 2015 FORMAT (1X,'C[',1X,'WITH BOUNDARY CONDITIONS : VALUE OF  Y(1)',   &
     &        T78,']*',/,1X,'C[',3X,'AT POINT  A = ',D10.3,'  :',D10.3, &
     &        T78,']*',/,1X,'C[',3X,'AT POINT  B = ',D10.3,'  :',D10.3, &
     &        T78,']*')                                                 
 2020 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',1X,'AT  A = ',D10.3,                              &
     &        ' THE FOLLOWING FUNCTION VALUES ARE PROVIDED:',T78,']*'   &
     &        ,/,1X,'C[',T78,']*',/,                                    &
     &        (1X,'C[',3X,'Y(',I2,') = ',D10.3,T78,']*'))               
 2030 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',1X,'REQUIRED PARAMETER:',T78,']*',/,              &
     &        1X,'C[',1X,19('-'),T78,']*',/,                            &
     &        1X,'C[',1X,'-PRECISION FOR THE IVP =',D20.10,T78,']*',/,  &
     &        1X,'C[',1X,'-PRECISION FOR THE BVP =',D20.10,T78,']*',/,  &
     &        1X,'C[',1X,'-STEP WIDTH FOR THE IVP=',D20.10,T78,']*',/,  &
     &        1X,'C[',1X,                                               &
     &        '-MAX. NUMBER OF F-EVALUATIONS FOR SOLUTION OF THE ',     &
     &        'IVP''S = ',I6,T78,']*',/,                                &
     &        1X,'C[',1X,                                               &
     &        '-MAX. NUMBER OF NEWTON-ITERATION STEPS               ',  &
     &        '   = ',I6,T78,']*',/,1X,'C[',T78,']*')                   
 2040 FORMAT (1X,'C[',1X,'SOLUTION OF THE PROBLEM:',T78,']*',/,         &
     &        1X,'C[',1X,24('-'),T78,']*',/,                            &
     &        1X,'C[',1X,'THE STARTING VALUE OF A SOLUTION OF THE BVP', &
     &        ' IS APPROXIMATED',T78,']*',/,                            &
     &        1X,'C[',1X,'AT A =',D10.3,' AS FOLLOWS:',T78,']*',/,      &
     &        1X,'C[',T78,']*',/,                                       &
     &        (1X,'C[',3X,'Y(',I2,') = ',D20.10,T78,']*'))              
 2050 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',1X,'THIS REQUIRES ',I5,' NEWTON-ITERATIONS!',     &
     &        T78,']*')                                                 
 2060 FORMAT (1X,'C[',T78,']*',/,1X,'C[',1X,A,T78,']*') 
      END PROGRAM TEST                              
!                                                                       
!                                                                       
      SUBROUTINE DGL (X, Y, N, F) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION Y (N), F (N) 
      F (1) = Y (2) 
      F (2) = - Y (1) **3 
      RETURN 
      END SUBROUTINE DGL                            
!                                                                       
!                                                                       
      SUBROUTINE RAND (YA, YB, N, R) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION YA (N), YB (N), R (N) 
      R (1) = YA (1) 
      R (2) = YB (1) 
      RETURN 
      END SUBROUTINE RAND                           
