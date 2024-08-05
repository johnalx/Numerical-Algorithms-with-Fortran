      PROGRAM TEST 
!                                                                       
!*********************************************************************  
!                                                                    *  
!  Testprogram for the  SUBROUTINE FRKFSY                            *  
!                                                                    *  
!  Test example:                                                     *  
!  We solve one initial value problem using the method of            *  
!  RUNGE-KUTTA-FEHLBERG.                                             *  
!                                                                    *  
!  Test problem:                                                     *  
!               Y1'=Y2                        Y1(0)=1   I=[0,20]     *  
!               Y2'=-Y1/(SQRT(Y1^2+Y3^2))^3   Y2(0)=0                *  
!               Y3'=Y4                        Y3(0)=0                *  
!               Y4'=-Y3/(SQRT(Y1^2+Y3^2))^3   Y4(0)=1                *  
!                                                                    *  
!                                                                    *  
!  Results:                                                          *  
!                                                                    *  
![                                                                  ]*  
![                                                                  ]*  
![ ABSERR= .00000000E+00    RELERR= .10000000E-05                   ]*  
![                                                                  ]*  
![                                                                  ]*  
![ IFEHL=  1     X=  20.000                                         ]*  
![                                                                  ]*  
![          SOLUTION:I APPROXIMATION        I EXACT                 ]*  
![          ---------+----------------------+---------------------- ]*  
![          Y1(X)    I  .40807628764144E+00 I  .40808206181339E+00  ]*  
![          Y2(X)    I -.91294791088434E+00 I -.91294525072763E+00  ]*  
![          Y3(X)    I  .91294748136748E+00 I  .91294525072763E+00  ]*  
![          Y4(X)    I  .40807650152524E+00 I  .40808206181339E+00  ]*  
!                                                                    *  
!*********************************************************************  
!                                                                       
! declarations                                                          
!                                                                       
      EXTERNAL DGL 
      DOUBLEPRECISION Y (4), YEX (4) 
      DOUBLEPRECISION X, BETA, H, HMX, ABSERR, RELERR 
!                                                                       
! initialize                                                            
!                                                                       
      ABSERR = 0.0D+00 
      RELERR = 1.0D-06 
      WRITE ( *, 2000) ABSERR, RELERR 
      N = 4 
      X = 0.0D+00 
      H = 0.01D+00 
      HMX = 0.1D+00 
      BETA = 20.D+00 
      Y (1) = 1.0D+00 
      Y (2) = 0.0D+00 
      Y (3) = 0.0D+00 
      Y (4) = 1.0D+00 
!                                                                       
! begin computations                                                    
!                                                                       
      CALL FRKFSY (X, BETA, N, Y, DGL, H, HMX, ABSERR, RELERR, IFEHL) 
      CALL EXAKT (X, N, YEX) 
      WRITE ( *, 2100) IFEHL, X 
!                                                                       
! put out results                                                       
!                                                                       
      DO 40 I = 1, N 
         WRITE ( *, 2200) I, Y (I), YEX (I) 
   40 END DO 
      STOP 
!                                                                       
! Format declarations                                                   
!                                                                       
 2000 FORMAT(1X,'C[',T70,']*',/,1X,'C[',T70,']*',/,                     &
     &       1X,'C[',1X,'ABSERR=',E14.8,'    RELERR=',E14.8,T70,']*',/, &
     &       1X,'C[',T70,']*')                                          
 2100 FORMAT(1X,'C[',T70,']*',/,                                        &
     &       1X,'C[',1X,'IFEHL=',I3,5X,'X=',F8.3,T70,']*',/,            &
     &       1X,'C[',T70,']*',/,                                        &
     &       1X,'C[',10X,'SOLUTION:I APPROXIMATION',8X,'I EXACT',       &
     &          T70,']*',/,                                             &
     &       1X,'C[',10X,9('-'),'+',22('-'),'+',22('-'),T70,']*')       
 2200 FORMAT(1X,'C[',10X,'Y',I1,'(X)',4X,'I ',E20.14,' I ',E20.14,      &
     &          T70,']*')                                               
      END PROGRAM TEST                              
!                                                                       
! Test problem                                                          
!                                                                       
      SUBROUTINE DGL (X, Y, F) 
      DOUBLEPRECISION Y (4), F (4), X 
      F (1) = Y (2) 
      F (2) = - Y (1) / (SQRT (Y (1) * Y (1) + Y (3) * Y (3) ) ) **     &
      3.D+00                                                            
      F (3) = Y (4) 
      F (4) = - Y (3) / (SQRT (Y (1) * Y (1) + Y (3) * Y (3) ) ) **     &
      3.D+00                                                            
      RETURN 
      END SUBROUTINE DGL                            
!                                                                       
      SUBROUTINE EXAKT (X, N, YEXAKT) 
      DOUBLEPRECISION X, YEXAKT (N) 
      YEXAKT (1) = COS (X) 
      YEXAKT (2) = - SIN (X) 
      YEXAKT (3) = SIN (X) 
      YEXAKT (4) = COS (X) 
      RETURN 
      END SUBROUTINE EXAKT                          
