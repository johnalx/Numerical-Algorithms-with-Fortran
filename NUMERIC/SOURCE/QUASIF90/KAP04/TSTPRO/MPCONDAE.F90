      PROGRAM TEST 
!                                                                       
!***********************************************************************
!                                                                       
!     Testprogram for the condition number subroutine CONDAE.           
!                                                                       
!     The test data produces the following output:                      
!                                                                       
![  EXAMPLE:                                                            
![  ========                                                            
![  MATRIX AO:                                                          
![                                                                      
![   .10000D+01   .50000D+00   .33333D+00   .25000D+00   .20000D+00     
![   .50000D+00   .33333D+00   .25000D+00   .20000D+00   .16667D+00     
![   .33333D+00   .25000D+00   .20000D+00   .16667D+00   .14286D+00     
![   .25000D+00   .20000D+00   .16667D+00   .14286D+00   .12500D+00     
![   .20000D+00   .16667D+00   .14286D+00   .12500D+00   .11111D+00     
![                                                                      
![  RIGHT HAND SIDE Y:                                                  
![                                                                      
![   .22833D+01   .14500D+01   .10929D+01   .88452D+00   .74563D+00     
![                                                                      
![                                                                      
![  MATRIX A:                                                           
![                                                                      
![   .10000D+01   .50000D+00   .33333D+00   .25000D+00   .20000D+00     
![   .20000D+00   .66667D-01   .76190D-01   .75000D-01   .71111D-01     
![   .50000D+00   .12500D+01  -.11905D-01  -.18750D-01  -.22222D-01     
![   .33333D+00   .12500D+01   .53333D+00  -.41667D-03  -.84656D-03     
![   .25000D+00   .11250D+01   .20000D+00   .64286D+00  -.11338D-04     
![                                                                      
![  SOLUTION VECTOR X:                                                  
![                                                                      
![   .10000D+01   .10000D+01   .10000D+01   .10000D+01   .10000D+01     
![                                                                      
![                                                                      
![  ESTIMATE OF THE CONDITION NUMBER OF A=-.40544090000000D+32 ;   N=  5
!                                                                       
!***********************************************************************
!                                                                       
      PARAMETER (N = 5, IA = N) 
      DOUBLEPRECISION A (IA, N), AO (IA, N), Y (N), X (N), D (N),       &
      Z (N), R (N)                                                      
      INTEGER IPIVOT (N) 
!                                                                       
      DO 10 I = 1, N 
         DO 10 K = 1, N 
            A (I, K) = 1.0D+00 / DBLE (I + K - 1) 
            AO (I, K) = A (I, K) 
   10 CONTINUE 
      DO 20 I = 1, N 
         Y (I) = 0.0D+00 
         DO 20 K = 1, N 
            Y (I) = Y (I) + A (I, K) 
   20 CONTINUE 
!                                                                       
! Output matrix AO                                                      
!                                                                       
      WRITE ( *, 2000) 
      DO 30 I = 1, 5 
         WRITE ( *, 2020) (AO (I, J), J = 1, 5) 
   30 END DO 
      WRITE ( *, 2030) (Y (I), I = 1, 5) 
!                                                                       
      CALL GAUSS (N, A, IA, Y, X, MARKE, D, IPIVOT) 
!                                                                       
! Output after SUBROUTINE GAUSS                                         
!                                                                       
      IF (MARKE.EQ.0) THEN 
         WRITE ( * , 2050) 'MATRIX SINGULAR', ' N=', N 
         STOP 
      ELSE 
         WRITE ( *, 2010) 
         DO 40 I = 1, 5 
            WRITE ( *, 2020) (A (I, J), J = 1, 5) 
   40    END DO 
         WRITE ( *, 2040) (X (I), I = 1, 5) 
      ENDIF 
!                                                                       
      CALL CONDAE (N, AO, A, IA, IPIVOT, Y, X, Z, R, CONDA) 
!                                                                       
      WRITE ( * , 2060) 'ESTIMATE OF THE CONDITION NUMBER OF A=', CONDA,&
     & ' ;   N=', N                                                     
      STOP 
 2000 FORMAT (1X,'C[',2X,'EXAMPLE:',T78,']*',/,                         &
     &        1X,'C[',2X,8('='),T78,']*',/,                             &
     &        1X,'C[',2X,'MATRIX A0:',T78,']*',/,                       &
     &        1X,'C[',T78,']*')                                         
 2010 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'MATRIX A:',T78,']*',/,                        &
     &        1X,'C[',T78,']*')                                         
 2020 FORMAT (1X,'C[',5(1X,D12.5),T78,']*') 
 2030 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'RIGHT HAND SIDE Y:',T78,']*',/,               &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',5(1X,D12.5),T78,']*',/,                           &
     &        1X,'C[',T78,']*')                                         
 2040 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'SOLUTION VECTOR X:',T78,']*',/,               &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',5(1X,D12.5),T78,']*',/,                           &
     &        1X,'C[',T78,']*',/,1X,'C[',T78,']*')                      
 2050 FORMAT (1X,'C[',2X,A,I3,T78,']*') 
 2060 FORMAT (1X,'C[',2X,A,D20.14,A,I3,T78,']*') 
      END PROGRAM TEST                              
