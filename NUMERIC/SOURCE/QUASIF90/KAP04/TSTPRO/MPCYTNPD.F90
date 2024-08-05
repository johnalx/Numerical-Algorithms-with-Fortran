      PROGRAM TEST 
!                                            ( Thomas Meuser )          
!***********************************************************************
!                                                                       
!     Testprogram for subroutine CYTNPD.                                
!     Solve a linear system  A*X = RS  with a symmetric strongly        
!     nonsingular  cyclic tridiagonal matrix A.                         
!                                                                       
!     The test example produces the following output:                   
!                                                                       
![  EXAMPLE:                                                           ]
![  ========                                                           ]
![  COEFFICIENT MATRIX A:                                              ]
![                                                                     ]
![   .50000D+01  -.10000D+01   .00000D+00  -.10000D+01                 ]
![  -.10000D+01   .50000D+01  -.10000D+01   .00000D+00                 ]
![   .00000D+00  -.10000D+01   .50000D+01  -.10000D+01                 ]
![  -.10000D+01   .00000D+00  -.10000D+01   .50000D+01                 ]
![                                                                     ]
![  RIGHT HAND SIDE                                                    ]
![                                                                     ]
![   .30000D+01   .30000D+01   .30000D+01   .30000D+01                 ]
![                                                                     ]
![  SOLUTION                                                           ]
![                                                                     ]
![   .10000D+01   .10000D+01   .10000D+01   .10000D+01                 ]
![  STOP. NO ERROR!                                                    ]
!                                                                       
!     Further tests with differing data are possible.                   
!                                                                       
!***********************************************************************
!                                                                       
      PARAMETER (N = 4, M = N * N) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION HD (N), OD (N), RSP (N), RS (N), X (N), A (N, N) 
!                                                                       
!   Initialize input; change data for other test problems               
!                                                                       
      DATA (HD (I), I = 1, N) / 5.D0, 5.D0, 5.D0, 5.D0 / 
      DATA (OD (I), I = 1, N) / - 1.D0, - 1.D0, - 1.D0, - 1.D0 / 
      DATA (RS (I), I = 1, N) / 3.D0, 3.D0, 3.D0, 3.D0 / 
!                                                                       
!   Output of test example in matrix form                               
!                                                                       
      DATA ( (A (J, I), I = 1, N), J = 1, N) / M * 0.0D0 / 
      A (1, 1) = HD (1) 
      A (1, 2) = OD (1) 
      A (1, N) = OD (N) 
      DO 110 I = 2, N - 1 
         A (I, I) = HD (I) 
         A (I, I + 1) = OD (I) 
         A (I, I - 1) = OD (I - 1) 
  110 END DO 
      A (N, 1) = OD (N) 
      A (N, N) = HD (N) 
      A (N, N - 1) = OD (N - 1) 
      WRITE ( *, 2010) 
      DO 120 I = 1, N 
  120 WRITE ( *, 2020) (A (I, J), J = 1, N) 
      WRITE ( *, 2030) (RS (I), I = 1, N) 
!                                                                       
      CALL CYTNPD (N, HD, OD, RSP, RS, X, MARKE) 
!                                                                       
      WRITE ( *, 2040) (X (I), I = 1, N) 
!                                                                       
!   output of error code  MARKE (-2/-1/0/1)                             
!                                                                       
      MARKE = MARKE+3 
      GOTO (10, 20, 30, 40), MARKE 
   10 WRITE ( * , 2050) 'ERROR:  N <= 2 !' 
      STOP 
   20 WRITE ( * , 2050) 'WARNING: MATRIX A IS NOT POSITIVE DEFINITE!' 
      STOP 
   30 WRITE ( * , 2050) 'ERROR: MATRIX A IS NOT STRONGLY NONSINGULAR!' 
      STOP 
   40 WRITE ( * , 2050) 'STOP. NO ERROR!' 
      STOP 
!                                                                       
 2010 FORMAT (1X,'C[',2X,'EXAMPLE:',T73,']*',/,                         &
     &        1X,'C[',2X,8('='),T73,']*',/,                             &
     &        1X,'C[',2X,'COEFFICIENT MATRIX A:',T73,']*',/,            &
     &        1X,'C[',T73,']*')                                         
 2020 FORMAT (1X,'C[',4(1X,D12.5),T73,']*') 
 2030 FORMAT (1X,'C[',T73,']*',/,                                       &
     &        1X,'C[',2X,'RIGHT HAND SIDE',T73,']*',/,                  &
     &        1X,'C[',T73,']*',/,                                       &
     &        1X,'C[',4(1X,D12.5),T73,']*')                             
 2040 FORMAT (1X,'C[',T73,']*',/,                                       &
     &        1X,'C[',2X,'SOLUTION',T73,']*',/,                         &
     &        1X,'C[',T73,']*',/,                                       &
     &        1X,'C[',4(1X,D12.5),T73,']*')                             
 2050 FORMAT (1X,'C[',2X,A,T73,']*') 
      END PROGRAM TEST                              
