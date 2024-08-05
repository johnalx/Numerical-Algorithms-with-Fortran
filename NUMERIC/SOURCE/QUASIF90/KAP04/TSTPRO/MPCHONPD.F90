      PROGRAM TEST 
!                                                ( Thomas Meuser )      
!***********************************************************************
!                                                                       
!     Testprogram for subroutine CHONPD.                                
!     Solve a linear system  A*X = RS  for a a symmetric strongly non-  
!     singular matrix A using the Cholesky method.                      
!                                                                       
!     The supplied input data produces the following results:           
!                                                                       
![  EXAMPLE:                                                           ]
![  ========                                                           ]
![  COEFFICIENT MATRIX A:                                              ]
![                                                                     ]
![   .50000D+01  -.10000D+01  -.10000D+01  -.10000D+01                 ]
![  -.10000D+01   .50000D+01  -.10000D+01  -.10000D+01                 ]
![  -.10000D+01  -.10000D+01   .50000D+01  -.10000D+01                 ]
![  -.10000D+01  -.10000D+01  -.10000D+01   .50000D+01                 ]
![                                                                     ]
![  RIGHT HAND SIDE                                                    ]
![                                                                     ]
![   .20000D+01   .20000D+01   .20000D+01   .20000D+01                 ]
![                                                                     ]
![  SOLUTION                                                           ]
![                                                                     ]
![   .10000D+01   .10000D+01   .10000D+01   .10000D+01                 ]
![  STOP. NO ERROR!                                                    ]
!                                                                       
!     Further tests possible for modified data.                         
!                                                                       
!***********************************************************************
!                                                                       
      PARAMETER (N = 4, IA = N) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A (IA, N), X (N), RS (N), Z (N) 
!                                                                       
!   Initialize; put in different data if desired                        
!                                                                       
      DATA ( (A (J, I), I = 1, N), J = 1, N) / 5.D0, - 1.D0, - 1.D0,    &
      - 1.D0, - 1.D0, 5.D0, - 1.D0, - 1.D0, - 1.D0, - 1.D0, 5.D0,       &
      - 1.D0, - 1.D0, - 1.D0, - 1.D0, 5.D0 /                            
      DATA (RS (I), I = 1, N) / 2.D0, 2.D0, 2.D0, 2.D0 / 
!                                                                       
!   Output of test example                                              
!                                                                       
      WRITE ( *, 2010) 
      DO 100 I = 1, N 
  100 WRITE ( *, 2020) (A (I, J), J = 1, N) 
      WRITE ( *, 2030) (RS (I), I = 1, N) 
!                                                                       
      CALL CHONPD (N, A, IA, RS, X, Z, MARKE) 
!                                                                       
      WRITE ( *, 2040) (X (I), I = 1, N) 
!                                                                       
!   Output of error code  MARKE (-1/0/1)                                
!                                                                       
      MARKE = MARKE+2 
      GOTO (10, 20, 30), MARKE 
   10 WRITE ( * , 2050) 'WARNING: MATRIX A IS NOT POSITIVE DEFINITE!' 
      STOP 
   20 WRITE ( * , 2050) 'ERROR: MATRIX A IS NOT STRONGLY NONSINGULAR!' 
      STOP 
   30 WRITE ( * , 2050) 'STOP. NO ERROR!' 
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
