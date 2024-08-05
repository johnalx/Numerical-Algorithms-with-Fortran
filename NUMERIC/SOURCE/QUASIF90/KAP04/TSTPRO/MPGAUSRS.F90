      PROGRAM TEST 
!                                                                       
!***********************************************************************
!                                                                       
!  Testprogram for the subroutines GAUSRS, GAUSSP, GAUSSS.              
!  Compute the inverse of a matrix A.                                   
!                                                                       
![  EXAMPLE:                                                           ]
![  ========                                                           ]
![  COEFFICIENT MATRIX A:                                              ]
![                                                                     ]
![   .30000D+01   .10000D+01   .20000D+01   .50000D+01                 ]
![   .20000D+01   .10000D+01   .30000D+01   .70000D+01                 ]
![   .30000D+01   .10000D+01   .20000D+01   .40000D+01                 ]
![   .40000D+01   .10000D+01   .30000D+01   .20000D+01                 ]
![                                                                     ]
![  N RIGHT HAND SIDES                                                 ]
![                                                                     ]
![   .10000D+01   .00000D+00   .00000D+00   .00000D+00                 ]
![   .00000D+00   .10000D+01   .00000D+00   .00000D+00                 ]
![   .00000D+00   .00000D+00   .10000D+01   .00000D+00                 ]
![   .00000D+00   .00000D+00   .00000D+00   .10000D+01                 ]
![                                                                     ]
![  MATRIX INVERSE                                                     ]
![                                                                     ]
![   .25000D+01  -.50000D+00  -.25000D+01   .50000D+00                 ]
![  -.10500D+02   .50000D+00   .13500D+02  -.25000D+01                 ]
![  -.50000D+00   .50000D+00  -.50000D+00   .50000D+00                 ]
![   .10000D+01   .00000D+00  -.10000D+01   .00000D+00                 ]
![  STOP. NO ERROR!                                                    ]
!                                                                       
!***********************************************************************
!                                                                       
      PARAMETER (N = 4) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A (1:N, 1:N), RS (1:N, 1:N), XL (1:N, 1:N), D (1:N) 
      INTEGER IPIVOT (1:N) 
      DATA A / 3.D0, 2.D0, 3.D0, 4.D0, 4 * 1.D0, 2.D0, 3.D0, 2.D0, 3.D0,&
      5.D0, 7.D0, 4.D0, 2.D0 /                                          
      DATA RS / 1.D0, 4 * 0.D0, 1.D0, 4 * 0.D0, 1.D0, 4 * 0.D0, 1.D0 / 
      WRITE ( *, 2000) 
      DO 10 I = 1, N 
         WRITE ( *, 2100) (A (I, K), K = 1, N) 
   10 END DO 
      WRITE ( *, 2200) 
      DO 20 I = 1, N 
         WRITE ( *, 2100) (RS (I, K), K = 1, N) 
   20 END DO 
      CALL GAUSRS (N, A, N, N, RS, XL, MARKE, D, IPIVOT) 
      WRITE ( *, 2300) 
      IF (MARKE.NE.0) THEN 
         DO 30 I = 1, N 
            WRITE ( *, 2100) (XL (I, K), K = 1, N) 
   30    END DO 
         WRITE ( * , 2400) 'STOP. NO ERROR!' 
      ELSE 
         WRITE ( * , 2400) 'STOP. ERROR --- MATRIX IS SINGULAR !' 
      ENDIF 
      STOP 
 2000 FORMAT (1X,'C[',2X,'EXAMPLE:',T73,']*',/,                         &
     &        1X,'C[',2X,8('='),T73,']*',/,                             &
     &        1X,'C[',2X,'COEFFICIENT MATRIX A:',T73,']*',/,            &
     &        1X,'C[',T73,']*')                                         
 2100 FORMAT (1X,'C[',4(1X,D12.5),T73,']*') 
 2200 FORMAT (1X,'C[',T73,']*',/,                                       &
     &        1X,'C[',2X,'N RIGHT HAND SIDES',T73,']*',/,               &
     &        1X,'C[',T73,']*')                                         
 2300 FORMAT (1X,'C[',T73,']*',/,                                       &
     &        1X,'C[',2X,'MATRIX INVERSE',T73,']*',/,                   &
     &        1X,'C[',T73,']*')                                         
 2400 FORMAT (1X,'C[',2X,A,T73,']*') 
      END PROGRAM TEST                              
