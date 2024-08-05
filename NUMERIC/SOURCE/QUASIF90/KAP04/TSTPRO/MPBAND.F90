      PROGRAM TEST 
!                                                ( Thomas Meuser )      
!***********************************************************************
!                                                                       
!     Testprogram for the subroutines BAND, BANDP, BANDS.               
!     Solve a system of linear equations  A*X = B with a banded matrix A
!                                                                       
!     For the included test examples the following results are computed:
!                                                                       
![  EXAMPLE:                                                            
![  ========                                                            
![  COEFFICIENT MATRIX A:                                               
![                                                                      
![   .50000D+01   .10000D+01   .10000D+01   .00000D+00   .00000D+00     
![   .20000D+01   .60000D+01   .20000D+01   .00000D+00   .00000D+00     
![   .20000D+01   .20000D+01   .70000D+01   .00000D+00   .20000D+01     
![   .00000D+00   .10000D+01   .30000D+01   .80000D+01   .20000D+01     
![   .00000D+00   .00000D+00   .30000D+01   .10000D+01   .90000D+01     
![                                                                      
![  RIGHT HAND SIDE                                                     
![                                                                      
![   .70000D+01   .10000D+02   .13000D+02   .14000D+02   .13000D+02     
![                                                                      
![                                                                      
![  MATRIX IN PACKED FORM:                                              
![                                                                      
![                                                                      
![   .00000D+00   .00000D+00   .50000D+01   .10000D+01   .10000D+01     
![   .00000D+00   .20000D+01   .60000D+01   .20000D+01   .00000D+00     
![   .20000D+01   .20000D+01   .70000D+01   .00000D+00   .20000D+01     
![   .10000D+01   .30000D+01   .80000D+01   .20000D+01   .00000D+00     
![   .30000D+01   .10000D+01   .90000D+01   .00000D+00   .00000D+00     
![                                                                      
![  SOLUTION                                                            
![                                                                      
![   .10000D+01   .10000D+01   .10000D+01   .10000D+01   .10000D+01     
![  STOP. NO ERROR !                                                    
!                                                                       
!     Further tests can be performed by changing the data.              
!                                                                       
!***********************************************************************
!                                                                       
      PARAMETER (N = 5, LDAP = N) 
      PARAMETER (MU = 2, MO = 2, MB = MU + MO + 1 + MU + 1) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION A (N, N), AP (LDAP, MB) 
      DOUBLEPRECISION B (N) 
      INTEGER IP (N) 
!                                                                       
!   Initialize; for other test examples, change the input data.         
!                                                                       
      DATA ( (A (J, I), I = 1, N), J = 1, N) / 5.D0, 1.D0, 1.D0, 0.D0,  &
      0.D0, 2.D0, 6.D0, 2.D0, 0.D0, 0.D0, 2.D0, 2.D0, 7.D0, 0.D0, 2.D0, &
      0.D0, 1.D0, 3.D0, 8.D0, 2.D0, 0.D0, 0.D0, 3.D0, 1.D0, 9.D0 /      
      DATA (B (I), I = 1, N) / 7.D0, 10.D0, 13.D0, 14.D0, 13.D0 / 
!                                                                       
!   Output of test example                                              
!                                                                       
      WRITE ( *, 2000) 
      DO 10 I = 1, N 
   10 WRITE ( *, 2010) (A (I, J), J = 1, N) 
      WRITE ( *, 2020) (B (I), I = 1, N) 
!                                                                       
!   Compute the packed matrix AP to use in BAND.                        
!                                                                       
      DO 20 I = 1, N 
         DO 20 K = MAX (1, I - MU), MIN (N, I + MO) 
            AP (I, MU + 1 + K - I) = A (I, K) 
   20 CONTINUE 
!                                                                       
!   Output of the condensed matrix AP                                   
!                                                                       
      WRITE ( *, 2023) 
      DO 30 I = 1, N 
         WRITE ( *, 2025) (0.D0, J = 1, MAX (0, MU + 1 - I) ), (AP (I,  &
         J), J = MAX (0, MU + 1 - I) + 1, MIN (MU + 1 + MO, N + MU + 1 -&
         I) ), (0.D0, J = MIN (MU + 1 + MO, N + MU + 1 - I) + 1, MU + 1 &
         + MO)                                                          
   30 END DO 
!                                                                       
      CALL BAND (AP, LDAP, MB, N, MU, MO, B, IFEHL, IP) 
!                                                                       
      IF (IFEHL.EQ.0) THEN 
!                                                                       
!   Output of a solution or an error message according to IFEHL         
!                                                                       
         WRITE ( *, 2030) (B (I), I = 1, N) 
         WRITE ( * , 2040) 'STOP. NO ERROR !' 
         STOP 
      ELSE 
         WRITE ( * , 2040) 'STOP. ERROR --- MATRIX IS SINGULAR !' 
         STOP 
      ENDIF 
!                                                                       
 2000 FORMAT (1X,'C[',2X,'EXAMPLE:',T78,']*',/,                         &
     &        1X,'C[',2X,8('='),T78,']*',/,                             &
     &        1X,'C[',2X,'COEFFICIENT MATRIX A:',T78,']*',/,            &
     &        1X,'C[',T78,']*')                                         
 2010 FORMAT (1X,'C[',5(1X,D12.5),T78,']*') 
 2020 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'RIGHT HAND SIDE',T78,']*',/,                  &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',5(1X,D12.5),T78,']*',/,                           &
     &        1X,'C[',T78,']*',/,1X,'C[',T78,']*')                      
 2023 FORMAT (1X,'C[',2X,'MATRIX IN PACKED FORM:',T78,']*',/,           &
     &        1X,'C[',T78,']*',/,1X,'C[',T78,']*')                      
 2025 FORMAT (1X,'C[',5(1X,D12.5),T78,']*') 
 2030 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'SOLUTION',T78,']*',/,                         &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',5(1X,D12.5),T78,']*')                             
 2040 FORMAT (1X,'C[',2X,A,T78,']*') 
      END PROGRAM TEST                              
