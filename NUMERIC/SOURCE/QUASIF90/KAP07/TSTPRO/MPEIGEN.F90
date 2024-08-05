      PROGRAM TEST 
!                                                                     18
!                                              (Thomas Meuser)          
!***********************************************************************
!                                                                       
!     Testprogram for the subroutine EIGEN.                             
!     Compute all eigenvalues and eigenvectors of a matrix A with the QR
!     method according to Martin, Parlett, Peters, Reinsch and Wilkinson
!                                                                       
!     The test gives the following results:                             
!                                                                       
![  TEST EXAMPLE:                                                       
![  =============                                                       
![  GIVEN MATRIX:                                                       
![                                                                      
![    -.20000D+01   .20000D+01   .20000D+01   .20000D+01                
![    -.30000D+01   .30000D+01   .20000D+01   .20000D+01                
![    -.20000D+01   .00000D+00   .40000D+01   .20000D+01                
![    -.10000D+01   .00000D+00   .00000D+00   .50000D+01                
![                                                                      
![  COMPUTED EIGENVALUES:                                               
![                                                                      
![            REAL PART   I IMAGINARY PART                              
![          ------------------------------                              
![           .10000D+01   I     .00000D+00                              
![           .20000D+01   I     .00000D+00                              
![           .30000D+01   I     .00000D+00                              
![           .40000D+01   I     .00000D+00                              
![                                                                      
![  NORMALIZED EIGENVECTORS:                                            
![                                                                      
![     .10000D+01   .10000D+01   .10000D+01   .10000D+01                
![     .75000D+00   .10000D+01   .10000D+01   .10000D+01                
![     .50000D+00   .66667D+00   .10000D+01   .10000D+01                
![     .25000D+00   .33333D+00   .50000D+00   .10000D+01                
![  STOP. NO ERROR !                                                    
!                                                                       
!     Other tests with different data are possible.                     
!                                                                       
!***********************************************************************
!                                                                       
      PARAMETER (N = 4, ND = N) 
      INTEGER CNT (N), EIGEN, BASIS, RES, LOW, HIGH 
      DOUBLEPRECISION MAT (ND, N), WERTR (N), WERTI (N), EIVEC (ND, N), &
      SKAL (N), D (N)                                                   
      LOGICAL VEC, ORTHO, EVNORM 
!                                                                       
!     Initialize; change if desired                                     
!                                                                       
      DATA MAT / - 2.0D+00, - 3.0D+00, - 2.0D+00, - 1.0D+00, 2.0D+00,   &
      3.0D+00, 0.0D+00, 0.0D+00, 2.0D+00, 2.0D+00, 4.0D+00, 0.0D+00,    &
      2.0D+00, 2.0D+00, 2.0D+00, 5.0D+00 /                              
      DATA BASIS / 2 / 
      DATA VEC, ORTHO, EVNORM / .TRUE., .FALSE., .TRUE. / 
!                                                                       
!     Print out test example                                            
!                                                                       
      WRITE ( *, 2000) 
      DO 10 I = 1, N 
   10 WRITE ( *, 2010) (MAT (I, J), J = 1, N) 
!                                                                       
!                                                                       
      RES = EIGEN (VEC, ORTHO, EVNORM, BASIS, ND, N, MAT, SKAL, D,      &
      EIVEC, WERTR, WERTI, CNT, LOW, HIGH)                              
!                                                                       
!                                                                       
      IF (RES.GT.0) THEN 
         RES = RES - 400 
!                                                                       
!     Output of error code in  RES (0/401/402/403)                      
!                                                                       
         GOTO (30, 40, 50), RES 
   30 WRITE ( * , 2050) 'ERROR: SIZE N OF THE MATRIX IS SMALLER THAN 1 !&
     &'                                                                 
         STOP 
   40    WRITE ( * , 2050) 'ERROR: MAT IS THE ZERO MATRIX !' 
         STOP 
   50 WRITE ( * , 2050) 'ERROR: MAX NUMBER OF STEPS FOR QR-METHOD IS REA&
     &CHED!'                                                            
         STOP 
      ELSE 
!                                                                       
!     Output                                                            
!                                                                       
         WRITE ( *, 2015) 
         WRITE ( *, 2020) (WERTR (I), WERTI (I), I = 1, N) 
         WRITE ( *, 2030) 
         DO 20 I = 1, N 
   20    WRITE ( *, 2040) (EIVEC (I, J), J = 1, N) 
         WRITE ( * , 2050) 'STOP. NO ERROR !' 
         STOP 
      ENDIF 
!                                                                       
!                                                                       
 2000 FORMAT(1X,'C[',2X,'TEST EXAMPLE:',T78,']*',/,                     &
     &       1X,'C[',2X,13('='),T78,']*',/,                             &
     &       1X,'C[',2X,'GIVEN MATRIX:',T78,']*',/,                     &
     &       1X,'C[',T78,']*')                                          
 2010 FORMAT(1X,'C[',2X,4(1X,D12.5),T78,']*') 
 2015 FORMAT(1X,'C[',T78,']*',/,                                        &
     &       1X,'C[',2X,'COMPUTED EIGENVALUES:',T78,']*',/,             &
     &       1X,'C[',T78,']*',/,                                        &
     &       1X,'C[',12X,'REAL PART   I ','IMAGINARY PART',T78,']*',/,  &
     &       1X,'C[',10X,30('-'),T78,']*')                              
 2020 FORMAT(1X,'C[',7X,D14.5,'   I ',D14.5,T78,']*') 
 2030 FORMAT(1X,'C[',T78,']*',/,                                        &
     &       1X,'C[',2X,'NORMALIZED EIGENVECTORS:',T78,']*',/,          &
     &       1X,'C[',T78,']*')                                          
 2040 FORMAT(1X,'C[',2X,4(1X,D12.5),T78,']*') 
 2050 FORMAT (1X,'C[',2X,A,T78,']*') 
      END PROGRAM TEST                              
