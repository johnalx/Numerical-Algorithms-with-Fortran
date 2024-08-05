!                                                                       
!***********************************************************************
!                                                                       
!     Testprogram for the subroutine  ADSOR.                            
!     Solve an inhomogeneous system of linear equations  A*X = B with   
!     det A  not zero iteratively via the SOR method.                   
!                                                                       
!     The output for the test example is as follows:                    
!                                                                       
![  EXAMPLE:                                                            
![  ========                                                            
![  COEFFICIENT MATRIX A:                                               
![                                                                      
![   .40000D+01  -.10000D+01  -.10000D+01   .00000D+00                  
![  -.10000D+01   .40000D+01   .00000D+00  -.10000D+01                  
![  -.10000D+01   .00000D+00   .40000D+01  -.10000D+01                  
![   .00000D+00  -.10000D+01  -.10000D+01   .40000D+01                  
![                                                                      
![  RIGHT HAND SIDE:                                                    
![                                                                      
![   .00000D+00   .00000D+00   .10000D+04   .10000D+04                  
![                                                                      
![  STARTING VECTOR X:                                                  
![                                                                      
![   .00000D+00   .00000D+00   .00000D+00   .00000D+00                  
![                                                                      
![                                                                      
![                                                                      
![                                                                      
![  INPUT PARAMETERS:                                                   
![    KADAPT =    4                                                     
![    EPS =  .10000000000000D-08                                        
![    KMAX =  100                                                       
![    IMETH  =    0                                                     
![    ISWITC =    0                                                     
![    OMEGA =  .10000000000000D+01                                      
![                                                                      
![                                                                      
![  OUTPUT PARAMETERS:                                                  
![    OMEGA =  .10739144783362D+01                                      
![    ITNUMB =   13                                                     
![    IERR =    1                                                       
![                                                                      
![                                                                      
![  SOLUTION VECTOR X:                                                  
![   .12500D+03   .12500D+03   .37500D+03   .37500D+03                  
![                                                                      
![  RESIDUAL VECTOR RES:                                                
![   .79888D-08  -.30867D-09  -.30867D-09  -.53348D-09                  
![                                                                      
![                                                                      
![  INPUT PARAMETERS:                                                   
![    KADAPT =    4                                                     
![    EPS =  .10000000000000D-08                                        
![    KMAX =  100                                                       
![    IMETH  =    1                                                     
![    ISWITC =    0                                                     
![    OMEGA =  .10739144783362D+01                                      
![                                                                      
![                                                                      
![  OUTPUT PARAMETERS:                                                  
![    OMEGA =  .10739144783362D+01                                      
![    ITNUMB =   10                                                     
![    IERR =    1                                                       
![                                                                      
![                                                                      
![  SOLUTION VECTOR X:                                                  
![   .12500D+03   .12500D+03   .37500D+03   .37500D+03                  
![                                                                      
![  RESIDUAL VECTOR RES:                                                
![   .15031D-07  -.56438D-08  -.44269D-08   .10584D-09                  
![                                                                      
![                                                                      
![  INPUT PARAMETERS:                                                   
![    KADAPT =    4                                                     
![    EPS =  .10000000000000D-08                                        
![    KMAX =  100                                                       
![    IMETH  =    2                                                     
![    ISWITC =    0                                                     
![    OMEGA =  .10739144783362D+01                                      
![                                                                      
![                                                                      
![  OUTPUT PARAMETERS:                                                  
![    OMEGA =  .10000000000000D+01                                      
![    ITNUMB =   17                                                     
![    IERR =    1                                                       
![                                                                      
![                                                                      
![  SOLUTION VECTOR X:                                                  
![   .12500D+03   .12500D+03   .37500D+03   .37500D+03                  
![                                                                      
![  RESIDUAL VECTOR RES:                                                
![   .43656D-07   .10914D-07   .10914D-07   .00000D+00                  
![                                                                      
![                                                                      
![  INPUT PARAMETERS:                                                   
![    KADAPT =    4                                                     
![    EPS =  .10000000000000D-08                                        
![    KMAX =  100                                                       
![    IMETH  =    0                                                     
![    ISWITC =    1                                                     
![    OMEGA =  .10000000000000D+01                                      
![                                                                      
![                                                                      
![  OUTPUT PARAMETERS:                                                  
![    OMEGA =  .10739144783362D+01                                      
![    ITNUMB =   13                                                     
![    IERR =    1                                                       
![                                                                      
![                                                                      
![  SOLUTION VECTOR X:                                                  
![   .12500D+03   .12500D+03   .37500D+03   .37500D+03                  
![                                                                      
![  RESIDUAL VECTOR RES:                                                
![   .79888D-08  -.30867D-09  -.30867D-09  -.53348D-09                  
![                                                                      
![                                                                      
![  INPUT PARAMETERS:                                                   
![    KADAPT =    4                                                     
![    EPS =  .10000000000000D-08                                        
![    KMAX =  100                                                       
![    IMETH  =    1                                                     
![    ISWITC =    1                                                     
![    OMEGA =  .10739144783362D+01                                      
![                                                                      
![                                                                      
![  OUTPUT PARAMETERS:                                                  
![    OMEGA =  .10739144783362D+01                                      
![    ITNUMB =   10                                                     
![    IERR =    1                                                       
![                                                                      
![                                                                      
![  SOLUTION VECTOR X:                                                  
![   .12500D+03   .12500D+03   .37500D+03   .37500D+03                  
![                                                                      
![  RESIDUAL VECTOR RES:                                                
![   .15031D-07  -.56438D-08  -.44269D-08   .10584D-09                  
![                                                                      
![                                                                      
![  INPUT PARAMETERS:                                                   
![    KADAPT =    4                                                     
![    EPS =  .10000000000000D-08                                        
![    KMAX =  100                                                       
![    IMETH  =    2                                                     
![    ISWITC =    1                                                     
![    OMEGA =  .10739144783362D+01                                      
![                                                                      
![                                                                      
![  OUTPUT PARAMETERS:                                                  
![    OMEGA =  .10000000000000D+01                                      
![    ITNUMB =   17                                                     
![    IERR =    1                                                       
![                                                                      
![                                                                      
![  SOLUTION VECTOR X:                                                  
![   .12500D+03   .12500D+03   .37500D+03   .37500D+03                  
![                                                                      
![  RESIDUAL VECTOR RES:                                                
![   .43656D-07   .10914D-07   .10914D-07   .00000D+00                  
![                                                                      
![                                                                      
![  INPUT PARAMETERS:                                                   
![    KADAPT =    4                                                     
![    EPS =  .10000000000000D-08                                        
![    KMAX =  100                                                       
![    IMETH  =    0                                                     
![    ISWITC =    2                                                     
![    OMEGA =  .10000000000000D+01                                      
![                                                                      
![                                                                      
![  OUTPUT PARAMETERS:                                                  
![    OMEGA =  .10739144783362D+01                                      
![    ITNUMB =   13                                                     
![    IERR =    1                                                       
![                                                                      
![                                                                      
![  SOLUTION VECTOR X:                                                  
![   .12500D+03   .12500D+03   .37500D+03   .37500D+03                  
![                                                                      
![  RESIDUAL VECTOR RES:                                                
![   .79888D-08  -.30867D-09  -.30867D-09  -.53348D-09                  
![                                                                      
![                                                                      
![  INPUT PARAMETERS:                                                   
![    KADAPT =    4                                                     
![    EPS =  .10000000000000D-08                                        
![    KMAX =  100                                                       
![    IMETH  =    1                                                     
![    ISWITC =    2                                                     
![    OMEGA =  .10739144783362D+01                                      
![                                                                      
![                                                                      
![  OUTPUT PARAMETERS:                                                  
![    OMEGA =  .10739144783362D+01                                      
![    ITNUMB =   10                                                     
![    IERR =    1                                                       
![                                                                      
![                                                                      
![  SOLUTION VECTOR X:                                                  
![   .12500D+03   .12500D+03   .37500D+03   .37500D+03                  
![                                                                      
      PROGRAM TEST 
![  RESIDUAL VECTOR RES:                                                
![   .15031D-07  -.56438D-08  -.44269D-08   .10584D-09                  
![                                                                      
![                                                                      
![  INPUT PARAMETERS:                                                   
![    KADAPT =    4                                                     
![    EPS =  .10000000000000D-08                                        
![    KMAX =  100                                                       
![    IMETH  =    2                                                     
![    ISWITC =    2                                                     
![    OMEGA =  .10739144783362D+01                                      
![                                                                      
![                                                                      
![  OUTPUT PARAMETERS:                                                  
![    OMEGA =  .10000000000000D+01                                      
![    ITNUMB =   17                                                     
![    IERR =    1                                                       
![                                                                      
![                                                                      
![  SOLUTION VECTOR X:                                                  
![   .12500D+03   .12500D+03   .37500D+03   .37500D+03                  
![                                                                      
![  RESIDUAL VECTOR RES:                                                
![   .43656D-07   .10914D-07   .10914D-07   .00000D+00                  
![                                                                      
![                                                                      
![  INPUT PARAMETERS:                                                   
![    KADAPT =    4                                                     
![    EPS =  .10000000000000D-08                                        
![    KMAX =  100                                                       
![    IMETH  =    0                                                     
![    ISWITC =    3                                                     
![    OMEGA =  .10000000000000D+01                                      
![                                                                      
![                                                                      
![  OUTPUT PARAMETERS:                                                  
![    OMEGA =  .10739144783362D+01                                      
![    ITNUMB =   13                                                     
![    IERR =    1                                                       
![                                                                      
![                                                                      
![  SOLUTION VECTOR X:                                                  
![   .12500D+03   .12500D+03   .37500D+03   .37500D+03                  
![                                                                      
![  RESIDUAL VECTOR RES:                                                
![   .79888D-08  -.30867D-09  -.30867D-09  -.53348D-09                  
![                                                                      
![                                                                      
![  INPUT PARAMETERS:                                                   
![    KADAPT =    4                                                     
![    EPS =  .10000000000000D-08                                        
![    KMAX =  100                                                       
![    IMETH  =    1                                                     
![    ISWITC =    3                                                     
![    OMEGA =  .10739144783362D+01                                      
![                                                                      
![                                                                      
![  OUTPUT PARAMETERS:                                                  
![    OMEGA =  .10739144783362D+01                                      
![    ITNUMB =   10                                                     
![    IERR =    1                                                       
![                                                                      
![                                                                      
![  SOLUTION VECTOR X:                                                  
![   .12500D+03   .12500D+03   .37500D+03   .37500D+03                  
![                                                                      
![  RESIDUAL VECTOR RES:                                                
![   .15031D-07  -.56438D-08  -.44269D-08   .10584D-09                  
![                                                                      
![                                                                      
![  INPUT PARAMETERS:                                                   
![    KADAPT =    4                                                     
![    EPS =  .10000000000000D-08                                        
![    KMAX =  100                                                       
![    IMETH  =    2                                                     
![    ISWITC =    3                                                     
![    OMEGA =  .10739144783362D+01                                      
![                                                                      
![                                                                      
![  OUTPUT PARAMETERS:                                                  
![    OMEGA =  .10000000000000D+01                                      
![    ITNUMB =   17                                                     
![    IERR =    1                                                       
![                                                                      
![                                                                      
![  SOLUTION VECTOR X:                                                  
![   .12500D+03   .12500D+03   .37500D+03   .37500D+03                  
![                                                                      
![  RESIDUAL VECTOR RES:                                                
![   .43656D-07   .10914D-07   .10914D-07   .00000D+00                  
!                                                                       
!***********************************************************************
!                                                                       
!                                                                       
      DOUBLEPRECISION A (1:4, 1:4), B (1:4), X (1:4), RES (1:4),        &
      WORK (1:4, 1:3), EPS, OMEGA                                       
      OMEGA = 1.0D+00 
!                                                                       
! output of the test example                                            
!                                                                       
      CALL INIT (A, B, X, KADAPT, EPS, KMAX) 
      WRITE ( *, 2000) 
      DO 10 I = 1, 4 
         WRITE ( *, 2010) (A (I, J), J = 1, 4) 
   10 END DO 
      WRITE ( *, 2020) (B (I), I = 1, 4) 
      WRITE ( *, 2030) (X (I), I = 1, 4) 
!                                                                       
      DO 100 ISWITC = 0, 3 
         DO 200 IMETH = 0, 2 
            CALL INIT (A, B, X, KADAPT, EPS, KMAX) 
!                                                                       
            WRITE ( *, 3000) KADAPT, EPS, KMAX, IMETH, ISWITC, OMEGA 
!                                                                       
            CALL ADSOR (A, 4, 4, B, X, KADAPT, EPS, KMAX, IMETH, ISWITC,&
            OMEGA, WORK, RES, ITNUMB, IERR)                             
!                                                                       
            WRITE ( *, 4000) OMEGA, ITNUMB, IERR 
!                                                                       
            WRITE ( *, 2040) (X (I), I = 1, 4) 
            WRITE ( *, 2050) (RES (I), I = 1, 4) 
  200    END DO 
  100 END DO 
      STOP 
 2000 FORMAT (1X,'C[',2X,'EXAMPLE:',T78,']*',/,                         &
     &        1X,'C[',2X,8('='),T78,']*',/,                             &
     &        1X,'C[',2X,'COEFFICIENT MATRIX A:',T78,']*',/,            &
     &        1X,'C[',T78,']*')                                         
 2010 FORMAT (1X,'C[',4(1X,D12.5),T78,']*') 
 2020 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'RIGHT HAND SIDE:',T78,']*',/,                 &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',4(1X,D12.5),T78,']*',/,                           &
     &        1X,'C[',T78,']*')                                         
 2030 FORMAT (1X,'C[',2X,'STARTING VECTOR X:',T78,']*',/,               &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',4(1X,D12.5),T78,']*',/,                           &
     &        1X,'C[',T78,']*',/,1X,'C[',T78,']*')                      
 2040 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'SOLUTION VECTOR X:',T78,']*',/,               &
     &        1X,'C[',4(1X,D12.5),T78,']*')                             
 2050 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'RESIDUAL VECTOR RES:',T78,']*',/,             &
     &        1X,'C[',4(1X,D12.5),T78,']*')                             
 3000 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'INPUT PARAMETERS:',T78,']*',/,                &
     &        1X,'C[',4X,'KADAPT = ',I4,T78,']*',/,                     &
     &        1X,'C[',4X,'EPS = ',D20.14,T78,']*',/,                    &
     &        1X,'C[',4X,'KMAX = ',I4,T78,']*',/,                       &
     &        1X,'C[',4X,'IMETH  = ',I4,T78,']*',/,                     &
     &        1X,'C[',4X,'ISWITC = ',I4,T78,']*',/,                     &
     &        1X,'C[',4X,'OMEGA = ',D20.14,T78,']*',/,                  &
     &        1X,'C[',T78,']*')                                         
 4000 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'OUTPUT PARAMETERS:',T78,']*',/,               &
     &        1X,'C[',4X,'OMEGA = ',D20.14,T78,']*',/,                  &
     &        1X,'C[',4X,'ITNUMB = ',I4,T78,']*',/,                     &
     &        1X,'C[',4X,'IERR = ',I4,T78,']*',/,                       &
     &        1X,'C[',T78,']*' )                                        
      END PROGRAM TEST                              
!                                                                       
!                                                                       
      SUBROUTINE INIT (A, B, X, KADAPT, EPS, KMAX) 
!                                                                       
!  Initialize system matrix A and right hand side B;                    
!  supply starting vector  X for the iteration                          
!                                                                       
      DOUBLEPRECISION A (1:4, 1:4), B (1:4), X (1:4), EPS 
      KADAPT = 4 
      EPS = 1.0D-09 
      KMAX = 100 
      DO 10 I = 1, 4 
         DO 20 J = 1, 4 
            IF (I.EQ.J) THEN 
               A (I, J) = 4.0D+00 
            ELSEIF ( (I + J) .EQ.5) THEN 
               A (I, J) = 0.0D+00 
            ELSE 
               A (I, J) = - 1.0D+00 
            ENDIF 
   20    END DO 
         IF (I.GT.2) THEN 
            B (I) = 1.0D+03 
         ELSE 
            B (I) = 0.0D+00 
         ENDIF 
         X (I) = 0.0D+00 
   10 END DO 
      RETURN 
      END SUBROUTINE INIT                           
