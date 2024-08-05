      PROGRAM TEST 
!                                                                       
!***********************************************************************
!                                                                       
!   Test program for the subroutines BICSP3 and FIBIC2      03/29/89    
!                                                                       
!     ( Elmar Pohl )                                                    
!                                                                       
!                                                                       
!                                                                       
!     Example for the use of  BICSP3 and FIBIC2:                        
!     We compute a bicubic spline surface with given normal vector and  
!     integrate the spline function on the rectangle R = [XU,XO] x [YU,Y
!                                                                       
!     (The output assumes fÅr M <= 10 )                                 
!                                                                       
![                                                                      
![ COMPUTE A BICUBIC SPLINE USING  BICSP3                               
![ ======================================                               
![                                                                      
![ X-COORDINATES OF THE NODES                                           
![ --------------------------                                           
![                                                                      
![    I     X(I)                                                        
![ -----------------                                                    
![    0    0.000D+00                                                    
![    1    1.000D+00                                                    
![    2    2.000D+00                                                    
![    3    3.400D+00                                                    
![    4    5.000D+00                                                    
![                                                                      
![ Y-COORDINATES OF THE NODES                                           
![ --------------------------                                           
![                                                                      
![    I     Y(I)                                                        
![ -----------------                                                    
![    0   -1.500D+00                                                    
![    1    0.000D+00                                                    
![    2    1.000D+00                                                    
![    3    1.500D+00                                                    
![                                                                      
![ GIVEN FUNCTION VALUES A(I,K,0,0) = S( X(I),Y(K) )                    
![ -------------------------------------------------                    
![            0           1           2           3                     
![ ---------------------------------------------------------------------
![    0   0.000D+00   0.000D+00   0.000D+00   0.000D+00                 
![    1   7.074D-02   1.000D+00   5.403D-01   7.074D-02                 
![    2   2.829D-01   4.000D+00   2.161D+00   2.829D-01                 
![    3   8.177D-01   1.156D+01   6.246D+00   8.177D-01                 
![    4   1.768D+00   2.500D+01   1.351D+01   1.768D+00                 
![                                                                      
![ SURFACE NORMAL FN(I,J)                                               
![ ----------------------                                               
![            0           1           2           3                     
![ ---------------------------------------------------------------------
![    0   0.000D+00   0.000D+00   0.000D+00   0.000D+00                 
![    0   0.000D+00   0.000D+00   0.000D+00   0.000D+00                 
![    0   1.000D+00   1.000D+00   1.000D+00   1.000D+00                 
![                                                                      
![    1  -1.415D-01  -2.000D+00  -1.081D+00  -1.415D-01                 
![    1  -9.975D-01   0.000D+00   8.415D-01   9.975D-01                 
![    1   1.000D+00   1.000D+00   1.000D+00   1.000D+00                 
![                                                                      
![    2  -2.829D-01  -4.000D+00  -2.161D+00  -2.829D-01                 
![    2  -3.990D+00   0.000D+00   3.366D+00   3.990D+00                 
![    2   1.000D+00   1.000D+00   1.000D+00   1.000D+00                 
![                                                                      
![    3  -4.810D-01  -6.800D+00  -3.674D+00  -4.810D-01                 
![    3  -1.153D+01   0.000D+00   9.727D+00   1.153D+01                 
![    3   1.000D+00   1.000D+00   1.000D+00   1.000D+00                 
![                                                                      
![    4  -7.074D-01  -1.000D+01  -5.403D+00  -7.074D-01                 
![    4  -2.494D+01   0.000D+00   2.104D+01   2.494D+01                 
![    4   1.000D+00   1.000D+00   1.000D+00   1.000D+00                 
![                                                                      
![                                                                      
![ EXAMPLE FOR APPLICATIONS OF FIBIC2:                                  
![ ===================================                                  
![ INTEGRATION OVER A RECTANGLE: R = [XU,XO] X [YU,YO]                  
![                                                                      
![ WHERE XU =  0.000D+00                                                
![       XO =  5.000D+00                                                
![       YU = -1.500D+00                                                
![       YO =  1.500D+00                                                
![                                                                      
![ RESULT OF INTEGRATION: I(S;R) =  8.28151D+01                         
!                                                                       
!***********************************************************************
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      PARAMETER (N = 4, M = 3, KDIM = 3, LDIM = 3) 
!                                                                       
                                                                        
!                                                                       
!     Variables used in format statements                               
!                                                                       
      CHARACTER FO1555 * 50, FO1560 * 50 
                                                                        
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), X (0:N), Y (0:M) 
      DIMENSION FN (0:N, 0:M, 3), F (80) 
!                                                                       
! assign mesh points                                                    
!                                                                       
      DATA (X (I), I = 0, N) / 0.0D0, 1.0D0, 2.0D0, 3.4D0, 5.0D0 / 
      DATA (Y (I), I = 0, M) / - 1.5D0, 0.0D0, 1.0D0, 1.5D0 / 
!                                                                       
! determine rectangular region for integration                          
!                                                                       
      DATA XU, XO, YU, YO / 0.0D0, 5.0D0, - 1.5D0, 1.5D0 / 
!                                                                       
! compute functional values at nodes                                    
!                                                                       
      DO 10 I = 0, N 
         DO 10 K = 0, M 
            A (I, K, 0, 0) = S (X (I), Y (K) ) 
   10 CONTINUE 
!                                                                       
! compute normal vector at nodes                                        
!                                                                       
      DO 20 I = 0, N 
         DO 20 J = 0, M 
            FN (I, J, 1) = - DSDX (X (I), Y (J) ) 
            FN (I, J, 2) = - DSDY (X (I), Y (J) ) 
            FN (I, J, 3) = 1.D0 
   20 CONTINUE 
!                                                                       
      WRITE ( *, 1520) 
      WRITE ( *, 1530) (I, X (I), I = 0, N) 
      WRITE ( *, 1540) (I, Y (I), I = 0, M) 
      WRITE ( *, 1550) 
      WRITE (FO1555, 1555) M + 1 
      WRITE ( *, FO1555) (I, I = 0, M) 
      WRITE ( *, 1556) 
      DO 30 I = 0, N 
         WRITE (FO1560, 1560) M + 1 
         WRITE ( *, FO1560) I, (A (I, K, 0, 0), K = 0, M) 
   30 END DO 
      WRITE ( *, 1570) 
      WRITE (FO1555, 1555) M + 1 
      WRITE ( *, FO1555) (I, I = 0, M) 
      WRITE ( *, 1556) 
      DO 40 I = 0, N 
         DO 50 K = 1, 3 
            WRITE (FO1560, 1560) M + 1 
            WRITE ( *, FO1560) I, (FN (I, J, K), J = 0, M) 
   50    END DO 
         WRITE ( *, 1565) 
   40 END DO 
      CALL BICSP3 (N, M, A, X, Y, FN, F, IFEHL) 
      IF (IFEHL.NE.0) THEN 
         WRITE ( *, 1500) IFEHL 
         STOP 
      ENDIF 
!                                                                       
!                                                                       
!     The spline coefficients have been computed and put out.           
!     Next we apply FIBIC2:                                             
!     -----------------------------------------------------             
!                                                                       
!                                                                       
      CALL FIBIC2 (N, M, A, X, Y, XU, YU, XO, YO, WERT, IFEHL) 
      IF (IFEHL.NE.0) THEN 
         WRITE ( *, 1620) IFEHL 
         STOP 
      ENDIF 
      WRITE ( *, 1610) XU, XO, YU, YO, WERT 
      STOP 
 1500 FORMAT(1X, 'C[', T76, ']*', /,                                    &
     &       1X, 'C[ *** ERROR IN BICSP3: IFEHL = ',I5,' ***',T76,']*') 
 1520 FORMAT(1X, 'C[', T76, ']*', /,                                    &
     &       1X, 'C[ COMPUTE A BICUBIC SPLINE USING  BICSP3',           &
     &           T76, ']*', /,                                          &
     &       1X, 'C[ ', 38('='), T76, ']*')                             
 1530 FORMAT(1X, 'C[', T76, ']*', /,                                    &
     &       1X, 'C[ X-COORDINATES OF THE NODES', T76, ']*', /,         &
     &       1X, 'C[ ', 26('-'), T76, ']*', /,                          &
     &       1X, 'C[', T76, ']*', /,                                    &
     &       1X, 'C[    I', 5X, 'X(I)', T76, ']*', /,                   &
     &       1X, 'C[ ', 17('-'), T76, ']*', /,                          &
     &      (1X, 'C[ ', 1P, I4, 3X, D10.3, T76, ']*'))                  
 1540 FORMAT(1X, 'C[', T76, ']*', /,                                    &
     &       1X, 'C[ Y-COORDINATES OF THE NODES', T76, ']*', /,         &
     &       1X, 'C[ ', 26('-'), T76, ']*', /,                          &
     &       1X, 'C[', T76, ']*', /,                                    &
     &       1X, 'C[    I', 5X, 'Y(I)', T76, ']*', /,                   &
     &       1X, 'C[ ', 17('-'), T76, ']*', /,                          &
     &      (1X, 'C[ ', 1P, I4, 3X, D10.3, T76, ']*'))                  
 1550 FORMAT(1X, 'C[', T76, ']*', /,                                    &
     &       1X, 'C[ GIVEN FUNCTION VALUES A(I,K,0,0) = S( X(I),',      &
     &           'Y(K) )', T76, ']*', /,                                &
     &       1X, 'C[ ', 49('-'),  T76, ']*')                            
 1555 FORMAT('(1X, ''C[ '',', I2, '(8X,I4), T76, '']*'')') 
 1556 FORMAT(1X, 'C[ ', 70('-'), T76, ']*') 
 1560 FORMAT('(1X, ''C[ '', I4, 1P, ', I2, '(2X,D10.3), T76, '']*'')') 
 1565 FORMAT(1X, 'C[', T76, ']*') 
 1570 FORMAT(1X, 'C[', T76, ']*', /,                                    &
     &       1X, 'C[ SURFACE NORMAL FN(I,J)', T76, ']*', /,             &
     &       1X, 'C[ ', 22('-'), T76, ']*')                             
 1610 FORMAT(1X, 'C[', T76, ']*', /,                                    &
     &       1X, 'C[ EXAMPLE FOR APPLICATIONS OF FIBIC2: ',             &
     &           T76, ']*', /,                                          &
     &       1X, 'C[ ', 35('='), T76, ']*', /,                          &
     &       1X, 'C[ INTEGRATION OVER A RECTANGLE: R = [XU,XO] X ',     &
     &           '[YU,YO]',T76, ']*', /,                                &
     &       1X, 'C[', T76, ']*', /,                                    &
     &       1X, 'C[ WHERE XU = ', 1PD10.3, T76, ']*', /,               &
     &       1X, 'C[', 7X, 'XO = ', D10.3, T76, ']*', /,                &
     &       1X, 'C[', 7X, 'YU = ', D10.3, T76, ']*', /,                &
     &       1X, 'C[', 7X, 'YO = ', D10.3, T76, ']*', /,                &
     &       1X, 'C[', T76, ']*', /,                                    &
     &       1X, 'C[ RESULT OF INTEGRATION: I(S;R) = ',D12.5,           &
     &           T76, ']*')                                             
 1620 FORMAT(1X, 'C[', T76, ']*', /,                                    &
     &       1X, 'C[ *** ERROR IN BICSP2: IFEHL = ',I5,' ***',T76,']*') 
      END PROGRAM TEST                              
!                                                                       
      DOUBLEPRECISION FUNCTION S (X, Y) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Compute functional value                                             
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION X, Y 
      S = X * X * DCOS (Y) 
      RETURN 
      END FUNCTION S                                
!                                                                       
      DOUBLEPRECISION FUNCTION DSDX (X, Y) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Compute partial derivative  DS/DX                            *       
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION X, Y 
      DSDX = 2.0D+00 * X * DCOS (Y) 
      RETURN 
      END FUNCTION DSDX                             
!                                                                       
      DOUBLEPRECISION FUNCTION DSDY (X, Y) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Compute partial derivative   DS/DY                            *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION X, Y 
      DSDY = X * X * ( - DSIN (Y) ) 
      RETURN 
      END FUNCTION DSDY                             
