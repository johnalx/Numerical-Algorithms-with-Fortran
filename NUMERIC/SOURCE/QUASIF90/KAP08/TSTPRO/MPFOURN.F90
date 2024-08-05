      PROGRAM TESFFN 
!                                                                       
!****************************************************************       
!                                                               *       
!  This program tests the subroutine  FOURN  with the values of *       
!  the function  F(X)  defined as                               *       
!        0.5          for  X = 0                                *       
!       (PI-X)/PI     for  X  in (0, PI]                        *       
!        0            otherwise,                                *       
!  whose (nonperiodic) Fourier transform  F^(V)  is given by    *       
!                                                               *       
!  (1-COS(2*PI*PI*V)+I*(SIN(2*PI*PI*V)-2*PI*PI*V)/(4*PI**3*V*V).*       
!                                                               *       
!  Here I denotes the imaginary unit (I**2 = -1).               *       
!  A second test function is  G(X)  defined as                  *       
!       SIN(X)        for  X  in [0, PI]                        *       
!        0            otherwise,                                *       
!  whose (nonperiodic) Fourier transform G^(V)  is              *       
!                                                               *       
!     (COS(2*PI*PI*V)+1-I*SIN(2*PI*PI*V))/(1-4*PI*PI*V*V).      *       
!                                                               *       
!  In both cases we use equidistant nodes. For increasing       *       
!  numbers of nodes, the accuracy must increase; specifically   *       
!  it must do so more rapidly for G(X) than for F(X).           *       
!                                                               *       
!****************************************************************       
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER M (7) 
      DIMENSION FUNKRE (0:1365), FUNKIM (0:1365) 
      DIMENSION GUNKRE (0:1365), GUNKIM (0:1365) 
      FEXRE (V) = (1.0D0 - COS (2.0D0 * PI2 * V) ) / (4.0D0 * PI2 * PI *&
      V * V)                                                            
      FEXIM (V) = (SIN (2.0D0 * PI2 * V) - 2.0D0 * PI2 * V) / (4.0D0 *  &
      PI2 * PI * V * V)                                                 
      GEXRE (V) = (COS (2.0D0 * PI2 * V) + 1.0D0) / (1.0D0 - 4.0D0 *    &
      PI2 * V * V)                                                      
      GEXIM (V) = - SIN (2.0D0 * PI2 * V) / (1.0D0 - 4.0D0 * PI2 * V *  &
      V)                                                                
      DATA M / 16, 64, 288, 512, 798, 1024, 1366 / 
!                                                                       
!     assign functional values                                          
!                                                                       
      PI = 4.0D0 * ATAN (1.0D0) 
      PI2 = PI * PI 
      DO 50 I = 1, 7 
         FAKT = PI / DBLE (M (I) ) 
         DO 10 J = 0, M (I) - 1 
            XJ = DBLE (J) * FAKT 
            FUNKRE (J) = (PI - XJ) / PI 
            FUNKIM (J) = 0.0D0 
            GUNKRE (J) = SIN (XJ) 
            GUNKIM (J) = 0.0D0 
   10    END DO 
         FUNKRE (0) = 0.5D0 
!                                                                       
!     Find Fourier transform                                            
!                                                                       
         CALL FOURN (M (I), FUNKRE, FUNKIM, 0.0D0, FAKT) 
         CALL FOURN (M (I), GUNKRE, GUNKIM, 0.0D0, FAKT) 
!                                                                       
!     Compare with theoretical result                                   
!                                                                       
         X = DBLE (M (I) ) * FAKT 
         FEHLRE = 0.0D0 
         FEHLIM = 0.0D0 
         GEHLRE = 0.0D0 
         GEHLIM = 0.0D0 
         DO 20 K = - M (I) / 2, - 1 
            TK = DBLE (K) / X 
            FEHLRE = MAX (FEHLRE, ABS (FEXRE (TK) - FUNKRE (K + M (I) ) &
            ) )                                                         
            FEHLIM = MAX (FEHLIM, ABS (FEXIM (TK) - FUNKIM (K + M (I) ) &
            ) )                                                         
            GEHLRE = MAX (GEHLRE, ABS (GEXRE (TK) - GUNKRE (K + M (I) ) &
            ) )                                                         
            GEHLIM = MAX (GEHLIM, ABS (GEXIM (TK) - GUNKIM (K + M (I) ) &
            ) )                                                         
   20    END DO 
         TK = 0.0D0 
         FEHLRE = MAX (FEHLRE, ABS (0.5D0 * PI - FUNKRE (0) ) ) 
         FEHLIM = MAX (FEHLIM, ABS (FUNKIM (0) ) ) 
         GEHLRE = MAX (GEHLRE, ABS (GEXRE (TK) - GUNKRE (0) ) ) 
         GEHLIM = MAX (GEHLIM, ABS (GEXIM (TK) - GUNKIM (0) ) ) 
         DO 30 K = 1, M (I) / 2 - 1 
            TK = DBLE (K) / X 
            FEHLRE = MAX (FEHLRE, ABS (FEXRE (TK) - FUNKRE (K) ) ) 
            FEHLIM = MAX (FEHLIM, ABS (FEXIM (TK) - FUNKIM (K) ) ) 
            GEHLRE = MAX (GEHLRE, ABS (GEXRE (TK) - GUNKRE (K) ) ) 
            GEHLIM = MAX (GEHLIM, ABS (GEXIM (TK) - GUNKIM (K) ) ) 
   30    END DO 
         WRITE ( * , 900) M (I) , 'F(X)', FEHLRE, FEHLIM 
         WRITE ( * , 900) M (I) , 'G(X)', GEHLRE, GEHLIM 
         WRITE ( *, * ) 
   50 END DO 
      STOP 
  900 FORMAT(1X,'M = ',I4,'  Error at  ',A,' : ',2(1PD10.3,2X)) 
      END PROGRAM TESFFN                            
