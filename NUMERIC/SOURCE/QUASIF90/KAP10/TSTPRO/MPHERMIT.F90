      PROGRAM TEST 
!                                         1.30.1993  ( Dubois Guido )   
!*********************************************************************  
!                                                                    *  
!  Test program for the subroutines  HERMIT, HMTVAL and HMTAB.       *  
!  We compute the coefficients of Hermit splines for several end     *  
!  point conditions:                                                 *  
!                                                                    *  
!     1. Natural Hermite spline;                                     *  
!     2. Hermite spline with given curvature radii and the ends.     *  
!                                                                    *  
![  NODES:                                                          ]*  
![                                                                  ]*  
![   I      X(I)       Y(I)      Y1(I)                              ]*  
![  -----------------------------------                             ]*  
![   1     2.0000     3.0000     1.0000                             ]*  
![   2     4.0000     4.0000      .0000                             ]*  
![   3     8.0000     3.0000     -.5000                             ]*  
![                                                                  ]*  
![                                                                  ]*  
![  END POINT CONDITIONS:  NATURAL SPLINE                           ]*  
![                    RB1 =   .000D+00     RBN =   .000D+00         ]*  
![                                                                  ]*  
![  COEFFICIENTS:                                                   ]*  
![                                                                  ]*  
![   I        A(I)           B(I)           C(I)                    ]*  
![  -----------------------------------------------                 ]*  
![   1    .300000D+01    .100000D+01    .000000D+00                 ]*  
![   2    .400000D+01    .000000D+00   -.125000D+00                 ]*  
![                                                                  ]*  
![   I        D(I)           E(I)           F(I)                    ]*  
![  -----------------------------------------------                 ]*  
![   1   -.312500D+00    .125000D+00   -.156250D-01                 ]*  
![   2    .625000D-01   -.195313D-01    .195313D-02                 ]*  
![                                                                  ]*  
![  FUNCTION VALUES AND DERIVATIVES AT THE NODES:                   ]*  
![                                                                  ]*  
![      X0           S             S(1)           S(2)              ]*  
![                  S(3)           S(4)           S(5)              ]*  
![  -----------------------------------------------------           ]*  
![    2.0000    .300000D+01    .100000D+01    .000000D+00           ]*  
![             -.187500D+01    .300000D+01   -.187500D+01           ]*  
![    4.0000    .400000D+01    .000000D+00   -.250000D+00           ]*  
![              .375000D+00   -.468750D+00    .234375D+00           ]*  
![    8.0000    .300000D+01   -.500000D+00    .000000D+00           ]*  
![              .375000D+00    .468750D+00    .234375D+00           ]*  
![                                                                  ]*  
![  TABLE OF VALUES:                                                ]*  
![                                                                  ]*  
![    I       X0         YTAB                                       ]*  
![  -----------------------------                                   ]*  
![    1     2.0000    .300000D+01                                   ]*  
![    2     2.2500    .324559D+01                                   ]*  
![    3     2.5000    .346826D+01                                   ]*  
![    4     2.7500    .365401D+01                                   ]*  
![    5     3.0000    .379688D+01                                   ]*  
![    6     3.2500    .389714D+01                                   ]*  
![    7     3.5000    .395947D+01                                   ]*  
![    8     3.7500    .399110D+01                                   ]*  
![    9     4.0000    .400000D+01                                   ]*  
![   10     4.2500    .399309D+01                                   ]*  
![   11     4.5000    .397540D+01                                   ]*  
![   12     4.7500    .395034D+01                                   ]*  
![   13     5.0000    .391992D+01                                   ]*  
![   14     5.2500    .388503D+01                                   ]*  
![   15     5.5000    .384564D+01                                   ]*  
![   16     5.7500    .380102D+01                                   ]*  
![   17     6.0000    .375000D+01                                   ]*  
![   18     6.2500    .369116D+01                                   ]*  
![   19     6.5000    .362311D+01                                   ]*  
![   20     6.7500    .354465D+01                                   ]*  
![   21     7.0000    .345508D+01                                   ]*  
![   22     7.2500    .335435D+01                                   ]*  
![   23     7.5000    .324335D+01                                   ]*  
![   24     7.7500    .312410D+01                                   ]*  
![   25     8.0000    .300000D+01                                   ]*  
![                                                                  ]*  
![                                                                  ]*  
![  END POINT CONDITIONS:  CURVATURE RADII                          ]*  
![                    RB1 =   .100D+01     RBN =   .100D+01         ]*  
![                                                                  ]*  
![  COEFFICIENTS:                                                   ]*  
![                                                                  ]*  
![   I        A(I)           B(I)           C(I)                    ]*  
![  -----------------------------------------------                 ]*  
![   1    .300000D+01    .100000D+01    .141421D+01                 ]*  
![   2    .400000D+01    .000000D+00    .266911D+00                 ]*  
![                                                                  ]*  
![   I        D(I)           E(I)           F(I)                    ]*  
![  -----------------------------------------------                 ]*  
![   1   -.223786D+01    .989705D+00   -.143413D+00                 ]*  
![   2   -.567404D-01   -.333944D-01    .674782D-02                 ]*  
![                                                                  ]*  
![  FUNCTION VALUES AND DERIVATIVES AT THE NODES:                   ]*  
![                                                                  ]*  
![      X0           S             S(1)           S(2)              ]*  
![                  S(3)           S(4)           S(5)              ]*  
![  -----------------------------------------------------           ]*  
![    2.0000    .300000D+01    .100000D+01    .282843D+01           ]*  
![             -.134272D+02    .237529D+02   -.172095D+02           ]*  
![    4.0000    .400000D+01    .000000D+00    .533822D+00           ]*  
![             -.340442D+00   -.801465D+00    .809738D+00           ]*  
![    8.0000    .300000D+01   -.500000D+00    .139754D+01           ]*  
![              .293160D+01    .243749D+01    .809738D+00           ]*  
![                                                                  ]*  
![  TABLE OF VALUES:                                                ]*  
![                                                                  ]*  
![    I       X0         YTAB                                       ]*  
![  -----------------------------                                   ]*  
![    1     2.0000    .300000D+01                                   ]*  
![    2     2.2500    .330715D+01                                   ]*  
![    3     2.5000    .363120D+01                                   ]*  
![    4     2.7500    .388051D+01                                   ]*  
![    5     3.0000    .402264D+01                                   ]*  
![    6     3.2500    .406749D+01                                   ]*  
![    7     3.5000    .405053D+01                                   ]*  
![    8     3.7500    .401597D+01                                   ]*  
![    9     4.0000    .400000D+01                                   ]*  
![   10     4.2500    .401567D+01                                   ]*  
![   11     4.5000    .405776D+01                                   ]*  
![   12     4.7500    .411724D+01                                   ]*  
![   13     5.0000    .418352D+01                                   ]*  
![   14     5.2500    .424529D+01                                   ]*  
![   15     5.5000    .429123D+01                                   ]*  
![   16     5.7500    .431087D+01                                   ]*  
![   17     6.0000    .429534D+01                                   ]*  
![   18     6.2500    .423818D+01                                   ]*  
![   19     6.5000    .413612D+01                                   ]*  
![   20     6.7500    .398989D+01                                   ]*  
![   21     7.0000    .380498D+01                                   ]*  
![   22     7.2500    .359246D+01                                   ]*  
![   23     7.5000    .336975D+01                                   ]*  
![   24     7.7500    .316143D+01                                   ]*  
![   25     8.0000    .300000D+01                                   ]*  
!                                                                    *  
!  The results are displayed on screen, but can also be sent to a    *  
!  file.                                                             *  
!                                                                    *  
!*********************************************************************  
!                                                                       
      INTEGER N, INTV, NTAB 
      PARAMETER (N = 3, INTV = 24, NTAB = INTV + N) 
!                                                                       
      DOUBLEPRECISION X (1:N), Y (1:N), Y1 (1:N), RB1, RBN, A (1:N),    &
      B (1:N), C (1:N), D (1:N), E (1:N), F (1:N), H (1:N), SUP (1:N),  &
      AINF (1:N), PRC (1:N), DXT (1:N), AR1 (1:N), AR2 (1:N), AR3 (1:N),&
      HMTVAL, AUSG (1:5), DELTX, S, XTAB (1:NTAB), YTAB (1:NTAB)        
      INTEGER MARG, IREP, MORS, I, J, LENTAB, IERR 
!                                                                       
!     initialize data; exchange if desired                              
!                                                                       
      DATA X / 2.0D0, 4.0D0, 8.0D0 / 
      DATA Y / 3.0D0, 4.0D0, 3.0D0 / 
      DATA Y1 / 1.0D0, 0.0D0, - 0.5D0 / 
!                                                                       
!     put out test example                                              
!                                                                       
      WRITE ( *, 900) 
      WRITE ( *, 910) (I, X (I), Y (I), Y1 (I), I = 1, N) 
!                                                                       
!     compute hermite spline coefficients:                              
!                                                                       
!     1. Natural Hermite spline                                         
!                                                                       
!     2. Hermite spline with given curvature radii                      
!        at end points                                                  
!                                                                       
      IREP = 0 
      DO 10 MARG = 2, 4, 2 
         IF (MARG.EQ.2) THEN 
            RB1 = 0.0D0 
            RBN = 0.0D0 
            WRITE ( *, 920) RB1, RBN 
         ELSEIF (MARG.EQ.4) THEN 
            RB1 = 1.0D0 
            RBN = 1.0D0 
            WRITE ( *, 930) RB1, RBN 
         ENDIF 
!                                                                       
         CALL HERMIT (N, MARG, X, Y, Y1, RB1, RBN, IREP, A, B, C, D, E, &
         F, MORS, H, SUP, AINF, PRC, DXT, AR1, AR2, AR3)                
         IREP = 1 
         IF (MORS.EQ.0) THEN 
!                                                                       
!     output of coefficients                                            
!                                                                       
            WRITE ( *, 940) 
            DO 20 I = 1, N - 1 
               WRITE ( *, 960) I, A (I), B (I), C (I) 
   20       END DO 
            WRITE ( *, 950) 
            DO 30 I = 1, N - 1 
               WRITE ( *, 960) I, D (I), E (I), F (I) 
   30       END DO 
!                                                                       
!     compute and print function values and those of derivatives        
!     at the nodes                                                      
!                                                                       
            WRITE ( *, 970) 
            DO 40 I = 1, N 
!                                                                       
               S = HMTVAL (N, X (I), A, B, C, D, E, F, X, AUSG) 
               WRITE ( *, 980) X (I), S, AUSG (1), AUSG (2) 
               WRITE ( *, 990) (AUSG (J), J = 3, 5) 
   40       END DO 
!                                                                       
!     table of values for                                               
!     X0 = X(1) (DELTX) X(N)                                            
!                                                                       
            DELTX = (X (N) - X (1) ) / INTV 
!                                                                       
            CALL HMTAB (N, NTAB, X (1), X (N), DELTX, X, A, B, C, D, E, &
            F, XTAB, YTAB, LENTAB, IERR)                                
            IF (IERR.EQ.0) THEN 
               WRITE ( *, 1000) 
               DO 50 I = 1, LENTAB 
                  WRITE ( *, 1010) I, XTAB (I), YTAB (I) 
   50          END DO 
            ELSE 
               WRITE ( *, 1020) IERR 
            ENDIF 
         ELSE 
            WRITE ( *, 1030) MORS 
         ENDIF 
   10 END DO 
      STOP 
!                                                                       
  900 FORMAT(1X,'C[  NODES:',T70,']*',/,                                &
     &       1X,'C[',T70,']*',/,                                        &
     &       1X,'C[   I',6X,'X(I)',7X,'Y(I)',6X,'Y1(I)',T70,']*',/,     &
     &       1X,'C[  ',35('-'),T70,']*')                                
  910 FORMAT(1X,'C[  ',I2,3X,F8.4,3X,F8.4,3X,F8.4,T70,']*') 
  920 FORMAT(1X,'C[',T70,']*',/,                                        &
     &       1X,'C[',T70,']*',/,                                        &
     &       1X,'C[  END POINT CONDITIONS:  NATURAL SPLINE',T70,']*',/, &
     &       1X,'C[',20X,'RB1 = ',D10.3,5X,'RBN = ',D10.3,T70,']*',/,   &
     &       1X,'C[',T70,']*')                                          
  930 FORMAT(1X,'C[',T70,']*',/,                                        &
     &       1X,'C[',T70,']*',/,                                        &
     &       1X,'C[  END POINT CONDITIONS:  CURVATURE RADII',T70,']*',/,&
     &       1X,'C[',20X,'RB1 = ',D10.3,5X,'RBN = ',D10.3,T70,']*',/,   &
     &       1X,'C[',T70,']*')                                          
  940 FORMAT(1X,'C[  COEFFICIENTS:',T70,']*',/,                         &
     &       1X,'C[',T70,']*',/,                                        &
     &       1X,'C[   I',8X,'A(I)',11X,'B(I)',11X,'C(I)',T70,']*',/,    &
     &       1X,'C[  ',47('-'),T70,']*')                                
  950 FORMAT(1X,'C[',T70,']*',/,                                        &
     &       1X,'C[   I',8X,'D(I)',11X,'E(I)',11X,'F(I)',T70,']*',/,    &
     &       1X,'C[  ',47('-'),T70,']*')                                
  960 FORMAT(1X,'C[  ',I2,3(3X,D12.6),T70,']*') 
  970 FORMAT(1X,'C[',T70,']*',/,                                        &
     &       1X,'C[  FUNCTION VALUES AND DERIVATIVES AT THE ',          &
     &          'NODES:',T70,']*',/,                                    &
     &       1X,'C[',T70,']*',/,                                        &
     &       1X,'C[',6X,'X0',11X,'S',13X,'S(1)',11X,'S(2)',T70,']*',/,  &
     &       1X,'C[',18X,'S(3)',11X,'S(4)',11X,'S(5)',T70,']*',/,       &
     &       1X,'C[  ',53('-'),T70,']*')                                
  980 FORMAT(1X,'C[  ',F8.4,3(3X,D12.6),T70,']*') 
  990 FORMAT(1X,'C[',10X,3(3X,D12.6),T70,']*') 
 1000 FORMAT(1X,'C[',T70,']*',/,                                        &
     &       1X,'C[  TABLE OF VALUES:',T70,']*',/,                      &
     &       1X,'C[',T70,']*',/,                                        &
     &       1X,'C[',4X,'I',7X,'X0',9X,'YTAB',T70,']*',/,               &
     &       1X,'C[  ',29('-'),T70,']*')                                
 1010 FORMAT(1X,'C[  ',I3,3X,F8.4,3X,D12.6,T70,']*') 
 1020 FORMAT(1X,'C[',T70,']*',/,                                        &
     &       1X,'C[  ERROR IN SUBROUTINE HMTAB, IERR = ',I4,T70,']*')   
 1030 FORMAT(1X,'C[',T70,']*',/,                                        &
     &       1X,'C[  ERROR IN SUBROUTINE HERMIT, MORS = ',I4,T70,']*')  
      END PROGRAM TEST                              
