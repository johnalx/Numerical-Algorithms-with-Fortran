      PROGRAM TEST 
!                                                             88/07/25  
!*********************************************************************  
!                                                                    *  
!  Test programm for the Adams-Bashforth-Moulton method              *  
!                                                                    *  
!  We perform the following tests:                                   *  
!                                                                    *  
![                                                                  ]*  
![ FUNCTION-NO.:  1                                                 ]*  
![ INITIAL VALUES:  X= .00000000E+00        Y= .10000000E+01        ]*  
![                                     RELEPS= .10000000E-07        ]*  
![ AFTER CALL:      X= .95000000E+00        Y= .25857096E+01        ]*  
![                                        YEX= .25857097E+01        ]*  
![              IFEHL=   1             ABSERR= .61702138E-07        ]*  
![             NEXT H= .50000000E-01   RELERR= .23862749E-07        ]*  
![ AFTER CALL:      X= .19000000E+01        Y= .66858941E+01        ]*  
![                                        YEX= .66858944E+01        ]*  
![              IFEHL=   1             ABSERR= .34458019E-06        ]*  
![             NEXT H= .50000000E-01   RELERR= .51538386E-07        ]*  
![ AT INTERVAL END:                                                 ]*  
![                0 UNSUCCESSFUL TRIES                              ]*  
![                0 TIMES THE ACCURACY WAS DECREASED                ]*  
![                                                                  ]*  
![ FUNCTION-NO.:  2                                                 ]*  
![ INITIAL VALUES:  X= .00000000E+00        Y= .10000000E+01        ]*  
![                                     RELEPS= .10000000E-07        ]*  
![ AFTER CALL:      X= .95000000E+00        Y= .38674101E+00        ]*  
![                                        YEX= .38674102E+00        ]*  
![              IFEHL=   1             ABSERR= .11224186E-07        ]*  
![             NEXT H= .50000000E-01   RELERR= .29022486E-07        ]*  
![ AFTER CALL:      X= .19000000E+01        Y= .14956861E+00        ]*  
![                                        YEX= .14956862E+00        ]*  
![              IFEHL=   1             ABSERR= .93764861E-08        ]*  
![             NEXT H= .50000000E-01   RELERR= .62690200E-07        ]*  
![ AT INTERVAL END:                                                 ]*  
![                0 UNSUCCESSFUL TRIES                              ]*  
![                0 TIMES THE ACCURACY WAS DECREASED                ]*  
![                                                                  ]*  
![ FUNCTION-NO.:  3                                                 ]*  
![ INITIAL VALUES:  X= .00000000E+00        Y= .10000000E+01        ]*  
![                                     RELEPS= .10000000E-07        ]*  
![ AFTER CALL:      X= .97500000E+00        Y= .18810957E+01        ]*  
![                                        YEX= .18810957E+01        ]*  
![              IFEHL=   1             ABSERR= .20232021E-08        ]*  
![             NEXT H= .25000000E-01   RELERR= .10755445E-08        ]*  
![ AFTER CALL:      X= .19625000E+01        Y= .30110410E+02        ]*  
![                                        YEX= .30110410E+02        ]*  
![              IFEHL=   1             ABSERR= .22762094E-08        ]*  
![             NEXT H= .25000000E-01   RELERR= .75595431E-10        ]*  
![ AT INTERVAL END:                                                 ]*  
![                0 UNSUCCESSFUL TRIES                              ]*  
![                0 TIMES THE ACCURACY WAS DECREASED                ]*  
![                                                                  ]*  
![ FUNCTION-NO.:  4                                                 ]*  
![ INITIAL VALUES:  X= .00000000E+00        Y= .10000000E+01        ]*  
![                                     RELEPS= .10000000E-07        ]*  
![ AFTER CALL:      X= .10000000E+01        Y= .60653066E+00        ]*  
![                                        YEX= .60653066E+00        ]*  
![              IFEHL=   0             ABSERR= .13253693E-08        ]*  
![             NEXT H= .41666667E-01   RELERR= .21851645E-08        ]*  
![ AT INTERVAL END:                                                 ]*  
![                0 UNSUCCESSFUL TRIES                              ]*  
![                0 TIMES THE ACCURACY WAS DECREASED                ]*  
![                                                                  ]*  
![ FUNCTION-NO.:  5                                                 ]*  
![ INITIAL VALUES:  X= .00000000E+00        Y= .10000000E+01        ]*  
![                                     RELEPS= .10000000E-07        ]*  
![ AFTER CALL:      X= .10000000E+01        Y= .12500000E+01        ]*  
![                                        YEX= .12500000E+01        ]*  
![              IFEHL=   0             ABSERR= .22204460E-15        ]*  
![             NEXT H= .33333333E+00   RELERR= .17763568E-15        ]*  
![ AT INTERVAL END:                                                 ]*  
![                0 UNSUCCESSFUL TRIES                              ]*  
![                0 TIMES THE ACCURACY WAS DECREASED                ]*  
![                                                                  ]*  
![ FUNCTION-NO.:  6                                                 ]*  
![ INITIAL VALUES:  X= .00000000E+00     Y(1)= .10000000E+01        ]*  
![                                       Y(2)= .00000000E+00        ]*  
![                                     RELEPS= .10000000E-07        ]*  
![ AFTER CALL:      X= .95000000E+00     Y(1)= .58168311E+00        ]*  
![                                       Y(2)= .81341552E+00        ]*  
![                                     YEX(1)= .58168309E+00        ]*  
![                                     YEX(2)= .81341550E+00        ]*  
![              IFEHL=   1             ABSERR= .19791918E-07        ]*  
![             NEXT H= .50000000E-01   RELERR= .24331866E-07        ]*  
![ AFTER CALL:      X= .19000000E+01     Y(1)=-.32328958E+00        ]*  
![                                       Y(2)= .94630014E+00        ]*  
![                                     YEX(1)=-.32328957E+00        ]*  
![                                     YEX(2)= .94630009E+00        ]*  
![              IFEHL=   1             ABSERR= .55283536E-07        ]*  
![             NEXT H= .50000000E-01   RELERR= .58420720E-07        ]*  
![ AT INTERVAL END:                                                 ]*  
![                0 UNSUCCESSFUL TRIES                              ]*  
![                0 TIMES THE ACCURACY WAS DECREASED                ]*  
!                                                                    *  
!*********************************************************************  
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      PARAMETER (NDIM = 2) 
!                                                                       
      DIMENSION Y (NDIM), FHILF (NDIM, 12), YEX (NDIM) 
      LOGICAL JUMP 
!                                                                       
      COMMON / FUNKTI / IFKT, XEND, JUMP 
!                                                                       
      JUMP = .TRUE. 
!.                                                                      
!     Interval end                                                      
!.                                                                      
      XEND = 1.0D0 
!.                                                                      
!     loop over all test functions                                      
!.                                                                      
      DO 60 IFKT = 1, 6 
         CALL BERECH (NDIM, Y, FHILF, YEX) 
   60 END DO 
      STOP 
!                                                                       
      END PROGRAM TEST                              
!                                                                       
!                                                                       
      SUBROUTINE BERECH (NDIM, Y, FHILF, YEX) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION Y (NDIM), FHILF (NDIM, 12), YEX (NDIM) 
      LOGICAL JUMP 
      COMMON / FUNKTI / IFKT, XEND, JUMP 
      EXTERNAL FKTN 
      RELEPS = 1.0D-08 
      ABSEPS = RELEPS 
      H = 0.1D0 
      HMAX = 10.0D0 * H 
      X = 0.0D0 
      Y (1) = 1.0D0 
      Y (2) = 0.0D0 
      N = 1 
      IF (IFKT.EQ.6) N = 2 
      IF (N.EQ.1) THEN 
         WRITE ( *, 3500) IFKT, X, Y (1), RELEPS 
      ELSEIF (N.EQ.2) THEN 
         WRITE ( *, 3600) IFKT, X, Y (1), Y (2), RELEPS 
      ENDIF 
      IFEHL = 0 
      IFZAH1 = 0 
      IFZAH2 = 0 
      IND = 1 
   20 XEND0 = X + 1.0D0 
!.                                                                      
!     Test function No. 5 has discontinuous first derivative at ABS(X)=0
!     integrate accordingly                                             
!.                                                                      
      IF (IFKT.EQ.5.AND.DABS (XEND0) .GT.0.5D0.AND.JUMP) THEN 
         H = 0.5D0 - X 
         IFEHL = 0 
         JUMP = .NOT.JUMP 
      ENDIF 
   30 CALL DESABM (X, Y, FKTN, N, XEND0, H, HMAX, ABSEPS, RELEPS, IND,  &
      IFEHL, FHILF, NDIM)                                               
      CALL LFKTN (X, YEX, N) 
      ABSERR = 0.0D0 
      DO 40 I = 1, N 
         ABSERR = DMAX1 (ABSERR, DABS (Y (I) - YEX (I) ) ) 
         RELERR = ABSERR / DABS (Y (I) ) 
   40 END DO 
      IF (IFEHL.LE.1) THEN 
         IF (N.EQ.1) THEN 
            WRITE ( *, 3700) X, Y (1), YEX (1), IFEHL, ABSERR, H,       &
            RELERR                                                      
         ELSEIF (N.EQ.2) THEN 
            WRITE ( *, 3800) X, Y (1), Y (2), YEX (1), YEX (2), IFEHL,  &
            ABSERR, H, RELERR                                           
         ENDIF 
      ELSEIF (IFEHL.EQ.2) THEN 
         IFZAHL = IFZAHL + 1 
         GOTO 30 
      ELSEIF (IFEHL.EQ.3) THEN 
         H = 0.1D0 
         RELEPS = RELEPS * 10.0D0 
         IFZAH2 = IFZAH2 + 1 
         IFZAHL = IFZAHL + 1 
         GOTO 30 
      ELSE 
         WRITE ( *, 3900) IFEHL 
         STOP 
      ENDIF 
      IF (X.LT.XEND) GOTO 20 
      WRITE ( *, 4000) IFZAH1, IFZAH2 
      RETURN 
 3500 FORMAT (1X,'C[',T70,']*',/,                                       &
     &        1X,'C[ FUNCTION NO.:',I3,T70,']*',/,                      &
     &        1X,'C[ INITIAL VALUES:  X=',E14.8,8X,'Y=',E14.8,          &
     &           T70,']*',/,                                            &
     &        1X,'C[',37X,'RELEPS=',E14.8,T70,']*')                     
 3600 FORMAT (1X,'C[',T70,']*',/,                                       &
     &        1X,'C[ FUNCTION NO.:',I3,T70,']*',/,                      &
     &        1X,'C[ INITIAL VALUES:  X=',E14.8,5X,'Y(1)=',E14.8,       &
     &           T70,']*',/,                                            &
     &        1X,'C[',39X,'Y(2)=',E14.8,T70,']*',/,                     &
     &        1X,'C[',37X,'RELEPS=',E14.8,T70,']*')                     
 3700 FORMAT (1X,'C[ AFTER CALL:      X=',E14.8,8X,'Y=',E14.8,          &
     &           T70,']*',/,                                            &
     &        1X,'C[',40X,'YEX=',E14.8,T70,']*',/,                      &
     &        1X,'C[',14X,'IFEHL= ',I3,13X,'ABSERR=',E14.8,T70,']*',/,  &
     &        1X,'C[',13X,'NEXT H=',E14.8,3X,'RELERR=',E14.8,T70,']*')  
 3800 FORMAT (1X,'C[ AFTER CALL:      X=',E14.8,5X,'Y(1)=',E14.8,       &
     &           T70,']*',/,                                            &
     &        1X,'C[',39X,'Y(2)=',E14.8,T70,']*',/,                     &
     &        1X,'C[',37X,'YEX(1)=',E14.8,T70,']*',/,                   &
     &        1X,'C[',37X,'YEX(2)=',E14.8,T70,']*',/,                   &
     &        1X,'C[',14X,'IFEHL= ',I3,13X,'ABSERR=',E14.8,T70,']*',/,  &
     &        1X,'C[',13X,'NEXT H=',E14.8,3X,'RELERR=',E14.8,T70,']*')  
 3900 FORMAT (1X,'C[ ERROR IN INPUT PARAMETERS, IFEHL=',I2,T70,']*') 
 4000 FORMAT (1X,'C[ AT INTERVAL END:',T70,']*',/,                      &
     &        1X,'C[',13X,I4,' UNSUCCESSFUL TRIES',T70,']*',/,          &
     &        1X,'C[',13X,I4,' TIMES THE ACCURACY WAS DECREASED',       &
     &           T70,']*')                                              
      END SUBROUTINE BERECH                         
!                                                                       
!                                                                       
!    various right hand sides of DEs (N-dimensional)                    
!                                                                       
      SUBROUTINE FKTN (X, Y, N, F) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION Y (N), F (N) 
      LOGICAL JUMP 
!                                                                       
      COMMON / FUNKTI / INUM, XEND, JUMP 
!                                                                       
      GOTO (10, 20, 30, 40, 50, 60) INUM 
   10 F (1) = Y (1) 
      RETURN 
!                                                                       
   20 F (1) = - Y (1) 
      RETURN 
!                                                                       
   30 F (1) = 5 * X**4 
      RETURN 
!                                                                       
   40 F (1) = - X * Y (1) 
      RETURN 
!                                                                       
   50 IF (DABS (X) .LT.0.5D0) THEN 
         F (1) = 0.0D0 
      ELSEIF (X.GE.0.5D0) THEN 
         F (1) = 2.0D0 * (X - 0.5D0) 
      ELSEIF (X.LE. - 0.5D0) THEN 
         F (1) = 2.0D0 * (X + 0.5D0) 
      ENDIF 
      RETURN 
!                                                                       
   60 F (1) = - Y (2) 
      F (2) = Y (1) 
      RETURN 
!                                                                       
      END SUBROUTINE FKTN                           
!                                                                       
!    solutions for the DEs                                              
!                                                                       
      SUBROUTINE LFKTN (X, YEX, N) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION YEX (N) 
      LOGICAL JUMP 
!                                                                       
      COMMON / FUNKTI / INUM, XEND, JUMP 
!                                                                       
      GOTO (10, 20, 30, 40, 50, 60) INUM 
!                                                                       
   10 YEX (1) = DEXP (X) 
      RETURN 
!                                                                       
   20 YEX (1) = DEXP ( - X) 
      RETURN 
!                                                                       
   30 YEX (1) = X**5 + 1.0D0 
      RETURN 
!                                                                       
   40 YEX (1) = DEXP ( - .5D0 * X * X) 
      RETURN 
!                                                                       
   50 IF (DABS (X) .LT.0.5D0) THEN 
         YEX (1) = 1.0D0 
      ELSEIF (X.GE.0.5D0) THEN 
         YEX (1) = (X - 0.5D0) **2 + 1.0D0 
      ELSEIF (X.LE. - 0.5D0) THEN 
         YEX (1) = (X + 0.5D0) **2 + 1.0D0 
      ENDIF 
      RETURN 
!                                                                       
   60 YEX (1) = DCOS (X) 
      YEX (2) = DSIN (X) 
      RETURN 
!                                                                       
      END SUBROUTINE LFKTN                          
