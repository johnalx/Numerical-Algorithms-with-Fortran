      PROGRAM TEST 
!                                                                       
!*********************************************************************  
!                                                                    *  
!     Test program for subroutine MULLRP.                            *  
!     Compute all zeros of a polynomial via Muller's method.         *  
!                                                                    *  
!     TEST EXAMPLE:                                                  *  
!     =============                                                  *  
!     FIND ZEROS OF THE POLYNOMIAL Y=3+7*X+8*X**2+8*X**3+5*X**4+X**5 *  
!                                                                    *  
!     MAXIMAL NUMBER OF ITERATIONS:                                  *  
!     100                                                            *  
!                                                                    *  
![ N = 5 POLYNOMIAL COEFFICIENTS =   .30D+01   .70D+01   .80D+01    ]*  
![    .80D+01   .50D+01   .10D+01                                   ]*  
![ NGEF = 5 ZERO =  -.10D+01   .00D+00   .43D-17   .10D+01          ]*  
![    .43D-17  -.10D+01  -.10D+01   .00D+00  -.30D+01   .00D+00     ]*  
![ ZERO ( 1)        = -.100000000038D+01  .000000000000D+00         ]*  
![     FUNCTION VALUE   = -.347267054831D-18  .000000000000D+00     ]*  
![     ERROR ESTIMATION =  .294209101873D-13  .000000000000D+00     ]*  
![ ZERO ( 2)        =  .432671294359D-17  .100000000000D+01         ]*  
![     FUNCTION VALUE   =  .112584757756D-33  .260208521397D-16     ]*  
![     ERROR ESTIMATION =  .299760216649D-14  .333066907388D-14     ]*  
![ ZERO ( 3)        =  .432671294359D-17 -.100000000000D+01         ]*  
![     FUNCTION VALUE   =  .112584757756D-33 -.260208521397D-16     ]*  
![     ERROR ESTIMATION =  .299760216649D-14  .333066907388D-14     ]*  
![ ZERO ( 4)        = -.999999999617D+00  .000000000000D+00         ]*  
![     FUNCTION VALUE   =  .221697338211D-15  .000000000000D+00     ]*  
![     ERROR ESTIMATION =  .294209101192D-13  .000000000000D+00     ]*  
![ ZERO ( 5)        = -.300000000000D+01  .000000000000D+00         ]*  
![     FUNCTION VALUE   =  .461852778244D-13  .000000000000D+00     ]*  
![     ERROR ESTIMATION =  .192079685490D-11  .000000000000D+00     ]*  
!                                                                    *  
!*********************************************************************  
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      INTEGER MXGRAD, NPOL, ITERMX 
      PARAMETER (MXGRAD = 100, ITERMX = 100, NPOL = 5) 
      DOUBLEPRECISION POLYNM (0:MXGRAD), WKARER (0:MXGRAD), ZERO (2, 10) 
      INTEGER NGEF, I, K 
      DOUBLEPRECISION WKARED (0:MXGRAD), XR, XI, FR, FI, XEB1, XEB0 
      DATA POLYNM / 3.D0, 7.D0, 8.D0, 8.D0, 5.D0, 1.D0, 95 * 0.D0 / 
      CALL MULLRP (NPOL, POLYNM, ITERMX, NGEF, ZERO, WKARER, WKARED) 
      WRITE ( * , 70) ' N =', NPOL, ' POLYNOMIAL COEFFICIENTS =',       &
      (POLYNM (I) , I = 0, NPOL)                                        
      WRITE ( * , 80) ' NGEF =', NGEF, ' ZERO =', ( (ZERO (K, I) , K =  &
      1, 2) , I = 1, NGEF)                                              
      XEPS = YEPS () 
      DO 10 I = 1, NGEF 
         XR = ZERO (1, I) 
         XI = ZERO (2, I) 
         CALL HORNCE (XEPS, FER, FEI, XEB1, XEB0, POLYNM, NPOL, XR, XI, &
         FR, FI, WKARED)                                                
         WRITE ( *, 90) I, XR, XI 
         WRITE ( *, 100) FR, FI 
         WRITE ( *, 110) FER, FEI 
   10 END DO 
      STOP 
   70 FORMAT(1X,'C[',A,I2,A,3D10.2,T70,']*',/,                          &
     &       1X,'C[',1X,3D10.2,T70,']*')                                
   80 FORMAT(1X,'C[',A,I2,A,4D10.2,T70,']*',/,                          &
     &       1X,'C[',1X,6D10.2,T70,']*')                                
   90 FORMAT(1X,'C[',1X,'ZERO (',I2,')        = ',                      &
     &          D18.12,1X,D18.12,T70,']*')                              
  100 FORMAT(1X,'C[',5X,'FUNCTION VALUE   = ',                          &
     &          D18.12,1X,D18.12,T70,']*')                              
  110 FORMAT(1X,'C[',5X,'ERROR ESTIMATION = ',                          &
     &          D18.12,1X,D18.12,T70,']*')                              
      END PROGRAM TEST                              
