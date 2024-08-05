      SUBROUTINE RKTRB (X, BETA, N, DES, Y, ABSERR, RELERR, IFLAG,      &
      WORK1, WORK2, IWORK, IERR)                                        
!                                                                       
!*****************************************************************      
!                                                                *      
! This program solves a system of at most 12 ordinary            *      
! differential equations of first order by using a RUNGE-KUTTA   *      
! embedding formula over the interval of integration             *      
! I= [X0,BETA].                                                  *      
! With the parameter IFLAG(1) one can chose the embedding formula*      
! and IFLAG(2) assists in step size control. When RKTRB is called*      
! first, we must have IFLAG(3)=0. RKTRB checks the input         *      
! parameters and determines the necessary constants.             *      
! The number of integration steps is adjusted so that maximally  *      
! 10000 function evaluations are performed. If IERR = -2 , then  *      
! the maximal number of steps was performed without reaching     *      
! BETA. RKTRB can then be repeated with IFLAG(3)=1. Other        *      
! parameters need not be adjusted.                               *      
! If several intermediate values for the solution in [X0,BETA]   *      
! are desired, then RKTRB must be started with IFLAG(3)=0.       *      
! All subsequent calls can be executed with IFLAG(3)=1 and BETA. *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! X       : DOUBLE PRECISION initial value for the integration:  *      
!           X=X0                                                 *      
! BETA    : DOUBLE PRECISION endpoint X=BETA at which we want to *      
!           find the solution                                    *      
! N       : number of differential equations in the system,      *      
!           or the size of Y:   0 < N < 13                       *      
! DES     : SUBROUTINE DES must be declared as EXTERNAL in the   *      
!           calling program. DES describes the system of         *      
!           differential equations and must have the following   *      
!           form:                                                *      
!                  SUBROUTINE DES(X,Y,N,YPUNKT)                  *      
!                  DOUBLE PRECISION Y(N),YPUNKT(N),X             *      
!                  YPUNKT(1)=....                                *      
!                  YPUNKT(2)=....                                *      
!                         .                                      *      
!                         .                                      *      
!                         .                                      *      
!                   YPUNKT(N)=....                               *      
!                   RETURN                                       *      
!                   END                                          *      
! Y       : DOUBLE PRECISION vector Y(1:N), the solution at X=X0 *      
! ABSERR  : DOUBLE PRECISION error bound for the absolute error  *      
!           (ABSERR >= 0). If ABSERR=0, then only the relative   *      
!           error is checked.                                    *      
! RELERR  : DOUBLE PRECISION error bound for the relative error  *      
!           (RELERR >= 0). If RELERR=0, then only the absolute   *      
!           error is checked.                                    *      
! IFLAG   : INTEGER vector IFLAG(1:3)                            *      
!           IFLAG(1) choses the embedding formula:               *      
!           IFLAG(1)=0 : RK3(2)                                  *      
!           IFLAG(1)=1 : RKF4(3)   (FSAL)                        *      
!           IFLAG(1)=2 : RKF5(4)                                 *      
!           IFLAG(1)=3 : RK5(4)6M                                *      
!           IFLAG(1)=4 : RKE5(4)                                 *      
!           IFLAG(1)=5 : HIHA5     (FSAL)                        *      
!           IFLAG(1)=6 : RK5(4)7S  (FSAL)                        *      
!           IFLAG(1)=7 : RK5(4)7M  (FSAL)                        *      
!           IFLAG(1)=8 : RK5(4)7C  (FSAL)                        *      
!           IFLAG(1)=9 : RK6(5)8M                                *      
!           IFLAG(1)=10: RK6(5)8S                                *      
!           IFLAG(1)=11: RK6(5)8C                                *      
!           IFLAG(1)=12: RKV6(5)                                 *      
!           IFLAG(1)=13: RKF6(5)A                                *      
!           IFLAG(1)=14: RKF6(5)B                                *      
!           IFLAG(1)=15: RKC6(5)   (FSAL)                        *      
!           IFLAG(1)=16: RKV6(5)9A (FSAL)                        *      
!           IFLAG(1)=17: RKV6(5)9B (FSAL)                        *      
!           IFLAG(1)=18: RKV7(6)                                 *      
!           IFLAG(1)=19: RK8(7)13M                               *      
!           IFLAG(1)=20: RKF8(7)                                 *      
!           IFLAG(1)=21: RKV8(7)                                 *      
!           IFLAG(1)=22: RKV9(8)                                 *      
!           IFLAG(2) choses the step size control method:        *      
!           IFLAG(2)=0 : HULL                                    *      
!           IFLAG(2)=1 : CIVPS                                   *      
!           IFLAG(3) describes the first or subsequent calls:    *      
!           IFLAG(3)=0 : first call                              *      
!           IFLAG(3)=1 : subsequent call                         *      
! WORK1   : DOUBLE PRECISION vector WORK1(1:4)                   *      
!           IFLAG(3)=0 : used for storage                        *      
!           IFLAG(3)=1 : WORK1(1) : 100 times the machine        *      
!                                   constant (EPS)               *      
!                        WORK1(2) : largest representable        *      
!                                   number for the computer      *      
!                                   (XZI=DLARGE/100)             *      
!                        WORK1(3) : global error order (QG)      *      
!                        WORK1(4) : initial step size (H)        *      
! WORK2   : 2-dim. DOUBLE PRECISION array WORK2(1:16,1:16)       *      
!           IFLAG(3)=0 : storage                                 *      
!           IFLAG(3)=1 : coefficients of the embedding formula   *      
! IWORK   : INTEGER vector IWORK(1:2)                            *      
!           IFLAG(3)=0 : storage                                 *      
!           IFLAG(3)=1 : IWORK(1) : level of the embedding       *      
!                                   formula                      *      
!                        IWORK(2) : maximally allowed number of  *      
!                                   integrations                 *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! X       : DOUBLE PRECISION value for X, where the integration  *      
!           has stopped (normally X=BETA)                        *      
! Y       : DOUBLE PRECISION solution vector Y(1:N) for X        *      
! WORK1   : DOUBLE PRECISION vector WORK1(1:4)                   *      
!           WORK1(1) : 100 times the machine constant (EPS)      *      
!           WORK1(2) : largest representable number for the      *      
!                      computer/100  (XZI)                       *      
!           WORK1(3) : global error order (QG)                   *      
!           WORK1(4) : last step size (H)                        *      
! WORK2   : 2-dim. DOUBLE PRECISION array WORK2(1:16,1:16), the  *      
!           coefficients of the embedding formula                *      
! IWORK   : INTEGER vector IWORK(1:2)                            *      
!           IWORK(1) : level of the embedding formula            *      
!           IWORK(2) : maximally allowed number of integrations  *      
! IERR    : error parameter:                                     *      
!               IERR=0     all is ok                             *      
!           Control of the transfer parameters:                  *      
!               IERR=-2   IFLAG(3) erroneous                     *      
!               If IERR is maximally 7 digits long:              *      
!               1st error   1st digit equal to 1                 *      
!                           wrong value for IFLAG(1)             *      
!                           IERR=IERR+1000000                    *      
!               2nd error   2nd digit equal to 1                 *      
!                           wrong value for IFLAG(2)             *      
!                           IERR=IERR+100000                     *      
!               3rd error   3rd digit equal to 1                 *      
!                           ABSERR and RELERR too small          *      
!                           IERR=IERR+10000                      *      
!               4th error   4th digit equal to 1                 *      
!                           Interval of integration too small    *      
!                           IERR=IERR+1000                       *      
!               5th error   5th digit equal to 1                 *      
!                           Interval of integration too small    *      
!                           X or BETA not representable in the   *      
!                           computer                             *      
!                           IERR=IERR+100                        *      
!               6th error   6th digit equal to 1                 *      
!                           N <= 0 or N > 12                     *      
!                           IERR=IERR+10                         *      
!               7th error   7th digit equal to 1                 *      
!                           initial condition not representable  *      
!                           in the computer                      *      
!                           IERR=IERR+1                          *      
!               EXAMPLE     If IERR=10101 then three errors      *      
!                           numbered 3, 5 and 7 occur            *      
!           Run time errors:                                     *      
!               IERR=-1    the desired relative accuracy is less *      
!                          than 100 times the machine constant   *      
!                          in certain parts of the integration   *      
!                          interval. In these regions we compute *      
!                          with 100 times the machine constant as*      
!                          an absolute error bound.              *      
!               IERR=-2    the nunber of maximally allowed steps *      
!                          has been reached.                     *      
!               IERR=-20   OVERFLOW, the program stops.          *      
!               IERR=-30   the computed step size is too small.  *      
!                          The program stops.                    *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! DSMALL  : DOUBLE PRECISION machine constant                    *      
! DLARGE  : largest computer DOUBLE PRECISION number             *      
! J       : loop variable                                        *      
! FSAL    : (LOGICAL) variable indicating whether the method     *      
!           FSAL (First Same As Last) is used by the RUNGE-KUTTA *      
!           embedding formula                                    *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: DMPREC, COEFFI, HSTART, HULL, CIVPS     *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Volker KrÅger                                      *      
!  Date     : 29.04.1993                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      EXTERNAL DES 
      DOUBLEPRECISION Y (N), WORK1 (4), WORK2 (16, 16) 
      INTEGER IFLAG (3), IWORK (2) 
      DOUBLEPRECISION X, BETA, ABSERR, RELERR, DSMALL, DLARGE 
      LOGICAL FSAL 
!                                                                       
! check input parameters, if erroneous,                                 
! return to calling program                                             
!                                                                       
      IERR = 0 
      IF (IFLAG (3) .LT.0.OR.IFLAG (3) .GT.1) THEN 
         IERR = 2 
         RETURN 
      ELSEIF (IFLAG (3) .EQ.0) THEN 
!                                                                       
!   determine the machine constant and the largest                      
!   representable number for the computer                               
!                                                                       
         CALL DMPREC (DSMALL, DLARGE) 
         WORK1 (1) = 100.0D0 * DSMALL 
         WORK1 (2) = DLARGE / 100.0D0 
      ENDIF 
!                                                                       
!   determine the appropriate error code                                
!                                                                       
      IF (IFLAG (1) .LT.0.OR.IFLAG (1) .GT.22) IERR = IERR + 1000000 
      IF (IFLAG (2) .LT.0.OR.IFLAG (2) .GT.1) IERR = IERR + 100000 
      IF (RELERR.LT.WORK1 (1) .AND.ABSERR.LT.WORK1 (1) ) IERR = IERR +  &
      10000                                                             
      IF (DABS (BETA - X) .LT.WORK1 (1) ) IERR = IERR + 1000 
      IF (DABS (X) .GT.WORK1 (2) .OR.DABS (BETA) .GT.WORK1 (2) ) IERR = &
      IERR + 100                                                        
      IF (N.LE.0.OR.N.GT.12) IERR = IERR + 10 
      DO 100 J = 1, N 
         IF (DABS (Y (J) ) .GT.WORK1 (2) ) THEN 
            IERR = IERR + 1 
            RETURN 
         ENDIF 
  100 END DO 
      IF (IERR.NE.0) RETURN 
!                                                                       
! End input check                                                       
!                                                                       
! Method FSAL is used by the RUNGE-KUTTA embedding formula              
!                                                                       
      FSAL = (IFLAG (1) .EQ.1.OR.IFLAG (1) .GE.5.AND.IFLAG (1)          &
      .LE.8.OR.IFLAG (1) .GE.15.AND.IFLAG (1) .LE.17)                   
      IF (IFLAG (3) .EQ.0) THEN 
!                                                                       
! Set the level of the embedding formula and initialize                 
! the maximally allowed number of integration steps                     
!                                                                       
         IF (IFLAG (1) .EQ.0) THEN 
            IWORK (1) = 3 
            IWORK (2) = 3330 
         ELSEIF (IFLAG (1) .EQ.1) THEN 
            IWORK (1) = 5 
            IWORK (2) = 2500 
         ELSEIF (IFLAG (1) .GE.2.AND.IFLAG (1) .LE.4) THEN 
            IWORK (1) = 6 
            IWORK (2) = 1670 
         ELSEIF (IFLAG (1) .GE.5.AND.IFLAG (1) .LE.8) THEN 
            IWORK (1) = 7 
            IWORK (2) = 1665 
         ELSEIF (IFLAG (1) .GE.9.AND.IFLAG (1) .LE.14) THEN 
            IWORK (1) = 8 
            IWORK (2) = 1250 
         ELSEIF (IFLAG (1) .GE.15.AND.IFLAG (1) .LE.17) THEN 
            IWORK (1) = 9 
            IWORK (2) = 1250 
         ELSEIF (IFLAG (1) .EQ.18) THEN 
            IWORK (1) = 10 
            IWORK (2) = 1000 
         ELSEIF (IFLAG (1) .GE.19.AND.IFLAG (1) .LE.21) THEN 
            IWORK (1) = 13 
            IWORK (2) = 770 
         ELSEIF (IFLAG (1) .EQ.22) THEN 
            IWORK (1) = 16 
            IWORK (2) = 625 
         ENDIF 
!                                                                       
! Find the coefficients and the global error order of the               
! low order method in the desired RUNGE-KUTTA pair                      
!                                                                       
         IF (IFLAG (1) .EQ.0) THEN 
            WORK1 (3) = 2.0D0 
         ELSE 
            CALL COEFFI (IWORK (1), IFLAG (1), WORK2, WORK1 (3) ) 
         ENDIF 
!                                                                       
! Initialize first step size                                            
!                                                                       
         CALL HSTART (DES, N, X, BETA, Y, RELERR, ABSERR, WORK1 (3),    &
         DSMALL, DLARGE, WORK1 (4) )                                    
      ENDIF 
!                                                                       
! Integrate using HULL or CIVPS                                         
!                                                                       
                                                                        
      IF (IFLAG (2) .EQ.0) THEN 
         CALL HULL (X, WORK1 (4), BETA, ABSERR, RELERR, N, FSAL, IWORK (&
         1), DES, Y, WORK1 (1), WORK1 (2), WORK1 (3), WORK2, IWORK (2), &
         IERR)                                                          
      ELSE 
         CALL CIVPS (X, WORK1 (4), BETA, ABSERR, RELERR, N, FSAL, IWORK &
         (1), DES, Y, WORK1 (1), WORK1 (2), WORK1 (3), WORK2, IWORK (2),&
         IERR)                                                          
      ENDIF 
      RETURN 
      END SUBROUTINE RKTRB                          
