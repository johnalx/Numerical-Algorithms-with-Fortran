C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 17.3.4.4}{Embedding Formulas}{Embedding Formulas}*)
C[LE*)
      SUBROUTINE RKTRB (X,BETA,N,DES,Y,ABSERR,RELERR,IFLAG,WORK1,
C[IX{RKTRB}*)
     +                  WORK2,IWORK,IERR)
C
C*****************************************************************
C                                                                *
C This program solves a system of at most 12 ordinary            *
C differential equations of first order by using a RUNGE-KUTTA   *
C embedding formula over the interval of integration             *
C I= [X0,BETA].                                                  *
C With the parameter IFLAG(1) one can chose the embedding formula*
C and IFLAG(2) assists in step size control.                     *
C[BE*)
C When RKTRB is called                                           *
C first, we must have IFLAG(3)=0. RKTRB checks the input         *
C parameters and determines the necessary constants.             *
C The number of integration steps is adjusted so that maximally  *
C 10000 function evaluations are performed. If IERR = -2 , then  *
C the maximal number of steps was performed without reaching     *
C BETA. RKTRB can then be repeated with IFLAG(3)=1. Other        *
C parameters need not be adjusted.                               *
C If several intermediate values for the solution in [X0,BETA]   *
C are desired, then RKTRB must be started with IFLAG(3)=0.       *
C All subsequent calls can be executed with IFLAG(3)=1 and BETA. *
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C X       : DOUBLE PRECISION initial value for the integration:  *
C           X=X0                                                 *
C BETA    : DOUBLE PRECISION endpoint X=BETA at which we want to *
C           find the solution                                    *
C N       : number of differential equations in the system,      *
C           or the size of Y:   0 < N < 13                       *
C DES     : SUBROUTINE DES must be declared as EXTERNAL in the   *
C           calling program. DES describes the system of         *
C           differential equations and must have the following   *
C           form:                                                *
C                  SUBROUTINE DES(X,Y,N,YPUNKT)                  *
C                  DOUBLE PRECISION Y(N),YPUNKT(N),X             *
C                  YPUNKT(1)=....                                *
C                  YPUNKT(2)=....                                *
C                         .                                      *
C                         .                                      *
C                         .                                      *
C                   YPUNKT(N)=....                               *
C                   RETURN                                       *
C                   END                                          *
C Y       : DOUBLE PRECISION vector Y(1:N), the solution at X=X0 *
C ABSERR  : DOUBLE PRECISION error bound for the absolute error  *
C           (ABSERR >= 0). If ABSERR=0, then only the relative   *
C           error is checked.                                    *
C RELERR  : DOUBLE PRECISION error bound for the relative error  *
C           (RELERR >= 0). If RELERR=0, then only the absolute   *
C           error is checked.                                    *
C IFLAG   : INTEGER vector IFLAG(1:3)                            *
C           IFLAG(1) choses the embedding formula:               *
C           IFLAG(1)=0 : RK3(2)                                  *
C           IFLAG(1)=1 : RKF4(3)   (FSAL)                        *
C           IFLAG(1)=2 : RKF5(4)                                 *
C           IFLAG(1)=3 : RK5(4)6M                                *
C           IFLAG(1)=4 : RKE5(4)                                 *
C           IFLAG(1)=5 : HIHA5     (FSAL)                        *
C           IFLAG(1)=6 : RK5(4)7S  (FSAL)                        *
C           IFLAG(1)=7 : RK5(4)7M  (FSAL)                        *
C           IFLAG(1)=8 : RK5(4)7C  (FSAL)                        *
C           IFLAG(1)=9 : RK6(5)8M                                *
C           IFLAG(1)=10: RK6(5)8S                                *
C           IFLAG(1)=11: RK6(5)8C                                *
C           IFLAG(1)=12: RKV6(5)                                 *
C           IFLAG(1)=13: RKF6(5)A                                *
C           IFLAG(1)=14: RKF6(5)B                                *
C           IFLAG(1)=15: RKC6(5)   (FSAL)                        *
C           IFLAG(1)=16: RKV6(5)9A (FSAL)                        *
C           IFLAG(1)=17: RKV6(5)9B (FSAL)                        *
C           IFLAG(1)=18: RKV7(6)                                 *
C           IFLAG(1)=19: RK8(7)13M                               *
C           IFLAG(1)=20: RKF8(7)                                 *
C           IFLAG(1)=21: RKV8(7)                                 *
C           IFLAG(1)=22: RKV9(8)                                 *
C           IFLAG(2) choses the step size control method:        *
C           IFLAG(2)=0 : HULL                                    *
C           IFLAG(2)=1 : CIVPS                                   *
C           IFLAG(3) describes the first or subsequent calls:    *
C           IFLAG(3)=0 : first call                              *
C           IFLAG(3)=1 : subsequent call                         *
C WORK1   : DOUBLE PRECISION vector WORK1(1:4)                   *
C           IFLAG(3)=0 : used for storage                        *
C           IFLAG(3)=1 : WORK1(1) : 100 times the machine        *
C                                   constant (EPS)               *
C                        WORK1(2) : largest representable        *
C                                   number for the computer      *
C                                   (XZI=DLARGE/100)             *
C                        WORK1(3) : global error order (QG)      *
C                        WORK1(4) : initial step size (H)        *
C WORK2   : 2-dim. DOUBLE PRECISION array WORK2(1:16,1:16)       *
C           IFLAG(3)=0 : storage                                 *
C           IFLAG(3)=1 : coefficients of the embedding formula   *
C IWORK   : INTEGER vector IWORK(1:2)                            *
C           IFLAG(3)=0 : storage                                 *
C           IFLAG(3)=1 : IWORK(1) : level of the embedding       *
C                                   formula                      *
C                        IWORK(2) : maximally allowed number of  *
C                                   integrations                 *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C X       : DOUBLE PRECISION value for X, where the integration  *
C           has stopped (normally X=BETA)                        *
C Y       : DOUBLE PRECISION solution vector Y(1:N) for X        *
C WORK1   : DOUBLE PRECISION vector WORK1(1:4)                   *
C           WORK1(1) : 100 times the machine constant (EPS)      *
C           WORK1(2) : largest representable number for the      *
C                      computer/100  (XZI)                       *
C           WORK1(3) : global error order (QG)                   *
C           WORK1(4) : last step size (H)                        *
C WORK2   : 2-dim. DOUBLE PRECISION array WORK2(1:16,1:16), the  *
C           coefficients of the embedding formula                *
C IWORK   : INTEGER vector IWORK(1:2)                            *
C           IWORK(1) : level of the embedding formula            *
C           IWORK(2) : maximally allowed number of integrations  *
C IERR    : error parameter:                                     *
C               IERR=0     all is ok                             *
C           Control of the transfer parameters:                  *
C               IERR=-2   IFLAG(3) erroneous                     *
C               If IERR is maximally 7 digits long:              *
C               1st error   1st digit equal to 1                 *
C                           wrong value for IFLAG(1)             *
C                           IERR=IERR+1000000                    *
C               2nd error   2nd digit equal to 1                 *
C                           wrong value for IFLAG(2)             *
C                           IERR=IERR+100000                     *
C               3rd error   3rd digit equal to 1                 *
C                           ABSERR and RELERR too small          *
C                           IERR=IERR+10000                      *
C               4th error   4th digit equal to 1                 *
C                           Interval of integration too small    *
C                           IERR=IERR+1000                       *
C               5th error   5th digit equal to 1                 *
C                           Interval of integration too small    *
C                           X or BETA not representable in the   *
C                           computer                             *
C                           IERR=IERR+100                        *
C               6th error   6th digit equal to 1                 *
C                           N <= 0 or N > 12                     *
C                           IERR=IERR+10                         *
C               7th error   7th digit equal to 1                 *
C                           initial condition not representable  *
C                           in the computer                      *
C                           IERR=IERR+1                          *
C               EXAMPLE     If IERR=10101 then three errors      *
C                           numbered 3, 5 and 7 occur            *
C           Run time errors:                                     *
C               IERR=-1    the desired relative accuracy is less *
C                          than 100 times the machine constant   *
C                          in certain parts of the integration   *
C                          interval. In these regions we compute *
C                          with 100 times the machine constant as*
C                          an absolute error bound.              *
C               IERR=-2    the nunber of maximally allowed steps *
C                          has been reached.                     *
C               IERR=-20   OVERFLOW, the program stops.          *
C               IERR=-30   the computed step size is too small.  *
C                          The program stops.                    *
C                                                                *
C                                                                *
C LOCAL VARIABLES:                                               *
C ================                                               *
C DSMALL  : DOUBLE PRECISION machine constant                    *
C DLARGE  : largest computer DOUBLE PRECISION number             *
C J       : loop variable                                        *
C FSAL    : (LOGICAL) variable indicating whether the method     *
C           FSAL (First Same As Last) is used by the RUNGE-KUTTA *
C           embedding formula                                    *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: DMPREC, COEFFI, HSTART, HULL, CIVPS     *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Volker KrÅger                                      *
C  Date     : 29.04.1993                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declarations
C
      EXTERNAL DES
      DOUBLE PRECISION Y(N),WORK1(4),WORK2(16,16)
      INTEGER IFLAG(3),IWORK(2)
      DOUBLE PRECISION X,BETA,ABSERR,RELERR,DSMALL,DLARGE
      LOGICAL FSAL
C
C check input parameters, if erroneous,
C return to calling program
C
      IERR=0
      IF(IFLAG(3) .LT. 0 .OR. IFLAG(3) .GT. 1) THEN
        IERR=2
        RETURN
      ELSEIF(IFLAG(3) .EQ. 0) THEN
C
C   determine the machine constant and the largest
C   representable number for the computer
C
        CALL DMPREC(DSMALL,DLARGE)
        WORK1(1)=100.0D0*DSMALL
        WORK1(2)=DLARGE/100.0D0
      ENDIF
C
C   determine the appropriate error code
C
      IF(IFLAG(1) .LT. 0 .OR. IFLAG(1) .GT. 22) IERR=IERR+1000000
      IF(IFLAG(2) .LT. 0 .OR. IFLAG(2) .GT. 1) IERR=IERR+100000
      IF(RELERR .LT. WORK1(1) .AND. ABSERR .LT. WORK1(1))
     +                                            IERR=IERR+10000
      IF(DABS(BETA-X) .LT. WORK1(1)) IERR=IERR+1000
      IF(DABS(X) .GT. WORK1(2) .OR. DABS(BETA) .GT. WORK1(2))
     +                                            IERR=IERR+100
      IF(N .LE. 0 .OR. N .GT. 12) IERR=IERR+10
      DO 100 J=1,N
         IF(DABS(Y(J)) .GT. WORK1(2)) THEN
           IERR=IERR+1
           RETURN
         ENDIF
100   CONTINUE
      IF(IERR .NE. 0) RETURN
C
C End input check
C
C Method FSAL is used by the RUNGE-KUTTA embedding formula
C
      FSAL=(IFLAG(1) .EQ. 1 .OR. IFLAG(1) .GE. 5 .AND.
     +      IFLAG(1) .LE. 8 .OR. IFLAG(1) .GE. 15 .AND.
     +      IFLAG(1) .LE. 17)
      IF(IFLAG(3) .EQ. 0) THEN
C
C Set the level of the embedding formula and initialize
C the maximally allowed number of integration steps
C
        IF(IFLAG(1) .EQ. 0) THEN
          IWORK(1)=3
          IWORK(2)=3330
        ELSEIF(IFLAG(1) .EQ. 1) THEN
          IWORK(1)=5
          IWORK(2)=2500
        ELSEIF(IFLAG(1) .GE. 2 .AND. IFLAG(1) .LE. 4) THEN
          IWORK(1)=6
          IWORK(2)=1670
        ELSEIF(IFLAG(1) .GE. 5 .AND. IFLAG(1) .LE. 8) THEN
          IWORK(1)=7
          IWORK(2)=1665
        ELSEIF(IFLAG(1) .GE. 9 .AND. IFLAG(1) .LE. 14) THEN
          IWORK(1)=8
          IWORK(2)=1250
        ELSEIF(IFLAG(1) .GE. 15 .AND. IFLAG(1) .LE. 17) THEN
          IWORK(1)=9
          IWORK(2)=1250
        ELSEIF(IFLAG(1) .EQ. 18) THEN
          IWORK(1)=10
          IWORK(2)=1000
        ELSEIF(IFLAG(1) .GE. 19 .AND. IFLAG(1) .LE. 21) THEN
          IWORK(1)=13
          IWORK(2)=770
        ELSEIF(IFLAG(1) .EQ. 22) THEN
          IWORK(1)=16
          IWORK(2)=625
        ENDIF
C
C Find the coefficients and the global error order of the
C low order method in the desired RUNGE-KUTTA pair
C
        IF(IFLAG(1) .EQ. 0) THEN
          WORK1(3)=2.0D0
        ELSE
          CALL COEFFI(IWORK(1),IFLAG(1),WORK2,WORK1(3))
        ENDIF
C
C Initialize first step size
C
        CALL HSTART(DES,N,X,BETA,Y,RELERR,ABSERR,WORK1(3),
     +              DSMALL,DLARGE,WORK1(4))
      ENDIF
C
C Integrate using HULL or CIVPS
C

      IF(IFLAG(2) .EQ. 0) THEN
        CALL HULL(X,WORK1(4),BETA,ABSERR,RELERR,N,FSAL,IWORK(1),
     +            DES,Y,WORK1(1),WORK1(2),WORK1(3),WORK2,IWORK(2),
     +            IERR)
      ELSE
        CALL CIVPS(X,WORK1(4),BETA,ABSERR,RELERR,N,FSAL,IWORK(1),
     +             DES,Y,WORK1(1),WORK1(2),WORK1(3),WORK2,IWORK(2),
     +             IERR)
      ENDIF
      RETURN
      END
