C[BA*)
C[LE*)
      SUBROUTINE HULL (X,H,BETA,ABSERR,RELERR,N,FSAL,M,DES,Y,EPS,
C[IX{HULL}*)
     +                 XZI,QG,COEFF,MAXSTP,IERR)
C
C*****************************************************************
C                                                                *
C This program solves a system of at most 12 ordinary            *
C differential equations of first order by using a RUNGE-KUTTA   *
C embedding formula over the interval of integration             *
C I= [X0,BETA]. The step size is automatically controlled using  *
C the method of HULL.                                            *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C X       : DOUBLE PRECISION initial value for the integration:  *
C           X=X0                                                 *
C H       : DOUBLE PRECISION initial step size                   *
C BETA    : DOUBLE PRECISION endpoint X=BETA at which we want to *
C           find the solution                                    *
C ABSERR  : DOUBLE PRECISION error bound for the absolute error  *
C           (ABSERR >= 0). If ABSERR=0, then only the relative   *
C           error is checked.                                    *
C RELERR  : DOUBLE PRECISION error bound for the relative error  *
C           (RELERR >= 0). If RELERR=0, then only the absolute   *
C           error is checked.                                    *
C N       : number of differential equations in the system,      *
C           or the size of Y:   0 < N < 13                       *
C FSAL    : (LOGICAL) variable indicating whether the method     *
C           FSAL (First Same As Last) is used by the RUNGE-KUTTA *
C           embedding formula                                    *
C M       : level of the embedding formula, also used for        *
C           dimensioning COEFF.                                  *
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
C EPS     : DOUBLE PRECISION 100 * machine constant              *
C XZI     : DOUBLE PRECISION largest representable number for    *
C           testing for OVERFLOW                                 *
C QG      : global error order of the low order RUNGE-KUTTA      *
C           method in use                                        *
C COEFF   : 2-dim. DOUBLE PRECISION array COEFF(1:16,1:M) which  *
C           contains the coefficients of the formula             *
C MAXSTP  : maximal number of allowed integration steps          *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C X       : DOUBLE PRECISION value for X, where the integration  *
C           has stopped (normally X=BETA)                        *
C H       : DOUBLE PRECISION last step size used                 *
C Y       : DOUBLE PRECISION solution vector Y(1:N) for X        *
C IERR    : error parameter:                                     *
C               IERR=0     all is ok                             *
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
C LSTSTP  : (LOGICAL)  LSTSTP=FALSE  continue integration        *
C                      LSTSTP=TRUE   stop integration            *
C ISTEP   : loop variable                                        *
C J       : loop variable                                        *
C YHILO   : 2-dim. DOUBLE PRECISION array YHILO(1:12,1:2)        *
C           see SUBROUTINE RKSTEP or STEP32                      *
C K       : 2-dim. DOUBLE PRECISION array K(1:12,1:16)           *
C           see SUBROUTINE RKSTEP or STEP32                      *
C YDIFF   : DOUBLE PRECISION auxiliary vector YDIFF(1:12)        *
C NOSTEP  : LOGICAL variable, see SUBROUTINE RKSTEP or STEP32    *
C FSALHP  : LOGICAL variable, auxiliary variable for FSAL        *
C DELTA   : DOUBLE PRECISION estimate of the local error         *
C EPSLON  : DOUBLE PRECISION tolerance for the maximal local     *
C           error                                                *
C TEMP    : DOUBLE PRECISION auxiliary variable, the last step   *
C           size                                                 *
C EXPO    : DOUBLE PRECISION variable                            *
C XDUMMY  : DOUBLE PRECISION variable                            *
C XEND    : DOUBLE PRECISION variable for checking the endpoint  *
C           of the interval                                      *
C JERR    : error parameter of the SUBROUTINE RKSTEP or STEP32   *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: RKSTEP, STEP32, VMNORM                  *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Volker Krger                                      *
C  Date     : 28.04.1993                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declarations
C
      EXTERNAL DES
      DOUBLE PRECISION Y(N),YHILO(12,2),K(12,16),YDIFF(12),
     +                 COEFF(16,M)
      DOUBLE PRECISION X,H,TEMP,EXPO,BETA,ABSERR,RELERR,VMNORM,
     +                 QG,DELTA,EPS,XZI,EPSLON,XDUMMY,XEND
      LOGICAL NOSTEP, FSAL, FSALHP, LSTSTP
C
C Initialize FSALHP, NOSTEP and LSTSTP
C
      FSALHP=.FALSE.
      NOSTEP=.FALSE.
      LSTSTP=.FALSE.
C
C Compute EXPO
C
      EXPO=1.0D0/(QG+1.0D0)
C
C Integrate:
C   the maximal number of steps is MAXSTP
C
      DO 20 ISTEP=1,MAXSTP
         IF(M .GT. 3) THEN
C
C   Integrate using the SUBROUTINE RKSTEP
C
           CALL RKSTEP(X,H,Y,N,M,K,DES,YHILO,COEFF,NOSTEP,FSALHP,
     +                 XZI,JERR)
         ELSE
C
C   Integrate using the SUBROUTINE STEP32
C
           CALL STEP32(X,H,Y,N,K,DES,YHILO,NOSTEP,XZI,JERR)
         ENDIF
C
C   after the first step we use the feature FSAL
C
         FSALHP=FSAL
C
C   If we encounter OVERFLOW:
C   return to calling program
C
         DO 200 J=1,N
            IF(JERR .EQ. 1 .OR. DABS(YHILO(J,1)) .GT. XZI .OR.
     +         DABS(YHILO(J,2)) .GT. XZI) THEN
              IERR=-20
              RETURN
            ENDIF
200      CONTINUE
C
C   Compute step size for next step
C
         DO 30 J=1,N
            YDIFF(J)=YHILO(J,1)-YHILO(J,2)
30       CONTINUE
         DELTA=VMNORM(YDIFF,N)
         EPSLON=ABSERR+RELERR*VMNORM(YHILO(1,1),N)
         IF(EPSLON .LT. EPS) THEN
           EPSLON=EPS
           IERR=-1
         ENDIF
C
C   DELTA exxceeds the tolerance EPSLON
C
         IF(DELTA .GT. EPSLON) THEN
           NOSTEP=.TRUE.
C
C   the new step size is limited to one fourth of the old step size.
C   If it then falls below the machine constant, we return to the
C   calling program
C
           H=H*DMAX1(2.5D-01,9.0D-01*(EPSLON/DELTA)**EXPO)
           IF(DABS(H) .LT. EPS) THEN
             IERR=-30
             RETURN
           ENDIF
C
C   Adjust LSTSTP if the last integration is repeated
C
           IF(LSTSTP) LSTSTP=.FALSE.
C
C   DELTA becomes smaller than the tolerance EPSLON
C
         ELSE
C
C   Reinitialize X and Y
C
           X=X+H
           DO 10 J=1,N
              Y(J)=YHILO(J,1)
10         CONTINUE
C
C   the last integration was successful, store the last step
C   size TEMP in H, go back to calling program
C
           IF(LSTSTP) THEN
             H=TEMP
             RETURN
C
C   prepare for a new integration
C
           ELSE
             NOSTEP=.FALSE.
C
C   the new step size is bounded by 4 times the old one
C
             IF(DELTA .LT. EPS .OR. (DELTA .LE. 1.0D+00 .AND.
     +          XZI*DELTA .LT. EPSLON)) THEN
C
C     Check prevent division by zero or OVERFLOW
C
               H=4.0D+00*H
             ELSE
               H=H*DMIN1(4.0D+00,9.0D-01*(EPSLON/DELTA)**EXPO)
C
C   If the step size becomes too small, return to calling program
C
               IF(DABS(H) .LT. EPS) THEN
                 IERR=-30
                 RETURN
               ENDIF
             ENDIF
C
C   Check whether the end of the interval has been reached
C
             XDUMMY=X+H
             XEND=BETA-0.1D0*H
             IF ((H .LT. 0.0D+00) .EQV. (XDUMMY .LT. XEND)) THEN
C
C   Adjust LSTSTP accordingly, store the last step size TEMP in H.
C   Compute H for the last integration step
C
               LSTSTP=.TRUE.
               TEMP=H
               H=BETA-X
             ENDIF
           ENDIF
         ENDIF
20    CONTINUE
C
C Integration is stopped after the maximally allowed
C number of integration steps
C
      IERR=-2
      RETURN
      END

