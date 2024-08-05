C[BA*)
C[LE*)
      SUBROUTINE K4RORI(USERF,A,B,IP,C,D,IQ,N,CREC,DIVIAT,WORK,
C[IX{K4RORI}*)
     +                  IERR,IUFCLL)
C
C*****************************************************************
C                                                                *
C Cubature over rectangular regions using the summed ROMBERG-    *
C RICHARDSON method.                                             *
C[BE*)
C                                                                *
C By using a simplified summed trapezoidal rule, we can          *
C approximate the integral of the FUNCTION USERF(X,Y) over       *
C the rectangle [A,B] x [C,D] for the ROMBERG sequence of step   *
C sizes. After this we can obtain a better approximation of the  *
C integral by RICHARDSON extrapolation. For H = B - A or D - C,  *
C respectively, we use the step sizes:                           *
C   H * (1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256 ...)      *
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C USERF   : user defined FUNCTION USERF(X,Y), whose integral is  *
C           to be computed.                                      *
C           The FUNCTION USERF must be declared as EXTERNAL in   *
C           the calling program.                                 *
C           The FUNCTION should have the following form:         *
C                  DOUBLE PRECISION FUNCTION USERF(X,Y)          *
C                  DOUBLE PRECISION X,Y                          *
C                         .                                      *
C                         .                                      *
C                         .                                      *
C                  USERF=F(X,Y)                                  *
C                         .                                      *
C                         .                                      *
C                         .                                      *
C                  RETURN                                        *
C                  STOP                                          *
C                                                                *
C A       : DOUBLE PRECISION left endpoint in X-direction        *
C B       : DOUBLE PRECISION right endpoint in X-direction       *
C IP      : INTEGER, the number of intervals in X-direction      *
C C       : DOUBLE PRECISION lower endpoint in Y-direction       *
C D       : DOUBLE PRECISION upper endpoint in Y-direction       *
C IQ      : INTEGER, the number of intervals in Y-direction      *
C N       : INTEGER, number of summed trapezoidal cubatures,     *
C           N > 1                                                *
C WORK    : DOUBLE PRECISION storage vector WORK(0:METHOD+2)     *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C CREC    : DOUBLE PRECISION approximate value for the integral  *
C DIVIAT  : DOUBLE PRECISION error estimate                      *
C IERR    : error parameter: IERR=0 all is ok                    *
C                            IERR=1 number of intervals in       *
C                                   X-direction erroneous        *
C                            IERR=2 number of intervals in       *
C                                   Y-direction erroneous        *
C                            IERR=3 N <= 1                       *
C                            IERR=4 integrating over an interval *
C                                   of length zero.              *
C IUFCLL  : INTEGER, the number of functional evaluations        *
C           performed.                                           *
C                                                                *
C                                                                *
C LOCAL VARIABLE:                                                *
C =================                                              *
C J   : INTEGER loop variable                                    *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  required subroutines: K4ROST, RORIEX                          *
C                                                                *
C*****************************************************************
C                                                                *
C  Author  : Volker Krger                                       *
C  Date    : 06.12.1991                                          *
C  Source  : FORTRAN 77                                          *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declarations
C
      DOUBLE PRECISION WORK(0:N-1,2),A,B,C,D,CREC,DIVIAT
      EXTERNAL USERF
C
C Initialize IUFCLL
C
      IUFCLL=0
C
C Check input data
C
      IF(IP .LT. 1) THEN
        IERR=1
        RETURN
      ELSEIF(IQ .LT. 1) THEN
        IERR=2
        RETURN
      ELSEIF(N .LT. 2) THEN
        IERR=3
        RETURN
      ELSEIF(A .EQ. B .OR. C .EQ. D) THEN
        IERR=4
        RETURN
      ELSE
        IERR=0
      ENDIF
C
C Execute trapezoidal cubatures for the ROMBERG sequence of step sizes
C
      DO 10 J=0,N-1
         CALL K4ROST(USERF,A,B,IP,C,D,IQ,J,WORK(J,1),IUFCLL)
10    CONTINUE
C
C Approximate and estimate the error via RICHARDSON extrapolation
C
      CALL RORIEX(WORK(0,1),WORK(0,2),N,2,CREC,DIVIAT)
C
C Return to calling program
C
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE RORIEX(ROMB,WORK,N,M,VAL,ERREST)
C[IX{RORIEX}*)
C
C*****************************************************************
C                                                                *
C RICHARDSON extrapolation for a given ROMBERG sequence.         *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C ROMB    : DOUBLE PRECISION vector ROMB(0:N-1), containing      *
C           the ROMBERG sequence                                 *
C WORK    : DOUBLE PRECISION storage vector WORK(0:N-1)          *
C N       : INTEGER, the number of elements in ROMB              *
C M       : INTEGER, the order of the method                     *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C VAL     : DOUBLE PRECISION, final value of the extrapolation   *
C ERREST  : DOUBLE PRECISION error estimate for VAL              *
C                                                                *
C                                                                *
C LOCALE VARIABLES:                                              *
C =================                                              *
C K,J     : INTEGER loop counters                                *
C S       : DOUBLE PRECISION ] auxiliary                         *
C P       : DOUBLE PRECISION ]      variables                    *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Volker Krger                                      *
C  Date     : 06.12.1991                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declarations
C
      DOUBLE PRECISION ROMB(0:N-1),WORK(0:N-1),VAL,ERREST,P,S
C
C Store ROMBERG sequence in auxiliary vector
C
      DO 10 K=0,N-1
         WORK(K)=ROMB(K)
10    CONTINUE
C
C Initialize auxiliary variables S and P
C
      S=2.0D0**(DBLE(M))
      P=S
C
C RICHARDSON extrapolation
C
      DO 20 K=1,N-1
         DO 30 J=0,N-K-1
            WORK(J)=(P*WORK(J+1)-WORK(J))/(P-1.0D0)
30       CONTINUE
      P=P*S
20    CONTINUE
C
C Set up values for the extrapolation
C
      VAL=WORK(0)
C
C Error estimation
C
      ERREST=DABS(WORK(0)-WORK(1))
C
C Return to calling program
C
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE K4ROST(USERF,A,B,IP,C,D,IQ,J,CROST,IUFCLL)
C[IX{K4ROST}*)
C
C*****************************************************************
C                                                                *
C Cubature over rectangular regions using the summed trapezoidal *
C rule in order to compute the Jth element of a ROMBERG sequence.*
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C USERF   : user defined FUNCTION USERF(X,Y), whose integral is  *
C           to be computed.                                      *
C           The FUNCTION USERF must be declared as EXTERNAL in   *
C           the calling program.                                 *
C           The FUNCTION should have the following form:         *
C                  DOUBLE PRECISION FUNCTION USERF(X,Y)          *
C                  DOUBLE PRECISION X,Y                          *
C                         .                                      *
C                         .                                      *
C                         .                                      *
C                  USERF=F(X,Y)                                  *
C                         .                                      *
C                         .                                      *
C                         .                                      *
C                  RETURN                                        *
C                  STOP                                          *
C                                                                *
C A       : DOUBLE PRECISION left endpoint in X-direction        *
C B       : DOUBLE PRECISION right endpoint in X-direction       *
C IP      : INTEGER, the number of intervals in X-direction      *
C C       : DOUBLE PRECISION lower endpoint in Y-direction       *
C D       : DOUBLE PRECISION upper endpoint in Y-direction       *
C IQ      : INTEGER, the number of intervals in Y-direction      *
C J       : INTEGER, the index of the element in the ROMBERG     *
C           sequence                                             *
C IUFCLL  : INTEGER, the number of previous function evaluations *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C CROST   : DOUBLE PRECISION approximate value for the integral  *
C IUFCLL  : INTEGER, the number of function evaluations          *
C                                                                *
C                                                                *
C LOKALE VARIABLEN:                                              *
C =================                                              *
C I,K     : INTEGER loop counters                                *
C HAB     : DOUBLE PRECISION step size in X-direction            *
C HCD     : DOUBLE PRECISION step size in X-direction            *
C FAC     : DOUBLE PRECISION weight for the node                 *
C ABNUM   : DOUBLE PRECISION value for the number of intervals   *
C           in X-direction                                       *
C CDNUM   : DOUBLE PRECISION value for the number of intervals   *
C           in Y-direction                                       *
C DJ      : DOUBLE PRECISION varible used to determine the number*
C           of intervals                                         *
C DBLEI   : DOUBLE PRECISION value for I                         *
C IABNUM  : INTEGER value for ABNUM                              *
C ICDNUM  : INTEGER value for CDNUM                              *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Volker Krger                                      *
C  Date     : 06.12.1991                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declarations
C
      DOUBLE PRECISION A,B,C,D,CROST,USERF,HAB,HCD,FAC,ABNUM,
     +                 CDNUM,DJ,DBLEI
C
C Factor to determine number of intervals
C
      DJ=2.0D0**(DBLE(J))
C
C Number of intervals in both X- and Y-direction
C
      ABNUM=DJ*DBLE(IP)
      CDNUM=DJ*DBLE(IQ)
      IABNUM=INT(ABNUM)
      ICDNUM=INT(CDNUM)
C
C Step size in  X- and Y-direction
C
      HAB=(B-A)/ABNUM
      HCD=(D-C)/CDNUM

C
C Initialize CROST
C
         CROST=0.0D0
C
C Find a value for the integral via the summed trapezoidal rule
C
         DO 10 I=0,IABNUM
            DBLEI=DBLE(I)
            DO 20 K=0,ICDNUM
C
C Determine weights
C
               FAC=1.0D0
               IF(I .GT. 0 .AND. I .LT. IABNUM) FAC=2.0D0
               IF(K .GT. 0 .AND. K .LT. ICDNUM) FAC=2.0D0*FAC
               CROST=CROST+FAC*USERF(A+HAB*DBLEI,C+HCD*DBLE(K))
               IUFCLL=IUFCLL+1
20          CONTINUE
10       CONTINUE
         CROST=0.25D0*HAB*HCD*CROST
C
C Return to calling program
C
      RETURN
      END
