C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 16.5}
C[  {Romberg Cubature for Rectangles and Triangles}
C[  {Romberg Cubature for Rectangles and Triangles}*)
C[LE*)
      SUBROUTINE K4BURI(USERF,A,B,IP,C,D,IQ,N,CREC,DIVIAT,WORK,
C[IX{K4BURI}*)
     +                  IERR,IUFCLL)
C
C*****************************************************************
C                                                                *
C Cubature over rectangular regions using the summed BULIRSCH-   *
C RICHARDSON method.                                             *
C[BE*)
C                                                                *
C We use a simplified summed two-dimensional trapezoidal rule    *
C for the BULIRSCH sequence of step sizes in order to find an    *
C approximation for the integral of the FUNCTION USERF(X,Y)      *
C over the rectangle [A,B] x [C,D]. Using RICHARDSON extrapola-  *
C tion we then compute an improved value for the integral.       *
C The step sizes follow the BULIRSCH sequences for the given     *
C interval length H and are:                                     *
C   H * (1/2, 1/3, 1/4, 1/6, 1/8, 1/12, 1/16, 1/24, 1/32 ...)    *
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
C A       : DOUBLE PRECISION left hand endpoint in X-direction   *
C B       : DOUBLE PRECISION right hand endpoint in X-direction  *
C IP      : INTEGER, the number of sub-intervals in X-direction  *
C C       : DOUBLE PRECISION lower endpoint in Y-direction       *
C D       : DOUBLE PRECISION upper endpoint in Y-direction       *
C IQ      : INTEGER, the number of sub-intervals in Y-direction  *
C N       : INTEGER, the number of summed trapezoidal cubatures  *
C           that are to be executed  ( N > 1 )                   *
C WORK    : 2-dimensional DOUBLE PRECISION array WORK(0:N-1,1:2) *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C CREC    : DOUBLE PRECISION value for the integral              *
C DIVIAT  : DOUBLE PRECISION error estimate                      *
C IERR    : error parameter: IERR=0 all is ok                    *
C                            IERR=1 number of intervals in       *
C                                   X-direction erroneous        *
C                            IERR=2 number of intervals in       *
C                                   Y-direction erroneous        *
C                            IERR=3 N <= 1                       *
C                            IERR=4 interval of integration has  *
C                                   length zero                  *
C IUFCLL  : INTEGER, the number of functional evaluations        *
C                                                                *
C                                                                *
C LOCAL VARIABLES:                                              *
C =================                                              *
C J       : INTEGER, loop counter                                *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: K4BUST, BURIEX                          *
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
C Perform cubature using trapezoidal rule for the
C BULIRSCH sequence of step sizes
C
      DO 10 J=0,N-1
         CALL K4BUST(USERF,A,B,IP,C,D,IQ,J,WORK(J,1),IUFCLL)
10    CONTINUE
C
C Find a better approximation and an error estimate
C for the integral using RICHARDSON extrapolation
C
      CALL BURIEX(WORK(0,1),WORK(0,2),N,2,CREC,DIVIAT)
C
C Return to calling program
C
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE K4BUST(USERF,A,B,IP,C,D,IQ,J,CBUST,IUFCLL)
C[IX{K4BUST}*)
C
C*****************************************************************
C                                                                *
C Cubature over rectangular regions using the summed trapezoidal *
C rule and computation of the Jth term of the BULIRSCH sequence. *
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
C A       : DOUBLE PRECISION left hand endpoint in X-direction   *
C B       : DOUBLE PRECISION right hand endpoint in X-direction  *
C IP      : INTEGER, the number of sub-intervals in X-direction  *
C C       : DOUBLE PRECISION lower endpoint in Y-direction       *
C D       : DOUBLE PRECISION upper endpoint in Y-direction       *
C IQ      : INTEGER, the number of sub-intervals in Y-direction  *
C J       : INTEGER, the index of the term of the BULIRSCH       *
C           sequence                                             *
C IUFCLL  : INTEGER, the number of functional evaluations of     *
C           previous calls                                       *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C CBUST   : DOUBLE PRECISION value for the integral              *
C IUFCLL  : INTEGER, the number of functional evaluations        *
C                                                                *
C                                                                *
C LOCAL VARIABLES:                                               *
C =================                                              *
C I,K     : INTEGERS, loop counters                              *
C HAB     : DOUBLE PRECISION step size in X-direction            *
C HCD     : DOUBLE PRECISION step size in Y-direction            *
C FAC     : DOUBLE PRECISION weights for the nodes               *
C ABNUM   : DOUBLE PRECISION number of intervals in X-direction  *
C CDNUM   : DOUBLE PRECISION number of intervals in Y-direction  *
C PQNUM   : DOUBLE PRECISION factor used to compute the number   *
C           of intervals                                         *
C DBLEI   : DOUBLE PRECISION value for I                         *
C IABNUM  : INTEGER value for ABNUM                              *
C ICDNUM  : INTEGER value for CDNUM                              *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines required: DENOM                                   *
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
      DOUBLE PRECISION A,B,C,D,CBUST,USERF,HAB,HCD,
     +                 FAC,DENOM,PQNUM,ABNUM,CDNUM,DBLEI
C
C Find the factor for the number of intervals
C
      PQNUM=DENOM(J)
C
C Number of intervals in X- and Y-directions
C
      ABNUM=PQNUM*DBLE(IP)
      CDNUM=PQNUM*DBLE(IQ)
      IABNUM=INT(ABNUM)
      ICDNUM=INT(CDNUM)
C
C Step sizes in X- and Y-directions
C
      HAB=(B-A)/ABNUM
      HCD=(D-C)/CDNUM
C
C Initialize CBUST
C
         CBUST=0.0D0
C
C Find approximate value for the integral by
C using the summed trapezoidal rule
C
         DO 10 I=0,IABNUM
            DBLEI=DBLE(I)
            DO 20 K=0,ICDNUM
C
C Determine weights for the nodes
C
               FAC=1.0D0
               IF(I .GT. 0 .AND. I .LT. IABNUM) FAC=2.0D0
               IF(K .GT. 0 .AND. K .LT. ICDNUM) FAC=2.0D0*FAC
               CBUST=CBUST+FAC*USERF(A+HAB*DBLEI,C+HCD*DBLE(K))
               IUFCLL=IUFCLL+1
20          CONTINUE
10       CONTINUE
         CBUST=0.25D0*HAB*HCD*CBUST
C
C Return to calling program
C
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BURIEX(BULI,WORK,N,M,VAL,ERREST)
C[IX{BURIEX}*)
C
C*****************************************************************
C                                                                *
C RICHARDSON extrapolation for a given BULIRSCH sequence         *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C BULI    : DOUBLE PRECISION vector BULI(0:N-1), containing the  *
C           BULIRSCH sequence                                    *
C WORK    : DOUBLE PRECISION vector WORK(0:N-1)                  *
C N       : INTEGER size of the vector BULI                      *
C M       : INTEGER order of the method                          *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C VAL     : DOUBLE PRECISION terminal value of the extrapolation *
C ERREST  : DOUBLE PRECISION error estimate for VAL              *
C                                                                *
C                                                                *
C LOCAL VARIABLES:                                               *
C ================                                               *
C K,J     : INTEGERS loop counters                               *
C P       : DOUBLE PRECISION auxiliary variable                  *
C DM      : DOUBLE PRECISION value for M                         *
C DK      : DOUBLE PRECISION value for K                         *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: DENOM                                   *
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
      DOUBLE PRECISION BULI(0:N-1),WORK(0:N-1),VAL,ERREST,P,
     +                 DENOM,DM,DK
C
C Change INTEGER value to DOUBLE PRECISION
C
      DM=DBLE(M)
C
C Store the BULIRSCH sequence in WORK
C
      DO 10 K=0,N-1
         WORK(K)=BULI(K)
10    CONTINUE
C
C RICHARDSON extrapolation
C
      DO 20 K=1,N-1
         DK=DBLE(K)
         DO 30 J=0,N-K-1
            P=(DENOM(J+K)/DENOM(J))**(DM*DK)
            WORK(J)=(P*WORK(J+1)-WORK(J))/(P-1.0D0)
30       CONTINUE
20    CONTINUE
C
C Store value of the extrapolation
C
      VAL=WORK(0)
C
C Error estimate
C
      ERREST=ABS(WORK(0)-WORK(1))
C
C Return to calling program
C
      RETURN
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION DENOM(J)
C[IX{DENOM}*)
C
C*****************************************************************
C                                                                *
C Determine the denominator for the step size in order to obtain *
C the Jth element of the BULIRSCH sequence                       *
C[BE*)
C                                                                *
C                                                                *
C Input PARAMETER:                                               *
C ================                                               *
C J       : INTEGER, the index of the element of the BULIRSCH    *
C           sequence                                             *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
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
C Determine denominator DENOM
C
      IF(J .EQ. 0) THEN
        DENOM=1.0D0
        RETURN
      ENDIF
      IF(MOD(J,2) .EQ. 0) THEN
C
C For even J
C
        DENOM=3.0D0*2.0D0**((DBLE(J)-2.0D0)*0.5D0)
      ELSE
C
C For odd J
C
        DENOM=2.0D0**((DBLE(J)+1.0D0)*0.5D0)
      ENDIF
C
C Return to calling program
C
      RETURN
      END
