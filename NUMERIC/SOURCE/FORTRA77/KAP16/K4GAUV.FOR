C[BA*)
C[LE*)
      SUBROUTINE K4GAUV(USERF,X,NX,Y,NY,METHOD,MOLD,CREC,
C[IX{K4GAUV}*)
     +                  ESTDIV,DIVIAT,WORK,IERR,IUFCLL)
C
C*****************************************************************
C                                                                *
C Cubature over rectangles using NEWTON-COTES formulas for       *
C GAUSSIAN nodes:                                                *
C[BE*)
C                                                                *
C The FUNCTION USERF(X,Y) shall be integrated using summed       *
C GAUSSIAN formulas for the rectangle [A,B] x [C,D].             *
C The sub-rectangles to be used are given via the vectors        *
C X and Y.                                                       *
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
C X       : DOUBLE PRECISION vector X(0:NX) containing the X     *
C           partition of [A,B].                                  *
C           A = X(0) < X(1) < ... < X(NX) = B                    *
C NX      : number of intervals in X-direction, NX > 0           *
C Y       : DOUBLE PRECISION vector Y(0:NY) containing the Y     *
C           partition of [C,D].                                  *
C           C = Y(0) < Y(1) < ... < Y(NY) = D                    *
C NY      : number of intervals in Y-direction, NY > 0           *
C METHOD  : INTEGER indicating method used, 0 <= METHOD <= 7     *
C MOLD    : INTEGER, the number in METHOD at the previous call.  *
C           Upon first call we must have: MOLD different from    *
C           METHOD                                               *
C           If K4NECN is called repeatedly with METHOD=MOLD the  *
C           internal initializing of parameters is skipped.      *
C ESTDIV  : LOGICAL variable, indicates whether error estimate   *
C           is to be computed (ESTDIV=TRUE) or not (ESTDIV=FALSE)*
C WORK    : 2-dimensional DOUBLE PRECISION array                 *
C           WORK(3,0:METHOD-1). If METHOD=MOLD this array must   *
C           contain the initializing parameters for the method.  *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C MOLD    : INTEGER indicating method used                       *
C CREC    : DOUBLE PRECISION value for the integral              *
C WORK    : 2-dimensional DOUBLE PRECISION array WORK(2,0:METHOD)*
C           containing the constants for the method              *
C DIVIAT  : DOUBLE PRECISION error estimate                      *
C           If ESTDIV=TRUE the error is estimated by one extra   *
C           cubature for the halved step size.                   *
C IERR    : error parameter: IERR=0 all is ok                    *
C                            IERR=1 X-interval of length zero    *
C                            IERR=2 Y-interval of length zero    *
C                            IERR=3 Number of method erroneous   *
C                            IERR=4 NX < 1 or NY < 1             *
C IUFCLL  : INTEGER, the number of functional evaluations        *
C                                                                *
C                                                                *
C LOCAL VARIABLES:                                               *
C =================                                              *
C I, J, K : INTEGER, loop counters                               *
C II, JJ  : INTEGER, loop counters                               *
C I1, J1  : INTEGER, loop counters                               *
C KMAX    : INTEGER number of cubature passes                    *
C DBLEX   : DOUBLE PRECISION value for K                         *
C HX      : DOUBLE PRECISION step size in X-direction            *
C HY      : DOUBLE PRECISION step size in Y-direction            *
C HXM     : DOUBLE PRECISION mid-point of X-interval             *
C HYM     : DOUBLE PRECISION mid-point of Y-interval             *
C CRECH   : DOUBLE PRECISION variable used for error estimation  *
C FAC     : DOUBLE PRECISION variable used for CREC              *
C HELPF   : DOUBLE PRECISION variable used for CREC              *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines required: K4GINI                                  *
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
      DOUBLE PRECISION WORK(2,0:METHOD),X(0:NX),Y(0:NY),CREC,
     +                 CRECH,DIVIAT,HX,HY,HXM,HYM,DBLEK,FAC,HELPF,
     +                 USERF
C
C LOGICAL variable ESTDIV
C
      LOGICAL ESTDIV
C
C Initialize IUFCLL
C
      IUFCLL=0
C
C Validate input data
C
C   Length of X-intervals
C
      DO 10 I=1,NX
         IF (X(I) .LE. X(I-1)) THEN
            IERR=1
            RETURN
         ENDIF
10    CONTINUE
C
C   Length of Y-intervals
C
      DO 20 I=1,NY
         IF (Y(I) .LE. Y(I-1)) THEN
            IERR=2
            RETURN
         ENDIF
20    CONTINUE
C
C   Check number of method
C
      IF (METHOD .LT. 0 .OR. METHOD .GT. 7) THEN
         IERR=3
         RETURN
C
C   Check number of sub-intervals
C
      ELSEIF (NX .LE. 0 .OR. NY .LE. 0) THEN
         IERR=4
         RETURN
      ELSE
         IERR=0
      ENDIF
C
C If necessary, check initial values
C
      IF (METHOD .NE. MOLD) THEN
        CALL K4GINI(METHOD,WORK)
        MOLD=METHOD
      ENDIF
C
C
      IF (ESTDIV) THEN
        KMAX=2
      ELSE
        KMAX=1
      ENDIF
C
C Loop over necessary cubature runs
C
      DO 30 K=1,KMAX
C
C Change K
C
         DBLEK=DBLE(K)
C
C Initialize CREC
C
         CREC=0.0D0
C
C Find approximation for the integral
C
         DO 40 I=0,NX-1
C
C Find step size in X-direction
C
            HX=(X(I+1)-X(I))/(2.0D0*DBLEK)
            DO 50 I1=1,2*K-1,2
               HXM=X(I)+DBLE(I1)*HX
               DO 60 J=0,NY-1
C
C Find step size in Y-direction
C
                  HY=(Y(I+1)-Y(I))/(2.0D0*DBLEK)
                  DO 70 J1=1,2*K-1,2
                     HYM=Y(J)+DBLE(J1)*HY
                     DO 80 II=0,METHOD
                        DO 90 JJ=0,METHOD
                           FAC=HX*HY*WORK(2,II)*WORK(2,JJ)
                           HELPF=USERF(HXM+HX*WORK(1,II),
     F                                 HYM+HY*WORK(1,JJ))
                           IUFCLL=IUFCLL+1
                           CREC=CREC+FAC*HELPF
90                      CONTINUE
80                   CONTINUE
70                CONTINUE
60             CONTINUE
50          CONTINUE
40       CONTINUE
C
C When estimating the error, store first integral value
C
         IF (ESTDIV .AND. K .EQ. 1) CRECH=CREC
30    CONTINUE
C
C Error estimation
C
      IF (ESTDIV) DIVIAT=(CREC-CRECH)/3.0D0
C
C Return to calling program
C
      RETURN
      END
