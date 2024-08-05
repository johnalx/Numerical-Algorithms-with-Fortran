C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 16.4}
C[  {Newton-Cotes Cubature Formulas for Triangles}
C[  {Newton-Cotes Cubature Formulas for Triangles}*)
C[LE*)
      SUBROUTINE K3NEC3(USERF,PX,PY,QX,QY,RX,RY,N,CTRI,IERR,
C[IX{K3NEC3}*)
     +                  IUFCLL)
C
C*****************************************************************
C                                                                *
C Cubature for triangular region using the three point NEWTON-   *
C COTES formulas:                                                *
C[BE*)
C                                                                *
C The FUNCTION USERF(X,Y) is integrated over the triangle PQR    *
C according to the summed NEWTON-COTES formulas using sub-tri-   *
C angles.                                                        *
C The dimensions of these sub-triangles are one Nth of those of  *
C the original triangle PQR.                                     *
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
C PX      : DOUBLE PRECISION X-coordinate of the vertex P        *
C PY      : DOUBLE PRECISION Y-coordinate of the vertex P        *
C QX      : DOUBLE PRECISION X-coordinate of the vertex Q        *
C QY      : DOUBLE PRECISION Y-coordinate of the vertex Q        *
C RX      : DOUBLE PRECISION X-coordinate of the vertex R        *
C RY      : DOUBLE PRECISION Y-coordinate of the vertex R        *
C N       : INTEGER, counting the number of sub-triangles formed *
C           along one edge of the triangle.                      *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C CTRI    : DOUBLE PRECISION approximate value for the integral  *
C IERR    : error parameter: IERR=0 all is ok                    *
C                            IERR=1 N is incorrect               *
C                            IERR=2 the vertices P Q and R are   *
C                                   collinear                    *
C IUFCLL  : INTEGER, the number of function evaluations performed*
C                                                                *
C                                                                *
C INTERMEDIATE VARIABLES:                                        *
C =======================                                        *
C I,J     : loop variables                                       *
C II,JJ   : loop variables                                       *
C DBLEN   : DOUBLE PRECISION version of N                        *
C DBLEI   : DOUBLE PRECISION version of I                        *
C DBLEJ   : DOUBLE PRECISION version of J                        *
C AREA    : DOUBLE PRECISION to check collinearity               *
C EPS     : DOUBLE PRECISION bound for collinearity check        *
C HPQX    : DOUBLE PRECISION ]   vectoriel representation of     *
C HPQY    : DOUBLE PRECISION ]   the steps taken along the edge  *
C HPRX    : DOUBLE PRECISION ]   PQ or PR, respectively          *
C HPRY    : DOUBLE PRECISION ]                                   *
C FAC     : DOUBLE PRECISION weight for the node                 *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author  : Volker KrÅger                                       *
C  Date    : 06.12.1991                                          *
C  Source  : FORTRAN 77                                          *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declarations
C
      DOUBLE PRECISION PX,PY,QX,QY,RX,RY,CTRI,AREA,EPS,HPQX,HPQY,
     +                 HPRX,HPRY,FAC,DBLEM,DBLEJ,DBLEI,USERF
C
C Initialize bound for collinearity test
C
      EPS=1.0D-06
C
C Initialize IUFCLL
C
      IUFCLL=0
C
C Check validity of N
C
      IF (N .LT. 1) THEN
         IERR=1
         RETURN
      ENDIF
C
C Test for collinearity
C
      AREA=PX*QY+QX*RY+RX*PY-PX*RY-QX*PY-RX*QY
      IF (AREA .LT. EPS) THEN
         IERR=2
         RETURN
      ELSE
         IERR=0
      ENDIF
C
C Number of halved triangular edges
C
      M=2*N
      DBLEM=DBLE(M)
C
C Vectoriel representation of the step sizes
C
      HPQX=(QX-PX)/DBLEM
      HPQY=(QY-PY)/DBLEM
      HPRX=(RX-PX)/DBLEM
      HPRY=(RY-PY)/DBLEM

C
C Initialize CTRI
C
         CTRI=0.0D0
C
C Compute approximate value for integral
C
         DO 10 J=0,M-1
            DBLEJ=DBLE(J)
            DO 20 I=0,M-J
               DBLEI=DBLE(I)
C
C Determine weights for the nodes
C
               IF (MOD(I,2) .NE. 0 .OR. MOD(J,2) .NE. 0) THEN
                  IF (I .EQ. 0 .OR. J .EQ. 0 .OR. I .EQ. M-J) THEN
                     FAC=1.0D0
                  ELSE
                     FAC=2.0D0
                  ENDIF
                  CTRI=CTRI+FAC*USERF(PX+HPQX*DBLEI+HPRX*
     F                 DBLEJ,PY+HPQY*DBLEI+HPRY*DBLEJ)
               IUFCLL=IUFCLL+1
               ENDIF
20          CONTINUE
10       CONTINUE
         CTRI=CTRI*AREA/(6.0D0*DBLE(N)**2.0D0)
C
C Return to the calling program
C
      RETURN
      END
