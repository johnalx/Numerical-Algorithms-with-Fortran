C[BA*)
C[LE*)
      SUBROUTINE K3RORI(USERF,PX,PY,QX,QY,RX,RY,N,WORK,CTRI,
C[IX{K3RORI}*)
     +                  DIVIAT,IERR,IUFCLL)
C
C*****************************************************************
C                                                                *
C Cubature over triangular regions using the summed 3-point      *
C formula and ROMBERG-RICHARDSON extrapolation:                  *
C[BE*)
C                                                                *
C Using the summed  3-point cubature formula of NEWTON-COTES     *
C we evaluate the integral of the Funktion USERF(X,Y) over the   *
C triangle P Q R by subdividing it into similar triangles.       *
C The number of cubature steps with halved sub-triangle edges is *
C designated by N.                                               *
C A RICHARDSON-extrapolation is used for an improved             *
C approximation of the integral.                                 *
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
C WORK    : 2-dimensional DOUBLE PRECISION array                 *
C           WORK(3,0:METHOD-1).                                  *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C CTRI    : DOUBLE PRECISION approximate value for the integral  *
C DIVIAT  : DOUBLE PRECISION error estimate                      *
C IERR    : error parameter: IERR=0 all is ok                    *
C                            IERR=1 N is incorrect               *
C                            IERR=2 the vertices P Q and R are   *
C                                   collinear                    *
C IUFCLL  : INTEGER, the number of function evaluations performed*
C                                                                *
C                                                                *
C INTERMEDIATE VARIABLES:                                        *
C =======================                                        *
C I       : loop variable                                        *
C IUFHLP  : auxiliary varialbe counting function evaluations     *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: K3GNEC3, RORIEX                         *
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
C declarations
C
      DOUBLE PRECISION WORK(0:N-1,2),PX,PY,QX,QY,RX,RY,CTRI,DIVIAT
      EXTERNAL USERF
C
C Initialize IUFCLL
C
      IUFCLL=0
C
C Check validity of N
C
      IF(N .LT. 2) THEN
        IERR=1
        RETURN
      ENDIF
C
C Perform N cubatures
C
      DO 10 I=0,N-1
         CALL K3NEC3(USERF,PX,PY,QX,QY,RX,RY,2**I,WORK(I,1),IERR,
     +               IUFHLP)
         IUFCLL=IUFCLL+IUFHLP
         IF (IERR .NE. 0) RETURN
10    CONTINUE
C
C Find an approximate integral value and an error estimate
C by using RICHARDSON-extrapolation
C
      CALL RORIEX(WORK(0,1),WORK(0,2),N,2,CTRI,DIVIAT)
C
C Return to calling program
C
      RETURN
      END
