C[BA*)
C[KA{F 16}{Numerical Cubature}{Numerical Cubature}*)
C[FE{F 16.3}
C[  {Newton-Cotes Cubature Formulas for Rectangles}
C[  {Newton-Cotes Cubature Formulas for Rectangles}*)
C[LE*)
      SUBROUTINE K4NECN(USERF,A,B,IP,C,D,IQ,METHOD,MOLD,CREC,
C[IX{K4NECN}*)
     +                  ESTDIV,DIVIAT,WORK,IERR,IUFCLL)
C
C*****************************************************************
C                                                                *
C Cubature for rectangular regions using the NEWTON-COTES        *
C formulas:                                                      *
C[BE*)
C                                                                *
C The FUNCTION USERF(X,Y) is integrated over the rectangle       *
C [A,B] x [C,D] using the summed NEWTON-COTES formulas.          *
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
C METHOD  : INTEGER designating the method used:                 *
C                     1: trapezoidal rule                        *
C                     2: SIMPSON rule                            *
C                     3: 3/8-rule                                *
C                     4: 4/90-rule                               *
C                     5: 5/288-rule                              *
C                     6: 6/840-rule                              *
C                     7: 7/17280-rule                            *
C MOLD    : INTEGER, the number of the method used at a previous *
C           call of this subroutine. On first call we must have  *
C           that MOLD does not equal METHOD !                    *
C           In a subsequent call of K4NECN with METHOD=MOLD      *
C           the internal initializing of parameters is skipped.  *
C ESTDIV  : LOGICAL variable; If ESTDIV=TRUE we compute an error *
C           estimate, if ESTDIV=FALSE we do not.                 *
C WORK    : DOUBLE PRECISION vector WORK(0:METHOD+2)             *
C           If METHOD=MOLD, WORK must contain the constants      *
C           initialized for the proper method.                   *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C MOLD    : INTEGER, the number of the method used               *
C CREC    : DOUBLE PRECISION value for the integral              *
C DIVIAT  : DOUBLE PRECISION error estimate                      *
C           If ESTDIV=TRUE, we perform an additional cubature    *
C           with halved stepsize for error estimation.           *
C WORK    : DOUBLE PRECISION vector WORK(0:METHOD+2),            *
C           contains the constants for the method used.          *
C IERR    : error parameter: IERR=0 all is ok                    *
C                            IERR=1 number of intervals in       *
C                                   X-direction erroneous        *
C                            IERR=2 number of intervals in       *
C                                   Y-direction erroneous        *
C                            IERR=3 method number erroneous      *
C                            IERR=4 integrating over an interval *
C                                   of length zero.              *
C IUFCLL  : INTEGER, the number of functional evaluations        *
C           performed.                                           *
C                                                                *
C                                                                *
C LOCAL VARIABLES:                                               *
C =================                                              *
C I,J,K   : INTEGER loop variables                               *
C KMAX    : INTEGER, the number of cubature passes desired       *
C IPX     : INTEGER, number of intervals in X-direction          *
C IPY     : INTEGER, number of intervals in Y-direction          *
C DBLEI   : DOUBLE PRECISION value for I                         *
C HX      : DOUBLE PRECISION step size in X-direction            *
C HY      : DOUBLE PRECISION step size in Y-direction            *
C CRECH   : DOUBLE PRECISION auxiliary variable                  *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  required subroutines: K4INIT, GRIDOT                          *
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
      DOUBLE PRECISION WORK(0:METHOD+2),A,B,C,D,CREC,CRECH,DIVIAT,
     +                 HX,HY,DBLEI,GRIDOT,USERF
C
C Initialize LOGICAL variable ESTDIV
C
      LOGICAL ESTDIV
C
C Initialize IUFCLL
C
      IUFCLL=0
C
C Check input data
C
      IF (IP .LT. 1) THEN
         IERR=1
         RETURN
      ELSEIF (IQ .LT. 1) THEN
         IERR=2
         RETURN
      ELSEIF (METHOD .LT. 1 .OR. METHOD .GT. 7) THEN
         IERR=3
         RETURN
      ELSEIF (A .EQ. B .OR. C .EQ. D) THEN
         IERR=4
         RETURN
      ELSE
         IERR=0
      ENDIF
C
C Initialize as needed
C
      IF (METHOD .NE. MOLD) THEN
        CALL K4INIT(METHOD,WORK)
        MOLD=METHOD
      ENDIF
C
C Set maximal number of cubature passes
C
      IF (ESTDIV) THEN
        KMAX=2
      ELSE
        KMAX=1
      ENDIF
C
C Determine number of X and Y sub-intervals
C
      DO 10 K=1,KMAX
         IPX=K*IP*METHOD
         IQY=K*IQ*METHOD
C
C Initialize CREC
C
         CREC=0.0D0
C
C Find step sizes in both X- and Y-directions
C
         HX=(B-A)/DBLE(IPX)
         HY=(D-C)/DBLE(IQY)
C
C Compute an approximate value of the integral
C
         DO 20 I=0,IPX
            DBLEI=DBLE(I)
            DO 30 J=0,IQY
               CREC=CREC+GRIDOT(I,J,WORK,METHOD,IPX,IQY)*
     +              USERF(A+DBLEI*HX,C+DBLE(J)*HY)
               IUFCLL=IUFCLL+1
30          CONTINUE
20       CONTINUE
         CREC=CREC*HX*HY*WORK(METHOD+1)
C
C In case error estimate is desired, store first value for the integral
C
         IF (ESTDIV .AND. K .EQ. 1) CRECH=CREC
10    CONTINUE
C
C Estimate the error
C
      IF (ESTDIV) DIVIAT=(CREC-CRECH)/(2.0D0**WORK(METHOD+2)-1.0D0)
C
C Return to calling program
C
      RETURN
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION GRIDOT(I,J,WORK,METHOD,IPX,IQY)
C[IX{GRIDOT}*)
C
C*****************************************************************
C                                                                *
C FUNCTION that determines the weights at the nodes.             *
C[BE*)
C                                                                *
C In summed Newton-Cotes cubature the functional values are      *
C given different weights that depend on their position on the   *
C boundary, center or at the join of two rectangles.             *
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C I       : INTEGER, the number of the node in X-direction       *
C J       : INTEGER, the number of the node in Y-direction       *
C WORK    : DOUBLE PRECISION vector WORK(0:METHOD), containing   *
C           the constants for the method                         *
C METHOD  : INTEGER designating the method chosen:               *
C                     1: trapezoidal rule                        *
C                     2: SIMPSON rule                            *
C                     3: 3/8-rule                                *
C                     4: 4/90-rule                               *
C                     5: 5/288-rule                              *
C                     6: 6/840-rule                              *
C                     7: 7/17280-rule                            *
C IPX     : INTEGER number of intervals used in X-direction      *
C IPY     : INTEGER number of intervals used in Y-direction      *
C                                                                *
C                                                                *
C LOCAL VARIABLE:                                                *
C ================                                               *
C K       : INTEGER auxiliary variable                           *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Volker KrÅger                                      *
C  Date     : 06.121991                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declarations
C
      DOUBLE PRECISION WORK(0:METHOD)
C
C Determine the weights for the nodes
C
      K=MOD(I,METHOD)
      GRIDOT=WORK(K)
      IF (K.EQ.0 .AND. I.GT.0 .AND. I.LT.IPX) GRIDOT=2.0D0*GRIDOT
      K=MOD(J,METHOD)
      GRIDOT=GRIDOT*WORK(K)
      IF (K.EQ.0 .AND. J.GT.0 .AND. J.LT.IQY) GRIDOT=2.0D0*GRIDOT
C
C Return to calling program
C
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE K4INIT(METHOD,WORK)
C[IX{K4INIT}*)
C
C*****************************************************************
C                                                                *
C Subroutine that initializes the constants for the various      *
C methods.                                                       *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C METHOD  : INTEGER, the number designating the method:          *
C                     1: trapezoidal rule                        *
C                     2: SIMPSON rule                            *
C                     3: 3/8 rule                                *
C                     4: 4/90 rule                               *
C                     5: 5/288 rule                              *
C                     6: 6/840 rule                              *
C                     7: 7/17280 rule                            *
C                                                                *
C                                                                *
C OUTPUT PARAMETER:                                              *
C =================                                              *
C WORK    : DOUBLE PRECISION vector WORK(0:METHOD+2),            *
C           containing the constants for the method              *
C                                                                *
C                                                                *
C LOCAL VARIABLES:                                               *
C ================                                               *
C I       : INTEGER loop counter                                 *
C J       : INTEGER auxiliary variable                           *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Volker KrÅger                                      *
C  Date     : 06.12.1991                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declarations
C
      DOUBLE PRECISION WORK(0:METHOD+2)
C
C Initialize upper half of WORK depending on method
C
      IF (METHOD .EQ. 1) THEN
         WORK(0)=1.0D0
         WORK(METHOD+1)=2.0D0
         WORK(METHOD+2)=2.0D0
         J=1
      ELSEIF (METHOD .EQ. 2) THEN
         WORK(0)=1.0D0
         WORK(1)=4.0D0
         WORK(METHOD+1)=6.0D0
         WORK(METHOD+2)=4.0D0
         J=2
      ELSEIF (METHOD .EQ. 3) THEN
         WORK(0)=1.0D0
         WORK(1)=3.0D0
         WORK(METHOD+1)=8.0D0
         WORK(METHOD+2)=4.0D0
         J=2
      ELSEIF (METHOD .EQ. 4) THEN
         WORK(0)=7.0D0
         WORK(1)=32.0D0
         WORK(2)=12.0D0
         WORK(METHOD+1)=90.0D0
         WORK(METHOD+2)=6.0D0
         J=3
      ELSEIF (METHOD .EQ. 5) THEN
         WORK(0)=19.0D0
         WORK(1)=75.0D0
         WORK(2)=50.0D0
         WORK(METHOD+1)=288.0D0
         WORK(METHOD+2)=6.0D0
         J=3
      ELSEIF (METHOD .EQ. 6) THEN
         WORK(0)=41.0D0
         WORK(1)=216.0D0
         WORK(2)=27.0D0
         WORK(3)=272.0D0
         WORK(METHOD+1)=840.0D0
         WORK(METHOD+2)=8.0D0
         J=4
      ELSEIF (METHOD .EQ. 7) THEN
         WORK(0)=751.0D0
         WORK(1)=3577.0D0
         WORK(2)=1323.0D0
         WORK(3)=2989.0D0
         WORK(METHOD+1)=17280.0D0
         WORK(METHOD+2)=8.0D0
         J=4
      ENDIF
C
C Determine lower half symmetrically
C
      DO 10 I=J,METHOD
         WORK(I)=WORK(METHOD-I)
10    CONTINUE
C
C Determine the multiplication factors for the summed values
C in the cubature formula
C
      WORK(METHOD+1)=DBLE(METHOD)/WORK(METHOD+1)
      WORK(METHOD+1)=WORK(METHOD+1)*WORK(METHOD+1)
C
C Return to calling program
C
      RETURN
      END
