C[BA*)
C[LE*)
      SUBROUTINE SHPLOK(X,Y,FX,FY,F,PPHI,W,R,N,DMUE,RR,PHI)
C[IX{SHPLOK}*)
C[LE*)
C********************************************************************
C[BE*)
C                                                                   *
C   Program name: SHPLOK                                            *
C                                                                   *
C********************************************************************
C                                                                   *
C[BA*)
C   This subroutine computes one functional value at (X,Y) for given*
C   nodes using the local Shepard method.                           *
C[BE*)
C                                                                   *
C********************************************************************
C                                                                   *
C                                                                   *
C   Input parameters:                                               *
C   =================                                               *
C   X    :  X value for which we want to interpolate the Z value    *
C   Y    :  Y value for which we want to interpolate the Z value    *
C   FX, FY, F   :  vectors ..(0:N) with X and Y coordinates of nodes*
C                  (FX,FY) and corresponding functional value F.    *
C   N    :  Index of last node                                      *
C   DMUE :  Exponent, 0 < DMUE < infinity, reasonable results can   *
C           be achieved for 2 < DMUE < 6. If on input DMUE <= 0, we *
C           set DMUE = 2 internally.                                *
C   RR   :  Radius around (X,Y) inside which all nodes are used to  *
C           interpolate at (X,Y).                                   *
C                                                                   *
C                                                                   *
C   AUX VECTORS:                                                    *
C   ============                                                    *
C   W, R, PPHI :  vectors ..(0:N)                                   *
C                                                                   *
C                                                                   *
C   Output parameters:                                              *
C   ==================                                              *
C   PHI  : Interpolated Z value at (X,Y)                            *
C                                                                   *
C                                                                   *
C********************************************************************
C                                                                   *
C   Required subroutines: none                                      *
C                                                                   *
C********************************************************************
C                                                                   *
C   Author      : Bjoern Terwege                                    *
C   Date        : 6.12.1995                                         *
C   Source code : FORTRAN 77                                        *
C                                                                   *
C[BA*)
C********************************************************************
C[BE*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(0:N),W(0:N),FX(0:N),FY(0:N),F(0:N),PPHI(0:N)
      IF(DMUE.LT.0.) DMUE=2.
      DO 20 I=0,N
          R(I)=DSQRT(((X-FX(I))*(X-FX(I)))+
     F         ((Y-FY(I))*(Y-FY(I))))
          IF (R(I).EQ.0) THEN
             PHI=F(I)
             GOTO 111
          ENDIF
20     CONTINUE
C
C  Compute the PPHI(I)
C
      DO 40 I=0,N
          IF(R(I).GE.RR) THEN
              PPHI(I)=0
          ELSE
              PPHI(I)=(RR/R(I))-1
          ENDIF
40    CONTINUE
C
C  Compute the denominators of the weights
C
      SUM=0
      DO 50 J=0,N
           IF(PPHI(J).EQ.0) GOTO 50
                SUM = SUM+1/(PPHI(J)**DMUE)
50    CONTINUE
C
C  Compute the weights
C
      DO 60 J=0,N
          IF(PPHI(J).EQ.0) THEN
               W(J)=0
               GOTO 60
          ELSE
               W(J)=1/(PPHI(J)**DMUE)*1/SUM
          ENDIF
60    CONTINUE
C
C  Compute the approximate function value Z at (X,Y)
C
      PHI=0
      DO 70 I=0,N
          PHI=PHI+W(I)*F(I)
70    CONTINUE
111   CONTINUE
      RETURN
      END
