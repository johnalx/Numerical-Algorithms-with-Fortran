C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 9.8.2}{Shepard Interpolation}{Shepard Interpolation}*)
C[LE*)
      SUBROUTINE SHPGLO(X,Y,FX,FY,F,W,R,N,DMUE,PHI)
C[IX{SHPGLO}*)
C[LE*)
C********************************************************************
C[BE*)
C                                                                   *
C   Program name: SHPGLO                                            *
C                                                                   *
C********************************************************************
C                                                                   *
C[BA*)
C   This subroutine computes one functional value at (X,Y) for given*
C   nodes using the global Shepard method.                          *
C[BE*)
C   The exponent dmue must be specified externally.                 *
C   DMUE should be chosen to lie betwen 2 and 6.                    *
C                                                                   *
C********************************************************************
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
C                                                                   *
C   AUX VECTORS:                                                    *
C   ============                                                    *
C   W, R :  vectors ..(0:N)                                         *
C                                                                   *
C                                                                   *
C   Output parameters:                                              *
C   ==================                                              *
C   PHI  : Interpolated Z value at (X,Y)                            *
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
C 
C[BA*)
C********************************************************************
C[BE*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(0:N),W(0:N),FX(0:N),FY(0:N),F(0:N)
C
C  Check that DMUE > 0, otherwise set DMUE = 2
C
      IF(DMUE.LE.0.) DMUE=2
C
C  Compute the R(I)
C
      DO 20 I=0,N
          R(I)=DSQRT(((X-FX(I))*(X-FX(I)))
     F         +((Y-FY(I))*(Y-FY(I))))
          IF (R(I).EQ.0) THEN
               PHI=F(I)
               GOTO 333
          ENDIF
20    CONTINUE
C
C  Compute the sum in the denominator of W(I)
C
      DO 30 J=0,N
           SUM=0
           DO 40 I=0,N
                SUM=SUM+1/(R(I)**DMUE)
40         CONTINUE
C       Berechnung der Gewichte W(I)
C
           W(J)=1/(R(J)**DMUE)*1/SUM
30    CONTINUE
      PHI=0
C
C  Compute the approximate function value at (X,Y)
C
      DO 50 J=0,N
           PHI=PHI+W(J)*F(J)
50    CONTINUE
333   CONTINUE
      RETURN
      END
