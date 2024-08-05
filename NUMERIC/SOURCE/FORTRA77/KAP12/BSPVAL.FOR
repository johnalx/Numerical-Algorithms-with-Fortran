C[BA*)
C[LE*)
      SUBROUTINE BSPVAL(N,M,A,X,Y,XX,YY,VALUE,IERR)
C[IX{BSPVAL}*)
C
C*****************************************************************
C                                                                *
C  Used for computing the functional value of a bicubic spline   *
C  at the point (XX,YY).                                         *
C[BE*)
C  IERR=1, if the point lies outside the domain of definition    *
C  of the spline.                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: XYINTV                                  *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Eberhard Heyne                                     *
C  date     : 02.15.1983                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (KDIM=3,LDIM=3)
C
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM),X(0:N),Y(0:M)
      DIMENSION XIP(0:3),ETAP(0:3)
      DATA XIP(0)/1.0D0/
      DATA ETAP(0)/1.0D0/
C
C*  determine intervals I, J and the relative coordinates XI, ETA
C
      CALL XYINTV(N,M,X,Y,I,J,XI,ETA,XX,YY,IERR)
      IF(IERR .NE. 0) RETURN
      S=0.0D0
      DO 101 K=1,3
        XIP(K)=XIP(K-1)*XI
        ETAP(K)=ETAP(K-1)*ETA
  101 CONTINUE
      DO 103 K=0,3
        DO 102 L=0,3
          S=S+A(I,J,K,L)*XIP(K)*ETAP(L)
  102   CONTINUE
  103 CONTINUE
      VALUE=S
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE XYINTV(N,M,X,Y,I,J,XI,ETA,XX,YY,IERR)
C[IX{XYINTV}*)
C
C*****************************************************************
C                                                                *
C  Determines the interval,in which the point (XX,YY) lies.      *
C[BE*)
C  The index I is determined so that  X(I).LE.XX .AND.           *
C  XX.LE.X(I+1), while J is determined so that  Y(J).LE.YY       *
C  .AND. YY.LE.Y(I+1). We use XI =  XX-X(I) and ETA = YY-Y(J).   *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Eberhard Heyne                                     *
C  date     : 02.15.1983                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(0:N),Y(0:M)
      IERR=1
      IF(XX .LT. X(0)) RETURN
      IF(XX .GT. X(N)) RETURN
      IF(YY .LT. Y(0)) RETURN
      IF(YY .GT. Y(M)) RETURN
      IERR=0
      LU=0
      LO=N
  100 L=(LU+LO)/2
      IF(XX .LT. X(L)) THEN
        LO=L
        GOTO 100
      ELSE IF(XX .GT. X(L+1)) THEN
        LU=L
        GOTO 100
      ELSE
        I=L
        XI=XX-X(L)
      ENDIF
      LU=0
      LO=M
  101 L=(LU+LO)/2
      IF(YY .LT. Y(L)) THEN
        LO=L
        GOTO 101
      ELSE IF(YY .GT. Y(L+1)) THEN
        LU=L
        GOTO 101
      ELSE
        J=L
        ETA=YY-Y(L)
      ENDIF
      RETURN
      END
