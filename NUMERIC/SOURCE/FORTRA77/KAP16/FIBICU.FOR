C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 16.8}
C[  {Riemann Double Integrals using Bicubic Splines}
C[  {Riemann Double Integrals using Bicubic Splines}*)
C[LE*)
      SUBROUTINE FIBICU (N,M,A,X,Y,VALUE)
C[IX{FIBICU}*)
C
C*****************************************************************
C                                                                *
C  This subroutine determines a double integral of a spline      *
C  function over the rectangle [X(0), X(N)] x [Y(0), Y(M)],      *
C  which is the complete domain of the spline.                   *
C[BE*)
C  This program is faster than FIBIC2.                           *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: FIBIC1                                  *
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
C
      S=0.0D0
      DO 102 I=0,N-1
        DO 101 J=0,M-1
          S=S+FIBIC1(N,M,A,I,J,X(I+1)-X(I),Y(J+1)-Y(J))
  101   CONTINUE
  102 CONTINUE
      VALUE=S
      RETURN
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION FIBIC1 (N,M,A,I,J,XI,ETA)
C[IX{FIBIC1}*)
C
C*****************************************************************
C                                                                *
C  The function FIBIC1 determines a double integral of a spline  *
C  function over the rectangle [0, XI] x [0, ETA]. XI and ETA    *
C  are relative coordinates in the rectangle I, J.               *
C                                                                *
C----------------------------------------------------------------*
C[BE*)
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  autor     : Eberhard Heyne                                    *
C  datum     : 02.15.1983                                        *
C  sourcecode: FORTRAN 77                                        *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (KDIM=3,LDIM=3)
C
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM)
      DIMENSION XIP(0:4),ETAP(0:4)
      DATA XIP(0) /1.0D0/
      DATA ETAP(0)/1.0D0/
      S=0.0D0
      DO 101 K=1,4
        XIP(K)=XIP(K-1)*XI
        ETAP(K)=ETAP(K-1)*ETA
  101 CONTINUE
      DO 103 K=0,3
        DO 102 L=0,3
          S=S+A(I,J,K,L)*XIP(K+1)*ETAP(L+1)/DBLE((K+1)*(L+1))
  102   CONTINUE
  103 CONTINUE
      FIBIC1=S
      RETURN
      END
