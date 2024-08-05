C[BA*)
C[LE*)
      SUBROUTINE BICSP2(N,M,A,X,Y,F,IERR)
C[IX{BICSP2}*)
C[BE*)
C
C*****************************************************************
C                                                                *
C  Determination of bicubic splines.                             *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N    : number of X-intervals                                  *
C  M    : number of Y-intervals                                  *
C  A    : 4-dimensional array A(0:N,0:M,0:KDIM,0:LDIM); contains *
C         the spline coefficients. On call, A(I,J,0,0) must con- *
C         tain the functional values U(I,J).                     *
C         BICSP1 determines all other A(I,J,K,L) for I=0 to N-1  *
C         and J=0 to M-1. Elements A(N,M,K,L), that are not      *
C         assigned any values on call, remain unassigned.        *
C  X    : (N+1)-vector X(0:N) contains the endpoints of the      *
C         X-intervals                                            *
C  Y    : (N+1)-vector Y(0:M) contains the endpoints of the      *
C         Y-intervals                                            *
C  F    : auxiliary vector F(1:9*MAX(N,M)-8)                     *
C  IERR : is initially set to 0. Will be changed to be nonzero,  *
C         if the algorithm detects an error. If errors occur,    *
C         the program does not complete the computations.        *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: BIC2S1, BIC2S2, BIC2S3, BICSP1          *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Eberhard Heyne                                     *
C  date     : 02.15.1983                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (KDIM=3,LDIM=3)
C
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM),X(0:N),Y(0:M),F(*)
C
C*  step 1
C
      CALL BIC2S1(N,M,A,X,F,JERR)
      IERR=JERR
      IF(JERR .NE. 0) RETURN
C
C*  step 2
C
      CALL BIC2S2(N,M,A,Y,F,JERR)
      IERR=JERR+1
      IF(JERR .NE. 0) RETURN
C
C*  step 3
C
      CALL BIC2S3(N,M,A,X,F,JERR)
      IERR=JERR+4
      IF(JERR .NE. 0) RETURN
C
C*  step 4 to 12
C
      CALL BICSP1(N,M,A,X,Y,F,JERR)
      IERR=JERR+4
      IF(JERR .NE. 0) RETURN
C
C*  algorithm has run successfully
C
      IERR=0
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BIC2S1(N,M,A,X,H,IERR)
C[IX{BIC2S1}*)
C[BE*)
C
C*****************************************************************
C                                                                *
C  step 1:                                                       *
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
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (KDIM=3,LDIM=3)
C
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM),X(0:N),H(0:N-1)
C
C*  determine H(0),H(1),H(N-2),H(N-1), test strict monotonicity
C
      DO 102 K=0,1
        DO 101 L=0,N-2,N-2
          I=K+L
          H(I)=X(I+1)-X(I)
          IF(H(I) .LE. 0.0D0) GOTO 900
  101   CONTINUE
  102 CONTINUE
C
C*  determine boundary values and store them in A
C
      DO 103 J=0,M
C
C*  converted formulas
C
        A(0,J,1,0)=
     +     (A(1,J,0,0)-A(0,J,0,0))*(1.0D0/H(0)+0.5D0/(H(0)+H(1)))
     +    -(A(2,J,0,0)-A(1,J,0,0))*H(0)/(H(1)*2.0D0*(H(0)+H(1)))
        A(N,J,1,0)=(A(N,J,0,0)-A(N-1,J,0,0))*
     +     (0.5D0/(H(N-2)+H(N-1))+1.0D0/H(N-1))
     +    -(A(N-1,J,0,0)-A(N-2,J,0,0))*H(N-1)/(2.0D0*(H(N-2)+
     +     H(N-1))*H(N-2))
C
  103 CONTINUE
      IERR=0
      RETURN
  900 IERR=2
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BIC2S2(N,M,A,Y,H,IERR)
C[IX{BIC2S2}*)
C[BE*)
C
C*****************************************************************
C                                                                *
C  step 2:                                                       *
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
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (KDIM=3,LDIM=3)
C
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM),Y(0:M),H(0:M-1)
C
C*  determine H(0),H(1),H(M-2),H(M-1), test strict monotonicity
C
      DO 102 K=0,1
        DO 101 L=0,M-2,M-2
          J=K+L
          H(J)=Y(J+1)-Y(J)
          IF(H(J) .LE. 0.0D0) GOTO 900
  101   CONTINUE
  102 CONTINUE
C
C*  determine boundary values and store them in A
C
      DO 103 I=0,N
C
C*  converted formulas
C
        A(I,0,0,1)=
     +     (A(I,1,0,0)-A(I,0,0,0))*(1.0D0/H(0)+0.5D0/(H(0)+H(1)))
     +    -(A(I,2,0,0)-A(I,1,0,0))*H(0)/(H(1)*2.0D0*(H(0)+H(1)))
        A(I,M,0,1)=(A(I,M,0,0)-A(I,M-1,0,0))*
     +     (0.5D0/(H(M-2)+H(M-1))+1.0D0/H(M-1))
     +    -(A(I,M-1,0,0)-A(I,M-2,0,0))*H(M-1)/(2.0D0*(H(M-2)+
     +     H(M-1))*H(M-2))
C
  103 CONTINUE
      IERR=0
      RETURN
  900 IERR=2
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BIC2S3(N,M,A,X,H,IERR)
C[IX{BIC2S3}*)
C[BE*)
C
C*****************************************************************
C                                                                *
C  step 3:                                                       *
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
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (KDIM=3,LDIM=3)
C
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM),X(0:N),H(0:N-1)
C
C*  determine H(0),H(1),H(N-2),H(N-1), test strict monotonicity
C
      DO 102 K=0,1
        DO 101 L=0,N-2,N-2
          I=K+L
          H(I)=X(I+1)-X(I)
          IF(H(I) .LE. 0.0D0) GOTO 900
  101   CONTINUE
  102 CONTINUE
C
C*  determine boundary values and store them in A
C
      DO 103 J=0,M,M
C
C*  converted formulas
C
        A(0,J,1,1)=
     +     (A(1,J,0,1)-A(0,J,0,1))*(1.0D0/H(0)+0.5D0/(H(0)+H(1)))
     +    -(A(2,J,0,1)-A(1,J,0,1))*H(0)/(H(1)*2.0D0*(H(0)+H(1)))
        A(N,J,1,1)=(A(N,J,0,1)-A(N-1,J,0,1))*
     +     (0.5D0/(H(N-2)+H(N-1))+1.0D0/H(N-1))
     +    -(A(N-1,J,0,1)-A(N-2,J,0,1))*H(N-1)/(2.0D0*(H(N-2)+
     +     H(N-1))*H(N-2))
C
  103 CONTINUE
      IERR=0
      RETURN
  900 IERR=2
      RETURN
      END
