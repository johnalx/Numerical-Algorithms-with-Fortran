C[BA*)
C[KA{F 12}
C[  {Two-Dim., B\'ezier, Surface, B-Splines}
C[  {Two-Dimensional Splines, Surface Splines,
C[   B\'ezier Splines, B-Splines}*)
C[FE{F 12.1}
C[  {Interpolating Two-Dimensional Cubic Splines}
C[  {Interpolating Two-Dimensional Cubic Splines for
C[   Constructing Smooth Surfaces}*)
C[LE*)
      SUBROUTINE BICSP1(N,M,A,X,Y,F,IERR)
C[IX{BICSP1}*)
C
C*****************************************************************
C                                                                *
C  Determination of bicubic splines.                             *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N    : number of X-intervals                                  *
C  M    : number of Y-intervals                                  *
C  A    : 4-dimensional array A(0:N,0:M,0:KDIM,0:LDIM); contains *
C                     the spline coefficients: on call,          *
C          A(I,J,0,0) must contain the functional values U(I,J)  *
C          A(I,J,1,0) must contain the derivatives P(I,J) for    *
C                     J=0 and J=M                                *
C          A(I,J,0,1) must contain the derivatives Q(I,J) for    *
C                     I=0 and I=N                                *
C          A(I,J,1,1) must contain the derivatives R(I,J) for    *
C                     I=0 and I=N, as well as for J=0 and J=M    *
C         BICSP1 determines all other A(I,J,K,L) for I=0 to N-1  *
C         and J=0 to M-1. Elements A(N,M,K,L), that are unas-    *
C         signed on call, will remain unassigned.                *
C  X    : (N+1)-vector X(0:N) containing the endpoints of the    *
C         X-intervals                                            *
C  Y    : (N+1)-vector Y(0:M) containing the endpoints of the    *
C         Y-intervals                                            *
C  F    : auxiliary vector F(1:9*MAX(N,M)-8)                     *
C  IERR : initially set to 0. Will be set different from zero,   *
C         if the algorithm detects an error. If errors occur the *
C         program does not complete the computations.            *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: BIC1S1, BIC1S2, BIC1S3, BIC1S4, BIC1S5, *
C                        BIC1S6, BIC1S7, BIC1S8, BIC1S9, TRIDIG  *
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
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM),X(0:N),Y(0:M),F(*)
      DIMENSION AA(0:3,0:3)
C
C*  splitting the vector F for finding X values
C
      IH=1
      ISA=IH+N
      ISB=ISA+N-1
      ISC=ISB+N-1
      ISD=ISC+N-1
      ISX=ISD+N-1
      ISGAMM=ISX+N-1
      ISALPH=ISGAMM+N-1
      ISG=ISALPH+N-1
C
C*  splitting the vector F for finding Y values
C
      JH=1
      JSA=JH+M
      JSB=JSA+M-1
      JSC=JSB+M-1
      JSD=JSC+M-1
      JSX=JSD+M-1
      JSGAMM=JSX+M-1
      JSALPH=JSGAMM+M-1
      JSG=JSALPH+M-1
C
C*  step 1
C
      CALL BIC1S1(N,M,A,X,
     +            F(IH),F(ISA),F(ISB),F(ISC),F(ISD),F(ISX),
     +            F(ISGAMM),F(ISALPH),F(ISG),JERR)
      IERR=JERR
      IF(JERR .NE. 0) RETURN
C
C*  step 2
C
      CALL BIC1S2(N,M,A,Y,
     +            F(JH),F(JSA),F(JSB),F(JSC),F(JSD),F(JSX),
     +            F(JSGAMM),F(JSALPH),F(JSG),JERR)
      IERR=JERR+2
      IF(JERR .NE. 0) RETURN
C
C*  step 3
C
      CALL BIC1S3(N,M,A,X,
     +            F(IH),F(ISA),F(ISB),F(ISC),F(ISD),F(ISX),
     +            F(ISGAMM),F(ISALPH),F(ISG),JERR)
      IERR=JERR+4
      IF(JERR .NE. 0) RETURN
C
C*  step 4
C
      CALL BIC1S4(N,M,A,Y,
     +            F(JH),F(JSA),F(JSB),F(JSC),F(JSD),F(JSX),
     +            F(JSGAMM),F(JSALPH),F(JSG),JERR)
      IERR=JERR+6
      IF(JERR .NE. 0) RETURN
C
C*  steps 5, 6, 7 are contained in step 8 ;
C*  looping over all X-Y values
C
      DO 112 I=0,N-1
        DO 111 J=0,M-1
          CALL BIC1S8(N,M,A,X,Y,I,J,AA)
C
C*  transfer  AA to A
C
          CALL BIC1S9(N,M,A,I,J,AA)
  111   CONTINUE
  112 CONTINUE
C
C*  all spline-coefficients have now been determined
C
      IERR=0
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BIC1S1(N,M,A,X,H,SA,SB,SC,SD,SX,SGAMM,
C[IX{BIC1S1}*)
     +                  SALPH,SG,IERR)
C[BE*)
C
C*****************************************************************
C                                                                *
C  step 1:                                                       *
C                                                                *
C  All vectors, that are needed for solving the linear system of *
C  equations start with the letter S. The notation is taken      *
C  directly from the formulas of the text section, chapter 12.1. *
C  For a description of the parameters see BICSP1.               *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: TRIDIG                                  *
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
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM),X(0:N)
      DIMENSION H(0:N-1)
      DIMENSION SX(N-1),SA(N-1),SB(N-1),SC(N-1),SG(N-1),SGAMM(N-1)
      DIMENSION SALPH(N-1),SD(N-1)
C
C*  determine H(I), test strict monotonicity of X(I)
C
      DO 100 I=0,N-1
        H(I)=X(I+1)-X(I)
        IF(H(I) .LE. 0) GOTO 900
  100 CONTINUE
C
C*  determine the columns of the tridiagonal system
C
      DO 101 I=1,N-1
        SB(I)=1.0D0/H(I-1)
        SC(I)=1.0D0/H(I)
        SD(I)=2.0D0*(SB(I)+SC(I))
  101 CONTINUE
C
C*  for the 1st system, TRIDIG is performed to the end
C
      IREP=0
C
C*  loop covering the M+1 systems of equations
C
      DO 104 J=0,M
C
C*  determine the right-hand side of each system
C
        DO 102 I=1,N-1
          SA(I)=3.0D0*((A(I,J,0,0)-A(I-1,J,0,0))/H(I-1)**2
     +               +(A(I+1,J,0,0)-A(I,J,0,0))/H(I)**2)
  102   CONTINUE
C
C*  corrections for the first and last equations only
C
        SA(1)=SA(1)-A(0,J,1,0)/H(0)
        SA(N-1)=SA(N-1)-A(N,J,1,0)/H(N-1)
C
C*  solve the system, SX is the solution vector
C
        CALL TRIDIG(N-1,SA,SB,SC,SD,IREP,SX,IERR,SGAMM,SALPH,SG)
C
C*  return with IERR=1, if the system is unsolvable
C
        IF(IERR .NE. 0) RETURN
C
C*  to increase speed from the 2nd system onwards (compare TRIDIG)
C
        IREP=1
C
C*  transfer solution vector to array A
C
        DO 103 I=1,N-1
          A(I,J,1,0)=SX(I)
  103   CONTINUE
  104 CONTINUE
C
C*  algorithm has run successfully
C
      RETURN
C
C*  return with IERR=2 if monotonicity of X(I) is not met
C
  900 IERR=2
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BIC1S2(N,M,A,Y,H,SA,SB,SC,SD,SX,SGAMM,
C[IX{BIC1S2}*)
     +                  SALPH,SG,IERR)
C[BE*)
C
C*****************************************************************
C                                                                *
C  Step 2:                                                       *
C                                                                *
C  All vectors, that are needed for solving the system of        *
C  equations start with the letter S. The notation is the same   *
C  as that of the formulas in chapter 12.1.                      *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: TRIDIG                                  *
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
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM),Y(0:M)
      DIMENSION H(0:M-1)
      DIMENSION SX(M-1),SA(M-1),SB(M-1),SC(M-1),SG(M-1),SGAMM(M-1)
      DIMENSION SALPH(M-1),SD(M-1)
C
C*  determine H(J), test strict monotonicity of Y(J)
C
      DO 100 J=0,M-1
        H(J)=Y(J+1)-Y(J)
        IF(H(J) .LE. 0) GOTO 900
  100 CONTINUE
C
C*  determine the columns of the tridiagonal system
C
      DO 101 J=1,M-1
        SB(J)=1.0D0/H(J-1)
        SC(J)=1.0D0/H(J)
        SD(J)=2.0D0*(SB(J)+SC(J))
  101 CONTINUE
C
C*  for the 1st system, TRIDIG is to be performed to the end
C
      IREP=0
C
C*  loop covering the N+1 systems of equations
C
      DO 104 I=0,N
C
C*  determine the right-hand side of each system
C
        DO 102 J=1,M-1
          SA(J)=3.0D0*((A(I,J,0,0)-A(I,J-1,0,0))/H(J-1)**2
     +               +(A(I,J+1,0,0)-A(I,J,0,0))/H(J)**2)
  102   CONTINUE
C
C*  corrections in the first and last equations only
C
        SA(1)=SA(1)-A(I,0,0,1)/H(0)
        SA(M-1)=SA(M-1)-A(I,M,0,1)/H(M-1)
C
C*  solve the system, SX is the solution vector
C
        CALL TRIDIG(M-1,SA,SB,SC,SD,IREP,SX,IERR,SGAMM,SALPH,SG)
C
C*  return with IERR=1, if the system is unsolvable
C
        IF(IERR .NE. 0) RETURN
C
C*  to increase the speed from the 2nd set onward (compare TRIDIG)
C
        IREP=1
C
C*  transfer solution vector to array A
C
        DO 103 J=1,M-1
          A(I,J,0,1)=SX(J)
  103   CONTINUE
  104 CONTINUE
C
C*  algorithm has run successfully
C
      RETURN
C
C*  return with IERR=2, if monotonicity of Y(J) is not true
C
  900 IERR=2
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BIC1S3(N,M,A,X,H,SA,SB,SC,SD,SX,SGAMM,
C[IX{BIC1S3}*)
     +                  SALPH,SG,IERR)
C[BE*)
C
C*****************************************************************
C                                                                *
C  Step 3:                                                       *
C                                                                *
C  All vectors, that are used for solving the system of equations*
C  start with the letter S. The notation follows directly from   *
C  the formulas in chapter 12.1.                                 *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: TRIDIG                                  *
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
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM),X(0:N)
      DIMENSION H(0:N-1)
      DIMENSION SX(N-1),SA(N-1),SB(N-1),SC(N-1),SG(N-1),SGAMM(N-1)
      DIMENSION SALPH(N-1),SD(N-1)
C
C*  determine H(I), test strict monotonicity of X(I)
C
      DO 100 I=0,N-1
        H(I)=X(I+1)-X(I)
        IF(H(I) .LE. 0) GOTO 900
  100 CONTINUE
C
C*  determine the columns of the tridiagonal system
C
      DO 101 I=1,N-1
        SB(I)=1.0D0/H(I-1)
        SC(I)=1.0D0/H(I)
        SD(I)=2.0D0*(SB(I)+SC(I))
  101 CONTINUE
C
C*  in the 1st system, TRIDIG is to be performed to the end
C
      IREP=0
C
C*  loop covering all the systems of equations
C
      DO 104 J=0,M,M
C
C*  determine the right-hand side of each system
C
        DO 102 I=1,N-1
          SA(I)=3.0D0*((A(I,J,0,1)-A(I-1,J,0,1))/H(I-1)**2
     +               +(A(I+1,J,0,1)-A(I,J,0,1))/H(I)**2)
  102   CONTINUE
C
C*  corrections in the first and last equations only
C
        SA(1)=SA(1)-A(0,J,1,1)/H(0)
        SA(N-1)=SA(N-1)-A(N,J,1,1)/H(N-1)
C
C*  solve the system, SX is the solution vector
C
        CALL TRIDIG(N-1,SA,SB,SC,SD,IREP,SX,IERR,SGAMM,SALPH,SG)
C
C*  return with IERR=1, if the system is unsolvable
C
        IF(IERR .NE. 0) RETURN
C
C*  to increase the speed from the 2nd set onwards (compare TRIDIG)
C
        IREP=1
C
C*  transfer solution vector to array A
C
        DO 103 I=1,N-1
          A(I,J,1,1)=SX(I)
  103   CONTINUE
  104 CONTINUE
C
C*  algorithm has run successfully
C
      RETURN
C
C*  return with IERR=2, if monotonicity of X(I) is not true
C
  900 IERR=2
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BIC1S4(N,M,A,Y,H,SA,SB,SC,SD,SX,SGAMM,
C[IX{BIC1S4}*)
     +                  SALPH,SG,IERR)
C[BE*)
C
C*****************************************************************
C                                                                *
C  Step 4:                                                       *
C                                                                *
C  All vectors, that are used for solving the systems of         *
C  equations start with the letter S. The notation is taken      *
C  directly from the formulas of chapter 12.1.                   *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: TRIDIG                                  *
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
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM),Y(0:M)
      DIMENSION H(0:M-1)
      DIMENSION SX(M-1),SA(M-1),SB(M-1),SC(M-1),SG(M-1),SGAMM(M-1)
      DIMENSION SALPH(M-1),SD(M-1)
C
C*  determine H(J), test strict monotonicity of Y(J)
C
      DO 100 J=0,M-1
        H(J)=Y(J+1)-Y(J)
        IF(H(J) .LE. 0) GOTO 900
  100 CONTINUE
C
C*  determine the columns of the tridiagonal system
C
      DO 101 J=1,M-1
        SB(J)=1.0D0/H(J-1)
        SC(J)=1.0D0/H(J)
        SD(J)=2.0D0*(SB(J)+SC(J))
  101 CONTINUE
C
C*  in the 1st set TRIDIG is performed completely
C
      IREP=0
C
C*  loop covering the N+1 system of equations
C
      DO 104 I=0,N
C
C*  determine the right-hand side of each system
C
        DO 102 J=1,M-1
          SA(J)=3.0D0*((A(I,J,1,0)-A(I,J-1,1,0))/H(J-1)**2
     +               +(A(I,J+1,1,0)-A(I,J,1,0))/H(J)**2)
  102   CONTINUE
C
C*  correction in the first and last equations
C
        SA(1)=SA(1)-A(I,0,1,1)/H(0)
        SA(M-1)=SA(M-1)-A(I,M,1,1)/H(M-1)
C
C*  solve the system, SX is the solution vector
C
        CALL TRIDIG(M-1,SA,SB,SC,SD,IREP,SX,IERR,SGAMM,SALPH,SG)
C
C*  return with IERR=1, if the system is unsolvable
C
        IF(IERR .NE. 0) RETURN
C
C*  to increase speed from the 2nd system onwards (compare TRIDIG)
C
        IREP=1
C
C*  transfer solution vector to array A
C
        DO 103 J=1,M-1
          A(I,J,1,1)=SX(J)
  103   CONTINUE
  104 CONTINUE
C
C*  algorithm has run successfully
C
      RETURN
C
C*  return with IERR=2, if monotonicity of Y(J) is not true
C
  900 IERR=2
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BIC1S5(I,G,N,X)
C[IX{BIC1S5}*)
C[BE*)
C
C*****************************************************************
C                                                                *
C  Step 5:                                                       *
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
      DIMENSION G(0:3,0:3),X(0:N),DAT(0:3,0:3),IPT(0:3,0:3),
     F          HPT(0:3)
C
C*  unchanging coefficients of the desired matrix
C
      DATA DAT/1.0D0,0.0D0,-3.0D0,2.0D0, 0.0D0,1.0D0,-2.0D0,1.0D0,
     F         0.0D0,0.0D0,3.0D0,-2.0D0, 0.0D0,0.0D0,-1.0D0,1.0D0/
C
C*  negative H-power of the desired matrix
C
      DATA IPT/0,0,2,3, 0,0,1,2, 0,0,2,3, 0,0,1,2/
C
C*  negative H-powers
C
      DATA HPT(0)/1.0D0/
      H=X(I+1)-X(I)
      DO 100 K=1,3
        HPT(K)=HPT(K-1)/H
  100 CONTINUE
C
C*  determine the matrix
C
      DO 102 K=0,3
        DO 101 L=0,3
          IPOT=IPT(K,L)
          G(K,L)=DAT(K,L)*HPT(IPOT)
  101   CONTINUE
  102 CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BIC1S6(J,G,M,Y)
C[IX{BIC1S6}*)
C[BE*)
C
C*****************************************************************
C                                                                *
C  Step 6:                                                       *
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
      DIMENSION G(0:3,0:3),Y(0:M),DAT(0:3,0:3),IPT(0:3,0:3),
     F          HPT(0:3)
C
C*  unchanging coefficients of the desired matrix
C
      DATA DAT/1.0D0,0.0D0,-3.0D0,2.0D0, 0.0D0,1.0D0,-2.0D0,1.0D0,
     F         0.0D0,0.0D0,3.0D0,-2.0D0, 0.0D0,0.0D0,-1.0D0,1.0D0/
C
C*  negative H-powers of the desired matrix
C
      DATA IPT/0,0,2,3, 0,0,1,2, 0,0,2,3, 0,0,1,2/
C
C*  negative H-powers
C
      DATA HPT(0)/1.0D0/
      H=Y(J+1)-Y(J)
      DO 100 K=1,3
        HPT(K)=HPT(K-1)/H
  100 CONTINUE
C
C*  determine the transposed matrix
C
      DO 102 K=0,3
        DO 101 L=0,3
          IPOT=IPT(K,L)
          G(L,K)=DAT(K,L)*HPT(IPOT)
  101   CONTINUE
  102 CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BIC1S7(N,M,A,I,J,W)
C[IX{BIC1S7}*)
C[BE*)
C
C*****************************************************************
C                                                                *
C  Step 7:                                                       *
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
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM),W(0:3,0:3)
      DO 102 L=0,1
        DO 101 K=0,1
          W(K,L)=A(I,J,K,L)
          W(K+2,L)=A(I+1,J,K,L)
          W(K,L+2)=A(I,J+1,K,L)
          W(K+2,L+2)=A(I+1,J+1,K,L)
  101   CONTINUE
  102 CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BIC1S8(N,M,A,X,Y,I,J,AA)
C[IX{BIC1S8}*)
C[BE*)
C
C*****************************************************************
C                                                                *
C  Step 8:                                                       *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: BIC1S5, BIC1S6, BIC1S7                  *
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
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM),X(0:N),Y(0:M)
      DIMENSION AA(0:3,0:3),GX(0:3,0:3),GYT(0:3,0:3),W(0:3,0:3)
      DIMENSION WGYT(0:3,0:3)
C
C*  determine G(I,J)(X) = GX
C
        CALL BIC1S5(I,GX,N,X)
C
C*  determine G(I,J)(Y) transposed = GY
C
      CALL BIC1S6(J,GYT,M,Y)
C
C*  determine W(I,J)
C
      CALL BIC1S7(N,M,A,I,J,W)
C
C*  W(I,J)*G(J)(Y) transposed = WGYT
C
      DO 103 K=0,3
        DO 102 L=0,3
          WGYT(K,L)=0.0D0
          DO 101 KL=0,3
            WGYT(K,L)=WGYT(K,L)+W(K,KL)*GYT(KL,L)
  101     CONTINUE
  102   CONTINUE
  103 CONTINUE
C
C*  determine A(I,J) = AA
C
      DO 106 K=0,3
        DO 105 L=0,3
        AA(K,L)=0.0D0
          DO 104 KL=0,3
            AA(K,L)=AA(K,L)+GX(K,KL)*WGYT(KL,L)
  104     CONTINUE
  105   CONTINUE
  106 CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BIC1S9(N,M,A,I,J,AA)
C[IX{BIC1S9}*)
C[BE*)
C
C*****************************************************************
C                                                                *
C  Step 9:                                                       *
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
      DIMENSION A(0:N,0:M,0:KDIM,0:LDIM)
      DIMENSION AA(0:3,0:3)
C
      DO 102 K=0,3
        DO 101 L=0,3
          A(I,J,K,L)=AA(K,L)
  101     CONTINUE
  102 CONTINUE
      RETURN
      END
