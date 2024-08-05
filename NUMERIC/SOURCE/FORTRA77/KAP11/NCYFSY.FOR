C[BA*)
C[LE*)
      SUBROUTINE NCYFSY (N,DM,DU1,DU2,RS,X,DML,DL1L,DL2L,RL2L,
C[IX{NCYFSY}*)
     +                   RL1L,DU1U,DU2U,RC2U,RC1U,IERR)
C
C*****************************************************************
C                                                                *
C  NCYFSY computes the solution X of a linear system of equations*
C  A*X=RS  for a symmetric, almost cyclic five-diagonal matrix A.*
C[BE*)
C  The matrix A is defined by the three vectors DM, DU1 and DU2  *
C  which give its diagonal and first and second co-diagonal      *
C  entries. The linear system has the form:                      *
C                                                                *
C  DM(1)*X(1) + DU1(1)*X(2) + DU2(1)*X(3) +                      *
C                       + DU2(N-1)*X(N-1) + DU1(N)*X(N) = RS(1)  *
C                                                                *
C  DU1(1)*X(1) + DM(2)*X(2) + DU1(2)*X(3)                        *
C                           + DU2(2)*X(4) + DU2(N)*X(N) = RS(2)  *
C                                                                *
C  DU2(I-2)*X(I-2) + DU1(I-1)*X(I-1) + DM(I)*X(I) +              *
C                       + DU1(I)*X(I+1) + DU2(I)*X(I+2) = RS(I), *
C                                 fÅr I=3(1)N-2                  *
C                                                                *
C  DU2(N-1)*X(1) + DU2(N-3)*X(N-3) +                             *
C   + DU1(N-2)*X(N-2) + DM(N-1)*X(N-1) + DU1(N-1)*X(N) = RS(N-1) *
C                                                                *
C  DU1(N)*X(1) + DU2(N)*X(2) +                                   *
C      + DU2(N-2)*X(N-2) + DU1(N-1)*X(N-1) + DM(N)*X(N) = RS(N)  *
C                                                                *
C                                                                *
C  ASSUMPTION:    N > 5                                          *
C  ===========                                                   *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N   :  Number of equations, size of the system matrix A       *
C                                                                *
C  DM  :]   vectors DM(1:N), DU1(1:N), DU2(1:N),                 *
C  DU1 :]   representing the diagonal, first and second          *
C  DU2 :]   co-diagonal of A, respectively                       *
C                                                                *
C  RS  : vector RS(1:N); the right hand side                     *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  DML  :  vector DML(1:N)    ]   describing the entries of the  *
C  DL1L :  vector DL1L(1:N)   ]   lower tringular matrix C in the*
C  DL2L :  vector DL2L(1:N)   ]   factorization A=C*B            *
C  RL2L :  vector RL2L(1:N-4) ]   (refer to NCYFSP)              *
C  RL1L :  vector RL1L(1:N-3) ]                                  *
C                                                                *
C  DU1U :  vector DU1U(1:N)    ]   describing the entries of the *
C  DU2U :  vector DU2U(1:N)    ]   upper tringular matrix B in   *
C  RC2U :  vector RC2U(1:N-4)  ]   the factorization A=C*B       *
C  RC1U :  vector RC1U(1:N-3)  ]   (refer to NCYFSP)             *
C                                                                *
C  X    :  vector X(1:N) ;  the solution                         *
C                                                                *
C  IERR :  error parameter                                       *
C          =  0 :  All is o.k.                                   *
C          = -1 :  N < 6                                         *
C          =  1 :  default in NCYFSP (system matrix singular)    *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: NCYFSP, NCYFSS                          *
C                                                                *
C*****************************************************************
C                                                                *
C  Author  : GÅnter Palm                                         *
C  Date    : 04.18.1988                                          *
C  Source  : FORTRAN 77                                          *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DM(N), DU1(N), DU2(N), RS(N), X(N), DML(N),
     +                 DL1L(N), DL2L(N), RL2L(N-4), RL1L(N-3),
     +                 DU1U(N), DU2U(N), RC2U(N-4), RC1U(N-3)
C
C-----Checking the input
C
      IERR = -1
      IF (N .LT. 6) RETURN
C
C-----Factoring the system matrix into a lower
C     times an upper triangular matrix
C
      CALL NCYFSP (N,DM,DU1,DU2,DML,DL1L,DL2L,RL2L,RL1L,DU1U,
     +             DU2U,RC2U,RC1U,IERR)
C
C-----If IERR = 0, we can update and backsubstitute
C
      IF (IERR .EQ. 0) THEN
         CALL NCYFSS (N,RS,X,DML,DL1L,DL2L,RL2L,RL1L,DU1U,DU2U,
     +                RC2U,RC1U)
      ENDIF
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE NCYFSP (N,DM,DU1,DU2,DML,DL1L,DL2L,RL2L,RL1L,
C[IX{NCYFSP}*)
     +                   DU1U,DU2U,RC2U,RC1U,IERR)
C
C*****************************************************************
C                                                                *
C  Factorization of a symmetric, almost cyclic five-diagonal     *
C  matrix A into a product of a lower triangular matrix C and an *
C  upper triangular matrix B.                                    *
C[BE*)
C  The matrix A is represented by                                *
C  three vectors DM, DU1 and DU2, as described in SUBROUTINE     *
C  NCYFSY.                                                       *
C                                                                *
C                                                                *
C  ASSUMPTIONS:    N > 5  (this is not checked again)            *
C  ============                                                  *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N   : size of the system matrix                               *
C                                                                *
C  DM  :]   vectors representing the diagonal and co-diagonals   *
C  DU1 :]   of A, dimensioned as ..(1:N)                         *
C  DU2 :]                                                        *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  DML  :  vector DML(1:N)    ]   describing the entries of the  *
C  DL1L :  vector DL1L(1:N)   ]   lower tringular matrix C in the*
C  DL2L :  vector DL2L(1:N)   ]   factorization A=C*B            *
C  RL2L :  vector RL2L(1:N-4) ]                                  *
C  RL1L :  vector RL1L(1:N-3) ]                                  *
C                                                                *
C  DU1U :  vector DU1U(1:N)    ]   describing the entries of the *
C  DU2U :  vector DU2U(1:N)    ]   upper tringular matrix B in   *
C  RC2U :  vector RC2U(1:N-4)  ]   the factorization A=C*B       *
C  RC1U :  vector RC1U(1:N-3)  ]                                 *
C                                                                *
C      In particular we use the following notation:              *
C      DML  :  main diagonal of C,      DML(I),  I=1, ..., N     *
C      DL1L :  first co-diagonal of C,  DL1L(I), I=2, ..., N     *
C      DL2L :  second co-diagonal of C, DL2L(I), I=3, ..., N     *
C      RL2L :  second to last row of C, except for the elements  *
C              labelled N-3, ..., N,    RL2L(I), I=1, ..., N-4   *
C      RL1L :  last row of C, except for elements labelled       *
C              N-2,...,N,               RL1L(I), I=1, ..., N-3   *
C      DU1U :  first co-diagonal of B,  DU1U(I), I=1, ..., N-1   *
C      DU2U :  second co-diagonal of B, DU2U(I), I=1, ..., N-2   *
C      RC2U :  last but one column of B, except for the elements *
C              labelled N-3,...,N,      RC2U(I), I=1, ..., N-4   *
C      RC1U :  last column of B, except for the elements labelled*
C              N-2,...,N,               RC1U(I), I=1, ..., N-3   *
C                                                                *
C  IERR :  error parameter                                       *
C          =  0 :  All is o.k.                                   *
C          =  1 :  default, because of intended division by an   *
C                  element DML(I), whose magnitude does not      *
C                  excced 4 * machine constant.                  *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: MACHPD                                  *
C                                                                *
C*****************************************************************
C                                                                *
C  Author  : GÅnter Palm                                         *
C  Date    : 04.18.1988                                          *
C  Source  : FORTRAN 77                                          *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DM(N), DU1(N), DU2(N), DML(N), DL1L(N),
     +                 DL2L(N), RL2L(N-4), RL1L(N-3), DU1U(N),
     +                 DU2U(N), RC2U(N-4), RC1U(N-3)
      LOGICAL FLAG
      SAVE FLAG, FMACHP
C
C-----Initializing
C
      DATA FLAG /.TRUE./
      IERR = 1
C
C-----Compute the machine constant (upon first call only)
C
      IF (FLAG) THEN
        FMACHP = 1.0D0
    5   FMACHP = 0.5D0*FMACHP
        IF (MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 5
        FMACHP = 8.0D0*FMACHP
        FLAG = .FALSE.
      ENDIF
C
C-----Factor the matrix A into a lower and an upper triangular
C     matrix. The computations are stopped if one of the diagonal
C     entries in the lower triangular factor C does not exceed
C     4 * machine constant in magnitude.
C
      DML(1)  = DM(1)
      IF (DML(1) .LE. FMACHP) RETURN
      DU1U(1) = DU1(1)/DML(1)
      RC2U(1) = DU2(N-1)/DML(1)
      RC1U(1) = DU1(N)/DML(1)
      DL1L(2) = DU1(1)
      DML(2)  = DM(2) - DL1L(2)*DU1U(1)
      IF (DML(2) .LE. FMACHP) RETURN
      RC2U(2) = -(RC2U(1)*DL1L(2))/DML(2)
      RC1U(2) = (DU2(N) - DL1L(2)*RC1U(1))/DML(2)
C
      DO 10 I=3,N-2,1
        K = I-1
        J = I-2
        DL2L(I) = DU2(I-2)
        DU2U(J) = DU2(J)/DML(J)
        DU1U(K) = (DU1(K) - DL1L(K)*DU2U(J))/DML(K)
        DL1L(I) = DU1(I-1) - DL2L(I)*DU1U(J)
        DML(I)  = DM(I) - DL1L(I)*DU1U(K) - DL2L(I)*DU2U(J)
        IF (DML(I) .LE. FMACHP) RETURN
   10 CONTINUE
C
      DO 20 I=3,N-4,1
        RC2U(I) = -(DL2L(I)*RC2U(I-2) + DL1L(I)*RC2U(I-1))/DML(I)
   20 CONTINUE
C
      DO 30 I=3,N-3,1
        RC1U(I) = -(DL2L(I)*RC1U(I-2) + DL1L(I)*RC1U(I-1))/DML(I)
   30 CONTINUE
C
      DU2U(N-3) = (DU2(N-3) - DL1L(N-3)*RC2U(N-4) -
     +             DL2L(N-3)*RC2U(N-5))/DML(N-3)
      DU2U(N-2) = (DU2(N-2) - DL1L(N-2)*RC1U(N-3) -
     +             DL2L(N-2)*RC1U(N-4))/DML(N-2)
      DU1U(N-2) = (DU1(N-2) - DL1L(N-2)*DU2U(N-3) -
     +             DL2L(N-2)*RC2U(N-4))/DML(N-2)
C
      RL2L(1) = DU2(N-1)
      RL2L(2) = -RL2L(1)*DU1U(1)
      DO 40 I=3,N-4,1
        RL2L(I) = -(RL2L(I-2)*DU2U(I-2) + RL2L(I-1)*DU1U(I-1))
   40 CONTINUE
      RL1L(1) = DU1(N)
      RL1L(2) = DU2(N) - RL1L(1)*DU1U(1)
C
      DO 50 I=3,N-3,1
        RL1L(I) = -(RL1L(I-2)*DU2U(I-2) + RL1L(I-1)*DU1U(I-1))
   50 CONTINUE
C
      DL2L(N-1) = DU2(N-3) - (RL2L(N-5)*DU2U(N-5) +
     +            RL2L(N-4)*DU1U(N-4))
      DL2L(N)   = DU2(N-2) - (RL1L(N-4)*DU2U(N-4) +
     +            RL1L(N-3)*DU1U(N-3))
      DL1L(N-1) = DU1(N-2) - (RL2L(N-4)*DU2U(N-4) +
     +            DL2L(N-1)*DU1U(N-3))
      DUMMY1 = 0.0D0
      DUMMY2 = 0.0D0
      DUMMY3 = 0.0D0
      DO 60 K=1,N-4,1
        DUMMY1 = DUMMY1 + RL1L(K)*RC2U(K)
        DUMMY2 = DUMMY2 + RL2L(K)*RC2U(K)
        DUMMY3 = DUMMY3 + RL2L(K)*RC1U(K)
   60 CONTINUE
C
      DL1L(N)   = DU1(N-1) - DUMMY1 - RL1L(N-3)*DU2U(N-3) -
     +            DL2L(N)*DU1U(N-2)
C
      DML(N-1)  = DM(N-1) - DUMMY2 - DL2L(N-1)*DU2U(N-3) -
     +            DL1L(N-1)*DU1U(N-2)
      IF (DML(N-1) .LE. FMACHP) RETURN
C
      DU1U(N-1) = (DU1(N-1) - DUMMY3 - DL2L(N-1)*RC1U(N-3) -
     +             DL1L(N-1)*DU2U(N-2))/DML(N-1)
C
      DUMMY1 = 0.0D0
      DO 70 K=1,N-3,1
        DUMMY1 = DUMMY1 + RL1L(K)*RC1U(K)
   70 CONTINUE
      DML(N) = DM(N) - DUMMY1 - DL2L(N)*DU2U(N-2) -
     +         DL1L(N)*DU1U(N-1)
      IF (DML(N) .LE. FMACHP) RETURN
C
      IERR = 0
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE NCYFSS (N,RS,X,DML,DL1L,DL2L,RL2L,RL1L,DU1U,
C[IX{NCYFSS}*)
     +                   DU2U,RC2U,RC1U)
C
C*****************************************************************
C                                                                *
C  Find the solution of a linear system of equations A*X=RS for  *
C  a symmetric, cyclic five-diagonal matrix A, once the factor   *
C  matrizes C and B have been determind in SUBROUTINE NCYFSP.    *
C[BE*)
C  C and B are input parameters and are represented by the       *
C  vectors DML, ... , RC1U as described in SUBROUTINE NCYFSY.    *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N    :  size of the matrix A                                  *
C  RS   :  vector RS(1:N); the right hand side                   *
C                                                                *
C  DML  :  vector DML(1:N)    ]   describing the entries of the  *
C  DL1L :  vector DL1L(1:N)   ]   lower tringular matrix C in the*
C  DL2L :  vector DL2L(1:N)   ]   factorization A=C*B            *
C  RL2L :  vector RL2L(1:N-4) ]   (refer to NCYFSP)              *
C  RL1L :  vector RL1L(1:N-3) ]                                  *
C                                                                *
C  DU1U :  vector DU1U(1:N)    ]   describing the entries of the *
C  DU2U :  vector DU2U(1:N)    ]   upper tringular matrix B in   *
C  RC2U :  vector RC2U(1:N-4)  ]   the factorization A=C*B       *
C  RC1U :  vector RC1U(1:N-3)  ]   (refer to NCYFSP)             *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  X :  vector  X(1:N); the solution vector                      *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author  : GÅnter Palm                                         *
C  Date    : 04.18.1988                                          *
C  Source  : FORTRAN 77                                          *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION RS(N), X(N), DML(N), DL1L(N), DL2L(N),
     +                 RL2L(N-4), RL1L(N-3), DU1U(N), DU2U(N),
     +                 RC2U(N-4), RC1U(N-3)
C
C-----Solving the system by updating and
C     backsubstitution
C
C     1. Updating the right hand side
C
      X(1) = RS(1)/DML(1)
      X(2) = (RS(2) -X(1)*DL1L(2))/DML(2)
      DO 80 I=3,N-2,1
        X(I) = (RS(I) - X(I-2)*DL2L(I) - X(I-1)*DL1L(I))/DML(I)
   80 CONTINUE
      DUMMY1 = 0.0D0
      DO 90 K=1,N-4,1
        DUMMY1 = DUMMY1 + X(K)*RL2L(K)
   90 CONTINUE
      X(N-1) = (RS(N-1)-DUMMY1-X(N-3)*DL2L(N-1)-X(N-2)*DL1L(N-1))
     +          / DML(N-1)
      DUMMY1 = 0.0D0
      DO 100 K=1,N-3,1
        DUMMY1 = DUMMY1 + X(K)*RL1L(K)
  100 CONTINUE
      X(N) = (RS(N) - DUMMY1 - X(N-2)*DL2L(N) - X(N-1)*DL1L(N))
     +        / DML(N)
C
C     2. Backsubstitution
C
      X(N-1) = X(N-1) - DU1U(N-1)*X(N)
      X(N-2) = X(N-2) - DU1U(N-2)*X(N-1) - DU2U(N-2)*X(N)
      X(N-3) = X(N-3) - DU1U(N-3)*X(N-2) - DU2U(N-3)*X(N-1) -
     +         RC1U(N-3)*X(N)
      DO 110 I=N-4,1,-1
        X(I) = X(I) - DU1U(I)*X(I+1) - DU2U(I)*X(I+2) -
     +         RC2U(I)*X(N-1) - RC1U(I)*X(N)
  110 CONTINUE
      RETURN
      END
