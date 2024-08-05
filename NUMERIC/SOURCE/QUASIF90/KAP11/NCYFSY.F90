      SUBROUTINE NCYFSY (N, DM, DU1, DU2, RS, X, DML, DL1L, DL2L, RL2L, &
      RL1L, DU1U, DU2U, RC2U, RC1U, IERR)                               
!                                                                       
!*****************************************************************      
!                                                                *      
!  NCYFSY computes the solution X of a linear system of equations*      
!  A*X=RS  for a symmetric, almost cyclic five-diagonal matrix A.*      
!  The matrix A is defined by the three vectors DM, DU1 and DU2  *      
!  which give its diagonal and first and second co-diagonal      *      
!  entries. The linear system has the form:                      *      
!                                                                *      
!  DM(1)*X(1) + DU1(1)*X(2) + DU2(1)*X(3) +                      *      
!                       + DU2(N-1)*X(N-1) + DU1(N)*X(N) = RS(1)  *      
!                                                                *      
!  DU1(1)*X(1) + DM(2)*X(2) + DU1(2)*X(3)                        *      
!                           + DU2(2)*X(4) + DU2(N)*X(N) = RS(2)  *      
!                                                                *      
!  DU2(I-2)*X(I-2) + DU1(I-1)*X(I-1) + DM(I)*X(I) +              *      
!                       + DU1(I)*X(I+1) + DU2(I)*X(I+2) = RS(I), *      
!                                 fÅr I=3(1)N-2                  *      
!                                                                *      
!  DU2(N-1)*X(1) + DU2(N-3)*X(N-3) +                             *      
!   + DU1(N-2)*X(N-2) + DM(N-1)*X(N-1) + DU1(N-1)*X(N) = RS(N-1) *      
!                                                                *      
!  DU1(N)*X(1) + DU2(N)*X(2) +                                   *      
!      + DU2(N-2)*X(N-2) + DU1(N-1)*X(N-1) + DM(N)*X(N) = RS(N)  *      
!                                                                *      
!                                                                *      
!  ASSUMPTION:    N > 5                                          *      
!  ===========                                                   *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N   :  Number of equations, size of the system matrix A       *      
!                                                                *      
!  DM  :]   vectors DM(1:N), DU1(1:N), DU2(1:N),                 *      
!  DU1 :]   representing the diagonal, first and second          *      
!  DU2 :]   co-diagonal of A, respectively                       *      
!                                                                *      
!  RS  : vector RS(1:N); the right hand side                     *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  DML  :  vector DML(1:N)    ]   describing the entries of the  *      
!  DL1L :  vector DL1L(1:N)   ]   lower tringular matrix C in the*      
!  DL2L :  vector DL2L(1:N)   ]   factorization A=C*B            *      
!  RL2L :  vector RL2L(1:N-4) ]   (refer to NCYFSP)              *      
!  RL1L :  vector RL1L(1:N-3) ]                                  *      
!                                                                *      
!  DU1U :  vector DU1U(1:N)    ]   describing the entries of the *      
!  DU2U :  vector DU2U(1:N)    ]   upper tringular matrix B in   *      
!  RC2U :  vector RC2U(1:N-4)  ]   the factorization A=C*B       *      
!  RC1U :  vector RC1U(1:N-3)  ]   (refer to NCYFSP)             *      
!                                                                *      
!  X    :  vector X(1:N) ;  the solution                         *      
!                                                                *      
!  IERR :  error parameter                                       *      
!          =  0 :  All is o.k.                                   *      
!          = -1 :  N < 6                                         *      
!          =  1 :  default in NCYFSP (system matrix singular)    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: NCYFSP, NCYFSS                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author  : GÅnter Palm                                         *      
!  Date    : 04.18.1988                                          *      
!  Source  : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DM (N), DU1 (N), DU2 (N), RS (N), X (N), DML (N), &
      DL1L (N), DL2L (N), RL2L (N - 4), RL1L (N - 3), DU1U (N), DU2U (N)&
      , RC2U (N - 4), RC1U (N - 3)                                      
!                                                                       
!-----Checking the input                                                
!                                                                       
      IERR = - 1 
      IF (N.LT.6) RETURN 
!                                                                       
!-----Factoring the system matrix into a lower                          
!     times an upper triangular matrix                                  
!                                                                       
      CALL NCYFSP (N, DM, DU1, DU2, DML, DL1L, DL2L, RL2L, RL1L, DU1U,  &
      DU2U, RC2U, RC1U, IERR)                                           
!                                                                       
!-----If IERR = 0, we can update and backsubstitute                     
!                                                                       
      IF (IERR.EQ.0) THEN 
         CALL NCYFSS (N, RS, X, DML, DL1L, DL2L, RL2L, RL1L, DU1U, DU2U,&
         RC2U, RC1U)                                                    
      ENDIF 
      RETURN 
      END SUBROUTINE NCYFSY                         
!                                                                       
!                                                                       
      SUBROUTINE NCYFSP (N, DM, DU1, DU2, DML, DL1L, DL2L, RL2L, RL1L,  &
      DU1U, DU2U, RC2U, RC1U, IERR)                                     
!                                                                       
!*****************************************************************      
!                                                                *      
!  Factorization of a symmetric, almost cyclic five-diagonal     *      
!  matrix A into a product of a lower triangular matrix C and an *      
!  upper triangular matrix B. The matrix A is represented by     *      
!  three vectors DM, DU1 and DU2, as described in SUBROUTINE     *      
!  NCYFSY.                                                       *      
!                                                                *      
!                                                                *      
!  ASSUMPTIONS:    N > 5  (this is not checked again)            *      
!  ============                                                  *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N   : size of the system matrix                               *      
!                                                                *      
!  DM  :]   vectors representing the diagonal and co-diagonals   *      
!  DU1 :]   of A, dimensioned as ..(1:N)                         *      
!  DU2 :]                                                        *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  DML  :  vector DML(1:N)    ]   describing the entries of the  *      
!  DL1L :  vector DL1L(1:N)   ]   lower tringular matrix C in the*      
!  DL2L :  vector DL2L(1:N)   ]   factorization A=C*B            *      
!  RL2L :  vector RL2L(1:N-4) ]                                  *      
!  RL1L :  vector RL1L(1:N-3) ]                                  *      
!                                                                *      
!  DU1U :  vector DU1U(1:N)    ]   describing the entries of the *      
!  DU2U :  vector DU2U(1:N)    ]   upper tringular matrix B in   *      
!  RC2U :  vector RC2U(1:N-4)  ]   the factorization A=C*B       *      
!  RC1U :  vector RC1U(1:N-3)  ]                                 *      
!                                                                *      
!      In particular we use the following notation:              *      
!      DML  :  main diagonal of C,      DML(I),  I=1, ..., N     *      
!      DL1L :  first co-diagonal of C,  DL1L(I), I=2, ..., N     *      
!      DL2L :  second co-diagonal of C, DL2L(I), I=3, ..., N     *      
!      RL2L :  second to last row of C, except for the elements  *      
!              labelled N-3, ..., N,    RL2L(I), I=1, ..., N-4   *      
!      RL1L :  last row of C, except for elements labelled       *      
!              N-2,...,N,               RL1L(I), I=1, ..., N-3   *      
!      DU1U :  first co-diagonal of B,  DU1U(I), I=1, ..., N-1   *      
!      DU2U :  second co-diagonal of B, DU2U(I), I=1, ..., N-2   *      
!      RC2U :  last but one column of B, except for the elements *      
!              labelled N-3,...,N,      RC2U(I), I=1, ..., N-4   *      
!      RC1U :  last column of B, except for the elements labelled*      
!              N-2,...,N,               RC1U(I), I=1, ..., N-3   *      
!                                                                *      
!  IERR :  error parameter                                       *      
!          =  0 :  All is o.k.                                   *      
!          =  1 :  default, because of intended division by an   *      
!                  element DML(I), whose magnitude does not      *      
!                  excced 4 * machine constant.                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author  : GÅnter Palm                                         *      
!  Date    : 04.18.1988                                          *      
!  Source  : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DM (N), DU1 (N), DU2 (N), DML (N), DL1L (N),      &
      DL2L (N), RL2L (N - 4), RL1L (N - 3), DU1U (N), DU2U (N), RC2U (N &
      - 4), RC1U (N - 3)                                                
      LOGICAL FLAG 
      SAVE FLAG, FMACHP 
!                                                                       
!-----Initializing                                                      
!                                                                       
      DATA FLAG / .TRUE. / 
      IERR = 1 
!                                                                       
!-----Compute the machine constant (upon first call only)               
!                                                                       
      IF (FLAG) THEN 
         FMACHP = 1.0D0 
    5    FMACHP = 0.5D0 * FMACHP 
         IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 5 
         FMACHP = 8.0D0 * FMACHP 
         FLAG = .FALSE. 
      ENDIF 
!                                                                       
!-----Factor the matrix A into a lower and an upper triangular          
!     matrix. The computations are stopped if one of the diagonal       
!     entries in the lower triangular factor C does not exceed          
!     4 * machine constant in magnitude.                                
!                                                                       
      DML (1) = DM (1) 
      IF (DML (1) .LE.FMACHP) RETURN 
      DU1U (1) = DU1 (1) / DML (1) 
      RC2U (1) = DU2 (N - 1) / DML (1) 
      RC1U (1) = DU1 (N) / DML (1) 
      DL1L (2) = DU1 (1) 
      DML (2) = DM (2) - DL1L (2) * DU1U (1) 
      IF (DML (2) .LE.FMACHP) RETURN 
      RC2U (2) = - (RC2U (1) * DL1L (2) ) / DML (2) 
      RC1U (2) = (DU2 (N) - DL1L (2) * RC1U (1) ) / DML (2) 
!                                                                       
      DO 10 I = 3, N - 2, 1 
         K = I - 1 
         J = I - 2 
         DL2L (I) = DU2 (I - 2) 
         DU2U (J) = DU2 (J) / DML (J) 
         DU1U (K) = (DU1 (K) - DL1L (K) * DU2U (J) ) / DML (K) 
         DL1L (I) = DU1 (I - 1) - DL2L (I) * DU1U (J) 
         DML (I) = DM (I) - DL1L (I) * DU1U (K) - DL2L (I) * DU2U (J) 
         IF (DML (I) .LE.FMACHP) RETURN 
   10 END DO 
!                                                                       
      DO 20 I = 3, N - 4, 1 
         RC2U (I) = - (DL2L (I) * RC2U (I - 2) + DL1L (I) * RC2U (I - 1)&
         ) / DML (I)                                                    
   20 END DO 
!                                                                       
      DO 30 I = 3, N - 3, 1 
         RC1U (I) = - (DL2L (I) * RC1U (I - 2) + DL1L (I) * RC1U (I - 1)&
         ) / DML (I)                                                    
   30 END DO 
!                                                                       
      DU2U (N - 3) = (DU2 (N - 3) - DL1L (N - 3) * RC2U (N - 4) - DL2L (&
      N - 3) * RC2U (N - 5) ) / DML (N - 3)                             
      DU2U (N - 2) = (DU2 (N - 2) - DL1L (N - 2) * RC1U (N - 3) - DL2L (&
      N - 2) * RC1U (N - 4) ) / DML (N - 2)                             
      DU1U (N - 2) = (DU1 (N - 2) - DL1L (N - 2) * DU2U (N - 3) - DL2L (&
      N - 2) * RC2U (N - 4) ) / DML (N - 2)                             
!                                                                       
      RL2L (1) = DU2 (N - 1) 
      RL2L (2) = - RL2L (1) * DU1U (1) 
      DO 40 I = 3, N - 4, 1 
         RL2L (I) = - (RL2L (I - 2) * DU2U (I - 2) + RL2L (I - 1)       &
         * DU1U (I - 1) )                                               
   40 END DO 
      RL1L (1) = DU1 (N) 
      RL1L (2) = DU2 (N) - RL1L (1) * DU1U (1) 
!                                                                       
      DO 50 I = 3, N - 3, 1 
         RL1L (I) = - (RL1L (I - 2) * DU2U (I - 2) + RL1L (I - 1)       &
         * DU1U (I - 1) )                                               
   50 END DO 
!                                                                       
      DL2L (N - 1) = DU2 (N - 3) - (RL2L (N - 5) * DU2U (N - 5) + RL2L (&
      N - 4) * DU1U (N - 4) )                                           
      DL2L (N) = DU2 (N - 2) - (RL1L (N - 4) * DU2U (N - 4) + RL1L (N - &
      3) * DU1U (N - 3) )                                               
      DL1L (N - 1) = DU1 (N - 2) - (RL2L (N - 4) * DU2U (N - 4) + DL2L (&
      N - 1) * DU1U (N - 3) )                                           
      DUMMY1 = 0.0D0 
      DUMMY2 = 0.0D0 
      DUMMY3 = 0.0D0 
      DO 60 K = 1, N - 4, 1 
         DUMMY1 = DUMMY1 + RL1L (K) * RC2U (K) 
         DUMMY2 = DUMMY2 + RL2L (K) * RC2U (K) 
         DUMMY3 = DUMMY3 + RL2L (K) * RC1U (K) 
   60 END DO 
!                                                                       
      DL1L (N) = DU1 (N - 1) - DUMMY1 - RL1L (N - 3) * DU2U (N - 3)     &
      - DL2L (N) * DU1U (N - 2)                                         
!                                                                       
      DML (N - 1) = DM (N - 1) - DUMMY2 - DL2L (N - 1) * DU2U (N - 3)   &
      - DL1L (N - 1) * DU1U (N - 2)                                     
      IF (DML (N - 1) .LE.FMACHP) RETURN 
!                                                                       
      DU1U (N - 1) = (DU1 (N - 1) - DUMMY3 - DL2L (N - 1) * RC1U (N - 3)&
      - DL1L (N - 1) * DU2U (N - 2) ) / DML (N - 1)                     
!                                                                       
      DUMMY1 = 0.0D0 
      DO 70 K = 1, N - 3, 1 
         DUMMY1 = DUMMY1 + RL1L (K) * RC1U (K) 
   70 END DO 
      DML (N) = DM (N) - DUMMY1 - DL2L (N) * DU2U (N - 2) - DL1L (N)    &
      * DU1U (N - 1)                                                    
      IF (DML (N) .LE.FMACHP) RETURN 
!                                                                       
      IERR = 0 
      RETURN 
      END SUBROUTINE NCYFSP                         
!                                                                       
!                                                                       
      SUBROUTINE NCYFSS (N, RS, X, DML, DL1L, DL2L, RL2L, RL1L, DU1U,   &
      DU2U, RC2U, RC1U)                                                 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Find the solution of a linear system of equations A*X=RS for  *      
!  a symmetric, cyclic five-diagonal matrix A, once the factor   *      
!  matrizes C and B have been determind in SUBROUTINE NCYFSP.    *      
!  C and B are input parameters and are represented by the       *      
!  vectors DML, ... , RC1U as described in SUBROUTINE NCYFSY.    *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N    :  size of the matrix A                                  *      
!  RS   :  vector RS(1:N); the right hand side                   *      
!                                                                *      
!  DML  :  vector DML(1:N)    ]   describing the entries of the  *      
!  DL1L :  vector DL1L(1:N)   ]   lower tringular matrix C in the*      
!  DL2L :  vector DL2L(1:N)   ]   factorization A=C*B            *      
!  RL2L :  vector RL2L(1:N-4) ]   (refer to NCYFSP)              *      
!  RL1L :  vector RL1L(1:N-3) ]                                  *      
!                                                                *      
!  DU1U :  vector DU1U(1:N)    ]   describing the entries of the *      
!  DU2U :  vector DU2U(1:N)    ]   upper tringular matrix B in   *      
!  RC2U :  vector RC2U(1:N-4)  ]   the factorization A=C*B       *      
!  RC1U :  vector RC1U(1:N-3)  ]   (refer to NCYFSP)             *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  X :  vector  X(1:N); the solution vector                      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author  : GÅnter Palm                                         *      
!  Date    : 04.18.1988                                          *      
!  Source  : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION RS (N), X (N), DML (N), DL1L (N), DL2L (N),       &
      RL2L (N - 4), RL1L (N - 3), DU1U (N), DU2U (N), RC2U (N - 4),     &
      RC1U (N - 3)                                                      
!                                                                       
!-----Solving the system by updating and                                
!     backsubstitution                                                  
!                                                                       
!     1. Updating the right hand side                                   
!                                                                       
      X (1) = RS (1) / DML (1) 
      X (2) = (RS (2) - X (1) * DL1L (2) ) / DML (2) 
      DO 80 I = 3, N - 2, 1 
         X (I) = (RS (I) - X (I - 2) * DL2L (I) - X (I - 1) * DL1L (I) )&
         / DML (I)                                                      
   80 END DO 
      DUMMY1 = 0.0D0 
      DO 90 K = 1, N - 4, 1 
         DUMMY1 = DUMMY1 + X (K) * RL2L (K) 
   90 END DO 
      X (N - 1) = (RS (N - 1) - DUMMY1 - X (N - 3) * DL2L (N - 1)       &
      - X (N - 2) * DL1L (N - 1) ) / DML (N - 1)                        
      DUMMY1 = 0.0D0 
      DO 100 K = 1, N - 3, 1 
         DUMMY1 = DUMMY1 + X (K) * RL1L (K) 
  100 END DO 
      X (N) = (RS (N) - DUMMY1 - X (N - 2) * DL2L (N) - X (N - 1)       &
      * DL1L (N) ) / DML (N)                                            
!                                                                       
!     2. Backsubstitution                                               
!                                                                       
      X (N - 1) = X (N - 1) - DU1U (N - 1) * X (N) 
      X (N - 2) = X (N - 2) - DU1U (N - 2) * X (N - 1) - DU2U (N - 2)   &
      * X (N)                                                           
      X (N - 3) = X (N - 3) - DU1U (N - 3) * X (N - 2) - DU2U (N - 3)   &
      * X (N - 1) - RC1U (N - 3) * X (N)                                
      DO 110 I = N - 4, 1, - 1 
         X (I) = X (I) - DU1U (I) * X (I + 1) - DU2U (I) * X (I + 2)    &
         - RC2U (I) * X (N - 1) - RC1U (I) * X (N)                      
  110 END DO 
      RETURN 
      END SUBROUTINE NCYFSS                         
