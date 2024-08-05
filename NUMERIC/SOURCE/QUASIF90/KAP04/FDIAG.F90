![            {Systems with Five--Diagonal Matrices}*)                  
      SUBROUTINE FDIAG (N, DL2, DL1, DM, DU1, DU2, RS, X, MARK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Solving a system of linear equations                       *      
!                      A * X = RS                                *      
!     with a five-diagonal, strongly nonsingular matrix A via    *      
!     Gauss algorithm without pivoting.                          *      
!     The matrix A is given as five N-vectors DL2, DL1, DM, DU1  *      
!     and DU2. The linear system has the form:                   *      
!                                                                *      
!     DM(1)*X(1)+DU1(1)*X(2)+DU2(1)*X(3)             = RS(1)     *      
!     DL1(2)*X(1)+DM(2)*X(2)+DU1(2)*X(3)+DU2(2)*X(4) = RS(2)     *      
!                                                                *      
!     DL2(I)*X(I-2)+DL1(I)*X(I-1)+                               *      
!           +DM(I)*X(I)+DU1(I)*X(I+1)+DU2(I)*X(I+2)  = RS(I)     *      
!            for I = 3, ..., N - 2, and                          *      
!                                                                *      
!     DL2(N-1)*X(N-3)+DL1(N-1)*X(N-2)+                           *      
!             +DM(N-1)*X(N-1)+DU1(N-1)+X(N)          = RS(N-1)   *      
!     DL2(N)*X(N-2)+DL1(N)*X(N-1)+DM(N)*X(N)         = RS(N)     *      
!                                                                *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N     : number of equations; N > 3                         *      
!     DL2   : N-vector DL2(1:N); second lower co-diagonal        *      
!             DL2(3), DL2(4), ... , DL2(N)                       *      
!     DL1   : N-vector DL1(1:N); lower co-diagonal               *      
!             DL1(2), DL1(3), ... , DL1(N)                       *      
!     DM    : N-vector DM(1:N); main diagonal                    *      
!             DM(1), DM(2), ... , DM(N)                          *      
!     DU1   : N-vector DU1(1:N); upper co-diagonal               *      
!             DU1(1), DU1(2), ... , DU1(N-1)                     *      
!     DU2   : N-vector DU2(1:N); second upper co-diagonal        *      
!             DU2(1), DU2(2), ... , DU2(N-2)                     *      
!     RS    : N-vector RS(1:N); the right hand side of the       *      
!             linear system                                      *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     DL2   :) overwritten with auxiliary vectors defining the   *      
!     DL1   :) factorization of the cyclically tridiagonal       *      
!     DM    :) matrix A                                          *      
!     DU1   :)                                                   *      
!     DU2   :)                                                   *      
!     X     : N-vector X(1:N); containing the solution of the    *      
!             the system of equations                            *      
!     MARK  : error parameter                                    *      
!             MARK=-1 : condition N > 3 is not satisfied         *      
!             MARK= 0 : numerically the matrix A is not strongly *      
!                       nonsingular                              *      
!             MARK= 1 : everything is o.k.                       *      
!                                                                *      
!     NOTE: if MARK = 1, the determinant of A is given by:       *      
!                DET A = DM(1) * DM(2) * ... * DM(N)             *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: FDIAGP, FDIAGS, MACHPD                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Gisela Engeln-Muellges                             *      
!  date     : 05.06.1988                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DL1 (1:N), DL2 (1:N), DM (1:N) 
      DOUBLEPRECISION DU1 (1:N), DU2 (1:N), RS (1:N), X (1:N) 
      MARK = - 1 
      IF (N.LT.4) RETURN 
!                                                                       
!  Factor the matrix A                                                  
!                                                                       
      CALL FDIAGP (N, DL2, DL1, DM, DU1, DU2, MARK) 
!                                                                       
!  if MARK = 1, update and bachsubstitute                               
!                                                                       
      IF (MARK.EQ.1) THEN 
         CALL FDIAGS (N, DL2, DL1, DM, DU1, DU2, RS, X) 
      ENDIF 
      RETURN 
      END SUBROUTINE FDIAG                          
!                                                                       
!                                                                       
      SUBROUTINE FDIAGP (N, DL2, DL1, DM, DU1, DU2, MARK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Factor a five-diagonal, strongly nonsingular matrix A      *      
!     that is defined by the five N-vectors DL2, DL1, DM, DU1    *      
!     and DU2, into its triangular factors  L * R  by applying   *      
!     Gaussian elimination specialized for five-diagonal matrices*      
!     (without pivoting).                                        *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N     : number of equations; N > 3                         *      
!     DL2   : N-vector DL2(1:N); second lower co-diagonal        *      
!             DL2(3), DL2(4), ... , DL2(N)                       *      
!     DL1   : N-vector DL1(1:N); lower co-diagonal               *      
!             DL1(2), DL1(3), ... , DL1(N)                       *      
!     DM    : N-vector DM(1:N); main diagonal                    *      
!             DM(1), DM(2), ... , DM(N)                          *      
!     DU1   : N-vector DU1(1:N); upper co-diagonal               *      
!             DU1(1), DU1(2), ... , DU1(N-1)                     *      
!     DU2   : N-vector DU2(1:N); second upper co-diagonal        *      
!             DU2(1), DU2(2), ... , DU2(N-2)                     *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     DL2   :) overwritten with auxiliary vectors that define    *      
!     DL1   :) the factors of the five-diagonal matrix A;        *      
!     DM    :) the three co-diagonals of the lower triangular    *      
!     DU1   :) matrix L are stored in the vectors DL2, DL1 and   *      
!     DU2   :) DM. The two co-diagonals of the unit upper        *      
!              triangular matrix R are stored in the vectors DU1 *      
!              and DU2, its diagonal elements each have the      *      
!              value  1.                                         *      
!     MARK  : error parameter                                    *      
!             MARK=-1 : condition N > 3 is violated              *      
!             MARK= 0 : numerically the matrix is not strongly   *      
!                       nonsingular                              *      
!             MARK= 1 : everything is o.k.                       *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Gisela Engeln-Muellges                             *      
!  date     : 05.06.1988                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DL2 (1:N), DL1 (1:N), DM (1:N), DU1 (1:N),        &
      DU2 (1:N)                                                         
!                                                                       
!  testing whether N > 3                                                
!                                                                       
      MARK = - 1 
      IF (N.LT.4) RETURN 
!                                                                       
!  calculating the machine constant                                     
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = FMACHP * 2.0D0 
!                                                                       
!  determining relative error bounds                                    
!                                                                       
      EPS = 4.0D0 * FMACHP 
!                                                                       
!  initializing the undefined vector components                         
!                                                                       
      DL2 (1) = 0.0D0 
      DL2 (2) = 0.0D0 
      DL1 (1) = 0.0D0 
      DU1 (N) = 0.0D0 
      DU2 (N - 1) = 0.0D0 
      DU2 (N) = 0.0D0 
!                                                                       
!  factoring the matrix A while checking for strong nonsingularity      
!  for N=1, 2                                                           
!                                                                       
      ROW = DABS (DM (1) ) + DABS (DU1 (1) ) + DABS (DU2 (1) ) 
      IF (ROW.EQ.0.0D0) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
      D = 1.0D0 / ROW 
      IF (DABS (DM (1) ) * D.LE.EPS) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
      DU1 (1) = DU1 (1) / DM (1) 
      DU2 (1) = DU2 (1) / DM (1) 
      ROW = DABS (DL1 (2) ) + DABS (DM (2) ) + DABS (DU1 (2) ) + DABS ( &
      DU2 (2) )                                                         
      IF (ROW.EQ.0.0D0) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
      D = 1.0D0 / ROW 
      DM (2) = DM (2) - DL1 (2) * DU1 (1) 
      IF (DABS (DM (2) ) * D.LE.EPS) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
      DU1 (2) = (DU1 (2) - DL1 (2) * DU2 (1) ) / DM (2) 
      DU2 (2) = DU2 (2) / DM (2) 
!                                                                       
!  factoring A while checking for strong nonsingularity of A            
!                                                                       
      DO 20 I = 3, N, 1 
         ROW = DABS (DL2 (I) ) + DABS (DL1 (I) ) + DABS (DM (I) )       &
         + DABS (DU1 (I) ) + DABS (DU2 (I) )                            
         IF (ROW.EQ.0.0D0) THEN 
            MARK = 0 
            RETURN 
         ENDIF 
         D = 1.0D0 / ROW 
         DL1 (I) = DL1 (I) - DL2 (I) * DU1 (I - 2) 
         DM (I) = DM (I) - DL2 (I) * DU2 (I - 2) - DL1 (I) * DU1 (I - 1) 
         IF (DABS (DM (I) ) * D.LE.EPS) THEN 
            MARK = 0 
            RETURN 
         ENDIF 
         IF (I.LT.N) THEN 
            DU1 (I) = (DU1 (I) - DL1 (I) * DU2 (I - 1) ) / DM (I) 
         ENDIF 
         IF (I.LT. (N - 1) ) THEN 
            DU2 (I) = DU2 (I) / DM (I) 
         ENDIF 
   20 END DO 
      MARK = 1 
      RETURN 
      END SUBROUTINE FDIAGP                         
!                                                                       
!                                                                       
      SUBROUTINE FDIAGS (N, DL2, DL1, DM, DU1, DU2, RS, X) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Solving a linear system of equations                       *      
!                A * X = RS                                      *      
!     for a five-diagonal, strongly nonsingular matrix A, once   *      
!     the factor matrices L * R have been calculated by          *      
!     SUBROUTINE FDIAGP. Here they are used as input arrays and  *      
!     they are stored in the five N-vectors DL2, DL1, DM, DU1    *      
!     and DU2.                                                   *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N     : number of equations; N > 3                         *      
!     DL2   : N-vector DL2(1:N); ) lower triangular matrix L     *      
!     DL1   : N-vector DL1(1:N); ) including the diagonal        *      
!     DM    : N-vector DM(1:N);  ) elements                      *      
!                                                                *      
!     DU1   : N-vector DU1(1:N); ) unit upper triangular matrix  *      
!     DU2   : N-vector DU2(1:N); ) R without its unit diagonal   *      
!                                   elements                     *      
!     RS    : N-vector RS1(1:N); right side of the linear system *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     X     : N-vector X(1:N); the solution of the linear system *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Gisela Engeln-Muellges                             *      
!  date     : 05.06.1988                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DL2 (1:N), DL1 (1:N), DM (1:N) 
      DOUBLEPRECISION DU1 (1:N), DU2 (1:N), RS (1:N), X (1:N) 
!                                                                       
!  updating                                                             
!                                                                       
      RS (1) = RS (1) / DM (1) 
      RS (2) = (RS (2) - DL1 (2) * RS (1) ) / DM (2) 
      DO 10 I = 3, N 
         RS (I) = (RS (I) - DL2 (I) * RS (I - 2) - DL1 (I) * RS (I - 1) &
         ) / DM (I)                                                     
   10 END DO 
!                                                                       
!  backsubstitution                                                     
!                                                                       
      X (N) = RS (N) 
      X (N - 1) = RS (N - 1) - DU1 (N - 1) * X (N) 
      DO 20 I = N - 2, 1, - 1 
         X (I) = RS (I) - DU1 (I) * X (I + 1) - DU2 (I) * X (I + 2) 
   20 END DO 
      RETURN 
      END SUBROUTINE FDIAGS                         
