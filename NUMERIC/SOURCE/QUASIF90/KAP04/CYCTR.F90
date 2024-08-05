![             Matrix}                                                  
![            {Systems with a Symmetric Cyclically Tridiagonal          
![             Matrix}*)                                                
      SUBROUTINE CYCTR (N, DL, DM, DU, RL, CR, RS, X, MARK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Solving a linear system of equations                       *      
!                   A * X = RS                                   *      
!     for a cyclically tridiagonal, strongly nonsingular matrix  *      
!     A. The system matrix A is defined via five N-vectors DL,   *      
!     DM, DU, RL and CR. The set of equations has the following  *      
!     form:                                                      *      
!                                                                *      
!     DM(1)*X(1)+DU(1)*X(2)+CR(1)*X(N)       =       RS(1)       *      
!                                                                *      
!     DL(I)*X(I-1)+DM(I)*X(I)+DU(I)*X(I+1)   =       RS(I)       *      
!            for I = 2, ..., N-1, and                            *      
!                                                                *      
!     RL(1)*X(1)+DL(N)*X(N-1)+DM(N)*X(N)     =       RS(N)       *      
!                                                                *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N   : number of equations; N > 2                           *      
!     DL  : N-vector DL(1:N); lower co-diagonal                  *      
!           DL(2), ... , DL(N)                                   *      
!     DM  : N-vector DM(1:N); main diagonal                      *      
!           DM(1), ... , DM(N)                                   *      
!     DU  : N-vector DU(1:N); upper co-diagonal                  *      
!           DU(1), ... , DU(N-1)                                 *      
!     RL  : N-vector RL(1:N); last row RL(1) with diagonal and   *      
!                             co-diagonal elements omitted       *      
!     CR  : N-vector CR(1:N); right most column RS(1) with       *      
!                             diagonal and codiagonal elements   *      
!                             omitted                            *      
!     RS  : N-vector RS(1:N); the right hand side                *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     DL   :) overwritten with auxiliary vectors that define the *      
!     DM   :) factorization matrices of the cyclically           *      
!     DU   :) tridiagonal matrix                                 *      
!     RL   :)                                                    *      
!     CR   :)                                                    *      
!     X    : N-vector X(1:N); the solution of the system of      *      
!            equations                                           *      
!     MARK : error parameter                                     *      
!            MARK=-1 : condition N > 2 is not satified           *      
!            MARK= 0 : numerically the matrix is not strongly    *      
!                      nonsingular                               *      
!            MARK= 1 : everything is o.k.                        *      
!                                                                *      
!                                                                *      
!     NOTE: If MARK = 1,  the determinant of A is given by:      *      
!              DET A = DM(1) * DM(2) * ... * DM(N).              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: CYCTRP, CYCTRS, MACHPD                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author    : Gisela Engeln-Muellges                            *      
!  date      : 05.05.1988                                        *      
!  source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DL (1:N), DM (1:N), DU (1:N), RL (1:N), CR (1:N), &
      RS (1:N), X (1:N)                                                 
      MARK = - 1 
      IF (N.LT.3) RETURN 
!                                                                       
!   factor the matrix A                                                 
!                                                                       
      CALL CYCTRP (N, DL, DM, DU, RL, CR, MARK) 
!                                                                       
!   if MARK = 1, update and bachsubstitute                              
!                                                                       
      IF (MARK.EQ.1) THEN 
         CALL CYCTRS (N, DL, DM, DU, RL, CR, RS, X) 
      ENDIF 
      RETURN 
      END SUBROUTINE CYCTR                          
!                                                                       
!                                                                       
      SUBROUTINE CYCTRP (N, DL, DM, DU, RL, CR, MARK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Factor a cyclically tridiagonal, strongly nonsingular      *      
!     matrix A, that is given by five N-vectors DL, DM, DU, RL   *      
!     and CR, into its factors  L * R  by applying a special     *      
!     Gaussian elimination method. The form of the  system of    *      
!     equations is as described in SUBROUTINE CYCTR.             *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N   : number of equations; N > 2                           *      
!     DL  : N-vector DL(1:N); lower co-diagonal                  *      
!           DL(2), ... , DL(N)                                   *      
!     DM  : N-vector DM(1:N); main diagonal                      *      
!           DM(1), ... , DM(N)                                   *      
!     DU  : N-vector DU(1:N); upper co-diagonal                  *      
!           DU(1), ... , DU(N-1)                                 *      
!     RL  : N-vector RL(1:N); last row RL(1) with diagonal and   *      
!                             co-diagonal elements omitted       *      
!     CR  : N-vector CR(1:N); right most column RS(1) with       *      
!                             diagonal and codiagonal elements   *      
!                             omitted                            *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     DL   :) overwritten with auxiliary vectors containing the  *      
!     DM   :) factors of the matrix A. The lower triangular      *      
!     DU   :) factor L is stored in DL, DM and RL, the unit      *      
!     RL   :) upper triangular matrix R is stored in DU and CR   *      
!     RS   :) with its main diagonal entries ( = 1) omitted.     *      
!     MARK : error parameter                                     *      
!            MARK=-1 : condition N > 2 is not satisfied          *      
!            MARK= 0 : numerically the matrix is not strongly    *      
!                      nonsingular                               *      
!            MARK= 1 : everything is o.k.                        *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author    : Gisela Engeln-Muellges                            *      
!  date      : 05.05.1988                                        *      
!  source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DL (1:N), DM (1:N), DU (1:N), RL (1:N), CR (1:N) 
!                                                                       
!   testing whether N > 2                                               
!                                                                       
      MARK = - 1 
      IF (N.LT.3) RETURN 
!                                                                       
!   computing the machine constant                                      
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = FMACHP * 2.0D0 
!                                                                       
!   determining bounds for the relative error                           
!                                                                       
      EPS = 4.0D0 * FMACHP 
!                                                                       
!   initializing the undefined vector components                        
!                                                                       
      DO 20 I = 2, N 
         RL (I) = 0.0D0 
         CR (I) = 0.0D0 
   20 END DO 
      DL (1) = 0.0D0 
      DU (N) = 0.0D0 
!                                                                       
!   checking for strong nonsingularity of the matrix for N=1            
!                                                                       
      ROW = DABS (DM (1) ) + DABS (DU (1) ) + DABS (CR (1) ) 
      IF (ROW.EQ.0.0D0) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
      D = 1.0D0 / ROW 
      IF (DABS (DM (1) ) * D.LE.EPS) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
!                                                                       
!   factoring the matrix A while checking for                           
!   strong nonsingularity of A                                          
!                                                                       
      DU (1) = DU (1) / DM (1) 
      CR (1) = CR (1) / DM (1) 
      DO 30 I = 2, N - 1 
         ROW = DABS (DL (I) ) + DABS (DM (I) ) + DABS (DU (I) ) 
         IF (ROW.EQ.0.0D0) THEN 
            MARK = 0 
            RETURN 
         ENDIF 
         D = 1.0D0 / ROW 
         DM (I) = DM (I) - DL (I) * DU (I - 1) 
         IF (DABS (DM (I) ) * D.LE.EPS) THEN 
            MARK = 0 
            RETURN 
         ENDIF 
         IF (I.LT. (N - 1) ) THEN 
            DU (I) = DU (I) / DM (I) 
            CR (I) = - DL (I) * CR (I - 1) / DM (I) 
         ENDIF 
   30 END DO 
      ROW = DABS (RL (1) ) + DABS (DL (N) ) + DABS (DM (N) ) 
      IF (ROW.EQ.0.0D0) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
      D = 1.0D0 / ROW 
      DO 40 K = 2, N - 2 
         RL (K) = - RL (K - 1) * DU (K - 1) 
   40 END DO 
      DL (N) = DL (N) - RL (N - 2) * DU (N - 2) 
      DU (N - 1) = (DU (N - 1) - DL (N - 1) * CR (N - 2) ) / DM (N - 1) 
      S = 0.0D0 
      DO 50 J = 1, N - 2 
         S = S + RL (J) * CR (J) 
   50 END DO 
      DM (N) = DM (N) - S - DL (N) * DU (N - 1) 
      IF (DABS (DM (N) ) * D.LE.EPS) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
      MARK = 1 
      RETURN 
      END SUBROUTINE CYCTRP                         
!                                                                       
!                                                                       
      SUBROUTINE CYCTRS (N, DL, DM, DU, RL, CR, RS, X) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Solving a linear system of equations                       *      
!               A * X = RS                                       *      
!     for a cyclically tridiagonal, strongly nonsingular matrix  *      
!     A, once its factors L, R have been calculated in           *      
!     SUBROUTINE CYCTRP.                                         *      
!     The elements of the lower triangular matrix L are stored   *      
!     in the vectors DL, DM, RL, while the elements of the unit  *      
!     upper triangular matrix R (except for the main diagonal)   *      
!     are stored in the vector DU and CR.                        *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N   : number of equations; N > 2                           *      
!     DL  : N-vector DL(1:N); ) these vectors DL, ... , CR       *      
!     DM  : N-vector DM(1:N); ) contain the factors of the       *      
!     RL  : N-vector RL(1:N); ) matrix A.                        *      
!     DU  : N-vector DU(1:N); ) (output vectors of CYCTRP)       *      
!     CR  : N-vector CR(1:N); )                                  *      
!     RS  : N-vector RS(1:N); right hand side                    *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     X   : N-vector X(1:N); the solution of the linear system   *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Gisela Engeln-Muellges                             *      
!  date     : 05.05.1988                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DL (1:N), DM (1:N), DU (1:N), RL (1:N), CR (1:N), &
      RS (1:N), X (1:N)                                                 
!                                                                       
!   updating                                                            
!                                                                       
      RS (1) = RS (1) / DM (1) 
      DO 10 I = 2, N - 1 
         RS (I) = (RS (I) - RS (I - 1) * DL (I) ) / DM (I) 
   10 END DO 
      S = 0.0D0 
      DO 20 J = 1, N - 2 
         S = S + RL (J) * RS (J) 
   20 END DO 
      RS (N) = (RS (N) - S - DL (N) * RS (N - 1) ) / DM (N) 
!                                                                       
!   backsubstitution                                                    
!                                                                       
      X (N) = RS (N) 
      X (N - 1) = RS (N - 1) - X (N) * DU (N - 1) 
      DO 30 I = N - 2, 1, - 1 
         X (I) = RS (I) - DU (I) * X (I + 1) - CR (I) * X (N) 
   30 END DO 
      RETURN 
      END SUBROUTINE CYCTRS                         
