![  {Systems with Trid.\ Symm.\ Strongly Nonsing.\ Matrices}            
![  {Systems with Tridiagonal Symmetric                                 
![   Strongly Nonsingular Matrices}*)                                   
      SUBROUTINE TRDSY (N, DM, DU, RS, X, MARK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Solving a linear system of equations                       *      
!                  A * X = RS                                    *      
!     for a tridiagonal, symmetric, positive definite matrix A.  *      
!     The matrix A is defined by the two N-vectors DM and DU.    *      
!     The system of equations is given as follows:               *      
!                                                                *      
!     DM(1) * X(1) + DU(1) * X(2)                      = RS(1)   *      
!                                                                *      
!     DU(I-1) * X(I-1) + DM(I) * X(I) + DU(I) * X(I+1) = RS(I)   *      
!            for I = 2, ... ,N-1, and                            *      
!                                                                *      
!     DU(N-1) * X(N-1) + DM(N) * X(N)                  = RS(N)   *      
!                                                                *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N    : number of equations, N > 2                          *      
!     DM   : N-vector DM(1:N); main diagonal of A                *      
!            DM(1), DM(2), ... , DM(N)                           *      
!     DU   : N-vector DU(1:N); co-diagonal of A                  *      
!            DU(1), DU(2), ... , DU(N-1)                         *      
!     RS   : N-vector X(1:N); the right hand side                *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     DM   :)                                                    *      
!     DU   :) overwritten with auxiliary vectors                 *      
!     RS   :)                                                    *      
!     X    : N-vector X(1:N) containing the solution of the      *      
!            system of equations                                 *      
!     MARK : error parameter                                     *      
!            MARK= 1 : ok                                        *      
!            MARK= 0 : numerically the matrix A is not strongly  *      
!                      nonsingular                               *      
!            MARK=-1 : A is not positive definite.               *      
!            MARK=-2 : condition N > 2 is not satisfied          *      
!                                                                *      
!     NOTE: If MARK = 1, then the determinant of A can be        *      
!           calculated as:                                       *      
!              DET A = DM(1) * DM(2) * ... * DM(N)               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: TRDSYP, TRDSYS, MACHPD                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Gisela Engeln-Muellges                           *      
!  date       : 25.04.1988                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DM (1:N), DU (1:N), RS (1:N), X (1:N) 
      MARK = - 2 
      IF (N.LT.3) RETURN 
!                                                                       
!  Factoring A                                                          
!                                                                       
      CALL TRDSYP (N, DM, DU, MARK) 
!                                                                       
!  if MARK = 1 update and backsubstitute                                
!                                                                       
      IF (MARK.EQ.1) THEN 
         CALL TRDSYS (N, DM, DU, RS, X) 
      ENDIF 
      RETURN 
      END SUBROUTINE TRDSY                          
!                                                                       
!                                                                       
      SUBROUTINE TRDSYP (N, DM, DU, MARK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Factoring a tridiagonal, symmetric, and positive definite  *      
!     matrix A, that is given by the two N-vectors DM and DU,    *      
!     into the product A = R(TRANSP) * D * R  for a unit upper   *      
!     triangular matrix R by applying the Cholesky-method for    *      
!     tridiagonal matrices. The form of the system matrix A is   *      
!     identical with the one desribed in SUBROUTINE TRDSY.       *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N    : number of equations, N > 2                          *      
!     DM   : N-vector DM(1:N); main diagonal of A                *      
!            DM(1), DM(2), ... , DM(N)                           *      
!     DU   : N-vector DU(1:N); co-diagonal of A                  *      
!            DU(1), DU(2), ... , DU(N-1);                        *      
!            due to symmetry of A its lower and upper            *      
!            co-diagonals coincide                               *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     DM   :) overwritten with auxiliary vectors containing the  *      
!     DU   :) factors of A. The co-diagonal of the unit upper    *      
!             triangular bidiagonal matrix R is stored in DU,    *      
!             while the diagonal matrix D is stored in DM.       *      
!     MARK : error parameter                                     *      
!            MARK= 1 : ok                                        *      
!            MARK= 0 : numerically the matrix A is not strongly  *      
!                      nonsingular.                              *      
!            MARK=-1 : A is not positive definite.               *      
!            MARK=-2 : condition N > 2 is not met.               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Gisela Engeln-Muellges                           *      
!  date       : 25.04.1988                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DM (1:N), DU (1:N) 
!                                                                       
!   calculating the machine constant                                    
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = FMACHP * 2.0D0 
!                                                                       
!   determining the relative error bound                                
!                                                                       
      EPS = 4.0D0 * FMACHP 
!                                                                       
!   checking whether N > 2                                              
!                                                                       
      MARK = - 2 
      IF (N.LT.3) RETURN 
      DU (N) = 0.0D0 
!                                                                       
!   testing for a positive definite, strong nonsingular matrix A        
!   for N=1                                                             
!                                                                       
      ROW = DABS (DM (1) ) + DABS (DU (1) ) 
      IF (ROW.EQ.0.0D0) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
      D = 1.0D0 / ROW 
      IF (DM (1) .LT.0.0D0) THEN 
         MARK = - 1 
         RETURN 
      ELSEIF (DABS (DM (1) ) * D.LE.EPS) THEN 
         MARK = 0 
         RETURN 
      ENDIF 
!                                                                       
!   factoring A while checking for a positive definite, strong          
!   nonsingular matrix A                                                
!                                                                       
      DUMMY = DU (1) 
      DU (1) = DU (1) / DM (1) 
      DO 20 I = 2, N, 1 
         ROW = (DABS (DM (I) ) + DABS (DU (I) ) + DABS (DUMMY) ) 
         IF (ROW.EQ.0.0D0) THEN 
            MARK = 0 
            RETURN 
         ENDIF 
         D = 1.0D0 / ROW 
         DM (I) = DM (I) - DUMMY * DU (I - 1) 
         IF (DM (I) .LT.0.0D0) THEN 
            MARK = - 1 
            RETURN 
         ELSEIF (DABS (DM (I) ) * D.LE.EPS) THEN 
            MARK = 0 
            RETURN 
         ENDIF 
         IF (I.LT.N) THEN 
            DUMMY = DU (I) 
            DU (I) = DU (I) / DM (I) 
         ENDIF 
   20 END DO 
      MARK = 1 
      RETURN 
      END SUBROUTINE TRDSYP                         
!                                                                       
!                                                                       
      SUBROUTINE TRDSYS (N, DM, DU, RS, X) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Solving a linear system of equations                       *      
!                  A * X = RS                                    *      
!     for a tridiagonal, symmetric, and positive definite matrix *      
!     A, whose tridiagonal factors have been calculated by the   *      
!     SUBROUTINE TRDSYP. Here the factoring matrices D and R are *      
!     used as input matrices and they are stored in the two      *      
!     N-vectors DM and DU, respectively.                         *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N    : number of equations, N > 2                          *      
!     DM   : N-vector DM(1:N); the diagonal matrix D             *      
!     DU   : N-vector DU(1:N); the upper co-diagonal entries of R*      
!     RS   : N-vector RS(1:N); the right hand side               *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETER:                                          *      
!     =================                                          *      
!     X    : N-vector X(1:N) containing the solution of the      *      
!            system of equations                                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: TRDSYP, TRDSYS                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Gisela Engeln-Muellges                           *      
!  date       : 25.04.1988                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DM (1:N), DU (1:N), RS (1:N), X (1:N) 
!                                                                       
!  updating                                                             
!                                                                       
      DUMMY = RS (1) 
      RS (1) = DUMMY / DM (1) 
      DO 10 I = 2, N, 1 
         DUMMY = RS (I) - DU (I - 1) * DUMMY 
         RS (I) = DUMMY / DM (I) 
   10 END DO 
!                                                                       
!  backsubstitution                                                     
!                                                                       
      X (N) = RS (N) 
      DO 20 I = N - 1, 1, - 1 
         X (I) = RS (I) - DU (I) * X (I + 1) 
   20 END DO 
      RETURN 
      END SUBROUTINE TRDSYS                         
