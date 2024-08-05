![            {Systems with Tridiagonal Matrices}*)                     
      SUBROUTINE TRDIG (N, DL, DM, DU, RS, X, MARK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     This routine solves a linear system of equations           *      
!                  A * X = RS                                    *      
!     for a tridiagonal, strongly nonsingular matrix A. The      *      
!     matrix is given by the three vectors DL, DM and DU which   *      
!     designate the lower co-diagonal, the diagonal and the      *      
!     upper co-diagonal elements of A, respectively.             *      
!     The system of equations has the form :                     *      
!                                                                *      
!     DM(1) * X(1) + DU(1) * X(2)                      = RS(1)   *      
!                                                                *      
!     DL(I) * X(I-1) + DM(I) * X(I) + DU(I) * X(I+1)   = RS(I)   *      
!            for I = 2, ..., N-1, and                            *      
!                                                                *      
!     DL(N) * X(N-1) + DM(N) * X(N)                    = RS(N)   *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N    : number of equations, N > 2                          *      
!     DL   : N-vector DL(1:N); lower co-diagonal of A            *      
!            DL(2), DL(3), ... ,DL(N)                            *      
!     DM   : N-vector DM(1:N); the diagonal of A                 *      
!            DM(1), DM(2), ... , DM(N)                           *      
!     DU   : N-vector DU(1:N); upper co-diagonal of A            *      
!            DU(1), DU(2), ... , DU(N-1)                         *      
!     RS   : N-vector RS(1:N); the right hand side               *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     DL   :)                                                    *      
!     DM   :)                                                    *      
!     DU   :) these are overwritten with auxiliary vectors       *      
!     RS   :)                                                    *      
!     X    : N-vector X(1:N), containing the solution of the     *      
!            system of equations                                 *      
!     MARK : error parameter                                     *      
!            MARK= 1 : everything is o.k.                        *      
!            MARK= 0 : the matrix A is not strongly nonsingular  *      
!            MARK=-1 : error on N: N <= 2                        *      
!                                                                *      
!     NOTE: if MARK = 1, the determinant of A can be calculated  *      
!           after this subroutine has run as follows:            *      
!              DET A = DM(1) * DM(2) * ... * DM(N)               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: TRDIGP, TRDIGS, MACHPD                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Gisela Engeln-Muellges                           *      
!  date       : 05.02.1988                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DL (1:N), DM (1:N), DU (1:N), RS (1:N), X (1:N) 
      MARK = - 1 
      IF (N.LT.3) RETURN 
!                                                                       
!  Factoring the matrix A                                               
!                                                                       
      CALL TRDIGP (N, DL, DM, DU, MARK) 
!                                                                       
!  If MARK = 1, update the right hand side and solve via backsubstitutio
!                                                                       
      IF (MARK.EQ.1) THEN 
         CALL TRDIGS (N, DL, DM, DU, RS, X) 
      ENDIF 
      RETURN 
      END SUBROUTINE TRDIG                          
!                                                                       
!                                                                       
      SUBROUTINE TRDIGP (N, DL, DM, DU, MARK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Factor the tridiagonal, strongly nonsingular matrix A,     *      
!     that is given by the vectors DL, DM and DU of co-diagonal  *      
!     and diagonal entries into the product  L * R  for a bi-    *      
!     diagonal lower triangular matrix  L  and a unit bidiagonal *      
!     upper triangular matrix R. The form of the system matrix A *      
!     is identical to that given in SUBROUTINE TRDIG.            *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N    : number of equations, N > 2                          *      
!     DL   : N-vector DL(1:N); lower co-diagonal of A            *      
!            DL(2), DL(3), ... ,DL(N)                            *      
!     DM   : N-vector DM(1:N); the diagonal of A                 *      
!            DM(1), DM(2), ... , DM(N)                           *      
!     DU   : N-vector DU(1:N); upper co-diagonal of A            *      
!            DU(1), DU(2), ... , DU(N-1)                         *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     DL   : the lower co-diagonal remains unchanged, since it   *      
!            forms the lower co-diagonal of L                    *      
!     DM   :) overwritten with auxiliary vectors that contain    *      
!     DU   :) the factors of A: The co-diagonal of the unit      *      
!             bidiagonal upper tridiagonal matrix R is stored in *      
!             DU. The diagonal elements of L are stored in DM    *      
!     MARK : error parameter                                     *      
!            MARK= 1 : everything is o.k.                        *      
!            MARK= 0 : the matrix is not strongly nonsingular    *      
!            MARK=-1 : condition N > 2 is not fulfilled          *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required subroutine: MACHPD                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Gisela Engeln-Muellges                           *      
!  date       : 05.02.1988                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DL (1:N), DM (1:N), DU (1:N) 
!                                                                       
!   testing whether N > 2                                               
!                                                                       
      MARK = - 1 
      IF (N.LT.3) RETURN 
!                                                                       
!   calculate the machine constant                                      
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = FMACHP * 2.0D0 
!                                                                       
!   determine the relative error bound                                  
!                                                                       
      EPS = 4.0D0 * FMACHP 
!                                                                       
!   checking for strong nonsingularity with N=1                         
!                                                                       
      ROW = DABS (DM (1) ) + DABS (DU (1) ) 
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
!   factoring A while checking for strong nonsingularity                
!                                                                       
      DL (1) = 0.0D0 
      DU (N) = 0.0D0 
      DU (1) = DU (1) / DM (1) 
      DO 20 I = 2, N, 1 
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
         IF (I.LT.N) THEN 
            DU (I) = DU (I) / DM (I) 
         ENDIF 
   20 END DO 
      MARK = 1 
      RETURN 
      END SUBROUTINE TRDIGP                         
!                                                                       
!                                                                       
      SUBROUTINE TRDIGS (N, DL, DM, DU, RS, X) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Solving a linear system of equations                       *      
!                  A * X = RS                                    *      
!     for a tridiagonal, strongly nonsingular matrix A, once the *      
!     tringular factorization has been calculated by SUBROUTINE  *      
!     TRDIGP. The factors are used as known input matrices and   *      
!     are stored in the three N-vectors DL, DM and DU.           *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N    : number of equations, N > 2                          *      
!     DL   : N-vector DL(1:N); ) the vectors DL, DM and DU con-  *      
!     DM   : N-vector DM(1:N); ) tain the factors of the matrix  *      
!     DU   : N-vector DU(1:N); ) A obtained as output of TRDIGP  *      
!     RS   : N-vector RS(1:N); right hand side                   *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     X    : N-vector X(1:N); containing the solution of the     *      
!            system of equations                                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Gisela Engeln-Muellges                           *      
!  date       : 05.02.1988                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION DL (1:N), DM (1:N), DU (1:N), RS (1:N), X (1:N) 
!                                                                       
!  updating                                                             
!                                                                       
      RS (1) = RS (1) / DM (1) 
      DO 10 I = 2, N, 1 
         RS (I) = (RS (I) - DL (I) * RS (I - 1) ) / DM (I) 
   10 END DO 
!                                                                       
!  backsubstitution                                                     
!                                                                       
      X (N) = RS (N) 
      DO 20 I = N - 1, 1, - 1 
         X (I) = RS (I) - DU (I) * X (I + 1) 
   20 END DO 
      RETURN 
      END SUBROUTINE TRDIGS                         
