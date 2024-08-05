![  {The Conjugate Gradient Method}                                     
![  {The Conjugate Gradient Method}*)                                   
      SUBROUTINE CG (A, N, IA, Y, X, IERR, D, G, AMULD) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE solves a linear system of equations AX = Y    *      
!  using the conjugate gradient method.                          *      
!                                                                *      
!  ASSUMPTION:                                                   *      
!  ===========                                                   *      
!          A must be a symmetric and positive definite N by N    *      
!          matrix                                                *      
!                                                                *      
!  Input PARAMETERS:                                             *      
!  =================                                             *      
!  A     : 2-dim. array  A(1:IA, 1:N), containing the NxN system *      
!          matrix A. Only the upper triangle of A shall be used  *      
!          and we do not check whether A is indeed symmetric.    *      
!  N     : order of the system                                   *      
!  IA    : leading dimension of A, as specified in the calling   *      
!          program                                               *      
!  Y     : N-vector Y(1:N), the right hand side                  *      
!                                                                *      
!  HILFSPARAMETER:                                               *      
!  ===============                                               *      
!  D     : N-vector D(1:N)                                       *      
!  G     : N-vector G(1:N)                                       *      
!  AMULD : N-vector AMULD(1:N) containing  A*D                   *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  X     : N-vector X(1:N), the solution of the linear system    *      
!  IERR  : error parameter:                                      *      
!            = 0, if the denominator of ALPHA vanishes           *      
!            = 1, all is ok, the solution has been found         *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Gisela Engeln-MÅllges                             *      
!  Date      : 02.12.1991                                        *      
!  Source    : FORTRAN  77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
!  Declarations                                                         
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION A (1:IA, 1:N), Y (1:N), X (1:N), D (1:N), G (1:N),&
      AMULD (1:N)                                                       
!                                                                       
!  Compute the machine constant                                         
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = FMACHP * 0.5D0 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = 8.0D0 * FMACHP 
!                                                                       
!  Initialize the auxiliary vectors  D and G and use                    
!  the zero vector for a start                                          
!                                                                       
      DO 20 I = 1, N 
         HELP = Y (I) 
         D (I) = HELP 
         G (I) = - HELP 
         X (I) = 0.0D0 
   20 END DO 
!                                                                       
!  Perform N conjugate gradient steps                                   
!                                                                       
      DO 30 K = 0, N - 1 
!                                                                       
!       Initialize numerator and denominator for ALPHA                  
!                                                                       
         XNUM = 0.0D0 
         DENOM = 0.0D0 
!                                                                       
!       update ALPHA according to:                                      
!       ALPHA = -(D(TRANSP)*G) /A*D(TRANSP)*(A*D))                      
!                                                                       
         DO 40 I = 1, N 
            XNUM = XNUM + D (I) * G (I) 
            HELP = 0.0D0 
            DO 50 J = 1, I - 1 
               HELP = HELP + A (J, I) * D (J) 
   50       END DO 
            DO 60 J = I, N 
               HELP = HELP + A (I, J) * D (J) 
   60       END DO 
            AMULD (I) = HELP 
            DENOM = DENOM + D (I) * HELP 
   40    END DO 
!                                                                       
!       check whether the denominator of  ALPHA  is zero                
!                                                                       
         IF (DABS (DENOM) .LT.FMACHP) THEN 
            IERR = 0 
            RETURN 
         ENDIF 
         ALPHA = - XNUM / DENOM 
!                                                                       
!       update  X := X + ALPHA * D                                      
!                                                                       
         DO 70 I = 1, N 
            X (I) = X (I) + ALPHA * D (I) 
   70    END DO 
                                                                        
!                                                                       
!       update  G := G + ALPHA * A * D                                  
!       and find its norm: NORM;                                        
!       we also check whether X is a good enough approximation          
!       of the solution so that computations can be stopped             
!       with less than  N  CG-steps.                                    
!                                                                       
         GNORM = 0.0D0 
         DO 80 I = 1, N 
            G (I) = G (I) + ALPHA * AMULD (I) 
            GNORM = GNORM + G (I) * G (I) 
   80    END DO 
         IF (GNORM.LT.FMACHP) THEN 
            IERR = 1 
            RETURN 
         ENDIF 
!                                                                       
!       Calculate a new  BETA :                                         
!                                                                       
!       BETA = (G(TRANSP)*(A*D)) / (D(TRANSP)*(A*D))                    
!                                                                       
         XNUM = 0.0D0 
         DO 90 I = 1, N 
            XNUM = XNUM + G (I) * AMULD (I) 
   90    END DO 
         BETA = XNUM / DENOM 
!                                                                       
!       update  D := -G + BETA * D                                      
!                                                                       
         DO 100 I = 1, N 
            D (I) = - G (I) + BETA * D (I) 
  100    END DO 
   30 END DO 
      IERR = 1 
!                                                                       
      RETURN 
      END SUBROUTINE CG                             
