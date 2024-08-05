      SUBROUTINE IRKCOE (MMAX, LUN, C, A, ALPHA, BETA) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the implicit   *      
! RUNGE-KUTTA methods (IRKM) of order 1 up to MMAX, as specified *      
! in the calling program.                                        *      
! The results are stored unformatted in an external file with    *      
! the logical number LUN. There they can be called up and pro-   *      
! cessed further by the SUBROUTINE IMRUKU.                       *      
! For each order the GAUSS-LEGENDRE nodes ALPHA(J), J=1, ..., M, *      
! of the interval of integration are determined first.           *      
! These are determined from the zeros of the LEGENDRE polynomials*      
! The coefficients BETA(I,J) and A(J), I,J=1, ..., M, are ob-    *      
! tained as the solution of a M*(M+1) linear system of equations.*      
! The solution of the linear system of equations can be deter-   *      
! mined by multiplying LAGRANGE polynomials.                     *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! MMAX   : maximum order up to which the coefficients of the     *      
!          IRKM are to be created                                *      
! LUN    : number of the output file, in which the coefficients  *      
!          are to be stored (unformatted)                        *      
! ALPHA  : vector ALPHA(1:MMAX);             ) the coefficients  *      
! BETA   : 2-dim. array BETA(1:MMAX,1:MMAX); ) of the IRKM       *      
! A      : vector A(1:MMAX);                 )                   *      
! C      : vector C(0:MMAX); auxiliary vector for GALE0 and      *      
!          vector used to store the coefficients of the          *      
!          LAGRANGE polynomials. Only used for storage space.    *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETER:                                              *      
! =================                                              *      
! n o n e                                                        *      
!                                                                *      
!                                                                *      
! All results, i.e., the orders and the coefficients of the IRKM,*      
! are saved unformatted in the file with the logic number LUN.   *      
! This file is to be read by SUBROUTINE IMRUKU.                  *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! M      : current order for which the coefficients are          *      
!          determined.                                           *      
! MM1    : auxiliary variable for M-1                            *      
! MM2    : auxiliary variable for M-2                            *      
! I,J,K  : control variables                                     *      
! JM1    : auxiliary variable for J-1                            *      
! JP1    : auxiliary variable for J+1                            *      
! NG     : counter for the number of factors                     *      
!          (ALPHA(K)-ALPHA(J)) (J constant,K=1,...,J-1,J+1,...,M)*      
!          of the LAGRANGE-polynomials that are multiplied.      *      
! ZJ     : counter of the J-th LAGRANGE-polynomial               *      
! BETAJK : ) auxiliary variables for determining the coefficients*      
! ALPHAK : ) BETA(J,K) or ALPHA(K) when multiplying  the         *      
!          ) LAGRANGE polynomials                                *      
! FLAG   : logic variable; input parameter of GALE0              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: GALE0                                   *      
!                                                                *      
!                                                                *      
!  source : 1. W. Glasmacher, D. Sommer, see [GLAS66].           *      
!                                                                *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Thomas Eul                                         *      
!  date     : 09.17.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
!     parameters                                                        
      INTEGER MMAX, LUN 
      DIMENSION C (0:MMAX), A (MMAX), ALPHA (MMAX), BETA (MMAX, MMAX) 
!                                                                       
!     local variables                                                   
      INTEGER M, MM1, MM2, I, J, JM1, JP1, K, NG 
      DOUBLEPRECISION ZJ, BETAJK, ALPHAK 
      LOGICAL FLAG 
!                                                                       
!*****************************************************************      
!* the weights A(J) and BETA(J,L), as well as the nodes          *      
!* ALPHA(J) are produced for the orders of 1 to MMAX and they    *      
!* are stored, unformatted, in the file numbered LUN.            *      
!*****************************************************************      
!                                                                       
      FLAG = .FALSE. 
      DO 2000 M = 1, MMAX 
         MM1 = M - 1 
         MM2 = M - 2 
!                                                                       
!*****************************************************************      
!*    G A U S S - L E G E N D R E   n o d e s                    *      
!*****************************************************************      
!        determine all zeros of the LEGENDRE polynomials, i.e.,         
!        the ALPHA(J). These are are all real and lie                   
!        symmetrically in the interval  -1. <= ALPHA(J) <= 1.           
!                                                                       
         IF (M.GT.1) THEN 
            CALL GALE0 (M, FLAG, ALPHA, C) 
         ELSE 
            ALPHA (1) = 0.0D0 
         ENDIF 
!                                                                       
!        transform the ALPHA(J) to the interval                         
!        0 <= alpha(j) <= 1                                             
!                                                                       
         DO 10 I = 1, M 
            ALPHA (I) = 0.5D0 * ALPHA (I) + 0.5D0 
   10    END DO 
!                                                                       
!*****************************************************************      
!*      determine the weights  BETA(J,K)  and  A(J)              *      
!*****************************************************************      
!        counter ZJ of the J-th LAGRANGE polynomial of degree M         
!                                                                       
         DO 1000 J = 1, M 
            JM1 = J - 1 
            JP1 = J + 1 
            ZJ = 1.0D0 
            DO 20 K = 1, JM1 
               ZJ = (ALPHA (J) - ALPHA (K) ) * ZJ 
   20       END DO 
            DO 30 K = JP1, M 
               ZJ = (ALPHA (J) - ALPHA (K) ) * ZJ 
   30       END DO 
!                                                                       
!           determine the coefficient of the J-th                       
!           LAGRANGE-polynomial of degree M                             
!                                                                       
            C (0) = 1.0D0 
            NG = 0 
            DO 60 K = 1, JM1 
               ALPHAK = - ALPHA (K) 
               DO 40 I = NG, 0, - 1 
                  C (I + 1) = C (I) 
   40          END DO 
               C (0) = ALPHAK * C (1) 
               DO 50 I = 1, NG 
                  C (I) = C (I) + ALPHAK * C (I + 1) 
   50          END DO 
               NG = NG + 1 
   60       END DO 
            DO 90 K = JP1, M 
               ALPHAK = - ALPHA (K) 
               DO 70 I = NG, 0, - 1 
                  C (I + 1) = C (I) 
   70          END DO 
               C (0) = ALPHAK * C (1) 
               DO 80 I = 1, NG 
                  C (I) = C (I) + ALPHAK * C (I + 1) 
   80          END DO 
               NG = NG + 1 
   90       END DO 
!                                                                       
!           determine the BETA(J,L) and all A(J)                        
!                                                                       
            ZJ = 1.0D0 / ZJ 
            AJ = 0.0D0 
            DO 110 K = 1, M 
               BETAJK = 0.0D0 
               DO 100 L = 1, M 
                  BETAJK = BETAJK + C (L - 1) * ALPHA (K) **L / DBLE (L) 
  100          END DO 
               BETA (J, K) = BETAJK * ZJ 
               AJ = AJ + C (K - 1) / DBLE (K) 
  110       END DO 
            A (J) = AJ * ZJ 
 1000    END DO 
!                                                                       
!        unformatted output of the order combined with the corres-      
!        ponding coefficients to the file numbered LUN                  
!                                                                       
         WRITE (LUN) M, (ALPHA (I), I = 1, M), ( (BETA (I, J), I = 1, M)&
         , J = 1, M), (A (I), I = 1, M)                                 
 2000 END DO 
      RETURN 
      END SUBROUTINE IRKCOE                         
