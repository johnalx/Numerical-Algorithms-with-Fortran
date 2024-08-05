      SUBROUTINE SESSPM (AZ, N, KPVT, B) 
!                                                                       
!*****************************************************************      
!     SESSPM solves the symmetric system of equations A*X = B    *      
!     by using the factorization of A produced by SUBROUTINE     *      
!     CEPSPM, ZSPMOK or ZSPMMK.                                  *      
!     Since the right-hand side B has already been updated in    *      
!     SUBROUTINE ZSPM.., only a backsubstitution is required     *      
!     here.                                                      *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     AZ    DOUBLE PRECISION vector AZ(1:N*(N+1)/2) containing   *      
!           the matrix A that was factored by CEPSPM or ZSPM.. . *      
!           If SESSPM is called after CEPSPM or ZSPMMK, AZ       *      
!           is the decomposed matrix WK from ZSPMMK.             *      
!           If SESSPM is called after ZSPMOK, AZ denotes the     *      
!           factored matrix AP from ZSPMOK.                      *      
!     N     dimension of the matrix A.                           *      
!     B     vector B(1:N) containing the right-hand side B of    *      
!           the system of equations  A*X = B.                    *      
!           Output of SUBROUTINE CEPSPM or ZSPM.. .              *      
!     KPVT  INTEGER vector KPVT(1:N) containing the PIVOT        *      
!           indices for the factorization.                       *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     B     DOUBLE PRECISION vector B(1:N), the solution         *      
!           vector X                                             *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!     subroutines required: SCAPRO                               *      
!                                                                *      
!                                                                *      
!     source  : Linpack User's Guide , SIAM Philadelphia, 1979   *      
!                                                                *      
!               the source, available in FORTRAN 4, was          *      
!               converted to FORTRAN 77. Some details were       *      
!               modified and adjusted for the requirements of    *      
!               our calling programs.                            *      
!               This program and the related subroutines are     *      
!               not compatible with the original ones from the   *      
!               Linpack User's Guide.                            *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     authors  :  Michael Groenheim, Ina Hinze                   *      
!     date     :  10.25.1989                                     *      
!     source   :  FORTRAN 77                                     *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER N, KPVT (N) 
      DOUBLEPRECISION AZ (1:N * (N + 1) / 2), B (N), D, D1, D2, TEMP,   &
      SCAPRO                                                            
      INTEGER IK, IKP1, K, KK, KP 
      K = 1 
      IK = 0 
      KK = 1 
   10 IF (K.GT.N) RETURN 
      IF (KPVT (K) .LT.0) THEN 
!                                                                       
!        2 x 2 PIVOT block                                              
!                                                                       
         IF (K.EQ.1) GOTO 20 
         B (K) = B (K) - SCAPRO (K - 1, AZ (IK + 1), B (1) ) 
         IKP1 = IK + K 
         B (K + 1) = B (K + 1) - SCAPRO (K - 1, AZ (IKP1 + 1), B (1) ) 
!                                                                       
!        determine the determinants                                     
!                                                                       
   20    D = AZ (KK) * AZ (KK + K + 1) - AZ (KK + K) * AZ (KK + K) 
         D1 = B (K) * AZ (KK + K + 1) - B (K + 1) * AZ (KK + K) 
         D2 = AZ (KK) * B (K + 1) - AZ (KK + K) * B (K) 
         B (K) = D1 / D 
         B (K + 1) = D2 / D 
         KP = IABS (KPVT (K) ) 
         IF (KP.NE.K) THEN 
!                                                                       
!           swap                                                        
!                                                                       
            TEMP = B (K) 
            B (K) = B (KP) 
            B (KP) = TEMP 
         ENDIF 
         IK = IK + K + K + 1 
         K = K + 2 
         KK = IK + K 
         GOTO 10 
      ELSE 
!                                                                       
!        1 x 1 PIVOT block                                              
!                                                                       
         IF (K.EQ.1) THEN 
            B (K) = B (K) / AZ (K) 
         ELSE 
            B (K) = (B (K) - SCAPRO (K - 1, AZ (IK + 1), B (1) ) )      &
            / AZ (KK)                                                   
            KP = KPVT (K) 
            IF (KP.NE.K) THEN 
!                                                                       
!              swap                                                     
!                                                                       
               TEMP = B (K) 
               B (K) = B (KP) 
               B (KP) = TEMP 
            ENDIF 
         ENDIF 
      ENDIF 
      IK = IK + K 
      K = K + 1 
      KK = IK + K 
      GOTO 10 
      END SUBROUTINE SESSPM                         
