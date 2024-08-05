      SUBROUTINE RKSTEP (X, H, Y, N, M, K, DES, YHILO, COEFF, NOSTEP,   &
      FSAL, XZI, IERR)                                                  
!                                                                       
!*****************************************************************      
!                                                                *      
! This program performs one integration step using the chosen    *      
! RUNGE-KUTTA embedding formula. The computed approximations for *      
! the solution Y of the system of differential equations at X + H*      
! are stored in the array YHILO. The first column of YHILO con-  *      
! tains the results for the high order RUNGE-KUTTA method, while *      
! those for the low order method appear in the second column.    *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! X       : DOUBLE PRECISION initial value for the integration   *      
!           step                                                 *      
! H       : DOUBLE PRECISION step size                           *      
! Y       : DOUBLE PRECISION vector Y(1:N), the initial condition*      
!           at X.                                                *      
! N       : number of differential equations in the system,      *      
!           or the size of Y:   0 < N < 13                       *      
! M       : level of the embedding formula, also used for        *      
!           dimensioning COEFF and K.                            *      
! K       : 2-dim. DOUBLE PRECISION array K(1:12,1:M):           *      
!           If an integation step is repeated for an x-value,    *      
!           i.e., if NOSTEP = .TRUE., then the values K(1,1),...,*      
!           K(N,1) must be supplied by the calling program.      *      
! DES     : SUBROUTINE DES must be declared as EXTERNAL in the   *      
!           calling program. DES describes the system of         *      
!           differential equations and must have the following   *      
!           form:                                                *      
!                  SUBROUTINE DES(X,Y,N,YPUNKT)                  *      
!                  DOUBLE PRECISION Y(N),YPUNKT(N),X             *      
!                  YPUNKT(1)=....                                *      
!                  YPUNKT(2)=....                                *      
!                         .                                      *      
!                         .                                      *      
!                         .                                      *      
!                   YPUNKT(N)=....                               *      
!                   RETURN                                       *      
!                   END                                          *      
! COEFF   : 2-dim. DOUBLE PRECISION array COEFF(1:16,1:M) which  *      
!           contains the coefficients of the formula             *      
! NOSTEP  : LOGICAL variable indicating whether a new step is    *      
!           performed (NOSTEP = .FALSE.) or the step is repeated *      
!           with decreased step size.                            *      
! FSAL    : (LOGICAL) variable indicating whether the method     *      
!           FSAL (First Same As Last) is used by the RUNGE-KUTTA *      
!           embedding formula                                    *      
! XZI     : DOUBLE PRECISION largest representable number for    *      
!           testing for OVERFLOW                                 *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! K       : 2-dim. DOUBLE PRECISION array K(1:12,1:M) containing *      
!           the K-values for the integration step                *      
! YHILO   : 2-dim. DOUBLE PRECISION array YHILO(1:12,1:2) con-   *      
!           taining the approximate solution at X + H            *      
! IERR    : error parameter: IERR=0  all is ok                   *      
!                            IERR=1  possible OVERFLOW           *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! I, J, I1: loop variables                                       *      
! XDUMMY  : DOUBLE PRECISION auxiliary variable                  *      
! YDUMMY  : DOUBLE PRECISION vector XDUMMY(1:12)                 *      
! SUM     : DOUBLE PRECISION vector SUM(1:2)                     *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Volker Krger                                      *      
!  Date     : 28.04.1993                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION Y (N), YHILO (12, 2), K (12, M), YDUMMY (12),     &
      COEFF (16, M), SUM (2), X, H, XDUMMY, XZI                         
!                                                                       
! Declare LOGICAL Variable NOSTEP and FSAL                              
!                                                                       
      LOGICAL NOSTEP, FSAL 
!                                                                       
! Initialize IERR                                                       
!                                                                       
      IERR = 0 
! for NOSTEP=.FALSE. we compute the K1 - values for the system          
! of differential equations.                                            
! for FSAL=.TRUE. we use the old K1 - values                            
! else:                                                                 
! We call SUBROUTINE DES, and store the K1 - values                     
! in the first column of K                                              
! In case of detected 'OVERFLOW :                                       
! return to calling program                                             
!                                                                       
      IF (.NOT.NOSTEP) THEN 
         IF (FSAL) THEN 
            DO 10 I = 1, N 
               K (I, 1) = K (I, M) 
   10       END DO 
         ELSE 
            CALL DES (X, Y, N, K (1, 1) ) 
            DO 100 I = 1, N 
               IF (DABS (K (I, 1) ) .GT.XZI) THEN 
                  IERR = 1 
                  RETURN 
               ENDIF 
  100       END DO 
         ENDIF 
      ENDIF 
!                                                                       
! the remaining K - values are computed.                                
! Call SUBROUTINE DES, the corresponding K - values                     
! are in the Ith column of K.                                           
! In case of OVERFLOW return to calling program.                        
!                                                                       
      DO 20 I = 2, M 
         XDUMMY = X + COEFF (I, 1) * H 
         DO 30 I1 = 1, N 
            SUM (1) = 0.0D0 
            DO 40 J = 2, I 
               SUM (1) = SUM (1) + COEFF (I, J) * K (I1, J - 1) 
   40       END DO 
            YDUMMY (I1) = Y (I1) + SUM (1) * H 
   30    END DO 
         CALL DES (XDUMMY, YDUMMY, N, K (1, I) ) 
         DO 110 J = 1, N 
            IF (DABS (K (J, I) ) .GT.XZI) THEN 
               IERR = 1 
               RETURN 
            ENDIF 
  110    END DO 
   20 END DO 
!                                                                       
! Determine approximate solution at X + H                               
!                                                                       
      DO 60 I1 = 1, N 
         SUM (2) = 0.0D0 
         DO 70 J = 1, M 
            SUM (2) = SUM (2) + COEFF (1, J) * K (I1, J) 
   70    END DO 
         SUM (1) = 0.0D0 
         DO 80 J = 3, M 
            SUM (1) = SUM (1) + COEFF (2, J) * K (I1, J - 2) 
   80    END DO 
         SUM (1) = SUM (1) + COEFF (3, 4) * K (I1, M - 1) 
         SUM (1) = SUM (1) + COEFF (3, 5) * K (I1, M) 
!                                                                       
! Approximate solutions of both orders at X + H                         
!                                                                       
         YHILO (I1, 1) = Y (I1) + SUM (1) * H 
         YHILO (I1, 2) = Y (I1) + SUM (2) * H 
   60 END DO 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE RKSTEP                         
