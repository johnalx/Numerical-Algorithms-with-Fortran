      SUBROUTINE STEP32 (X, H, Y, N, K, DES, YHILO, NOSTEP, XZI, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine performs one integration step using the        *      
! RUNGE-KUTTA embedding formula composed of the RUNGE-KUTTA      *      
! method of third order and the improved EULER-CAUCHY method of  *      
! of second order. The computed approximate solutions for the    *      
! ordinary differential equation system at X + H are stored in   *      
! YHILO: Its first column contains the solution for the third    *      
! order method, while the second one has the second order solu-  *      
! tion.                                                          *      
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
! K       : 2-dim. DOUBLE PRECISION array K(1:12,1:3):           *      
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
! NOSTEP  : LOGICAL variable indicating whether a new step is    *      
!           performed (NOSTEP = .FALSE.) or the step is repeated *      
!           with decreased step size.                            *      
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
! I       : loop variable                                        *      
! XDUMMY  : DOUBLE PRECISION auxiliary variable                  *      
! YDUMMY  : DOUBLE PRECISION vector XDUMMY(1:12)                 *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Volker Krger                                      *      
!  Date     : 07.08.1990                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION Y (N), YHILO (12, 2), K (12, 3), YDUMMY (12),     &
      X, H, XDUMMY, XZI                                                 
!                                                                       
! Initialize LOGICAL Variable NOSTEP                                    
!                                                                       
      LOGICAL NOSTEP 
!                                                                       
! Initialize IERR                                                       
!                                                                       
      IERR = 0 
!                                                                       
! IF NOSTEP=.FALSE., we compute the K1 - values for the                 
! system of differential equations.                                     
! Call SUBROUTINE DES, the K1 - values are in the first column          
! of K.                                                                 
! In case of OVERFLOW return to calling program.                        
!                                                                       
      IF (.NOT.NOSTEP) THEN 
         CALL DES (X, Y, N, K (1, 1) ) 
         DO 100 I = 1, N 
            IF (DABS (K (I, 1) ) .GT.XZI) THEN 
               IERR = 1 
               RETURN 
            ENDIF 
  100    END DO 
      ENDIF 
!                                                                       
! the K2 - values are now being computed.                               
! Call SUBROUTINE DES, the K2 - values are in column two of K.          
! In case of detected OVERFLOW return to calling program.               
!                                                                       
      XDUMMY = X + 0.5D0 * H 
      DO 10 I = 1, N 
         YDUMMY (I) = Y (I) + H * 0.5D0 * K (I, 1) 
   10 END DO 
      CALL DES (XDUMMY, YDUMMY, N, K (1, 2) ) 
      DO 110 I = 1, N 
         IF (DABS (K (I, 2) ) .GT.XZI) THEN 
            IERR = 1 
            RETURN 
         ENDIF 
  110 END DO 
!                                                                       
! Compute the K3 - values.                                              
! Call SUBROUTINE DES, the K3 - values will be in column three          
! of K.                                                                 
! In case of detected OVERFLOW we return to the calling program.        
!                                                                       
      XDUMMY = X + H 
      DO 20 I = 1, N 
         YDUMMY (I) = Y (I) + H * (2.0D0 * K (I, 2) - K (I, 1) ) 
   20 END DO 
      CALL DES (XDUMMY, YDUMMY, N, K (1, 3) ) 
      DO 120 I = 1, N 
         IF (DABS (K (I, 3) ) .GT.XZI) THEN 
            IERR = 1 
            RETURN 
         ENDIF 
  120 END DO 
!                                                                       
! Compute the approximate solutions of third and                        
! second order at X + H.                                                
!                                                                       
      DO 30 I = 1, N 
         YHILO (I, 1) = Y (I) + H * ( (K (I, 1) + K (I, 3) ) / 6.0D0 +  &
         2.0D0 * K (I, 2) / 3.0D0)                                      
         YHILO (I, 2) = Y (I) + H * K (I, 2) 
   30 END DO 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE STEP32                         
