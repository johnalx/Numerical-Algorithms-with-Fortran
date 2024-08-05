![          {Differentiation by the Romberg Method}*)                   
      SUBROUTINE DIFROM (FCT, X0, EPS, N, H, RES, ERREST, NEND, HEND,   &
      IERR, D)                                                          
!                                                                       
!*****************************************************************      
!                                                                *      
!  The subroutine DIFROM approximately computes the first        *      
!  derivative of a given function FCT at X0 by using the         *      
!  ROMBERG-method.                                               *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  FCT      : name of the function that is to be differentiated. *      
!             It must be provided by the user in the form:       *      
!                  DOUBLE PRECISION FUNCTION FCT(X)              *      
!                    ...                                         *      
!                  RETURN                                        *      
!                  END                                           *      
!             in the calling program. It has to be defined as    *      
!             EXTERNAL.                                          *      
!  X0       : location where the first derivative is wanted.     *      
!  EPS      : desired accuracy for the derivative.               *      
!  N        : maximum number of rows or columns of the ROMBERG   *      
!             scheme - 1 ; N has to be > 0.                      *      
!  H        : initial step size; H has to be >= 4 * machine      *      
!             constant.                                          *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  RES      : approximate value for the first derivative of the  *      
!             function FCT at X0.                                *      
!  ERREST   : error estimate for the approximate value RES.      *      
!  NEND     : number of rows or columns of the ROMBERG scheme    *      
!             that was actually used.                            *      
!  HEND     : terminal step size.                                *      
!  EPS      : error bound actually used.                         *      
!  IERR     : error parameter.                                   *      
!             = 0, no error.                                     *      
!             = 1, error in the input data.                      *      
!             = 2, desired accuracy was not reached after N steps*      
!                  ERREST >= EPS.                                *      
!             = 3, step size for the difference quotient became  *      
!                  too small ( < 4 * machine constant).          *      
!                                                                *      
!                                                                *      
!  AUXILIARY  PARAMETERS:                                        *      
!  ======================                                        *      
!  D        : (N+1)-vector D(0:N); storage of the current row    *      
!             in the ROMBERG scheme.                             *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Elmar Pohl                                         *      
!  editor   : Guido Dubois                                       *      
!  date     : 04.25.88                                           *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      DOUBLEPRECISION D (0:N) 
!                                                                       
!  local storage of the minimal error bound EPSMIN in                   
!  case that the subroutine is called repeatedly                        
!                                                                       
      SAVE EPSMIN, IFLAG 
      DATA IFLAG / 0 / 
!                                                                       
!  determine the machine constant and initialize the                    
!  minimal error bound                                                  
!                                                                       
      IF (IFLAG.EQ.0) THEN 
         IFLAG = 1 
         FMACHP = 1.0D0 
   10    FMACHP = 0.5D0 * FMACHP 
         IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
         EPSMIN = 8.0D0 * FMACHP 
      ENDIF 
!                                                                       
!  test the input data                                                  
!                                                                       
      IF (EPS.LT.EPSMIN) EPS = EPSMIN 
      IERR = 1 
      IF (N.LT.1.OR.H.LT.EPSMIN) RETURN 
      IERR = 0 
!                                                                       
!  determine the first central difference quotient                      
!                                                                       
      H2 = 2.0D0 * H 
      H0 = H 
      D (0) = (FCT (X0 + H0) - FCT (X0 - H0) ) / H2 
      K = 0 
   30 K = K + 1 
      D (K) = 0.0D0 
      D0 = D (0) 
      H2 = H0 
      H0 = 0.5D0 * H0 
      IF (H0.LT.EPSMIN) THEN 
         IERR = 3 
      ELSE 
!                                                                       
!  determine the central difference quotient                            
!                                                                       
         D (0) = (FCT (X0 + H0) - FCT (X0 - H0) ) / H2 
         IK = 1 
!                                                                       
!  determine the linear combinations                                    
!                                                                       
         DO 20 J = 1, K 
            IK = IK * 4 
            D1 = D (J) 
            D (J) = (IK * D (J - 1) - D0) / (IK - 1) 
            D0 = D1 
   20    END DO 
!                                                                       
!  check the break-off criteria                                         
!                                                                       
         ERREST = DABS (D (K) - D (K - 1) ) 
         IF (ERREST.GE.EPS) THEN 
            IF (K.LE.N) GOTO 30 
            IERR = 2 
         ENDIF 
      ENDIF 
      RES = D (K) 
      NEND = K + 1 
      HEND = H0 
      RETURN 
      END SUBROUTINE DIFROM                         
