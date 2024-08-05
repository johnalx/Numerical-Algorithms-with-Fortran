      SUBROUTINE K4RORI (USERF, A, B, IP, C, D, IQ, N, CREC, DIVIAT,    &
      WORK, IERR, IUFCLL)                                               
!                                                                       
!*****************************************************************      
!                                                                *      
! Cubature over rectangular regions using the summed ROMBERG-    *      
! RICHARDSON method.                                             *      
!                                                                *      
! By using a simplified summed trapezoidal rule, we can          *      
! approximate the integral of the FUNCTION USERF(X,Y) over       *      
! the rectangle [A,B] x [C,D] for the ROMBERG sequence of step   *      
! sizes. After this we can obtain a better approximation of the  *      
! integral by RICHARDSON extrapolation. For H = B - A or D - C,  *      
! respectively, we use the step sizes:                           *      
!   H * (1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256 ...)      *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! USERF   : user defined FUNCTION USERF(X,Y), whose integral is  *      
!           to be computed.                                      *      
!           The FUNCTION USERF must be declared as EXTERNAL in   *      
!           the calling program.                                 *      
!           The FUNCTION should have the following form:         *      
!                  DOUBLE PRECISION FUNCTION USERF(X,Y)          *      
!                  DOUBLE PRECISION X,Y                          *      
!                         .                                      *      
!                         .                                      *      
!                         .                                      *      
!                  USERF=F(X,Y)                                  *      
!                         .                                      *      
!                         .                                      *      
!                         .                                      *      
!                  RETURN                                        *      
!                  STOP                                          *      
!                                                                *      
! A       : DOUBLE PRECISION left endpoint in X-direction        *      
! B       : DOUBLE PRECISION right endpoint in X-direction       *      
! IP      : INTEGER, the number of intervals in X-direction      *      
! C       : DOUBLE PRECISION lower endpoint in Y-direction       *      
! D       : DOUBLE PRECISION upper endpoint in Y-direction       *      
! IQ      : INTEGER, the number of intervals in Y-direction      *      
! N       : INTEGER, number of summed trapezoidal cubatures,     *      
!           N > 1                                                *      
! WORK    : DOUBLE PRECISION storage vector WORK(0:METHOD+2)     *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! CREC    : DOUBLE PRECISION approximate value for the integral  *      
! DIVIAT  : DOUBLE PRECISION error estimate                      *      
! IERR    : error parameter: IERR=0 all is ok                    *      
!                            IERR=1 number of intervals in       *      
!                                   X-direction erroneous        *      
!                            IERR=2 number of intervals in       *      
!                                   Y-direction erroneous        *      
!                            IERR=3 N <= 1                       *      
!                            IERR=4 integrating over an interval *      
!                                   of length zero.              *      
! IUFCLL  : INTEGER, the number of functional evaluations        *      
!           performed.                                           *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLE:                                                *      
! =================                                              *      
! J   : INTEGER loop variable                                    *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required subroutines: K4ROST, RORIEX                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author  : Volker Krger                                       *      
!  Date    : 06.12.1991                                          *      
!  Source  : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION WORK (0:N - 1, 2), A, B, C, D, CREC, DIVIAT 
      EXTERNAL USERF 
!                                                                       
! Initialize IUFCLL                                                     
!                                                                       
      IUFCLL = 0 
!                                                                       
! Check input data                                                      
!                                                                       
      IF (IP.LT.1) THEN 
         IERR = 1 
         RETURN 
      ELSEIF (IQ.LT.1) THEN 
         IERR = 2 
         RETURN 
      ELSEIF (N.LT.2) THEN 
         IERR = 3 
         RETURN 
      ELSEIF (A.EQ.B.OR.C.EQ.D) THEN 
         IERR = 4 
         RETURN 
      ELSE 
         IERR = 0 
      ENDIF 
!                                                                       
! Execute trapezoidal cubatures for the ROMBERG sequence of step sizes  
!                                                                       
      DO 10 J = 0, N - 1 
         CALL K4ROST (USERF, A, B, IP, C, D, IQ, J, WORK (J, 1),        &
         IUFCLL)                                                        
   10 END DO 
!                                                                       
! Approximate and estimate the error via RICHARDSON extrapolation       
!                                                                       
      CALL RORIEX (WORK (0, 1), WORK (0, 2), N, 2, CREC, DIVIAT) 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE K4RORI                         
!                                                                       
!                                                                       
      SUBROUTINE RORIEX (ROMB, WORK, N, M, VAL, ERREST) 
!                                                                       
!*****************************************************************      
!                                                                *      
! RICHARDSON extrapolation for a given ROMBERG sequence.         *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! ROMB    : DOUBLE PRECISION vector ROMB(0:N-1), containing      *      
!           the ROMBERG sequence                                 *      
! WORK    : DOUBLE PRECISION storage vector WORK(0:N-1)          *      
! N       : INTEGER, the number of elements in ROMB              *      
! M       : INTEGER, the order of the method                     *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! VAL     : DOUBLE PRECISION, final value of the extrapolation   *      
! ERREST  : DOUBLE PRECISION error estimate for VAL              *      
!                                                                *      
!                                                                *      
! LOCALE VARIABLES:                                              *      
! =================                                              *      
! K,J     : INTEGER loop counters                                *      
! S       : DOUBLE PRECISION ] auxiliary                         *      
! P       : DOUBLE PRECISION ]      variables                    *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Volker Krger                                      *      
!  Date     : 06.12.1991                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION ROMB (0:N - 1), WORK (0:N - 1), VAL, ERREST, P, S 
!                                                                       
! Store ROMBERG sequence in auxiliary vector                            
!                                                                       
      DO 10 K = 0, N - 1 
         WORK (K) = ROMB (K) 
   10 END DO 
!                                                                       
! Initialize auxiliary variables S and P                                
!                                                                       
      S = 2.0D0** (DBLE (M) ) 
      P = S 
!                                                                       
! RICHARDSON extrapolation                                              
!                                                                       
      DO 20 K = 1, N - 1 
         DO 30 J = 0, N - K - 1 
            WORK (J) = (P * WORK (J + 1) - WORK (J) ) / (P - 1.0D0) 
   30    END DO 
         P = P * S 
   20 END DO 
!                                                                       
! Set up values for the extrapolation                                   
!                                                                       
      VAL = WORK (0) 
!                                                                       
! Error estimation                                                      
!                                                                       
      ERREST = DABS (WORK (0) - WORK (1) ) 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE RORIEX                         
!                                                                       
!                                                                       
      SUBROUTINE K4ROST (USERF, A, B, IP, C, D, IQ, J, CROST, IUFCLL) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Cubature over rectangular regions using the summed trapezoidal *      
! rule in order to compute the Jth element of a ROMBERG sequence.*      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! USERF   : user defined FUNCTION USERF(X,Y), whose integral is  *      
!           to be computed.                                      *      
!           The FUNCTION USERF must be declared as EXTERNAL in   *      
!           the calling program.                                 *      
!           The FUNCTION should have the following form:         *      
!                  DOUBLE PRECISION FUNCTION USERF(X,Y)          *      
!                  DOUBLE PRECISION X,Y                          *      
!                         .                                      *      
!                         .                                      *      
!                         .                                      *      
!                  USERF=F(X,Y)                                  *      
!                         .                                      *      
!                         .                                      *      
!                         .                                      *      
!                  RETURN                                        *      
!                  STOP                                          *      
!                                                                *      
! A       : DOUBLE PRECISION left endpoint in X-direction        *      
! B       : DOUBLE PRECISION right endpoint in X-direction       *      
! IP      : INTEGER, the number of intervals in X-direction      *      
! C       : DOUBLE PRECISION lower endpoint in Y-direction       *      
! D       : DOUBLE PRECISION upper endpoint in Y-direction       *      
! IQ      : INTEGER, the number of intervals in Y-direction      *      
! J       : INTEGER, the index of the element in the ROMBERG     *      
!           sequence                                             *      
! IUFCLL  : INTEGER, the number of previous function evaluations *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! CROST   : DOUBLE PRECISION approximate value for the integral  *      
! IUFCLL  : INTEGER, the number of function evaluations          *      
!                                                                *      
!                                                                *      
! LOKALE VARIABLEN:                                              *      
! =================                                              *      
! I,K     : INTEGER loop counters                                *      
! HAB     : DOUBLE PRECISION step size in X-direction            *      
! HCD     : DOUBLE PRECISION step size in X-direction            *      
! FAC     : DOUBLE PRECISION weight for the node                 *      
! ABNUM   : DOUBLE PRECISION value for the number of intervals   *      
!           in X-direction                                       *      
! CDNUM   : DOUBLE PRECISION value for the number of intervals   *      
!           in Y-direction                                       *      
! DJ      : DOUBLE PRECISION varible used to determine the number*      
!           of intervals                                         *      
! DBLEI   : DOUBLE PRECISION value for I                         *      
! IABNUM  : INTEGER value for ABNUM                              *      
! ICDNUM  : INTEGER value for CDNUM                              *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Volker Krger                                      *      
!  Date     : 06.12.1991                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION A, B, C, D, CROST, USERF, HAB, HCD, FAC, ABNUM,   &
      CDNUM, DJ, DBLEI                                                  
!                                                                       
! Factor to determine number of intervals                               
!                                                                       
      DJ = 2.0D0** (DBLE (J) ) 
!                                                                       
! Number of intervals in both X- and Y-direction                        
!                                                                       
      ABNUM = DJ * DBLE (IP) 
      CDNUM = DJ * DBLE (IQ) 
      IABNUM = INT (ABNUM) 
      ICDNUM = INT (CDNUM) 
!                                                                       
! Step size in  X- and Y-direction                                      
!                                                                       
      HAB = (B - A) / ABNUM 
      HCD = (D-C) / CDNUM 
                                                                        
!                                                                       
! Initialize CROST                                                      
!                                                                       
      CROST = 0.0D0 
!                                                                       
! Find a value for the integral via the summed trapezoidal rule         
!                                                                       
      DO 10 I = 0, IABNUM 
         DBLEI = DBLE (I) 
         DO 20 K = 0, ICDNUM 
!                                                                       
! Determine weights                                                     
!                                                                       
            FAC = 1.0D0 
            IF (I.GT.0.AND.I.LT.IABNUM) FAC = 2.0D0 
            IF (K.GT.0.AND.K.LT.ICDNUM) FAC = 2.0D0 * FAC 
            CROST = CROST + FAC * USERF (A + HAB * DBLEI, C + HCD *     &
            DBLE (K) )                                                  
            IUFCLL = IUFCLL + 1 
   20    END DO 
   10 END DO 
      CROST = 0.25D0 * HAB * HCD * CROST 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE K4ROST                         
