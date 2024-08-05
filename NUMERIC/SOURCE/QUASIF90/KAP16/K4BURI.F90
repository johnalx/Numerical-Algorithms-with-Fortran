![  {Romberg Cubature for Rectangles and Triangles}                     
![  {Romberg Cubature for Rectangles and Triangles}*)                   
      SUBROUTINE K4BURI (USERF, A, B, IP, C, D, IQ, N, CREC, DIVIAT,    &
      WORK, IERR, IUFCLL)                                               
!                                                                       
!*****************************************************************      
!                                                                *      
! Cubature over rectangular regions using the summed BULIRSCH-   *      
! RICHARDSON method.                                             *      
!                                                                *      
! We use a simplified summed two-dimensional trapezoidal rule    *      
! for the BULIRSCH sequence of step sizes in order to find an    *      
! approximation for the integral of the FUNCTION USERF(X,Y)      *      
! over the rectangle [A,B] x [C,D]. Using RICHARDSON extrapola-  *      
! tion we then compute an improved value for the integral.       *      
! The step sizes follow the BULIRSCH sequences for the given     *      
! interval length H and are:                                     *      
!   H * (1/2, 1/3, 1/4, 1/6, 1/8, 1/12, 1/16, 1/24, 1/32 ...)    *      
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
! A       : DOUBLE PRECISION left hand endpoint in X-direction   *      
! B       : DOUBLE PRECISION right hand endpoint in X-direction  *      
! IP      : INTEGER, the number of sub-intervals in X-direction  *      
! C       : DOUBLE PRECISION lower endpoint in Y-direction       *      
! D       : DOUBLE PRECISION upper endpoint in Y-direction       *      
! IQ      : INTEGER, the number of sub-intervals in Y-direction  *      
! N       : INTEGER, the number of summed trapezoidal cubatures  *      
!           that are to be executed  ( N > 1 )                   *      
! WORK    : 2-dimensional DOUBLE PRECISION array WORK(0:N-1,1:2) *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! CREC    : DOUBLE PRECISION value for the integral              *      
! DIVIAT  : DOUBLE PRECISION error estimate                      *      
! IERR    : error parameter: IERR=0 all is ok                    *      
!                            IERR=1 number of intervals in       *      
!                                   X-direction erroneous        *      
!                            IERR=2 number of intervals in       *      
!                                   Y-direction erroneous        *      
!                            IERR=3 N <= 1                       *      
!                            IERR=4 interval of integration has  *      
!                                   length zero                  *      
! IUFCLL  : INTEGER, the number of functional evaluations        *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                              *       
! =================                                              *      
! J       : INTEGER, loop counter                                *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: K4BUST, BURIEX                          *      
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
! Perform cubature using trapezoidal rule for the                       
! BULIRSCH sequence of step sizes                                       
!                                                                       
      DO 10 J = 0, N - 1 
         CALL K4BUST (USERF, A, B, IP, C, D, IQ, J, WORK (J, 1),        &
         IUFCLL)                                                        
   10 END DO 
!                                                                       
! Find a better approximation and an error estimate                     
! for the integral using RICHARDSON extrapolation                       
!                                                                       
      CALL BURIEX (WORK (0, 1), WORK (0, 2), N, 2, CREC, DIVIAT) 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE K4BURI                         
!                                                                       
!                                                                       
      SUBROUTINE K4BUST (USERF, A, B, IP, C, D, IQ, J, CBUST, IUFCLL) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Cubature over rectangular regions using the summed trapezoidal *      
! rule and computation of the Jth term of the BULIRSCH sequence. *      
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
! A       : DOUBLE PRECISION left hand endpoint in X-direction   *      
! B       : DOUBLE PRECISION right hand endpoint in X-direction  *      
! IP      : INTEGER, the number of sub-intervals in X-direction  *      
! C       : DOUBLE PRECISION lower endpoint in Y-direction       *      
! D       : DOUBLE PRECISION upper endpoint in Y-direction       *      
! IQ      : INTEGER, the number of sub-intervals in Y-direction  *      
! J       : INTEGER, the index of the term of the BULIRSCH       *      
!           sequence                                             *      
! IUFCLL  : INTEGER, the number of functional evaluations of     *      
!           previous calls                                       *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! CBUST   : DOUBLE PRECISION value for the integral              *      
! IUFCLL  : INTEGER, the number of functional evaluations        *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! =================                                              *      
! I,K     : INTEGERS, loop counters                              *      
! HAB     : DOUBLE PRECISION step size in X-direction            *      
! HCD     : DOUBLE PRECISION step size in Y-direction            *      
! FAC     : DOUBLE PRECISION weights for the nodes               *      
! ABNUM   : DOUBLE PRECISION number of intervals in X-direction  *      
! CDNUM   : DOUBLE PRECISION number of intervals in Y-direction  *      
! PQNUM   : DOUBLE PRECISION factor used to compute the number   *      
!           of intervals                                         *      
! DBLEI   : DOUBLE PRECISION value for I                         *      
! IABNUM  : INTEGER value for ABNUM                              *      
! ICDNUM  : INTEGER value for CDNUM                              *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required: DENOM                                   *      
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
      DOUBLEPRECISION A, B, C, D, CBUST, USERF, HAB, HCD, FAC, DENOM,   &
      PQNUM, ABNUM, CDNUM, DBLEI                                        
!                                                                       
! Find the factor for the number of intervals                           
!                                                                       
      PQNUM = DENOM (J) 
!                                                                       
! Number of intervals in X- and Y-directions                            
!                                                                       
      ABNUM = PQNUM * DBLE (IP) 
      CDNUM = PQNUM * DBLE (IQ) 
      IABNUM = INT (ABNUM) 
      ICDNUM = INT (CDNUM) 
!                                                                       
! Step sizes in X- and Y-directions                                     
!                                                                       
      HAB = (B - A) / ABNUM 
      HCD = (D-C) / CDNUM 
!                                                                       
! Initialize CBUST                                                      
!                                                                       
      CBUST = 0.0D0 
!                                                                       
! Find approximate value for the integral by                            
! using the summed trapezoidal rule                                     
!                                                                       
      DO 10 I = 0, IABNUM 
         DBLEI = DBLE (I) 
         DO 20 K = 0, ICDNUM 
!                                                                       
! Determine weights for the nodes                                       
!                                                                       
            FAC = 1.0D0 
            IF (I.GT.0.AND.I.LT.IABNUM) FAC = 2.0D0 
            IF (K.GT.0.AND.K.LT.ICDNUM) FAC = 2.0D0 * FAC 
            CBUST = CBUST + FAC * USERF (A + HAB * DBLEI, C + HCD *     &
            DBLE (K) )                                                  
            IUFCLL = IUFCLL + 1 
   20    END DO 
   10 END DO 
      CBUST = 0.25D0 * HAB * HCD * CBUST 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE K4BUST                         
!                                                                       
!                                                                       
      SUBROUTINE BURIEX (BULI, WORK, N, M, VAL, ERREST) 
!                                                                       
!*****************************************************************      
!                                                                *      
! RICHARDSON extrapolation for a given BULIRSCH sequence         *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! BULI    : DOUBLE PRECISION vector BULI(0:N-1), containing the  *      
!           BULIRSCH sequence                                    *      
! WORK    : DOUBLE PRECISION vector WORK(0:N-1)                  *      
! N       : INTEGER size of the vector BULI                      *      
! M       : INTEGER order of the method                          *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! VAL     : DOUBLE PRECISION terminal value of the extrapolation *      
! ERREST  : DOUBLE PRECISION error estimate for VAL              *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! K,J     : INTEGERS loop counters                               *      
! P       : DOUBLE PRECISION auxiliary variable                  *      
! DM      : DOUBLE PRECISION value for M                         *      
! DK      : DOUBLE PRECISION value for K                         *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: DENOM                                   *      
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
      DOUBLEPRECISION BULI (0:N - 1), WORK (0:N - 1), VAL, ERREST, P,   &
      DENOM, DM, DK                                                     
!                                                                       
! Change INTEGER value to DOUBLE PRECISION                              
!                                                                       
      DM = DBLE (M) 
!                                                                       
! Store the BULIRSCH sequence in WORK                                   
!                                                                       
      DO 10 K = 0, N - 1 
         WORK (K) = BULI (K) 
   10 END DO 
!                                                                       
! RICHARDSON extrapolation                                              
!                                                                       
      DO 20 K = 1, N - 1 
         DK = DBLE (K) 
         DO 30 J = 0, N - K - 1 
            P = (DENOM (J + K) / DENOM (J) ) ** (DM * DK) 
            WORK (J) = (P * WORK (J + 1) - WORK (J) ) / (P - 1.0D0) 
   30    END DO 
   20 END DO 
!                                                                       
! Store value of the extrapolation                                      
!                                                                       
      VAL = WORK (0) 
!                                                                       
! Error estimate                                                        
!                                                                       
      ERREST = ABS (WORK (0) - WORK (1) ) 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE BURIEX                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION DENOM (J) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Determine the denominator for the step size in order to obtain *      
! the Jth element of the BULIRSCH sequence                       *      
!                                                                *      
!                                                                *      
! Input PARAMETER:                                               *      
! ================                                               *      
! J       : INTEGER, the index of the element of the BULIRSCH    *      
!           sequence                                             *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author  : Volker Krger                                       *      
!  Date    : 06.12.1991                                          *      
!  Source  : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Determine denominator DENOM                                           
!                                                                       
      IF (J.EQ.0) THEN 
         DENOM = 1.0D0 
         RETURN 
      ENDIF 
      IF (MOD (J, 2) .EQ.0) THEN 
!                                                                       
! For even J                                                            
!                                                                       
         DENOM = 3.0D0 * 2.0D0** ( (DBLE (J) - 2.0D0) * 0.5D0) 
      ELSE 
!                                                                       
! For odd J                                                             
!                                                                       
         DENOM = 2.0D0** ( (DBLE (J) + 1.0D0) * 0.5D0) 
      ENDIF 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END FUNCTION DENOM                            
