![KA{P 16}{Numerical Cubature}{Numerical Cubature}*)                    
![  {Newton--Cotes Cubature Formulas for Rectangles}                    
![  {Newton--Cotes Cubature Formulas for Rectangles}*)                  
      SUBROUTINE K4NECN (USERF, A, B, IP, C, D, IQ, METHOD, MOLD, CREC, &
      ESTDIV, DIVIAT, WORK, IERR, IUFCLL)                               
!                                                                       
!*****************************************************************      
!                                                                *      
! Cubature for rectangular regions using the NEWTON-COTES        *      
! formulas:                                                      *      
!                                                                *      
! The FUNCTION USERF(X,Y) is integrated over the rectangle       *      
! [A,B] x [C,D] using the summed NEWTON-COTES formulas.          *      
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
! METHOD  : INTEGER designating the method used:                 *      
!                     1: trapezoidal rule                        *      
!                     2: SIMPSON rule                            *      
!                     3: 3/8-rule                                *      
!                     4: 4/90-rule                               *      
!                     5: 5/288-rule                              *      
!                     6: 6/840-rule                              *      
!                     7: 7/17280-rule                            *      
! MOLD    : INTEGER, the number of the method used at a previous *      
!           call of this subroutine. On first call we must have  *      
!           that MOLD does not equal METHOD !                    *      
!           In a subsequent call of K4NECN with METHOD=MOLD      *      
!           the internal initializing of parameters is skipped.  *      
! ESTDIV  : LOGICAL variable; If ESTDIV=TRUE we compute an error *      
!           estimate, if ESTDIV=FALSE we do not.                 *      
! WORK    : DOUBLE PRECISION vector WORK(0:METHOD+2)             *      
!           If METHOD=MOLD, WORK must contain the constants      *      
!           initialized for the proper method.                   *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! MOLD    : INTEGER, the number of the method used               *      
! CREC    : DOUBLE PRECISION value for the integral              *      
! DIVIAT  : DOUBLE PRECISION error estimate                      *      
!           If ESTDIV=TRUE, we perform an additional cubature    *      
!           with halved stepsize for error estimation.           *      
! WORK    : DOUBLE PRECISION vector WORK(0:METHOD+2),            *      
!           contains the constants for the method used.          *      
! IERR    : error parameter: IERR=0 all is ok                    *      
!                            IERR=1 number of intervals in       *      
!                                   X-direction erroneous        *      
!                            IERR=2 number of intervals in       *      
!                                   Y-direction erroneous        *      
!                            IERR=3 method number erroneous      *      
!                            IERR=4 integrating over an interval *      
!                                   of length zero.              *      
! IUFCLL  : INTEGER, the number of functional evaluations        *      
!           performed.                                           *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! =================                                              *      
! I,J,K   : INTEGER loop variables                               *      
! KMAX    : INTEGER, the number of cubature passes desired       *      
! IPX     : INTEGER, number of intervals in X-direction          *      
! IPY     : INTEGER, number of intervals in Y-direction          *      
! DBLEI   : DOUBLE PRECISION value for I                         *      
! HX      : DOUBLE PRECISION step size in X-direction            *      
! HY      : DOUBLE PRECISION step size in Y-direction            *      
! CRECH   : DOUBLE PRECISION auxiliary variable                  *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required subroutines: K4INIT, GRIDOT                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author  : Volker KrÅger                                       *      
!  Date    : 06.12.1991                                          *      
!  Source  : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION WORK (0:METHOD+2), A, B, C, D, CREC, CRECH,       &
      DIVIAT, HX, HY, DBLEI, GRIDOT, USERF                              
!                                                                       
! Initialize LOGICAL variable ESTDIV                                    
!                                                                       
      LOGICAL ESTDIV 
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
      ELSEIF (METHOD.LT.1.OR.METHOD.GT.7) THEN 
         IERR = 3 
         RETURN 
      ELSEIF (A.EQ.B.OR.C.EQ.D) THEN 
         IERR = 4 
         RETURN 
      ELSE 
         IERR = 0 
      ENDIF 
!                                                                       
! Initialize as needed                                                  
!                                                                       
      IF (METHOD.NE.MOLD) THEN 
         CALL K4INIT (METHOD, WORK) 
         MOLD = METHOD 
      ENDIF 
!                                                                       
! Set maximal number of cubature passes                                 
!                                                                       
      IF (ESTDIV) THEN 
         KMAX = 2 
      ELSE 
         KMAX = 1 
      ENDIF 
!                                                                       
! Determine number of X and Y sub-intervals                             
!                                                                       
      DO 10 K = 1, KMAX 
         IPX = K * IP * METHOD 
         IQY = K * IQ * METHOD 
!                                                                       
! Initialize CREC                                                       
!                                                                       
         CREC = 0.0D0 
!                                                                       
! Find step sizes in both X- and Y-directions                           
!                                                                       
         HX = (B - A) / DBLE (IPX) 
         HY = (D-C) / DBLE (IQY) 
!                                                                       
! Compute an approximate value of the integral                          
!                                                                       
         DO 20 I = 0, IPX 
            DBLEI = DBLE (I) 
            DO 30 J = 0, IQY 
               CREC = CREC + GRIDOT (I, J, WORK, METHOD, IPX, IQY)      &
               * USERF (A + DBLEI * HX, C + DBLE (J) * HY)              
               IUFCLL = IUFCLL + 1 
   30       END DO 
   20    END DO 
         CREC = CREC * HX * HY * WORK (METHOD+1) 
!                                                                       
! In case error estimate is desired, store first value for the integral 
!                                                                       
         IF (ESTDIV.AND.K.EQ.1) CRECH = CREC 
   10 END DO 
!                                                                       
! Estimate the error                                                    
!                                                                       
      IF (ESTDIV) DIVIAT = (CREC - CRECH) / (2.0D0**WORK (METHOD+2)     &
      - 1.0D0)                                                          
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE K4NECN                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION GRIDOT (I, J, WORK, METHOD, IPX, IQY) 
!                                                                       
!*****************************************************************      
!                                                                *      
! FUNCTION that determines the weights at the nodes.             *      
!                                                                *      
! In summed Newton-Cotes cubature the functional values are      *      
! given different weights that depend on their position on the   *      
! boundary, center or at the join of two rectangles.             *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! I       : INTEGER, the number of the node in X-direction       *      
! J       : INTEGER, the number of the node in Y-direction       *      
! WORK    : DOUBLE PRECISION vector WORK(0:METHOD), containing   *      
!           the constants for the method                         *      
! METHOD  : INTEGER designating the method chosen:               *      
!                     1: trapezoidal rule                        *      
!                     2: SIMPSON rule                            *      
!                     3: 3/8-rule                                *      
!                     4: 4/90-rule                               *      
!                     5: 5/288-rule                              *      
!                     6: 6/840-rule                              *      
!                     7: 7/17280-rule                            *      
! IPX     : INTEGER number of intervals used in X-direction      *      
! IPY     : INTEGER number of intervals used in Y-direction      *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLE:                                                *      
! ================                                               *      
! K       : INTEGER auxiliary variable                           *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Volker KrÅger                                      *      
!  Date     : 06.121991                                         *       
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION WORK (0:METHOD) 
!                                                                       
! Determine the weights for the nodes                                   
!                                                                       
      K = MOD (I, METHOD) 
      GRIDOT = WORK (K) 
      IF (K.EQ.0.AND.I.GT.0.AND.I.LT.IPX) GRIDOT = 2.0D0 * GRIDOT 
      K = MOD (J, METHOD) 
      GRIDOT = GRIDOT * WORK (K) 
      IF (K.EQ.0.AND.J.GT.0.AND.J.LT.IQY) GRIDOT = 2.0D0 * GRIDOT 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END FUNCTION GRIDOT                           
!                                                                       
!                                                                       
      SUBROUTINE K4INIT (METHOD, WORK) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Subroutine that initializes the constants for the various      *      
! methods.                                                       *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! METHOD  : INTEGER, the number designating the method:          *      
!                     1: trapezoidal rule                        *      
!                     2: SIMPSON rule                            *      
!                     3: 3/8 rule                                *      
!                     4: 4/90 rule                               *      
!                     5: 5/288 rule                              *      
!                     6: 6/840 rule                              *      
!                     7: 7/17280 rule                            *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETER:                                              *      
! =================                                              *      
! WORK    : DOUBLE PRECISION vector WORK(0:METHOD+2),            *      
!           containing the constants for the method              *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! I       : INTEGER loop counter                                 *      
! J       : INTEGER auxiliary variable                           *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Volker KrÅger                                      *      
!  Date     : 06.12.1991                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION WORK (0:METHOD+2) 
!                                                                       
! Initialize upper half of WORK depending on method                     
!                                                                       
      IF (METHOD.EQ.1) THEN 
         WORK (0) = 1.0D0 
         WORK (METHOD+1) = 2.0D0 
         WORK (METHOD+2) = 2.0D0 
         J = 1 
      ELSEIF (METHOD.EQ.2) THEN 
         WORK (0) = 1.0D0 
         WORK (1) = 4.0D0 
         WORK (METHOD+1) = 6.0D0 
         WORK (METHOD+2) = 4.0D0 
         J = 2 
      ELSEIF (METHOD.EQ.3) THEN 
         WORK (0) = 1.0D0 
         WORK (1) = 3.0D0 
         WORK (METHOD+1) = 8.0D0 
         WORK (METHOD+2) = 4.0D0 
         J = 2 
      ELSEIF (METHOD.EQ.4) THEN 
         WORK (0) = 7.0D0 
         WORK (1) = 32.0D0 
         WORK (2) = 12.0D0 
         WORK (METHOD+1) = 90.0D0 
         WORK (METHOD+2) = 6.0D0 
         J = 3 
      ELSEIF (METHOD.EQ.5) THEN 
         WORK (0) = 19.0D0 
         WORK (1) = 75.0D0 
         WORK (2) = 50.0D0 
         WORK (METHOD+1) = 288.0D0 
         WORK (METHOD+2) = 6.0D0 
         J = 3 
      ELSEIF (METHOD.EQ.6) THEN 
         WORK (0) = 41.0D0 
         WORK (1) = 216.0D0 
         WORK (2) = 27.0D0 
         WORK (3) = 272.0D0 
         WORK (METHOD+1) = 840.0D0 
         WORK (METHOD+2) = 8.0D0 
         J = 4 
      ELSEIF (METHOD.EQ.7) THEN 
         WORK (0) = 751.0D0 
         WORK (1) = 3577.0D0 
         WORK (2) = 1323.0D0 
         WORK (3) = 2989.0D0 
         WORK (METHOD+1) = 17280.0D0 
         WORK (METHOD+2) = 8.0D0 
         J = 4 
      ENDIF 
!                                                                       
! Determine lower half symmetrically                                    
!                                                                       
      DO 10 I = J, METHOD 
         WORK (I) = WORK (METHOD-I) 
   10 END DO 
!                                                                       
! Determine the multiplication factors for the summed values            
! in the cubature formula                                               
!                                                                       
      WORK (METHOD+1) = DBLE (METHOD) / WORK (METHOD+1) 
      WORK (METHOD+1) = WORK (METHOD+1) * WORK (METHOD+1) 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE K4INIT                         
