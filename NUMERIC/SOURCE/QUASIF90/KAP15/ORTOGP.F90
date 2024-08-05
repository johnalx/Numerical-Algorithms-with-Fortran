![           Formulas}*)                                                
      SUBROUTINE ORTOGP (N, AINT, X, W, IERR, WK, WKD, IWK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This subroutine determines the nodes and weights of the       *      
!  generalized GAUSSIAN quadrature formula.                      *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N   :  number of nodes                                        *      
!  AINT:  vector AINT(0:2*N-1) containing the values AINT(I) of  *      
!         the integral of the function (X**I) * G(X) over        *      
!         specified intervals.                                   *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  X   : the vector of nodes for the integration formula.        *      
!  W   : the vector of weights of the integration formula for    *      
!        the above nodes.                                        *      
!  IERR: error parameter:                                        *      
!           if IERR is different from zero, no useable values    *      
!           could be determined, since the linear system of      *      
!           equations for the nodes and weights is too ill-con-  *      
!           ditioned or because the zero finding algorithm could *      
!           not find the zeros reliably.                         *      
!                                                                *      
!  AUXILIARY VARIABLES:                                          *      
!  ====================                                          *      
!  WK    vector WK(1:N**2+4*N+1)                                 *      
!  WKD   DOUBLE PRECISION vector of length N+1                   *      
!  IWK   INTEGER vector of length N                              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: ORTPOL, GAUSS, MULLRP                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Eberhard Heyne                                     *      
!  date     : 05.25.1988                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION AINT (0:2 * N - 1), X (N), W (N), WK ( * ), IWK ( * ) 
      DIMENSION WKD (0: * ) 
!                                                                       
!  WK is a vector of length N**2+4*N+1.                                 
!  Here this vector is split into several subvectors, that              
!  are partially overlapping, as these are used                         
!  successively.                                                        
!  The DOUBLE PRECISION vector WKD is entered externally, since         
!  it was proven that not all compilers accept a type-new-              
!  declaration of preexisting variables                                 
!                                                                       
      CALL ORTPOL (N, AINT, X, W, IERR, WK (1), WK (N + 1), WK (2 * N + &
      1), IWK, WK (4 * N + 1), WK (N + 2), WK (2 * N + 3), WK (3 * N +  &
      4), WKD)                                                          
      RETURN 
      END SUBROUTINE ORTOGP                         
!                                                                       
!                                                                       
      SUBROUTINE ORTPOL (N, AINT, X, W, IERR, Q, RS, DFG, IPFG, GL, WKR,&
      ETA, Z, WKD)                                                      
!                                                                       
!*****************************************************************      
!                                                                *      
!  Auxiliary routine for ORTOGP.                                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : GAUSS, MULLRP                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Eberhard Heyne                                     *      
!  date     : 05.25.1988                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION AINT (0:2 * N - 1), Q (0:N), ETA (0:N), W (N), X (N) 
      DIMENSION GL (0:N - 1, 0:N - 1), RS (0:N) 
      DIMENSION DFG (N), IPFG (N), WKR (0:N), Z (0:1, N), WKD (0:N) 
!                                                                       
!  System matrix GL and the right hand side RS of the linear system     
!  that determines the polynomial coefficients                          
!                                                                       
      DO 10 I = 0, N - 1 
         DO 12 K = 0, N - 1 
            GL (I, K) = - AINT (I + K) 
   12    END DO 
         RS (I) = AINT (I + N) 
   10 END DO 
!                                                                       
!  solve the system of equations, obtain the polynomial coefficients Q  
!                                                                       
      CALL GAUSS (N, GL, N, RS, Q, MARK, DFG, IPFG) 
      IERR = 1 
      IF (MARK.EQ.0) RETURN 
      Q (N) = 1.0D0 
!                                                                       
!  determine zeros of the polynomial Q                                  
!                                                                       
      CALL MULLRP (N, Q, 200, NFND, Z, WKR, WKD) 
      IERR = 2 
      IF (NFND.NE.N) RETURN 
      IERR = 0 
!                                                                       
!  from the theory we know that all zeros must be                       
!  real. They are moved from Z to X.                                    
!                                                                       
      DO 1002 K = 1, N 
 1002 X (K) = Z (0, K) 
!                                                                       
!  loop over all zeros                                                  
!                                                                       
      DO 14 I = 1, N 
!                                                                       
!  ETA/F is the LAGRANGE interpolation polynomial for X(I)              
!                                                                       
         F = 1.0D0 
         ETA (N - 1) = 1.0D0 
         DO 16 K = N - 2, 0, - 1 
            ETA (K) = Q (K + 1) + ETA (K + 1) * X (I) 
            F = ETA (K) + F * X (I) 
   16    END DO 
!                                                                       
!  determine weights W(I) of the quadrature formula                     
!                                                                       
         W (I) = 0.0D0 
         DO 18 K = 0, N - 1 
            W (I) = W (I) + ETA (K) * AINT (K) / F 
   18    END DO 
   14 END DO 
      RETURN 
      END SUBROUTINE ORTPOL                         
