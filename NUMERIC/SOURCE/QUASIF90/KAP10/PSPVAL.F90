      SUBROUTINE PSPVAL (TV, N, T, AX, BX, CX, DX, AY, BY, CY, DY, SX,  &
      SY)                                                               
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE evaluates a parametric cubic spline given by  *      
!  its component functions SX(T), SY(T) in the form:             *      
!                                                                *      
!  SX := SX(T) = AX(I) + BX(I)(T-T(I)) + CX(I)(T-T(I))**2 +      *      
!                                      + DX(I)(T-T(I))**3        *      
!                                                                *      
!  SY := SY(T) = AY(I) + BY(I)(T-T(I)) + CY(I)(T-T(I))**2 +      *      
!                                      + DY(I)(T-T(I))**3        *      
!                                                                *      
!  for T in the interval [T(I),T(I+1)], I=0,1,...,N-1.           *      
!                                                                *      
!  The program computes the values of SX(T) and SY(T) at T=TV.   *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  TV :  value where the component functions SX(T) and SY(T)     *      
!        are to be evaluated                                     *      
!  N  :  index of the last node                                  *      
!  T  :  vector T(0:N); the nodes T(I), I=0,1,...,N              *      
!  AX :  ] N+1-vectors ..(0:N);                                  *      
!  BX :  ] their elements 0, ..., N-1 describe the component     *      
!  CX :  ] function SX(T).                                       *      
!  DX :  ]                                                       *      
!                                                                *      
!  AY :  ] N+1-vectors ..(0:N);                                  *      
!  BY :  ] their elements 0, ..., N-1 describe the component     *      
!  CY :  ] function SY(T).                                       *      
!  DY :  ]                                                       *      
!          The elements BX(N), CX(N), DX(N), BY(N), CY(N) and    *      
!          DY(N) are used for auxiliary purposes.                *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  SX :  Functional value for SX(T) at T=TV                      *      
!  SY :  Functional value for SY(T) at T=TV                      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required: SPVAL                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : GÅnter Palm                                        *      
!  Date     : 06.01.1991                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!-----Declarations                                                      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION T (0:N), AX (0:N), BX (0:N), CX (0:N), DX (0:N),  &
      AY (0:N), BY (0:N), CY (0:N), DY (0:N)                            
!                                                                       
!-----Compute the functional values of the component splines            
!     SX(T) and SY(T)                                                   
!                                                                       
      SX = SPVAL (TV, N, T, AX, BX, CX, DX) 
      SY = SPVAL (TV, N, T, AY, BY, CY, DY) 
!                                                                       
      RETURN 
      END SUBROUTINE PSPVAL                         
