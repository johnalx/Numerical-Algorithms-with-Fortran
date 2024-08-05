      SUBROUTINE PSPANA (TV, N, T, AX, BX, CX, DX, AY, BY, CY, DY, SX,  &
      SY, XDVT, YDVT, C1, C2, IERR)                                     
!                                                                       
!*****************************************************************      
!                                                                *      
!  A program to evaluate parametric cubic splines with component *      
!  functions SX(T), SY(T) given in the following form:           *      
!                                                                *      
!  SX := SX(T) = AX(I) + BX(I)(T-T(I)) + CX(I)(T-T(I))**2 +      *      
!                                      + DX(I)(T-T(I))**3        *      
!                                                                *      
!  SY := SY(T) = AY(I) + BY(I)(T-T(I)) + CY(I)(T-T(I))**2 +      *      
!                                      + DY(I)(T-T(I))**3        *      
!                                                                *      
!  for T in the interval [T(I),T(I+1)], I=0,1,...,N-1.           *      
!                                                                *      
!  This program determines the function value and that of the    *      
!  1st, 2nd and 3rd derivative of the component functions SX(T)  *      
!  and SY(T) at T=TV, as well as the 1st and 2nd derivative of   *      
!  the curve K that is described in the plane by of the component*      
!  functions SX and SY.                                          *      
!                                                                *      
!                                                                *      
!  NOTE:  This evaluation routine is not well suited for making  *      
!  =====  a table of values of SX(T) and SY(T).                  *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  TV :  Location where the spline functions SX(T) and SY(T) are *      
!        to be evaluated.                                        *      
!  N  :  Index of the final node T(N)                            *      
!  T  :  vector T(0:N) ; the nodes T(I), I=0,1,...,N.            *      
!  AX :  ] N+1-vectors ..(0:N);                                  *      
!  BX :  ] the components in positions 0 to N-1 contain the      *      
!  CX :  ] spline coefficients for the component function SX(T)  *      
!  DX :  ]                                                       *      
!  AY :  ] N+1-vectors ..(0:N);                                  *      
!  BY :  ] the components in positions 0 to N-1 contain the      *      
!  CY :  ] spline coefficients for the component function SY(T)  *      
!  DY :  ]                                                       *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  SX :  Function value of SX(T) at T=TV                         *      
!  SY :  Function value of SY(T) at T=TV                         *      
!  XDVT : ]  vectors ..(1:3);                                    *      
!  YDVT : ]  the element in positin I contains the Ith derivative*      
!            of the function SX(T) or SY(T), respectively, at    *      
!            T=TV for I = 1, 2, 3.                               *      
!  C1 :  1st derivative of the curve K at T=TV.                  *      
!        It is obtained from the following equation:             *      
!        C1 = YDVT(1)/XDVT(1)                                    *      
!  C2 :  2nd derivative of the curve K at T=TV.                  *      
!        It is obtained from the following equation:             *      
!        C2 = (XDVT(1)*YDVT(2) - XDVT(2)*YDVT(1))/(XDVT(1)**3)   *      
!                                                                *      
!  IERR :  Error parameter                                       *      
!          = 0 : Everything o.k.                                 *      
!          = 1 : XDVT(1) = 0; values for C1 and C2 were not      *      
!                determined                                      *      
!          = 2 : The magnitude of the 1st derivative of SX(TV),  *      
!                ABS(XDVT(1)), is positive but less than four    *      
!                times the machine constant. The values of C1 and*      
!                C2 can therefore not be determined accurately.  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: SPLFVD, MACHPD                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Guenter Palm                                       *      
!  date     : 04.15.1988                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!-----declarations------------------------------------------------      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION T (0:N), AX (0:N), BX (0:N), CX (0:N), DX (0:N),  &
      AY (0:N), BY (0:N), CY (0:N), DY (0:N), XDVT (3), YDVT (3)        
      LOGICAL FLAG 
      SAVE FMACHP, FLAG 
!                                                                       
!-----initializing------------------------------------------------      
!                                                                       
      DATA FLAG / .TRUE. / 
      IERR = 0 
!                                                                       
!-----determine the machine constant (only on 1st call)-----------      
!                                                                       
      IF (FLAG) THEN 
         FMACHP = 1.0D0 
   10    FMACHP = 0.5D0 * FMACHP 
         IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
         FMACHP = 8.0D0 * FMACHP 
         FLAG = .FALSE. 
      ENDIF 
!                                                                       
!-----determine the functional values and the derivatives----           
!     of the component functions SX(T) and SY(T)                        
!                                                                       
      CALL SPLFVD (TV, N, T, AX, BX, CX, DX, SX, XDVT (1), XDVT (2),    &
      XDVT (3) )                                                        
      CALL SPLFVD (TV, N, T, AY, BY, CY, DY, SY, YDVT (1), YDVT (2),    &
      YDVT (3) )                                                        
!                                                                       
!-----determine the 1st and 2nd derivatives of the curve K--------      
!                                                                       
      IF (XDVT (1) .EQ.0.0D0) THEN 
         IERR = 1 
      ELSE 
         IF (DABS (XDVT (1) ) .LE.FMACHP) IERR = 2 
         C1 = YDVT (1) / XDVT (1) 
         C2 = (XDVT (1) * YDVT (2) - XDVT (2) * YDVT (1) ) / (XDVT (1)  &
         **3)                                                           
      ENDIF 
      RETURN 
      END SUBROUTINE PSPANA                         
