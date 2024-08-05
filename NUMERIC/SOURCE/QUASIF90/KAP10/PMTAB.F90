      SUBROUTINE PMTAB (N, NTAB, TBEG, TEND, DELT, T, AX, BX, CX, DX,   &
      EX, FX, AY, BY, CY, DY, EY, FY, XTAB, YTAB, LENTAB, IERR)         
!                                                                       
!*****************************************************************      
!                                                                *      
!  Tabulates a parametric Hermite spline. The nodes of the       *      
!  spline are retained, as long as they lie in the table interval*      
!  [TBEG, TEND]. This allows the program to create input data    *      
!  for graphics subroutines.                                     *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N        : Number of nodes for the spline                     *      
!  NTAB     : Maximal length of the table. NTAB should be at     *      
!             least   (TEND-TBEG)/DELT+N.                        *      
!  TBEG     : ) Parameter interval, where the spline is tabulated*      
!  TEND     : ) Necessary condition:                             *      
!                     T(1) <= TBEG <= TEND <= T(N)               *      
!  DELT     : Step size. The table of values is created for      *      
!             T = TBEG, TBEG + DELT, ..., TEND                   *      
!  T        : N-vector T(1:N); the nodes of the splineparameter T*      
!  AX, BX   : ) N-vectors ..(1:N); the coefficients of the spline*      
!  CX, DX   : ) component SX(T)                                  *      
!  EX, FX   : )                                                  *      
!  AY, BY   : ) N-vectors ..(1:N); the coefficients of the Spline*      
!  CY, DY   : ) component SY(T)                                  *      
!  EY, FY   : )                                                  *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  =================                                             *      
!  XTAB     : ) NTAB-vectors ..(1:NTAB); the table of values     *      
!  YTAB     : )                                                  *      
!  LENTAB   : length of the table                                *      
!  IERR     : = 0, no error                                      *      
!             = 1, TBEG > TEND .OR. TBEG < T(1) .OR. TEND > T(N) *      
!             = 2, DELT <= 0.                                    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Guido Dubois                                    *      
!  Date        : 4.18.1993                                      *       
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER N, NTAB, LENTAB, IERR, I, J, K, M, IBEG, IEND, LBEG,      &
      IFLAG, MACHPD, IBP1, IEM1                                         
      DOUBLEPRECISION T (1:N), AX (1:N), BX (1:N), CX (1:N), DX (1:N),  &
      EX (1:N), FX (1:N), AY (1:N), BY (1:N), CY (1:N), DY (1:N),       &
      EY (1:N), FY (1:N), XTAB (1:NTAB), YTAB (1:NTAB), TBEG, TEND,     &
      DELT, FMACHP, EPS, T0, T1                                         
!                                                                       
!  Local storage of the error EPS in case that this subroutine          
!  is called repeatedly                                                 
!                                                                       
      SAVE EPS, IFLAG 
      DATA IFLAG / 0 / 
      IERR = 0 
!                                                                       
!  Check input parameters                                               
!                                                                       
      IF (TBEG.GT.TEND.OR.TBEG.LT.T (1) .OR.TEND.GT.T (N) ) THEN 
         IERR = 1 
         RETURN 
      ENDIF 
      IF (DELT.LE.0.0D0) THEN 
         IERR = 2 
         RETURN 
      ENDIF 
!                                                                       
!  Compute the machine constant                                         
!                                                                       
      IF (IFLAG.EQ.0) THEN 
         IFLAG = 1 
         FMACHP = 1.0D0 
    5    FMACHP = 0.5D0 * FMACHP 
         IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 5 
         FMACHP = 2.0D0 * FMACHP 
         EPS = 1000.0D0 * FMACHP 
      ENDIF 
!                                                                       
!  Determine the starting and terminal intervals of computation         
!                                                                       
      LENTAB = 0 
      I = 1 
      K = N 
   10 M = (I + K) / 2 
      IF (M.NE.I) THEN 
         IF (TBEG.GE.T (M) ) THEN 
            I = M 
         ELSE 
            K = M 
         ENDIF 
         GOTO 10 
      ENDIF 
      IBEG = I 
      K = N 
   20 M = (I + K) / 2 
      IF (M.NE.I) THEN 
         IF (TEND.GT.T (M) ) THEN 
            I = M 
         ELSE 
            K = M 
         ENDIF 
         GOTO 20 
      ENDIF 
      IEND = I 
!                                                                       
      T0 = TBEG 
      T1 = T0 - T (IBEG) 
      IF (IBEG.NE.IEND) THEN 
!                                                                       
!  First interval                                                       
!                                                                       
         LENTAB = INT ( (T (IBEG + 1) - TBEG + EPS) / DELT) + 1 
         DO 30 J = 1, LENTAB 
            XTAB (J) = ( ( ( (FX (IBEG) * T1 + EX (IBEG) ) * T1 + DX (  &
            IBEG) ) * T1 + CX (IBEG) ) * T1 + BX (IBEG) ) * T1 + AX (   &
            IBEG)                                                       
            YTAB (J) = ( ( ( (FY (IBEG) * T1 + EY (IBEG) ) * T1 + DY (  &
            IBEG) ) * T1 + CY (IBEG) ) * T1 + BY (IBEG) ) * T1 + AY (   &
            IBEG)                                                       
            T0 = T0 + DELT 
            T1 = T1 + DELT 
   30    END DO 
!                                                                       
!  Second to (N-1)st interval                                           
!                                                                       
         IF ( (IEND-IBEG) .NE.1) THEN 
            IBP1 = IBEG + 1 
            IEM1 = IEND-1 
            DO 40 I = IBP1, IEM1 
               IF (DABS (T0 - DELT - T (I) ) .GT.EPS) THEN 
                  LENTAB = LENTAB + 1 
                  XTAB (LENTAB) = AX (I) 
                  YTAB (LENTAB) = AY (I) 
               ENDIF 
               LBEG = LENTAB + 1 
               LENTAB = LENTAB + INT ( (T (I + 1) - T0 + EPS) / DELT)   &
               + 1                                                      
               T1 = T0 - T (I) 
               DO 50 J = LBEG, LENTAB 
                  XTAB (J) = ( ( ( (FX (I) * T1 + EX (I) ) * T1 + DX (I)&
                  ) * T1 + CX (I) ) * T1 + BX (I) ) * T1 + AX (I)       
                  YTAB (J) = ( ( ( (FY (I) * T1 + EY (I) ) * T1 + DY (I)&
                  ) * T1 + CY (I) ) * T1 + BY (I) ) * T1 + AY (I)       
                  T0 = T0 + DELT 
                  T1 = T1 + DELT 
   50          END DO 
   40       END DO 
         ENDIF 
      ELSE 
         LENTAB = LENTAB + 1 
         XTAB (LENTAB) = ( ( ( (FX (IBEG) * T1 + EX (IBEG) ) * T1 + DX (&
         IBEG) ) * T1 + CX (IBEG) ) * T1 + BX (IBEG) ) * T1 + AX (IBEG) 
         YTAB (LENTAB) = ( ( ( (FY (IBEG) * T1 + EY (IBEG) ) * T1 + DY (&
         IBEG) ) * T1 + CY (IBEG) ) * T1 + BY (IBEG) ) * T1 + AY (IBEG) 
         T0 = T0 + DELT 
         T1 = T1 + DELT 
      ENDIF 
!X(IBEG))*T1+BX(IBEG))*T1+AX(IBEG)                                      
      IF (DABS (T0 - DELT - T (IEND) ) .GT.EPS.AND.T (IEND) .GT.TBEG)   &
      THEN                                                              
         LENTAB = LENTAB + 1 
         XTAB (LENTAB) = AX (IEND) 
         YTAB (LENTAB) = AY (IEND) 
      ENDIF 
      LBEG = LENTAB + 1 
      LENTAB = LENTAB + INT ( (TEND-T0 + EPS) / DELT) + 1 
      T1 = T0 - T (IEND) 
      IF (LENTAB.GE.LBEG) THEN 
         DO 60 J = LBEG, LENTAB 
            XTAB (J) = ( ( ( (FX (I) * T1 + EX (I) ) * T1 + DX (I) )    &
            * T1 + CX (I) ) * T1 + BX (I) ) * T1 + AX (I)               
            YTAB (J) = ( ( ( (FY (I) * T1 + EY (I) ) * T1 + DY (I) )    &
            * T1 + CY (I) ) * T1 + BY (I) ) * T1 + AY (I)               
            T0 = T0 + DELT 
            T1 = T1 + DELT 
   60    END DO 
      ENDIF 
      IF (DABS (T0 - DELT - TEND) .GT.EPS) THEN 
         LENTAB = LENTAB + 1 
         T0 = TEND 
         T1 = T0 - T (IEND) 
         XTAB (LENTAB) = ( ( ( (FX (I) * T1 + EX (I) ) * T1 + DX (I) )  &
         * T1 + CX (I) ) * T1 + BX (I) ) * T1 + AX (I)                  
         YTAB (LENTAB) = ( ( ( (FY (I) * T1 + EY (I) ) * T1 + DY (I) )  &
         * T1 + CY (I) ) * T1 + BY (I) ) * T1 + AY (I)                  
      ENDIF 
      RETURN 
      END SUBROUTINE PMTAB                          
