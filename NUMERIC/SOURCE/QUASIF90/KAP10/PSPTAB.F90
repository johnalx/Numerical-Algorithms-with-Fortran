      SUBROUTINE PSPTAB (N, NL, TBEG, TEND, T, AX, BX, CX, DX, AY, BY,  &
      CY, DY, NT, XTAB, YTAB, IERR)                                     
!                                                                       
!*****************************************************************      
!                                                                *      
!  Program to create a table of values for parametric cubic      *      
!  splines with component functions SX(T), SY(T), given in the   *      
!  following form:                                               *      
!                                                                *      
!  SX := SX(T) = AX(I) + BX(I)(T-T(I)) + CX(I)(T-T(I))**2 +      *      
!                                      + DX(I)(T-T(I))**3        *      
!                                                                *      
!  SY := SY(T) = AY(I) + BY(I)(T-T(I)) + CY(I)(T-T(I))**2 +      *      
!                                      + DY(I)(T-T(I))**3        *      
!                                                                *      
!  for T in the interval [T(I),T(I+1)], I=0,1,...,N-1.           *      
!                                                                *      
!                                                                *      
!  This program creates a table of function values XTAB = SX(TW) *      
!  and YTAB = SY(TW) where TW lies in [TBEG,TEND]. We use the    *      
!  following conventions:                                        *      
!   - if TBEG < T(0), the end polynomial P(0) will be evaluated  *      
!     for all values XTAB < T(0)                                 *      
!   - if TEND > T(N), the end polynomial P(N-1) will be evaluated*      
!     for all values XTAB > T(N)                                 *      
!   - in every table the interval end points TBEG and TEND and   *      
!     all nodes T(I) in between will occur in the table          *      
!   - in each subinterval [T(I),T(I+1)] the table is created for *      
!     equidistant steps of size H. Thus H will always depend on  *      
!     the length of the given interval and on the length NL of   *      
!     the table.                                                 *      
!   - the input parameter NL presets an approximate table length;*      
!     the actual table length is NT+1 (NT denotes the final index*      
!     of the table of values). We must have  0 < NT < NL+N+3.    *      
!                                                                *      
!                                                                *      
!  ASSUMPTIONS:   TBEG <  TEND                                   *      
!  ============   NL   >= 0                                      *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N    :  Index of the final knot T(N)                          *      
!  NL   :  Given table length for dimensioning of vectors XTAB   *      
!          and YTAB                                              *      
!  TBEG :  starting value of the table                           *      
!  TEND :  ending value of the table                             *      
!  T    :  N+1-vector T(0:N); the nodes T(I), I=0,1,...,N        *      
!  AX :  ] N+1-vectors ..(0:N);                                  *      
!  BX :  ] the components in positions 0 to N-1 contain the      *      
!  CX :  ] spline coefficients for the component function SX(T)  *      
!  DX :  ]                                                       *      
!                                                                *      
!  AY :  ] N+1-vectors ..(0:N);                                  *      
!  BY :  ] the components in positions 0 to N-1 contain the      *      
!  CY :  ] spline coefficients for the component function SY(T)  *      
!  DY :  ]                                                       *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  NT   :  final table index; equals actual table length - 1     *      
!  XTAB :  vector XTAB(0:NL+N+2) ] the elements in positions     *      
!  YTAB :  vector YTAB(0:NL+N+2) ] 0 to NT contain the table     *      
!  IERR :  Error parameter                                       *      
!          = 0 : Everything o.k.                                 *      
!          = 1 : Stop caused by TBEG >= TEND                     *      
!          = 2 : Stop caused by NL < 0                           *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Guenter Palm                                       *      
!  date     : 03.28.1989                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!-----declarations------------------------------------------------      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION XTAB (0:NL + N + 2), YTAB (0:NL + N + 2), T (0:N),&
      AX (0:N), BX (0:N), CX (0:N), DX (0:N), AY (0:N), BY (0:N),       &
      CY (0:N), DY (0:N)                                                
!                                                                       
!-----checking the input------------------------------------------      
!                                                                       
      IF (TEND.LE.TBEG) THEN 
         IERR = 1 
         RETURN 
      ELSEIF (NL.LT.0) THEN 
         IERR = 2 
         RETURN 
      ENDIF 
      IERR = 0 
!                                                                       
!-----determine the interval [T(I),T(I+1)]-------------------------     
!     that contains TBEG; index it by IBEG                              
!                                                                       
      I = 0 
      K = N 
   10 M = (I + K) / 2 
      IF (TBEG.LT.T (M) ) THEN 
         K = M 
      ELSE 
         I = M 
      ENDIF 
      IF (K.GT.I + 1) GOTO 10 
      IBEG = I 
!                                                                       
!-----determine the interval [T(I),T(I+1)]------------------------      
!     that contains TEND; index it by IEND                              
!                                                                       
      K = N 
   20 M = (I + K) / 2 
      IF (TEND.LT.T (M) ) THEN 
         K = M 
      ELSE 
         I = M 
      ENDIF 
      IF (K.GT.I + 1) GOTO 20 
      IEND = I 
!                                                                       
!-----determine the table values XTAB(I), YTAB(I), I=0,1,...,NT---      
!                                                                       
!     Initialize                                                        
!                                                                       
      HP = TEND-TBEG 
      FC = NL / HP 
      NT = 0 
      TW = TBEG 
!                                                                       
      IF (IBEG.NE.IEND) THEN 
!                                                                       
         IF (TBEG.LT.T (0) ) THEN 
            IP = 0 
         ELSE 
            IP = 1 
         ENDIF 
!                                                                       
         IF (TEND.GT.T (N) ) THEN 
            IM = 0 
         ELSE 
            IM = 1 
         ENDIF 
!                                                                       
!     determine the values for the table at TBEG to T(IBEG+IP)          
!                                                                       
         I = IBEG 
         TD = TW - T (I) 
         XTAB (NT) = ( (DX (I) * TD+CX (I) ) * TD+BX (I) ) * TD+AX (I) 
         YTAB (NT) = ( (DY (I) * TD+CY (I) ) * TD+BY (I) ) * TD+AY (I) 
         DIF = T (I + IP) - TBEG 
         TIV = DIF * FC 
         ITV = INT (TIV) 
         IF ( (TIV - ITV) .GT.0.0D0) ITV = ITV + 1 
         IF (ITV.GT.0) H = DIF / ITV 
         DO 30 J = 1, ITV - 1, 1 
            NT = NT + 1 
            TW = TW + H 
            TD = TW - T (I) 
            XTAB (NT) = ( (DX (I) * TD+CX (I) ) * TD+BX (I) ) * TD+AX ( &
            I)                                                          
            YTAB (NT) = ( (DY (I) * TD+CY (I) ) * TD+BY (I) ) * TD+AY ( &
            I)                                                          
   30    END DO 
         NT = NT + 1 
         IF ( (IEND-IBEG) .NE.1) THEN 
!                                                                       
!           determine the table values at T(IBEG+IP) to                 
!           T(IEND-IM+1)                                                
!                                                                       
            IBP = IBEG + IP 
            IEM = IEND-IM 
            DO 40 I = IBP, IEM, 1 
               TW = T (I) 
               XTAB (NT) = AX (I) 
               YTAB (NT) = AY (I) 
               DIF = T (I + 1) - T (I) 
               TIV = DIF * FC 
               ITV = INT (TIV) 
               IF ( (TIV - ITV) .GT.0.0D0) ITV = ITV + 1 
               IF (ITV.GT.0) H = DIF / ITV 
               DO 50 J = 1, ITV - 1, 1 
                  NT = NT + 1 
                  TW = TW + H 
                  TD = TW - T (I) 
                  XTAB (NT) = ( (DX (I) * TD+CX (I) ) * TD+BX (I) )     &
                  * TD+AX (I)                                           
                  YTAB (NT) = ( (DY (I) * TD+CY (I) ) * TD+BY (I) )     &
                  * TD+AY (I)                                           
   50          END DO 
               NT = NT + 1 
   40       END DO 
         ENDIF 
         TW = T (IEND-IM + 1) 
      ENDIF 
!                                                                       
!     determine the table values from the location                      
!     tabulated last to TEND                                            
!                                                                       
      I = IEND 
      TD = TW - T (I) 
      XTAB (NT) = ( (DX (I) * TD+CX (I) ) * TD+BX (I) ) * TD+AX (I) 
      YTAB (NT) = ( (DY (I) * TD+CY (I) ) * TD+BY (I) ) * TD+AY (I) 
      DIF = TEND-TW 
      TIV = DIF * FC 
      ITV = INT (TIV) 
      IF ( (TIV - ITV) .GT.0.0D0) ITV = ITV + 1 
      IF (ITV.GT.0) H = DIF / ITV 
      DO 60 J = 1, ITV - 1, 1 
         NT = NT + 1 
         TW = TW + H 
         TD = TW - T (I) 
         XTAB (NT) = ( (DX (I) * TD+CX (I) ) * TD+BX (I) ) * TD+AX (I) 
         YTAB (NT) = ( (DY (I) * TD+CY (I) ) * TD+BY (I) ) * TD+AY (I) 
   60 END DO 
      NT = NT + 1 
      TD = TEND-T (I) 
      XTAB (NT) = ( (DX (I) * TD+CX (I) ) * TD+BX (I) ) * TD+AX (I) 
      YTAB (NT) = ( (DY (I) * TD+CY (I) ) * TD+BY (I) ) * TD+AY (I) 
      RETURN 
      END SUBROUTINE PSPTAB                         
