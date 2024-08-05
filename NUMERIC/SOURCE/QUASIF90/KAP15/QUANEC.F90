![KA{P 15}{Numerical Integration}{Numerical Integration}*)              
      SUBROUTINE QUANEC (A, B, N, NRP, FCT, VAL, IEG, ERREST, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This subroutine integrates a given function FCT(X) over the   *      
!  interval [A,B] using the summed NEWTON-COTES formulas.        *      
!  N indicates how many subdivisions are used for the NC-formula.*      
!  Thus the step size becomes H1=(B-A)/(N*VN) and for this we    *      
!  obtain the approximate value QH1 for the integral. In case           
!  N > 1 the process is repeated with step size                         
!  H2=(B-A)/((N-1)*VN). A different H2 may be chosen as well; the*      
!  program can be adjusted accordingly. An error estimate is com-*      
!  puted for QH1 by comparison with the second approximation QH2.*      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  A    : left endpoint of the interval.                         *      
!  B    : right endpoint of the interval.                        *      
!  N    : number N of the sub-intervals in the summed quadrature *      
!         formula.                                               *      
!  NRP  : index for the specific NEWTON-COTES formula to be used:*      
!          = 1: trapezoidal rule                                 *      
!          = 2: SIMPSON rule                                     *      
!          = 3: 3/8 - formula                                    *      
!          = 4: 4/90 - rule                                      *      
!          = 5: 5/288 - rule                                     *      
!          = 6: 6/840 - rule                                     *      
!          = 7: 7/17280 - rule                                   *      
!  FCT  : name of the FUNCTION FCT(X) that we want to integrate. *      
!         It has to be provided by the user in the following     *      
!         format:                                                *      
!              DOUBLE PRECISION FUNCTION FCT (X).                *      
!         The function has to be defined as EXTERNAL in the      *      
!         calling program.                                       *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  VAL    : computed approximate value of the integral.          *      
!  IEG    : error order of the method chosen.                    *      
!  ERREST : error estimate for N > 1.                            *      
!  IERR   : error parameter:                                     *      
!           = 0: everything o.k.                                 *      
!           = 1: false input parameter.                          *      
!           = 2: A > B, VAL is given the appropriate sign        *      
!           = 3: A = B: VAL = 0.0                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : none                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Hermann-Josef Rheinbach                            *      
!  date     : 08.19.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION CFORML (7, 8), FRML (7), Q (2), H (2), IRANGE (7) 
      INTEGER IEG 
!                                                                       
      DATA CFORML / 1.0D0, 1.0D0, 1.0D0, 7.0D0, 19.0D0, 41.0D0, 751.0D0,&
      1.0D0, 4.0D0, 3.0D0, 32.0D0, 75.0D0, 216.0D0, 3577.0D0, 0.0D0,    &
      1.0D0, 3.0D0, 12.0D0, 50.0D0, 27.0D0, 1323.0D0, 0.0D0, 0.0D0,     &
      1.0D0, 32.0D0, 50.0D0, 272.0D0, 2989.0D0, 0.0D0, 0.0D0, 0.0D0,    &
      7.0D0, 75.0D0, 27.0D0, 2989.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,      &
      19.0D0, 216.0D0, 1323.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,     &
      41.0D0, 3577.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,       &
      751.0D0 /                                                         
!                                                                       
      DATA FRML / 2.0D0, 6.0D0, 8.0D0, 90.0D0, 288.0D0, 840.0D0,        &
      17280.0D0 /                                                       
      DATA IRANGE / 2, 4, 4, 6, 6, 8, 8 / 
      DATA IFLAG / 0 / 
!                                                                       
      IF (IFLAG.EQ.0) THEN 
         IFLAG = 1 
         DO 10 K = 1, 8 
            DO 50 I = 1, 7 
               CFORML (I, K) = CFORML (I, K) / FRML (I) * DBLE (I) 
   50       END DO 
   10    END DO 
      ENDIF 
!                                                                       
!*  test the input parameters                                           
!                                                                       
      IERR = 1 
      IF (NRP.LE.0.OR.NRP.GE.8..OR.N.LE.0) RETURN 
!                                                                       
      IF (A.EQ.B) THEN 
         IERR = 3 
         VAL = 0.0D0 
         RETURN 
      ENDIF 
!                                                                       
!*  check the interval endpoints                                        
!                                                                       
      IERR = 0 
      IF (A.GT.B) THEN 
         IERR = 2 
         BL = B 
         BU = A 
         FC = - 1.0D0 
      ELSE 
         BL = A 
         BU = B 
         FC = 1.0D0 
      ENDIF 
!                                                                       
      IEG = IRANGE (NRP) 
      DO 20 I = 1, 2 
         KEND = N + 1 - I 
         H (I) = (BU - BL) / DBLE (KEND * NRP) 
         Q (I) = 0.0D0 
         DO 30 K = 1, KEND 
            XN = BL + DBLE ( (K - 1) * NRP) * H (I) 
            DO 40 L = 1, NRP + 1 
               Q (I) = Q (I) + CFORML (NRP, L) * FCT (XN + DBLE (L - 1) &
               * H (I) )                                                
   40       END DO 
   30    END DO 
         Q (I) = Q (I) * H (I) 
         IF (N.EQ.1) THEN 
            VAL = Q (1) * FC 
            RETURN 
         ENDIF 
   20 END DO 
!                                                                       
      ERREST = DABS ( (Q (1) - Q (2) ) / ( (H (2) / H (1) ) **IEG -     &
      1.0D0) )                                                          
      VAL = Q (1) * FC 
!                                                                       
      RETURN 
      END SUBROUTINE QUANEC                         
