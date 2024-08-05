      SUBROUTINE HMTAB (N, NTAB, XBEG, XEND, DELTX, X, A, B, C, D, E, F,&
      XTAB, YTAB, LENTAB, IERR)                                         
!                                                                       
!*****************************************************************      
!                                                                *      
!  Constructs a value table for a hermitian polynomial spline of *      
!  degree five over an arbitrary interval inside the interval of *      
!  definition (X(1),X(N)). The nodes for the spline inside the   *      
!  subinterval are also tabulated.                               *      
!  This allows the program to be used to create input data for   *      
!  graphics subroutines.                                         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N        : Number of nodes for the spline                     *      
!  NTAB     : Maximal length of the table. NTAB should be at     *      
!             least  (XEND-XBEG)/DELTX+N                         *      
!  XBEG     : ) Interval, where the spline is to be tabulated    *      
!  XEND     : ) with the necessary inclusion condition:          *      
!                   X(1) <= XBEG <= XEND <= X(N)                 *      
!  DELTX    : Step size. The values are created for x-ordinates  *      
!             X = XBEG, XBEG + DELTX, ..., XEND                  *      
!  X        : N-vector X(1:N); the nodes for the spline          *      
!  A, B, C  : ) N-vectors ..(1:N); the spline coefficients       *      
!  D, E, F  : )                                                  *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  =================                                             *      
!  XTAB     : ) NTAB-vectors ..(1:NTAB); the value table         *      
!  YTAB     : ) Specifically, YTAB(I) = S(XTAB(I)) for           *      
!             )          I = 1, ..., LENTAB                      *      
!  LENTAB   : size of the table                                  *      
!  IERR     : = 0, no error                                      *      
!             = 1, XBEG > XEND .OR. XBEG < X(1) .OR. XEND > X(N) *      
!             = 2, DELTX <= 0.                                   *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Guido Dubois                                    *      
!  Date        : 1.30.1993                                       *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER N, NTAB, LENTAB, IERR, I, J, K, M, IBEG, IEND, LBEG,      &
      IFLAG, MACHPD, IBP1, IEM1                                         
      DOUBLEPRECISION X (1:N), A (1:N), B (1:N), C (1:N), D (1:N),      &
      E (1:N), F (1:N), XTAB (1:NTAB), YTAB (1:NTAB), XBEG, XEND, DELTX,&
      X0, X1, FMACHP, EPS                                               
!                                                                       
!  Local storage of the error EPS in case this subroutine is            
!  called repeatedly.                                                   
!                                                                       
      SAVE EPS, IFLAG 
      DATA IFLAG / 0 / 
      IERR = 0 
!                                                                       
!  Check input parameters                                               
!                                                                       
      IF (XBEG.GT.XEND.OR.XBEG.LT.X (1) .OR.XEND.GT.X (N) ) THEN 
         IERR = 1 
         RETURN 
      ENDIF 
      IF (DELTX.LE.0.0D0) THEN 
         IERR = 2 
         RETURN 
      ENDIF 
!                                                                       
!  Find the machine constant                                            
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
!  Determine the initial and terminal intervals for the computations    
!                                                                       
      LENTAB = 0 
      I = 1 
      K = N 
   10 M = (I + K) / 2 
      IF (M.NE.I) THEN 
         IF (XBEG.GE.X (M) ) THEN 
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
         IF (XEND.GT.X (M) ) THEN 
            I = M 
         ELSE 
            K = M 
         ENDIF 
         GOTO 20 
      ENDIF 
      IEND = I 
!                                                                       
      X0 = XBEG 
      X1 = X0 - X (IBEG) 
      IF (IBEG.NE.IEND) THEN 
!                                                                       
!  First interval                                                       
!                                                                       
         LENTAB = INT ( (X (IBEG + 1) - XBEG + EPS) / DELTX) + 1 
         DO 30 J = 1, LENTAB 
            XTAB (J) = X0 
            YTAB (J) = ( ( ( (F (IBEG) * X1 + E (IBEG) ) * X1 + D (IBEG)&
            ) * X1 + C (IBEG) ) * X1 + B (IBEG) ) * X1 + A (IBEG)       
            X0 = X0 + DELTX 
            X1 = X1 + DELTX 
   30    END DO 
!                                                                       
!  Second to (N-1)st interval                                           
!                                                                       
         IF ( (IEND-IBEG) .NE.1) THEN 
            IBP1 = IBEG + 1 
            IEM1 = IEND-1 
            DO 40 I = IBP1, IEM1 
               IF (DABS (X0 - DELTX - X (I) ) .GT.EPS) THEN 
                  LENTAB = LENTAB + 1 
                  XTAB (LENTAB) = X (I) 
                  YTAB (LENTAB) = A (I) 
               ENDIF 
               LBEG = LENTAB + 1 
               LENTAB = LENTAB + INT ( (X (I + 1) - X0 + EPS) / DELTX)  &
               + 1                                                      
               X1 = X0 - X (I) 
               DO 50 J = LBEG, LENTAB 
                  XTAB (J) = X0 
                  YTAB (J) = ( ( ( (F (I) * X1 + E (I) ) * X1 + D (I) ) &
                  * X1 + C (I) ) * X1 + B (I) ) * X1 + A (I)            
                  X0 = X0 + DELTX 
                  X1 = X1 + DELTX 
   50          END DO 
   40       END DO 
         ENDIF 
      ELSE 
         LENTAB = LENTAB + 1 
         XTAB (LENTAB) = X0 
         YTAB (LENTAB) = ( ( ( (F (IBEG) * X1 + E (IBEG) ) * X1 + D (   &
         IBEG) ) * X1 + C (IBEG) ) * X1 + B (IBEG) ) * X1 + A (IBEG)    
         X0 = X0 + DELTX 
         X1 = X1 + DELTX 
      ENDIF 
!                                                                       
!  Nth interval                                                         
!                                                                       
      IF (DABS (X0 - DELTX - X (IEND) ) .GT.EPS.AND.X (IEND) .GT.XBEG)  &
      THEN                                                              
         LENTAB = LENTAB + 1 
         XTAB (LENTAB) = X (IEND) 
         YTAB (LENTAB) = A (IEND) 
      ENDIF 
      LBEG = LENTAB + 1 
      LENTAB = LENTAB + INT ( (XEND-X0 + EPS) / DELTX) + 1 
      X1 = X0 - X (IEND) 
      IF (LENTAB.GE.LBEG) THEN 
         DO 60 J = LBEG, LENTAB 
            XTAB (J) = X0 
            YTAB (J) = ( ( ( (F (IEND) * X1 + E (IEND) ) * X1 + D (IEND)&
            ) * X1 + C (IEND) ) * X1 + B (IEND) ) * X1 + A (IEND)       
            X0 = X0 + DELTX 
            X1 = X1 + DELTX 
   60    END DO 
      ENDIF 
      IF (DABS (X0 - DELTX - XEND) .GT.EPS) THEN 
         LENTAB = LENTAB + 1 
         X0 = XEND 
         X1 = X0 - X (IEND) 
         XTAB (LENTAB) = X0 
         YTAB (LENTAB) = ( ( ( (F (IEND) * X1 + E (IEND) ) * X1 + D (   &
         IEND) ) * X1 + C (IEND) ) * X1 + B (IEND) ) * X1 + A (IEND)    
      ENDIF 
      RETURN 
      END SUBROUTINE HMTAB                          
