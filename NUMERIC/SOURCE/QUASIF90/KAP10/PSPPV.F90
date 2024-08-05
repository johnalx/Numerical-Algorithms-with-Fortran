      SUBROUTINE PSPPV (N, XN, FN, T, MT, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  PSPPV computes the parameter values T(I), I=0,1,...,N, for    *      
!  parametric splines. By using the parameter MT we can specify  *      
!  whether PSPPV determines the T(I) from the chordal length or  *      
!  from the arc length of the curve.                             *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N  :  Index of final node                                     *      
!  XN :  vector XN(0:N); the nodes XN(I), I = 0,1,..,N           *      
!  FN :  vector FN(0:N); the function values at the nodes FN(I) =*      
!        = FN( XN(I)).                                           *      
!  MT :  Indicates the method of determining the parameter:      *      
!        MT =  1: The parameter values are determined from the   *      
!                 chordal length                                 *      
!        MT <> 1: The parameter values are determined using the  *      
!                 arc length                                     *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  T    :  vector T(0:N); the parameter values T(I)              *      
!  IERR :  Error parameter                                       *      
!          = 0: Everything o.k.                                  *      
!          = 1: The parameter values T(I) are not monotonic,     *      
!               T(I) >= T(I+1) for some I between 0 and N-1      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required: none                                    *      
!                                                                *      
!                                                                *      
!================================================================*      
!                                                                *      
!  author   : Guenter Palm                                       *      
!  date     : 10.11.1989                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!-----declarations------------------------------------------------      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION XN (0:N), FN (0:N), T (0:N) 
!                                                                       
!-----initializing------------------------------------------------      
!                                                                       
      IERR = 1 
!                                                                       
!-----determine the parameter values ... -------------------------      
!                                                                       
      IF (MT.EQ.1) THEN 
!                                                                       
!        ... via the chordal length                                     
!                                                                       
         T (0) = 0.0D0 
         DO 10 I = 1, N, 1 
            DELTX = XN (I) - XN (I - 1) 
            DELTY = FN (I) - FN (I - 1) 
            DELTA = DELTX * DELTX + DELTY * DELTY 
            IF (DELTA.LE.0.0D0) RETURN 
            T (I) = T (I - 1) + DSQRT (DELTA) 
   10    END DO 
      ELSE 
!                                                                       
!        ... or using the arc length                                    
!                                                                       
         T (0) = 0.0D0 
         DO 20 I = 0, N - 2 
            A = XN (I + 1) - XN (I) 
            B = FN (I + 1) - FN (I) 
            C = XN (I + 2) - XN (I + 1) 
            D = FN (I + 2) - FN (I + 1) 
            E = XN (I + 2) - XN (I) 
            F = FN (I + 2) - FN (I) 
            DN = A * D-B * C 
            IF (DN.EQ.0.0D0) THEN 
               G = 1.0D0 
            ELSE 
               DZ = C * E+D * F 
               IF (DZ.EQ.0.0D0) THEN 
                  G = 1.57D0 
               ELSE 
                  DZ = DZ / DN 
                  G = DSQRT (1.0D0 + DZ * DZ) * DATAN (1.0D0 / DABS (DZ)&
                  )                                                     
               ENDIF 
            ENDIF 
            DT = G * DSQRT (A * A + B * B) 
            IF (DT.LE.0.0D0) RETURN 
            T (I + 1) = T (I) + DT 
   20    END DO 
         G = A 
         A = - C 
         C = - G 
         G = B 
         B = - D 
         D = - G 
         E = - E 
         F = - F 
         DN = A * D-B * C 
         IF (DN.EQ.0.0D0) THEN 
            G = 1.0D0 
         ELSE 
            DZ = C * E+D * F 
            IF (DZ.EQ.0.0D0) THEN 
               G = 1.57D0 
            ELSE 
               DZ = DZ / DN 
               G = DSQRT (1.0D0 + DZ * DZ) * DATAN (1.0D0 / DABS (DZ) ) 
            ENDIF 
         ENDIF 
         DT = G * DSQRT (A * A + B * B) 
         IF (DT.LE.0.0D0) RETURN 
         T (N) = T (N - 1) + DT 
      ENDIF 
!                                                                       
      IERR = 0 
      RETURN 
      END SUBROUTINE PSPPV                          
