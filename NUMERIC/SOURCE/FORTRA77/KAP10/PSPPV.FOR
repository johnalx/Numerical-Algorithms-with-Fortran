C[BA*)
C[LE*)
      SUBROUTINE PSPPV (N, XN, FN, T, MT, IERR)
C[IX{PSPPV}*)
C
C*****************************************************************
C                                                                *
C  PSPPV computes the parameter values T(I), I=0,1,...,N, for    *
C  parametric splines. By using the parameter MT we can specify  *
C  whether PSPPV determines the T(I) from the chordal length or  *
C  from the arc length of the curve.                             *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N  :  Index of final node                                     *
C  XN :  vector XN(0:N); the nodes XN(I), I = 0,1,..,N           *
C  FN :  vector FN(0:N); the function values at the nodes FN(I) =*
C        = FN( XN(I)).                                           *
C  MT :  Indicates the method of determining the parameter:      *
C        MT =  1: The parameter values are determined from the   *
C                 chordal length                                 *
C        MT <> 1: The parameter values are determined using the  *
C                 arc length                                     *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  T    :  vector T(0:N); the parameter values T(I)              *
C  IERR :  Error parameter                                       *
C          = 0: Everything o.k.                                  *
C          = 1: The parameter values T(I) are not monotonic,     *
C               T(I) >= T(I+1) for some I between 0 and N-1      *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines required: none                                    *
C                                                                *
C                                                                *
C================================================================*
C                                                                *
C  author   : Guenter Palm                                       *
C  date     : 10.11.1989                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C-----declarations------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION XN(0:N), FN(0:N), T(0:N)
C
C-----initializing------------------------------------------------
C
      IERR = 1
C
C-----determine the parameter values ... -------------------------
C
      IF (MT .EQ. 1) THEN
C
C        ... via the chordal length
C
        T(0) = 0.0D0
        DO 10 I=1,N,1
          DELTX = XN(I) - XN(I-1)
          DELTY = FN(I) - FN(I-1)
          DELTA = DELTX*DELTX + DELTY*DELTY
          IF (DELTA .LE. 0.0D0) RETURN
          T(I) = T(I-1) + DSQRT(DELTA)
   10   CONTINUE
      ELSE
C
C        ... or using the arc length
C
        T(0) = 0.0D0
        DO 20 I=0,N-2
          A = XN(I+1) - XN(I)
          B = FN(I+1) - FN(I)
          C = XN(I+2) - XN(I+1)
          D = FN(I+2) - FN(I+1)
          E = XN(I+2) - XN(I)
          F = FN(I+2) - FN(I)
          DN = A*D - B*C
          IF (DN .EQ. 0.0D0) THEN
            G = 1.0D0
          ELSE
            DZ = C*E + D*F
            IF (DZ .EQ. 0.0D0) THEN
              G = 1.57D0
            ELSE
              DZ = DZ/DN
              G = DSQRT(1.0D0+DZ*DZ) * DATAN(1.0D0/DABS(DZ))
            ENDIF
          ENDIF
          DT = G * DSQRT(A*A + B*B)
          IF (DT .LE. 0.0D0) RETURN
          T(I+1) = T(I) + DT
   20   CONTINUE
        G =  A
        A = -C
        C = -G
        G =  B
        B = -D
        D = -G
        E = -E
        F = -F
        DN = A*D - B*C
        IF (DN .EQ. 0.0D0) THEN
          G = 1.0D0
        ELSE
          DZ = C*E + D*F
          IF (DZ .EQ. 0.0D0) THEN
            G = 1.57D0
          ELSE
            DZ = DZ/DN
            G = DSQRT(1.0D0+DZ*DZ) * DATAN(1.0D0/DABS(DZ))
          ENDIF
        ENDIF
        DT = G * DSQRT(A*A + B*B)
        IF (DT .LE. 0.0D0) RETURN
        T(N) = T(N-1) + DT
      ENDIF
C
      IERR = 0
      RETURN
      END
