C[KA{F 0}{Auxiliary Library}{Auxiliary Library}*)
      INTEGER FUNCTION MACHPD(X)
C[IX{MACHPD}*)
      DOUBLE PRECISION X
      MACHPD=0
      IF (1.0D0 .LT. X) MACHPD=1
      RETURN
      END
