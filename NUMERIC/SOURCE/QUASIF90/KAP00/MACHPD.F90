      INTEGER FUNCTION MACHPD (X) 
      DOUBLEPRECISION X 
      MACHPD = 0 
      IF (1.0D0.LT.X) MACHPD = 1 
      RETURN 
      END FUNCTION MACHPD                           
