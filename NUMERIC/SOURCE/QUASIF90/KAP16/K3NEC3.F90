![  {Newton--Cotes Cubature Formulas for Triangles}                     
![  {Newton--Cotes Cubature Formulas for Triangles}*)                   
      SUBROUTINE K3NEC3 (USERF, PX, PY, QX, QY, RX, RY, N, CTRI, IERR,  &
      IUFCLL)                                                           
!                                                                       
!*****************************************************************      
!                                                                *      
! Cubature for triangular region using the three point NEWTON-   *      
! COTES formulas:                                                *      
!                                                                *      
! The FUNCTION USERF(X,Y) is integrated over the triangle PQR    *      
! according to the summed NEWTON-COTES formulas using sub-tri-   *      
! angles.                                                        *      
! The dimensions of these sub-triangles are one Nth of those of  *      
! the original triangle PQR.                                     *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! USERF   : user defined FUNCTION USERF(X,Y), whose integral is  *      
!           to be computed.                                      *      
!           The FUNCTION USERF must be declared as EXTERNAL in   *      
!           the calling program.                                 *      
!           The FUNCTION should have the following form:         *      
!                  DOUBLE PRECISION FUNCTION USERF(X,Y)          *      
!                  DOUBLE PRECISION X,Y                          *      
!                         .                                      *      
!                         .                                      *      
!                         .                                      *      
!                  USERF=F(X,Y)                                  *      
!                         .                                      *      
!                         .                                      *      
!                         .                                      *      
!                  RETURN                                        *      
!                  STOP                                          *      
!                                                                *      
! PX      : DOUBLE PRECISION X-coordinate of the vertex P        *      
! PY      : DOUBLE PRECISION Y-coordinate of the vertex P        *      
! QX      : DOUBLE PRECISION X-coordinate of the vertex Q        *      
! QY      : DOUBLE PRECISION Y-coordinate of the vertex Q        *      
! RX      : DOUBLE PRECISION X-coordinate of the vertex R        *      
! RY      : DOUBLE PRECISION Y-coordinate of the vertex R        *      
! N       : INTEGER, counting the number of sub-triangles formed *      
!           along one edge of the triangle.                      *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! CTRI    : DOUBLE PRECISION approximate value for the integral  *      
! IERR    : error parameter: IERR=0 all is ok                    *      
!                            IERR=1 N is incorrect               *      
!                            IERR=2 the vertices P Q and R are   *      
!                                   collinear                    *      
! IUFCLL  : INTEGER, the number of function evaluations performed*      
!                                                                *      
!                                                                *      
! INTERMEDIATE VARIABLES:                                        *      
! =======================                                        *      
! I,J     : loop variables                                       *      
! II,JJ   : loop variables                                       *      
! DBLEN   : DOUBLE PRECISION version of N                        *      
! DBLEI   : DOUBLE PRECISION version of I                        *      
! DBLEJ   : DOUBLE PRECISION version of J                        *      
! AREA    : DOUBLE PRECISION to check collinearity               *      
! EPS     : DOUBLE PRECISION bound for collinearity check        *      
! HPQX    : DOUBLE PRECISION ]   vectoriel representation of     *      
! HPQY    : DOUBLE PRECISION ]   the steps taken along the edge  *      
! HPRX    : DOUBLE PRECISION ]   PQ or PR, respectively          *      
! HPRY    : DOUBLE PRECISION ]                                   *      
! FAC     : DOUBLE PRECISION weight for the node                 *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author  : Volker KrÅger                                       *      
!  Date    : 06.12.1991                                          *      
!  Source  : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION PX, PY, QX, QY, RX, RY, CTRI, AREA, EPS, HPQX,    &
      HPQY, HPRX, HPRY, FAC, DBLEM, DBLEJ, DBLEI, USERF                 
!                                                                       
! Initialize bound for collinearity test                                
!                                                                       
      EPS = 1.0D-06 
!                                                                       
! Initialize IUFCLL                                                     
!                                                                       
      IUFCLL = 0 
!                                                                       
! Check validity of N                                                   
!                                                                       
      IF (N.LT.1) THEN 
         IERR = 1 
         RETURN 
      ENDIF 
!                                                                       
! Test for collinearity                                                 
!                                                                       
      AREA = PX * QY + QX * RY + RX * PY - PX * RY - QX * PY - RX * QY 
      IF (AREA.LT.EPS) THEN 
         IERR = 2 
         RETURN 
      ELSE 
         IERR = 0 
      ENDIF 
!                                                                       
! Number of halved triangular edges                                     
!                                                                       
      M = 2 * N 
      DBLEM = DBLE (M) 
!                                                                       
! Vectoriel representation of the step sizes                            
!                                                                       
      HPQX = (QX - PX) / DBLEM 
      HPQY = (QY - PY) / DBLEM 
      HPRX = (RX - PX) / DBLEM 
      HPRY = (RY - PY) / DBLEM 
                                                                        
!                                                                       
! Initialize CTRI                                                       
!                                                                       
      CTRI = 0.0D0 
!                                                                       
! Compute approximate value for integral                                
!                                                                       
      DO 10 J = 0, M - 1 
         DBLEJ = DBLE (J) 
         DO 20 I = 0, M - J 
            DBLEI = DBLE (I) 
!                                                                       
! Determine weights for the nodes                                       
!                                                                       
            IF (MOD (I, 2) .NE.0.OR.MOD (J, 2) .NE.0) THEN 
               IF (I.EQ.0.OR.J.EQ.0.OR.I.EQ.M - J) THEN 
                  FAC = 1.0D0 
               ELSE 
                  FAC = 2.0D0 
               ENDIF 
               CTRI = CTRI + FAC * USERF (PX + HPQX * DBLEI + HPRX *    &
               DBLEJ, PY + HPQY * DBLEI + HPRY * DBLEJ)                 
               IUFCLL = IUFCLL + 1 
            ENDIF 
   20    END DO 
   10 END DO 
      CTRI = CTRI * AREA / (6.0D0 * DBLE (N) **2.0D0) 
!                                                                       
! Return to the calling program                                         
!                                                                       
      RETURN 
      END SUBROUTINE K3NEC3                         
