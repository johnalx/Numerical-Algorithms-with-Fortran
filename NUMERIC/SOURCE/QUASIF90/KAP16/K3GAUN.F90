![          {Gau"s Cubature Formulas for Triangles}*)                   
      SUBROUTINE K3GAUN (USERF, PX, PY, QX, QY, RX, RY, N, METHOD, MOLD,&
      CTRI, WORK, IERR, IUFCLL)                                         
!                                                                       
!*****************************************************************      
!                                                                *      
! Gaussian cubature over triangular regions:                     *      
!                                                                *      
! The FUNCTION USERF(X,Y) is integrated over the triangle PQR    *      
! according to the summed N point gaussian formula using N*N     *      
! sub-triangles.                                                 *      
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
! METHOD  : INTEGER, designating the method: If METHOD=1,2,3 or 7*      
!           the 1, 2, 3 or 7 point Gauss formula is chosen.      *      
! MOLD    : INTEGER, the number in METHOD at the previous call.  *      
!           Upon first call we must have: MOLD different from    *      
!           METHOD                                               *      
!           If K3GAUN is called repeatedly with METHOD=MOLD the  *      
!           internal initializing of parameters is skipped.      *      
! WORK    : 2-dimensional DOUBLE PRECISION array                 *      
!           WORK(3,0:METHOD-1). If METHOD=MOLD this array must   *      
!           contain the initializing parameters for the method.  *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! MOLD    : INTEGER indicating the number of points used in the  *      
!           Gauss method.                                        *      
! CTRI    : DOUBLE PRECISION approximate value for the integral  *      
! WORK    : 2-dimensional DOUBLE PRECISION array                 *      
!           WORK(3,0:METHOD-1) containing the constants needed   *      
!           for the specified method                             *      
! IERR    : error parameter: IERR=0 all is ok                    *      
!                            IERR=1 N is incorrect               *      
!                            IERR=2 the vertices P Q and R are   *      
!                                   collinear                    *      
!                            IERR=3 invalid Number for the method*      
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
! FAC     : DOUBLE PRECISION number, indicates type of triangle  *      
!                              FAC=1.0  (not symmetric) or       *      
!                              FAC=-1.0 (reflection symmetric)   *      
! X       : DOUBLE PRECISION ]   coordinates of the top vertex   *      
! Y       : DOUBLE PRECISION ]   of the sub-triangle in use.     *      
!                                                                *      
! XX      : DOUBLE PRECISION ]   These are auxiliary variables   *      
! YY      : DOUBLE PRECISION ]   that determine the weights.     *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: K3GINI                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author  : Volker KrÅger                                       *      
!  Date    : 06.12.1991                                          *      
!  Source  : FORTRAN 77                                          *      
!                                                                *      
!*****************************************************************      
!                                                                       
! declarations                                                          
!                                                                       
      DOUBLEPRECISION WORK (3, 0:METHOD-1), PX, PY, QX, QY, RX, RY,     &
      CTRI, AREA, EPS, HPQX, HPQY, HPRX, HPRY, FAC, DBLEN, DBLEJ, DBLEI,&
      X, XX, Y, YY, USERF                                               
!                                                                       
! initialize bound for collinearity check                               
!                                                                       
      EPS = 1.0D-06 
!                                                                       
! check N                                                               
!                                                                       
      IF (N.LT.1) THEN 
         IERR = 1 
         RETURN 
      ENDIF 
!                                                                       
!   test for collinearity                                               
!                                                                       
      AREA = PX * QY + QX * RY + RX * PY - PX * RY - QX * PY - RX * QY 
      IF (AREA.LT.EPS) THEN 
         IERR = 2 
         RETURN 
!                                                                       
!   check validity of method number                                     
!                                                                       
      ELSEIF (METHOD.LT.0.OR.METHOD.GT.7.OR.METHOD.GT.3.AND.METHOD.LT.7)&
      THEN                                                              
         IERR = 3 
         RETURN 
      ELSE 
         IERR = 0 
      ENDIF 
!                                                                       
! Initialize if necessary                                               
!                                                                       
      IF (METHOD.NE.MOLD) THEN 
         CALL K3GINI (METHOD, WORK) 
         MOLD = METHOD 
      ENDIF 
!                                                                       
! Initialize IUFCLL                                                     
!                                                                       
      IUFCLL = 0 
!                                                                       
! find twice the area                                                   
!                                                                       
      DBLEN = DBLE (N) 
      AREA = AREA / (DBLEN * DBLEN) 
!                                                                       
! vektorize the step size                                               
!                                                                       
      HPQX = (QX - PX) / DBLEN 
      HPQY = (QY - PY) / DBLEN 
      HPRX = (RX - PX) / DBLEN 
      HPRY = (RY - PY) / DBLEN 
                                                                        
!                                                                       
! Initialize CTRI                                                       
!                                                                       
      CTRI = 0.0D0 
!                                                                       
! Approximate the integral                                              
!                                                                       
      DO 10 JJ = 0, 1 
!                                                                       
!   triangle reflection symmetric or not                                
!                                                                       
         IF (JJ.EQ.0) THEN 
            FAC = 1.0D0 
         ELSE 
            FAC = - 1.0D0 
         ENDIF 
!                                                                       
!   loop along the edge PR                                              
!                                                                       
         DO 20 J = JJ, N - 1 
            DBLEJ = DBLE (J) 
!                                                                       
!   loop along the edge PQ                                              
!                                                                       
            DO 30 I = JJ, N - 1 - J + JJ 
               DBLEI = DBLE (I) 
!                                                                       
!   find the coordinates of the top vertex of the sub-triangle          
!                                                                       
               X = PX + HPQX * DBLEI + HPRX * DBLEJ 
               Y = PY + HPQY * DBLEI + HPRY * DBLEJ 
!                                                                       
!   sum the weighted functional values                                  
!                                                                       
               DO 40 II = 0, METHOD-1 
                  XX = HPQX * WORK (2, II) + HPRX * WORK (3, II) 
                  YY = HPQY * WORK (2, II) + HPRY * WORK (3, II) 
                  CTRI = CTRI + WORK (1, II) * USERF (X + FAC * XX, Y + &
                  FAC * YY)                                             
!                                                                       
!   count number of functional evaluations                              
!                                                                       
                  IUFCLL = IUFCLL + 1 
   40          END DO 
   30       END DO 
   20    END DO 
   10 END DO 
!                                                                       
! Multiply by the area                                                  
!                                                                       
      CTRI = CTRI * AREA 
!                                                                       
! return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE K3GAUN                         
!                                                                       
!                                                                       
      SUBROUTINE K3GINI (METHOD, WORK) 
!                                                                       
!*****************************************************************      
!                                                                *      
! SUBROUTINE that initializes the constants in accordance with   *      
! the method.                                                    *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! METHOD  : INTEGER designating the method: METHOD= 1,2,3 or 7.  *      
!           This indicates the number of points used.            *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! WORK    : 2-dimensional DOUBLE PRECISION array                 *      
!           WORK(3,0:METHOD-1) containing the constants for the  *      
!           method used.                                         *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Volker KrÅger                                      *      
!  Date     : 06.12.1991                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION WORK (3, 0:METHOD-1) 
!                                                                       
!  Initialize the array WORK depending on specified method              
!                                                                       
      IF (METHOD.EQ.1) THEN 
         WORK (1, 0) = 0.5D0 
         WORK (2, 0) = 0.3333333333333333D0 
         WORK (3, 0) = 0.3333333333333333D0 
      ELSEIF (METHOD.EQ.2) THEN 
         WORK (1, 0) = 0.25D0 
         WORK (2, 0) = 0.1666666666666667D0 
         WORK (3, 0) = 0.5D0 
         WORK (1, 1) = 0.25D0 
         WORK (2, 1) = 0.5D0 
         WORK (3, 1) = 0.1666666666666667D0 
      ELSEIF (METHOD.EQ.3) THEN 
         WORK (1, 0) = 0.16666666666666667D0 
         WORK (2, 0) = 0.1666666666666667D0 
         WORK (3, 0) = 0.1666666666666667D0 
         WORK (1, 1) = 0.16666666666666667D0 
         WORK (2, 1) = 0.6666666666666667D0 
         WORK (3, 1) = 0.1666666666666667D0 
         WORK (1, 2) = 0.16666666666666667D0 
         WORK (2, 2) = 0.1666666666666667D0 
         WORK (3, 2) = 0.6666666666666667D0 
      ELSEIF (METHOD.EQ.7) THEN 
         WORK (1, 0) = 0.1125D0 
         WORK (2, 0) = 0.3333333333333333D0 
         WORK (3, 0) = 0.3333333333333333D0 
         WORK (1, 1) = 0.0661970763942531D0 
         WORK (2, 1) = 0.4701420641051151D0 
         WORK (3, 1) = 0.4701420641051151D0 
         WORK (1, 2) = 0.0661970763942531D0 
         WORK (2, 2) = 0.05971587178976981D0 
         WORK (3, 2) = 0.4701420641051151D0 
         WORK (1, 3) = 0.0661970763942531D0 
         WORK (2, 3) = 0.4701420641051151D0 
         WORK (3, 3) = 0.05971587178976981D0 
         WORK (1, 4) = 0.06296959027241357D0 
         WORK (2, 4) = 0.1012865073234563D0 
         WORK (3, 4) = 0.1012865073234563D0 
         WORK (1, 5) = 0.06296959027241357D0 
         WORK (2, 5) = 0.7974269853530873D0 
         WORK (3, 5) = 0.1012865073234563D0 
         WORK (1, 6) = 0.06296959027241357D0 
         WORK (2, 6) = 0.1012865073234563D0 
         WORK (3, 6) = 0.7974269853530873D0 
      ENDIF 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE K3GINI                         
