      SUBROUTINE K3RORI (USERF, PX, PY, QX, QY, RX, RY, N, WORK, CTRI,  &
      DIVIAT, IERR, IUFCLL)                                             
!                                                                       
!*****************************************************************      
!                                                                *      
! Cubature over triangular regions using the summed 3-point      *      
! formula and ROMBERG-RICHARDSON extrapolation:                  *      
!                                                                *      
! Using the summed  3-point cubature formula of NEWTON-COTES     *      
! we evaluate the integral of the Funktion USERF(X,Y) over the   *      
! triangle P Q R by subdividing it into similar triangles.       *      
! The number of cubature steps with halved sub-triangle edges is *      
! designated by N.                                               *      
! A RICHARDSON-extrapolation is used for an improved             *      
! approximation of the integral.                                 *      
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
! WORK    : 2-dimensional DOUBLE PRECISION array                 *      
!           WORK(3,0:METHOD-1).                                  *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! CTRI    : DOUBLE PRECISION approximate value for the integral  *      
! DIVIAT  : DOUBLE PRECISION error estimate                      *      
! IERR    : error parameter: IERR=0 all is ok                    *      
!                            IERR=1 N is incorrect               *      
!                            IERR=2 the vertices P Q and R are   *      
!                                   collinear                    *      
! IUFCLL  : INTEGER, the number of function evaluations performed*      
!                                                                *      
!                                                                *      
! INTERMEDIATE VARIABLES:                                        *      
! =======================                                        *      
! I       : loop variable                                        *      
! IUFHLP  : auxiliary varialbe counting function evaluations     *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: K3GNEC3, RORIEX                         *      
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
      DOUBLEPRECISION WORK (0:N - 1, 2), PX, PY, QX, QY, RX, RY, CTRI,  &
      DIVIAT                                                            
      EXTERNAL USERF 
!                                                                       
! Initialize IUFCLL                                                     
!                                                                       
      IUFCLL = 0 
!                                                                       
! Check validity of N                                                   
!                                                                       
      IF (N.LT.2) THEN 
         IERR = 1 
         RETURN 
      ENDIF 
!                                                                       
! Perform N cubatures                                                   
!                                                                       
      DO 10 I = 0, N - 1 
         CALL K3NEC3 (USERF, PX, PY, QX, QY, RX, RY, 2**I, WORK (I, 1), &
         IERR, IUFHLP)                                                  
         IUFCLL = IUFCLL + IUFHLP 
         IF (IERR.NE.0) RETURN 
   10 END DO 
!                                                                       
! Find an approximate integral value and an error estimate              
! by using RICHARDSON-extrapolation                                     
!                                                                       
      CALL RORIEX (WORK (0, 1), WORK (0, 2), N, 2, CTRI, DIVIAT) 
!                                                                       
! Return to calling program                                             
!                                                                       
      RETURN 
      END SUBROUTINE K3RORI                         
