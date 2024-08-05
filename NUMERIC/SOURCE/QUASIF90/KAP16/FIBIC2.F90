      SUBROUTINE FIBIC2 (N, M, A, X, Y, XU, YU, XO, YO, VALUE, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  The subroutine determines a double integral of a spline       *      
!  function over a rectangle [XU, XO] x [YU, YO].                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: XYINTV, FIBIC1                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Eberhard Heyne                                     *      
!  date     : 02.15.1983                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      PARAMETER (KDIM = 3, LDIM = 3) 
!                                                                       
      DIMENSION A (0:N, 0:M, 0:KDIM, 0:LDIM), X (0:N), Y (0:M) 
!                                                                       
!*  determine intervals and relative coordinates                        
!                                                                       
      CALL XYINTV (N, M, X, Y, IU, JU, XIU, ETAU, XU, YU, JERR) 
      IERR = JERR 
      IF (JERR.NE.0) RETURN 
      CALL XYINTV (N, M, X, Y, IO, JO, XIO, ETAO, XO, YO, JERR) 
      IERR = JERR + 1 
      IF (JERR.NE.0) RETURN 
      FACTOR = 1.0D0 
      IF (IU.GT.IO) THEN 
         I = IO 
         IO = IU 
         IU = I 
         XI = XIO 
         XIO = XIU 
         XIU = XI 
         FACTOR = - FACTOR 
      ENDIF 
      IF (JU.GT.JO) THEN 
         J = JO 
         JO = JU 
         JU = J 
         ETA = ETAO 
         ETAO = ETAU 
         ETAU = ETA 
         FACTOR = - FACTOR 
      ENDIF 
      S = FIBIC1 (N, M, A, IU, JU, XIU, ETAU) 
      S = S - FIBIC1 (N, M, A, IU, JO, XIU, ETAO) 
      S = S - FIBIC1 (N, M, A, IO, JU, XIO, ETAU) 
      S = S + FIBIC1 (N, M, A, IO, JO, XIO, ETAO) 
      DO 101 I = IU, IO - 1 
         S = S + FIBIC1 (N, M, A, I, JO, X (I + 1) - X (I), ETAO) 
         S = S - FIBIC1 (N, M, A, I, JU, X (I + 1) - X (I), ETAU) 
  101 END DO 
      DO 103 J = JU, JO - 1 
         S = S + FIBIC1 (N, M, A, IO, J, XIO, Y (J + 1) - Y (J) ) 
         S = S - FIBIC1 (N, M, A, IU, J, XIU, Y (J + 1) - Y (J) ) 
         DO 102 I = IU, IO - 1 
            S = S + FIBIC1 (N, M, A, I, J, X (I + 1) - X (I), Y (J + 1) &
            - Y (J) )                                                   
  102    END DO 
  103 END DO 
      IERR = 0 
      VALUE = FACTOR * S 
      RETURN 
      END SUBROUTINE FIBIC2                         
