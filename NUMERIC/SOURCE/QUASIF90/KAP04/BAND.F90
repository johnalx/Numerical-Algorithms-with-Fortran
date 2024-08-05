![          {Linear Systems with Band Matrices using Pivots}*)          
      SUBROUTINE BAND (AP, LDAP, MB, N, ML, MU, B, IERR, IP) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     BAND solves a linear system of equations for a banded      *      
!     system matrix AP and a right hand side B. This program     *      
!     should not be used if several systems with the same system *      
!     matrix AP are to be solved. In this case it is better to   *      
!     call BANDP below once and BANDS separately for each right  *      
!     hand side.                                                 *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     AP     : 2-dimensional array AP(1:LDAP,1:MB); the system   *      
!              matrix in condensed form. If A denotes the        *      
!              original system matrix, then                      *      
!                 A(I,K)=AP(I,ML+1+K-I)                          *      
!              holds for all indices I, K inside the diagonal    *      
!              band of A. Thus the co-diagonals of A become      *      
!              columns of AP; its rows are condensed to the left.*      
!     LDAP   : leading dimension of AP as defined in the calling *      
!              program                                           *      
!     MB     : number of columns of AP. It is at least           *      
!                 M + MIN(ML,MU) + 1                             *      
!              where  M = ML+MU+1  denotes the band width of AP  *      
!     N      : size of the matrix AP                             *      
!     ML     : number of lower co-diagonals                      *      
!     MU     : number of upper co-diagonals                      *      
!     B      : N-vector  B(1:N); the right hand side             *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     B      : the solution (if IERR .NE. 0)                     *      
!     IERR   : error parameter;                                  *      
!                 IERR = 0 : no error                            *      
!                 IERR.NE.0: matrix AP is numerically singular   *      
!     AP     : as above; tranformed band matrix and transforma-  *      
!              tion coefficients. May be used as input for BANDS *      
!              in order to find solutions for additional right   *      
!              hand sides.                                       *      
!     IP     : N-vector IP(1:N); this vector contains pivoting   *      
!              information that is needed by BANDS.              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: BANDP, BANDS, MACHPD                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Elmar Pohl                                         *      
!  date     : 01.16.1992                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      DIMENSION AP (LDAP, MB), B (N), IP (N) 
!                                                                       
      CALL BANDP (AP, LDAP, MB, N, ML, MU, IP, ISIG, IERR) 
      IF (IERR.NE.0) RETURN 
      CALL BANDS (AP, LDAP, MB, N, ML, MU, B, IP) 
      RETURN 
      END SUBROUTINE BAND                           
!                                                                       
!                                                                       
      SUBROUTINE BANDP (AP, LDAP, MB, N, ML, MU, IP, ISIG, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     BANDP transforms a banded matrix AP to upper or lower      *      
!     triangular form by applying a Gauss algorithm with         *      
!     partial pivot search. The triangular matrix and the        *      
!     transformation matrix are overwritten onto the original    *      
!     array. If IERR=0, the subroutine BANDS may be applied in   *      
!     order to solve several linear systems of equations.        *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     AP     : 2-dimensional array AP(1:LDAP,1:MB); the band     *      
!              matrix in condensed form. The rows of this array  *      
!              contain the rows of the matrix AP, shifted to the *      
!              left. The diagonals of AP, however, appear as     *      
!              columns of AP:                                    *      
!                    A(I,K) = AP(I,ML+1+K-I)                     *      
!              where A denotes the matrix in uncondensed form.   *      
!     LDAP   : leading dimension of AP as defined in the calling *      
!              program                                           *      
!     MB     : column dimension of AP. May differ from that used *      
!              in the calling program; however, it must be at    *      
!              least M + MIN(ML,MU) + 1, where M = MU+ML+1 is    *      
!              the band width of AP.                             *      
!     N      : order of the matrix AP                            *      
!     ML     : number of lower co-diagonals                      *      
!     MU     : number of upper co-diagonals                      *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     AP     : as above; triangular matrix and transforming      *      
!              matrix                                            *      
!     IP     : N-vector IP(1:N); pivoting information which is   *      
!              needed for BANDS                                  *      
!     ISIG   : = (-1)**IZ, where IZ denotes the number of row    *      
!              permutations; to be used for calculating the      *      
!              determinant:                                      *      
!                                    N                           *      
!                 DET(A) = ISIG * PRODUCT AP(I,ML+1)             *      
!                                   I=1                          *      
!                                                                *      
!     IERR   : error parameter;                                  *      
!                 IERR = 0 : no error                            *      
!                 IERR.GT.0: in the IERR-th step, the pivot      *      
!                            element becomes zero, the matrix AP *      
!                            is numerically singular             *      
!                 IERR.LT.0: row IERR contains only zeros,       *      
!                            matrix is singular                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Elmar Pohl                                         *      
!  date     : 01.16.1992                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      DIMENSION AP (LDAP, MB), IP (N) 
      LOGICAL UT 
      M1 = ML + 1 
      M = MU + M1 
      MMB = M + MIN0 (ML, MU) 
!                                                                       
!     calculate the machine constant                                    
!                                                                       
      FMACHP = 1.0D0 
   80 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 80 
      FMACHP = FMACHP * 2.0D0 
!                                                                       
!     determinaing the relative error bound                             
      EPS = 4.0D0 * FMACHP 
!                                                                       
!     Due to the pivot search the band width may increase. When         
!     transforming to upper triangular form, ML additional              
!     upper co-diagonals will appear, while transforming                
!     to lower triangular form creates MU additional lower              
!     co-diagonals. The program automatically selects the               
!     most advantageous triangularization.                              
!                                                                       
      UT = ML.LT.MU 
!                                                                       
!     UT = .TRUE.   : transforming to upper triangular form             
!     UT = .FALSE.  : transforming to lower triangular form             
!                                                                       
!     the additional co-diagonals are stored in the (M+1)th             
!     to MMB-th column of AP. At first they are initialized             
!     to zero.                                                          
!     Simultaneously we compute and store the summed absolute           
!     values of each row in the MMB + 1st column. (To be used           
!     for checking nonsingularity)                                      
!                                                                       
      DO 10 I = 1, N 
         S = 0.0D0 
         DO 90 K = MAX0 (1, I - ML), MIN0 (N, I + MU) 
            S = S + DABS (AP (I, M1 + K - I) ) 
   90    END DO 
         IF (S.EQ.0.0D0) THEN 
            IERR = - I 
            RETURN 
         ENDIF 
         AP (I, MMB + 1) = S 
         DO 20 K = M + 1, MMB 
            AP (I, K) = 0.0D0 
   20    END DO 
   10 END DO 
!                                                                       
      ISIG = 1 
      IERR = 0 
!                                                                       
!     the index bounds and the direction of the sweep in the            
!     following loops are set depending on UT in such a manner          
!     that an upper or a lower triangular form is obtained.             
!                                                                       
      IF (UT) THEN 
         IBEG = 1 
         IEND = N - 1 
         ISTEP = 1 
         KBEG = 1 
         JBEG = 1 
      ELSE 
         IBEG = N 
         IEND = 2 
         ISTEP = - 1 
         KBEG = - 1 
         JBEG = - 1 
      ENDIF 
!                                                                       
!     loop over all rows                                                
!                                                                       
      DO 30 I = IBEG, IEND, ISTEP 
         IF (UT) THEN 
            KEND = MIN0 (ML, N - I) 
         ELSE 
            KEND = MAX0 (1 - I, - MU) 
         ENDIF 
!                                                                       
!        pivot search: searching for the element of largest             
!        magnitude in the I-th column; store its index                  
!        (relative to I) in IV.                                         
!                                                                       
         IV = 0 
         AM = DABS (AP (I, M1) ) 
         DO 40 K = KBEG, KEND, ISTEP 
            IF (DABS (AP (K + I, M1 - K) ) .GT.AM) THEN 
               AM = DABS (AP (K + I, M1 - K) ) 
               IV = K 
            ENDIF 
   40    END DO 
         IF (AM / AP (I + IV, MMB + 1) .LE.EPS) THEN 
!                                                                       
!           element of largest magnitude of the                         
!           column is zero, thus the matrix is singular                 
!                                                                       
            IERR = I 
            RETURN 
         ENDIF 
         IP (I) = IV 
         IF (UT) THEN 
            KJEND = MIN0 (IV + MU, N - I) 
         ELSE 
            KJEND = MAX0 (1 - I, IV - ML) 
         ENDIF 
         IF (IV.NE.0) THEN 
!                                                                       
!           interchange row I with row I+IV                             
!                                                                       
            ISIG = - ISIG 
            DO 50 K = 0, KJEND, ISTEP 
               KM = K + M1 
               IF (KM.LE.0) KM = KM + MMB 
               H = AP (I, KM) 
               AP (I, KM) = AP (I + IV, K + M1 - IV) 
               AP (I + IV, K + M1 - IV) = H 
   50       END DO 
            H = AP (I, MMB + 1) 
            AP (I, MMB + 1) = AP (I + IV, MMB + 1) 
            AP (I + IV, MMB + 1) = H 
         ENDIF 
!                                                                       
!        elimination;                                                   
!        loop including all rows below row I:                           
!                                                                       
         DO 60 K = KBEG, KEND, ISTEP 
            AP (K + I, M1 - K) = AP (K + I, M1 - K) / AP (I, M1) 
!                                                                       
!           loop including all columns to the right of column I         
!                                                                       
            DO 70 J = JBEG, KJEND, ISTEP 
               JK = J + M1 - K 
               JM = J + M1 
!                                                                       
!              With our condensing scheme, additional upper             
!              co-diagonals would have to be stored to the right        
!              of the M-th column and additional lower ones to the      
!              left of the first column.                                
!              The latter would mean that we would need access to       
!              storage space in front of AP(1,1).                       
!              In order to prevent this, MMB is added to all            
!              column subscripts < 1. Thus additional lower             
!              co-diagonals are also stored to the right of AP.         
!                                                                       
               IF (JK.LE.0) JK = JK + MMB 
               IF (JM.LE.0) JM = JM + MMB 
               AP (K + I, JK) = AP (K + I, JK) - AP (K + I, M1 - K)     &
               * AP (I, JM)                                             
   70       END DO 
   60    END DO 
   30 END DO 
      IP (IEND+ISTEP) = 0 
      RETURN 
      END SUBROUTINE BANDP                          
!                                                                       
!                                                                       
      SUBROUTINE BANDS (AP, LDAP, MB, N, ML, MU, B, IP) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     After calling BANDP, the subroutine BANDS solves the linear*      
!     system of equations consisting of a banded, condensed      *      
!     coefficient matrix AP and a right hand side B.             *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     AP     : 2-dimensional array AP(1:LDAP,1:MB); the trans-   *      
!              formed coefficient matrix as put out by BANDP     *      
!     LDAP   : leading dimension of AP as defined in the calling *      
!              program                                           *      
!     MB     : column dimension of AP which must be at least     *      
!                         M + MIN(ML,MU) + 1                     *      
!     N      : order of the matrix                               *      
!     ML     : number of lower co-diagonals                      *      
!     MU     : number of upper co-diagonals                      *      
!     B      : N-vector B(1:N); the right hand side              *      
!     IP     : N-vector IP(1:N); vector containing pivot infor-  *      
!              mation. Should be used unchanged from BANDP.      *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETER:                                          *      
!     =================                                          *      
!     B            as above;  the solution                       *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Elmar Pohl                                         *      
!  date     : 01.16.1992                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      DIMENSION AP (LDAP, MB), IP (N), B (N) 
      LOGICAL UT 
!                                                                       
      M1 = ML + 1 
      M = MU + M1 
      MMB = M + MIN0 (ML, MU) 
      UT = ML.LT.MU 
      IF (UT) THEN 
         IBEG = 1 
         IEND = N - 1 
         ISTEP = 1 
         KBEG = 1 
      ELSE 
         IBEG = N 
         IEND = 2 
         ISTEP = - 1 
         KBEG = - 1 
      ENDIF 
      DO 10 I = IBEG, IEND, ISTEP 
         IF (IP (I) .NE.0) THEN 
            H = B (I) 
            B (I) = B (I + IP (I) ) 
            B (I + IP (I) ) = H 
         ENDIF 
         IF (UT) THEN 
            KEND = MIN0 (ML, N - I) 
         ELSE 
            KEND = MAX0 (1 - I, - MU) 
         ENDIF 
         DO 20 K = KBEG, KEND, ISTEP 
            B (K + I) = B (K + I) - AP (K + I, M1 - K) * B (I) 
   20    END DO 
   10 END DO 
      B (IEND+ISTEP) = B (IEND+ISTEP) / AP (IEND+ISTEP, M1) 
      DO 30 I = IEND, IBEG, - ISTEP 
         IF (UT) THEN 
            KEND = MIN0 (N - I, MU + ML) 
         ELSE 
            KEND = MAX0 (1 - I, - MU - ML) 
         ENDIF 
         DO 40 K = KBEG, KEND, ISTEP 
            KM = K + M1 
            IF (KM.LE.0) KM = KM + MMB 
            B (I) = B (I) - AP (I, KM) * B (K + I) 
   40    END DO 
         B (I) = B (I) / AP (I, M1) 
   30 END DO 
      RETURN 
      END SUBROUTINE BANDS                          
