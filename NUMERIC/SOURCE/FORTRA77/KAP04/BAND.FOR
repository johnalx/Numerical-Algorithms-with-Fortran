C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 4.13}
C[  {Linear Systems with Band Matrices}
C[  {Linear Systems with Band Matrices}*)
C[LE*)
C[FE{}
C[  {Linear Systems with Band Matrices using Pivots}
C[  {Linear Systems with Band Matrices using Pivots}*)
C[LE*)
      SUBROUTINE BAND(AP,LDAP,MB,N,ML,MU,B,IERR,IP)
C[IX{BAND}*)
C
C*****************************************************************
C                                                                *
C     BAND solves a linear system of equations for a banded      *
C     system matrix AP and a right hand side B. This program     *
C     should not be used if several systems with the same system *
C     matrix AP are to be solved. In this case it is better to   *
C     call BANDP below once and BANDS separately for each right  *
C     hand side.                                                 *
C[BE*)
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     AP     : 2-dimensional array AP(1:LDAP,1:MB); the system   *
C              matrix in condensed form. If A denotes the        *
C              original system matrix, then                      *
C                 A(I,K)=AP(I,ML+1+K-I)                          *
C              holds for all indices I, K inside the diagonal    *
C              band of A. Thus the co-diagonals of A become      *
C              columns of AP; its rows are condensed to the left.*
C     LDAP   : leading dimension of AP as defined in the calling *
C              program                                           *
C     MB     : number of columns of AP. It is at least           *
C                 M + MIN(ML,MU) + 1                             *
C              where  M = ML+MU+1  denotes the band width of AP  *
C     N      : size of the matrix AP                             *
C     ML     : number of lower co-diagonals                      *
C     MU     : number of upper co-diagonals                      *
C     B      : N-vector  B(1:N); the right hand side             *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     B      : the solution (if IERR .NE. 0)                     *
C     IERR   : error parameter;                                  *
C                 IERR = 0 : no error                            *
C                 IERR.NE.0: matrix AP is numerically singular   *
C     AP     : as above; tranformed band matrix and transforma-  *
C              tion coefficients. May be used as input for BANDS *
C              in order to find solutions for additional right   *
C              hand sides.                                       *
C     IP     : N-vector IP(1:N); this vector contains pivoting   *
C              information that is needed by BANDS.              *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: BANDP, BANDS, MACHPD                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Elmar Pohl                                         *
C  date     : 01.16.1992                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION AP(LDAP,MB),B(N),IP(N)
C
      CALL BANDP(AP,LDAP,MB,N,ML,MU,IP,ISIG,IERR)
      IF(IERR .NE. 0) RETURN
      CALL BANDS(AP,LDAP,MB,N,ML,MU,B,IP)
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BANDP(AP,LDAP,MB,N,ML,MU,IP,ISIG,IERR)
C[IX{BANDP}*)
C
C*****************************************************************
C                                                                *
C     BANDP transforms a banded matrix AP to upper or lower      *
C     triangular form by applying a Gauss algorithm with         *
C     partial pivot search.                                      *
C[BE*)
C     The triangular matrix and the                              *
C     transformation matrix are overwritten onto the original    *
C     array. If IERR=0, the subroutine BANDS may be applied in   *
C     order to solve several linear systems of equations.        *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     AP     : 2-dimensional array AP(1:LDAP,1:MB); the band     *
C              matrix in condensed form. The rows of this array  *
C              contain the rows of the matrix AP, shifted to the *
C              left. The diagonals of AP, however, appear as     *
C              columns of AP:                                    *
C                    A(I,K) = AP(I,ML+1+K-I)                     *
C              where A denotes the matrix in uncondensed form.   *
C     LDAP   : leading dimension of AP as defined in the calling *
C              program                                           *
C     MB     : column dimension of AP. May differ from that used *
C              in the calling program; however, it must be at    *
C              least M + MIN(ML,MU) + 1, where M = MU+ML+1 is    *
C              the band width of AP.                             *
C     N      : order of the matrix AP                            *
C     ML     : number of lower co-diagonals                      *
C     MU     : number of upper co-diagonals                      *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     AP     : as above; triangular matrix and transforming      *
C              matrix                                            *
C     IP     : N-vector IP(1:N); pivoting information which is   *
C              needed for BANDS                                  *
C     ISIG   : = (-1)**IZ, where IZ denotes the number of row    *
C              permutations; to be used for calculating the      *
C              determinant:                                      *
C                                    N                           *
C                 DET(A) = ISIG * PRODUCT AP(I,ML+1)             *
C                                   I=1                          *
C                                                                *
C     IERR   : error parameter;                                  *
C                 IERR = 0 : no error                            *
C                 IERR.GT.0: in the IERR-th step, the pivot      *
C                            element becomes zero, the matrix AP *
C                            is numerically singular             *
C                 IERR.LT.0: row IERR contains only zeros,       *
C                            matrix is singular                  *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: MACHPD                                  *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Elmar Pohl                                         *
C  date     : 01.16.1992                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION AP(LDAP,MB),IP(N)
      LOGICAL UT
      M1 = ML + 1
      M = MU + M1
      MMB = M + MIN0(ML,MU)
C
C     calculate the machine constant
C
      FMACHP = 1.0D0
   80 FMACHP = 0.5D0 * FMACHP
      IF (MACHPD(1.0D0 + FMACHP) .EQ. 1) GOTO 80
      FMACHP = FMACHP * 2.0D0
C
C     determinaing the relative error bound
      EPS = 4.0D0 * FMACHP
C
C     Due to the pivot search the band width may increase. When
C     transforming to upper triangular form, ML additional
C     upper co-diagonals will appear, while transforming
C     to lower triangular form creates MU additional lower
C     co-diagonals. The program automatically selects the
C     most advantageous triangularization.
C
      UT = ML .LT. MU
C
C     UT = .TRUE.   : transforming to upper triangular form
C     UT = .FALSE.  : transforming to lower triangular form
C
C     the additional co-diagonals are stored in the (M+1)th
C     to MMB-th column of AP. At first they are initialized
C     to zero.
C     Simultaneously we compute and store the summed absolute
C     values of each row in the MMB + 1st column. (To be used
C     for checking nonsingularity)
C
      DO 10 I=1,N
         S = 0.0D0
         DO 90 K=MAX0(1,I-ML),MIN0(N,I+MU)
            S = S + DABS(AP(I,M1+K-I))
   90    CONTINUE
         IF (S .EQ. 0.0D0) THEN
            IERR = -I
            RETURN
         ENDIF
         AP(I,MMB+1) = S
         DO 20 K=M+1,MMB
            AP(I,K) = 0.0D0
   20    CONTINUE
   10 CONTINUE
C
      ISIG = 1
      IERR = 0
C
C     the index bounds and the direction of the sweep in the
C     following loops are set depending on UT in such a manner
C     that an upper or a lower triangular form is obtained.
C
      IF (UT) THEN
         IBEG = 1
         IEND = N-1
         ISTEP = 1
         KBEG = 1
         JBEG = 1
      ELSE
         IBEG = N
         IEND = 2
         ISTEP = -1
         KBEG = -1
         JBEG = -1
      ENDIF
C
C     loop over all rows
C
      DO 30 I=IBEG,IEND,ISTEP
         IF (UT) THEN
            KEND = MIN0(ML,N-I)
         ELSE
            KEND = MAX0(1-I,-MU)
         ENDIF
C
C        pivot search: searching for the element of largest
C        magnitude in the I-th column; store its index
C        (relative to I) in IV.
C
         IV = 0
         AM = DABS(AP(I,M1))
         DO 40 K=KBEG,KEND,ISTEP
            IF (DABS(AP(K+I,M1-K)) .GT. AM) THEN
               AM = DABS(AP(K+I,M1-K))
               IV = K
            ENDIF
   40    CONTINUE
         IF (AM/AP(I+IV,MMB+1) .LE. EPS) THEN
C
C           element of largest magnitude of the
C           column is zero, thus the matrix is singular
C
            IERR = I
            RETURN
         ENDIF
         IP(I) = IV
         IF (UT) THEN
            KJEND = MIN0(IV+MU,N-I)
         ELSE
            KJEND = MAX0(1-I,IV-ML)
         ENDIF
         IF (IV .NE. 0) THEN
C
C           interchange row I with row I+IV
C
            ISIG = -ISIG
            DO 50 K=0,KJEND,ISTEP
               KM = K+M1
               IF (KM .LE. 0) KM = KM+MMB
               H = AP(I,KM)
               AP(I,KM) = AP(I+IV,K+M1-IV)
               AP(I+IV,K+M1-IV) = H
   50       CONTINUE
            H = AP(I,MMB+1)
            AP(I,MMB+1) = AP(I+IV,MMB+1)
            AP(I+IV,MMB+1) = H
         ENDIF
C
C        elimination;
C        loop including all rows below row I:
C
         DO 60 K=KBEG,KEND,ISTEP
            AP(K+I,M1-K) = AP(K+I,M1-K) / AP(I,M1)
C
C           loop including all columns to the right of column I
C
            DO 70 J=JBEG,KJEND,ISTEP
               JK = J+M1-K
               JM = J+M1
C
C              With our condensing scheme, additional upper
C              co-diagonals would have to be stored to the right
C              of the M-th column and additional lower ones to the
C              left of the first column.
C              The latter would mean that we would need access to
C              storage space in front of AP(1,1).
C              In order to prevent this, MMB is added to all
C              column subscripts < 1. Thus additional lower
C              co-diagonals are also stored to the right of AP.
C
               IF (JK .LE. 0) JK = JK+MMB
               IF (JM .LE. 0) JM = JM+MMB
               AP(K+I,JK) = AP(K+I,JK)-AP(K+I,M1-K)*AP(I,JM)
   70       CONTINUE
   60    CONTINUE
   30 CONTINUE
      IP(IEND+ISTEP) = 0
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BANDS(AP,LDAP,MB,N,ML,MU,B,IP)
C[IX{BANDS}*)
C
C*****************************************************************
C                                                                *
C     After calling BANDP, the subroutine BANDS solves the linear*
C     system of equations consisting of a banded, condensed      *
C     coefficient matrix AP and a right hand side B.             *
C[BE*)
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     AP     : 2-dimensional array AP(1:LDAP,1:MB); the trans-   *
C              formed coefficient matrix as put out by BANDP     *
C     LDAP   : leading dimension of AP as defined in the calling *
C              program                                           *
C     MB     : column dimension of AP which must be at least     *
C                         M + MIN(ML,MU) + 1                     *
C     N      : order of the matrix                               *
C     ML     : number of lower co-diagonals                      *
C     MU     : number of upper co-diagonals                      *
C     B      : N-vector B(1:N); the right hand side              *
C     IP     : N-vector IP(1:N); vector containing pivot infor-  *
C              mation. Should be used unchanged from BANDP.      *
C                                                                *
C                                                                *
C     OUTPUT PARAMETER:                                          *
C     =================                                          *
C     B            as above;  the solution                       *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Elmar Pohl                                         *
C  date     : 01.16.1992                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION AP(LDAP,MB),IP(N),B(N)
      LOGICAL UT
C
      M1 = ML+1
      M = MU+M1
      MMB = M+MIN0(ML,MU)
      UT = ML .LT. MU
      IF (UT) THEN
         IBEG = 1
         IEND = N-1
         ISTEP = 1
         KBEG = 1
      ELSE
         IBEG = N
         IEND = 2
         ISTEP = -1
         KBEG = -1
      ENDIF
      DO 10 I=IBEG,IEND,ISTEP
         IF (IP(I) .NE. 0) THEN
            H = B(I)
            B(I) = B(I+IP(I))
            B(I+IP(I)) = H
         ENDIF
         IF (UT) THEN
            KEND = MIN0(ML,N-I)
         ELSE
            KEND = MAX0(1-I,-MU)
         ENDIF
         DO 20 K=KBEG,KEND,ISTEP
            B(K+I) = B(K+I)-AP(K+I,M1-K)*B(I)
   20    CONTINUE
   10 CONTINUE
      B(IEND+ISTEP) = B(IEND+ISTEP) / AP(IEND+ISTEP,M1)
      DO 30 I=IEND,IBEG,-ISTEP
         IF (UT) THEN
            KEND = MIN0(N-I,MU+ML)
         ELSE
            KEND = MAX0(1-I,-MU-ML)
         ENDIF
         DO 40 K=KBEG,KEND,ISTEP
            KM = K+M1
            IF (KM .LE. 0) KM = KM+MMB
            B(I) = B(I)-AP(I,KM)*B(K+I)
   40    CONTINUE
         B(I) = B(I)/AP(I,M1)
   30 CONTINUE
      RETURN
      END
