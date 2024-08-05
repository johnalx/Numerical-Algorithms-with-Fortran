C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 6.2.4}
C[  {Brown's Method for Nonlinear Systems}
C[  {Brown's Method for Nonlinear Systems}*)
C[LE*)
      SUBROUTINE BROWN (N,FCT,X0,EPS,EPSM,IOUT,MAXIT,IHF,LDIHF,
C[IX{BROWN}*)
     1                  HF,LDHF,X1,NUMIT,IERR)
C
C*****************************************************************
C                                                                *
C SUBROUTINE BROWN determines the zero of a nonlinear system of  *
C equations consisting of N equations in N unknowns by applying  *
C the Brown-method.                                              *
C[BE*)
C                                                                *
C In order to perfornm one iteration step with Brown's algorithm *
C the SUBROUTINE ITER4 is called. The iteration is continued     *
C until the given maximal number of iteration steps is reached.  *
C The iteration is stopped if one of the following break-off     *
C criteria is met:                                               *
C - if the relative change of two successive iterative solutions *
C   is smaller than EPS.                                         *
C - if the functional value is smaller or equal to EPSM.         *
C - if the limiting accuracy is reached.                         *
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C N       : number of equations and unknowns                     *
C FCT     : function subroutine, given as follows:               *
C                 SUBROUTINE FCT (K,X,F)                         *
C                 DOUBLE PRECISION X(2)                          *
C                 GO TO (1,2) K                                  *
C               1 F = F1(X(1),X(2))                              *
C                 RETURN                                         *
C               2 F = F2(X(1),X(2))                              *
C                 RETURN                                         *
C                 END                                            *
C X0      : N-vector X0(1:N); starting vector with N components  *
C EPS     : desired accuracy                                     *
C EPSM    : machine constant                                     *
C IOUT    : IOUT=1 after each iteration step the program puts    *
C                  out the difference of the last approximate    *
C                  solutions, the current approximate solution   *
C                  and the functional value                      *
C           IOUT=0 no output is provided, except upon            *
C                  convergence                                   *
C MAXIT   : maximum number of iteration steps allowed            *
C IHF     : 2-dim. array IHF(1:LDIHF,1:N+1); ) auxiliary arrays  *
C HF      : 2-dim. array   HF(1:LDHF,1:N+3); ) used as           *
C                                            ) intermediate      *
C                                            ) storage space     *
C LDIHF   : leading dimension of IHF as defined in the calling   *
C           program. LDIHF >= N                                  *
C LDHF    : leading dimension of HF as defined in the calling    *
C           program. LDHF >= N                                   *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C X1      : N-vector X1(1:N); the approximate solution vector    *
C NUMIT   : number of iterations executed                        *
C IERR    : error parameter                                      *
C           =0; no error, zero has been found                    *
C           =1; after MAXIT steps the desired accuracy was not   *
C               reached                                          *
C           =2; singular matrix                                  *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: ITER4, SUBST                            *
C                                                                *
C                                                                *
C  sources : Brown K.M.: A quadratically convergent Newton-like  *
C            method based upon Gaussian elimination,             *
C            Siam J.Numer.Anal.vol 6 (1969), 560 - 569           *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Johannes Karfusehr                                 *
C  Editor   : Thomas Eul                                         *
C  Date     : 08.21.1985                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C   declarations
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL FCT
      INTEGER N, IERR, LDIHF, LDHF
      DIMENSION X0(N),X1(N)
      DIMENSION IHF(LDIHF,N+1),HF(LDHF,N+3)
C
C   local variables
C
      LOGICAL SING
C
C   initializing
C
      SING=.FALSE.
      ICRIT=0
      DELTA0=1.0D-2
      DO 10 J=1,N
         HF(J,N+2)=X0(J)
   10 CONTINUE
      IF(IOUT .EQ. 1) WRITE(*,1000)
C
C   iteration
C
      DO 20 M=1,MAXIT
         CALL ITER4(N,FCT,EPSM,IHF,LDIHF,HF,LDHF,X1,SING)
         NUMIT=M
C
C   after each iteration step we output the current
C   starting vector and the next iterate
C
         IF(IOUT .EQ. 1) THEN
            WRITE(*,3000) M
            IF(.NOT. SING) THEN
C
               WRITE(*,4000)
               DO 30 J=1,N
                  CALL FCT(J,X1,FVAL)
                  WRITE(*,5000) X1(J)-HF(J,N+2),J,X1(J),FVAL
   30          CONTINUE
            ELSE
               WRITE(*,7000)
            ENDIF
         ENDIF
C
         IF(.NOT. SING) THEN
C
C   test the break-off criterion
C
C   test the relative change of the iterates
C
            DO 40 I=1,N
               RELF=(X1(I)-HF(I,N+2))/(HF(I,N+2)+EPS)
               IF(DABS(RELF) .GE. EPS) THEN
                  GOTO 41
               ENDIF
   40       CONTINUE
            ICRIT=1
            GOTO 21
   41       CONTINUE
C
C   test the functional value
C
            DO 50 I=1,N
               CALL FCT(I,X1,FVAL)
               IF(DABS(FVAL) .GT. EPSM) GOTO 51
   50       CONTINUE
            ICRIT=2
            GOTO 21
   51       CONTINUE
C
C   test the limiting accuracy
C
            DELTA1=DABS(X1(1)-HF(1,N+2))
            DO 60 I=2,N
               DELTA1=DMAX1(DELTA1,DABS(X1(I)-HF(I,N+2)))
   60       CONTINUE
            IF(DELTA1 .LE. 1.0D-3) THEN
               IF(DELTA0 .LE. DELTA1) THEN
                  ICRIT=3
                  GOTO 21
               ENDIF
            ENDIF
            DELTA0=DELTA1
            IF(NUMIT .LT. MAXIT) THEN
               DO 70 I=1,N
                  HF(I,N+2)=X1(I)
   70          CONTINUE
            ENDIF
         ELSE
            GOTO 21
         ENDIF
   20 CONTINUE
   21 CONTINUE
      IF(SING) THEN
         IERR=2
      ELSE
         IF(ICRIT .EQ. 0) THEN
            IERR=1
         ELSE
            IERR=0
         ENDIF
      ENDIF
C
 1000 FORMAT ('1')
 3000 FORMAT (1X,I4,'-TH ITERATION STEP :')
 4000 FORMAT (1X,6X,'DIFFERENCE',7X,2X,'KOMP',2X,6X,'APPROXIMATION',7X,
     1        2X,4X,'FUNCTIONAL VALUE')
 5000 FORMAT (1X,D22.15,2X,I4,2X,D22.15,2X,D22.15)
 7000 FORMAT(1X,'THE JACOBI MATRIX IS SINGULAR')
C
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE ITER4 (N,FCT,EPSM,IHF,LDIHF,HF,LDHF,X1,SING)
C[IX{ITER4}*)
C
C*****************************************************************
C                                                                *
C SUBROUTINE ITER4 computes an approximation by Brown's          *
C algorithm.                                                     *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETER:                                               *
C ================                                               *
C N       : number of equations or unknowns                      *
C FCT     : function subroutine (compare with SUBROUTINE BROWN)  *
C EPSM    : machine constant                                     *
C IHF     : 2-dim. array IHF(1:LDIHF,1:N+1); ) auxiliary array   *
C HF      : 2-dim. array   HF(1:LDHF,1:N+3); ) providing storage *
C                                            ) space             *
C LDIHF   : leading dimension of IHF as defined in the calling   *
C           program. LDIHF >= N                                  *
C LDHF    : leading dimension of HF as defined in the calling    *
C           program. LDHF >= N                                   *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C X1      : N-vector X1(1:N); the approximation vector           *
C SING    : error parameter; indicates whether the Jacobi matrix *
C           is singular or not                                   *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: SUBST                                   *
C                                                                *
C                                                                *
C  sources : Brown K.M.: A quadratically convergent Newton-like  *
C            method based upon Gaussian elimination,             *
C            Siam J.Numer.Anal.,vol 6 (1969), 560 - 569          *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Johannes Karfusehr                                 *
C  editor   : Thomas Eul                                         *
C  date     : 08.21.1985                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C     declarations
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SING
      INTEGER N, LDIHF, LDHF
      DIMENSION IHF(LDIHF,N+1), HF(LDHF,N+3), X1(N)
C
C   initializing
C
      DO 10 J=1,N
         IHF(1,J)=J
         X1(J)=HF(J,N+2)
   10 CONTINUE
C
C   linearization of the K-th coordinate function
C
      DO 20 K=1,N
         ICOUNT=0
         FACTOR=0.001D0
         DO 30 J=1,3
            IF(K .GT. 1) CALL SUBST(K,N,IHF,LDIHF,HF,LDHF,X1)
            CALL FCT(K,X1,F)
C
C   determining the I-th discretization stepsize
C   and the I-th difference quotient
C
            DO 40 I=K,N
               ITEMP=IHF(K,I)
               HOLD=X1(ITEMP)
               H=FACTOR*HOLD
               IF(DABS(H) .LE. EPSM) H=0.001D0
               X1(ITEMP)=HOLD+H
               IF(K .GT. 1) CALL SUBST(K,N,IHF,LDIHF,HF,LDHF,X1)
               CALL FCT(K,X1,FPLUS)
               X1(ITEMP)=HOLD
               HF(ITEMP,N+3)=(FPLUS-F)/H
               IF(DABS(HF(ITEMP,N+3)) .LE. EPSM) THEN
                  ICOUNT=ICOUNT+1
               ELSE
                  IF(DABS(F/HF(ITEMP,N+3)) .GE. 1.0D20) ICOUNT=ICOUNT+1
               ENDIF
   40       CONTINUE
            IF(ICOUNT .LE. N-K) THEN
               SING=.FALSE.
               GOTO 31
            ELSE
               SING=.TRUE.
               FACTOR=FACTOR*10.0D0
               ICOUNT=0
            ENDIF
   30    CONTINUE
   31    CONTINUE
C
         IF(.NOT. SING) THEN
            IF(K .LT. N) THEN
               KMAX=IHF(K,K)
C
C   determining the difference quotient of largest magnitude
C
               DERMAX=DABS(HF(KMAX,N+3))
               KPLUS=K+1
               DO 50 I=KPLUS,N
                  JSUB=IHF(K,I)
                  TEST=DABS(HF(JSUB,N+3))
                  IF(TEST .LT. DERMAX) THEN
                     IHF(KPLUS,I)=JSUB
                  ELSE
                     IHF(KPLUS,I)=KMAX
                     KMAX=JSUB
                  ENDIF
   50          CONTINUE
               IF(DABS(HF(KMAX,N+3)) .LE. EPSM) SING=.TRUE.
               IHF(K,N+1)=KMAX
               IF(.NOT. SING) THEN
                  HF(K,N+1)=0.0D0
C
C   solving the K-th equation for XMAX
C
                  DO 60 J=KPLUS,N
                     JSUB=IHF(KPLUS,J)
                     HF(K,JSUB)=-HF(JSUB,N+3)/HF(KMAX,N+3)
                     HF(K,N+1)=HF(K,N+1)+HF(JSUB,N+3)*X1(JSUB)
   60             CONTINUE
                  HF(K,N+1)=(HF(K,N+1)-F)/HF(KMAX,N+3)+X1(KMAX)
               ELSE
                  GOTO 21
               ENDIF
            ELSE
C
C   solving the N-th coordinate function by use of the
C   discrete Newton-method for one variable
C
               IF(DABS(HF(ITEMP,N+3)) .LE. EPSM) THEN
                  SING=.TRUE.
               ELSE
                  HF(K,N+1)=0.0D0
                  KMAX=ITEMP
                  HF(K,N+1)=(HF(K,N+1)-F)/HF(KMAX,N+3)+X1(KMAX)
               ENDIF
            ENDIF
         ELSE
            GOTO 21
         ENDIF
   20 CONTINUE
   21 CONTINUE
      IF(.NOT. SING) THEN
C
C   determining of the approximate solution by backsubstitution
C
         X1(KMAX)=HF(N,N+1)
         IF(N .GT. 1) CALL SUBST(N,N,IHF,LDIHF,HF,LDHF,X1)
      ENDIF
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE SUBST (K,N,IHF,LDIHF,HF,LDHF,X1)
C[IX{SUBST}*)
C
C*****************************************************************
C                                                                *
C SUBROUTINE SUBST solves a linear system of equations.          *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C K       : index of the coordinate function                     *
C N       : total number of equations and unknowns               *
C IHF     : 2-dim. array IHF(1:LDIHF,1:N+1); ) auxiliary arrays, *
C HF      : 2-dim. array   HF(1:LDHF,1:N+3); ) providing storage *
C                                            ) space             *
C LDIHF   : leading dimension of IHF as defined in the calling   *
C           program. LDIHF >= N                                  *
C LDHF    : leading dimension of HF as defined in the calling    *
C           program. LDHF >= N                                   *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C X1      : N-vector X1(1:N); approximate solution vector        *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C                                                                *
C  source : Brown K.M.: A quadratically convergent Newton-like   *
C           method based upon Gaussian elimination,              *
C           Siam J. Numer. Anal., vol 6 (1969), 560 - 569        *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Johannes Karfusehr                                 *
C  editor   : Thomas Eul                                         *
C  date     : 08.21.1985                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C   declarations
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER K, N, LDIHF, LDHF
      DIMENSION IHF(LDIHF,N+1), HF(LDHF,N+3), X1(N)
C
      DO 10 KM=K,2,-1
         KMAX=IHF(KM-1,N+1)
         X1(KMAX)=0.0D0
         DO 20 J=KM,N
            JSUB=IHF(KM,J)
            X1(KMAX)=X1(KMAX)+HF(KM-1,JSUB)*X1(JSUB)
   20    CONTINUE
         X1(KMAX)=X1(KMAX)+HF(KM-1,N+1)
   10 CONTINUE
      RETURN
      END
