      PROGRAM TEST
C
C*****************************************************************
C                                                                *
C  Testprogram for SUBROUTINE SHOUSE.                            *
C                                                                *
C  Solve the linear system  A * X = Y  with the Hilbert matrix   *
C                                                                *
C  A(I,K)=1./(I+K-1) ,  I,K = 1, .., N  for  N = 1, .., 10       *
C                                                                *
C  and  Y(I)=A(I,1)+A(I,2)+...+A(I,N) ,  I=1, ..., N.            *
C                                                                *
C  Solution:  X(I) = 1. , I = 1, .., N.                          *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(0:10,0:10),X(0:10),D(0:10)
      CHARACTER FORM*30,TEX1*2
      DO 10 N=1,10
         DO 20 I=0,N-1
            SUM=0.D0
            DO 30 K=0,N-1
               A(I,K)=1.D0/(I+K+1.D0)
               SUM=SUM+A(I,K)
   30       CONTINUE
            A(I,N)=SUM
   20    CONTINUE
         WRITE(TEX1,'(I2)')N
         FORM='(1X,'//TEX1//'(F5.3,1X),2X,F5.3)'
         DO 40 I=0,N-1
            WRITE(*,FORM)(A(I,K),K=0,N)
   40    CONTINUE
         CALL SHOUSE(A,10,N-1,N-1,D,X,IFEHL)
         WRITE(*,'(//,1X,''SOLUTION:'')')
         WRITE(*,*)'*********'
         WRITE(*,*)
         WRITE(*,'(/,1X,''IERR = '',I3,//)')IFEHL
         IF(IFEHL.EQ.0)THEN
            DO 50 I=0,N-1
               WRITE(*,'(1X,D20.14,10X,D20.14)')X(I),DABS(X(I)-1.)
   50       CONTINUE
         END IF
   10 CONTINUE
      STOP
      END
