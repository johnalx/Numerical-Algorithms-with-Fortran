      PROGRAM TEST
C
C*****************************************************************
C                                                                *
C  Testprogram for subroutines GAUSRS, GAUSSP, GAUSSS, HACOND,   *
C  and POSTIT.                                                   *
C                                                                *
C  Solve the linear system  A * XL = RS  with the Hilbert matrix *
C                                                                *
C  A(I,K)=1./(I+K-1) ,  I,K = 1, .., N  for  N = 1, ..., 10      *
C                                                                *
C   and the right hand sides                                     *
C                                                                *
C  RS(I,K)=K*(A(I,1)+A(I,2)+...+A(I,N)),  I = 1,..,N, K = 1,..,5.*
C                                                                *
C  Solution:  XL(I,K) = K , I = 1, .., N, K = 1, .. , 5.         *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1:10,1:10),A0(1:10,1:10),RS(1:10,1:5),XL(1:10,1:5),
     +     D(1:10),Z(1:10),RES(1:10),R(1:10)
      INTEGER IPIVOT(1:10)
      CHARACTER FORM*40,TEX1*2
      EPS=0.1D-12
      DO 10 N=1,10
         DO 20 I=1,N
            SUM=0.D0
            DO 30 K=1,N
               A(I,K)=1.D0/(I+K-1.D0)
               A0(I,K)=A(I,K)
               SUM=SUM+A(I,K)
   30       CONTINUE
            RS(I,1)=SUM
   20    CONTINUE
         DO 40 K=2,5
            DO 50 I=1,N
               RS(I,K)=K*RS(I,1)
   50       CONTINUE
   40    CONTINUE
         WRITE(TEX1,'(I2)')N
         FORM='(1X,'//TEX1//'(F5.3,1X),5(1X,F6.3))'
         DO 60 I=1,N
            WRITE(*,FORM)(A(I,K),K=1,N),(RS(I,K),K=1,5)
   60    CONTINUE
         CALL GAUSRS(N,A,10,5,RS,XL,MARKE,D,IPIVOT)
         WRITE(*,'(//,1X,''SOLUTION:'')')
         WRITE(*,*)'*********'
         WRITE(*,'(/,1X,''MARK = '',I3)')MARKE
         IF(MARKE.NE.0)THEN
            DO 70 K=1,5
               WRITE(*,'(/)')
               DO 80 I=1,N
                  WRITE(*,'(1X,D20.14,10X,D20.14)')XL(I,K),
     +                                             DABS(XL(I,K)-K)
   80          CONTINUE
   70       CONTINUE
            CALL HACOND(N,A0,A,10,MARKE,HKOND)
            WRITE(*,'(//,1X,''CONDITION NUMBER: '',E20.14)')HKOND
            DO 90 K=1,5
               CALL POSTIT(N,A0,A,10,IPIVOT,RS(1,K),XL(1,K),EPS,10,
     +                     ITANZ,IFEHL,Z,R,RES)
               WRITE(*,'(//,1X,''BREAK-OFF CONDITION: '',I3)')IFEHL
               IF(IFEHL.EQ.0)THEN
                  WRITE(*,900)ITANZ,EPS
                  DO 100 I=1,N
                     WRITE(*,'(1X,D20.14,10X,D20.14)')XL(I,K),
     +                                                DABS(XL(I,K)-K)
  100             CONTINUE
               END IF
   90       CONTINUE
         END IF
         WRITE(*,'(1H1)')
   10 CONTINUE
      STOP
  900 FORMAT(/,1X,'NUMBER OF ITERATIONS: ',I3,10X,'ERRONEOUS BOUND: ',
     +       D20.14,//)
      END
