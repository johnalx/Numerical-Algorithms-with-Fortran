      PROGRAM TEST 
!                                            4.24.1993  ( Dubois Guido )
!***********************************************************************
!                                                                       
!  Test program for subroutine DIFROM.                                  
!                                                                       
!  Test example:                                                        
!  =============                                                        
!                                                                       
!  Function:  F(X) = 3. * SIN(X),  for  X = 0., PI/20., ..., PI         
!                                                                       
!                                                                       
![       X         RESULT     ERROR EST.   ACTUAL ERR. NEND    HEND    ]
![  .000000D+00  .300000D+01  .218181D-11  .888178D-15  4  .196350D-01 ]
![  .157080D+00  .296307D+01  .215516D-11  .888178D-15  4  .196350D-01 ]
![  .314159D+00  .285317D+01  .207523D-11  .133227D-14  4  .196350D-01 ]
![  .471239D+00  .267302D+01  .194422D-11  .444089D-15  4  .196350D-01 ]
![  .628319D+00  .242705D+01  .176525D-11  .133227D-14  4  .196350D-01 ]
![  .785398D+00  .212132D+01  .154277D-11  .444089D-15  4  .196350D-01 ]
![  .942478D+00  .176336D+01  .128253D-11  .444089D-15  4  .196350D-01 ]
![  .109956D+01  .136197D+01  .990985D-12  .168754D-13  4  .196350D-01 ]
![  .125664D+01  .927051D+00  .674238D-12  .344169D-14  4  .196350D-01 ]
![  .141372D+01  .469303D+00  .341505D-12  .971445D-14  4  .196350D-01 ]
![  .157080D+01  .000000D+00  .000000D+00  .183691D-15  2  .785398D-01 ]
![  .172788D+01 -.469303D+00  .341505D-12  .100475D-13  4  .196350D-01 ]
![  .188496D+01 -.927051D+00  .674238D-12  .388578D-14  4  .196350D-01 ]
![  .204204D+01 -.136197D+01  .990541D-12  .177636D-14  4  .196350D-01 ]
![  .219911D+01 -.176336D+01  .128253D-11  .131006D-13  4  .196350D-01 ]
![  .235619D+01 -.212132D+01  .154277D-11  .932587D-14  4  .196350D-01 ]
![  .251327D+01 -.242705D+01  .176525D-11  .932587D-14  4  .196350D-01 ]
![  .267035D+01 -.267302D+01  .194422D-11  .186517D-13  4  .196350D-01 ]
![  .282743D+01 -.285317D+01  .207523D-11  .159872D-13  4  .196350D-01 ]
![  .298451D+01 -.296307D+01  .215516D-11  .155431D-13  4  .196350D-01 ]
![  .314159D+01 -.300000D+01  .218181D-11  .155431D-13  4  .196350D-01 ]
!                                                                       
!***********************************************************************
!                                                                       
      INTEGER NMAX 
      PARAMETER (NMAX = 100) 
!                                                                       
      DOUBLEPRECISION D (1:NMAX), A, B, H, EPS, X, FKT, RES, SCHETZ,    &
      HEND, FSTR, EXF                                                   
      INTEGER N, IANZ, NEND, IFEHL 
!                                                                       
      EXTERNAL FKT 
!                                                                       
      DATA N / 10 / 
      DATA IANZ / 20 / 
      DATA EPS / 1.0D-10 / 
!                                                                       
      A = 0.0D0 
      B = 4.0D0 * DATAN (1.0D0) 
      H = (B - A) / IANZ 
      WRITE ( *, 900) 
      DO 10 X = A, B + .5D0 * H, H 
!                                                                       
         CALL DIFROM (FKT, X, EPS, N, H, RES, SCHETZ, NEND, HEND, IFEHL,&
         D)                                                             
         IF (IFEHL.EQ.1) THEN 
            WRITE ( * , * ) 'ERROR: IFEHL = 1' 
            GOTO 10 
         ENDIF 
         EXF = DABS (RES - FSTR (X) ) 
         WRITE ( *, 910) X, RES, SCHETZ, EXF, NEND, HEND 
   10 END DO 
      STOP 
  900 FORMAT(1X,'C[       X         RESULT     ERROR',                  &
     &          ' EST.   ACTUAL ERR. NEND    HEND',T73,']*')            
  910 FORMAT(1X,'C[',4(1X,D12.6),1X,I2,1X,D12.6,T73,']*') 
      END PROGRAM TEST                              
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION FKT (X) 
      DOUBLEPRECISION X 
      FKT = 3.0D0 * DSIN (X) 
      RETURN 
      END FUNCTION FKT                              
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION FSTR (X) 
      DOUBLEPRECISION X 
      FSTR = 3.0D0 * DCOS (X) 
      RETURN 
      END FUNCTION FSTR                             
