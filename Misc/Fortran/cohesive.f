C     -----------------------------------------------------
C     PROGRAM COHESIVE:  SEP. 9 1998, JU LI (MIT)
C     LOCATE REACTION PATH FROM SUITABLE REFERENCE STATES
C     AND MAKE BEST USE OF RAW LDA TOTAL ENERGY DATA WITH 
C     THE UNIVERSAL BINDING FUNCTION OR TAYLOR EXPANSION 
C     TO CALCULATE THE COMPOSITE COHESIVE ENERGY CURVE.
C     f77 -o cohesive cohesive.f; cohesive < uon
C     -----------------------------------------------------
      
      PROGRAM COHESIVE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      PARAMETER (MOP=10,MMT=5,MMS=1000,RY_IN_EV=13.605826,
     A       AU_IN_A=0.529177249,EV_IN_JOULE=1.602192D-19,
     A       STRESS_IN_MBAR=EV_IN_JOULE*1D30*1D-11)
      COMMON /A/ VOLUME_LDA(MMS),ENERGY_LDA(MMS),NUM_LDA,IFT,NOP
      DIMENSION NAME_ATOM(MMT),NUM_ATOM(MMT),NAME_REF(MMT),
     A          NUM_REF_ATOM(MMT,MMT),TOT_REF_ATOM(MMT),
     A          REF_LDA(MMT),REF_EXP(MMT),VOLUME_FIT(MMS),
     A          ENERGY_FIT(MMS),RATIO(MMT),IPVT(MMT),Z(MMT),
     A          BOUND(MOP,2),GUESS(MOP)
      DOUBLE PRECISION NUM_REF_ATOM
      CHARACTER LABEL_LDA*45,NAME_ATOM*2,NAME_REF*50,NAME*20
      EXTERNAL FE
      DATA LP /20/
      
      READ (*,'("LDA data set label:")')
      READ (*,'(A45)') LABEL_LDA
      
      READ (*,'(/,"Number of atomic species: ",I)') NUM_SPECIES
      IF (NUM_SPECIES.GT.MMT) THEN 
      PRINT *,'cohesive(): MMT reached.'
      STOP
      ENDIF
      TOT_ATOM = 0.
      DO I = 1,NUM_SPECIES
      READ (*,'("Species: ",A2,I)') NAME_ATOM(I),NUM_ATOM(I)
      RATIO(I) = NUM_ATOM(I)
      TOT_ATOM = TOT_ATOM + NUM_ATOM(I)
      ENDDO

      DO I = 1,NUM_SPECIES
      READ (*,'(/,"Reference state: ",A50)') NAME_REF(I)
      READ (*,'("Composition: ",5F3.0)')
     A     (NUM_REF_ATOM(J,I),J=1,NUM_SPECIES)
      TOT_REF_ATOM(I) = 0.
      DO J = 1,NUM_SPECIES
      TOT_REF_ATOM(I) = TOT_REF_ATOM(I) + NUM_REF_ATOM(J,I)
      ENDDO
      READ (*,'("LDA total energy (eV/atom): ",F)') REF_LDA(I)
      READ (*,'("Exp''t cohesive energy (eV/atom): ",F)') REF_EXP(I)
      ENDDO
      
C     LU factorize the composition matrix. Warning: 
C     NUM_REF_ATOM will be changed after calling.
      CALL DGECO (NUM_REF_ATOM,MMT,NUM_SPECIES,IPVT,RCOND,Z)
      IF (1.0+RCOND.EQ.1.0) THEN 
      PRINT *, 'cohesive(): matrix is singular.'
      STOP
      ENDIF
C     Solve the linear system of equations:
      CALL DGESL (NUM_REF_ATOM,MMT,NUM_SPECIES,IPVT,RATIO,0)
      WRITE (*,'(/,
     A "We can infer target''s cohesive energy from the reaction")')
      WRITE (*,'("target <=", F9.5," * (",A7,")",
     A      4(A2,F9.5," * (",A7,")"))') RATIO(1),NAME_REF(1),
     A      (' +',RATIO(I),NAME_REF(I),I=2,NUM_SPECIES)
      
C     Calculate the reaction energy + reference cohesive energy
      SHIFT = 0.D0
      DO J = 1,NUM_SPECIES
      SHIFT = SHIFT + RATIO(J) * TOT_REF_ATOM(J) / 
     A        TOT_ATOM * (REF_EXP(J) - REF_LDA(J))
      ENDDO
      
C     Load in raw LDA total energy data:
      READ (*,'(/,"LDA total energy (Ry) and cell volume (au^3):")')
      NUM_LDA = 0
 10   READ *, VOLUME, ENERGY
      IF (VOLUME.NE.-9999) THEN 
      NUM_LDA = NUM_LDA + 1
      IF (NUM_LDA.GT.MMS) THEN 
      PRINT *,'cohesive(): MMS reached.'
      STOP
      ENDIF
C     From cell volume (au^3) to volume/atom (A^3):
      VOLUME_LDA(NUM_LDA) = VOLUME * AU_IN_A**3 / TOT_ATOM 
C     From cell energy (Ry) to energy/atom (eV):
      ENERGY_LDA(NUM_LDA) = ENERGY * RY_IN_EV / TOT_ATOM 
C     Calculate the composite cohesive energy targets from 
C     experimental cohesive energy of reference states and
C     LDA predicted reaction energy per particle:
      ENERGY_LDA(NUM_LDA) = ENERGY_LDA(NUM_LDA) + SHIFT
      GOTO 10
      ENDIF
      IF (NUM_LDA.LT.1) THEN 
      PRINT *,'cohesive(): there should be at least one valid data.'
      STOP
      ENDIF
      
      READ (*, '(/,
     A "Fit to (-1: Universal binding  >=0: Taylor expansion):",/,I)')
     A IFT
      IF (IFT.LT.0) THEN 
      NOP = 3
      ELSEIF (IFT.EQ.1) THEN
      PRINT *, 'cohesive(): IFT cannot be 1' 
      STOP
      ELSE
      NOP = IFT+1
      ENDIF
      
C     Find the energy minimum and volume bounds of the raw data:
      A = ENERGY_LDA(1)
      B = VOLUME_LDA(1)
      C = VOLUME_LDA(1)
      D = VOLUME_LDA(1)
      DO J = 2,NUM_LDA
      IF (A.GT.ENERGY_LDA(J)) THEN 
      A = ENERGY_LDA(J)
      B = VOLUME_LDA(J)
      ENDIF
      IF (C.GT.VOLUME_LDA(J)) C = VOLUME_LDA(J)
      IF (D.LT.VOLUME_LDA(J)) D = VOLUME_LDA(J)
      ENDDO
      
      WRITE (*,'(/,
     A "Raw data volume range = [",F9.6,",",F9.6,"] A^3/atom,")') C,D
      WRITE (*,'("         minimum energy is", F11.5,
     A             " Ry/cell,   or",F11.6," eV/atom,")')
     A       A*TOT_ATOM/RY_IN_EV, A
      WRITE (*,'("     equilibrium volume is",F11.5,
     A           " au^3/cell, or",F11.6," A^3/atom.")')
     A      B*TOT_ATOM/AU_IN_A**3, B
      
      READ (*, '(/,"Parameters bound (-9999 means default):")')
      READ (*, '("Lowest Total Cell Energy (Ry):")')
      READ *, BOUND(1,2),GUESS(1),BOUND(1,1)
C     Assign defaults for the bound and guess:
      IF (GUESS(1).EQ.-9999) THEN 
      GUESS(1) = A
      ELSE
      GUESS(1) = GUESS(1) * RY_IN_EV / TOT_ATOM
      ENDIF
      IF (BOUND(1,2).EQ.-9999) THEN 
      BOUND(1,2) = -0.1
      ELSE
      BOUND(1,2) = BOUND(1,2) * RY_IN_EV / TOT_ATOM
      ENDIF
      IF (BOUND(1,1).EQ.-9999) THEN 
      BOUND(1,1) = A - 2.5
      ELSE
      BOUND(1,1) = BOUND(1,1) * RY_IN_EV / TOT_ATOM
      ENDIF
      
      READ (*, '("Equilibrium Cell Volume (au^3):")')
      READ *, BOUND(2,2),GUESS(2),BOUND(2,1)
      IF (GUESS(2).EQ.-9999) THEN 
      GUESS(2) = B
      ELSE
      GUESS(2) = GUESS(2) * AU_IN_A**3 / TOT_ATOM
      ENDIF
      IF (BOUND(2,2).EQ.-9999) THEN 
      BOUND(2,2) = B + 3.
      ELSE
      BOUND(2,2) = BOUND(2,2) * AU_IN_A**3 / TOT_ATOM
      ENDIF
      IF (BOUND(2,1).EQ.-9999) THEN 
      BOUND(2,1) = MAX(B-3.,2.)
      ELSE
      BOUND(2,1) = BOUND(2,1) * AU_IN_A**3 / TOT_ATOM
      ENDIF
      
      READ (*, '("Bulk Modulus (MBar):")')
      READ *, BOUND(3,2),GUESS(3),BOUND(3,1)
      IF (GUESS(3).EQ.-9999) THEN
      GUESS(3) = 2. / STRESS_IN_MBAR
      ELSE
      GUESS(3) = GUESS(3) / STRESS_IN_MBAR
      ENDIF
      IF (BOUND(3,2).EQ.-9999) THEN 
      BOUND(3,2) = 7. / STRESS_IN_MBAR
      ELSE
      BOUND(3,2) = ABS(BOUND(3,2)) / STRESS_IN_MBAR
      ENDIF
      IF (BOUND(3,1).EQ.-9999) THEN 
      BOUND(3,1) = 0.1 / STRESS_IN_MBAR
      ELSE
      BOUND(3,1) = BOUND(3,1) / STRESS_IN_MBAR
      ENDIF
      
C     Nonlinearities in the Taylor expansion:
C     Should not have too many terms either.
      DO I = 4,NOP
      GUESS(I) = 0.
      BOUND(I,2) = 100.
      BOUND(I,1) = -100.
      ENDDO
      
      READ (*,'(/,"DEL: ",F)') DEL
      READ (*,'("NDIV: ",I)') NDIV
      READ (*,'("Save composite energies on files (*.fit,*.m): ",
     A           A20)') NAME
      READ (*,'("Number of mesh points: ",I)') MESH
      IF (MESH.EQ.-9999) MESH=100
      READ (*,'("Minimum volume/atom (A^3): ",F)') VOL_MIN
      IF (VOL_MIN.EQ.-9999) VOL_MIN=C*0.98
      READ (*,'("Maximum volume/atom (A^3): ",F)') VOL_MAX
      IF (VOL_MAX.EQ.-9999) VOL_MAX=D*1.02
      
      WRITE (*,'(/,"Optimization starts:",/)')
      CALL MINA (FE, NOP, NDIV, DEL, BOUND, GUESS, GUESS, EMIN, IERR)
      
      CC = FE(GUESS)
      WRITE (*,'(/,"The average error for parameter set")')
      WRITE (*,'(9F12.6)') (GUESS(I),I=1,NOP)
      WRITE (*,'("is", F12.6," eV.")') CC

      WRITE (*,'(/,"Cohesive energy minimum:",F12.6," eV/atom,")')
     A   GUESS(1)
      WRITE (*,'("Equilibrium volume:",F12.6," A^3/atom,")') GUESS(2)
      WRITE (*,'("Bulk modulus:",F12.6," MBar.",/)')
     A   GUESS(3)*STRESS_IN_MBAR
      
      PRINT *,'Saving cohesive energy in ',
     A NAME(1:INDEX(NAME,' ')-1)//'.fit'//' ...'
      OPEN (UNIT=LP, FILE=NAME(1:INDEX(NAME,' ')-1)//'.fit',
     A      STATUS='UNKNOWN',FORM='FORMATTED')
      WRITE (LP,'(3I11)') 2+MESH,2+MESH,2+MESH
      VOLUME_FIT(1) = GUESS(2)
      DO I = 0,MESH
      VOLUME_FIT(2+I) = VOL_MIN + I*(VOL_MAX-VOL_MIN)/MESH
      ENDDO
      DO I = 1,MESH+2
      IF (IFT.LT.0) THEN 
      ENERGY_FIT(I) = UNIV_COHESION(GUESS,VOLUME_FIT(I))
      ELSE
      ENERGY_FIT(I) = TAYLOR_COHESION(IFT,GUESS,VOLUME_FIT(I))
      ENDIF
      WRITE (LP,'(3F12.6)') VOLUME_FIT(I), ENERGY_FIT(I), 
     A                      ENERGY_FIT(I)-SHIFT
      ENDDO
      WRITE (LP,'(A45)') LABEL_LDA
      WRITE (LP,'("Bulk modulus:",F9.6," MBar, error:",F9.6," eV.")')
     A   GUESS(3)*STRESS_IN_MBAR, EMIN
      CLOSE (LP)
      
C     Matlab script to plot out the data:
      PRINT *,'Writing Matlab script in ',
     A  NAME(1:INDEX(NAME,' ')-1)//'.m'//' ...'
      PRINT *,' '
      OPEN (UNIT=LP,FILE=NAME(1:INDEX(NAME,' ')-1)//'.m',
     A      STATUS='UNKNOWN',FORM='FORMATTED')
      WRITE (LP,'("clear; fit=[")')
      WRITE (LP,'(2F12.6)') (VOLUME_fit(I),ENERGY_FIT(I),I=2,MESH+2)
      WRITE (LP,'("];clf;plot(fit(:,1),fit(:,2));hold on;")')
      WRITE (LP,'("raw=[")')
      WRITE (LP,'(2F12.6)') (VOLUME_LDA(I),ENERGY_LDA(I),I=1,NUM_LDA)
      WRITE (LP,'("];plot(raw(:,1),raw(:,2),''ro'');")')
      WRITE (LP,'("title(''",A45,"(",3F7.3,") error=",F6.5," eV.'');")')
     A LABEL_LDA,GUESS(1),GUESS(2),GUESS(3)*STRESS_IN_MBAR,EMIN
      WRITE (LP,'("xlabel(''Volume/atom [A^3]'');")')
      WRITE (LP,'("ylabel(''Composite cohesive energy/atom [eV]'');")')
      CLOSE (LP)
      
      STOP
      END
      
      
      DOUBLE PRECISION FUNCTION FE (GUESS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MMS=1000,EXPONENT=2.D0)
      COMMON /A/ VOLUME_LDA(MMS),ENERGY_LDA(MMS),NUM_LDA,IFT,NOP
      DIMENSION GUESS(*)
      FE = 0.D0
      DO I = 1,NUM_LDA
      IF (IFT.LT.0) THEN 
      A = UNIV_COHESION(GUESS,VOLUME_LDA(I))
      ELSE
      A = TAYLOR_COHESION(IFT,GUESS,VOLUME_LDA(I))
      ENDIF
      FE = FE + (A-ENERGY_LDA(I))**2 / NUM_LDA
      ENDDO
      FE = SQRT(FE)
      WRITE (*,'("fe=", F8.4," eV, ", 5F11.6)')
     A  FE,(GUESS(I),I=1,NOP)
      RETURN
      END
      
      
      DOUBLE PRECISION FUNCTION TAYLOR_COHESION (N,PARAM,VOLUME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PARAM(*)
      TAYLOR_COHESION = PARAM(1)
      IFAC = 1
      DO I = 2,N
      IFAC = IFAC * I
      TAYLOR_COHESION = TAYLOR_COHESION + PARAM(I+1) *
     A (VOLUME-PARAM(2))**I / PARAM(2)**(I-1) / IFAC
      ENDDO
      RETURN
      END
      
      
      DOUBLE PRECISION FUNCTION UNIV_COHESION (PARAM,VOLUME)
C     -------------------------------------------------------
C     GIVEN COHESION WELL-DEPTH, EQUILIBRIUM VOLUME AND BULK 
C     MODULUS, COMPUTE THE COHESIVE ENERGY AT ANY VOLUME.
C     -------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (PI=3.14159265358979323846D0)
      DIMENSION PARAM(*)
C     UNIVERSAL BINDING CURVE FOR METALS AFTER SCALING:
C     WIGNER-SEITZ RADIUS
      R  = (3.*VOLUME/4./PI)**(1./3.)
      R0 = (3.*PARAM(2)/4./PI)**(1./3.)
C     INTRINSIC LENGTHSCALE GIVEN BY BULK MODULUS
      RL = SQRT(ABS(PARAM(1))/12./PI/R0/PARAM(3))
C     DIMENSIONLESS LENGTH:
      X = (R-R0) / RL
      UNIV_COHESION = PARAM(1)*(1+X+0.05*X**3)*EXP(-X)
      RETURN
      END
      
      
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
C***BEGIN PROLOGUE  DGESL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D2A1
C***KEYWORDS  DOUBLE PRECISION,LINEAR ALGEBRA,LINPACK,MATRIX,SOLVE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Solves the double precision system  A*X=B or  TRANS(A)*X=B
C            using the factors computed by DGECO or DGEFA.
C***DESCRIPTION
C
C     DGESL solves the double precision system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGECO or DGEFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DGECO or DGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGECO or DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B  where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGECO has set RCOND .GT. 0.0
C        or DGEFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DDOT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DDOT
C***END PROLOGUE  DGESL
      INTEGER LDA,N,IPVT(1),JOB
      DOUBLE PRECISION A(LDA,1),B(1)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
C***FIRST EXECUTABLE STATEMENT  DGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C***BEGIN PROLOGUE  DAXPY
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D1A7
C***KEYWORDS  BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,TRIAD,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P computation y = a*x + y
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY), where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N
C       and LY is defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DAXPY
C
      DOUBLE PRECISION DX(1),DY(1),DA
C***FIRST EXECUTABLE STATEMENT  DAXPY
      IF(N.LE.0.OR.DA.EQ.0.D0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DA*DX(I) + DY(I)
   70     CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C***BEGIN PROLOGUE  DDOT
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D1A4
C***KEYWORDS  BLAS,DOUBLE PRECISION,INNER PRODUCT,LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P. inner product of d.p. vectors
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C     DDOT  double precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of double precision DX and DY.
C     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY)
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DDOT
C
      DOUBLE PRECISION DX(1),DY(1)
C***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     1   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
      RETURN
C
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DDOT = DDOT + DX(I)*DY(I)
   70     CONTINUE
      RETURN
      END
      

      SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)
C***BEGIN PROLOGUE  DGECO
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D2A1
C***KEYWORDS  CONDITION,DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,
C             MATRIX
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a double precision matrix by Gaussian elimination
C            and estimates the condition of the matrix.
C***DESCRIPTION
C
C     DGECO factors a double precision matrix by Gaussian elimination
C     and estimates the condition of the matrix.
C
C     If  RCOND  is not needed, DGEFA is slightly faster.
C     To solve  A*X = B , follow DGECO by DGESL.
C     To compute  INVERSE(A)*C , follow DGECO by DGESL.
C     To compute  DETERMINANT(A) , follow DGECO by DGEDI.
C     To compute  INVERSE(A) , follow DGECO by DGEDI.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an INTEGER vector of pivot indices.
C
C        RCOND   DOUBLE PRECISION
C                an estimate of the reciprocal condition of  A .
C                For the system  A*X = B , relative perturbations
C                in  A  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  A  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C
C        Z       DOUBLE PRECISION(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     LINPACK DGEFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     Fortran DABS,DMAX1,DSIGN
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DASUM,DAXPY,DDOT,DGEFA,DSCAL
C***END PROLOGUE  DGECO
      INTEGER LDA,N,IPVT(1)
      DOUBLE PRECISION A(LDA,1),Z(1)
      DOUBLE PRECISION RCOND
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
C
C     COMPUTE 1-NORM OF A
C
C***FIRST EXECUTABLE STATEMENT  DGECO
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL DGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30
            S = DABS(A(K,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + DABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150
            S = DABS(A(K,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
C***BEGIN PROLOGUE  DGEFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D2A1
C***KEYWORDS  DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a double precision matrix by Gaussian elimination.
C***DESCRIPTION
C
C     DGEFA factors a double precision matrix by Gaussian elimination.
C
C     DGEFA is usually called by DGECO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGESL or DGEDI will divide by zero
C                     if called.  Use  RCOND  in DGECO for a reliable
C                     indication of singularity.
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DSCAL,IDAMAX
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DSCAL,IDAMAX
C***END PROLOGUE  DGEFA
      INTEGER LDA,N,IPVT(1),INFO
      DOUBLE PRECISION A(LDA,1)
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
C***FIRST EXECUTABLE STATEMENT  DGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
      SUBROUTINE DSCAL(N,DA,DX,INCX)
C***BEGIN PROLOGUE  DSCAL
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D1A6
C***KEYWORDS  BLAS,LINEAR ALGEBRA,SCALE,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P. vector scale x = a*x
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scale factor
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  double precision result (unchanged if N.LE.0)
C
C     Replace double precision DX by double precision DA*DX.
C     For I = 0 to N-1, replace DX(1+I*INCX) with  DA * DX(1+I*INCX)
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DSCAL
C
      DOUBLE PRECISION DA,DX(1)
C***FIRST EXECUTABLE STATEMENT  DSCAL
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          DX(I) = DA*DX(I)
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
C***BEGIN PROLOGUE  IDAMAX
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D1A2
C***KEYWORDS  BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,MAXIMUM COMPONENT,
C             VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Find largest component of d.p. vector
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C   IDAMAX  smallest index (zero if N .LE. 0)
C
C     Find smallest index of maximum magnitude of double precision DX.
C     IDAMAX =  first I, I = 1 to N, to minimize  ABS(DX(1-INCX+I*INCX)
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  IDAMAX
C
      DOUBLE PRECISION DX(1),DMAX,XMAG
C***FIRST EXECUTABLE STATEMENT  IDAMAX
      IDAMAX = 0
      IF(N.LE.0) RETURN
      IDAMAX = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      DMAX = DABS(DX(1))
      NS = N*INCX
      II = 1
          DO 10 I = 1,NS,INCX
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 5
          IDAMAX = II
          DMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 30
          IDAMAX = I
          DMAX = XMAG
   30 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
C***BEGIN PROLOGUE  DASUM
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D1A3A
C***KEYWORDS  ADD,BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,MAGNITUDE,SUM,
C             VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Sum of magnitudes of d.p. vector components
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    DASUM  double precision result (zero if N .LE. 0)
C
C     Returns sum of magnitudes of double precision DX.
C     DASUM = sum from 0 to N-1 of DABS(DX(1+I*INCX))
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DASUM
C
      DOUBLE PRECISION DX(1)
C***FIRST EXECUTABLE STATEMENT  DASUM
      DASUM = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I=1,NS,INCX
          DASUM = DASUM + DABS(DX(I))
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DASUM = DASUM + DABS(DX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
         DASUM = DASUM + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2))
     1   + DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
   50 CONTINUE
      RETURN
      END

      
      
      SUBROUTINE MINA(FN,NV,NDIV,DEL,A,GUESS,X,FOFX,IERR)               MINA.2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MOP=10)
c                                                                       ADDRESS.2
c     SANDIA MATHEMATICAL PROGRAM LIBRARY                               ADDRESS.3
c     APPLIED MATHEMATICS DIVISION 2613                                 ADDRESS.
c     SANDIA LABORATORIES                                               ADDRESS.5
c     ALBUQUERQUE, NEW MEXICO  87185                                    JUN0278.1
c     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978                     JUN0278.2
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.1
c                    ISSUED BY SANDIA LABORATORIES                     *MAR1378.2
c  *                   A PRIME CONTRACTOR TO THE                       *MAR1378.3
c  *                UNITED STATES DEPARTMENT OF ENERGY                 *MAR1378.4
c  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * * *MAR1378.5
c  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE   *MAR1378.6
c  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE      *MAR1378.7
c  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,    *MAR1378.8
c  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES  *MAR1378.9
c  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL      *MAR1378.10
c  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR     *MAR1378.11
c  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS      *MAR1378.12
c  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE          *MAR1378.13
c  * OWNED RIGHTS.                                                     *MAR1378.14
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.15
c  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS     *MAR1378.16
c  * PART IS SAND77-1441.                                              *MAR1378.17
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.18
c                                                                       ADDRESS.9
c     ORIGINAL ROUTINE WAS H2 SAND MIN, BY Z. BEISINGER AND S. BELL     MINA.4
c     PRESENT VERSION BY R E JONES                                      MINA.5
c                                                                       MINA.6
c     ABSTRACT                                                          MINA.7
c        MINA FINDS AN APPROXIMATE MINIMUM OF A REAL FUNCTION OF        MINA.8
c        NV VARIABLES, GIVEN AN INITIAL ESTIMATE OF THE POSITION OF     MINA.9
c        THE MINIMUM AND RANGES FOR EACH OF THE VARIABLES.              MINA.10
c        MINA USES A SELECTIVE DIRECTED SEARCH OF A SURROUNDING         MINA.11
c        NV-DIMENSIONAL GRID OF POINTS TO FIND A DIRECTION IN WHICH     MINA.12
c        THE FUNCTION DECREASES.  IT THEN PROCEEDS IN THIS DIRECTION    MINA.13
c        AS FAR AS THE FUNCTION DECREASES, THEN DETERMINES A NEW        MINA.14
c        DIRECTION TO TRAVEL.  WHEN NO SUCH DIRECTION IS FOUND THE      MINA.15
c        SEARCH INCREMENT FACTOR IS DECREASED AND THE PROCESS           MINA.16
c        IS REPEATED.                                                   MINA.17
c                                                                       MINA.18
c     DESCRIPTION OF ARGUMENTS                                          MINA.19
c        THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST  MINA.20
c              A(NV,2), GUESS(NV), X(NV)                                MINA.21
c                                                                       MINA.22
c        INPUT--                                                        MINA.23
c        FN   - NAME OF FUNCTION OF NV VARIABLES TO BE MINIMIZED.       MINA.24
c               (THIS NAME MUST APPEAR IN AN EXTERNAL STATEMENT.)       MINA.25
c               FORM OF THE CALLING SEQUENCE MUST BE FUNCTION FN(X),    MINA.26
c               WHERE X IS AN ARRAY OF NV VARIABLE VALUES. THE          MINA.27
c               ORDERING OF THE VARIABLES IS ARBITRARY, EXCEPT          MINA.28
c               THAT IT MUST AGREE WITH THE ORDERING USED IN            MINA.29
c               ARRAYS A AND GUESS.                                     MINA.30
c        NV   - NUMBER OF VARIABLES.  (NV .GE. 1)                       MINA.31
c        NDIV - NUMBER OF REFINEMENTS OF THE SEARCH INCREMENTS TO USE.  MINA.32
c               AT EACH REFINEMENT, THE INCREMENT IN EACH DIMENSION     MINA.33
c               IS DIVIDED BY 10.  (USUALLY NDIV IS ABOUT 3 OR 4.)      MINA.34
c        DEL  - FRACTION OF VARIABLE RANGE (IN EACH DIMENSION) TO USE   MINA.35
c               AS THE INITIAL INCREMENT (IN THAT DIMENSION)            MINA.36
c        A    - ARRAY OF SEARCH BOUNDS, DIMENSIONED NV BY 2.            MINA.37
c               A(I,1) SHOULD BE THE LOWER BOUND OF THE I-TH VARIABLE.  MINA.38
c               A(I,2) SHOULD BE THE UPPER BOUND OF THE I-TH VARIABLE.  MINA.39
c        GUESS- ARRAY OF NV INITIAL VALUES.  GUESS(I) SHOULD BE THE     MINA.40
c               INITIAL VALUE TO USE FOR THE I-TH VARIABLE.             MINA.41
c                                                                       MINA.42
c        OUTPUT--                                                       MINA.43
c        X    - ARRAY (DIMENSIONED NV) GIVING THE VALUES OF THE         MINA.44
c               VARIABLES AT THE MINIMUM.  X(I) WILL BE THE VALUE       MINA.45
c               OF THE I-TH VARIABLE.                                   MINA.46
c        FOFX - FUNCTION VALUE AT THE MINIMUM                           MINA.47
c        IERR - A STATUS CODE                                           MINA.48
c              -NORMAL CODE                                             MINA.49
c               =1 MEANS THE SEARCH FOR A MINIMUM PROCEEDED FOR THE     MINA.50
c                  SPECIFIED NUMBER OF REFINEMENTS.                     MINA.51
c              -ABNORMAL CODES                                          MINA.52
c               =2 MEANS NV IS GREATER THAN 50                          MINA.53
c               =3 MEANS A RANGE MINIMUM IS GREATER THAN THE            MINA.54
c                  CORRESPONDING MAXIMUM                                MINA.55
c                                                                       MINA.57
      DIMENSION A(MOP,2),GUESS(MOP),X(MOP)                              MINA.58
      DIMENSION XNOW(MOP),XNEW(MOP),R(MOP)                              MINA.59
      EXTERNAL FN
c                                                                       MESS.2
c     CALL MLMESS(10HSANDIA7.2 )                                        MESS3.1
c                                                                       MESS.4
c                                                                       MINA.60
c     INITIALIZE                                                        MINA.61
c                                                                       MINA.62
      IERR = 1                                                          MINA.64
      IF (NV.LE.MOP) GO TO 2                                             MINA.65
cc      CALL ERRCHK(33,33HIN MINA  , NV IS GREATER THAN 100.)            MINA.66
      write(6,*)'nv is greater than 100'
      IERR = 2                                                          MINA.67
      RETURN                                                            MINA.68
    2 NX = NV                                                           MINA.69
      IDIV = 0                                                          MINA.70
      DO 5 I=1,NX                                                       MINA.71
      XNOW(I) = GUESS(I)                                                MINA.72
      IF (XNOW(I).LT.A(I,1)) XNOW(I) = A(I,1)                           MINA.73
      IF (XNOW(I).GT.A(I,2)) XNOW(I) = A(I,2)                           MINA.74
      IF (A(I,1)-A(I,2)) 5,5,4                                          MINA.75
    4  write(6,*) 'range minimum greater than maximum' 
       print *,A(i,1),A(i,2)
c    4 CALL ERRCHK(46,46HIN MINA  , RANGE MINIMUM GREATER THAN MAXIMUM.) MINA.76
      IERR = 3                                                          MINA.77
      RETURN                                                            MINA.78
    5 R(I) = A(I,2)-A(I,1)                                              MINA.79
      DELTA = DEL                                                       MINA.80
      IF (DELTA.LE.0.0) DELTA = 0.1                                     MINA.81
      FNOW = FN(XNOW)                                                   MINA.82
c                                                                       MINA.83
c     FIND NEW DIRECTION                                                MINA.84
c                                                                       MINA.85
    7 DO 8 I=1,NX                                                       MINA.86
    8 XNEW(I) = XNOW(I)                                                 MINA.87
      FOLD = FNOW                                                       MINA.88
   10 DO 40 I=1,NX                                                      MINA.89
      IF (XNOW(I).GE.A(I,2)) GO TO 20                                   MINA.90
      XNEW(I) = MIN(XNOW(I)+DELTA*R(I),A(I,2))
      FNEW = FN(XNEW)                                                   MINA.92
      IF (FNEW.LT.FNOW) GO TO 30                                        MINA.93
   20 IF (XNOW(I) .LE. A(I,1)) GO TO 25                                 MINA.94
      XNEW(I) = MAX(XNOW(I)-DELTA*R(I),A(I,1))
      FNEW = FN(XNEW)                                                   MINA.96
      IF (FNEW.LT.FNOW) GO TO 30                                        MINA.97
   25 XNEW(I) = XNOW(I)                                                 MINA.98
      GO TO 40                                                          MINA.99
   30 FNOW = FNEW                                                       MINA.100
   40 CONTINUE                                                          MINA.101
      ISTEP = 1                                                         MINA.102
c                                                                       MINA.103
c     REFINE IF NEEDED                                                  MINA.104
c                                                                       MINA.105
      IF (FNOW.LT.FOLD) GO TO 50                                        MINA.106
      IF (IDIV.GE.NDIV) GO TO 100                                       MINA.107
      DELTA = DELTA*0.1                                                 MINA.108
      IDIV = IDIV+1                                                     MINA.109
      GO TO 10                                                          MINA.110
c                                                                       MINA.111
c     TRY TO CONTINUE IN CHOSEN DIRECTION                               MINA.112
c                                                                       MINA.113
   50 ICHNG = 0                                                         MINA.114
      FAc = 1.0                                                         MINA.115
      IF ((ISTEP/10)*10.EQ.ISTEP) FAc = 2.0                             MINA.116
      DO 60 I=1,NX                                                      MINA.117
      DX = (XNEW(I)-XNOW(I))*FAc                                        MINA.118
      XNOW(I) = XNEW(I)                                                 MINA.119
      IF (DX) 52,54,56                                                  MINA.120
   52 XNEW(I) = MAX(XNOW(I)+DX,A(I,1))
      IF (XNEW(I).LT.XNOW(I)) ICHNG = 1                                 MINA.122
      GO TO 60                                                          MINA.123
   54 XNEW(I) = XNOW(I)                                                 MINA.124
      GO TO 60                                                          MINA.125
   56 XNEW(I) = MIN(XNOW(I)+DX,A(I,2))
      IF (XNEW(I).GT.XNOW(I)) ICHNG = 1                                 MINA.127
   60 CONTINUE                                                          MINA.128
      IF (ICHNG.EQ.0) GO TO 7                                           MINA.129
      FNEW = FN(XNEW)                                                   MINA.130
      IF (FNEW.GE.FNOW) GO TO 7                                         MINA.131
      FNOW = FNEW                                                       MINA.132
      ISTEP = ISTEP+1                                                   MINA.133
      GO TO 50                                                          MINA.134
c                                                                       MINA.135
c     RETURN ANSWERS                                                    MINA.136
c                                                                       MINA.137
  100 FOFX = FOLD                                                       MINA.138
      DO 110 I=1,NX                                                     MINA.139
  110 X(I) = XNOW(I)                                                    MINA.140
      RETURN                                                            MINA.141
      END                                                               MINA.142
      SUBROUTINE SIMIN (F,K,EPS,ANS,S,NEV,ICONT,Y)                      SIMIN.2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c                                                                       ADDRESS.2
c     SANDIA MATHEMATICAL PROGRAM LIBRARY                               ADDRESS.3
c     APPLIED MATHEMATICS DIVISION 2613                                 ADDRESS.4
c     SANDIA LABORATORIES                                               ADDRESS.5
c     ALBUQUERQUE, NEW MEXICO  87185                                    JUN0278.1
c     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978                     JUN0278.2
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.1
c                    ISSUED BY SANDIA LABORATORIES                     *MAR1378.2
c  *                   A PRIME CONTRACTOR TO THE                       *MAR1378.3
c  *                UNITED STATES DEPARTMENT OF ENERGY                 *MAR1378.4
c  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * * *MAR1378.5
c  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE   *MAR1378.6
c  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE      *MAR1378.7
c  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,    *MAR1378.8
c  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES  *MAR1378.9
c  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL      *MAR1378.10
c  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR     *MAR1378.11
c  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS      *MAR1378.12
c  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE          *MAR1378.13
c  * OWNED RIGHTS.                                                     *MAR1378.14
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.15
c  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS     *MAR1378.16
c  * PART IS SAND77-1441.                                              *MAR1378.17
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *MAR1378.18
c                                                                       ADDRESS.9
c     ORIGINAL ROUTINE BY L F SHAMPINE, AS DESCRIBED IN REF.1 BELOW.    SIMIN.4
c     PREPARATION FOR MATH LIBRARY BY R E JONES.                        SIMIN.5
c                                                                       SIMIN.6
c     ABSTRACT                                                          SIMIN.7
c        SIMIN FINDS AN APPROXIMATE MINIMUM OF A REAL FUNCTION OF K     SIMIN.8
c        VARIABLES, GIVEN AN INITIAL ESTIMATE OF THE POSITION OF THE    SIMIN.9
c        MINIMUM. THE SIMPLEX METHOD IS USED. SEE REFERENCE 1 BELOW     SIMIN.10
c        FOR A FULL EXPLANATION OF THIS METHOD.  BRIEFLY, A SET OF      SIMIN.11
c        K+1 POINTS IN K-DIMENSIONAL SPACE IS CALLED A SIMPLEX.         SIMIN.12
c        THE MINIMIZATION PROCESS ITERATES BY REPLACING THE POINT       SIMIN.13
c        WITH THE LARGEST FUNCTION VALUE BY A NEW POINT WITH A          SIMIN.14
c        SMALLER FUNCTION VALUE.  ITERATION CONTINUES UNTIL ALL THE     SIMIN.15
c        POINTS CLUSTER SUFFICIENTLY CLOSE TO A MINIMUM.                SIMIN.16
c                                                                       SIMIN.17
c     REFERENCES                                                        SIMIN.18
c      1. L F SHAMPINE, A ROUTINE FOR UNCONSTRAINED OPTIMIZATION,       SIMIN.19
c         SC-TM-72130  OR  SC-RR-720657                                 SIMIN.20
c      2. J A NELDER AND R MEAD, A SIMPLEX METHOD FOR FUNCTION          SIMIN.21
c         MINIMIZATION, COMPUTER JOURNAL, 7(1965) 308-313               SIMIN.22
c                                                                       SIMIN.23
c     DESCRIPTION OF PARAMETERS                                         SIMIN.24
c      --INPUT--                                                        SIMIN.25
c        F  - NAME OF FUNCTION OF K VARIABLES TO BE MINIMIZED.          SIMIN.26
c             (THIS NAME MUST APPEAR IN AN EXTERNAL STATEMENT.)         SIMIN.27
c             FORM OF THE CALLING SEQUENCE MUST BE FUNCTION F(X),       SIMIN.28
c             WHERE X IS AN ARRAY OF K VARIABLES.                       SIMIN.29
c        K  - THE NUMBER OF VARIABLES.  K MUST BE AT LEAST 2.           SIMIN.30
c             NORMALLY K SHOULD BE LESS THAN ABOUT 10, AS SIMIN         SIMIN.31
c             BECOMES LESS EFFECTIVE FOR LARGER VALUES OF K.            SIMIN.32
c        EPS- THE CONVERGENCE CRITERION.  LET YAVG BE THE AVERAGE       SIMIN.33
c             VALUE OF THE FUNCTION F AT THE K+1 POINTS OF THE          SIMIN.34
c             SIMPLEX, AND LET R BE THEIR STANDARD ERROR.  (THAT IS,    SIMIN.35
c             THE ROOT-MEAN-SQUARE OF THE SET OF VALUES (Y(I)-YAVG),    SIMIN.36
c             WHERE Y(I) IS THE FUNCTION VALUE AT THE I-TH POINT OF     SIMIN.37
c             THE SIMPLEX.)  THEN--                                     SIMIN.38
c             IF EPS.GT.0, CONVERGENCE IS OBTAINED IF  R.LE.EPS.        SIMIN.39
c             IF EPS.LT.0, CONVERGENCE IS IF  R.LE.ABS(EPS*YAVG).       SIMIN.40
c             IF EPS=0, THE PROCESS WILL NOT CONVERGE BUT INSTEAD WILL  SIMIN.41
c             QUIT WHEN NEV FUNCTION EVALUATIONS HAVE BEEN USED.        SIMIN.42
c        ANS- AN ARRAY OF LENGTH K CONTAINING A GUESS FOR THE LOCATION  SIMIN.43
c             OF A MINIMUM OF F.                                        SIMIN.44
c        S  - A SCALE PARAMETER, WHICH MAY BE A SIMPLE VARIABLE OR AN   SIMIN.45
c             ARRAY OF LENGTH K.  USE OF AN ARRAY IS SIGNALLED BY       SIMIN.46
c             SETTING S(1) NEGATIVE.                                    SIMIN.47
c             -SIMPLE VARIABLE CASE.  HERE S IS THE LENGTH OF EACH      SIMIN.48
c             SIDE OF THE INITIAL SIMPLEX.  THUS, THE INITIAL SEARCH    SIMIN.49
c             RANGE IS THE SAME FOR ALL THE VARIABLES.                  SIMIN.50
c             -ARRAY CASE.  HERE THE LENGTH OF SIDE I OF THE INITIAL    SIMIN.51
c             SIMPLEX IS ABS(S(I)).  THUS, THE INITIAL SEARCH RANGE     SIMIN.52
c             MAY BE DIFFERENT FOR DIFFERENT VARIABLES.                 SIMIN.53
c             NOTE-- THE VALUE(S) USED FOR S ARE NOT VERY CRITICAL.     SIMIN.54
c             ANY REASONABLE GUESS SHOULD DO O.K.                       SIMIN.55
c        NEV- THE MAXIMUM NUMBER OF FUNCTION EVALUATIONS TO BE USED.    SIMIN.56
c             (THE ACTUAL NUMBER USED MAY EXCEED THIS SLIGHTLY SO THE   SIMIN.57
c             LAST SEARCH ITERATION MAY BE COMPLETED.)                  SIMIN.58
c        ICONT - ICONT SHOULD BE ZERO ON ANY CALL TO SIMIN WHICH        SIMIN.59
c             IS NOT A CONTINUATION OF A PREVIOUS CALL.                 SIMIN.60
c             IF ICONT=1 THE PROBLEM WILL BE CONTINUED.  IN THIS        SIMIN.61
c             CASE THE WORK ARRAY Y MUST BE THE SAME ARRAY THAT WAS     SIMIN.62
c             USED IN THE CALL THAT IS BEING CONTINUED (AND THE VALUES  SIMIN.63
c             IN IT MUST BE UNCHANGED).  THE REASON FOR THIS IS THAT    SIMIN.64
c             IF ICONT=1 THEN THE ARGUMENT S IS IGNORED AND THE SIMPLEX SIMIN.65
c             AND RELATED FUNCTION VALUES THAT WERE STORED IN ARRAY Y   SIMIN.66
c             DURING A PREVIOUS EXECUTION ARE USED TO CONTINUE THAT     SIMIN.67
c             PREVIOUS PROBLEM.                                         SIMIN.68
c        Y  - A WORK ARRAY CONTAINING AT LEAST K*K + 5*K + 1 WORDS.     SIMIN.69
c             IF ICONT=1 THIS MUST BE THE SAME ARRAY USED IN THE CALL   SIMIN.70
c             THAT IS BEING CONTINUED.                                  SIMIN.71
c      --OUTPUT--                                                       SIMIN.72
c        ANS- ANS WILL CONTAIN THE LOCATION OF THE POINT WITH THE       SIMIN.73
c             SMALLEST VALUE OF THE FUNCTION THAT WAS FOUND.            SIMIN.74
c        S  - IN THE SIMPLE VARIABLE CASE S WILL BE RETURNED AS THE     SIMIN.75
c             AVERAGE DISTANCE FROM THE VERTICES TO THE CENTROID OF     SIMIN.76
c             THE SIMPLEX.                                              SIMIN.77
c             IN THE ARRAY CASE S(I) WILL BE RETURNED AS THE AVERAGE    SIMIN.78
c             DISTANCE IN THE I-TH DIMENSION OF VERTICES FROM           SIMIN.79
c             THE CENTROID.  (S(1) WILL BE NEGATED.)                    SIMIN.80
c             NOTE-- THE VALUE(S) RETURNED IN S ARE USEFUL FOR          SIMIN.81
c             ASSESSING THE FLATNESS OF THE FUNCTION NEAR THE           SIMIN.82
c             MINIMUM.  THE LARGER THE VALUE OF S (FOR A GIVEN          SIMIN.83
c             VALUE OF EPS), THE FLATTER THE FUNCTION.                  SIMIN.84
c        NEV- NEV WILL BE THE COUNT OF THE ACTUAL NUMBER OF FUNCTION    SIMIN.85
c             EVALUATIONS USED.                                         SIMIN.86
c        Y  - WILL CONTAIN ALL DATA NEEDED TO CONTINUE THE MINIMIZATION SIMIN.87
c             SEARCH EFFICIENTLY IN A SUBSEQUENT CALL.                  SIMIN.88
c             NOTE -- THE FIRST K+1 ELEMENTS OF Y WILL CONTAIN THE      SIMIN.89
c             FUNCTION VALUES AT THE K+1 POINTS OF THE LATEST SIMPLEX.  SIMIN.90
c             THE NEXT K*(K+1) ELEMENTS OF Y WILL BE THE K+1 POINTS     SIMIN.91
c             OF THE SIMPLEX (IN EXACT CORRESPONDENSE TO THE ARRAY      SIMIN.92
c             P DISCUSSED IN REFERENCE 1 ABOVE).  THE REMAINING 3*K     SIMIN.93
c             WORDS ARE TEMPORARY WORKING STORAGE ONLY.                 SIMIN.94
c                                                                       SIMIN.96
      DIMENSION ANS(K),S(K),Y(1)                                        SIMIN.97
      EXTERNAL F                                                        SIMIN.98
c                                                                       MESS.2
c     CALL MLMESS(10HSANDIA7.2 )                                        MESS3.1
c                                                                       MESS.4
c                                                                       SIMIN.99
      IF(K.GE.2 .AND. S(1).NE.0.) GO TO 10                              SIMIN.101
      write(6,*) 's(1)=0 or k is less than 2'
c      CALL ERRCHK(39,39HIN SIMIN , S(1)=0. OR K IS LESS THAN 2.)        SIMIN.102
      RETURN                                                            SIMIN.103
c   10 IF (K.GT.100) CALL ERRCHK(31,31HIN SIMIN , K IS LARGER THAN 100.)   SIMIN.104
   10 write(6,*) 'k is large than 100'
c                                                                       SIMIN.105
      IP = K+2                                                          SIMIN.106
      Ic = IP+K*(K+1)                                                   SIMIN.107
      IR = IC+K                                                         SIMIN.108
      IRR = IR+K                                                        SIMIN.109
      CALL SIMINA(F,K,EPS,ANS,S,NEV,ICONT,Y,Y(IP),Y(IC),Y(IR),Y(IRR))   SIMIN.110
      RETURN                                                            SIMIN.111
      END                                                               SIMIN.112
      SUBROUTINE SIMINA(F,K,EPS,ANS,S,NEV,ICONT,Y,P,PC,PR,PRR)          SIMIN.113
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ANS(K),S(K),Y(3),P(K,3),PC(K),PR(K),PRR(K)              SIMIN.114
      DATA ALPHA,BETA,GAMMA/1.0,0.5,2.0/                                SIMIN.115
      EXTERNAL F
c                                                                       SIMIN.116
c     SIMINA IS A CORE MINIMIZATION ROUTINE CALLED ONLY BY SIMIN.       SIMIN.117
c                                                                       SIMIN.118
c     INITIALIZATION                                                    SIMIN.119
c                                                                       SIMIN.120
      KOUNT = 0                                                         SIMIN.122
      KK = K                                                            SIMIN.123
      IF (KK.LE.1) GO TO 99                                             SIMIN.124
      ONEK = 1.0/FLOAT(KK)                                              SIMIN.125
      KP1 = KK+1                                                        SIMIN.126
      ONEKP1 = 1.0/FLOAT(KP1)                                           SIMIN.127
      TOL = FLOAT(KP1)*EPS**2                                           SIMIN.128
c                                                                       SIMIN.129
c     INITIAL SIMPLEX                                                   SIMIN.130
c                                                                       SIMIN.131
      IF (ICONT.GE.1) GO TO 10                                          SIMIN.132
      IF (S(1)) 4,99,1                                                  SIMIN.133
    1 SKP1 = S(1)*ONEKP1                                                SIMIN.134
      DO 2 I=1,KP1                                                      SIMIN.135
      DO 2 J=1,KK                                                       SIMIN.136
    2 P(J,I) = ANS(J) - SKP1                                            SIMIN.137
      DO 3 J=1,KK                                                       SIMIN.138
    3 P(J,J+1) = P(J,J+1) + S(1)                                        SIMIN.139
      GO TO 7                                                           SIMIN.140
c                                                                       SIMIN.141
    4 DO 5 I=1,KP1                                                      SIMIN.142
      DO 5 J=1,KK                                                       SIMIN.143
    5 P(J,I) = ANS(J) - ABS(S(J))*ONEKP1                                SIMIN.144
      DO 6 J=1,KK                                                       SIMIN.145
    6 P(J,J+1) = P(J,J+1) + ABS(S(J))                                   SIMIN.146
c                                                                       SIMIN.147
c     FUNCTION VALUES FOR INITIAL SIMPLEX                               SIMIN.148
c                                                                       SIMIN.149
    7 I1 = 1                                                            SIMIN.150
      DO 8 I=1,KP1                                                      SIMIN.151
      Y(I) = F(P(1,I))                                                  SIMIN.152
      IF (Y(I).GT.Y(I1)) I1 = I                                         SIMIN.153
    8 CONTINUE                                                          SIMIN.154
      YANS = F(ANS)                                                     SIMIN.155
      KOUNT = KP1+1                                                     SIMIN.156
      IF (YANS.GE.Y(I1)) GO TO 10                                       SIMIN.157
      Y(I1) = YANS                                                      SIMIN.158
      DO 9 J=1,KK                                                       SIMIN.159
    9 P(J,I1) = ANS(J)                                                  SIMIN.160
c                                                                       SIMIN.161
c     RE-START / NEXT ITERATION                                         SIMIN.162
c     IF K.LT.0 VALUES IN THE P AND Y ARRAYS (AND ONLY THESE VALUES)    SIMIN.163
c     WILL NOT HAVE BEEN DEFINED IN THIS CALL.  THIS IS NON-ANSI USAGE. SIMIN.164
c                                                                       SIMIN.165
c     FIRST FIND LARGEST, SECOND LARGEST, AND SMALLEST FUNCTION VALUES. SIMIN.166
c                                                                       SIMIN.167
   10 I1 = 1                                                            SIMIN.168
      IL = 1                                                            SIMIN.169
      DO 12 I=2,KP1                                                     SIMIN.170
      IF (Y(I).LT.Y(IL)) IL = I                                         SIMIN.171
      IF (Y(I).GT.Y(I1)) I1 = I                                         SIMIN.172
   12 CONTINUE                                                          SIMIN.173
      I2 = IL                                                           SIMIN.174
      DO 13 I=1,KP1                                                     SIMIN.175
      IF (I.EQ.I1) GO TO 13                                             SIMIN.176
      IF (Y(I).GT.Y(I2)) I2 = I                                         SIMIN.177
   13 CONTINUE                                                          SIMIN.178
c                                                                       SIMIN.179
c     COMPUTE CENTROID, LEAVING OUT P(*,I1)                             SIMIN.180
c                                                                       SIMIN.181
      DO 15 J=1,KK                                                      SIMIN.182
      SUM = 0.0                                                         SIMIN.183
      DO 14 I=1,KP1                                                     SIMIN.184
      IF (I.EQ.I1) GO TO 14                                             SIMIN.185
      SUM = SUM + P(J,I)                                                SIMIN.186
   14 CONTINUE                                                          SIMIN.187
   15 PC(J) = SUM*ONEK                                                  SIMIN.188
c                                                                       SIMIN.189
c     FORM REFLECTED POINT AND TEST                                     SIMIN.190
c                                                                       SIMIN.191
      DO 20 J=1,KK                                                      SIMIN.192
   20 PR(J) = PC(J) + ALPHA*(PC(J)-P(J,I1))                             SIMIN.193
      YR = F(PR)                                                        SIMIN.194
      KOUNT = KOUNT+1                                                   SIMIN.195
      IF (YR.LT.Y(IL)) GO TO 30                                         SIMIN.196
      IF (YR.GE.Y(I2)) GO TO 40                                         SIMIN.197
c                                                                       SIMIN.198
c     ACCEPT REFLECTED POINT                                            SIMIN.199
c                                                                       SIMIN.200
   21 Y(I1) = YR                                                        SIMIN.201
      DO 22 J=1,KK                                                      SIMIN.202
   22 P(J,I1) = PR(J)                                                   SIMIN.203
      GO TO 60                                                          SIMIN.204
c                                                                       SIMIN.205
c     EXPAND IN FAVORABLE DIRECTION AND TEST                            SIMIN.206
c                                                                       SIMIN.207
   30 DO 31 J=1,KK                                                      SIMIN.208
   31 PRR(J) = PR(J) + GAMMA*(PR(J)-PC(J))                              SIMIN.209
      YRR = F(PRR)                                                      SIMIN.210
      KOUNT = KOUNT+1                                                   SIMIN.211
      IF (YRR.GE.YR) GO TO 21                                           SIMIN.212
c                                                                       SIMIN.213
c     ACCEPT EXPANDED POINT                                             SIMIN.214
c                                                                       SIMIN.215
      Y(I1) = YRR                                                       SIMIN.216
      DO 32 J=1,KK                                                      SIMIN.217
   32 P(J,I1) = PRR(J)                                                  SIMIN.218
      GO TO 60                                                          SIMIN.219
c                                                                       SIMIN.220
c     DECIDE WHETHER TO ACCEPT REFLECTED POINT.                         SIMIN.221
c                                                                       SIMIN.222
   40 IF (YR.GE.Y(I1)) GO TO 42                                         SIMIN.223
      Y(I1) = YR                                                        SIMIN.224
      DO 41 J=1,KK                                                      SIMIN.225
   41 P(J,I1) = PR(J)                                                   SIMIN.226
c                                                                       SIMIN.227
c     TRY CONTRACTION.                                                  SIMIN.228
c                                                                       SIMIN.229
   42 DO 43 J=1,KK                                                      SIMIN.230
   43 PR(J) = PC(J) + BETA*(P(J,I1)-PC(J))                              SIMIN.231
      YCT = F(PR)                                                       SIMIN.232
      KOUNT = KOUNT+1                                                   SIMIN.233
      IF (YCT.GT.Y(I1)) GO TO 50                                        SIMIN.234
      Y(I1) = YCT                                                       SIMIN.235
      DO 44 J=1,KK                                                      SIMIN.236
   44 P(J,I1) = PR(J)                                                   SIMIN.237
      GO TO 60                                                          SIMIN.238
c                                                                       SIMIN.239
c     ALL EFFORTS FAILED.  SHRINK THE SIMPLEX ABOUT BEST POINT.         SIMIN.240
c                                                                       SIMIN.241
   50 DO 52 I=1,KP1                                                     SIMIN.242
      IF (I.EQ.IL) GO TO 52                                             SIMIN.243
      DO 51 J=1,KK                                                      SIMIN.244
   51 P(J,I) = 0.5*(P(J,I)+P(J,IL))                                     SIMIN.245
      Y(I) = F(P(1,I))                                                  SIMIN.246
   52 CONTINUE                                                          SIMIN.247
      KOUNT = KOUNT+KP1                                                 SIMIN.248
c                                                                       SIMIN.249
c     CHECK FOR CONVERGENCE                                             SIMIN.250
c                                                                       SIMIN.251
   60 IF (KOUNT.GE.NEV) GO TO 65                                        SIMIN.252
      IF (EPS.EQ.0.0) GO TO 10                                          SIMIN.253
      SUM = 0.0                                                         SIMIN.254
      DO 61 I=1,KP1                                                     SIMIN.255
   61 SUM = SUM + Y(I)                                                  SIMIN.256
      YAVG = SUM*ONEKP1                                                 SIMIN.257
      SUM = 0.0                                                         SIMIN.258
      DO 62 I=1,KP1                                                     SIMIN.259
   62 SUM = SUM + (Y(I)-YAVG)**2                                        SIMIN.260
      IF (EPS) 64,63,63                                                 SIMIN.261
   63 IF (SUM-TOL) 65,65,10                                             SIMIN.262
   64 IF (SUM-TOL*ABS(YAVG)) 65,65,10                                   SIMIN.263
c                                                                       SIMIN.264
c     CONVERGENCE OBTAINED.                                             SIMIN.265
c                                                                       SIMIN.266
c     COMPUTE CENTROID                                                  SIMIN.267
c                                                                       SIMIN.268
   65 DO 68 J=1,KK                                                      SIMIN.269
      SUM = 0.0                                                         SIMIN.270
      DO 67 I=1,KP1                                                     SIMIN.271
   67 SUM = SUM+P(J,I)                                                  SIMIN.272
   68 PC(J) = SUM*ONEKP1                                                SIMIN.273
      IF (S(1)) 73,69,69                                                SIMIN.274
c                                                                       SIMIN.275
c     COMPUTE S(1) AS AVERAGE DISTANCE OF VERTICES FROM CENTROID.       SIMIN.276
c                                                                       SIMIN.277
   69 DIST = 0.0                                                        SIMIN.278
      DO 71 I=1,KP1                                                     SIMIN.279
      SUM = 0.0                                                         SIMIN.280
      DO 70 J=1,KK                                                      SIMIN.281
   70 SUM = SUM + (P(J,I)-PC(J))**2                                     SIMIN.282
   71 DIST = DIST + SQRT(SUM)                                           SIMIN.283
      S(1) = DIST*ONEKP1                                                SIMIN.284
      GO TO 80                                                          SIMIN.285
c                                                                       SIMIN.286
c     COMPUTE S(J) AS AVERAGE DISTANCE IN J-TH DIMENSION OF             SIMIN.287
c     VERTICES FROM THE CENTROID.                                       SIMIN.288
c                                                                       SIMIN.289
   73 DO 75 J=1,KK                                                      SIMIN.290
      SUM = 0.0                                                         SIMIN.291
      DO 74 I=1,KP1                                                     SIMIN.292
   74 SUM = SUM + ABS(P(J,I)-PC(J))                                     SIMIN.293
   75 S(J) = SUM*ONEKP1                                                 SIMIN.294
      S(1) = -S(1)                                                      SIMIN.295
c                                                                       SIMIN.296
c     RETURN P(*,IL) AS ANSWER                                          SIMIN.297
c                                                                       SIMIN.298
   80 IL = 1                                                            SIMIN.299
      DO 82 I=2,KP1                                                     SIMIN.300
      IF (Y(I).LT.Y(IL)) IL = I                                         SIMIN.301
   82 CONTINUE                                                          SIMIN.302
      DO 84 J=1,KK                                                      SIMIN.303
   84 ANS(J) = P(J,IL)                                                  SIMIN.304
      NEV = KOUNT                                                       SIMIN.305
      RETURN                                                            SIMIN.306
c                                                                       SIMIN.307
c     ERROR MESSAGE                                                     SIMIN.308
c                                                                       SIMIN.309
c     99 CALL ERRCHK(39,39HIN SIMINA, S(1)=0. OR K IS LESS THAN 2.)        SIMIN.310
   99 write(6,*) 's(1)=0. or k is less than 2.'
      RETURN                                                            SIMIN.311
      END                                                               SIMIN.312



