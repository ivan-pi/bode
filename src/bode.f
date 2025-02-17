      SUBROUTINE BODE(XIN,XOUT,N,YN,YMIN,EMAX,XSTEP,MONIT,IMN,M1,IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YN(N), YMIN(N), YOLD1(75), ARR(75,20), TL(75,20), X1(75)
C***********************************************************************
C   THIS SUBROUTINE ATTEMPTS THE SOLUTION OF A SYSTEM OF N
C  FIRST ORDER DIFFERENTIAL EQUATIONS OF THE FROM
C
C                 B*DY(X)/DX = F(X,Y)
C
C   AT X=XOUT WHERE B IS A BANDED NXN MATRIX
C
C   THE USER NEEDS TO PROVIDE TWO SUBROUTINE TO DEFINE
C   THE SYSTEM (A)O-
C
C   1) SUBROUTINE PMULT(P,N,Q)
C      DIMENSION P(N),Q(N)
C
C   ON EXIT FROM PMULT THE ARRAY Q SHOULD CONTAIN THE
C   RESULT OF MULTIPLYING THE GIVEN VECTOR OF VALUES, P,
C   BY THE MATRIX B.
C
C   2) SUBROUTINE DERIV(Y,N,Q,X,H)
C      DIMENSION Y(N), Q(N)
C
C   ON EXIT FROM DERIV THE ARRAY Q SHOULD CONTAIN THE
C   DERIVATIVE VALUES F(X,Y)*H GIVEN THE VALUES OF THE
C   INDEPENDENT VARIABLE X AND THE ARRAY Y CONTAINING
C   THE N VALUES OF THE DEPENDENT VARIABLES
C
C   ARGUMENT LIST
C   -------------
C
C    XIN   - VALUE OF X AT WHICH THE INITIAL VALUES ARE GIVEN
C    XOUT  - VALUE OF X AT WHICH INTEGRATION IS TO BE TERMINATED
C           IF FAILURE OCCURS THIS CONTAINS, ON EXIT, THE VALEU
C           OF X AT WHICH THE LAST ACCEPTED APPROXIMATION WAS
C           OBTAINED PRIOR TO FAILURE.
C    N     - NUMBER OF EQUATIONS IN THE SYSTEM (>=3)
C    YN    - ARRAY OF LENTH N ON ENTRY CONTAINS THE VALUES OF
C           THE INITIAL CONDITIONS. ON EXIT CONTAINS THE COMPUTED
C           VALUES AT X=XOUT
C    YMIN  - ARRAY OF LENGTH N CONTAINING USER SUPPLIED VALUES OF
C           THE MAGNITUDES OF THE DEPENDENT VARIABLE BELOW WHICH
C           AN ABSOLUTE ERROR TEST IS TO BE USED FOR THE ESTIMATED
C           ERROR IN ANY PARTICULAR INTEGRATION STEP. THE ABSOLUTE
C           ERROR LIMIT USED IS EMAX*YMIN.
C    EMAX  - RELATIVE ACCURACY REQUIRED PER TIME STEP
C    XSTEP - THE SIZE OF THE FIRST INITIAL TIME STEP. IF THIS IS SET
C           ZERO THE INITIAL STEP IS SET TO ABS(XOUT-XIN)*0.25
C    MONIT - A USER SUPPLIED ROUTINE OF THE FORM
C
C           SUBROUTINE MONIT(Y,N,X)
C           DIMENSION Y(N)
C
C           WHICH ALLOWS THE USER TO MONITOR THE PROGRESS OF THE
C           INTEGRATION. THE ARRAY Y CONTAINS THE COMPUTED SOLUTION
C           AT TIME X.
C
C    IMN   - THE ROUTINE MONIT IS CALLED AFTER EVERY IMN SUCCESFUL
C           TIME STEPS UNLESS IMN<=0 WHEN IT IS NEVER CALLED
C
C    M1    - A MEASURE OF THE BANDWIDTH OF THE MATRIX B
C           DEFINED SUCH THATO-
C           M1 = MAX(C) S.T. B(I,J)=0 FOR MOD(I,J)=C
C           FOR ALL I,J=1(1)N.
C    IFAIL - ON EXIT IFAIL MAY TAKE THE VALUEO-
C                  0 - FOR SUCCESFUL INTEGRATION TO X=XOUT
C                  1 - IF THE TIME STEP HAS BEEN HALVED 20 TIMES
C                      SUCCESIVELY
C                  2 - IF THE STEP IN THE X-DIRECTION WOULD BE
C                      LESS THAN THE VALUE TMIN SET IN THE DATA
C                      STATEMENT IN ODE
C***********************************************************************
      DIMENSION DEL(75), FN(75), YN1(75), YOLD(75), PYP(75), FN2(75),
     * V1(75)
      DIMENSION PYC(75), INT(75), HJAC(75)
C***********************************************************************
C   DATA STATEMENT
C   --------------
C
C   TMIN  - SMALLEST TIME STEP SUCH THAT X+TMIN AND X
C           ARE DIFFERENT WITHIN THE MACHINE
C***********************************************************************
      DATA TMIN/1.0D-10/
C***********************************************************************
C INITIALIZE THE VARIABLES
C***********************************************************************
      IMON=IMN
      M2=M1+1
      M=2*M1+1
      ALPHA=0.55D0
      A=1.0D0-ALPHA
C***********************************************************************
C INITIALIZE THE OTHER REQUIRED VARIABLES
C***********************************************************************
      MON=0
      IDOHA=0
      NHALF=0
      RAT=0.0D0
      LSTEP=0
      IM=0
      IFSTEP=0
      IJAC=0
      B=1.0D0/6.0D0
      IF (XSTEP.EQ.0.0D0) XSTEP=ABS(XOUT-XIN)*0.25D0
C***********************************************************************
C SET UP INITIAL TIME STEP
C***********************************************************************
      H=XSTEP
C***********************************************************************
C TEST IF ONLY ONE STEP IS REQUIRED AND IF SO SET LSTEP
C***********************************************************************
      IF (ABS(XOUT-XIN)/H.GT.1.1D0) GO TO 10
      LSTEP=1
   10 X=XIN+H
C***********************************************************************
C PREPARE REQUIRED ARRAYS FOR FIRST STEP
C***********************************************************************
      CALL PMULT (YN,N,YN1)
   20 CALL DERIV (YN,N,FN,X,H)
C***********************************************************************
C USE INITIAL VALUES AS FIRST PREDICTED VALEUS
C***********************************************************************
      DO 30 I=1,N
      YOLD(I)=YN(I)
      HJAC(I)=1.0D-3
   30 PYP(I)=YN(I)
C***********************************************************************
C SOLUTION OF THE CORRECTOR STEP IN THE FORM OF A SYSTEM OF
C NONLINEAR EQUATIONS
C***********************************************************************
   40 DO 50 I=1,N
      X1(I)=PYP(I)
   50 DEL(I)=YN1(I)+A*FN(I)
      H1=H*ALPHA
      IFAIL=0
      CALL NONLIN(X1,N,DEL,IJAC,HJAC,X,ARR,TL,INT,M,M1,M2,EMAX,H1,IFAIL)
      DO 60 I=1,N
   60 PYC(I)=X1(I)
C***********************************************************************
C TEST FOR ERROR 0-
C    IFAIL=1  - FAILURE IN FUN EVALUATION
C    IFAIL=2  - NO CONVERGENCE  HALVE STEP SIZE
C    IFAIL=3  - DIVERGENCE IN N-R ITERATION
C***********************************************************************
      IF (IFAIL.NE.0) GO TO 210
C***********************************************************************
C CALCULATE REQUIRED ARRAYS FROM CORRECTED VALUES
C***********************************************************************
      CALL DERIV (PYC,N,FN2,X,H)
C***********************************************************************
C CALCULATE VECTORS REQUIRED FOR LOCAL ERROR ESTIMATES
C***********************************************************************
      DO 70 I=1,N
   70 PYP(I)=FN2(I)-FN(I)
      CALL TSOL (ARR,TL,M1,M2,M,N,INT,PYP,DEL)
C***********************************************************************
C IF FIRST STEP OMIT THESE CALCULATIONS FOR LOCAL ERROR
C***********************************************************************
      IF (IFSTEP.EQ.0) GO TO 110
      IF (RAT .GE. 0.5D0) GO TO 80
      CON=B*RAT/(1.0D0+RAT)
      GO TO 90
   80 CON=RAT*(-B+ALPHA*(1.0D0-ALPHA))/(1.0D0+2.0D0*ALPHA*(RAT-1.0D0))
   90 DO 100 I=1,N
  100 PYP(I)=CON*(DEL(I)-V1(I)*RAT)
      GO TO 130
C***********************************************************************
C ZERO PYP ARRAY IF FIRST STEP
C***********************************************************************
  110 DO 120 I =1,N
  120 PYP(I)=0.0D0
  130 REL=0.0D0
C***********************************************************************
C CALCULATE MAXIMUM RELATIVE LOCAL TRUNCATION ERROR ESTIMATES
C***********************************************************************
      DO 140 I=1,N
      DEL1=ABS(PYP(I)+(ALPHA-0.5D0)*DEL(I))/(ABS(PYC(I))+YMIN(I))
      IF (DEL1.GT.REL) REL=DEL1
  140 CONTINUE
      REL=REL/EMAX
C***********************************************************************
C IF REL > 1  HALVE STEP SIZE
C***********************************************************************
      IF (REL.GT.1.0D0) GO TO 210
C***********************************************************************
C SUCCESFUL STEP TEST FOR LAST STEP
C***********************************************************************
      IF (LSTEP.EQ.1) GO TO 300
      RAT=1.0D0
C***********************************************************************
C DECIDE FROM VALUE OF REL WHETHER STEP SIZE SHOULD REMAIN
C UNCHANGED OR BE DOUBLED
C***********************************************************************
      IDOHA=0
      IF (REL.GT.0.2D0) GO TO 150
      RAT=2.0D0
      IDOHA=1
      IJAC=0
C***********************************************************************
C TEST IF NEXT STEP WILL BE LAST STEP
C***********************************************************************
  150 IF (ABS(XOUT-X)/(H*RAT).GT.1.1D0) GO TO 160
      LSTEP=1
      IJAC=0
      XSTEP=H*RAT
      RAT=(XOUT-X)/H
  160 H=H*RAT
      IFSTEP=1
      NHALF=0
      X=X+H
      CON=1.0D0+ALPHA*(RAT-1.0D0)
C***********************************************************************
C REST ARRAYS FOR NEXT STEP
C***********************************************************************
      DO 170 I=1,N
      V1(I)=RAT*DEL(I)
      YOLD1(I)=PYC(I)-YOLD(I)
      PYP(I)=RAT*YOLD(I)+PYC(I)+CON*V1(I)
      YOLD(I)=PYC(I)
      YN(I)=YN1(I)
  170 FN(I)=RAT*FN2(I)
      CALL PMULT (PYC,N,YN1)
C***********************************************************************
C OUTPUT TO ROUTINE MONIT IF REQURIED
C***********************************************************************
      IF (IMON.LT.1) GO TO 190
      MON=MON+1
      IF (IMON.NE.MON) GO TO 190
      MON=0
      IM=-1
      T=X-H
      DO 180 I=1,N
  180 YN(I)=PYC(I)
      IM=0
      CALL MONIT (YN,N,T)
C***********************************************************************
C CALCULATE NEW PREDICTOR
C***********************************************************************
  190 DO 200 I=1,N
  200 PYP(I)=(1.0D0+RAT)*PYC(I)-RAT*YOLD(I)+CON*V1(I)
      GO TO 40
C***********************************************************************
C STEP REJECTED
C***********************************************************************
  210 NHALF=NHALF+1
C***********************************************************************
C IF STEP SIZE HAS BEEN DOUBLES TO H (IDOHA SET 1) AND THE
C NEXT STEP ATTEMPTS TO HALVE AGAIN NEW STEP SET 3*H/4
C***********************************************************************
      IF (IDOHA.EQ.1) GO TO 220
      RAT1=0.5D0
      GO TO 230
  220 RAT1=0.75D0
      IDOHA=0
C***********************************************************************
C MAXIMUM NUMBER OF SUCCESIVE HALVINGS ALLOWED IS 20
C***********************************************************************
  230 IF (NHALF.LT.20) GO TO 240
      IFAIL=1
      GO TO 280
  240 X=X-(1.0D0-RAT1)*H
      H=H*RAT1
      IF (H.GT.TMIN) GO TO 250
      IFAIL=2
      GO TO 280
  250 IJAC=0
C***********************************************************************
C RESET ARRAYS FOR HALVED STEP SIZE AND CARRY ON
C***********************************************************************
      RAT=RAT*RAT1
      DO 260 I=1,N
  260 FN(I)=FN(I)*RAT1
      IF (IFSTEP.EQ.0) GO TO 20
      CON=1.0D0+ALPHA*(RAT-1.0D0)
      DO 270 I=1,N
      V1(I)=V1(I)*RAT1
  270 PYP(I)=YOLD(I)-YOLD1(I)*RAT+CON*V1(I)
      GO TO 40
C***********************************************************************
C NON SUCCESFUL EXIT
C***********************************************************************
  280 CONTINUE
      XOUT=X-H
      DO 290 I=1,N
  290 YN(I)=YOLD(I)
      RETURN
C***********************************************************************
C SUCCESSFUL EXIT  RETURN VALUES
C  PUT RESULTS INTO YN ARRAY FOR OUTPUT
C***********************************************************************
  300 XOUT=X
      DO 310 I=1,N
  310 YN(I)=PYC(I)
      IFAIL=0
      RETURN
      END
      SUBROUTINE NONLIN(X1,N,DEL,IJAC,HJAC,X,A,TL,INT,M,M1,M2,EMAX,
     *  H,IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***********************************************************************
C SUBROUTINE NONLIN SOLVES THE SYSTEM OF NON-LINEAR
C EQUATIONS GENERATED BY THE CORRECTOR STEP USING
C THE MODIFIED NEWTON-RAPHSON ITERATION. IT ALSO PRODUCES A
C FACTORIZED FORM OF THE JACOBIAN IN THE ARRAYS A,TL.
C ACCOUNT IS TAKEN OF THE BANDED STRUCTURE OF THE
C JACOBIAN TO MINIMIZE THE NUMBER OF FUNCTION EVALUATIONS
C REQUIRED TO APPROXIMATE THE DERIVATIVES
C***********************************************************************
      DIMENSION X1(N),INT(N),DEL(N),HJAC(N),A(75,M),TL(75,M2)
      DIMENSION Y(75),F(75),GG(75),V2(75)
      DATA HMN,ETA/1.0D-5,1.0D-4/
C***********************************************************************
C DATA STATEMENT
C --------------
C
C THE MAGNITUDE OF THE DIFFERENCE H USED IN THE APPROXIMATION
C OF THE DERIVATIVE USINGO-
C              (F(X+H)-F(X))/H
C IS   MAX(HMN,ETA*MOD(X))
C
C THE MODIFIED NEWTON-RAPHSON ITERATION IS ASSUMED TO HAVE
C DIVERGED IF THE MAGNITUE OF ANY COMPONENT OF THE
C INCREMENTAL VECTOR IS GREATER THAN BIG.
C***********************************************************************
      DATA BIG/1.0D10/
C***********************************************************************
C MAXF - MAXIMUM NUMBER OF ITERATIONS ALLOWED BEFORE NON-
C        CONVERGENCE IS ASSUMED
C***********************************************************************
      MAXF=N+10
C***********************************************************************
C IF STEP-SIZE HAS BEEN CHANGED REAPPROXIMATE INVERSE JACOBIAN
C***********************************************************************
      IF (IJAC.EQ.1) GO TO 120
      CALL PMULT (X1,N,Y)
      CALL DERIV (X1,N,GG,X,H)
      DO 10 I=1,N
      F(I)=Y(I)-GG(I)-DEL(I)
      Y(I)=X1(I)
      HJAC(I)=MIN(1.0D0/HMN,MAX(ETA*ABS(X1(I)),HMN))
      DO 10 J=1,M
      A(I,J)=0.0D0
   10 CONTINUE
      K=1
C***********************************************************************
C INCREMENT RELEVANT VARIABLES - USE BANDED STRUCTURE
C***********************************************************************
   20 J=K
   30 IF (J.GT.N) GO TO 40
      Y(J)=X1(J)+HJAC(J)
      J=J+M
      GO TO 30
   40 CALL PMULT (Y,N,V2)
      CALL DERIV (Y,N,GG,X,H)
      DO 50 I=1,N
   50 GG(I)=V2(I)-GG(I)-DEL(I)
      J=K
   60 IF (J.GT.N) GO TO 70
      Y(J)=X1(J)
      J=J+M
      GO TO 60
   70 J=K
C***********************************************************************
C APPROXIMATE DERIVATIVES BY DIFFERENCING
C***********************************************************************
   80 JJ=M
      IMIN=-M1
      IF (J.GE.(M1+1)) GO TO 90
      IMIN=-J+1
      JJ=M1+J
   90 IMAX=M1
      IF (J.GT.(N-M1)) IMAX=-J+N
      I=IMIN
  100 IJ=I+J
      A(IJ,JJ)=(GG(IJ)-F(IJ))/HJAC(J)
      JJ=JJ-1
      I=I+1
      IF (I.LE.IMAX) GO TO 100
      J=J+M
      IF (J.LE.N) GO TO 80
      K=K+1
      IF (M.GT.N.AND.K.GT.N) GO TO 110
      IF (K.LE.M) GO TO 20
  110 CONTINUE
C***********************************************************************
C OBTAIN L-U DECOMPOSITION
C***********************************************************************
      CALL LTRI (A,TL,N,INT,M,M1,M2,IFAIL)
      IF (IFAIL.NE.0) RETURN
      IF (IJAC.GE.0) IJAC=1
  120 DO 160 II=1,MAXF
      SUM1=0.0D0
      IFLAG=0
C***********************************************************************
C DO NOT RECALCULATE F IF  JACOBIAN HAS BEEN REAPPROXIMATED
C***********************************************************************
      CALL PMULT (X1,N,V2)
      CALL DERIV (X1,N,GG,X,H)
      DO 130 I=1,N
  130 F(I)=V2(I)-GG(I)-DEL(I)
C***********************************************************************
C SOLVE FOR INCREMENTAL VECTOR
C***********************************************************************
      CALL TSOL (A,TL,M1,M2,M,N,INT,F,GG)
      DO 150 I=1,N
      GGI=GG(I)
      IF (ABS(GGI).LT.BIG) GO TO 140
      IFAIL=3
      GO TO 170
  140 X1(I)=X1(I)-GGI
      IF (ABS(GGI).GT.(ABS(X1(I))+1.0D-6)*EMAX) IFLAG=1
  150 CONTINUE
      IF (IFLAG.EQ.0) GO TO 170
  160 CONTINUE
      IFAIL=2
  170 CONTINUE
      RETURN
      END
      SUBROUTINE LTRI (A,TL,N,INT,M21,M1,M3,IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***********************************************************************
C FORTRAN VERSION OF ROUTINE BANDET1
C MARTIN A WILKINSON NUM. MATH. VOL 9 P279-301 (1967)
C***********************************************************************
      DIMENSION A(75,M21), TL(75,M3), INT(N)
      L=M1
      DO 20 I=1,M1
      M4=M1+I
      DO 10 J=1,M4
      JL=J+L
   10 A(I,J)=A(I,JL)
      L=L-1
      M21L=M21-L
      DO 20 J=M21L,M21
   20 A(I,J)=0.0D0
      L=M1
      DO 100 K=1,N
      X=A(K,1)
      I=K
      IF (L.LT.N) L=L+1
      K1=K+1
      IF (K1.GT.L) GO TO 40
      DO 30 J=K1,L
      IF (ABS(A(J,1)).LT.ABS(X)) GO TO 30
      X=A(J,1)
      I=J
   30 CONTINUE
   40 INT(K)=I
      IF (X.NE.0.0D0) GO TO 50
      IFAIL=1
      RETURN
   50 IF (I.EQ.K) GO TO 70
      DO 60 J=1,M21
      X=A(K,J)
      A(K,J)=A(I,J)
   60 A(I,J)=X
   70 K1=K+1
      IF (K1.GT.L) GO TO 100
      DO 90 I=K1,L
      X=A(I,1)/A(K,1)
      IK=I-K
      TL(K,IK)=X
      DO 80 J=2,M21
   80 A(I,J-1)=A(I,J)-X*A(K,J)
   90 A(I,M21)=0.0D0
  100 CONTINUE
      RETURN
      END
      SUBROUTINE TSOL(A,TL,M1,M3,M21,N,INT,XIN,XOUT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***********************************************************************
C   FORTRAN VERSION OF ROUTINE BANSOL1
C   MARTIN A WILKINSION NUM. MATH. VOL 9 P279-301 (1967)
C***********************************************************************
      DIMENSION A(75,M21), TL(75,M3), INT(N), XIN(N), XOUT(N)
      DO 10 I=1,N
   10 XOUT(I)=XIN(I)
      L=M1
      DO 40 K=1,N
      I=INT(K)
      IF (I.EQ.K) GO TO 20
      X=XOUT(K)
      XOUT(K)=XOUT(I)
      XOUT(I)=X
   20 IF (L.LT.N) L=L+1
      K1=K+1
      IF (K1.GT.L) GO TO 40
      DO 30 I=K1,L
      IK=I-K
      X=TL(K,IK)
   30 XOUT(I)=XOUT(I)-X*XOUT(K)
   40 CONTINUE
      L=1
      DO 70 II=1,N
      I=N-II+1
      X=XOUT(I)
      IF (L.EQ.1) GO TO 60
      DO 50 K=2,L
      K1=K+I-1
   50 X=X-A(I,K)*XOUT(K1)
   60 XOUT(I)=X/A(I,1)
      IF (L.LT.M21) L=L+1
   70 CONTINUE
      RETURN
      END
