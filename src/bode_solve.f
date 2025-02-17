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


      SUBROUTINE BODE(XIN,XOUT,N,YN,YMIN,EMAX,XSTEP,MONIT,IMN,M1,IFAIL)
      use bode_mod, only: wp, ld, tsol
      implicit none
      real(wp), intent(in) :: xin
      real(wp), intent(inout) :: xout
      integer, intent(in) :: n
      real(wp), intent(inout) :: yn(n)
      real(wp), intent(in) :: ymin(n)
      real(wp), intent(in) :: emax
      real(wp), intent(inout) :: xstep
      procedure() :: monit
      integer, intent(in) :: imn, m1
      integer, intent(out) :: ifail

      interface
      SUBROUTINE NONLIN(X1,N,DEL,IJAC,HJAC,X,A,TL,LD,IPIV,M,M1,M2,
     *EMAX,H,IFAIL)
        import wp
         implicit none
         integer, intent(in) :: n, ld, m, m1, m2
         real(wp), intent(inout) :: x1(n)
         real(wp), intent(in) :: del(n)
         integer, intent(inout) :: ijac
         real(wp), intent(inout) :: hjac(n)
         real(wp), intent(in) :: x, h   ! scalars passed to routine deriv

         ! storage for the factorized Jacobian
         real(wp), intent(inout) :: a(ld,m), tl(ld,m2)
         integer, intent(inout) :: ipiv(n)

         real(wp), intent(in) :: emax
         integer, intent(inout) :: ifail
      end subroutine
      end interface

C
C LOCAL VARIABLES
C
      real(wp) :: yold1(ld), arr(ld,20), tl(ld,20), x1(ld)
      real(wp) :: del(ld), fn(ld), yn1(ld), yold(ld), pyp(ld), fn2(ld)
      real(wp) :: v1(ld), pyc(ld), hjac(ld)
      integer :: ipiv(ld)
      logical :: laststep, firststep, halve_step
C
      real(wp) :: con, del1, h, h1, t, rel, rat, rat1, x

      integer :: i, idoha, ijac, im, imon, mon, m, m2, nhalf

C   TMIN  - SMALLEST TIME STEP SUCH THAT X+TMIN AND X
C           ARE DIFFERENT WITHIN THE MACHINE
C
      real(wp), parameter :: tmin = 1.0e-10_wp
      real(wp), parameter :: B = 1.0_wp/6.0_wp
      real(wp), parameter :: alpha = 0.55_wp
      real(wp), parameter :: a = 1.0_wp - alpha
C***********************************************************************
C INITIALIZE THE VARIABLES
C***********************************************************************
      IMON = IMN
      M2 = M1 + 1
      M = 2*M1 + 1
      if (m > 20) error stop "Bandwidth exceeds 20."
C***********************************************************************
C INITIALIZE THE OTHER REQUIRED VARIABLES
C***********************************************************************
      MON = 0
      IDOHA = 0
      NHALF = 0
      RAT = 0.0_wp
      firststep = .true.
      IM = 0
      IJAC = 0
C***********************************************************************
C SET UP INITIAL TIME STEP
C***********************************************************************
      if (xstep == 0.0_wp) xstep = abs(xout - xin)*0.25_wp
      h = xstep
C***********************************************************************
C TEST IF ONLY ONE STEP IS REQUIRED AND IF SO SET LASTSTEP
C***********************************************************************
      laststep = abs(xout - xin)/h <= 1.1_wp
      x = xin + h
C***********************************************************************
C PREPARE REQUIRED ARRAYS FOR FIRST STEP
C***********************************************************************
      call pmult(yn,n,yn1)
      call deriv(yn,n,fn,x,h)
C***********************************************************************
C USE INITIAL VALUES AS FIRST PREDICTED VALEUS
C***********************************************************************
      do i = 1, n
        yold(i) = yn(i)
        hjac(i) = 1.0e-3_wp
        pyp(i) = yn(i)
      end do

      solve: do
!
! solution of the corrector step in the form of a system of
! nonlinear equations
!
        attempt_step: block
          do i =1,n
            x1(i) = pyp(i) ! y
            del(i) = yn1(i) + a*fn(i) ! B*y + (1 - gamma)*(h*f(y))
          end do
          h1 = h*alpha
          ifail = 0
          call nonlin(x1,n,del,ijac,hjac,x,arr,tl,ld,ipiv,m,m1,m2,emax,
     *      h1,ifail)
          do i = 1, n
            pyc(i)=x1(i)
          end do

          halve_step = IFAIL /= 0
          if (halve_step) exit attempt_step

C***********************************************************************
C CALCULATE REQUIRED ARRAYS FROM CORRECTED VALUES
C***********************************************************************
          CALL DERIV(PYC,N,FN2,X,H) ! fn2 := h*f(x,pyc)
C***********************************************************************
C CALCULATE VECTORS REQUIRED FOR LOCAL ERROR ESTIMATES
C***********************************************************************
          do i = 1, n
            pyp(i) = fn2(i) - fn(i)
          end do
          call tsol(arr,tl,ld,m1,m2,m,n,ipiv,pyp,del) ! del :=
          if (firststep) then
            ! first step
            do i = 1, n
                pyp(i) = 0.0_wp
            end do
          else
              ! CALCULATIONS FOR LOCAL ERROR, omitted if first step
              if (rat >= 0.5_wp) then
                con = rat*(-b + alpha*(1 - alpha)) /
     *                   (1 + 2*alpha*(rat - 1))
              else
                con = b*rat/(1 + rat)
              end if
              do i = 1, n
                pyp(i) = con*(del(i) - v1(i)*rat)
              end do
          end if

C***********************************************************************
C CALCULATE MAXIMUM RELATIVE LOCAL TRUNCATION ERROR ESTIMATES
C***********************************************************************
          rel = 0.0_wp
          do i = 1, n
            del1 = abs(pyp(i) + (alpha - 0.5_wp)*del(i)) /
     *             (abs(pyc(i)) + ymin(i))
            rel = max(rel,del1)
          end do
          rel = rel/emax
C***********************************************************************
C IF REL > 1  HALVE STEP SIZE
C***********************************************************************
          halve_step = REL > 1.0_wp
          if (halve_step) exit attempt_step
C***********************************************************************
C SUCCESFUL STEP TEST FOR LAST STEP
C***********************************************************************
          IF (laststep) then
        ! last step, normal exit
C***********************************************************************
C SUCCESSFUL EXIT  RETURN VALUES
C  PUT RESULTS INTO YN ARRAY FOR OUTPUT
C***********************************************************************
            xout = x
            do i = 1, n
              yn(i) = pyc(i)
            end do
            ifail = 0
            return
          end if

          RAT = 1.0_wp
C***********************************************************************
C DECIDE FROM VALUE OF REL WHETHER STEP SIZE SHOULD REMAIN
C UNCHANGED OR BE DOUBLED
C***********************************************************************
          idoha = 0
          if (rel <= 0.2_wp) then
            rat = 2.0_wp
            idoha = 1
            ijac = 0
          end if
C***********************************************************************
C TEST IF NEXT STEP WILL BE LAST STEP
C***********************************************************************
          if (abs(xout-x)/(h*rat) <= 1.1_wp) then
            laststep = .true.
            ijac = 0
            xstep = h*rat
            rat = (xout - x)/h
          end if
          h = h*rat
          firststep = .false.
          nhalf = 0
          x = x + h
          con = 1 + alpha*(rat - 1)
C***********************************************************************
C RESET ARRAYS FOR NEXT STEP
C***********************************************************************
          do i = 1, n
            v1(i) = rat*del(i)
            yold1(i) = pyc(i) - yold(i)
            pyp(i) = rat*yold(i) + pyc(i) + con*v1(i)
            yold(i) = pyc(i)
            yn(i) = yn1(i)
            fn(i) = rat*fn2(i)
          end do
          call pmult(pyc,n,yn1) ! yn1 := B pyc
c***********************************************************************
c output to routine monit if requried
c***********************************************************************
          if (imon > 0) then
            mon = mon + 1
            ! output only every imon steps
            if (imon == mon) then
              mon = 0
              t = x - h
              do i = 1, n
                yn(i) = pyc(i)
              end do
              call monit(yn,n,t)
            end if
          end if
C***********************************************************************
C CALCULATE NEQ PREDICTOR
C***********************************************************************
          do i = 1, n
            pyp(i) = (1 + rat)*pyc(i) - rat*yold(i) + con*v1(i)
          end do

        end block attempt_step

        if (halve_step) then
C***********************************************************************
C STEP REJECTED
C***********************************************************************
          nhalf = nhalf + 1
C***********************************************************************
C IF STEP SIZE HAS BEEN DOUBLES TO H (IDOHA SET 1) AND THE
C NEXT STEP ATTEMPTS TO HALVE AGAIN NEW STEP SET 3*H/4
C***********************************************************************
          if (idoha == 1) then
            rat1 = 0.75_wp
            idoha = 0
          else
            rat1 = 0.5_wp
          end if
C***********************************************************************
C MAXIMUM NUMBER OF SUCCESIVE HALVINGS ALLOWED IS 20
C***********************************************************************
          if (nhalf >= 20) then
            ifail = 1
            exit solve
          end if

          ! calculate new time and stepsize
          x = x - (1 - rat1)*h
          h = h*rat1

          if (h <= tmin) then
            ifail = 2
            exit solve
          end if

          ijac = 0 ! Request Jacobian
!
! reset arrays for halved step size and carry on
!
          rat = rat*rat1
          do i = 1, n
            fn(i) = fn(i)*rat1
          end do
          if (firststep) then
            ! halving required on first step
            call deriv(yn,n,fn,x,h)
            do i = 1, n
              yold(i) = yn(i)
              hjac(i) = 1.0e-3_wp
              pyp(i) = yn(i)
            end do
          else
            con = 1 + alpha*(rat - 1)
            do i = 1, n
              v1(i) = v1(i)*rat1
              pyp(i) = yold(i) - yold1(i)*rat + con*v1(i)
            end do
          end if
        end if ! halve_step

      end do solve
!
! Non-succesful exit
!
      failed: block
        xout = x - h
        do i = 1, n
          yn(i) = yold(i)
        end do
      end block failed

      return
      end