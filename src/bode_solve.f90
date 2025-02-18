!
!  THIS SUBROUTINE ATTEMPTS THE SOLUTION OF A SYSTEM OF N
!  FIRST ORDER DIFFERENTIAL EQUATIONS OF THE FROM
!
!                 B*DY(X)/DX = F(X,Y)
!
!   AT X=XOUT WHERE B IS A BANDED NXN MATRIX
!
!   THE USER NEEDS TO PROVIDE TWO SUBROUTINE TO DEFINE
!   THE SYSTEM (A)O-
!
!   1) SUBROUTINE PMULT(P,N,Q)
!      DIMENSION P(N),Q(N)
!
!   ON EXIT FROM PMULT THE ARRAY Q SHOULD CONTAIN THE
!   RESULT OF MULTIPLYING THE GIVEN VECTOR OF VALUES, P,
!   BY THE MATRIX B.
!
!   2) SUBROUTINE DERIV(Y,N,Q,X,H)
!      DIMENSION Y(N), Q(N)
!
!   ON EXIT FROM DERIV THE ARRAY Q SHOULD CONTAIN THE
!   DERIVATIVE VALUES F(X,Y)*H GIVEN THE VALUES OF THE
!   INDEPENDENT VARIABLE X AND THE ARRAY Y CONTAINING
!   THE N VALUES OF THE DEPENDENT VARIABLES
!
!   ARGUMENT LIST
!   -------------
!
!    XIN   - VALUE OF X AT WHICH THE INITIAL VALUES ARE GIVEN
!    XOUT  - VALUE OF X AT WHICH INTEGRATION IS TO BE TERMINATED
!           IF FAILURE OCCURS THIS CONTAINS, ON EXIT, THE VALEU
!           OF X AT WHICH THE LAST ACCEPTED APPROXIMATION WAS
!           OBTAINED PRIOR TO FAILURE.
!    N     - NUMBER OF EQUATIONS IN THE SYSTEM (>=3)
!    YN    - ARRAY OF LENTH N ON ENTRY CONTAINS THE VALUES OF
!           THE INITIAL CONDITIONS. ON EXIT CONTAINS THE COMPUTED
!           VALUES AT X=XOUT
!    YMIN  - ARRAY OF LENGTH N CONTAINING USER SUPPLIED VALUES OF
!           THE MAGNITUDES OF THE DEPENDENT VARIABLE BELOW WHICH
!           AN ABSOLUTE ERROR TEST IS TO BE USED FOR THE ESTIMATED
!           ERROR IN ANY PARTICULAR INTEGRATION STEP. THE ABSOLUTE
!           ERROR LIMIT USED IS EMAX*YMIN.
!    EMAX  - RELATIVE ACCURACY REQUIRED PER TIME STEP
!    XSTEP - THE SIZE OF THE FIRST INITIAL TIME STEP. IF THIS IS SET
!           ZERO THE INITIAL STEP IS SET TO ABS(XOUT-XIN)*0.25
!    MONIT - A USER SUPPLIED ROUTINE OF THE FORM
!
!           SUBROUTINE MONIT(Y,N,X,NHALF,R)
!           DIMENSION Y(N)
!
!           WHICH ALLOWS THE USER TO MONITOR THE PROGRESS OF THE
!           INTEGRATION. THE ARRAY Y CONTAINS THE COMPUTED SOLUTION
!           AT TIME X.
!
!    IMN   - THE ROUTINE MONIT IS CALLED AFTER EVERY IMN SUCCESFUL
!           TIME STEPS UNLESS IMN<=0 WHEN IT IS NEVER CALLED
!
!    IOPT  - INTEGER ARRAY OF LENGTH 3 THAT SPECIFIES PROPERTIES OF
!            THE PROBLEM
!
!    IOPT(1): A MEASURE OF THE BANDWIDTH OF THE MATRIX B
!             DEFINED SUCH THATO-
!             M1 = MAX(C) S.T. B(I,J)=0 FOR MOD(I,J)=C
!             FOR ALL I,J=1(1)N.
!
!             M1 MUST REMAIN THE FIRST ELEMENT OF IOPT FOR BACKWARD
!             COMPATIBILITY
!
!    IOPT(2): IF EQUAL TO 1, A USER-PROVIDED JACOBIAN SUBROUTINE WILL
!             BE USED. THE USER SUBROUTINE IS OF THE FORM
!
!             SUBROUTINE PJACB(X,X1,N,A,LD,M1,H)
!               INTEGER N, LD, M1
!               DOUBLE PRECISION X, X1(N), A(LD,2*M1+1), H
!             END SUBROUTINE
!
!             THE PROCEDURE SHOULD EVALUATE THE MATRIX
!               B - h J
!             AND STORE THE VALUES IN ARRAY A USING BANDED
!             COLUMN STORAGE, ASSUMING A STRUCTURALLY SYMMETRIC
!             BANDED MATRIX WITH BANDWIDTH ML = MU = M1.
!
!    IOPT(3): USE DENSE OR BANDED JACOBIAN
!             = 0, dense jacobian is used
!             = 1, banded jacobian is used (structurally-symmetric)
!
!    IFAIL - ON EXIT IFAIL MAY TAKE THE VALUEO-
!                  0 - FOR SUCCESFUL INTEGRATION TO X=XOUT
!                  1 - IF THE TIME STEP HAS BEEN HALVED 20 TIMES
!                      SUCCESIVELY
!                  2 - IF THE STEP IN THE X-DIRECTION WOULD BE
!                      LESS THAN THE VALUE TMIN SET IN THE DATA
!                      STATEMENT IN ODE
!
subroutine bode(xin,xout,n,yn,ymin,emax,xstep,monit,imn,iopt,ifail)
  use bode_mod, only: wp, ld, tsol, pmonit, nfev, njev, nlu, nbsol
  implicit none
  real(wp), intent(in) :: xin
  real(wp), intent(inout) :: xout
  integer, intent(in) :: n
  real(wp), intent(inout) :: yn(n)
  real(wp), intent(in) :: ymin(n)
  real(wp), intent(in) :: emax
  real(wp), intent(inout) :: xstep
  procedure(pmonit) :: monit
  integer, intent(in) :: imn
  integer, intent(in) :: iopt(3)
  integer, intent(out) :: ifail

  interface
    subroutine nonlin(x1,n,del,ijac,hjac,user_jac,dense_jac,x,a,tl,ld,ipiv,&
                      m,m1,m2,emax,h,ifail)
      import wp
      implicit none
      integer, intent(in) :: n, ld, m, m1, m2
      real(wp), intent(inout) :: x1(n)
      real(wp), intent(in) :: del(n)
      integer, intent(inout) :: ijac
      real(wp), intent(inout) :: hjac(n)
      logical, intent(in) :: user_jac, dense_jac
      real(wp), intent(in) :: x, h   ! scalars passed to routine deriv

      ! storage for the factorized jacobian
      real(wp), intent(inout) :: a(ld,m), tl(ld,m2)
      integer, intent(inout) :: ipiv(n)

      real(wp), intent(in) :: emax
      integer, intent(inout) :: ifail
    end subroutine
  end interface

  ! LAPACK procedures
  external :: dgetrs
!
! local variables
!
  real(wp) :: yold1(ld), arr(ld,20), tl(ld,20), x1(ld)
  real(wp) :: del(ld), fn(ld), yn1(ld), yold(ld), pyp(ld), fn2(ld)
  real(wp) :: v1(ld), pyc(ld), hjac(ld)
  integer :: ipiv(ld)
  logical :: laststep, firststep, halve_step

  real(wp) :: con, h, h1, t, rel, rat, rat1, x

  integer :: i, idoha, ijac, im, imon, mon, m, m1, m2, nhalf, la_info
  logical :: user_jac, dense_jac
!
! tmin  - smallest time step such that x+tmin and x
!         are different within the machine
  real(wp), parameter :: tmin = 1.0e-10_wp

  real(wp), parameter :: b = 1.0_wp/6.0_wp
  real(wp), parameter :: alpha = 0.55_wp
  real(wp), parameter :: a = 1.0_wp - alpha

  m1 = iopt(1)
  user_jac = iopt(2) == 1
  dense_jac = iopt(3) == 0

  if (dense_jac .neqv. user_jac) error stop "FD Jacobian only supported with banded solver"

  nfev = 0
  njev = 0
  nlu = 0
  nbsol = 0
!
! initialize the variables
!
  imon = imn
  m2 = m1 + 1
  m = 2*m1 + 1
  if (m > 20) error stop "bandwidth exceeds limit of 20"
  if (n > 20) error stop "size exceeds workspace of arr"
!
! initialize the other required variables
!
  mon = 0
  idoha = 0
  nhalf = 0
  rat = 0.0_wp
  firststep = .true.
  im = 0
  ijac = 0
!
! set up initial time step
!
  if (xstep == 0.0_wp) xstep = abs(xout - xin)*0.25_wp
  h = xstep
!
! test if only one step is required and if so set laststep
!
  laststep = abs(xout - xin)/h <= 1.1_wp
  x = xin + h
!
! prepare required arrays for first step
!
  call pmult(yn,n,yn1)
  call deriv(yn,n,fn,x,h)
  nfev = nfev + 1
!
! use initial values as first predicted valeus
!
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
      call nonlin (x1,n,del,ijac,hjac,user_jac,dense_jac,x,arr,tl,ld,ipiv,&
                   m,m1,m2,emax,h1,ifail)

      do i = 1, n
        pyc(i)=x1(i)
      end do

      halve_step = IFAIL /= 0
      if (halve_step) exit attempt_step

    !
    ! calculate required arrays from corrected values
    !
      call deriv(pyc,n,fn2,x,h) ! fn2 := h*f(x,pyc)
      nfev = nfev + 1
    !
    ! calculate vectors required for local error estimates
    !
      do i = 1, n
        pyp(i) = fn2(i) - fn(i)
      end do

      if (dense_jac) then
        del(1:n) = pyp(1:n)
        call dgetrs('N',n,1,arr,ld,ipiv,del,ld,la_info)
        if (la_info /= 0) then
          write(*,'("DGETRF INFO = ", I0)') la_info
          error stop
        end if
      else
        call tsol(arr,tl,ld,m1,m2,m,n,ipiv,pyp,del) ! del :=
      end if
      nbsol = nbsol + 1

      if (firststep) then
        ! first step
        do i = 1, n
            pyp(i) = 0.0_wp
        end do
      else
          ! CALCULATIONS FOR LOCAL ERROR, omitted if first step
          if (rat >= 0.5_wp) then
            con = rat*(-b + alpha*(1 - alpha)) / (1 + 2*alpha*(rat - 1))
          else
            con = b*rat/(1 + rat)
          end if
          do i = 1, n
            pyp(i) = con*(del(i) - v1(i)*rat)
          end do
      end if
!
! calculate maximum relative local truncation error estimates
!
      rel = maxval(abs(pyp(1:n) + (alpha - 0.5_wp)*del(1:n)) / &
                   (abs(pyc(1:n)) + ymin)) / emax
!
! if rel > 1  halve step size
!
      halve_step = rel > 1.0_wp
      if (halve_step) exit attempt_step
!
! succesful step test for last step
!
      if (laststep) then
      ! last step, normal exit
!
! successful exit  return values
!  put results into yn array for output
!
        xout = x
        do i = 1, n
          yn(i) = pyc(i)
        end do
        ifail = 0
        return
      end if

      rat = 1.0_wp
!
! decide from value of rel whether step size should remain
! unchanged or be doubled
!
      idoha = 0
      if (rel <= 0.2_wp) then
        rat = 2.0_wp
        idoha = 1
        ijac = 0
      end if
!
! test if next step will be last step
!
      if (abs(xout-x)/(h*rat) <= 1.1_wp) then
        laststep = .true.
        ijac = 0
        xstep = h*rat
        rat = (xout - x)/h
      end if
      h = h*rat
      firststep = .false.

      x = x + h
      con = 1 + alpha*(rat - 1)
!
! reset arrays for next step
!
      do i = 1, n
        v1(i) = rat*del(i)
        yold1(i) = pyc(i) - yold(i)
        pyp(i) = rat*yold(i) + pyc(i) + con*v1(i)
        yold(i) = pyc(i)
        yn(i) = yn1(i)
        fn(i) = rat*fn2(i)
      end do
      call pmult(pyc,n,yn1) ! yn1 := B pyc
!
! output to routine monit if requried
!
      if (imon > 0) then
        mon = mon + 1
        ! output only every imon steps
        if (imon == mon) then
          mon = 0
          t = x - h
          do i = 1, n
            yn(i) = pyc(i)
          end do
          call monit(yn,n,t,nhalf,rat)
        end if
      end if
      nhalf = 0
!
! calculate neq predictor
!
      do i = 1, n
        pyp(i) = (1 + rat)*pyc(i) - rat*yold(i) + con*v1(i)
      end do


    end block attempt_step

    if (halve_step) then
!
! step rejected
!
      nhalf = nhalf + 1
!
! if step size has been doubles to h (idoha set 1) and the
! next step attempts to halve again new step set 3*h/4
!
      if (idoha == 1) then
        rat1 = 0.75_wp
        idoha = 0
      else
        rat1 = 0.5_wp
      end if
!
! maximum number of succesive halvings allowed is 20
!
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
        nfev = nfev + 1
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

end subroutine