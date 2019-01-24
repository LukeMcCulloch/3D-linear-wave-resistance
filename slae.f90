MODULE SLAE
!
! Module for the solution of systems of linear equations.
! Last update: 09.01.2001 by Lothar Birk
!  
  USE precise, ONLY: DEFAULTP
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: WP=DEFAULTP


CONTAINS


  SUBROUTINE gauss(A,B,X,ERROR)

! Gauss elimination with partial Pivot element search (row switching)
! Solves the system of linear equations
!     Ax = b
! A and b are changed by the routine!
!
  REAL(WP), INTENT(IN OUT) :: A(:,:)
  REAL(WP), INTENT(IN OUT) :: B(:)

  INTEGER, INTENT(IN OUT) :: ERROR
  REAL(WP), INTENT(OUT) :: X(:)
  REAL(WP) :: am, P
  INTEGER :: I, J, K, ITMP, NEQ, l
  INTEGER, ALLOCATABLE :: IDX(:)

  ! Initialize some local variables
  ERROR = 0
  NEQ = SIZE(B)
  ALLOCATE(IDX(NEQ))
  
  ! Initialize index array (the rows are not actually switched. Only
  ! the inices are interchanged
  DO i = 1, NEQ
     IDX(i) = i
  END DO

  ! Loop over the equations from first to next to last Eq
  DO K = 1, NEQ-1
     J = K
     ! Search for Pivot element (largest number in column below
     ! current position
     DO I = K+1, NEQ
        IF ( ABS(A(IDX(I),K)) > ABS(A(IDX(J),K))) J = I
     END DO
     ! Make row with Pivot element the current row
     ITMP = IDX(J)
     IDX(J) = IDX(K)
     IDX(K) = ITMP

     ! Return error condition if pivot element is too small
     IF (ABS(A(IDX(K),K)) < 0.00000000000001) THEN
        PRINT*, IDX(K), K, A(IDX(K),K)
        ERROR=-1
        RETURN
     ENDIF
     ! Gauss elimination procedure
     DO I = K+1, NEQ

        ! Divide rows by pivot element
        am = A(IDX(I),K) / A(IDX(K),K)

        ! Subtract pivot element row from rows below
        !DO J = K+1, NEQ
        DO J = 1, NEQ
           A(IDX(I),J) = A(IDX(I),J) - am*A(IDX(K),J)
        END DO
        B(IDX(I)) = B(IDX(I)) - am*B(IDX(K))
     !!! Poor man debugger; uncomment if you want to see the matrix changes
     !!do l = 1, neq
     !!   WRITE(6,'(5G12.6)') (A(IDX(l),J), J = 1, neq), B(IDX(l))
     !!END DO
     !!WRITE(6,*)
     END DO
  END DO

  ! Back substitution part
  ! Start with last unknown
  ! Walk your way back up from next to last to first row (unknown)
  DO K = NEQ, 1, -1
     X(K) = B(IDX(K)) 
     DO I = K+1, NEQ
        X(K) = X(K) - A(IDX(K),I)*X(I)
     END DO
     ! Compute unknown
     X(K) =  X(K) / A(IDX(K), K)
  END DO

  ! Done
  RETURN
  END SUBROUTINE gauss


  SUBROUTINE gauss_simple(A,B,X,error)
! Bare bone Gauss elimination
! Solves the system of linear equations
!     Ax = b
! A and b are changed by the routine!
! Whenever a diagonal element is zero this routine will fail to produce
! a solution

  REAL(WP), INTENT(IN OUT) :: A(:,:)
  REAL(WP), INTENT(IN OUT) :: B(:)
 
  INTEGER, INTENT(IN OUT) :: ERROR
  REAL(WP), INTENT(OUT) :: X(:)

  INTEGER :: i, k, J, N
  REAL(wp) :: md

  N = SIZE(B) 
  ERROR = 0
  
  ! Gauss elimination (creation of upper triangular matrix)
  DO k = 1, N-1
     IF (ABS(a(k,k)) < 0.0000001) THEN
        error = -1
        RETURN
     ENDIF 
     DO i = k+1, N
        md = a(i,k)/a(k,k)
        DO j = k+1, n
           a(i,j) = a(i,j) - md*a(k,j)
        END DO
        b(i) = b(i) - md*b(k)
     END DO 
     !!! Poor man debugger; uncomment if you want to see the matrix changes
     !!DO i = 1, n
     !!   WRITE(6,'(5G12.6)') (A(i,J), J = 1, n), B(I)
     !!END do
  END DO 
  
  ! Back substitution part
  DO k = n, 1, -1
     x(k) = b(k)
     DO I = K+1, N
        x(k) = x(k) - A(k,i)*X(I)
     END DO
     x(k) = x(k) / a(k,k)
  END DO 
  RETURN 
  END SUBROUTINE gauss_simple




  SUBROUTINE GSIT(A,B,EPS,OMEGA,MAXIT,X,ERROR)
!
! Last update: 14.07.1999 by Lothar Birk
!
! Gauss-Seidel over-relaxation algorithm for the iterative solution of
! systems of linear equations:
!                        A*X = B
!
! The matrix A is of size (NEQ,NEQ) and has to be diagonal dominant.
! B is the right hand side vector of size (NEQ).
! eps is a small real number specifying the required accuracy of
! the solution.
! OMEGA is the relaxation parameter (real). Must be in the open interval
! (0,2). If you do not know the problem dependend proper value, let
! omega=1.0 (no relaxation, pure Gauss-Seidel iteration).
! Start testing with omega=1.5, if you are of the adventurous kind.
! MaxIt is the maximum number of iterations allowd for solving the LGS.
! X contains the initial solution when the procedure is called and
! on return provides the solution vector X.
! Error returns the status of the subroutine:
!    error = -1 : Critical error condition. One of the main diagonal 
!                 elements of the matrix is (almost) zero. The
!                 system can not be solved with this routine.
!    error = -2 : The solution did not converge to the required accuracy.
!                 The parameter MaxIt iterations was reached.
!    error > 0 :  The system is successfully solved and ERROR iterations
!                 have been used for the solution.
!
  INTEGER, INTENT(IN) :: MAXIT
  REAL(WP), INTENT(IN) :: A(:,:)
  REAL(WP), INTENT(IN) :: B(:)
  REAL(WP), INTENT(IN) :: EPS,OMEGA
 
  INTEGER, INTENT(IN OUT) :: ERROR
  REAL(WP), INTENT(OUT) :: X(:)

  INTEGER :: I, J
  INTEGER :: NEQ
  REAL(WP) :: Xold(SIZE(X))

  ! Initialize variables
  ERROR = 0
  NEQ = SIZE(B)  ! Number of equations

  ! Simple error check. None of the elements on the main diagonal
  ! of A is allowed to be zero
  DO i = 1, NEQ
    IF (ABS(A(i,i)) < 0.0000001) THEN
      ERROR = -1
      RETURN
    END IF
  END DO

  DO J = 1, MAXIT   ! LOOP OVER ITERATIONS
    Xold = X  ! Store last solution for convergence test
    ! The following loop contains the Gauss-Seidel over relaxation algorithm
    DO I = 1, NEQ
      X(i) = (1.-omega)*X(i) &
                  - omega*SUM(A(i,1:i-1)*X(1:i-1))/A(i,i) &
                      - omega*SUM(A(i,i+1:NEQ)*X(i+1:NEQ))/A(i,i) &
                          + omega*B(i)/A(i,i)
    END DO 
    ! Check for convergence. If the absolute difference between
    ! the next to last (Xold) and the last iteration step is smaller
    ! than eps exit iteration loop
    IF (SQRT(SUM((X-Xold)**2))/SQRT(SUM(Xold**2)) < eps ) THEN
      error = j
      EXIT 
    END IF
    ! This line should not be reached if the matrix is properly
    ! scaled
    error = -2  ! Return value for bad convergence, i.e. j=MaxIt
  END DO

  RETURN
  END SUBROUTINE gsit  


   SUBROUTINE SIMQIT(A,R,N,NMAX,SING,AV,DX,X,iter)
! from Soeding/Bertram
! Iterative solution of a non-singular, REAL linear system of equations
! A X = R, WHERE A, X and R are (n by n), (n by 1) and (n by 1) matrices,
! respectively. (n by c) means: n rows, c columns.

! Method: Incomplete Gauss elimination: Those matrix coefficients 
! the absolute 
! value of which exceeds the column pivot element times the value 
! AVERH are 
! eliminated before a Gauss-Seidel iteration WITH column pivoting 
! is applied. 
! IF no convergence is encountered, another elimination step 
! is performed
! WITH half the previous value of AVERH. The starting value of 
! AVERH is AV.

! PARAMETER A (input; changed): Two-index array of size 1:NMAX >= n,
! 1: >= n + 2) containing the elements of the matrix A. FIRST INDEX 
! COUNTS COLUMNS, LAST INDEX COUNTS ROWS.
! PARAMETER R (input; changed): One-index array containing the elements 
!  of the   matrix R
! PARAMETER N (input): Number of unknowns
! PARAMETER NMAX (input): Maximum defined value of first index of A
! PARAMETER SING (input): Lower limit of absolute value of pivot
! elements. IF
!     a division by an element the absolute value of which is smaller 
!     than SING
!     would be required, elimination is stopped WITH an error message 
!     indicating
!     the indices of those equations which are responsible for the 
!     difficulty.
!     It is recommended to generate the equation system such that 
!     maximum 
!     absolute values of the elements of A in each row and column DO
!      not differ 
!     by orders of magnitude. THEN, for SING a value of 10**-5 times 
!     that 
!     maximum absolute value may be recommended.
! PARAMETER AV (input): See above under method. Recommended value: 
!     SQRT(1/N).
! PARAMETER DX: The iteration stops IF the maximum absolute change of 
!     any 
!     unknown within one iteration step is smaller than DX.
! PARAMETER X (input and RESULT): array of size 1: >= n containing 
!     an estimated
!     solution when calling the SUBROUTINE (IF no estimation is 
!     available, 0 is
!     recommended), and the solution vector X after RETURN
! Parameter iter (output): integer, number of iterations needed
!

   USE precise, ONLY: DEFAULTP
   IMPLICIT NONE
   INTEGER, PARAMETER :: WP=DEFAULTP

   INTEGER IS, ITER, IZ, IZPIV, IZ1, N, NMAX
   REAL(wp) AA,AV, AVERH, ALIMIT, BMAX, DX, DXMAX, DXMAXV, SAVEAIJ, SING
   !REAL(wp) A(NMAX,N+2), R(N), X(N)
   REAL(wp), DIMENSION(:,:) :: A
   REAL(wp), DIMENSION(:) :: R, X
! Initial values for R and (n+2)nd column of A (equation numbers)
   DO IS=1,N
      A(IS,N+2)=IS+0.5
   END DO 
   DO IZ=1,N
      DO IS=1,N
         R(IZ)=R(IZ)-X(IS)*A(IS,IZ)
      END DO 
   END DO 
   AVERH=AV*1.4
   iter=1
3 CONTINUE
      AVERH=AVERH/1.4
      iter=iter+1
      !PRINT*, "simqit: iter=", iter
!     For all rows
      !!PRINT*, A(N,:)
      DO IZ = 1,N
!     Search for pivot element 
         !PRINT*, IZ
	 BMAX=0.
	 DO  IZ1 = IZ,N
	    IF(ABS(A(IZ,IZ1)).GT.BMAX)THEN
	       IZPIV = IZ1
	       BMAX = ABS(A(IZ,IZ1))
	    ENDIF
         END DO 
         !!PRINT*, BMAX
	 IF (BMAX.LT.SING) THEN
	    WRITE(*,*)'Equation system is singular: Bmax=',BMAX
	    STOP
	 ENDIF
!        Exchange of pivot and IZ rows; normalization of pivot row
	 ALIMIT = BMAX*AVERH
	 BMAX=SIGN(BMAX,A(IZ,IZPIV))
	 IF (IZPIV.EQ.IZ) THEN
	    DO IS=1,N
               A(IS,IZ)=A(IS,IZ)/BMAX
            END DO 
	 ELSE
	    DO IS = 1,N
               SAVEAIJ=A(IS,IZPIV)
	       A(IS,IZPIV)=A(IS,IZ)
               A(IS,IZ)=SAVEAIJ/BMAX
            END DO 
	 ENDIF
	 SAVEAIJ=A(IZPIV,N+2)
	 A(IZPIV,N+2)=A(IZ,N+2)
	 A(IZ,N+2)=SAVEAIJ
	 SAVEAIJ=R(IZPIV)
	 R(IZPIV)=R(IZ)
	 R(IZ)=SAVEAIJ/BMAX
!     Elimination of large coefficients
	 DO IZ1=IZ+1,N
	    AA=-A(IZ,IZ1)
	    IF (ABS(AA).GT.ALIMIT) THEN
	       DO IS=1,N
                  A(IS,IZ1)=A(IS,IZ1)+AA*A(IS,IZ)
               END DO 
	       R(IZ1)   =R(IZ1)   +AA*R(IZ)
	    ENDIF
         END DO 
      END DO 
!     Determination of changes of X stored in A(IZ,N+1)
      DXMAXV=1.E30
60    CONTINUE 
      DXMAX=0
      DO IZ=N,1,-1
	 A(IZ,N+1)=R(IZ)
	 DO IS=IZ+1,N
            A(IZ,N+1)=A(IZ,N+1)-A(IS,IZ)*A(IS,N+1)
         END DO 
	 DXMAX=MAX(DXMAX,ABS(A(IZ,N+1)))
         X(IZ)=X(IZ)+A(IZ,N+1)
      END DO 
      IF(DXMAX.LT.DX) RETURN
!     Update of absolute value vector R
      DO IZ=1,N
	 R(IZ)=0.
	 DO IS=1,IZ-1
            R(IZ)=R(IZ)-A(IS,N+1)*A(IS,IZ)
         END DO 
      END DO 
      IF(DXMAX.GT.0.8*DXMAXV) GOTO 3
      DXMAXV=DXMAX
      GOTO 60
   END SUBROUTINE simqit



END MODULE SLAE 
