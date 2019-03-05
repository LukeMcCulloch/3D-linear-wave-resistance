!!
!! Luke McCulloch
!! 6160 Hw 1
!! 2-21-10
!! Simple program uses subroutines to read input file (fifi.dat)
!! and store output in a seperate file  (out.dat)
!! The input file is fifi
!!
!!
!!
PROGRAM hw1    ! Every Fortran program starts with a program statement

  !! modules are program units that themselves contain variable
  !! declarations and subroutines and functions.
  !! The contents of the module is made availbe in other program units
  !! with the USE statement

  use precise, only : defaultp
  use io
!  USE fourthconst  ! incorporates everything from the module
  !USE thirdconst, ONLY : pi ! would only incorporate pi and ignore the rest
!  USE sphere ! a second module
!  use vtkxmlmod

!! variable declarations

  IMPLICIT NONE    ! makes sure the compiler complains about 
                   ! undeclared variables
  integer, parameter :: WP=defaultp
  
  INTEGER :: narg      ! Integer variables
  INTEGER :: i, j, n
  INTEGER :: npoints, npanels

  ! An allocatable array. Its actual size will be defined during
  ! execution of the program.
  ! The array has two dimensions (:,:), i.e. number of comma separated 
  ! colons determines dimensions.
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: panels        !-----------------------------------------------

  CHARACTER(len=24) :: inputfile   ! file name can be max. 24 char. long
  CHARACTER(len=24) :: outputfile  ! file name can be max. 24 char. long
  CHARACTER(len=30) :: title


  LOGICAL :: flexists  ! a logical variable. I .true. or .false.


  ! An allocatable array for real values
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: points        !-----------------------------------------------
  
  real(wp), allocatable, dimension(:,:) :: op
!  real(wp), allocatable, dimension(:,:) :: opi
  
 
  ! Variable declarations go before any executable statement

!! the actual program starts here

  ! Simple start up message
  ! The character A is used for formating character constants and
  ! variables
  WRITE(6,'(A)') 'HW1, v @.0, 2-20-10, Main Program, Luke McCulloch '

  ! Count the number command line arguments (actually Fortran 2003)
  ! If your compiler does not support this Fortran 2003 feature look
  ! for the common extension iargc()
  ! narg = iargc()
  narg = command_argument_count()
  WRITE(6, FMT='(AI3A)') 'we have ', narg, ' command line arguments '
  
  ! Now read the name of the input file and output file from the command line

  IF ( narg > 1 ) THEN
     ! If there are 2 (or more) arguments we will retrieve them

     ! We use the first argument as name for the input file to
     ! the program
     CALL get_command_argument(1, inputfile)

     ! We use the second argument as name for the output file to
     ! the program
     CALL get_command_argument(2, outputfile)

  ELSE
     ! Only 1 or no arguments have been provided on the command
     ! line. Report the error and print a simple help.
     WRITE(6,'(A)') ' Output file name and/or sphere diameter missing!'
     WRITE(6,'(A)') ' Usage:>  third <outputfile> <diameter> <nl> <ng>'
     ! Of course <inputfile> and <outputfile> have to be replaced
     ! with file names and not typed literally
     STOP ! abort program
  ENDIF

  INQUIRE(file=inputfile, exist=flexists)

  IF (flexists) THEN
  
     !Here are two different ways to write the output to the screen.
     WRITE(6,'(AA)') ' input file,  ', inputfile
     WRITE(6,*) 'output file, ', outputfile
  
!! Call a subroutine to read the input file back into the main program
!! and output the input file into the out.dat file as a check.
!!----------------------------------------------------------------------
     CALL input( inputfile, outputfile, title, &
                         npanels, npoints, &
                         panels, points,i,j)
!!----------------------------------------------------------------------

  Else
     ! If not, the program stops
     WRITE(6,'(A)') ' Input file does not exist!'
     STOP ! abort program
  ENDIF
  
  ! We read fifi in the "input" subroutine
  ! Write fifi to the screen
  ! write(6,'(A)') title
  write(6,*) title
  write(6,*) i,j
  write(6,*) npanels, npoints
  write(6,'(4I4)') panels
  write(6,'(3D16.8)') points

  ! Manipulate the data (the "interpret" part of the homework)
  ! Compute the "Matrix INNER Product"
  ! i.e. there are 3 column vectors.  Each is to be dotted with itself
  ! And every other column vector to form
  !a 3X3 matrix of inner products.
  !  This matrix will be called "op"

  ALLOCATE( op(3,3) )
  write(6,*)
  write(6,*) 'The main program computes the matrix inner product'
  write(6,*) 'Of points'
  
  do n=1,3
    do i=1,3
       do j=1, npoints
          op(i,n) = op(i,n)+points(i,j)*points(n,j)
       end do
    end do
  end do

  write(6,'(3D16.8)') op
  
  ! Easy to check only if FIFI contains a small set of points!

  
  

!! Call a subroutine to Send the matrix of dot products
!!  to the file "out.dat".
!!----------------------------------------------------------------------
  CALL output( outputfile, op )
!!----------------------------------------------------------------------

  ! De-allocate all arrays
  IF (ALLOCATED(points)) DEALLOCATE(points)
  IF (ALLOCATED(panels)) DEALLOCATE(panels)
  IF (ALLOCATED(op)) DEALLOCATE(op)
!  IF (ALLOCATED(opi)) DEALLOCATE(opi)

 
  ! Every program ends with end program
END PROGRAM hw1
! end of file

