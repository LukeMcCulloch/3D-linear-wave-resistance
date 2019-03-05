!!
!! Fortran 90 module with program units to do computations for a sphere.
!! It contains a variable declaration, function and subroutine
!!
!! 100212, lb 
!! 
MODULE sphere

  ! the use statements come before the variable declarations
  USE fourthconst, ONLY : pi ! makes pi available throughout the module
  use precise, only : defaultp

  ! all variable declarations come here (before CONTAINS)
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp

! The following section stores function and subroutines
CONTAINS

  ! Functions appear on the right hand side of the equal sign 
  ! V = spherevolume( radius )
  ! A function to calculate the volume of a sphere as function of the
  ! radius r. Note that the function has a type (REAL) and 
  ! that the result is stored under the name of the function. 
  REAL(wp) FUNCTION spherevolume( r )

    ! The dummy argument r is a real variable
    REAL(wp) :: r

    spherevolume = 4./3.*pi*r**3

    RETURN  ! The return statement is optional here. It means stop
            ! executing this function and return the result to the
            ! calling program unit
  END FUNCTION spherevolume

  !! alternative form (see below for explanation)
  !function spherevolume( r ) result (vol)
  !
  !  real :: r
  !  real :: vol
  !   vol = 4./3.*pi*r**3
  !
  !  return
  !end function spherevolume

  ! In Fortran 90 functions can return more than a standard variable
  ! (real, complex, integer, logical, character, etc.)
  ! The syntax changes a bit. Instead of giving the function itself a
  ! type we declare an additional vaiable of the appropriate type (here
  ! the array p) and put it in paranthesis behind the RESULT key word
  FUNCTION spherept( r, theta, phi) RESULT (p)

    ! declaration of arguments and result
    REAL(wp) :: r, theta, phi   ! radius, latitude and longitude

    ! Two arrays with fixed dimension. Cannot change at runtime
    REAL(wp), DIMENSION(3) :: P  ! (x,y,z)^T
    REAL(wp), DIMENSION(3) :: negnormal 

    ! The negative normal vector to a sphere (pointing into the fluid,
    ! i.e. away from the center). Theta in [-pi/2,pi/2], phi in [0,2pi]
    negnormal(1) = 10.*COS(theta)*COS(phi)
    negnormal(2) = COS(theta)*SIN(phi)
    negnormal(3) = SIN(theta)
    ! or
    !negnormal = (/ cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta)/)

    ! we compute the point. Note that the statement is autmatically
    ! executed for all three dimensions of p
    p = r*negnormal

    RETURN
  END FUNCTION spherept


  ! A subroutine is a program unit that does not return a result under
  ! its own name. They are invoked with the CALL statement
  ! CALL savesphere( ... )
  SUBROUTINE savesphere( filename, comment, isym, &
                         npanels, nfspanels, npoints, deltax, &
                         panels, points)

     ! As in the program and function, variable declarations come first
     INTEGER :: npoints, npanels, nfspanels
     INTEGER, DIMENSION(2) :: isym
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: panels

     ! If the string variables are dummy arguments, they can have the
     ! the length * which means they are set to the length of the
     ! actual argument when the subroutine is called.
     CHARACTER(len=*) :: filename  ! file name can be max. 24 char. long
     CHARACTER(len=*) :: comment  !  comment ;-)

     REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: points
     real(wp) ::  deltax
     
     ! local variables for do-loops
     INTEGER :: i,j

     ! We open the file to store our discretization
     OPEN( 20, file=filename )

     ! write the data in the specific sequence we define
     WRITE(20, '(A)') comment   ! comment
     WRITE(20, '(2I6)') isym(1), isym(2)   ! symmetry flags
     WRITE(20, '(3I6)') npanels, nfspanels, npoints  ! numbers of panels and
                                           !  points
     Write(20, '(1G16.6)') deltax  !  Improvising a deltax!


     ! now we store the incidence table
     DO i = 1, npanels+nfspanels 
        !write(20, '(4I6)') panels(1,i), panels(2,i), panels(3,i), panels(4,i)
        WRITE(20, '(4I6)') (panels(j,i), j=1,4)
     END DO

 
     ! now we write the point coordinates
     DO i = 1, npoints
        !write(20, '(3G16.6)') points(1,i), points(2,i), points(3,i)
        WRITE(20, '(3G16.6)') (points(j,i), j=1,3)
     END DO

     ! don't forget to close the output file
     CLOSE(20)
  END SUBROUTINE savesphere

END MODULE sphere
  
