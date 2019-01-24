!!
!! Program to compute geometric properties of a quadrilateral panel.
!! 100226 lb
!!
!!
PROGRAM fifth    ! Every Fortran program starts with a program statement

  !! modules are program units that themselves contain variable
  !! declarations and subroutines and functions.
  !! The contents of the module is made availbe in other program units
  !! with the USE statement

  USE precise, ONLY : defaultp
  USE fifthpanel

!! variable declarations

  IMPLICIT NONE    ! makes sure the compiler complains about 
                   ! undeclared variables
  INTEGER, PARAMETER :: WP=defaultp
  
  INTEGER :: i, j, k

  INTEGER, PARAMETER :: npanels=1

  ! An array for panel corner coordinates
  INTEGER, DIMENSION(4,npanels) :: panels

  ! An array for point coordinates
  REAL(WP), DIMENSION(3,4) :: points
  REAL(WP), DIMENSION(3,4) :: corners  ! temporary storage for panel corners

  REAL(wp), DIMENSION(3,npanels) :: center
  REAL(wp), DIMENSION(3,3,npanels) :: coordsys
  REAL(wp), DIMENSION(2,4,npanels) :: cornerslocal
  REAL(wp), DIMENSION(npanels) :: area

  ! Variable declarations go before any executable statement

!! the actual program starts here

  ! Data for a test panel
  panels = RESHAPE( (/ 1, 2, 3, 4 /), (/4,1/))
  points = RESHAPE( (/ -1., -1., 0., &
                        1., -1., 0., &
                        1.,  1., 0., &
                       -1.,  1., 0. /), &
                        (/3,4/) )

  ! Simple start up message
  ! The character A is used for formating character constants and 
  ! variables 
  WRITE(6,'(A)') 'fifth, v 1.0, 100307, lb'
  WRITE(6,*)

  ! process geometry for each panel (here it is only one)
  DO i = 1, npanels
     WRITE(6,'(A)') 'Panel point numbers:'
     WRITE(6,'(4I6)') (panels(j,i), j=1,4)

     WRITE(6,'(A)') 'Panel point coordinates:'
     DO j = 1, 4
        WRITE(6,'(I5'': ''3(G12.6,2X))') panels(j,i), &
             (points(k,j), k = 1,3)
     END DO
     
     ! we store the panel corner points in a 3 by 4 array
     DO j=1, 4
        corners(:,j) = points(:,panels(j,i))
     END DO
 
     ! We call the panel geometry function in the module fifthpanel
     ! Output data is collected in arrays with additional dimension
     ! for number of panels
     CALL panelgeometry( corners, coordsys(:,:,i), cornerslocal(:,:,i), &
                      area(i), center(:,i) )

     ! We pick some output to check.
     WRITE(6,'(''nv: '',3(G12.6,2X))') coordsys(:,3,i)
     WRITE(6,'(''cl: '',8(G12.6,2X))') cornerslocal(:,:,i)
     WRITE(6,'(''A: '',G12.6)') area(i)

  END DO 
! Every program ends with end program
END PROGRAM fifth
! end of file

