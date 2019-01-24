!!
!! Simple program to generate panels for a quarter sphere
!! The sphere is also stored in a VTK file for Paraview
!! 100212, lb
!!
!!
PROGRAM fourth    ! Every Fortran program starts with a program statement

  !! modules are program units that themselves contain variable
  !! declarations and subroutines and functions.
  !! The contents of the module is made availbe in other program units
  !! with the USE statement

  use precise, only : defaultp
  USE fourthconst  ! incorporates everything from the module
  !USE thirdconst, ONLY : pi ! would only incorporate pi and ignore the rest
  USE sphere ! a second module
  use vtkxmlmod

!! variable declarations

  IMPLICIT NONE    ! makes sure the compiler complains about 
                   ! undeclared variables
  integer, parameter :: WP=defaultp
  
  INTEGER :: narg      ! Integer variables
  INTEGER :: i, j, n
  INTEGER :: ipt, ipn
  INTEGER :: nl, ng   ! controls for panel distribution
  INTEGER :: iounit
  
  INTEGER :: npoints, npanels

  ! An allocatable array. Its actual size will be defined during
  ! execution of the program.
  ! The array has two dimensions (:,:), i.e. number of comma separated 
  ! colons determines dimensions.
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: panels

  CHARACTER(len=24) :: outputfile  ! file name can be max. 24 char. long
  CHARACTER(len=24) :: datafile1  ! file name can be max. 24 char. long
  CHARACTER(len=24) :: datafile2  ! file name can be max. 24 char. long
  ! if you need longer file names increase the value of len, e.g. len=35
  CHARACTER(len=16) :: argstr  ! string for command line argument


  REAL(WP) :: diameter
  REAL(WP) :: theta, phi
  ! An allocatable array for real values
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: points

  ! an array for data we make up
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: celldata1
 
  ! Variable declarations go before any executable statement

!! the actual program starts here

  ! Simple start up message
  ! The character A is used for formating character constants and 
  ! variables 
  WRITE(6,'(A)') 'fourth, v 1.0, 100212, lb'

  ! Count the number command line arguments (actually Fortran 2003)
  ! If your compiler does not support this Fortran 2003 feature look
  ! for the common extension iargc()
  ! narg = iargc()
  narg = command_argument_count()

  IF ( narg > 3 ) THEN
     ! If there are 2 (or more) arguments we will retrieve them
     ! We use the first argument as name for the output file to
     ! the program
     CALL get_command_argument(1, outputfile)

     ! The second argument is used as the name of the output file
     ! name
     CALL get_command_argument(2, argstr)
     ! Convert argstr into a REAL number (USE of internal files)
     READ(argstr,*) diameter
     WRITE(6,*) argstr

     ! thrid argument is the number of panels lengthwise
     CALL get_command_argument(3, argstr)
     ! Convert argstr into an INTEGER number (USE of internal files)
     READ(argstr,*) nl

     ! thrid argument is the number of panels girthwise
     CALL get_command_argument(4, argstr)
     ! Convert argstr into an INTEGER number (USE of internal files)
     READ(argstr,*) ng

  ELSE
     ! Only 1 or no arguments have been provided on the command 
     ! line. Report the error and print a simple help.
     WRITE(6,'(A)') ' Output file name and/or sphere diameter missing!'
     WRITE(6,'(A)') ' Usage:>  third <outputfile> <diameter> <nl> <ng>'
     ! Of course <inputfile> and <outputfile> have to be replaced
     ! with file names and not typed literally
     STOP ! abort program
  ENDIF 
  
  WRITE(6,*) outputfile
  WRITE(6,*) diameter
  WRITE(6,*) nl
  WRITE(6,*) ng
  WRITE(6,'(G20.14)') pi

  
  ! We first determine a set of points on the sphere. The points
  ! are organized in ng+1 rows with n+1 points each.

  ! Allocate memory to store point x,y,z coordinates
  npoints = (nl+1)*(ng+1)  ! number of points (points per WL x number
                           ! of WLs)
  ! this statement reserves memery space (if available) for the two
  ! -dimensional  array points. First dimension is 3 (for x,y,z) and
  ! second dimension counts number of points
  ALLOCATE( points(3,npoints) )

  ! We define a loop over the waterlines (starting at the South pole 
  ! of the sphere)
  DO j = 1, ng+1
     ! the angle of latidude (-pi/2 up to 0 (WL))
     !MODIFIED********************************************************
      !theta = 0.5*pi*((j-1)/(0.5*REAL(ng)) -1.)!Whole sphere
      theta = 0.5*pi*((j-1)/REAL(ng) -1.)  ! original

     ! The loop over stations (points on waterline)
     DO i = 1, nl+1
        ! the angle og longitude (0 (bow) up to pi (stern))
        !Modified*****************************************************
        ! phi = 2.0*(i-1)*pi / REAL(nl)   ! whole sphere
         phi = 1.0*(i-1)*pi / REAL(nl) !1/4 sphere (original)

        ! the index of the point
        ipt = i + (j-1)*(nl+1)
        ! point coordinates. Note that the function spherept (module
        ! sphere) returns an array!
        points(:,ipt) = spherept( 0.5*diameter, theta, phi )

        WRITE(6,*) ipt, points(:,ipt)
     END DO ! end of loop over stations
  END DO ! end of loop over waterlines

  ! We now organize the points into panels. Panels are small
  ! quadrilateral elements of the body surface. Each panel is
  ! defined by its four corner points. We will only store the
  ! number of the corner points for each panel. The actual coordinates
  ! of the points are stored in points. Such a table is also
  ! known as an incidence table.

  ! Allocate memory for the incidence table
  npanels = nl*ng    ! number of panels
  ALLOCATE( panels(4,npanels) ) ! 4 corners per panel
  
  ! Loop over waterlines
  DO j = 1, ng
     ! loop over stations
     DO i = 1, nl
        ! panel number
        ipn = i + (j-1)*nl

        ! point numbers are stored in clockwise sequence if the
        ! panel is looked at from the flow region (from outside of 
        ! the body). This ensures that the normal vector computed 
        ! later will point out of the flow (inside the body).
        panels(1,ipn) = i + (j-1)*(nl+1)
        panels(2,ipn) = i + (j)*(nl+1)
        panels(3,ipn) = i + 1 + (j)*(nl+1)
        panels(4,ipn) = i + 1 + (j-1)*(nl+1)
        WRITE(6,'(4I6)') panels(:,ipn)
     END DO
  END DO 

  ! Finally we save the sphere discretization in the output file.
  ! The ampersand & can be used to continue a line.
  CALL savesphere( outputfile, 'A sphere discretization', (/1,1/), &
                   npanels, npoints, panels, points )


  ! Store data in a VTK files

  ! First we store data valid for a panel (cell in VTK jargon)

  ! Create a file name with the appropriate extension
  i = INDEX(outputfile, ".", BACK=.TRUE.)
  IF (i == 0) THEN
     WRITE(datafile1, '(A,A)')  TRIM(outputfile), '.vtp'
  ELSE
     IF (i > 20) i = 20 ! To make sure the filename does not exceed 24 char
     WRITE(datafile1, '(A,A)')  outputfile(1:i-1), '.vtp'
  ENDIF 
  PRINT*, datafile1
  
  ! In lieu of real data we create an array which we 
  ! fill with the average x component of the panel corners
  ALLOCATE(celldata1(npanels))
  
  DO i = 1, npanels
     celldata1(i) = 0.25*( points(1,panels(1,i)) &
                             + points(1,panels(2,i)) &
                               + points(1,panels(3,i)) &
                                 + points(1,panels(4,i)) )
  END DO 
  !
  ! Note that we have to reduce the point numbers by 1 (panels-1) 
  ! as VTK starts counting from zero.
  !
  iounit = 10
  CALL VtkXmlPolyDataCellScalar( iounit, datafile1, npoints, &
        0, 0, 0, npanels, &
        RESHAPE(points, (/3*npoints/)), &
        RESHAPE((panels-1), (/4*npanels/)), &
        scalar1=celldata1, &
        namescalar1='Average x')

  ! Now we store data valid for each corner point of the panels

  ! Create a file name with the appropriate extension
  i = INDEX(outputfile, ".", BACK=.TRUE.)
  IF (i == 0) THEN
     WRITE(datafile2, '(A,A,A)')  TRIM(outputfile), 'pt', '.vtp'
  ELSE
     IF (i > 18) i = 18 ! To make sure the filename does not exceed 24 char
     WRITE(datafile2, '(A,A,A)')  outputfile(1:i-1), 'pt', '.vtp'
  ENDIF 
  PRINT*, datafile2
  
  ! In lieu of real data we just store the z-component of the
  ! point coordinates
  !
  ! Note that we have to reduce the point numbers by 1 (panels-1) 
  ! as VTK starts counting from zero.
  !
  iounit = 10
  CALL VtkXmlPolyDataPtScalar( iounit, datafile2, npoints, &
        0, 0, 0, npanels, &
        RESHAPE(points, (/3*npoints/)), &
        RESHAPE((panels-1), (/4*npanels/)), &
        scalar1=points(3,:), &
        namescalar1='z')

  ! It is good practice to deallocate arrays before the program ends
  ! or when you do not need them anymore
  IF (ALLOCATED(celldata1)) DEALLOCATE(celldata1)
  IF (ALLOCATED(points)) DEALLOCATE(points)
  IF (ALLOCATED(panels)) DEALLOCATE(panels)

 
  ! Every program ends with end program
END PROGRAM fourth
! end of file

