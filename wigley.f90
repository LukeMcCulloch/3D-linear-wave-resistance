!!
!! program to discretize a Wigley hull and adjacent free surface
!! into quadrilateral panels
!!
PROGRAM wigley

  ! load module to adjust float precision
  USE precise, ONLY: DEFAULTP
  USE vtkxmlmod
  USE sphere
  IMPLICIT NONE
  INTEGER, PARAMETER :: WP=DEFAULTP

  REAL(wp), PARAMETER :: zero=REAL(0.0, kind=wp)

  REAL(wp) :: L, B, T
  REAL(wp) :: Lfs, Bfs
  REAL(wp) :: zfspanels

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: panels 
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: points 

  INTEGER :: i, j, ipt, ipan, io, narg

  INTEGER :: nhl, nhg, nhf
  INTEGER :: nfsl, nfst

  INTEGER :: npanels, nfspanels, npoints 
  integer :: nhpt

  REAL(wp) :: x, y, yfs, z, dx, dy, dz, deltax

!!**********************************************************************
!!
  CHARACTER(len=24) :: outputfile  ! file name can be max. 24 char. long
 ! CHARACTER(len=24) :: datafile1  ! file name can be max. 24 char. long
 ! CHARACTER(len=24) :: datafile2  ! file name can be max. 24 char. long
  ! if you need longer file names increase the value of len, e.g. len=35
 ! CHARACTER(len=16) :: argstr  ! string for command line argument
!!
!!***********************************************************************


!!**********************************************************************
!!
  ! Count the number command line arguments (actually Fortran 2003)
  ! If your compiler does not support this Fortran 2003 feature look
  ! for the common extension iargc()
  ! narg = iargc()
  narg = command_argument_count()

  IF ( narg > 0 ) THEN
     ! If there are 2 (or more) arguments we will retrieve them
     ! We use the first argument as name for the output file to
     ! the program
     CALL get_command_argument(1, outputfile)

!!$     ! The second argument is used as the name of the output file
!!$     ! name
!!$     CALL get_command_argument(2, argstr)
!!$     ! Convert argstr into a REAL number (USE of internal files)
!!$     READ(argstr,*) diameter
!!$     WRITE(6,*) argstr
!!$
!!$     ! thrid argument is the number of panels lengthwise
!!$     CALL get_command_argument(3, argstr)
!!$     ! Convert argstr into an INTEGER number (USE of internal files)
!!$     READ(argstr,*) nl
!!$
!!$     ! thrid argument is the number of panels girthwise
!!$     CALL get_command_argument(4, argstr)
!!$     ! Convert argstr into an INTEGER number (USE of internal files)
!!$     READ(argstr,*) ng

  ELSE
     ! no arguments have been provided on the command 
     ! line. Report the error and print a simple help.
     WRITE(6,'(A)') ' Output file name  missing!'
     WRITE(6,'(A)') ' Usage:>  wigley  <outputfile>'! <diameter> <nl> <ng>'
     ! Of course <inputfile> and <outputfile> have to be replaced
     ! with file names and not typed literally
     STOP ! abort program
  ENDIF 
!!
!!***********************************************************************


  !! Definitions
  L = 1.0   ! length
  B = 0.1   ! beam
  T = 0.0625   ! draft


  ! number of hull panels and points (small numbers to test program)
  nhl = 32  ! lengthwise 
  nhg = 8   ! girthwise below WL  ex: for Fn 0.3 wavelength hull of length 1, I need # panels
  nhf = 1   ! girthwise above WL

  ! number of hull panels and points (small numbers to test program)
!!$  nhl = 16  ! lengthwise 
!!$  nhg = 8   ! girthwise below WL  ex: for Fn 0.3 wavelength hull of length 1, I need # panels
!!$  nhf = 1   ! girthwise above WL
  
  npanels = nhl*(nhg+nhf)
  npoints = (nhl+1)*(nhg+nhf+1)

  !! number of panels on free surface
  
  ! We discretize a free surface area of 4*L lengthwise
  ! That is larger than but we want to see the waves!

  Lfs = 4.*L

  nfsl = 4*nhl
  nfst = INT(0.75*nhl) ! may need to be adjusted

  nfspanels = nfsl*nfst
  npoints = npoints + (nfsl+1)*(nfst+1)

  ! position of free surface panels above z=0
  ! H. Raven suggests: about 0.8 of panel length
  !     zfspanels = 0.8 * L/ nhl
  zfspanels = 0.65*L / real(nhl)


  !! Allocate memory
  ALLOCATE( panels(4, npanels+nfspanels) )
  ALLOCATE( points(3, npoints) )

  !! Compute hull offsets
  !! We go from bow to stern and from keel to WL
  ipt = 1  ! point counter

  dz = T / REAL(nhg)
 
  DO j = 1, nhg + nhf + 1  ! girthwise (vertical) subdivision
     ! a true subdivision of the girth would be nice but too
     ! complicated for now.

     IF (j <= nhg+1) THEN 
        z = REAL(j-1)*dz - T ! waterlines
     ELSE 
        z = zfspanels/REAL(nhf) * REAL(j-nhg-1) 
     ENDIF 

     DO i = 1, nhl+1  ! lengthwise subdivision
        dx = L / REAL(nhl)
      
        x = 0.5*L - REAL(i-1)*dx ! stations

        points(:, ipt ) = (/ x, halfbreadth(L,B,T,x,z), z /)
        WRITE(6,'(I4,2X,3(G16.8,2X))') ipt, points(:, ipt)

        ipt = ipt + 1

     END DO 
  END DO 


 
  ! number of points on the hull
  nhpt = ipt-1



  !! Compute hull panels
  !! We go from bow to stern and from keel to WL
  ipan = 1  ! panel counter

  DO j = 1, nhg + nhf ! girthwise (vertical) 
     DO i = 1, nhl  ! lengthwise 
        
        panels(1,ipan) = i + (j-1)*(nhl+1)
        panels(2,ipan) = i + j*(nhl+1)
        panels(3,ipan) = i + 1 + j*(nhl+1)
        panels(4,ipan) = i + 1+ (j-1)*(nhl+1)
        WRITE(6,'(I4'': '',4(I5,2X))') ipan, panels(:, ipan)

        ipan = ipan+1

     END DO 
  END DO 

  !! Compute free surface offsets
  !! We go from bow to stern and from midships to max. width

  dx = Lfs/REAL(nfsl)
!!******************************* Improvising an x shift!*************************
  deltax = .8*dx
!!*******************************Improvising an x shift!**************************
  dy = 0.02 * L ! width of first free surface strip midships
  yfs = zero  !start at midships

  ! width of free surface
  Bfs = 0.
  DO j = 0, nfst-1
     Bfs = Bfs + 1.1**j
  END DO
  Bfs = dy*Bfs
  
  DO j = 1, nfst + 1 ! transverse

     DO i = 1, nfsl + 1 ! lengthwise

        x = L - REAL(i-1)*dx

        y = (Bfs-halfbreadth(L,B,T,x,zfspanels))/Bfs*yfs &
                + halfbreadth(L,B,T,x,zfspanels)

        points(:, ipt ) = (/ x, y, zfspanels /)
        WRITE(6,'(I4,2X,3(G16.8,2X))') ipt, points(:, ipt)

        ipt = ipt + 1

     END DO 

     yfs = yfs + dy*1.1**(j-1)  

  END DO 







!! Compute the free surface panels

  do j = 1, nfst     ! nfs transverse loop
     do i = 1, nfsl  ! nfs lengthwise loop
        

        panels(1,ipan) = nhpt + i + (j-1)*(nfsl+1)
        panels(2,ipan) = nhpt + i + j*(nfsl+1)
        panels(3,ipan) = nhpt + i + 1 + j*(nfsl+1)
        panels(4,ipan) = nhpt + i + 1 + (j-1)*(nfsl+1)

        WRITE(6,'(I4'': '',4(I5,2X))') ipan, panels(:, ipan)
        
        ipan = ipan+1

        end do

  end do


  io = 10
  open(unit=io,file='wigley.pan')

  write(io,'(A)') 'Wigley hull'
  write(io,'(A)') '0  1'
  write(io,'(4(I5,2X))') npanels, nfspanels, npoints
  
  do i = 1, npanels+nfspanels
     write(io,'(4(I5,2X))') panels(:,i)
  end do

  do i = 1, npoints
     write(io,'(3(G16.8,2X))') points(:, i)
  end do

  close(io)




!!
!!***********************************************************************
!!

  ! Finally we save the hull + free surface discretization in the output file.
  ! The ampersand & can be used to continue a line.
  CALL savesphere( outputfile, 'A Wigley discretization', (/0,1/), &
                   npanels, nfspanels, npoints, deltax, panels, points )


!!
!!**********************************************************************
!!







 ! iounit = 10
  CALL VtkXmlPolyDataCellScalar( io, 'wig.vtp', npoints, &
        0, 0, 0, npanels+nfspanels, &
        RESHAPE(points, (/3*npoints/)), &
        RESHAPE((panels-1), (/4*(npanels+nfspanels)/)), &
        scalar1=real(panels(1,:),kind=wp), &
        namescalar1='x')




  ! program ends here


  DEALLOCATE( points )
  DEALLOCATE( panels )


CONTAINS

  FUNCTION halfbreadth( L, B, T, x, z ) RESULT (y) 

    !! Function to compute Wigley hull halfbreadth
    !! Standard Wigley hull has 
    !!   B/L = 0.1
    !!   T/L = 0.0625
    !!   CB = 0.444, CP = 0.667, S/L^2=0.1487

    REAL(wp) :: y
    REAL(wp) :: L, B, T
    REAL(wp) :: x,z

    !! The following simple quadratic equation refers to
    !! L. Landweber (1979). Wigley parabolic hull group discussion
    !! Proceedings of the workshop on ship wave-resistance computations.
    !! DTNSRDC, p.51, Bethesda, MY.
    !! 
    !! Origin is midships in the waterline. z-axis is pointing upwards,
    !! x-axis pointing forward, and y-axis points to port.
    !!

    IF ((ABS(x) > 0.5*L).OR.(z < -T)) THEN ! return zero if outside of hull
       y = 0.                              ! limits
    ELSE IF (z > 0.) THEN ! return maximum beam above WL
       y = 0.5*B*(1-(2.*x/L)**2)
    ELSE ! actual halfbreadth
       y = 0.5*B*((1.-(z/T)**2)*(1-(2.*x/L)**2))
    ENDIF 

   END FUNCTION halfbreadth

END PROGRAM wigley
