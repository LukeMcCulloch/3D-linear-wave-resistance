!!
!! Luke McCulloch
!! p1.f90
!! 6160 CFDpart3
!! 4-17-10
!! Program will read in the hull and free surface data
!! Solve for the source strengths
!! Find the flow properties
!! and forces on the hull


PROGRAM p3

  !! modules are program units that 
  !! themselves contain variable declaerations
  !! and subroutines and functions.
  !! The contents of the module is made availbe 
  !! in other program units
  !! with the USE statement

  use precise, only : defaultp  ! Double Precision
  use io                        ! Input Output (read in fifi)
  use fifthpanel                ! Get the Panel Geometry vn, vt1, vt2
  use A                         ! Compute Dij coefficients in the "A" matrix
  use vtkxmlmod                 ! Create Paraview files
  use slae                      ! invert the Dij "A" matrix

!  USE fourthconst  ! incorporates everything from the module
!  USE thirdconst, ONLY : pi ! would only incorporate pi and ignore the rest
!  USE sphere ! a second module

!! variable declarations

  IMPLICIT NONE    ! makes sure the compiler complains about 
                   ! undeclared variables
  integer, parameter :: WP=defaultp
  
  INTEGER :: narg      ! Integer variables
  INTEGER :: i, j, k, error, Row, Col
  INTEGER :: npoints, npanels, nfspanels
  INTEGER :: iounit

  CHARACTER(len=24) :: inputfile   ! file name can be max. 24 char. long
  CHARACTER(len=24) :: outputfile  ! file name can be max. 24 char. long
  CHARACTER(len=24) :: datafile1
  CHARACTER(len=24) :: datafile2
  CHARACTER(len=24) :: datafile3
  CHARACTER(len=24) :: datafile4
  CHARACTER(len=24) :: datafile5

  CHARACTER(len=30) :: title

  LOGICAL :: flexists  ! a logical variable. I .true. or .false.

  ! An allocatable array. Its actual size will be defined during
  ! execution of the program.
  ! The array has two dimensions (:,:), i.e. number of comma separated 
  ! colons determines dimensions.
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: panels  
  ! allocated in the io module
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: points 
  !allocated in the io module


  ! Arrays for the fifth.f90 panel geometry property computations
  REAL(WP), DIMENSION(3,4) :: corners  ! temporary storage for panel corners
  REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: center
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: coordsys
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: cornerslocal
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: area
  
  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: vel
  REAL(WP), DIMENSION(3) :: vp
  REAL(WP), DIMENSION(3) :: vpp
  REAL(WP), DIMENSION(3) :: vppp
  REAL(WP), DIMENSION(3) :: fieldpoint, fieldpointp, fieldpointpp
  REAL(WP), DIMENSION(3) :: fieldpointppp


  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: am
  !!REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: asim  !! dont need for simquit?
  

  REAL(WP), ALLOCATABLE, DIMENSION(:) :: b
  REAL(WP), DIMENSION(3) :: vinf

  REAL(WP), ALLOCATABLE, DIMENSION(:) :: sigma
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: vtotal
  REAL(WP), ALLOCATABLE, DIMENSION(:) :: vn
  REAL(WP), ALLOCATABLE, DIMENSION(:) :: vt1
  REAL(WP), ALLOCATABLE, DIMENSION(:) :: vt2
  REAL(WP), ALLOCATABLE, DIMENSION(:) :: cp

  REAL(WP),DIMENSION(3) :: force

  REAL(WP), DIMENSION(3) :: vtottest
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: velt
  REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: vtest  
  REAL(WP), DIMENSION(3) :: testpoint
  REAL(WP), DIMENSION(3) :: vtp
  REAL(WP), DIMENSION(3) :: vtpp
  REAL(WP), DIMENSION(3) :: vtppp
  REAL(WP), DIMENSION(3) :: testpointp, testpointpp
  REAL(WP), DIMENSION(3) :: testpointppp  
  REAL(WP), DIMENSION(1) :: cpt
  !!REAL(WP) :: testi, testj, testk

  real(wp) :: deltax
  !!simquit stuff
  REAL(wp) :: Amatmax, sing, av, dx
  INTEGER :: nmax, iter
  !!end simquit stuff
  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:,:) :: hessglobal
  REAL(WP), DIMENSION(3,3) :: hp
  REAL(wp) :: cw, S, sourcesink
  REAL(WP), ALLOCATABLE,DIMENSION(:) :: zeta
  REAL(WP), DIMENSION(2,128) :: profilewave
  

!  REAL(wp), ALLOCATABLE, DIMENSION(:) :: celldata1
!  REAL(wp), ALLOCATABLE, DIMENSION(:) :: celldata2
!  REAL(wp), DIMENSION(3) :: paneloffsets

  ! Variable declarations go before any executable statement

!! the actual program starts here

  ! Simple start up message
  ! The character A is used for formating character constants and
  ! variables
  WRITE(6,'(A)') 'Linear Free Surface Panel Method, v @.0, 4-17-10, Main Program, Luke McCulloch '

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
!! This is from the "io" MODULE
!! ALLOCATE 2 ARRAYS : "points" and "panels"
!!
!! and output the input file into the out.dat file as a check.
!!----------------------------------------------------------------------
     CALL input( inputfile, outputfile, title, &
                         npanels, nfspanels, npoints, &
                         panels, points,i,j,deltax)
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
  write(6,*) npanels, nfspanels, npoints, deltax
  write(6,'(4I4)') panels
  write(6,'(3D16.8)') points

!!  Test - input the "test point" from the dos prompt.
!   Need to move this
!   Should be computed after source strengths are found
!
!!$Write(6,*) 'choose test point i value'
!!$read(*,*) testi
!!$Write(6,*) 'choose test point j value'
!!$read(*,*) testj
!!$Write(6,*) 'choose test point k  value'
!!$read(*,*) testk
!
!
!!
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!!  compute panel properties as in fifth.f90
  

 ALLOCATE( center(3,npanels+nfspanels) )
 ALLOCATE( coordsys(3,3,npanels+nfspanels) )
 ALLOCATE( cornerslocal(2,4,npanels+nfspanels) )
 ALLOCATE( area(npanels+nfspanels) )



  DO i = 1, npanels+nfspanels
 !    WRITE(6,'(A)') 'Panel point numbers:'
 !    WRITE(6,'(4I6)') (panels(j,i), j=1,4)

 !    WRITE(6,'(A)') 'Panel point coordinates:'
 !    DO j = 1, 4
 !       WRITE(6,'(I5'': ''3(G12.6,2X))') panels(j,i), &
 !            (points(k,panels(j,i)), k = 1,3)
 !    END DO
 !    
 !    ! we store the panel corner points in a 3 by 4 array
     DO j=1, 4
        corners(:,j) = points(:,panels(j,i))
     END DO
 
     ! We call the panel geometry function IN THE  MODULE "fifthpanel"
     ! Output data is collected in arrays with additional dimension
     ! for number of panels
     CALL panelgeometry( corners, coordsys(:,:,i), cornerslocal(:,:,i), &
                      area(i), center(:,i) )

     ! We pick some output to check.
     !WRITE(6,'(''nv: '',3(G12.6,2X))') coordsys(:,3,i)
     !WRITE(6,'(''cl: '',8(G12.6,2X))') cornerslocal(:,:,i)
     !WRITE(6,'(''A: '',G12.6)') area(i)
     !WRITE(6,'(''center: '',3(G12.6,2X))') center(:,i)     

  END DO 

 
!! end of fifth.f90 computations
!!--------------------------------------------------------------------

  

!! Compute the A matrix
!!---------------------------------------------------------------------
!!---------------------------------------------------------------------
    !!
    !! procedure to compute the A matrix of v dot a 
    !!
    !! Input:  the locations of panel centers pi and qj 
    !!         the normal vector of each panel center
    !!         the area of each panel
    !!         
    !!
    !!
    !! Output:  the velocity at panel i induced by panel j
    !!          the A matrix
    !!  
    !!
    !!
    !! Note we will use the following  from fifthpanel Output:
    !!     coordsys: a 3x3 matrix containing the two tangent vectors
    !!               and the normal vector
    !!                  t1x   t2x   nx
    !!                  t1y   t2y   ny
    !!                  t1z   t2z   nz
    !!               i.e. coordsys(:,3) is the normal vector
    !!     area: the area of the panel
    !!     center: the 3D coordinates of the panel center
    !!

  !! Allocate needed memory
  ALLOCATE(vel(3,npanels+nfspanels,npanels+nfspanels))
  ALLOCATE(hessglobal(3,3,npanels+nfspanels,npanels+nfspanels))
  ALLOCATE(am(npanels+nfspanels+2,npanels+nfspanels+2))
  am = REAL(0.,kind=wp)
  ALLOCATE(b(npanels+nfspanels))
  ALLOCATE(sigma(npanels+nfspanels))



  !!  Define Vinfinity
  vinf = (/1.503404137,0.,0./) !Corresponding to Fr=.350

  !!  Double Nested Double Do Loop to compute the A matrix 
  !!  PART 1:  Body B.C.
  !!**********************************************************************************************************************************************

    !! First get the image system of mirrored points  (enforce the centerplane!)
    !! in lieu of mirroring panels
    DO i=1,npanels
       
       Row=i

       fieldpoint(:)=center(:,Row)

!       !! First mirror the fieldpoints about the Z-X plane
       fieldpointp(:)=fieldpoint(:)
       fieldpointp(2)=-fieldpoint(2)

    
       !! Compute the Influence of the Body Panels on hull fieldpoints
       !!****************************************************************************
       Do j=1,npanels
          
          Col = j
!
!          !First get the singular points
          IF(Row==Col)THEN
             vel(:,Row,Col)=-0.5*coordsys(:,3,Row)
!
!          !Then get the nonsingular points
          ELSE
             vel(:,Row,Col)=hsinfluence(fieldpoint(:),center(:,Col),area(Col),coordsys(:,:,Col),cornerslocal(:,:,Col))
          END IF

          !Now get the influence of the panel 
          !on the mirrored points  (must enforce no flow through the centerplane of the ship!)
          vp(:)=hsinfluence(fieldpointp(:),center(:,Col),area(Col),coordsys(:,:,Col),cornerslocal(:,:,Col))
!

          !Get Influence of panel sj 
          !and its mirrors @ point pi

          !Use no mirrors if modeling the whole sphere.
          
          !Velocities for the half ship (flow is symetric)
          vel(1,Row,Col) = vel(1,Row,Col)+vp(1)
          vel(2,Row,Col) = vel(2,Row,Col)-vp(2)
          vel(3,Row,Col) = vel(3,Row,Col)+vp(3)
 


          !Finally, make the a(i,j) matrix components satisfying no flow through the body surface
          am(Row,Col)= DOT_PRODUCT(coordsys(:,3,Row),vel(:,Row,Col))

       END DO

       !! Compute the Influence of the Free Surface Panels on hull fieldpoints
       !! ***********************************************************************************
       Do j=1,nfspanels
          
          Col = j+npanels
!
!          !There are no singular points.  We put the sources above the free surface!
          IF(Row==Col)THEN
             vel(:,Row,Col)=-0.5*coordsys(:,3,Row)
!
!          ! Get the nonsingular points
          ELSE
             vel(:,Row,Col)=hsinfluence(fieldpoint(:),center(:,Col),area(Col),coordsys(:,:,Col),cornerslocal(:,:,Col))
          END IF

          !Now get the influence of the panel 
          !on the mirrored points
          vp(:)=hsinfluence(fieldpointp(:),center(:,Col),area(Col),coordsys(:,:,Col),cornerslocal(:,:,Col))
!

          !Get Influence of panel sj 
          !and its mirror @ point pi

          
          !Velocities for the half ship (flow is symetric)
          vel(1,Row,Col) = vel(1,Row,Col)+vp(1)
          vel(2,Row,Col) = vel(2,Row,Col)-vp(2)
          vel(3,Row,Col) = vel(3,Row,Col)+vp(3)
 


          !Finally, make the a(i,j) matrix components satisfying no flow through the body surface
          am(Row,Col)= DOT_PRODUCT(coordsys(:,3,Row),vel(:,Row,Col))

       END DO

    END DO


!!***********************************************************************************************************************************************
!!
!!  Double Nested Double Do Loop to compute the A matrix 
!!  PART 2:  Free Surface B.C.

    Do i=1,nfspanels

       !Get the field points on the first order free surface, include x offset 
     
       Row = npanels + i

       fieldpoint(:)=center(:,Row)

!      !! No mirrors, but need offsets in i and k
       fieldpoint(1)=fieldpoint(1)+deltax
       fieldpoint(3)=0.

       fieldpointp(:) = fieldpoint(:)
       fieldpointp(2) = -fieldpoint(2)

       !! Calculate the contribution of body panels to the free surface condition
       !! ***********************************************************************
       Do j=1,npanels
          
          Col=j

          !  Call what kind of influence function?  Use the Linear free surface condition
          vel(:,Row,Col)=hsinfluence(fieldpoint(:),center(:,Col),area(Col),coordsys(:,:,Col),cornerslocal(:,:,Col))
          vp(:)=hsinfluence(fieldpointp(:),center(:,Col),area(Col),coordsys(:,:,Col),cornerslocal(:,:,Col))

          !Velocities for the half ship (flow is symetric)
          vel(1,Row,Col) = vel(1,Row,Col)+vp(1)
          vel(2,Row,Col) = vel(2,Row,Col)-vp(2)
          vel(3,Row,Col) = vel(3,Row,Col)+vp(3)


          hessglobal(:,:,Row,Col)=phixxinfluence(fieldpoint(:),center(:,Col),area(Col),coordsys(:,:,Col),cornerslocal(:,:,Col))
          hp(:,:)=phixxinfluence(fieldpointp(:),center(:,Col),area(Col),coordsys(:,:,Col),cornerslocal(:,:,Col))

          hessglobal(1,1,Row,Col)=hessglobal(1,1,Row,Col)+hp(1,1)

          !Finally, make the a(i,j) matrix components satisfying the FS BC
          am(Row,Col)= (vinf(1)**2.)*hessglobal(1,1,Row,Col) + 9.81*vel(3,Row,Col)

       END DO          

       !! Calculate the contribution of free surface panels to the free surface condition
       !!********************************************************************************
       Do j=1,nfspanels
          
          Col=j+npanels

          !   Enforce PHIz 
          vel(:,Row,Col)=hsinfluence(fieldpoint(:),center(:,Col),area(Col),coordsys(:,:,Col),cornerslocal(:,:,Col))
          vp(:)=hsinfluence(fieldpointp(:),center(:,Col),area(Col),coordsys(:,:,Col),cornerslocal(:,:,Col))

             !Velocities for the half ship (flow is symetric)
          vel(1,Row,Col) = vel(1,Row,Col)+vp(1)
          vel(2,Row,Col) = vel(2,Row,Col)-vp(2)
          vel(3,Row,Col) = vel(3,Row,Col)+vp(3)
          
          ! Enforce PHIxx squared
          hessglobal(:,:,Row,Col)=phixxinfluence(fieldpoint(:),center(:,Col),area(Col),coordsys(:,:,Col),cornerslocal(:,:,Col))
          hp(:,:)=phixxinfluence(fieldpointp(:),center(:,Col),area(Col),coordsys(:,:,Col),cornerslocal(:,:,Col))

          hessglobal(1,1,Row,Col)=hessglobal(1,1,Row,Col)+hp(1,1)

          !Finally, make the a(i,j) matrix components satisfying the FS BC
          am(Row,Col)= (vinf(1)**2.)*hessglobal(1,1,Row,Col) + 9.81*vel(3,Row,Col)

       END DO

    END DO
!!************************************************************************************************





!  write(6,*) 'the velocity on each panel center'
!  write(6,'(3D16.8)') vel
!  write(6,*)''
!  write(6,*) 'the a matrix'
!  write(6,'(6D16.8)') am
!  write(6,*)''

!!  End of a(i,j) matrix computations
!!---------------------------------------------------------------------


!!---------------------------------------------------------------------
!!  Compute the RHS


  DO i=1, npanels

     b(i)= DOT_PRODUCT(coordsys(:,3,i),vinf)  !From MODULE "SLAE"

  END DO

  Do i=1, nfspanels
     row=i+npanels
     b(row)=0.
  end do

!!  RHS done
!!---------------------------------------------------------------------


!! Additional Control Parameters for the Simqit matrix inverter
  Amatmax = maxval(am)                  !returns the largest value in am
  sing = REAL(1.e-6,kind=wp)*Amatmax    !make the matrix singular if any value exceeds Amatmax
  av = sqrt(1./REAL(npanels+nfspanels))
  dx = 0.000001
  nmax = npanels+nfspanels+2


  !!  Invert the a matrix and solve for the Source Strength
  call SIMQIT(TRANSPOSE(am), b, npanels+nfspanels, nmax, &
       sing, av, dx, sigma, iter)

  write(6,*) 'simquit iter', iter

!Dealocate a and b matrices!

  IF (ALLOCATED(am)) DEALLOCATE(am)
  IF (ALLOCATED(b)) DEALLOCATE(b)

! Allocate space for flow data
  ALLOCATE(vtotal(3,npanels+nfspanels))
  ALLOCATE(vn(npanels+nfspanels))
  ALLOCATE(vt1(npanels+nfspanels))
  ALLOCATE(vt2(npanels+nfspanels))
  ALLOCATE(cp(npanels+nfspanels))
  ALLOCATE(zeta(npanels+nfspanels))
  !ALLOCATE(profilewave(nfspanels/32))

  !!Evaluate Flow Properties
  !!************************************************************
  Do i=1,npanels+nfspanels
     vtotal(:,i) = -vinf(:)
     Do j=1, npanels+nfspanels
        vtotal(:,i)=vtotal(:,i)+sigma(j)*vel(:,i,j)
      End Do

     !! Get the NORMAL VELOCITY on each panel:
     vn(i)=Dot_Product(coordsys(:,3,i),vtotal(:,i))
     !Should be 0 for all panels

     !write(6,*) vn

     !! Get the tangent components on a panel
     vt1(i)=dot_product(coordsys(:,1,i),vtotal(:,i))
     vt2(i)=dot_product(coordsys(:,2,i),vtotal(:,i))

     ! Get the pressure on each panel from the steady Bernoulli Eq
    ! cp(i)=1.-((sum(vtotal(:,i)**2))/(sum(vinf**2)))

  !write(6,*) 'the normal velocity on each panel'
  !write(6,'(1D16.8)') vn
!  write(6,*)''
!  write(6,*) 'cp', cp

  End Do

  Do i=1,npanels

     ! Get the pressure on each panel from the steady Bernoulli Eq
     cp(i)=1.-((sum(vtotal(:,i)**2))/(sum(vinf**2)))

  End do

  sourcesink=0.
  Do i=1,npanels+nfspanels
     sourcesink = sourcesink+sigma(i)*area(i)
  end do

  write(6,*) 'sourcesink sum (should be near 0.)',sourcesink
     
!!************************************************************
!! Compute the wave height!
  zeta=0.
  Do i=1,nfspanels
     !zeta=0.
!     vtotal(:,i) = -vinf(:)
     !Do j=1, nfspanels
        Row=i+npanels
       ! Col=j+npanels
        !zeta(Row) = zeta(Row)+ ((vinf(1)/9.81)*sigma(Col)*vel(1,Row,Col))     ! Free surface wave elevation
        !zeta(Row) = zeta(Row)+ ((vinf(1)/9.81)*vtotal(1,Row))   !4-26-10 old
        zeta(Row) = (vinf(1)/9.81)*(vtotal(1,Row) + vinf(1))
      !End Do
   End Do


!! Coompute the forces on the hull

     force(:)=(/0.,0.,0./)
     
     Do i=1,npanels  ! Hull summation only

        if(center(3,i)<0.)Then

        force(:)=force(:)+ area(i)*( .5*1025.*sum(vinf**2.)*cp(i)-1025.*9.81*center(3,i) )*coordsys(:,3,i)
        Else 

        End if

  end do
  !force(1)=2.*force(1)

!! Wave Resistance Coefficient

  S=0.
  Do i=1,npanels
     if(center(3,i)<0.)Then
        S=S+area(i)
        Else
        End if
     End do

  cw = -force(1)/(.5*1000.*(vinf(1)**2.)*S)



  write(6,*) 'cw', cw


!!$! Write a seperate loop to compute vtotals and cp's at field points specified by you!
!!$! 1st, an Example: For the point (5,0,0) find the velocity induced by the sphere.
!!$
!!$
!!$!! Allocate needed memory
!!$ALLOCATE(vtest(3,npanels))
!!$ALLOCATE(velt(3,npanels))
!!$
!!$       testpoint = (/ testi, testj, testk /)
!!$       !testpoint=center(:,1)-.001*coordsys(:,3,i)
!!$       !! First mirror the fieldpoints about the Z-X plane
!!$       testpointp(:)=testpoint(:)
!!$       testpointp(2)=-testpoint(2)
!!$
!!$       ! Now mirror original fieldpoints about the Y-X plane
!!$       testpointpp(:)=testpoint(:)
!!$       testpointpp(3)=-testpointpp(3)
!!$
!!$       ! Now mirror the first mirrored fieldpoints about the Y-X plane
!!$       testpointppp(:)=testpointp(:)
!!$       testpointppp(3)=-testpointppp(3)
!!$
!!$
!!$
!!$     vtottest=(/0.,0.,0./)
!!$     Do j=1, npanels
!!$        vtest(:,j) = panelinfluence(testpoint(:),center(:,j),area(j),coordsys(:,:,j),cornerslocal(:,:,j))
!!$
!!$        !Now get the influence of the panel 
!!$        !on the mirrored points
!!$        vtp(:)=panelinfluence(testpointp(:),center(:,j),area(j),coordsys(:,:,j),cornerslocal(:,:,j))
!!$        vtpp(:)=panelinfluence(testpointpp(:),center(:,j),area(j),coordsys(:,:,j),cornerslocal(:,:,j))
!!$        vtppp(:)=panelinfluence(testpointppp(:),center(:,j),area(j),coordsys(:,:,j),cornerslocal(:,:,j))
!!$
!!$        !Velocities for the quarter body (in "double body" flow)
!!$        velt(1,j) = vtest(1,j)+vtp(1)+vtpp(1)+vtppp(1)
!!$        velt(2,j) = vtest(2,j)-vtp(2)+vtpp(2)-vtppp(2)
!!$        velt(3,j) = vtest(3,j)+vtp(3)-vtpp(3)-vtppp(3)
!!$
!!$
!!$         vtottest(:)=vinf(:)+sigma(j)*velt(:,j)
!!$        cpt=1.-((sum(vtottest(:)**2))/(sum(vinf**2)))
!!$!        vtottest(:)=vinf(:)+sigma(j)*vtest(:,j)
!!$!        vtottest= vtottest(:)+vtest(:,j)
!!$!        cpt=1.-((sum(vtottest(:)**2))/(sum(vinf**2)))
!!$
!!$     end do   
!!$
!!$  write(6,*) 'v test'
!!$  write(6,'(3D16.8)') vtottest
!!$   write(6,*) 'cpt', cpt
   write(6,*) 'force', force

   do i=npanels,npanels+128
      row=(i-npanels)
      profilewave(:,row) = (/center(1,i)+deltax, zeta(i)/)
   end do


   !end do

!!Write some important data to the output file
!!----------------------------------------------------------------------
     CALL output( outputfile, iter, sourcesink, cw, profilewave )
!!----------------------------------------------------------------------

  !! Lets try to store the waves on each free surface panel...
  !!***************************************************************

    ! First we store data valid for a panel (cell in VTK jargon)

  ! Create a file name with the appropriate extension
   i = INDEX(outputfile, ".", BACK=.TRUE.)
   IF (i == 0) THEN
      WRITE(datafile1, '(A,A,A)')  TRIM(outputfile), 'wave', '.vtp'
   ELSE
      IF (i > 20) i = 20 ! To make sure the filename does not exceed 24 char
      WRITE(datafile1, '(A,A,A)')  outputfile(1:i-1), 'wave', '.vtp'
   ENDIF 
   PRINT*, datafile1


   iounit = 10
   CALL VtkXmlPolyDataCellScalar( iounit, datafile1, npoints, &
         0, 0, 0, npanels+nfspanels, &
         RESHAPE(points, (/3*npoints/)), &
         RESHAPE((panels-1), (/4*(npanels+nfspanels)/)), &
         scalar1=zeta, &
         scalar2=sigma, &
         namescalar1='zeta', &
         namescalar2='sigma')





  !! Lets try to store the waves on each just the y=0 free surface panels...
  !!***************************************************************

    ! First we store data valid for a panel (cell in VTK jargon)

  ! Create a file name with the appropriate extension
   i = INDEX(outputfile, ".", BACK=.TRUE.)
   IF (i == 0) THEN
      WRITE(datafile2, '(A,A,A)')  TRIM(outputfile), 'waveprofile', '.vtp'
   ELSE
      IF (i > 20) i = 20 ! To make sure the filename does not exceed 24 char
      WRITE(datafile2, '(A,A,A)')  outputfile(1:i-1), 'waveprofile', '.vtp'
   ENDIF 
   PRINT*, datafile2


   iounit = 10
   CALL VtkXmlPolyDataCellScalar( iounit, datafile2, npoints, &
         0, 0, 0, npanels+128, &
         RESHAPE(points, (/3*npoints/)), &
         RESHAPE((panels-1), (/4*(npanels+128)/)), &
         scalar1=zeta, &
         namescalar1='zeta')







  !! Lets try to store the cp on each hull panel...
  !!***************************************************************

    ! First we store data valid for a panel (cell in VTK jargon)

  ! Create a file name with the appropriate extension
   i = INDEX(outputfile, ".", BACK=.TRUE.)
   IF (i == 0) THEN
      WRITE(datafile3, '(A,A,A)')  TRIM(outputfile), 'scalar', '.vtp'
   ELSE
      IF (i > 20) i = 20 ! To make sure the filename does not exceed 24 char
      WRITE(datafile3, '(A,A,A)')  outputfile(1:i-1), 'scalar', '.vtp'
   ENDIF 
   PRINT*, datafile3


   iounit = 10
   CALL VtkXmlPolyDataCellScalar( iounit, datafile3, npoints, &
         0, 0, 0, npanels, &
         RESHAPE(points, (/3*npoints/)), &
         RESHAPE((panels-1), (/4*(npanels)/)), &
         scalar1=cp, &
         namescalar1='cp')       



  !!******************************************************************

  !! Now the total velocity vector on each panel...
  !!***************************************************************
   ! Create a file name with the appropriate extension
  i = INDEX(outputfile, ".", BACK=.TRUE.)
  IF (i == 0) THEN
     WRITE(datafile4, '(A,A,A)')  TRIM(outputfile), 'vtotal', '.vtp'
  ELSE
     IF (i > 20) i = 20 ! To make sure the filename does not exceed 24 char
     WRITE(datafile4, '(A,A,A)')  outputfile(1:i-1), 'vtotal', '.vtp'
  ENDIF 
  PRINT*, datafile4
  

  iounit = 10
  CALL VtkXmlSaveVectorFlowData( iounit, datafile4, npanels+nfspanels, center, &
        vtotal) 

! This method probably won't work for later versions of the program!

! vt1 and vt2 must be re-converted into vectors to display as vectors
! multiply coordsys()times the scalar vt1 or vt2 to scale the coordsys
! vector.  Then they should work.

  !!******************************************************************


  !! Now the vector data on each panel...
  !!***************************************************************
   ! Create a file name with the appropriate extension
  i = INDEX(outputfile, ".", BACK=.TRUE.)
  IF (i == 0) THEN
     WRITE(datafile5, '(A,A,A)')  TRIM(outputfile), 'components', '.vtp'
  ELSE
     IF (i > 20) i = 20 ! To make sure the filename does not exceed 24 char
     WRITE(datafile5, '(A,A,A)')  outputfile(1:i-1), 'components', '.vtp'
  ENDIF 
  PRINT*, datafile5
  

  iounit = 10
  CALL VtkXmlPolyDataCellScalar( iounit, datafile5, npoints, &
        0, 0, 0, npanels+nfspanels, &
        RESHAPE(points, (/3*npoints/)), &
        RESHAPE((panels-1), (/4*(npanels+nfspanels)/)), &
        scalar1=vn, &
        scalar2=vt1, &
        scalar3=vt2, &
        namescalar1='vn', &
        namescalar2='vt1', &
        namescalar3='vt2')
        


! Face specific
! Scalar components of the velocity vector
! correct format
  !!******************************************************************


  ! De-allocate all arrays
  IF (ALLOCATED(points)) DEALLOCATE(points)
  IF (ALLOCATED(panels)) DEALLOCATE(panels)


  IF (ALLOCATED(center)) DEALLOCATE(center)
  IF (ALLOCATED(coordsys)) DEALLOCATE(coordsys)
  IF (ALLOCATED(cornerslocal)) DEALLOCATE(cornerslocal)
  IF (ALLOCATED(area)) DEALLOCATE(area)

  IF (ALLOCATED(vel)) DEALLOCATE(vel)
!  IF (ALLOCATED(am)) DEALLOCATE(am)
!  IF (ALLOCATED(b)) DEALLOCATE(b)

  IF (ALLOCATED(sigma)) DEALLOCATE(sigma)

  IF (ALLOCATED(vtotal)) DEALLOCATE(vtotal)
  IF (ALLOCATED(vn)) DEALLOCATE(vn)
  IF (ALLOCATED(vt1)) DEALLOCATE(vt1)
  IF (ALLOCATED(vt2)) DEALLOCATE(vt2)
  IF (ALLOCATED(cp)) DEALLOCATE(cp)
  IF (ALLOCATED(vtest)) DEALLOCATE(vtest)
  IF (ALLOCATED(velt)) DEALLOCATE(velt)
  IF (ALLOCATED(zeta)) DEALLOCATE(zeta)

!  IF (ALLOCATED(celldata1)) DEALLOCATE(celldata1)
!  IF (ALLOCATED(celldata2)) DEALLOCATE(celldata2)

  ! Every program ends with end program
END PROGRAM p3
! end of file

