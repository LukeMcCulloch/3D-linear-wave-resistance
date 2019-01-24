!!
!! Simple program to read fifi and return the array to hw1
!!
!!
MODULE io    ! Every Fortran program starts with a program statement



  use precise, only : defaultp
  ! all variable declarations come here (before CONTAINS)
  ! Except subroutines
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp

CONTAINS

!! Subroutine to read the input file-----------------------------------
!! Subroutint to read the input file-----------------------------------
  SUBROUTINE input( inputfile, outputfile, title, &
                         npanels, nfspanels, npoints, &
                         panels, points,i,j,deltax )

!     INTEGER :: narg      ! Integer variables
     INTEGER :: i, j, n
     INTEGER :: npoints, npanels, nfspanels

     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: panels

     CHARACTER(len=*) :: inputfile   ! file name
     CHARACTER(len=*) :: outputfile  ! file name
     CHARACTER(len=*) :: title

     LOGICAL :: flexists  ! a logical variable. I .true. or .false.

     REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: points
     real(wp) :: deltax

     ! Variable declarations go before any executable statement
     
     OPEN(10,file=inputfile) ! open the file
     
     ! Read the data-------------------------------------------------
     read(10,'(AAAAA)') title
     read(10,*) i,j
     read(10,*) npanels, nfspanels, npoints
     read(10,*) deltax
     ALLOCATE( panels(4,npanels+nfspanels) )
     ALLOCATE( points(3,npoints) )
     read(10,*) panels
     read(10,*) points
     close(10)
     ! Finished Reading the data-------------------------------------
  
     ! Write the fifi data (just because I want to)------------------
     OPEN(20,file=outputfile)
     Write(20,*) ' Luke McCulloch, NAME 6160 Project, Outputfile '
     write(20,*)
     write(20,*) '(1.) Read the fifi data file and copy to "out.dat" '
     write(20,*) ' using the subroutine "input" via module "io"'
     write(20,*) '---------------------------------------------------------'
     write(20,*) title
     WRITE(20, FMT='(AI3A)') ' we have ', npanels+nfspanels, ' panels '
     WRITE(20, FMT='(AI3A)') ' we have ', npoints, ' points '
     write(20,'(4I4)') panels
     write(20,'(3D16.8)') points  ! Use Implicit do loop array formating
     write(20,*) '--------------------------------------------------------'

     close(20)
  
     write(*,*)  ! Blank space
     !  IF (ALLOCATED(points)) DEALLOCATE(points)
     !  IF (ALLOCATED(panels)) DEALLOCATE(panels)
     !  IF (ALLOCATED(panels)) DEALLOCATE(op)

  END SUBROUTINE input
!! End of Subroutine to read the data ----------------------------------
!! End of Subroutine to read the data ----------------------------------



!! Subroutine writes to file the data manipulated by the main program---
!! Subroutine writes to file the data manipulated by the main program---

  SUBROUTINE output( outputfile, iter, sourcesink, cw, profilewave )
  
     CHARACTER(len=*) :: outputfile
     !real(wp), allocatable, dimension(:,:) :: op
     real(wp) :: sourcesink, cw
     integer :: iter
     real(wp), dimension(2,128) :: profilewave

     Open(20,file= outputfile, Access= 'append', status='old')
     write(20,*)
     Write(20,*) ' Results From the Wigley Hull Linear Free Surface Panel Code '
     write(20,*) '---------------------------------------------------------------'
     write(20,*) ' (2.) Output data using the subroutine "output" in module "io" '
     write(20,*)
     write(20,*) 'source sink summation over all panels'
     write(20,'(3D16.8)') sourcesink
     write(20,*)
     write(20,*) 'Wave Resistance'
     write(20,'(3D16.8)') cw
     write(20,*)
     write(20,*) 'iter count from simquit'
     write(20,'(I4)') iter
     write(20,*)
     write(20,*)
     write(20,*) 'now output the wave profile on the centerline + ship hull'
     write(20,'(2D16.8)') profilewave
     write(20,*) '----------------------------------------------------------------'
     close(20)
     
  END SUBROUTINE output
     
!! End of Subroutine to write the data ---------------------------------
!! End of Subroutine to write the data ---------------------------------

END MODULE io


