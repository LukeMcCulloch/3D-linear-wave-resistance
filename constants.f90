!!
!! Simple Fortran 90 module with program constants
!! We could have delcared pi directly in the program units that
!! need pi. However the purpose is to show how to use modules 
!! to make pieces of the code reusable and easier to maintain.
!! 
MODULE constants

use precise, only : defaultp
IMPLICIT NONE

INTEGER, PARAMETER, PRIVATE :: WP=defaultp

! This defines the number pi. 
! Declared as a parameter to makes sure
! the number pi is not overridden during program execution.
REAL(WP), PARAMETER :: one=1.0  
REAL(WP), PARAMETER :: four=4.0  
REAL(WP), PARAMETER :: pi=Four*ATAN(one)     ! using a function to define the
                                             ! parameter value is a Fortran 2003
                                             ! extension.
!real(wp), parameter :: pi=3.141592653589793 ! alternative form
                                             ! for older compilers 
REAL(WP), PARAMETER :: gravity=9.80665       ! metric
REAL(WP), PARAMETER :: length=1.0            ! normalized by length

END MODULE constants
