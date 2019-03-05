MODULE geometry

  USE precise, ONLY : defaultp



  IMPLICIT NONE    ! makes sure the compiler complains about 
                   ! undeclared variables
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp

CONTAINS

  FUNCTION cross_product( a, b ) RESULT (c)
    !
    ! Function to compute the cross product of two 3D vectors a and b.
    ! It will return a 3D vector.
    ! Usage: c = cross_product(a,b)
    !
    REAL(wp), DIMENSION(3) :: a, b, c

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  END FUNCTION cross_product


  SUBROUTINE panelgeometry(corners, coordsys, cornerslocal, area, center )
    !!
    !! procedure to compute the geometric properties of a
    !! quadrilateral panel. The actual quadrilateral is replaced by
    !! a flat panel.
    !! Input are the four corner points (3D). The points must form
    !! an anti-clockwise sequence if viewed from inside the hull. 
    !! The normal vector will then point inside the hull.
    !!
    !! Output:
    !!     coordsys: a 3x3 matrix containing the two tangent vectors
    !!               and the normal vector
    !!                  t1x   t2x   nx
    !!                  t1y   t2y   ny
    !!                  t1z   t2z   nz
    !!               i.e. coordsys(:,3) is the normal vector
    !!     cornerlocal: 2D - coordinates of the panel corners with
    !!                  respect to the panel center
    !!     area: the area of the panel
    !!     center: the 3D coordinates of the panel center
    !!
    INTEGER :: i

    INTEGER, DIMENSION(4), PARAMETER :: ip1 = (/ 2, 3, 4, 1 /)

    REAL(wp) :: area, dxi, deta
    REAL(wp) :: momentxi, momenteta, xic, etac

    REAL(wp), DIMENSION(3,4) :: corners      ! input, coordinates of 
                                             ! the four corner points
    REAL(wp), DIMENSION(3) :: m1, m2, m3, m4 ! midpoints
    REAL(wp), DIMENSION(3) :: iv, jv, jvbar, nv ! temporary vectors
    REAL(wp), DIMENSION(3) :: center            ! center of panel
    REAL(wp), DIMENSION(3,3) :: coordsys        ! local coordinate system
    REAL(wp), DIMENSION(2,4) :: cornerslocal    ! local coordinate system

    ! first compute the midpoints of all four panel edges. The points
    ! m1 through m4 are situated on a plane even if the corners are not.
    m1(:) = 0.5*(corners(:,1)+corners(:,2))
    m2(:) = 0.5*(corners(:,2)+corners(:,3))
    m3(:) = 0.5*(corners(:,3)+corners(:,4))
    m4(:) = 0.5*(corners(:,4)+corners(:,1))

    ! second compute a coordinate system local to the panel
    iv = m2-m4                  ! first local coordinate vector
    iv = iv / SQRT(SUM(iv**2))  ! normalize to unit length

    jvbar = m3-m1 ! second temporary local coordinate vector

    ! Normal vector
    nv = cross_product( iv, jvbar )
    nv = nv / SQRT(SUM(nv**2))  ! normalize to unit length
    
    ! find second tangent vector which is perpendicular to nv and iv
    jv = cross_product( nv, iv )

    ! store the local coordinate system vectors in the matrix coordsys
    coordsys(:,1) = iv
    coordsys(:,2) = jv
    coordsys(:,3) = nv
 
    ! temporary origin
    center = 0.25*(m1+m2+m3+m4)

    ! corners in local (plane) coordinate system
    DO i = 1, 4
       ! xi component
       cornerslocal(1,i) = DOT_PRODUCT( (corners(:,i)-center), iv) 
       ! eta component
       cornerslocal(2,i) = DOT_PRODUCT( (corners(:,i)-center), jv) 
    END DO 
    
    ! area and moment of panel. The equations are derived by means
    ! of Stokes theorem.
    area = 0.
    momentxi = 0.  ! the moments are with respect to the panel tangent 
                   ! vectors
    momenteta = 0.
    DO i = 1, 4
       dxi = cornerslocal(1,ip1(i))-cornerslocal(1,i) ! difference in xi
       deta = cornerslocal(2,ip1(i))-cornerslocal(2,i)! difference in eta
     
       area = area + deta*(cornerslocal(1,ip1(i))+cornerslocal(1,i))

       momentxi = momentxi + dxi * ( &
            cornerslocal(2,ip1(i))**2 &
               + cornerslocal(2,i)*cornerslocal(2,ip1(i)) &
                   + cornerslocal(2,i)**2 &
                       )
       momenteta = momenteta + deta * ( &
            cornerslocal(1,ip1(i))**2 &
               + cornerslocal(1,i)*cornerslocal(1,ip1(i)) &
                   + cornerslocal(1,i)**2 &
                       )
    END DO 

    area = 0.5*area
    momentxi = -1./6.*momentxi
    momenteta = 1./6.*momenteta

    ! centroid
    IF (area > 0.0000000001) THEN
       ! We make sure that the panel is not degenerated with zero area
       xic = momenteta / area
       etac = momentxi / area
    ELSE
       xic = 0.
       etac = 0.
    ENDIF

    ! Correct the local coordinates of the corners to be measured from
    ! the true panel centroid.
    DO i = 1, 4
       cornerslocal(1,i) = cornerslocal(1,i) - xic
       cornerslocal(2,i) = cornerslocal(2,i) - etac
    END DO 

    ! Finally move the centroid to its exact position
    center = center + xic*iv + etac*jv

  END SUBROUTINE panelgeometry

END MODULE geometry
