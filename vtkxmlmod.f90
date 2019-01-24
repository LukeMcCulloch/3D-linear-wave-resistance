MODULE vtkxmlmod

   !! Lothar Birk, 100211
  
   !! A module with procedures to write VTK XML datafiles
   !! This is not a general module. The subroutines are narrowly defined
   !! for the context of a panel method (BEM potential flow)

   !! Needs much more documentation

   USE precise, ONLY: DEFAULTP
   IMPLICIT NONE
   INTEGER, PARAMETER, PRIVATE :: WP=DEFAULTP

   INTERFACE vtkxmlDataArray
      MODULE PROCEDURE vtkxmlDataArrayI,  vtkxmlDataArrayR 
   END INTERFACE 


CONTAINS

   SUBROUTINE VtkXmlPolyDataCellScalar(io, file, &
           & npoints, nverts, nlines, nstrips, npolys, &
           & points, polys, &
           & scalar1, scalar2, scalar3, &
           & namescalar1, namescalar2, namescalar3 )
      !!
      !! Subroutine stores up to three scalar data values per cell.
      !! Cells are all quadrilaterals.
      !!
      !! Length of scalar? is equal to number of polygons (cells)! 
      !!
      CHARACTER(len=*) :: file
      INTEGER :: io, npoints, nverts, nlines, nstrips, npolys
      REAL(wp), DIMENSION(:) :: points
      REAL(wp), DIMENSION(:), OPTIONAL :: scalar1, scalar2, scalar3
      CHARACTER(len=*), OPTIONAL :: namescalar1, namescalar2, namescalar3

      INTEGER, DIMENSION(:) :: polys
      INTEGER, DIMENSION(:), ALLOCATABLE :: iarray

      INTEGER :: i
      CHARACTER(Len=255) :: tag

      ! open VTK XML file 
      CALL vtkxmlopen(io, file, 'PolyData')

      ! Store panels as polygons (all quadrilaterals)
      CALL vtkxmlwritetag(io,'<PolyData>')
      WRITE(tag,'(''<Piece NumberOfPoints="'',I5, &
                & ''" NumberOfVerts="'',I5, &
                & ''" NumberOfLines="'',I5, &
                & ''" NumberOfStrips="'',I5, &
                & ''" NumberOfPolys="'',I5,''">'')') &
                & npoints, nverts, nlines, nstrips, npolys
      CALL vtkxmlwritetag(io, TRIM(tag))
      CALL vtkxmlwritetag(io,'<Points>')
      CALL vtkxmldataarray(io, points, &
                               'type="Float32"', &
                               'NumberOfComponents="3"'&
                                )
      CALL vtkxmlwritetag(io, '</Points>')
      CALL vtkxmlwritetag(io, '<Polys>')

      CALL vtkxmldataarray(io, polys, &
                               'type="Int32"', &
                               'Name="connectivity"'&
                                )
      ! This assumes all polygons are quadrilaterals!
      ALLOCATE(iarray(npolys))
      DO i = 1, npolys
         iarray(i) = 4*i
      END DO 
      CALL vtkxmldataarray(io, iarray, &
                               'type="Int32"', &
                               'Name="offsets"'&
                                )

      CALL vtkxmlwritetag(io, '</Polys>')

      CALL vtkxmlwritetag(io, '<CellData>')
         IF (PRESENT(scalar1)) THEN 
            IF (PRESENT(namescalar1)) THEN
               tag = ''
               WRITE(tag, '(''Name="'',A,''"'')') TRIM(namescalar1)
            ELSE
               tag = 'Name="scalar1"'
            ENDIF  
            CALL vtkxmldataarray(io, scalar1, &
                               'type="Float32"', &
                               tag)
         END IF 
         IF (PRESENT(scalar2)) THEN 
            IF (PRESENT(namescalar2)) THEN
               tag = ''
               WRITE(tag, '(''Name="'',A,''"'')') TRIM(namescalar2)
            ELSE
               tag = 'Name="scalar2"'
            ENDIF  
            CALL vtkxmldataarray(io, scalar2, &
                               'type="Float32"', &
                               tag)
         END IF 
         IF (PRESENT(scalar3)) THEN 
            IF (PRESENT(namescalar3)) THEN
               tag = ''
               WRITE(tag, '(''Name="'',A,''"'')') TRIM(namescalar3)
            ELSE
               tag = 'Name="scalar3"'
            ENDIF  
            CALL vtkxmldataarray(io, scalar3, &
                               'type="Float32"', &
                               tag)
         END IF 
      CALL vtkxmlwritetag(io, '</CellData>')
      
      
      CALL vtkxmlwritetag(io, '</Piece>')
      CALL vtkxmlwritetag(io, '</PolyData>')
      CALL vtkxmlclose(io)

      DEALLOCATE(iarray)
   END SUBROUTINE VtkXmlPolyDataCellScalar

   SUBROUTINE VtkXmlPolyDataCellVector(io, file, &
           & npoints, nverts, nlines, nstrips, npolys, &
           & points, polys, &
           & vector1, vector2, vector3, &
           & namevector1, namevector2, namevector3 )
      !!
      !! Subroutine stores up to three vectors per cell.
      !! Cells are all quadrilaterals.
      !!
      CHARACTER(len=*) :: file
      INTEGER :: io, npoints, nverts, nlines, nstrips, npolys
      REAL(wp), DIMENSION(:) :: points
      REAL(wp), DIMENSION(:), OPTIONAL :: vector1, vector2, vector3
      CHARACTER(len=*), OPTIONAL :: namevector1, namevector2, namevector3

      INTEGER, DIMENSION(:) :: polys
      INTEGER, DIMENSION(:), ALLOCATABLE :: iarray

      INTEGER :: i
      CHARACTER(Len=255) :: tag

      ! open VTK XML file 
      CALL vtkxmlopen(io, file, 'PolyData')

      ! Store panels as polygons (quadrilaterals)
      CALL vtkxmlwritetag(io,'<PolyData>')
      WRITE(tag,'(''<Piece NumberOfPoints="'',I5, &
                & ''" NumberOfVerts="'',I5, &
                & ''" NumberOfLines="'',I5, &
                & ''" NumberOfStrips="'',I5, &
                & ''" NumberOfPolys="'',I5,''">'')') &
                & npoints, nverts, nlines, nstrips, npolys
      CALL vtkxmlwritetag(io, TRIM(tag))
      CALL vtkxmlwritetag(io,'<Points>')
      CALL vtkxmldataarray(io, points, &
                               'type="Float32"', &
                               'NumberOfComponents="3"'&
                                )
      CALL vtkxmlwritetag(io, '</Points>')
      CALL vtkxmlwritetag(io, '<Polys>')

      CALL vtkxmldataarray(io, polys, &
                               'type="Int32"', &
                               'Name="connectivity"'&
                                )
      ALLOCATE(iarray(npolys))
      DO i = 1, npolys
         iarray(i) = 4*i
      END DO 
      CALL vtkxmldataarray(io, iarray, &
                               'type="Int32"', &
                               'Name="offsets"'&
                                )

      CALL vtkxmlwritetag(io, '</Polys>')

      CALL vtkxmlwritetag(io, '<CellData>')
         IF (PRESENT(vector1)) THEN 
            IF (PRESENT(namevector1)) THEN
               tag = ''
               WRITE(tag, '(''Name="'',A,''" NumberOfComponents=3'')') &
                  & TRIM(namevector1)
            ELSE
               tag = 'Name="vector1"'
            ENDIF  
            CALL vtkxmldataarray(io, vector1, &
                               'type="Float32"', &
                               tag)
         END IF 
         IF (PRESENT(vector2)) THEN 
            IF (PRESENT(namevector2)) THEN
               tag = ''
               WRITE(tag, '(''Name="'',A,''" NumberOfComponents=3'')') &
                  &  TRIM(namevector2)
            ELSE
               tag = 'Name="vector2"'
            ENDIF  
            CALL vtkxmldataarray(io, vector2, &
                               'type="Float32"', &
                               tag)
         END IF 
         IF (PRESENT(vector3)) THEN 
            IF (PRESENT(namevector3)) THEN
               tag = ''
               WRITE(tag, '(''Name="'',A,''" NumberOfComponents=3'')') & 
                  & TRIM(namevector3)
            ELSE
               tag = 'Name="vector3"'
            ENDIF  
            CALL vtkxmldataarray(io, vector3, &
                               'type="Float32"', &
                               tag)
         END IF 
     CALL vtkxmlwritetag(io, '</CellData>')
      
      
      CALL vtkxmlwritetag(io, '</Piece>')
      CALL vtkxmlwritetag(io, '</PolyData>')
      CALL vtkxmlclose(io)

      DEALLOCATE(iarray)
   END SUBROUTINE VtkXmlPolyDataCellVector


   SUBROUTINE VtkXmlPolyDataPtScalar(io, file, &
           & npoints, nverts, nlines, nstrips, npolys, &
           & points, polys, &
           & scalar1, scalar2, scalar3, &
           & namescalar1, namescalar2, namescalar3 )
      !!
      !! Subroutine stores up to three scalar data values per vertex
      !! of the cell.
      !! Cells are all quadrilaterals.
      !!
      !! Length of scalar? is equal to number of points! 
      !!
      CHARACTER(len=*) :: file
      INTEGER :: io, npoints, nverts, nlines, nstrips, npolys
      REAL(wp), DIMENSION(:) :: points
      REAL(wp), DIMENSION(:), OPTIONAL :: scalar1, scalar2, scalar3
      CHARACTER(len=*), OPTIONAL :: namescalar1, namescalar2, namescalar3

      INTEGER, DIMENSION(:) :: polys
      INTEGER, DIMENSION(:), ALLOCATABLE :: iarray

      INTEGER :: i
      CHARACTER(Len=255) :: tag

      ! open VTK XML file 
      CALL vtkxmlopen(io, file, 'PolyData')

      ! Store panels as polygons (all quadrilaterals)

      CALL vtkxmlwritetag(io,'<PolyData>')
      WRITE(tag,'(''<Piece NumberOfPoints="'',I5, &
                & ''" NumberOfVerts="'',I5, &
                & ''" NumberOfLines="'',I5, &
                & ''" NumberOfStrips="'',I5, &
                & ''" NumberOfPolys="'',I5,''">'')') &
                & npoints, nverts, nlines, nstrips, npolys
      CALL vtkxmlwritetag(io, TRIM(tag))
      CALL vtkxmlwritetag(io,'<Points>')
      CALL vtkxmldataarray(io, points, &
                               'type="Float32"', &
                               'NumberOfComponents="3"'&
                                )
      CALL vtkxmlwritetag(io, '</Points>')
      CALL vtkxmlwritetag(io, '<Polys>')

      CALL vtkxmldataarray(io, polys, &
                               'type="Int32"', &
                               'Name="connectivity"'&
                                )
      ! Assumes all cells are quadrilaterals)
      ALLOCATE(iarray(npolys))
      DO i = 1, npolys
         iarray(i) = 4*i
      END DO 
      CALL vtkxmldataarray(io, iarray, &
                               'type="Int32"', &
                               'Name="offsets"'&
                                )

      CALL vtkxmlwritetag(io, '</Polys>')

      CALL vtkxmlwritetag(io, '<PointData>')
         IF (PRESENT(scalar1)) THEN 
            IF (PRESENT(namescalar1)) THEN
               tag = ''
               WRITE(tag, '(''Name="'',A,''"'')') TRIM(namescalar1)
            ELSE
               tag = 'Name="scalar1"'
            ENDIF  
            CALL vtkxmldataarray(io, scalar1, &
                               'type="Float32"', &
                               tag)
         END IF 
         IF (PRESENT(scalar2)) THEN 
            IF (PRESENT(namescalar2)) THEN
               tag = ''
               WRITE(tag, '(''Name="'',A,''"'')') TRIM(namescalar2)
            ELSE
               tag = 'Name="scalar2"'
            ENDIF  
            CALL vtkxmldataarray(io, scalar2, &
                               'type="Float32"', &
                               tag)
         END IF 
         IF (PRESENT(scalar3)) THEN 
            IF (PRESENT(namescalar3)) THEN
               tag = ''
               WRITE(tag, '(''Name="'',A,''"'')') TRIM(namescalar3)
            ELSE
               tag = 'Name="scalar3"'
            ENDIF  
            CALL vtkxmldataarray(io, scalar3, &
                               'type="Float32"', &
                               tag)
         END IF 
      CALL vtkxmlwritetag(io, '</PointData>')
      
      
      CALL vtkxmlwritetag(io, '</Piece>')
      CALL vtkxmlwritetag(io, '</PolyData>')
      CALL vtkxmlclose(io)

      DEALLOCATE(iarray)
   END SUBROUTINE VtkXmlPolyDataPtScalar

!---------------------------------------------------------------------
!---------------------------------------------------------------------

   SUBROUTINE VtkXmlSaveScalarFlowData(io, file, npanel, paneloffsets, &
             sigma, cp, vnv, zeta)
      CHARACTER(len=*) :: file
      INTEGER :: io, npanel
      REAL(wp), DIMENSION(:,:,:) :: paneloffsets
      REAL(wp), DIMENSION(:) :: sigma, cp
      REAL(wp), DIMENSION(:), OPTIONAL :: vnv, zeta

      INTEGER, DIMENSION(:,:), ALLOCATABLE :: iele
      INTEGER, DIMENSION(:), ALLOCATABLE :: iarray

      INTEGER :: i
      CHARACTER(Len=79) :: tag

      ! open VTK XML file 
      CALL vtkxmlopen(io, file, 'UnstructuredGrid')

      ! Store panels as Unstructured grid of quadrilaterals
      CALL vtkxmlwritetag(io,'<UnstructuredGrid>')
      WRITE(tag,'(''<Piece NumberOfPoints="'',I5,''" NumberOfCells="''&
           &,I5,''">'')') 4*npanel, npanel
      CALL vtkxmlwritetag(io, tag)
      CALL vtkxmlwritetag(io,'<Points>')
      CALL vtkxmldataarray(io, RESHAPE(paneloffsets, (/3*4*npanel/)), &
                               'type="Float32"', &
                               'NumberOfComponents="3"'&
                                )
      CALL vtkxmlwritetag(io, '</Points>')
      CALL vtkxmlwritetag(io, '<Cells>')
      ! create iele
      ALLOCATE(iele(4,npanel))
      DO i = 1, npanel
         iele(1,i) = (i-1)*4
         iele(2,i) = (i-1)*4 + 1
         iele(3,i) = (i-1)*4 + 2
         iele(4,i) = (i-1)*4 + 3
      END DO 
      CALL vtkxmldataarray(io, RESHAPE(iele, (/4*npanel/)), &
                               'type="Int32"', &
                               'Name="connectivity"',&
                               'NumberOfComponents="1"'&
                                )
      ALLOCATE(iarray(npanel))
      DO i = 1, npanel
         iarray(i) = 4*i
      END DO 
      CALL vtkxmldataarray(io, iarray, &
                               'type="Int32"', &
                               'Name="offsets"',&
                               'NumberOfComponents="1"'&
                                )
      iarray = 9  ! all quadrilateral panels
      CALL vtkxmldataarray(io, iarray, &
                               'type="Int32"', &
                               'Name="types"',&
                               'NumberOfComponents="1"'&
                                )
      CALL vtkxmlwritetag(io, '</Cells>')

      WRITE(tag,'(''<CellData Scalars="Cp">'')')
      CALL vtkxmlwritetag(io, tag)
      CALL vtkxmldataarray(io, sigma, &
                               'type="Float32"', &
                               'Name="Sigma"')
      CALL vtkxmldataarray(io, Cp, &
                               'type="Float32"', &
                               'Name="Cp"')
      IF (PRESENT(vnv)) THEN 
         CALL vtkxmldataarray(io, vnv, &
                               'type="Float32"', &
                               'Name="Vnormal"')
      END IF 
      IF (PRESENT(zeta)) THEN 
         CALL vtkxmldataarray(io, zeta, &
                               'type="Float32"', &
                               'Name="Zeta"')
      END IF 
      
      CALL vtkxmlwritetag(io, '</CellData>')
      
      CALL vtkxmlwritetag(io, '</Piece>')
      CALL vtkxmlwritetag(io, '</UnstructuredGrid>')
      CALL vtkxmlclose(io)

      DEALLOCATE(iarray, iele)
   END SUBROUTINE VtkXmlSaveScalarFlowData

   SUBROUTINE VtkXmlSaveVectorFlowData(io, file, npanel, panelcenters, &
             velocity, vt1, vt2)
      CHARACTER(len=*) :: file
      INTEGER :: io, npanel
      REAL(wp), DIMENSION(:,:) :: panelcenters
      REAL(wp), DIMENSION(:,:) :: velocity
      REAL(wp), DIMENSION(:,:), OPTIONAL :: vt1, vt2

      INTEGER, DIMENSION(:), ALLOCATABLE :: iele
      INTEGER, DIMENSION(:), ALLOCATABLE :: iarray

      INTEGER :: i, dim
      CHARACTER(Len=79) :: tag

      ! open VTK XML file 
      dim = SIZE(panelcenters,1)

      CALL vtkxmlopen(io, file, 'UnstructuredGrid')

      ! Store panels as Unstructured grid of quadrilaterals
      CALL vtkxmlwritetag(io,'<UnstructuredGrid>')
      WRITE(tag,'(''<Piece NumberOfPoints="'',I5,''" NumberOfCells="''&
           &,I5,''">'')') npanel, npanel
      CALL vtkxmlwritetag(io, tag)
      CALL vtkxmlwritetag(io,'<Points>')
      WRITE(tag,'(''NumberOfComponents="'',I1,''"'')') dim
      CALL vtkxmldataarray(io, RESHAPE(panelcenters, (/dim*npanel/)), &
                               'type="Float32"', &
                               tag &
                                )
      CALL vtkxmlwritetag(io, '</Points>')
      CALL vtkxmlwritetag(io, '<Cells>')
      ! create iele
      ALLOCATE(iele(npanel))
      DO i = 1, npanel
         iele(i) = (i-1)
      END DO 
      CALL vtkxmldataarray(io, iele, &
                               'type="Int32"', &
                               'Name="connectivity"',&
                               'NumberOfComponents="1"'&
                                )
      ALLOCATE(iarray(npanel))
      DO i = 1, npanel
         iarray(i) = i
      END DO 
      CALL vtkxmldataarray(io, iarray, &
                               'type="Int32"', &
                               'Name="offsets"',&
                               'NumberOfComponents="1"'&
                                )
      iarray = 1  ! all panels are points
      CALL vtkxmldataarray(io, iarray, &
                               'type="Int32"', &
                               'Name="types"',&
                               'NumberOfComponents="1"'&
                                )
      CALL vtkxmlwritetag(io, '</Cells>')

      WRITE(tag,'(''<PointData Vectors="Velocity">'')')
      !!WRITE(tag,'(''<CellData Vectors="Velocity">'')') !does not plot correctly
      CALL vtkxmlwritetag(io, tag)
      
      dim = SIZE(velocity,1)
      WRITE(tag,'(''NumberOfComponents="'',I1,''"'')') dim
      CALL vtkxmldataarray(io,  RESHAPE(velocity, (/dim*npanel/)), &
                               'type="Float32"', &
                               'Name="velocity"', &
                               tag)
!       IF (PRESENT(vn)) THEN
!          dim = SIZE(vn,1)
!          WRITE(tag,'(''NumberOfComponents="'',I1,''"'')') dim
!          CALL vtkxmldataarray(io,  RESHAPE(vn, (/dim*npanel/)), &
!                                'type="Float32"', &
!                                'Name="normal_velocity"', &
!                                tag)
!       ENDIF 
      IF (PRESENT(vt1)) THEN
         dim = SIZE(vt1,1)
         WRITE(tag,'(''NumberOfComponents="'',I1,''"'')') dim
         CALL vtkxmldataarray(io,  RESHAPE(vt1, (/dim*npanel/)), &
                               'type="Float32"', &
                               'Name="tangential_velocity_1"', &
                               tag)
      ENDIF 
      IF (PRESENT(vt2)) THEN
         dim = SIZE(vt2,1)
         WRITE(tag,'(''NumberOfComponents="'',I1,''"'')') dim
         CALL vtkxmldataarray(io,  RESHAPE(vt2, (/dim*npanel/)), &
                               'type="Float32"', &
                               'Name="tangential_velocity_2"', &
                               tag)
      ENDIF 
      CALL vtkxmlwritetag(io, '</PointData>')
      !!!!CALL vtkxmlwritetag(io, '</CellData>')! Does not plot properly
      
      CALL vtkxmlwritetag(io, '</Piece>')
      CALL vtkxmlwritetag(io, '</UnstructuredGrid>')
      CALL vtkxmlclose(io)

      DEALLOCATE(iarray, iele)
   END SUBROUTINE VtkXmlSaveVectorFlowData


   SUBROUTINE vtkxmlwritetag(io, tag)
     INTEGER :: io
     CHARACTER(len=*) :: tag

     WRITE(io,'(A)') tag

   END SUBROUTINE vtkxmlwritetag

   SUBROUTINE vtkxmlDataArrayI(io, dataarray, typestr, namestr, tag)
     ! only ASCII format for now
     INTEGER :: io
     CHARACTER(len=*) :: typestr, namestr
     CHARACTER(len=*), OPTIONAL :: tag
     INTEGER, DIMENSION(:) :: dataarray
     INTEGER :: i,j,n
     n = SIZE(dataarray)
     IF (PRESENT(tag))THEN
        WRITE(io,'(8A)') '<DataArray ',TRIM(typestr), ' ',TRIM(namestr),&
                         ' ', TRIM(tag), ' ', &
                         ' format="ascii">'
     ELSE
        WRITE(io,'(5A)') '<DataArray ',TRIM(typestr), ' ',TRIM(namestr), &
                         ' format="ascii">'
     ENDIF
     DO i = 1,n
        WRITE(io,'(I8,2X)', Advance='No') dataarray(i)
     END DO 
     WRITE(io,*)
     WRITE(io,'(A)') '</DataArray>'

   END SUBROUTINE vtkxmlDataArrayI

   SUBROUTINE vtkxmlDataArrayR(io, dataarray, typestr, namestr, tag)
     ! only ASCII format for now
     INTEGER :: io
     CHARACTER(len=*) :: typestr, namestr
     CHARACTER(len=*), OPTIONAL :: tag
     REAL(wp), DIMENSION(:) :: dataarray
     INTEGER :: i,j,n

     n = SIZE(dataarray)
     IF (PRESENT(tag))THEN
        WRITE(io,'(8A)') '<DataArray ',TRIM(typestr), ' ',TRIM(namestr),&
                         ' ', TRIM(tag), ' ', &
                         ' format="ascii">'
     ELSE
        WRITE(io,'(5A)') '<DataArray ',TRIM(typestr), ' ',TRIM(namestr), &
                         ' format="ascii">'
     ENDIF
     DO i = 1,n
        WRITE(io,'(G16.6,1X)', Advance='No') dataarray(i)
     END DO 
     WRITE(io,*)
     WRITE(io,'(A)') '</DataArray>'

   END SUBROUTINE vtkxmlDataArrayR

   SUBROUTINE vtkxmlopen(io, filename, TYPE)

     INTEGER :: io
     CHARACTER(len=*) :: filename
     CHARACTER(len=*), OPTIONAL :: TYPE
     CHARACTER(len=79) :: locdatatype
     IF (.NOT.(PRESENT(TYPE))) THEN
        PRINT*, "fifi"
        locdatatype = "UnstructuredGrid"
     ELSE 
        locdatatype = TRIM(TYPE)
     ENDIF
     !PRINT*, locdatatype

     OPEN(io, file=filename)
     WRITE(io,'(A,A,A)') '<VTKFile type="',TRIM(locdatatype),'">'
 
     RETURN 
   END SUBROUTINE vtkxmlopen

   SUBROUTINE vtkxmlclose(io)

     INTEGER io
     WRITE(io,'(A)') '</VTKFile>'
     CLOSE(io)

     RETURN
   END SUBROUTINE vtkxmlclose



END MODULE vtkxmlmod
