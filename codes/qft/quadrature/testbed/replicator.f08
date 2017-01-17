program replicator

    use mReplicate, only : duplicate
    use mStride,    only : stride_up, stride_dn

    implicit none

    integer, parameter :: nDim = 3

    !rank 1
    integer, dimension ( 1 : 1000 ) :: myList = 1
    integer, dimension ( 1 : nDim ) :: shape_cells ( 1 : nDim ) = [ 2, 3, 4 ], shape_verts ( 1 : nDim ) = 0
    integer, dimension ( 1 : nDim ) :: stride_up_cells = 0, stride_dn_cells = 0, stride_up_verts = 0, stride_dn_verts = 0
    !rank 0
    integer :: nCells = 0, nVerts = 0, j = 0

        shape_verts = shape_cells ( : ) + 1

        call stride_up ( shape_cells, stride_up_cells, nCells )
        ! call stride_dn ( shape_cells, stride_dn_cells, nCells )

        call stride_up ( shape_verts, stride_up_verts, nVerts )
        ! call stride_dn ( shape_verts, stride_dn_verts, nVerts )

        myList ( 1 ) = 1 ! nucleate

        write ( *, 100 ) "shape_cells: ", shape_cells ( : )
        write ( *, 100 ) "shape_verts: ", shape_verts ( : )

        write ( *, 100 ) "number of cells    = ", nCells
        write ( *, 100 ) "number of vertices = ", nVerts

        write ( *, 100 ) "stride_up_cells = ", stride_up_cells
        write ( *, 100 ) "stride_up_verts = ", stride_up_verts

        call duplicate ( anchors = myList, from = 1, to = 2,  length = 1, shift = 1 )

        call duplicate ( anchors = myList, from = 1, to = 3,  length = 2, shift = 3 )
        call duplicate ( anchors = myList, from = 3, to = 5,  length = 2, shift = 3 )

        call duplicate ( anchors = myList, from = 1, to = 7,  length = 6, shift = 12 )
        call duplicate ( anchors = myList, from = 7, to = 13, length = 6, shift = 12 )
        !call duplicate ( anchors = myList, from = 7, to = 12, length = 6, shift = 12 )
        !call replicate ( anchors = myList, length = 2, position = p , copies = 2, shift = 3 )

        do j = 1, 25
            write ( *, 100 ) "cell ", j, " anchor vertex = ", myList ( j )
        end do

    stop 'successful completion for replicator . . .'

    100 format ( * ( g0, ' ' ) )

end program replicator
