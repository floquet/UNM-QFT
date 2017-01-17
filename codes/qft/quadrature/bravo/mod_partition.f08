module mPartition

    use mConstants,                     only : stdout
    use mParameters,                    only : nDim
    use mSetPrecision,                  only : ip, rp

    implicit none

    type :: partitions
        real    ( rp ), allocatable :: mesh ( : )
        integer ( ip )              :: nCell
    contains
        private
        procedure, public :: allocate_mesh => allocate_mesh_sub
    end type partitions

    type :: GrandPartition
        type ( partitions ) :: partition ( 1 : nDim )
        integer ( ip )      :: nTotalCells
    end type GrandPartition

contains

!   =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    function allocate_mesh_sub ( me ) result ( status )

        class ( partitions ), target :: me

        ! in and out
        integer :: status
        ! local variables
        integer                 :: alloc_status = 0
        character ( len = 256 ) :: alloc_msg    = ''

            status = 0

            if ( allocated ( me % mesh ) ) then
                status = -1
                write ( stdout, 100 )
                write ( stdout, 110 )
                write ( stdout, 200 ) status
                stop 'Fatal error'
            end if

            allocate ( me % mesh ( 0 : me % nCell ), stat = alloc_status, errmsg = alloc_msg )
            if ( alloc_status /= 0 ) then
                write ( stdout, 100 )
                write ( stdout, 120 ) trim ( alloc_msg ), alloc_status
                write ( stdout, 200 ) status
                status = -2
                stop 'Fatal error'
            end if

            status = 0  ! success

        return

    100 format ( /, 'Error!', /, 'module: mPartition', /, 'function:  allocate_mesh_sub', / )
    110 format ( 'partitions % mesh is already allocated' )
    120 format ( 'allocation errmsg = ', g0, /, 'allocation alloc_status = ', g0 )

    200 format ( 'status = ', g0, / )

    end function allocate_mesh_sub

end module mPartition
