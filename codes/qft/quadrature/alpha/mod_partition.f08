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
                write ( stdout, '( "array is allocated" )' )
                status = -1
                return
            end if

            allocate ( me % mesh ( 0 : me % nCell ), stat = alloc_status, errmsg = alloc_msg )


    end function allocate_mesh_sub

end module mPartition
