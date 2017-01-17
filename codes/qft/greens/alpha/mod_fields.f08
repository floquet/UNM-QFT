module mFields

    use mConstants,                     only : one, stdout
    use mExtents,                       only : extents
    use mSetPrecision,                  only : ip, rp

    implicit none

    type :: fields
        ! rank 4
        real ( rp ), allocatable, dimension ( : , : , : , : ) :: f
        ! rank 3
        real ( rp ), allocatable, dimension ( : , : , : )     :: A, C
        ! rank 1
        real ( rp ), allocatable, dimension ( : )             :: phi, gphi, dphi
        ! spatial, temporal extents
        type ( extents ) :: myExtents
    contains
        private
        procedure, public :: allocator       => allocator_fcn
        procedure, public :: allocate_rank_1 => allocate_rank_1_fcn
    end type fields

    ! local variables
    integer , private                :: alloc_status  = -1
    character ( len = 512 ), private :: alloc_message = 'null'
    character ( len = * ),   private, parameter :: error_fatal   = 'Program halting in module mFields due to fatal error.'

    private :: allocator_fcn
    private :: allocate_rank_1_fcn

contains

    function allocator_fcn ( me ) result ( fcn_success )

        class ( fields ), target        :: me

        integer :: fcn_success

            fcn_success = 0
            ! me % allocate_rank_4_fcn ( me, me % f, extents % Ns, extents % Ns, extents % Ns, extents % Nt )
            !
            ! status_alloc = status_alloc + allocate_rank_3_fcn ( me, me % A, extents % Nphi + 1, extents % Nphi + 1, extents % Nphi + 1 )
            ! status_alloc = status_alloc + allocate_rank_3_fcn ( me, me % C, extents % Nphi + 1, extents % Nphi + 1, extents % Nphi + 1 )

            fcn_success = allocate_rank_1_fcn ( me, me % phi,  me % myExtents % Nphi + 1 )
            fcn_success = allocate_rank_1_fcn ( me, me % gphi, me % myExtents % Nphi + 1 )
            fcn_success = allocate_rank_1_fcn ( me, me % dphi, me % myExtents % Nphi + 1 )

    end function allocator_fcn

    function allocate_rank_1_fcn ( me, array, length ) result  ( fcn_success )

        class ( fields ), target                   :: me

        real ( rp ), allocatable, intent ( inout ) :: array ( : )
        integer ( ip ),           intent ( in )    :: length
        integer                                    :: fcn_success

            fcn_success = -1

            ! deallocate if needed
            if ( allocated ( array ) ) then
                write ( stdout, 100 ) 'Warning: deallocating rank 1 array...'
                fcn_success = 1
                deallocate ( array, stat = alloc_status, errmsg = alloc_message )
                if ( alloc_status /= 0 ) then
                    write ( stdout, 100 ) 'Error deallocating rank 1 array of ', length, ' elements, type real ( rp )'
                    write ( stdout, 100 ) 'error message: ', trim ( alloc_message ), '.'
                    write ( stdout, 100 ) 'error number:  ', alloc_status
                    stop error_fatal
                end if
            end if

            ! allocate array
            allocate ( array ( 1 : length ), stat = alloc_status, errmsg = alloc_message )
            if ( alloc_status /= 0 ) then
                write ( stdout, 100 ) 'Error allocating rank 1 array of ', length, ' elements, type real ( rp )'
                write ( stdout, 100 ) 'error message: ', trim ( alloc_message ), '.'
                write ( stdout, 100 ) 'error number:  ', alloc_status
                stop error_fatal
            end if

            fcn_success = 0 ! success

        return

    100 format ( * ( g0 ) )

    end function allocate_rank_1_fcn

end module mFields
