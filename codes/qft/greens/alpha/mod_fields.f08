module mFields

    use mConstants,                     only : one, stdout
    use mExtents,                       only : extents
    use mSetPrecision,                  only : ip, rp

    implicit none

    type :: fields
        ! rank 4
        real ( rp ),    allocatable, dimension ( : , : , : , : ) :: f
        ! rank 3
        real ( rp ),    allocatable, dimension ( : , : , : )     :: A, C
        ! rank 1
        real ( rp ),    allocatable, dimension ( : )             :: phi, gphi, dphi
        integer ( ip ), allocatable, dimension ( : )             :: ups, dns, upt, dnt
        ! spatial, temporal extents
        type ( extents ) :: myExtents
    ! contains
    !     private
    !     procedure, public :: housekeeping => housekeeping_sub
    end type fields

    ! local variables
    integer , private                :: alloc_status  = -1
    character ( len = 512 ), private :: alloc_message = 'null'
    character ( len = * ),   private, parameter :: error_fatal = 'Program halting in module mFields due to fatal error.'

    private :: allocator_sub
    private :: allocate_rank_1_rp_sub
    private :: allocate_rank_1_ip_sub

contains

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine housekeeping_sub ( me )

        class ( fields ), target :: me

            call allocator_sub ( me ) ! allocate all arrays
            call neighbors_sub ( me ) ! populate pointers to neighbors

    end subroutine housekeeping_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine allocator_sub ( me )

        class ( fields ), target :: me

            ! me % allocate_rank_4_sub ( me, me % f, extents % Ns, extents % Ns, extents % Ns, extents % Nt )
            !
            ! status_alloc = status_alloc + allocate_rank_3_sub ( me, me % A, extents % Nphi + 1, extents % Nphi + 1, extents % Nphi + 1 )
            ! status_alloc = status_alloc + allocate_rank_3_sub ( me, me % C, extents % Nphi + 1, extents % Nphi + 1, extents % Nphi + 1 )

            ! fcn_success = allocate_rank_1_rp_sub ( me, me % phi,  me % myExtents % Nphi + 1 )
            ! fcn_success = allocate_rank_1_rp_sub ( me, me % gphi, me % myExtents % Nphi + 1 )
            ! fcn_success = allocate_rank_1_rp_sub ( me, me % dphi, me % myExtents % Nphi + 1 )
            call allocate_rank_1_rp_sub ( me, me % phi,  me % myExtents % Nphi + 1 )
            call allocate_rank_1_rp_sub ( me, me % gphi, me % myExtents % Nphi + 1 )
            call allocate_rank_1_rp_sub ( me, me % dphi, me % myExtents % Nphi + 1 )
            ! pbc
            call allocate_rank_1_ip_sub ( me, me % ups,  me % myExtents % Ns )
            call allocate_rank_1_ip_sub ( me, me % dns,  me % myExtents % Ns )
            call allocate_rank_1_ip_sub ( me, me % upt,  me % myExtents % Nt )
            call allocate_rank_1_ip_sub ( me, me % dnt,  me % myExtents % Nt )

    end subroutine allocator_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine allocate_rank_1_rp_sub ( me, array, length )
    !subroutine allocate_rank_1_rp_sub ( me, array, length ) result  ( fcn_success )

        class ( fields ), target :: me

        real ( rp ), allocatable, intent ( inout ) :: array ( : )
        integer ( ip ),           intent ( in )    :: length

            ! deallocate if needed
            if ( allocated ( array ) ) then
                write ( stdout, 100 ) 'Warning: deallocating rank 1 array...'
                deallocate ( array, stat = alloc_status, errmsg = alloc_message )
                if ( alloc_status /= 0 ) then
                    write ( stdout, 100 ) 'Error deallocating rank 1 array of ', length, ' elements, type integer ( ip )'
                    write ( stdout, 100 ) 'error message: ', trim ( alloc_message ), '.'
                    write ( stdout, 100 ) 'error number:  ', alloc_status
                    stop error_fatal
                end if
            end if

            ! allocate array
            allocate ( array ( 1 : length ), stat = alloc_status, errmsg = alloc_message )
            if ( alloc_status /= 0 ) then
                write ( stdout, 100 ) 'Error allocating rank 1 array of ', length, ' elements, type integer ( ip )'
                write ( stdout, 100 ) 'error message: ', trim ( alloc_message ), '.'
                write ( stdout, 100 ) 'error number:  ', alloc_status
                stop error_fatal
            end if

        return

    100 format ( * ( g0 ) )

    end subroutine allocate_rank_1_rp_sub
    !end subroutine allocate_rank_1_rp_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine allocate_rank_1_ip_sub ( me, array, length )

        class ( fields ), target :: me

        integer ( ip ), allocatable, intent ( inout ) :: array ( : )
        integer ( ip ),              intent ( in )    :: length

            ! deallocate if needed
            if ( allocated ( array ) ) then
                write ( stdout, 100 ) 'Warning: deallocating rank 1 array...'
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

        return

    100 format ( * ( g0 ) )

    end subroutine allocate_rank_1_ip_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine neighbors_sub ( me )

        class ( fields ), target :: me
        ! local
        integer ( ip ) :: i

            ! make tables for periodic boundary conditions
            ! time
            me % upt ( me % myExtents % Nt ) = 1
            me % dnt ( 1 )  = me % myExtents % Nt
            do i = 1, me % myExtents % Nt - 1
                me % upt ( i ) = i + 1
            end do
            do i = 2, me % myExtents % Nt
                me % dnt ( i ) = i - 1
            end do
            ! space
            me % ups ( me % myExtents % Ns ) = 1
            me % dns ( 1 )  = me % myExtents % Ns
            do i = 1, me % myExtents % Ns - 1
                me % ups ( i ) = i + 1
            end do
            do i = 2, me % myExtents % Ns
                me % dns ( i ) = i - 1
            end do

    end subroutine neighbors_sub

end module mFields
