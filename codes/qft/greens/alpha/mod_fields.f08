module mFields

    use mConstants,                     only : one, zero, stdout
    use mExtents,                       only : extents
    use mMasses,                        only : masses
    use mSetPrecision,                  only : ip, rp

    implicit none

    type :: fields
        ! rank 4
        real ( rp ),    allocatable, dimension ( : , : , : , : ) :: f
        ! rank 3
        real ( rp ),    allocatable, dimension ( : , : , : )     :: A, C
        ! rank 1
        real ( rp ),    allocatable, dimension ( : )             :: phi, gphi, dphi, E_0
        integer ( ip ), allocatable, dimension ( : )             :: ups, dns, upt, dnt
        real ( rp )                                              :: G ( 0 : 3 )
        ! rank 0
        ! spatial, temporal extents
        type ( extents ) :: myExtents
        type ( masses )  :: myMasses
     contains
    !     private
    !     procedure, public :: housekeeping => housekeeping_sub
        !procedure, private, nopass :: allocate_rank_1_rp_sub
    end type fields

    ! local variables
    integer,                 private :: alloc_status  = -1
    character ( len = 512 ), private :: alloc_message = 'null'
    character ( len = * ),   private, parameter :: error_fatal = 'Program halting in module mFields due to fatal error.'

    private :: allocator_sub
    private :: allocate_rank_1_rp_sub
    private :: allocate_rank_1_ip_sub

contains

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine greens_two_point ( me )

        class ( fields ), target :: me

        !real ( rp ), intent ( in ) :: as, at
        ! locals
        integer ( ip ) :: i, j, k, l, iu, ju, ku, lu, luu, luuu
        !real ( rp ) :: Esweep

            !Esweep = zero

            do i = 1, me % myExtents % Ns
                iu = me % ups ( i )
                do j = 1, me % myExtents % Ns
                    ju = me % ups ( j )
                    do k = 1, me % myExtents % Ns
                        ku = me % ups ( k )
                        do l = 1, me % myExtents % Nt ! time
                            lu = me % upt ( l ); luu = me % upt ( lu ); luuu = me % upt ( luu )
                            me % G ( 0 ) = me % G ( 0 ) + me % f ( i, j, k, l    ) * me % f ( i, j, k, l )
                            me % G ( 1 ) = me % G ( 1 ) + me % f ( i, j, k, lu   ) * me % f ( i, j, k, l )
                            me % G ( 2 ) = me % G ( 2 ) + me % f ( i, j, k, luu  ) * me % f ( i, j, k, l )
                            me % G ( 3 ) = me % G ( 3 ) + me % f ( i, j, k, luuu ) * me % f ( i, j, k, l )
                        end do
                    end do
                end do
            end do
            !  % E_0 ( sweep ) = me % myExtents % volume_ip / me % myExtents % avolume / me % volume_rp

    end subroutine greens_two_point

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine housekeeping_sub ( me )

        class ( fields ), target :: me

            call allocator_sub ( me ) ! allocate all arrays
            call neighbors_sub ( me ) ! populate pointers to neighbors
            me % G ( : ) = zero

    end subroutine housekeeping_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine allocator_sub ( me )

        class ( fields ), target :: me

            ! rank 4
            call allocate_rank_4_rp_sub ( me % f, me % myExtents % Ns, &
                                                  me % myExtents % Ns, &
                                                  me % myExtents % Ns, &
                                                  me % myExtents % Nt )
            ! rank 3
            call allocate_rank_3_rp_sub ( me % A, me % myExtents % Nphi  + 1, &
                                                  me % myExtents % Ngphi + 1, &
                                                  me % myExtents % Ndphi + 1 )
            call allocate_rank_3_rp_sub ( me % C, me % myExtents % Nphi  + 1, &
                                                  me % myExtents % Ngphi + 1, &
                                                  me % myExtents % Ndphi + 1 )
            ! rank 1
            call allocate_rank_1_rp_sub ( me % phi,  me % myExtents % Nphi + 1 )
            call allocate_rank_1_rp_sub ( me % gphi, me % myExtents % Nphi + 1 )
            call allocate_rank_1_rp_sub ( me % dphi, me % myExtents % Nphi + 1 )
            call allocate_rank_1_rp_sub ( me % E_0,  me % myExtents % Nsweeps )
            ! pbc
            call allocate_rank_1_ip_sub ( me % ups,  me % myExtents % Ns )
            call allocate_rank_1_ip_sub ( me % dns,  me % myExtents % Ns )
            call allocate_rank_1_ip_sub ( me % upt,  me % myExtents % Nt )
            call allocate_rank_1_ip_sub ( me % dnt,  me % myExtents % Nt )

    end subroutine allocator_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine allocate_rank_4_rp_sub ( array, length1, length2, length3, length4 )

        real ( rp ), allocatable, intent ( inout ) :: array ( : , : , : , : )
        integer ( ip ),           intent ( in )    :: length1, length2, length3, length4

            ! deallocate if needed
            if ( allocated ( array ) ) then
                write ( stdout, 100 ) 'Warning: deallocating rank 4 array...'
                deallocate ( array, stat = alloc_status, errmsg = alloc_message )
                if ( alloc_status /= 0 ) then
                    write ( stdout, 100 ) 'Error deallocating rank 4 array of ', length1, ' x ', &
                                                                                 length2, ' x ', &
                                                                                 length3, ' x ', &
                                                                                 length4, ' elements, type real ( rp )'
                    write ( stdout, 100 ) 'error message: ', trim ( alloc_message ), '.'
                    write ( stdout, 100 ) 'error number:  ', alloc_status
                    stop error_fatal
                end if
            end if

            ! allocate array
            allocate ( array ( 1 : length1, 1 : length2, 1 : length3, 1 : length4 ), stat = alloc_status, errmsg = alloc_message )
            if ( alloc_status /= 0 ) then
                write ( stdout, 100 ) 'Error allocating rank 4 array of ', length1, ' x ', &
                                                                           length2, ' x ', &
                                                                           length3, ' x ', &
                                                                           length4, ' elements, type real ( rp )'
                write ( stdout, 100 ) 'error message: ', trim ( alloc_message ), '.'
                write ( stdout, 100 ) 'error number:  ', alloc_status
                stop error_fatal
            end if

        return

    100 format ( * ( g0 ) )

    end subroutine allocate_rank_4_rp_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine allocate_rank_3_rp_sub ( array, length1, length2, length3 )

        real ( rp ), allocatable, intent ( inout ) :: array ( : , : , : )
        integer ( ip ),           intent ( in )    :: length1, length2, length3

            ! deallocate if needed
            if ( allocated ( array ) ) then
                write ( stdout, 100 ) 'Warning: deallocating rank 1 array...'
                deallocate ( array, stat = alloc_status, errmsg = alloc_message )
                if ( alloc_status /= 0 ) then
                write ( stdout, 100 ) 'Error deallocating rank 4 array of ', length1, ' x ', &
                                                                             length2, ' x ', &
                                                                             length3, ' elements, type real ( rp )'
                    write ( stdout, 100 ) 'error message: ', trim ( alloc_message ), '.'
                    write ( stdout, 100 ) 'error number:  ', alloc_status
                    stop error_fatal
                end if
            end if

            ! allocate array
            allocate ( array ( 1 : length1, 1 : length2, 1 : length3 ), stat = alloc_status, errmsg = alloc_message )
            if ( alloc_status /= 0 ) then
                write ( stdout, 100 ) 'Error allocating rank 4 array of ', length1, ' x ', &
                                                                           length2, ' x ', &
                                                                           length3, ' elements, type real ( rp )'
                write ( stdout, 100 ) 'error message: ', trim ( alloc_message ), '.'
                write ( stdout, 100 ) 'error number:  ', alloc_status
                stop error_fatal
            end if

        return

    100 format ( * ( g0 ) )

    end subroutine allocate_rank_3_rp_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine allocate_rank_1_rp_sub ( array, length )

        real ( rp ), allocatable, intent ( inout ) :: array ( : )
        integer ( ip ),           intent ( in )    :: length

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
                write ( stdout, 100 ) 'Error allocating rank 1 array of ', length, ' elements, type integer ( ip )'
                write ( stdout, 100 ) 'error message: ', trim ( alloc_message ), '.'
                write ( stdout, 100 ) 'error number:  ', alloc_status
                stop error_fatal
            end if

        return

    100 format ( * ( g0 ) )

    end subroutine allocate_rank_1_rp_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine allocate_rank_1_ip_sub ( array, length )

        integer ( ip ), allocatable, intent ( inout ) :: array ( : )
        integer ( ip ),              intent ( in )    :: length

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
