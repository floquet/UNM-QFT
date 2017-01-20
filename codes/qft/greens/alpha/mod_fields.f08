module mFields

    use mConstants,                     only : one, zero, half, stdout, mille
    use mExtents,                       only : extents
    use mFileHandling,                  only : safeopen_readonly
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
        real ( rp ),                 dimension ( 0 : 3 )         :: G
        integer ( ip ), allocatable, dimension ( : )             :: ups, dns, upt, dnt
        ! rank 0
        integer ( ip ) :: naccept, nreject
        ! spatial, temporal extents
        type ( extents ) :: myExtents
        type ( masses )  :: myMasses
     contains
         private
         !procedure :: aA
         procedure, public :: housekeeping => housekeeping_sub
         procedure, public :: thermalize   => thermalize_sub
         !procedure, public :: update_f     => update_f_fcn
        !procedure, private, nopass :: allocate_rank_1_rp_sub
    end type fields

    ! local variables
    integer,                 private :: alloc_status  = -1
    character ( len = 512 ), private :: alloc_message = 'null'
    character ( len = * ),   private, parameter :: error_fatal = 'Program halting in module mFields due to fatal error.'

    private :: aA
    private :: allocate_rank_1_rp_sub, allocate_rank_3_rp_sub, allocate_rank_4_rp_sub
    private :: allocate_rank_1_ip_sub
    private :: housekeeping_sub, thermalize_sub

contains

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine thermalize_sub ( me, temp, farray )

        class ( fields ), target :: me

        character ( len = 64 ), intent ( in ) :: farray
        character ( len =  * ), intent ( in ) :: temp
        ! locals
        real ( rp ) :: random
        integer ( ip ) :: i, j, k, l
        integer :: io_in_farray

            if ( temp == 'hot' ) then
                write ( stdout, 100 ) 'The temp is ', temp
                io_in_farray = safeopen_readonly ( farray )
                read ( io_in_farray, * ) me % f
            else if ( temp == 'cold') then
                write ( stdout, 100 ) 'The temp is ', temp
                do i = 1, me % myExtents % Ns
                    do j = 1, me % myExtents % Ns
                        do k = 1, me % myExtents % Ns
                            do l = 1, me % myExtents % Nt
                                call random_number ( random )
                                me % f ( i, j, k, l ) = ( one - random ) * mille
                            end do ! l
                        end do ! k
                    end do ! j
                end do ! i
                else
                    write ( stdout, 100 ) 'Unrecognized temperature: should be "hot" or "cold". Input value was ', temp, '.'
                    stop 'Fatal error - I need to know the temperature.'
            end if
        100 format ( * ( g0 ) )

    end subroutine thermalize_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    function aA ( iu, ju, ku, lu, i, j, k, l ) result ( fcn_result )

        class ( fields ), target :: me

        integer ( ip ),  intent ( in ) :: iu, ju, ku, lu, i, j, k, l
        real ( rp ) :: fcn_result
        ! locals
        real ( rp ) :: phil, gphil, dphil, V

            phil = abs      ( me % f (  i,  j,  k,  l ) )
            gphil = sqrt( ( ( me % f ( iu,  j,  k,  l ) - me % f ( i, j, k, l ) )**2 &
                         +  ( me % f (  i, ju,  k,  l ) - me % f ( i, j, k, l ) )**2 &
                         +  ( me % f (  i,  j, ku,  l ) - me % f ( i, j, k, l ) )**2 ) )
            dphil = abs     ( me % f (  i,  j,  k, lu ) - me % f ( i, j, k, l ) )
            V = ( gphil / me % myExtents % as )**2 + ( me % myMasses % m * phil )**2
            fcn_result  = me % myExtents % at * sqrt ( V ) / &
                        ( me % myExtents % as**3 * ( me % myExtents % at**2 * V + dphil**2 ) )
     end function aA

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine update_f_fcn ( me )

        class ( fields ), target :: me

        type ( extents ), pointer :: ex

        real ( rp ) :: oldf, oldP, newP, rdn, df

        integer ( ip ) :: sweep
        integer ( ip ) :: i,  j,  k,  l,  &
                          iu, ju, ku, lu, &
                          id, jd, kd, ld

            ex => me % myExtents

            sweepdo: do sweep = 1, ex % Nsweeps
                ido: do i = 1, ex % Ns
                    iu = me % ups ( i )
                    id = me % dns ( i )
                    jdo: do j = 1, ex % Ns
                        ju = me % ups ( j )
                        jd = me % dns ( j )
                        kdo: do k = 1, ex % Ns
                            ku = me % ups ( k )
                            kd = me % dns ( k )
                            ldo: do l = 1, ex % Nt ! time
                                lu = me % upt ( l )
                                ld = me % dnt ( l )
                                oldf = me % f ( i, j, k, l )
                                oldP = aA ( iu, ju, ku, lu,  i,  j,  k,  l ) &
                                     * aA ( iu, ju, ku,  l,  i,  j,  k, ld ) &
                                     * aA (  i, ju, ku, lu, id,  j,  k,  l ) &
                                     * aA ( iu,  j, ku, lu,  i, jd,  k,  l ) &
                                     * aA ( iu, ju,  k, lu,  i,  j, kd,  l )
                                call random_number ( rdn )
                                me % f ( i, j, k, l ) = oldf + df * ( rdn - half )
                                newP = aA ( iu, ju, ku, lu,  i,  j,  k,  l ) &
                                     * aA ( iu, ju, ku,  l,  i,  j,  k, ld ) &
                                     * aA (  i, ju, ku, lu, id,  j,  k,  l ) &
                                     * aA ( iu,  j, ku, lu,  i, jd,  k,  l ) &
                                     * aA ( iu, ju,  k, lu,  i,  j, kd,  l )
                                if ( newP >= oldP ) then ! accept
                                    me % naccept = me % naccept + 1
                                else
                                    call random_number ( rdn )
                                    if ( newP/oldP >= rdn ) then ! accept
                                        me % naccept = me % naccept + 1
                                    else ! reject
                                        me % nreject = me % nreject + 1
                                        me % f ( i, j, k, l )  = oldf
                                    end if
                                end if
                            end do ldo
                        end do kdo
                    end do jdo
                end do ido
            end do sweepdo

            ex => null ( )

    end subroutine update_f_fcn

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
            me % naccept = 0
            me % nreject = 0

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
