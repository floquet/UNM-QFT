! 3456789 123456789 223456789 323456789 423456789 523456789 623456789 723456789 823456789 923456789 023456789 123456789 223456789 32
module mFields

    use mConstants,                     only : one, half, stdout, mille, fmt_generic
    use mExtents,                       only : extents
    use mFileHandling,                  only : safeopen_readonly, safeopen_writenew
    use mInputs,                        only : inputs
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
        real ( rp ) :: meanE, sigma
        real ( rp ) :: A_max, A_min, C_max, C_min
        real ( rp ) :: ratio_accept, ratio_reject
        real ( rp ) :: highphi, highgphi, highdphi

        integer ( ip ) :: naccept, nreject

        ! derived types
        type ( extents ) :: myExtents ! spatial, temporal extents
        type ( masses )  :: myMasses

     contains
         private
         procedure, private :: aA_fcn
         procedure, public  :: writer           => writer_sub
         procedure, public  :: extrema          => extrema_sub
         procedure, public  :: update_f         => update_f_sub
         procedure, public  :: neighbors        => neighbors_sub
         procedure, public  :: allocator        => allocator_sub
         procedure, public  :: thermalize       => thermalize_sub
         procedure, public  :: housekeeping     => housekeeping_sub
         procedure, public  :: compute_sigma    => compute_sigma_sub
         procedure, public  :: greens_two_point => greens_two_point_sub
    end type fields

    ! local variables
    character ( len = * ),   private, parameter :: error_fatal = 'Program halting in module mFields due to fatal error.'

    ! alphabetical
    private :: allocator_sub
    private :: compute_sigma_sub
    private :: extrema_sub
    private :: greens_two_point_sub
    private :: housekeeping_sub
    private :: neighbors_sub
    private :: thermalize_sub
    private :: update_f_sub
    private :: writer_sub
    private :: write_f_sub

    interface

        module subroutine allocator_sub ( me )
            class ( fields ), target :: me
        end subroutine allocator_sub

        module subroutine housekeeping_sub ( me )
            class ( fields ), target :: me
        end subroutine housekeeping_sub

        module subroutine extrema_sub ( me )
            class ( fields ), target :: me
        end subroutine extrema_sub

        module subroutine compute_sigma_sub ( me )
            class ( fields ), target :: me
        end subroutine compute_sigma_sub

        module subroutine neighbors_sub ( me )
            class ( fields ), target :: me
        end subroutine neighbors_sub

        module subroutine write_f_sub ( me, myInputs )
            class ( fields ), target        :: me
            type  ( inputs ), intent ( in ) :: myInputs
        end subroutine write_f_sub

        module subroutine writer_sub ( me, io_output_handle, myInputs )
            class ( fields ), target :: me
            type  ( inputs ), intent ( in ), target :: myInputs
            integer,          intent ( in )         :: io_output_handle
        end subroutine writer_sub

    end interface

contains
    ! aA_fcn
    ! greens_two_point_sub
    ! thermalize_sub
    ! update_f_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    function aA_fcn ( me, iu, ju, ku, lu, i, j, k, l ) result ( fcn_result )

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
     end function aA_fcn

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine greens_two_point_sub ( me )

        class ( fields ), target :: me

        ! locals
        integer ( ip ) :: i, j, k, l, iu, ju, ku, lu, luu, luuu

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

    end subroutine greens_two_point_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine thermalize_sub ( me, temp, farray )

        class ( fields ), target :: me

        character ( len = * ), intent ( in ) :: farray
        character ( len = * ), intent ( in ) :: temp
        ! locals
        real ( rp ) :: random
        integer ( ip ) :: i, j, k, l
        integer :: io_in_farray

            ! read in a thermalized array or heat up a cold array
            if ( temp == 'hot' ) then
                write ( stdout, fmt_generic ) 'The temp is ', temp
                io_in_farray = safeopen_readonly ( farray )
                read ( io_in_farray, * ) me % f
            else if ( temp == 'cold') then
                write ( stdout, fmt_generic ) 'The temp is ', temp
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
                    write ( stdout, fmt_generic ) 'Unrecognized temperature: should be "hot" or "cold". Input value was ', temp, '.'
                    stop 'Fatal error - I need to know the temperature.'
            end if

    end subroutine thermalize_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    subroutine update_f_sub ( me )

        class ( fields ), target :: me

        type ( extents ), pointer :: ex

        real ( rp ) :: oldf, oldP, newP, rdn
        real ( rp ) :: total

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
                                oldP = me % aA_fcn ( iu, ju, ku, lu,  i,  j,  k,  l ) &
                                     * me % aA_fcn ( iu, ju, ku,  l,  i,  j,  k, ld ) &
                                     * me % aA_fcn (  i, ju, ku, lu, id,  j,  k,  l ) &
                                     * me % aA_fcn ( iu,  j, ku, lu,  i, jd,  k,  l ) &
                                     * me % aA_fcn ( iu, ju,  k, lu,  i,  j, kd,  l )
                                call random_number ( rdn )
                                me % f ( i, j, k, l ) = oldf + ex % df * ( rdn - half )
                                newP = me % aA_fcn ( iu, ju, ku, lu,  i,  j,  k,  l ) &
                                     * me % aA_fcn ( iu, ju, ku,  l,  i,  j,  k, ld ) &
                                     * me % aA_fcn (  i, ju, ku, lu, id,  j,  k,  l ) &
                                     * me % aA_fcn ( iu,  j, ku, lu,  i, jd,  k,  l ) &
                                     * me % aA_fcn ( iu, ju,  k, lu,  i,  j, kd,  l )
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
                call greens_two_point_sub ( me )
                me % E_0 ( sweep ) = ex % avolume
            end do sweepdo

            ex => null ( )

            total = real ( me % naccept + me % nreject, rp )
            me % ratio_accept = real ( me % naccept, rp ) / total
            me % ratio_reject = real ( me % nreject, rp ) / total

    end subroutine update_f_sub

end module mFields
