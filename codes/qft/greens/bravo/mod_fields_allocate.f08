! 3456789 123456789 223456789 323456789 423456789 523456789 623456789 723456789 823456789 923456789 023456789 123456789 223456789 32
submodule ( mFields ) smFieldsAllocate

    use mConstants,                     only : zero, stdout

    integer                 :: alloc_status  = -1
    character ( len = 512 ) :: alloc_message = 'null'

    contains
        ! housekeeping_sub
        ! allocator_sub
        ! allocate_rank_4_rp_sub
        ! allocate_rank_3_rp_sub
        ! allocate_rank_2_rp_sub
        ! allocate_rank_1_rp_sub
        ! allocate_rank_1_ip_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    module subroutine housekeeping_sub ( me )

        class ( fields ), target :: me

            call allocator_sub ( me ) ! allocate all arrays
            call neighbors_sub ( me ) ! populate pointers to neighbors

            me % naccept = 0
            me % nreject = 0

            me % highphi  = zero
            me % highgphi = zero
            me % highdphi = zero

            me % myExtents % outoftable = zero

            me % G ( : ) = zero

    end subroutine housekeeping_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    module subroutine allocator_sub ( me )

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
            call allocate_rank_1_rp_sub ( me % phi,         me % myExtents % Nphi + 1 )
            call allocate_rank_1_rp_sub ( me % gphi,        me % myExtents % Nphi + 1 )
            call allocate_rank_1_rp_sub ( me % dphi,        me % myExtents % Nphi + 1 )
            call allocate_rank_1_rp_sub ( me % E_0,         me % myExtents % Nsweeps )
            ! pbc
            call allocate_rank_1_ip_sub ( me % ups,  me % myExtents % Ns )
            call allocate_rank_1_ip_sub ( me % dns,  me % myExtents % Ns )
            call allocate_rank_1_ip_sub ( me % upt,  me % myExtents % Nt )
            call allocate_rank_1_ip_sub ( me % dnt,  me % myExtents % Nt )

    end subroutine allocator_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    module subroutine allocate_rank_4_rp_sub ( array, length1, length2, length3, length4 )

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

            array ( : , : , : , : ) = zero

        return

    100 format ( * ( g0 ) )

    end subroutine allocate_rank_4_rp_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    module subroutine allocate_rank_3_rp_sub ( array, length1, length2, length3 )

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

            array ( : , : , : ) = zero

        return

    100 format ( * ( g0 ) )

    end subroutine allocate_rank_3_rp_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    module subroutine allocate_rank_1_rp_sub ( array, length )

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

            array ( : ) = zero

        return

    100 format ( * ( g0 ) )

    end subroutine allocate_rank_1_rp_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    module subroutine allocate_rank_1_ip_sub ( array, length )

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

            array ( : ) = 0

        return

    100 format ( * ( g0 ) )

    end subroutine allocate_rank_1_ip_sub

end submodule smFieldsAllocate
