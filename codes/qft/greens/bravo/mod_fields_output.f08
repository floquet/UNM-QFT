! 3456789 123456789 223456789 323456789 423456789 523456789 623456789 723456789 823456789 923456789 023456789 123456789 223456789 32
submodule ( mFields ) smFieldsOutput

    use mConstants,                     only : stdout, biggest
    use mFileHandling,                  only : safeopen_readonly, find_IU_info

    contains
        ! write_f_sub
        ! writer_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    module subroutine write_f_sub ( me, myInputs )

        class ( fields ), target :: me

        type  ( inputs ), intent ( in ) :: myInputs
        ! locals
        integer :: io_output, io_stat
        character ( len = 256 ) :: io_msg

            io_output = safeopen_readonly ( myInputs % farray )
            write ( io_output, *, iostat = io_stat, iomsg = io_msg ) me % f
            if ( io_stat /= 0 ) then
                write ( stdout, fmt_generic ) 'Error attempting to write the rank 4 array f:'
                write ( stdout, 110 ) 'shape ( f ) =', shape ( me % f )
                write ( stdout, fmt_generic ) 'number of elements in f = ', product ( shape ( me % f ) )
                write ( stdout, fmt_generic ) 'iomsg = ', trim ( io_msg ), '.'
                write ( stdout, fmt_generic ) 'iostat = ', io_stat
                call find_IU_info ( io_output )
                write ( stdout, fmt_generic ) 'Execution nervously continues...'
            endif
            close ( io_output )

        return
        110 format ( * ( g0, ' ' ) )

    end subroutine write_f_sub

    !  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    module subroutine writer_sub ( me, io_output_handle, myInputs )

        class ( fields ), target :: me

        type  ( inputs ), intent ( in ), target :: myInputs
        integer,          intent ( in )         :: io_output_handle
        !locals
        type ( extents ), pointer :: ex
        type ( inputs ),  pointer :: in
        integer :: i, j, k

            ex => me % myExtents
            in => myInputs

            write ( io_output_handle, fmt_generic ) 'maxPhi,    maxGphi,    maxDphi'
            write ( io_output_handle, 200 ) ex % maxPhi, ex % maxGphi, ex % maxDphi

            write ( io_output_handle, fmt_generic ) 'Nphi,   Ngphi,    Ndphi'
            write ( io_output_handle, 210 ) ex % Nphi, ex % Ngphi, ex % Ndphi

            write ( io_output_handle, fmt_generic ) 'as,   at,   Mass,   m,   df'
            write ( io_output_handle, 220 ) ex % as, ex % at, me % myMasses % Mass, me % myMasses % m, ex % df

            write ( io_output_handle, fmt_generic ) 'tablename,     temp,    root,    farray,    index'
            write ( io_output_handle, 230 ) trim ( in % tablename ), ' hot ', &
                                            trim ( in % root ),      '  ',    &
                                            trim ( in % farray ),    in % index + 1

            write ( io_output_handle, fmt_generic ) 'Nsweeps,   Ns,   Nt'
            write ( io_output_handle, 210 ) ex % Nsweeps, ex % Ns, ex % Nt

            write ( io_output_handle, fmt_generic ) 'Results from run ', in % index
            write ( io_output_handle, fmt_generic ) 'E_0 = ', me % meanE, ', sigma = ', me % sigma

            do k = 0, 3
                write ( io_output_handle, fmt_generic ) 'G( ', k, ' ) = ', me % G ( k ) / real ( ex % kount, rp )
            end do

            write ( io_output_handle, fmt_generic ) 'minimum of A = ', me % A_min
            write ( io_output_handle, fmt_generic ) 'maximum of A = ', me % A_max

            write ( io_output_handle, fmt_generic ) 'minimum of C = ', me % C_min
            write ( io_output_handle, fmt_generic ) 'maximum of C = ', me % C_max

            write ( io_output_handle, fmt_generic ) 'highphi  = ', me % highphi
            write ( io_output_handle, fmt_generic ) 'highgphi = ', me % highgphi
            write ( io_output_handle, fmt_generic ) 'highdphi = ', me % highdphi

            i = int ( me % highphi  / ex % phistep  ) + 1
            j = int ( me % highgphi / ex % gphistep ) + 1
            k = int ( me % highdphi / ex % dphistep ) + 1

            write ( io_output_handle, fmt_generic ) 'A ( highphi = ', i, ', highgphi = ', j, ', highdphi = ', k, ' ) = ', &
                                                     me % A ( i, j, k )

            write ( io_output_handle, fmt_generic ) 'out of table = ', ex % outoftable
            write ( io_output_handle, fmt_generic ) 'minA = ', biggest

            ex => null ( )
            in => null ( )

        return

    200 format ( 3 ( f12.3 ) )
    210 format ( 3 ( I10 ) )
    220 format ( 5 ( f12.6 ) )
    230 format ( 5 ( a ), I10 )

    end subroutine writer_sub

end submodule smFieldsOutput
