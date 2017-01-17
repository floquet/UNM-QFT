! 3456789 123456789 223456789 323456789 423456789 523456789 623456789 723456789 823456789 923456789 023456789 123456789 223456789 32
program makesAC

    use, intrinsic :: iso_fortran_env,  only : compiler_version, compiler_options, stdin => input_unit, stdout => output_unit

    use mFileHandling,                  only : safeopen_writereplace
    use mMPI_SGI,                       only : mpi_comm_world, mpi_wtime, mpi_char, mpi_double, mpi_long, mpi_max, mpi_min  ! SGI
    use mLibraryFunctions,              only : H, D, F, X
    use mSetPrecision,                  only : ip, rp

    implicit none

    external :: MPI_INIT, MPI_COMM_RANK, MPI_FINALIZE, MPI_COMM_SIZE, MPI_BCAST, MPI_BARRIER, MPI_REDUCE!, MPI_Wtime

    ! parameters
    real ( rp ), parameter :: zero = 0.0_rp, one = 1.0_rp, two = 2.0_rp
    real ( rp ), parameter :: half = 0.5_rp, athird = one / 3.0_rp, tthird = two * athird, sixth = half * athird

    real ( rp ), parameter :: eps = 10.0_rp * epsilon ( one ), threshold = huge ( one )

    ! variables
    real ( rp ) :: at = zero, as = zero, pi = zero, phi = zero, dphi = zero, gphi = zero, Mass = zero, m = zero
    real ( rp ) :: pistep = zero, phistep = zero, dphistep = zero, gphistep = zero
    real ( rp ) :: maxPhi = zero, maxDphi = zero, maxGphi = zero, maxPi = zero, thismaxPi = zero
    real ( rp ) :: A = zero, C = zero, test = zero

    real ( rp ) :: Amax  = zero, Cmax  = zero, Amin  = zero, Cmin  = zero
    real ( rp ) :: Amax0 = zero, Cmax0 = zero, Amin0 = zero, Cmin0 = zero
    real ( rp ) :: timei = zero, timef = zero, percent = zero

    integer ( ip ) :: i = 0, j = 0, k = 0, l = 0, ijk = 0, Npi = 0, Nphi = 0, Ndphi = 0, Ngphi = 0
    integer ( ip ) :: Nd = 0, Ng = 0, N = 0, Npoints = 0, thisNpi = 0
    integer ( ip ) :: proc = 0, numproc
    integer        :: ierror = 0, myproc = 0, mynumproc = 0, io_proc_table = 0

    character ( len =  50 ) :: tablename = ''
    character ( len =  55 ) :: proctablename = ''
    character ( len =   5 ) :: och5 = ''

        call mpi_init ( ierror )
        call mpi_comm_rank ( mpi_comm_world, myproc,    ierror )
        call mpi_comm_size ( mpi_comm_world, mynumproc, ierror )

        ! convert to long integer for integration loop: allows Npoints > 2**16
        proc = myproc
        numproc = mynumproc
        if ( proc == 0 ) timei = mpi_wtime ()

        Amax = -threshold
        Cmax = -threshold
        Amin =  threshold
        Cmin =  threshold

        if ( proc == 0 ) then
            read ( stdin, * ); read ( stdin, * ) maxPhi, maxGphi, maxDphi, maxPi
            read ( stdin, * ); read ( stdin, * ) Nphi, Ndphi, Ngphi, Npi
            read ( stdin, * ); read ( stdin, * ) as, at, Mass, m
            read ( stdin, * ); read ( stdin, * ) tablename
        end if

        call mpi_barrier ( mpi_comm_world, ierror )

        call mpi_bcast ( maxPhi , 1, MPI_DOUBLE, 0, mpi_comm_world, ierror )
        call mpi_bcast ( maxDphi, 1, MPI_DOUBLE, 0, mpi_comm_world, ierror )
        call mpi_bcast ( maxGphi, 1, MPI_DOUBLE, 0, mpi_comm_world, ierror )
        call mpi_bcast ( maxPi  , 1, MPI_DOUBLE, 0, mpi_comm_world, ierror )

        call mpi_bcast ( Nphi,    1, MPI_LONG,   0, mpi_comm_world, ierror )
        call mpi_bcast ( Ndphi,   1, MPI_LONG,   0, mpi_comm_world, ierror )
        call mpi_bcast ( Ngphi,   1, MPI_LONG,   0, mpi_comm_world, ierror )
        call mpi_bcast ( Npi,     1, MPI_LONG,   0, mpi_comm_world, ierror )

        call mpi_bcast ( at,      1, MPI_DOUBLE, 0, mpi_comm_world, ierror )
        call mpi_bcast ( as,      1, MPI_DOUBLE, 0, mpi_comm_world, ierror )
        call mpi_bcast ( Mass,    1, MPI_DOUBLE, 0, mpi_comm_world, ierror )
        call mpi_bcast ( m,       1, MPI_DOUBLE, 0, mpi_comm_world, ierror )

        call mpi_bcast ( tablename, 50, MPI_CHAR, 0, mpi_comm_world, ierror )

        call mpi_barrier ( mpi_comm_world, ierror )

        write ( och5, '( I5.5 )' ) proc

        proctablename = trim ( tablename ) // och5

        call mpi_barrier ( mpi_comm_world, ierror )

        io_proc_table = safeopen_writereplace ( trim ( proctablename ) )

        phistep  = maxPhi  / real ( Nphi,  rp )
        gphistep = maxGphi / real ( Ngphi, rp )
        dphistep = maxDphi / real ( Ndphi, rp )
        pistep   = maxPi   / real ( Npi,   rp )

        Nd = Ndphi + 1
        Ng = Ngphi + 1
        N  = Nphi  + 1

        Npoints = Nd * Ng * N

        percent = 5.0_rp

        do ijk = proc + 1, Npoints, numproc

            k = ( ijk - 1 ) / ( Ng * N ) + 1
            j = ( ijk - 1 - ( k - 1 ) * Ng * N ) / N + 1
            i = ijk - N * ( j - 1 + Ng * ( k - 1 ) )

            dphi = real ( k - 1, rp ) * dphistep
            gphi = real ( j - 1, rp ) * gphistep
            phi  = real ( i - 1, rp ) *  phistep
            pistep = maxPi / real ( Npi, rp )

            thismaxPi = min ( maxPi, sqrt ( Mass**4 + X ( phi, gphi, as, m ) ) )
            thisNpi   = int ( thismaxPi / pistep, rp )

            pi   = zero

            A = - sixth * pistep * cos ( as**3 * F ( phi, gphi, dphi, pi, as, m, Mass ) )
            A =  A * exp ( -at * as**3 * H ( phi, gphi, pi, as, m, Mass ) ) * D ( phi, gphi, pi, as, m, Mass )

            C = A * H ( phi, gphi, pi, as, m, Mass )

            pido: do l = 1, thisNpi

                A = A + &
                    athird * pistep * cos ( as ** 3 * F ( phi, gphi, dphi, pi, as, m, Mass ) ) &
                    * exp ( -at * as ** 3 * H ( phi, gphi, pi, as, m, Mass ) ) &
                    * D ( phi, gphi, pi, as, m, Mass )
                    A = A + tthird * pistep * cos ( as ** 3 * F ( phi, gphi, dphi, pi + half * pistep, as, m, Mass ) ) &
                    * exp ( -at * as ** 3 * H ( phi, gphi, pi + half * pistep, as, m, Mass ) ) &
                    * D ( phi, gphi, pi + half * pistep, as, m, Mass )

                C = C + &
                    athird * pistep * cos ( as ** 3 * F ( phi, gphi, dphi, pi, as, m, Mass ) ) &
                    * exp ( -at * as ** 3 * H ( phi, gphi, pi, as, m, Mass) ) &
                    * H ( phi, gphi, pi, as, m, Mass ) &
                    * D ( phi, gphi, pi, as, m, Mass )
                C = C + tthird * pistep * cos ( as ** 3 * F ( phi, gphi, dphi, pi + half * pistep, as, m, Mass ) ) &
                    * exp ( -at * as ** 3 * H ( phi, gphi, pi + half * pistep, as, m, Mass ) ) &
                    * H ( phi, gphi, pi + half * pistep, as, m, Mass ) &
                    * D ( phi, gphi, pi + half * pistep, as, m, Mass )

                test = pistep * exp ( -at * as ** 3 * H ( phi, gphi, pi, as, m, Mass ) ) * D ( phi, gphi, pi, as, m, Mass )
                if ( test < eps ) exit

                pi = pi + real ( l - 1, rp ) * pistep ! avoid accumlation errors

            end do pido

            write ( io_proc_table, '( 3 (f16.5), 2 (es20.11) )' ) phi, gphi, dphi, A, C

            Amax = max ( A, Amax )
            Cmax = max ( C, Cmax )
            Amin = min ( A, Amin )
            Cmin = min ( C, Cmin )

            if ( proc == 0 .and. 100.0_rp * real ( ijk, rp ) / real ( Npoints, rp ) > percent ) then
                write ( stdout, '( f10.4, g0 )' ) percent, ' % done'
                percent = percent + 5.0_rp
            end if

        end do

        close ( io_proc_table )

        call mpi_barrier ( mpi_comm_world, ierror )

        call mpi_reduce ( Amin, Amin0, 1, MPI_DOUBLE, MPI_MIN, 0, mpi_comm_world, ierror )
        call mpi_reduce ( Cmin, Cmin0, 1, MPI_DOUBLE, MPI_MIN, 0, mpi_comm_world, ierror )
        call mpi_reduce ( Amax, Amax0, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm_world, ierror )
        call mpi_reduce ( Cmax, Cmax0, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm_world, ierror )

        call mpi_barrier(mpi_comm_world, ierror )

        if (proc == 0) then
            write ( stdout, * ) 'maxPhi, maxGphi, maxDphi, maxPi'
            write ( stdout, * ) maxPhi,maxGphi,  maxDphi, maxPi
            write ( stdout, * ) 'Nphi, Ngphi,  Ndphi, Npi'
            write ( stdout, * ) Nphi, Ngphi, Ndphi, Npi
            write ( stdout, * ) 'at, as, m, Mass, m'
            write ( stdout, * ) at, as, m, Mass, m
            write ( stdout, * ) 'tablename'
            write ( stdout, * ) tablename
            write ( stdout, * ) 'minimum of A =', Amin0
            write ( stdout, * ) 'maximum of A =', Amax0
            write ( stdout, * ) 'minimum of C =', Cmin0
            write ( stdout, * ) 'maximum of C =', Cmax0
            write ( stdout, '( /, "Fortran compiler version:    ", g0    )' ) compiler_version ( )
            write ( stdout, '(    "Fortran compilation options: ", g0, / )' ) compiler_options ( )
            timef = mpi_wtime()
            write ( stdout, '( /, g0, ": cpu time, seconds")' ) timef - timei
        end if

        call mpi_finalize ( ierror )

end program makesAC
