! 3456789 123456789 223456789 323456789 423456789 523456789 623456789 723456789 823456789 23456789 023456789 123456789 223456789 32
program makesAC

    use, intrinsic :: iso_fortran_env,  only : compiler_version, compiler_options, stdin => input_unit, stdout => output_unit

    use mSetPrecision,                  only : ip, rp
    use mFileHandling,                  only : safeopen_writereplace
    use mMPI_SGI,                       only : MPI_COMM_WORLD  ! SGI

    implicit none

    external :: MPI_INIT, MPI_COMM_RANK, MPI_FINALIZE

    character ( len =  50 ) :: tablename = ''
    character ( len =  55 ) :: proctablename = ''
    character ( len =   5 ) :: och5 = ''
    character ( len = 512 ) :: io_msg = ''

    ! parameters
    real ( rp ), parameter :: zero = 0.0_rp, one = 1.0_rp, two = 2.0_rp
    real ( rp ), parameter :: half = 0.5_rp, athird = one / 3.0_rp, tthird = two / 3.0_rp, sixth = half * athird

    ! real ( rp ), parameter :: twopi = 2.0_rp * acos ( -one )
    real ( rp ), parameter :: eps = 10.0_rp * epsilon ( one )

    ! variables
    integer( ip ) :: i = 0, j = 0, k = 0, l = 0, ijk = 0, Npi = 0, Nphi = 0, Ndphi = 0, Ngphi = 0
    integer( ip ) :: Nd = 0, Ng = 0, N = 0, Npoints = 0, thisNpi = 0
    integer       :: ierror = 0, proc = 0, numproc = 0, io_status = 0, io_proc_table = 0

    real ( rp ) :: at = zero, as = zero, pi = zero, phi = zero, dphi = zero, gphi = zero, Mass = zero, m = zero
    real ( rp ) :: pistep = zero, phistep = zero, dphistep = zero, gphistep = zero
    real ( rp ) :: maxPhi = zero, maxDphi = zero, maxGphi = zero, maxPi = zero, thismaxPi = zero
    real ( rp ) :: A = zero, C = zero

    real ( rp ) :: Amax = zero, Cmax = zero, Amin = zero, Cmin = zero
    real ( rp ) :: Amax0 = zero, Cmax0 = zero, Amin0 = zero, Cmin0 = zero
    real ( rp ) :: timei = zero, timef = zero, percent = zero

        call mpi_init ( ierror )
        call mpi_comm_rank ( mpi_comm_world, proc,    ierror )
        call mpi_comm_size ( mpi_comm_world, numproc, ierror )

        timei = zero
        timef = zero
        if ( proc == 0 ) timei = mpi_wtime()

        Amax = real (-1.e6, rp )
        Cmax = real (-1.e6, rp )
        Amin = real ( 1.e6, rp )
        Cmin = real ( 1.e6, rp )

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

        call mpi_bcast ( Nphi,  1, MPI_LONG, 0, mpi_comm_world, ierror )
        call mpi_bcast ( Ndphi, 1, MPI_LONG, 0, mpi_comm_world, ierror )
        call mpi_bcast ( Ngphi, 1, MPI_LONG, 0, mpi_comm_world, ierror )
        call mpi_bcast ( Npi,   1, MPI_LONG, 0, mpi_comm_world, ierror )

        call mpi_bcast ( at,   1, MPI_DOUBLE, 0, mpi_comm_world, ierror )
        call mpi_bcast ( as,   1, MPI_DOUBLE, 0, mpi_comm_world, ierror )
        call mpi_bcast ( Mass, 1, MPI_DOUBLE, 0, mpi_comm_world, ierror )
        call mpi_bcast ( m,    1, MPI_DOUBLE, 0, mpi_comm_world, ierror )

        call mpi_bcast ( tablename, 50, MPI_CHAR, 0, mpi_comm_world, ierror )

        call mpi_barrier ( mpi_comm_world, ierror )

        write ( och5, '( I5.5 )' ) proc

        proctablename = trim ( tablename ) // och5

        call mpi_barrier ( mpi_comm_world, ierror )

        ! ( newunit = io_proc_table, file = trim ( proctablename ), iostat = io_status, iomsg = io_msg )
        io_proc_table = safeopen_writereplace ( trim ( proctablename ) )

        phistep  = maxPhi  /real ( Nphi,  rp )
        gphistep = maxGphi /real ( Ngphi, rp )
        dphistep = maxDphi /real ( Ndphi, rp )
        pistep   = maxPi   /real ( Npi,   rp )
        Nd = Ndphi + 1
        Ng = Ngphi + 1
        N  = Nphi  + 1

        Npoints = Nd * Ng * N

        percent = 5.0_rp

        do ijk = proc+1, Npoints, numproc

            k = (ijk - 1)/(Ng * N) + 1
            j = (ijk - 1 - (k - 1) * Ng * N ) / N + 1
            i = ijk - N*(j-1 + Ng*(k-1))

            dphi = real(k-1, rp )*dphistep
            gphi = real(j-1, rp )*gphistep
            phi  = real(i-1, rp )*phistep
            pistep = maxPi/real(Npi, rp )

            thismaxPi = min(maxPi, sqrt(Mass**4 + (gphi/as)**2 + (m*phi)**2))
            thisNpi = int(thismaxPi/pistep, rp )

            pi   = zero

            A = - sixth*pistep*cos(as**3*F(phi, gphi, dphi, pi)) &
          *exp(-at*as**3*H(phi,gphi,pi))*D(phi, gphi, pi)
          C = - sixth*pistep*cos(as**3*F(phi, gphi, dphi, pi)) &
          *exp(-at*as**3*H(phi,gphi,pi))*H(phi,gphi,pi) &
          *D(phi, gphi, pi)

     pido: do l = 1, thisNpi

        A = A + &
             athird*pistep*cos(as**3*F(phi, gphi, dphi, pi)) &
             *exp(-at*as**3*H(phi,gphi,pi))*D(phi, gphi, pi) &
             + tthird*pistep*cos(as**3*F(phi, gphi, dphi, pi+half*pistep)) &
             *exp(-at*as**3*H(phi,gphi,pi+half*pistep)) &
             *D(phi, gphi, pi+half*pistep)

        C = C + &
             athird*pistep*cos(as**3*F(phi, gphi, dphi, pi)) &
             *exp(-at*as**3*H(phi,gphi,pi))*H(phi,gphi,pi) &
             *D(phi, gphi, pi) &
             + tthird*pistep*cos(as**3*F(phi, gphi, dphi, pi+half*pistep)) &
             *exp(-at*as**3*H(phi,gphi,pi+half*pistep)) &
             *H(phi,gphi,pi+half*pistep)*D(phi, gphi, pi+half*pistep)

        if ( pistep*exp(-at*as**3*H(phi,gphi,pi))*D(phi, gphi, pi) < eps ) then
           exit
        end if

        pi = pi + pistep

     end do pido

     write (io_proc_table,'(3(f16.5),2(es20.11))') phi,gphi,dphi,A,C

     Amax = max( A, Amax)
     Cmax = max(C,Cmax)
     Amin = min( A, Amin)
     Cmin = min(C,Cmin)

     if (proc == 0 .and. 100.0*real(ijk)/real(Npoints) > percent) then
        write(6,'(f10.4,a)') percent,' % done'
        percent = percent + 5.0
     end if

  end do

  close (io_proc_table)

  call mpi_barrier(mpi_comm_world, ierror )

  call mpi_reduce( Amin, Amin0, 1, MPI_DOUBLE,MPI_MIN, 0, mpi_comm_world, ierror )
  call mpi_reduce(Cmin,Cmin0, 1, MPI_DOUBLE,MPI_MIN, 0, mpi_comm_world, ierror )
  call mpi_reduce( Amax, Amax0, 1, MPI_DOUBLE,MPI_MAX, 0, mpi_comm_world, ierror )
  call mpi_reduce(Cmax,Cmax0, 1, MPI_DOUBLE,MPI_MAX, 0, mpi_comm_world, ierror )

  call mpi_barrier(mpi_comm_world, ierror )

  if (proc == 0) then
     write(6,*) 'maxPhi, maxGphi, maxDphi, maxPi'
     write(6,*) maxPhi,maxGphi,  maxDphi, maxPi
     write(6,*) 'Nphi, Ngphi,  Ndphi, Npi'
     write(6,*) Nphi, Ngphi, Ndphi, Npi
     write(6,*) 'at, as, Mass, m'
     write(6,*) at, as, Mass, m
     write(6,*) 'tablename'
     write(6,*) tablename
     write(6,*) 'minimum of A =', Amin0
     write(6,*) 'maximum of A =', Amax0
     write(6,*) 'minimum of C =', Cmin0
     write(6,*) 'maximum of C =', Cmax0
     timef = mpi_wtime()
     print *, timef - timei," seconds"
  end if

  call mpi_finalize(ierror)

contains
  function H(phi, gphi, pi)
    real ( rp ) :: H, phi, gphi, pi
    H = ( two*pi**2 - (gphi/as)**2 - (m*phi)**2 - Mass**4 ) &
         /( ( one - ( pi**2 - (gphi/as)**2 - (m*phi)**2 )/Mass**4 ) &
         * sqrt( one - ( pi**2 - (gphi/as)**2 - (m*phi)**2 )/Mass**4 ) ) &
         + Mass**4
  end function H
  function D(phi, gphi, pi)
    real ( rp ) :: D, phi, gphi, pi
    D = ( one + ( two*pi**2 +  (gphi/as)**2 + (m*phi)**2 )/Mass**4 ) &
         /( ( one - ( pi**2 -  (gphi/as)**2 - (m*phi)**2 )/Mass**4 )**2 &
         *sqrt( ( one - ( pi**2 -  (gphi/as)**2 - (m*phi)**2 )/Mass**4 ) ) )
  end function D
  function F(phi, gphi, dphi, pi)
    real ( rp ) :: F, phi, gphi, dphi, pi
    F = dphi*pi  &
         /( ( one - ( pi**2 -  (gphi/as)**2 - (m*phi)**2 )/Mass**4 ) &
         *sqrt( ( one - ( pi**2 -  (gphi/as)**2 - (m*phi)**2 )/Mass**4 ) ) )
  end function F
end program makesAC
