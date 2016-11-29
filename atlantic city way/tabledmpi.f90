program tabledmpi ! keep Nphi > numproc
     use, intrinsic :: iso_fortran_env, only : compiler_version, compiler_options
 implicit none

  include 'mpif.h'

  character(len=50):: tablename
  character(len=55):: proctablename

  character(len=1)::och1
  character(len=2)::och2
  character(len=3)::och3
  character(len=4)::och4
  character(len=5)::och5

  integer(kind=8):: i, j, k, l, Npi, Nphi, Ndphi, Ngphi
  real(kind=8):: at, a, pi, phi, dphi, gphi, Mass, m
  real(kind=8):: pistep, phistep, dphistep, gphistep
  real(kind=8):: maxPhi, maxDphi, maxGphi, maxPi
  real(kind=8),parameter:: zero=0.0d0, one=1.0d0, half=0.5d0
  real(kind=8),parameter:: athird=one/3.0d0, tthird=2.0d0/3.0d0
  real(kind=8),parameter:: twopi = 8.0d0*atan(1.0d0), sixth=one/6.0d0
  real(kind=8):: Integral, HIntegral, timei, timef
  real(kind=8),parameter:: eps=10.d0*epsilon(1.0d0)

  integer(kind=4):: ierror,proc,numproc
  real(kind=8):: Intmax,Hintmax,Intmin,Hintmin
  real(kind=8):: Intmax0,Hintmax0,Intmin0,Hintmin0

  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world,proc,ierror)
  call mpi_comm_size(mpi_comm_world,numproc,ierror)

  if( proc == 0) timei = mpi_wtime()

   Intmax = real(-1.e6,kind=8)
  Hintmax = real(-1.e6,kind=8)
   Intmin = real( 1.e6,kind=8)
  Hintmin = real( 1.e6,kind=8)

  if (proc == 0) then
     read(5,*); read(5,*) maxPhi, maxDphi, maxGphi, maxPi
     read(5,*); read(5,*) Nphi, Ndphi, Ngphi, Npi
     read(5,*); read(5,*) at, a, Mass, m
     read(5,*); read(5,*) tablename
  end if

  call mpi_barrier(mpi_comm_world,ierror)

  call mpi_bcast(maxPhi ,1,MPI_DOUBLE,0,mpi_comm_world,ierror)
  call mpi_bcast(maxDphi,1,MPI_DOUBLE,0,mpi_comm_world,ierror)
  call mpi_bcast(maxGphi,1,MPI_DOUBLE,0,mpi_comm_world,ierror)
  call mpi_bcast(maxPi  ,1,MPI_DOUBLE,0,mpi_comm_world,ierror)

  call mpi_bcast(Nphi ,1,MPI_LONG,0,mpi_comm_world,ierror)
  call mpi_bcast(Ndphi,1,MPI_LONG,0,mpi_comm_world,ierror)
  call mpi_bcast(Ngphi,1,MPI_LONG,0,mpi_comm_world,ierror)
  call mpi_bcast(Npi  ,1,MPI_LONG,0,mpi_comm_world,ierror)

  call mpi_bcast(at  ,1,MPI_DOUBLE,0,mpi_comm_world,ierror)
  call mpi_bcast(a   ,1,MPI_DOUBLE,0,mpi_comm_world,ierror)
  call mpi_bcast(Mass,1,MPI_DOUBLE,0,mpi_comm_world,ierror)
  call mpi_bcast(m   ,1,MPI_DOUBLE,0,mpi_comm_world,ierror)

  call mpi_bcast(tablename,50,MPI_CHAR,0,mpi_comm_world,ierror)

  call mpi_barrier(mpi_comm_world,ierror)

  if (proc < 10) then
      write (och1,'(I1)') proc
      och5 = '0000'//och1
  else if (proc >= 10 .and. proc < 100) then
      write (och2,'(I2)') proc
      och5 = '000'//och2
  else if (proc >= 100 .and. proc < 1000) then
      write (och3,'(I3)') proc
      och5 = '00'//och3
  else if (proc >= 1000 .and. proc < 10000) then
      write (och4,'(I4)') proc
      och5 = '0'//och4
  else
      write (och5,'(I5)') proc
  end if

  proctablename = trim(tablename)//och5

  call mpi_barrier(mpi_comm_world,ierror)

  open(45,file=trim(proctablename))

  phistep = maxPhi/real(Nphi,kind=8)
  gphistep = maxGphi/real(Ngphi,kind=8)
  dphistep = maxDphi/real(Ndphi,kind=8)
  pistep = maxPi/real(Npi,kind=8)

  dphi = zero
  dphido: do k = 1, Ndphi + 1
    gphi = zero
    gphido: do j= 1, Ngphi + 1
      phido: do i = proc+1, Nphi + 1, numproc
        phi = real(i-1,kind=8)*phistep
        pi = zero

        Integral = - sixth*pistep*cos(a**3*dphi*pi) &
                    *exp(-at*a**3*H(phi,gphi,pi))
        HIntegral = - sixth*pistep*cos(a**3*dphi*pi) &
                    *exp(-at*a**3*H(phi,gphi,pi))*H(phi,gphi,pi)

        pido: do l = 1, Npi + 1

           Integral = Integral + &
                      athird*pistep*cos(a**3*dphi*pi) &
                      *exp(-at*a**3*H(phi,gphi,pi)) + &
                      tthird*pistep*cos(a**3*dphi*(pi+half*pistep)) &
                      *exp(-at*a**3*H(phi,gphi,pi+half*pistep))

           HIntegral = HIntegral + &
                       athird*pistep*cos(a**3*dphi*pi) &
                       *exp(-at*a**3*H(phi,gphi,pi))*H(phi,gphi,pi) + &
                       tthird*pistep*cos(a**3*dphi*(pi+half*pistep)) &
                       *exp(-at*a**3*H(phi,gphi,pi+half*pistep)) &
                       *H(phi,gphi,pi+half*pistep)

          pi = pi + pistep

          if (pistep*exp(-at*a**3*H(phi,gphi,pi+half*pistep)) < eps ) then
            exit
          end if

        end do pido

        write (45,'(3(f12.5),2(es20.11))') phi,gphi,dphi,Integral,Hintegral

         Intmax = max( Integral, Intmax)
        HIntmax = max(HIntegral,HIntmax)
         Intmin = min( Integral, Intmin)
        HIntmin = min(HIntegral,HIntmin)

      end do phido
      gphi = gphi + gphistep
    end do gphido
    if (proc == 0) write(6,'(f10.4,a)') real(100*k)/real(Ndphi),' % done'
    dphi = dphi + dphistep
  end do dphido

  close (45)

  call mpi_barrier(mpi_comm_world,ierror)

  call mpi_reduce( Intmin, Intmin0,1,MPI_DOUBLE,MPI_MIN,0,mpi_comm_world,ierror)
  call mpi_reduce(Hintmin,Hintmin0,1,MPI_DOUBLE,MPI_MIN,0,mpi_comm_world,ierror)
  call mpi_reduce( Intmax, Intmax0,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm_world,ierror)
  call mpi_reduce(Hintmax,Hintmax0,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm_world,ierror)

  call mpi_barrier(mpi_comm_world,ierror)

  if (proc == 0) then
     write(6,*) 'maxPhi, maxDphi, maxGphi, maxPi'
     write(6,*) maxPhi, maxDphi, maxGphi, maxPi
     write(6,*) 'Nphi, Ndphi, Ngphi, Npi'
     write(6,*) Nphi, Ndphi, Ngphi, Npi
     write(6,*) 'at, a, Mass, m'
     write(6,*) at, a, Mass, m
     write(6,*) 'tablename'
     write(6,*) tablename
     write(6,*) 'minimum of Integral =', Intmin0
     write(6,*) 'maximum of Integral =', Intmax0
     write(6,*) 'minimum of HIntegral =', Hintmin0
     write(6,*) 'maximum of HIntegral =', Hintmax0
  end if
  if( proc == 0) then
      timef = mpi_wtime()
      print *, timef - timei, 'seconds'
      write ( *, '( "Fortran compiler version: ", g0 )'    ) compiler_version ( )
      write ( *, '( "Fortran compiler options: ", g0, / )' ) compiler_options ( )
  endif
  call mpi_finalize(ierror)

contains
  function H(phi, gphi, pi)
    real(kind=8):: H, phi, gphi, pi
    H = sqrt((Mass**4 + pi**2)*(Mass**4 + gphi**2 + (m*phi)**2)) - Mass**4
  end function H
end program tabledmpi
