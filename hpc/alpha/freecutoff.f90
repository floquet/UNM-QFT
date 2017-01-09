program freecutoff ! free field
  implicit none
  ! gphi = sqrt( (f(iu,j,k,l) - f(i,j,k,l))**2 + ... )     (without /as**2)
  ! order is phi, gphi, dphi, pi
  ! |phi| < Mass
  include 'mpif.h'

  character(len=50):: tablename, infile
  character(len=55):: proctablename

  character(len=5)::och5

  integer(kind=8):: i, j, k, l, ijk, Npi, Nphi, Ndphi, Ngphi
  integer(kind=8):: Nd, Ng, N, Npoints
  real(kind=8):: at, as, pi, phi, dphi, gphi, Mass, m
  real(kind=8):: pistep, phistep, dphistep, gphistep
  real(kind=8):: maxPhi, maxGphi, maxDphi, maxPi
  real(kind=8),parameter:: zero=0.0d0, one=1.0d0, half=0.5d0
  real(kind=8),parameter:: athird=one/3.0d0, tthird=2.0d0/3.0d0
  real(kind=8),parameter:: twopi = 8.0d0*atan(1.0d0), sixth=one/6.0d0
  real(kind=8):: A, C
  real(kind=8),parameter:: eps=10.d0*epsilon(1.0d0)

  integer(kind=4):: ierror,proc,numproc
  real(kind=8):: Amax,Cmax,Amin,Cmin
  real(kind=8):: Amax0,Cmax0,Amin0,Cmin0
  real(kind=8):: timei, timef, percent

  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world,proc,ierror)
  call mpi_comm_size(mpi_comm_world,numproc,ierror)

  timei = zero
  timef = zero
  if (proc == 0) timei = mpi_wtime()

  Amax = real(-1.e6,kind=8)
  Cmax = real(-1.e6,kind=8)
  Amin = real( 1.e6,kind=8)
  Cmin = real( 1.e6,kind=8)

  if (proc == 0) then
     call getarg(1,infile)
     open(5,file=infile)
     read(5,*); read(5,*) maxPhi, maxGphi, maxDphi, maxPi
     read(5,*); read(5,*) Nphi, Ngphi, Ndphi, Npi
     read(5,*); read(5,*) as, at, Mass, m
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
  call mpi_bcast(as   ,1,MPI_DOUBLE,0,mpi_comm_world,ierror)
  call mpi_bcast(Mass,1,MPI_DOUBLE,0,mpi_comm_world,ierror)
  call mpi_bcast(m   ,1,MPI_DOUBLE,0,mpi_comm_world,ierror)

  call mpi_bcast(tablename,50,MPI_CHAR,0,mpi_comm_world,ierror)

  call mpi_barrier(mpi_comm_world,ierror)

  write (och5,'(I5.5)') proc

  proctablename = trim(tablename)//och5

  call mpi_barrier(mpi_comm_world,ierror)

  open(45,file=trim(proctablename))

  phistep = maxPhi/real(Nphi,kind=8)
  gphistep = maxGphi/real(Ngphi,kind=8)
  dphistep = maxDphi/real(Ndphi,kind=8)
  pistep = maxPi/real(Npi,kind=8)

  Nd = Ndphi + 1
  Ng = Ngphi + 1
  N  = Nphi  ! avoid maxPhi

  Npoints = Nd*Ng*N

  percent = 5.0
  do ijk = proc+1, Npoints, numproc

     k = (ijk - 1)/(Ng*N) + 1
     j = (ijk - 1 - (k-1)*Ng*N)/N + 1
     i = ijk - N*(j-1 + Ng*(k-1))

     dphi = real(k-1,8)*dphistep
     gphi = real(j-1,8)*gphistep
     phi  = real(i-1,8)*phistep
     pi   = zero

     A = - sixth*pistep*cos(as**3*dphi*pi) &
          *exp(-at*as**3*H(phi,gphi,pi))
     C = - sixth*pistep*cos(as**3*dphi*pi) &
          *exp(-at*as**3*H(phi,gphi,pi))*H(phi,gphi,pi)

     pido: do l = 1, Npi + 1

        A = A + &
             athird*pistep*cos(as**3*dphi*pi) &
             *exp(-at*as**3*H(phi,gphi,pi)) + &
             tthird*pistep*cos(as**3*dphi*(pi+half*pistep)) &
             *exp(-at*as**3*H(phi,gphi,pi+half*pistep))

        C = C + &
             athird*pistep*cos(as**3*dphi*pi) &
             *exp(-at*as**3*H(phi,gphi,pi))*H(phi,gphi,pi) + &
             tthird*pistep*cos(as**3*dphi*(pi+half*pistep)) &
             *exp(-at*as**3*H(phi,gphi,pi+half*pistep)) &
             *H(phi,gphi,pi+half*pistep)

        pi = pi + pistep

        if (pistep*exp(-at*as**3*H(phi,gphi,pi+half*pistep)) < eps ) then
           exit
        end if

     end do pido

     write (45,'(3(f10.5),2(es20.11))') phi,gphi,dphi,A,C

     Amax = max( A, Amax)
     Cmax = max(C,Cmax)
     Amin = min( A, Amin)
     Cmin = min(C,Cmin)

     if (proc == 0 .and. 100.0*real(ijk)/real(Npoints) > percent) then
        write(6,'(f10.4,a)') percent,' % done'
        percent = percent + 5.0
     end if

  end do

  close (45)

  call mpi_barrier(mpi_comm_world,ierror)

  call mpi_reduce( Amin, Amin0,1,MPI_DOUBLE,MPI_MIN,0,mpi_comm_world,ierror)
  call mpi_reduce(Cmin,Cmin0,1,MPI_DOUBLE,MPI_MIN,0,mpi_comm_world,ierror)
  call mpi_reduce( Amax, Amax0,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm_world,ierror)
  call mpi_reduce(Cmax,Cmax0,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm_world,ierror)

  call mpi_barrier(mpi_comm_world,ierror)

  if (proc == 0) timef = mpi_wtime()

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
  end if

  if (proc == 0) print *, timef - timei," seconds"

  call mpi_finalize(ierror)

contains
  function H(phi, gphi, pi)
    real(kind=8):: H, phi, gphi, pi
    H = half*( pi**2 +  (gphi/as)**2 ) + m**4/sqrt(one - (phi/Mass)**2) - m**4
  end function H
end program freecutoff
