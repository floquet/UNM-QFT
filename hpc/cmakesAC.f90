program makesAC
  implicit none
  ! gphi = sqrt( (f(iu,j,k,l) - f(i,j,k,l))**2 + ... )     (without /as**2)
  ! order is phi, gphi, dphi, pi
  ! PI STANDS FOR DOT PSI
  include 'mpif.h'

  character(len=50):: tablename
  character(len=55):: proctablename

  character(len=5)::och5

  integer(8):: i, j, k, l, ijk, Npi, Nphi, Ndphi, Ngphi
  integer(8):: Nd, Ng, N, Npoints, thisNpi
  real(8):: at, as, pi, phi, dphi, gphi, Mass, m
  real(8):: pistep, phistep, dphistep, gphistep
  real(8):: maxPhi, maxDphi, maxGphi, maxPi, thismaxPi
  real(8),parameter:: zero=0.0d0, one=1.0d0
  real(8),parameter:: half=0.5d0, two = 2.0d0
  real(8),parameter:: athird=one/3.0d0, tthird=2.0d0/3.0d0
  real(8),parameter:: twopi = 8.0d0*atan(1.0d0), sixth=one/6.0d0
  real(8):: A, C
  real(8),parameter:: eps=10.d0*epsilon(1.0d0)

  integer(kind=4):: ierror,proc,numproc
  real(8):: Amax,Cmax,Amin,Cmin
  real(8):: Amax0,Cmax0,Amin0,Cmin0
  real(8):: timei, timef, percent

  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world,proc,ierror)
  call mpi_comm_size(mpi_comm_world,numproc,ierror)

  timei = zero
  timef = zero
  if (proc == 0) timei = mpi_wtime()

  Amax = real(-1.e6,8)
  Cmax = real(-1.e6,8)
  Amin = real( 1.e6,8)
  Cmin = real( 1.e6,8)

  if (proc == 0) then
     open(5,file='incm1a10')
     read(5,*); read(5,*) maxPhi, maxGphi, maxDphi, maxPi
     read(5,*); read(5,*) Nphi, Ndphi, Ngphi, Npi
     read(5,*); read(5,*) as, at, Mass, m
     read(5,*); read(5,*) tablename
     close(5)
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

  phistep = maxPhi/real(Nphi,8)
  gphistep = maxGphi/real(Ngphi,8)
  dphistep = maxDphi/real(Ndphi,8)
  pistep = maxPi/real(Npi,8)

  Nd = Ndphi + 1
  Ng = Ngphi + 1
  N  = Nphi + 1

  Npoints = Nd*Ng*N

  percent = 5.0
  do ijk = proc+1, Npoints, numproc

     k = (ijk - 1)/(Ng*N) + 1
     j = (ijk - 1 - (k-1)*Ng*N)/N + 1
     i = ijk - N*(j-1 + Ng*(k-1))

     dphi = real(k-1,8)*dphistep
     gphi = real(j-1,8)*gphistep
     phi  = real(i-1,8)*phistep
     pistep = maxPi/real(Npi,8)

     thismaxPi = min(maxPi, sqrt(Mass**4 + (gphi/as)**2 + (m*phi)**2))
     thisNpi = int(thismaxPi/pistep,8)
     
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

     write (45,'(3(f16.5),2(es20.11))') phi,gphi,dphi,A,C

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
    real(8):: H, phi, gphi, pi
    H = ( two*pi**2 - (gphi/as)**2 - (m*phi)**2 - Mass**4 ) &
         /( ( one - ( pi**2 - (gphi/as)**2 - (m*phi)**2 )/Mass**4 ) &
         * sqrt( one - ( pi**2 - (gphi/as)**2 - (m*phi)**2 )/Mass**4 ) ) &
         + Mass**4
  end function H
  function D(phi, gphi, pi) 
    real(8):: D, phi, gphi, pi
    D = ( one + ( two*pi**2 +  (gphi/as)**2 + (m*phi)**2 )/Mass**4 ) &
         /( ( one - ( pi**2 -  (gphi/as)**2 - (m*phi)**2 )/Mass**4 )**2 &
         *sqrt( ( one - ( pi**2 -  (gphi/as)**2 - (m*phi)**2 )/Mass**4 ) ) )
  end function D
  function F(phi, gphi, dphi, pi) 
    real(8):: F, phi, gphi, dphi, pi
    F = dphi*pi  &
         /( ( one - ( pi**2 -  (gphi/as)**2 - (m*phi)**2 )/Mass**4 ) &
         *sqrt( ( one - ( pi**2 -  (gphi/as)**2 - (m*phi)**2 )/Mass**4 ) ) )
  end function F
end program makesAC
