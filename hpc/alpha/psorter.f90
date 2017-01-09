program psorter

  implicit none

  character(len=50):: tablename
  character(len=55):: proctablename
  character(len=5)::och5

  integer(kind=4):: i, j, k, Nphi, Ndphi, Ngphi
  integer(8):: Npi
  integer(kind=4):: ijk, Nijk
  real(kind=8):: at, as, Mass, m
  real(kind=8):: maxPi, maxPhi, maxGphi, maxDphi
  real(8):: phistep, dphistep, gphistep
  real(kind=8),parameter:: zero=0.0d0, one=1.0d0, half=0.5d0
  real(kind=8),parameter:: athird=one/3.0d0, tthird=2.0d0/3.0d0
  real(kind=8),parameter:: twopi = 8.0d0*atan(1.0d0), sixth=one/6.0d0
  real(kind=8),parameter:: two = 2.0_8, three = 3.0_8
  real(kind=8),parameter:: eps=10.d0*epsilon(1.0d0)
  real(kind=8):: phi, gphi, dphi
  real(kind=8), allocatable:: A(:,:,:), C(:,:,:)
  integer(kind=4)::numproc, proc
  integer(kind=4)::N, Ng, Nd
  

  read(5,*); read(5,*) maxPhi, maxGphi, maxDphi, maxPi
  read(5,*); read(5,*) Nphi, Ngphi, Ndphi, Npi, numproc
  read(5,*); read(5,*) as, at, Mass, m
  read(5,*); read(5,*) tablename

  N  = Nphi 
  Ng = Ngphi + 1
  Nd = Ndphi + 1

  Nijk = N*Ng*Nd
  allocate(A(N,Ng,Nd))
  allocate(C(N,Ng,Nd))

  A = zero; C = zero

  do proc = 0, numproc-1

     write(och5,'(I5.5)') proc

     proctablename = trim(tablename)//och5

     open(38,file=trim(proctablename),form='formatted',status='unknown')

     do ijk = proc+1, Nijk, numproc
        k = (ijk - 1)/(N*Ng) + 1
        j = (ijk - 1 - (k-1)*N*Ng)/N + 1
        i = ijk - N*(j - 1 + Ng*(k - 1))
        read (38,*) phi,gphi,dphi,A(i,j,k),C(i,j,k)
     end do
     close(38)

  end do

  phistep = maxPhi/real(Nphi,8)
  gphistep = maxGphi/real(Ngphi,8)
  dphistep = maxDphi/real(Ndphi,8)

  open(8,file=trim(tablename))
  dphiw: do k = 1, Nd
     gphiw: do j= 1, Ng
        phiw: do i = 1, N
           phi = real(i-1,8)*phistep
           gphi = real(j-1,8)*gphistep
           dphi = real(k-1,8)*dphistep
           write (8,'(3(f16.5,1x),2(es20.11,1x))') phi, gphi, dphi, &
                A(i,j,k), C(i,j,k)
        end do phiw
     end do gphiw
  end do dphiw
  close(8)

  write(6,*)'minimum of A =', minval(A)
  write(6,*) 'maximum of A =', maxval(A)
  write(6,*) 'minimum of C =', minval(C)
  write(6,*) 'maximum of C =', maxval(C)

end program psorter
