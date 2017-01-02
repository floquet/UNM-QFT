program psorter

  implicit none

  character(len=50):: tablename
  character(len=55):: proctablename
  character(len=5)::och5

  integer(kind=8):: i, j, k, Npi, Nphi, Ndphi, Ngphi
  integer(kind=8):: ijk, Nijk
  real(kind=8):: at, as, Mass, m
  real(kind=8):: maxPhi, maxDphi, maxGphi, maxPi
  real(kind=8),parameter:: zero=0.0d0, one=1.0d0, half=0.5d0
  real(kind=8),parameter:: athird=one/3.0d0, tthird=2.0d0/3.0d0
  real(kind=8),parameter:: twopi = 8.0d0*atan(1.0d0), sixth=one/6.0d0
  real(kind=8),parameter:: eps=10.d0*epsilon(1.0d0)
  real(kind=8), allocatable:: phia(:), gphia(:), dphia(:)
  real(kind=8), allocatable:: A(:,:,:), C(:,:,:)
  integer(kind=4)::numproc, proc

  read(5,*); read(5,*) maxPhi, maxGphi, maxDphi, maxPi
  read(5,*); read(5,*) Nphi, Ngphi, Ndphi, Npi, numproc
  read(5,*); read(5,*) as, at, Mass, m
  read(5,*); read(5,*) tablename

  Nijk = (Nphi + 1)*(Ngphi + 1)*(Ndphi + 1)
  allocate(phia(Nphi+1), gphia(Ngphi+1), dphia(Ndphi+1))
  allocate(A(Nphi+1,Ngphi+1,Ndphi+1))
  allocate(C(Nphi+1,Ngphi+1,Ndphi+1))

  phia = zero; gphia = zero; dphia = zero
  A = zero; C = zero

  do proc = 0, numproc-1

    write(och5,'(I5.5)') proc

    proctablename = trim(tablename)//och5

    open(38,file=trim(proctablename),form='formatted',status='unknown')


    do ijk = proc+1, Nijk, numproc
      k = (ijk - 1)/((Nphi + 1)*(Ngphi + 1)) + 1
      j = (ijk - 1 - (k-1)*(Nphi + 1)*(Ngphi + 1))/(Nphi + 1) + 1
      i = ijk - (Nphi + 1)*(j - 1 + (Ngphi + 1)*(k - 1))
      read (38,*) phia(i),gphia(j),dphia(k),A(i,j,k),C(i,j,k)
    end do
    close(38)

  end do

  open(8,file=trim(tablename))
  dphiw: do k = 1, Ndphi+1
    gphiw: do j= 1, Ngphi+1
      phiw: do i = 1, Nphi+1
        write (8,'(3(f16.5,1x),2(es20.11,1x))') phia(i), gphia(j), dphia(k), &
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
