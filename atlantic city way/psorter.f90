program psorter

  implicit none

  character(len=50):: tablename
  character(len=55):: proctablename

  character(len=1)::och1
  character(len=2)::och2
  character(len=3)::och3
  character(len=4)::och4
  character(len=5)::och5

  integer(kind=8):: i, j, k, Npi, Nphi, Ndphi, Ngphi
  real(kind=8):: at, a, Mass, m
  real(kind=8):: maxPhi, maxDphi, maxGphi, maxPi
  real(kind=8),parameter:: zero=0.0d0, one=1.0d0, half=0.5d0
  real(kind=8),parameter:: athird=one/3.0d0, tthird=2.0d0/3.0d0
  real(kind=8),parameter:: twopi = 8.0d0*atan(1.0d0), sixth=one/6.0d0
  real(kind=8),parameter:: eps=10.d0*epsilon(1.0d0)
  real(kind=8), allocatable:: phia(:), gphia(:), dphia(:)
  real(kind=8), allocatable:: Integral(:,:,:), HIntegral(:,:,:)
  integer(kind=4)::numproc, proc

  read(5,*); read(5,*) maxPhi, maxDphi, maxGphi, maxPi
  read(5,*); read(5,*) Nphi, Ndphi, Ngphi, Npi
  read(5,*); read(5,*) at, a, Mass, m
  read(5,*); read(5,*) tablename, numproc

  allocate(phia(Nphi+1), gphia(Ngphi+1), dphia(Ndphi+1))
  allocate(Integral(Nphi+1,Ngphi+1,Ndphi+1))
  allocate(HIntegral(Nphi+1,Ngphi+1,Ndphi+1))

  phia = zero; gphia = zero; dphia = zero
  Integral = zero; HIntegral = zero

  do proc = 0,numproc-1

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

     open(38,file=trim(proctablename),form='formatted',status='unknown')

     do k = 1, Ndphi + 1
        do j = 1, Ngphi + 1
           do i = proc+1, Nphi + 1, numproc
              read (38,*) phia(i),gphia(j),dphia(k), &
                          Integral(i,j,k),Hintegral(i,j,k)
           end do
        end do
     end do

     close(38)

  end do

  open(8,file=trim(tablename))
  dphiw: do k = 1, Ndphi+1
    gphiw: do j= 1, Ngphi+1
      phiw: do i = 1, Nphi+1
        write (8,'(3(f12.5),2(es20.11))') phia(i), gphia(j), dphia(k), &
        Integral(i,j,k), HIntegral(i,j,k)
      end do phiw
    end do gphiw
  end do dphiw
  close(8)

  write(6,*)'minimum of Integral =', minval(Integral)
  write(6,*) 'maximum of Integral =', maxval(Integral)
  write(6,*) 'minimum of HIntegral =', minval(HIntegral)
  write(6,*) 'maximum of HIntegral =', maxval(HIntegral)

end program psorter
