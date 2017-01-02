program ammonte
  implicit none
  ! gphi = sqrt( (f(iu,j,k,l) - f(i,j,k,l))**2 + ... )     (without /as**2)
  character(len=2):: och2
  character(len=50)::tablename, temp, farray, root, infile
  integer(8):: i, j, k, l, iu, ju, ku, lu, id, jd, kd, ld, Nsweeps, sweep
  integer(8):: index, Nphi, Ndphi, Ngphi, Ns, Nt
  integer(8), allocatable::ups(:), dns(:), upt(:), dnt(:)
  real(8),allocatable:: f(:,:,:,:)
  real(8):: as, at, rdn, Mass, m, maxPhi, maxDphi, maxGphi
  real(8):: df, count, E, meanE, oldf
  real(8):: naccept, nreject, Esweep, minA
  real(8):: newP, oldP, sigma, highphi, highgphi, highdphi
  real(8):: phistep, gphistep, dphistep, outoftable
  real(8),allocatable:: A(:,:,:), C(:,:,:), E_0(:)
  real(8),allocatable:: phi(:), gphi(:), dphi(:)
  real(8),parameter:: zero=0.0d0, half = 0.5d0, one=1.0d0

  read(5,*); read(5,*) maxPhi, maxGphi, maxDphi
  read(5,*); read(5,*) Nphi, Ngphi, Ndphi
  read(5,*); read(5,*) as, at, Mass, m, df
  read(5,*); read(5,*) tablename, temp, root, farray, index
  read(5,*); read(5,*) Nsweeps, Ns, Nt

  allocate (f(Ns,Ns,Ns,Nt))
  allocate(A(Nphi+1, Ngphi+1, Ndphi+1), C(Nphi+1, Ngphi+1, Ndphi+1))
  allocate(phi(Nphi+1), gphi(Ngphi+1), dphi(Ndphi+1))
  allocate(ups(Ns), dns(Ns), upt(Nt), dnt(Nt))
  allocate(E_0(Nsweeps))

  highphi = zero; highgphi= zero; highdphi = zero
  minA = 1.0D200; outoftable = zero

  upt(Nt) = 1
  dnt(1) = Nt
  do i = 1, Nt-1
     upt(i) = i+1 ! make tables for periodic boundary conditions
  end do
  do i = 2, Nt
     dnt(i) = i - 1
  end do
  ups(Ns) = 1
  dns(1) = Ns
  do i = 1, Ns-1
     ups(i) = i+1 ! make tables for periodic boundary conditions
  end do
  do i = 2, Ns
     dns(i) = i - 1
  end do

  phistep = maxPhi/real(Nphi,8)
  gphistep = maxGphi/real(Ngphi,8)
  dphistep = maxDphi/real(Ndphi,8)

  open(4, file=trim(tablename))
  do k = 1, Ndphi+1
     do j = 1, Ngphi+1
        do i = 1, Nphi+1
           read(4,*) phi(i), gphi(j), dphi(k), A(i,j,k), C(i,j,k)
           C(i,j,k) = C(i,j,k)/A(i,j,k)
        end do
     end do
  end do
  close(4)

  call init_random_seed() ! set new seed
  if (temp == 'hot') then
     open(4,file=farray)
     read(4,*) f
     write(6,*) 'The temp is ', temp
  else if (temp == 'cold') then
     f = zero
     write(6,*) 'The temp is ', temp
  else
     write(6,*) 'I need to know the temperature.'
     stop
  end if

  E_0 = zero
  E = zero; count = zero; naccept = zero; nreject = zero

  sweepdo: do sweep = 1, Nsweeps
     do i = 1, Ns
        iu=ups(i); id=dns(i)
        do j = 1, Ns
           ju=ups(j); jd=dns(j)
           do k = 1, Ns
              ku=ups(k); kd=dns(k)
              do l = 1, Nt ! time
                 lu=upt(l); ld=dnt(l)
                 oldf = f(i,j,k,l)
                 oldP = aA(iu,ju,ku,lu,i,j,k,l)*aA(iu,ju,ku,l,i,j,k,ld) &
                      *aA(i,ju,ku,lu,id,j,k,l)*aA(iu,j,ku,lu,i,jd,k,l) &
                      *aA(iu,ju,k,lu,i,j,kd,l)
                 call random_number(rdn)
                 f(i,j,k,l) = oldf + df*( rdn - half )
                 newP = aA(iu,ju,ku,lu,i,j,k,l)*aA(iu,ju,ku,l,i,j,k,ld) &
                      *aA(i,ju,ku,lu,id,j,k,l)*aA(iu,j,ku,lu,i,jd,k,l) &
                      *aA(iu,ju,k,lu,i,j,kd,l)
                 if ( newP >= oldP ) then ! accept
                    naccept = naccept + one
                 else
                    call random_number(rdn)
                    if ( newP/oldP >= rdn ) then ! accept
                       naccept = naccept + one
                    else ! reject
                       nreject = nreject + one
                       f(i,j,k,l) = oldf
                    end if
                 end if
              end do
           end do
        end do
     end do

     Esweep = zero
     do i = 1, Ns
        iu=ups(i)
        do j = 1, Ns
           ju=ups(j)
           do k = 1, Ns
              ku=ups(k)
              do l = 1, Nt ! time
                 lu=upt(l)
                 Esweep = Esweep + aC(iu,ju,ku,lu,i,j,k,l)
              end do
           end do
        end do
     end do
     E_0(sweep) = Esweep/real(Nt*Ns**3,8)

  end do sweepdo

  meanE = sum(E_0)/real(Nsweeps,8)
  sigma = zero
  do i = 1, Nsweeps
     sigma = sigma + (E_0(i) - meanE)**2
  end do
  sigma = sigma/real(Nsweeps,8)
  sigma = sqrt(sigma/real(Nsweeps-1,8))

  write(och2,"(i2.2)") index
  farray = trim('array'//trim(root)//och2)
  infile = trim('in'//trim(root)//och2)
  open(8,file=trim(infile))
  write(8,*) 'maxPhi,    maxGphi,    maxDphi'
  write(8,'(3(f12.3))') maxPhi, maxGphi, maxDphi
  write(8,*) 'Nphi,   Ngphi,    Ndphi'
  write(8,'(3(I10))') Nphi, Ngphi, Ndphi
  write(8,*) 'as,   at,   Mass,   m,   df'
  write(8,'(5(f12.6))') as, at, Mass, m, df
  write(8,*) 'tablename,     temp,    root,    farray,    index'
  write(8,'(5(a),I10)') trim(tablename), ' hot ', trim(root), '  ', trim(farray), index+1
  write(8,*) 'Nsweeps,   Ns,   Nt'
  write(8,'(3(I10))') Nsweeps, Ns, Nt
  write(8,*) 'Results from run', index
  write(6,*) 'Results from run', index
  write(8,*)  'E_0 =', real(meanE,4), 'sigma =', sigma
  write(6,*)  'E_0 =', real(meanE,4), 'sigma =', sigma
  write(8,*) 'accept ratio = ', real(naccept/(naccept+nreject), 4)
  write(6,*) 'accept ratio = ', real(naccept/(naccept+nreject), 4)
  write(6,*) 'minimum of A =', minval(A)
  write(6,*) 'maximum of A =', maxval(A)
  write(6,*) 'minimum of C =', minval(C)
  write(6,*) 'maximum of C =', maxval(C)
  write(8,*)'minimum of A =', minval(A)
  write(8,*) 'maximum of A =', maxval(A)
  write(8,*) 'minimum of C =', minval(C)
  write(8,*) 'maximum of C =', maxval(C)
  write(6,*) 'highphi =', highphi
  write(6,*) 'highgphi =', highgphi
  write(6,*) 'highdphi =', highdphi
  write(8,*) 'highphi =', highphi
  write(8,*) 'highgphi =', highgphi
  write(8,*) 'highdphi =', highdphi
  i = int(highphi/phistep) + 1
  j = int(highgphi/gphistep) + 1
  k = int(highdphi/dphistep) + 1
  write(8,*) 'A(highphi,highgphi,highdphi) =', A(i,j,k)
  write(6,*) 'A(highphi,highgphi,highdphi) =', A(i,j,k)
  write(8,*) 'minimum A =', minA
  write(6,*) 'minimum A =', minA
  write(8,*) 'out of table =', outoftable
  write(6,*) 'out of table =', outoftable
  open(7,file=trim(farray))
  write(7,*) f
  close(7)

contains

  function aA(iu,ju,ku,lu,i,j,k,l)
    integer(8):: iphi, igphi, idphi, iu,ju,ku,lu,i,j,k,l
    real(8):: aA, phil, gphil, dphil, u, v, w
    phil = abs(f(i,j,k,l))
    gphil = sqrt( ( (f(iu,j,k,l) - f(i,j,k,l))**2 &
         +  (f(i,ju,k,l) - f(i,j,k,l))**2 &
         +  (f(i,j,ku,l) - f(i,j,k,l))**2 ) )
    dphil = abs(f(i,j,k,lu) - f(i,j,k,l))
    highphi = max(phil, highphi); highgphi = max(gphil, highgphi)
    highdphi = max(dphil, highdphi)
    if ( phil < maxPhi .and. gphil < maxGphi .and. dphil < maxDphi ) then
       iphi = int(phil/phistep) + 1
       igphi = int(gphil/gphistep) + 1
       idphi = int(dphil/dphistep) + 1
       minA = min(minA, A(iphi,igphi,idphi))
       u = (phil - phi(iphi))/(phi(iphi+1) - phi(iphi))
       v = (gphil - gphi(igphi))/(gphi(igphi+1) - gphi(igphi))
       w = (dphil - dphi(idphi))/(dphi(idphi+1) - dphi(idphi))
       aA = A(iphi,igphi,idphi)*(one-u)*(one-v)*(one-w) &
            + A(iphi+1,igphi,idphi)*(u)*(one-v)*(one-w) &
            + A(iphi,igphi+1,idphi)*(one-u)*(v)*(one-w) &
            + A(iphi,igphi,idphi+1)*(one-u)*(one-v)*(w) &
            + A(iphi+1,igphi+1,idphi)*(u)*(v)*(one-w) &
            + A(iphi+1,igphi,idphi+1)*(u)*(one-v)*(w) &
            + A(iphi,igphi+1,idphi+1)*(one-u)*(v)*(w) &
            + A(iphi+1,igphi+1,idphi+1)*(u)*(v)*(w)
    else
       aA = zero
       outoftable = outoftable + one
    end if
  end function aA

  function aC(iu,ju,ku,lu,i,j,k,l)
    integer(8):: iphi, igphi, idphi, iu,ju,ku,lu,i,j,k,l
    real(8):: aC, phil, gphil, dphil, u, v, w
    phil = abs(f(i,j,k,l))
    gphil = sqrt( ((f(iu,j,k,l) - f(i,j,k,l))**2 &
         +  (f(i,ju,k,l) - f(i,j,k,l))**2 &
         +  (f(i,j,ku,l) - f(i,j,k,l))**2) )
    dphil = abs(f(i,j,k,lu) - f(i,j,k,l))
    if ( phil < maxPhi .and. gphil < maxGphi .and. dphil < maxDphi ) then
       iphi = int(phil/phistep) + 1
       igphi = int(gphil/gphistep) + 1
       idphi = int(dphil/dphistep) + 1
       u = (phil - phi(iphi))/(phi(iphi+1) - phi(iphi))
       v = (gphil - gphi(igphi))/(gphi(igphi+1) - gphi(igphi))
       w = (dphil - dphi(idphi))/(dphi(idphi+1) - dphi(idphi))
       aC = C(iphi,igphi,idphi)*(one-u)*(one-v)*(one-w) &
            + C(iphi+1,igphi,idphi)*(u)*(one-v)*(one-w) &
            + C(iphi,igphi+1,idphi)*(one-u)*(v)*(one-w) &
            + C(iphi,igphi,idphi+1)*(one-u)*(one-v)*(w) &
            + C(iphi+1,igphi+1,idphi)*(u)*(v)*(one-w) &
            + C(iphi+1,igphi,idphi+1)*(u)*(one-v)*(w) &
            + C(iphi,igphi+1,idphi+1)*(one-u)*(v)*(w) &
            + C(iphi+1,igphi+1,idphi+1)*(u)*(v)*(w) 
    else
       aC = zero
    end if
  end function aC

  subroutine init_random_seed()
    implicit none
    integer ir, nr, clock
    integer, dimension(:), allocatable :: seed
    call random_seed(size = nr) ! find size of seed
    allocate(seed(nr))
    call system_clock(count=clock) ! get time of processor clock
    seed = clock + 37 * (/ (ir-1, ir=1, nr) /) ! make seed
    call random_seed(put=seed) ! set seed
    deallocate(seed)
  end subroutine init_random_seed
end program ammonte
