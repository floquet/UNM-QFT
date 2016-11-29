program table
  implicit none
  character(len=50):: tablename
  integer(kind=8):: i, j, k, l, Npi, Nphi, Ndphi, Ngphi
  real(kind=8):: at, a, pi, phi, dphi, gphi, Mass, m
  real(kind=8):: pistep, phistep, dphistep, gphistep
  real(kind=8):: maxPhi, maxDphi, maxGphi, maxPi
  real(kind=8),parameter:: zero=0.0d0, one=1.0d0, half=0.5d0
  real(kind=8),parameter:: athird=one/3.0d0, tthird=2.0d0/3.0d0
  real(kind=8),parameter:: twopi = 8.0d0*atan(1.0d0), sixth=one/6.0d0
  real(kind=8), allocatable:: phia(:), gphia(:), dphia(:)
  real(kind=8), allocatable:: Integral(:,:,:), HIntegral(:,:,:)
  real(kind=8),parameter:: eps=10.d0*epsilon(1.0d0)

  read(5,*); read(5,*) maxPhi, maxDphi, maxGphi, maxPi
  read(5,*); read(5,*) Nphi, Ndphi, Ngphi, Npi
  read(5,*); read(5,*) at, a, Mass, m
  read(5,*); read(5,*) tablename

  allocate(phia(Nphi+1), gphia(Ngphi+1), dphia(Ndphi+1))
  allocate(Integral(Nphi+1,Ngphi+1,Ndphi+1))
  allocate(HIntegral(Nphi+1,Ngphi+1,Ndphi+1))

  Integral=zero; HIntegral=zero
  phistep = maxPhi/real(Nphi,kind=8)
  gphistep = maxGphi/real(Ngphi,kind=8)
  dphistep = maxDphi/real(Ndphi,kind=8)
  pistep = maxPi/real(Npi,kind=8)

  dphi = zero
  dphido: do k = 1, Ndphi + 1
    dphia(k) = dphi
    gphi = zero
    gphido: do j= 1, Ngphi + 1
      gphia(j) = gphi
      phi = zero
      phido: do i = 1, Nphi + 1
        phia(i) = phi
        pi = zero
        Integral(i,j,k) = - sixth*pistep*cos(a**3*dphi*pi) &
        *exp(-at*a**3*H(phi,gphi,pi))
        HIntegral(i,j,k) = - sixth*pistep*cos(a**3*dphi*pi) &
        *exp(-at*a**3*H(phi,gphi,pi))*H(phi,gphi,pi)
        pido: do l = 1, Npi + 1
          Integral(i,j,k) = Integral(i,j,k) + &
          athird*pistep*cos(a**3*dphi*pi) &
          *exp(-at*a**3*H(phi,gphi,pi)) + &
          tthird*pistep*cos(a**3*dphi*(pi+half*pistep)) &
          *exp(-at*a**3*H(phi,gphi,pi+half*pistep))
          HIntegral(i,j,k) = HIntegral(i,j,k) + &
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
        phi = phi + phistep
      end do phido
      gphi = gphi + gphistep
    end do gphido
    write(6,'(f10.4,a)') real(100*k)/real(Ndphi),' % done'
    dphi = dphi + dphistep
  end do dphido

  open(8,file=trim(tablename))
  dphiw: do k = 1, Ndphi+1
    gphiw: do j= 1, Ngphi+1
      phiw: do i = 1, Nphi+1
        write(8,*) phia(i), gphia(j), dphia(k), Integral(i,j,k), HIntegral(i,j,k)
      end do phiw
    end do gphiw
  end do dphiw
  close(8)

  write(6,*) 'maxPhi, maxDphi, maxGphi, maxPi'
  write(6,*) maxPhi, maxDphi, maxGphi, maxPi
  write(6,*) 'Nphi, Ndphi, Ngphi, Npi'
  write(6,*) Nphi, Ndphi, Ngphi, Npi
  write(6,*) 'at, a, Mass, m'
  write(6,*) at, a, Mass, m
  write(6,*) 'tablename'
  write(6,*) tablename
  write(6,*) 'minimum of Integral =', minval(Integral)
  write(6,*) 'maximum of Integral =', maxval(Integral)
  write(6,*) 'minimum of HIntegral =', minval(HIntegral)
  write(6,*) 'maximum of HIntegral =', maxval(HIntegral)

contains
  function H(phi, gphi, pi)
    real(kind=8):: H, phi, gphi, pi
    H = sqrt((Mass**4 + pi**2)*(Mass**4 + gphi**2 + (m*phi)**2)) - Mass**4
  end function H
end program
