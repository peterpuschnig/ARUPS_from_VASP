!!$******************************* ARUPS ************************************
!!$*
!!$*   input the WAVECAR file in binary format and the IBZKPT file from VASP
!!$*   output calculated ARUPS intensities for 3 different modes specified in
!!$*   the input file. 
!!$*
!!$*   Compile with gfortran or ifort. Flag "-assume byterecl" is required
!!$*   for ifort.
!!$*
!!$*   with the help of the program Wavetrans version 2.0 - July 5, 2012 - by R. M. Feenstra and M. Widom
!!$*
!!$*
!!$*
!!$*   format of the input file which has to be named ARUPS.in
!!$*
!!$*   "string"                      ! path and name of output file 
!!$*   integerN                      ! number of k-points
!!$*   integerN integerN             ! first and last bandindex for summation over initial states
!!$*   realN                         ! Fermi-level (energy of HOMO) in [eV]
!!$*   realN                         ! Temperature in [eV] for Fermi-Dirac smearing, negative value means no FD smearing
!!$*   realN                         ! work function in [eV]
!!$*   realN realN                   ! energy range of incident photon [eV], only first number used for Energy vs k and kx-ky plot
!!$*   realN                         ! energy broadening [eV]
!!$*   realN realN                   ! broadening in k-space [1/A]
!!$*   integerN                      ! 0 = Energy versus k, 1 = kx-ky plot , 2 = CIS scan
!!$*   realN realN realN             ! binding energy range and step size in [eV], kx-ky plot and CIS scan use first binding energy only (set binding energy below fermi energy to negative value)
!!$*   realN realN realN             ! range of k-vectors and step size in [1/A], Energy vs k and CIS scan use first vector only
!!$*   realN realN realN             ! unit vector for normal emission
!!$*   realN realN realN             ! unit vector for grazing emission
!!$*   realN realN                   ! damping parameter and z-position of surface, negative damping parameter means no damping
!!$*   integerN                      ! 0: LNONCOLLINEAR=.FALSE. , 1: LNONCOLLINEAR=.TRUE.
!!$*
!!$*
!!$*   note that the energy eigenvalues are complex, as provided in the
!!$*   WAVECAR file, but the imaginary part is zero (at least for all cases
!!$*   investigated thus far)
!!$*     

implicit none
real*8                  ::      c, wtime, deltak, ekinmax
integer                 ::      iost, iplane
real*8                  ::      pi, two_pi
real*8                  ::      starttime, endtime
real*8                  ::      epsilon, ekintok
complex*16              ::      imag

!!$* INPUT-file
character (1000)        ::      outfile
integer                 ::      khi
integer                 ::      bandlow, bandhi
real*8                  ::      efermi
real*8                  ::      temperature
real*8                  ::      phiwork
real*8                  ::      ephotonmin, ephotonmax
real*8                  ::      ebroad
real*8                  ::      kxbroad, kybroad
integer                 ::      arupstype
real*8                  ::      eklo,ekhi,ekstep
real*8                  ::      kxlo,kxhi,kxstep
real*8                  ::      vecz(3)
real*8                  ::      vecx(3)
real*8                  ::      gamma,z0 
integer                 ::      LNONCOLLINEAR

!!$* arrays
real*8, allocatable     ::      photoint_sum(:,:),photoint2_sum(:,:),photointspin_sum(:,:),photointspin2_sum(:,:) 
real*8, allocatable     ::      photoint_pol(:,:),photoint2_pol(:,:),photointspin_pol(:,:),photointspin2_pol(:,:) 
real*8, allocatable     ::      photoint_nup(:,:),photoint2_nup(:,:),photointspin_nup(:,:),photointspin2_nup(:,:) 
real*8, allocatable     ::      photoint_ndo(:,:),photoint2_ndo(:,:),photointspin_ndo(:,:),photointspin2_ndo(:,:) 

!!$*Wavetrans
real*8                  ::      xnwk, xnband, ecut
integer                 ::      nwk, nband
real*8                  ::      a1(3),a2(3),a3(3),b1(3),b2(3),b3(3), vtmp(3), wk(3), sumkg(3)
real*8                  ::      b1mag,b2mag,b3mag
real*8                  ::      Vcell, ccell
real*8                  ::      xnrecl, xnspin, xnprec
integer                 ::      nrecl, nspin, nprec, irec
complex*16, allocatable ::      cener(:)
real*8, allocatable     ::      occ(:)
real*8, allocatable     ::      igall(:,:)
integer, allocatable    ::      hkl(:,:)
complex*8, allocatable  ::      coeff_u(:,:,:), coeff_d(:,:,:), coeff(:,:,:)
real*8                  ::      phi12, sinphi123, phi13, phi23, phi123, vmag
integer                 ::      nb1maxA, nb2maxA, nb3maxA, npmaxA
integer                 ::      nb1maxB, nb2maxB, nb3maxB, npmaxB
integer                 ::      nb1maxC, nb2maxC, nb3maxC, npmaxC
integer                 ::      nb1max, nb2max, nb3max, npmax
real*8                  ::      xnplane
integer                 ::      iband, nplane, nplane2
integer                 ::      ncnt
!integer                 ::      npmaxA, npmaxB, npmaxC

integer                 ::      i, j, isp, iwk, bandcount, iek, jp, isymmetry, jgz, ikx, iky, jg, ieph
integer                 ::      ig3, ig3p, ig2, ig2p, ig1, ig1p
integer                 ::      nkx, nek, nephk, nsymmetry
logical                 ::      pureplanewave
integer                 ::      iwindow
real*8                  ::      ekinmin
real*8                  ::      vecy(3), qplusG(3)
real*8                  ::      gtot, etot
real*8                  ::      eigval, ebind, fermifac
real*8                  ::      efinal
integer                 ::      iekstart, iekend
real*8                  ::      ekarr, kfinal, deltae, l1, l2
integer                 ::      hmil, kmil, lmil
real*8                  ::      qplusGlen, qplusGx, qplusGy, qplusGz
integer                 ::      ikxend, ikxstart, ikyend, ikystart
real*8                  ::      cgcabs_s, cgcabs_p, cgcabs_u, cgcabs_d
real*8                  ::      cgcabs2_s, cgcabs2_p, cgcabs2_u, cgcabs2_d
real*8                  ::      kxarr, kyarr, kpar, kperp
integer                 ::      ephotonstart, ephotonend, hmilz, kmilz, lmilz
complex*16              ::      fgz, cgccomplex_s, cgccomplex_p, cgccomplex_u, cgccomplex_d
real*8                  ::      gz, dnorm, ephoton


!!$*   constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 (value is
!!$*   adjusted in final decimal places to agree with VASP value)

data c/0.262465831d0/

!data c/0.26246582250210965422d0/ 
pi=4.*atan(1.)
two_pi    =  2.0d0*3.14159265d0

starttime = wtime ( ); ! get wall time

OPEN(UNIT=81,FILE='ARUPS.in', STATUS='UNKNOWN')

read(81,*) outfile               ! name of output file
read(81,*) khi                   !  dummy !! not used anymore, this line is only kept for backward compatibility
read(81,*) bandlow,bandhi        ! first and last bandindex for summation over initial states
read(81,*) efermi                ! Fermi-level (energy of HOMO) in [eV]
read(81,*) temperature           ! Temperature in [eV]
read(81,*) phiwork               ! work function in [eV]
read(81,*) ephotonmin,ephotonmax ! energy range of incident photon [eV]
read(81,*) ebroad                ! energy broadening [eV]
read(81,*) kxbroad,kybroad       ! broadening in k-space [1/A]
read(81,*) arupstype             ! 0 = Energy versus k, 1 = kx-ky plot at constant energy
read(81,*) eklo,ekhi,ekstep      ! kinetic energy range (1 - 2) and step size (3) in [eV] 
read(81,*) kxlo,kxhi,kxstep      ! range of final state k-vectors (1 - 2) and step (3) in [1/A]
read(81,*) vecz                  ! unit vector for normal emission
read(81,*) vecx                  ! unit vector for grazing emission
read(81,*) gamma,z0              ! damping parameter and z-position of surface in unit cell
                                 ! used for damping the plane wave final state
!read(81,*) LNONCOLLINEAR         ! determine type of WAVECAR, 0->normal, 1->spinor-form
LNONCOLLINEAR = 1
close(UNIT=81)


OPEN(UNIT=80,FILE='ARUPS.log',STATUS='UNKNOWN')
OPEN(UNIT=82,FILE='ARUPS_s.out',STATUS='UNKNOWN')
OPEN(UNIT=83,FILE='ARUPS_p.out',STATUS='UNKNOWN')
OPEN(UNIT=84,FILE='ARUPS_u.out',STATUS='UNKNOWN')
OPEN(UNIT=85,FILE='ARUPS_d.out',STATUS='UNKNOWN')

WRITE(80,*) 'CALCULATION OF ARUPS SPECTRA FROM WAVECAR FILE' 


!!$* allocate photointesity array
nkx = (kxhi - kxlo)/kxstep + 1
nek = (ekhi - eklo)/ekstep + 1 
nephk = (ephotonmax - ephotonmin)/ekstep + 1

if (arupstype .eq. 0) then
    allocate(photoint_sum(nkx,nek))
    allocate(photointspin_sum(nkx,nek))
    allocate(photoint_pol(nkx,nek))
    allocate(photointspin_pol(nkx,nek))
    allocate(photoint_nup(nkx,nek))
    allocate(photointspin_nup(nkx,nek))
    allocate(photoint_ndo(nkx,nek))
    allocate(photointspin_ndo(nkx,nek))
elseif (arupstype .eq. 1) then
    allocate(photoint_sum(nkx,nkx))
    allocate(photointspin_sum(nkx,nkx))
    allocate(photoint_pol(nkx,nkx))
    allocate(photointspin_pol(nkx,nkx))
    allocate(photoint_nup(nkx,nkx))
    allocate(photointspin_nup(nkx,nkx))
    allocate(photoint_ndo(nkx,nkx))
    allocate(photointspin_ndo(nkx,nkx))
else
    allocate(photoint_sum(nkx,nephk))
    allocate(photointspin_sum(nkx,nephk))
    allocate(photoint_pol(nkx,nephk))
    allocate(photointspin_pol(nkx,nephk))
    allocate(photoint_nup(nkx,nephk))
    allocate(photointspin_nup(nkx,nephk))
    allocate(photoint_ndo(nkx,nephk))
    allocate(photointspin_ndo(nkx,nephk))
end if
photoint_sum(:,:) = 0.0
photointspin_sum(:,:) = 0.0
photoint_pol(:,:) = 0.0
photointspin_pol(:,:) = 0.0
photoint_nup(:,:) = 0.0
photointspin_nup(:,:) = 0.0
photoint_ndo(:,:) = 0.0
photointspin_ndo(:,:) = 0.0


!!$* open WAVECAR and read header
nrecl = 24
open(unit=10, file='WAVECAR',access='direct',recl=nrecl, &
    iostat=iost,status='old')
if (iost.ne.0) then
    write(80,*) 'WAVECAR open error - iostat =',iost            
    stop
endif

read(unit=10,rec=1) xnrecl,xnspin,xnprec
close(unit=10)

nrecl=nint(xnrecl)
nspin=nint(xnspin)
nprec=nint(xnprec)

if(nprec.eq.45210) then
   write(80,*) '*** error - WAVECAR_double requires complex*16'
   stop
endif

if(nspin.ne.1) then
    write(80,*) 'program requires VASP calculation LNONCOLLINEAR=.TRUE.'
    stop
endif

write(80,*) 
write(80,*) 'record length  =',nrecl,' spins =',nspin, &
     ' prec flag ',nprec

!!$* reopen WAVECAR
open(unit=10,file='WAVECAR',access='direct',recl=nrecl, &
     iostat=iost,status='old')
if (iost.ne.0) then
    write(80,*) 'WAVECAR reopen error - iostat =',iost            
    stop
endif
read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3), &
     (a3(j),j=1,3)

nwk=nint(xnwk)
nband=nint(xnband)
allocate(occ(nband))
allocate(cener(nband))

write(80,*) 'no. k points =',nwk
write(80,*) 'no. bands =',nband
write(80,*) 'max. energy =',sngl(ecut)
write(80,*) 'real space lattice vectors:'
write(80,*) 'a1 =',(sngl(a1(j)),j=1,3)
write(80,*) 'a2 =',(sngl(a2(j)),j=1,3)
write(80,*) 'a3 =',(sngl(a3(j)),j=1,3)
write(80,*) ' '

call vcross(vtmp,a2,a3)
Vcell=a1(1)*vtmp(1)+a1(2)*vtmp(2)+a1(3)*vtmp(3)
write(80,*) 'unit cell volume =',sngl(Vcell)
call vcross(b1,a2,a3)
call vcross(b2,a3,a1)
call vcross(b3,a1,a2)
do j=1,3
   b1(j)=2.*pi*b1(j)/Vcell
   b2(j)=2.*pi*b2(j)/Vcell
   b3(j)=2.*pi*b3(j)/Vcell
enddo
b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)
write(80,*) 'reciprocal lattice vectors:'
write(80,*) 'b1 =',(sngl(b1(j)),j=1,3)
write(80,*) 'b2 =',(sngl(b2(j)),j=1,3)
write(80,*) 'b3 =',(sngl(b3(j)),j=1,3)

phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag))
call vcross(vtmp,b1,b2)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(vmag*b3mag)
nb1maxA=(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1
nb2maxA=(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
nb3maxA=(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)
      
phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
call vcross(vtmp,b1,b3)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(vmag*b2mag)
phi123=abs(asin(sinphi123))
nb1maxB=(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
nb2maxB=(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
nb3maxB=(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
      
phi23=acos((b2(1)*b3(1)+b2(2)*b3(2)+b2(3)*b3(3))/(b2mag*b3mag))
call vcross(vtmp,b2,b3)
vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(vmag*b1mag)
phi123=abs(asin(sinphi123))
nb1maxC=(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
nb2maxC=(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
nb3maxC=(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1 
npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)

nb1max=max0(nb1maxA,nb1maxB,nb1maxC)
nb2max=max0(nb2maxA,nb2maxB,nb2maxC)
nb3max=max0(nb3maxA,nb3maxB,nb3maxC)
npmax=min0(npmaxA,npmaxB,npmaxC)

allocate (igall(3,npmax))
allocate (hkl(3,npmax))

allocate (coeff_u(npmax,nband,nspin))
allocate (coeff_d(npmax,nband,nspin))
allocate (coeff(2*npmax,nband,nspin))

write(80,*) ' '
WRITE(80,*)
WRITE(80,*) 'read header from WAVECAR file'

call vcross(vecy,vecz,vecx)

if (gamma .le. 0.0d0) then ! normal plane wave as final state
    pureplanewave = .true.
else                       ! treat final state as z-damped plane wave
    pureplanewave = .false.
    ccell         = a3(3)
    OPEN(UNIT=86,FILE='ARUPS_s_PWFS.out',STATUS='UNKNOWN')
    OPEN(UNIT=87,FILE='ARUPS_p_PWFS.out',STATUS='UNKNOWN')
    OPEN(UNIT=88,FILE='ARUPS_u_PWFS.out',STATUS='UNKNOWN')
    OPEN(UNIT=89,FILE='ARUPS_d_PWFS.out',STATUS='UNKNOWN')
if (arupstype .eq. 0) then
    allocate(photoint2_sum(nkx,nek))
    allocate(photointspin2_sum(nkx,nek))
    allocate(photoint2_pol(nkx,nek))
    allocate(photointspin2_pol(nkx,nek))
    allocate(photoint2_nup(nkx,nek))
    allocate(photointspin2_nup(nkx,nek))
    allocate(photoint2_ndo(nkx,nek))
    allocate(photointspin2_ndo(nkx,nek))
elseif (arupstype .eq. 1) then
    allocate(photoint2_sum(nkx,nkx))
    allocate(photointspin2_sum(nkx,nkx))
    allocate(photoint2_pol(nkx,nkx))
    allocate(photointspin2_pol(nkx,nkx))
    allocate(photoint2_nup(nkx,nkx))
    allocate(photointspin2_nup(nkx,nkx))
    allocate(photoint2_ndo(nkx,nkx))
    allocate(photointspin2_ndo(nkx,nkx))
else
    allocate(photoint2_sum(nkx,nephk))
    allocate(photointspin2_sum(nkx,nephk))
    allocate(photoint2_pol(nkx,nephk))
    allocate(photointspin2_pol(nkx,nephk))
    allocate(photoint2_nup(nkx,nephk))
    allocate(photointspin2_nup(nkx,nephk))
    allocate(photoint2_ndo(nkx,nephk))
    allocate(photointspin2_ndo(nkx,nephk))
end if
photoint2_sum(:,:) = 0.0
photointspin2_sum(:,:) = 0.0
photoint2_pol(:,:) = 0.0
photointspin2_pol(:,:) = 0.0
photoint2_nup(:,:) = 0.0
photointspin2_nup(:,:) = 0.0
photoint2_ndo(:,:) = 0.0
photointspin2_ndo(:,:) = 0.0        
end if

epsilon = 1.0d-8
ekintok = 0.512316807
imag    = (0.d0,1.0d0)

iwindow = 4
ekinmin = eklo - iwindow*ebroad
      
      
if (arupstype .eq. 0) then                             ! 0 = Energy versus k
   ekinmax = ekhi + iwindow*ebroad  
   write(80,*) 'Entering energy vs. k plot mode ...'    
elseif (arupstype .eq. 1) then                         ! kx-ky plot
   ekinmax = eklo + iwindow*ebroad
   write(80,*) 'Entering kx-ky plot mode ...'
else                                                   ! CIS scan
   ekinmax = eklo + iwindow*ebroad
   write(80,*) 'Entering CIS plot mode ...'
end if

write(80,*)
write(80,*) 'Total number of k-points in WAVECAR : ',nwk
write(80,*) '                Summation runs over : ',1, nwk
write(80,*)

!!$* beginn loop over: spin, k
irec=2
do isp=1,nspin
    do iwk=1,nwk                                   ! loop over q points in BZ
        !write(6,*) iwk,' of ', nwk
        irec=irec+1
        read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
            (cener(iband),occ(iband),iband=1,nband)
        
        nplane  = nint(xnplane)
        nplane2 = nplane / 2
        
        write(80,*) 'k point #',iwk,'  input no. of plane waves =', &
            nplane2
        write(80,*) 'k value =',(sngl(wk(j)),j=1,3)
        do iband=1,nband
            irec=irec+1
            !read(unit=10,rec=irec) (coeff_u(iplane,iband,isp), &
            !    iplane=1,nplane2), (coeff_d(iplane,iband,isp), &
            !    iplane=1+nplane2,nplane)
            read(unit=10,rec=irec) (coeff(iplane,iband,isp), &
                iplane=1,nplane)
        enddo
        
        ncnt=0
        do ig3=0,2*nb3max
            ig3p=ig3
            if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1
            do ig2=0,2*nb2max
                ig2p=ig2
                if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
                do ig1=0,2*nb1max
                    ig1p=ig1
                    if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
                    do j=1,3
                        sumkg(j)=(wk(1)+ig1p)*b1(j)+ &
                            (wk(2)+ig2p)*b2(j)+(wk(3)+ig3p)*b3(j)
                    enddo
                    gtot=sqrt(sumkg(1)**2+sumkg(2)**2+sumkg(3)**2)
                    etot=gtot**2/c
                    if (etot.lt.ecut) then
                        ncnt=ncnt+1
                        igall(1,ncnt)=sumkg(1)
                        igall(2,ncnt)=sumkg(2)
                        igall(3,ncnt)=sumkg(3)
                        hkl(1,ncnt)=ig1p
                        hkl(2,ncnt)=ig2p
                        hkl(3,ncnt)=ig3p
                    end if
                enddo
            enddo
        enddo
        if ((2*ncnt).ne.nplane) then
            write(80,*) '*** error - computed no. ',ncnt,' != input no. ',nplane
            stop
        endif
        if (ncnt.gt.npmax) then
            write(80,*) '*** error - plane wave count exceeds estimate'
            stop
        endif

        nsymmetry = 0

        do bandcount = bandlow, bandhi
            eigval  = real(cener(bandcount))
            ebind   = -efermi + eigval
            
            if ( temperature .le. 0.0d0 ) then
                fermifac = 1.0d0
            else
                fermifac = 1.0d0/(exp( (eigval - efermi)/temperature) + 1.0d0) 
            endif

            if     (arupstype .eq. 0 ) then     ! bandmap

                efinal = ephotonmin - ( phiwork - ebind )
                if (ebind.gt.ekinmin .and. ebind.lt.ekinmax) then   ! binding energy in desired range
                    iekstart = max(1  ,1+floor((ebind - iwindow*ebroad - eklo)/ekstep))
                    iekend   = min(nek,1+floor((ebind + iwindow*ebroad - eklo)/ekstep))
                    
                    if (iekend .ge. iekstart) then
                        write(80,'(a3,i3,a6,i1,a4,i4,a4,f9.4,a6,f9.4,a10,2i4,a9,f9.4)') &
                            'iq=',iwk,' qsym=',nsymmetry,' ib=',bandcount,' Ei=',ebind,&
                            ' Ekin=',efinal,' iekrange=',iekstart,iekend, ' fermifac', fermifac 
                        do iek = iekstart, iekend
                            ekarr  = ephotonmin - ( phiwork - (eklo + (iek - 1)*ekstep))    ! kinetic energy in (eV)
                            kfinal = ekintok*sqrt(ekarr)       ! final wave vector in (1/A) 
                            deltae = efinal - ekarr
                            l1     = (1/(ebroad*sqrt(two_pi)))*exp(-deltae*deltae/(2*ebroad*ebroad))      ! gauss broadening
                            
                            do jg = 1, nplane2   ! loop over G-vectors
                                do j = 1,3
                                    qplusG(j) =  igall(j,jg)   
                                end do 
                                
                                hmil = hkl(1,jg)
                                kmil = hkl(2,jg)
                                lmil = hkl(3,jg)         
                                qplusGlen=sqrt(qplusG(1)**2+qplusG(2)**2+qplusG(3)**2)
                                if ( (qplusGlen.ge.(kfinal - iwindow*kxbroad)) .and. &
                                     (qplusGlen.le.(kfinal + iwindow*kxbroad)) ) then 
                                    
                                    do isymmetry = 0, nsymmetry
                                        
                                        qplusGx = 0.0d0
                                        qplusGy = 0.0d0
                                        qplusGz = 0.0d0
                                        do j = 1,3
                                            qplusGx = qplusGx+qplusG(j)*vecx(j)
                                            qplusGy = qplusGy+qplusG(j)*vecy(j)
                                            qplusGz = qplusGz+qplusG(j)*vecz(j)
                                        enddo
                                        
                                        if ((abs(qplusGy) .lt. iwindow*kybroad) .and. (qplusGz.gt.0.0d0)) then ! restriction of ky-values
                                            ikxstart = max(1  ,1+floor((qplusGx - iwindow*kxbroad - kxlo)/kxstep))
                                            ikxend   = min(nkx,1+floor((qplusGx + iwindow*kxbroad - kxlo)/kxstep))
                                            
                                            if (ikxend .ge. ikxstart) then   ! ikxend .ge. ikxstart
                                                if (pureplanewave) then    ! final state is a plane wave
                                                    cgcabs_s = abs(coeff(jg,bandcount,isp) + coeff(jg+nplane2,bandcount,isp))**2.0d0
                                                    cgcabs_p = abs(coeff(jg,bandcount,isp))**2.0d0 - &
                                                               abs(coeff(jg+nplane2,bandcount,isp))**2.0d0
                                                    cgcabs_u = abs(coeff(jg,bandcount,isp))**2.0d0
                                                    cgcabs_d = abs(coeff(jg+nplane2,bandcount,isp))**2.0d0
                                                 else
                                                    cgccomplex_s = (0.0d0,0.0d0)
                                                    cgccomplex_p = (0.0d0,0.0d0)
                                                    cgccomplex_u = (0.0d0,0.0d0)
                                                    cgccomplex_d = (0.0d0,0.0d0)
                                                    
                                                    do jgz = 1,nplane   ! loop over Gz-vectors
                                                        hmilz = hkl(1,jgz)
                                                        kmilz = hkl(2,jgz)
                                                        lmilz = hkl(3,jgz)
                                                        
                                                        if ( (hmil .eq. hmilz) .and. (kmil .eq. kmilz)) then
                                                            lmilz = lmilz - lmil
                                                            gz    = lmilz*b3(3) 
                                                            dnorm = 1/sqrt(Vcell/ccell) * sqrt(1/(ccell-z0+(1/(2*gamma)) &
                                                                - (exp(-2*gamma*z0)/(2*gamma)))) 

                                                            if (lmilz .eq. 0) then ! Gz .eq. 0 
                                                                fgz = dnorm* (ccell - (gamma*z0 + exp(-gamma*z0) - 1.0d0) &
                                                                    /gamma)
                                                            else                   ! Gz .neq. 0
                                                                fgz = dnorm*( (exp(-imag*gz*z0)    - exp(-gamma*z0)) &
                                                                    / (gamma - imag*gz) + (exp(-imag*gz*ccell) & 
                                                                    - exp(-imag*gz*z0))/(-imag*gz))
                                                            end if 
                        cgccomplex_s = cgccomplex_s + fgz*conjg(coeff(jgz,bandcount,isp) + coeff(jgz+nplane2,bandcount,isp))
                        cgccomplex_p = cgccomplex_p + fgz*conjg(coeff(jgz,bandcount,isp) - coeff(jgz+nplane2,bandcount,isp))
                        cgccomplex_u = cgccomplex_u + fgz*conjg(coeff(jgz,bandcount,isp))
                        cgccomplex_d = cgccomplex_d + fgz*conjg(coeff(jgz+nplane2,bandcount,isp))
                                                        endif
                                                        
                                                enddo   ! loop over Gz-vectors
                        cgcabs_s = abs(cgccomplex_s)**2.0d0
                        cgcabs_p = abs(cgccomplex_u)**2.0d0-abs(cgccomplex_d)**2.0d0
                        cgcabs_u = abs(cgccomplex_u)**2.0d0
                        cgcabs_d = abs(cgccomplex_d)**2.0d0
                        cgcabs2_s = abs(coeff(jg,bandcount,isp) + coeff(jg+nplane2,bandcount,isp))**2.0d0
                        cgcabs2_p = abs(coeff(jg,bandcount,isp))**2.0d0 - abs(coeff(jg+nplane2,bandcount,isp))**2.0d0
                        cgcabs2_u = abs(coeff(jg,bandcount,isp))**2.0d0
                        cgcabs2_d = abs(coeff(jg+nplane2,bandcount,isp))**2.0d0
                                                    
                                                endif   ! final state is a plane wave
                                                
	                                            do ikx = ikxstart, ikxend  ! loop over kx values  
	                                                kxarr    = kxlo + (ikx - 1)*kxstep
	                                                kperp    = dsqrt(kfinal*kfinal - kxarr*kxarr)
	                                                
	                                                if (kxarr .le. kfinal .and. kxarr .ge. -kfinal) then   ! k range
	                                                    !deltak   = sqrt((qplusGx - kxarr)*(qplusGx - kxarr) +&
	                                                    !   (qplusGy - 0)    *(qplusGy - 0))
	                                                    deltak   = dsqrt((qplusGx - kxarr)*(qplusGx - kxarr) +&
	                                                        (qplusGy - 0)    *(qplusGy - 0) +&
	                                                        (qplusGz - kperp)*(qplusGz - kperp))
	                                                                                                        
	                                                    l2      = (1/(kxbroad*sqrt(two_pi)))*exp(-deltak*deltak/(2*kxbroad*kxbroad))
	                                                    photoint_sum(ikx,iek) = photoint_sum(ikx,iek) + fermifac*l1*l2*cgcabs_s
	                                                    photoint_pol(ikx,iek) = photoint_pol(ikx,iek) + fermifac*l1*l2*cgcabs_p
	                                                    photoint_nup(ikx,iek) = photoint_nup(ikx,iek) + fermifac*l1*l2*cgcabs_u
	                                                    photoint_ndo(ikx,iek) = photoint_ndo(ikx,iek) + fermifac*l1*l2*cgcabs_d
	                                                    if (.not. pureplanewave) then
	                                                        photoint2_sum(ikx,iek)= photoint2_sum(ikx,iek)+ fermifac*l1*l2*cgcabs2_s
	                                                        photoint2_pol(ikx,iek)= photoint2_pol(ikx,iek)+ fermifac*l1*l2*cgcabs2_p
                                                            photoint2_nup(ikx,iek)= photoint2_nup(ikx,iek)+ fermifac*l1*l2*cgcabs2_u
	                                                        photoint2_ndo(ikx,iek)= photoint2_ndo(ikx,iek)+ fermifac*l1*l2*cgcabs2_d
	                                                    end if  ! final state is a plane wave                                 
	                                                end if  ! k range               
	                                            end do  ! end loop over kx values 
                                            endif   ! ikxend .ge. ikxstart
                                        endif   ! restriction of ky-values
                                    enddo   !nsymmetry
                                endif   !qplusGlen
                            enddo   !jp
                        enddo   !iek
                    endif   !iekend
                endif   ! binding energy in desired range
                    
            elseif (arupstype .eq. 1 ) then     ! momentum map
                if (ebind.gt.ekinmin .and. ebind.lt.ekinmax) then  ! kinetic energy in desired range
                    efinal   = ephotonmin - (phiwork - ebind)
                    ekarr  = eklo + ephotonmin - phiwork
                    kfinal = ekintok*sqrt(ekarr)
                    deltae = efinal - ekarr
                    l1     = (1/(ebroad*sqrt(two_pi)))*exp(-deltae*deltae/(2*ebroad*ebroad))      ! gauss broadening
                    
                    write(80,'(a3,i3,a6,i1,a4,i4,a4,f9.4,a6,f9.4,a10,f9.4)') &
                        'iq=',iwk,' qsym=',nsymmetry,' ib=',bandcount,' Ei=',ebind,&
                        ' Ekin=',efinal,' Ekin_set=',ekarr
                    
                    do  jg = 1, nplane2   ! loop over G-vectors
                        do j = 1,3
                            qplusG(j) =  igall(j,jg)   
                        end do 
                        hmil = hkl(1,jg)
                        kmil = hkl(2,jg)
                        lmil = hkl(3,jg)         
                        qplusGlen=sqrt(qplusG(1)**2+qplusG(2)**2+qplusG(3)**2)
                        
                        if ((qplusGlen.ge.(kfinal - iwindow*kxbroad)) .and. &
                                (qplusGlen.le.(kfinal + iwindow*kxbroad)) ) then     ! qplusGlen
                            do isymmetry = 0, nsymmetry   ! loop to symmetrize +/- q-values
                                qplusGx = 0
                                qplusGy = 0
                                qplusGz = 0
                                do j = 1,3
                                    qplusGx = qplusGx+qplusG(j)*vecx(j)
                                    qplusGy = qplusGy+qplusG(j)*vecy(j)
                                    qplusGz = qplusGz+qplusG(j)*vecz(j)
                                enddo
                                
                                if (qplusGz.gt.0.0d0) then ! restriction to positive kz-values
                                    ikxstart = max(1  ,1+floor((qplusGx - iwindow*kxbroad - kxlo)/kxstep))
                                    ikxend   = min(nkx,1+floor((qplusGx + iwindow*kxbroad - kxlo)/kxstep))
                                    ikystart = max(1  ,1+floor((qplusGy - iwindow*kybroad - kxlo)/kxstep))
                                    ikyend   = min(nkx,1+floor((qplusGy + iwindow*kybroad - kxlo)/kxstep))
                                    
                                                                        
                                                                        
                                    if ((ikxend .ge. ikxstart) .and. (ikyend .ge. ikystart)) then  !  ikx iky cond.   
                                        if (pureplanewave) then    ! final state is a plane wave
                                            cgcabs_s = abs(coeff(jg,bandcount,isp) + coeff(jg+nplane2,bandcount,isp))**2.0d0
                                            cgcabs_p = abs(coeff(jg,bandcount,isp))**2.0d0 - &
                                                       abs(coeff(jg+nplane2,bandcount,isp))**2.0d0
                                            cgcabs_u = abs(coeff(jg,bandcount,isp))**2.0d0
                                            cgcabs_d = abs(coeff(jg+nplane2,bandcount,isp))**2.0d0
                                        else
                                            cgccomplex_s = (0.0d0,0.0d0)
                                            cgccomplex_p = (0.0d0,0.0d0)
                                            cgccomplex_u = (0.0d0,0.0d0)
                                            cgccomplex_d = (0.0d0,0.0d0)
                                            
                                            do jgz = 1,nplane   ! loop over Gz-vectors
                                                hmilz = hkl(1,jgz)
                                                kmilz = hkl(2,jgz)
                                                lmilz = hkl(3,jgz)
                                                
                    if ( (hmil .eq. hmilz) .and. (kmil .eq. kmilz)) then
                        lmilz = lmilz - lmil
                        gz    = lmilz*b3(3) 
                        dnorm = 1/sqrt(Vcell/ccell) * sqrt(1/(ccell-z0+(1/(2*gamma))-(exp(-2*gamma*z0)/(2*gamma)))) 

                        if (lmilz .eq. 0) then ! Gz .eq. 0 
                            fgz = dnorm* (ccell - (gamma*z0 + exp(-gamma*z0) - 1.0d0)/gamma)
                        else                   ! Gz .neq. 0
                            fgz = dnorm*((exp(-imag*gz*z0)    - exp(-gamma*z0)) /(gamma - imag*gz) &
                                + (exp(-imag*gz*ccell) - exp(-imag*gz*z0))/(       -imag*gz))
                        end if 
                        cgccomplex_s = cgccomplex_s + fgz*conjg(coeff(jgz,bandcount,isp) + coeff(jgz+nplane2,bandcount,isp))
                        cgccomplex_p = cgccomplex_p + fgz*conjg(coeff(jgz,bandcount,isp) - coeff(jgz+nplane2,bandcount,isp))
                        cgccomplex_u = cgccomplex_u + fgz*conjg(coeff(jgz,bandcount,isp))
                        cgccomplex_d = cgccomplex_d + fgz*conjg(coeff(jgz+nplane2,bandcount,isp))
                    endif
                                                
                                            enddo   ! loop over Gz-vectors
                                            cgcabs_s = abs(cgccomplex_s)**2.0d0
                                            cgcabs_p = abs(cgccomplex_u)**2.0d0-abs(cgccomplex_d)**2.0d0
                                            cgcabs_u = abs(cgccomplex_u)**2.0d0
                                            cgcabs_d = abs(cgccomplex_d)**2.0d0
                                            cgcabs2_s = abs(coeff(jg,bandcount,isp) + coeff(jg+nplane2,bandcount,isp))**2.0d0
                                            cgcabs2_p = abs(coeff(jg,bandcount,isp))**2.0d0 - &
                                                        abs(coeff(jg+nplane2,bandcount,isp))**2.0d0
                                            cgcabs2_u = abs(coeff(jg,bandcount,isp))**2.0d0
                                            cgcabs2_d = abs(coeff(jg+nplane2,bandcount,isp))**2.0d0
                                                    
                                        endif   !final state is a plane wave
                                        
                                        do ikx = ikxstart, ikxend  ! loop over kx values
                                            
                                            kxarr    = kxlo + (ikx - 1)*kxstep
                                            
                                            do iky = ikystart, ikyend  ! loop over ky values  
                                                
                                                kyarr    = kxlo + (iky - 1)*kxstep                     
                                                kpar     = dsqrt(kxarr*kxarr + kyarr*kyarr)
                                                kperp    = dsqrt(kfinal*kfinal - kpar*kpar)
                                                
                                                if (kpar .le. kfinal) then
                                                    
                                                    deltak   = dsqrt((qplusGx - kxarr)*(qplusGx - kxarr) +&
                                                                (qplusGy - kyarr)*(qplusGy - kyarr) +&
                                                                (qplusGz - kperp)*(qplusGz - kperp))
                                                    
                                                    l2       = (1/(kxbroad*sqrt(two_pi)))*exp(-deltak*deltak/(2*kxbroad*kxbroad))  
                                                    
                                                    !l2      = (1/(kxbroad*kybroad*two_pi))* & 
                                                    !               exp(-(qplusGx - kxarr)*(qplusGx - kxarr)/(2*kxbroad*kxbroad) ) &
                                                    !               * exp(-(qplusGy - kyarr)*(qplusGy - kyarr)/(2*kybroad*kybroad) )    
                                                    
                                                    photoint_sum(ikx,iky) = photoint_sum(ikx,iky) + fermifac*l1*l2*cgcabs_s
                                                    photoint_pol(ikx,iky) = photoint_pol(ikx,iky) + fermifac*l1*l2*cgcabs_p
                                                    photoint_nup(ikx,iky) = photoint_nup(ikx,iky) + fermifac*l1*l2*cgcabs_u
                                                    photoint_ndo(ikx,iky) = photoint_ndo(ikx,iky) + fermifac*l1*l2*cgcabs_d
                                                    if (.not. pureplanewave) then
                                                        photoint2_sum(ikx,iky)= photoint2_sum(ikx,iky)+ fermifac*l1*l2*cgcabs2_s
                                                        photoint2_pol(ikx,iky)= photoint2_pol(ikx,iky)+ fermifac*l1*l2*cgcabs2_p
                                                        photoint2_nup(ikx,iky)= photoint2_nup(ikx,iky)+ fermifac*l1*l2*cgcabs2_u
                                                        photoint2_ndo(ikx,iky)= photoint2_ndo(ikx,iky)+ fermifac*l1*l2*cgcabs2_d
                                                    end if                                        
                                                end if                 
                                            end do  ! end loop over ky values       
                                        end do  ! end loop over kx values  
                                    end if !  ikx iky cond.
                                endif   !restriction to positive kz-values
                            enddo   !loop to symmetrize +/- q-values
                        endif   ! qplusGlen
                    enddo    ! loop over G-vectors
                endif   ! kinetic energy in desired range
                
            elseif( arupstype .eq. 2 ) then       !arupstype = 2, default ! constant initial state
                if (ebind.gt.ekinmin .and. ebind.lt.ekinmax) then  ! kinetic energy in desired range
                    ephotonstart = 1       ! max(1  ,1+floor((efinal - iwindow*ebroad - eklo)/ekstep))
                    ephotonend   = nephk    ! min(nek,1+floor((efinal + iwindow*ebroad - eklo)/ekstep))
                    do ieph = ephotonstart, ephotonend  ! loop over photon energies
                        ephoton = ephotonmin + (ieph-1) * ekstep
                        efinal   = ephoton - (phiwork - ebind)
                        ekarr  = eklo + ephoton - phiwork
                        kfinal = ekintok*sqrt(ekarr)       ! final wave vector in (1/A)
                        deltae = efinal - ekarr
                        l1     = (1/(ebroad*sqrt(two_pi)))*exp(-deltae*deltae/(2*ebroad*ebroad))  ! gauss broadening
                        
                        write(80,'(a3,i3,a3,f8.4,a4,i4,a5,f9.4,a6,f9.4,a10,f9.4)') &
                            'iq=',iwk,' ei=',ebind,' ib=',bandcount,' Eph=',ephoton,&
                            ' Ekin=',efinal,' Ekin_set=',ekarr
                        
                        do  jg = 1, nplane2   ! loop over G-vectors
                            do j = 1,3
                               qplusG(j) =  igall(j,jg)   
                            end do 
                            hmil = hkl(1,jg)
                            kmil = hkl(2,jg)
                            lmil = hkl(3,jg)         
                            qplusGlen=sqrt(qplusG(1)**2+qplusG(2)**2+qplusG(3)**2)
                        
                            if ( (qplusGlen.ge.(kfinal - iwindow*kxbroad)) .and. &
                                 (qplusGlen.le.(kfinal + iwindow*kxbroad)) ) then   !qplusGlen
                                do isymmetry = 0, nsymmetry   ! loop to symmetrize +/- q-values
                                    qplusGx = 0
                                    qplusGy = 0
                                    qplusGz = 0
                                    do j = 1,3
                                        qplusGx = qplusGx+qplusG(j)*vecx(j)
                                        qplusGy = qplusGy+qplusG(j)*vecy(j)
                                        qplusGz = qplusGz+qplusG(j)*vecz(j)
                                    enddo
                                    
                                    if ((abs(qplusGy) .lt. iwindow*kybroad) .and. (qplusGz.gt.0.0d0)) then ! restriction of ky-values
                                        ikxstart = max(1  ,1+floor((qplusGx - iwindow*kxbroad - kxlo)/kxstep))
                                        ikxend   = min(nkx,1+floor((qplusGx + iwindow*kxbroad - kxlo)/kxstep))
                                        
                                        if ((ikxend .ge. ikxstart) .and. (ikyend .ge. ikystart)) then  !  ikx iky cond.   
                                            if (pureplanewave) then    ! final state is a plane wave
                                                cgcabs_s = abs(coeff(jg,bandcount,isp) + coeff(jg+nplane2,bandcount,isp))**2.0d0
                                                cgcabs_p = abs(coeff(jg,bandcount,isp))**2.0d0 - &
                                                           abs(coeff(jg+nplane2,bandcount,isp))**2.0d0
                                                cgcabs_u = abs(coeff(jg,bandcount,isp))**2.0d0
                                                cgcabs_d = abs(coeff(jg+nplane2,bandcount,isp))**2.0d0
                                            else
                                                cgccomplex_s = (0.0d0,0.0d0)
                                                cgccomplex_p = (0.0d0,0.0d0)
                                                cgccomplex_u = (0.0d0,0.0d0)
                                                cgccomplex_d = (0.0d0,0.0d0)
                                                
                                                do jgz = 1,nplane   ! loop over Gz-vectors
                                                    hmilz = hkl(1,jgz)
                                                    kmilz = hkl(2,jgz)
                                                    lmilz = hkl(3,jgz)
                                                    
                                                    if ( (hmil .eq. hmilz) .and. (kmil .eq. kmilz)) then
                                                        lmilz = lmilz - lmil
                                                        gz    = lmilz*b3(3) 
                                                        dnorm = sqrt(1/(ccell-z0+(1/(2*gamma))-(exp(-2*gamma*z0)/(2*gamma)))) 

                                                        if (lmilz .eq. 0) then ! Gz .eq. 0 
                                                            fgz = dnorm* (ccell - (gamma*z0 + exp(-gamma*z0) - 1.0d0)/gamma)
                                                        else                   ! Gz .neq. 0
                                                            fgz = dnorm*((exp(-imag*gz*z0)    - exp(-gamma*z0)) /(gamma - imag*gz) &
                                                                + (exp(-imag*gz*ccell) - exp(-imag*gz*z0))/(       -imag*gz))
                                                        end if 
                                cgccomplex_s = cgccomplex_s + fgz*conjg(coeff(jgz,bandcount,isp) + coeff(jgz+nplane2,bandcount,isp))
                                cgccomplex_p = cgccomplex_p + fgz*conjg(coeff(jgz,bandcount,isp) - coeff(jgz+nplane2,bandcount,isp))
                                cgccomplex_u = cgccomplex_u + fgz*conjg(coeff(jgz,bandcount,isp))
                                cgccomplex_d = cgccomplex_d + fgz*conjg(coeff(jgz+nplane2,bandcount,isp))
                                                    endif
                                                    
                                                enddo   ! loop over Gz-vectors
                                                cgcabs_s = abs(cgccomplex_s)**2.0d0
                                                cgcabs_p = abs(cgccomplex_u)**2.0d0-abs(cgccomplex_d)**2.0d0
                                                cgcabs_u = abs(cgccomplex_u)**2.0d0
                                                cgcabs_d = abs(cgccomplex_d)**2.0d0
                                                cgcabs2_s = abs(coeff(jg,bandcount,isp) + coeff(jg+nplane2,bandcount,isp))**2.0d0 
                                                cgcabs2_p = abs(coeff(jg,bandcount,isp))**2.0d0 - &
                                                            abs(coeff(jg+nplane2,bandcount,isp))**2.0d0
                                                cgcabs2_u = abs(coeff(jg,bandcount,isp))**2.0d0
                                                cgcabs2_d = abs(coeff(jg+nplane2,bandcount,isp))**2.0d0
                                                    
                                            endif   !final state is a plane wave
                                            
                                            do ikx = ikxstart, ikxend  ! loop over kx values  
                                                
                                                kxarr    = kxlo + (ikx - 1)*kxstep
                                                kperp    = dsqrt(kfinal*kfinal - kxarr*kxarr)
                                                
                                                if (kxarr .le. kfinal .and. kxarr .ge. -kfinal) then
                                                    
                                                    !deltak   = sqrt((qplusGx - kxarr)*(qplusGx - kxarr) +&
                                                    !                (qplusGy - 0)    *(qplusGy - 0))
                                                    deltak   = dsqrt((qplusGx - kxarr)*(qplusGx - kxarr) +&
                                                                     (qplusGy - 0)    *(qplusGy - 0) +&
                                                                     (qplusGz - kperp)*(qplusGz - kperp))
                                                    
                                                    l2 = (1/(kxbroad*sqrt(two_pi)))*exp(-deltak*deltak/(2*kxbroad*kxbroad))
                                                     
                                                    photoint_sum(ikx,ieph) = photoint_sum(ikx,ieph) + fermifac*l1*l2*cgcabs_s
                                                    photoint_pol(ikx,ieph) = photoint_pol(ikx,ieph) + fermifac*l1*l2*cgcabs_p
                                                    photoint_nup(ikx,ieph) = photoint_nup(ikx,ieph) + fermifac*l1*l2*cgcabs_u
                                                    photoint_ndo(ikx,ieph) = photoint_ndo(ikx,ieph) + fermifac*l1*l2*cgcabs_d
                                                    if (.not. pureplanewave) then
                                                        photoint2_sum(ikx,ieph)= photoint2_sum(ikx,ieph)+ fermifac*l1*l2*cgcabs2_s
                                                        photoint2_pol(ikx,ieph)= photoint2_pol(ikx,ieph)+ fermifac*l1*l2*cgcabs2_p
                                                        photoint2_nup(ikx,ieph)= photoint2_nup(ikx,ieph)+ fermifac*l1*l2*cgcabs2_u
                                                        photoint2_ndo(ikx,ieph)= photoint2_ndo(ikx,ieph)+ fermifac*l1*l2*cgcabs2_d
                                                    end if                                      
                                                end if                       
                                            end do  ! end loop over kx values
                                        endif   !  ikx iky cond.
                                    endif   ! restriction of ky-values
                                enddo ! loop to symmetrize
                            endif !qplusGlen
                        enddo    ! loop over G-vectors
                    enddo   ! loop over photon energies
                endif ! kinetic energy in desired range
            endif       !arupstype
            
        enddo   !bandcount
    enddo   !iwk

    photointspin_sum = photointspin_sum + photoint_sum
    photointspin_pol = photointspin_pol + photoint_pol 
    photointspin_ndo = photointspin_ndo + photoint_ndo
    photointspin_nup = photointspin_nup + photoint_nup
    if (.not. pureplanewave) then 
        photointspin2_sum = photointspin2_sum + photoint2_sum
        photointspin2_pol = photointspin2_pol + photoint2_pol 
        photointspin2_ndo = photointspin2_ndo + photoint2_ndo
        photointspin2_nup = photointspin2_nup + photoint2_nup
    end if

enddo   !isp

do i = 82,89
    if ( (pureplanewave) .and. (i.ge.86)) then
        EXIT
    endif
    write(i,'(I5,t35,A)') nwk,           '! number of k-points'
    write(i,'(2I5,t35,A)') bandlow,bandhi,'! first and last bandindex for summation over initial states'
    write(i,'(f10.5,t35,A)') efermi,      '! Fermi-level (energy of HOMO) in [eV]'
    write(i,'(f10.5,t35,A)') temperature, '! Temperature in [eV]'
    write(i,'(f10.5,t35,A)') phiwork,     '! work function in [eV]'
    
    if (arupstype .eq.1) then
        write(i,'(2f10.5,t35,A)') ephotonmin, eklo, '! energy of incident photon [eV], binding energy [eV]'
    else
        write(i,'(f10.5,t35,A)') ephotonmin, '! energy of incident photon [eV], '
    end if
    
    write(i,'(f10.5,t35,A)') ebroad,      '! energy broadening [eV]'
    write(i,'(2f10.5,t35,A)') kxbroad,kybroad,'! broadening in k-space [A^-1]'
    write(i,'(i10,t35,A)') arupstype,     '! 0 = Energy versus k, 1 = kx-ky plot at constant energy, 2 = CIS scan'
    write(i,'(3f10.5,t35,A)') vecz,       '! unit vector for normal emission (z-axis)'
    write(i,'(3f10.5,t35,A)') vecx,       '! unit vector for grazing emission (x-axis)' 

    if (arupstype .eq.0) then
        write(i,'(2I5)')    nek,nkx
        write(i,'(3F15.9)') eklo,ekhi,ekstep 
        write(i,'(3F15.9)') kxlo,kxhi,kxstep 
    elseif (arupstype .eq.1) then
        write(i,'(2I5)')    nkx,nkx
        write(i,'(3F15.9)') kxlo,kxhi,kxstep       
        write(i,'(3F15.9)') kxlo,kxhi,kxstep    
    else
        write(i,'(2I5)')    nephk,nkx
        write(i,'(3F15.9)') ephotonmin,ephotonmax,ekstep       
        write(i,'(3F15.9)') kxlo,kxhi,kxstep  
    end if
enddo   !i

write(82,'(6(f16.10,2X))') photointspin_sum
write(83,'(6(f16.10,2X))') photointspin_pol
write(84,'(6(f16.10,2X))') photointspin_nup
write(85,'(6(f16.10,2X))') photointspin_ndo

if (.not. pureplanewave) then
    write(86,'(6(f16.10,2X))') photointspin2_sum
    write(87,'(6(f16.10,2X))') photointspin2_pol
    write(88,'(6(f16.10,2X))') photointspin2_nup
    write(89,'(6(f16.10,2X))') photointspin2_ndo
endif

CLOSE(UNIT=80)
CLOSE(UNIT=82)
CLOSE(UNIT=83)
CLOSE(UNIT=84)
CLOSE(UNIT=85)
CLOSE(UNIT=86)
CLOSE(UNIT=87)
CLOSE(UNIT=88)
CLOSE(UNIT=89)

endtime = wtime ( ); ! get wall time
!write(6,*) 'ELAPSED WALL TIME = ',endtime - starttime,' seconds'

WRITE(80,*)
WRITE(80,*) 'ELAPSED WALL TIME = ',endtime - starttime,' seconds'


WRITE(80,*)
WRITE(80,*) 'DONE WITH ARUPS'
WRITE(80,*) '(for now)'



end program

subroutine vcross(a,b,c)
  implicit real*8(a-h,o-z)
  dimension a(3),b(3),c(3)
  
  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)
  return
end subroutine vcross

function wtime ( )

!*****************************************************************************80
!
!! WTIME returns a reading of the wall clock time.
!
!  Discussion:
!
!    To get the elapsed wall clock time, call WTIME before and after a given
!    operation, and subtract the first reading from the second.
!
!    This function is meant to suggest the similar routines:
!
!      "omp_get_wtime ( )" in OpenMP,
!      "MPI_Wtime ( )" in MPI,
!      and "tic" and "toc" in MATLAB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) WTIME, the wall clock reading, in seconds.
!
  implicit none

  integer ( kind = 4 ) clock_max
  integer ( kind = 4 ) clock_rate
  integer ( kind = 4 ) clock_reading
  real ( kind = 8 ) wtime

  call system_clock ( clock_reading, clock_rate, clock_max )

  wtime = real ( clock_reading, kind = 8 ) &
        / real ( clock_rate, kind = 8 )

  return
end
