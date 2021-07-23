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
!!$*
!!$*
!!$*   note that the energy eigenvalues are complex, as provided in the
!!$*   WAVECAR file, but the imaginary part is zero (at least for all cases
!!$*   investigated thus far)
!!$*     
program arups_MPI  
                                                                    
implicit real*8 (a-h, o-z)    

! Start: MPI configuration
include 'mpif.h'                                    
INTEGER :: numcpu, rank, ierr
! End MPI

complex*8, allocatable :: coeff(:,:,:)
complex*16, allocatable :: cener(:)
real*8, allocatable :: occ(:),igall(:,:)
integer, allocatable :: hkl(:,:) !,weights(:)
dimension a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),vtmp(3),sumkg(3),wk(3)


character (1000) :: outfile
integer              :: kcount,nparwavecar,khi !,numk
integer              :: bandcount,bandlow,bandhi
integer              :: nkx,nek,ikx,iky,iek,iekstart,iekend
integer              :: iwindow,ephotonstart,ephotonend,nephk
integer              :: jg , hmil , kmil,  lmil , line
integer              :: jgz, hmilz, kmilz, lmilz      
integer              :: arupstype,isymmetry,nsymmetry
integer              :: ikxstart,ikxend,ikystart,ikyend
real*8               :: two_pi,ekintok,epsilon
real*8               :: efermi,temperature,fermifac,phiwork
real*8               :: ephotonmin,ephotonmax,ebroad,kxbroad,kybroad,eigval
real*8               :: photonangle,evector(3),qplusG2(3)
real*8               :: vecz(3),vecx(3),vecy(3)  
real*8               :: kxlo,kxhi,kxstep
real*8               :: eklo,ekhi,ekstep,ekinmin,ekinmax
real*8  ,allocatable :: photoint(:,:),photoint2(:,:),photointspin(:,:),photointspin2(:,:)  
real*8  ,allocatable :: photoint_rank(:,:)
real*8  ,allocatable :: photoint2_rank(:,:)
real*8               :: qvec(3),qweight,qplusG(3),qplusGx,qplusGy,qplusGz,qplusGlen
real*8               :: kfvec(3)
real*8               :: bvec(3,3)
real*8               :: deltae,deltak,ebind,efinal,kfinal,kperp
real*8               :: cgcabs,cgcabs2
real*8               :: gz,gamma,z0,ccell
complex*16           :: fgz,cgccomplex(3),imag
real*8               :: l1,l2,kxarr,kyarr,kpar,ekarr      
logical              :: pureplanewave
real*8               :: starttime, endtime

     
!!$*   constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 (value is
!!$*   adjusted in final decimal places to agree with VASP value)

data c/0.262465831d0/
pi=4.*atan(1.)
      
!!$*   input

! We parallelize the loop over G-vectors
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numcpu, ierr)

if (rank .eq. 0) then
    starttime = wtime ( ); ! get wall time
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! read weights from IBZKPT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!OPEN(UNIT=85,FILE='IBZKPT',STATUS='old',action='read', iostat=iost) 
!if (iost.ne.0) write(80,*) 'IBZKPT open error - iostat =',iost
!read(85,*) 
!read(85,*) numk
!read(85,*)
!allocate (weights(numk))
!do line = 1,numk
!   read(85,*) dummy, dummy, dummy, weights(line)
!end do
!numk = sum(weights)
!close(unit=85)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


OPEN(UNIT=81,FILE='ARUPS.in', STATUS='UNKNOWN')

      
two_pi    =  2.0d0*3.14159265d0

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

if (rank .eq. 0) then
    OPEN(UNIT=80,FILE='ARUPS.log',STATUS='UNKNOWN')
    OPEN(UNIT=82,FILE=TRIM(ADJUSTL(outfile)),STATUS='UNKNOWN')
    WRITE(80,*) 'CALCULATION OF ARUPS SPECTRA FROM WAVECAR FILE' 
    if     (arupstype .eq. 0) then
        WRITE(80,*) 'We are using ',numcpu,' CPUs to parallelize the G-vector loop'
    elseif (arupstype .eq. 1) then
        WRITE(80,*) 'We are using ',numcpu,' CPUs to parallelize the G-vector loop'
    else
        WRITE(80,*) 'We are using ',numcpu,' CPUs to parallelize loop over photon energies'
    end if
end if

nkx = (kxhi - kxlo)/kxstep + 1
nek = (ekhi - eklo)/ekstep + 1 
nephk = (ephotonmax - ephotonmin)/ekstep + 1

if (arupstype .eq. 0) then
    allocate(photoint(nkx,nek))
    allocate(photoint_rank(nkx,nek))
    allocate(photointspin(nkx,nek))
elseif (arupstype .eq. 1) then
    allocate(photoint(nkx,nkx))
    allocate(photoint_rank(nkx,nkx))
    allocate(photointspin(nkx,nkx))
else
    allocate(photoint(nkx,nephk))
    allocate(photoint_rank(nkx,nephk))
    allocate(photointspin(nkx,nephk))
end if
photoint(:,:) = 0.0
photoint_rank(:,:) = 0.0
photointspin(:,:) = 0.0



nrecl=24
open(unit=10,file='WAVECAR',access='direct',recl=nrecl, &
     iostat=iost,status='old')
if (iost.ne.0) write(80,*) 'WAVECAR open error - iostat =',iost            
read(unit=10,rec=1) xnrecl,xnspin,xnprec
close(unit=10)

nrecl=nint(xnrecl)
nspin=nint(xnspin)
nprec=nint(xnprec)


if (rank .eq. 0) then
    if(nprec.eq.45210) then
        write(80,*) '*** error - WAVECAR_double requires complex*16'
        stop
    endif

    write(80,*) 
    write(80,*) 'record length  =',nrecl,' spins =',nspin, &
        ' prec flag ',nprec
end if

open(unit=10,file='WAVECAR',access='direct',recl=nrecl, &
     iostat=iost,status='old')

if (rank .eq. 0) then
    if (iost.ne.0) write(80,*) 'open error - iostat =',iost  
end if

read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3), &
     (a3(j),j=1,3)

nwk=nint(xnwk)
nband=nint(xnband)
allocate(occ(nband))
allocate(cener(nband))

if (rank .eq. 0) then
    write(80,*) 'no. k points =',nwk
    write(80,*) 'no. bands =',nband
    write(80,*) 'max. energy =',sngl(ecut)
    write(80,*) 'real space lattice vectors:'
    write(80,*) 'a1 =',(sngl(a1(j)),j=1,3)
    write(80,*) 'a2 =',(sngl(a2(j)),j=1,3)
    write(80,*) 'a3 =',(sngl(a3(j)),j=1,3)
    write(80,*) ' '
end if


call vcross(vtmp,a2,a3)
Vcell=a1(1)*vtmp(1)+a1(2)*vtmp(2)+a1(3)*vtmp(3)

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

if (rank .eq. 0) then
    write(80,*) 'unit cell volume =',sngl(Vcell)
    write(80,*) 'reciprocal lattice vectors:'
    write(80,*) 'b1 =',(sngl(b1(j)),j=1,3)
    write(80,*) 'b2 =',(sngl(b2(j)),j=1,3)
    write(80,*) 'b3 =',(sngl(b3(j)),j=1,3)
end if

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


allocate (coeff(npmax,nband,nspin))


if (rank .eq. 0) then
    write(80,*) ' '
    WRITE(80,*)
    WRITE(80,*) 'read header from WAVECAR file'
end if


call vcross(vecy,vecz,vecx)


if (gamma .le. 0.0d0) then ! normal plane wave as final state
    pureplanewave = .true.
else                       ! treat final state as z-damped plane wave
    pureplanewave = .false.
    ccell         = a3(3)
    if (rank .eq. 0) then
        OPEN(UNIT=83,FILE='ARUPS_PWFS.out',STATUS='UNKNOWN')
    end if
    if (arupstype .eq. 0) then
        allocate(photoint2(nkx,nek))
        allocate(photoint2_rank(nkx,nek))
        allocate(photointspin2(nkx,nek))
    elseif (arupstype .eq. 1) then
        allocate(photoint2(nkx,nkx))
        allocate(photoint2_rank(nkx,nkx))
        allocate(photointspin2(nkx,nkx))
    else
        allocate(photoint2(nkx,nephk))
        allocate(photoint2_rank(nkx,nephk))
        allocate(photointspin2(nkx,nephk)) 
    end if
    photoint2(:,:) = 0.0
    photoint2_rank(:,:) = 0.0
    photointspin2(:,:) = 0.0          
end if

epsilon = 1.0d-8
ekintok = 0.512316807
imag    = (0.d0,1.0d0)

iwindow = 4
ekinmin = eklo - iwindow*ebroad
      
      
if (arupstype .eq. 0) then                             ! 0 = Energy versus k
    ekinmax = ekhi + iwindow*ebroad  
    if (rank .eq. 0) then
        write(80,*) 'Entering energy vs. k plot mode ...'    
    end if
elseif (arupstype .eq. 1) then                         ! kx-ky plot
    ekinmax = eklo + iwindow*ebroad
    if (rank .eq. 0) then
        write(80,*) 'Entering kx-ky plot mode ...'
    end if
else                                                   ! CIS scan
    ekinmax = eklo + iwindow*ebroad
    if (rank .eq. 0) then
        write(80,*) 'Entering CIS plot mode ...'
    end if
end if

if (rank .eq. 0) then
    write(80,*)
    write(80,*) 'Total number of k-points in WAVECAR : ',nwk
    write(80,*) '                Summation runs over : ',1, nwk
    write(80,*)
end if


irec=2
do isp=1,nspin
   do iwk=1,nwk                                   ! loop over q points in BZ
      irec=irec+1
      read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
           (cener(iband),occ(iband),iband=1,nband)
      nplane=nint(xnplane)

      if (rank .eq. 0) then  
        write(80,*) 'k point #',iwk,'  input no. of plane waves =', nplane
        write(80,*) 'k value =',(sngl(wk(j)),j=1,3)
      end if

      do iband=1,nband
         irec=irec+1
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
      if (ncnt.ne.nplane) then
        if (rank .eq. 0) then
         write(6,*) '*** warning - computed no. != input no.'
         write(6,*) 'iwk = ',iwk,' ncnt = ',ncnt,' nplane = ',nplane
	 nplane = min(ncnt,nplane)
         ! stop
        endif
      endif
      if (ncnt.gt.npmax) then
         write(6,*) '*** error - plane wave count exceeds estimate'
         stop
      endif
      
      if ( wk(1) .eq. 0 .and. wk(2).eq. 0  .and. wk(3).eq. 0) then
         nsymmetry = 0
      else
         nsymmetry = 1
      end if    
      
      do bandcount = bandlow, bandhi ! loop over bands
        
        eigval   = real(cener(bandcount))
        ebind    = -efermi + eigval

        if ( temperature .le. 0.0d0) then
          fermifac = 1.0d0
        else
          fermifac = 1.0d0/(exp( (eigval - efermi)/temperature) + 1.0d0) 
        end if       

        if (arupstype .eq. 0) then           ! 0 = Energy versus k

          efinal   = ephotonmin - (phiwork - ebind)
          if (ebind.gt.ekinmin .and. ebind.lt.ekinmax) then                ! binding energy in desired range

            iekstart = max(1  ,1+floor((ebind - iwindow*ebroad - eklo)/ekstep))
            iekend   = min(nek,1+floor((ebind + iwindow*ebroad - eklo)/ekstep))


            if (iekend .ge. iekstart) then        

                if (rank .eq. 0) then
                    write(80,'(a3,i3,a6,i1,a4,i4,a4,f9.4,a6,f9.4,a10,2i4,a9,f9.4)') &
                     'iq=',iwk,' qsym=',nsymmetry,' ib=',bandcount,' Ei=',ebind,&
                     ' Ekin=',efinal,' iekrange=',iekstart,iekend, ' fermifac', fermifac 
                end if
   
              do iek = iekstart, iekend   ! loop over kinetic energies           

!                efinal = ephotonmin - (phiwork - ebind)
                ekarr  = ephotonmin - ( phiwork - (eklo + (iek - 1)*ekstep))    ! kinetic energy in (eV)
                            


                kfinal = ekintok*sqrt(ekarr)       ! final wave vector in (1/A)           
                deltae = efinal - ekarr
                l1     = (1/(ebroad*sqrt(two_pi)))*exp(-deltae*deltae/(2*ebroad*ebroad))      ! gauss broadening


!                write(80,'(a4,i4,a4,f9.4,a7,f9.4,a6,f9.4,a5,i4,a9,f9.4,a7,f9.4)') &
!                     ' ib=',bandcount,' Ei=',ebind, ' ekarr=',ekarr,&
!                     ' Ekin=',efinal,' iek=',iek, ' fermifac', fermifac, ' gauss=', l1


                do  jg = rank+1, nplane, numcpu   ! distribute loop over G-vectors to numcpu processes
                  ! write(*,*) rank,jg
              
                  do j = 1,3
                    qplusG(j) =  igall(j,jg)   
                  end do 

                  hmil = hkl(1,jg)
                  kmil = hkl(2,jg)
                  lmil = hkl(3,jg)         
                  qplusGlen=sqrt(qplusG(1)**2+qplusG(2)**2+qplusG(3)**2)
              
                  if ( (qplusGlen.ge.(kfinal - iwindow*kxbroad)) .and. &
                       (qplusGlen.le.(kfinal + iwindow*kxbroad)) ) then 
     
                    do isymmetry = 0, nsymmetry   ! loop to symmetrize +/- q-values
      
                      if (isymmetry .eq. 0) then
                        qplusG =   qplusG          
                      else
                        qplusG =   -qplusG                          
                      end if
                                                
                      qplusGy   = qplusG(1)*vecy(1) + &
                                  qplusG(2)*vecy(2) + &
                                  qplusG(3)*vecy(3)
                      qplusGz   = qplusG(1)*vecz(1) + &
                                  qplusG(2)*vecz(2) + &
                                  qplusG(3)*vecz(3)
         
                      if ((abs(qplusGy) .lt. iwindow*kybroad) .and. (qplusGz.gt.0.0d0)) then ! restriction of ky-values

                        qplusGx   = qplusG(1)*vecx(1) + &
                                    qplusG(2)*vecx(2) + &
                                    qplusG(3)*vecx(3)

                        ikxstart = max(1  ,1+floor((qplusGx - iwindow*kxbroad - kxlo)/kxstep))
                        ikxend   = min(nkx,1+floor((qplusGx + iwindow*kxbroad - kxlo)/kxstep))
 
                        if (ikxend .ge. ikxstart) then 
                
                          if (pureplanewave) then    ! final state is a plane wave
                
                            cgcabs = abs(coeff(jg,bandcount,isp))**2.0d0
                     
                          else     ! final state is a damped plane wave (exp(gamma*z) factor)
                
                            cgccomplex = (0.0d0,0.0d0)
                     
                            do  jgz = 1, nplane   ! loop over Gz-vectors
                              hmilz = hkl(1,jgz)
                              kmilz = hkl(2,jgz)
                              lmilz = hkl(3,jgz)
                              if ( (hmil .eq. hmilz) .and. (kmil .eq. kmilz)) then
                     
                                lmilz = lmilz - lmil
                                gz    = lmilz*b3(3) 
                                dnorm = 1/sqrt(Vcell/ccell)  &
                                        * sqrt(1/(ccell-z0+(1/(2*gamma))-(exp(-2*gamma*z0)/(2*gamma)))) 
                        
                                if (lmilz .eq. 0) then ! Gz .eq. 0 
                        
                                  fgz = dnorm* (ccell - (gamma*z0 + exp(-gamma*z0) - 1.0d0)/gamma)
                        
                                else                   ! Gz .neq. 0
                            
                                  fgz = dnorm*((exp(-imag*gz*z0)    - exp(-gamma*z0))  /(gamma - imag*gz) &
                                            + (exp(-imag*gz*ccell) - exp(-imag*gz*z0))/(       -imag*gz))
                            
                                end if 
                     
                                cgccomplex(1) = cgccomplex(1) &
                                                + fgz*conjg(coeff(jgz,bandcount,isp))
                       
                              end if
                            end do                               ! end loop over Gz-vectors
                
                            cgcabs = abs(cgccomplex(1))**2.0d0
                            cgcabs2= abs(coeff(jg,bandcount,isp))**2.0d0
                          
                          end if   ! decision on type of final state                

                          do ikx = ikxstart, ikxend  ! loop over kx values  
                          
                            kxarr    = kxlo + (ikx - 1)*kxstep
                            kperp    = dsqrt(kfinal*kfinal - kxarr*kxarr)

                            if (kxarr .le. kfinal .and. kxarr .ge. -kfinal) then   ! k range
                     
!                              deltak   = sqrt((qplusGx - kxarr)*(qplusGx - kxarr) +&
!                                              (qplusGy - 0)    *(qplusGy - 0))
                             deltak   = dsqrt((qplusGx - kxarr)*(qplusGx - kxarr) +&
                                            (qplusGy - 0)    *(qplusGy - 0) +&
                                            (qplusGz - kperp)*(qplusGz - kperp))
                                                                                                                                     
                              l2                = (1/(kxbroad*sqrt(two_pi)))*exp(-deltak*deltak/(2*kxbroad*kxbroad))
                              photoint_rank(ikx,iek) = photoint_rank(ikx,iek) + fermifac*l1*l2*cgcabs
                              if (.not. pureplanewave) then
                                photoint2_rank(ikx,iek)= photoint2_rank(ikx,iek)+ fermifac*l1*l2*cgcabs2
                              end if                                     
                            end if  ! k range               
                          end do  ! end loop over kx values                
                        end if  ! ikxend .ge. ikxstart
                      end if  ! restriction of ky-values
                    end do  ! loop to symmetrize +/- q-values
                  end if          
                end do  ! end loop over G-vectors  
                                 
              end do ! end loop over kinetic energies
            end if ! iekend .ge. iekstart
          end if  ! kinetic energy in desired range

        elseif (arupstype .eq. 1) then 

          if (ebind.gt.ekinmin .and. ebind.lt.ekinmax) then  ! kinetic energy in desired range

            efinal   = ephotonmin - (phiwork - ebind)
            ekarr  = eklo + ephotonmin - phiwork
            kfinal = ekintok*sqrt(ekarr)
            deltae = efinal - ekarr
            l1     = (1/(ebroad*sqrt(two_pi)))*exp(-deltae*deltae/(2*ebroad*ebroad))      ! gauss broadening
          
            if (rank .eq. 0) then
                write(80,'(a3,i3,a6,i1,a4,i4,a4,f9.4,a6,f9.4,a10,f9.4)') &
                  'iq=',iwk,' qsym=',nsymmetry,' ib=',bandcount,' Ei=',ebind,&
                  ' Ekin=',efinal,' Ekin_set=',ekarr
            end if

            do  jg = rank+1, nplane, numcpu   ! distribute loop over G-vectors to numcpu processes

              do j = 1,3
                qplusG(j) =  igall(j,jg)   
              end do 
              hmil = hkl(1,jg)
              kmil = hkl(2,jg)
              lmil = hkl(3,jg)         
              qplusGlen=sqrt(qplusG(1)**2+qplusG(2)**2+qplusG(3)**2)

              if ( (qplusGlen.ge.(kfinal - iwindow*kxbroad)) .and. &
                   (qplusGlen.le.(kfinal + iwindow*kxbroad)) ) then     ! qplusGlen

                do isymmetry = 0, nsymmetry   ! loop to symmetrize +/- q-values

                  if (isymmetry .eq. 0) then
                    qplusG =   qplusG                
                  else 
                    qplusG =   -qplusG                                   
                  end if

                  qplusGz   = qplusG(1)*vecz(1) + &
                              qplusG(2)*vecz(2) + &
                              qplusG(3)*vecz(3)

                  if (qplusGz.gt.0.0d0) then ! restriction to positive kz-values


                    qplusGx   = qplusG(1)*vecx(1) + &
                                qplusG(2)*vecx(2) + &
                                qplusG(3)*vecx(3)
                    qplusGy   = qplusG(1)*vecy(1) + &
                                qplusG(2)*vecy(2) + &
                                qplusG(3)*vecy(3)

                    ikxstart = max(1  ,1+floor((qplusGx - iwindow*kxbroad - kxlo)/kxstep))
                    ikxend   = min(nkx,1+floor((qplusGx + iwindow*kxbroad - kxlo)/kxstep))
                    ikystart = max(1  ,1+floor((qplusGy - iwindow*kybroad - kxlo)/kxstep))
                    ikyend   = min(nkx,1+floor((qplusGy + iwindow*kybroad - kxlo)/kxstep))

                    if ((ikxend .ge. ikxstart) .and. (ikyend .ge. ikystart)) then  !  ikx iky cond.

                      if (pureplanewave) then    ! final state is a plane wave

                        cgcabs = abs(coeff(jg,bandcount,isp))**2.0d0


                      else     ! final state is a damped plane wave (exp(gamma*z) factor)

                        cgccomplex = (0.0d0,0.0d0)

                        do  jgz = 1, nplane   ! loop over Gz-vectors
                          hmilz = hkl(1,jgz)
                          kmilz = hkl(2,jgz)
                          lmilz = hkl(3,jgz)
                          if ( (hmil .eq. hmilz) .and. (kmil .eq. kmilz)) then
                     
                            lmilz = lmilz - lmil
                            gz    = lmilz*b3(3) 
                            dnorm = 1/sqrt(Vcell/ccell)  &
                                    * sqrt(1/(ccell-z0+(1/(2*gamma))-(exp(-2*gamma*z0)/(2*gamma))))                        

                            if (lmilz .eq. 0) then ! Gz .eq. 0 
                        
                              fgz = dnorm*(ccell - (gamma*z0 + exp(-gamma*z0) - 1.0d0)/gamma)
                        
                            else                   ! Gz .neq. 0
                            
                              fgz =dnorm* ((exp(-imag*gz*z0)    - exp(-gamma*z0))  /(gamma - imag*gz) &
                                         + (exp(-imag*gz*ccell) - exp(-imag*gz*z0))/(      -imag*gz))
                             
                            end if 
                     
                            cgccomplex(1) = cgccomplex(1) &
                                           + fgz*conjg(coeff(jgz,bandcount,isp))
          
                          end if
                        end do                               ! end loop over Gz-vectors
                
                        cgcabs =  abs(cgccomplex(1))**2.0d0
                        cgcabs2=   abs(coeff(jg,bandcount,isp))**2.0d0
  
                      end if ! final state is a plane wave

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
                                   
!                            l2                = (1/(kxbroad*kybroad*two_pi))* & 
!                                                exp(-(qplusGx - kxarr)*(qplusGx - kxarr)/(2*kxbroad*kxbroad) ) &
!                                               * exp(-(qplusGy - kyarr)*(qplusGy - kyarr)/(2*kybroad*kybroad) )    
  
                            photoint_rank(ikx,iky) = photoint_rank(ikx,iky) + fermifac*l1*l2*cgcabs
                            if (.not. pureplanewave) then
                              photoint2_rank(ikx,iky)= photoint2_rank(ikx,iky)+ fermifac*l1*l2*cgcabs2
                            end if                                        
                          end if                 
                        end do  ! end loop over ky values       
                      end do  ! end loop over kx values  
                    end if !  ikx iky cond.
                  end if ! restriction to positive kz-values
                end do ! loop to symmetrize +/- q-values
              end if ! qplusGlen
            end do ! loop over G-vectors
          end if ! kinetic energy in desired range
!        photoint = photoint**2.0d0
        elseif (arupstype .eq. 2) then   ! 2  cis scan

          if (ebind.gt.ekinmin .and. ebind.lt.ekinmax) then  ! kinetic energy in desired range

            do ieph = 1+rank, nephk, numcpu  ! distribute loop over photon-energies to numcpu processes

              ephoton = ephotonmin + (ieph-1) * ekstep
              efinal   = ephoton - (phiwork - ebind)
              ekarr  = eklo + ephoton - phiwork
              kfinal = ekintok*sqrt(ekarr)       ! final wave vector in (1/A)
              deltae = efinal - ekarr
              l1     = (1/(ebroad*sqrt(two_pi)))*exp(-deltae*deltae/(2*ebroad*ebroad))      ! gauss broadening

              if (rank .eq. 0) then
                write(80,'(a3,i3,a3,f8.4,a4,i4,a4,f9.4,a6,f9.4,a10,f9.4)') &
                   'iq=',iwk,' ei=',ebind,' ib=',bandcount,' Eph=',ephoton,&
                   ' Ekin=',efinal,' Ekin_set=',ekarr
              end if

              do  jg = 1, nplane
              
                do j = 1,3
                  qplusG(j) =  igall(j,jg)   
                end do 
                hmil = hkl(1,jg)
                kmil = hkl(2,jg)
                lmil = hkl(3,jg)         
                qplusGlen=sqrt(qplusG(1)**2+qplusG(2)**2+qplusG(3)**2)
                         
                if ( (qplusGlen.ge.(kfinal - iwindow*kxbroad)) .and. &
                     (qplusGlen.le.(kfinal + iwindow*kxbroad)) ) then
     
                  do isymmetry = 0, nsymmetry   ! loop to symmetrize +/- q-values
           
                    if (isymmetry .eq. 0) then
                      qplusG =   qplusG              
                    else
                      qplusG =   -qplusG                          
                    end if
     
                                           
                    qplusGy   = qplusG(1)*vecy(1) + &
                                qplusG(2)*vecy(2) + &
                                qplusG(3)*vecy(3)
                    qplusGz   = qplusG(1)*vecz(1) + &
                                qplusG(2)*vecz(2) + &
                                qplusG(3)*vecz(3)
         
                    if ((abs(qplusGy) .lt. iwindow*kybroad) .and. (qplusGz.gt.0.0d0)) then ! restriction of ky-values

                      qplusGx   = qplusG(1)*vecx(1) + &
                                  qplusG(2)*vecx(2) + &
                                  qplusG(3)*vecx(3)

                      ikxstart = max(1  ,1+floor((qplusGx - iwindow*kxbroad - kxlo)/kxstep))
                      ikxend   = min(nkx,1+floor((qplusGx + iwindow*kxbroad - kxlo)/kxstep))
 
                      if (ikxend .ge. ikxstart) then 
                
                        if (pureplanewave) then    ! final state is a plane wave
                
                          cgcabs = abs(coeff(jg,bandcount,isp))**2.0d0
                     
                        else     ! final state is a damped plane wave (exp(gamma*z) factor)
                
                          cgccomplex = (0.0d0,0.0d0)
                     
                          do  jgz = 1, nplane   ! loop over Gz-vectors
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
                     
                              cgccomplex(1) = cgccomplex(1) &
                                              + fgz*conjg(coeff(jgz,bandcount,isp))
         
                            end if
                          end do                               ! end loop over Gz-vectors
                
                          cgcabs = abs(cgccomplex(1))**2.0d0
                          cgcabs2= abs(coeff(jg,bandcount,isp))**2.0d0
                          
                        end if   ! decision on type of final state

                        do ikx = ikxstart, ikxend  ! loop over kx values  

                          kxarr    = kxlo + (ikx - 1)*kxstep
                          kperp    = dsqrt(kfinal*kfinal - kxarr*kxarr)

                          if (kxarr .le. kfinal .and. kxarr .ge. -kfinal) then

!                            deltak   = sqrt((qplusGx - kxarr)*(qplusGx - kxarr) +&
!                                            (qplusGy - 0)    *(qplusGy - 0))
                            deltak   = dsqrt((qplusGx - kxarr)*(qplusGx - kxarr) +&
                                            (qplusGy - 0)    *(qplusGy - 0) +&
                                            (qplusGz - kperp)*(qplusGz - kperp))
                                                                                                                                     
                            l2 = (1/(kxbroad*sqrt(two_pi)))*exp(-deltak*deltak/(2*kxbroad*kxbroad))
                            photoint_rank(ikx,ieph) = photoint_rank(ikx,ieph) + fermifac*l1*l2*cgcabs
                            if (.not. pureplanewave) then
                              photoint2_rank(ikx,ieph)= photoint2_rank(ikx,ieph)+ fermifac*l1*l2*cgcabs2
                            end if                                       
                          end if                       
                        end do  ! end loop over kx values  
                      end if  ! ikxend .ge. ikxstart
                    end if  ! restriction of ky-values
                  end do  ! loop to symmetrize +/- q-values
                end if          
              end do  ! end loop over G-vectors         
            end do !end loop over photon energies
          end if ! kinetic energy in desired range

        end if ! decision 0/1/2 E-vs k / kx-ky plot / cis
      end do ! loop over bands
   end do! loop over q points in BZ

   if (arupstype .eq. 0) then
        call MPI_REDUCE(photoint_rank,photoint,nkx*nek,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (.not. pureplanewave) then
            call MPI_REDUCE(photoint2_rank,photoint2,nkx*nek,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        end if
   elseif (arupstype .eq. 1) then
        call MPI_REDUCE(photoint_rank,photoint,nkx*nkx,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (.not. pureplanewave) then
            call MPI_REDUCE(photoint2_rank,photoint2,nkx*nkx,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        end if
   else
        call MPI_REDUCE(photoint_rank,photoint,nkx*nephk,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (.not. pureplanewave) then
            call MPI_REDUCE(photoint2_rank,photoint2,nkx*nephk,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        end if
   end if

   if (rank .eq. 0) then
       photointspin = photointspin + photoint
       if (.not. pureplanewave) then 
          photointspin2 = photointspin2 + photoint2
       end if
   end if
end do ! loop over spin

if (rank .eq. 0) then
   photointspin = photointspin/(nwk*nspin)
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (rank .eq. 0) then 

      write(82,'(I5,t35,A)') nwk,           '! number of k-points'
      write(82,'(2I5,t35,A)') bandlow,bandhi,'! first and last bandindex for summation over initial states'
      write(82,'(f10.5,t35,A)') efermi,      '! Fermi-level (energy of HOMO) in [eV]'
      write(82,'(f10.5,t35,A)') temperature, '! Temperature in [eV]'
      write(82,'(f10.5,t35,A)') phiwork,     '! work function in [eV]'
      if (arupstype .eq.1) then
        write(82,'(2f10.5,t35,A)') ephotonmin, eklo, '! energy of incident photon [eV], binding energy [eV]'
      else
        write(82,'(f10.5,t35,A)') ephotonmin, '! energy of incident photon [eV], '
      end if
      write(82,'(f10.5,t35,A)') ebroad,      '! energy broadening [eV]'
      write(82,'(2f10.5,t35,A)') kxbroad,kybroad,'! broadening in k-space [A^-1]'
      write(82,'(i10,t35,A)') arupstype,     '! 0 = Energy versus k, 1 = kx-ky plot at constant energy, 2 = CIS scan'
      write(82,'(3f10.5,t35,A)') vecz,       '! unit vector for normal emission (z-axis)'
      write(82,'(3f10.5,t35,A)') vecx,       '! unit vector for grazing emission (x-axis)' 

      if (arupstype .eq.0) then
          write(82,'(2I5)')    nek,nkx
          write(82,'(3F15.9)') eklo,ekhi,ekstep 
          write(82,'(3F15.9)') kxlo,kxhi,kxstep 
      elseif (arupstype .eq.1) then
          write(82,'(2I5)')    nkx,nkx
          write(82,'(3F15.9)') kxlo,kxhi,kxstep       
          write(82,'(3F15.9)') kxlo,kxhi,kxstep    
      else
          write(82,'(2I5)')    nephk,nkx
          write(82,'(3F15.9)') ephotonmin,ephotonmax,ekstep       
          write(82,'(3F15.9)') kxlo,kxhi,kxstep  
      end if
          
      write(82,'(6f16.10)') photointspin


      if (.not. pureplanewave) then
         write(83,'(I5,t35,A)') nwk,            '! number of k-points'
         write(83,'(2I5,t35,A)') bandlow,bandhi,'! first and last bandindex for summation over initial states'
         write(83,'(f10.5,t35,A)') efermi,      '! Fermi-level (energy of HOMO) in [eV]'
         write(83,'(f10.5,t35,A)') temperature, '! Temperature in [eV]'
         write(83,'(f10.5,t35,A)') phiwork,     '! work function in [eV]'
         if (arupstype .eq.1) then
           write(83,'(2f10.5,t35,A)') ephotonmin, eklo, '! energy of incident photon [eV], binding energy [eV]'
         else
           write(83,'(f10.5,t35,A)') ephotonmin, '! energy of incident photon [eV], '
         end if
         write(83,'(f10.5,t35,A)') ebroad,      '! energy broadening [eV]'
         write(83,'(2f10.5,t35,A)') kxbroad,kybroad,'! broadening in k-space [A^-1]'
         write(83,'(i10,t35,A)') arupstype,     '! 0 = Energy versus k, 1 = kx-ky plot at constant energy, 2 = CIS scan'
         write(83,'(3f10.5,t35,A)') vecz,       '! unit vector for normal emission (z-axis)'
         write(83,'(3f10.5,t35,A)') vecx,       '! unit vector for grazing emission (x-axis)' 
  
         if (arupstype .eq.0) then
             write(83,'(2I5)')    nek,nkx
             write(83,'(3F15.9)') eklo,ekhi,ekstep 
             write(83,'(3F15.9)') kxlo,kxhi,kxstep 
         elseif (arupstype .eq.1) then
             write(83,'(2I5)')    nkx,nkx
             write(83,'(3F15.9)') kxlo,kxhi,kxstep       
             write(83,'(3F15.9)') kxlo,kxhi,kxstep 
         else
             write(83,'(2I5)')    nek,nkx
             write(83,'(3F15.9)') ephotonmin,ephotonmax,ekstep       
             write(83,'(3F15.9)') kxlo,kxhi,kxstep        
         end if
         
         write(83,'(6f16.10)') photointspin2
         
         close(83)
      
      end if

      endtime = wtime ( ); ! get wall time 
      WRITE(80,*)
      WRITE(80,*) 'ELAPSED WALL TIME = ',endtime - starttime,' seconds'
       

      WRITE(80,*)
      WRITE(80,*) 'DONE WITH ARUPS'

end if ! rank = 0

call mpi_finalize(ierr)

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

