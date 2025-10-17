module wave_data_pw

    use comms

    public ::  read_wavecar

contains



subroutine read_wavecar(KKK)

    implicit real*8 (a-h, o-z)

    integer,intent(in)::KKK

    complex*8 , allocatable :: coeff(:)
    complex*16, allocatable :: cener(:)
    real*8    , allocatable :: occ  (:)

    real*8   :: omeg, n_x,n_y,n_z

    dimension  a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),a2xa3(3)
    dimension  sumkg(3),vtmp(3)

    character*75  filename

    PARAMETER (PI=3.141592653589793d0)      


!-----read WAVECAR HEAD------------
    filename="WAVECAR"
    data c/0.26246582250210965422d0/ 
    data TOLK/1E-4/ 
    data TOLPH/1E-10/ 


!----- set nrecl ------------------
    nrecl=24

!-END-read WAVECAR HEAD------------
    open(unit=10,file=filename,access='direct',recl=nrecl, &
         iostat=iost,status='old')
    if (iost.ne.0) write(6,*) 'open error - iostat =',iost            

    read(unit=10,rec=1) xnrecl,xispin,xnprec
    nrecl=nint(xnrecl)
    ispin=nint(xispin)
    nprec=nint(xnprec)
    if(nprec.eq.45210) stop  '*** error - WAVECAR_double requires complex*16'
    if(ispin.eq.2) then
       !write(0,*) '*** NOTE - COMPLETE for FM  case, ISPIN =',ispin
    endif
    close(unit=10)
!-END-read WAVECAR HEAD------------

!-----open FILE--------
    open(unit=10,file=filename,access='direct',recl=nrecl, &
         iostat=iost,status='old')
    if (iost.ne.0) write(6,*) 'open error - iostat =',iost

    read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3), &
                                         (a2(j),j=1,3), &
                                         (a3(j),j=1,3)
    nwk=nint(xnwk)
    nband=nint(xnband)
    if ( NSPIN .ne. ispin) then
        write(0,*)   '*** error - selected k=',NSPIN,' > max k=',ispin
        stop
    endif
    if ( num_k .ne. nwk) then
        write(0,*)   '*** error - selected k=',num_k,' > max k=',nwk
        stop
    endif
    if ( num_bands .ne. nband ) then
        write(0,*)   '*** error - selected band=',num_bands,' > max band=',nband
        stop
    endif
      


    allocate(cener(nband));cener=0.d0
    allocate(occ  (nband));  occ=0.d0

!---------tmp for 624
    if(num_k<10.and.kkk==1) then
        write(624,'(I3)') NSPIN*num_k
        do iK=1,NSPIN*num_k
            irec=3+(iK-1)*(nband+1)
            read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
                                   (cener(iband),occ(iband),iband=1,nband)
            write(624,'(3F12.6)') (wk(i),i=1,3)
         enddo
    endif
!---------tmp for 624


!------Find desired wavefunction----------
    irec=3+(KKK-1)*(nband+1)
    read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
                           (cener(iband),occ(iband),iband=1,nband)
    nplane=nint(xnplane)
    EE=0.d0;EE(1:nband)=real(cener(1:nband))

!-----
    do i=1,3
        IF(abs(3.d0*wk(i)-1.d0) .lt. 0.003d0)  wk(i)= 1.d0/3.d0
        IF(abs(3.d0*wk(i)+1.d0) .lt. 0.003d0)  wk(i)=-1.d0/3.d0

        IF(abs(3.d0*wk(i)-2.d0) .lt. 0.003d0)  wk(i)= 2.d0/3.d0
        IF(abs(3.d0*wk(i)+2.d0) .lt. 0.003d0)  wk(i)=-2.d0/3.d0
    enddo

!.....output
    !WRITE(6,500)
    !WRITE(6,509) KKK
    !WRITE(6,510) (WK(I),I=1,3)

 500  FORMAT(/,80('*'))
 509  FORMAT(/,/,'knum =',I3,4X,'kname= ',A10)
 510  FORMAT('k =',3F9.6,/)


      
!----------------  FROM ECUT  -----------------
!-----compute reciprocal properties------------
    call vcross(a2xa3,a2,a3)
    Vcell=dot_product(a1,a2xa3)
    a3mag=dsqrt(dot_product(a3,a3))
    call vcross(b1,a2,a3)
    call vcross(b2,a3,a1)
    call vcross(b3,a1,a2)
    b1=2.*pi*b1/Vcell
    b2=2.*pi*b2/Vcell
    b3=2.*pi*b3/Vcell
    b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
    b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
    b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)

      
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

    npmax=2*(1+2*nb1max)*(1+2*nb2max)*(1+2*nb3max)

    igall=0
    KV=0.d0
 
    igall = 0
!-----Calculate plane waves---------------------
!-----FOR a special K point -------------------
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
             gtot=sqrt(dot_product(sumkg,sumkg))
             etot=gtot**2/c
             if (etot.lt.ecut) then
                ncnt=ncnt+1
                igall(1,ncnt)=ig1p
                igall(2,ncnt)=ig2p
                igall(3,ncnt)=ig3p
                KV(:,ncnt)=DBLE(igall(:,ncnt))+WK(:)
             end if
          enddo
       enddo
    enddo

!--------judge SOC----------
      !! 2* to handle two component spinors
    if (2*ncnt.eq.nplane) then
        !write(6,*) 'Spin-orbit Wavefunctions(INCOMPLETE): just for parity'
    elseif ( ncnt .eq. nplane) then
        !write(6,*) 'No spin-orbit Wavefunctions'
    else
        write(0,*) '*** error - computed 2*ncnt=',2*ncnt, &
                   ' != input nplane=',nplane
        stop
    endif
!---------------- ----------
!
    allocate (coeff(2*ncnt )); coeff=0.0
    coeffa=0.d0
    coeffb=0.d0
    L=0
    DO iband=1,nband
        !------Read desired wavefunction----------
        coeff=0.0
        irec=3+( KKK -1)*(nband+1)+iband
        read(unit=10,rec=irec) (coeff(iplane), iplane=1,nplane)

        !WRITE(547,500)
        !WRITE(547,507) KKK,iband
        !WRITE(547,510) (WK(I),I=1,3)
        !do iplane=1,nplane
        !    write(547,511) iplane,igall(:,iplane),coeff(iplane),abs(coeff(iplane))!coeff(iplane+ncnt)
        !enddo

        coeffa(1:ncnt,iband)=coeff(1:ncnt)
        coeffb(1:ncnt,iband)=coeff(ncnt+1:2*ncnt)
    ENDDO

 507  FORMAT(/,/,'knum =',I3,4X,'nband= ',I3)
 511  FORMAT(4I5,2X,' (',2F10.6,' )',2X,2F12.6)

    deallocate(coeff)
    deallocate(cener)
    deallocate(occ)
    close(10)

    return
end subroutine read_wavecar 


subroutine bubble_sort_with_index(nband, A, idx)
    integer, intent(in) :: nband
      real(dp)   , intent(inout) :: A(nband)
      integer, intent(out) :: idx(nband)

      integer :: n, i, j, tmp_i
      real    :: tmp_r

      n = size(A)
      idx = [(i, i=1, n)]  

      do i = 1, n-1
         do j = 1, n-i
            if (A(idx(j)) > A(idx(j+1))) then
               tmp_i      = idx(j)
               idx(j)     = idx(j+1)
               idx(j+1)   = tmp_i
            end if
         end do
      end do
      A = A(idx)
   end subroutine bubble_sort_with_index


subroutine read_wavecar_spin_polarized(kkk)
    integer,intent(in)::KKK
    EE = 0.d0
    coeffa=0.d0
    coeffb=0.d0
    call read_wavecar_up(KKK)
    call read_wavecar_dn(KKK)
    call bubble_sort_with_index(num_bands*2,EE,order)
    coeffa(:,:) = coeffa(:,order)
    coeffb(:,:) = coeffb(:,order)

endsubroutine

subroutine read_wavecar_up(KKK)

    implicit real*8 (a-h, o-z)

    integer,intent(in)::KKK

    complex*8 , allocatable :: coeff(:)
    complex*16, allocatable :: cener(:)
    real*8    , allocatable :: occ  (:)

    real*8   :: omeg, n_x,n_y,n_z

    dimension  a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),a2xa3(3)
    dimension  sumkg(3),vtmp(3)

    character*75  filename

    PARAMETER (PI=3.141592653589793d0)      


!-----read WAVECAR HEAD------------
    filename="WAVECAR"
    data c/0.26246582250210965422d0/ 
    data TOLK/1E-4/ 
    data TOLPH/1E-10/ 


!----- set nrecl ------------------
    nrecl=24

!-END-read WAVECAR HEAD------------
    open(unit=10,file=filename,access='direct',recl=nrecl, &
         iostat=iost,status='old')
    if (iost.ne.0) write(6,*) 'open error - iostat =',iost            

    read(unit=10,rec=1) xnrecl,xispin,xnprec
    nrecl=nint(xnrecl)
    ispin=nint(xispin)
    nprec=nint(xnprec)
    if(nprec.eq.45210) stop  '*** error - WAVECAR_double requires complex*16'
    if(ispin.eq.2) then
       !write(0,*) '*** NOTE - COMPLETE for FM  case, ISPIN =',ispin
    endif
    close(unit=10)
!-END-read WAVECAR HEAD------------

!-----open FILE--------
    open(unit=10,file=filename,access='direct',recl=nrecl, &
         iostat=iost,status='old')
    if (iost.ne.0) write(6,*) 'open error - iostat =',iost

    read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3), &
                                         (a2(j),j=1,3), &
                                         (a3(j),j=1,3)
    nwk=nint(xnwk)
    nband=nint(xnband)
    if ( NSPIN .ne. ispin) then
        write(0,*)   '*** error - selected k=',NSPIN,' > max k=',ispin
        stop
    endif
    if ( num_k .ne. nwk) then
        write(0,*)   '*** error - selected k=',num_k,' > max k=',nwk
        stop
    endif
    if ( num_bands .ne. nband ) then
        write(0,*)   '*** error - selected band=',num_bands,' > max band=',nband
        stop
    endif
      


    allocate(cener(nband));cener=0.d0
    allocate(occ  (nband));  occ=0.d0

!---------tmp for 624
    if(num_k<10.and.kkk==1) then
        write(624,'(I3)') NSPIN*num_k
        do iK=1,NSPIN*num_k
            irec=3+(iK-1)*(nband+1)
            read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
                                   (cener(iband),occ(iband),iband=1,nband)
            write(624,'(3F12.6)') (wk(i),i=1,3)
         enddo
    endif
!---------tmp for 624


!------Find desired wavefunction----------
    irec=3+(KKK-1)*(nband+1)
    read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
                           (cener(iband),occ(iband),iband=1,nband)
    nplane=nint(xnplane)
    EE(1:nband)=real(cener(1:nband))

!-----
    do i=1,3
        IF(abs(3.d0*wk(i)-1.d0) .lt. 0.003d0)  wk(i)= 1.d0/3.d0
        IF(abs(3.d0*wk(i)+1.d0) .lt. 0.003d0)  wk(i)=-1.d0/3.d0

        IF(abs(3.d0*wk(i)-2.d0) .lt. 0.003d0)  wk(i)= 2.d0/3.d0
        IF(abs(3.d0*wk(i)+2.d0) .lt. 0.003d0)  wk(i)=-2.d0/3.d0
    enddo
      
!----------------  FROM ECUT  -----------------
!-----compute reciprocal properties------------
    call vcross(a2xa3,a2,a3)
    Vcell=dot_product(a1,a2xa3)
    a3mag=dsqrt(dot_product(a3,a3))
    call vcross(b1,a2,a3)
    call vcross(b2,a3,a1)
    call vcross(b3,a1,a2)
    b1=2.*pi*b1/Vcell
    b2=2.*pi*b2/Vcell
    b3=2.*pi*b3/Vcell
    b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
    b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
    b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)

      
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

    npmax=2*(1+2*nb1max)*(1+2*nb2max)*(1+2*nb3max)

    igall=0
    KV=0.d0
 
    igall = 0
!-----Calculate plane waves---------------------
!-----FOR a special K point -------------------
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
             gtot=sqrt(dot_product(sumkg,sumkg))
             etot=gtot**2/c
             if (etot.lt.ecut) then
                ncnt=ncnt+1
                igall(1,ncnt)=ig1p
                igall(2,ncnt)=ig2p
                igall(3,ncnt)=ig3p
                KV(:,ncnt)=DBLE(igall(:,ncnt))+WK(:)
             end if
          enddo
       enddo
    enddo

    allocate (coeff(2*ncnt )); coeff=0.0
    L=0
    DO iband=1,nband
        !------Read desired wavefunction----------
        coeff=0.0
        irec=3+( KKK -1)*(nband+1)+iband
        read(unit=10,rec=irec) (coeff(iplane), iplane=1,nplane)

        coeffa(1:ncnt,iband)=coeff(1:ncnt)
        coeffb(1:ncnt,iband)=coeff(ncnt+1:2*ncnt)
    ENDDO

    deallocate(coeff)
    deallocate(cener)
    deallocate(occ)
    close(10)

    return
end subroutine read_wavecar_up



subroutine read_wavecar_dn(KKK_)

    implicit real*8 (a-h, o-z)

    integer,intent(in)::KKK_

    complex*8 , allocatable :: coeff(:)
    complex*16, allocatable :: cener(:)
    real*8    , allocatable :: occ  (:)

    real*8   :: omeg, n_x,n_y,n_z

    dimension  a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),a2xa3(3)
    dimension  sumkg(3),vtmp(3)

    character*75  filename

    PARAMETER (PI=3.141592653589793d0)     
    integer :: KKK

    KKK = KKK_ + num_k 


!-----read WAVECAR HEAD------------
    filename="WAVECAR"
    data c/0.26246582250210965422d0/ 
    data TOLK/1E-4/ 
    data TOLPH/1E-10/ 


!----- set nrecl ------------------
    nrecl=24

!-END-read WAVECAR HEAD------------
    open(unit=10,file=filename,access='direct',recl=nrecl, &
         iostat=iost,status='old')
    if (iost.ne.0) write(6,*) 'open error - iostat =',iost            

    read(unit=10,rec=1) xnrecl,xispin,xnprec
    nrecl=nint(xnrecl)
    ispin=nint(xispin)
    nprec=nint(xnprec)
    if(nprec.eq.45210) stop  '*** error - WAVECAR_double requires complex*16'
    if(ispin.eq.2) then
       !write(0,*) '*** NOTE - COMPLETE for FM  case, ISPIN =',ispin
    endif
    close(unit=10)
!-END-read WAVECAR HEAD------------

!-----open FILE--------
    open(unit=10,file=filename,access='direct',recl=nrecl, &
         iostat=iost,status='old')
    if (iost.ne.0) write(6,*) 'open error - iostat =',iost

    read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3), &
                                         (a2(j),j=1,3), &
                                         (a3(j),j=1,3)
    nwk=nint(xnwk)
    nband=nint(xnband)
    if ( NSPIN .ne. ispin) then
        write(0,*)   '*** error - selected k=',NSPIN,' > max k=',ispin
        stop
    endif
    if ( num_k .ne. nwk) then
        write(0,*)   '*** error - selected k=',num_k,' > max k=',nwk
        stop
    endif
    if ( num_bands .ne. nband ) then
        write(0,*)   '*** error - selected band=',num_bands,' > max band=',nband
        stop
    endif
      


    allocate(cener(nband));cener=0.d0
    allocate(occ  (nband));  occ=0.d0

!---------tmp for 624
    if(num_k<10.and.kkk==1) then
        write(624,'(I3)') NSPIN*num_k
        do iK=1,NSPIN*num_k
            irec=3+(iK-1)*(nband+1)
            read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
                                   (cener(iband),occ(iband),iband=1,nband)
            write(624,'(3F12.6)') (wk(i),i=1,3)
         enddo
    endif
!---------tmp for 624


!------Find desired wavefunction----------
    irec=3+(KKK-1)*(nband+1)
    read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
                           (cener(iband),occ(iband),iband=1,nband)
    nplane=nint(xnplane)
    EE(nband+1:2*nband)=real(cener(1:nband))

!-----
    do i=1,3
        IF(abs(3.d0*wk(i)-1.d0) .lt. 0.003d0)  wk(i)= 1.d0/3.d0
        IF(abs(3.d0*wk(i)+1.d0) .lt. 0.003d0)  wk(i)=-1.d0/3.d0

        IF(abs(3.d0*wk(i)-2.d0) .lt. 0.003d0)  wk(i)= 2.d0/3.d0
        IF(abs(3.d0*wk(i)+2.d0) .lt. 0.003d0)  wk(i)=-2.d0/3.d0
    enddo
      
!----------------  FROM ECUT  -----------------
!-----compute reciprocal properties------------
    call vcross(a2xa3,a2,a3)
    Vcell=dot_product(a1,a2xa3)
    a3mag=dsqrt(dot_product(a3,a3))
    call vcross(b1,a2,a3)
    call vcross(b2,a3,a1)
    call vcross(b3,a1,a2)
    b1=2.*pi*b1/Vcell
    b2=2.*pi*b2/Vcell
    b3=2.*pi*b3/Vcell
    b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
    b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
    b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)

      
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

    npmax=2*(1+2*nb1max)*(1+2*nb2max)*(1+2*nb3max)

    igall=0
    KV=0.d0
 
    igall = 0
!-----Calculate plane waves---------------------
!-----FOR a special K point -------------------
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
             gtot=sqrt(dot_product(sumkg,sumkg))
             etot=gtot**2/c
             if (etot.lt.ecut) then
                ncnt=ncnt+1
                igall(1,ncnt)=ig1p
                igall(2,ncnt)=ig2p
                igall(3,ncnt)=ig3p
                KV(:,ncnt)=DBLE(igall(:,ncnt))+WK(:)
             end if
          enddo
       enddo
    enddo

    allocate (coeff(2*ncnt )); coeff=0.0
    L=0
    DO iband=1,nband
        !------Read desired wavefunction----------
        coeff=0.0
        irec=3+( KKK -1)*(nband+1)+iband
        read(unit=10,rec=irec) (coeff(iplane), iplane=1,nplane)

        coeffb(1:ncnt,iband+nband)=coeff(1:ncnt)
        coeffa(1:ncnt,iband+nband)=coeff(ncnt+1:2*ncnt)
    ENDDO

    deallocate(coeff)
    deallocate(cener)
    deallocate(occ)
    close(10)

    return
end subroutine read_wavecar_dn


!!$*   routine for computing vector cross-product
subroutine vcross(a,b,c)

    implicit real*8(a-h,o-z)
  
    dimension a(3),b(3),c(3)
  
    a(1)=b(2)*c(3)-b(3)*c(2)
    a(2)=b(3)*c(1)-b(1)*c(3)
    a(3)=b(1)*c(2)-b(2)*c(1)
    return
end subroutine vcross      

end module wave_data_pw 
