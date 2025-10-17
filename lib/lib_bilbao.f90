module bilbao 

    use lib_comms
    implicit none 

    real(dp),   public,  save :: Kc2p(3,3), p2cR(3,3)

contains


subroutine bilbao_read(sgn)

    !!! This subroutine read symmetry operations and 
    !!! character tables from Bilbao.
    !!! Compare the input symmetry operations with Bilbao. 

    integer,          intent(in) :: sgn

    character(len=3)   :: csgn
    character(len=180) :: spgpath, spgfile
    character(len=10)  :: symbol_sg

    real(dp)           :: Df(2,4), abcde(5), ktmp(3), ttmp(1,3)
    real(dp)           :: tmp33(3,3) 

    integer            :: Numk, antiss, tnir, wi 

    integer            :: i, j, itmp, jtmp 
    integer            :: iir, nele, ikt 
    character(len=5)   :: irtmp 
    character(len=2)   :: nametmp, nametmp2 
    character(len=15)  :: ckpoint 
    character(len=40)  :: ListIrrep

    ! bilbao table file
    integer,     parameter   :: bb = 11
    ! output file : operations in conventional basis
    integer,     parameter   :: op = 9
    !! output file : special k points
    !integer,     parameter   :: sk = 8

    if     (sgn < 10)  then; write(csgn,'(I1)') sgn
    elseif (sgn < 100) then; write(csgn,'(I2)') sgn
    else                   ; write(csgn,'(I3)') sgn
    endif 
! Prefer a compile-time absolute path if provided
#ifdef IRSSGDATA_PATH
    spgpath = IRSSGDATA_PATH
#elif defined(IRVSPDATA)
    call get_environment_variable('IRVSPDATA',spgpath)
#else
    write(6,*) "Environment variable 'IRVSPDATA' must be provided "
    write(6,*) "Please run the following commands to make the library:"
    write(6,*) "./configure.sh"
    write(6,*) "source ~/.bashrc"
    write(6,*) "make lib"
    stop
#endif 
    spgfile = trim(spgpath)//'/kLittleGroups/kLG_'//trim(csgn)//'.data'
    write(*,*) "SPGFILE :", trim(adjustl(spgfile)) 

    open(unit=bb, file=spgfile, status='old', form='unformatted')
    !open(unit=sk, file='SGklist_'//trim(csgn)//'.cht', status='unknown')

    read(bb) num_doub_sym, symbol_sg

    if     (symbol_sg(1:1) == 'P') then
        Kc2p = Pabc
    elseif (symbol_sg(1:1) == 'C') then 
        Kc2p = Cabc
        if (sgn == 68) Kc2p = Cabc68    ! added by szhang on 17/4/2024
    elseif (symbol_sg(1:1) == 'B') then
        Kc2p = Babc
    elseif (symbol_sg(1:1) == 'A') then
        Kc2p = Aabc
    elseif (symbol_sg(1:1) == 'R') then
        Kc2p = Rabc
    elseif (symbol_sg(1:1) == 'F') then
        Kc2p = Fabc
    elseif (symbol_sg(1:1) == 'I') then
        Kc2p = Iabc
    else
        stop "Error in space-group-symbol"
    endif 
    call invreal33(Kc2p, p2cR)

    rot_bilbao = 0
    tau_bilbao = 0.d0
    SU2_bilbao = 0.d0
    do i = 1, num_doub_sym 
        read(bb) rot_bilbao(:,:,i), tau_bilbao(:,i), Df(:,:)
        SU2_bilbao(1,1,i)=cmplx(Df(1,1)*dcos(PI*Df(2,1)),Df(1,1)*dsin(PI*Df(2,1)),dp)
        SU2_bilbao(1,2,i)=cmplx(Df(1,2)*dcos(PI*Df(2,2)),Df(1,2)*dsin(PI*Df(2,2)),dp)
        SU2_bilbao(2,1,i)=cmplx(Df(1,3)*dcos(PI*Df(2,3)),Df(1,3)*dsin(PI*Df(2,3)),dp)
        SU2_bilbao(2,2,i)=cmplx(Df(1,4)*dcos(PI*Df(2,4)),Df(1,4)*dsin(PI*Df(2,4)),dp)
    enddo 

    ! conventional cell -> primitive cell
    invrot_bilbao = 0
    do i = 1, num_doub_sym 
        tmp33(:,:) = dble(transpose(rot_bilbao(:,:,i)))
        rot_bilbao(:,:,i) = nint(matmul(matmul(p2cR, tmp33), Kc2p(:,:))) 
        tau_bilbao(:,i) = matmul(p2cR, tau_bilbao(:,i))
        call invmati(rot_bilbao(:,:,i), invrot_bilbao(:,:,i))
    enddo 


    ! read character tables
    nirreps(:)=0
    sirreps(:)=0
    antisym(:)=0
    Herringrule(:,:)=-2

    labels=0
    iir=0;ikt=0
    tableTraces(:,:,:)=cmplx_0 
    chartTraces(:,:,:)=cmplx_0 
    coeff_uvw(:,:,:,:)=cmplx_0
    factsTraces(:,:,:)='            '
    chkpoint(:)='***************'
    nametmp='  '
    read(11)  Numk,tnir

    DO wi=1,tnir
      read(11) ListIrrep
      nametmp2=nametmp
      read(ListIrrep,*) ktmp(:),antiss,irtmp, itmp,itmp,nametmp,jtmp
      IF(nametmp/=nametmp2) THEN 
         IF(ikt>0) THEN
            nirreps(ikt)=iir
           !call dumptableofIrs(ikt,9)
            call Kreal2string(samplek(:,ikt),ckpoint) 
            ttmp(1,1:3)=samplek(1:3,ikt)
         ENDIF
        !------out
         ikt=ikt+1; iir=0
         antisym(ikt)=antiss
      ENDIF
!
      samplek(:,ikt)=ktmp(:)
      samplekname(ikt)=nametmp
      iir=iir+1
      nele=0
      irk:DO j=1,num_doub_sym
          read(11) itmp,itmp
          IF(itmp==1) THEN
             nele=nele+1
             labels(1,j,iir,ikt)=itmp
             read(11) itmp;labels(2,j,iir,ikt)=itmp
             abcde(:)=0._dp
             IF(itmp==1) THEN
                read(11) abcde(1:2)
             ELSEIF(itmp==2) THEN
                read(11) abcde(1:5)
                coeff_uvw(:,j,iir,ikt)=abcde(3:5)
                write(factsTraces(j,iir,ikt),"(3F4.1)") abcde(3:5)
             ELSE
                STOP "Error!" 
             ENDIF
                chartTraces(j,iir,ikt)=cmplx(abcde(1)*dcos(PI*abcde(2)) &
                                            ,abcde(1)*dsin(PI*abcde(2)),dp)
          ELSEIF(itmp==0) THEN
          ELSE
           STOP "Error!!" 
          ENDIF
      END DO irk
!
      nelelittle (    ikt) = nele
      Irrepsname (iir,ikt) = irtmp
      Herringrule(iir,ikt) = jtmp
    END DO

            nirreps(ikt)=iir

            call Kreal2string(samplek(:,ikt),ckpoint) 
            ttmp(1,1:3)=samplek(1:3,ikt)

    IF(ikt/=Numk) STOP"ERROR in little groups of k-points"
    num_ktype=ikt



end subroutine bilbao_read 


subroutine bilbao_getkid(k_input, &
                         kname,k_frac_symbol)

    real(dp)        , intent(inout) :: k_input(3)
    character(len=15),intent(out)   :: k_frac_symbol
    integer            :: ind_ktype
    character(len=3), intent(out)   :: kname
    integer         :: ind_rot 
    logical         :: timerev_k

    integer     :: ind_ktype_tmp 
    integer     :: nkt(num_doub_sym), nvar(num_doub_sym)

    integer     :: iR,j,ivar,iw
    real(dp)    :: k_trial(3), k_input_conv(3)

    integer     :: iir
    real(dp)    :: AGW
    complex(dp) :: PHW

    timerev_k = .false.

    ind_ktype  = 0
    nkt(:)     = 0
    nvar(:)    = 0
    k_trial(:) = 0._dp
    do iR=1, num_doub_sym/2 + 1
        if (iR < num_doub_sym/2 + 1) then 
            k_trial(:) = matmul(k_input(:),invrot_bilbao(:,:,iR))
            call getkid2(k_trial,ind_ktype_tmp, ivar)
        else 
            k_trial(:) = -k_input(:)
            call getkid2(k_trial, ind_ktype_tmp, ivar)
        endif 
        nkt(iR)  = ind_ktype_tmp
        nvar(iR) = ivar
        if (ivar == 0) exit 
    enddo 

    if (iR == num_doub_sym/2 + 1) then 
        timerev_k = .true.
    endif 

    if (iR == num_doub_sym/2 + 2) then 
        zj:DO iw = 1,3 
            do iR = 1, num_doub_sym/2
                ivar = nvar(iR)
                if (ivar == iw) exit zj
            enddo
        ENDDO zj
    endif 
    
    ind_rot   = iR
    ind_ktype = nkt(iR)
    ivar      = nvar(iR)

    if (ind_ktype > num_ktype) stop " Nonsymmorphic kpoint is NOT found."
    kname = samplekname(ind_ktype)//' '
    call Kreal2string(samplek(:,ind_ktype),k_frac_symbol) 

    k_trial(1) = dot_product(k_input(:), invrot_bilbao(:,1,iR))
    k_trial(2) = dot_product(k_input(:), invrot_bilbao(:,2,iR))
    k_trial(3) = dot_product(k_input(:), invrot_bilbao(:,3,iR))
    if (timerev_k) k_trial(:) = -k_input(:)
    k_input(:) = k_trial(:)
     
    k_input_conv(:) = matmul(k_input(:), p2cR(:,:))
    do iir = 1,nirreps(ind_ktype)
        do j = 1,num_doub_sym

            if(labels(1,j,iir,ind_ktype) == 1) then
                
                if(labels(2,j,iir,ind_ktype) == 1) then
                    tableTraces(j,iir,ind_ktype)=chartTraces(j,iir,ind_ktype)
                elseif(labels(2,j,iir,ind_ktype) == 2) then
                    AGW=0._DP
                    AGW=PI*DOT_PRODUCT(coeff_uvw(1:3,j,iir,ind_ktype),k_input_conv(:))
                    PHW=CMPLX(DCOS(AGW),DSIN(AGW),DP)
                    tableTraces(j,iir,ind_ktype)=chartTraces(j,iir,ind_ktype)*PHW
                else
                    stop "Error" 
                endif
            endif

        enddo
    enddo

end subroutine bilbao_getkid 


      subroutine getkid2(k_input,ik,ivar)
      real(dp)   , intent(in ) :: k_input(3)
      integer    , intent(out) :: ik
      integer    , intent(out) :: ivar

      real(dp)    :: uvw(3)
      integer     :: ikt
      character*15 , allocatable  :: ckpoint(:)
      logical :: is_variable(3)
      integer :: j
      integer :: num_var_tmp, ref_varnum
      integer, allocatable :: num_var(:)
      integer, allocatable :: ind_u(:), ind_v(:), ind_w(:)
      
      real(dp) :: uvw_tmp(3)
      real(dp) :: k_input_conv(3), refkpoint(3), refkpointp(3)
      real(dp) :: diff, diff1, diff2, diff3
      real(dp), parameter :: tol=1e-5

      character*15 :: chr_tmp
      ivar=-1

      allocate(ckpoint(num_ktype))
      do ikt=1, num_ktype
         call Kreal2string(samplek(:,ikt),ckpoint(ikt)) 
      enddo

      

      allocate(num_var(num_ktype))
      num_var = 0
      allocate(ind_u(num_ktype))
      allocate(ind_v(num_ktype))
      allocate(ind_w(num_ktype))
      ind_u=0; ind_v=0; ind_w=0

     !print*,"================================================"
     !   write(*,'(20X,A5,10X,A6,20X,A5)') 'conv', 'symbol', 'prim'
     !do ikt=1,num_ktype
     !   write(*,'(I5,3F9.3,2X,A15,3F12.6)') ikt,samplek(:,ikt),ckpoint(ikt),matmul(samplek(:,ikt),Kc2p(:,:))
     !enddo
     !   write(*,*) 'Caution!!! conv  : 0.333 -> 0.3333333333'
     !   write(*,*) 'Caution!!! symbol: 0.33  -> 0.3333333333'
     ! !--- add something here

      !--- add something below-------
      ! get the number of variables for the reference kpoints
      do ikt = 1, num_ktype
         num_var_tmp = 0
         ind_u(ikt)   = index(ckpoint(ikt), "  u  ") 
         ind_v(ikt)   = index(ckpoint(ikt), "  v  ") 
         ind_w(ikt)   = index(ckpoint(ikt), "  w  ") 
         if (ind_u(ikt).ne.0) num_var_tmp = num_var_tmp+1
         if (ind_v(ikt).ne.0) num_var_tmp = num_var_tmp+1
         if (ind_w(ikt).ne.0) num_var_tmp = num_var_tmp+1
         num_var(ikt) = num_var_tmp 
      enddo 

      k_input_conv = matmul(k_input(:), p2cR(:,:))

      uvw = 0d0 
      ref_varnum = 9999
      do ikt = 1, num_ktype
         is_variable = .false.
         refkpoint = 9999d0
         uvw_tmp = 0d0 
         if(ckpoint(ikt)(1:5)=='  u  ') then 
            is_variable(1) = .true.
            uvw_tmp(1) = k_input_conv(1)
            refkpoint(1) = k_input_conv(1)
            if(ckpoint(ikt)(6:10)=='  u  ') refkpoint(2) = refkpoint(1)
            if(ckpoint(ikt)(6:10)==' 1-u ') refkpoint(2) = 1d0-refkpoint(1)
            if(ckpoint(ikt)(6:10)=='-2u  ') refkpoint(2) = -2d0*refkpoint(1)
            if(ckpoint(ikt)(6:10)==' -u  ') refkpoint(2) = -refkpoint(1)
            if(ckpoint(ikt)(6:10)==' 1+u ') refkpoint(2) = 1d0+refkpoint(1)
            if(ckpoint(ikt)(11:15)=='  u  ') refkpoint(3) = refkpoint(1)
            if(ckpoint(ikt)(11:15)==' 1-u ') refkpoint(3) = 1d0-refkpoint(1)
            if(ckpoint(ikt)(11:15)=='-2u  ') refkpoint(3) = -2d0*refkpoint(1)
            if(ckpoint(ikt)(11:15)==' -u  ') refkpoint(3) = -refkpoint(1)
            if(ckpoint(ikt)(11:15)==' 1+u ') refkpoint(3) = 1d0+refkpoint(1)
         endif 
         if(ckpoint(ikt)(6:10)=='  v  ') then 
            is_variable(2) = .true.
            uvw_tmp(2) = k_input_conv(2)
            refkpoint(2) = k_input_conv(2)
            if (ckpoint(ikt)(1:5)==' 1-v ') refkpoint(1) = 1d0-refkpoint(2)
            if (ckpoint(ikt)(11:15)==' 1-v ') refkpoint(3) = 1d0-refkpoint(2)
            if (ckpoint(ikt)(11:15)=='  v  ') refkpoint(3) = refkpoint(2)
         endif 
         if(ckpoint(ikt)(11:15)=='  w  ') then 
             is_variable(3) = .true.
             uvw_tmp(3) = k_input_conv(3)
             refkpoint(3) = k_input_conv(3)
         endif 

         ! some exceptions
         ! 0.5 u 0.0  : 195 198 200 201 205
         ! 1+u 1-u 0.0: 197 199 204 206 211 214 217 220 229 230
         if(ckpoint(ikt)(1:5).ne.'  u  '.and.ckpoint(ikt)(6:10).eq.'  u  ') then 
            is_variable(1) = .true.
            uvw_tmp(1) = k_input_conv(2)
            refkpoint(2) = k_input_conv(2)
         endif 
         if(ckpoint(ikt)(1:5).eq.' 1+u '.and.ckpoint(ikt)(6:10).eq.' 1-u ') then 
            is_variable(1) = .true.
            uvw_tmp(1) = k_input_conv(1) - 1
            refkpoint(1) = k_input_conv(1)
            refkpoint(2) = 2-refkpoint(1)
         endif 
         ! end exceptions
         
         do j = 1,3
            if (abs(refkpoint(j)-9999d0) < tol) then 
               chr_tmp = adjustr(ckpoint(ikt)((j-1)*5+1:j*5))
               read(chr_tmp, *) refkpoint(j)
               if (abs(refkpoint(j)-0.33) < tol) refkpoint(j) = 1d0/3d0
            endif 
         enddo 

         refkpointp = matmul(refkpoint,kc2p)

         diff = abs((refkpointp(1)-floor(refkpointp(1)))-(k_input(1)-floor(k_input(1)))) + &
                abs((refkpointp(2)-floor(refkpointp(2)))-(k_input(2)-floor(k_input(2)))) + &
                abs((refkpointp(3)-floor(refkpointp(3)))-(k_input(3)-floor(k_input(3))))    

         if (diff < tol) then 
             if (num_var(ikt) < ref_varnum) then 
                 ref_varnum = num_var(ikt)
                 ik = ikt 
                 uvw = 0d0
                 do j = 1, 3
                    if(is_variable(j)) uvw(j) = uvw_tmp(j)
                 enddo 
             endif 
         endif 

      enddo 
      ivar= num_var(ik)
     !write(*,*) "The id of reference kpoint for the input kpoint is:", ik,num_var(ik)
     !write(*,"(A10,3F12.6)") "uvw(:)=",uvw(:)
     !print*,"================================================"
      deallocate(ckpoint)
      deallocate(num_var)

      deallocate(ind_u)
      deallocate(ind_v)
      deallocate(ind_w)
      end subroutine getkid2


      SUBROUTINE Kreal2string(coorkp,kstr) 
      real(dp), intent(in)::  coorkp(3) ! Coordinates 
      character*15,intent(out):: kstr

!real(dp),    parameter  :: su=0.123_dp,s1_u=0.877_dp,s_2u=-0.246_dp
!real(dp),    parameter  :: sv=0.313_dp,s1_v=0.687_dp
!real(dp),    parameter  :: sw=0.427_dp,s_u=-0.123_dp,s1u = 1.123_dp
      ! u     ->  0.123  : su
      ! v     ->  0.313  : sv
      ! w     ->  0.427  : sw
      ! 1 - u ->  0.877  : s1_u
      ! -2 u  -> -0.246  : s_2u
      ! -u    -> -0.123  : s_u
      ! 1 + u ->  1.123  : s1u
      ! 1 - v ->  0.687  : s1_v
      !------------codes for u,v,w

      if     ( abs( coorkp(1)-su   )<epsil) then;kstr(1:5)="  u  "
      elseif ( abs( coorkp(1)-sv   )<epsil) then;kstr(1:5)="  v  "
      elseif ( abs( coorkp(1)-sw   )<epsil) then;kstr(1:5)="  w  "
      elseif ( abs( coorkp(1)-s1_u )<epsil) then;kstr(1:5)=" 1-u "
      elseif ( abs( coorkp(1)-s_2u )<epsil) then;kstr(1:5)="-2u  "
      elseif ( abs( coorkp(1)-s_u  )<epsil) then;kstr(1:5)=" -u  "
      elseif ( abs( coorkp(1)-s1u  )<epsil) then;kstr(1:5)=" 1+u "
      elseif ( abs( coorkp(1)-s1_v )<epsil) then;kstr(1:5)=" 1-v "
      else
       write(kstr(1:5),501) coorkp(1)
      endif

      if     ( abs( coorkp(2)-su   )<epsil) then;kstr(6:10)="  u  "
      elseif ( abs( coorkp(2)-sv   )<epsil) then;kstr(6:10)="  v  "
      elseif ( abs( coorkp(2)-sw   )<epsil) then;kstr(6:10)="  w  "
      elseif ( abs( coorkp(2)-s1_u )<epsil) then;kstr(6:10)=" 1-u "
      elseif ( abs( coorkp(2)-s_2u )<epsil) then;kstr(6:10)="-2u  "
      elseif ( abs( coorkp(2)-s_u  )<epsil) then;kstr(6:10)=" -u  "
      elseif ( abs( coorkp(2)-s1u  )<epsil) then;kstr(6:10)=" 1+u "
      elseif ( abs( coorkp(2)-s1_v )<epsil) then;kstr(6:10)=" 1-v "
      else
       write(kstr(6:10),501) coorkp(2)
      endif

      if     ( abs( coorkp(3)-su   )<epsil) then;kstr(11:15)="  u  "
      elseif ( abs( coorkp(3)-sv   )<epsil) then;kstr(11:15)="  v  "
      elseif ( abs( coorkp(3)-sw   )<epsil) then;kstr(11:15)="  w  "
      elseif ( abs( coorkp(3)-s1_u )<epsil) then;kstr(11:15)=" 1-u "
      elseif ( abs( coorkp(3)-s_2u )<epsil) then;kstr(11:15)="-2u  "
      elseif ( abs( coorkp(3)-s_u  )<epsil) then;kstr(11:15)=" -u  "
      elseif ( abs( coorkp(3)-s1u  )<epsil) then;kstr(11:15)=" 1+u "
      elseif ( abs( coorkp(3)-s1_v )<epsil) then;kstr(11:15)=" 1-v "
      else
       write(kstr(11:15),501) coorkp(3)
      endif

 501  FORMAT(F5.2)
      END SUBROUTINE Kreal2string

end module bilbao
